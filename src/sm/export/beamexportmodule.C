/*
*
*                 #####    #####   ######  ######  ###   ###
*               ##   ##  ##   ##  ##      ##      ## ### ##
*              ##   ##  ##   ##  ####    ####    ##  #  ##
*             ##   ##  ##   ##  ##      ##      ##     ##
*            ##   ##  ##   ##  ##      ##      ##     ##
*            #####    #####   ##      ######  ##     ##
*
*
*             OOFEM : Object Oriented Finite Element Code
*
*               Copyright (C) 1993 - 2013   Borek Patzak
*
*
*
*       Czech Technical University, Faculty of Civil Engineering,
*   Department of Structural Mechanics, 166 29 Prague, Czech Republic
*
*  This library is free software; you can redistribute it and/or
*  modify it under the terms of the GNU Lesser General Public
*  License as published by the Free Software Foundation; either
*  version 2.1 of the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*  Lesser General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General Public
*  License along with this library; if not, write to the Free Software
*  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "beamexportmodule.h"
#include "timestep.h"
#include "element.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/Beams/beambaseelement.h"
#include "../sm/Elements/Beams/beam3d.h"
#include "../sm/Elements/Beams/beam2d.h"
#include "gausspoint.h"
#include "engngm.h"
#include "material.h"
#include "classfactory.h"
#include "generalboundarycondition.h"
#include "constantedgeload.h"
#include "fei3dlinelin.h"
#include "inputrecord.h"
#include "../sm/EngineeringModels/responsespectrum.h"
#include <math.h>

using namespace std;

namespace oofem {
	REGISTER_ExportModule(BeamExportModule)

		BeamExportModule::BeamExportModule(int n, EngngModel *e) : ExportModule(n, e) { }

	BeamExportModule :: ~BeamExportModule() { }

	IRResultType
		BeamExportModule::initializeFrom(InputRecord *ir)
	{
		IRResultType result;                 // Required by IR_GIVE_FIELD macro
		isRespSpec = ir->hasField(_IFT_BeamExportModule_isrespspec);

		if (isRespSpec) {
			const char* name = this->emodel->giveClassName();
			if (!strcmp(name, "ResponseSpectrum") == 0) OOFEM_ERROR("Using rspec mode without a ResponseSpectrum engineering model");
			rs = dynamic_cast<ResponseSpectrum*>(this->emodel);
			if (!rs) OOFEM_ERROR("Error retrieving engmodel.");
		}

		return ExportModule::initializeFrom(ir);
	}

	void
		addComponents(FloatArray &dst, FloatArray &src, double pos, double len, bool momentOnly = false)
	{
		if (!momentOnly) {
			// axial force. linear interpolation
			dst.at(1) += src.at(1) * (len / 2 - pos);

			// shear forces. linear interpolation
			dst.at(3) += src.at(3) * (len / 2 - pos);
			dst.at(2) += src.at(2) * (len / 2 - pos);
		}

		// only these two parabolic. linear interpolation should be enough for the other directions
		dst.at(5) += src.at(3) * (pos) / 2 * (len - pos);
		dst.at(6) -= src.at(2) * (pos) / 2 * (len - pos);
	}

	//int checkValidType(const char* name)
	//{
	//	return (strcmp(name, "Beam3d") == 0) || (strcmp(name, "Beam2d") == 0) || (strcmp(name, "beam3d") == 0) || (strcmp(name, "beam2d") == 0);
	//}

	void
		BeamExportModule::doOutput(TimeStep *tStep, bool forcedOutput)
	{
		if (!(testTimeStepOutput(tStep) || forcedOutput)) {
			return;
		}

		vector< int >beamIDs;
		map<int, FloatArray >BeamLoads;
		IntArray temp;
		// loop through the beam elements
		Domain *d = emodel->giveDomain(1);
		if ((!this->isRespSpec && tStep->giveIntrinsicTime() == 0) || (tStep->giveIntrinsicTime() != 0)) {

			for (auto &elem : d->giveElements()) {
				if (this->checkValidType(elem->giveClassName())) {   // check if elem is beam (LIbeam?)

					int elNum;
					elNum = elem->giveNumber();
					//elNum = elem->giveLabel();

					// store IDs of known beams
					beamIDs.push_back(elNum);
					//beamIDs.push_back(elNum);

					BeamBaseElement *SElem;

					SElem = static_cast<BeamBaseElement *>(elem.get());

					double ksi, l = elem->computeLength();
					FloatArray Fl, loadEndForces;

					// get end forces considering the reduction due to winkler soil presence
					if (strcmp(SElem->giveClassName(), "Beam3d") == 0 || strcmp(SElem->giveClassName(), "beam3d") == 0){
						Beam3d *B3d = static_cast<Beam3d *>(SElem);
						B3d->giveEndForcesVector(Fl, tStep);
					}
					else {
						Beam2d *B2d = static_cast<Beam2d *>(SElem);
						B2d->giveEndForcesVector(Fl, tStep);
					}

					map< double, FloatArray >ForceDict;
					FloatArray I, E, Diff, dI, dE;

					I.resize(6);

					temp.resize(6);
					temp.at(1) = 1;
					temp.at(2) = 2;
					temp.at(3) = 3;
					temp.at(4) = 4;
					temp.at(5) = 5;
					temp.at(6) = 6;
					I.beSubArrayOf(Fl, temp);
					I *= -1; // invert sign of I end

					E.resize(6);
					for (int i = 1; i < 7; i++) {
						temp.at(i) = temp.at(i) + 6;
					}
					E.beSubArrayOf(Fl, temp);

					Diff.beDifferenceOf(E, I);

					FloatArray FinalLoads;
					FinalLoads.resize(6);
					FinalLoads.zero();

					FloatMatrix T;
					elem->computeGtoLRotationMatrix(T);
					T.resizeWithData(6, 6);

					temp.resize(6);
					temp.at(1) = 1;
					temp.at(2) = 2;
					temp.at(3) = 3;
					temp.at(4) = 4;
					temp.at(5) = 5;
					temp.at(6) = 6;

					IntArray *loads = elem->giveBoundaryLoadArray();
					for (auto &loadNum : *loads)
					{
						GeneralBoundaryCondition *bc = d->giveBc(loadNum);
						if (strcmp(bc->giveClassName(), "ConstantEdgeLoad") == 0) {
							ConstantEdgeLoad *CLoad = static_cast<ConstantEdgeLoad *>(bc);
							FloatArray compArr;
							const FloatArray coords;

							// CLoad->computeValues(compArr, tStep, NULL, temp, VM_Total);
							CLoad->computeValues(compArr, tStep, coords, temp, VM_Total);

							// transform to local coordinates
							compArr.rotatedWith(T, 'n');
							FinalLoads.add(compArr);
						}
					}

					addComponents(I, FinalLoads, 0.0, l, true);
					ForceDict[0.0] = I;

					for (GaussPoint *gp : *elem->giveDefaultIntegrationRulePtr()) {
						//double dV = elem->computeVolumeAround(gp);
						FloatArray ipState;
						double pos;

						ksi = 0.5 + 0.5 * gp->giveNaturalCoordinate(1);
						pos = ksi*l;

						// can't use this until beam is fixed?
						// elem->giveGlobalIPValue(ipState, gp, (InternalStateType)1, tStep); // IST_StressTensor
						ipState.zero();
						ipState.beScaled(ksi, Diff);
						ipState.add(I);

						addComponents(ipState, FinalLoads, pos, l, true);

						ForceDict[pos] = ipState;
					}

					addComponents(E, FinalLoads, l, l, true);
					ForceDict[l] = E;

					BeamForces[elem->giveNumber()] = ForceDict;
					//BeamForces[elem->giveLabel()] = ForceDict;

					//pair <double, double> loadPair;
					//loadPair.first = FinalLoads.at(2);
					//loadPair.second = FinalLoads.at(3);

					// save loads
					BeamLoads[elem->giveNumber()] = FinalLoads;

					//elem->giveBodyLoadArray
				}
			}

			//vector< unique_ptr< GeneralBoundaryCondition > > BCs = d->giveBcs();

			// tamper with stuff only if sets are defined.
			if (d->giveNumberOfSets()) {
				temp.resize(6);
				temp.at(1) = 1;
				temp.at(2) = 2;
				temp.at(3) = 3;
				temp.at(4) = 4;
				temp.at(5) = 5;
				temp.at(6) = 6;
				// loop through the loads
				for (auto &bc : d->giveBcs()) {
					// int bType = bc->giveBCValType(); // UNUSED: ConstantEdgeLoad is never == 2, they're all == 0 == unknown
					//if (bc->giveBCValType() == ForceLoadBVT) {
					if (strcmp(bc->giveClassName(), "ConstantEdgeLoad") == 0) {
						ConstantEdgeLoad *CLoad = static_cast<ConstantEdgeLoad *>(bc.get());

						// is it in a set?
						int nSet = CLoad->giveSetNumber();
						if (nSet) {
							Set *mySet = d->giveSet(nSet);
							// contains any of our beams?
							const IntArray &EdgeList = mySet->giveEdgeList();
							int numEdges = EdgeList.giveSize();

							int c = 1;
							while (c <= numEdges) {
								FloatArray compArr, tempArr;
								FloatMatrix T;
								int elNum = EdgeList.at(c);
								c += 2; // increment counter in interleaved array

								Element* ele = d->giveElement(elNum);
								if (!this->checkValidType(ele->giveClassName())) continue;

								const FloatArray coords;

								// CLoad->computeValues(compArr, tStep, NULL, temp, VM_Total);
								CLoad->computeValues(compArr, tStep, coords, temp, VM_Total);
								//d->giveElement(elNum)->computeBoundaryEdgeLoadVector(compArr, CLoad, edgeNum, ExternalForcesVector, VM_Total, tStep); // always vm_total???

								// transform to local coordinates
								d->giveElement(elNum)->computeGtoLRotationMatrix(T);
								T.resizeWithData(6, 6);
								compArr.rotatedWith(T, 'n');

								// add loads to our map
								BeamLoads[elNum] += compArr;
								//BeamLoads[elNum].first += compArr.at(2);
								//BeamLoads[elNum].second += compArr.at(3);

								// compute contribution to internal forces
								map< double, FloatArray >Dst = BeamForces[elNum];
								for (auto &PointVals : Dst) {
									const double &pos = PointVals.first;
									FloatArray Vals = PointVals.second;

									// tamper with values?
									addComponents(Vals, compArr, pos, d->giveElement(elNum)->computeLength());

									// update in point-forces map
									PointVals.second = Vals;
								}
								// update in beam-forces map
								BeamForces[elNum] = Dst;
							}
						}
					}
					//}
				}
			}

			for (auto beamPair : BeamForces)
			{
				int elNum = beamPair.first;
				Element *elem = d->giveElement(elNum);
				FloatArray rl, dI, dE; // used to store element end displacements
				FloatArray dNI, dNE; // used to store nodal displacements - may be different from the previous because of releases.
				FloatArray ddN; // used to store the difference between the ends.
				map< double, FloatArray >DispDict;
				double l = elem->computeLength();
				double l_2 = l*l;
				double l_3 = l_2*l;
				double l_4 = l_2*l_2;
				double ksi;
				//FloatMatrix shapeFunctions(2, 12);
				bool calc = false;
				//double phiy, phiz;

				elem->computeVectorOf(VM_Total, tStep, rl);
				temp.resize(6);
				temp.at(1) = 1;
				temp.at(2) = 2;
				temp.at(3) = 3;
				temp.at(4) = 4;
				temp.at(5) = 5;
				temp.at(6) = 6;

				FloatMatrix T;
				elem->computeGtoLRotationMatrix(T);
				T.resizeWithData(6, 6);

				dI.beSubArrayOf(rl, temp);
				//dI.rotatedWith(T, 'n');  // no need?

				DofManager *dofMan = elem->giveDofManager(1);
				dofMan->giveCompleteUnknownVector(dNI, VM_Total, tStep);
				FloatMatrix N;
				if (dofMan->computeL2GTransformation(N, NULL))
				{
					dNI.rotatedWith(N, 'n'); // rotate to global c.s.
				}
				dNI.rotatedWith(T, 'n');	// rotate to element c.s.

				dofMan = elem->giveDofManager(2);
				dofMan->giveCompleteUnknownVector(dNE, VM_Total, tStep);
				if (dofMan->computeL2GTransformation(N, NULL))
				{
					dNE.rotatedWith(N, 'n'); // rotate to global c.s.
				}
				dNE.rotatedWith(T, 'n');	// rotate to element c.s.

				ddN = dNE - dNI;

				// increment id array
				for (int i = 1; i <= 6; i++) temp.at(i) += 6;
				dE.beSubArrayOf(rl, temp);
				//dE.rotatedWith(T, 'n');

				//ddN = dE - dI;

				DispDict[0.0] = dI; // -dNI;

				CrossSection *Sect = elem->giveCrossSection();
				StructuralCrossSection *SCSect = static_cast<StructuralCrossSection *>(Sect);
				FloatMatrix MatStiffness;

				double EJyy, EJzz, EA, GJ, GKyAy, GKzAz;
				double psi_y, psi_z;
				double ay, by, cy, dy, fy;
				double anx, bnx, cnx;
				double az, bz, cz, dz, fz;
				double atx, btx, ctx;

				for (GaussPoint *gp : *elem->giveDefaultIntegrationRulePtr()) {
					FloatArray ipState;
					double pos, pos_2, pos_3, pos_4;

					ksi = 0.5 + 0.5 * gp->giveNaturalCoordinate(1);
					pos = ksi*l;

					pos_2 = pos*pos;
					pos_3 = pos_2*pos;
					pos_4 = pos_2*pos_2;

					// calculate this stuff on the first pass. Section are supposed prismatic (constant section)
					if (!calc){
						SCSect->give3dBeamStiffMtrx(MatStiffness, ElasticStiffness, gp, tStep);

						EA = MatStiffness.at(1, 1);
						GJ = MatStiffness.at(4, 4);
						EJzz = MatStiffness.at(6, 6);
						EJyy = MatStiffness.at(5, 5);
						GKyAy = MatStiffness.at(2, 2);
						GKzAz = MatStiffness.at(3, 3);

						psi_y = EJzz / GKyAy;
						psi_z = EJyy / GKzAz;

						double vy_0, vy_l, vz_0, vz_l;			// transversal displacements
						double phiy_0, phiy_l, phiz_0, phiz_l;	// rotation / first derivatives
						//double By, Ay, Bz, Az;
						double dx_0, dx_l;						// axial displacements
						double tx_0, tx_l;						// torsional rotations

						//map<double, FloatArray> &td = BeamDisplacements[elNum];
						FloatArray &bl = BeamLoads[elNum];
						FloatArray *disps = &dI;
						vy_0 = disps->at(2);
						vz_0 = disps->at(3);
						phiy_0 = -disps->at(5); // inverted signs for angles about y. in this case phi is used as first derivative
						phiz_0 = disps->at(6);
						dx_0 = disps->at(1);
						tx_0 = disps->at(4);

						disps = &dE;
						vy_l = disps->at(2);
						vz_l = disps->at(3);
						phiy_l = -disps->at(5); // inverted signs for angles about y. in this case phi is used as first derivative
						phiz_l = disps->at(6);
						dx_l = disps->at(1);
						tx_l = disps->at(4);


						// euler-bernoulli formulation
						//ey = vy_0;
						//dy = phiz_0;
						//ay = bl.at(2) / 24 / EJzz;

						//Ay = (vy_l - vy_0) / l_2 - ay*l_2 - dy / l;
						//By = (phiz_l - phiz_0) / l - ay*l_2 * 4;

						//by = (By - 2 * Ay) / l;
						//cy = 3 * Ay - By;


						//ez = vz_0;
						//dz = -phiy_0; // inverted signs for angles
						//az = bl.at(3) / 24 / EJyy;

						//Az = (vz_l - vz_0) / l_2 - az*l_2 - dz / l;
						//Bz = (-phiy_l + phiy_0) / l - az*l_2 * 4; // inverted signs for angles

						//bz = (Bz - 2 * Az) / l;
						//cz = 3 * Az - Bz;



						// timoshenko formulation for transversal displacements



						ay = bl.at(2) / 24 / EJzz;
						dy = phiz_0;
						fy = vy_0;
						by = (2 / l*(vy_0 - vy_l) + (phiz_l + phiz_0)) / (l_2 + 12 * psi_y) - 2 * ay*l;
						cy = -(72 * EJzz*GKyAy*l* (vy_0 - vy_l) + GKyAy*l_2*(24 * EJzz*(2 * phiz_0 + phiz_l) - l_3*bl.at(2)) + 12 * EJzz*(12 * EJzz*(phiz_0 - phiz_l) - l_3*bl.at(2))) / (24 * EJzz*l* (GKyAy*l_2 + 12 * EJzz));



						az = bl.at(3) / 24 / EJyy;
						dz = phiy_0;
						fz = vz_0;
						bz = (2 / l*(vz_0 - vz_l) + (phiy_l + phiy_0)) / (l_2 + 12 * psi_z) - 2 * az*l;
						cz = -(72 * EJyy*GKzAz*l* (vz_0 - vz_l) + GKzAz*l_2*(24 * EJyy*(2 * phiy_0 + phiy_l) - l_3*bl.at(3)) + 12 * EJyy*(12 * EJyy*(phiy_0 - phiy_l) - l_3*bl.at(3))) / (24 * EJyy*l* (GKzAz*l_2 + 12 * EJyy));

						// axial displacements
						cnx = dx_0;
						anx = -bl.at(1) / 2 / EA;
						bnx = (dx_l - dx_0) / l - anx*l;

						// torsional rotations
						ctx = tx_0;
						atx = -bl.at(4) / 2 / GJ;
						btx = (tx_l - tx_0) / l - atx*l;

						calc = true;
					}

					FloatArray disps(6);
					disps.at(1) = anx*pos_2 + bnx*pos + cnx;
					disps.at(2) = ay*pos_4 + by*pos_3 + (cy - 12 * ay*psi_y)*pos_2 + (dy - 6 * by*psi_y)*pos + fy;
					disps.at(3) = az*pos_4 + bz*pos_3 + (cz - 12 * az*psi_z)*pos_2 + (dz - 6 * bz*psi_z)*pos + fz;
					disps.at(4) = atx*pos_2 + btx*pos + ctx;
					disps.at(5) = -(4 * az*pos_3 + 3 * bz*pos_2 + 2 * cz*pos + dz);  // inverted signs for rotations about y.
					disps.at(6) = 4 * ay*pos_3 + 3 * by*pos_2 + 2 * cy*pos + dy;

					//disps -= (dI+ddN*ksi);

					DispDict[pos] = disps;

					//ipDisp.beProductOf(shapeFunctions, rl);
				}

				DispDict[l] = dE; // -dNE;

				// save the displacements
				BeamDisplacements[elem->giveNumber()] = DispDict;
				//BeamDisplacements[elem->giveLabel()] = DispDict;

			}
		}

		//for (auto &set : d->giveSets()) {
		//	IntArray &ElEdges = set->giveEdgeList();
		//}

		//	d->giveSets or d->giveLoad ?

		if (this->isRespSpec && tStep->giveIntrinsicTime()!=0){
			// square and save
			BeamDisplacementsList.push_back(BeamDisplacements);
			BeamForcesList.push_back(BeamForces);
			//BeamExportModule::addMultiply(combBeamDisplacements, BeamDisplacements,,1.0);
			//BeamExportModule::addMultiply(combBeamForces, BeamForces,,1.0);

		} else {

			if (this->isRespSpec && tStep->giveIntrinsicTime() == 0) {

				if (rs->giveComboType() == RSC_SRSS) {
					this->SRSS();
				} else {
					this->CQC();
				}

				BeamDisplacements = combBeamDisplacements;
				BeamForces = combBeamForces;

				combBeamDisplacements.clear();
				combBeamForces.clear();
			}

			double curTime = tStep->giveTargetTime();
			map<int, map<double, FloatArray>>::iterator BForces_it = BeamForces.begin();
			map<int, map<double, FloatArray>>::iterator BDisps_it = BeamDisplacements.begin();
			for (;
				BForces_it != BeamForces.end();
				++BForces_it, ++BDisps_it)
			{
				map< double, FloatArray > &BForces = BForces_it->second;
				map< double, FloatArray > &BDisps = BDisps_it->second;
				Element* elem = d->giveElement(BForces_it->first);
				int ID = elem->giveLabel();

				map<double, FloatArray >::iterator forces_it = BForces.begin();
				map<double, FloatArray >::iterator disps_it = BDisps.begin();

				for (;
					forces_it != BForces.end();
					++forces_it, ++disps_it)
				{
					double pos = forces_it->first;
					FloatArray forces = forces_it->second;
					FloatArray disps = disps_it->second;
					fprintf(this->stream, "%10.3e;%d;%10.3e;", curTime, ID, pos);

					for (auto &val : forces) {
						fprintf(this->stream, "%10.3e;", val);
					}
					for (auto &val : disps) {
						fprintf(this->stream, "%10.3e;", val);
					}
					fprintf(this->stream, "\n");
				}

			}

			//for (auto &bForces : BeamForces) {
			//	map< double, FloatArray >pForces = bForces.second;
			//	int ID = bForces.first;
			//	for (auto &vals : pForces) {
			//		double pos = vals.first;
			//		FloatArray forces = vals.second;
			//		fprintf(this->stream, "%10.3e;%d;%10.3e;", curTime, ID, pos);
			//		for (auto &val : forces) {
			//			fprintf(this->stream, "%10.3e;", val);
			//		}
			//		fprintf(this->stream, "\n");
			//	}
			//}

			// write file in the format:
			// elementNumber distanceFromIend N_x T_z T_y M_x M_y M_z
			// if 3 Gauss points are used, there would be 5 lines per beam (at distances 0, 0.1127*L, 0.5*L, 0.8873*L, L), ->>> to check

			//fprintf(this->stream, "%d ", avgState.giveSize());
			//for ( auto s: avgState ) {
			//    fprintf(this->stream, "%e ", s);
			//}
			//fprintf(this->stream, "    ");

			BeamDisplacements.clear();
			BeamForces.clear();

			fflush(this->stream);
		}
	}

	void BeamExportModule::populateElResults(map<int, map<double, FloatArray>> &answer, map<int, map<double, FloatArray>> &src)
	{

		map<int, map<double, FloatArray>>::iterator srcElem_it = src.begin();
		for (; srcElem_it != src.end(); ++srcElem_it)
		{
			map<double, FloatArray> *destBRespMap = new map<double, FloatArray>;
			map<double, FloatArray > &srcBRespMap = srcElem_it->second;

			map<double, FloatArray>::iterator srcBRespMap_it = srcBRespMap.begin();
			for (; srcBRespMap_it != srcBRespMap.end(); ++srcBRespMap_it)
			{
				FloatArray &srcRespArray = srcBRespMap_it->second;
				FloatArray *destRespArray = new FloatArray(srcRespArray.giveSize());

				destBRespMap->operator[](srcBRespMap_it->first) = *destRespArray;
			}

			answer[srcElem_it->first] = *destBRespMap;
		}
	}

	void BeamExportModule::addMultiply(map<int, map<double, FloatArray>> &answer, map<int, map<double, FloatArray>> &src, map<int, map<double, FloatArray>> &src2, double fact)
	{
		if (answer.size() == 0) {
			populateElResults(answer, src);
		}

		// awful iteration
		map<int, map<double, FloatArray>>::iterator destElem_it = answer.begin();
		map<int, map<double, FloatArray>>::iterator srcElem_it = src.begin();
		map<int, map<double, FloatArray>>::iterator srcElem_it2 = src2.begin();
		for (; destElem_it != answer.end(); ++destElem_it, ++srcElem_it, ++srcElem_it2)
		{
			map<double, FloatArray> &destRespMap = destElem_it->second;
			map<double, FloatArray> &srcRespMap = srcElem_it->second;
			map<double, FloatArray> &srcRespMap2 = srcElem_it2->second;

			map<double, FloatArray>::iterator destRespMap_it = destRespMap.begin();
			map<double, FloatArray>::iterator srcRespMap_it = srcRespMap.begin();
			map<double, FloatArray>::iterator srcRespMap_it2 = srcRespMap2.begin();
			for (; destRespMap_it != destRespMap.end(); ++destRespMap_it, ++srcRespMap_it, ++srcRespMap_it2)
			{
				FloatArray &destRespArray = destRespMap_it->second;
				FloatArray &srcRespArray = srcRespMap_it->second;
				FloatArray &srcRespArray2 = srcRespMap_it2->second;

				for (int i = 1; i <= srcRespArray.giveSize(); i++)
				{
					// square it and add it
					destRespArray.at(i) += srcRespArray.at(i)*srcRespArray2.at(i)*fact;
				}
			}
		}
	}

	void BeamExportModule::calcRoot(map<int, map<double, FloatArray>> &answer)
	{
		// another awful iteration
		map<int, map<double, FloatArray>>::iterator destElem_it = answer.begin();
		for (; destElem_it != answer.end(); ++destElem_it)
		{
			map<double, FloatArray> &destRespMap = destElem_it->second;

			map<double, FloatArray>::iterator destRespMap_it = destRespMap.begin();
			for (; destRespMap_it != destRespMap.end(); ++destRespMap_it)
			{
				FloatArray &destRespArray = destRespMap_it->second;

				for (int i = 1; i <= destRespArray.giveSize(); i++)
				{
					// square it and add it
					destRespArray.at(i) = sqrt(destRespArray.at(i));
				}
			}
		}
	}

	void BeamExportModule::SRSS(){
		list<map<int, map<double, FloatArray>>>::iterator disps_it = BeamDisplacementsList.begin();
		for (; disps_it != BeamDisplacementsList.end(); ++disps_it)
		{
			addMultiply(combBeamDisplacements, *disps_it, *disps_it);
		}
		calcRoot(combBeamDisplacements);

		list<map<int, map<double, FloatArray>>>::iterator forces_it = BeamForcesList.begin();
		for (; forces_it != BeamForcesList.end(); ++forces_it)
		{
			addMultiply(combBeamForces, *forces_it, *forces_it);  // mult by 1.0
		}
		calcRoot(combBeamForces);

	}

	void BeamExportModule::CQC(){
		FloatMatrix rhos;
		rs->giveRhos(rhos);

		list<map<int, map<double, FloatArray>>>::iterator disps_it = BeamDisplacementsList.begin();

		for (int i = 1; disps_it != BeamDisplacementsList.end(); ++disps_it, i++)
		{
			list<map<int, map<double, FloatArray>>>::iterator disps_it2 = BeamDisplacementsList.begin();
			for (int j = 1; disps_it2 != BeamDisplacementsList.end(); ++disps_it2, j++)
			{
				addMultiply(combBeamDisplacements, *disps_it, *disps_it2, rhos.at(i, j));
			}
		}
		calcRoot(combBeamDisplacements);

		list<map<int, map<double, FloatArray>>>::iterator forces_it = BeamForcesList.begin();

		for (int i = 1; forces_it != BeamForcesList.end(); ++forces_it, i++)
		{
			list<map<int, map<double, FloatArray>>>::iterator forces_it2 = BeamForcesList.begin();
			for (int j = 1; forces_it2 != BeamForcesList.end(); ++forces_it2, j++)
			{
				addMultiply(combBeamForces, *forces_it, *forces_it2, rhos.at(i, j));
			}
		}
		calcRoot(combBeamForces);
	}

	void
		BeamExportModule::initialize()
	{
		string fileName = emodel->giveOutputBaseFileName() + ".bem";
		if ((this->stream = fopen(fileName.c_str(), "w")) == NULL) {
			OOFEM_ERROR("failed to open file %s", fileName.c_str());
		}
		// ";" as separator
		fprintf(this->stream, "#Time;BeamNo;DistanceFromI;N_x;T_y;T_z;M_x;M_y;M_z;dx;dy;dz;rx;ry;rz;");
		//for ( int var: this->ists ) {
		//    fprintf(this->stream, "%s    ", __InternalStateTypeToString( ( InternalStateType ) var) );
		//}
		fprintf(this->stream, "\n");
		fflush(this->stream);
	}

	void
		BeamExportModule::terminate()
	{
		fclose(this->stream);
	}
} // end namespace oofem
