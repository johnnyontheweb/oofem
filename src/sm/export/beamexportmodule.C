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
#include "linearedgeload.h"
#include "fei3dlinelin.h"
#include "inputrecord.h"
#include "../sm/EngineeringModels/responsespectrum.h"
#include <math.h>
#include <functional>

#ifdef MEMSTR
#include <io.h>
#include <fcntl.h>
#endif

using namespace std;

namespace oofem {
	REGISTER_ExportModule(BeamExportModule)

		BeamExportModule::BeamExportModule(int n, EngngModel *e) : ExportModule(n, e) { }

	BeamExportModule :: ~BeamExportModule() { }

	IRResultType
		BeamExportModule::initializeFrom(InputRecord *ir)
	{
		// IRResultType result;                 // Required by IR_GIVE_FIELD macro
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
		addComponents(FloatArray &dst, std::pair<FloatArray,FloatArray> &src, double pos, double len, bool momentOnly = false)
	{
		FloatArray &qi = src.first, &qf = src.second;
		double tot1, tot2, tot3;
		tot1 = (qi.at(1) + qf.at(1))*len / 2;
		tot2 = (qi.at(2) + qf.at(2))*len / 2;
		tot3 = (qi.at(3) + qf.at(3))*len / 2;
		double center1 = len / 2, center2 = len / 2, center3 = len / 2;
		if (tot1 != 0.0) center1 = (qi.at(1) + 2 * qf.at(1))*len / 3 / (qi.at(1) + qf.at(1));
		if (tot2 != 0.0) center2 = (qi.at(2) + 2 * qf.at(2))*len / 3 / (qi.at(2) + qf.at(2));
		if (tot3 != 0.0) center3 = (qi.at(3) + 2 * qf.at(3))*len / 3 / (qi.at(3) + qf.at(3));
		double V1I, V2I, V3I;
		V1I = tot1*(len-center1) / len;
		V2I = tot2*(len-center2) / len;
		V3I = tot3*(len-center3) / len;
		double V1E, V2E, V3E;
		V1E = V1I-tot1;
		V2E = V2I-tot2;
		V3E = V3I-tot3;
		double pos2 = pos*pos, pos3 = pos2*pos;

		if (!momentOnly) {
			// axial force.
			dst.at(1) -= (qi.at(1)*pos + (qf.at(1) - qi.at(1)) / 2 / len*pos*pos);

			// shear forces.
			dst.at(2) -= (qi.at(2)*pos + (qf.at(2) - qi.at(2)) / 2 / len*pos*pos);
			dst.at(3) -= (qi.at(3)*pos + (qf.at(3) - qi.at(3)) / 2 / len*pos*pos);
		}

		// moments.
		dst.at(5) += (V3E + V3I) / len / len*pos3 - (2 * V3I + V3E) / len*pos2 + V3I*pos;
		dst.at(6) -= (V2E + V2I) / len / len*pos3 - (2 * V2I + V2E) / len*pos2 + V2I*pos;
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
		map<int, std::pair<FloatArray,FloatArray> >BeamLoads;
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
					FloatArray Fl, Fl2, loadEndForces;

					// get end forces considering the reduction due to winkler soil presence
					if (strcmp(SElem->giveClassName(), "Beam3d") == 0 || strcmp(SElem->giveClassName(), "beam3d") == 0){
						Beam3d *B3d = static_cast<Beam3d *>(SElem);
						B3d->giveEndForcesVector(Fl, tStep, true);
					}
					else {
						Beam2d *B2d = static_cast<Beam2d *>(SElem);
						B2d->giveEndForcesVector(Fl, tStep, true);
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

					std::pair<FloatArray, FloatArray> FinalLoads;
					FinalLoads.first.resize(6);
					FinalLoads.first.zero();
					FinalLoads.second.resize(6);
					FinalLoads.second.zero();

					//FloatArray FinalLoads;
					//FinalLoads.resize(6);
					//FinalLoads.zero();

					//FloatArray FinalLoadsL;
					//FinalLoadsL.resize(12); // 6 by 2 position (start & end)
					//FinalLoadsL.zero();

					// temporary stuff for winker
					//double wy, wz;
					//wy = wz = 0.0;
					//double qIy, qIz, qEy, qEz;
					//qIy = qIz = qEy = qEz = 0.0;

					// get end forces without the reduction due to winkler soil presence
					//if (strcmp(SElem->giveClassName(), "Beam3d") == 0 || strcmp(SElem->giveClassName(), "beam3d") == 0){
					//	Beam3d *B3d = static_cast<Beam3d *>(SElem);
					//	B3d->giveEndForcesVector(Fl2, tStep, false);
					//	FloatArray compArr, diff;
					//	diff.beDifferenceOf(Fl, Fl2);
					//	compArr.resize(6);

					//	wy = compArr.at(2) = -(diff.at(2) + diff.at(8)) / l;
					//	wz = compArr.at(3) = -(diff.at(3) + diff.at(9)) / l;
					//	FinalLoads.add(compArr);

					//	qIy = -(3 * diff.at(2) - diff.at(8)) / l;
					//	qIz = -(3 * diff.at(3) - diff.at(9)) / l;
					//	qEy = -(3 * diff.at(8) - diff.at(2)) / l;
					//	qEz = -(3 * diff.at(9) - diff.at(3)) / l;
					//}
					//else {
					//	Beam2d *B2d = static_cast<Beam2d *>(SElem);
					//	B2d->giveEndForcesVector(Fl2, tStep, false);
					//	FloatArray compArr, diff;
					//	diff.beDifferenceOf(Fl, Fl2);
					//	compArr.resize(6);

					//	wy = compArr.at(2) = -(diff.at(2) + diff.at(8)) / l;  // @TODO: check array positions for 2d beams
					//	wz = compArr.at(3) = -(diff.at(3) + diff.at(9)) / l;
					//	FinalLoads.add(compArr);

					//	qIy = -(3 * diff.at(2) - diff.at(8)) / l;
					//	qIz = -(3 * diff.at(3) - diff.at(9)) / l;
					//	qEy = -(3 * diff.at(8) - diff.at(2)) / l;
					//	qEz = -(3 * diff.at(9) - diff.at(3)) / l;
					//}

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

					std::list<double> midPoints;
					//midPoints.push_back(0.125);
					midPoints.push_back(0.25);
					//midPoints.push_back(0.375);
					midPoints.push_back(0.5);
					//midPoints.push_back(0.625);
					midPoints.push_back(0.75);
					//midPoints.push_back(0.875);

					IntArray *loads = elem->giveBoundaryLoadArray();
					int nBoundaryLoads = loads->giveSize() / 2;
					for (int i = 1; i <= nBoundaryLoads; i++)
					{
						int loadNum = loads->at(1 + (i - 1) * 2);
						GeneralBoundaryCondition *bc = d->giveBc(loadNum);
						if (strcmp(bc->giveClassName(), "ConstantEdgeLoad") == 0) {
							ConstantEdgeLoad *CLoad = static_cast<ConstantEdgeLoad *>(bc);
							FloatArray compArr;
							FloatArray coords;

							// CLoad->computeValues(compArr, tStep, NULL, temp, VM_Total);
							CLoad->computeValues(compArr, tStep, coords, temp, VM_Total);

							// transform to local coordinates
							if (CLoad->giveCoordSystMode() == Load::CoordSystType::CST_Global)	compArr.rotatedWith(T, 'n');
							// they're the same
							FinalLoads.first.add(compArr);
							FinalLoads.second.add(compArr);
						}
						else if (strcmp(bc->giveClassName(), "LinearEdgeLoad") == 0) {
							LinearEdgeLoad *CLoad = static_cast<LinearEdgeLoad *>(bc);
							FloatArray compArr;
							FloatArray coords(3);

							coords.at(1) = -1;  // 1st node param. coord

							// CLoad->computeValues(compArr, tStep, NULL, temp, VM_Total);
							CLoad->computeValues(compArr, tStep, coords, temp, VM_Total);

							// transform to local coordinates
							if (CLoad->giveCoordSystMode() == Load::CoordSystType::CST_Global)	compArr.rotatedWith(T, 'n');
							FinalLoads.first.add(compArr);
							
							coords.at(1) = 1;  // 2nd node param. coord

							// CLoad->computeValues(compArr, tStep, NULL, temp, VM_Total);
							CLoad->computeValues(compArr, tStep, coords, temp, VM_Total);

							// transform to local coordinates
							if (CLoad->giveCoordSystMode() == Load::CoordSystType::CST_Global)	compArr.rotatedWith(T, 'n');
							FinalLoads.second.add(compArr);

							//// add stations inside load
							//midPoints.push_back(CLoad->startLocal);
							//double halfL = CLoad->endLocal - CLoad->startLocal;
							//midPoints.push_back(CLoad->startLocal + 0.25 * halfL);
							//midPoints.push_back(CLoad->startLocal + 0.5 * halfL);
							//midPoints.push_back(CLoad->startLocal + 0.75 * halfL);
							//midPoints.push_back(CLoad->endLocal);
						}
					}

					addComponents(I, FinalLoads, 0.0, l, true);
					ForceDict[0.0] = I;

					// temporary stuff for winkler
					//FloatArray wI, wE;
					//wI.resize(2); wE.resize(2);
					//wI.at(1) = qIy; wI.at(2) = qIz;
					//wE.at(1) = qEy; wE.at(2) = qEz;
					//winkDict[0.0] = wI;


					//double tempmidpoint = 0;
					//int count = 0;

					//for (GaussPoint *gp : *elem->giveDefaultIntegrationRulePtr()) {
					//	//double dV = elem->computeVolumeAround(gp);
					//	FloatArray ipState;
					//	FloatArray winkState;
					//	double pos;

					//	ksi = 0.5 + 0.5 * gp->giveNaturalCoordinate(1);
					//	pos = ksi*l;
					//	//if (count) {
					//	//	tempmidpoint += ksi / 2;
					//	//	midPoints.push_back(tempmidpoint);
					//	//}

					//	// can't use this until beam is fixed?
					//	// elem->giveGlobalIPValue(ipState, gp, (InternalStateType)1, tStep); // IST_StressTensor
					//	ipState.zero();
					//	ipState.beScaled(ksi, Diff);
					//	ipState.add(I);

					//	addComponents(ipState, FinalLoads, pos, l, true);

					//	ForceDict[pos] = ipState;

					//	//winkState.resize(2);
					//	//winkState.at(1) = qIy * (l - pos) / l + qEy * pos / l;
					//	//winkState.at(2) = qIz * (l - pos) / l + qEz * pos / l;
					//	//winkDict[pos] = winkState;
					//	//tempmidpoint = ksi / 2;
					//	//count += 1;
					//}

					for (double midP : midPoints) {
						FloatArray ipState;
						FloatArray winkState;
						double pos;

						pos = midP*l;

						// can't use this until beam is fixed?
						// elem->giveGlobalIPValue(ipState, gp, (InternalStateType)1, tStep); // IST_StressTensor
						ipState.resize(6); ipState.zero();
						// ipState.beScaled(midP, Diff);
						// moments only
						ipState.at(5) = Diff.at(5) * midP;
						ipState.at(6) = Diff.at(6) * midP;
						ipState.add(I);

						addComponents(ipState, FinalLoads, pos, l, false);

						ForceDict[pos] = ipState;
					}

					addComponents(E, FinalLoads, l, l, true);
					ForceDict[l] = E;
					/*winkDict[l] = wE;*/

					BeamForces[elem->giveNumber()] = ForceDict;
					//BeamWinkler[elem->giveNumber()] = winkDict;
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

			// only if sets are defined. actually not used, wrong end forces in oofem for beam edge loads.
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

								FloatArray coords;

								// CLoad->computeValues(compArr, tStep, NULL, temp, VM_Total);
								CLoad->computeValues(compArr, tStep, coords, temp, VM_Total);
								//d->giveElement(elNum)->computeBoundaryEdgeLoadVector(compArr, CLoad, edgeNum, ExternalForcesVector, VM_Total, tStep); // always vm_total???

								// transform to local coordinates
								d->giveElement(elNum)->computeGtoLRotationMatrix(T);
								T.resizeWithData(6, 6);
								if (CLoad->giveCoordSystMode() == Load::CoordSystType::CST_Global)	compArr.rotatedWith(T, 'n');

								// add loads to our map
								BeamLoads[elNum].first += compArr;

								std::pair<FloatArray, FloatArray> localpair;
								localpair.first = compArr;

								// add loads to our map
								BeamLoads[elNum].second += compArr;
								localpair.second = compArr;

								double ll = d->giveElement(elNum)->computeLength();
								FloatArray Diff; Diff.beDifferenceOf(BeamForces[elNum][ll], BeamForces[elNum][0.0]);

								// compute contribution to internal forces
								map< double, FloatArray >Dst = BeamForces[elNum];
								for (auto &PointVals : Dst) {
									const double &pos = PointVals.first;
									FloatArray Vals = PointVals.second;

									// tamper with values?
									// addComponents(Vals, localpair, pos, d->giveElement(elNum)->computeLength());
									if (pos == 0.0 || pos == ll) {
										addComponents(Vals, localpair, pos, ll, true);
									} else {
										// moments only
										Vals.at(5) = Diff.at(5) * pos / ll;
										Vals.at(6) = Diff.at(6) * pos / ll;
										Vals.add(BeamForces[elNum][0.0]);
										addComponents(Vals, localpair, pos, ll, false);
									}
									
									// update in point-forces map
									PointVals.second = Vals;
								}
								// update in beam-forces map
								BeamForces[elNum] = Dst;
							}
						}
					} else if (strcmp(bc->giveClassName(), "LinearEdgeLoad") == 0) {
						LinearEdgeLoad *CLoad = static_cast<LinearEdgeLoad *>(bc.get());

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

								FloatArray coords(3);
								coords.at(1) = -1;

								// CLoad->computeValues(compArr, tStep, NULL, temp, VM_Total);
								CLoad->computeValues(compArr, tStep, coords, temp, VM_Total);
								//d->giveElement(elNum)->computeBoundaryEdgeLoadVector(compArr, CLoad, edgeNum, ExternalForcesVector, VM_Total, tStep); // always vm_total???

								// transform to local coordinates
								d->giveElement(elNum)->computeGtoLRotationMatrix(T);
								T.resizeWithData(6, 6);
								if (CLoad->giveCoordSystMode() == Load::CoordSystType::CST_Global)	compArr.rotatedWith(T, 'n');

								// add loads to our map
								BeamLoads[elNum].first += compArr;

								std::pair<FloatArray, FloatArray> localpair;
								localpair.first = compArr;

								coords.at(1) = 1;
								CLoad->computeValues(compArr, tStep, coords, temp, VM_Total);
								// transform to local coordinates
								d->giveElement(elNum)->computeGtoLRotationMatrix(T);
								T.resizeWithData(6, 6);
								if (CLoad->giveCoordSystMode() == Load::CoordSystType::CST_Global)	compArr.rotatedWith(T, 'n');

								// add loads to our map
								BeamLoads[elNum].second += compArr;
								localpair.second = compArr;

								double ll = d->giveElement(elNum)->computeLength();
								FloatArray Diff; Diff.beDifferenceOf(BeamForces[elNum][ll], BeamForces[elNum][0.0]);
								
								// compute contribution to internal forces
								map< double, FloatArray >Dst = BeamForces[elNum];
								for (auto &PointVals : Dst) {
									const double &pos = PointVals.first;
									FloatArray Vals = PointVals.second;

									// tamper with values?
									// addComponents(Vals, localpair, pos, d->giveElement(elNum)->computeLength());
									if (pos == 0.0 || pos == ll) {
										addComponents(Vals, localpair, pos, ll, true);
									}
									else {
										// moments only
										Vals.at(5) = Diff.at(5) * pos / ll;
										Vals.at(6) = Diff.at(6) * pos / ll;
										Vals.add(BeamForces[elNum][0.0]);
										addComponents(Vals, localpair, pos, ll, false);
									}

									// update in point-forces map
									PointVals.second = Vals;
								}
								// update in beam-forces map
								BeamForces[elNum] = Dst;
							}
						}
					}

				}
			}

			// in the next section all deflections are calculated.
			// beam on soil deflections and forces are calculated by directly solving the 4th order ODE for winkler formulation v(IV) + 4*lambda^4*v = q.
			// Boundary conditions considered are those relative to displacements and rotation (relative to v and v(I)).
			// Shears and moments are calculated using the closed form derivatives.
			// For beam on soil formulation, v(IV) - alpha*v(II) + lambda*v = q.

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
				// double l_4 = l_2*l_2;
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

				double EJyy = 0, EJzz = 0, EA = 0, GJ = 0, GKyAy = 0, GKzAz = 0;
				double psi_y = 0, psi_z = 0;
				double ay = 0, by = 0, cy = 0, dy = 0, ey=0, fy = 0;
				double anx = 0, bnx = 0, cnx = 0, dnx=0;
				double az = 0, bz = 0, cz = 0, dz = 0, ez=0, fz = 0;
				double atx = 0, btx = 0, ctx = 0, dtx=0;
				bool hasWinklerY, hasWinklerZ;
				hasWinklerY = hasWinklerZ = false;
				double wy = 0, wz = 0;  // winkler stiffness
				double lambdaY=0.0, lambdaZ=0.0;
				FloatArray W0(6), WL(6);
				// values for winkler+timoshenko
				double deltaY = 0, deltaZ = 0;
				double alphaY = 0, alphaZ = 0;
				double lambdaY1 = 0, lambdaY2 = 0;
				double lambdaZ1 = 0, lambdaZ2 = 0;

				// saving winkler reaction for each gp
				map< double, FloatArray > WinkDict;

				FloatArray *disps = &dI;

				// shorthands for the loads
				FloatArray &qi = BeamLoads[elNum].first, &qf = BeamLoads[elNum].second;
				
				std::map< double, FloatArray > &Forces = BeamForces[elem->giveNumber()];

				for (auto &force : Forces) {
				//for (GaussPoint *gp : *elem->giveDefaultIntegrationRulePtr()) {
					//FloatArray ipState;
					double pos = force.first;
					double pos_2, pos_3, pos_4, pos_5;

					//ksi = 0.5 + 0.5 * gp->giveNaturalCoordinate(1);
					//pos = ksi*l;

					pos_2 = pos*pos;
					pos_3 = pos_2*pos;
					pos_4 = pos_2*pos_2;
					pos_5 = pos_2*pos_3;

					// calculate this stuff on the first pass. Constant section along the length
					if (!calc){
						GaussPoint *gp = *(elem->giveDefaultIntegrationRulePtr()->begin());
						SCSect->give3dBeamStiffMtrx(MatStiffness, ElasticStiffness, gp, tStep);

						EA = MatStiffness.at(1, 1);
						GJ = MatStiffness.at(4, 4);
						EJzz = MatStiffness.at(6, 6);
						EJyy = MatStiffness.at(5, 5);
						GKyAy = MatStiffness.at(2, 2);
						GKzAz = MatStiffness.at(3, 3);

						if (GKyAy < 1e-6) {
							psi_y = 0;
						}
						else {
							psi_y = EJzz / GKyAy;
						}

						if (GKzAz < 1e-6) {
							psi_z = 0;
						}
						else {
							psi_z = EJyy / GKzAz;
						}

						double vy_0, vy_l, vz_0, vz_l;			// transversal displacements
						double phiy_0, phiy_l, phiz_0, phiz_l;	// rotation / first derivatives
						//double By, Ay, Bz, Az;
						double dx_0, dx_l;						// axial displacements
						double tx_0, tx_l;						// torsional rotations

						//map<double, FloatArray> &td = BeamDisplacements[elNum];

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

						Beam3d *B3d = static_cast<Beam3d *>(elem);
						Material* subSoilMatTemp = B3d->giveSubSoilMaterial();
						if (subSoilMatTemp != NULL)
						{
							StructuralMaterial* subSoilMat = (StructuralMaterial*)subSoilMatTemp;
							FloatMatrix subMat;
							subSoilMat->give3dBeamSubSoilStiffMtrx(subMat, MatResponseMode::TangentStiffness, gp, tStep);
							wy = subMat.at(2, 2);
							if (wy > 1e-6) hasWinklerY = true;

							W0.at(2) = -dI.at(2)*wy;
							WL.at(2) = -dE.at(2)*wy;


							wz = subMat.at(3, 3);
							if (wz > 1e-6) hasWinklerZ = true;

							W0.at(3) = -dI.at(3)*wz;
							WL.at(3) = -dE.at(3)*wz;
						}

						WinkDict[0.0] = W0;

						if (hasWinklerY) {
							FloatMatrix odeMtrx(4, 4);
							FloatArray rhs(4);
							FloatArray abcd(4);

							if (psi_y == 0.0) { // Euler Bernoulli formulation
								lambdaY = sqrt(sqrt(wy / 4 / EJzz));

								odeMtrx.at(1, 1) = 1;
								odeMtrx.at(1, 3) = 1;
								odeMtrx.at(2, 1) = cos(l*lambdaY) * exp(l*lambdaY);
								odeMtrx.at(2, 2) = sin(l*lambdaY) * exp(l*lambdaY);
								odeMtrx.at(2, 3) = cos(l*lambdaY) / exp(l*lambdaY);
								odeMtrx.at(2, 4) = sin(l*lambdaY) / exp(l*lambdaY);
								odeMtrx.at(3, 1) = lambdaY;
								odeMtrx.at(3, 2) = lambdaY;
								odeMtrx.at(3, 3) = -lambdaY;
								odeMtrx.at(3, 4) = lambdaY;
								odeMtrx.at(4, 1) = lambdaY * (cos(l*lambdaY) - sin(l*lambdaY)) * exp(l*lambdaY);
								odeMtrx.at(4, 2) = lambdaY * (cos(l*lambdaY) + sin(l*lambdaY)) * exp(l*lambdaY);
								odeMtrx.at(4, 3) = -lambdaY * (cos(l*lambdaY) + sin(l*lambdaY)) / exp(l*lambdaY);
								odeMtrx.at(4, 4) = lambdaY * (cos(l*lambdaY) - sin(l*lambdaY)) / exp(l*lambdaY);
							}
							else {  // Timoshenko formulation
								alphaY = wy / GKyAy;
								lambdaY = wy / EJzz;
								deltaY = alphaY*alphaY - 4 * lambdaY;

								if (deltaY > 0) {
									lambdaY1 = sqrt(alphaY / 2 + 0.5*sqrt(deltaY));
									lambdaY2 = sqrt(alphaY / 2 - 0.5*sqrt(deltaY));

									odeMtrx.at(1, 1) = 1;
									odeMtrx.at(1, 2) = 1;
									odeMtrx.at(1, 3) = 1;
									odeMtrx.at(1, 4) = 1;
									odeMtrx.at(2, 1) = exp(l*lambdaY1);
									odeMtrx.at(2, 2) = exp(-l*lambdaY1);
									odeMtrx.at(2, 3) = exp(l*lambdaY2);
									odeMtrx.at(2, 4) = exp(-l*lambdaY2);
									odeMtrx.at(3, 1) =  lambdaY1*(lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1);
									odeMtrx.at(3, 2) = -lambdaY1*(lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1);
									odeMtrx.at(3, 3) =  lambdaY2*(lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1);
									odeMtrx.at(3, 4) = -lambdaY2*(lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1);
									odeMtrx.at(4, 1) =  lambdaY1*(lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1)*exp(l*lambdaY1);
									odeMtrx.at(4, 2) = -lambdaY1*(lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1)/exp(l*lambdaY1);
									odeMtrx.at(4, 3) =  lambdaY2*(lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1)*exp(l*lambdaY2);
									odeMtrx.at(4, 4) = -lambdaY2*(lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1)/exp(l*lambdaY2);
								}
								else if (deltaY == 0) {
									lambdaY1 = sqrt(alphaY / 2);
									lambdaY2 = -sqrt(alphaY / 2);

									odeMtrx.at(1, 1) = 1;
									odeMtrx.at(1, 2) = 1;
									odeMtrx.at(2, 1) = exp(l*lambdaY1);
									odeMtrx.at(2, 2) = exp(-l*lambdaY1);
									odeMtrx.at(2, 3) = l*exp(l*lambdaY1);
									odeMtrx.at(2, 4) = l*exp(-l*lambdaY1);
									odeMtrx.at(3, 1) =  lambdaY1*(lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1);
									odeMtrx.at(3, 2) = -lambdaY1*(lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1);
									odeMtrx.at(3, 3) =  lambdaY1*(3*lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1);
									odeMtrx.at(3, 4) =  lambdaY1*(3*lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1);
									odeMtrx.at(4, 1) =  lambdaY1*(lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1)*exp(l*lambdaY1);
									odeMtrx.at(4, 2) = -lambdaY1*(lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1)/exp(l*lambdaY1);
									odeMtrx.at(4, 3) =  (l*lambdaY1*(lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1) + 3 * lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1)*exp(l*lambdaY1);
									odeMtrx.at(4, 4) =  (-l*lambdaY1*(lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1) + 3 * lambdaY1*lambdaY1*psi_y - alphaY*psi_y + 1)/exp(l*lambdaY1);
								}
								else {
									lambdaY1 = sqrt(alphaY / 4 + 0.5*sqrt(lambdaY));
									lambdaY2 = sqrt(-alphaY / 4 + 0.5*sqrt(lambdaY));

									odeMtrx.at(1, 1) = 1;
									odeMtrx.at(1, 2) = 1;
									odeMtrx.at(2, 1) = cos(l*lambdaY2) * exp(l*lambdaY1);
									odeMtrx.at(2, 2) = cos(l*lambdaY2) / exp(l*lambdaY1);
									odeMtrx.at(2, 3) = sin(l*lambdaY2) * exp(l*lambdaY1);
									odeMtrx.at(2, 4) = sin(l*lambdaY2) / exp(l*lambdaY1);
									odeMtrx.at(3, 1) =  lambdaY1*(lambdaY1*lambdaY1*psi_y - 3 * lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1);
									odeMtrx.at(3, 2) = -lambdaY1*(lambdaY1*lambdaY1*psi_y - 3 * lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1);
									odeMtrx.at(3, 3) =  lambdaY2*(3 * lambdaY1*lambdaY1*psi_y - lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1);
									odeMtrx.at(3, 4) =  lambdaY2*(3 * lambdaY1*lambdaY1*psi_y - lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1);
									odeMtrx.at(4, 1) =  ( lambdaY1*(lambdaY1*lambdaY1*psi_y - 3 * lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1)*cos(l*lambdaY2) - lambdaY2*(3 * lambdaY1*lambdaY1*psi_y - lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1)*sin(l*lambdaY2)) * exp(l*lambdaY1);
									odeMtrx.at(4, 2) = -( lambdaY1*(lambdaY1*lambdaY1*psi_y - 3 * lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1)*cos(l*lambdaY2) + lambdaY2*(3 * lambdaY1*lambdaY1*psi_y - lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1)*sin(l*lambdaY2)) / exp(l*lambdaY1);
									odeMtrx.at(4, 3) =  ( lambdaY2*(3 * lambdaY1*lambdaY1*psi_y - lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1)*cos(l*lambdaY2) + lambdaY1*(lambdaY1*lambdaY1*psi_y - 3 * lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1)*sin(l*lambdaY2)) * exp(l*lambdaY1);
									odeMtrx.at(4, 4) = -(-lambdaY2*(3 * lambdaY1*lambdaY1*psi_y - lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1)*cos(l*lambdaY2) + lambdaY1*(lambdaY1*lambdaY1*psi_y - 3 * lambdaY2*lambdaY2*psi_y - alphaY*psi_y + 1)*sin(l*lambdaY2)) / exp(l*lambdaY1);
								}
							}

							//rhs.at(1) = dI.at(2) - bl.at(2) / wy;
							//rhs.at(2) = dE.at(2) - bl.at(2) / wy;
							//rhs.at(3) = dI.at(6);
							//rhs.at(4) = dE.at(6);

							rhs.at(1) = dI.at(2) - qi.at(2) / wy;
							rhs.at(2) = dE.at(2) - qf.at(2) / wy;
							rhs.at(3) = dI.at(6) + (qf.at(2) - qi.at(2)) / (wy*l);
							rhs.at(4) = dE.at(6) + (qf.at(2) - qi.at(2)) / (wy*l);

							//odeMtrx.computeReciprocalCondition('1');
							odeMtrx.solveForRhs(rhs, abcd);

#ifdef DEBUG
							// odeMtrx.writeCSV("ODE.csv");
#endif

							ay = abcd.at(1);
							by = abcd.at(2);
							cy = abcd.at(3);
							dy = abcd.at(4);
						}
						else {
							if (psi_y == 0.0) {  // Euler Bernoulli formulation
								//ay = bl.at(2) / 24 / EJzz;
								//dy = phiz_0;
								//fy = vy_0;
								//by = (2 * (vy_0 - vy_l) + l*(phiz_l + phiz_0)) / (l_3)-2 * ay*l;
								//cy = -(3 * (vy_0 - vy_l) + l*(2 * phiz_0 + phiz_l)) / l_2 + ay*l_2;

								ay = (qf.at(2) - qi.at(2)) / 120 / EJzz / l;
								by = qi.at(2) / 24 / EJzz;
								cy = (240 * EJzz*(vy_0 - vy_l) + l*(120 * EJzz*(phiz_0 + phiz_l) - l_3*(3 * qf.at(2) + 7 * qi.at(2)))) / (120 * EJzz*l_3);
								dy = -(360 * EJzz*(vy_0 - vy_l) + l*(120 * EJzz*(2 * phiz_0 + phiz_l) - l_3*(2 * qf.at(2) + 3 * qi.at(2)))) / (120 * EJzz*l_2);
								ey = phiz_0;
								fy = vy_0;
							}
							else { // timoshenko formulation for transversal displacements
								//ay = bl.at(2) / 24 / EJzz;
								//dy = phiz_0;
								//fy = vy_0;
								//by = (2 / l*(vy_0 - vy_l) + (phiz_l + phiz_0)) / (l_2 + 12 * psi_y) - 2 * ay*l;
								//cy = -(72 * EJzz*GKyAy*l* (vy_0 - vy_l) + GKyAy*l_2*(24 * EJzz*(2 * phiz_0 + phiz_l) - l_3*bl.at(2)) + 12 * EJzz*(12 * EJzz*(phiz_0 - phiz_l) - l_3*bl.at(2))) / (24 * EJzz*l* (GKyAy*l_2 + 12 * EJzz));
								
								ay = (qf.at(2) - qi.at(2)) / 120 / EJzz / l;
								by = qi.at(2) / 24 / EJzz;
								cy = (240 * EJzz*GKyAy*(vy_0 - vy_l) + l*(GKyAy*(120 * EJzz*(phiz_0 + phiz_l) - l_3*(3 * qf.at(2) + 7 * qi.at(2))) - 40 * EJzz*l*(qf.at(2) + 2 * qi.at(2)))) / (120 * EJzz*l*(GKyAy*l_2 + 12 * EJzz));
								dy = -(360 * EJzz*GKyAy*l*(vy_0 - vy_l) + GKyAy*l_2*(120 * EJzz*(2 * phiz_0 + phiz_l) - l_3*(2 * qf.at(2) + 3 * qi.at(2))) + 30 * EJzz*(24 * EJzz*(phiz_0 - phiz_l) - l_3*(qf.at(2) + qi.at(2)))) / (120 * EJzz*l*(GKyAy*l_2 + 12 * EJzz));
								ey = phiz_0;
								fy = vy_0;
							}
						}

						if (hasWinklerZ) {
							FloatMatrix odeMtrx(4, 4);
							FloatArray rhs(4);
							FloatArray abcd(4);

							if (psi_z == 0.0) { // Euler Bernoulli formulation
								lambdaZ = sqrt(sqrt(wz / 4 / EJyy));

								odeMtrx.at(1, 1) = 1;
								odeMtrx.at(1, 3) = 1;
								odeMtrx.at(2, 1) = cos(l*lambdaZ) * exp(l*lambdaZ);
								odeMtrx.at(2, 2) = sin(l*lambdaZ) * exp(l*lambdaZ);
								odeMtrx.at(2, 3) = cos(l*lambdaZ) / exp(l*lambdaZ);
								odeMtrx.at(2, 4) = sin(l*lambdaZ) / exp(l*lambdaZ);
								odeMtrx.at(3, 1) = lambdaZ;
								odeMtrx.at(3, 2) = lambdaZ;
								odeMtrx.at(3, 3) = -lambdaZ;
								odeMtrx.at(3, 4) = lambdaZ;
								odeMtrx.at(4, 1) = lambdaZ * (cos(l*lambdaZ) - sin(l*lambdaZ)) * exp(l*lambdaZ);
								odeMtrx.at(4, 2) = lambdaZ * (cos(l*lambdaZ) + sin(l*lambdaZ)) * exp(l*lambdaZ);
								odeMtrx.at(4, 3) = -lambdaZ * (cos(l*lambdaZ) + sin(l*lambdaZ)) / exp(l*lambdaZ);
								odeMtrx.at(4, 4) = lambdaZ * (cos(l*lambdaZ) - sin(l*lambdaZ)) / exp(l*lambdaZ);
							}
							else {  // Timoshenko formulation
								alphaZ = wz / GKzAz;
								lambdaZ = wz / EJyy;
								deltaZ = alphaZ*alphaZ - 4 * lambdaZ;

								if (deltaZ > 0) {
									lambdaZ1 = sqrt(alphaZ / 2 + 0.5*sqrt(deltaZ));
									lambdaZ2 = sqrt(alphaZ / 2 - 0.5*sqrt(deltaZ));

									odeMtrx.at(1, 1) = 1;
									odeMtrx.at(1, 2) = 1;
									odeMtrx.at(1, 3) = 1;
									odeMtrx.at(1, 4) = 1;
									odeMtrx.at(2, 1) = exp(l*lambdaZ1);
									odeMtrx.at(2, 2) = exp(-l*lambdaZ1);
									odeMtrx.at(2, 3) = exp(l*lambdaZ2);
									odeMtrx.at(2, 4) = exp(-l*lambdaZ2);
									odeMtrx.at(3, 1) = lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1);
									odeMtrx.at(3, 2) = -lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1);
									odeMtrx.at(3, 3) = lambdaZ2*(lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1);
									odeMtrx.at(3, 4) = -lambdaZ2*(lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1);
									odeMtrx.at(4, 1) = lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1)*exp(l*lambdaZ1);
									odeMtrx.at(4, 2) = -lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1) / exp(l*lambdaZ1);
									odeMtrx.at(4, 3) = lambdaZ2*(lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1)*exp(l*lambdaZ2);
									odeMtrx.at(4, 4) = -lambdaZ2*(lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1) / exp(l*lambdaZ2);
								}
								else if (deltaZ == 0) {
									lambdaZ1 = sqrt(alphaZ / 2);
									lambdaZ2 = -sqrt(alphaZ / 2);

									odeMtrx.at(1, 1) = 1;
									odeMtrx.at(1, 2) = 1;
									odeMtrx.at(2, 1) = exp(l*lambdaZ1);
									odeMtrx.at(2, 2) = exp(-l*lambdaZ1);
									odeMtrx.at(2, 3) = l*exp(l*lambdaZ1);
									odeMtrx.at(2, 4) = l*exp(-l*lambdaZ1);
									odeMtrx.at(3, 1) = lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1);
									odeMtrx.at(3, 2) = -lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1);
									odeMtrx.at(3, 3) = lambdaZ1*(3 * lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1);
									odeMtrx.at(3, 4) = lambdaZ1*(3 * lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1);
									odeMtrx.at(4, 1) = lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1)*exp(l*lambdaZ1);
									odeMtrx.at(4, 2) = -lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1) / exp(l*lambdaZ1);
									odeMtrx.at(4, 3) = (l*lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1) + 3 * lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1)*exp(l*lambdaZ1);
									odeMtrx.at(4, 4) = (-l*lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1) + 3 * lambdaZ1*lambdaZ1*psi_z - alphaZ*psi_z + 1) / exp(l*lambdaZ1);
								}
								else {
									lambdaZ1 = sqrt(alphaZ / 4 + 0.5*sqrt(lambdaZ));
									lambdaZ2 = sqrt(-alphaZ / 4 + 0.5*sqrt(lambdaZ));

									odeMtrx.at(1, 1) = 1;
									odeMtrx.at(1, 2) = 1;
									odeMtrx.at(2, 1) = cos(l*lambdaZ2) * exp(l*lambdaZ1);
									odeMtrx.at(2, 2) = cos(l*lambdaZ2) / exp(l*lambdaZ1);
									odeMtrx.at(2, 3) = sin(l*lambdaZ2) * exp(l*lambdaZ1);
									odeMtrx.at(2, 4) = sin(l*lambdaZ2) / exp(l*lambdaZ1);
									odeMtrx.at(3, 1) = lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - 3 * lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1);
									odeMtrx.at(3, 2) = -lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - 3 * lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1);
									odeMtrx.at(3, 3) = lambdaZ2*(3 * lambdaZ1*lambdaZ1*psi_z - lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1);
									odeMtrx.at(3, 4) = lambdaZ2*(3 * lambdaZ1*lambdaZ1*psi_z - lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1);
									odeMtrx.at(4, 1) = (lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - 3 * lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1)*cos(l*lambdaZ2) - lambdaZ2*(3 * lambdaZ1*lambdaZ1*psi_z - lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1)*sin(l*lambdaZ2)) * exp(l*lambdaZ1);
									odeMtrx.at(4, 2) = -(lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - 3 * lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1)*cos(l*lambdaZ2) + lambdaZ2*(3 * lambdaZ1*lambdaZ1*psi_z - lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1)*sin(l*lambdaZ2)) / exp(l*lambdaZ1);
									odeMtrx.at(4, 3) = (lambdaZ2*(3 * lambdaZ1*lambdaZ1*psi_z - lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1)*cos(l*lambdaZ2) + lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - 3 * lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1)*sin(l*lambdaZ2)) * exp(l*lambdaZ1);
									odeMtrx.at(4, 4) = -(-lambdaZ2*(3 * lambdaZ1*lambdaZ1*psi_z - lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1)*cos(l*lambdaZ2) + lambdaZ1*(lambdaZ1*lambdaZ1*psi_z - 3 * lambdaZ2*lambdaZ2*psi_z - alphaZ*psi_z + 1)*sin(l*lambdaZ2)) / exp(l*lambdaZ1);
								}
							}

							//rhs.at(1) = dI.at(3) - bl.at(3) / wz;
							//rhs.at(2) = dE.at(3) - bl.at(3) / wz;
							//rhs.at(3) = -dI.at(5);
							//rhs.at(4) = -dE.at(5);

							rhs.at(1) = dI.at(3) - qi.at(3) / wy;
							rhs.at(2) = dE.at(3) - qf.at(3) / wy;
							rhs.at(3) = -dI.at(5) - (qf.at(3) - qi.at(3)) / (wz*l); // - o + per il secondo addendo?
							rhs.at(4) = -dE.at(5) - (qf.at(3) - qi.at(3)) / (wz*l);

							odeMtrx.solveForRhs(rhs, abcd);

							az = abcd.at(1);
							bz = abcd.at(2);
							cz = abcd.at(3);
							dz = abcd.at(4);
						}
						else {
							if (psi_z == 0.0) {  // Euler Bernoulli formulation
								//az = bl.at(3) / 24 / EJyy;
								//dz = phiy_0;
								//fz = vz_0;
								//bz = (2 * (vz_0 - vz_l) + l*(phiy_l + phiy_0)) / (l_3)-2 * az*l;
								//cz = -(3 * (vz_0 - vz_l) + l*(2 * phiy_0 + phiy_l)) / l_2 + az*l_2;

								az = (qf.at(3) - qi.at(3)) / 120 / EJyy / l;
								bz = qi.at(3) / 24 / EJyy;
								cz = (240 * EJyy*(vy_0 - vy_l) + l*(120 * EJyy*(phiy_0 + phiy_l) - l_3*(3 * qf.at(3) + 7 * qi.at(3)))) / (120 * EJyy*l_3);
								dz = -(360 * EJyy*(vy_0 - vy_l) + l*(120 * EJyy*(2 * phiy_0 + phiy_l) - l_3*(2 * qf.at(3) + 3 * qi.at(3)))) / (120 * EJyy*l_2);
								ez = phiy_0;
								fz = vz_0;
							}
							else { // timoshenko formulation for transversal displacements
								//az = bl.at(3) / 24 / EJyy;
								//dz = phiy_0;
								//fz = vz_0;
								//bz = (2 / l*(vz_0 - vz_l) + (phiy_l + phiy_0)) / (l_2 + 12 * psi_z) - 2 * az*l;
								//cz = -(72 * EJyy*GKzAz*l* (vz_0 - vz_l) + GKzAz*l_2*(24 * EJyy*(2 * phiy_0 + phiy_l) - l_3*bl.at(3)) + 12 * EJyy*(12 * EJyy*(phiy_0 - phiy_l) - l_3*bl.at(3))) / (24 * EJyy*l* (GKzAz*l_2 + 12 * EJyy));
							
								az = (qf.at(3) - qi.at(3)) / 120 / EJyy / l;
								bz = qi.at(3) / 24 / EJyy;
								cz = (240 * EJyy*GKzAz*(vz_0 - vz_l) + l*(GKzAz*(120 * EJyy*(phiy_0 + phiy_l) - l_3*(3 * qf.at(3) + 7 * qi.at(3))) - 40 * EJyy*l*(qf.at(3) + 2 * qi.at(3)))) / (120 * EJyy*l*(GKzAz*l_2 + 12 * EJyy));
								dz = -(360 * EJyy*GKzAz*l*(vz_0 - vz_l) + GKzAz*l_2*(120 * EJyy*(2 * phiy_0 + phiy_l) - l_3*(2 * qf.at(3) + 3 * qi.at(3))) + 30 * EJyy*(24 * EJyy*(phiy_0 - phiy_l) - l_3*(qf.at(3) + qi.at(3)))) / (120 * EJyy*l*(GKzAz*l_2 + 12 * EJyy));
								ez = phiy_0;
								fz = vz_0;
							}
						}

						// axial displacements
						//anx = -bl.at(1) / 2 / EA;
						//bnx = (dx_l - dx_0) / l - anx*l;
						//cnx = dx_0;

						anx = (qf.at(1) - qi.at(1)) / (6*EA*l);
						bnx = qi.at(1) / (2*EA);
						cnx = -(6 * EA*(dx_0 - dx_l) + l_2*(qf.at(1) + 2 * qi.at(1))) / (6 * EA*l);
						dnx = dx_0;

						// torsional rotations
						//atx = -bl.at(4) / 2 / GJ;
						//btx = (tx_l - tx_0) / l - atx*l;
						//ctx = tx_0;

						atx = (qf.at(1) - qi.at(1)) / (6 * GJ*l);
						btx = qi.at(1) / (2 * GJ);
						ctx = -(6 * GJ*(tx_0 - tx_l) + l_2*(qf.at(1) + 2 * qi.at(1))) / (6 * GJ*l);
						dtx = tx_0;
						
						calc = true;
					}

					FloatArray disps(6);
					FloatArray wink(6); // winkler reactions
					disps.at(1) = anx*pos_3 + bnx*pos_2 + cnx*pos + dnx;
					double lamxY, lamxZ, lamxY1, lamxY2, lamxZ1, lamxZ2;
					lamxY = lambdaY*pos;
					lamxZ = lambdaZ*pos;
					lamxY1 = lambdaY1*pos;
					lamxZ1 = lambdaZ1*pos;
					lamxY2 = lambdaY2*pos;
					lamxZ2 = lambdaZ2*pos;

					if (hasWinklerY)
					{
						if (psi_y == 0.0) {
							// displacement
							// disps.at(2) = exp(lamxY)*(ay*cos(lamxY) + by*sin(lamxY)) + (cy*cos(lamxY) + dy*sin(lamxY)) / exp(lamxY) + bl.at(2) / wy;
							disps.at(2) = exp(lamxY)*(ay*cos(lamxY) + by*sin(lamxY)) + (cy*cos(lamxY) + dy*sin(lamxY)) / exp(lamxY) + (qi.at(2)+(qf.at(2)-qi.at(2))*pos/l) / wy;
							// rotation
							disps.at(6) = exp(lamxY)*(lambdaY*(ay + by)*cos(lamxY) + lambdaY*(by - ay)*sin(lamxY)) - (lambdaY*(cy - dy)*cos(lamxY) + lambdaY*(cy + dy)*sin(lamxY)) / (exp(lamxY));
							// adjust the diagrams
							BeamForces[elem->giveNumber()].at(pos).at(6) = 2 * lambdaY*lambdaY*EJzz* (exp(lamxY)*(by*cos(lamxY) - ay*sin(lamxY)) + (-dy*cos(lamxY) + cy*sin(lamxY)) / exp(lamxY));
							BeamForces[elem->giveNumber()].at(pos).at(2) = -2 * lambdaY*lambdaY*lambdaY*EJzz* (-exp(lamxY)*((ay - by)*cos(lamxY) + (by + ay)*sin(lamxY)) + ((cy + dy)*cos(lamxY) + (-cy + dy)*sin(lamxY)) / (exp(lamxY)));
						}
						else {
							if (deltaY > 0) {
								// displacement
								disps.at(2) = ay*exp(lamxY1) + by / exp(lamxY1) + cy * exp(lamxY2) + dy / exp(lamxY2) + (qi.at(2) + (qf.at(2) - qi.at(2))*pos / l) / wy;
								// rotation
								disps.at(6) = psi_y*(ay*lambdaY1*lambdaY1*lambdaY1*exp(lamxY1)-by*lambdaY1*lambdaY1*lambdaY1*exp(-lamxY1) + cy*lambdaY2*lambdaY2*lambdaY2*exp(lamxY2)-dy*lambdaY2*lambdaY2*lambdaY2*exp(-lamxY2)) - (alphaY*psi_y - 1)*(ay*lambdaY1*exp(lamxY1)-by*lambdaY1*exp(-lamxY1) + cy*lambdaY2*exp(lamxY2)-dy*lambdaY2*exp(-lamxY2));
								// adjust the diagrams
								BeamForces[elem->giveNumber()].at(pos).at(6) = EJzz* (ay*lambdaY1*lambdaY1*exp(lamxY1) + by*lambdaY1*lambdaY1*exp(-lamxY1) + cy*lambdaY2*lambdaY2*exp(lamxY2) + dy*lambdaY2*lambdaY2*exp(-lamxY2) - alphaY*(ay*exp(lamxY1) + by*exp(-lamxY1) + cy*exp(lamxY2) + dy*exp(-lamxY2)));
								BeamForces[elem->giveNumber()].at(pos).at(2) = -EJzz* (ay*lambdaY1*lambdaY1*lambdaY1*exp(lamxY1) - by*lambdaY1*lambdaY1*lambdaY1*exp(-lamxY1) + cy*lambdaY2*lambdaY2*lambdaY2*exp(lamxY2) - dy*lambdaY2*lambdaY2*lambdaY2*exp(-lamxY2) + alphaY*(ay*lambdaY1*exp(lamxY1) - by*lambdaY1*exp(-lamxY1) + cy*lambdaY2*exp(lamxY2) - dy*lambdaY2*exp(-lamxY2)));
							}
							else if (deltaY == 0) { 
								// displacement
								disps.at(2) = ay*exp(lamxY1) + by / exp(lamxY1) + pos*(cy * exp(lamxY1) + dy / exp(lamxY1)) + (qi.at(2) + (qf.at(2) - qi.at(2))*pos / l) / wy;
								// rotation
								disps.at(6) = psi_y*(lambdaY1*lambdaY1*exp(lamxY1)*(cy*lamxY1 + ay*lambdaY1 + 3*cy) - lambdaY1*lambdaY1*exp(-lamxY1)*(dy*lamxY1 + by*lambdaY1 - 3*dy)) - (alphaY*psi_y - 1)*(exp(lamxY1)*(cy*lamxY1 + ay*lambdaY1 + cy) - exp(-lamxY1)*(dy*lamxY1 + by*lambdaY1 - dy));
								// adjust the diagrams
								BeamForces[elem->giveNumber()].at(pos).at(6) = EJzz*(lambdaY1*exp(lamxY1)*(cy*lamxY1 + ay*lambdaY1 + 2 * cy) + lambdaY1*exp(-lamxY1)*(dy*lamxY1 + by*lambdaY1 - 2 * dy) - alphaY*(ay*exp(lamxY1) + by*exp(-lamxY1) + cy*pos*exp(lamxY1) + dy*pos*exp(-lamxY1)));
								BeamForces[elem->giveNumber()].at(pos).at(2) = -EJzz*(lambdaY1*lambdaY1*exp(lamxY1)*(cy*lamxY1 + ay*lambdaY1 + 3 * cy) - lambdaY1*lambdaY1*exp(-lamxY1)*(dy*lamxY1 + by*lambdaY1 - 3 * dy) + alphaY*(exp(lamxY1)*(cy*lamxY1 + ay*lambdaY1 + cy) - exp(-lamxY1)*(dy*lamxY1 + by*lambdaY1 - dy)));
							}
							else {
								// displacement
								disps.at(2) = exp(lamxY1)*(ay*cos(lamxY2) + cy*sin(lamxY2)) + (by*cos(lamxY2) + dy*sin(lamxY2)) / exp(lamxY1) + (qi.at(2) + (qf.at(2) - qi.at(2))*pos / l) / wy;
								// rotation
								disps.at(6) = psi_y*(exp(lamxY1)*((ay*lambdaY1*(lambdaY1*lambdaY1 - 3*lambdaY2*lambdaY2) + cy*lambdaY2*(3*lambdaY1*lambdaY1 - lambdaY2*lambdaY2))*cos(lamxY2) - (ay*lambdaY2*(3*lambdaY1*lambdaY1 - lambdaY2*lambdaY2) + cy*lambdaY1*(3*lambdaY2*lambdaY2 - lambdaY1*lambdaY1))*sin(lamxY2)) - exp(-lamxY1)*((by*lambdaY1*(lambdaY1*lambdaY1 - 3*lambdaY2*lambdaY2) + dy*lambdaY2*(lambdaY2*lambdaY2 - 3*lambdaY1*lambdaY1))*cos(lamxY2) + (by*lambdaY2*(3*lambdaY1*lambdaY1 - lambdaY2*lambdaY2) + dy*lambdaY1*(lambdaY1*lambdaY1 - 3*lambdaY2*lambdaY2))*sin(lamxY2))) - (alphaY*psi_y - 1)*(exp(lamxY1)*((ay*lambdaY1 + cy*lambdaY2)*cos(lamxY2) + (cy*lambdaY1 - ay*lambdaY2)*sin(lamxY2)) - exp(-lamxY1)*((by*lambdaY1 - dy*lambdaY2)*cos(lamxY2) + (by*lambdaY2 + dy*lambdaY1)*sin(lamxY2)));
								// adjust the diagrams
								BeamForces[elem->giveNumber()].at(pos).at(6) = EJzz *(exp(lamxY1)*((ay*(lambdaY1*lambdaY1 - lambdaY2*lambdaY2) + 2 * cy*lambdaY1*lambdaY2)*cos(lamxY2) - (2 * ay*lambdaY1*lambdaY2 + cy*(lambdaY2*lambdaY2 - lambdaY1*lambdaY1))*sin(lamxY2)) + exp(-lamxY1)*((by*(lambdaY1*lambdaY1 - lambdaY2*lambdaY2) - 2 * dy*lambdaY1*lambdaY2)*cos(lamxY2) + (2 * by*lambdaY1*lambdaY2 + dy*(lambdaY1*lambdaY1 - lambdaY2*lambdaY2))*sin(lamxY2)) - alphaY*(ay*exp(lamxY1)*cos(lamxY2) + by*exp(-lamxY1)*cos(lamxY2) + cy*exp(lamxY1)*sin(lamxY2) + dy*exp(-lamxY1)*sin(lamxY2)));
								BeamForces[elem->giveNumber()].at(pos).at(2) = -EJzz *(exp(lamxY1)*((ay*lambdaY1*(lambdaY1*lambdaY1 - 3 * lambdaY2*lambdaY2) + cy*lambdaY2*(3 * lambdaY1*lambdaY1 - lambdaY2*lambdaY2))*cos(lamxY2) - (ay*lambdaY2*(3 * lambdaY1*lambdaY1 - lambdaY2*lambdaY2) + cy*lambdaY1*(3 * lambdaY2*lambdaY2 - lambdaY1*lambdaY1))*sin(lamxY2)) - exp(-lamxY1)*((by*lambdaY1*(lambdaY1*lambdaY1 - 3 * lambdaY2*lambdaY2) + dy*lambdaY2*(lambdaY2*lambdaY2 - 3 * lambdaY1*lambdaY1))*cos(lamxY2) + (by*lambdaY2*(3 * lambdaY1*lambdaY1 - lambdaY2*lambdaY2) + dy*lambdaY1*(lambdaY1*lambdaY1 - 3 * lambdaY2*lambdaY2))*sin(lamxY2)) + alphaY*(exp(lamxY1)*((ay*lambdaY1 + cy*lambdaY2)*cos(lamxY2) + (cy*lambdaY1 - ay*lambdaY2)*sin(lamxY2)) - exp(-lamxY1)*((by*lambdaY1 - dy*lambdaY2)*cos(lamxY2) + (by*lambdaY2 + dy*lambdaY1)*sin(lamxY2))));
							}
						}

						wink.at(2) = -disps.at(2)*wy;
					}
					else{
						// displacement
						if (psi_y == 0.0) {
							disps.at(2) = ay*pos_5 + by*pos_4 + cy*pos_3 + dy*pos_2 + ey*pos + fy;
						}
						else {
							disps.at(2) = ay*pos_5 + by*pos_4 + (cy - 20 * ay*psi_y)*pos_3 + (dy - 12 * by*psi_y)*pos_2 + (ey-6*cy*psi_y)*pos + fy;
						}
						// rotation
						//disps.at(6) = 4 * ay*pos_3 + 3 * by*pos_2 + 2 * cy*pos + dy;
						disps.at(6) = 5 * ay*pos_4 + 4 * by*pos_3 + 3 * cy*pos_2 + 2*dy*pos + ey;
					}

					if (hasWinklerZ) {
						if (psi_z == 0.0) {
							// displacement
							disps.at(3) = exp(lamxZ)*(az*cos(lamxZ) + bz*sin(lamxZ)) + (cz*cos(lamxZ) + dz*sin(lamxZ)) / exp(lamxZ) + (qi.at(3) + (qf.at(3) - qi.at(3))*pos / l) / wz;
							// rotation
							disps.at(5) = -(exp(lamxZ)*(lambdaZ*(az + bz)*cos(lamxZ) + lambdaZ*(bz - az)*sin(lamxZ)) - (lambdaZ*(cz - dz)*cos(lamxZ) + lambdaZ*(cz + dz)*sin(lamxZ)) / (exp(lamxZ)));
							// adjust the diagrams
							BeamForces[elem->giveNumber()].at(pos).at(5) = -2 * lambdaZ*lambdaZ*EJyy* (exp(lamxZ)*(bz*cos(lamxZ) - az*sin(lamxZ)) + (-dz*cos(lamxZ) + cz*sin(lamxZ)) / exp(lamxZ));
							BeamForces[elem->giveNumber()].at(pos).at(3) = -2 * lambdaZ*lambdaZ*lambdaZ*EJyy* (-exp(lamxZ)*((az - bz)*cos(lamxZ) + (bz + az)*sin(lamxZ)) + ((cz + dz)*cos(lamxZ) + (-cz + dz)*sin(lamxZ)) / (exp(lamxZ)));
						}
						else {
							if (deltaZ > 0) {
								// displacement
								disps.at(3) = az*exp(lamxZ1) + bz / exp(lamxZ1) + cz * exp(lamxZ2) + dz / exp(lamxZ2) + (qi.at(3) + (qf.at(3) - qi.at(3))*pos / l) / wz;
								// rotation
								disps.at(5) = -psi_z*(az*lambdaZ1*lambdaZ1*lambdaZ1*exp(lamxZ1) - bz*lambdaZ1*lambdaZ1*lambdaZ1*exp(-lamxZ1) + cz*lambdaZ2*lambdaZ2*lambdaZ2*exp(lamxZ2) - dz*lambdaZ2*lambdaZ2*lambdaZ2*exp(-lamxZ2)) + (alphaZ*psi_z - 1)*(az*lambdaZ1*exp(lamxZ1) - bz*lambdaZ1*exp(-lamxZ1) + cz*lambdaZ2*exp(lamxZ2) - dz*lambdaZ2*exp(-lamxZ2));
								// adjust the diagrams
								BeamForces[elem->giveNumber()].at(pos).at(5) = -EJyy* (az*lambdaZ1*lambdaZ1*exp(lamxZ1) + bz*lambdaZ1*lambdaZ1*exp(-lamxZ1) + cz*lambdaZ2*lambdaZ2*exp(lamxZ2) + dz*lambdaZ2*lambdaZ2*exp(-lamxZ2) - alphaZ*(az*exp(lamxZ1) + bz*exp(-lamxZ1) + cz*exp(lamxZ2) + dz*exp(-lamxZ2)));
								BeamForces[elem->giveNumber()].at(pos).at(3) = -EJyy* (az*lambdaZ1*lambdaZ1*lambdaZ1*exp(lamxZ1) - bz*lambdaZ1*lambdaZ1*lambdaZ1*exp(-lamxZ1) + cz*lambdaZ2*lambdaZ2*lambdaZ2*exp(lamxZ2) - dz*lambdaZ2*lambdaZ2*lambdaZ2*exp(-lamxZ2) + alphaZ*(az*lambdaZ1*exp(lamxZ1) - bz*lambdaZ1*exp(-lamxZ1) + cz*lambdaZ2*exp(lamxZ2) - dz*lambdaZ2*exp(-lamxZ2)));
							}
							else if (deltaZ == 0) {
								// displacement
								disps.at(3) = az*exp(lamxZ1) + bz / exp(lamxZ1) + pos*(cz * exp(lamxZ1) + dz / exp(lamxZ1)) + (qi.at(3) + (qf.at(3) - qi.at(3))*pos / l) / wz;
								// rotation
								disps.at(5) = -psi_z*(lambdaZ1*lambdaZ1*exp(lamxZ1)*(cz*lamxZ1 + az*lambdaZ1 + 3 * cz) - lambdaZ1*lambdaZ1*exp(-lamxZ1)*(dz*lamxZ1 + bz*lambdaZ1 - 3 * dz)) + (alphaZ*psi_z - 1)*(exp(lamxZ1)*(cz*lamxZ1 + az*lambdaZ1 + cz) - exp(-lamxZ1)*(dz*lamxZ1 + bz*lambdaZ1 - dz));
								// adjust the diagrams
								BeamForces[elem->giveNumber()].at(pos).at(5) = -EJyy*(lambdaZ1*exp(lamxZ1)*(cz*lamxZ1 + az*lambdaZ1 + 2 * cz) + lambdaZ1*exp(-lamxZ1)*(dz*lamxZ1 + bz*lambdaZ1 - 2 * dz) - alphaZ*(az*exp(lamxZ1) + bz*exp(-lamxZ1) + cz*pos*exp(lamxZ1) + dz*pos*exp(-lamxZ1)));
								BeamForces[elem->giveNumber()].at(pos).at(3) = -EJyy*(lambdaZ1*lambdaZ1*exp(lamxZ1)*(cz*lamxZ1 + az*lambdaZ1 + 3 * cz) - lambdaZ1*lambdaZ1*exp(-lamxZ1)*(dz*lamxZ1 + bz*lambdaZ1 - 3 * dz) + alphaZ*(exp(lamxZ1)*(cz*lamxZ1 + az*lambdaZ1 + cz) - exp(-lamxZ1)*(dz*lamxZ1 + bz*lambdaZ1 - dz)));
							}
							else {
								// displacement
								disps.at(3) = exp(lamxZ1)*(az*cos(lamxZ2) + cz*sin(lamxZ2)) + (bz*cos(lamxZ2) + dz*sin(lamxZ2)) / exp(lamxZ1) + (qi.at(3) + (qf.at(3) - qi.at(3))*pos / l) / wz;
								// rotation
								disps.at(5) = -psi_z*(exp(lamxZ1)*((az*lambdaZ1*(lambdaZ1*lambdaZ1 - 3 * lambdaZ2*lambdaZ2) + cz*lambdaZ2*(3 * lambdaZ1*lambdaZ1 - lambdaZ2*lambdaZ2))*cos(lamxZ2) - (az*lambdaZ2*(3 * lambdaZ1*lambdaZ1 - lambdaZ2*lambdaZ2) + cz*lambdaZ1*(3 * lambdaZ2*lambdaZ2 - lambdaZ1*lambdaZ1))*sin(lamxZ2)) - exp(-lamxZ1)*((bz*lambdaZ1*(lambdaZ1*lambdaZ1 - 3 * lambdaZ2*lambdaZ2) + dz*lambdaZ2*(lambdaZ2*lambdaZ2 - 3 * lambdaZ1*lambdaZ1))*cos(lamxZ2) + (bz*lambdaZ2*(3 * lambdaZ1*lambdaZ1 - lambdaZ2*lambdaZ2) + dz*lambdaZ1*(lambdaZ1*lambdaZ1 - 3 * lambdaZ2*lambdaZ2))*sin(lamxZ2))) + (alphaZ*psi_z - 1)*(exp(lamxZ1)*((az*lambdaZ1 + cz*lambdaZ2)*cos(lamxZ2) + (cz*lambdaZ1 - az*lambdaZ2)*sin(lamxZ2)) - exp(-lamxZ1)*((bz*lambdaZ1 - dz*lambdaZ2)*cos(lamxZ2) + (bz*lambdaZ2 + dz*lambdaZ1)*sin(lamxZ2)));
								// adjust the diagrams
								BeamForces[elem->giveNumber()].at(pos).at(5) = -EJyy *(exp(lamxZ1)*((az*(lambdaZ1*lambdaZ1 - lambdaZ2*lambdaZ2) + 2 * cz*lambdaZ1*lambdaZ2)*cos(lamxZ2) - (2 * az*lambdaZ1*lambdaZ2 + cz*(lambdaZ2*lambdaZ2 - lambdaZ1*lambdaZ1))*sin(lamxZ2)) + exp(-lamxZ1)*((bz*(lambdaZ1*lambdaZ1 - lambdaZ2*lambdaZ2) - 2 * dz*lambdaZ1*lambdaZ2)*cos(lamxZ2) + (2 * bz*lambdaZ1*lambdaZ2 + dz*(lambdaZ1*lambdaZ1 - lambdaZ2*lambdaZ2))*sin(lamxZ2)) - alphaZ*(az*exp(lamxZ1)*cos(lamxZ2) + bz*exp(-lamxZ1)*cos(lamxZ2) + cz*exp(lamxZ1)*sin(lamxZ2) + dz*exp(-lamxZ1)*sin(lamxZ2)));
								BeamForces[elem->giveNumber()].at(pos).at(3) = -EJyy *(exp(lamxZ1)*((az*lambdaZ1*(lambdaZ1*lambdaZ1 - 3 * lambdaZ2*lambdaZ2) + cz*lambdaZ2*(3 * lambdaZ1*lambdaZ1 - lambdaZ2*lambdaZ2))*cos(lamxZ2) - (az*lambdaZ2*(3 * lambdaZ1*lambdaZ1 - lambdaZ2*lambdaZ2) + cz*lambdaZ1*(3 * lambdaZ2*lambdaZ2 - lambdaZ1*lambdaZ1))*sin(lamxZ2)) - exp(-lamxZ1)*((bz*lambdaZ1*(lambdaZ1*lambdaZ1 - 3 * lambdaZ2*lambdaZ2) + dz*lambdaZ2*(lambdaZ2*lambdaZ2 - 3 * lambdaZ1*lambdaZ1))*cos(lamxZ2) + (bz*lambdaZ2*(3 * lambdaZ1*lambdaZ1 - lambdaZ2*lambdaZ2) + dz*lambdaZ1*(lambdaZ1*lambdaZ1 - 3 * lambdaZ2*lambdaZ2))*sin(lamxZ2)) + alphaZ*(exp(lamxZ1)*((az*lambdaZ1 + cz*lambdaZ2)*cos(lamxZ2) + (cz*lambdaZ1 - az*lambdaZ2)*sin(lamxZ2)) - exp(-lamxZ1)*((bz*lambdaZ1 - dz*lambdaZ2)*cos(lamxZ2) + (bz*lambdaZ2 + dz*lambdaZ1)*sin(lamxZ2))));
							}
						}
						wink.at(3) = -disps.at(3)*wz;
					}
					else {
						// displacement
						if (psi_z == 0) {
							disps.at(3) = az*pos_5 + bz*pos_4 + cz*pos_3 + dz*pos_2 + ez*pos + fz;
						}
						else {
							disps.at(3) = az*pos_5 + bz*pos_4 + (cz - 20 * az*psi_z)*pos_3 + (dz - 12 * bz*psi_z)*pos_2 + (ez - 6 * cz*psi_z) * pos + fz;
						}
						// rotation
						// disps.at(5) = -(4 * az*pos_3 + 3 * bz*pos_2 + 2 * cz*pos + dz);  // inverted signs for rotations about y.
						disps.at(5) = -(5 * az*pos_4 + 4 * bz*pos_3 + 3 * cz*pos_2 + 2*dz*pos + ez);  // inverted signs for rotations about y.
					}
					
					disps.at(4) = atx*pos_3 + btx*pos_2 + ctx*pos + dtx;

					//disps -= (dI+ddN*ksi);

					DispDict[pos] = disps;
					WinkDict[pos] = wink;

					//ipDisp.beProductOf(shapeFunctions, rl);
				}

				DispDict[l] = dE; // -dNE;
				WinkDict[l] = WL;

				// save the displacements
				BeamDisplacements[elem->giveNumber()] = DispDict;
				BeamWinkler[elem->giveNumber()] = WinkDict;
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
			BeamWinklerList.push_back(BeamWinkler);
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
				BeamWinkler = combBeamWinkler;

#ifdef DEBUG
				map<double, FloatArray> elem = combBeamDisplacements[7];
				OOFEM_WARNING("Stop %e", elem.at(4).at(1));
#endif

				combBeamDisplacements.clear();
				combBeamForces.clear();
				combBeamWinkler.clear();
			}

			double curTime = tStep->giveTargetTime();
			map<int, map<double, FloatArray>>::iterator BForces_it = BeamForces.begin();
			map<int, map<double, FloatArray>>::iterator BDisps_it = BeamDisplacements.begin();
			map<int, map<double, FloatArray>>::iterator BWinks_it = BeamWinkler.begin();
			for (;
				BForces_it != BeamForces.end();
				++BForces_it, ++BDisps_it, ++BWinks_it)
			{
				map< double, FloatArray > &BForces = BForces_it->second;
				map< double, FloatArray > &BDisps = BDisps_it->second;
				map< double, FloatArray > &BWinks = BWinks_it->second;
				Element* elem = d->giveElement(BForces_it->first);
				int ID = elem->giveLabel();

				map<double, FloatArray >::iterator forces_it = BForces.begin();
				map<double, FloatArray >::iterator disps_it = BDisps.begin();
				map<double, FloatArray >::iterator winks_it = BWinks.begin();

				for (;
					forces_it != BForces.end();
					++forces_it, ++disps_it, ++winks_it)
				{
					double pos = forces_it->first;
					FloatArray forces = forces_it->second;
					FloatArray disps = disps_it->second;
					FloatArray winks = winks_it->second;
					fprintf(this->stream, "%10.5e;%d;%10.5e;", curTime, ID, pos);

					for (auto &val : forces) {
						fprintf(this->stream, "%10.5e;", val);
					}
					for (auto &val : disps) {
						fprintf(this->stream, "%10.5e;", val);
					}
					fprintf(this->stream, "%10.5e;", winks.at(2));
					fprintf(this->stream, "%10.5e;", winks.at(3));

					fprintf(this->stream, "\n");
				}

			}

#ifdef MEMSTR
			fprintf(this->stream, "endStep\n");
			//if (usestream) fprintf(this->stream, "endStep\n");
#endif
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
			BeamWinkler.clear();

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
					destRespArray.at(i) += fabs(srcRespArray.at(i)*srcRespArray2.at(i)*fact);
				}
			}
		}
	}

	void BeamExportModule::calcRoot(map<int, map<double, FloatArray>> &answer)
	{
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


	void BeamExportModule::correctSigns(map<int, map<double, FloatArray>> &answer, map<int, map<double, FloatArray>> &src, bool use1stmode)
	{
		map<int, map<double, FloatArray>>::iterator destElem_it = answer.begin();
		map<int, map<double, FloatArray>>::iterator srcElem_it = src.begin();
		for (; destElem_it != answer.end(); ++destElem_it, ++srcElem_it)
		{
			map<double, FloatArray> &destRespMap = destElem_it->second;
			map<double, FloatArray> &srcRespMap = srcElem_it->second; // I'm using this from the dominant mode only

			map<double, FloatArray>::iterator destRespMap_it = destRespMap.begin();
			map<double, FloatArray>::iterator srcRespMap_it = srcRespMap.begin();
			for (; destRespMap_it != destRespMap.end(); ++destRespMap_it, ++srcRespMap_it)
			{
				FloatArray &destRespArray = destRespMap_it->second;
				FloatArray &srcRespArray = srcRespMap_it->second;

				for (int i = 1; i <= destRespArray.giveSize(); i++)
				{
					// square it and add it
					double res = destRespArray.at(i);
					if (use1stmode) res *= signbit(srcRespArray.at(i)) ? -1 : 1;
					destRespArray.at(i) = res;
				}
			}
		}
	}

	void BeamExportModule::SRSS(){
		int dominantMode;
		rs->giveDominantMode(dominantMode);

		list<map<int, map<double, FloatArray>>>::iterator disps_it = BeamDisplacementsList.begin();
		for (; disps_it != BeamDisplacementsList.end(); ++disps_it)
		{
			addMultiply(combBeamDisplacements, *disps_it, *disps_it);
		}
		calcRoot(combBeamDisplacements);
		disps_it = std::next(BeamDisplacementsList.begin(), dominantMode - 1);  // dominant mode only
		correctSigns(combBeamDisplacements, *disps_it, true);

		list<map<int, map<double, FloatArray>>>::iterator forces_it = BeamForcesList.begin();
		for (; forces_it != BeamForcesList.end(); ++forces_it)
		{
			addMultiply(combBeamForces, *forces_it, *forces_it);  // mult by 1.0
		}
		calcRoot(combBeamForces);

		list<map<int, map<double, FloatArray>>>::iterator winks_it = BeamWinklerList.begin();
		for (; winks_it != BeamWinklerList.end(); ++winks_it)
		{
			addMultiply(combBeamWinkler, *winks_it, *winks_it);  // mult by 1.0
		}
		calcRoot(combBeamWinkler);
		winks_it = std::next(BeamWinklerList.begin(), dominantMode - 1);	// dominant mode only
		correctSigns(combBeamWinkler, *winks_it, true);

	}

	void BeamExportModule::CQC(){
		FloatMatrix rhos;
		rs->giveRhos(rhos);
		int dominantMode;
		rs->giveDominantMode(dominantMode);

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
		disps_it = std::next(BeamDisplacementsList.begin(), dominantMode - 1);  // dominant mode only
		correctSigns(combBeamDisplacements, *disps_it, true);

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

		list<map<int, map<double, FloatArray>>>::iterator winks_it = BeamWinklerList.begin();

		for (int i = 1; winks_it != BeamWinklerList.end(); ++winks_it, i++)
		{
			list<map<int, map<double, FloatArray>>>::iterator winks_it2 = BeamWinklerList.begin();
			for (int j = 1; winks_it2 != BeamWinklerList.end(); ++winks_it2, j++)
			{
				addMultiply(combBeamWinkler, *winks_it, *winks_it2, rhos.at(i, j));
			}
		}
		calcRoot(combBeamWinkler);
		winks_it = std::next(BeamWinklerList.begin(), dominantMode - 1);	// dominant mode only
		correctSigns(combBeamWinkler, *winks_it, true);
	}

	void
		BeamExportModule::initialize()
	{
#ifdef MEMSTR
		this->stream = nullptr;
		FILE *source = classFactory.giveMemoryStream("bem");
		int sourceFD = _open_osfhandle((intptr_t)source, _O_APPEND);
		if (sourceFD != -1) {
			this->stream = _fdopen(sourceFD, "a");
		}
		if (!(this->stream)) {  // if not, write to file
#endif
			string fileName = emodel->giveOutputBaseFileName() + ".bem";
			if ((this->stream = fopen(fileName.c_str(), "w")) == NULL) {
				OOFEM_ERROR("failed to open file %s", fileName.c_str());
			}
#ifdef MEMSTR
			usestream = false;
		}
#endif
		// ";" as separator
		fprintf(this->stream, "#Time;BeamNo;DistanceFromI;N_x;T_y;T_z;M_x;M_y;M_z;dx;dy;dz;rx;ry;rz;q_winky;q_winkz;");
		//for ( int var: this->ists ) {
		//    fprintf(this->stream, "%s    ", __InternalStateTypeToString( ( InternalStateType ) var) );
		//}
		fprintf(this->stream, "\n");
		fflush(this->stream);
	}

	void
		BeamExportModule::terminate()
	{
#ifdef MEMSTR
		fprintf(this->stream, "strTerm\n");
		//if (usestream) fprintf(this->stream, "strTerm\n");
#endif
		fflush(this->stream);
// #ifndef MEMSTR
		fclose(this->stream);
// #endif
	}
} // end namespace oofem
