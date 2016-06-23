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

#include "NodalRecoveryModule.h"
#include "timestep.h"
#include "element.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/Beams/beam3d.h"
#include "gausspoint.h"
#include "engngm.h"
#include "material.h"
#include "classfactory.h"
#include "generalboundarycondition.h"
#include "constantedgeload.h"
#include "fei3dlinelin.h"
#include "inputrecord.h"
#include "../sm/engineeringmodels/responseSpectrum.h"
#include "cltypes.h"
#include "materialinterface.h"

#include <vector>

using namespace std;

namespace oofem {
	REGISTER_ExportModule(NodalRecoveryModule)

		NodalRecoveryModule::NodalRecoveryModule(int n, EngngModel *e) : ExportModule(n, e) { }

	NodalRecoveryModule :: ~NodalRecoveryModule() { }

	IRResultType
		NodalRecoveryModule::initializeFrom(InputRecord *ir)
	{
		IRResultType result;                 // Required by IR_GIVE_FIELD macro
		int val;

		val = NodalRecoveryModel::NRM_NodalAveraging;
		IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_NodalRecoveryModule_stype);
		stype = (NodalRecoveryModel::NodalRecoveryModelType) val;

		IR_GIVE_OPTIONAL_FIELD(ir, internalVarsToExport, _IFT_NodalRecoveryModule_rtypes);
		if (internalVarsToExport.giveSize() == 0) OOFEM_ERROR("NodalRecoveryModule - No response types defined");
		for (int i = 0; internalVarsToExport.giveSize(); i++) {
			InternalStateType ist = static_cast<InternalStateType>(internalVarsToExport.at(i));
			if (ist==InternalStateType::IST_Undefined) OOFEM_ERROR("NodalRecoveryModule - Invalid response type defined at pos %d", i);
		}

		isRespSpec = ir->hasField(_IFT_NodalRecoveryModule_isrespspec);

		if (isRespSpec) {
			const char* name = this->emodel->giveClassName();
			if (!strcmp(name, "ResponseSpectrum") == 0) OOFEM_ERROR("NodalRecoveryModule - Using rspec mode without a ResponseSpectrum engineering model");
			rs = dynamic_cast<ResponseSpectrum*>(this->emodel);
			if (!rs) OOFEM_ERROR("NodalRecoveryModule - Error retrieving engmodel.");
		}

		return ExportModule::initializeFrom(ir);
	}

	void
		NodalRecoveryModule::exportIntVars(FILE *stream, TimeStep *tStep)
	{
		// should be performed over regions

		int i, n = internalVarsToExport.giveSize();
		//int nnodes;
		//Domain *d = emodel->giveDomain(1);
		InternalStateType type;

		if (n == 0) {
			return;
		}

		for (i = 1; i <= n; i++) {
			type = (InternalStateType)internalVarsToExport.at(i);
			InternalStateValueType iType = giveInternalStateValueType(type);
			this->exportIntVarAs(type, iType, stream, tStep);
		}
	}

	void
		NodalRecoveryModule::exportIntVarAs(InternalStateType valID, InternalStateValueType type, FILE *stream, TimeStep *tStep)
	{
		Domain *d = emodel->giveDomain(1);
		int ireg;
		int nnodes = d->giveNumberOfDofManagers(), inode;
		int j, jsize;
		FloatArray iVal(3);
		FloatMatrix t(3, 3);
		const FloatArray *val = NULL;

		//this->giveSmoother();

		int nindx = giveInternalStateTypeSize(type);

		for (int indx = 1; indx <= nindx; indx++) {
			// print header
			if (type == ISVT_SCALAR) {
				fprintf(stream, "SCALARS %s double 1\n", __InternalStateTypeToString(valID));
			}
			else if (type == ISVT_VECTOR) {
				fprintf(stream, "VECTORS %s double\n", __InternalStateTypeToString(valID));
			}
			else if ((type == ISVT_TENSOR_S3) || (type == ISVT_TENSOR_S3E)) {
				fprintf(stream, "TENSORS %s double\n", __InternalStateTypeToString(valID));
			}
			else if (type == ISVT_TENSOR_G) {
				fprintf(stream, "SCALARS %s_%d double 1\n", __InternalStateTypeToString(valID), indx);
			}
			else {
				fprintf(stderr, "VTKExportModule :: exportIntVarAs: unsupported variable type %s\n", __InternalStateTypeToString(valID));
			}

			if ((type == ISVT_SCALAR) || (type == ISVT_TENSOR_G)) {
				fprintf(stream, "LOOKUP_TABLE default\n");
			}

			if (!((valID == IST_DisplacementVector) || (valID == IST_MaterialInterfaceVal))) {
				this->smoother->recoverValues(*elemSet, valID, tStep);
			}

			IntArray regionNodalNumbers(nnodes);
			int regionDofMans = 0, offset = 0;
			ireg = -1;
			int defaultSize = 0;

			this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, offset, d, ireg, 1);
			if (!((valID == IST_DisplacementVector) || (valID == IST_MaterialInterfaceVal))) {
				// assemble local->global map
				defaultSize = giveInternalStateTypeSize(type);
			}
			else {
				regionDofMans = nnodes;
			}

			for (inode = 1; inode <= regionDofMans; inode++) {
				if (valID == IST_DisplacementVector) {
					iVal.resize(3);
					val = &iVal;
					for (j = 1; j <= 3; j++) {
						iVal.at(j) = d->giveNode(regionNodalNumbers.at(inode))->giveUpdatedCoordinate(j, tStep, 1.0) -
							d->giveNode(regionNodalNumbers.at(inode))->giveCoordinate(j);
					}
				}
				else if (valID == IST_MaterialInterfaceVal) {
					MaterialInterface *mi = emodel->giveMaterialInterface(1);
					if (mi) {
						iVal.resize(1);
						val = &iVal;
						iVal.at(1) = mi->giveNodalScalarRepresentation(regionNodalNumbers.at(inode));
					}
				}
				else {
					this->smoother->giveNodalVector(val, regionNodalNumbers.at(inode));
				}

				if (val == NULL) {
					iVal.resize(defaultSize);
					iVal.zero();
					val = &iVal;
					//OOFEM_ERROR("internal error: invalid dofman data");
				}
			}


			if (type == ISVT_SCALAR) {
				if (val->giveSize()) {
					fprintf(stream, "%e ", val->at(1));
				}
				else {
					fprintf(stream, "%e ", 0.0);
				}
			}
			else if (type == ISVT_VECTOR) {
				jsize = min(3, val->giveSize());
				for (j = 1; j <= jsize; j++) {
					fprintf(stream, "%e ", val->at(j));
				}

				for (j = jsize + 1; j <= 3; j++) {
					fprintf(stream, "0.0 ");
				}
			}
			else if (type == ISVT_TENSOR_S3) {
				t.zero();
				for (int ii = 1; ii <= 6; ii++) {
					if (ii == 1) {
						t.at(1, 1) = val->at(ii);
					}
					else if (ii == 2) {
						t.at(2, 2) = val->at(ii);
					}
					else if (ii == 3) {
						t.at(3, 3) = val->at(ii);
					}
					else if (ii == 4) {
						t.at(2, 3) = val->at(ii);
						t.at(3, 2) = val->at(ii);
					}
					else if (ii == 5) {
						t.at(1, 3) = val->at(ii);
						t.at(3, 1) = val->at(ii);
					}
					else if (ii == 6) {
						t.at(1, 2) = val->at(ii);
						t.at(2, 1) = val->at(ii);
					}
				}

				for (int ii = 1; ii <= 3; ii++) {
					for (int jj = 1; jj <= 3; jj++) {
						fprintf(stream, "%e ", t.at(ii, jj));
					}

					fprintf(stream, "\n");
				}
			}
			else if (type == ISVT_TENSOR_G) { // export general tensor values as scalars
				fprintf(stream, "%e ", val->at(indx));
			}

			fprintf(stream, "\n");
		}
	}


	int
		NodalRecoveryModule::initRegionNodeNumbering(IntArray &regionNodalNumbers, int &regionDofMans,
		int offset, Domain *domain, int reg, int mode)
	{
		// if mode == 0 then regionNodalNumbers is array with mapping from global numbering to local region numbering.
		// The i-th value contains the corresponding local region number (or zero, if global numbar is not in region).

		// if mode == 1 then regionNodalNumbers is array with mapping from local to global numbering.
		// The i-th value contains the corresponding global node number.


		int nnodes = domain->giveNumberOfDofManagers();
		int elemNodes;
		int currOffset = offset + 1;

		regionNodalNumbers.resize(nnodes);
		regionNodalNumbers.zero();
		regionDofMans = 0;

		for (auto &element : domain->giveElements()) {

			if (element->giveParallelMode() != Element_local) {
				continue;
			}

			elemNodes = element->giveNumberOfNodes();
			//  elemSides = element->giveNumberOfSides();

			// determine local region node numbering
			for (int elementNode = 1; elementNode <= elemNodes; elementNode++) {
				int node = element->giveNode(elementNode)->giveNumber();
				if (regionNodalNumbers.at(node) == 0) { // assign new number
					/* mark for assignement. This is done later, as it allows to preserve
					* natural node numbering.
					*/
					regionNodalNumbers.at(node) = 1;
					regionDofMans++;
				}
			}
		}

		if (mode == 1) {
			IntArray answer(nnodes);
			for (int i = 1; i <= nnodes; i++) {
				if (regionNodalNumbers.at(i)) {
					regionNodalNumbers.at(i) = currOffset++;
					answer.at(regionNodalNumbers.at(i)) = i;
				}
			}

			regionNodalNumbers = answer;
		}
		else {
			for (int i = 1; i <= nnodes; i++) {
				if (regionNodalNumbers.at(i)) {
					regionNodalNumbers.at(i) = currOffset++;
				}
			}
		}

		return 1;
	}


	void
		NodalRecoveryModule::doOutput(TimeStep *tStep, bool forcedOutput)
	{
		if (!(testTimeStepOutput(tStep) || forcedOutput)) {
			return;
		}

		// gather stuff
		Domain *d = emodel->giveDomain(1);
		for (auto &elem : d->giveElements()) {
			if (this->checkValidType(elem->giveClassName())) {

			}

		}

		if (this->isRespSpec && tStep->giveIntrinsicTime()!=0){
			// square and save
			BeamDisplacementsList.push_back(BeamDisplacements);
			BeamForcesList.push_back(BeamForces);

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

			BeamDisplacements.clear();
			BeamForces.clear();

			fflush(this->stream);
		}
	}

	void NodalRecoveryModule::populateElResults(map<int, map<double, FloatArray>> &answer, map<int, map<double, FloatArray>> &src)
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

	void NodalRecoveryModule::addMultiply(map<int, map<double, FloatArray>> &answer, map<int, map<double, FloatArray>> &src, map<int, map<double, FloatArray>> &src2, double fact)
	{
		if (answer.size() == 0) {
			populateElResults(answer, src);
		}

		// awful iteration
		map<int, map<double, FloatArray>>::iterator destElem_it = answer.begin();
		map<int, map<double, FloatArray>>::iterator srcElem_it = src.begin();
		map<int, map<double, FloatArray>>::iterator srcElem_it2 = src.begin();
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

	void NodalRecoveryModule::calcRoot(map<int, map<double, FloatArray>> &answer)
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

	void NodalRecoveryModule::SRSS(){
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

	void NodalRecoveryModule::CQC(){
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
		NodalRecoveryModule::initialize()
	{
		Domain *d = emodel->giveDomain(1);

		if (this->smoother == NULL) {
			this->smoother = classFactory.createNodalRecoveryModel(this->stype, d);
		}

		// create a new set containing all elements
		elemSet = new Set(0, d);
		elemSet->addAllElements();

		for (int i = 1; i <= internalVarsToExport.giveSize(); i++) {
			InternalStateType type = (InternalStateType)internalVarsToExport.at(i);
			valueTypesStr.push_back(__InternalStateTypeToString(type));
		}

		IntArray map(d->giveNumberOfDofManagers());

		// asemble local->global region map
		this->initRegionNodeNumbering(map, regionDofMans, 0, d, -1, 1);

		string fileName = emodel->giveOutputBaseFileName() + ".nrm";
		if ((this->stream = fopen(fileName.c_str(), "w")) == NULL) {
			OOFEM_ERROR("failed to open file %s", fileName.c_str());
		}
		// ";" as separator
		//fprintf(this->stream, "#Time;ElemNo;ValNo;N_x;T_y;T_z;M_x;M_y;M_z;dx;dy;dz;rx;ry;rz;");
		//for ( int var: this->ists ) {
		//    fprintf(this->stream, "%s    ", __InternalStateTypeToString( ( InternalStateType ) var) );
		//}
		fprintf(this->stream, "\n");
		fflush(this->stream);
	}

	void
		NodalRecoveryModule::terminate()
	{
		fclose(this->stream);
	}
} // end namespace oofem