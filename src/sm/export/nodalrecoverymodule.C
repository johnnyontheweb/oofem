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

#include "nodalrecoverymodule.h"
#include "timestep.h"
#include "element.h"
#include "gausspoint.h"
#include "engngm.h"
#include "node.h"
#include "dof.h"
#include "material.h"
#include "classfactory.h"
#include "generalboundarycondition.h"
#include "inputrecord.h"
#include "../sm/EngineeringModels/responsespectrum.h"
#include "cltypes.h"
#include "materialinterface.h"
#include "nodalaveragingrecoverymodel.h"
#include <math.h>
#include <vector>

#ifdef MEMSTR
#include <io.h>
#include <fcntl.h>
#endif

using namespace std;

namespace oofem {
	REGISTER_ExportModule(NodalRecoveryModule)

	IntArray NodalRecoveryModule::redToFull = {
		1, 5, 9, 8, 7, 4, 6, 3, 2
	};

	NodalRecoveryModule::NodalRecoveryModule(int n, EngngModel *e) : ExportModule(n, e) { }

	NodalRecoveryModule :: ~NodalRecoveryModule()
	{
		delete elemSet;
	}

	IRResultType
		NodalRecoveryModule::initializeFrom(InputRecord *ir)
	{
		IRResultType result;                 // Required by IR_GIVE_FIELD macro
		int val;

		val = NodalRecoveryModel::NRM_NodalAveraging;
		IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_NodalRecoveryModule_stype);
		stype = (NodalRecoveryModel::NodalRecoveryModelType) val;

		IR_GIVE_OPTIONAL_FIELD(ir, internalVarsToExport, _IFT_NodalRecoveryModule_vars);
		if (internalVarsToExport.giveSize() == 0) OOFEM_ERROR("NodalRecoveryModule - No response types defined");
		for (int i = 1; i<=internalVarsToExport.giveSize(); i++) {
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
		NodalRecoveryModule::exportIntVars(TimeStep *tStep)
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
			map< int, FloatArray >  nodVec;
			type = (InternalStateType)internalVarsToExport.at(i);
			InternalStateValueType iType = giveInternalStateValueType(type);
			this->exportIntVarAs(nodVec, type, iType, tStep);
			nodalValues[i] = nodVec;
		}
	}

	void
		NodalRecoveryModule::exportIntVarAs(map< int, FloatArray > &answer, InternalStateType valID, InternalStateValueType type, TimeStep *tStep)
	{
		Domain *d = emodel->giveDomain(1);
		int ireg;
		int nnodes = d->giveNumberOfDofManagers(), inode;
		int j;
		FloatArray iVal(3);
		FloatArray t(9);
		const FloatArray *val = NULL;

		//this->giveSmoother();

		if (!((valID == IST_DisplacementVector) || (valID == IST_MaterialInterfaceVal))) {
			this->smoother->recoverValues(*this->elemSet, valID, tStep);
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
			//regionDofMans = nnodes;
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

#ifdef DEBUG
			if (type == ISVT_SCALAR || type == ISVT_VECTOR || type == ISVT_TENSOR_S3 || type == ISVT_TENSOR_G) {
				FloatArray t;
				this->makeFullTensorForm(t, *val, type);
				answer[regionNodalNumbers.at(inode)] = t;

			}
			else {
#endif
				answer[regionNodalNumbers.at(inode)] = *val;
#ifdef DEBUG
			}
#endif

		}



	}

	void
		NodalRecoveryModule::makeFullTensorForm(FloatArray &answer, const FloatArray &reducedForm, InternalStateValueType vtype)
	{
		answer.resize(9);
		answer.zero();

		for (int i = 1; i <= reducedForm.giveSize(); i++) {
			answer.at(redToFull.at(i)) = reducedForm.at(i);
		}

		if (vtype == ISVT_TENSOR_S3E) {
			answer.at(4) *= 0.5;
			answer.at(7) *= 0.5;
			answer.at(8) *= 0.5;
		}

		// Symmetrize if needed
		if (vtype != ISVT_TENSOR_G) {
			answer.at(2) = answer.at(4);
			answer.at(3) = answer.at(7);
			answer.at(6) = answer.at(8);
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

		if ((!this->isRespSpec && tStep->giveIntrinsicTime() == 0) || (tStep->giveIntrinsicTime() != 0)){
			// gather stuff
			//Domain *d = emodel->giveDomain(1);
			//for (auto &elem : d->giveElements()) {
				//if (this->checkValidType(elem->giveClassName())) { // ALL?
					this->exportIntVars(tStep); // save to map.
				//}
			//}
		}

		// is this responseSpectrum? then pile up stuff to combine later
		if (this->isRespSpec && tStep->giveIntrinsicTime()!=0){
			// square and save
			combNodalValuesList.push_back(nodalValues);
		} else {

			// is RespSpec and magic timeStep? then time to combine results
			if (this->isRespSpec && tStep->giveIntrinsicTime() == 0) {

				if (rs->giveComboType() == RSC_SRSS) {
					this->SRSS();
				} else {
					this->CQC();
				}

				nodalValues = combNodalValues;

				combNodalValues.clear();
			}


			// export all the stuff!
			Domain *d = emodel->giveDomain(1);
			double curTime = tStep->giveTargetTime();
			map<int, map<int, FloatArray>>::iterator values_it = nodalValues.begin();
			list<string>::iterator names_it = valueTypesStr.begin();
			for (;
				values_it != nodalValues.end();
				++values_it, ++names_it)
			{
				map< int, FloatArray > &nodes = values_it->second;
				const string resp = *names_it;

				map<int, FloatArray >::iterator nodes_it = nodes.begin();

				for (;
					nodes_it != nodes.end();
					++nodes_it)
				{
					int node = nodes_it->first;
					Node* nd = d->giveNode(node);
					FloatArray resps = nodes_it->second;
					double chRes = resps.computeNorm();
					if (resps.giveSize() && chRes != 0 && isnan(chRes) == false) {
						fprintf(this->stream, "%10.3e;%s;%d;", curTime, resp.c_str(), nd->giveLabel());

						for (auto &val : resps) {
							fprintf(this->stream, "%10.3e;", val);
						}
						fprintf(this->stream, "\n");
					}
				}

			}

			// clean up
			nodalValues.clear();
			fflush(this->stream);
		}
	}

	void NodalRecoveryModule::populateElResults(map< int, map< int, FloatArray > > &answer, map< int, map< int, FloatArray > > &src)
	{

		map< int, map< int, FloatArray > >::iterator srcValue_it = src.begin();
		for (; srcValue_it != src.end(); ++srcValue_it)
		{
			map< int, FloatArray > *destNRespMap = new map< int, FloatArray >;
			map< int, FloatArray > &srcNRespMap = srcValue_it->second;

			map< int, FloatArray >::iterator srcNRespMap_it = srcNRespMap.begin();
			for (; srcNRespMap_it != srcNRespMap.end(); ++srcNRespMap_it)
			{
				FloatArray &srcRespArray = srcNRespMap_it->second;
				FloatArray *destRespArray = new FloatArray(srcRespArray.giveSize());

				destNRespMap->operator[](srcNRespMap_it->first) = *destRespArray;
			}

			answer[srcValue_it->first] = *destNRespMap;
		}
	}

	void NodalRecoveryModule::addMultiply(std::map< int, std::map< int, FloatArray > > &answer, std::map< int, std::map< int, FloatArray > > &src, std::map< int, std::map< int, FloatArray > > &src2, double fact)
	{
		if (answer.size() == 0) {
			populateElResults(answer, src);
		}

		// awful iteration
		std::map< int, std::map< int, FloatArray > >::iterator destValue_it = answer.begin();
		std::map< int, std::map< int, FloatArray > >::iterator srcValue_it = src.begin();
		std::map< int, std::map< int, FloatArray > >::iterator srcValue_it2 = src2.begin();
		for (; destValue_it != answer.end(); ++destValue_it, ++srcValue_it, ++srcValue_it2)
		{
			map<int, FloatArray> &destRespMap = destValue_it->second;
			map<int, FloatArray> &srcRespMap = srcValue_it->second;
			map<int, FloatArray> &srcRespMap2 = srcValue_it2->second;

			map<int, FloatArray>::iterator destRespMap_it = destRespMap.begin();
			map<int, FloatArray>::iterator srcRespMap_it = srcRespMap.begin();
			map<int, FloatArray>::iterator srcRespMap_it2 = srcRespMap2.begin();
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

	void NodalRecoveryModule::calcRoot(std::map< int, std::map< int, FloatArray > > &answer)
	{
		// another awful iteration
		std::map< int, std::map< int, FloatArray > >::iterator destValue_it = answer.begin();
		for (; destValue_it != answer.end(); ++destValue_it)
		{
			map<int, FloatArray> &destRespMap = destValue_it->second;

			map<int, FloatArray>::iterator destRespMap_it = destRespMap.begin();
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
		list<map<int, map<int, FloatArray>>>::iterator values_it = combNodalValuesList.begin();
		for (; values_it != combNodalValuesList.end(); ++values_it)
		{
			addMultiply(combNodalValues, *values_it, *values_it);
		}
		calcRoot(combNodalValues);

	}

	void NodalRecoveryModule::CQC(){
		FloatMatrix rhos;
		rs->giveRhos(rhos);

		list<map<int, map<int, FloatArray>>>::iterator disps_it = combNodalValuesList.begin();

		for (int i = 1; disps_it != combNodalValuesList.end(); ++disps_it, i++)
		{
			list<map<int, map<int, FloatArray>>>::iterator disps_it2 = combNodalValuesList.begin();
			for (int j = 1; disps_it2 != combNodalValuesList.end(); ++disps_it2, j++)
			{
				addMultiply(combNodalValues, *disps_it, *disps_it2, rhos.at(i, j));
			}
		}
		calcRoot(combNodalValues);
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

		// save result type strings
		for (int i = 1; i <= internalVarsToExport.giveSize(); i++) {
			InternalStateType type = (InternalStateType)internalVarsToExport.at(i);
			valueTypesStr.push_back(__InternalStateTypeToString(type));
		}

#ifdef MEMSTR
		this->stream = nullptr;
		FILE *source = classFactory.giveMemoryStream("nrm");
		int sourceFD = _open_osfhandle((intptr_t)source, _O_APPEND);
		if (sourceFD != -1) {
			this->stream = _fdopen(sourceFD, "a");
		}
		if (!(this->stream)) {  // if not, write to file
#endif
			string fileName = emodel->giveOutputBaseFileName() + ".nrm";
			if ((this->stream = fopen(fileName.c_str(), "w")) == NULL) {
				OOFEM_ERROR("failed to open file %s", fileName.c_str());
			}
#ifdef MEMSTR
		}
#endif
		// ";" as separator
		fprintf(this->stream, "#Time;ISType;NodeNo;vars;");
		//for ( int var: this->ists ) {
		//    fprintf(this->stream, "%s    ", __InternalStateTypeToString( ( InternalStateType ) var) );
		//}
		fprintf(this->stream, "\n");
		fflush(this->stream);
	}

	void
		NodalRecoveryModule::terminate()
	{
#ifndef MEMSTR
		fclose(this->stream);
#endif
	}
} // end namespace oofem
