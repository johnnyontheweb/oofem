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

#include "../sm/EngineeringModels/responsespectrum.h"
#include "timestep.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "exportmodulemanager.h"
#include "verbose.h"
#include "classfactory.h"
#include "datastream.h"
#include "geneigvalsolvertype.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dofmanager.h"
#include "dof.h"
#include "domain.h"
#include "element.h"
#include "node.h"
#include "unknownnumberingscheme.h"
#include "function.h"
#include "activebc.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_EngngModel(ResponseSpectrum);

NumericalMethod *ResponseSpectrum :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
        nMethod.reset( classFactory.createGeneralizedEigenValueSolver(solverType, this->giveDomain(1), this) );
        if ( !nMethod ) {
            OOFEM_ERROR("solver creation failed");
        }
    }

    return nMethod.get();
}

IRResultType
ResponseSpectrum::initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    //EngngModel::instanciateFrom (ir);

	IR_GIVE_FIELD(ir, numberOfRequiredEigenValues, _IFT_ResponseSpectrum_nroot);

    // numberOfSteps set artificially to numberOfRequiredEigenValues
    // in order to allow
    // use restoreContext function for different eigenValues
    // numberOfSteps = numberOfRequiredEigenValues;
    numberOfSteps = 1;

	IR_GIVE_FIELD(ir, rtolv, _IFT_ResponseSpectrum_rtolv);
    if ( rtolv < 1.e-12 ) {
        rtolv =  1.e-12;
    }

    if ( rtolv > 0.01 ) {
        rtolv =  0.01;
    }

    int val = 0;
	IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_ResponseSpectrum_stype);
    solverType = ( GenEigvalSolverType ) val;

    val = 0; //Default Skyline
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

	IR_GIVE_FIELD(ir, val, _IFT_ResponseSpectrum_func);
	func = (int)val; // we'll check in postInitialize whether this id exists or not

	IR_GIVE_FIELD(ir, val, _IFT_ResponseSpectrum_dir);
	dir = (int)val;
	if (dir < D_u && dir > D_w){
		OOFEM_ERROR("Invalid direction. Currently 1 to 3 are supported");
	}

    return IRRT_OK;
}


void ResponseSpectrum::postInitialize()
{
	EngngModel::postInitialize();

	// we check whether the spectrum function exists or not.
	Domain *d = this->giveDomain(1);
	Function *f = d->giveFunction(func);

	if (f == NULL) OOFEM_ERROR("Invalid function given");
}


double ResponseSpectrum::giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, eigenvalue.
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
#ifdef DEBUG
    if ( eq == 0 ) {
        OOFEM_ERROR("invalid equation number");
    }
#endif

    switch ( mode ) {
    case VM_Total:  // EigenVector
    case VM_Incremental:
        return eigVec.at( eq, ( int ) tStep->giveTargetTime() );

    default:
        OOFEM_ERROR("Unknown is of undefined type for this problem");
    }

    return 0.;
}


TimeStep *ResponseSpectrum::giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    StateCounterType counter = 1;

    if ( currentStep ) {
        istep = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = std :: move(currentStep);
    currentStep.reset( new TimeStep(istep, this, 1, ( double ) istep, 0., counter) );

    return currentStep.get();
}


// gets the Spectral acceleration from the function given in input
double ResponseSpectrum::calcSpectrumOrdinate(double period)
{
	Function *f = this->giveDomain(1)->giveFunction(this->func);
	return f->evaluateAtTime(period);
}


void ResponseSpectrum::solveYourselfAt(TimeStep *tStep)
{
    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling stiffness and mass matrices\n");
#endif

    if ( tStep->giveNumber() == 1 ) {
        //
        // first step  assemble stiffness Matrix
        //

        stiffnessMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

        massMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        massMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

        this->assemble( *stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        this->assemble( *massMatrix, tStep, MassMatrixAssembler(),
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //
        // create resulting objects eigVec and eigVal
        //
        eigVec.resize(this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ), numberOfRequiredEigenValues);
        eigVec.zero();
        eigVal.resize(numberOfRequiredEigenValues);
        eigVal.zero();
    }

    //
    // set-up numerical model
    //
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );

    //
    // call numerical model to solve arised problem
    //
#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving ...\n");
#endif

    nMethod->solve(*stiffnessMatrix, *massMatrix, eigVal, eigVec, rtolv, numberOfRequiredEigenValues);

		FloatMatrix *unitDisp = new FloatMatrix();
	FloatArray *tempCol = new FloatArray();
	FloatArray *tempCol2 = new FloatArray();

	

	Domain *domain = this->giveDomain(1);
	IntArray dofIDArry, loc;
	dofIDArry = domain->giveDefaultNodeDofIDArry();
	int nelem = domain->giveNumberOfElements();

	// matrix and array initialization
	totMass.resize(6);		// 3 are the translational dofs - enough for the moment
	centroid.resize(3);
	partFact.resize(numberOfRequiredEigenValues, 6);
	massPart.resize(numberOfRequiredEigenValues, 6);
	unitDisp->resize(this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()), 6);
	tempCol->resize(this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()));
	tempCol2->resize(this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()));
	periods.resize(numberOfRequiredEigenValues);

	// mass normalization
	for (int i = 1; i <= numberOfRequiredEigenValues; i++)
	{
		tempCol->beColumnOf(eigVec, i);
		massMatrix->times(*tempCol, *tempCol2);
		double m = tempCol->dotProduct(*tempCol2);
		if (m != 0.0) m = 1 / sqrt(m);
		tempCol->times(m);
		eigVec.setColumn(*tempCol, i);
		periods.at(i)= 2 * M_PI / sqrt(eigVal.at(i));
	}
	// eigVec has been normalized

	IntArray masterDofIDs, nodalArray, ids;
	IntArray locationArray;
	IntArray *dofIdArray = new IntArray();
	dofIdArray->clear();

	FloatMatrix tempMat, tempMat2;
	FloatArray tempCoord, coordArray;
	tempMat.resize(this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()), 6);
	tempMat2.resize(this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()), 6);

	//
	// create unit displacement vectors
	//
	// first from nodes themselves
	for (std::unique_ptr<DofManager> &node : domain->giveDofManagers()) {
		//node->giveLocationArray(dofIDArry, loc, EModelDefaultEquationNumbering());
		if (!node->giveNumberOfDofs()) continue;

		for (int dType = D_u; dType <= D_w; dType++)
		{
			auto myDof = node->findDofWithDofId((DofIDItem)dType);
			if (myDof == node->end()) {
				//OOFEM_ERROR("incompatible dof (%d) requested", dType);
				continue;
			}

			int eqN = EModelDefaultEquationNumbering().giveDofEquationNumber(*myDof);

			// save unit displacement and coordinate.
			// TODO consider the fact that nodes may have own UCS.
			if (eqN)
			{
				unitDisp->at(eqN, dType) = 1.0;
				tempMat2.at(eqN, dType) = node->giveCoordinate(dType);
			}
		}
	}  // end of search among nodes

	// then from internaldof managers
	for (int ielem = 1; ielem <= nelem; ielem++) {
		Element *element = domain->giveElement(ielem);

		// the following may be simplified.
		// retrieve internal dof managers and location array
		for (int i = 1; i <= element->giveNumberOfInternalDofManagers(); i++) {
			DofManager *intDofMan = element->giveInternalDofManager(i);

			if (!intDofMan) continue; // you may never know...

			locationArray.clear();
			tempCoord.clear();

			element->giveInternalDofManDofIDMask(i, ids);
			intDofMan->giveLocationArray(ids, nodalArray, EModelDefaultEquationNumbering());
			locationArray.followedBy(nodalArray);

			intDofMan->giveMasterDofIDArray(ids, masterDofIDs);
			dofIdArray->followedBy(masterDofIDs);

			coordArray.resize(masterDofIDs.giveSize());
			int c = 1;
			for (int dof : masterDofIDs)
			{
				if (intDofMan->hasCoordinates()) {
					coordArray.at(c) = intDofMan->giveCoordinate(dof);
				}
				else
				{
					// get the number. ghostNode FEMComponent number is changed on purpose.
					int tempN = intDofMan->giveNumber();
					//element->giveDofManager(tempN)->requiresTransformation()
					coordArray.at(c) = element->giveDofManager(tempN)->giveCoordinate(dof);
				}
				c++; // Increment the counter in any case. The position index must match the index in masterDofIDs. We'll sort the dofs needed hereafter
			}
			tempCoord.append(coordArray);

			int partialDofCount = locationArray.giveSize();
			if (partialDofCount) {
				// search for our dofs in there
				for (int myDofIndex = 1; myDofIndex <= partialDofCount; myDofIndex++)
				{
					int dType = dofIdArray->at(myDofIndex);
					int eqN = locationArray.at(myDofIndex);

					if ((dType >= D_u) && (dType <= D_w) && eqN) {
						// save unit displacement and coordinate

						unitDisp->at(eqN, dType) = 1.0;
						tempMat2.at(eqN, dType) = tempCoord.at(myDofIndex);
					}
				}
			}
		}
	}  // end of search among internal dof managers
	// end of creation of translational unit displacement vectors


	for (int i = 1; i <= 3; i++)
	{
		tempCol->beColumnOf(*unitDisp, i);
		massMatrix->times(*tempCol, *tempCol2);  // now tempCol2 has only the masses pertaining the i-th direction
		tempMat.setColumn(*tempCol2, i);
		totMass.at(i) = tempCol->dotProduct(*tempCol2);  // total mass for direction i-th direction
		tempCol->beColumnOf(tempMat2, i);	// fetch coordinates in i-th direction
		if (totMass.at(i) != 0.0) centroid.at(i) = tempCol->dotProduct(*tempCol2) / totMass.at(i);  // dot multiply to get first moment, then divide by total mass in i-th direction to get i-th coordinate of the centroid
	}


	// we have the centroid. we can now calculate rotational components. first from nodes.
	for (std::unique_ptr<DofManager> &node : domain->giveDofManagers()) {
		//node->giveLocationArray(dofIDArry, loc, EModelDefaultEquationNumbering());
		if (!node->giveNumberOfDofs()) continue;

		FloatArray* nodeCoords = node->giveCoordinates();
		if (nodeCoords){
			FloatArray vk(3);
			IntArray eq(3);

			// TODO consider own UCS if present
			for (int dType = D_u; dType <= D_w; dType++)
			{
				auto myDof = node->findDofWithDofId((DofIDItem)dType);
				if (myDof == node->end()){
					vk.at(dType) = 0.0;
					eq.at(dType) = 0;
					continue;
				}
				vk.at(dType) = node->giveCoordinate(dType) - centroid.at(dType);
				eq.at(dType) = EModelDefaultEquationNumbering().giveDofEquationNumber(*myDof);
			}

			// set mixed contribution due to rotation about centroid
			if (eq.at(1)){
				unitDisp->at(eq.at(1), 5) = vk.at(3);
				unitDisp->at(eq.at(1), 6) = -vk.at(2);
			}

			if (eq.at(2)){
				unitDisp->at(eq.at(2), 4) = -vk.at(3);
				unitDisp->at(eq.at(2), 6) = vk.at(1);
			}

			if (eq.at(3)){
				unitDisp->at(eq.at(3), 4) = vk.at(2);
				unitDisp->at(eq.at(3), 5) = -vk.at(1);
			}

			// set pure rotational contribution
			for (int dType = R_u; dType <= R_w; dType++)
			{
				auto myDof = node->findDofWithDofId((DofIDItem)dType);
				if (myDof == node->end()) {
					//OOFEM_ERROR("incompatible dof (%d) requested", dType);
					continue;
				}

				int eqN = EModelDefaultEquationNumbering().giveDofEquationNumber(*myDof);

				// save unit displacement and coordinate
				// TODO consider own UCS if present
				if (eqN)
				{
					unitDisp->at(eqN, dType) = 1.0;
					tempMat2.at(eqN, dType) = node->giveCoordinate(dType);
				}

			}
		}
	} // end of search among nodes

	// then from internaldof managers
	for (int ielem = 1; ielem <= nelem; ielem++) {
		Element *element = domain->giveElement(ielem);

		// the following may be simplified
		//	retrieve internal dof managers and location array
		for (int i = 1; i <= element->giveNumberOfInternalDofManagers(); i++) {
			DofManager *intDofMan = element->giveInternalDofManager(i);

			if (!intDofMan) continue; // you may never know...

			locationArray.clear();
			tempCoord.clear();

			element->giveInternalDofManDofIDMask(i, ids);
			intDofMan->giveLocationArray(ids, nodalArray, EModelDefaultEquationNumbering());
			locationArray.followedBy(nodalArray);

			intDofMan->giveMasterDofIDArray(ids, masterDofIDs);
			dofIdArray->followedBy(masterDofIDs);

			coordArray.resize(masterDofIDs.giveSize());
			int c = 1;
			for (int dof : masterDofIDs)
			{
				if (intDofMan->hasCoordinates()) {
					coordArray.at(c) = intDofMan->giveCoordinate(dof);
				}
				else
				{
					// get the number. ghostNode FEMComponent number is changed on purpose.
					int tempN = intDofMan->giveNumber();
					//element->giveDofManager(tempN)->requiresTransformation()
					coordArray.at(c) = element->giveDofManager(tempN)->giveCoordinate(dof);
				}
				c++; // Increment the counter in any case. The position index must match the index in masterDofIDs. We'll sort the dofs needed hereafter
			}
			tempCoord.append(coordArray);

			// TODO should element CS be considered even with internal dof managers
			if (locationArray.giveSize()){

				FloatArray vk(3);
				IntArray eq(3);
				for (int dType = D_u; dType <= D_w; dType++)
				{
					int myDof = dofIdArray->findFirstIndexOf((DofIDItem)dType);
					if (myDof == 0){
						vk.at(dType) = 0.0;
						eq.at(dType) = 0;
						continue;
					}
					vk.at(dType) = tempCoord.at(myDof) - centroid.at(dType);
					eq.at(dType) = locationArray.at(myDof);
				}

				// set mixed contribution due to rotation about centroid
				if (eq.at(1)){
					unitDisp->at(eq.at(1), 5) = vk.at(3);
					unitDisp->at(eq.at(1), 6) = -vk.at(2);
				}

				if (eq.at(2)){
					unitDisp->at(eq.at(2), 4) = -vk.at(3);
					unitDisp->at(eq.at(2), 6) = vk.at(1);
				}

				if (eq.at(3)){
					unitDisp->at(eq.at(3), 4) = vk.at(2);
					unitDisp->at(eq.at(3), 5) = -vk.at(1);
				}


				// search for our dofs in there
				for (int dType = R_u; dType <= R_w; dType++)
				{
					int myDof = dofIdArray->findFirstIndexOf(dType);
					if (myDof == 0) continue;

					int eqN = locationArray.at(myDof);

					// save unit displacement and coordinate
					if (eqN)
					{
						unitDisp->at(eqN, dType) = 1.0;
						tempMat2.at(eqN, dType) = tempCoord.at(myDof);
					}
				}
			}
		}  // end of search among internal dof managers
		}
	// end of creation of translational unit displacement vectors

	for (int i = 4; i <= 6; i++)
	{
		tempCol->beColumnOf(*unitDisp, i);
		massMatrix->times(*tempCol, *tempCol2);  // now tempCol2 has only the masses pertaining the i-th direction
		tempMat.setColumn(*tempCol2, i);
		totMass.at(i) = tempCol->dotProduct(*tempCol2);  // total mass for direction i-th direction
		tempCol->beColumnOf(tempMat2, i);	// fetch coordinates in i-th direction
	}

	//
	// calculate participation factors and mass participation
	//

	// participation factors
	partFact.beTProductOf(eigVec, tempMat);

	// mass participation ratios
	for (int i = 1; i <= numberOfRequiredEigenValues; i++)
	{
		for (int j = 1; j <= 6; j++)
		{
			if (totMass.at(j) != 0.0)
			{
				massPart.at(i, j) = pow(partFact.at(i, j), 2) / totMass.at(j);
			}
			//else
			//{
			//	massPart.at(i, j) = 0.0;
			//}
		}
	}

	//
	// start creating loaded models
	//
#ifdef VERBOSE
	OOFEM_LOG_INFO("Starting creation of loaded models ...\n");
#endif


	// create linear solver
	std::unique_ptr< SparseLinearSystemNM > nLinMethod; // = NULL;

	if (!nLinMethod) {
		if (isParallel()) {
			if ((solverType == ST_Petsc) || (solverType == ST_Feti)) {
				nLinMethod.reset(classFactory.createSparseLinSolver(ST_Direct , this->giveDomain(1), this));
			}
		}
		else {
			nLinMethod.reset(classFactory.createSparseLinSolver(ST_Direct , this->giveDomain(1), this));
		}
		if (!nLinMethod) {
			OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
		}
	}
	
	for (int dN = 1; dN <= numberOfRequiredEigenValues; dN++)
	{
		OOFEM_LOG_INFO("Creation of loaded model %d...\n", dN);

		double sAcc = calcSpectrumOrdinate(periods.at(dN));

		loadVector.beColumnOf(eigVec, dN);
		loadVector *= (sAcc * partFact.at(dN, dir)); // scaled forces

		// this shouldn't be needed...
		//NM_Status s = nLinMethod->solve(*stiffnessMatrix, loadVector, dummyDisps);  // solve linear system
		//if (!(s & NM_Success)) {
		//	OOFEM_ERROR("No success in solving system.");
		//}

		FloatArray reactions;
		IntArray dofManMap, dofidMap, eqnMap;
		this->buildReactionTable(dofManMap, dofidMap, eqnMap, tStep, 1);

		// compute reaction forces
		this->computeReaction(reactions, tStep, 1);

		//std::stringstream outName;
		//FILE *outputContext;

		//outName << this->giveOutputBaseFileName().c_str() << dN;

		//Domain *d2 = domain->Clone();
		

		//if ((outputContext = fopen(outName.str().c_str(), "w")) == NULL) {
		//	OOFEM_ERROR("Can't open output file %s", outName.str().c_str());
		//}

		//FileDataStream outputContextStream(outputContext);

		//d2->saveContext(outputContextStream, CM_Definition);
	}

	//
	// zero matrix
	//
    stiffnessMatrix.reset(NULL);
    massMatrix.reset(NULL);

	// dispose the rest of the stuff
	delete unitDisp;
	delete tempCol;
	delete tempCol2;
	delete dofIdArray;
}


void ResponseSpectrum::buildReactionTable(IntArray &restrDofMans, IntArray &restrDofs,
	IntArray &eqn, TimeStep *tStep, int di)
{
	// determine number of restrained dofs
	Domain *domain = this->giveDomain(di);
	int numRestrDofs = this->giveNumberOfDomainEquations(di, EModelDefaultPrescribedEquationNumbering());
	int ndofMan = domain->giveNumberOfDofManagers();
	int rindex, count = 0;

	// initialize corresponding dofManagers and dofs for each restrained dof
	restrDofMans.resize(numRestrDofs);
	restrDofs.resize(numRestrDofs);
	eqn.resize(numRestrDofs);

	for (int i = 1; i <= ndofMan; i++) {
		DofManager *inode = domain->giveDofManager(i);
		for (Dof *jdof : *inode) {
			if (jdof->isPrimaryDof() && (jdof->hasBc(tStep))) { // skip slave dofs
				rindex = jdof->__givePrescribedEquationNumber();
				if (rindex) {
					count++;
					restrDofMans.at(count) = i;
					restrDofs.at(count) = jdof->giveDofID();
					eqn.at(count) = rindex;
				}
				else {
					// NullDof has no equation number and no prescribed equation number
					//_error("No prescribed equation number assigned to supported DOF");
				}
			}
		}
	}
	// Trim to size.
	restrDofMans.resizeWithValues(count);
	restrDofs.resizeWithValues(count);
	eqn.resizeWithValues(count);
}

void
ResponseSpectrum::computeReaction(FloatArray &answer, TimeStep *tStep, int di)
{
	//FloatArray contribution;

	answer.resize(this->giveNumberOfDomainEquations(di, EModelDefaultPrescribedEquationNumbering()));
	answer.zero();

	// Add internal forces
	this->assembleVector(answer, tStep, LastEquilibratedInternalForceAssembler(), VM_Total,
		EModelDefaultPrescribedEquationNumbering(), this->giveDomain(di));
	// Subtract external loading
	///@todo All engineering models should be using this (for consistency)
	//this->assembleVector( answer, tStep, ExternalForceAssembler(), VM_Total,
	//                    EModelDefaultPrescribedEquationNumbering(), this->giveDomain(di) );
	///@todo This method is overloaded in some functions, it needs to be generalized.
	//this->computeExternalLoadReactionContribution(contribution, tStep, di); changed to loadVector
	answer.subtract(loadVector);
	this->updateSharedDofManagers(answer, EModelDefaultPrescribedEquationNumbering(), ReactionExchangeTag);
}


void
ResponseSpectrum::updateInternalState(TimeStep *tStep)
{
	for (auto &domain : domainList) {
		if (requiresUnknownsDictionaryUpdate()) {
			for (auto &dman : domain->giveDofManagers()) {
				this->updateDofUnknownsDictionary(dman.get(), tStep);
			}
		}

		for (auto &bc : domain->giveBcs()) {
			ActiveBoundaryCondition *abc;

			if ((abc = dynamic_cast< ActiveBoundaryCondition * >(bc.get()))) {
				int ndman = abc->giveNumberOfInternalDofManagers();
				for (int j = 1; j <= ndman; j++) {
					this->updateDofUnknownsDictionary(abc->giveInternalDofManager(j), tStep);
				}
			}
		}

		if (true) { //internalVarUpdateStamp != tStep->giveSolutionStateCounter()
			for (auto &elem : domain->giveElements()) {
				elem->updateInternalState(tStep);
			}

			//internalVarUpdateStamp = tStep->giveSolutionStateCounter();
		}
	}
}


void ResponseSpectrum::updateYourself(TimeStep *tStep)
{
	this->updateInternalState(tStep);
	EngngModel::updateYourself(tStep);
}


void ResponseSpectrum::terminate(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    FILE *outputStream = this->giveOutputStream();

    // print loadcase header
    fprintf(outputStream, "\nOutput for time %.3e \n\n", 1.0);
    // print eigen values on output
    fprintf(outputStream, "\n\nEigen Values (Omega^2) are:\n-----------------\n");

    for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        fprintf( outputStream, "%15.8e ", eigVal.at(i) );
        if ( ( i % 5 ) == 0 ) {
            fprintf(outputStream, "\n");
        }
    }
	
	
	fprintf(outputStream, "\n\n\nCentroid Coordinates are:\n-----------------\n\tX\t|\tY\t|\tZ\n");
	for (int i = 1; i <= centroid.giveSize(); ++i) {
		fprintf(outputStream, "%10.3e ", centroid.at(i));
	}

	fprintf(outputStream, "\n");
	
	fprintf(outputStream, "\n\nParticipation Factors are:\n-----------------\n\tDx\t|\tDy\t|\tDz\t|\tRx\t|\tRy\t|\tRz\n");
	for (int i = 1; i <= partFact.giveNumberOfRows(); ++i) {
		for (int j = 1; j <= partFact.giveNumberOfColumns(); ++j) {
			fprintf(outputStream, "%10.3e ", partFact.at(i,j));
		}
		fprintf(outputStream, "\n");
	}

	fprintf(outputStream, "\n\nTotal Masses are:\n-----------------\n\tDx\t|\tDy\t|\tDz\t|\tRx\t|\tRy\t|\tRz\n");
	for (int i = 1; i <= totMass.giveSize(); ++i) {
		fprintf(outputStream, "%10.3e ", totMass.at(i));
	}

	fprintf(outputStream, "\n");

	fprintf(outputStream, "\n\nMass Ratios are:\n-----------------\n\tDx\t|\tDy\t|\tDz\t|\tRx\t|\tRy\t|\tRz\n");
	for (int i = 1; i <= massPart.giveNumberOfRows(); ++i) {
		for (int j = 1; j <= massPart.giveNumberOfColumns(); ++j) {
			fprintf(outputStream, "%10.3e ", massPart.at(i,j));
		}
		fprintf(outputStream, "\n");
	}

    fprintf(outputStream, "\n\n");

    for ( int i = 1; i <=  numberOfRequiredEigenValues; i++ ) {
        fprintf(outputStream, "\nOutput for eigen value no.  %.3e \n", ( double ) i);
        fprintf( outputStream,
                "Printing eigen vector no. %d, corresponding eigen value is %15.8e\n\n",
                i, eigVal.at(i) );
        tStep->setTime( ( double ) i ); // we use time as intrinsic eigen value index

        if ( this->requiresUnknownsDictionaryUpdate() ) {
            for ( auto &dman : domain->giveDofManagers() ) {
                this->updateDofUnknownsDictionary(dman.get(), tStep);
            }
        }

        for ( auto &dman : domain->giveDofManagers() ) {
            dman->updateYourself(tStep);
            dman->printOutputAt(outputStream, tStep);
        }
    }

    for ( int i = 1; i <=  numberOfRequiredEigenValues; i++ ) {
        // export using export manager
        tStep->setTime( ( double ) i ); // we use time as intrinsic eigen value index
        tStep->setNumber(i);
        exportModuleManager->doOutput(tStep);
    }
    fflush( this->giveOutputStream() );
    this->saveStepContext(tStep);
}


contextIOResultType ResponseSpectrum::saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves state variable - displacement vector
//
{
    int closeFlag = 0;
    contextIOResultType iores;
    FILE *file = NULL;

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, this->giveCurrentStep()->giveNumber(),
                                    this->giveCurrentStep()->giveVersion(), contextMode_write) ) {
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = EngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = eigVal.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( eigVec.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    } // ensure consistent records

    return CIO_OK;
}


contextIOResultType ResponseSpectrum::restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore state variable - displacement vector
//
{
    int closeFlag = 0;
    int activeVector = this->resolveCorrespondingEigenStepNumber(obj);
    int istep = 1, iversion = 0;
    contextIOResultType iores;
    FILE *file = NULL;

    if ( restoreFlag == 0 ) { // not restored before
        if ( stream == NULL ) {
            if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
                THROW_CIOERR(CIO_IOERR); // override
            }

            stream = new FileDataStream(file);
            closeFlag = 1;
        }

        // save element context

        if ( ( iores = EngngModel :: restoreContext(stream, mode, ( void * ) & istep) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = eigVal.restoreYourself(*stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = eigVec.restoreYourself(*stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( closeFlag ) {
            fclose(file);
            delete stream;
            stream = NULL;
        } // ensure consistent records
    }

    if ( activeVector > numberOfRequiredEigenValues ) {
        activeVector = numberOfRequiredEigenValues;
    }

    OOFEM_LOG_INFO( "Restoring - corresponding index is %d, EigenValue is %f\n", activeVector, eigVal.at(activeVector) );
    this->giveCurrentStep()->setTime( ( double ) activeVector );
    this->restoreFlag = 1;

    return CIO_OK;
}


int ResponseSpectrum::resolveCorrespondingEigenStepNumber(void *obj)
{
    //
    // returns corresponding eigen step number
    //
    if ( obj == NULL ) {
        return 1;
    }

    int *istep = ( int * ) obj;

    if ( * istep > numberOfRequiredEigenValues ) {
        return numberOfRequiredEigenValues;
    }

    if ( * istep <= 0 ) {
        return 1;
    }

    return * istep;
}

} // end namespace oofem