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

// please activate or de-activate next line
//#define LIN_STAB_COMPATIBILITY_MODE

#include "sm/EngineeringModels/linearstability.h"
#include "timestep.h"
#include "element.h"
#include "contextioerr.h"
#include "floatmatrix.h"
#include "verbose.h"
#include "floatarray.h"
#include "classfactory.h"
#include "datastream.h"
#include "exportmodulemanager.h"
#include "dofmanager.h"
#include "dof.h"
#include "unknownnumberingscheme.h"
#include "eigenvectorprimaryfield.h"


#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_EngngModel(LinearStability);

LinearStability :: LinearStability(int i, EngngModel *master) : StructuralEngngModel(i, master),
    numberOfRequiredEigenValues(1),
    rtolv(1e-6),
    solverType(GES_SubspaceIt)
{
    numberOfSteps = 1;
    ndomains = 1;
}


NumericalMethod *LinearStability :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
        nMethod = classFactory.createGeneralizedEigenValueSolver(solverType, this->giveDomain(1), this);
        if ( !nMethod ) {
            OOFEM_ERROR("solver creation failed");
        }
    }

    return nMethod.get();
}


SparseLinearSystemNM *LinearStability :: giveNumericalMethodForLinStaticProblem(TimeStep *tStep)
{
    if ( !nMethodLS ) {
        nMethodLS = classFactory.createSparseLinSolver(ST_Direct, this->giveDomain(1), this); ///@todo Support other solvers
        if ( !nMethodLS ) {
            OOFEM_ERROR("solver creation failed");
        }
    }

    return nMethodLS.get();
}


void
LinearStability :: initializeFrom(InputRecord &ir)
{
    //StructuralEngngModel::instanciateFrom(ir);
    IR_GIVE_FIELD(ir, numberOfRequiredEigenValues, _IFT_LinearStability_nroot);
    // numberOfSteps set artifficially to numberOfRequiredEigenValues
    // in order to allow
    // use restoreContext function for different eigenValues
    numberOfSteps = numberOfRequiredEigenValues;
    this->field = std::make_unique<EigenVectorPrimaryField>(this, 1, FT_Displacements, numberOfRequiredEigenValues + 1); // +1 for eq. solution

    IR_GIVE_FIELD(ir, rtolv, _IFT_LinearStability_rtolv);
    if ( rtolv < 1.e-12 ) {
        rtolv =  1.e-12;
    } else if ( rtolv > 0.01 ) {
        rtolv =  0.01;
    }

    int val = 1; // inverseit
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_LinearStability_stype);
    solverType = ( GenEigvalSolverType ) val;

	sparseMtrxType= SparseMtrxType::SMT_Skyline; //Default Skyline
	linStype = LinSystSolverType::ST_Direct;

    nMetaSteps = 0;

    suppressOutput = ir.hasField(_IFT_EngngModel_suppressOutput);

    if (suppressOutput) {
        printf("Suppressing output.\n");
    } else {
    nMetaSteps = 0;

<<<<<<< .mine
    return IRRT_OK;



=======
        fprintf(outputStream, "%s", PRG_HEADER);
        fprintf(outputStream, "\nStarting analysis on: %s\n", ctime(& this->startTime) );
        fprintf(outputStream, "%s\n", simulationDescription.c_str());
    }
>>>>>>> .theirs
}


int LinearStability :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
{
    return tStep->giveNumber() % (this->numberOfRequiredEigenValues + 1); // +1 for eq. solution 
    }

<<<<<<< .mine
    int activeVector = ( int ) tStep->giveTargetTime();
    switch ( mode ) {
	case VM_Incremental:
    case VM_Total: // EigenVector
        if ( activeVector ) {
            return eigVec.at(eq, activeVector);
        }
=======







>>>>>>> .theirs

<<<<<<< .mine
        return displacementVector.at(eq);
	case VM_Velocity:
		return 0.;
	case VM_Acceleration:
		return 0.;
    default:
        OOFEM_ERROR("Unknown is of undefined type for this problem");
=======
double LinearStability :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    return field->giveUnknownValue(dof, mode, tStep);




>>>>>>> .theirs
    }


TimeStep *LinearStability :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    StateCounterType counter = 1;

    if ( currentStep ) {
        istep =  currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = std :: move(currentStep);
    currentStep = std::make_unique<TimeStep>(istep, this, 1, 0., 0., counter);

    return currentStep.get();
}


void LinearStability :: solveYourself()
{
    this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);
    // update state according to new meta step
    this->giveNextStep();
    this->updateAttributes( this->giveCurrentMetaStep() );
    this->solveYourselfAt( this->giveCurrentStep() );
    this->terminate( this->giveCurrentStep() );
}


void LinearStability :: solveYourselfAt(TimeStep *tStep)
{
    tStep->setNumber(0);
    tStep->setTime(0.0);

    // creates system of governing eq's and solves them at given time step
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );
    this->giveNumericalMethodForLinStaticProblem(tStep);

    // first assemble problem at current time step
    if ( !stiffnessMatrix ) {
        //
        // first step - solve linear static problem
        //
<<<<<<< .mine
		stiffnessMatrix.reset(classFactory.createSparseMtrx(sparseMtrxType));
=======
        stiffnessMatrix = classFactory.createSparseMtrx(SMT_Skyline); ///@todo Don't hardcode skyline matrix only
>>>>>>> .theirs
        stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
    }

#ifndef LIN_STAB_COMPATIBILITY_MODE
    OOFEM_LOG_INFO("Assembling stiffness matrix\n");
    stiffnessMatrix->zero();
    this->assemble( *stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
                   EModelDefaultEquationNumbering(), this->giveDomain(1) );
#endif

<<<<<<< .mine
#ifdef VERBOSE
=======

>>>>>>> .theirs
    OOFEM_LOG_INFO("Assembling load\n");
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );
    FloatArray displacementVector(neq), loadVector(neq);

    // Internal forces first, negated;
    field->update(VM_Total, tStep, displacementVector, EModelDefaultEquationNumbering());
    this->assembleVector( loadVector, tStep, InternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    loadVector.negated();

    this->assembleVector( loadVector, tStep, ExternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->updateSharedDofManagers(loadVector, EModelDefaultEquationNumbering(), ReactionExchangeTag);

    OOFEM_LOG_INFO("Solving linear static problem\n");
    nMethodLS->solve(*stiffnessMatrix, loadVector, displacementVector);
    // Initial displacements are stored at position 0; this is a bit of a hack. In the future, a cleaner approach of handling fields could be suitable,
    // but currently, it all converges down to the same giveUnknownComponent, so this is the easisest approach.
    field->update(VM_Total, tStep, displacementVector, EModelDefaultEquationNumbering());
    // terminate linear static computation (necessary, in order to compute stresses in elements).
    // Recompute for updated state:
    this->assembleVector( loadVector, tStep, InternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->terminateLinStatic( tStep );

    // Normal forces already known, proceed with linear stability
    stiffnessMatrix->zero();
    if ( !initialStressMatrix ) {
        initialStressMatrix = stiffnessMatrix->clone();
    } else {
        initialStressMatrix->zero();
    }

    OOFEM_LOG_INFO("Assembling stiffness matrix\n");
    this->assemble( *stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
                   EModelDefaultEquationNumbering(), this->giveDomain(1) );

    OOFEM_LOG_INFO("Assembling initial stress matrix\n");
    this->assemble( *initialStressMatrix, tStep, InitialStressMatrixAssembler(),
                   EModelDefaultEquationNumbering(), this->giveDomain(1) );
    initialStressMatrix->times(-1.0);

	//FloatMatrix stiff, initstr;
	//stiffnessMatrix->writeToFile("stiff.txt");
	//initialStressMatrix->writeToFile("istr.txt");

    //  stiffnessMatrix->printYourself();
    //  initialStressMatrix->printYourself();

    FloatMatrix eigVec(neq, numberOfRequiredEigenValues);
    eigVal.resize(numberOfRequiredEigenValues);
    eigVal.zero();

    OOFEM_LOG_INFO("Solving ...\n");
#ifdef DEBUG
	//stiffnessMatrix->writeToFile("K.dat");
	//initialStressMatrix->writeToFile("M.dat");
#endif

    nMethod->solve(*stiffnessMatrix, *initialStressMatrix, eigVal, eigVec, rtolv, numberOfRequiredEigenValues);
<<<<<<< .mine

=======
    this->field->updateAll(eigVec, EModelDefaultEquationNumbering());
>>>>>>> .theirs
}


void LinearStability :: updateYourself(TimeStep *tStep)
{ }


void
LinearStability :: terminateLinStatic(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
<<<<<<< .mine
    // FILE *File = this->giveOutputStream();
    tStep->setTime(0.);
=======


>>>>>>> .theirs

<<<<<<< .mine
    //fprintf(File, "\nOutput for time %.3e \n\n", tStep->giveTargetTime() ); // FP: removed output
    //fprintf(File, "Linear static:\n\n"); // FP: removed output

    if ( requiresUnknownsDictionaryUpdate() ) {
=======




>>>>>>> .theirs
        for ( auto &dman : domain->giveDofManagers() ) {
        dman->updateYourself(tStep);
        //dman->printOutputAt(File, tStep); // FP: removed output
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Updated nodes ", domain->giveNumberOfDofManagers())
#  endif

    for ( auto &elem : domain->giveElements() ) {
        elem->updateInternalState(tStep);
        elem->updateYourself(tStep);
        //elem->printOutputAt(File, tStep); // FP: removed output
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Updated Elements ", domain->giveNumberOfElements())
#  endif
	// fprintf(File, "\n");
    /*
     * // save context if required
     * // default - save only if ALWAYS is set ( see cltypes.h )
     *
     * if ((domain->giveContextOutputMode() == COM_Always ) ||
     * (domain->giveContextOutputMode() == COM_Required )) {
     * this->saveContext(NULL);
     * }
     * else if (domain->giveContextOutputMode() == COM_UserDefined ) {
     * if (tStep->giveNumber()%domain->giveContextOutputStep() == 0)
     * this->saveContext(NULL);
     * }
     */

    //this->printReactionForces(tStep, 1);  // FP: removed output
}


void LinearStability :: terminate(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
<<<<<<< .mine
    FILE *outputStream = this->giveOutputStream();




=======
    // i = 0  represents the linear solution, which is followed by the eigen vectors starting at i = 1
    for ( int i = 0; i <= numberOfRequiredEigenValues; i++ ) {
        TimeStep step = *tStep;
        step.setTime( ( double ) i );
        step.setNumber(i);
>>>>>>> .theirs

<<<<<<< .mine
    // print eigen values on output
    fprintf(outputStream, "\nLinear Stability:");
    fprintf(outputStream, "\nEigen Values are:\n-----------------\n");
=======
        for ( auto &dman : domain->giveDofManagers() ) {
            dman->updateYourself(&step);
        }
>>>>>>> .theirs

<<<<<<< .mine
	for (int i = 1; i <= eigVal.giveSize(); i++) { // numberOfRequiredEigenValues
        fprintf( outputStream, "%15.8e ", eigVal.at(i) );
        if ( ( i % 5 ) == 0 ) {
            fprintf(outputStream, "\n");
=======
        exportModuleManager.doOutput(&step);



>>>>>>> .theirs
        }
    }

    fprintf(outputStream, "\n\n");

    int nnodes = domain->giveNumberOfDofManagers();

<<<<<<< .mine
	for (int i = 1; i <= eigVal.giveSize(); i++) { // numberOfRequiredEigenValues
        fprintf(outputStream, "\nOutput for eigen value no.  %.3e \n", ( double ) i);
        fprintf( outputStream,
                "Printing eigen vector no. %d, corresponding eigen value is %15.8e\n\n",
                i, eigVal.at(i) );

















=======
    fprintf(file, "\nLinear Stability:");
    fprintf(file, "\nEigen Values are:\n-----------------\n");

    for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        fprintf(file, "%15.8e ", eigVal.at(i) );
        if ( ( i % 5 ) == 0 ) {
            fprintf(file, "\n");
        }
    }

    fprintf(file, "\n\n");

    for ( int i = 0; i <= numberOfRequiredEigenValues; i++ ) {
        TimeStep step = *tStep;
        step.setTime( ( double ) i );
        step.setNumber(i);

        if ( i == 0 ) {
            fprintf(file, "\nLinear solution\n\n");
        } else {
            fprintf(file, "\nEigen vector no. %d, corresponding eigen value is %15.8e\n\n", i, eigVal.at(i));
        }
>>>>>>> .theirs

            for ( auto &dman : domain->giveDofManagers() ) {
            dman->updateYourself(&step);
            dman->printOutputAt(file, &step);
            }
		
<<<<<<< .mine
        for ( auto &dman : domain->giveDofManagers() ) {
            dman->updateYourself(tStep);
            dman->printOutputAt(outputStream, tStep);
        }

        tStep->setNumber(i);
        exportModuleManager->doOutput(tStep);
    }

=======









>>>>>>> .theirs
<<<<<<< .mine
#  ifdef VERBOSE
    VERBOSE_PRINT0("Updated nodes & sides ", nnodes)
#  endif
    fflush( this->giveOutputStream() );
    // save context if required
    this->saveStepContext(tStep);
=======
        if ( i == 0 ) {
            for ( auto &elem : domain->giveElements() ) {
                elem->printOutputAt(file, &step);



>>>>>>> .theirs
}
<<<<<<< .mine




=======
            this->printReactionForces(&step, 1., file);
        }
    }
}
>>>>>>> .theirs


void LinearStability :: setActiveVector(int activeVector)
//
// saves state variable - displacement vector
//
{
    this->giveCurrentStep()->setTime( ( double ) activeVector );
    this->giveCurrentStep()->setNumber( ( double ) activeVector );
        }


void LinearStability :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralEngngModel :: saveContext(stream, mode);

    field->saveContext(stream);

    contextIOResultType iores;
    if ( ( iores = eigVal.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    }


void LinearStability :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralEngngModel :: restoreContext(stream, mode);

    field->restoreContext(stream);

<<<<<<< .mine

contextIOResultType LinearStability :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore state variable - displacement vector
//
{
    int activeVector, version;
    int istep = 1, iversion = 1;
    int closeFlag = 0;
=======









>>>>>>> .theirs
    contextIOResultType iores;
    if ( ( iores = eigVal.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        }


} // end namespace oofem
