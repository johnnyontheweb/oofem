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

#include "sm/EngineeringModels/varlinearstability.h"
#include "oofemcfg.h"
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
#include "outputmanager.h"
#include "eigenvectorprimaryfield.h"


#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

#ifdef MEMSTR
    #include <io.h>
    #include <fcntl.h>
#endif

namespace oofem {
REGISTER_EngngModel(VarLinearStability);

VarLinearStability::VarLinearStability(int i, EngngModel *master) : StructuralEngngModel(i, master),
    numberOfRequiredEigenValues(1),
    rtolv(1e-6),
    solverType(GES_SubspaceIt)
{
    numberOfSteps = 1;
    ndomains = 1;
}


NumericalMethod *VarLinearStability :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
        nMethod = classFactory.createGeneralizedEigenValueSolver(solverType, this->giveDomain(1), this);
        if ( !nMethod ) {
            OOFEM_ERROR("solver creation failed");
        }
    }

    return nMethod.get();
}


SparseLinearSystemNM *VarLinearStability :: giveNumericalMethodForLinStaticProblem(TimeStep *tStep)
{
    if ( !nMethodLS ) {
        nMethodLS = classFactory.createSparseLinSolver( linStype, this->giveDomain( 1 ), this ); ///@todo Support other solvers
        if ( !nMethodLS ) {
            OOFEM_ERROR("solver creation failed");
        }
    }

    return nMethodLS.get();
}


void
VarLinearStability :: initializeFrom(InputRecord &ir)
{
    //StructuralEngngModel::instanciateFrom(ir);
    IR_GIVE_FIELD(ir, numberOfRequiredEigenValues, _IFT_VarLinearStability_nroot);
    // numberOfSteps set artifficially to numberOfRequiredEigenValues
    // in order to allow
    // use restoreContext function for different eigenValues
    numberOfSteps = numberOfRequiredEigenValues;
    this->field = std::make_unique<EigenVectorPrimaryField>(this, 1, FT_Displacements, numberOfRequiredEigenValues + 1); // +1 for eq. solution

    IR_GIVE_FIELD(ir, rtolv, _IFT_VarLinearStability_rtolv);
    if ( rtolv < 1.e-12 ) {
        rtolv =  1.e-12;
    } else if ( rtolv > 0.01 ) {
        rtolv =  0.01;
    }

    int val = 1; // inverseit
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_VarLinearStability_stype);
    solverType = ( GenEigvalSolverType ) val;

    sparseMtrxType = SparseMtrxType::SMT_Skyline; //Default Skyline
    linStype = LinSystSolverType::ST_Direct;

    if (solverType == GenEigvalSolverType::GES_Eigen){
	    sparseMtrxType = SparseMtrxType::SMT_EigenSparse; // linStype = ST_Spooles;
	    linStype = LinSystSolverType::ST_EigenLib;
    }

    IR_GIVE_OPTIONAL_FIELD( ir, flexkg, _IFT_VarLinearStability_flexkg );
    
    nMetaSteps = 0;

    suppressOutput = ir.hasField(_IFT_EngngModel_suppressOutput);

    if (suppressOutput) {
        printf("Suppressing output.\n");
    } else {
#ifdef MEMSTR
        outputStream = nullptr;
        FILE* source = classFactory.giveMemoryStream("out");
        int sourceFD = _open_osfhandle((intptr_t)source, _O_APPEND);
        if (sourceFD != -1) { outputStream = _fdopen(sourceFD, "a"); }
        if (!(outputStream)) {
            // if not, write to file
#endif
            if ((outputStream = fopen(this->dataOutputFileName.c_str(), "w")) == NULL) { OOFEM_ERROR("Can't open output file %s", this->dataOutputFileName.c_str()); }
#ifdef MEMSTR
            usestream = false;
        }
#endif

        fprintf(outputStream, "%s", PRG_HEADER);
        fprintf(outputStream, "\nStarting analysis on: %s\n", ctime(& this->startTime) );
        fprintf(outputStream, "%s\n", simulationDescription.c_str());
    }
}


int VarLinearStability :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
{
    return tStep->giveNumber() % (this->numberOfRequiredEigenValues + 1); // +1 for eq. solution 
}


double VarLinearStability :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    return field->giveUnknownValue(dof, mode, tStep);
}


TimeStep *VarLinearStability :: giveNextStep()
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


void VarLinearStability :: solveYourself()
{
    this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);
    // update state according to new meta step
    this->giveNextStep();
    this->updateAttributes( this->giveCurrentMetaStep() );
    this->solveYourselfAt( this->giveCurrentStep() );
    this->terminate( this->giveCurrentStep() );
}


void VarLinearStability :: solveYourselfAt(TimeStep *tStep)
{
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );
    FloatArray displacementVector( neq ), loadVector( neq ), loadVector2( neq );

    tStep->setNumber(0); // variable loading a t=0
    tStep->setTime(0.0);

    // creates system of governing eq's and solves them at given time step
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );
    this->giveNumericalMethodForLinStaticProblem(tStep);

    // first assemble problem at current time step
    if ( !stiffnessMatrix ) {
        //
        // first step - solve linear static problem
        //
        stiffnessMatrix = classFactory.createSparseMtrx( sparseMtrxType );
        stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
    }

#ifndef LIN_STAB_COMPATIBILITY_MODE
    OOFEM_LOG_INFO("Assembling stiffness matrix\n");
    stiffnessMatrix->zero();
    this->assemble( *stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
                   EModelDefaultEquationNumbering(), this->giveDomain(1) );
#endif

    // Normal forces already known, proceed with linear stability
    // stiffnessMatrix->zero();
    if ( !initialStressMatrix ) {
        initialStressMatrix = stiffnessMatrix->clone(); // copy matrix and its values
    } else {
        initialStressMatrix->zero();
    }

    // OOFEM_LOG_INFO("Assembling stiffness matrix\n"); // is this here for changes in stiffness after first resolution? if no, use initial
    // this->assemble( *stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
    //                EModelDefaultEquationNumbering(), this->giveDomain(1) );

    // -----------------------------------------------------------------
    // create the actual initial stiffness matrix using only time = 1
    TimeStep tStep1( *tStep );
    tStep1.setNumber( 1 );
    tStep1.setTime( 1.0 ); // constant loading at t=1
    TimeStep *tStep1Ptr = &tStep1;

    //  Internal forces first, negated;
    field->update( VM_Total, tStep1Ptr, displacementVector, EModelDefaultEquationNumbering() ); // not really needed here
    this->assembleVector( loadVector2, tStep1Ptr, InternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain( 1 ) );
    loadVector2.negated();
    // external forces
    this->assembleVector( loadVector2, tStep1Ptr, ExternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain( 1 ) );
    this->updateSharedDofManagers( loadVector2, EModelDefaultEquationNumbering(), ReactionExchangeTag );

    if ( loadVector2.computeNorm() > 1.e-10 ) { // constant loading
        //loadVector.add( loadVector2 ); // sum to be consistent with deformed shape?
        OOFEM_LOG_INFO( "Solving linear static problem\n" );
        nMethodLS->solve( *stiffnessMatrix, loadVector2, displacementVector );
        // Initial displacements are stored at position 0; this is a bit of a hack. In the future, a cleaner approach of handling fields could be suitable,
        // but currently, it all converges down to the same giveUnknownComponent, so this is the easiest approach.
        field->update( VM_Total, tStep1Ptr, displacementVector, EModelDefaultEquationNumbering() );
        // terminate linear static computation (necessary, in order to compute stresses in elements).
        // Recompute for updated state:
        this->assembleVector(loadVector2, tStep1Ptr, InternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1));
        this->terminateLinStatic( tStep1Ptr );

        initialStressMatrix->zero(); // rewrite initialStressMatrix
        OOFEM_LOG_INFO( "Assembling initial stress matrix for constant loading\n" );
        this->assemble( *initialStressMatrix, tStep1Ptr, InitialStressMatrixAssembler(), EModelDefaultEquationNumbering(), this->giveDomain( 1 ) );
        // add constant loads to update the stiffness matrix
        stiffnessMatrix->add( 1.0, *initialStressMatrix );
        displacementVector.zero(); // don't consider deformed shape for multipliers of variable loading, otherwise it will add constant loads again
    } // constant loading
    // -----------------------------------------------------------------

    initialStressMatrix->zero(); // rewrite initialStressMatrix
    OOFEM_LOG_INFO("Assembling variable loading\n");
    // Internal forces first, negated;
    field->update(VM_Total, tStep, displacementVector, EModelDefaultEquationNumbering());
    this->assembleVector( loadVector, tStep, InternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    loadVector.negated();
    // external forces
    this->assembleVector( loadVector, tStep, ExternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->updateSharedDofManagers(loadVector, EModelDefaultEquationNumbering(), ReactionExchangeTag);

    OOFEM_LOG_INFO("Solving linear static problem\n");
    nMethodLS->solve(*stiffnessMatrix, loadVector, displacementVector);
    // Initial displacements are stored at position 0; this is a bit of a hack. In the future, a cleaner approach of handling fields could be suitable,
    // but currently, it all converges down to the same giveUnknownComponent, so this is the easiest approach.
    field->update(VM_Total, tStep, displacementVector, EModelDefaultEquationNumbering());
    // terminate linear static computation (necessary, in order to compute stresses in elements).
    // Recompute for updated state:
    this->assembleVector( loadVector, tStep, InternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->terminateLinStatic( tStep );

    OOFEM_LOG_INFO("Assembling initial stress matrix for variable loading\n");
    this->assemble( *initialStressMatrix, tStep, InitialStressMatrixAssembler(),
                   EModelDefaultEquationNumbering(), this->giveDomain(1) );
    initialStressMatrix->times(-1.0);

    FloatMatrix eigVec(neq, numberOfRequiredEigenValues);
    eigVal.resize(numberOfRequiredEigenValues);
    eigVal.zero();

    OOFEM_LOG_INFO("Solving ...\n");

#ifdef DEBUG
    //  stiffnessMatrix->printYourself();
    //  initialStressMatrix->printYourself();
    
    //stiffnessMatrix->writeToFile("K.dat");
    //initialStressMatrix->writeToFile("M.dat");
#endif

    auto cr = nMethod->solve(*stiffnessMatrix, *initialStressMatrix, eigVal, eigVec, rtolv, numberOfRequiredEigenValues);
    this->field->updateAll(eigVec, EModelDefaultEquationNumbering());
    if ( cr != CR_CONVERGED ) {
        OOFEM_ERROR( "Buckling solver couldn't find a solution." );
    }
}


void VarLinearStability :: updateYourself(TimeStep *tStep)
{ }


void
VarLinearStability :: terminateLinStatic(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);

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
#if 0
    // save context if required
    // default - save only if ALWAYS is set ( see cltypes.h )

    if ((domain->giveContextOutputMode() == COM_Always ) ||
        (domain->giveContextOutputMode() == COM_Required )) {
        this->saveContext(NULL);
    } else if (domain->giveContextOutputMode() == COM_UserDefined ) {
        if (tStep->giveNumber()%domain->giveContextOutputStep() == 0)
            this->saveContext(NULL);
    }
#endif
}


void VarLinearStability :: doStepOutput(TimeStep *tStep)
{
    if ( !suppressOutput ) {
        this->printOutputAt(this->giveOutputStream(), tStep);
        fflush( this->giveOutputStream() );
    }

    //Domain *domain = this->giveDomain(1);
    //// i = 0  represents the linear solution, which is followed by the eigen vectors starting at i = 1
    //for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
    //    TimeStep step = *tStep;
    //    step.setTime( ( double ) i );
    //    step.setNumber(i);

    //    for ( auto &dman : domain->giveDofManagers() ) {
    //        dman->updateYourself(&step);
    //    }
    //    for ( auto &elem : domain->giveElements() ) {
    //        elem->updateYourself(&step);
    //    }

    //    exportModuleManager.doOutput(&step);
    //}
}


void VarLinearStability :: printOutputAt(FILE *file, TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    if ( !domain->giveOutputManager()->testTimeStepOutput(tStep) ) {
        return;
    }

    fprintf(file, "\nLinear Stability:");
    fprintf(file, "\nEigen Values are:\n-----------------\n");

    for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        fprintf(file, "%15.8e ", eigVal.at(i) );
        if ( ( i % 5 ) == 0 ) {
            fprintf(file, "\n");
        }
    }

    fprintf(file, "\n\n");

    for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        TimeStep step = *tStep;
        step.setTime( ( double ) i );
        step.setNumber(i);
        // added for output
        fprintf( file, "\nOutput for eigen value no.  %.3e \n", (double)i );
        fprintf( file,"Printing eigen vector no. %d, corresponding eigen value is %15.8e\n\n", i, eigVal.at( i ) );
        // end added for output
        
        //if ( i == 0 ) {
        //    fprintf(file, "\nLinear solution\n\n");
        //} else {
        //    fprintf(file, "\nEigen vector no. %d, corresponding eigen value is %15.8e\n\n", i, eigVal.at(i));
        //}

        for ( auto &dman : domain->giveDofManagers() ) {
            dman->updateYourself(&step);
            dman->printOutputAt(file, &step);
        }

        //if ( i == 0 ) {
            //for ( auto &elem : domain->giveElements() ) {
            //    elem->updateYourself( &step );
            //    elem->printOutputAt(file, &step);
            //}
        //this->printReactionForces(&step, 1., file);
        //}

        for ( auto &elem : domain->giveElements() ) {
            elem->updateYourself( &step );
        }

        exportModuleManager.doOutput( &step );
    }
}


void VarLinearStability :: setActiveVector(int activeVector)
{
    this->giveCurrentStep()->setTime( ( double ) activeVector );
    this->giveCurrentStep()->setNumber( ( double ) activeVector );
}


void VarLinearStability :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralEngngModel :: saveContext(stream, mode);

    field->saveContext(stream);

    contextIOResultType iores;
    if ( ( iores = eigVal.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


void VarLinearStability :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralEngngModel :: restoreContext(stream, mode);

    field->restoreContext(stream);

    contextIOResultType iores;
    if ( ( iores = eigVal.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


} // end namespace oofem
