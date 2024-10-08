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

#include "sm/EngineeringModels/adaptnlinearstatic.h"
#include "mathfem.h"
#include "verbose.h"
#include "timer.h"
#include "metastep.h"
#include "timestep.h"
#include "nummet.h"
#include "element.h"
#include "node.h"
#include "domain.h"
#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "remeshingcrit.h"
#include "mesherinterface.h"
#include "dof.h"
#include "eleminterpunknownmapper.h"
#include "errorestimator.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "oofem_terminate.h"
#include "unknownnumberingscheme.h"

#ifdef __MPI_PARALLEL_MODE
 #include "parallelcontext.h"
 #include "loadbalancer.h"
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

#include <cstdlib>

namespace oofem {
REGISTER_EngngModel(AdaptiveNonLinearStatic);

AdaptiveNonLinearStatic :: AdaptiveNonLinearStatic(int i, EngngModel *_master) : NonLinearStatic(i, _master),
    d2_totalDisplacement(), d2_incrementOfDisplacement(), timeStepLoadLevels()
{
    meshPackage = MPT_T3D;
    equilibrateMappedConfigurationFlag = 0;

#ifdef __MPI_PARALLEL_MODE
    this->preMappingLoadBalancingFlag = false;
#endif
}


AdaptiveNonLinearStatic :: ~AdaptiveNonLinearStatic()
{ }


void
AdaptiveNonLinearStatic :: initializeFrom(InputRecord &ir)
{
    NonLinearStatic :: initializeFrom(ir);

    int _val;

    int meshPackageId = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, meshPackageId, _IFT_AdaptiveNonLinearStatic_meshpackage);
    meshPackage = ( MeshPackageType ) meshPackageId;

    equilibrateMappedConfigurationFlag =  0;
    IR_GIVE_OPTIONAL_FIELD(ir, equilibrateMappedConfigurationFlag, _IFT_AdaptiveNonLinearStatic_equilmc);
    _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_AdaptiveNonLinearStatic_preMappingLoadBalancingFlag);
    preMappingLoadBalancingFlag = _val > 0;


    // check if error estimator initioalized
    if (this->defaultErrEstimator == NULL) {
      OOFEM_ERROR ("AdaptiveNonLinearStatic :: initializeFrom: Error estimator not defined [eetype missing]");
    }
}

void
AdaptiveNonLinearStatic :: solveYourselfAt(TimeStep *tStep)
{
    proceedStep(1, tStep);
    this->updateYourself(tStep);

#ifdef __OOFEG
    ESIEventLoop( YES, const_cast< char * >("AdaptiveNonLinearStatic: Solution finished; Press Ctrl-p to continue") );
#endif

    //this->terminate( this->giveCurrentStep() );

#ifdef __MPI_PARALLEL_MODE
    if ( preMappingLoadBalancingFlag ) {
        this->balanceLoad( this->giveCurrentStep() );
    }

#endif

    // evaluate error of the reached solution
    this->defaultErrEstimator->estimateError( equilibratedEM, this->giveCurrentStep() );
    //this->defaultErrEstimator->estimateError( temporaryEM, this->giveCurrentStep() );
    this->defaultErrEstimator->giveRemeshingCrit()->estimateMeshDensities( this->giveCurrentStep() );
    RemeshingStrategy strategy = this->defaultErrEstimator->giveRemeshingCrit()->giveRemeshingStrategy( this->giveCurrentStep() );

    // if ((strategy == RemeshingFromCurrentState_RS) && (this->giveDomain(1)->giveSerialNumber() == 0))
    //  strategy = RemeshingFromPreviousState_RS;

    if ( strategy == NoRemeshing_RS ) {
        //
    } else if ( ( strategy == RemeshingFromCurrentState_RS ) || ( strategy == RemeshingFromPreviousState_RS ) ) {

        this->terminate( this->giveCurrentStep() ); // make output 

        // do remeshing
        auto mesher = classFactory.createMesherInterface( meshPackage, this->giveDomain(1) );

        Domain *newDomain;
        MesherInterface :: returnCode result = mesher->createMesh(this->giveCurrentStep(), 1,
                                                                  this->giveDomain(1)->giveSerialNumber() + 1, & newDomain);

        if ( result == MesherInterface :: MI_OK ) {
            this->initFlag = 1;
            this->adaptiveRemap(newDomain);
        } else if ( result == MesherInterface :: MI_NEEDS_EXTERNAL_ACTION ) {
            if ( strategy == RemeshingFromCurrentState_RS ) {
                // ensure the updating the step
                this->setContextOutputMode(COM_Always);
                //this->terminate (this->giveCurrentStep());
            } else {
                // save previous step (because update not called)
            }

            this->terminateAnalysis();
            throw OOFEM_Terminate();
        } else {
            OOFEM_ERROR("createMesh failed");
        }
    }
}

void
AdaptiveNonLinearStatic :: updateYourself(TimeStep *tStep)
{
    if ( timeStepLoadLevels.isEmpty() ) {
        timeStepLoadLevels.resize( this->giveNumberOfSteps() );
    }

    // in case of adaptive restart from given timestep
    // a time step with same number and incremented version will be generated.
    // Then load level reached on new discretization will overwrite the old one obtained for
    // old discretization for this step. But this is consistent, since when initialLoadVector
    // is requested to be recovered the reference load vectors are assembled
    // on actual discretization.
    timeStepLoadLevels.at( tStep->giveNumber() ) = loadLevel;

    NonLinearStatic :: updateYourself(tStep);
}


double AdaptiveNonLinearStatic :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    int eq = dof->__giveEquationNumber();
#ifdef DEBUG
    if ( eq == 0 ) {
        OOFEM_ERROR("invalid equation number");
    }
#endif

    if ( tStep != this->giveCurrentStep() ) {
        OOFEM_ERROR("unknown time step encountered");
    }

    if ( d->giveNumber() == 2 ) {
        switch ( mode ) {
        case VM_Incremental:
            if ( d2_incrementOfDisplacement.isNotEmpty() ) {
                return d2_incrementOfDisplacement.at(eq);
            } else {
                return 0.;
            }

        case VM_Total:
            if ( d2_totalDisplacement.isNotEmpty() ) {
                return d2_totalDisplacement.at(eq);
            } else {
                return 0.;
            }

        default:
            OOFEM_ERROR("Unknown is of undefined ValueModeType for this problem");
        }
    } else {
        return NonLinearStatic :: giveUnknownComponent(mode, tStep, d, dof);
    }

    // return 0.0;
}


int
AdaptiveNonLinearStatic :: initializeAdaptiveFrom(EngngModel *sourceProblem)
{
    int result = 1;

    // measure time consumed by mapping
    Timer timer;
    double mc1, mc2, mc3;
    timer.startTimer();

    if ( dynamic_cast< AdaptiveNonLinearStatic * >(sourceProblem) ) {
        OOFEM_ERROR("source problem must also be AdaptiveNonlinearStatic.");
    }

    this->currentStep = std::make_unique<TimeStep>( * ( sourceProblem->giveCurrentStep() ) );
    if ( sourceProblem->givePreviousStep() ) {
        this->previousStep = std::make_unique<TimeStep>( * ( sourceProblem->givePreviousStep() ) );
    }

    // map primary unknowns
    EIPrimaryUnknownMapper mapper;

    totalDisplacement.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    incrementOfDisplacement.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    totalDisplacement.zero();
    incrementOfDisplacement.zero();

    result &= mapper.mapAndUpdate( totalDisplacement, VM_Total,
                                  sourceProblem->giveDomain(1), this->giveDomain(1), sourceProblem->giveCurrentStep() );

    result &= mapper.mapAndUpdate( incrementOfDisplacement, VM_Incremental,
                                  sourceProblem->giveDomain(1), this->giveDomain(1), sourceProblem->giveCurrentStep() );

    timer.stopTimer();
    mc1 = timer.getUtime();
    timer.startTimer();

    // map internal ip state
    for ( auto &e : this->giveDomain(1)->giveElements() ) {
        result &= e->adaptiveMap( sourceProblem->giveDomain(1), sourceProblem->giveCurrentStep() );
    }

    timer.stopTimer();
    mc2 = timer.getUtime();
    timer.startTimer();

    // computes the stresses and calls updateYourself to mapped state
    for ( auto &e : this->giveDomain(1)->giveElements() ) {
        result &= e->adaptiveUpdate(currentStep.get());
    }

    // finish mapping process
    for ( auto &e : this->giveDomain(1)->giveElements() ) {
        result &= e->adaptiveFinish(currentStep.get());
    }


    // increment time step if mapped state will be considered as new solution stepL
    /*
     * this->giveNextStep();
     * if (equilibrateMappedConfigurationFlag) {
     * // we need to  assemble the load vector in same time as the restarted step,
     * // so new time step is generated with same intrincic time as has the
     * // previous step if equilibrateMappedConfigurationFlag is set.
     * // this allows to equlibrate the previously reached state
     * TimeStep* cts = this->giveCurrentStep();
     * cts->setTime(cts->giveTime()-cts->giveTimeIncrement());
     * cts = this->givePreviousStep();
     * cts->setTime(cts->giveTime()-cts->giveTimeIncrement());
     * }
     *
     *
     * if (this->giveCurrentStep()->giveNumber() ==
     * this->giveMetaStep(this->giveCurrentStep()->giveMetaStepNumber())->giveFirstStepNumber()) {
     * this->updateAttributes (this->giveCurrentStep());
     * }
     */
    this->updateAttributes( this->giveCurrentMetaStep() );

    // increment solution state counter - not needed, IPs are updated by adaptiveUpdate previously called
    // and there is no change in primary vars.
    // this->giveCurrentStep()->incrementStateCounter();

    // assemble new initial load for new discretization
    this->assembleInitialLoadVector( initialLoadVector, initialLoadVectorOfPrescribed,
                                    static_cast< AdaptiveNonLinearStatic * >(sourceProblem), 1, this->giveCurrentStep() );
    // assemble new total load for new discretization
    // this->assembleCurrentTotalLoadVector (totalLoadVector, totalLoadVectorOfPrescribed, this->giveCurrentStep());
    // set bcloadVector to zero (no increment within same step)

    timer.stopTimer();
    mc3 = timer.getUtime();

    // compute processor time used by the program
    OOFEM_LOG_INFO("user time consumed by primary mapping: %.2fs\n", mc1);
    OOFEM_LOG_INFO("user time consumed by ip mapping:      %.2fs\n", mc2);
    OOFEM_LOG_INFO("user time consumed by ip update:       %.2fs\n", mc3);
    OOFEM_LOG_INFO("user time consumed by mapping:         %.2fs\n", mc1 + mc2 + mc3);

    //

    //
    // bring mapped configuration into equilibrium
    //
    if ( equilibrateMappedConfigurationFlag ) {
        // use secant stiffness to resrore equilibrium
        NonLinearStatic_stiffnessMode oldStiffMode = this->stiffMode;
        stiffMode = nls_secantStiffness;



        if ( initFlag ) {
            stiffnessMatrix = classFactory.createSparseMtrx(sparseMtrxType);
            if ( !stiffnessMatrix ) {
                OOFEM_ERROR("sparse matrix creation failed");
            }

            if ( nonlocalStiffnessFlag ) {
                if ( !stiffnessMatrix->isAsymmetric() ) {
                    OOFEM_ERROR("stiffnessMatrix does not support asymmetric storage");
                }
            }

            stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
            stiffnessMatrix->zero(); // zero stiffness matrix
            this->assemble( *stiffnessMatrix, this->giveCurrentStep(), TangentAssembler(SecantStiffness),
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            initFlag = 0;
        }

        // updateYourself() not necessary - the adaptiveUpdate previously called does the job
        //this->updateYourself(this->giveCurrentStep());
#ifdef VERBOSE
        OOFEM_LOG_INFO( "Equilibrating mapped configuration [step number %5d.%d]\n",
                       this->giveCurrentStep()->giveNumber(), this->giveCurrentStep()->giveVersion() );
#endif

        //double deltaL = nMethod->giveUnknownComponent (StepLength, 0);
        double deltaL = nMethod->giveCurrentStepLength();
        this->assembleIncrementalReferenceLoadVectors( incrementalLoadVector, incrementalLoadVectorOfPrescribed,
                                                      refLoadInputMode, this->giveDomain(1), this->giveCurrentStep() );
        //
        // call numerical model to solve arised problem
        //
#ifdef VERBOSE
        OOFEM_LOG_RELEVANT( "Solving [step number %5d.%d]\n",
                           this->giveCurrentStep()->giveNumber(), this->giveCurrentStep()->giveVersion() );
#endif

        //nMethod -> solveYourselfAt(this->giveCurrentStep()) ;
        nMethod->setStepLength(deltaL / 5.0);
        if ( initialLoadVector.isNotEmpty() ) {
          numMetStatus = nMethod->solve( *stiffnessMatrix, incrementalLoadVector, & initialLoadVector,
                                          totalDisplacement, incrementOfDisplacement, internalForces,
                                          internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, this->giveCurrentStep() );
        } else {
          numMetStatus = nMethod->solve( *stiffnessMatrix, incrementalLoadVector, NULL, 
                                          totalDisplacement, incrementOfDisplacement, internalForces,
                                          internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, this->giveCurrentStep() );
        }


        loadVector.zero();

        this->updateYourself( this->giveCurrentStep() );
        this->terminate( this->giveCurrentStep() );
        this->updateLoadVectors( this->giveCurrentStep() );

        // restore old step length
        nMethod->setStepLength(deltaL);

        stiffMode = oldStiffMode;
    } else {
        // comment this, if output for mapped configuration (not equilibrated) not wanted
        this->printOutputAt( this->giveOutputStream(), this->giveCurrentStep() );
    }

    return result;
}


int
AdaptiveNonLinearStatic :: initializeAdaptive(int tStepNumber)
{
    try {
        FileDataStream stream(this->giveContextFileName(tStepNumber, 0), false);
        this->restoreContext(stream, CM_State);
    } catch(ContextIOERR & c) {
        c.print();
        exit(1);
    }

    this->initStepIncrements();

    int sernum = this->giveDomain(1)->giveSerialNumber();
    OOFEM_LOG_INFO("restoring domain %d.%d\n", 1, sernum + 1);
    Domain *dNew = new Domain(2, sernum + 1, this);
    OOFEMTXTDataReader domainDr(this->giveDomainFileName(1, sernum + 1));
    if ( !dNew->instanciateYourself(domainDr) ) {
        OOFEM_ERROR("domain Instanciation failed");
    }

    // remap solution to new domain
    return this->adaptiveRemap(dNew);
}


int
AdaptiveNonLinearStatic :: adaptiveRemap(Domain *dNew)
{
    int result = 1;

    this->initStepIncrements();

    this->ndomains = 2;
    this->domainNeqs.resize(2);
    this->domainPrescribedNeqs.resize(2);
    this->domainNeqs.at(2) = 0;
    this->domainPrescribedNeqs.at(2) = 0;
    this->domainList.emplace(domainList.begin() + 1, dNew);
    this->parallelContextList.emplace(parallelContextList.begin() + 1, this);

    // init equation numbering
    //this->forceEquationNumbering(2);
    this->forceEquationNumbering();

    // measure time consumed by mapping
    Timer timer;
    double mc1, mc2, mc3;

    timer.startTimer();

    // map primary unknowns
    EIPrimaryUnknownMapper mapper;

    d2_totalDisplacement.resize( this->giveNumberOfDomainEquations( 2, EModelDefaultEquationNumbering() ) );
    d2_incrementOfDisplacement.resize( this->giveNumberOfDomainEquations( 2, EModelDefaultEquationNumbering() ) );
    d2_totalDisplacement.zero();
    d2_incrementOfDisplacement.zero();

    result &= mapper.mapAndUpdate( d2_totalDisplacement, VM_Total,
                                  this->giveDomain(1), this->giveDomain(2), this->giveCurrentStep() );

    result &= mapper.mapAndUpdate( d2_incrementOfDisplacement, VM_Incremental,
                                  this->giveDomain(1), this->giveDomain(2), this->giveCurrentStep() );

    timer.stopTimer();
    mc1 = timer.getUtime();
    timer.startTimer();

    // map internal ip state
    for ( auto &e : this->giveDomain(2)->giveElements() ) {
        /* HUHU CHEATING */

        if ( e->giveParallelMode() == Element_remote ) {
            continue;
        }

        result &= e->adaptiveMap( this->giveDomain(1), this->giveCurrentStep() );
    }

    /* replace domains */
    OOFEM_LOG_DEBUG("deleting old domain\n");
    //delete domainList->at(1);
    //domainList->put(1, dNew);
    //dNew->setNumber(1);
    //domainList->put(2, NULL);

    //domainList[0] = std :: move(domainList[1]);
    domainList.erase(domainList.begin());
    domainList[0]->setNumber(1);

    parallelContextList = {parallelContextList[1]};

    // keep equation numbering of new domain
    this->numberOfEquations = this->domainNeqs.at(1) = this->domainNeqs.at(2);
    this->numberOfPrescribedEquations = this->domainPrescribedNeqs.at(1) = this->domainPrescribedNeqs.at(2);
    this->equationNumberingCompleted = 1;

    // update solution
    totalDisplacement = d2_totalDisplacement;
    incrementOfDisplacement = d2_incrementOfDisplacement;


    this->ndomains = 1;
    // init equation numbering
    // this->forceEquationNumbering();
    this->updateDomainLinks();

#ifdef __MPI_PARALLEL_MODE
    if ( isParallel() ) {
        // set up communication patterns
        this->initializeCommMaps(true);
        this->exchangeRemoteElementData(RemoteElementExchangeTag);
    }

#endif

    timer.stopTimer();
    mc2 = timer.getUtime();
    timer.startTimer();

    // computes the stresses and calls updateYourself to mapped state
    for ( auto &e : this->giveDomain(1)->giveElements() ) {
        /* HUHU CHEATING */
        if ( e->giveParallelMode() == Element_remote ) {
            continue;
        }

        result &= e->adaptiveUpdate( this->giveCurrentStep() );
    }

    // finish mapping process
    for ( auto &e : this->giveDomain(1)->giveElements() ) {
        /* HUHU CHEATING */
        if ( e->giveParallelMode() == Element_remote ) {
            continue;
        }

        result &= e->adaptiveFinish( this->giveCurrentStep() );
    }

    nMethod->reinitialize();


    // increment time step if mapped state will be considered as new solution stepL
    // this->giveNextStep();
    if ( equilibrateMappedConfigurationFlag ) {
        // we need to  assemble the load vector in same time as the restarted step,
        // so new time step is generated with same intrincic time as has the
        // previous step if equilibrateMappedConfigurationFlag is set.
        // this allows to equlibrate the previously reached state
        TimeStep *cts = this->giveCurrentStep();
        // increment version of solution step
        cts->incrementVersion();

        //cts->setTime(cts->giveTime()-cts->giveTimeIncrement());
        //cts = this->givePreviousStep();
        //cts->setTime(cts->giveTime()-cts->giveTimeIncrement());
    }

    if ( this->giveCurrentStep()->giveNumber() == this->giveCurrentMetaStep()->giveFirstStepNumber() ) {
        this->updateAttributes( this->giveCurrentMetaStep() );
    }

    // increment solution state counter - not needed, IPs are updated by adaptiveUpdate previously called
    // and there is no change in primary vars.
    // this->giveCurrentStep()->incrementStateCounter();

    // assemble new initial load for new discretization
    this->assembleInitialLoadVector( initialLoadVector, initialLoadVectorOfPrescribed,
                                    this, 1, this->giveCurrentStep() );
    this->assembleIncrementalReferenceLoadVectors( incrementalLoadVector, incrementalLoadVectorOfPrescribed,
                                                  refLoadInputMode, this->giveDomain(1),
                                                  this->giveCurrentStep() );

    // assemble new total load for new discretization
    // this->assembleCurrentTotalLoadVector (totalLoadVector, totalLoadVectorOfPrescribed, this->giveCurrentStep());
    // set bcloadVector to zero (no increment within same step)

    timer.stopTimer();
    mc3 = timer.getUtime();

    // compute processor time used by the program
    OOFEM_LOG_INFO("user time consumed by primary mapping: %.2fs\n", mc1);
    OOFEM_LOG_INFO("user time consumed by ip mapping:      %.2fs\n", mc2);
    OOFEM_LOG_INFO("user time consumed by ip update:       %.2fs\n", mc3);
    OOFEM_LOG_INFO("user time consumed by mapping:         %.2fs\n", mc1 + mc2 + mc3);

    //

    /********
     * #if 0
     * {
     *
     * // evaluate the force error of mapped configuration
     * this->updateComponent (this->giveCurrentStep(), InternalRhs);
     * FloatArray rhs;
     *
     * loadVector.resize (this->numberOfEquations);
     * loadVector.zero();
     * this->assemble (loadVector, this->giveCurrentStep(), ExternalForcesVector_Total, this->giveDomain(1)) ;
     * this->assemble (loadVector, this->giveCurrentStep(), ExternalForcesVector_Total, this->giveDomain(1));
     *
     * rhs =  loadVector;
     * rhs.times(loadLevel);
     * if (initialLoadVector.isNotEmpty()) rhs.add(initialLoadVector);
     * rhs.subtract(internalForces);
     *
     * //
     * // compute forceError
     * //
     * // err is relative error of unbalanced forces
     * double RR, RR0, forceErr = dotProduct (rhs.givePointer(),rhs.givePointer(),rhs.giveSize());
     * if (initialLoadVector.isNotEmpty())
     * RR0 = dotProduct (initialLoadVector.givePointer(), initialLoadVector.givePointer(), initialLoadVector.giveSize());
     * else
     * RR0 = 0.0;
     * RR = dotProduct(loadVector.givePointer(),loadVector.givePointer(),loadVector.giveSize());
     * // we compute a relative error norm
     * if ((RR0 + RR * loadLevel * loadLevel) < calm_SMALL_NUM) forceErr = 0.;
     * else forceErr = sqrt (forceErr / (RR0+RR * loadLevel * loadLevel));
     *
     * printf ("Relative Force Error of Mapped Configuration is %-15e\n", forceErr);
     *
     * }
     **#endif
     *************/

#ifdef __OOFEG
    ESIEventLoop( YES, const_cast< char * >("AdaptiveRemap: Press Ctrl-p to continue") );
#endif

    //
    // bring mapped configuration into equilibrium
    //
    if ( equilibrateMappedConfigurationFlag ) {
        // use secant stiffness to resrore equilibrium
        NonLinearStatic_stiffnessMode oldStiffMode = this->stiffMode;
        stiffMode = nls_secantStiffness;



        if ( initFlag ) {
            if ( !stiffnessMatrix ) {
                stiffnessMatrix = classFactory.createSparseMtrx(sparseMtrxType);
                if ( !stiffnessMatrix ) {
                    OOFEM_ERROR("sparse matrix creation failed");
                }
            }

            if ( nonlocalStiffnessFlag ) {
                if ( !stiffnessMatrix->isAsymmetric() ) {
                    OOFEM_ERROR("stiffnessMatrix does not support asymmetric storage");
                }
            }

            stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
            stiffnessMatrix->zero(); // zero stiffness matrix
            this->assemble( *stiffnessMatrix, this->giveCurrentStep(), TangentAssembler(SecantStiffness),
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );

            initFlag = 0;
        }

        // updateYourself() not necessary - the adaptiveUpdate previously called does the job
        //this->updateYourself(this->giveCurrentStep());
#ifdef VERBOSE
        OOFEM_LOG_INFO( "Equilibrating mapped configuration [step number %5d.%d]\n",
                       this->giveCurrentStep()->giveNumber(), this->giveCurrentStep()->giveVersion() );
#endif

        //double deltaL = nMethod->giveUnknownComponent (StepLength, 0);
        double deltaL = nMethod->giveCurrentStepLength();
        //
        // call numerical model to solve arised problem
        //
#ifdef VERBOSE
        OOFEM_LOG_RELEVANT( "Solving [step number %5d.%d]\n",
                           this->giveCurrentStep()->giveNumber(), this->giveCurrentStep()->giveVersion() );
#endif

        //nMethod -> solveYourselfAt(this->giveCurrentStep()) ;
        nMethod->setStepLength(deltaL / 5.0);
        if ( initialLoadVector.isNotEmpty() ) {
          numMetStatus = nMethod->solve( *stiffnessMatrix, incrementalLoadVector, & initialLoadVector,
                                          totalDisplacement, incrementOfDisplacement, internalForces,
                                          internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, this->giveCurrentStep() );
        } else {
          numMetStatus = nMethod->solve( *stiffnessMatrix, incrementalLoadVector, NULL,
                                          totalDisplacement, incrementOfDisplacement, internalForces,
                                          internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, this->giveCurrentStep() );
        }


        loadVector.zero();

        this->updateYourself( this->giveCurrentStep() );
        this->terminate( this->giveCurrentStep() );
        // this->updateLoadVectors (this->giveCurrentStep()); // already in terminate

        // restore old step length
        nMethod->setStepLength(deltaL);

        stiffMode = oldStiffMode;
    } else {
        // comment this, if output for mapped configuration (not equilibrated) not wanted
        this->printOutputAt( this->giveOutputStream(), this->giveCurrentStep() );
    }

    return result;
}


void
AdaptiveNonLinearStatic :: saveContext(DataStream &stream, ContextMode mode)
{
    NonLinearStatic :: saveContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = timeStepLoadLevels.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}

void
AdaptiveNonLinearStatic :: restoreContext(DataStream &stream, ContextMode mode)
{
    NonLinearStatic :: restoreContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = timeStepLoadLevels.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


void
AdaptiveNonLinearStatic :: updateDomainLinks()
{
    NonLinearStatic :: updateDomainLinks();
    this->defaultErrEstimator->setDomain( this->giveDomain(1) );
}


void
AdaptiveNonLinearStatic :: assembleInitialLoadVector(FloatArray &loadVector, FloatArray &loadVectorOfPrescribed,
                                                     AdaptiveNonLinearStatic *sourceProblem, int domainIndx,
                                                     TimeStep *tStep)
{
    int mStepNum = tStep->giveMetaStepNumber();
    int hasfixed, mode;
    FloatArray _incrementalLoadVector, _incrementalLoadVectorOfPrescribed;
    SparseNonLinearSystemNM :: referenceLoadInputModeType rlm;
    //Domain* sourceDomain = sourceProblem->giveDomain(domainIndx);

    loadVector.resize( this->giveNumberOfDomainEquations( domainIndx, EModelDefaultEquationNumbering() ) );
    loadVectorOfPrescribed.resize( this->giveNumberOfDomainEquations( domainIndx, EModelDefaultPrescribedEquationNumbering() ) );
    loadVector.zero();
    loadVectorOfPrescribed.zero();
    _incrementalLoadVector.resize( this->giveNumberOfDomainEquations( domainIndx, EModelDefaultEquationNumbering() ) );
    _incrementalLoadVectorOfPrescribed.resize( this->giveNumberOfDomainEquations( domainIndx, EModelDefaultPrescribedEquationNumbering() ) );
    _incrementalLoadVector.zero();
    _incrementalLoadVectorOfPrescribed.zero();

    for ( int imstep = 1; imstep < mStepNum; imstep++ ) {
        auto iMStep = this->giveMetaStep(imstep);
        auto &ir = iMStep->giveAttributesRecord();
        //hasfixed = ir.hasField("fixload");
        hasfixed = 1;
        if ( hasfixed ) {
            // test for control mode
            // here the algorithm works only for direct load control.
            // Direct displacement control requires to know the quasi-rections, and the controlled nodes
            // should have corresponding node on new mesh -> not supported
            // Indirect control -> the load level from prevous steps is required, currently nt supported.

            // additional problem: direct load control supports the reduction of step legth if convergence fails
            // if this happens, this implementation does not work correctly.
            // But there is NO WAY HOW TO TEST IF THIS HAPPEN

            mode = 0;
            IR_GIVE_OPTIONAL_FIELD(ir, mode, _IFT_AdaptiveNonLinearStatic_controlmode);

            // check if displacement control takes place
            if ( ir.hasField(_IFT_AdaptiveNonLinearStatic_ddm) ) {
                OOFEM_ERROR("fixload recovery not supported for direct displacement control");
            }

            int firststep = iMStep->giveFirstStepNumber();
            int laststep  = iMStep->giveLastStepNumber();

            int _val = 0;
            IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_AdaptiveNonLinearStatic_refloadmode);
            rlm = ( SparseNonLinearSystemNM :: referenceLoadInputModeType ) _val;

            if ( mode == ( int ) nls_directControl ) { // and only load control
                for ( int istep = firststep; istep <= laststep; istep++ ) {
                    // bad practise here
                    ///@todo Likely memory leak here with new TimeStep; Check.
                    TimeStep *old = new TimeStep(istep, this, imstep, istep - 1.0, deltaT, 0);
                    this->assembleIncrementalReferenceLoadVectors(_incrementalLoadVector, _incrementalLoadVectorOfPrescribed,
                                                                  rlm, this->giveDomain(domainIndx), old);

                    _incrementalLoadVector.times( sourceProblem->giveTimeStepLoadLevel(istep) );
                    loadVector.add(_incrementalLoadVector);
                    loadVectorOfPrescribed.add(_incrementalLoadVectorOfPrescribed);
                }
            } else if ( mode == ( int ) nls_indirectControl ) {
                // bad practise here
                if ( !ir.hasField(_IFT_NonLinearStatic_donotfixload) ) {
                    TimeStep *old = new TimeStep(firststep, this, imstep, firststep - 1.0, deltaT, 0);
                    this->assembleIncrementalReferenceLoadVectors(_incrementalLoadVector, _incrementalLoadVectorOfPrescribed,
                                                                  rlm, this->giveDomain(domainIndx), old);

                    _incrementalLoadVector.times( sourceProblem->giveTimeStepLoadLevel(laststep) );
                    loadVector.add(_incrementalLoadVector);
                    loadVectorOfPrescribed.add(_incrementalLoadVectorOfPrescribed);
                }
            } else {
                OOFEM_ERROR("fixload recovery not supported");
            }
        }
    } // end loop over meta-steps

    /* if direct control; add to initial load also previous steps in same metestep */
    auto iMStep = this->giveMetaStep(mStepNum);
    auto &ir = iMStep->giveAttributesRecord();
    mode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, mode, _IFT_AdaptiveNonLinearStatic_controlmode);
    int firststep = iMStep->giveFirstStepNumber();
    int laststep  = tStep->giveNumber();
    int _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_AdaptiveNonLinearStatic_refloadmode);
    rlm = ( SparseNonLinearSystemNM :: referenceLoadInputModeType ) _val;

    if ( mode == ( int ) nls_directControl ) { // and only load control
        for ( int istep = firststep; istep <= laststep; istep++ ) {
            // bad practise here
            TimeStep *old = new TimeStep(istep, this, mStepNum, istep - 1.0, deltaT, 0);
            this->assembleIncrementalReferenceLoadVectors(_incrementalLoadVector, _incrementalLoadVectorOfPrescribed,
                                                          rlm, this->giveDomain(domainIndx), old);

            _incrementalLoadVector.times( sourceProblem->giveTimeStepLoadLevel(istep) );
            loadVector.add(_incrementalLoadVector);
            loadVectorOfPrescribed.add(_incrementalLoadVectorOfPrescribed);
        }
    }
}

/*
 * void
 * AdaptiveNonLinearStatic::assembleCurrentTotalLoadVector (FloatArray& loadVector,
 *                           FloatArray& loadVectorOfPrescribed,
 *                           AdaptiveNonLinearStatic* sourceProblem, int domainIndx,
 *                           TimeStep* tStep)
 * {
 * int mStepNum = tStep->giveMetaStepNumber() ;
 * int mode;
 * InputRecord* ir;
 * MetaStep* mStep = sourceProblem->giveMetaStep(mStepNum);
 * FloatArray _incrementalLoadVector, _incrementalLoadVectorOfPrescribed;
 * SparseNonLinearSystemNM::referenceLoadInputModeType rlm;
 * //Domain* sourceDomain = sourceProblem->giveDomain(domainIndx);
 *
 * loadVector.resize(this->giveNumberOfDomainEquations());
 * loadVectorOfPrescribed.resize(this->giveNumberOfPrescribedDomainEquations());
 * loadVector.zero(); loadVectorOfPrescribed.zero();
 * _incrementalLoadVector.resize(this->giveNumberOfDomainEquations());
 * _incrementalLoadVectorOfPrescribed.resize(this->giveNumberOfPrescribedDomainEquations());
 * _incrementalLoadVector.zero(); _incrementalLoadVectorOfPrescribed.zero();
 *
 * // ASK CURRENT MSTEP FOR ITS RECORD
 * ir = mStep->giveAttributesRecord();
 *
 * mode = 0;
 * IR_GIVE_OPTIONAL_FIELD (ir, mode, _IFT_AdaptiveNonLinearStatic_controlmode, "controlmode");
 *
 * // check if displacement control takes place
 * if (ir->hasField(_IFT_AdaptiveNonLinearStatic_ddm, "ddm"))
 * OOFEM_ERROR("fixload recovery not supported for direct displacement control");
 * int _val = 0;
 * IR_GIVE_OPTIONAL_FIELD (ir, _val, _IFT_AdaptiveNonLinearStatic_refloadmode, "refloadmode");
 *
 * int firststep = mStep->giveFirstStepNumber();
 * int laststep  = tStep->giveNumber()-1;
 *
 * rlm = (SparseNonLinearSystemNM::referenceLoadInputModeType) _val;
 *
 * if (mode == (int)nls_directControl) { // and only load control
 * for (int istep = firststep; istep<=laststep; istep++) {
 * // bad practise here
 * TimeStep* old = new TimeStep (istep, this, mStepNum, istep-1.0, deltaT, 0);
 * this->assembleIncrementalReferenceLoadVectors (_incrementalLoadVector, _incrementalLoadVectorOfPrescribed,
 *                         rlm, this->giveDomain(domainIndx), old);
 *
 * _incrementalLoadVector.times(sourceProblem->giveTimeStepLoadLevel(istep));
 * loadVector.add(_incrementalLoadVector);
 * loadVectorOfPrescribed.add(_incrementalLoadVectorOfPrescribed);
 * }
 * } else if (mode == (int)nls_indirectControl) {
 * // bad practise here
 * TimeStep* old = new TimeStep (firststep, this, mStepNum, firststep-1.0, deltaT, 0);
 * this->assembleIncrementalReferenceLoadVectors (_incrementalLoadVector, _incrementalLoadVectorOfPrescribed,
 *                        rlm, this->giveDomain(domainIndx), old);
 *
 * _incrementalLoadVector.times(sourceProblem->giveTimeStepLoadLevel(laststep));
 * loadVector.add(_incrementalLoadVector);
 * loadVectorOfPrescribed.add(_incrementalLoadVectorOfPrescribed);
 * } else {
 * OOFEM_ERROR("totalload recovery not supported");
 * }
 *
 *
 * }
 */


double
AdaptiveNonLinearStatic :: giveTimeStepLoadLevel(int istep)
{
    if ( ( istep >= this->giveNumberOfFirstStep() ) && ( istep <= this->giveNumberOfSteps() ) ) {
        return timeStepLoadLevels.at(istep);
    } else {
        OOFEM_ERROR("solution step out of range");
    }

    // return 0.0; // to make compiler happy
}


#ifdef __MPI_PARALLEL_MODE
LoadBalancer *
AdaptiveNonLinearStatic :: giveLoadBalancer()
{
    if ( lb ) {
        return lb.get();
    }

    if ( loadBalancingFlag || preMappingLoadBalancingFlag ) {
        lb = classFactory.createLoadBalancer( "parmetis", this->giveDomain(1) );
        return lb.get();
    } else {
        return nullptr;
    }
}
LoadBalancerMonitor *
AdaptiveNonLinearStatic :: giveLoadBalancerMonitor()
{
    if ( lbm ) {
        return lbm.get();
    }

    if ( loadBalancingFlag || preMappingLoadBalancingFlag ) {
        lbm = classFactory.createLoadBalancerMonitor( "wallclock", this);
        return lbm.get();
    } else {
        return nullptr;
    }
}
#endif
} // end namespace oofem
