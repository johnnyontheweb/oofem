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

#include "sm/EngineeringModels/pdeltastatic.h"
#include "sm/Elements/structuralelement.h"
#include "sm/Elements/structuralelementevaluator.h"
#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "dof.h"
#include "sparsemtrx.h"
#include "verbose.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "unknownnumberingscheme.h"
#include "geneigvalsolvertype.h"
#include "sparsegeneigenvalsystemnm.h"

#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
 #include "communicator.h"
#endif

#include <typeinfo>

namespace oofem {
REGISTER_EngngModel(PDeltaStatic);

PDeltaStatic :: PDeltaStatic(int i, EngngModel *_master) : StructuralEngngModel(i, _master), loadVector(), displacementVector()
{
    ndomains = 1;
    initFlag = 1;
    solverType = ST_Direct;
}


PDeltaStatic :: ~PDeltaStatic() { }


NumericalMethod *PDeltaStatic :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
        if ( isParallel() ) {
            if ( ( solverType == ST_Petsc ) || ( solverType == ST_Feti ) ) {
                nMethod = classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this);
            }
        } else {
            nMethod = classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this);
        }
        if ( !nMethod ) {
            OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
        }
    }


    return nMethod.get();
}

void
PDeltaStatic :: initializeFrom(InputRecord &ir)
{
    StructuralEngngModel :: initializeFrom(ir);

    // int val = 0;
    // IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    // solverType = ( LinSystSolverType ) val;
	solverType = ST_EigenLib;

    // val = 0;
    // IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    // sparseMtrxType = ( SparseMtrxType ) val;
	sparseMtrxType = SMT_EigenSparse;

	IR_GIVE_FIELD(ir, rtolv, _IFT_PDeltaStatic_rtolv);

#ifdef __PARALLEL_MODE
    if ( isParallel() ) {
        commBuff = new CommunicatorBuff( this->giveNumberOfProcesses() );
        communicator = new NodeCommunicator(this, commBuff, this->giveRank(),
                                            this->giveNumberOfProcesses());
    }

#endif
}


double PDeltaStatic :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation eq
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
#ifdef DEBUG
    if ( eq == 0 ) {
        OOFEM_ERROR("invalid equation number");
    }
#endif

    if ( tStep != this->giveCurrentStep() ) {
        OOFEM_ERROR("unknown time step encountered");
        return 0.;
    }

    switch ( mode ) {
    case VM_Total:
    case VM_Incremental:
        if ( displacementVector.isNotEmpty() ) {
            return displacementVector.at(eq);
        } else {
            return 0.;
        }
	case VM_Velocity:
		return 0.;
	case VM_Acceleration:
		return 0.;

    default:
        OOFEM_ERROR("Unknown is of undefined type for this problem");
    }

    return 0.;
}


TimeStep *PDeltaStatic :: giveNextStep()
{
    if ( !currentStep ) {
        // first step -> generate initial step
        //currentStep.reset( new TimeStep(*giveSolutionStepWhenIcApply()) );
        currentStep.reset( new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 1, 0., 1., 0) );
    }
    previousStep = std :: move(currentStep);
    currentStep.reset( new TimeStep(*previousStep, 1.) );

    return currentStep.get();
}

void PDeltaStatic :: solveYourself()
{
    if ( this->isParallel() ) {
#ifdef __VERBOSE_PARALLEL
        // force equation numbering before setting up comm maps
        int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering());
        OOFEM_LOG_INFO("[process rank %d] neq is %d\n", this->giveRank(), neq);
#endif

        this->initializeCommMaps();
    }

    StructuralEngngModel :: solveYourself();
}


void PDeltaStatic :: solveYourselfAt(TimeStep *tStep)
{
    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step

    if ( initFlag ) {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling stiffness matrix\n");
#endif

        //
        // first step  assemble stiffness Matrix
        //
        stiffnessMatrix = classFactory.createSparseMtrx(sparseMtrxType);
        if ( !stiffnessMatrix ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

        this->assemble( *stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );

        initFlag = 0;

		initialStressMatrix = classFactory.createSparseMtrx(sparseMtrxType);
		initialStressMatrix->buildInternalStructure(this, 1, EModelDefaultEquationNumbering());
    }

#ifdef VERBOSE
    OOFEM_LOG_DEBUG("Assembling load\n");
#endif

    //
    // allocate space for displacementVector
    //
    displacementVector.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    displacementVector.zero();

    //
    // assembling the load vector
    //
    loadVector.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    loadVector.zero();
    this->assembleVector( loadVector, tStep, ExternalForceAssembler(), VM_Total,
                         EModelDefaultEquationNumbering(), this->giveDomain(1) );

    //
    // internal forces (from Dirichlet b.c's, or thermal expansion, etc.)
    //
    FloatArray internalForces( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    internalForces.zero();
    this->assembleVector( internalForces, tStep, InternalForceAssembler(), VM_Total,
                         EModelDefaultEquationNumbering(), this->giveDomain(1) );

    loadVector.subtract(internalForces);

    this->updateSharedDofManagers(loadVector, EModelDefaultEquationNumbering(), ReactionExchangeTag);

    //
    // set-up numerical model
    //
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );

    //
    // call numerical model to solve arose problem
    //
#ifdef VERBOSE
    OOFEM_LOG_INFO("\n\nSolving initial ...\n\n");
#endif
	//NM_Status s = 
	nMethod->solve(*stiffnessMatrix, loadVector, displacementVector);
	// norm of previous displ. vector
	double oldNorm = displacementVector.computeSquaredNorm(); double newNorm = 0;
	bool escape = false; int maxIter = 0;
	//initialStressMatrix.reset(classFactory.createSparseMtrx(sparseMtrxType)); // stresses are in the model now
	//initialStressMatrix->buildInternalStructure(this, 1, EModelDefaultEquationNumbering());
    FloatArray rhs = loadVector;
	do {
		// PDELTA approx solution with iterations - maximum 10 iterations
		if (newNorm!=0) oldNorm = newNorm;
		maxIter += 1;

        if ( true ) {
                
		// terminate linear static computation (necessary, in order to compute stresses in elements).
		this->updateAfterStatic(this->giveCurrentStep(), this->giveDomain(1)); // not needed for beam - conservatively left (shells?)
#ifdef VERBOSE
	OOFEM_LOG_INFO("Assembling initial stress matrix\n");
#endif
		// initialStressMatrix.reset(stiffnessMatrix->GiveCopy());
		initialStressMatrix->zero();

		this->assemble(*initialStressMatrix, tStep, InitialStressMatrixAssembler(),	EModelDefaultEquationNumbering(), this->giveDomain(1));
		//initialStressMatrix->times(-1.0);

//#ifdef DEBUG
//	stiffnessMatrix->writeToFile("K-pds.dat");
//	initialStressMatrix->writeToFile("KG-pds.dat");
//#endif

		//int numEigv = 1;
		//FloatMatrix eigVec; eigVec.resize(this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()), numEigv);
		//FloatArray eigVal; eigVal.resize(numEigv);
		//GenEigvalSolverType eigSolver = GES_Eigen;
		//std::unique_ptr< SparseGeneralEigenValueSystemNM > nMethodST;
		//nMethodST.reset(classFactory.createGeneralizedEigenValueSolver(eigSolver, this->giveDomain(1), this));
		//nMethodST->solve(*stiffnessMatrix, *initialStressMatrix, eigVal, eigVec, rtolv, numEigv);
		//fprintf(outputStream, "# Eigenvalue found: %.4e \n", abs(eigVal.at(1)));
		//stiffnessMatrix->add(1 / abs(eigVal.at(1)), *initialStressMatrix); // without "-" if ->times(-1.0) is not used.

		//stiffnessMatrix->add(1, *initialStressMatrix); // without "-" if ->times(-1.0) is not used.
		std::unique_ptr< SparseMtrx > Kiter;
		Kiter = stiffnessMatrix->clone();
		Kiter->add(1, *initialStressMatrix);

//#ifdef DEBUG
//		stiffnessMatrix->writeToFile("Ke.dat");
//		Kiter->writeToFile("Kiter.dat");
//#endif
		// initialStressMatrix->times(eigVal.at(1));
		// initialStressMatrix->add(1.0, *stiffnessMatrix);
		// solve again
#ifdef VERBOSE
	OOFEM_LOG_INFO("\nSolving iteration %d ...\n", maxIter);
#endif
		// displacementVector.zero(); // not needed
		nMethod->solve(*Kiter, loadVector, displacementVector);

        } else {

        FloatArray feq( displacementVector.giveSize() );
        this->assembleVector( feq, tStep, MatrixProductAssembler( InitialStressMatrixAssembler() ),
            VM_Total, EModelDefaultEquationNumbering(), this->giveDomain( 1 ) );
        rhs.subtract( feq );
#ifdef VERBOSE
        OOFEM_LOG_INFO( "\nSolving iteration %d ...\n", maxIter );
#endif
        nMethod->solve( *stiffnessMatrix, rhs, displacementVector );
        }

		// check convergence on DISPLACEMENTS: ( u(i)^2 - u(i-1)^2 ) / u(i)^2
		newNorm = displacementVector.computeSquaredNorm();
		double toll = abs((newNorm - oldNorm) / newNorm);
		if (toll <= rtolv || maxIter >= 20) escape = true;
#ifdef VERBOSE
		OOFEM_LOG_INFO("\nCurrent displ. residual: %.2e \n\n", toll);
#endif
		// PDELTA end p-delta stiffness
	} while (escape == false);

    tStep->incrementStateCounter();            // update solution state counter
}


//SparseMtrx* add(SparseMtrx &n, SparseMtrx &m)
//{
//	SparseMtrx* res;
//
//	// for square matrices
//	if (n.giveNumberOfColumns() != m.giveNumberOfColumns()) {
//		OOFEM_ERROR("dimension of 'n' and 'm' mismatch");
//	}
//
//	res->add(1, m);
//	res->add(1, n);
//	return res;
//}

void PDeltaStatic :: saveContext(DataStream &stream, ContextMode mode)
//
// saves state variable - displacement vector
//
{
    contextIOResultType iores;
    //FILE *file = NULL;

    StructuralEngngModel::saveContext(stream, mode);

    if ( ( iores = displacementVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


void PDeltaStatic :: restoreContext(DataStream &stream, ContextMode mode)
//
// restore state variable - displacement vector
//
{
    contextIOResultType iores;

    StructuralEngngModel::restoreContext(stream, mode);

    if ( ( iores = displacementVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}

void
PDeltaStatic::terminate(TimeStep *tStep)
{
	StructuralEngngModel::terminate(tStep);
	//this->printReactionForces(tStep, 1, this->giveOutputStream());
	fflush(this->giveOutputStream());
}

void
PDeltaStatic::updateAfterStatic(TimeStep *tStep, Domain *domain)
{
	// Domain *domain = this->giveDomain(1);
	//tStep->setTime(0.);

	//if (requiresUnknownsDictionaryUpdate()) {
	//	for (auto &dman : domain->giveDofManagers()) {
	//		this->updateDofUnknownsDictionary(dman.get(), tStep);
	//	}
	//}

	for (auto &dman : domain->giveDofManagers()) {
		dman->updateYourself(tStep);
	}

	for (auto &elem : domain->giveElements()) {
		elem->updateInternalState(tStep);
		elem->updateYourself(tStep);
	}
}

void
PDeltaStatic :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}


int
PDeltaStatic :: estimateMaxPackSize(IntArray &commMap, DataStream &buff, int packUnpackType)
{
    int count = 0, pcount = 0;
    Domain *domain = this->giveDomain(1);

    if ( packUnpackType == 0 ) { ///@todo Fix this old ProblemCommMode__NODE_CUT value
        for ( int map: commMap ) {
            DofManager *dman = domain->giveDofManager( map );
            for ( Dof *dof: *dman ) {
                if ( dof->isPrimaryDof() && ( dof->__giveEquationNumber() ) ) {
                    count++;
                } else {
                    pcount++;
                }
            }
        }

        // --------------------------------------------------------------------------------
        // only pcount is relevant here, since only prescribed components are exchanged !!!!
        // --------------------------------------------------------------------------------

        return ( buff.givePackSizeOfDouble(1) * pcount );
    } else if ( packUnpackType == 1 ) {
        for ( int map: commMap ) {
            count += domain->giveElement( map )->estimatePackSize(buff);
        }

        return count;
    }

    return 0;
}

} // end namespace oofem
