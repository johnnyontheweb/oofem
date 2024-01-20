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

#include "sm/EngineeringModels/eigenvaluedynamic.h"
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
#include "simpleslavedof.h"
#include "slavedof.h"
#include "domain.h"
#include "element.h"
#include "node.h"
#include "rigidarmnode.h"
#include "unknownnumberingscheme.h"
#include "sm/Elements/lumpedmasselement.h"
#include "sm/Materials/structuralmaterial.h"
#include <math.h>

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#endif

#ifdef MEMSTR
#include <io.h>
#include <fcntl.h>
#endif

using namespace std;

namespace oofem {

/* this class is meant to give the plain numbering size for the matrix*/
class DummyNumbering : public UnknownNumberingScheme
{
protected:
    int neqs = 0;
    map<int, IntArray> dofs;

public:
    DummyNumbering( int n, map<int, IntArray> dofs ) :
        UnknownNumberingScheme(), neqs( n ), dofs( dofs ){};

    bool isDefault() const override { return false; }
    int giveRequiredNumberOfDomainEquation() const override { return neqs; }
    int giveDofEquationNumber( Dof *dof ) const override
    {
        int num          = dof->giveDofManager()->giveNumber();
        const auto dType = dof->giveDofType();
        return dofs.at( num ).at( dType );
    }
};


REGISTER_EngngModel( EigenValueDynamic );


EigenValueDynamic ::EigenValueDynamic( int i, EngngModel *master ) :
    EngngModel( i, master )
{
    numberOfSteps = 1;
    ndomains      = 1;
}


NumericalMethod *EigenValueDynamic ::giveNumericalMethod( MetaStep *mStep )
{
    if ( !nMethod ) {
        nMethod = classFactory.createGeneralizedEigenValueSolver( solverType, this->giveDomain( 1 ), this );
        if ( !nMethod ) {
            OOFEM_ERROR( "solver creation failed" );
        }
    }

    return nMethod.get();
}


void EigenValueDynamic ::initializeFrom( InputRecord &ir )
{
    //EngngModel::instanciateFrom (ir);

    IR_GIVE_FIELD( ir, numberOfRequiredEigenValues, _IFT_EigenValueDynamic_nroot );
    this->field = std::make_unique<EigenVectorPrimaryField>( this, 1, FT_Displacements, numberOfRequiredEigenValues );

    // numberOfSteps set artificially to numberOfRequiredEigenValues
    // in order to allow
    // use restoreContext function for different eigenValues
    // numberOfSteps = numberOfRequiredEigenValues;
    numberOfSteps = 1;

    IR_GIVE_FIELD( ir, rtolv, _IFT_EigenValueDynamic_rtolv );
    if ( rtolv < 1.e-12 ) {
        rtolv = 1.e-12;
    } else if ( rtolv > 0.01 ) {
        rtolv = 0.01;
    }

    int val = 1; // inverseit
    IR_GIVE_OPTIONAL_FIELD( ir, val, _IFT_EigenValueDynamic_stype );
    solverType = (GenEigvalSolverType)val;

    val = 0; //Default Skyline
    IR_GIVE_OPTIONAL_FIELD( ir, val, _IFT_EngngModel_smtype );
    sparseMtrxType = (SparseMtrxType)val;

    if ( solverType == GenEigvalSolverType::GES_Eigen )
        sparseMtrxType = SparseMtrxType::SMT_EigenSparse;
    suppressOutput = ir.hasField( _IFT_EngngModel_suppressOutput );

    if ( suppressOutput ) {
        printf( "Suppressing output.\n" );
    } else {
#ifdef MEMSTR
        outputStream = nullptr;
        FILE *source = classFactory.giveMemoryStream( "out" );
        int sourceFD = _open_osfhandle( (intptr_t)source, _O_APPEND );
        if ( sourceFD != -1 ) {
            outputStream = _fdopen( sourceFD, "a" );
        }
        if ( !( outputStream ) ) {
            // if not, write to file
#endif
            if ( ( outputStream = fopen( this->dataOutputFileName.c_str(), "w" ) ) == NULL ) {
                OOFEM_ERROR( "Can't open output file %s", this->dataOutputFileName.c_str() );
            }
#ifdef MEMSTR
            usestream = false;
        }
#endif

        fprintf( outputStream, "%s", PRG_HEADER );
        fprintf( outputStream, "\nStarting analysis on: %s\n", ctime( &this->startTime ) );
        fprintf( outputStream, "%s\n", simulationDescription.c_str() );
    }
}


int EigenValueDynamic ::giveUnknownDictHashIndx( ValueModeType mode, TimeStep *tStep )
{
    return tStep->giveNumber() % this->numberOfRequiredEigenValues;
}


double EigenValueDynamic ::giveUnknownComponent( ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof )
{
    return field->giveUnknownValue( dof, mode, tStep );
}


TimeStep *EigenValueDynamic ::giveNextStep()
{
    int istep                = giveNumberOfFirstStep();
    StateCounterType counter = 1;

    if ( currentStep ) {
        istep   = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = std ::move( currentStep );
    currentStep  = std::make_unique<TimeStep>( istep, this, 1, (double)istep, 0., counter );

    return currentStep.get();
}


void EigenValueDynamic ::solveYourself()
{
    this->timer.startTimer( EngngModelTimer ::EMTT_AnalysisTimer );

    TimeStep *tStep = this->giveNextStep();
    this->updateAttributes( this->giveCurrentMetaStep() );

    OOFEM_LOG_INFO( "Assembling stiffness and mass matrices\n" );

    FloatMatrix eigVec;
    std ::unique_ptr<SparseMtrx> stiffnessMatrix;
    std ::unique_ptr<SparseMtrx> massMatrix;

    const EModelDefaultEquationNumbering defNumbering;

    stiffnessMatrix = classFactory.createSparseMtrx( sparseMtrxType );
    stiffnessMatrix->buildInternalStructure( this, 1, defNumbering );

    massMatrix = classFactory.createSparseMtrx( sparseMtrxType );
    massMatrix->buildInternalStructure( this, 1, defNumbering );

    this->assemble( *stiffnessMatrix, tStep, TangentAssembler( TangentStiffness ), defNumbering, this->giveDomain( 1 ) );
    this->assemble( *massMatrix, tStep, MassMatrixAssembler(), defNumbering, this->giveDomain( 1 ) );

    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );
    OOFEM_LOG_INFO( "Solving ...\n" );
    nMethod->solve( *stiffnessMatrix, *massMatrix, eigVal, eigVec, rtolv, numberOfRequiredEigenValues );

    // cut off, but output at least one
    int nValid = 1;
    for ( int i = 1; i < eigVal.giveSize(); ++i ) {
        eigVal( i ) = abs( eigVal( i ) ); // force positive results even for ill-conditioned matrices
        if ( eigVal( i ) > 1e10 && i > 0 ) {
            nValid = i;
            break;
        }
        nValid = i + 1;
    }
    numberOfRequiredEigenValues = nValid; // better having a dedicated field
    eigVal.resize( numberOfRequiredEigenValues );
    eigVec.resizeWithData( eigVec.giveNumberOfRows(), numberOfRequiredEigenValues );

    Domain *domain  = this->giveDomain( 1 );
    const int nelem = domain->giveNumberOfElements();

    // matrix and array initialization
    totMass.resize( 6 ); // 3 translational masses + 3 rotational inertias
    centroid.resize( 3 );
    partFact.resize( numberOfRequiredEigenValues, 6 );
    massPart.resize( numberOfRequiredEigenValues, 6 );
    massPart.zero();
    // matrix of unit displacements for the 6 dofs
    FloatMatrix unitDisp( this->giveNumberOfDomainEquations( 1, defNumbering ), 6 );
    // auxillary vectors for mass normalization
    FloatArray tempCol( this->giveNumberOfDomainEquations( 1, defNumbering ) );
    FloatArray tempCol2( this->giveNumberOfDomainEquations( 1, defNumbering ) );

    // mass normalization
    for ( int i = 1; i <= numberOfRequiredEigenValues; ++i ) {
        tempCol.beColumnOf( eigVec, i );
        massMatrix->times( tempCol, tempCol2 );
        double m = tempCol.dotProduct( tempCol2 );
        if ( m != 0.0 ) m = 1 / sqrt( m );
        tempCol.times( m );
        eigVec.setColumn( tempCol, i );
        // if ( isnan( periods.at( i ) ) ) {
        //     printf( "stop" );
        // }
    }
    // eigVec has been normalized, now can be set into dofs' dictionaries
    this->field->updateAll( eigVec, defNumbering );

    // Temporary variables used in loops
    IntArray masterDofIDs, nodalEqArray, ids;

    static const DofIDItem dofIDs[] = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };

    bool warn          = false;
    const int nDofMans = domain->giveNumberOfDofManagers();

    FloatMatrix nodeCoords( 3, nDofMans ); // cache the coordinates
    std::map<int, IntArray> nodeActive; // cache which dofs are active

    FloatArray geomCenter( 3 );
    IntArray geomWeight( 3 );
    IntArray coordFilter( 3 ); // if no mass or free dof is defined in a direction, then centroid in that direction is undefined
    // cache dofs and node coordinates, unit displacements for translational dofs
    for ( std::unique_ptr<DofManager> &node : domain->giveDofManagers() ) {
        // no support yet for LCS
        Node *actualNode = dynamic_cast<Node *>( node.get() );
        if ( actualNode && actualNode->hasLocalCS() ) {
            if ( !warn ) {
                OOFEM_WARNING( "Nodes with LCS are unsupported, use at own risk." );
                warn = true;
            }
            continue;
        }
        if ( !node->giveNumberOfDofs() ) {
            continue;
        }

        IntArray activeDofs( 3 );
        if ( strcmp( node->giveClassName(), "Node" ) == 0 || strcmp( node->giveClassName(), "RigidArmNode" ) == 0 ) {
            for ( int iDof = 0; iDof < 3; ++iDof ) {
                DofIDItem dType = dofIDs[iDof];
                auto pos        = node->findDofWithDofId( dType );
                if ( pos == node->end() ) {
                    continue;
                }

                int eqN = 0;
                if ( ( *pos )->isPrimaryDof() ) {
                    eqN = ( *pos )->giveEquationNumber( defNumbering );
                } else {
                    // find the master dofmanager from the slave, get the equations from it.
                    // Plain Nodes have master and simpleslave dofs, rigid arm nodes have master and slave dofs
                    IntArray masterDofMans;
                    ( *pos )->giveMasterDofManArray( masterDofMans );

                    auto *masterMan = domain->giveDofManager( masterDofMans.at( 1 ) ); // there's only 1
                    auto masterDof  = masterMan->findDofWithDofId( dType );
                    eqN             = ( *masterDof )->giveEquationNumber( defNumbering );
                }

                activeDofs.at( dType ) = eqN > 0;
            }
        }

        // save coordinate and increment the weight
        const int num      = node->giveNumber();
        const auto &coords = node->giveCoordinates();
        nodeCoords.setColumn( coords, num );
        for ( int iDof = 1; iDof <= 3; ++iDof ) {
            if ( activeDofs.at( iDof ) ) {
                geomCenter.at( iDof ) += coords.at( iDof );
                geomWeight.at( iDof ) += activeDofs.at( iDof );
            }
        }
        nodeActive[num] = std::move( activeDofs );
    }

    FloatArray massPos( 3 ); // stores the product between mass and position
    // get mass contribution to the center of mass from elements
    for ( int ielem = 1; ielem <= nelem; ++ielem ) {
        Element *element = domain->giveElement( ielem );

        // we only support masses from lumpedmasselements.
        if ( strcmp( element->giveClassName(), "LumpedMassElement" ) != 0 ) {
            if ( element->giveMaterialNumber() && element->giveIntegrationRulesArray().size() && element->giveMaterial()->give( 'd', element->giveIntegrationRulesArray()[0]->getIntegrationPoint( 0 ) ) != 0 ) {
                if ( !warn ) {
                    OOFEM_WARNING( "Only masses from LumpedMassElements are supported." );
                    warn = true;
                }
            }
            continue;
        }

        LumpedMassElement *massElement = dynamic_cast<LumpedMassElement *>( element );
        if ( !massElement ) {
            continue;
        }
        const int dofMan = massElement->giveDofManArray().at( 1 );
        IntArray dofMask;
        massElement->giveElementDofIDMask( dofMask );
        FloatMatrix mat;
        massElement->computeLumpedMassMatrix( mat, tStep );

        if ( mat.isNotEmpty() ) {
            for ( auto iDof : dofMask ) {
                if ( iDof >= D_u && iDof <= D_w ) {
                    const int idx = dofMask.findFirstIndexOf( iDof );
                    if ( idx && nodeActive[dofMan].at( iDof ) ) {
                        const double m = mat.at( idx, idx );
                        massPos.at( iDof ) += nodeCoords.at( iDof, dofMan ) * m;
                        totMass.at( iDof ) += m;
                    }
                }
            }
        }
    }

    // center of geometry
    for ( int i = 1; i <= 3; ++i ) {
        if ( geomWeight.at( i ) ) {
            geomCenter.at( i ) /= geomWeight.at( i );
            coordFilter.at( i ) |= 1;
        } else {
            OOFEM_WARNING( "No free dof in direction %d, results may be affected", i );
        }
    }
    // center of mass
    for ( int i = 1; i <= 3; ++i ) {
        if ( totMass.at( i ) > 0 ) {
            centroid.at( i ) = massPos.at( i ) / totMass.at( i );
            coordFilter.at( i ) |= 1;
        } else {
            // OOFEM_WARNING( "No mass in direction %d, results may be affected", i ); // no warning for this
            centroid.at( i ) = geomCenter.at( i );
        }
    }

    // create plain numbering system
    std::map<int, IntArray> dofMansEqns; // association between dofmanagers and plain numbering. BC dofs have equation set to 0.
    std::map<int, IntArray> intDofMansEqns; // equation numbers with default numbering for internal dofs
    int totalDofs = this->giveNumberOfDomainEquations( 1, defNumbering ); // initialize the counting from the default numbering

    // fill the arrays with dof numbers, equations etc
    // also internal dofs are parsed, so the equation numbers will match for both numbering systems
    // for dofs defined in both. the plain numbering will have additional dofs (the slave ones) but
    // their numbering is the continuation
    for ( std::unique_ptr<DofManager> &node : domain->giveDofManagers() ) {
        // no support yet for LCS
        Node *actualNode = dynamic_cast<Node *>( node.get() );
        if ( actualNode && actualNode->hasLocalCS() ) {
            continue;
        }

        IntArray eqns( 6 );
        if ( strcmp( node->giveClassName(), "Node" ) == 0 || strcmp( node->giveClassName(), "RigidArmNode" ) == 0 ) {
            for ( int iDof = 0; iDof < 6; ++iDof ) {
                DofIDItem dType = dofIDs[iDof];
                auto pos        = node->findDofWithDofId( dType );
                if ( pos == node->end() ) {
                    continue;
                }

                if ( ( *pos )->isPrimaryDof() ) {
                    // use the default numbering for the master dof
                    eqns.at( dType ) = ( *pos )->giveEquationNumber( defNumbering );
                } else {
                    // add a dummy equation number for the slave dof
                    eqns.at( dType ) = ++totalDofs;
                }
            }
        }
        dofMansEqns[actualNode->giveNumber()] = std::move( eqns );
    }
    // do it for internal dofs as well
    warn = false;
    for ( int ielem = 1; ielem <= nelem; ++ielem ) {
        Element *element = domain->giveElement( ielem );
        if ( element->giveNumberOfInternalDofManagers() == 0 ) {
            continue;
        }

        // we only support end releases on 3d beams for the moment. other internal dof managers should be skipped.
        if ( strcmp( element->giveClassName(), "Beam3d" ) != 0 ) {
            if ( !warn ) {
                OOFEM_WARNING( "Only internal dof managers of Beam3d elements are supported. Skipping." );
                warn = true;
            }
            continue;
        }

        IntArray eqArray;
        // retrieve internal dof managers and location array
        for ( int i = 1; i <= element->giveNumberOfInternalDofManagers(); ++i ) {
            DofManager *intDofMan = element->giveInternalDofManager( i );
            if ( !intDofMan ) continue; // you may never know...

            element->giveInternalDofManDofIDMask( i, ids );
            intDofMan->giveLocationArray( ids, nodalEqArray, defNumbering );
            eqArray.followedBy( nodalEqArray );
            intDofMan->giveMasterDofIDArray( ids, masterDofIDs );
        }
        intDofMansEqns[element->giveNumber()] = std::move( eqArray );
    }

    FloatMatrix plainUnitDisps( totalDofs, 6 ); // vector for unit displacements expressed in plain numbering
    FloatMatrix plainEigVecs( totalDofs, numberOfRequiredEigenValues ); // vector for eigen vectors expressed in plain numbering

    IntArray dofArray{ 1, 2, 3, 4, 5, 6 };

    warn = false;
    // assemble into the plain numbering eigenvectors
    // and plain numbering unit displacement vectors
    for ( std::unique_ptr<DofManager> &node : domain->giveDofManagers() ) {
        // no support yet for LCS
        Node *actualNode = dynamic_cast<Node *>( node.get() );
        if ( actualNode && actualNode->hasLocalCS() ) {
            if ( !warn ) {
                OOFEM_WARNING( "Nodes with LCS are unsupported, use at own risk." );
                warn = true;
            }
            continue;
        }
        if ( !node->giveNumberOfDofs() ) {
            continue;
        }

        const int num   = node->giveNumber();
        const auto &loc = dofMansEqns[num];
        for ( int dN = 1; dN <= numberOfRequiredEigenValues; ++dN ) {
            this->setActiveVector( dN ); // we use time as intrinsic eigen value index

            FloatArray charVec( 6 );
            for ( auto dof : *node ) {
                // fetch the unknown from the eigen field. see giveUnknownComponent
                charVec.at( dof->giveDofID() ) = dof->giveUnknown( VM_Total, tStep );
            }
            plainEigVecs.assemble( charVec, loc, { dN } );
        }

        const FloatArray diff = node->giveCoordinates() - centroid;

        FloatMatrix dofTransf( 6, 6 );
        for ( int iDof = 1; iDof <= 6; ++iDof ) {
            dofTransf.at( iDof, iDof ) = 1.0;

            // translational components due to rotation about the center of mass of the model
            switch ( iDof ) {
            case 4:
                dofTransf.at( 2, 4 ) = -diff.at( 3 );
                dofTransf.at( 3, 4 ) = diff.at( 2 );
                break;
            case 5:
                dofTransf.at( 1, 5 ) = diff.at( 3 );
                dofTransf.at( 3, 5 ) = -diff.at( 1 );
                break;
            case 6:
                dofTransf.at( 1, 6 ) = -diff.at( 2 );
                dofTransf.at( 2, 6 ) = diff.at( 1 );
                break;
            default:
                break;
            }
        }
        plainUnitDisps.assemble( dofTransf, loc, dofArray );
    }

    // mass matrix for plain numbering
    std::unique_ptr<SparseMtrx> plainMassMatrix = classFactory.createSparseMtrx( sparseMtrxType );
    DummyNumbering dummyNumbering( totalDofs, dofMansEqns );
    plainMassMatrix->buildInternalStructure( this, 1, dummyNumbering );

    // create mass matrix with plain numbering
    for ( int ielem = 1; ielem <= nelem; ++ielem ) {
        Element *element = domain->giveElement( ielem );

        LumpedMassElement *massElement = dynamic_cast<LumpedMassElement *>( element );
        if ( !massElement ) {
            continue;
        }
        const int dofMan = massElement->giveDofManagerNumber( 1 );
        IntArray dofMask;
        massElement->giveElementDofIDMask( dofMask );
        FloatMatrix mat;
        massElement->computeLumpedMassMatrix( mat, tStep );

        if ( mat.isNotEmpty() ) {
            IntArray loc( 6 );
            const IntArray &eqns = dofMansEqns[dofMan];
            for ( int iLoc = 1; iLoc <= dofMask.giveSize(); ++iLoc ) {
                const int iDof = dofMask.at( iLoc );
                loc.at( iDof ) = eqns.at( iDof );
            }

            if ( plainMassMatrix->assemble( loc, mat ) == 0 ) {
                OOFEM_ERROR( "sparse matrix assemble error" );
            }
        }
    }

    // force the creation of the entries from the triplets
    plainMassMatrix->assembleBegin();
    plainMassMatrix->assembleEnd();

    FloatMatrix tempMat( totalDofs, 6 ); // matrix for temporary results

    for ( int i = 1; i <= 6; ++i ) {
        tempCol.beColumnOf( plainUnitDisps, i );
        plainMassMatrix->times( tempCol, tempCol2 ); // now tempCol2 has only the masses pertaining the i-th direction
        tempMat.setColumn( tempCol2, i );
        totMass.at( i ) = tempCol.dotProduct( tempCol2 );
    }

    //
    // calculate participation factors and mass participation
    //

    // Mass renormalization
    // In presence of rigid arms in the model, the mass matrix will contain inertias wrt the master nodes as this is how dof transformation works.
    // If master nodes aren't placed barycentrically among the related nodes, it means the inertias will be higher than with barycentric placement.
    // Our calculation of the centroid is instead based on the actual geometric placement of lumped masses, being this less misleading and easier to do.
    // This in turn will compromise the computation of participation factors and masses when unit rotation vectors are determined about our centroid
    // while eigenvector refer to the master node inertias.
    // Eigenvectors are therefore renormalized on the mass matrix done with the plain numbering, to adjust the results of the coming calculations.
    // The recommendation is however to place master nodes in the barycenter of the nodes they constrain, as this is only a hack.
    // If the model has barycentric master nodes, then this renormalization has no effect.
    for ( int i = 1; i <= numberOfRequiredEigenValues; ++i ) {
        tempCol.beColumnOf( plainEigVecs, i );
        plainMassMatrix->times( tempCol, tempCol2 );
        double m = tempCol.dotProduct( tempCol2 );
        if ( m != 0.0 ) {
            m = 1 / sqrt( m );
        }
        tempCol.times( m );
        plainEigVecs.setColumn( tempCol, i );
    }

    // participation factors
    partFact.beTProductOf( plainEigVecs, tempMat );

    // mass participation ratios
    for ( int i = 1; i <= numberOfRequiredEigenValues; ++i ) {
        for ( int j = 1; j <= 6; ++j ) {
            if ( totMass.at( j ) > 1e-10 ) {
                massPart.at( i, j ) = pow( partFact.at( i, j ), 2 ) / totMass.at( j );
            }
        }
    }

    //
    // zero matrix
    //
    stiffnessMatrix.reset( NULL );
    massMatrix.reset( NULL );

    this->terminate( tStep );

    double steptime = this->giveSolutionStepTime();
    OOFEM_LOG_INFO( "EngngModel info: user time consumed by solution: %.2fs\n", steptime );
}


void EigenValueDynamic ::updateYourself( TimeStep *tStep )
{
}


void EigenValueDynamic ::doStepOutput( TimeStep *tStep )
{
    if ( !suppressOutput ) {
        this->printOutputAt( this->giveOutputStream(), tStep );
        fflush( this->giveOutputStream() );
    }

    for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        // export using export manager
        tStep->setTime( (double)i ); // we use time as intrinsic eigen value index
        tStep->setNumber( i );
        exportModuleManager.doOutput( tStep );
    }
}


void EigenValueDynamic ::printOutputAt( FILE *file, TimeStep *tStep )
{
    Domain *domain     = this->giveDomain( 1 );
    FILE *outputStream = this->giveOutputStream();

    // print loadcase header
    fprintf( file, "\nOutput for Eigen analysis \n\n", 1.0 );
    // print eigen values on output
    fprintf( file, "\n\nEigen Values (Omega^2) are:\n-----------------\n" );

    for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        fprintf( outputStream, "%15.8e ", eigVal.at( i ) );
        if ( ( i % 5 ) == 0 ) {
            fprintf( outputStream, "\n" );
        }
    }


    fprintf( outputStream, "\n\n\nCenter of Mass:\n-----------------\n\tX\t|\tY\t|\tZ\n" );
    for ( int i = 1; i <= centroid.giveSize(); ++i ) {
        fprintf( outputStream, "%15.8e ", centroid.at( i ) );
    }

    fprintf( outputStream, "\n" );

    fprintf( outputStream, "\n\nParticipation Factors are:\n-----------------\n\tDx\t|\tDy\t|\tDz\t|\tRx\t|\tRy\t|\tRz\n" );
    for ( int i = 1; i <= partFact.giveNumberOfRows(); ++i ) {
        for ( int j = 1; j <= partFact.giveNumberOfColumns(); ++j ) {
            fprintf( outputStream, "%15.8e ", partFact.at( i, j ) );
        }
        fprintf( outputStream, "\n" );
    }

    fprintf( outputStream, "\n\nTotal Masses are:\n-----------------\n\tDx\t|\tDy\t|\tDz\t|\tRx\t|\tRy\t|\tRz\n" );
    for ( int i = 1; i <= totMass.giveSize(); ++i ) {
        fprintf( outputStream, "%15.8e ", totMass.at( i ) );
    }

    fprintf( outputStream, "\n" );

    fprintf( outputStream, "\n\nMass Ratios are:\n-----------------\n\tDx\t|\tDy\t|\tDz\t|\tRx\t|\tRy\t|\tRz\n" );
    for ( int i = 1; i <= massPart.giveNumberOfRows(); ++i ) {
        for ( int j = 1; j <= massPart.giveNumberOfColumns(); ++j ) {
            fprintf( outputStream, "%15.8e ", massPart.at( i, j ) );
        }
        fprintf( outputStream, "\n" );
    }

    fprintf( outputStream, "\n\n" );

    for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        fprintf( outputStream, "\nOutput for eigen value no.  %.3e \n", (double)i );
        fprintf( outputStream,
            "Printing eigen vector no. %d, corresponding eigen value is %15.8e\n\n",
            i, eigVal.at( i ) );
        tStep->setTime( (double)i ); // we use time as intrinsic eigen value indextStep->setNumber(i);
        tStep->setNumber( i );

        if ( this->requiresUnknownsDictionaryUpdate() ) {
            for ( auto &dman : domain->giveDofManagers() ) {
                this->updateDofUnknownsDictionary( dman.get(), tStep );
            }
        }


        for ( auto &dman : domain->giveDofManagers() ) {
            dman->updateYourself( tStep );
            dman->printOutputAt( outputStream, tStep );
        }
    }

    fflush( this->giveOutputStream() );

    double utsec = this->timer.getUtime( EngngModelTimer::EMTT_AnalysisTimer );
    fprintf( file, "\nUser time consumed by solution step: %.3f [s]\n\n", utsec );
}


void EigenValueDynamic ::saveContext( DataStream &stream, ContextMode mode )
{
    EngngModel ::saveContext( stream, mode );

    contextIOResultType iores;
    if ( ( iores = eigVal.storeYourself( stream ) ) != CIO_OK ) {
        THROW_CIOERR( iores );
    }

    this->field->saveContext( stream );
}


void EigenValueDynamic ::restoreContext( DataStream &stream, ContextMode mode )
{
    EngngModel ::restoreContext( stream, mode );

    contextIOResultType iores;
    if ( ( iores = eigVal.restoreYourself( stream ) ) != CIO_OK ) {
        THROW_CIOERR( iores );
    }

    this->field->restoreContext( stream );
}


void EigenValueDynamic ::setActiveVector( int i )
{
    this->activeVector = i;
    if ( activeVector > numberOfRequiredEigenValues ) {
        activeVector = numberOfRequiredEigenValues;
    }

    this->giveCurrentStep()->setNumber( activeVector );
    this->giveCurrentStep()->setTime( (double)activeVector );
}

} // end namespace oofem
