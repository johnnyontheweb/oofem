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

#include "sm/EngineeringModels/responsespectrum.h"
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
#include "slavedof.h"
#include "simpleslavedof.h"
#include "domain.h"
#include "element.h"
#include "node.h"
#include "unknownnumberingscheme.h"
#include "function.h"
#include "activebc.h"
#include "integrationrule.h"
#include "gausspoint.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Elements/lumpedmasselement.h"
#include "outputmanager.h"
#include "dynamicdatareader.h"
#include "dynamicinputrecord.h"
#include "inputrecord.h"
#include "sm/Elements/structuralelement.h"
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
    vector<IntArray> dofs;

public:
    DummyNumbering( int n, vector<IntArray> dofs ) :
        UnknownNumberingScheme(), neqs( n ), dofs(dofs){};

    bool isDefault() const override { return false; }
    int giveRequiredNumberOfDomainEquation() const override { return neqs; }
    int giveDofEquationNumber( Dof *dof ) const override { 
        int num = dof->giveDofManager()->giveNumber() - 1;  // make it zero based for std container
        const auto dType = dof->giveDofType();
        return dofs.at( num ).at( dType );
    }
};


REGISTER_EngngModel( ResponseSpectrum );

NumericalMethod *ResponseSpectrum::giveNumericalMethod( MetaStep *mStep )
{
    if ( !nMethod ) {
        nMethod = classFactory.createGeneralizedEigenValueSolver( solverType, this->giveDomain( 1 ), this );
        if ( !nMethod ) {
            OOFEM_ERROR( "solver creation failed" );
        }
    }

    return nMethod.get();
}

void ResponseSpectrum::initializeFrom( InputRecord &ir )
{
    // EngngModel::instanciateFrom (ir);

    IR_GIVE_FIELD( ir, numberOfRequiredEigenValues, _IFT_ResponseSpectrum_nroot );
    this->field = std::make_unique<EigenVectorPrimaryField>( this, 1, FT_Displacements, numberOfRequiredEigenValues );

    // numberOfSteps set artificially to numberOfRequiredEigenValues
    // in order to allow
    // use restoreContext function for different eigenValues
    // numberOfSteps = numberOfRequiredEigenValues;
    numberOfSteps = 1;

    IR_GIVE_FIELD( ir, rtolv, _IFT_ResponseSpectrum_rtolv );
    if ( rtolv < 1.e-12 ) {
        rtolv = 1.e-12;
    }

    if ( rtolv > 0.01 ) {
        rtolv = 0.01;
    }

    int val = 1;
    IR_GIVE_OPTIONAL_FIELD( ir, val, _IFT_ResponseSpectrum_stype );
    solverType = (GenEigvalSolverType)val;

    val = 0; // Default Skyline
    IR_GIVE_OPTIONAL_FIELD( ir, val, _IFT_EngngModel_smtype );
    sparseMtrxType = (SparseMtrxType)val;

    if ( solverType == GenEigvalSolverType::GES_Eigen ) {
        sparseMtrxType = SparseMtrxType::SMT_EigenSparse;
        linStype       = LinSystSolverType::ST_EigenLib;
    }

    IR_GIVE_FIELD( ir, val, _IFT_ResponseSpectrum_func );
    func = (int)val; // we'll check in postInitialize whether this id exists or not

    IR_GIVE_FIELD( ir, dir, _IFT_ResponseSpectrum_dir );
    // dir = (int)val;
    if ( !( dir.giveSize() ) ) {
        OOFEM_ERROR( "No direction vector set." );
    }
    if ( dir.giveSize() > 3 ) {
        OOFEM_WARNING( "more than 3 vector components set. Trimming direction vector" );
        dir.resizeWithValues( 3 );
    } else if ( dir.giveSize() < 3 ) {
        OOFEM_WARNING( "less than 3 vector components set. Setting the remaining to zero" );
        dir.resizeWithValues( 3 );
    }
    dir.normalize();

    val = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, val, _IFT_ResponseSpectrum_modalCombo );
    modalCombo = (RSpecComboType)val;

    double damp = 0.05; // default damping ratio
    IR_GIVE_OPTIONAL_FIELD( ir, damp, _IFT_ResponseSpectrum_damp );
    csi = damp;

    suppressOutput = ir.hasField( _IFT_EngngModel_suppressOutput );

    if ( suppressOutput ) {
        // printf("Suppressing output.\n");
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


void ResponseSpectrum::postInitialize()
{
    EngngModel::postInitialize();

    // we check whether the spectrum function exists or not.
    Domain *d   = this->giveDomain( 1 );
    Function *f = d->giveFunction( func );

    if ( f == NULL ) OOFEM_ERROR( "Invalid function given" );
}


double ResponseSpectrum::giveUnknownComponent( ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof )
// returns unknown quantity like displacement, eigenvalue.
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
#ifdef DEBUG
    if ( eq == 0 ) {
        OOFEM_ERROR( "invalid equation number" );
    }
#endif

    switch ( mode ) {
    case VM_Total: // EigenVector
    case VM_Incremental:
        if ( tStep->giveIntrinsicTime() == 0.0 ) {
            return combDisps.at( eq ); // , (int)tStep->giveTargetTime());
        } else {
            return dummyDisps.at( eq );
        }
    case VM_Velocity:
        return 0.;
    case VM_Acceleration:
        return 0.;
    default:
        OOFEM_WARNING( "Unknown is of undefined type for this problem" );
    }

    return 0.;
}


TimeStep *ResponseSpectrum::giveNextStep()
{
    int istep                = giveNumberOfFirstStep();
    StateCounterType counter = 1;

    if ( currentStep ) {
        istep   = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = std::move( currentStep );
    currentStep.reset( new TimeStep( istep, this, 1, (double)istep, 0., counter ) );

    return currentStep.get();
}


// gets the Spectral acceleration from the function given in input
double ResponseSpectrum::calcSpectrumOrdinate( double period )
{
    Function *f = this->giveDomain( 1 )->giveFunction( this->func );
    return f->evaluateAtTime( period );
}

// forward declaration
void addMultiply( map<int, map<int, map<int, map<string, FloatArray> > > > &answer, map<int, map<int, map<int, map<string, FloatArray> > > > &src, map<int, map<int, map<int, map<string, FloatArray> > > > &src2, double fact = 1.0 );
void calcRoot( map<int, map<int, map<int, map<string, FloatArray> > > > &answer );
void addMultiply( map<int, map<string, FloatArray> > &answer, map<int, map<string, FloatArray> > &src, map<int, map<string, FloatArray> > &src2, double fact = 1.0 );
void calcRoot( map<int, map<string, FloatArray> > &answer );

void ResponseSpectrum::solveYourself()
{
    this->timer.startTimer( EngngModelTimer ::EMTT_AnalysisTimer );

    TimeStep *tStep = this->giveNextStep();
    this->updateAttributes( this->giveCurrentMetaStep() );

    OOFEM_LOG_INFO( "Assembling stiffness and mass matrices\n" );

    stiffnessMatrix = classFactory.createSparseMtrx( sparseMtrxType );
    stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

    massMatrix = classFactory.createSparseMtrx( sparseMtrxType );
    massMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

    this->assemble( *stiffnessMatrix, tStep, TangentAssembler( TangentStiffness ), EModelDefaultEquationNumbering(), this->giveDomain( 1 ) );
    this->assemble( *massMatrix, tStep, MassMatrixAssembler(), EModelDefaultEquationNumbering(), this->giveDomain( 1 ) );

    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );
    OOFEM_LOG_INFO( "Solving ...\n" );
    nMethod->solve( *stiffnessMatrix, *massMatrix, eigVal, eigVec, rtolv, numberOfRequiredEigenValues );

    // cut off, but output at least one
    int nValid = 1;
    for ( int i = 1; i < eigVal.giveSize(); ++i ) {
        if ( eigVal( i ) > 1e10 && i > 0 ) {
            nValid = i;
            break;
        }
        nValid = i + 1;
    }
    numberOfRequiredEigenValues = nValid; // better having a dedicated field
    eigVal.resize( numberOfRequiredEigenValues );
    eigVec.resizeWithData( eigVec.giveNumberOfRows(), numberOfRequiredEigenValues );

    this->field->updateAll( eigVec, EModelDefaultEquationNumbering() );

    FloatMatrix unitDisp;
    FloatArray tempCol;
    FloatArray tempCol2;

    Domain *domain = this->giveDomain( 1 );
    IntArray dofIDArry, loc;
    dofIDArry = domain->giveDefaultNodeDofIDArry();
    int nelem = domain->giveNumberOfElements();

    // matrix and array initialization
    totMass.resize( 6 ); // 3 translational masses + 3 rotational inertias
    centroid.resize( 3 );
    partFact.resize( numberOfRequiredEigenValues, 6 );
    massPart.resize( numberOfRequiredEigenValues, 6 );
    // matrix of unit displacements for the 6 dofs
    unitDisp.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ), 6 );
    // auxillary vectors for mass normalization
    tempCol.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    tempCol2.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    periods.resize( numberOfRequiredEigenValues );
    combReactions.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultPrescribedEquationNumbering() ) );
    ;
    combDisps.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    rhos.resize( numberOfRequiredEigenValues, numberOfRequiredEigenValues );

    // mass normalization
    for ( int i = 1; i <= numberOfRequiredEigenValues; ++i ) {
        tempCol.beColumnOf( eigVec, i );
        massMatrix->times( tempCol, tempCol2 );
        double m = tempCol.dotProduct( tempCol2 );
        if ( m != 0.0 ) m = 1 / sqrt( m );
        tempCol.times( m );
        eigVec.setColumn( tempCol, i );
        periods.at( i ) = 2 * M_PI / sqrt( eigVal.at( i ) );
        // if ( isnan( periods.at( i ) ) ) {
        //     printf( "stop" );
        // }
    }
    // eigVec has been normalized

    for ( int i = 1; i <= numberOfRequiredEigenValues; ++i )
        for ( int j = 1; j <= numberOfRequiredEigenValues; ++j ) {
            double beta     = periods.at( i ) / periods.at( j );
            rhos.at( i, j ) = 8 * pow( this->csi, 2.0 ) * pow( beta, 1.5 ) / ( 1.0 + beta ) / ( pow( 1 - beta, 2.0 ) + 4 * pow( this->csi, 2.0 ) * beta );
            // if(isnan( rhos.at( i, j ))) {
            //     printf( "stop" );
            // }
        }

    IntArray masterDofIDs, nodalEqArray, ids;
    IntArray eqArray;
    IntArray dofManArray;

    FloatMatrix tempMat, tempMat2;
    tempMat.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ), 6 );
    tempMat2.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ), 6 );

    static const DofIDItem dofids[] = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };

    bool warn = false;
    const EModelDefaultEquationNumbering defNumbering;
    const int nDofMans = domain->giveNumberOfDofManagers();
    FloatMatrix nodeCoords( 3, nDofMans );
    std::vector<IntArray> nodeActive( nDofMans );

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

        const int num = node->giveNumber();

        const auto &coords = node->giveCoordinates();
        IntArray activeDofs( 3 );
        if ( strcmp( node->giveClassName(), "Node" ) == 0 || strcmp( node->giveClassName(), "RigidArmNode" ) == 0 ) {
            for ( int iDof = 0; iDof < 3; ++iDof ) {
                DofIDItem dType = dofids[iDof];
                auto pos        = node->findDofWithDofId( dType );
                if ( pos == node->end() ) {
                    continue;
                }

                int eqN = 0;
                if ( ( *pos )->isPrimaryDof() ) {
                    eqN = ( *pos )->giveEquationNumber( defNumbering );
                    // if ( eqN ) {
                    //     unitDisp.at( eqN, dType ) = 1.0;
                    // }
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
        nodeActive[num - 1] = activeDofs;

        // save coordinate and increment the weight
        nodeCoords.setColumn( coords, num );
        for ( int iDof = 1; iDof <= 3; ++iDof ) {
            if ( activeDofs.at( iDof ) ) {
                geomCenter.at( iDof ) += coords.at( iDof );
                geomWeight.at( iDof ) += activeDofs.at( iDof );
            }
        }
    }

    warn = false;
    LumpedMassVectorAssembler lma;
    FloatMatrix R;
    FloatArray massPos( 3 ); // stores the product between mass and position
    // then from internaldof managers
    for ( int ielem = 1; ielem <= nelem; ++ielem ) {
        Element *element = domain->giveElement( ielem );

        // we only support masses from lumpedmasselements.
        if ( strcmp( element->giveClassName(), "LumpedMassElement" ) != 0 ) {
            if ( !warn ) {
                OOFEM_WARNING( "Only masses from LumpedMassElements are suppported." );
                warn = true;
            }
            continue;
        }

        LumpedMassElement *massElement = dynamic_cast<LumpedMassElement *>( element );
        FloatArray m;
        const int dofMan = massElement->giveDofManArray().at( 1 );
        lma.vectorFromElement( m, *massElement, tStep, VM_Unknown );

        IntArray dofMask;
        massElement->giveElementDofIDMask( dofMask );

        for ( auto iDof : dofMask ) {
            if ( iDof >= D_u && iDof <= D_w ) {
                const int idx = dofMask.findFirstIndexOf( iDof );
                if ( idx && nodeActive[dofMan - 1].at( iDof ) ) {
                    massPos.at( iDof ) += nodeCoords.at( iDof, dofMan ) * m.at( idx );
                    totMass.at( iDof ) += m.at( idx );
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
            OOFEM_WARNING( "No mass in direction %d, results may be affected", i );
            centroid.at( i ) = geomCenter.at( i );
        }
    }

    // create plain numbering system
    std::vector<IntArray> dofMansEqns( nDofMans );
    std::map<int, IntArray> intDofMansEqns;
    int totalDofs = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

    for ( std::unique_ptr<DofManager> &node : domain->giveDofManagers() ) {
        // no support yet for LCS
        Node *actualNode = dynamic_cast<Node *>( node.get() );
        if ( actualNode && actualNode->hasLocalCS() ) {
            continue;
        }

        IntArray eqns( 6 );
        if ( strcmp( node->giveClassName(), "Node" ) == 0 || strcmp( node->giveClassName(), "RigidArmNode" ) == 0 ) {
            for ( int iDof = 0; iDof < 6; ++iDof ) {
                DofIDItem dType = dofids[iDof];
                auto pos        = node->findDofWithDofId( dType );
                if ( pos == node->end() ) {
                    continue;
                }

                int eqN = 0;
                if ( ( *pos )->isPrimaryDof() ) {
                    eqN              = ( *pos )->giveEquationNumber( defNumbering );
                    eqns.at( dType ) = eqN;
                } else {
                    // add a dummy equation number for the slave dof
                    eqns.at( dType ) = ++totalDofs;
                }
            }
        }
        dofMansEqns[actualNode->giveNumber() - 1] = std::move( eqns );
    }
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

        eqArray.clear();
        // retrieve internal dof managers and location array
        for ( int i = 1; i <= element->giveNumberOfInternalDofManagers(); ++i ) {
            DofManager *intDofMan = element->giveInternalDofManager( i );
            if ( !intDofMan ) continue; // you may never know...

            element->giveInternalDofManDofIDMask( i, ids );
            intDofMan->giveLocationArray( ids, nodalEqArray, EModelDefaultEquationNumbering() );
            eqArray.followedBy( nodalEqArray );
            intDofMan->giveMasterDofIDArray( ids, masterDofIDs );
        }
        intDofMansEqns[element->giveNumber()] = std::move( eqArray );
    }

    std::unique_ptr<SparseMtrx> plainMassMatrix = classFactory.createSparseMtrx( sparseMtrxType );
    DummyNumbering dummyNumbering(totalDofs, dofMansEqns);
    plainMassMatrix->buildInternalStructure( this, 1, dummyNumbering);

    // todo: clean up all duplication wtf
    LumpedMassVectorAssembler lma;
    FloatMatrix R;
    FloatArray massPos(3); // stores the product between mass and position
    // then from internaldof managers
    for (int ielem = 1; ielem <= nelem; ++ielem) {
        Element* element = domain->giveElement(ielem);

        // we only support masses from lumpedmasselements.
        if (strcmp(element->giveClassName(), "LumpedMassElement") != 0) {
            if (!warn) {
                OOFEM_WARNING("Only masses from LumpedMassElements are suppported.");
                warn = true;
            }
            continue;
        }

        LumpedMassElement* massElement = dynamic_cast<LumpedMassElement*>(element);
        FloatArray m;
        const int dofMan = massElement->giveDofManagerNumber(1);
        lma.vectorFromElement(m, *massElement, tStep, VM_Unknown);

        IntArray dofMask;
        massElement->giveElementDofIDMask(dofMask);
        FloatMatrix mat, R;
        mat.beDiagonal(m);

        if (mat.isNotEmpty()) {
            const IntArray loc = dofMansEqns[dofMan];
            ///@todo This rotation matrix is not flexible enough.. it can only work with full size matrices and doesn't allow for flexibility in the matrixassembler.
            if (element->giveRotationMatrix(R)) {
                mat.rotatedWith(R);
            }

            if (plainMassMatrix->assemble(loc, mat) == 0) {
                OOFEM_ERROR("sparse matrix assemble error");
            }
        }
    }

    //
    // create unit displacement vectors
    //
    // first from nodes themselves, and get the coordinates to compute the center of mass
    for ( std::unique_ptr<DofManager> &node : domain->giveDofManagers() ) {
        // no support yet for LCS
        Node *actualNode = dynamic_cast<Node *>( node.get() );
        if ( actualNode && actualNode->hasLocalCS() ) {
            continue;
        }
        if ( !node->giveNumberOfDofs() ) continue;

        IntArray mstrDofs, locArr;
        node->givePrimaryDofs( mstrDofs );
        node->giveLocationArray( mstrDofs, locArr, EModelDefaultEquationNumbering() );

        // search for our dofs in there
        for ( int myDofIndex = 1; myDofIndex <= locArr.giveSize(); ++myDofIndex ) {
            int dType = mstrDofs.at( myDofIndex );
            int eqN   = locArr.at( myDofIndex );

            if ( ( dType >= D_u ) && ( dType <= D_w ) && eqN ) {
                // save unit displacement and coordinate
                unitDisp.at( eqN, dType ) = 1.0;
                // tempMat2.at( eqN, dType ) = node->giveCoordinate( dType );
            }
        }
    } // end of search among nodes

    warn = false;
    // unit displacement vectors from nodes were created before for translational dofs
    // now from internal dof managers
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

        // use the lcs transform to work out which ghost dofs are affected by global displacement
        FloatMatrix GtoL;
        element->computeGtoLRotationMatrix( GtoL );
        FloatMatrix lcs;
        element->giveLocalCoordinateSystem( lcs );

        eqArray.clear();

        // the following may be simplified.
        // retrieve internal dof managers and location array
        for ( int i = 1; i <= element->giveNumberOfInternalDofManagers(); ++i ) {
            DofManager *intDofMan = element->giveInternalDofManager( i );
            if ( !intDofMan ) continue; // you may never know...

            element->giveInternalDofManDofIDMask( i, ids );
            intDofMan->giveLocationArray( ids, nodalEqArray, EModelDefaultEquationNumbering() );
            eqArray.followedBy( nodalEqArray );
            intDofMan->giveMasterDofIDArray( ids, masterDofIDs );

            int tempN          = intDofMan->giveNumber();
            const auto &coords = element->giveDofManager( tempN )->giveCoordinates();
            // for (auto eqN : nodalEqArray) {
            //     if (eqN == 0) {
            //         continue;
            //     }
            //     tempMat2.addSubVectorRow(coords, eqN, 1);
            // }
        }

        const int nBaseDofs = element->giveNumberOfDofManagers() * 6;
        for ( int iDof = 0; iDof < 3; ++iDof ) {
            const auto dType = dofids[iDof];
            FloatArray localVector( nBaseDofs );
            FloatArray globalVector;
            FloatArray dirVector;
            dirVector.beColumnOf( lcs, iDof % 3 + 1 );
            // get unit vectors in local dofs
            for ( int inode = 0; inode < element->giveNumberOfDofManagers(); ++inode ) {
                localVector.addSubVector( dirVector, ( iDof / 3 ) * 3 + inode * 6 + 1 );
            }
            // work it back to global dofs
            globalVector.beTProductOf( GtoL, localVector );

            // apply the components to the condensed dofs on the global matrix
            for ( int i = 1; i <= eqArray.giveSize(); ++i ) {
                const auto eqN            = eqArray.at( i );
                const auto val            = globalVector( i + nBaseDofs - 1 );
                unitDisp.at( eqN, dType ) = val;
            }
        }
    } // end of translational unit vectors among internal dof managers

    for ( int i = 1; i <= 3; ++i ) {
        tempCol.beColumnOf( unitDisp, i );
        massMatrix->times( tempCol, tempCol2 ); // now tempCol2 has only the masses pertaining the i-th direction
        tempMat.setColumn( tempCol2, i );
        totMass.at( i ) = tempCol.dotProduct( tempCol2 );
    }

    // we have the centroid. we can now calculate rotational components. first from nodes.
    for ( std::unique_ptr<DofManager> &node : domain->giveDofManagers() ) {
        if ( !node->giveNumberOfDofs() ) continue;

        FloatArray vk( 3 );
        IntArray eq( 3 );

        // TODO consider own UCS if present

        IntArray mstrDofs, locArr;
        node->givePrimaryDofs( mstrDofs );
        node->giveLocationArray( mstrDofs, locArr, EModelDefaultEquationNumbering() );

        int partialDofCount = locArr.giveSize();
        if ( partialDofCount ) {
            // search for our dofs in there
            for ( int myDofIndex = 1; myDofIndex <= partialDofCount; myDofIndex++ ) {
                int dType = mstrDofs.at( myDofIndex );
                int eqN   = locArr.at( myDofIndex );

                if ( ( dType >= D_u ) && ( dType <= D_w ) && eqN ) {
                    // save unit displacement and coordinate
                    vk.at( dType ) = coordFilter.at( dType ) * ( node->giveCoordinate( dType ) - centroid.at( dType ) );
                    eq.at( dType ) = eqN;
                }
            }
        }

        // set mixed contribution due to rotation about centroid
        if ( eq.at( 1 ) ) {
            unitDisp.at( eq.at( 1 ), R_v ) = vk.at( 3 );
            unitDisp.at( eq.at( 1 ), R_w ) = -vk.at( 2 );
        }

        if ( eq.at( 2 ) ) {
            unitDisp.at( eq.at( 2 ), R_u ) = -vk.at( 3 );
            unitDisp.at( eq.at( 2 ), R_w ) = vk.at( 1 );
        }

        if ( eq.at( 3 ) ) {
            unitDisp.at( eq.at( 3 ), R_u ) = vk.at( 2 );
            unitDisp.at( eq.at( 3 ), R_v ) = -vk.at( 1 );
        }

        if ( partialDofCount ) {
            // search for our dofs in there
            for ( int myDofIndex = 1; myDofIndex <= partialDofCount; myDofIndex++ ) {
                int dType = mstrDofs.at( myDofIndex );
                int eqN   = locArr.at( myDofIndex );

                if ( ( dType >= R_u ) && ( dType <= R_w ) && eqN ) {
                    // save unit displacement and coordinate
                    unitDisp.at( eqN, dType ) = 1.0;
                }
            }
        }
    } // end of search among nodes

    // then from internaldof managers
    for ( int ielem = 1; ielem <= nelem; ++ielem ) {
        Element *element = domain->giveElement( ielem );
        if ( !element->giveNumberOfInternalDofManagers() ) {
            continue;
        }

        // we only support end releases on 3d beams for the moment. other internal dof managers should be skipped.
        if ( strcmp( element->giveClassName(), "Beam3d" ) != 0 ) {
            continue;
        }

        eqArray.clear();

        // use the lcs transform to work out which ghost dofs are affected by global displacement
        FloatMatrix GtoL;
        element->computeGtoLRotationMatrix( GtoL );
        FloatMatrix lcs;
        element->giveLocalCoordinateSystem( lcs );

        FloatMatrix componentMatrix( 3, 3 );
        //	retrieve internal dof managers and location array
        for ( int i = 1; i <= element->giveNumberOfInternalDofManagers(); ++i ) {
            DofManager *intDofMan = element->giveInternalDofManager( i );
            if ( !intDofMan ) continue; // you may never know...

            element->giveInternalDofManDofIDMask( i, ids );
            intDofMan->giveLocationArray( ids, nodalEqArray, EModelDefaultEquationNumbering() );
            eqArray.followedBy( nodalEqArray );
            intDofMan->giveMasterDofIDArray( ids, masterDofIDs );

            int tempN          = intDofMan->giveNumber();
            const auto &coords = element->giveDofManager( tempN )->giveCoordinates();
            componentMatrix.setColumn( coords - centroid, 2 );

            for ( int iDof = 1; iDof <= masterDofIDs.giveSize(); ++iDof ) {
                const int myDof = masterDofIDs.at( iDof );
                // now we're only interested about the effect of rotation on translational dofs
                if ( myDof > D_w ) {
                    continue;
                }
                const int eqN = nodalEqArray.at( iDof );
                for ( int j = 1; j <= 3; ++j ) {
                    componentMatrix.at( j, 3 ) = lcs.at( myDof, j );
                }
                // project the displacement component due to rotation onto the local direction of the dof
                for ( int rotDof = R_u; rotDof <= R_w; ++rotDof ) {
                    FloatArray rotationVector( 3 );
                    rotationVector.at( rotDof - 3 ) = 1.0;
                    componentMatrix.setColumn( rotationVector, 1 );
                    const double displacement  = componentMatrix.giveDeterminant(); // triple product: rotation ^ (node-centroid) . dofVector
                    unitDisp.at( eqN, rotDof ) = displacement;
                }
            }
        } // end of search among internal dof managers

        const int nBaseDofs = element->giveNumberOfDofManagers() * 6;
        for ( int iDof = 0; iDof < 6; ++iDof ) {
            const auto dType = dofids[iDof];
            if ( dType < R_u || dType > R_w ) {
                continue;
            }

            FloatArray localVector( nBaseDofs );
            FloatArray globalVector;
            FloatArray dirVector;
            dirVector.beColumnOf( lcs, iDof % 3 + 1 );
            // get unit vectors in local dofs
            for ( int inode = 0; inode < element->giveNumberOfDofManagers(); ++inode ) {
                localVector.addSubVector( dirVector, ( iDof / 3 ) * 3 + inode * 6 + 1 );
            }
            // work it back to global dofs
            globalVector.beTProductOf( GtoL, localVector );

            // apply the components to the condensed dofs on the global matrix
            for ( int i = 1; i <= eqArray.giveSize(); ++i ) {
                const auto eqN            = eqArray.at( i );
                const auto val            = globalVector( i + nBaseDofs - 1 );
                unitDisp.at( eqN, dType ) = val;
            }
        }
    }
    // end of creation of translational unit displacement vectors

    for ( int i = 4; i <= 6; ++i ) {
        tempCol.beColumnOf( unitDisp, i );
        massMatrix->times( tempCol, tempCol2 ); // now tempCol2 has only the masses pertaining the i-th direction
        tempMat.setColumn( tempCol2, i );
        totMass.at( i ) = tempCol.dotProduct( tempCol2 ); // total mass for i-th direction
    }

    std::unique_ptr<EigenVectorPrimaryField> unitDispField = std::make_unique<EigenVectorPrimaryField>( this, 1, FT_Displacements, 6 ); // 6 dofs
    unitDispField->updateAll( unitDisp, defNumbering );

    //
    // calculate participation factors and mass participation
    //

    // participation factors
    partFact.beTProductOf( eigVec, tempMat );

    // mass participation ratios
    for ( int i = 1; i <= numberOfRequiredEigenValues; ++i ) {
        for ( int j = 1; j <= 6; ++j ) {
            if ( totMass.at( j ) > 1e-10 ) {
                massPart.at( i, j ) = pow( partFact.at( i, j ), 2 ) / totMass.at( j );
            }
            // else
            //{
            //	massPart.at(i, j) = 0.0;
            // }
        }
    }

    //
    // start creating loaded models
    //
#ifdef VERBOSE
    OOFEM_LOG_INFO( "Starting analysis for each mode ...\n" );
#endif


    // create linear solver
    std::unique_ptr<SparseLinearSystemNM> nLinMethod; // = NULL;

    if ( !nLinMethod ) {
        if ( isParallel() ) {
            if ( ( solverType == ST_Petsc ) || ( solverType == ST_Feti ) ) {
                nLinMethod = classFactory.createSparseLinSolver( ST_Direct, this->giveDomain( 1 ), this );
            }
        } else {
            nLinMethod = classFactory.createSparseLinSolver( linStype, this->giveDomain( 1 ), this );
        }
        if ( !nLinMethod ) {
            OOFEM_ERROR( "linear solver creation failed for lstype %d", solverType );
        }
    }

    // determine dominat mode in requested direction
    FloatArray dirVect( dir );
    FloatArray dirFactors( numberOfRequiredEigenValues );
    dirVect.normalize();
    for ( int dN = 1; dN <= numberOfRequiredEigenValues; ++dN ) {
        for ( int nDir = 1; nDir <= 3; ++nDir )
            dirFactors.at( dN ) += partFact.at( dN, nDir ) * dirVect.at( nDir );
        dirFactors.at( dN ) = fabs( dirFactors.at( dN ) );
    }
    // mode with highes participation factor in requested direction
    dominantMode = std::distance( dirFactors.begin(), std::max_element( dirFactors.begin(), dirFactors.end() ) ) + 1;

    for ( int dN = 1; dN <= numberOfRequiredEigenValues; ++dN ) {
        OOFEM_LOG_INFO( "Analyzing mode %d...\n", dN );

        double sAcc = calcSpectrumOrdinate( periods.at( dN ) );

        loadVector.clear();

        // TODO: recreate the load vector from the unitdisp field
        // iterate on dofmans, create load given displacement field
        for ( int nDir = 1; nDir <= 3; ++nDir ) {
            if ( dir.at( nDir ) != 0 ) {
                tempCol.beColumnOf( eigVec, dN );
                massMatrix->times( tempCol, tempCol2 );
                tempCol2 *= ( sAcc * partFact.at( dN, nDir ) * dir.at( nDir ) ); // scaled forces
                loadVector.add( tempCol2 );
            }
        }

        // solve linear system
        NM_Status s = nLinMethod->solve( *stiffnessMatrix, loadVector, dummyDisps ); // solve linear system
        if ( !( s & NM_Success ) ) {
            OOFEM_ERROR( "No success in solving system." );
        }

        // update elements to we can get internal state!!!
        this->updateInternalState( tStep );
        EngngModel::updateYourself( tStep );

        FloatArray reactions;
        this->buildReactionTable( dofManMap, dofidMap, eqnMap, tStep, 1 );

        // compute reaction forces
        this->computeReaction( reactions, tStep, 1 );

        reactionsList.push_back( reactions );
        dispList.push_back( dummyDisps );

        map<int, map<int, map<int, map<string, FloatArray> > > > elemResponse;
        map<int, map<string, FloatArray> > beamResponse;

        for ( auto &elem : domain->giveElements() ) {
            // test for remote element in parallel mode
            if ( elem->giveParallelMode() == Element_remote ) {
                continue;
            }

            // for (int i = 1; i <= elem->giveNumberOfInternalDofManagers(); ++i) {
            //	DofManager *dman = elem->giveInternalDofManager(i);
            //	dman->printOutputAt(outputStream, tStep);
            // }

            map<int, map<int, map<string, FloatArray> > > *eir = new map<int, map<int, map<string, FloatArray> > >;

            for ( int i = 0; i < elem->giveNumberOfIntegrationRules(); ++i ) {
                map<int, map<string, FloatArray> > *ir = NULL;
                this->getIntRuleOutputAt( elem->giveIntegrationRule( i ), tStep, ir );
                if ( ir ) eir->operator[]( i ) = *ir;
            }

            elemResponse[elem->giveNumber()] = *eir;

            const char *name = elem->giveClassName();

            if ( ( strcmp( name, "Beam3d" ) == 0 ) || ( strcmp( name, "Beam2d" ) == 0 ) || ( strcmp( name, "beam3d" ) == 0 ) || ( strcmp( name, "beam2d" ) == 0 ) ) {
                map<string, FloatArray> *b = new map<string, FloatArray>;

                FloatArray rl, Fl;
                elem->computeVectorOf( VM_Total, tStep, rl );
                // ask for global element end forces vector
                StructuralElement *SElem;

                b->operator[]( "enddisp" ) = rl;

                SElem = static_cast<StructuralElement *>( elem.get() );
                SElem->giveInternalForcesVector( Fl, tStep );

                // FloatArray loadEndForces;

                //// add exact end forces due to nonnodal loading
                // SElem->computeForceLoadVector(loadEndForces, tStep, VM_Total);
                // if (loadEndForces.giveSize()) {
                //	Fl.subtract(loadEndForces);
                // }

                b->operator[]( "endforces" ) = Fl;

                // #ifdef DEBUG
                //					if (SElem->giveNumber() == 7) {
                //							FloatArray tempnd;
                //							DofManager *dofMan = SElem->giveDofManager(2);
                //							dofMan->giveCompleteUnknownVector(tempnd, VM_Total, tStep);
                //							OOFEM_WARNING("Stop %e",tempnd.at(1));
                //							}
                // #endif
                beamResponse[elem->giveNumber()] = *b;
            }
        }

        beamResponseList.push_back( beamResponse );
        elemResponseList.push_back( elemResponse );

        // bem to call addSquared
        exportModuleManager.doOutput( tStep );
    }

    if ( modalCombo == RSpecComboType::RSC_SRSS ) {
        this->SRSS();
    } else {
        this->CQC();
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

void ResponseSpectrum::SRSS()
{
    list<map<int, map<int, map<int, map<string, FloatArray> > > > >::iterator elem_it = elemResponseList.begin();
    list<map<int, map<string, FloatArray> > >::iterator beam_it                       = beamResponseList.begin();
    for ( ; elem_it != elemResponseList.end(); ++elem_it ) {
        addMultiply( combElemResponse, *elem_it, *elem_it );
    }

    for ( ; beam_it != beamResponseList.end(); ++beam_it ) {
        addMultiply( combBeamResponse, *beam_it, *beam_it );
    }

    list<FloatArray>::iterator reac_it = reactionsList.begin();
    list<FloatArray>::iterator disp_it = dispList.begin();

    for ( ; reac_it != reactionsList.end(); ++reac_it ) {
        FloatArray &reactions = *reac_it;
        for ( int z = 1; z <= reactions.giveSize(); ++z ) {
            combReactions.at( z ) += pow( reactions.at( z ), 2 );
        }
    }

    for ( ; disp_it != dispList.end(); ++disp_it ) {
        FloatArray &disps = *disp_it;
        for ( int z = 1; z <= disps.giveSize(); ++z ) {
            combDisps.at( z ) += pow( disps.at( z ), 2 );
        }
    }

    for ( int z = 1; z <= combReactions.giveSize(); ++z ) {
        combReactions.at( z ) = sqrt( combReactions.at( z ) );
    }

    disp_it           = std::next( dispList.begin(), dominantMode - 1 );
    FloatArray &disps = *disp_it; // dominant mode
    for ( int z = 1; z <= combDisps.giveSize(); ++z ) {
        double res = sqrt( combDisps.at( z ) );
        res *= signbit( disps.at( z ) ) ? -1 : 1;
        combDisps.at( z ) = res;
    }

    calcRoot( combElemResponse );
    calcRoot( combBeamResponse );
}


void ResponseSpectrum::CQC()
{
    list<map<int, map<int, map<int, map<string, FloatArray> > > > >::iterator elem_it = elemResponseList.begin();

    for ( int i = 1; elem_it != elemResponseList.end(); ++elem_it, ++i ) {
        list<map<int, map<int, map<int, map<string, FloatArray> > > > >::iterator elem_it2 = elemResponseList.begin();
        for ( int j = 1; elem_it2 != elemResponseList.end(); ++elem_it2, ++j ) {
            addMultiply( combElemResponse, *elem_it, *elem_it2, rhos.at( i, j ) );
        }
    }
    calcRoot( combElemResponse );

    list<map<int, map<string, FloatArray> > >::iterator beam_it = beamResponseList.begin();
    for ( int i = 1; beam_it != beamResponseList.end(); ++beam_it, ++i ) {
        list<map<int, map<string, FloatArray> > >::iterator beam_it2 = beamResponseList.begin();
        for ( int j = 1; beam_it2 != beamResponseList.end(); ++beam_it2, ++j ) {
            addMultiply( combBeamResponse, *beam_it, *beam_it2, rhos.at( i, j ) );
        }
    }
    calcRoot( combBeamResponse );

    list<FloatArray>::iterator reac_it = reactionsList.begin();
    for ( int i = 1; reac_it != reactionsList.end(); ++reac_it, ++i ) {
        FloatArray &reactions               = *reac_it;
        list<FloatArray>::iterator reac_it2 = reactionsList.begin();
        for ( int j = 1; reac_it2 != reactionsList.end(); ++reac_it2, ++j ) {
            FloatArray &reactions2 = *reac_it2;
            for ( int z = 1; z <= reactions.giveSize(); ++z ) {
                combReactions.at( z ) += fabs( reactions.at( z ) * reactions2.at( z ) * rhos.at( i, j ) );
            }
        }
    }
    for ( int z = 1; z <= combReactions.giveSize(); ++z ) {
        combReactions.at( z ) = sqrt( combReactions.at( z ) );
    }

    list<FloatArray>::iterator disp_it = dispList.begin();
    for ( int i = 1; disp_it != dispList.end(); ++disp_it, ++i ) {
        FloatArray &disps                   = *disp_it;
        list<FloatArray>::iterator disp_it2 = dispList.begin();
        for ( int j = 1; disp_it2 != dispList.end(); ++disp_it2, ++j ) {
            FloatArray &disps2 = *disp_it2;
            for ( int z = 1; z <= disps.giveSize(); ++z ) {
                combDisps.at( z ) += fabs( disps.at( z ) * disps2.at( z ) * rhos.at( i, j ) );
            }
        }
    }

    disp_it           = std::next( dispList.begin(), dominantMode - 1 );
    FloatArray &disps = *disp_it; // dominant mode
    for ( int z = 1; z <= combDisps.giveSize(); ++z ) {
        double res = sqrt( combDisps.at( z ) );
        res *= signbit( disps.at( z ) ) ? -1 : 1;
        combDisps.at( z ) = res;
    }
}


void ResponseSpectrum::getGPOutputAt( GaussPoint *gp, TimeStep *tStep, std::map<std::string, FloatArray> *&ips )
{
    // int iruleNumber = 0;

    // if (gp->irule) {
    //	iruleNumber = irule->giveNumber();
    // }

    // fprintf(File, "  GP %2d.%-2d :", iruleNumber, number);

    // invoke printOutputAt method for all managed statuses
    IntegrationPointStatus *status = gp->giveMaterialStatus();
    if ( status ) {
        this->getIntPointStatusOutputAt( status, tStep, gp->giveMaterialMode(), ips );
    }

    // if (gp->gaussPoints.size() != 0) { // layered material
    //	fprintf(File, "Layers report \n{\n");
    //	for (GaussPoint *gp : gaussPoints) {
    //		gp->printOutputAt(File, tStep);
    //	}

    //	fprintf(File, "} end layers report\n");
    //}
}

void ResponseSpectrum::giveRhos( FloatMatrix &ans ) { ans = rhos; }

void ResponseSpectrum::giveDominantMode( int &mode ) { mode = dominantMode; }

RSpecComboType ResponseSpectrum::giveComboType() { return modalCombo; }

void ResponseSpectrum::getIntRuleOutputAt( IntegrationRule *iRule, TimeStep *tStep, map<int, map<string, FloatArray> > *&ir )
{
    map<string, FloatArray> *igp = NULL;
    ir                           = new map<int, map<string, FloatArray> >;
    for ( GaussPoint *gp : *iRule ) {
        this->getGPOutputAt( gp, tStep, igp );
        if ( igp ) ir->operator[]( gp->giveNumber() ) = *igp;
    }
}

void ResponseSpectrum::getIntPointStatusOutputAt( IntegrationPointStatus *iStatus, TimeStep *tStep, MaterialMode materialMode, map<string, FloatArray> *&ir )
{
    ir                              = NULL;
    StructuralMaterialStatus *strMS = dynamic_cast<StructuralMaterialStatus *>( iStatus );
    if ( strMS ) {
        FloatArray helpVec;
        ir = new map<string, FloatArray>;

        StructuralMaterial::giveFullSymVectorForm( helpVec, strMS->giveStrainVector(), materialMode );
        // for (auto &var : helpVec) {
        //	fprintf(File, " %.4e", var);
        // }
        ir->operator[]( "strains" ) = helpVec;

        // fprintf(File, "\n              stresses");
        StructuralMaterial::giveFullSymVectorForm( helpVec, strMS->giveStressVector(), materialMode );
        ir->operator[]( "stresses" ) = helpVec;

        // for (auto &var : helpVec) {
        //	fprintf(File, " %.4e", var);
        // }
        // fprintf(File, "\n");
    }
}

void populateElResults( map<int, map<int, map<int, map<string, FloatArray> > > > &answer, map<int, map<int, map<int, map<string, FloatArray> > > > &src )
{
    map<int, map<int, map<int, map<string, FloatArray> > > >::iterator srcElem_it = src.begin();
    for ( ; srcElem_it != src.end(); ++srcElem_it ) {
        map<int, map<int, map<string, FloatArray> > > *destElIntRuleMap = new map<int, map<int, map<string, FloatArray> > >;
        map<int, map<int, map<string, FloatArray> > > &srcElIntRuleMap  = srcElem_it->second;

        map<int, map<int, map<string, FloatArray> > >::iterator srcElIntRuleMap_it = srcElIntRuleMap.begin();
        for ( ; srcElIntRuleMap_it != srcElIntRuleMap.end(); ++srcElIntRuleMap_it ) {
            map<int, map<string, FloatArray> > *destGPMap = new map<int, map<string, FloatArray> >;
            map<int, map<string, FloatArray> > &srcGPMap  = srcElIntRuleMap_it->second;

            map<int, map<string, FloatArray> >::iterator srcGPMap_it = srcGPMap.begin();
            for ( ; srcGPMap_it != srcGPMap.end(); ++srcGPMap_it ) {
                map<string, FloatArray> *destRespMap = new map<string, FloatArray>;
                map<string, FloatArray> &srcRespMap  = srcGPMap_it->second;

                map<string, FloatArray>::iterator srcRespMap_it = srcRespMap.begin();
                for ( ; srcRespMap_it != srcRespMap.end(); ++srcRespMap_it ) {
                    FloatArray &srcRespArray  = srcRespMap_it->second;
                    FloatArray *destRespArray = new FloatArray( srcRespArray.giveSize() );

                    destRespMap->operator[]( srcRespMap_it->first ) = *destRespArray;
                }

                destGPMap->operator[]( srcGPMap_it->first ) = *destRespMap;
            }

            destElIntRuleMap->operator[]( srcElIntRuleMap_it->first ) = *destGPMap;
        }

        answer[srcElem_it->first] = *destElIntRuleMap;
    }
}

void addMultiply( map<int, map<int, map<int, map<string, FloatArray> > > > &answer, map<int, map<int, map<int, map<string, FloatArray> > > > &src, map<int, map<int, map<int, map<string, FloatArray> > > > &src2, double fact )
{
    if ( answer.size() == 0 ) {
        populateElResults( answer, src );
    }

    map<int, map<int, map<int, map<string, FloatArray> > > >::iterator destElem_it = answer.begin();
    map<int, map<int, map<int, map<string, FloatArray> > > >::iterator srcElem_it  = src.begin();
    map<int, map<int, map<int, map<string, FloatArray> > > >::iterator srcElem_it2 = src2.begin();
    for ( ; destElem_it != answer.end(); ++destElem_it, ++srcElem_it, ++srcElem_it2 ) {
        map<int, map<int, map<string, FloatArray> > > &destElIntRuleMap = destElem_it->second;
        map<int, map<int, map<string, FloatArray> > > &srcElIntRuleMap  = srcElem_it->second;
        map<int, map<int, map<string, FloatArray> > > &srcElIntRuleMap2 = srcElem_it2->second;

        map<int, map<int, map<string, FloatArray> > >::iterator destElIntRuleMap_it = destElIntRuleMap.begin();
        map<int, map<int, map<string, FloatArray> > >::iterator srcElIntRuleMap_it  = srcElIntRuleMap.begin();
        map<int, map<int, map<string, FloatArray> > >::iterator srcElIntRuleMap_it2 = srcElIntRuleMap2.begin();
        for ( ; destElIntRuleMap_it != destElIntRuleMap.end(); ++destElIntRuleMap_it, ++srcElIntRuleMap_it, ++srcElIntRuleMap_it2 ) {
            map<int, map<string, FloatArray> > &destGPMap = destElIntRuleMap_it->second;
            map<int, map<string, FloatArray> > &srcGPMap  = srcElIntRuleMap_it->second;
            map<int, map<string, FloatArray> > &srcGPMap2 = srcElIntRuleMap_it2->second;

            map<int, map<string, FloatArray> >::iterator destGPMap_it = destGPMap.begin();
            map<int, map<string, FloatArray> >::iterator srcGPMap_it  = srcGPMap.begin();
            map<int, map<string, FloatArray> >::iterator srcGPMap_it2 = srcGPMap2.begin();
            for ( ; destGPMap_it != destGPMap.end(); ++destGPMap_it, ++srcGPMap_it, ++srcGPMap_it2 ) {
                map<string, FloatArray> &destRespMap = destGPMap_it->second;
                map<string, FloatArray> &srcRespMap  = srcGPMap_it->second;
                map<string, FloatArray> &srcRespMap2 = srcGPMap_it2->second;

                map<string, FloatArray>::iterator destRespMap_it = destRespMap.begin();
                map<string, FloatArray>::iterator srcRespMap_it  = srcRespMap.begin();
                map<string, FloatArray>::iterator srcRespMap_it2 = srcRespMap2.begin();
                for ( ; destRespMap_it != destRespMap.end(); ++destRespMap_it, ++srcRespMap_it, ++srcRespMap_it2 ) {
                    FloatArray &destRespArray = destRespMap_it->second;
                    FloatArray &srcRespArray  = srcRespMap_it->second;
                    FloatArray &srcRespArray2 = srcRespMap_it2->second;

                    for ( int i = 1; i <= srcRespArray.giveSize(); ++i ) {
                        // square it and add it
                        destRespArray.at( i ) += fabs( srcRespArray.at( i ) * srcRespArray2.at( i ) * fact );
                    }
                }
            }
        }
    }
}

void calcRoot( map<int, map<int, map<int, map<string, FloatArray> > > > &answer )
{
    map<int, map<int, map<int, map<string, FloatArray> > > >::iterator destElem_it = answer.begin();
    for ( ; destElem_it != answer.end(); ++destElem_it ) {
        map<int, map<int, map<string, FloatArray> > > &destElIntRuleMap = destElem_it->second;

        map<int, map<int, map<string, FloatArray> > >::iterator destElIntRuleMap_it = destElIntRuleMap.begin();
        for ( ; destElIntRuleMap_it != destElIntRuleMap.end(); ++destElIntRuleMap_it ) {
            map<int, map<string, FloatArray> > &destGPMap = destElIntRuleMap_it->second;

            map<int, map<string, FloatArray> >::iterator destGPMap_it = destGPMap.begin();
            for ( ; destGPMap_it != destGPMap.end(); ++destGPMap_it ) {
                map<string, FloatArray> &destRespMap = destGPMap_it->second;

                map<string, FloatArray>::iterator destRespMap_it = destRespMap.begin();
                for ( ; destRespMap_it != destRespMap.end(); ++destRespMap_it ) {
                    FloatArray &destRespArray = destRespMap_it->second;

                    for ( int i = 1; i <= destRespArray.giveSize(); ++i ) {
                        // square it and add it
                        destRespArray.at( i ) = sqrt( destRespArray.at( i ) );
                    }
                }
            }
        }
    }
}

void populateElResults( map<int, map<string, FloatArray> > &answer, map<int, map<string, FloatArray> > &src )
{
    map<int, map<string, FloatArray> >::iterator srcElem_it = src.begin();
    for ( ; srcElem_it != src.end(); ++srcElem_it ) {
        map<string, FloatArray> *destBRespMap = new map<string, FloatArray>;
        map<string, FloatArray> &srcBRespMap  = srcElem_it->second;

        map<string, FloatArray>::iterator srcBRespMap_it = srcBRespMap.begin();
        for ( ; srcBRespMap_it != srcBRespMap.end(); ++srcBRespMap_it ) {
            FloatArray &srcRespArray  = srcBRespMap_it->second;
            FloatArray *destRespArray = new FloatArray( srcRespArray.giveSize() );

            destBRespMap->operator[]( srcBRespMap_it->first ) = *destRespArray;
        }

        answer[srcElem_it->first] = *destBRespMap;
    }
}

void addMultiply( map<int, map<string, FloatArray> > &answer, map<int, map<string, FloatArray> > &src, map<int, map<string, FloatArray> > &src2, double fact )
{
    if ( answer.size() == 0 ) {
        populateElResults( answer, src );
    }

    map<int, map<string, FloatArray> >::iterator destElem_it = answer.begin();
    map<int, map<string, FloatArray> >::iterator srcElem_it  = src.begin();
    map<int, map<string, FloatArray> >::iterator srcElem_it2 = src2.begin();
    for ( ; destElem_it != answer.end(); ++destElem_it, ++srcElem_it, ++srcElem_it2 ) {
        map<string, FloatArray> &destRespMap = destElem_it->second;
        map<string, FloatArray> &srcRespMap  = srcElem_it->second;
        map<string, FloatArray> &srcRespMap2 = srcElem_it2->second;

        map<string, FloatArray>::iterator destRespMap_it = destRespMap.begin();
        map<string, FloatArray>::iterator srcRespMap_it  = srcRespMap.begin();
        map<string, FloatArray>::iterator srcRespMap_it2 = srcRespMap2.begin();
        for ( ; destRespMap_it != destRespMap.end(); ++destRespMap_it, ++srcRespMap_it, ++srcRespMap_it2 ) {
            FloatArray &destRespArray = destRespMap_it->second;
            FloatArray &srcRespArray  = srcRespMap_it->second;
            FloatArray &srcRespArray2 = srcRespMap_it2->second;

            for ( int i = 1; i <= srcRespArray.giveSize(); ++i ) {
                // square it and add it
                destRespArray.at( i ) += fabs( srcRespArray.at( i ) * srcRespArray2.at( i ) * fact );
            }
        }
    }
}

void calcRoot( map<int, map<string, FloatArray> > &answer )
{
    map<int, map<string, FloatArray> >::iterator destElem_it = answer.begin();
    for ( ; destElem_it != answer.end(); ++destElem_it ) {
        map<string, FloatArray> &destRespMap = destElem_it->second;

        map<string, FloatArray>::iterator destRespMap_it = destRespMap.begin();
        for ( ; destRespMap_it != destRespMap.end(); ++destRespMap_it ) {
            FloatArray &destRespArray = destRespMap_it->second;

            for ( int i = 1; i <= destRespArray.giveSize(); ++i ) {
                // square it and add it
                destRespArray.at( i ) = sqrt( destRespArray.at( i ) );
            }
        }
    }
}

void ResponseSpectrum::computeExternalLoadReactionContribution( FloatArray &reactions, TimeStep *tStep, int di )
{
    reactions.resize( this->giveNumberOfDomainEquations( di, EModelDefaultPrescribedEquationNumbering() ) );
    reactions.zero();
    this->assembleVector( reactions, tStep, ExternalForceAssembler(), VM_Total,
        EModelDefaultPrescribedEquationNumbering(), this->giveDomain( di ) );
}


void ResponseSpectrum::buildReactionTable( IntArray &restrDofMans, IntArray &restrDofs,
    IntArray &eqn, TimeStep *tStep, int di )
{
    // determine number of restrained dofs
    Domain *domain   = this->giveDomain( di );
    int numRestrDofs = this->giveNumberOfDomainEquations( di, EModelDefaultPrescribedEquationNumbering() );
    int ndofMan      = domain->giveNumberOfDofManagers();
    int rindex, count = 0;

    // initialize corresponding dofManagers and dofs for each restrained dof
    restrDofMans.resize( numRestrDofs );
    restrDofs.resize( numRestrDofs );
    eqn.resize( numRestrDofs );

    for ( int i = 1; i <= ndofMan; ++i ) {
        DofManager *inode = domain->giveDofManager( i );
        for ( Dof *jdof : *inode ) {
            if ( jdof->isPrimaryDof() && ( jdof->hasBc( tStep ) ) ) {
                // skip slave dofs
                rindex = jdof->__givePrescribedEquationNumber();
                if ( rindex ) {
                    count++;
                    restrDofMans.at( count ) = i;
                    restrDofs.at( count )    = jdof->giveDofID();
                    eqn.at( count )          = rindex;
                } else {
                    // NullDof has no equation number and no prescribed equation number
                    //_error("No prescribed equation number assigned to supported DOF");
                }
            }
        }
    }
    // Trim to size.
    restrDofMans.resizeWithValues( count );
    restrDofs.resizeWithValues( count );
    eqn.resizeWithValues( count );
}

void ResponseSpectrum::computeReaction( FloatArray &answer, TimeStep *tStep, int di )
{
    FloatArray contribution;

    answer.resize( this->giveNumberOfDomainEquations( di, EModelDefaultPrescribedEquationNumbering() ) );
    answer.zero();

    // Add internal forces
    this->assembleVector( answer, tStep, LastEquilibratedInternalForceAssembler(), VM_Total,
        EModelDefaultPrescribedEquationNumbering(), this->giveDomain( di ) );
    // Subtract external loading
    ///@todo All engineering models should be using this (for consistency)
    // this->assembleVector( answer, tStep, ExternalForceAssembler(), VM_Total,
    //                     EModelDefaultPrescribedEquationNumbering(), this->giveDomain(di) );
    ///@todo This method is overloaded in some functions, it needs to be generalized.
    this->computeExternalLoadReactionContribution( contribution, tStep, di );
    answer.subtract( contribution );
    this->updateSharedDofManagers( answer, EModelDefaultPrescribedEquationNumbering(), ReactionExchangeTag );
}


void ResponseSpectrum::updateInternalState( TimeStep *tStep )
{
    for ( auto &domain : domainList ) {
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( auto &dman : domain->giveDofManagers() ) {
                this->updateDofUnknownsDictionary( dman.get(), tStep );
            }
        }

        for ( auto &bc : domain->giveBcs() ) {
            ActiveBoundaryCondition *abc;

            if ( ( abc = dynamic_cast<ActiveBoundaryCondition *>( bc.get() ) ) ) {
                int ndman = abc->giveNumberOfInternalDofManagers();
                for ( int j = 1; j <= ndman; ++j ) {
                    this->updateDofUnknownsDictionary( abc->giveInternalDofManager( j ), tStep );
                }
            }
        }

        if ( true ) {
            // internalVarUpdateStamp != tStep->giveSolutionStateCounter()
            for ( auto &elem : domain->giveElements() ) {
                elem->updateInternalState( tStep );
            }

            // internalVarUpdateStamp = tStep->giveSolutionStateCounter();
        }
    }
}


void ResponseSpectrum::updateYourself( TimeStep *tStep )
{
    // this->updateInternalState(tStep);
    // EngngModel::updateYourself(tStep);
}


void ResponseSpectrum::terminate( TimeStep *tStep )
{
    Domain *domain     = this->giveDomain( 1 );
    FILE *outputStream = this->giveOutputStream();

    // print loadcase header
    fprintf( outputStream, "\nOutput for Response Spectrum analysis \n\n" );
    // print eigen values on output
    fprintf( outputStream, "\n\nEigen Values (Omega^2) are:\n-----------------\n" );

    for ( int i = 1; i <= numberOfRequiredEigenValues; ++i ) {
        fprintf( outputStream, "%15.8e ", eigVal.at( i ) );
        if ( ( i % 5 ) == 0 ) {
            fprintf( outputStream, "\n" );
        }
    }


    fprintf( outputStream, "\n\n\nCentroid Coordinates are:\n-----------------\n\tX\t|\tY\t|\tZ\n" );
    for ( int i = 1; i <= centroid.giveSize(); ++i ) {
        fprintf( outputStream, "%10.3e ", centroid.at( i ) );
    }

    fprintf( outputStream, "\n" );

    fprintf( outputStream, "\n\nParticipation Factors are:\n-----------------\n\tDx\t|\tDy\t|\tDz\t|\tRx\t|\tRy\t|\tRz\n" );
    for ( int i = 1; i <= partFact.giveNumberOfRows(); ++i ) {
        for ( int j = 1; j <= partFact.giveNumberOfColumns(); ++j ) {
            fprintf( outputStream, "%10.3e ", partFact.at( i, j ) );
        }
        fprintf( outputStream, "\n" );
    }

    fprintf( outputStream, "\n\nTotal Masses are:\n-----------------\n\tDx\t|\tDy\t|\tDz\t|\tRx\t|\tRy\t|\tRz\n" );
    for ( int i = 1; i <= totMass.giveSize(); ++i ) {
        fprintf( outputStream, "%10.3e ", totMass.at( i ) );
    }

    fprintf( outputStream, "\n" );

    fprintf( outputStream, "\n\nMass Ratios are:\n-----------------\n\tDx\t|\tDy\t|\tDz\t|\tRx\t|\tRy\t|\tRz\n" );
    for ( int i = 1; i <= massPart.giveNumberOfRows(); ++i ) {
        for ( int j = 1; j <= massPart.giveNumberOfColumns(); ++j ) {
            fprintf( outputStream, "%10.3e ", massPart.at( i, j ) );
        }
        fprintf( outputStream, "\n" );
    }

    fprintf( outputStream, "\n\n" );

    fprintf( outputStream, "\n\nDofManager output:\n------------------\n" );

    // change tStep to 0 to allow final displacements extraction
    tStep->setIntrinsicTime( 0.0 );

    for ( auto &dman : domain->giveDofManagers() ) {
        dman->updateYourself( tStep );
        fprintf( outputStream, "%-8s%8d (%8d):\n", dman->giveClassName(), dman->giveLabel(), dman->giveNumber() );
        for ( Dof *dof : *dman ) {
            this->printDofOutputAt( outputStream, dof, tStep );
        }
    }

    // for ( int i = 1; i <=  numberOfRequiredEigenValues; ++i ) {
    //     fprintf(outputStream, "\nOutput for eigen value no.  %.3e \n", ( double ) i);
    //     fprintf( outputStream,
    //             "Printing eigen vector no. %d, corresponding eigen value is %15.8e\n\n",
    //             i, eigVal.at(i) );
    //     tStep->setTime( ( double ) i ); // we use time as intrinsic eigen value index

    //    if ( this->requiresUnknownsDictionaryUpdate() ) {
    //        for ( auto &dman : domain->giveDofManagers() ) {
    //            this->updateDofUnknownsDictionary(dman.get(), tStep);
    //        }
    //    }

#if 0
    fprintf( outputStream, "\n\nElement output:\n---------------\n" );

    for ( auto &elem : domain->giveElements() ) {
        // test for remote element in parallel mode
        if ( elem->giveParallelMode() == Element_remote ) { continue; }

        map<int, map<int, map<string, FloatArray>>> &destElIntRuleMap = combElemResponse[elem->giveNumber()];

        map<int, map<string, FloatArray>>::iterator it = combBeamResponse.find( elem->giveNumber() );
        //beam found;
        if ( it != combBeamResponse.end() ) {
            fprintf( outputStream, "beam element %d (%8d) :\n", elem->giveLabel(), elem->giveNumber() );
            fprintf( outputStream, "  local displacements " );
            for ( auto &val : ( it->second )["enddisp"] ) { fprintf( outputStream, " %.4e", val ); }

            fprintf( outputStream, "\n  local end forces    " );
            for ( auto &val : ( it->second )["endforces"] ) { fprintf( outputStream, " %.4e", val ); }
            fprintf( outputStream, "\n" );
        } else { fprintf( outputStream, "element %d (%8d) :\n", elem->giveLabel(), elem->giveNumber() ); }

        for ( int i = 1; i <= elem->giveNumberOfInternalDofManagers(); ++i ) {
            DofManager *dman = elem->giveInternalDofManager( i );
            //dman->printOutputAt(file, tStep);
        }

        for ( int i = 0; i < elem->giveNumberOfIntegrationRules(); ++i ) {
            IntegrationRule *iRule = elem->giveIntegrationRule( i );

            map<int, map<string, FloatArray>> &destGPMap = destElIntRuleMap[i];

            for ( GaussPoint *gp : *iRule ) {
                fprintf( outputStream, "  GP %2d.%-2d :", iRule->giveNumber(), gp->giveNumber() );
                IntegrationPointStatus *status = gp->giveMaterialStatus();
                if ( status ) {
                    StructuralMaterialStatus *strMS = dynamic_cast<StructuralMaterialStatus *>(status);
                    if ( strMS ) {
                        map<string, FloatArray> &destRespMap = destGPMap[gp->giveNumber()];

                        fprintf( outputStream, "  strains " );
                        //StructuralMaterial::giveFullSymVectorForm(helpVec, strMS->giveStrainVector(), gp->giveMaterialMode());
                        for ( auto &var : destRespMap["strains"] ) { fprintf( outputStream, " %.4e", var ); }

                        fprintf( outputStream, "\n              stresses" );
                        //StructuralMaterial::giveFullSymVectorForm(helpVec, strMS->giveStressVector(), gp->giveMaterialMode());

                        for ( auto &var : destRespMap["stresses"] ) { fprintf( outputStream, " %.4e", var ); }
                        fprintf( outputStream, "\n" );
                    }
                }
            }
        }
    }

#endif

    fprintf( outputStream, "\n\n\tR E A C T I O N S  O U T P U T:\n\t_______________________________\n\n\n" );

    for ( int i = 1; i <= dofManMap.giveSize(); ++i ) {
        if ( domain->giveOutputManager()->testDofManOutput( dofManMap.at( i ), tStep ) ) {
            fprintf( outputStream, "\tNode %8d iDof %2d reaction % .4e    [bc-id: %d]\n",
                domain->giveDofManager( dofManMap.at( i ) )->giveLabel(),
                dofidMap.at( i ), combReactions.at( eqnMap.at( i ) ),
                domain->giveDofManager( dofManMap.at( i ) )->giveDofWithID( dofidMap.at( i ) )->giveBcId() );
        }
    }

    // for ( int i = 1; i <=  numberOfRequiredEigenValues; ++i ) {
    //     // export using export manager
    //     tStep->setTime( ( double ) i ); // we use time as intrinsic eigen value index
    //     tStep->setNumber(i);
    exportModuleManager.doOutput( tStep ); // forcing bem with intrinsicTime=0 to compute the square root for SRSS and print results
    //}
    fprintf( outputStream, "strTerm\n" );
    fflush( this->giveOutputStream() );
    this->saveStepContext( tStep, CM_State | CM_Definition );
}


void ResponseSpectrum::saveContext( DataStream &stream, ContextMode mode )
//
// saves state variable - displacement vector
//
{
    contextIOResultType iores;
    FILE *file = NULL;

    EngngModel::saveContext( stream, mode );

    if ( ( iores = eigVal.storeYourself( stream ) ) != CIO_OK ) {
        THROW_CIOERR( iores );
    }

    if ( ( eigVec.storeYourself( stream ) ) != CIO_OK ) {
        THROW_CIOERR( iores );
    }
}


void ResponseSpectrum::restoreContext( DataStream &stream, ContextMode mode )
//
// restore state variable - displacement vector
//
{
    contextIOResultType iores;

    if ( restoreFlag == 0 ) {
        // save element context

        EngngModel::restoreContext( stream, mode );

        if ( ( iores = eigVal.restoreYourself( stream ) ) != CIO_OK ) {
            THROW_CIOERR( iores );
        }

        if ( ( iores = eigVec.restoreYourself( stream ) ) != CIO_OK ) {
            THROW_CIOERR( iores );
        }
    }

    if ( activeVector > numberOfRequiredEigenValues ) {
        activeVector = numberOfRequiredEigenValues;
    }

    this->restoreFlag = 1;
}

} // end namespace oofem
