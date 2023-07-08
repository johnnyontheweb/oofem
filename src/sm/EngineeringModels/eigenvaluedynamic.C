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
#include "domain.h"
#include "element.h"
#include "node.h"
#include "unknownnumberingscheme.h"
#include <math.h>

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

#ifdef MEMSTR
    #include <io.h>
    #include <fcntl.h>
#endif

namespace oofem {
REGISTER_EngngModel(EigenValueDynamic);


EigenValueDynamic :: EigenValueDynamic(int i, EngngModel *master) : EngngModel(i, master)
{
    numberOfSteps = 1;
    ndomains = 1;
}


NumericalMethod *EigenValueDynamic :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
        nMethod = classFactory.createGeneralizedEigenValueSolver(solverType, this->giveDomain(1), this);
        if ( !nMethod ) {
            OOFEM_ERROR("solver creation failed");
        }
    }

    return nMethod.get();
}


void
EigenValueDynamic :: initializeFrom(InputRecord &ir)
{
    //EngngModel::instanciateFrom (ir);

    IR_GIVE_FIELD(ir, numberOfRequiredEigenValues, _IFT_EigenValueDynamic_nroot);
    this->field = std::make_unique<EigenVectorPrimaryField>(this, 1, FT_Displacements, numberOfRequiredEigenValues);

    // numberOfSteps set artificially to numberOfRequiredEigenValues
    // in order to allow
    // use restoreContext function for different eigenValues
    // numberOfSteps = numberOfRequiredEigenValues;
    numberOfSteps = 1;

    IR_GIVE_FIELD(ir, rtolv, _IFT_EigenValueDynamic_rtolv);
    if ( rtolv < 1.e-12 ) {
        rtolv =  1.e-12;
    } else if ( rtolv > 0.01 ) {
        rtolv =  0.01;
    }

    int val = 1; // inverseit
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EigenValueDynamic_stype);
    solverType = ( GenEigvalSolverType ) val;

    val = 0; //Default Skyline
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = (SparseMtrxType)val;

    if (solverType == GenEigvalSolverType::GES_Eigen)
	sparseMtrxType = SparseMtrxType::SMT_EigenSparse;
    suppressOutput = ir.hasField(_IFT_EngngModel_suppressOutput);

    if ( suppressOutput ) {
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


int EigenValueDynamic :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
{
    return tStep->giveNumber() % this->numberOfRequiredEigenValues; 
}


double EigenValueDynamic :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    return field->giveUnknownValue(dof, mode, tStep);
}


TimeStep *EigenValueDynamic :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    StateCounterType counter = 1;

    if ( currentStep ) {
        istep = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = std :: move(currentStep);
    currentStep = std::make_unique<TimeStep>(istep, this, 1, ( double ) istep, 0., counter);

    return currentStep.get();
}


void EigenValueDynamic :: solveYourself()
{
    this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);

    TimeStep *tStep = this->giveNextStep();
    this->updateAttributes( this->giveCurrentMetaStep() );

    OOFEM_LOG_INFO("Assembling stiffness and mass matrices\n");

    FloatMatrix eigVec;
    std :: unique_ptr< SparseMtrx > stiffnessMatrix;
    std :: unique_ptr< SparseMtrx > massMatrix;

    //if ( tStep->giveNumber() == 1 ) {
    stiffnessMatrix = classFactory.createSparseMtrx(sparseMtrxType);
    stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

    massMatrix = classFactory.createSparseMtrx(sparseMtrxType);
    massMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

    this->assemble( *stiffnessMatrix, tStep, TangentAssembler(TangentStiffness), EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->assemble( *massMatrix, tStep, MassMatrixAssembler(), EModelDefaultEquationNumbering(), this->giveDomain(1) );

    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );
    OOFEM_LOG_INFO("Solving ...\n");
    nMethod->solve(*stiffnessMatrix, *massMatrix, eigVal, eigVec, rtolv, numberOfRequiredEigenValues);
    //}

    // cut off, but output at least one
    int nValid = 1;
    for (int i = 1; i < eigVal.giveSize(); ++i) {
        if (eigVal(i) > 1e10 && i > 0) { nValid = i; break; }
    }
    numberOfRequiredEigenValues = nValid; // better having a dedicated field
    eigVal.resize(numberOfRequiredEigenValues);
    eigVec.resizeWithData(eigVec.giveNumberOfRows(), numberOfRequiredEigenValues);

    this->field->updateAll(eigVec, EModelDefaultEquationNumbering());

    // custom code for output
    FloatMatrix *unitDisp = new FloatMatrix();
    FloatArray *tempCol   = new FloatArray();
    FloatArray *tempCol2  = new FloatArray();

    Domain *domain = this->giveDomain( 1 );
    IntArray dofIDArry, loc;
    dofIDArry = domain->giveDefaultNodeDofIDArry();
    int nelem = domain->giveNumberOfElements();

    // matrix and array initialization
    totMass.resize( 6 ); // 3 are the translational dofs - enough for the moment
    centroid.resize( 3 );
    partFact.resize( numberOfRequiredEigenValues, 6 );
    massPart.resize( numberOfRequiredEigenValues, 6 );
    unitDisp->resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ), 6 );
    tempCol->resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    tempCol2->resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );

    // mass normalization
    for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        tempCol->beColumnOf( eigVec, i );
        massMatrix->times( *tempCol, *tempCol2 );
        double m = tempCol->dotProduct( *tempCol2 );
        if ( m != 0.0 ) m = 1 / sqrt( m );
        tempCol->times( m );
        eigVec.setColumn( *tempCol, i );
    }
    // eigVec has been normalized

    IntArray masterDofIDs, nodalArray, ids;
    IntArray locationArray;
    IntArray *dofIdArray = new IntArray();
    dofIdArray->clear();

    FloatMatrix tempMat, tempMat2;
    FloatArray tempCoord, coordArray;
    tempMat.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ), 6 );
    tempMat2.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ), 6 );

    //
    // create unit displacement vectors
    //
    // first from nodes themselves
    for ( std::unique_ptr<DofManager> &node : domain->giveDofManagers() ) {
        //node->giveLocationArray(dofIDArry, loc, EModelDefaultEquationNumbering());
        if ( !node->giveNumberOfDofs() ) continue;

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

                    unitDisp->at( eqN, dType ) = 1.0;
                    tempMat2.at( eqN, dType )  = node->giveCoordinate( dType );
                }
            }
        }
    } // end of search among nodes

    // then from internaldof managers
    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        Element *element = domain->giveElement( ielem );

        // the following may be simplified.
        // retrieve internal dof managers and location array
        for ( int i = 1; i <= element->giveNumberOfInternalDofManagers(); i++ ) {
            DofManager *intDofMan = element->giveInternalDofManager( i );

            if ( !intDofMan ) continue; // you may never know...

            locationArray.clear();
            tempCoord.clear();

            element->giveInternalDofManDofIDMask( i, ids );
            intDofMan->giveLocationArray( ids, nodalArray, EModelDefaultEquationNumbering() );
            locationArray.followedBy( nodalArray );

            intDofMan->giveMasterDofIDArray( ids, masterDofIDs );
            dofIdArray->followedBy( masterDofIDs );

            coordArray.resize( masterDofIDs.giveSize() );
            int c = 1;
            for ( int dof : masterDofIDs ) {
                if ( intDofMan->giveCoordinates().giveSize() ) {
                    coordArray.at( c ) = intDofMan->giveCoordinate( dof );
                } else {
                    // get the number. ghostNode FEMComponent number is changed on purpose.
                    int tempN = intDofMan->giveNumber();
                    //element->giveDofManager(tempN)->requiresTransformation()
                    coordArray.at( c ) = element->giveDofManager( tempN )->giveCoordinate( dof );
                }
                c++; // Increment the counter in any case. The position index must match the index in masterDofIDs. We'll sort the dofs needed hereafter
            }
            tempCoord.append( coordArray );

            int partialDofCount = locationArray.giveSize();
            if ( partialDofCount ) {
                // search for our dofs in there
                for ( int myDofIndex = 1; myDofIndex <= partialDofCount; myDofIndex++ ) {
                    int dType = dofIdArray->at( myDofIndex );
                    int eqN   = locationArray.at( myDofIndex );

                    if ( ( dType >= D_u ) && ( dType <= D_w ) && eqN ) {
                        // save unit displacement and coordinate

                        unitDisp->at( eqN, dType ) = 1.0;
                        tempMat2.at( eqN, dType )  = tempCoord.at( myDofIndex );
                    }
                }
            }
        }
    } // end of search among internal dof managers
    // end of creation of translational unit displacement vectors

    // if no mass defined in a direction, then centroid in that direction is undefined
    int contr[3];

    for ( int i = 1; i <= 3; i++ ) {
        tempCol->beColumnOf( *unitDisp, i );
        massMatrix->times( *tempCol, *tempCol2 ); // now tempCol2 has only the masses pertaining the i-th direction
        tempMat.setColumn( *tempCol2, i );
        totMass.at( i ) = tempCol->dotProduct( *tempCol2 ); // total mass for direction i-th direction
        tempCol->beColumnOf( tempMat2, i ); // fetch coordinates in i-th direction
        if ( totMass.at( i ) != 0.0 ) {
            centroid.at( i ) = tempCol->dotProduct( *tempCol2 ) / totMass.at( i ); // dot multiply to get first moment, then divide by total mass in i-th direction to get i-th coordinate of the centroid
            contr[i - 1]     = 1;
        } else { contr[i - 1] = 0; }
    }

    // we have the centroid. we can now calculate rotational components. first from nodes.
    for ( std::unique_ptr<DofManager> &node : domain->giveDofManagers() ) {
        //node->giveLocationArray(dofIDArry, loc, EModelDefaultEquationNumbering());
        if ( !node->giveNumberOfDofs() ) continue;

        const FloatArray &nodeCoords = node->giveCoordinates();
        FloatArray vk( 3 );
        IntArray eq( 3 );

        // TODO consider own UCS if present
        //for (int dType = D_u; dType <= D_w; dType++)
        //{
        //	auto myDof = node->findDofWithDofId((DofIDItem)dType);
        //	if (myDof == node->end()){
        //		vk.at(dType) = 0.0;
        //		eq.at(dType) = 0;
        //		continue;
        //	}
        //	vk.at(dType) = node->giveCoordinate(dType) - centroid.at(dType);
        //	eq.at(dType) = EModelDefaultEquationNumbering().giveDofEquationNumber(*myDof);
        //}

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

                    vk.at( dType ) = contr[dType - 1] * ( node->giveCoordinate( dType ) - centroid.at( dType ) );
                    eq.at( dType ) = eqN;
                }
            }
        }

        // set mixed contribution due to rotation about centroid
        if ( eq.at( 1 ) ) {
            unitDisp->at( eq.at( 1 ), 5 ) = vk.at( 3 );
            unitDisp->at( eq.at( 1 ), 6 ) = -vk.at( 2 );
        }

        if ( eq.at( 2 ) ) {
            unitDisp->at( eq.at( 2 ), 4 ) = -vk.at( 3 );
            unitDisp->at( eq.at( 2 ), 6 ) = vk.at( 1 );
        }

        if ( eq.at( 3 ) ) {
            unitDisp->at( eq.at( 3 ), 4 ) = vk.at( 2 );
            unitDisp->at( eq.at( 3 ), 5 ) = -vk.at( 1 );
        }

        // set pure rotational contribution
        //for (int dType = R_u; dType <= R_w; dType++)
        //{
        //	auto myDof = node->findDofWithDofId((DofIDItem)dType);
        //	if (myDof == node->end()) {
        //		//OOFEM_ERROR("incompatible dof (%d) requested", dType);
        //		continue;
        //	}

        //	int eqN = EModelDefaultEquationNumbering().giveDofEquationNumber(*myDof);

        //	// save unit displacement and coordinate
        //	// TODO consider own UCS if present
        //	if (eqN)
        //	{
        //		unitDisp->at(eqN, dType) = 1.0;
        //		tempMat2.at(eqN, dType) = node->giveCoordinate(dType);
        //	}

        //}

        if ( partialDofCount ) {
            // search for our dofs in there
            for ( int myDofIndex = 1; myDofIndex <= partialDofCount; myDofIndex++ ) {
                int dType = mstrDofs.at( myDofIndex );
                int eqN   = locArr.at( myDofIndex );

                if ( ( dType >= R_u ) && ( dType <= R_w ) && eqN ) {
                    // save unit displacement and coordinate

                    unitDisp->at( eqN, dType ) = 1.0;
                    tempMat2.at( eqN, dType )  = node->giveCoordinate( dType );
                }
            }
        }
    } // end of search among nodes

    // then from internaldof managers
    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        Element *element = domain->giveElement( ielem );

        // the following may be simplified
        //	retrieve internal dof managers and location array
        for ( int i = 1; i <= element->giveNumberOfInternalDofManagers(); i++ ) {
            DofManager *intDofMan = element->giveInternalDofManager( i );

            if ( !intDofMan ) continue; // you may never know...

            locationArray.clear();
            tempCoord.clear();

            element->giveInternalDofManDofIDMask( i, ids );
            intDofMan->giveLocationArray( ids, nodalArray, EModelDefaultEquationNumbering() );
            locationArray.followedBy( nodalArray );

            intDofMan->giveMasterDofIDArray( ids, masterDofIDs );
            dofIdArray->followedBy( masterDofIDs );

            coordArray.resize( masterDofIDs.giveSize() );
            int c = 1;
            for ( int dof : masterDofIDs ) {
                if ( intDofMan->giveCoordinates().giveSize() ) {
                    coordArray.at( c ) = intDofMan->giveCoordinate( dof );
                } else {
                    // get the number. ghostNode FEMComponent number is changed on purpose.
                    int tempN = intDofMan->giveNumber();
                    //element->giveDofManager(tempN)->requiresTransformation()
                    coordArray.at( c ) = element->giveDofManager( tempN )->giveCoordinate( dof );
                }
                c++; // Increment the counter in any case. The position index must match the index in masterDofIDs. We'll sort the dofs needed hereafter
            }
            tempCoord.append( coordArray );

            // TODO should element CS be considered even with internal dof managers
            if ( locationArray.giveSize() ) {
                FloatArray vk( 3 );
                IntArray eq( 3 );
                for ( int dType = D_u; dType <= D_w; dType++ ) {
                    int myDof = dofIdArray->findFirstIndexOf( (DofIDItem)dType );
                    if ( myDof == 0 ) {
                        vk.at( dType ) = 0.0;
                        eq.at( dType ) = 0;
                        continue;
                    }
                    vk.at( dType ) = contr[dType - 1] * ( tempCoord.at( myDof ) - centroid.at( dType ) );
                    eq.at( dType ) = locationArray.at( myDof );
                }

                // set mixed contribution due to rotation about centroid
                if ( eq.at( 1 ) ) {
                    unitDisp->at( eq.at( 1 ), 5 ) = vk.at( 3 );
                    unitDisp->at( eq.at( 1 ), 6 ) = -vk.at( 2 );
                }

                if ( eq.at( 2 ) ) {
                    unitDisp->at( eq.at( 2 ), 4 ) = -vk.at( 3 );
                    unitDisp->at( eq.at( 2 ), 6 ) = vk.at( 1 );
                }

                if ( eq.at( 3 ) ) {
                    unitDisp->at( eq.at( 3 ), 4 ) = vk.at( 2 );
                    unitDisp->at( eq.at( 3 ), 5 ) = -vk.at( 1 );
                }


                // search for our dofs in there
                for ( int dType = R_u; dType <= R_w; dType++ ) {
                    int myDof = dofIdArray->findFirstIndexOf( dType );
                    if ( myDof == 0 ) continue;

                    int eqN = locationArray.at( myDof );

                    // save unit displacement and coordinate
                    if ( eqN ) {
                        unitDisp->at( eqN, dType ) = 1.0;
                        tempMat2.at( eqN, dType )  = tempCoord.at( myDof );
                    }
                }
            }
        } // end of search among internal dof managers
    }
    // end of creation of translational unit displacement vectors

    for ( int i = 4; i <= 6; i++ ) {
        tempCol->beColumnOf( *unitDisp, i );
        massMatrix->times( *tempCol, *tempCol2 ); // now tempCol2 has only the masses pertaining the i-th direction
        tempMat.setColumn( *tempCol2, i );
        totMass.at( i ) = tempCol->dotProduct( *tempCol2 ); // total mass for direction i-th direction
        tempCol->beColumnOf( tempMat2, i ); // fetch coordinates in i-th direction
    }

    //
    // calculate participation factors and mass participation
    //

    // participation factors
    partFact.beTProductOf( eigVec, tempMat );

    // mass participation ratios
    for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        for ( int j = 1; j <= 6; j++ ) {
            if ( totMass.at( j ) > 1e-10 ) { massPart.at( i, j ) = pow( partFact.at( i, j ), 2 ) / totMass.at( j ); }
            //else
            //{
            //	massPart.at(i, j) = 0.0;
            //}
        }
    }

    //
    // zero matrix
    //
    stiffnessMatrix.reset( NULL );
    massMatrix.reset( NULL );

    // dispose the rest of the stuff
    delete unitDisp;
    delete tempCol;
    delete tempCol2;
    delete dofIdArray;

    this->terminate( tStep );

    double steptime = this->giveSolutionStepTime();
    OOFEM_LOG_INFO( "EngngModel info: user time consumed by solution: %.2fs\n", steptime );
}


void EigenValueDynamic :: updateYourself(TimeStep *tStep)
{ }


void EigenValueDynamic :: doStepOutput(TimeStep *tStep)
{
    if ( !suppressOutput ) {
        this->printOutputAt(this->giveOutputStream(), tStep);
        fflush( this->giveOutputStream() );
    }

    for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        // export using export manager
        tStep->setTime( ( double ) i ); // we use time as intrinsic eigen value index
        tStep->setNumber(i);
        exportModuleManager.doOutput(tStep);
    }
}


void EigenValueDynamic :: printOutputAt(FILE *file, TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    FILE *outputStream = this->giveOutputStream();

    // print loadcase header
    fprintf( file, "\nOutput for Eigen analysis \n\n", 1.0 );
    // print eigen values on output
    fprintf( file, "\n\nEigen Values (Omega^2) are:\n-----------------\n" );

    for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        fprintf( outputStream, "%15.8e ", eigVal.at( i ) );
        if ( ( i % 5 ) == 0 ) { fprintf( outputStream, "\n" ); }
    }


    fprintf( outputStream, "\n\n\nCentroid Coordinates are:\n-----------------\n\tX\t|\tY\t|\tZ\n" );
    for ( int i = 1; i <= centroid.giveSize(); ++i ) { fprintf( outputStream, "%10.3e ", centroid.at( i ) ); }

    fprintf( outputStream, "\n" );

    fprintf( outputStream, "\n\nParticipation Factors are:\n-----------------\n\tDx\t|\tDy\t|\tDz\t|\tRx\t|\tRy\t|\tRz\n" );
    for ( int i = 1; i <= partFact.giveNumberOfRows(); ++i ) {
        for ( int j = 1; j <= partFact.giveNumberOfColumns(); ++j ) { fprintf( outputStream, "%10.3e ", partFact.at( i, j ) ); }
        fprintf( outputStream, "\n" );
    }

    fprintf( outputStream, "\n\nTotal Masses are:\n-----------------\n\tDx\t|\tDy\t|\tDz\t|\tRx\t|\tRy\t|\tRz\n" );
    for ( int i = 1; i <= totMass.giveSize(); ++i ) { fprintf( outputStream, "%10.3e ", totMass.at( i ) ); }

    fprintf( outputStream, "\n" );

    fprintf( outputStream, "\n\nMass Ratios are:\n-----------------\n\tDx\t|\tDy\t|\tDz\t|\tRx\t|\tRy\t|\tRz\n" );
    for ( int i = 1; i <= massPart.giveNumberOfRows(); ++i ) {
        for ( int j = 1; j <= massPart.giveNumberOfColumns(); ++j ) { fprintf( outputStream, "%10.3e ", massPart.at( i, j ) ); }
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

        if ( this->requiresUnknownsDictionaryUpdate() ) { for ( auto &dman : domain->giveDofManagers() ) { this->updateDofUnknownsDictionary( dman.get(), tStep ); } }


        for ( auto &dman : domain->giveDofManagers() ) {
            dman->updateYourself( tStep );
            dman->printOutputAt( outputStream, tStep );
        }
    }

    fflush( this->giveOutputStream() );

    double utsec = this->timer.getUtime( EngngModelTimer::EMTT_AnalysisTimer );
    fprintf( file, "\nUser time consumed by solution step: %.3f [s]\n\n", utsec );
}


void EigenValueDynamic :: saveContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: saveContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = eigVal.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    this->field->saveContext(stream);
}


void EigenValueDynamic :: restoreContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: restoreContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = eigVal.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    this->field->restoreContext(stream);
}


void EigenValueDynamic :: setActiveVector(int i)
{
    this->activeVector = i;
    if ( activeVector > numberOfRequiredEigenValues ) {
        activeVector = numberOfRequiredEigenValues;
    }

    this->giveCurrentStep()->setNumber( activeVector );
    this->giveCurrentStep()->setTime( ( double ) activeVector );
}

} // end namespace oofem
