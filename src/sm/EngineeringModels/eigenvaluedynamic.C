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
        nValid = i+1;
    }
    numberOfRequiredEigenValues = nValid; // better having a dedicated field
    eigVal.resize(numberOfRequiredEigenValues);
    eigVec.resizeWithData(eigVec.giveNumberOfRows(), numberOfRequiredEigenValues);

    this->field->updateAll(eigVec, EModelDefaultEquationNumbering());

    // custom code for output
    FloatMatrix unitDisp;
    FloatArray tempCol;
    FloatArray tempCol2;

    Domain *domain = this->giveDomain( 1 );
    IntArray dofIDArry, loc;
    dofIDArry = domain->giveDefaultNodeDofIDArry();
    int nelem = domain->giveNumberOfElements();

    // matrix and array initialization
    totMass.resize( 6 ); // 3 are the translational dofs - enough for the moment
    centroid.resize( 3 );
    partFact.resize( numberOfRequiredEigenValues, 6 );
    massPart.resize( numberOfRequiredEigenValues, 6 );
    // matrix of unit displacements for the 6 dofs
    unitDisp.resize(this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()), 6);
    // auxillary vectors for mass normalization
    tempCol.resize(this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()));
    tempCol2.resize(this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()));

    // mass normalization
    for (int i = 1; i <= numberOfRequiredEigenValues; ++i) {
        tempCol.beColumnOf(eigVec, i);
        massMatrix->times(tempCol, tempCol2);
        double m = tempCol.dotProduct(tempCol2);
        if (m != 0.0) m = 1 / sqrt(m);
        tempCol.times(m);
        eigVec.setColumn(tempCol, i);
    }
    // eigVec has been normalized

    IntArray masterDofIDs, nodalEqArray, ids;
    IntArray eqArray;
    IntArray dofManArray;

    FloatMatrix tempMat, tempMat2;
    tempMat.resize(this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()), 6);
    tempMat2.resize(this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()), 6);

    static const DofIDItem dofids[] = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };

    const EModelDefaultEquationNumbering defNumbering;
    const int nDofMans = domain->giveNumberOfDofManagers();
    bool warn = false;
    FloatMatrix nodeCoords(3, nDofMans);
    std::vector<IntArray> nodeActive(nDofMans);

    FloatArray geomCenter(3);
    IntArray geomWeight(3);
    IntArray coordFilter(3); // if no mass or free dof is defined in a direction, then centroid in that direction is undefined
    // cache dofs and node coordinates, unit displacements for translational dofs
    for (std::unique_ptr<DofManager>& node : domain->giveDofManagers()) {
        // no support yet for LCS
        Node* actualNode = dynamic_cast<Node*>(node.get());
        if (actualNode && actualNode->hasLocalCS()) {
            if (!warn) {
                OOFEM_WARNING("Nodes with LCS are unsupported, use at own risk.");
                warn = true;
            }
            continue;
        }
        if (!node->giveNumberOfDofs()) {
            continue;
        }

        const int num = node->giveNumber();

        const auto& coords = node->giveCoordinates();
        IntArray activeDofs(3);
        if (strcmp(node->giveClassName(), "Node") == 0 || strcmp(node->giveClassName(), "RigidArmNode") == 0) {
            for (int iDof = 0; iDof < 3; ++iDof) {
                DofIDItem dType = dofids[iDof];
                auto pos = node->findDofWithDofId(dType);
                if (pos == node->end()) {
                    continue;
                }

                int eqN = 0;
                if ((*pos)->isPrimaryDof()) {
                    eqN = (*pos)->giveEquationNumber(defNumbering);
                    //if ( eqN ) {
                    //    unitDisp.at( eqN, dType ) = 1.0;
                    //}
                }
                else {
                    // find the master dofmanager from the slave, get the equations from it.
                    // Plain Nodes have master and simpleslave dofs, rigid arm nodes have master and slave dofs
                    IntArray masterDofMans;
                    (*pos)->giveMasterDofManArray(masterDofMans);

                    auto* masterMan = domain->giveDofManager(masterDofMans.at(1)); // there's only 1
                    auto masterDof = masterMan->findDofWithDofId(dType);
                    eqN = (*masterDof)->giveEquationNumber(defNumbering);
                }

                activeDofs.at(dType) = eqN > 0;
            }
        }
        nodeActive[num - 1] = activeDofs;

        // save coordinate and increment the weight
        nodeCoords.setColumn(coords, num);
        for (int iDof = 1; iDof <= 3; ++iDof) {
            if (activeDofs.at(iDof)) {
                geomCenter.at(iDof) += coords.at(iDof);
                geomWeight.at(iDof) += activeDofs.at(iDof);
            }
        }
    }

    warn = false;
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
        const int dofMan = massElement->giveDofManArray().at(1);
        const IntArray loc;
        lma.vectorFromElement(m, *massElement, tStep, VM_Unknown);

        IntArray dofMask;
        massElement->giveElementDofIDMask(dofMask);

        for (auto iDof : dofMask) {
            if (iDof >= D_u && iDof <= D_w) {
                const int idx = dofMask.findFirstIndexOf(iDof);
                if (idx && nodeActive[dofMan - 1].at(iDof)) {
                    massPos.at(iDof) += nodeCoords.at(iDof, dofMan) * m.at(idx);
                    totMass.at(iDof) += m.at(idx);
                }
            }
        }
    }

    // center of geometry
    for (int i = 1; i <= 3; ++i) {
        if (geomWeight.at(i)) {
            geomCenter.at(i) /= geomWeight.at(i);
            coordFilter.at(i) |= 1;
        }
        else {
            OOFEM_WARNING("No free dof in direction %d, results may be affected", i);
        }
    }
    // center of mass
    for (int i = 1; i <= 3; ++i) {
        if (totMass.at(i) > 0) {
            centroid.at(i) = massPos.at(i) / totMass.at(i);
            coordFilter.at(i) |= 1;
        }
        else {
            OOFEM_WARNING("No mass in direction %d, results may be affected", i);
            centroid.at(i) = geomCenter.at(i);
        }
    }

    //
    // create unit displacement vectors
    //
    // first from nodes themselves, and get the coordinates to compute the center of mass
    for (std::unique_ptr<DofManager>& node : domain->giveDofManagers()) {
        // no support yet for LCS
        Node* actualNode = dynamic_cast<Node*>(node.get());
        if (actualNode && actualNode->hasLocalCS()) {
            continue;
        }
        if (!node->giveNumberOfDofs()) continue;

        IntArray mstrDofs, locArr;
        node->givePrimaryDofs(mstrDofs);
        node->giveLocationArray(mstrDofs, locArr, EModelDefaultEquationNumbering());

        // search for our dofs in there
        for (int myDofIndex = 1; myDofIndex <= locArr.giveSize(); ++myDofIndex) {
            int dType = mstrDofs.at(myDofIndex);
            int eqN = locArr.at(myDofIndex);

            if ((dType >= D_u) && (dType <= D_w) && eqN) {
                // save unit displacement and coordinate
                unitDisp.at(eqN, dType) = 1.0;
                //tempMat2.at(eqN, dType) = node->giveCoordinate(dType);
            }
        }
    } // end of search among nodes

    warn = false;
    // then from internaldof managers
    for (int ielem = 1; ielem <= nelem; ++ielem) {
        Element* element = domain->giveElement(ielem);
        if (element->giveNumberOfInternalDofManagers() == 0) {
            continue;
        }

        // we only support end releases on 3d beams for the moment. other internal dof managers should be skipped.
        if (strcmp(element->giveClassName(), "Beam3d") != 0) {
            if (!warn) {
                OOFEM_WARNING("Only internal dof managers of Beam3d elements are supported. Skipping.");
                warn = true;
            }
            continue;
        }

        // use the lcs transform to work out which ghost dofs are affected by global displacement
        FloatMatrix GtoL;
        element->computeGtoLRotationMatrix(GtoL);
        FloatMatrix lcs;
        element->giveLocalCoordinateSystem(lcs);

        eqArray.clear();

        // the following may be simplified.
        // retrieve internal dof managers and location array
        for (int i = 1; i <= element->giveNumberOfInternalDofManagers(); ++i) {
            DofManager* intDofMan = element->giveInternalDofManager(i);
            if (!intDofMan) continue; // you may never know...

            element->giveInternalDofManDofIDMask(i, ids);
            intDofMan->giveLocationArray(ids, nodalEqArray, EModelDefaultEquationNumbering());
            eqArray.followedBy(nodalEqArray);
            intDofMan->giveMasterDofIDArray(ids, masterDofIDs);

            int tempN = intDofMan->giveNumber();
            const auto& coords = element->giveDofManager(tempN)->giveCoordinates();
            //for (auto eqN : nodalEqArray) {
            //    if (eqN == 0) {
            //        continue;
            //    }
            //    tempMat2.addSubVectorRow(coords, eqN, 1);
            //}
        }

        const int nBaseDofs = element->giveNumberOfDofManagers() * 6;
        for (int iDof = 0; iDof < 3; ++iDof) {
            const auto dType = dofids[iDof];
            FloatArray localVector(nBaseDofs);
            FloatArray globalVector;
            FloatArray dirVector;
            dirVector.beColumnOf(lcs, iDof % 3 + 1);
            // get unit vectors in local dofs
            for (int inode = 0; inode < element->giveNumberOfDofManagers(); ++inode) {
                localVector.addSubVector(dirVector, (iDof / 3) * 3 + inode * 6 + 1);
            }
            // work it back to global dofs
            globalVector.beTProductOf(GtoL, localVector);

            // apply the components to the condensed dofs on the global matrix
            for (int i = 1; i <= eqArray.giveSize(); ++i) {
                const auto eqN = eqArray.at(i);
                const auto val = globalVector(i + nBaseDofs - 1);
                unitDisp.at(eqN, dType) = val;
            }
        }


    } // end of search among internal dof managers
    // end of translational unit vectors among internal dof managers

    for (int i = 1; i <= 3; ++i) {
        tempCol.beColumnOf(unitDisp, i);
        massMatrix->times(tempCol, tempCol2); // now tempCol2 has only the masses pertaining the i-th direction
        tempMat.setColumn(tempCol2, i);
        totMass.at(i) = tempCol.dotProduct(tempCol2);
    }

    // we have the centroid. we can now calculate rotational components. first from nodes.
    for (std::unique_ptr<DofManager>& node : domain->giveDofManagers()) {
        if (!node->giveNumberOfDofs()) continue;

        FloatArray vk(3);
        IntArray eq(3);

        // TODO consider own UCS if present

        IntArray mstrDofs, locArr;
        node->givePrimaryDofs(mstrDofs);
        node->giveLocationArray(mstrDofs, locArr, EModelDefaultEquationNumbering());

        int partialDofCount = locArr.giveSize();
        if (partialDofCount) {
            // search for our dofs in there
            for (int myDofIndex = 1; myDofIndex <= partialDofCount; myDofIndex++) {
                int dType = mstrDofs.at(myDofIndex);
                int eqN = locArr.at(myDofIndex);

                if ((dType >= D_u) && (dType <= D_w) && eqN) {
                    // save unit displacement and coordinate
                    vk.at( dType ) = coordFilter( dType ) * ( node->giveCoordinate( dType ) - centroid.at( dType ) );
                    eq.at( dType ) = eqN;
                }
            }
        }

        // set mixed contribution due to rotation about centroid
        if (eq.at(1)) {
            unitDisp.at(eq.at(1), R_v) = vk.at(3);
            unitDisp.at(eq.at(1), R_w) = -vk.at(2);
        }

        if (eq.at(2)) {
            unitDisp.at(eq.at(2), R_u) = -vk.at(3);
            unitDisp.at(eq.at(2), R_w) = vk.at(1);
        }

        if (eq.at(3)) {
            unitDisp.at(eq.at(3), R_u) = vk.at(2);
            unitDisp.at(eq.at(3), R_v) = -vk.at(1);
        }

        if (partialDofCount) {
            // search for our dofs in there
            for (int myDofIndex = 1; myDofIndex <= partialDofCount; myDofIndex++) {
                int dType = mstrDofs.at(myDofIndex);
                int eqN = locArr.at(myDofIndex);

                if ((dType >= R_u) && (dType <= R_w) && eqN) {
                    // save unit displacement and coordinate
                    unitDisp.at(eqN, dType) = 1.0;
                }
            }
        }
    } // end of search among nodes

    // then from internaldof managers
    for (int ielem = 1; ielem <= nelem; ++ielem) {
        Element* element = domain->giveElement(ielem);
        if (!element->giveNumberOfInternalDofManagers()) {
            continue;
        }

        // we only support end releases on 3d beams for the moment. other internal dof managers should be skipped.
        if (strcmp(element->giveClassName(), "Beam3d") != 0) {
            continue;
        }

        eqArray.clear();

        // use the lcs transform to work out which ghost dofs are affected by global displacement
        FloatMatrix GtoL;
        element->computeGtoLRotationMatrix(GtoL);
        FloatMatrix lcs;
        element->giveLocalCoordinateSystem(lcs);

        FloatMatrix componentMatrix(3, 3);
        //	retrieve internal dof managers and location array
        for (int i = 1; i <= element->giveNumberOfInternalDofManagers(); ++i) {
            DofManager* intDofMan = element->giveInternalDofManager(i);
            if (!intDofMan) continue; // you may never know...

            element->giveInternalDofManDofIDMask(i, ids);
            intDofMan->giveLocationArray(ids, nodalEqArray, EModelDefaultEquationNumbering());
            eqArray.followedBy(nodalEqArray);
            intDofMan->giveMasterDofIDArray(ids, masterDofIDs);

            int tempN = intDofMan->giveNumber();
            const auto& coords = element->giveDofManager(tempN)->giveCoordinates();
            componentMatrix.setColumn(coords - centroid, 2);

            for (int iDof = 1; iDof <= masterDofIDs.giveSize(); ++iDof) {
                const int myDof = masterDofIDs.at(iDof);
                // now we're only interested about the effect of rotation on translational dofs
                if (myDof > D_w) { continue; }
                const int eqN = nodalEqArray.at(iDof);
                for (int j = 1; j <= 3; ++j) {
                    componentMatrix.at(j, 3) = lcs.at(myDof, j);
                }
                // project the displacement component due to rotation onto the local direction of the dof
                for (int rotDof = R_u; rotDof <= R_w; ++rotDof) {
                    FloatArray rotationVector(3);
                    rotationVector.at(rotDof - 3) = 1.0;
                    componentMatrix.setColumn(rotationVector, 1);
                    const double displacement = componentMatrix.giveDeterminant(); // triple product: rotation ^ (node-centroid) . dofVector
                    unitDisp.at(eqN, rotDof) = displacement;
                }
            }
        } // end of search among internal dof managers

        const int nBaseDofs = element->giveNumberOfDofManagers() * 6;
        for (int iDof = 0; iDof < 6; ++iDof) {
            const auto dType = dofids[iDof];
            if (dType < R_u || dType > R_w) {
                continue;
            }

            FloatArray localVector(nBaseDofs);
            FloatArray globalVector;
            FloatArray dirVector;
            dirVector.beColumnOf(lcs, iDof % 3 + 1);
            // get unit vectors in local dofs
            for (int inode = 0; inode < element->giveNumberOfDofManagers(); ++inode) {
                localVector.addSubVector(dirVector, (iDof / 3) * 3 + inode * 6 + 1);
            }
            // work it back to global dofs
            globalVector.beTProductOf(GtoL, localVector);

            // apply the components to the condensed dofs on the global matrix
            for (int i = 1; i <= eqArray.giveSize(); ++i) {
                const auto eqN = eqArray.at(i);
                const auto val = globalVector(i + nBaseDofs - 1);
                unitDisp.at(eqN, dType) = val;
            }
        }
    }
    // end of creation of translational unit displacement vectors

    for (int i = 4; i <= 6; ++i) {
        tempCol.beColumnOf(unitDisp, i);
        massMatrix->times(tempCol, tempCol2); // now tempCol2 has only the masses pertaining the i-th direction
        tempMat.setColumn(tempCol2, i);
        totMass.at(i) = tempCol.dotProduct(tempCol2); // total mass for direction i-th direction
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
