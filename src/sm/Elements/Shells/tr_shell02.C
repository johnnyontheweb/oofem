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

#include "sm/Elements/Shells/tr_shell02.h"
#include "fei2dtrlin.h"
#include "contextioerr.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "node.h"
#include "load.h"

#ifdef __OOFEG
 #include "node.h"
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif


namespace oofem {
REGISTER_Element(TR_SHELL02);

IntArray TR_SHELL02 :: loc_plate = {3, 4, 5, 9, 10, 11, 15, 16, 17};
IntArray TR_SHELL02 :: loc_membrane = {1, 2, 6, 7, 8, 12, 13, 14, 18};

TR_SHELL02 :: TR_SHELL02(int n, Domain *aDomain) : StructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this), ZZErrorEstimatorInterface(this), SpatialLocalizerInterface(this),
    plate(std::make_unique<DKTPlate3d>(-1, aDomain)),
    membrane(std::make_unique<TrPlanestressRotAllman3d>(-1, aDomain))
{
    numberOfDofMans = 3;
}


void
TR_SHELL02 :: initializeFrom(InputRecord &ir)
{
    // proc tady neni return = this...   ??? termitovo
    StructuralElement :: initializeFrom(ir);

    // optional record for 1st local axes
    la1.resize(3);
    la1.at(1) = 0; la1.at(2) = 0; la1.at(3) = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->la1, _IFT_TR_SHELL02_FirstLocalAxis);

    this->macroElem = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->macroElem, _IFT_TR_SHELL02_macroElem);

    plate->initializeFrom(ir);
    plate->la1 = la1;
    membrane->initializeFrom(ir);
    membrane->la1 = la1;
}

void
TR_SHELL02 :: postInitialize()
{
    StructuralElement :: postInitialize();
    this->numberOfGaussPoints = plate->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
    if ( plate->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints() != membrane->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints() ) {
        OOFEM_ERROR("incompatible integration rules detected");
    }
}

void
TR_SHELL02 :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    StructuralElement :: updateLocalNumbering(f);
    plate->updateLocalNumbering(f);
    membrane->updateLocalNumbering(f);
}

void TR_SHELL02 :: setCrossSection(int csIndx)
{
    StructuralElement :: setCrossSection(csIndx);
    plate->setCrossSection(csIndx);
    membrane->setCrossSection(csIndx);
}



void
TR_SHELL02 :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep)
//
// returns characteristics vector of receiver accordind to mtrx
//
{
    FloatArray aux;

    answer.resize(18);
    answer.zero();

    plate->giveCharacteristicVector(aux, mtrx, mode, tStep);
    if ( !aux.isEmpty() ) answer.assemble(aux, loc_plate);

    membrane->giveCharacteristicVector(aux, mtrx, mode, tStep);
    if ( !aux.isEmpty() ) answer.assemble(aux, loc_membrane);
}

void
TR_SHELL02 :: giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver accordind to mtrx
//
{
    FloatMatrix aux;

    answer.resize(18, 18);
    answer.zero();

    plate->giveCharacteristicMatrix(aux, mtrx, tStep);
    if ( aux.isNotEmpty() ) answer.assemble(aux, loc_plate);

    membrane->giveCharacteristicMatrix(aux, mtrx, tStep);
    if ( aux.isNotEmpty() ) answer.assemble(aux, loc_membrane);
}

void
TR_SHELL02::giveNodeCoordinates(FloatArray& nc1, FloatArray& nc2, FloatArray& nc3)
{
    nc1.resize(3);
    nc2.resize(3);
    nc3.resize(3);

	this->giveLocalCoordinates(nc1, this->giveNode(1)->giveCoordinates());
	this->giveLocalCoordinates(nc2, this->giveNode(2)->giveCoordinates());
	this->giveLocalCoordinates(nc3, this->giveNode(3)->giveCoordinates());
}

void
TR_SHELL02::giveLocalCoordinates(FloatArray &answer, const FloatArray &global)
// Returns global coordinates given in global vector
// transformed into local coordinate system of the
// receiver
{
	FloatArray offset;
	// test the parametr
	if (global.giveSize() != 3) {
		OOFEM_ERROR("cannot transform coordinates - size mismatch");
		exit(1);
	}

	//this->computeGtoLRotationMatrix();

	offset = global;
	offset.subtract(this->giveNode(1)->giveCoordinates());
	answer.beProductOf(*this->plate->computeGtoLRotationMatrix(), offset);
}

bool
TR_SHELL02 :: giveRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix aux1, aux2;
    int ncol;

    bool t1 = plate->giveRotationMatrix(aux1);
    bool t2 =  membrane->giveRotationMatrix(aux2);

    if ( t1 != t2 ) {
        OOFEM_ERROR("Transformation demand mismatch");
    }

    if ( t1 ) {
        ncol = aux1.giveNumberOfColumns();
        answer.resize(18, ncol);

        for ( int i = 1; i <= 9; i++ ) { // row index
            for ( int j = 1; j <= ncol; j++ ) {
                answer.at(loc_plate.at(i), j) = aux1.at(i, j);
            }
        }

        for ( int i = 1; i <= 9; i++ ) { // row index
            for ( int j = 1; j <= ncol; j++ ) {
                answer.at(loc_membrane.at(i), j) = aux2.at(i, j);
            }
        }
    }

    return t1;
}

void
TR_SHELL02::computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
{
	answer.resize(18, 18);
	answer.zero();

	// matrix assembly vectors
	IntArray asmz{ 3, 9, 15 };  // local z

	// stress vector
	FloatArray str, str2;
	//this->giveCharacteristicVector(str, InternalForcesVector, VM_Total, tStep);
    double wsum = 0;
	for (GaussPoint *gp : *this->giveDefaultIntegrationRulePtr()) {
		this->giveIPValue(str2, gp, IST_ShellForceTensor, tStep);
        str2.times( gp->giveWeight() ); wsum += gp->giveWeight();
		str.add(str2);
	}
	str.times(1.0 / wsum); // the weights sum up to 0.5 for a tria
	// this needs to be transformed to local
	FloatMatrix strmat{ 3, 3 };
	strmat.at(1, 1) = str.at(1); strmat.at(2, 2) = str.at(2); strmat.at(3, 3) = str.at(3);
	strmat.at(1, 2) = str.at(6); strmat.at(1, 3) = str.at(5); strmat.at(2, 3) = str.at(4);
	strmat.symmetrized();
	strmat.rotatedWith(*this->plate->computeGtoLRotationMatrix(), 't'); // back to local

	// if above are local and Forces by unit length
	// first average them
	double sx = 0, sy = 0, sxy = 0;
	//for (auto it : asm1) sx += str.at(it);
	//for (auto it : asm2) sy += str.at(it);
	//for (auto it : asm3) sxy += str.at(it);
	//sx = str.at(1);
	//sy = str.at(2);
	//sxy = str.at(6);
	sx = strmat.at(1, 1);
	sy = strmat.at(2, 2);
	sxy = strmat.at(1, 2);
	// then divide by 4*A
	double area = plate->computeArea();
	sx /= (4 * area);
	sy /= (4 * area);
	sxy /= (4 * area);

	// partial matrices
	FloatMatrix Kgx{ 3, 3 }, Kgy{ 3, 3 }, Kgxy{ 3, 3 };

	// calculate the matrices - hp constant thickness
	FloatArray n1, n2, n3;
	this->giveNodeCoordinates(n1, n2, n3);

	double x12, x23, x13, y12, y23, y13;
	x12 = n1.at(1) - n2.at(1);
	y12 = n1.at(2) - n2.at(2);
	x13 = n1.at(1) - n3.at(1);
	y13 = n1.at(2) - n3.at(2);
	x23 = n2.at(1) - n3.at(1);
	y23 = n2.at(2) - n3.at(2);

	Kgx.at( 1, 1 ) = y23 * y23;
        Kgx.at( 2, 2 ) = y13 * y13;
        Kgx.at( 3, 3 ) = y12 * y12;
        Kgx.at( 1, 2 ) = y23 * -y13;
        Kgx.at( 1, 3 ) = y23 * y12;
        Kgx.at( 2, 3 ) = -y13 * y12;

        Kgy.at( 1, 1 ) = x23 * x23;
        Kgy.at( 2, 2 ) = x13 * x13;
        Kgy.at( 3, 3 ) = x12 * x12;
        Kgy.at( 1, 2 ) = -x23 * x13;
        Kgy.at( 1, 3 ) = x23 * x12;
        Kgy.at( 2, 3 ) = x13 * -x12;

        Kgxy.at( 1, 1 ) = 2 * y23 * -x23;
        Kgxy.at( 2, 2 ) = 2 * -y13 * x13;
        Kgxy.at( 3, 3 ) = 2 * y12 * -x12;
        Kgxy.at( 1, 2 ) = y23 * x13 + x23 * y13;
        Kgxy.at( 1, 3 ) = y23 * -x12 + -x23 * y12;
        Kgxy.at( 2, 3 ) = y13 * x12 + x13 * y12;

	Kgx.symmetrized(); Kgy.symmetrized(); Kgxy.symmetrized();

	// assemble them
	Kgx.times(sx);
	Kgx.add(sy, Kgy);
	Kgx.add(sxy, Kgxy);
	// once for local z
	answer.assemble(Kgx, asmz);
	// done
}

void
TR_SHELL02 :: updateInternalState(TimeStep *tStep)
// Updates the receiver at end of step.
{
    plate->updateInternalState(tStep);
    membrane->updateInternalState(tStep);
}

void
TR_SHELL02 :: updateYourself(TimeStep *tStep)
{
    StructuralElement :: updateYourself(tStep);

    plate->updateYourself(tStep);
    membrane->updateYourself(tStep);
}


Interface *
TR_SHELL02 :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    } else if ( interface == ZZErrorEstimatorInterfaceType ) {
        return static_cast< ZZErrorEstimatorInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    }


    return NULL;
}

double
TR_SHELL02 :: computeVolumeAround(GaussPoint *gp)
{
    return plate->computeVolumeAround(gp);
}

void
TR_SHELL02 :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
{
    // OOFEM_ERROR("This function is not implemented yet.");
    FloatArray aux;

    answer.resize(18);
    answer.zero();

    plate->computeBodyLoadVectorAt(aux, forLoad, tStep, mode);
    if (!aux.isEmpty()) answer.assemble(aux, loc_plate);

    membrane->computeBodyLoadVectorAt(aux, forLoad, tStep, mode);
    if (!aux.isEmpty()) answer.assemble(aux, loc_membrane);
}

int
TR_SHELL02 :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_ShellForceTensor || type == IST_ShellStrainTensor ||
        type == IST_ShellMomentTensor || type == IST_CurvatureTensor ) {
        FloatArray aux;
        GaussPoint *membraneGP = membrane->giveDefaultIntegrationRulePtr()->getIntegrationPoint(gp->giveNumber() - 1);
        GaussPoint *plateGP = plate->giveDefaultIntegrationRulePtr()->getIntegrationPoint(gp->giveNumber() - 1);

        plate->giveIPValue(answer, plateGP, type, tStep);
	if (type == IST_ShellMomentTensor || type == IST_CurvatureTensor) { answer.negated(); } // tension at bottom is positive bending
        membrane->giveIPValue(aux, membraneGP, type, tStep);
        answer.add(aux);
        return 1;
    } else if (type == IST_StrainTensor || type == IST_StressTensor) {
	    answer.resize(6);
	    answer.zero();
	    return 1;
    } else {
        return StructuralElement :: giveIPValue(answer, gp, type, tStep);
    }
}


//
// The element interface required by NodalAveragingRecoveryModel
//
void
TR_SHELL02 :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                         InternalStateType type, TimeStep *tStep)
{
    //if (node < 1 || node>3) {
    //	OOFEM_ERROR("unsupported node");
    //}

    //int size = 0;
    //FloatArray val;
    //FloatArray Nweights;
    //FEIElementGeometryWrapper wrap(this);
    //for (int i = 0; i < plate->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints(); i++) {
    //	GaussPoint *gpP = plate->giveDefaultIntegrationRulePtr()->getIntegrationPoint(i);
    //	//GaussPoint *gpM = membrane->giveDefaultIntegrationRulePtr()->getIntegrationPoint(i);
    //	giveIPValue(val, gpP, type, tStep);
    //	if (size == 0) {
    //		size = val.giveSize();
    //		answer.resize(size);
    //	}
    //	plate->giveInterpolation()->evalN(Nweights, gpP->giveNaturalCoordinates(), wrap);
    //	//membrane->giveInterpolation()->evalN(Nweights, gpM->giveNaturalCoordinates(), wrap);
    //	answer.add(Nweights.at(node), val);
    //}

    // choose least squares or not
    if (this->numberOfGaussPoints == 1) {
	this->giveIPValue(answer, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0), type, tStep);
    }else{
	double x1 = 0.0, x2 = 0.0, y = 0.0;
	FloatMatrix A(3, 3);
	FloatMatrix b, r;
	FloatArray val;
	double u, v;

	int size = 0;

	//for (GaussPoint *gp : *integrationRulesArray[0]) {
	for (int i = 0; i < plate->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints(); i++) {
	    GaussPoint *gp = plate->giveDefaultIntegrationRulePtr()->getIntegrationPoint(i);
	    giveIPValue(val, gp, type, tStep);
	    if (size == 0) {
		size = val.giveSize();
		b.resize(3, size);
		r.resize(3, size);
		A.zero();
		r.zero();
	    }

	    const FloatArray &coord = gp->giveNaturalCoordinates();
	    u = coord.at(1);
	    v = coord.at(2);

	    A.at(1, 1) += 1;
	    A.at(1, 2) += u;
	    A.at(1, 3) += v;
	    A.at(2, 1) += u;
	    A.at(2, 2) += u * u;
	    A.at(2, 3) += u * v;
	    A.at(3, 1) += v;
	    A.at(3, 2) += v * u;
	    A.at(3, 3) += v * v;

	    for (int j = 1; j <= size; j++) {
		y = val.at(j);
		r.at(1, j) += y;
		r.at(2, j) += y * u;
		r.at(3, j) += y * v;
	    }
	}

	A.solveForRhs(r, b);

	switch (node) {
	case 1:
	    x1 = 0.0;
	    x2 = 0.0;
	    break;
	case 2:
	    x1 = 1.0;
	    x2 = 0.0;
	    break;
	case 3:
	    x1 = 0.0;
	    x2 = 1.0;
	    break;
	default:
	    OOFEM_ERROR("unsupported node");
	}

	answer.resize(size);
	for (int j = 1; j <= size; j++) {
	    answer.at(j) = b.at(1, j) + x1 *b.at(2, j) + x2 *b.at(3, j);
	}
    }
}


void
TR_SHELL02 :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    fprintf(file, "element %d (%8d) macroelem %d :\n", this->giveLabel(), this->giveNumber(), this->macroElem);

//#ifdef DEBUG
//    FloatArray v, aux;
//    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
//    //fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );
//
//    for ( GaussPoint *gp: *iRule ) {
//        fprintf( file, "  GP %2d.%-2d :", iRule->giveNumber(), gp->giveNumber() );
//        GaussPoint *membraneGP = membrane->giveDefaultIntegrationRulePtr()->getIntegrationPoint(gp->giveNumber() - 1);
//        // Strain - Curvature
//        plate->giveIPValue(v, gp, IST_ShellStrainTensor, tStep);
//        membrane->giveIPValue(aux, membraneGP, IST_ShellStrainTensor, tStep);
//        v.add(aux);
//
//        fprintf(file, "  strains    ");
//        for ( auto &val : v ) fprintf(file, " %.4e", val);
//
//        plate->giveIPValue(v, gp, IST_CurvatureTensor, tStep);
//        membrane->giveIPValue(aux, membraneGP, IST_CurvatureTensor, tStep);
//        v.add(aux);
//
//        fprintf(file, "\n              curvatures ");
//        for ( auto &val : v ) fprintf(file, " %.4e", val);
//
//        // Forces - Moments
//        plate->giveIPValue(v, gp, IST_ShellForceTensor, tStep);
//        membrane->giveIPValue(aux, membraneGP, IST_ShellForceTensor, tStep);
//        v.add(aux);
//
//        fprintf(file, "\n              stresses   ");
//        for ( auto &val : v ) fprintf(file, " %.4e", val);
//
//        plate->giveIPValue(v, gp, IST_ShellMomentTensor, tStep);
//        membrane->giveIPValue(aux, membraneGP, IST_ShellMomentTensor, tStep);
//        v.add(aux);
//
//        fprintf(file, "\n              moments    ");
//        for ( auto &val : v ) fprintf(file, " %.4e", val);
//
//        fprintf(file, "\n");
//    }
//#endif
}

void
TR_SHELL02::computeEdgeNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords)
{
    FloatArray n_vec;
    this->giveInterpolation()->boundaryEdgeEvalN(n_vec, boundaryID, lcoords, FEIElementGeometryWrapper(this));
    answer.beNMatrixOf(n_vec, 6);
}

void
TR_SHELL02::computeSurfaceNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords)
{
    FloatArray n_vec;
    this->giveInterpolation()->boundarySurfaceEvalN(n_vec, boundaryID, lcoords, FEIElementGeometryWrapper(this));
    answer.beNMatrixOf(n_vec, 6);
}

double
TR_SHELL02::computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    return this->computeVolumeAround(gp);
}

double
TR_SHELL02::computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    std::vector< FloatArray >lc = {
	FloatArray(3), FloatArray(3), FloatArray(3)
    };
    this->giveNodeCoordinates(lc[0], lc[1], lc[2]);

    double detJ = this->plate->interp_lin.edgeGiveTransformationJacobian(iEdge, gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(lc));
    return detJ * gp->giveWeight();
}


void
TR_SHELL02 :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralElement :: saveContext(stream, mode);
    this->plate->saveContext(stream, mode);
    this->membrane->saveContext(stream, mode);
}


void
TR_SHELL02 :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralElement :: restoreContext(stream, mode);
    this->plate->restoreContext(stream, mode);
    this->membrane->restoreContext(stream, mode);
}


IntegrationRule *
TR_SHELL02 :: ZZErrorEstimatorI_giveIntegrationRule()
{
    if ( !this->compositeIR ) {
        this->compositeIR = std::make_unique<GaussIntegrationRule>(1, this, 1, 12);
        this->compositeIR->SetUpPointsOnTriangle(plate->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints(), _3dShell);
    }
    return this->compositeIR.get();
}

void
TR_SHELL02 :: ZZErrorEstimatorI_computeLocalStress(FloatArray &answer, FloatArray &sig)
{
    // sig is global ShellForceMomentTensor
    FloatMatrix globTensor(3, 3);
    const FloatMatrix *GtoLRotationMatrix = plate->computeGtoLRotationMatrix();
    FloatMatrix LtoGRotationMatrix;

    answer.resize(8); // reduced, local form
    LtoGRotationMatrix.beTranspositionOf(* GtoLRotationMatrix);

    // Forces
    globTensor.at(1, 1) = sig.at(1);  //sxForce
    globTensor.at(1, 2) = sig.at(6);  //qxyForce
    globTensor.at(1, 3) = sig.at(5);  //qxzForce

    globTensor.at(2, 1) = sig.at(6);  //qxyForce
    globTensor.at(2, 2) = sig.at(2);  //syForce
    globTensor.at(2, 3) = sig.at(4);  //syzForce

    globTensor.at(3, 1) = sig.at(5);  //qxzForce
    globTensor.at(3, 2) = sig.at(4);  //syzForce
    globTensor.at(3, 3) = sig.at(3);  //szForce

    globTensor.rotatedWith(LtoGRotationMatrix);
    // Forces: now globTensoris transformed into local c.s

    // answer should be in reduced, local  form
    answer.at(1) = globTensor.at(1, 1); //sxForce
    answer.at(2) = globTensor.at(2, 2); //syForce
    answer.at(3) = globTensor.at(1, 2); //qxyForce
    answer.at(7) = globTensor.at(2, 3); //syzForce
    answer.at(8) = globTensor.at(1, 3); //qxzForce


    // Moments:
    globTensor.at(1, 1) = sig.at(7);  //mxForce
    globTensor.at(1, 2) = sig.at(12); //mxyForce
    globTensor.at(1, 3) = sig.at(11); //mxzForce

    globTensor.at(2, 1) = sig.at(12); //mxyForce
    globTensor.at(2, 2) = sig.at(8);  //myForce
    globTensor.at(2, 3) = sig.at(10); //myzForce

    globTensor.at(3, 1) = sig.at(11); //mxzForce
    globTensor.at(3, 2) = sig.at(10); //myzForce
    globTensor.at(3, 3) = sig.at(9);  //mzForce

    globTensor.rotatedWith(LtoGRotationMatrix);
    // now globTensoris transformed into local c.s

    answer.at(4)  = globTensor.at(1, 1); //mxForce
    answer.at(5)  = globTensor.at(2, 2); //myForce
    answer.at(6) = globTensor.at(1, 2); //mxyForce
}


void
TR_SHELL02 :: SpatialLocalizerI_giveBBox(FloatArray &bb0, FloatArray &bb1)
{
    FloatArray lt3, gt3; // global vector in the element thickness direction of lenght thickeness/2
    const FloatMatrix *GtoLRotationMatrix = plate->computeGtoLRotationMatrix();

    // setup vector in the element local cs. perpendicular to element plane of thickness/2 length
    lt3 = {0., 0., 1.}; //this->giveCrossSection()->give(CS_Thickness)/2.0; // HUHU
    // transform it to globa cs
    gt3.beTProductOf(* GtoLRotationMatrix, lt3);

    // use gt3 to construct element bounding box respecting true element volume

    FloatArray _c;

    for ( int i = 1; i <= this->giveNumberOfNodes(); ++i ) {
        const auto &coordinates = this->giveNode(i)->giveCoordinates();

        _c = coordinates;
        _c.add(gt3);
        if ( i == 1 ) {
            bb0 = bb1 = _c;
        } else {
            bb0.beMinOf(bb0, _c);
            bb1.beMaxOf(bb1, _c);
        }

        _c = coordinates;
        _c.subtract(gt3);
        bb0.beMinOf(bb0, _c);
        bb1.beMaxOf(bb1, _c);
    }
}

//
// io routines
//
#ifdef __OOFEG

void
TR_SHELL02 :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 3 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( this->giveMaterial()->isActivated(tStep) ) {
        EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
        EASValsSetColor( gc.getElementColor() );
        EASValsSetEdgeColor( gc.getElementEdgeColor() );
        EASValsSetEdgeFlag(true);
        EASValsSetFillStyle(FILL_SOLID);
        EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
        p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
        p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
        p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveCoordinate(3);

        go =  CreateTriangle3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
        EGAttachObject(go, ( EObjectP ) this);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}

void
TR_SHELL02 :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    WCRec p [ 3 ];
    GraphicObj *go;
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( this->giveMaterial()->isActivated(tStep) ) {
        EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
        EASValsSetColor( gc.getDeformedElementColor() );
        EASValsSetEdgeColor( gc.getElementEdgeColor() );
        EASValsSetEdgeFlag(true);
        EASValsSetFillStyle(FILL_SOLID);
        EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale);
        p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
        p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
        p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(3, tStep, defScale);

        go =  CreateTriangle3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}

void
TR_SHELL02 :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    FloatArray v1, v2, v3;
    double s [ 3 ], defScale;
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( !this->giveMaterial()->isActivated(tStep) ) {
        return;
    }

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, gc.giveIntVarType(), gc.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, gc.giveIntVarType(), gc.giveIntVarMode(), 2, tStep);
        result += this->giveInternalStateAtNode(v3, gc.giveIntVarType(), gc.giveIntVarMode(), 3, tStep);
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        double tot_w = 0.;
        FloatArray a, v;
        for ( GaussPoint *gp: *plate->giveDefaultIntegrationRulePtr() ) {
            this->giveIPValue(a, gp, IST_ShellMomentTensor, tStep);
            v.add(gp->giveWeight(), a);
            tot_w += gp->giveWeight();
        }
        v.times(1. / tot_w);
        v1 = v;
        v2 = v;
        v3 = v;
    }

    indx = gc.giveIntVarIndx();

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);
    s [ 2 ] = v3.at(indx);
    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 3; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = gc.getDefScale();
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(3, tStep, defScale);
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
            }
        }
        //     //EASValsSetColor(gc.getYieldPlotColor(ratio));
        gc.updateFringeTableMinMax(s, 3);
        tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

#endif
} // end namespace oofem
