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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#include "sm/Elements/tria1platesubsoil.h"
#include "sm/Materials/structuralms.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "fei2dtrlin.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "load.h"
#include "mathfem.h"
#include "classfactory.h"
#include "angle.h"

namespace oofem {
REGISTER_Element(Tria1PlateSubSoil);

FEI2dTrLin Tria1PlateSubSoil :: interp_lin(1, 2);

Tria1PlateSubSoil :: Tria1PlateSubSoil(int n, Domain *aDomain) :
    StructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this),
    SPRNodalRecoveryModelInterface()  // , NodalAveragingRecoveryModelInterface()
{
    numberOfGaussPoints = 1;
    numberOfDofMans = 3;
}


FEInterpolation *
Tria1PlateSubSoil :: giveInterpolation(DofIDItem id) const
{
    return & interp_lin;
}


FEInterpolation *
Tria1PlateSubSoil :: giveInterpolation() const { return & interp_lin; }


void
Tria1PlateSubSoil :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 5);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
Tria1PlateSubSoil :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
{
  OOFEM_ERROR("Body load not supported, use surface load instead");
}


void
Tria1PlateSubSoil :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [3x3] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    FloatArray n;
    FloatMatrix dn;

    this->interp_lin.evaldNdx( dn, gp->giveNaturalCoordinates(),  FEIElementGeometryWrapper(this) );
    this->interp_lin.evalN( n, gp->giveNaturalCoordinates(),  FEIElementGeometryWrapper(this) );

    answer.resize(3, 3);
    answer.zero();

    ///@todo Check sign here
    for ( int i = 0; i < 3; ++i ) {
      answer(0, i) = n(i); // eps_z
      answer(1, i) = dn(i, 0); // gamma_xz
      answer(2, i) = dn(i, 1); // gamma_yz
    }
}


void
Tria1PlateSubSoil :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->giveGeneralizedStress_PlateSubSoil(strain, gp, tStep);
}


void
Tria1PlateSubSoil :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->give2dPlateSubSoilStiffMtrx(rMode, gp, tStep);
}


void
Tria1PlateSubSoil :: initializeFrom(InputRecord &ir)
{
    this->numberOfGaussPoints = 1;
    StructuralElement::initializeFrom(ir);
    // optional record for 1st local axes - here it is not used, unuseful
    la1.resize(3);
    la1.at(1) = 0; la1.at(2) = 0; la1.at(3) = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->la1, _IFT_Tria1PlateSubSoil_lcs);

    this->macroElem = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->macroElem, _IFT_Tria1PlateSubSoil_macroelem);
}


void
Tria1PlateSubSoil :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_w};
}


void
Tria1PlateSubSoil :: computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp)
{
    FloatArray u, v;
    u.beDifferenceOf( this->giveNode(2)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );
    v.beDifferenceOf( this->giveNode(3)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );

    answer.beVectorProductOf(u, v);
    answer.normalize();
}


double
Tria1PlateSubSoil :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
//
// returns receiver's characteristic length for crack band models
// for a crack formed in the plane with normal normalToCrackPlane.
//
{
    return this->giveCharacteristicLengthForPlaneElements(normalToCrackPlane);
}


double
Tria1PlateSubSoil :: computeVolumeAround(GaussPoint *gp)
{
    double detJ, weight;

    weight = gp->giveWeight();
    detJ = fabs( this->interp_lin.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    return detJ * weight;
}


void
Tria1PlateSubSoil :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
  // OOFEM_ERROR("Mass matrix not provided");
}

void 
Tria1PlateSubSoil ::computeGtoLmatrix( FloatMatrix &lc )
{
    FloatArray e1, e2, e3, help;
    // compute e1' = [N2-N1]  and  help = [N3-N1]
    e1.beDifferenceOf( this->giveNode( 2 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );
    help.beDifferenceOf( this->giveNode( 3 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );

    // let us normalize e1'
    e1.normalize();

    // compute e3' : vector product of e1' x help
    e3.beVectorProductOf( e1, help );
    // let us normalize
    e3.normalize();

    // now from e3' x e1' compute e2'
    e2.beVectorProductOf( e3, e1 );

    // rotate as to have the 1st local axis equal to la1
    if ( la1.computeNorm() != 0 ) {
        double ang = -Angle::giveAngleIn3Dplane( la1, e1, e3 ); // radians
        e1         = Angle::rotate( e1, e3, ang );
        e2         = Angle::rotate( e2, e3, ang );
    }
    // rot. matrix
    lc.resize( 3, 3 );
    // 1st axis is local z
    for ( int i = 1; i <= 3; i++ ) {
        lc.at( 3, i ) = e1.at( i );
        lc.at( 2, i ) = e2.at( i );
        lc.at( 1, i ) = e3.at( i );
    }
}

FloatMatrixF<3, 3>
Tria1PlateSubSoil ::P3SSMI_getUnknownsGtoLRotationMatrix() const
// Returns the rotation matrix for element unknowns
{
    FloatMatrixF<3, 3> answer;
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at( i, j ) = this->lcs.at( i, j );
        }
    }
    return answer;
}

int
Tria1PlateSubSoil :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
  return StructuralElement :: giveIPValue(answer, gp, type, tStep);
}

void
Tria1PlateSubSoil::NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
InternalStateType type, TimeStep *tStep)
{
	this->giveIPValue(answer, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0), type, tStep);
}

Interface *
Tria1PlateSubSoil :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if (interface == NodalAveragingRecoveryModelInterfaceType) {
    //	return static_cast< NodalAveragingRecoveryModelInterface * >(this);
        return NULL; // avoid results
    } else if ( interface == Plate3dSubsoilMaterialInterfaceType ) {
        if ( this->lcs.giveNumberOfRows() == 0 ) this->computeGtoLmatrix( this->lcs );
        return static_cast<Plate3dSubsoilMaterialInterface *>( this );
    }

    return NULL;
}

void
Tria1PlateSubSoil :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    for ( int i = 1; i < 4; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
Tria1PlateSubSoil :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i < 4; i++ ) {
        if ( pap == this->giveNode(i)->giveNumber() ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("node unknown");
    }
}

void 
Tria1PlateSubSoil::computeInitialStressMatrix( FloatMatrix &answer, TimeStep *tStep )
{
    answer.resize( 3, 3 );
    answer.zero();
}

void
Tria1PlateSubSoil ::computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the [1x3] displacement interpolation matrix {N}
{
    FloatArray N(3);
    giveInterpolation()->evalN(N, iLocCoord, FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(N, 1);
}


void
Tria1PlateSubSoil :: computeSurfaceNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords)
{
    if (boundaryID == 1) {
        this->computeNmatrixAt(lcoords, answer);
    } else {
        OOFEM_ERROR("computeSurfaceNMatrix: Only one surface is supported with id=1");
    }
}


void
Tria1PlateSubSoil :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
	fprintf(file, "element %d (%8d) macroelem %d :\n", this->giveLabel(), number, this->macroElem);
}

} // end namespace oofem
