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

#include "sm/Elements/PlaneStrain/trplanestrain.h"
#include "../sm/Materials/structuralmaterial.h"
#include "fei2dtrlin.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"
#include "angle.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
#endif

namespace oofem {
REGISTER_Element(TrPlaneStrain);

FEI2dTrLin TrPlaneStrain :: interp(1, 2);

TrPlaneStrain :: TrPlaneStrain(int n, Domain *aDomain) :
    PlaneStrainElement(n, aDomain), ZZNodalRecoveryModelInterface(this), NodalAveragingRecoveryModelInterface(),
    SPRNodalRecoveryModelInterface(), SpatialLocalizerInterface(this),
    ZZErrorEstimatorInterface(this),
    HuertaErrorEstimatorInterface()
    // Constructor.
{
    numberOfDofMans  = 3;
    area = -1;
    numberOfGaussPoints = 1;
}


FEInterpolation *TrPlaneStrain :: giveInterpolation() const { return & interp; }


Interface *
TrPlaneStrain :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    } else if ( interface == ZZErrorEstimatorInterfaceType ) {
        return static_cast< ZZErrorEstimatorInterface * >(this);
    } else if ( interface == HuertaErrorEstimatorInterfaceType ) {
        return static_cast< HuertaErrorEstimatorInterface * >(this);
    }

    return NULL;
}

void
TrPlaneStrain :: initializeFrom(InputRecord &ir)
{
    numberOfGaussPoints = 1;
    PlaneStrainElement :: initializeFrom(ir);

	// optional record for 1st local axes
	la1.resize(3);
	la1.at(1) = 0; la1.at(2) = 0; la1.at(3) = 0;
	IR_GIVE_OPTIONAL_FIELD(ir, this->la1, _IFT_TrPlaneStrain_FirstLocalAxis);

    if ( numberOfGaussPoints != 1 ) {
        numberOfGaussPoints = 1;
    }
}

const FloatMatrix *
TrPlaneStrain::computeGtoLRotationMatrix()
// Returns the rotation matrix of the receiver of the size [3,3]
// coords(local) = T * coords(global)
//
// local coordinate (described by vector triplet e1',e2',e3') is defined as follows:
//
// e1'    : [N2-N1]    Ni - means i - th node
// help   : [N3-N1]
// e3'    : e1' x help
// e2'    : e3' x e1'
{
	if (!GtoLRotationMatrix.isNotEmpty()) {
		FloatArray e1, e2, e3, help;

		// compute e1' = [N2-N1]  and  help = [N3-N1]
		e1.beDifferenceOf(*this->giveNode(2)->giveCoordinates(), *this->giveNode(1)->giveCoordinates());
		help.beDifferenceOf(*this->giveNode(3)->giveCoordinates(), *this->giveNode(1)->giveCoordinates());

		// let us normalize e1'
		e1.normalize();

		// compute e3' : vector product of e1' x help
		e3.beVectorProductOf(e1, help);
		// let us normalize
		e3.normalize();

		//if (la1.computeNorm() != 0) {
		//	// custom local axes
		//	e1 = la1;
		//}

		// now from e3' x e1' compute e2'
		e2.beVectorProductOf(e3, e1);

		// rotate as to have the 1st local axis equal to la1
		if (la1.computeNorm() != 0) {
			double ang = -Angle::giveAngleIn3Dplane(la1, e1, e3); // radians
			e1 = Angle::rotate(e1, e3, ang);
			e2 = Angle::rotate(e2, e3, ang);
		}
		// rot. matrix
		GtoLRotationMatrix.resize(3, 3);

		for (int i = 1; i <= 3; i++) {
			GtoLRotationMatrix.at(1, i) = e1.at(i);
			GtoLRotationMatrix.at(2, i) = e2.at(i);
			GtoLRotationMatrix.at(3, i) = e3.at(i);
		}
	}

	return &GtoLRotationMatrix;
}


// #define POINT_TOL 1.e-3
//int
//TrPlaneStrain :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
//{
//	double detJ, x1, x2, x3, y1, y2, y3;
//
//	x1 = cellgeo.giveVertexCoordinates(1)->at(xind);
//	x2 = cellgeo.giveVertexCoordinates(2)->at(xind);
//	x3 = cellgeo.giveVertexCoordinates(3)->at(xind);
//
//	y1 = cellgeo.giveVertexCoordinates(1)->at(yind);
//	y2 = cellgeo.giveVertexCoordinates(2)->at(yind);
//	y3 = cellgeo.giveVertexCoordinates(3)->at(yind);
//
//	detJ = x1 * (y2 - y3) + x2 * (-y1 + y3) + x3 * (y1 - y2);
//
//	answer.resize(3);
//	answer.at(1) = ((x2 * y3 - x3 * y2) + (y2 - y3) * coords.at(xind) + (x3 - x2) * coords.at(yind)) / detJ;
//	answer.at(2) = ((x3 * y1 - x1 * y3) + (y3 - y1) * coords.at(xind) + (x1 - x3) * coords.at(yind)) / detJ;
//
//	// check if point is inside
//	bool inside = true;
//	for (int i = 1; i <= 2; i++) {
//		if (answer.at(i) < (0. - POINT_TOL)) {
//			answer.at(i) = 0.;
//			inside = false;
//		}
//		else if (answer.at(i) > (1. + POINT_TOL)) {
//			answer.at(i) = 1.;
//			inside = false;
//		}
//	}
//
//	if ((answer.at(1) + answer.at(2)) > 1.0) {
//		const double temp = 0.5*(answer.at(1) + answer.at(2) - 1.);
//		answer.at(1) -= temp;
//		answer.at(2) -= temp;
//		inside = false;
//	}
//
//	answer.at(3) = 1. - answer.at(1) - answer.at(2);
//
//	return inside;
//}

void
TrPlaneStrain :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep)
{
    //GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    //this->giveIPValue(answer, gp, type, tStep);
	answer.zero();
	GaussPoint *gp = integrationRulesArray[0]->getIntegrationPoint(0);
	FloatArray localStress, localStrain;
	this->computeStrainVector(localStrain, gp, tStep);
	this->computeStressVector(localStress, localStrain, gp, tStep);

	if (type == IST_ShellForceTensor) {
		double thickness, w;
		FloatArray mLocal;
		mLocal.resize(6);
		mLocal.zero();

		thickness = this->giveCrossSection()->give(CS_Thickness, gp->giveGlobalCoordinates(), this, false);
		w = gp->giveWeight() * thickness;
		// mLocal.add(w, localStress);
		mLocal.at(1) = localStress.at(1) * w;
		mLocal.at(2) = localStress.at(2) * w;
		mLocal.at(6) = localStress.at(3) * w;

		// local to global
		this->computeGtoLRotationMatrix();
		StructuralMaterial::transformStressVectorTo(answer, GtoLRotationMatrix, mLocal, false);
	}
	else if (type == IST_ShellMomentTensor || type == IST_CurvatureTensor) {
		// nothing
	}
	else if (type == IST_ShellStrainTensor) {
		FloatArray mLocal;
		mLocal.resize(6);
		mLocal.zero();
		mLocal.at(1) = localStrain.at(1);
		mLocal.at(2) = localStrain.at(2);
		mLocal.at(6) = localStrain.at(3);

		// local to global
		this->computeGtoLRotationMatrix();
		StructuralMaterial::transformStressVectorTo(answer, GtoLRotationMatrix, mLocal, false);
	}
}


void
TrPlaneStrain :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}


void
TrPlaneStrain :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    answer.resize(1);
    if ( ( pap == this->giveNode(1)->giveNumber() ) ||
        ( pap == this->giveNode(2)->giveNumber() ) ||
        ( pap == this->giveNode(3)->giveNumber() ) ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("node unknown");
    }
}


int
TrPlaneStrain :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return 1;
}


SPRPatchType
TrPlaneStrain :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}


void
TrPlaneStrain :: HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                                  IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                                  HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                                  int &localNodeId, int &localElemId, int &localBcId,
                                                                  IntArray &controlNode, IntArray &controlDof,
                                                                  HuertaErrorEstimator :: AnalysisMode aMode)
{
    int nodes = 3, sides = 3;
    FloatArray corner [ 3 ], midSide [ 3 ], midNode, cor [ 3 ];
    double x = 0.0, y = 0.0;

    static int sideNode [ 3 ] [ 2 ] = { { 1, 2 }, { 2, 3 }, { 3, 1 } };

    if ( sMode == HuertaErrorEstimatorInterface :: NodeMode ||
        ( sMode == HuertaErrorEstimatorInterface :: BCMode && aMode == HuertaErrorEstimator :: HEE_linear ) ) {
        for ( int inode = 0; inode < nodes; inode++ ) {
            corner [ inode ] = this->giveNode(inode + 1)->giveCoordinates();
            if ( corner [ inode ].giveSize() != 3 ) {
                cor [ inode ].resize(3);
                cor [ inode ].at(1) = corner [ inode ].at(1);
                cor [ inode ].at(2) = corner [ inode ].at(2);
                cor [ inode ].at(3) = 0.0;

                corner [ inode ] = cor [ inode ];
            }

            x += corner [ inode ].at(1);
            y += corner [ inode ].at(2);
        }

        for ( int iside = 0; iside < sides; iside++ ) {
            midSide [ iside ].resize(3);

            int nd1 = sideNode [ iside ] [ 0 ] - 1;
            int nd2 = sideNode [ iside ] [ 1 ] - 1;

            midSide [ iside ].at(1) = ( corner [ nd1 ].at(1) + corner [ nd2 ].at(1) ) / 2.0;
            midSide [ iside ].at(2) = ( corner [ nd1 ].at(2) + corner [ nd2 ].at(2) ) / 2.0;
            midSide [ iside ].at(3) = 0.0;
        }

        midNode.resize(3);

        midNode.at(1) = x / nodes;
        midNode.at(2) = y / nodes;
        midNode.at(3) = 0.0;
    }

    this->setupRefinedElementProblem2D(this, refinedElement, level, nodeId, localNodeIdArray, globalNodeIdArray,
                                       sMode, tStep, nodes, corner, midSide, midNode,
                                       localNodeId, localElemId, localBcId,
                                       controlNode, controlDof, aMode, "Quad1PlaneStrain");
}


void TrPlaneStrain :: HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    computeNmatrixAt(gp->giveSubPatchCoordinates(), answer);
}



#ifdef __OOFEG
 #define TR_LENGHT_REDUCT 0.3333

void TrPlaneStrain :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 3 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = 0.;

    go =  CreateTriangle3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void TrPlaneStrain :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    WCRec p [ 3 ];
    GraphicObj *go;
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 2 ].z = 0.;

    go =  CreateTriangle3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void TrPlaneStrain :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    FloatArray v1, v2, v3;
    double s [ 3 ], defScale;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, gc.giveIntVarType(), gc.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, gc.giveIntVarType(), gc.giveIntVarMode(), 2, tStep);
        result += this->giveInternalStateAtNode(v3, gc.giveIntVarType(), gc.giveIntVarMode(), 3, tStep);
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        result += giveIPValue(v1, gp, gc.giveIntVarType(), tStep);
        v2 = v1;
        v3 = v1;
        result *= 3;
    }

    if ( result != 3 ) {
        return;
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
                p [ i ].z = 0.;
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                p [ i ].z = 0.;
            }
        }

        //EASValsSetColor(gc.getYieldPlotColor(ratio));
        gc.updateFringeTableMinMax(s, 3);
        tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    } else if ( ( gc.getScalarAlgo() == SA_ZPROFILE ) || ( gc.getScalarAlgo() == SA_COLORZPROFILE ) ) {
        double landScale = gc.getLandScale();

        for ( i = 0; i < 3; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = gc.getDefScale();
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                p [ i ].z = s [ i ] * landScale;
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                p [ i ].z = s [ i ] * landScale;
            }
        }

        if ( gc.getScalarAlgo() == SA_ZPROFILE ) {
            EASValsSetColor( gc.getDeformedElementColor() );
            EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangle3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | FILL_MASK | LAYER_MASK, tr);
        } else {
            gc.updateFringeTableMinMax(s, 3);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
            EGWithMaskChangeAttributes(FILL_MASK | LAYER_MASK, tr);
        }

        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

#endif
} // end namespace oofem
