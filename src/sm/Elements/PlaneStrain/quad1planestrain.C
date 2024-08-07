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

#include "sm/Elements/PlaneStrain/quad1planestrain.h"
#include "fei2dquadlin.h"
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
 #include "sm/Materials/rcm2.h"
#endif

namespace oofem {
REGISTER_Element(Quad1PlaneStrain);

FEI2dQuadLin Quad1PlaneStrain :: interp(1, 2);

Quad1PlaneStrain :: Quad1PlaneStrain(int n, Domain *aDomain) :
    PlaneStrainElement(n, aDomain), ZZNodalRecoveryModelInterface(this), SPRNodalRecoveryModelInterface(),
    SpatialLocalizerInterface(this),
    HuertaErrorEstimatorInterface()
{
    numberOfDofMans  = 4;
    numberOfGaussPoints = 4;
}


Quad1PlaneStrain :: ~Quad1PlaneStrain()
{ }


FEInterpolation *Quad1PlaneStrain :: giveInterpolation() const { return & interp; }

const FloatMatrix *
Quad1PlaneStrain::computeGtoLRotationMatrix()
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
		e1.beDifferenceOf(this->giveNode(2)->giveCoordinates(), this->giveNode(1)->giveCoordinates());
		help.beDifferenceOf(this->giveNode(3)->giveCoordinates(), this->giveNode(1)->giveCoordinates());

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


void
Quad1PlaneStrain :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns the [4x8] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,epsilon_z,gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3,u4,v4)
{
    FloatMatrix dN;

    this->interp.evaldNdx( dN, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    // Reshape
    answer.resize(4, 8);
    answer.zero();

#ifdef Quad1PlaneStrain_reducedVolumetricIntegration
    FloatMatrix dN_red;
    FloatArray lcoords_red(2);
    lcoords_red.zero();
    this->interp.evaldNdx( dN_red, lcoords_red, FEIElementGeometryWrapper(this) );

    for ( int i = 1; i <= 4; i++ ) {
      answer.at(1, 2 * i - 1) = ( dN.at(i, 1) + dN_red.at(i, 1) ) / 2.;
      answer.at(1, 2 * i - 0) = (-dN.at(i, 2) + dN_red.at(i, 2) ) / 2.;
      answer.at(2, 2 * i - 1) = (-dN.at(i, 1) + dN_red.at(i, 1) ) / 2.;
      answer.at(2, 2 * i - 0) = ( dN.at(i, 2) + dN_red.at(i, 2) ) / 2.;
      answer.at(4, 2 * i - 1) = dN.at(i, 2);
      answer.at(4, 2 * i - 0) = dN.at(i, 1);
    }
#else
#ifdef Quad1PlaneStrain_reducedShearIntegration
    FloatMatrix dN_red;
    FloatArray lcoords_red(2);
    lcoords_red.zero();
    this->interp.evaldNdx( dN_red, lcoords_red, FEIElementGeometryWrapper(this) );

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dN.at(i, 1);
        answer.at(2, 2 * i - 0) = dN.at(i, 2);
        answer.at(4, 2 * i - 1) = dN_red.at(i, 2);
        answer.at(4, 2 * i - 0) = dN_red.at(i, 1);
    }

#else
    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dN.at(i, 1);
        answer.at(2, 2 * i - 0) = dN.at(i, 2);
        answer.at(4, 2 * i - 1) = dN.at(i, 2);
        answer.at(4, 2 * i - 0) = dN.at(i, 1);
    }

#endif
#endif
}

void
Quad1PlaneStrain :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [5x8] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
{
    FloatMatrix dnx;
    this->interp.evaldNdx( dnx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(5, 8);
    answer.zero();
    // 3rd row is zero -> dw/dz = 0
    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);     // du/dx -1
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);     // dv/dy -2
        answer.at(4, 2 * i - 1) = dnx.at(i, 2);     // du/dy -6
        answer.at(5, 2 * i - 0) = dnx.at(i, 1);     // dv/dx -9
    }
}


void
Quad1PlaneStrain :: initializeFrom(InputRecord &ir)
{
    numberOfGaussPoints = 4;
    PlaneStrainElement :: initializeFrom(ir);

    if ( !( ( numberOfGaussPoints == 4 ) ||
            ( numberOfGaussPoints == 1 ) ||
            ( numberOfGaussPoints == 9 ) ||
            ( numberOfGaussPoints == 16 ) ) ) {
        numberOfGaussPoints = 4;
    }

    // optional record for 1st local axes
    la1.resize(3);
    la1.at(1) = 0; la1.at(2) = 0; la1.at(3) = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->la1, _IFT_Quad1PlaneStrain_FirstLocalAxis);
}



Interface *
Quad1PlaneStrain :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    } else if ( interface == HuertaErrorEstimatorInterfaceType ) {
        return static_cast< HuertaErrorEstimatorInterface * >(this);
    }

    return NULL;
}


void
Quad1PlaneStrain :: HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                                     IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                                     HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                                     int &localNodeId, int &localElemId, int &localBcId,
                                                                     IntArray &controlNode, IntArray &controlDof,
                                                                     HuertaErrorEstimator :: AnalysisMode aMode)
{
    int nodes = 4, sides = 4;
    double x = 0.0, y = 0.0;

    static int sideNode [ 4 ] [ 2 ] = { { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 1 } };

    FloatArray corner [ 4 ], midSide [ 4 ], midNode, cor [ 4 ];
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


void Quad1PlaneStrain :: HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    computeNmatrixAt(gp->giveSubPatchCoordinates(), answer);
}



#ifdef __OOFEG
 #define TR_LENGHT_REDUCT 0.3333

void Quad1PlaneStrain :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 4 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_HOLLOW);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveCoordinate(1);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveCoordinate(2);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Quad1PlaneStrain :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    WCRec p [ 4 ];
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
    EASValsSetFillStyle(FILL_HOLLOW);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}



void Quad1PlaneStrain :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 4 ];
    GraphicObj *tr;
    FloatArray v [ 4 ];
    double s [ 4 ], defScale;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( gc.giveIntVarMode() == ISM_recovered ) {
        for ( i = 1; i <= 4; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], gc.giveIntVarType(), gc.giveIntVarMode(), i, tStep);
        }

        if ( result != 4 ) {
            return;
        }

        indx = gc.giveIntVarIndx();

        for ( i = 1; i <= 4; i++ ) {
            s [ i - 1 ] = v [ i - 1 ].at(indx);
        }

        if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
            for ( i = 0; i < 4; i++ ) {
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
            gc.updateFringeTableMinMax(s, 4);
            tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        } else if ( ( gc.getScalarAlgo() == SA_ZPROFILE ) || ( gc.getScalarAlgo() == SA_COLORZPROFILE ) ) {
            double landScale = gc.getLandScale();

            for ( i = 0; i < 4; i++ ) {
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

                // this fixes a bug in ELIXIR
                if ( fabs(s [ i ]) < 1.0e-6 ) {
                    s [ i ] = 1.0e-6;
                }
            }

            if ( gc.getScalarAlgo() == SA_ZPROFILE ) {
                EASValsSetColor( gc.getDeformedElementColor() );
                EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
                tr =  CreateQuad3D(p);
                EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
            } else {
                gc.updateFringeTableMinMax(s, 4);
                tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
                EGWithMaskChangeAttributes(LAYER_MASK, tr);
            }

            EMAddGraphicsToModel(ESIModel(), tr);
        }
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        if ( numberOfGaussPoints != 4 ) {
            return;
        }

        IntArray ind(4);
        WCRec pp [ 9 ];

        for ( i = 0; i < 4; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = gc.getDefScale();
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                pp [ i ].z = 0.;
            } else {
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                pp [ i ].z = 0.;
            }
        }

        for ( i = 0; i < 3; i++ ) {
            pp [ i + 4 ].x = 0.5 * ( pp [ i ].x + pp [ i + 1 ].x );
            pp [ i + 4 ].y = 0.5 * ( pp [ i ].y + pp [ i + 1 ].y );
            pp [ i + 4 ].z = 0.5 * ( pp [ i ].z + pp [ i + 1 ].z );
        }

        pp [ 7 ].x = 0.5 * ( pp [ 3 ].x + pp [ 0 ].x );
        pp [ 7 ].y = 0.5 * ( pp [ 3 ].y + pp [ 0 ].y );
        pp [ 7 ].z = 0.5 * ( pp [ 3 ].z + pp [ 0 ].z );

        pp [ 8 ].x = 0.25 * ( pp [ 0 ].x + pp [ 1 ].x + pp [ 2 ].x + pp [ 3 ].x );
        pp [ 8 ].y = 0.25 * ( pp [ 0 ].y + pp [ 1 ].y + pp [ 2 ].y + pp [ 3 ].y );
        pp [ 8 ].z = 0.25 * ( pp [ 0 ].z + pp [ 1 ].z + pp [ 2 ].z + pp [ 3 ].z );

        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
            const FloatArray& gpCoords = gp->giveNaturalCoordinates();
            if ( ( gpCoords.at(1) > 0. ) && ( gpCoords.at(2) > 0. ) ) {
                ind.at(1) = 0;
                ind.at(2) = 4;
                ind.at(3) = 8;
                ind.at(4) = 7;
            } else if ( ( gpCoords.at(1) < 0. ) && ( gpCoords.at(2) > 0. ) ) {
                ind.at(1) = 4;
                ind.at(2) = 1;
                ind.at(3) = 5;
                ind.at(4) = 8;
            } else if ( ( gpCoords.at(1) < 0. ) && ( gpCoords.at(2) < 0. ) ) {
                ind.at(1) = 5;
                ind.at(2) = 2;
                ind.at(3) = 6;
                ind.at(4) = 8;
            } else {
                ind.at(1) = 6;
                ind.at(2) = 3;
                ind.at(3) = 7;
                ind.at(4) = 8;
            }

            if ( giveIPValue(v [ 0 ], gp, gc.giveIntVarType(), tStep) == 0 ) {
                return;
            }

            indx = gc.giveIntVarIndx();

            for ( i = 1; i <= 4; i++ ) {
                s [ i - 1 ] = v [ 0 ].at(indx);
            }

            for ( i = 0; i < 4; i++ ) {
                p [ i ].x = pp [ ind.at(i + 1) ].x;
                p [ i ].y = pp [ ind.at(i + 1) ].y;
                p [ i ].z = pp [ ind.at(i + 1) ].z;
            }

            gc.updateFringeTableMinMax(s, 4);
            tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        }
    }
}


void
Quad1PlaneStrain :: drawSpecial(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i;
    WCRec l [ 2 ];
    GraphicObj *tr;
    double defScale = gc.getDefScale();
    FloatArray crackStatuses, cf;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarType() == IST_CrackState ) {
        // ask if any active crack exist
        int crackStatus;
        double ax, ay, bx, by, norm, xc, yc, length;
        FloatArray crackDir;
        FloatArray gpglobalcoords;

        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {

            if ( this->giveIPValue(cf, gp, IST_CrackedFlag, tStep) == 0 ) {
                return;
            }

            if ( ( int ) cf.at(1) == 0 ) {
                return;
            }

            if ( this->giveIPValue(crackDir, gp, IST_CrackDirs, tStep) ) {
                this->giveIPValue(crackStatuses, gp, IST_CrackStatuses, tStep);
                for ( i = 1; i <= 3; i++ ) {
                    crackStatus = ( int ) crackStatuses.at(i);
                    if ( ( crackStatus != pscm_NONE ) && ( crackStatus != pscm_CLOSED ) ) {
                        // draw a crack
                        // this element is 2d element in x-y plane
                        //
                        // compute perpendicular line to normal in xy plane
                        ax = crackDir.at(i);
                        ay = crackDir.at(3 + i);
                        if ( fabs(ax) > 1.e-6 ) {
                            by = 1.;
                            bx = -ay * by / ax;
                            norm = sqrt(bx * bx + by * by);
                            bx = bx / norm; // normalize to obtain unit vector
                            by = by / norm;
                        } else {
                            by = 0.0;
                            bx = 1.0;
                        }

                        // obtain gp global coordinates
                        if ( gc.getInternalVarsDefGeoFlag() ) {
                            double ksi, eta, n1, n2, n3, n4;
                            ksi = gp->giveNaturalCoordinate(1);
                            eta = gp->giveNaturalCoordinate(2);

                            n1 = ( 1. + ksi ) * ( 1. + eta ) * 0.25;
                            n2 = ( 1. - ksi ) * ( 1. + eta ) * 0.25;
                            n3 = ( 1. - ksi ) * ( 1. - eta ) * 0.25;
                            n4 = ( 1. + ksi ) * ( 1. - eta ) * 0.25;

                            gpglobalcoords.resize(2);

                            gpglobalcoords.at(1) = ( n1 * this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale) +
                                                     n2 * this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale) +
                                                     n3 * this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale) +
                                                     n4 * this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale) );
                            gpglobalcoords.at(2) = ( n1 * this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale) +
                                                     n2 * this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale) +
                                                     n3 * this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale) +
                                                     n4 * this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale) );
                        } else {
                            computeGlobalCoordinates( gpglobalcoords, gp->giveNaturalCoordinates() );
                        }

                        xc = gpglobalcoords.at(1);
                        yc = gpglobalcoords.at(2);
                        length = sqrt( computeVolumeAround(gp) / this->giveCrossSection()->give(CS_Thickness, gp) ) / 2.;
                        l [ 0 ].x = ( FPNum ) xc + bx * length;
                        l [ 0 ].y = ( FPNum ) yc + by * length;
                        l [ 0 ].z = 0.;
                        l [ 1 ].x = ( FPNum ) xc - bx * length;
                        l [ 1 ].y = ( FPNum ) yc - by * length;
                        l [ 1 ].z = 0.;
                        EASValsSetLayer(OOFEG_CRACK_PATTERN_LAYER);
                        EASValsSetLineWidth(OOFEG_CRACK_PATTERN_WIDTH);
                        if ( ( crackStatus == pscm_SOFTENING ) || ( crackStatus == pscm_OPEN ) ) {
                            EASValsSetColor( gc.getActiveCrackColor() );
                        } else {
                            EASValsSetColor( gc.getCrackPatternColor() );
                        }

                        tr = CreateLine3D(l);
                        EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
                        EMAddGraphicsToModel(ESIModel(), tr);
                    }
                }
            }
        }
    }
}

#endif


void
Quad1PlaneStrain :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(4);
    for ( int i = 1; i < 5; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}


void
Quad1PlaneStrain :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i < 5; i++ ) {
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


int
Quad1PlaneStrain :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return this->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
}


SPRPatchType
Quad1PlaneStrain :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}

} // end namespace oofem
