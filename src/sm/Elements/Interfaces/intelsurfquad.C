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

#include "../sm/Elements/Interfaces/IntElSurfQuad.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"

 #include <Emarkwd3d.h>
#endif

namespace oofem {
REGISTER_Element(IntElSurfQuad);

FEI3dHexaLin IntElSurfQuad :: interpolation;

IntElSurfQuad :: IntElSurfQuad(int n, Domain *aDomain) :
    StructuralInterfaceElement(n, aDomain)
{
    numberOfDofMans = 8;
}


void
IntElSurfQuad :: computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer)
{
    // Returns the modified N-matrix which multiplied with u give the spatial jump.

    FloatArray n;
    this->interpolation.evalN( n, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(3, 24);
    answer.zero();

	answer.at(1, 13) = answer.at(2, 14) = answer.at(3, 15) = n.at(1);
	answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = -n.at(1);

	answer.at(1, 16) = answer.at(2, 17) = answer.at(3, 18) = n.at(2);
	answer.at(1, 4) = answer.at(2, 5) = answer.at(3, 6) = -n.at(2);

	answer.at(1, 19) = answer.at(2, 20) = answer.at(3, 21) = n.at(3);
	answer.at(1, 7) = answer.at(2, 8) = answer.at(3, 9) = -n.at(3);

	answer.at(1, 22) = answer.at(2, 23) = answer.at(3, 24) = n.at(4);
	answer.at(1, 10) = answer.at(2, 11) = answer.at(3, 12) = -n.at(4);
}


void
IntElSurfQuad :: computeGaussPoints()
{
    // Sets up the array of Gauss Points of the receiver.
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        //integrationRulesArray[0] = new LobattoIntegrationRule (1,domain, 1, 2);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) ); 
		integrationRulesArray[0]->setUpIntegrationPoints(_Square, 4, _3dInterface); ///@todo add parameter for ngp
    }
}


void
IntElSurfQuad :: computeCovarBaseVectorsAt(IntegrationPoint *ip, FloatArray &G1, FloatArray &G2)
{
    FloatMatrix dNdxi;
    this->interpolation.evaldNdxi( dNdxi, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    G1.resize(4);
    G2.resize(4);
    G1.zero();
    G2.zero();

    FloatArray meanNode;
    int numNodes = this->giveNumberOfNodes();
    for ( int i = 1; i <= dNdxi.giveNumberOfRows(); i++ ) {
        meanNode = 0.5 * ( *this->giveNode(i)->giveCoordinates() + *this->giveNode(i + numNodes / 2)->giveCoordinates() );
        G1 += dNdxi.at(i, 1) * meanNode;
        G2 += dNdxi.at(i, 2) * meanNode;
    }
}


double
IntElSurfQuad :: computeAreaAround(IntegrationPoint *ip)
{
    FloatArray G1, G2, G3;
    this->computeCovarBaseVectorsAt(ip, G1, G2);
    double weight  = ip->giveWeight();
    G3.beVectorProductOf(G1, G2);
    return G3.computeNorm() * weight; // 0.5 *
}


IRResultType
IntElSurfQuad :: initializeFrom(InputRecord *ir)
{
    return StructuralInterfaceElement :: initializeFrom(ir);
}


void
IntElSurfQuad :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = { D_u, D_v, D_w };
}

void
IntElSurfQuad :: computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // Transformation matrix to the local coordinate system
    FloatArray G1, G2, Normal;
    this->computeCovarBaseVectorsAt(gp, G1, G2);
    Normal.beVectorProductOf(G1, G2);
    Normal.normalize();
    answer.beLocalCoordSys(Normal);
    
}



int
IntElSurfQuad :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    FloatArray N, meanNode;
    this->interpolation.evalN( N, lcoords, FEIElementGeometryWrapper(this) );
    answer.resize(3);
    answer.zero();
    for ( int i = 1; i <= 3; i++ ) {
        meanNode = 0.5 * ( *this->giveNode(i)->giveCoordinates() + *this->giveNode(i + 3)->giveCoordinates() );
        answer += N.at(i) * meanNode;
    }

    return 1;
}


bool
IntElSurfQuad :: computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords)
{
    OOFEM_ERROR("Not implemented");
    return false;
}

///@todo this code not tested, onlu copied from interfaceelem3dtrlin.C //JB
#ifdef __OOFEG
void IntElSurfQuad :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p [ 3 ]; /* triangle */
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
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveCoordinate(3);

    go =  CreateTriangle3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void IntElSurfQuad :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{ }


void IntElSurfQuad :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{ }

#endif
} // end namespace oofem
