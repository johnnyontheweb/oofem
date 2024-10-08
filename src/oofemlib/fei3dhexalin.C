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

#include "fei3dhexalin.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"

namespace oofem {

FloatArrayF<8>
FEI3dHexaLin :: evalN(const FloatArrayF<3> &lcoords)
{
    ///auto [x, y, z] = lcoords; // structured bindings C++17 would be nice.
    double x = lcoords[0];
    double y = lcoords[1];
    double z = lcoords[2];

    return {
        0.125 * ( 1. - x ) * ( 1. - y ) * ( 1. + z ),
        0.125 * ( 1. - x ) * ( 1. + y ) * ( 1. + z ),
        0.125 * ( 1. + x ) * ( 1. + y ) * ( 1. + z ),
        0.125 * ( 1. + x ) * ( 1. - y ) * ( 1. + z ),
        0.125 * ( 1. - x ) * ( 1. - y ) * ( 1. - z ),
        0.125 * ( 1. - x ) * ( 1. + y ) * ( 1. - z ),
        0.125 * ( 1. + x ) * ( 1. + y ) * ( 1. - z ),
        0.125 * ( 1. + x ) * ( 1. - y ) * ( 1. - z )
    };
}


void
FEI3dHexaLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
#if 0
    answer = evalN(lcoords);
#else
    double x = lcoords[0];
    double y = lcoords[1];
    double z = lcoords[2];

    answer = {
        0.125 * ( 1. - x ) * ( 1. - y ) * ( 1. + z ),
        0.125 * ( 1. - x ) * ( 1. + y ) * ( 1. + z ),
        0.125 * ( 1. + x ) * ( 1. + y ) * ( 1. + z ),
        0.125 * ( 1. + x ) * ( 1. - y ) * ( 1. + z ),
        0.125 * ( 1. - x ) * ( 1. - y ) * ( 1. - z ),
        0.125 * ( 1. - x ) * ( 1. + y ) * ( 1. - z ),
        0.125 * ( 1. + x ) * ( 1. + y ) * ( 1. - z ),
        0.125 * ( 1. + x ) * ( 1. - y ) * ( 1. - z )
    };
#endif
}


FloatMatrixF<3,8>
FEI3dHexaLin :: evaldNdxi(const FloatArrayF<3> &lcoords) 
{
    //auto [u, v, w] = lcoords;
    double u = lcoords.at(1);
    double v = lcoords.at(2);
    double w = lcoords.at(3);

    return {
        -0.125 * ( 1. - v ) * ( 1. + w ),
        -0.125 * ( 1. - u ) * ( 1. + w ),
        0.125 * ( 1. - u ) * ( 1. - v ),
        -0.125 * ( 1. + v ) * ( 1. + w ),
        0.125 * ( 1. - u ) * ( 1. + w ),
        0.125 * ( 1. - u ) * ( 1. + v ),
        0.125 * ( 1. + v ) * ( 1. + w ),
        0.125 * ( 1. + u ) * ( 1. + w ),
        0.125 * ( 1. + u ) * ( 1. + v ),
        0.125 * ( 1. - v ) * ( 1. + w ),
        -0.125 * ( 1. + u ) * ( 1. + w ),
        0.125 * ( 1. + u ) * ( 1. - v ),
        -0.125 * ( 1. - v ) * ( 1. - w ),
        -0.125 * ( 1. - u ) * ( 1. - w ),
        -0.125 * ( 1. - u ) * ( 1. - v ),
        -0.125 * ( 1. + v ) * ( 1. - w ),
        0.125 * ( 1. - u ) * ( 1. - w ),
        -0.125 * ( 1. - u ) * ( 1. + v ),
        0.125 * ( 1. + v ) * ( 1. - w ),
        0.125 * ( 1. + u ) * ( 1. - w ),
        -0.125 * ( 1. + u ) * ( 1. + v ),
        0.125 * ( 1. - v ) * ( 1. - w ),
        -0.125 * ( 1. + u ) * ( 1. - w ),
        -0.125 * ( 1. + u ) * ( 1. - v ),
    };
}

void
FEI3dHexaLin :: evaldNdxi(FloatMatrix &dN, const FloatArray &lcoords, const FEICellGeometry &) const
{
#if 0
    dN = evaldNdxi(lcoords);
#else
    double u, v, w;
    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);

    dN.resize(8, 3);

    dN.at(1, 1) = -0.125 * ( 1. - v ) * ( 1. + w );
    dN.at(2, 1) = -0.125 * ( 1. + v ) * ( 1. + w );
    dN.at(3, 1) =  0.125 * ( 1. + v ) * ( 1. + w );
    dN.at(4, 1) =  0.125 * ( 1. - v ) * ( 1. + w );
    dN.at(5, 1) = -0.125 * ( 1. - v ) * ( 1. - w );
    dN.at(6, 1) = -0.125 * ( 1. + v ) * ( 1. - w );
    dN.at(7, 1) =  0.125 * ( 1. + v ) * ( 1. - w );
    dN.at(8, 1) =  0.125 * ( 1. - v ) * ( 1. - w );

    dN.at(1, 2) = -0.125 * ( 1. - u ) * ( 1. + w );
    dN.at(2, 2) =  0.125 * ( 1. - u ) * ( 1. + w );
    dN.at(3, 2) =  0.125 * ( 1. + u ) * ( 1. + w );
    dN.at(4, 2) = -0.125 * ( 1. + u ) * ( 1. + w );
    dN.at(5, 2) = -0.125 * ( 1. - u ) * ( 1. - w );
    dN.at(6, 2) =  0.125 * ( 1. - u ) * ( 1. - w );
    dN.at(7, 2) =  0.125 * ( 1. + u ) * ( 1. - w );
    dN.at(8, 2) = -0.125 * ( 1. + u ) * ( 1. - w );

    dN.at(1, 3) =  0.125 * ( 1. - u ) * ( 1. - v );
    dN.at(2, 3) =  0.125 * ( 1. - u ) * ( 1. + v );
    dN.at(3, 3) =  0.125 * ( 1. + u ) * ( 1. + v );
    dN.at(4, 3) =  0.125 * ( 1. + u ) * ( 1. - v );
    dN.at(5, 3) = -0.125 * ( 1. - u ) * ( 1. - v );
    dN.at(6, 3) = -0.125 * ( 1. - u ) * ( 1. + v );
    dN.at(7, 3) = -0.125 * ( 1. + u ) * ( 1. + v );
    dN.at(8, 3) = -0.125 * ( 1. + u ) * ( 1. - v );
#endif
}


std::pair<double, FloatMatrixF<3,8>>
FEI3dHexaLin :: evaldNdx(const FloatArrayF<3> &lcoords, const FEICellGeometry &cellgeo) 
{
    auto dNduvw = evaldNdxi(lcoords);
    FloatMatrixF<3,8> coords;
    for ( int i = 0; i < 8; i++ ) {
        ///@todo cellgeo should give a FloatArrayF<3>, this will add a "costly" construction now:
        coords.setColumn(cellgeo.giveVertexCoordinates(i+1), i);
    }
    auto jacT = dotT(dNduvw, coords);
    return {det(jacT), dot(inv(jacT), dNduvw)};
}


double
FEI3dHexaLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
#if 0
    auto tmp = evaldNdx(lcoords, cellgeo);
    //auto [det, dndx] = evaldNdx(lcoords, cellgeo);
    answer = tmp.second;
    return tmp.first;
#else
    FloatMatrix jacobianMatrix, inv, dNduvw, coords;

    this->evaldNdxi(dNduvw, lcoords, cellgeo);
    coords.resize(3, 8);
    for ( int i = 1; i <= 8; i++ ) {
        coords.setColumn(cellgeo.giveVertexCoordinates(i), i);
    }
    jacobianMatrix.beProductOf(coords, dNduvw);
    inv.beInverseOf(jacobianMatrix);

    answer.beProductOf(dNduvw, inv);
    return jacobianMatrix.giveDeterminant();
#endif
}

void
FEI3dHexaLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    FloatArray n;
    this->evalN(n, lcoords, cellgeo);

    answer.clear();
    for ( int i = 1; i <= 8; i++ ) {
        answer.add( n.at(i), cellgeo.giveVertexCoordinates(i) );
    }
}

#define POINT_TOL 1.e-3

int
FEI3dHexaLin :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo) const
{
    double x1, x2, x3, x4, x5, x6, x7, x8, a1, a2, a3, a4, a5, a6, a7, a8;
    double y1, y2, y3, y4, y5, y6, y7, y8, b1, b2, b3, b4, b5, b6, b7, b8;
    double z1, z2, z3, z4, z5, z6, z7, z8, c1, c2, c3, c4, c5, c6, c7, c8;
    double xp, yp, zp, u, v, w;
    FloatMatrix p(3, 3);
    FloatArray r(3), delta;
    int nite = 0;

    x1 = cellgeo.giveVertexCoordinates(1).at(1);
    x2 = cellgeo.giveVertexCoordinates(2).at(1);
    x3 = cellgeo.giveVertexCoordinates(3).at(1);
    x4 = cellgeo.giveVertexCoordinates(4).at(1);
    x5 = cellgeo.giveVertexCoordinates(5).at(1);
    x6 = cellgeo.giveVertexCoordinates(6).at(1);
    x7 = cellgeo.giveVertexCoordinates(7).at(1);
    x8 = cellgeo.giveVertexCoordinates(8).at(1);

    y1 = cellgeo.giveVertexCoordinates(1).at(2);
    y2 = cellgeo.giveVertexCoordinates(2).at(2);
    y3 = cellgeo.giveVertexCoordinates(3).at(2);
    y4 = cellgeo.giveVertexCoordinates(4).at(2);
    y5 = cellgeo.giveVertexCoordinates(5).at(2);
    y6 = cellgeo.giveVertexCoordinates(6).at(2);
    y7 = cellgeo.giveVertexCoordinates(7).at(2);
    y8 = cellgeo.giveVertexCoordinates(8).at(2);

    z1 = cellgeo.giveVertexCoordinates(1).at(3);
    z2 = cellgeo.giveVertexCoordinates(2).at(3);
    z3 = cellgeo.giveVertexCoordinates(3).at(3);
    z4 = cellgeo.giveVertexCoordinates(4).at(3);
    z5 = cellgeo.giveVertexCoordinates(5).at(3);
    z6 = cellgeo.giveVertexCoordinates(6).at(3);
    z7 = cellgeo.giveVertexCoordinates(7).at(3);
    z8 = cellgeo.giveVertexCoordinates(8).at(3);

    xp = coords.at(1);
    yp = coords.at(2);
    zp = coords.at(3);

    a1 =  x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8;
    a2 = -x1 - x2 + x3 + x4 - x5 - x6 + x7 + x8;
    a3 = -x1 + x2 + x3 - x4 - x5 + x6 + x7 - x8;
    a4 =  x1 + x2 + x3 + x4 - x5 - x6 - x7 - x8;
    a5 =  x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8;
    a6 = -x1 - x2 + x3 + x4 + x5 + x6 - x7 - x8;
    a7 = -x1 + x2 + x3 - x4 + x5 - x6 - x7 + x8;
    a8 =  x1 - x2 + x3 - x4 - x5 + x6 - x7 + x8;

    b1 =  y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8;
    b2 = -y1 - y2 + y3 + y4 - y5 - y6 + y7 + y8;
    b3 = -y1 + y2 + y3 - y4 - y5 + y6 + y7 - y8;
    b4 =  y1 + y2 + y3 + y4 - y5 - y6 - y7 - y8;
    b5 =  y1 - y2 + y3 - y4 + y5 - y6 + y7 - y8;
    b6 = -y1 - y2 + y3 + y4 + y5 + y6 - y7 - y8;
    b7 = -y1 + y2 + y3 - y4 + y5 - y6 - y7 + y8;
    b8 =  y1 - y2 + y3 - y4 - y5 + y6 - y7 + y8;

    c1 =  z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8;
    c2 = -z1 - z2 + z3 + z4 - z5 - z6 + z7 + z8;
    c3 = -z1 + z2 + z3 - z4 - z5 + z6 + z7 - z8;
    c4 =  z1 + z2 + z3 + z4 - z5 - z6 - z7 - z8;
    c5 =  z1 - z2 + z3 - z4 + z5 - z6 + z7 - z8;
    c6 = -z1 - z2 + z3 + z4 + z5 + z6 - z7 - z8;
    c7 = -z1 + z2 + z3 - z4 + z5 - z6 - z7 + z8;
    c8 =  z1 - z2 + z3 - z4 - z5 + z6 - z7 + z8;

    // setup initial guess
    answer.resize(3);
    answer.zero();

    // apply Newton-Raphson to solve the problem
    for ( ;; ) {
        if ( ( ++nite ) > 10 ) {
            answer.zero();
            return false;
        }

        u = answer.at(1);
        v = answer.at(2);
        w = answer.at(3);

        // compute the residual
        r.at(1) = a1 + u * a2 + v * a3 + w * a4 + u * v * a5 + u * w * a6 + v * w * a7 + u * v * w * a8 - 8.0 * xp;
        r.at(2) = b1 + u * b2 + v * b3 + w * b4 + u * v * b5 + u * w * b6 + v * w * b7 + u * v * w * b8 - 8.0 * yp;
        r.at(3) = c1 + u * c2 + v * c3 + w * c4 + u * v * c5 + u * w * c6 + v * w * c7 + u * v * w * c8 - 8.0 * zp;

        // check for convergence
        if ( r.computeSquaredNorm() < 1.e-20 ) {
            break;                                  // sqrt(1.e-20) = 1.e-10
        }

        p.at(1, 1) = a2 + v * a5 + w * a6 + v * w * a8;
        p.at(1, 2) = a3 + u * a5 + w * a7 + u * w * a8;
        p.at(1, 3) = a4 + u * a6 + v * a7 + u * v * a8;

        p.at(2, 1) = b2 + v * b5 + w * b6 + v * w * b8;
        p.at(2, 2) = b3 + u * b5 + w * b7 + u * w * b8;
        p.at(2, 3) = b4 + u * b6 + v * b7 + u * v * b8;

        p.at(3, 1) = c2 + v * c5 + w * c6 + v * w * c8;
        p.at(3, 2) = c3 + u * c5 + w * c7 + u * w * c8;
        p.at(3, 3) = c4 + u * c6 + v * c7 + u * v * c8;

        // solve for corrections
        p.solveForRhs(r, delta);

        // update guess
        answer.subtract(delta);
    }

    // test if inside
    bool inside = true;
    for ( int i = 1; i <= 3; i++ ) {
        if ( answer.at(i) < ( -1. - POINT_TOL ) ) {
            answer.at(i) = -1.;
            inside = false;
        } else if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            answer.at(i) = 1.;
            inside = false;
        }
    }

    return inside;
}

void
FEI3dHexaLin :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double ksi = lcoords.at(1);
    answer.resize(2);

    answer.at(1) = ( 1. - ksi ) * 0.5;
    answer.at(2) = ( 1. + ksi ) * 0.5;
}

void
FEI3dHexaLin :: edgeEvaldNdx(FloatMatrix &answer, int iedge,
                             const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    const auto &edgeNodes = this->computeLocalEdgeMapping(iedge);
    double l = this->edgeComputeLength(edgeNodes, cellgeo);

    answer.resize(2, 1);
    answer.at(1, 1) = -1.0 / l;
    answer.at(2, 1) =  1.0 / l;
}


void
FEI3dHexaLin :: edgeEvaldNdxi(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    answer.resize(2);
    answer(0) = -0.5;
    answer(1) =  0.5;
}


void
FEI3dHexaLin :: edgeLocal2global(FloatArray &answer, int iedge,
                                 const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    FloatArray n;
    const auto &edgeNodes = this->computeLocalEdgeMapping(iedge);
    this->edgeEvalN(n, iedge, lcoords, cellgeo);

    answer.resize(3);
    answer.at(1) = n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(1) +
                   n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(1);
    answer.at(2) = n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(2) +
                   n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(2);
    answer.at(3) = n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(3) +
                   n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(3);
}


double
FEI3dHexaLin :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    const auto &edgeNodes = this->computeLocalEdgeMapping(iedge);
    return 0.5 * this->edgeComputeLength(edgeNodes, cellgeo);
}


IntArray
FEI3dHexaLin :: computeLocalEdgeMapping(int iedge) const
{
    if ( iedge == 1 ) { // edge between nodes 1 2
        return {1, 2};
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        return {2, 3};
    } else if ( iedge == 3 ) { // edge between nodes 3 4
        return {3, 4};
    } else if ( iedge == 4 ) { // edge between nodes 4 1
        return {4, 1};
    } else if ( iedge == 5 ) { // edge between nodes 1 5
        return {1, 5};
    } else if ( iedge == 6 ) { // edge between nodes 2 6
        return {2, 6};
    } else if ( iedge == 7 ) { // edge between nodes 3 7
        return {3, 7};
    } else if ( iedge == 8 ) { // edge between nodes 4 8
        return {4, 8};
    } else if ( iedge == 9 ) { // edge between nodes 5 6
        return {5, 6};
    } else if ( iedge == 10 ) { // edge between nodes 6 7
        return {6, 7};
    } else if ( iedge == 11 ) { // edge between nodes 7 8
        return {7, 8};
    } else if ( iedge == 12 ) { // edge between nodes 8 5
        return {8, 5};
    } else {
        throw std::range_error("invalid edge number");
        //return {};
    }
}

double
FEI3dHexaLin :: edgeComputeLength(const IntArray &edgeNodes, const FEICellGeometry &cellgeo) const
{
    return distance(cellgeo.giveVertexCoordinates( edgeNodes.at(2) ), cellgeo.giveVertexCoordinates( edgeNodes.at(1) ));
}

void
FEI3dHexaLin :: surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double ksi = lcoords.at(1);
    double eta = lcoords.at(2);

    answer.resize(4);
    answer.at(1) = ( 1. + ksi ) * ( 1. + eta ) * 0.25;
    answer.at(2) = ( 1. - ksi ) * ( 1. + eta ) * 0.25;
    answer.at(3) = ( 1. - ksi ) * ( 1. - eta ) * 0.25;
    answer.at(4) = ( 1. + ksi ) * ( 1. - eta ) * 0.25;
}

void
FEI3dHexaLin :: surfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    // Note, this must be in correct order, not just the correct nodes, therefore we must use snodes;
    const auto &snodes = this->computeLocalSurfaceMapping(isurf);

    FloatArray lcoords_hex;

    ///@note Nodal, surface, edge ordering on this class is a mess. No consistency or rules. Have to convert surface->volume coords manually:
#if 1
    if ( isurf == 1 ) { // surface 1 - nodes 1 4 3 2
        lcoords_hex = {-lcoords.at(1), -lcoords.at(2), 1};
    } else if ( isurf == 2 ) { // surface 2 - nodes 5 6 7 8
        lcoords_hex = {-lcoords.at(2), -lcoords.at(1), -1};
    } else if ( isurf == 3 ) { // surface 3 - nodes 1 2 6 5
        lcoords_hex = {-1, -lcoords.at(1), lcoords.at(2)};
    } else if ( isurf == 4 ) { // surface 4 - nodes 2 3 7 6
        lcoords_hex = {-lcoords.at(1), 1, lcoords.at(2)};
    } else if ( isurf == 5 ) { // surface 5 - nodes 3 4 8 7
        lcoords_hex = {1, lcoords.at(1), lcoords.at(2)};
    } else if ( isurf == 6 ) { // surface 6 - nodes 4 1 5 8
        lcoords_hex = {lcoords.at(1), -1, lcoords.at(2)};
    } else {
        OOFEM_ERROR("wrong surface number (%d)", isurf);
    }
#else
    ///@note This would be somewhat consistent at least.
    if ( isurf == 1 ) { // surface 1 - nodes 3 4 8 7
        lcoords_hex = {-1, lcoords.at(1), lcoords.at(2)};
    } else if ( isurf == 2 ) { // surface 2 - nodes 2 1 5 6
        lcoords_hex = {1, lcoords.at(1), lcoords.at(2)};
    } else if ( isurf == 3 ) { // surface 3 - nodes 3 7 6 2
        lcoords_hex = {lcoords.at(1), -1, lcoords.at(2)};
    } else if ( isurf == 4 ) { // surface 4 - nodes 4 8 5 1
        lcoords_hex = {lcoords.at(1), 1, lcoords.at(2)};
    } else if ( isurf == 5 ) { // surface 5 - nodes 3 2 1 4
        lcoords_hex = {lcoords.at(1), lcoords.at(2), -1};
    } else if ( isurf == 6 ) { // surface 6 - nodes 7 6 5 8
        lcoords_hex = {lcoords.at(1), lcoords.at(2), 1};
    } else {
        OOFEM_ERROR("wrong surface number (%d)", isurf);
    }
#endif

    FloatMatrix fullB;
    this->evaldNdx(fullB, lcoords_hex, cellgeo);
    answer.resize(snodes.giveSize(), 3);
    for ( int i = 1; i <= snodes.giveSize(); ++i ) {
        for ( int j = 1; j <= 3; ++j ) {
            answer.at(i, j) = fullB.at(snodes.at(i), j);
        }
    }
}

void
FEI3dHexaLin :: surfaceLocal2global(FloatArray &answer, int iedge,
                                    const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    FloatArray n;

    const auto &nodes = this->computeLocalSurfaceMapping(iedge);
    this->surfaceEvalN(n, iedge, lcoords, cellgeo);

    answer.resize(3);
    answer.at(1) = n.at(1) * cellgeo.giveVertexCoordinates( nodes.at(1) ).at(1) + n.at(2) * cellgeo.giveVertexCoordinates( nodes.at(2) ).at(1) +
                   n.at(3) * cellgeo.giveVertexCoordinates( nodes.at(3) ).at(1) + n.at(4) * cellgeo.giveVertexCoordinates( nodes.at(4) ).at(1);
    answer.at(2) = n.at(1) * cellgeo.giveVertexCoordinates( nodes.at(1) ).at(2) + n.at(2) * cellgeo.giveVertexCoordinates( nodes.at(2) ).at(2) +
                   n.at(3) * cellgeo.giveVertexCoordinates( nodes.at(3) ).at(2) + n.at(4) * cellgeo.giveVertexCoordinates( nodes.at(4) ).at(2);
    answer.at(3) = n.at(1) * cellgeo.giveVertexCoordinates( nodes.at(1) ).at(3) + n.at(2) * cellgeo.giveVertexCoordinates( nodes.at(2) ).at(3) +
                   n.at(3) * cellgeo.giveVertexCoordinates( nodes.at(3) ).at(3) + n.at(4) * cellgeo.giveVertexCoordinates( nodes.at(4) ).at(3);
}

double
FEI3dHexaLin :: surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    FloatArray a, b, dNdksi(4), dNdeta(4);
    const auto &snodes = this->computeLocalSurfaceMapping(isurf);

    double ksi = lcoords.at(1);
    double eta = lcoords.at(2);

    // No need to divide by 1/4, we'll normalize anyway;
    dNdksi.at(1) =  ( 1. + eta );
    dNdksi.at(2) = -( 1. + eta );
    dNdksi.at(3) = -( 1. - eta );
    dNdksi.at(4) =  ( 1. - eta );

    dNdeta.at(1) =  ( 1. + ksi );
    dNdeta.at(2) =  ( 1. - ksi );
    dNdeta.at(3) = -( 1. - ksi );
    dNdeta.at(4) = -( 1. + ksi );

    for ( int i = 1; i <= 4; ++i ) {
        a.add( dNdksi.at(i), cellgeo.giveVertexCoordinates( snodes.at(i) ) );
        b.add( dNdeta.at(i), cellgeo.giveVertexCoordinates( snodes.at(i) ) );
    }

    answer.beVectorProductOf(a, b);
    return answer.normalize() * 0.0625;
}

double
FEI3dHexaLin :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                  const FEICellGeometry &cellgeo) const
{
    FloatArray normal;
    return this->surfaceEvalNormal(normal, isurf, lcoords, cellgeo);
}

IntArray
FEI3dHexaLin :: computeLocalSurfaceMapping(int isurf) const
{
    if ( isurf == 1 ) { // surface 1 - nodes 1 4 3 2
        return  {1, 4, 3, 2};
    } else if ( isurf == 2 ) { // surface 2 - nodes 5 6 7 8
        return  {5, 6, 7, 8};
    } else if ( isurf == 3 ) { // surface 3  - nodes 1 2 6 5
        return  {1, 2, 6, 5};
    } else if ( isurf == 4 ) { // surface 4 - nodes 2 3 7 6
        return  {2, 3, 7, 6};
    } else if ( isurf == 5 ) { // surface 5 - nodes 3 4 8 7
        return {3, 4, 8, 7};
    } else if ( isurf == 6 ) { // surface 6 - nodes 4 1 5 8
        return {4, 1, 5, 8};
    } else {
        throw std::runtime_error("invalid surface number");
    }
}

void
FEI3dHexaLin :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
// Returns the jacobian matrix  J (x,y,z)/(ksi,eta,dzeta)  of the receiver.
{
    FloatMatrix dNduvw, coords;
    this->evaldNdxi(dNduvw, lcoords, cellgeo);
    coords.resize(3, 8);
    for ( int i = 1; i <= 8; i++ ) {
        coords.setColumn(cellgeo.giveVertexCoordinates(i), i);
    }
    jacobianMatrix.beProductOf(coords, dNduvw);
}


double
FEI3dHexaLin :: evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo) const
{
    const auto &fNodes = this->computeLocalSurfaceMapping(iEdge);

    const auto &c1 = cellgeo.giveVertexCoordinates( fNodes.at(1) );
    const auto &c2 = cellgeo.giveVertexCoordinates( fNodes.at(2) );
    const auto &c3 = cellgeo.giveVertexCoordinates( fNodes.at(3) );
    const auto &c4 = cellgeo.giveVertexCoordinates( fNodes.at(4) );

    return (
        c4(2) * ( c1(1) * ( -c2(0) - c3(0) ) + c2(1) * ( c1(0) - c3(0) ) + c3(1) * ( c1(0) + c2(0) ) ) +
        c3(2) * ( c1(1) * ( -c2(0) + c4(0) ) + c2(1) * ( c1(0) + c4(0) ) +                          c4(1) * ( -c1(0) - c2(0) ) ) +
        c2(2) * ( c1(1) * ( c3(0) + c4(0) ) +                          c3(1) * ( -c1(0) - c4(0) ) + c4(1) * ( -c1(0) + c3(0) ) ) +
        c1(2) * ( c2(1) * ( -c3(0) - c4(0) ) + c3(1) * ( c2(0) - c4(0) ) + c4(1) * ( c2(0) + c3(0) ) ) ) * 0.25;
}

std::unique_ptr<IntegrationRule>
FEI3dHexaLin :: giveIntegrationRule(int order, Element_Geometry_Type egt) const
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Cube, order + 6);
    iRule->SetUpPointsOnCube(points, _Unknown);
    return std::move(iRule);
}

std::unique_ptr<IntegrationRule>
FEI3dHexaLin :: giveBoundaryIntegrationRule(int order, int boundary, Element_Geometry_Type egt) const
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Square, order + 2);
    iRule->SetUpPointsOnSquare(points, _Unknown);
    return std::move(iRule);
}
} // end namespace oofem
