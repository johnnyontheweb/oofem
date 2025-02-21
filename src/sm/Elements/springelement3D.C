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
#include <cmath>
#include "springelement3D.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "classfactory.h"
#include "node.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element( SpringElement3D );

SpringElement3D::SpringElement3D( int n, Domain *aDomain ) :
    StructuralElement( n, aDomain )
{
    numberOfDofMans = 2;
    springC1        = springC2 = springC3 = springC4 = springC5 = springC6 = 0.0;
}

void
SpringElement3D::computeStiffnessMatrix( FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep )
{
    /* spring stiffness matrix in local coordinate system (along orientation axis) */
    double dr = d / 2;

    answer.resize( 12, 12 );
    answer.at( 1, 1 ) = answer.at( 7, 7 ) = this->springC1;
    answer.at( 1, 7 ) = answer.at( 7, 1 ) = -this->springC1;

    answer.at( 2, 2 ) = answer.at( 8, 8 ) = this->springC2;
    answer.at( 2, 8 ) = answer.at( 8, 2 ) = -this->springC2;

    answer.at( 3, 3 ) = answer.at( 9, 9 ) = this->springC3;
    answer.at( 3, 9 ) = answer.at( 9, 3 ) = -this->springC3;

    answer.at( 4, 4 )  = answer.at( 10, 10 ) = this->springC4;
    answer.at( 4, 10 ) = answer.at( 10, 4 )  = -this->springC4;

    answer.at( 5, 5 )  = answer.at( 11, 11 ) = this->springC5 + this->springC3 * ( dr * dr );
    answer.at( 5, 11 ) = answer.at( 11, 5 ) = this->springC3 * ( dr * dr ) - this->springC5;

    answer.at( 6, 6 )   = answer.at( 12, 12 ) = this->springC6 + this->springC2 * ( dr * dr );
    answer.at( 6, 12 ) = answer.at( 12, 6 ) = this->springC2 * ( dr * dr ) - this->springC6;

    // rigid link transport terms ------------------------------------------------------
    answer.at( 12, 2 ) = answer.at( 2, 12 ) = -dr * this->springC2;
    answer.at( 8, 12 ) = answer.at( 12, 8 ) = dr * this->springC2;

    answer.at( 11, 3 ) = answer.at( 3, 11 ) = dr * this->springC3;
    answer.at( 11, 9 ) = answer.at( 9, 11 ) = -dr * this->springC3;

    answer.at( 2, 6 ) = answer.at( 6, 2 ) = -dr * this->springC2;
    answer.at( 8, 6 ) = answer.at( 6, 8 ) = dr * this->springC2;

    answer.at( 3, 5 ) = answer.at( 5, 3 ) = dr * this->springC3;
    answer.at( 5, 9 ) = answer.at( 9, 5 ) = -dr * this->springC3;
}


void
SpringElement3D::giveInternalForcesVector( FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord )
{
    answer = this->computeSpringInternalForce( tStep );
    //FloatArray f = this->computeSpringInternalForce( tStep );
    //answer.resize( 12 );
    //for ( int i = 1; i <= 6; i++ ) {
    //    answer.at( i ) = -f.at( i );
    //    answer.at( i+6 ) = f.at( i );
    //}
}

bool
SpringElement3D::computeGtoLRotationMatrix( FloatMatrix &answer )
{
    /*
     * Spring3D is defined as 3D element along orientation axis
     */
    FloatMatrix lcs; // local axes
    int ndofs = 12;
    answer.resize( ndofs, ndofs );
    answer.zero();

    this->giveLocalCoordinateSystem( lcs );
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at( i, j )         = lcs.at( i, j );
            answer.at( i + 3, j + 3 ) = lcs.at( i, j );
            answer.at( i + 6, j + 6 ) = lcs.at( i, j );
            answer.at( i + 9, j + 9 ) = lcs.at( i, j );
        }
    }
    return 1;
}

int
SpringElement3D::giveLocalCoordinateSystem( FloatMatrix &answer )
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    FloatArray lx, ly, lz, help( 3 );
    help.zero();
    lx = dir;
    FloatMatrix rot( 3, 3 );
    double theta = referenceAngle * M_PI / 180.0;

    rot.at( 1, 1 ) = cos( theta ) + pow( lx.at( 1 ), 2 ) * ( 1 - cos( theta ) );
    rot.at( 1, 2 ) = lx.at( 1 ) * lx.at( 2 ) * ( 1 - cos( theta ) ) - lx.at( 3 ) * sin( theta );
    rot.at( 1, 3 ) = lx.at( 1 ) * lx.at( 3 ) * ( 1 - cos( theta ) ) + lx.at( 2 ) * sin( theta );

    rot.at( 2, 1 ) = lx.at( 2 ) * lx.at( 1 ) * ( 1 - cos( theta ) ) + lx.at( 3 ) * sin( theta );
    rot.at( 2, 2 ) = cos( theta ) + pow( lx.at( 2 ), 2 ) * ( 1 - cos( theta ) );
    rot.at( 2, 3 ) = lx.at( 2 ) * lx.at( 3 ) * ( 1 - cos( theta ) ) - lx.at( 1 ) * sin( theta );

    rot.at( 3, 1 ) = lx.at( 3 ) * lx.at( 1 ) * ( 1 - cos( theta ) ) - lx.at( 2 ) * sin( theta );
    rot.at( 3, 2 ) = lx.at( 3 ) * lx.at( 2 ) * ( 1 - cos( theta ) ) + lx.at( 1 ) * sin( theta );
    rot.at( 3, 3 ) = cos( theta ) + pow( lx.at( 3 ), 2 ) * ( 1 - cos( theta ) );

    help.at( 3 ) = 1.0; // up-vector
    // here is ly is used as a temp var
    // double prvect = acos(lx.dotProduct(help));
    // if (prvect < 0.001 || prvect > M_PI - 0.001) { // Check if it is vertical
    if ( fabs( lx.dotProduct( help ) ) > 0.999 ) {
        // Check if it is vertical
        lz = { 1., 0., 0. };
    } else {
        ly.beVectorProductOf( lx, help );
        lz.beVectorProductOf( ly, lx );
    }
    ly.beProductOf( rot, lz );
    ly.normalize();
    lz.beVectorProductOf( lx, ly );
    lz.normalize();

    answer.resize( 3, 3 );
    answer.zero();
    for ( int i = 1; i <= 3; i++ ) {
        answer.at( 1, i ) = lx.at( i );
        answer.at( 2, i ) = ly.at( i );
        answer.at( 3, i ) = lz.at( i );
    }

    return 1;
}

void
SpringElement3D::giveDofManDofIDMask( int inode, IntArray &answer ) const { answer = { D_u, D_v, D_w, R_u, R_v, R_w }; }

FloatArray
SpringElement3D::computeSpringInternalForce( TimeStep *tStep )
{
    FloatArray u;
    this->computeVectorOf( VM_Total, tStep, u );
    FloatArray res;
    FloatMatrix k;
    this->computeStiffnessMatrix( k, TangentStiffness, tStep );
    res.beProductOf( k, u );
    return res;
    //FloatArray ans; ans.resize( 6 );
    //for ( int i = 1; i <= 6; i++ ) {
    //    ans.at( i ) = res.at(i+6)-res.at(i);
    //} 
    //return ans;
}

void
SpringElement3D::computeLumpedMassMatrix( FloatMatrix &answer, TimeStep *tStep )
{
    answer.resize( 12, 12 );
    answer.zero(); // translational dofs only
    answer.at( 1, 1 ) = answer.at( 2, 2 ) = answer.at( 3, 3 ) = this->mass / 2.0;
    answer.at( 7, 7 ) = answer.at( 8, 8 ) = answer.at( 9, 9 ) = this->mass / 2.0;
}

int
SpringElement3D::computeNumberOfGlobalDofs() { return 12; }

void
SpringElement3D::initializeFrom( InputRecord &ir )
{
    StructuralElement::initializeFrom( ir );

    IR_GIVE_FIELD( ir, springC1, _IFT_SpringElement3D_springC1 );
    IR_GIVE_FIELD( ir, springC2, _IFT_SpringElement3D_springC2 );
    IR_GIVE_FIELD( ir, springC3, _IFT_SpringElement3D_springC3 );
    IR_GIVE_FIELD( ir, springC4, _IFT_SpringElement3D_springC4 );
    IR_GIVE_FIELD( ir, springC5, _IFT_SpringElement3D_springC5 );
    IR_GIVE_FIELD( ir, springC6, _IFT_SpringElement3D_springC6 );

    this->mass = 0.0;
    IR_GIVE_OPTIONAL_FIELD( ir, this->mass, _IFT_SpringElement3D_mass );

    IR_GIVE_FIELD( ir, this->dir, _IFT_SpringElement3D_orientation );
    this->dir.normalize();

    this->macroElem = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, this->macroElem, _IFT_SpringElement3D_macroElem );

    this->referenceAngle = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, this->referenceAngle, _IFT_SpringElement3D_refangle );

    this->d = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, this->d, _IFT_SpringElement3D_actAsRigidLink );
}

void SpringElement3D::postInitialize() { 
    this->d = this->computeLength(); // always recalculate length!
}

void SpringElement3D::printOutputAt( FILE *File, TimeStep *tStep )
{
    if ( this->macroElem != 0 ) {
        FloatArray u;
        this->computeVectorOf( VM_Total, tStep, u );
        FloatArray res;
        res.resize( 6 );
        res.at( 1 ) = ( u.at( 7 ) - u.at( 1 ) );
        res.at( 2 ) = ( u.at( 8 ) - u.at( 2 ) );
        res.at( 3 ) = ( u.at( 9 ) - u.at( 3 ) );
        res.at( 4 ) = ( u.at( 10 ) - u.at( 4 ) );
        res.at( 5 ) = ( u.at( 11 ) - u.at( 5 ) );
        res.at( 6 ) = ( u.at( 12 ) - u.at( 6 ) );
        FloatArray resF = this->computeSpringInternalForce( tStep );
        fprintf( File, "SpringElement3D %d dir 3 %.4e %.4e %.4e refangle %.4e disp 6 %.4e %.4e %.4e %.4e %.4e %.4e macroelem %d : %.4e %.4e %.4e %.4e %.4e %.4e\n", this->giveLabel(), 
            this->dir.at( 1 ), this->dir.at( 2 ), this->dir.at( 3 ), this->referenceAngle, res.at( 1 ), res.at( 2 ), res.at( 3 ), res.at( 4 ), res.at( 5 ), res.at( 6 ), this->macroElem, 
            resF.at( 1 ), resF.at( 2 ), resF.at( 3 ), resF.at( 4 ), resF.at( 5 ), resF.at( 6 ) ); // print the first 6 items
    } else { 
        FloatArray res = this->computeSpringInternalForce( tStep );
        fprintf( File, "SpringElement3D %d :%.4e %.4e %.4e %.4e %.4e %.4e\n", this->giveLabel(), res.at(1), res.at(2), res.at(3), res.at(4), res.at(5), res.at(6) ); // print the first 6 items
    }
}

double
SpringElement3D::computeLength()
// Returns the length of the receiver.
{
    double dx, dy, dz, length, leng2;
    Node *nodeA, *nodeB;

    nodeA  = this->giveNode( 1 );
    nodeB  = this->giveNode( 2 );
    dx     = nodeB->giveCoordinate( 1 ) - nodeA->giveCoordinate( 1 );
    dy     = nodeB->giveCoordinate( 2 ) - nodeA->giveCoordinate( 2 );
    dz     = nodeB->giveCoordinate( 3 ) - nodeA->giveCoordinate( 3 );
    leng2  = dx * dx + dy * dy + dz * dz;
    if ( leng2 > 0 )
        length = sqrt( leng2 );
    else
        length = 0;

    return length;
}

void SpringElement3D::computeInitialStressMatrix( FloatMatrix &answer, TimeStep *tStep )
{
    // computes initial stress matrix of receiver (or geometric stiffness matrix)
    answer.resize( 12, 12 );
    answer.zero();

    double l = this->computeLength();
    if ( l > 0 && d==0 ) {
        double N;

        answer.at( 2, 2 ) = 1;
        answer.at( 2, 8 ) = -1;

        answer.at( 3, 3 ) = 1;
        answer.at( 3, 9 ) = -1;

        answer.at( 8, 2 ) = -1;
        answer.at( 8, 8 ) = 1;

        answer.at( 9, 3 ) = -1;
        answer.at( 9, 9 ) = 1;

        FloatMatrix lcs;
        this->giveLocalCoordinateSystem( lcs );

        FloatMatrix transf( 12, 12 );
        this->computeGtoLRotationMatrix( transf );

        answer.rotatedWith( transf, 'n' );
        // ask end forces in g.c.s
        FloatArray endForces;
        this->giveInternalForcesVector( endForces, tStep );

        FloatArray N1, N2;
        IntArray ind( { 1, 2, 3 } );
        IntArray ind2( { 7, 8, 9 } );
        N1.beSubArrayOf( endForces, ind );
        N2.beSubArrayOf( endForces, ind2 );

        FloatArray lx;
        lx.beDifferenceOf( this->giveNode( 2 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );
        lx.normalize();

        // sign of N
        N = ( -N1.dotProduct( lx ) + N2.dotProduct( lx ) ) / 2.;
        answer.times( N / l );
    }
}

} // end namespace oofem