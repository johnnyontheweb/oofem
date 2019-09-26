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

#include "../sm/Elements/SpringElement3D.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(SpringElement3D);

SpringElement3D :: SpringElement3D(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{
    numberOfDofMans = 12;
    springC1 = springC2 = springC3 = springC4 = springC5 = springC6 = 0.0;
}

void
SpringElement3D :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    /* spring stiffness matrix in local coordinate system (along orientation axis) */
    answer.resize(12, 12);
    answer.at(1, 1) = answer.at(7, 7) = this->springC1;
    answer.at(1, 7) = answer.at(7, 1) = -this->springC1;

	answer.at(2, 2) = answer.at(8, 8) = this->springC2;
	answer.at(2, 8) = answer.at(8, 2) = -this->springC2;

	answer.at(3, 3) = answer.at(9, 9) = this->springC3;
	answer.at(3, 9) = answer.at(9, 3) = -this->springC3;

	answer.at(4, 4) = answer.at(10, 10) = this->springC4;
	answer.at(4, 10) = answer.at(10, 4) = -this->springC4;

	answer.at(5, 5) = answer.at(11, 11) = this->springC5;
	answer.at(5, 11) = answer.at(11, 5) = -this->springC5;

	answer.at(6, 6) = answer.at(12, 12) = this->springC6;
	answer.at(6, 12) = answer.at(12, 6) = -this->springC6;
}


void
SpringElement3D :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatArray f = this->computeSpringInternalForce(tStep);
    answer.resize(12);
	answer.at(1) = -f.at(1); answer.at(2) = -f.at(2); answer.at(3) = -f.at(3); 	answer.at(6) = -f.at(4); answer.at(5) = -f.at(5); answer.at(6) = -f.at(6);
	answer.at(7) = f.at(1); answer.at(8) = f.at(2); answer.at(9) = f.at(3);	answer.at(10) = f.at(4); answer.at(11) = f.at(6); answer.at(12) = f.at(6);
}


bool
SpringElement3D :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    /*
     * Spring3D is defined as 3D element along orientation axis
     */
	FloatMatrix lcs; // local axes
	int ndofs = 12;
	answer.resize(ndofs, ndofs);
	answer.zero();

	this->giveLocalCoordinateSystem(lcs);
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			answer.at(i, j) = lcs.at(i, j);
			answer.at(i + 3, j + 3) = lcs.at(i, j);
			answer.at(i + 6, j + 6) = lcs.at(i, j);
			answer.at(i + 9, j + 9) = lcs.at(i, j);
		}
	}
    return 1;
}

int
SpringElement3D::giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
	FloatArray lx, ly, lz, help(3);
	lx = dir;
	help.at(3) = 1.0;         // up-vector
	// here is ly is used as a temp var
	if (fabs(lx.dotProduct(help)) > 0.999) { // Check if it is vertical
		lz = { 1., 0., 0. };
	}
	else {
		ly.beVectorProductOf(lx, help);
		lz.beVectorProductOf(ly, lx);
	}
	// ly.beProductOf(rot, lz); // for theta angle if necessary
	ly.normalize();
	lz.beVectorProductOf(lx, ly);
	lz.normalize();

	answer.resize(3, 3);
	answer.zero();
	for (int i = 1; i <= 3; i++) {
		answer.at(1, i) = lx.at(i);
		answer.at(2, i) = ly.at(i);
		answer.at(3, i) = lz.at(i);
	}

	return 1;
}

void
SpringElement3D :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
	answer = { D_u, D_v, D_w, R_u, R_v, R_w };
}

FloatArray
SpringElement3D :: computeSpringInternalForce(TimeStep *tStep)
{
    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u);
	FloatArray res; res.resize(6);
	res.at(1) = this->springC1 * (u.at(7) - u.at(1));
	res.at(2) = this->springC2 * (u.at(8) - u.at(2));
	res.at(3) = this->springC3 * (u.at(9) - u.at(3));
	res.at(4) = this->springC4 * (u.at(10) - u.at(4));
	res.at(5) = this->springC5 * (u.at(11) - u.at(5));
	res.at(6) = this->springC6 * (u.at(12) - u.at(6));
	return res;
}

void 
SpringElement3D::computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
{
	answer.resize(12, 12); answer.zero(); // translational dofs only
	answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = this->mass / 2.0;
	answer.at(7, 7) = answer.at(8, 8) = answer.at(9, 9) = this->mass / 2.0;
}

int
SpringElement3D :: computeNumberOfGlobalDofs()
{
    return 12;
}

IRResultType
SpringElement3D :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

	IR_GIVE_FIELD(ir, springC1, _IFT_SpringElement3D_springC1);
	IR_GIVE_FIELD(ir, springC2, _IFT_SpringElement3D_springC2);
	IR_GIVE_FIELD(ir, springC3, _IFT_SpringElement3D_springC3);
	IR_GIVE_FIELD(ir, springC4, _IFT_SpringElement3D_springC4);
	IR_GIVE_FIELD(ir, springC5, _IFT_SpringElement3D_springC5);
	IR_GIVE_FIELD(ir, springC6, _IFT_SpringElement3D_springC6);
    
	this->mass = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->mass, _IFT_SpringElement3D_mass);

	IR_GIVE_FIELD(ir, this->dir, _IFT_SpringElement3D_orientation);
	this->dir.normalize();

	this->macroElem = 0;
	IR_GIVE_OPTIONAL_FIELD(ir, this->macroElem, _IFT_SpringElement3D_macroElem);

    return StructuralElement :: initializeFrom(ir);
}

void SpringElement3D :: printOutputAt(FILE *File, TimeStep *tStep)
{
	if (this->macroElem != 0) {
		FloatArray u;
		FloatArray res; res.resize(6);
		res.at(1) = (u.at(7) - u.at(1));
		res.at(2) = (u.at(8) - u.at(2));
		res.at(3) = (u.at(9) - u.at(3));
		res.at(4) = (u.at(10) - u.at(4));
		res.at(5) = (u.at(11) - u.at(5));
		res.at(6) = (u.at(12) - u.at(6));
		this->computeVectorOf(VM_Total, tStep, u);
		fprintf(File, "SpringElement3D %d dir 3 %.4e %.4e %.4e disp 6 %.4e %.4e %.4e %.4e %.4e %.4e macroelem %d : %.4e %.4e %.4e %.4e %.4e %.4e\n", this->giveLabel(), this->dir.at(1), this->dir.at(2), this->dir.at(3), res.at(1),res.at(2),res.at(3), res.at(4), res.at(5), res.at(6), this->macroElem, this->computeSpringInternalForce(tStep));
	} else {
		fprintf(File, "SpringElement3D %d :%.4e %.4e %.4e %.4e %.4e %.4e\n", this->giveLabel(), this->computeSpringInternalForce(tStep));
	}
}
} // end namespace oofem