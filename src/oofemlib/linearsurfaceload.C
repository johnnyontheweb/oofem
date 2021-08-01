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

#include "linearsurfaceload.h"
#include "function.h"
#include "floatarray.h"
#include "timestep.h"
#include "classfactory.h"

namespace oofem {
REGISTER_BoundaryCondition(LinearSurfaceLoad);

void
LinearSurfaceLoad :: initializeFrom(InputRecord &ir)
{
    // support for normalized values for loading magnitude
    IR_GIVE_OPTIONAL_FIELD(ir, normVals, _IFT_LinearSurfaceLoad_normVals);
    if (normVals.giveSize() == 0) {
	    normVals.resize(4); // maximum supported values
	    for (int i = 1; i <= 4; i++) {
		    normVals.at(i) = 1.;
	    }
    }
    return BoundaryLoad :: initializeFrom(ir);
}

void 
LinearSurfaceLoad::computeNArray(FloatArray &answer, const FloatArray &coords) const
{
	// compute local isoparametric coordinates of given point
	answer.resize(coords.giveSize());
	answer.zero();
}

void
LinearSurfaceLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
{
    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        OOFEM_ERROR("mode not supported");
    }

    double factor = this->giveTimeFunction()->evaluate(tStep, mode);
    answer.beScaled(factor, componentArray);

	double fplane = loadPlane(coords);
	answer.beScaled(fplane, componentArray);
}

double
LinearSurfaceLoad::loadPlane(const FloatArray &coords)
{
	int n = normVals.giveSize();
	// load plane
	FloatArray xc; xc.resize(n); 
	FloatArray yc; yc.resize(n); 
	if (n == 4) {
		xc.at(1) = -1; xc.at(2) = -1; xc.at(3) = 1; xc.at(4) = 1;
		yc.at(1) = 1; yc.at(2) = -1; yc.at(3) = -1; yc.at(4) = 1;
	} else {
		xc.at(1) = 0; xc.at(2) = 1; xc.at(3) = 0;
		yc.at(1) = 0; yc.at(2) = 0; yc.at(3) = 1;
	}
	FloatMatrix A(3, 3); A.zero();
	FloatArray b(3), x; b.zero(); x.zero();
	A.at(3, 3) = n;
	for (int i = 1; i <= n; i++) {
		//  sum_i x[i] * x[i], sum_i x[i] * y[i], sum_i x[i]
		//	sum_i x[i] * y[i], sum_i y[i] * y[i], sum_i y[i]
		//	sum_i x[i], sum_i y[i], n
		A.at(1, 1) += xc.at(i)*xc.at(i);
		A.at(2, 1) += xc.at(i)*yc.at(i);
		A.at(3, 1) += xc.at(i);

		A.at(1, 2) += xc.at(i)*yc.at(i);
		A.at(2, 2) += yc.at(i)*yc.at(i);
		A.at(3, 2) += yc.at(i);

		A.at(1, 3) += xc.at(i);
		A.at(2, 3) += yc.at(i);

		//{sum_i x[i]*z[i],   sum_i y[i]*z[i],    sum_i z[i]}
		b.at(1) += xc.at(i)*normVals.at(i);
		b.at(2) += yc.at(i)*normVals.at(i);
		b.at(3) += normVals.at(i);
	}
	// solve
	A.solveForRhs(b, x);

	double fplane = x.at(1)*coords.at(1) + x.at(2)*coords.at(2) + x.at(3);
	return fplane;
}

} // end namespace oofem
