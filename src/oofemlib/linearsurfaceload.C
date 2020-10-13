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

IRResultType
LinearSurfaceLoad :: initializeFrom(InputRecord *ir)
{
	IRResultType result;                // Required by IR_GIVE_FIELD macro

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
	//// compute local isoparametric coordinates of given point
	//FloatArray ksi; ksi.resize(2); ksi.zero();
	//double w = 1. / coords.giveSize();
	answer.resize(coords.giveSize());
	answer.zero();

	//ksi = coords.at(1);
	//if ((ksi < -1.0) || (ksi > 1.0)) {
	//	OOFEM_WARNING("point out of receiver, skipped", 1);
	//	return;
	//}

	//for (int j=1; j <= normVals.giveSize(); j++) {
	//	int i=j-1;
	//	if (i == 0) j = normVals.giveSize();

	//	answer.at(j) = (1. - ksi) * w;
	//}
}

void
LinearSurfaceLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
{
    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        OOFEM_ERROR("mode not supported");
    }

    double factor = this->giveTimeFunction()->evaluate(tStep, mode);
    answer.beScaled(factor, componentArray);
}
} // end namespace oofem
