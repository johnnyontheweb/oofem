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

#ifndef angle_h
#define angle_h

#include "floatmatrix.h"
#include "floatarray.h"

namespace oofem {
	class Angle
	{
	public:
	static double giveAngleIn3Dplane(FloatArray v1, FloatArray v2, FloatArray vn)
	{
		v1.normalize(); v2.normalize(); vn.normalize();
		double dot = v1.at(1) * v2.at(1) + v1.at(2) * v2.at(2) + v1.at(3) * v2.at(3);
		double det = v1.at(1) * v2.at(2) * vn.at(3) + v2.at(1) * vn.at(2) * v1.at(3) + vn.at(1) * v1.at(2) * v2.at(3) - v1.at(3) * v2.at(2) * vn.at(1) - v2.at(3) * vn.at(2) * v1.at(1) - vn.at(3) * v1.at(2) * v2.at(1);
		return atan2(det, dot);
	}

	static FloatArray rotate(FloatArray vec, FloatArray dir, double theta)
	{
		FloatArray res(3); //res.resize(3);
		FloatArray orig; orig.resize(3); orig.zero(); // always around origin

		double x = vec.at(1); double y = vec.at(2); double z = vec.at(3);
        //if (abs(x) < 1e-7) x = 0; if (abs(y) < 1e-7) y = 0; if (abs(z) < 1e-7) z = 0;
		double a = orig.at(1); // origin of rot. axis
		double b = orig.at(2); double c = orig.at(3);
		double u = dir.at(1); // direction of rot. axis
		double v = dir.at(2); double w = dir.at(3);
		double cosTh = cos(theta); double senTh = sin(theta);
        //if (abs(cosTh) < 1e-7) cosTh = 0; if (abs(senTh) < 1e-7) senTh = 0;

		res.at(1) = (a * (pow(v, 2) + pow(w, 2)) - u * (b * v + c * w - u * x - v * y - w * z)) * (1 - cosTh) + x * cosTh + (-c * v + b * w - w * y + v * z) * senTh;
		res.at(2) = (b * (pow(u, 2) + pow(w, 2)) - v * (a * u + c * w - u * x - v * y - w * z)) * (1 - cosTh) + y * cosTh + (c * u - a * w + w * x - u * z) * senTh;
		res.at(3) = (c * (pow(u, 2) + pow(v, 2)) - w * (a * u + b * v - u * x - v * y - w * z)) * (1 - cosTh) + z * cosTh + (-b * u + a * v - v * x + u * y) * senTh;

		return res;
	}

	}; // end class
} // end namespace oofem
#endif