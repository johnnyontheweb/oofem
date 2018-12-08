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

#ifndef isoheatadvmat_h
#define isoheatadvmat_h

#include "transportmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "scalarfunction.h"

///@name Input fields for IsotropicHeatAdvTransferMaterial
//@{
#define _IFT_IsotropicHeatAdvTransferMaterial_Name "isoheatadv"
#define _IFT_IsotropicHeatAdvTransferMaterial_k "k" ///< Conductivity
#define _IFT_IsotropicHeatAdvTransferMaterial_c "c" ///< Specific heat
#define _IFT_IsotropicHeatAdvTransferMaterial_maturityT0 "maturityt0" ///< Baseline for maturity method
#define _IFT_IsotropicHeatAdvTransferMaterial_d "d"
#define _IFT_IsotropicHeatAdvTransferMaterial_dFunc "funcd"
#define _IFT_IsotropicHeatAdvTransferMaterial_kFunc "funck"
#define _IFT_IsotropicHeatAdvTransferMaterial_cFunc "funcc"
//@}

namespace oofem {

/**
 * This class implements an isotropic linear heat  material. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 */
class IsotropicHeatAdvTransferMaterial : public TransportMaterial
{
protected:
    ScalarFunction conductivity; ///< Conductivity (k in input file).
    ScalarFunction capacity;     ///< Capacity (c in input file).
    ScalarFunction density;      ///< Density (td in input file).
    double maturityT0;           ///< Baseline for maturity mathod
	int funcD=0,funcK=0,funcC=0;

public:
    IsotropicHeatAdvTransferMaterial(int n, Domain * d);
    virtual ~IsotropicHeatAdvTransferMaterial();

    virtual void giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep);

    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual double giveIsotropicConductivity(GaussPoint *gp, TimeStep *tStep);

    virtual double giveCharacteristicValue(MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep);

    virtual double  giveMaturityT0() { return maturityT0; }

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual const char *giveInputRecordName() const { return _IFT_IsotropicHeatAdvTransferMaterial_Name; }
    virtual const char *giveClassName() const { return "IsotropicHeatAdvTransferMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual double give(int aProperty, GaussPoint *gp, TimeStep *tStep);
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
    double giveTemperature(GaussPoint *gp);
	virtual void postInitialize();
};

class IsotropicHeatAdvTransferMaterialStatus : public TransportMaterialStatus
{
public:
    IsotropicHeatAdvTransferMaterialStatus(int n, Domain * d, GaussPoint * g);
    virtual ~IsotropicHeatAdvTransferMaterialStatus();
    virtual void updateYourself(TimeStep *tStep);
};

} // end namespace oofem
#endif // isoheatadvmat_h
