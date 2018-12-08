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

#include "isoheatadvmat.h"
#include "floatmatrix.h"
#include "function.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "engngm.h"

namespace oofem {
REGISTER_Material(IsotropicHeatAdvTransferMaterial);

IsotropicHeatAdvTransferMaterial :: IsotropicHeatAdvTransferMaterial(int n, Domain *d) : TransportMaterial(n, d)
{
    // constructor
    maturityT0 = 0.;
}

IsotropicHeatAdvTransferMaterial :: ~IsotropicHeatAdvTransferMaterial() {
    // destructor
}

IRResultType
IsotropicHeatAdvTransferMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, conductivity, _IFT_IsotropicHeatAdvTransferMaterial_k);
    IR_GIVE_FIELD(ir, capacity, _IFT_IsotropicHeatAdvTransferMaterial_c);
    IR_GIVE_OPTIONAL_FIELD(ir, maturityT0, _IFT_IsotropicHeatAdvTransferMaterial_maturityT0);
    IR_GIVE_OPTIONAL_FIELD(ir, density, _IFT_IsotropicHeatAdvTransferMaterial_d); // old td - reverted to d

	IR_GIVE_OPTIONAL_FIELD(ir, funcD, _IFT_IsotropicHeatAdvTransferMaterial_dFunc);
	if (funcD) density.setMultiplierReference(funcD);
	IR_GIVE_OPTIONAL_FIELD(ir, funcK, _IFT_IsotropicHeatAdvTransferMaterial_kFunc);
	if (funcK) conductivity.setMultiplierReference(funcK);
	IR_GIVE_OPTIONAL_FIELD(ir, funcC, _IFT_IsotropicHeatAdvTransferMaterial_cFunc);
	if (funcC) capacity.setMultiplierReference(funcC);
	
    return Material :: initializeFrom(ir);
}

void IsotropicHeatAdvTransferMaterial::postInitialize() // not called
{
	// we check whether the function exists or not.
	if (funcD > 0) {
		Function *fD = this->giveDomain()->giveFunction(funcD);
		if (fD == NULL) OOFEM_ERROR("Invalid function given");
	}
	if (funcK > 0) {
		Function *fK = this->giveDomain()->giveFunction(funcK);
		if (fK == NULL) OOFEM_ERROR("Invalid function given");
	}
	if (funcC > 0) {
		Function *fC = this->giveDomain()->giveFunction(funcC);
		if (fC == NULL) OOFEM_ERROR("Invalid function given");
	}
}

double
IsotropicHeatAdvTransferMaterial :: give(int aProperty, GaussPoint *gp, TimeStep *tStep)
//
// Returns the value of the property aProperty (e.g. 'k' the conductivity of the receiver).
//
{
    if ( aProperty == 'k' ) { //thermal conductivity [W/m/K]   
        return conductivity.eval( { { "te", giveTemperature(gp) }, { "t", tStep->giveIntrinsicTime() } }, this->giveDomain(), gp, giveTemperature(gp) );
    } else if ( aProperty == 'c' ) { //mass-specific heat capacity [J/kg/K]
        return capacity.eval( { { "te", giveTemperature(gp) }, { "t", tStep->giveIntrinsicTime() } }, this->giveDomain(), gp, giveTemperature(gp) );
    } else if ( aProperty == 'd' && density.isDefined() ) { //density [kg/m3]
        return density.eval( { { "te", giveTemperature(gp) }, { "t", tStep->giveIntrinsicTime() } }, this->giveDomain(), gp, giveTemperature(gp) );
    } else if ( aProperty == HeatCapaCoeff ) { //volume-specific heat capacity [J/m3/K]
        return ( this->give('c', gp, tStep) * this->give('d', gp, tStep) );
    }

    return this->Material :: give(aProperty, gp);
}


void
IsotropicHeatAdvTransferMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );

    ms->setTempField(field);
    ms->setTempGradient(grad);

    ///@todo Shouldn't the conductivity typically depend on the primary field and/or its gradient?
    answer.beScaled(-this->giveIsotropicConductivity(gp, tStep), grad);

    ms->setTempFlux(answer);
}


void
IsotropicHeatAdvTransferMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                          MatResponseMode mode,
                                                          GaussPoint *gp,
                                                          TimeStep *tStep)
{
    /*
     * returns constitutive (conductivity) matrix of receiver
     */
    MaterialMode mMode = gp->giveMaterialMode();
    double cond = this->giveIsotropicConductivity(gp, tStep);

    switch  ( mMode ) {
    case _1dHeat:
        answer.resize(1, 1);
        answer.at(1, 1) = cond;
    case _2dHeat:
        answer.resize(2, 2);
        answer.at(1, 1) = cond;
        answer.at(2, 2) = cond;
        return;

    case _3dHeat:
        answer.resize(3, 3);
        answer.at(1, 1) = cond;
        answer.at(2, 2) = cond;
        answer.at(3, 3) = cond;
        return;

    default:
        OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

double
IsotropicHeatAdvTransferMaterial :: giveIsotropicConductivity(GaussPoint *gp, TimeStep *tStep) {
    return give('k', gp, tStep);
}

double
IsotropicHeatAdvTransferMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                                         GaussPoint *gp,
                                                         TimeStep *tStep)
{
    if ( mode == Capacity ) {
        return ( this->give('c', gp, tStep) * this->give('d', gp, tStep) );
    } else {
        OOFEM_ERROR("unknown mode (%s)", __MatResponseModeToString(mode) );
    }

    return 0.;
}


int
IsotropicHeatAdvTransferMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if (  type == IST_HydrationDegree ) {
        answer.resize(1);
        answer.at(1) = 0.;
        return 1;
    }

    if ( type == IST_Temperature ) {
        answer = FloatArray{ this->giveTemperature(gp) };
        return 1;
    } else if ( type == IST_Density ) {
        answer = FloatArray{ this->give('d', gp, tStep) };
        return 1;
    } else if ( type == IST_HeatCapacity ) {
        answer = FloatArray{ this->give('c', gp, tStep) };
        return 1;
    } else if ( type == IST_ThermalConductivityIsotropic ) {
        answer = FloatArray{ this->give('k', gp, tStep) };
        return 1;
    }

    return TransportMaterial :: giveIPValue(answer, gp, type, tStep);
}






MaterialStatus *
IsotropicHeatAdvTransferMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new IsotropicHeatAdvTransferMaterialStatus(1, domain, gp);
}

IsotropicHeatAdvTransferMaterialStatus :: IsotropicHeatAdvTransferMaterialStatus(int n, Domain *d, GaussPoint *g) : TransportMaterialStatus(n, d, g)
{
    //constructor
}


IsotropicHeatAdvTransferMaterialStatus :: ~IsotropicHeatAdvTransferMaterialStatus()
{
    //destructor
}

void
IsotropicHeatAdvTransferMaterialStatus :: updateYourself(TimeStep *tStep)
{
    TransportMaterialStatus :: updateYourself(tStep);
}

double IsotropicHeatAdvTransferMaterial :: giveTemperature(GaussPoint *gp)
{
    IsotropicHeatAdvTransferMaterialStatus *ms = static_cast< IsotropicHeatAdvTransferMaterialStatus * >( this->giveStatus(gp) );
    return ms->giveTempField().at(1);
}

} // end namespace oofem
