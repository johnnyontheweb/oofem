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

#include "misesmatgrad.h"
#include "stressvector.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Material(MisesMatGrad);

/////////////////////////////////////////////////////////////////
//gradient regularization of Mises plasticity coupled with isotropic damage
/////////////////////////////////////////////////////////////////

MisesMatGrad :: MisesMatGrad(int n, Domain *d) : MisesMat(n, d), GradientDamageMaterialExtensionInterface(d)
{}


bool
MisesMatGrad :: hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _1dMat || mode == _PlaneStrain || mode == _3dMat;
}


void
MisesMatGrad :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}



void
MisesMatGrad :: giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _1dMat:
        answer = give1dStressStiffMtrx(mode, gp, tStep);
        break;
    case _PlaneStrain:
        answer = givePlaneStrainStiffMtrx(mode, gp, tStep);
        break;
    case _3dMat:
        answer = give3dMaterialStiffnessMatrix(mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
MisesMatGrad :: giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _1dMat:
        give1dKappaMatrix(answer, mode, gp, tStep);
        break;
    case _PlaneStrain:
        givePlaneStrainKappaMatrix(answer, mode, gp, tStep);
        break;
    case _3dMat:
        give3dKappaMatrix(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
MisesMatGrad :: giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _1dMat:
        give1dGprime(answer, mode, gp, tStep);
        break;
    case _PlaneStrain:
        givePlaneStrainGprime(answer, mode, gp, tStep);
        break;
    case _3dMat:
        give3dGprime(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
MisesMatGrad :: giveGradientDamageStiffnessMatrix_dd_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _1dMat:
        giveInternalLength(answer, mode, gp, tStep);
        break;
    case _PlaneStrain:
        giveInternalLength(answer, mode, gp, tStep);
        break;
    case _3dMat:
        giveInternalLength(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}



FloatMatrixF<1,1>
MisesMatGrad :: give1dStressStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    double E = this->linearElasticMaterial.give('E', gp);
    if ( mode != TangentStiffness ) {
        return {E};
    }

    MisesMatGradStatus *status = static_cast< MisesMatGradStatus * >( this->giveStatus(gp) );
    double tempKappa = status->giveTempCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double dKappa = ( tempKappa - status->giveCumulativePlasticStrain() );
    /********************************************************************/
    double tempDamage = status->giveTempDamage();
    double damage = status->giveDamage();
    /*********************************************************************/
    double nlKappa =  status->giveNonlocalCumulatedStrain();
    double kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
    if ( dKappa <= 0.0 ) {
        return {( 1. - tempDamage ) * E};
    }


    // === plastic loading ===
    const FloatArray &stressVector = status->giveTempEffectiveStress();
    double stress = stressVector.at(1);

    auto tangent = ( 1. - tempDamage ) * E * H / ( E + H );
    if ( tempDamage > damage ) {
        tangent -= ( 1. - mParam ) * computeDamageParamPrime(kappa) * E / ( E + H ) * signum(stress) * stress;
    }
    return {tangent};
}


FloatMatrixF<4,4>
MisesMatGrad :: givePlaneStrainStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto d = this->linearElasticMaterial.givePlaneStrainStiffMtrx(mode, gp, tStep);
    if ( mode != TangentStiffness ) {
        return d;
    }

    MisesMatGradStatus *status = static_cast< MisesMatGradStatus * >( this->giveStatus(gp) );
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double kappa = status->giveCumulativePlasticStrain();
    double dKappa = tempKappa - kappa;
    double tempDamage = status->giveTempDamage();
    double damage = status->giveDamage();

    if ( dKappa <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        return d;
    }

    // === plastic loading ===
    // yield stress at the beginning of the step
    double sigmaY = this->giveS( gp, tStep) + H * kappa;
    // trial deviatoric stress and its norm
    StressVector trialStressDev(status->giveTrialStressDev(), _PlaneStrain);
    double trialS = trialStressDev.computeStressNorm();
    // volumetric stress
    //double trialStressVol = status->giveTrialStressVol();
    // one correction term
    FloatMatrix stiffnessCorrection(4, 4);
    stiffnessCorrection.beDyadicProductOf(trialStressDev, trialStressDev);
    double factor = -2. * sqrt(6.) * G * G / trialS;
    double factor1 = factor * sigmaY / ( ( H + 3. * G ) * trialS * trialS );
    stiffnessCorrection.times(factor1);
    d += FloatMatrixF<4,4>(stiffnessCorrection);
    // another correction term
    stiffnessCorrection.zero();
    stiffnessCorrection.at(1, 1) = stiffnessCorrection.at(2, 2) = stiffnessCorrection.at(3, 3) = 2. / 3.;
    stiffnessCorrection.at(1, 2) = stiffnessCorrection.at(1, 3) = stiffnessCorrection.at(2, 1) = -1. / 3.;
    stiffnessCorrection.at(2, 3) = stiffnessCorrection.at(3, 1) = stiffnessCorrection.at(3, 2) = -1. / 3.;
    stiffnessCorrection.at(4, 4) = 0.5;
    double factor2 = factor * dKappa;
    stiffnessCorrection.times(factor2);
    d += FloatMatrixF<4,4>(stiffnessCorrection);
    //influence of damage
    d *= 1 - tempDamage;
    if ( tempDamage > damage ) {
        const FloatArray &effStress = status->giveTempEffectiveStress();
        double nlKappa = status->giveNonlocalCumulatedStrain();
        kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
        double omegaPrime = computeDamageParamPrime(kappa);
        double scalar = -omegaPrime *sqrt(6.) * G / ( 3. * G + H ) / trialS;
        stiffnessCorrection.beDyadicProductOf(effStress, trialStressDev);
        stiffnessCorrection.times( scalar * ( 1. - mParam ) );
        d += FloatMatrixF<4,4>(stiffnessCorrection);
    }
    return d;
}


FloatMatrixF<6,6>
MisesMatGrad :: give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< MisesMatGradStatus * >( this->giveStatus(gp) );
    auto d = this->linearElasticMaterial.give3dMaterialStiffnessMatrix(mode, gp, tStep);
    // start from the elastic stiffness
    if ( mode != TangentStiffness ) {
        return d;
    }

    double tempKappa = status->giveTempCumulativePlasticStrain();
    double kappa = status->giveCumulativePlasticStrain();
    double dKappa = tempKappa - kappa;

    if ( dKappa > 0.0 ) {
        double tempDamage = status->giveTempDamage();
        double damage = status->giveDamage();
        double sigmaY = this->giveS(gp, tStep) + H * kappa;
        // trial deviatoric stress and its norm
        const FloatArrayF<6> trialStressDev = status->giveTrialStressDev();
        /*****************************************************/
        //double trialStressVol = status->giveTrialStressVol();
        /****************************************************/
        double trialS = computeStressNorm(trialStressDev);
        // one correction term
        double factor = -2. * sqrt(6.) * G * G / trialS;
        double factor1 = factor * sigmaY / ( ( H + 3. * G ) * trialS * trialS );
        d += factor1 * dyad(trialStressDev, trialStressDev);
        // another correction term
        d += factor * dKappa * I_dev6;
        //influence of damage
        d *= (1. - tempDamage);
        if ( tempDamage > damage ) {
            const FloatArrayF<6> effStress = status->giveTempEffectiveStress();
            double nlKappa = status->giveNonlocalCumulatedStrain();
            kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
            double omegaPrime = computeDamageParamPrime(kappa);
            double scalar = -omegaPrime *sqrt(6.) * G / ( 3. * G + H ) / trialS;
            d += scalar * ( 1. - mParam ) * dyad(effStress, trialStressDev);
        }
    }
    return d;
}


void
MisesMatGrad :: give1dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(1, 1);
    answer.zero();
    double E = linearElasticMaterial.give('E', gp);
    MisesMatGradStatus *status = static_cast< MisesMatGradStatus * >( this->giveStatus(gp) );
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double dKappa = tempKappa - status->giveCumulativePlasticStrain();
    const FloatArray &effStress = status->giveTempEffectiveStress();
    double stress = effStress.at(1);
    if ( dKappa > 0 ) {
        double trialS = signum(stress);
        double factor = trialS * E / ( E + H );
        answer.at(1, 1) = factor;
    }
}


void
MisesMatGrad :: givePlaneStrainKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MisesMatGradStatus *status = static_cast< MisesMatGradStatus * >( this->giveStatus(gp) );
    StressVector trialStressDev(status->giveTrialStressDev(), _PlaneStrain);
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double dKappa = tempKappa - status->giveCumulativePlasticStrain();
    answer.resize(1, 4);
    answer.zero();
    if ( dKappa > 0 ) {
        double trialS = trialStressDev.computeStressNorm();
        answer.at(1, 1) = trialStressDev.at(1);
        answer.at(1, 2) = trialStressDev.at(2);
        answer.at(1, 3) = trialStressDev.at(3);
        answer.at(1, 4) = trialStressDev.at(4);
        double factor = sqrt(6.) * G / ( 3. * G + H ) / trialS;
        answer.times(factor);
    }
}


void
MisesMatGrad :: give3dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MisesMatGradStatus *status = static_cast< MisesMatGradStatus * >( this->giveStatus(gp) );
    const FloatArray &trialStressDev = status->giveTrialStressDev();
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double dKappa = tempKappa - status->giveCumulativePlasticStrain();
    answer.resize(1, 6);
    answer.zero();
    if ( dKappa > 0 ) {
        double trialS = computeStressNorm(trialStressDev);
        for ( int i = 1; i <= 6; i++ ) {
            answer.at(1, i) = trialStressDev.at(i);
        }

        double factor = sqrt(6.) * G / ( 3. * G + H ) / trialS;
        answer.times(factor);
    }
}


void
MisesMatGrad :: give1dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MisesMatGradStatus *status = static_cast< MisesMatGradStatus * >( this->giveStatus(gp) );
    double damage, tempDamage;
    double nlKappa, kappa;
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double gPrime;
    answer.resize(1, 1);
    damage = status->giveDamage();
    tempDamage = status->giveTempDamage();
    nlKappa =  status->giveNonlocalCumulatedStrain();
    kappa = mParam * nlKappa + ( 1 - mParam ) * tempKappa;
    if ( ( tempDamage - damage ) > 0 ) {
        const FloatArray &tempEffStress = status->giveTempEffectiveStress();
        answer.at(1, 1) = tempEffStress.at(1);
        gPrime = computeDamageParamPrime(kappa);
        answer.times(gPrime * mParam);
    } else {
        answer.zero();
    }
}

void
MisesMatGrad :: givePlaneStrainGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MisesMatGradStatus *status = static_cast< MisesMatGradStatus * >( this->giveStatus(gp) );
    double damage, tempDamage;
    double nlKappa, kappa;
    answer.resize(4, 1);
    answer.zero();
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double gPrime;
    damage = status->giveDamage();
    tempDamage = status->giveTempDamage();
    nlKappa =  status->giveNonlocalCumulatedStrain();
    kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
    if ( ( tempDamage - damage ) > 0 ) {
        const FloatArray &tempEffStress = status->giveTempEffectiveStress();
        answer.at(1, 1) = tempEffStress.at(1);
        answer.at(2, 1) = tempEffStress.at(2);
        answer.at(3, 1) = tempEffStress.at(3);
        answer.at(4, 1) = tempEffStress.at(4);
        gPrime = computeDamageParamPrime(kappa);
        answer.times(gPrime * mParam);
    } else {
        answer.zero();
    }
}

void
MisesMatGrad :: give3dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MisesMatGradStatus *status = static_cast< MisesMatGradStatus * >( this->giveStatus(gp) );
    answer.resize(6, 1);
    answer.zero();
    double damage, tempDamage;
    double nlKappa, kappa;
    double gPrime;
    double tempKappa = status->giveTempCumulativePlasticStrain();
    damage = status->giveDamage();
    tempDamage = status->giveTempDamage();
    nlKappa =  status->giveNonlocalCumulatedStrain();
    kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
    if ( ( tempDamage - damage ) > 0 ) {
        const FloatArray &tempEffStress = status->giveTempEffectiveStress();
        for ( int i = 1; i <= 6; i++ ) {
            answer.at(i, 1) = tempEffStress.at(i);
        }

        gPrime = computeDamageParamPrime(kappa);
        answer.times(gPrime * mParam);
    } else {
        answer.zero();
    }
}

void
MisesMatGrad :: giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(1, 1);
    answer.at(1, 1) = L;
}




void
MisesMatGrad :: giveRealStressVectorGradientDamage(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep)
{
    auto status = static_cast< MisesMatGradStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    MisesMat :: performPlasticityReturn(totalStrain, gp, tStep);
    status->letTempStrainVectorBe(totalStrain);
    double tempDamage = computeDamage(gp, tStep);
    const auto &tempEffStress = status->giveTempEffectiveStress();
    answer1.beScaled(1.0 - tempDamage, tempEffStress);
    answer2 = status->giveTempCumulativePlasticStrain();

    status->setNonlocalCumulatedStrain(nonlocalCumulatedStrain);
    status->setTempDamage(tempDamage);
    status->letTempEffectiveStressBe(tempEffStress);
    status->letTempStressVectorBe(answer1);
}


double
MisesMatGrad :: computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< MisesMatGradStatus * >( this->giveStatus(gp) );
    double localCumPlastStrain = status->giveTempCumulativePlasticStrain();
    double nlCumPlastStrain = status->giveNonlocalCumulatedStrain();

    return mParam * nlCumPlastStrain + ( 1 - mParam ) * localCumPlastStrain;
}


void
MisesMatGrad :: giveNonlocalInternalForces_N_factor(double &answer, double nlDamageDrivingVariable, GaussPoint *gp, TimeStep *tStep)
{
    answer = nlDamageDrivingVariable;
}

void
MisesMatGrad :: giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nlDamageDrivingVariable_grad, GaussPoint *gp, TimeStep *tStep)
{
    answer = nlDamageDrivingVariable_grad;
    answer.times(internalLength * internalLength);
}




void
MisesMatGrad :: computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep)
{
    MisesMatGradStatus *status = static_cast< MisesMatGradStatus * >( this->giveStatus(gp) );
    answer =  status->giveTempCumulativePlasticStrain();
}



void
MisesMatGrad :: initializeFrom(InputRecord &ir)
{
    MisesMat :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, L, _IFT_MisesMatGrad_l);
    if ( L < 0.0 ) {
        L = 0.0;
    }

    mParam = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, _IFT_MisesMatGrad_m);
}


MisesMatGradStatus :: MisesMatGradStatus(GaussPoint *g) :
    MisesMatStatus(g)
{}


void
MisesMatGradStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status {");
    fprintf(file, "kappa %f, damage %f ", this->kappa, this->damage);
    fprintf(file, "}\n");
}


void
MisesMatGradStatus :: initTempStatus()
{
    // MisesMatStatus::initTempStatus();
    StructuralMaterialStatus :: initTempStatus();

    if ( plasticStrain.giveSize() == 0 ) {
        if ( gp->giveMaterialMode() == _1dMat ) {
            plasticStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector(_1dMat) );
        } else if ( gp->giveMaterialMode() == _PlaneStrain ) {
            plasticStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector(_PlaneStrain) );
        } else if ( gp->giveMaterialMode() == _PlaneStress ) {
            plasticStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector(_PlaneStress) );
        } else if ( gp->giveMaterialMode() == _3dMat ) {
            plasticStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector(_3dMat) );
        }

        plasticStrain.zero();
    }

    tempDamage = damage;
    tempPlasticStrain = plasticStrain;
    tempKappa = kappa;
    trialStressD.clear();
}


void
MisesMatGradStatus :: updateYourself(TimeStep *tStep)
{
    MisesMatStatus :: updateYourself(tStep);
}
} // end namespace oofem
