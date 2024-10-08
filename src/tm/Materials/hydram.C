/* TODO
 *
 * could only one HydrationModel instance be used during staggered analysis for both tm and sm analysis?
 * if initialized in one of the slave problems, it can hardly be accessed from the other domain,
 * if initialized in master staggered problem, it could be assigned to both of the slave problems (setHydrationModel(*)), but it's not very clean
 * in the future, when more types of hydration models are implemented, a separate input file for the hydration model could be used.
 *
 * the staggered problem should be extended to enable
 * 1) different time step lengths - done by dtf
 * 2) alternate stepping - several time steps performed in one analysis, while only one in the other
 *    e.g. plastic evolution - short time x tm analysis without load changes
 *
 * should provide a function to compute we (water consumption) based on chemical cement composition and/or concrete mixture
 * w/c and a/c or [c/V .. result of wc,ac,rho] ratio.
 *
 * UpdateInternalState:
 *   maybe hydration model should ensure that the hydration degree is NOT recalculated in isothermal analysis for each gp,
 *  when each gp refers to the material level status. Maybe it will be necessary to save time stamp in status?
 *  but even then, there is no way to check if the state vector isn't changed.
 * ==> keep it at caller (only HellMat is now the case)
 * in case of the  planned several hydration model statuses + interpolation, internal state update must be done for each gp,
 *  and HydrationModel must keep track of the status variables extents, to select the appropriate range for interpolation.
 *  ??? / choose 'typical' integration points at beginning of analysis, compute hydration there
 \ use extreme computed values, which may be at different places in time
 *
 */

#include "tm/Materials/hydram.h"
#include "gausspoint.h"
#include "datastream.h"
#include "mathfem.h"
#include "contextioerr.h"

#include <cstdio>

namespace oofem {
// ======= class HydrationModelStatus implementation =======
HydrationModelStatus :: HydrationModelStatus(GaussPoint *g) : MaterialStatus(g) { }

void
HydrationModelStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    fprintf(file, " ksi %.5f", hydrationDegree);
}

void
HydrationModelStatus :: initTempStatus()
// initializes temp status - resets temp variables to saved equilibrium state
// Normally, this is used to clear temp status variables during equilibrium iteration
// (e.g. plastic multiplier, active yield surface, ...)
// In HellmichMaterial, there is number of auxiliary status values that only need to be evaluated once for particular
// time, and have to recalculated only if the iteration is restarted.

// In the hydration model, this means resetting to hydration degree at start of step.
// This should only be done at step restart or restoring context (InitStepIncrements->initGpForNewStep)

// Because hydration model initTempStatus is never called directly, it should be determined in master status
// initTempStatus(), whether hydration model init is necessary.
{
    tempHydrationDegree = hydrationDegree;
}

void
HydrationModelStatus :: updateYourself(TimeStep *tStep)
// update after new equilibrium state reached - prepare status for saving context and for next time increment
{
    hydrationDegree = tempHydrationDegree;
}

//modified_zh 25.6.2004 obj = NULL is set in .h
void
HydrationModelStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    if ( !stream.write(hydrationDegree) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}

//modified_zh 25.6.2004 obj = NULL is set in .h
void
HydrationModelStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    if ( !stream.read(hydrationDegree) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    tempHydrationDegree = hydrationDegree; // for oofeg
}


// ======= class HydrationModel implementation =======
// default constructor - Lafarge mixture
HydrationModel :: HydrationModel() : Material(0, nullptr)
{
    setMixture(mtLafarge);
    useFindRoot = frMixed;
}

HydrationModel :: HydrationModel(MixtureType mix, FindRootMethod usefr) : Material(0, nullptr)
{
    setMixture(mix);
    if ( ( usefr ) && ( usefr <= frMixed ) ) {
        useFindRoot = usefr;
    } else {
        OOFEM_ERROR("unknown FindRootMethod");
    }
}

void
HydrationModel :: initializeFrom(InputRecord &ir)
{
    double value;

    //hydration>0  ->  initial hydration degree
    initialHydrationDegree = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, initialHydrationDegree, _IFT_HydrationModel_hydration);
    if ( initialHydrationDegree >= 0. ) {
        OOFEM_LOG_INFO("HydrationModel: Hydration from %.2f.", initialHydrationDegree);
    } else {
        throw ValueInputException(ir, _IFT_HydrationModel_hydration, "must be between 0 and 1");
    }

    if ( ir.hasField(_IFT_HydrationModel_c60mix) ) {
        OOFEM_LOG_INFO("HydrationModel: Model parameters for Skanska C60/75 mixture.");
        setMixture(mtC60);
    }

    timeScale = 1.;
    value = -1.;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_HydrationModel_timeScale);
    if ( value >= 0. ) {
        timeScale = value;
        OOFEM_LOG_INFO("HydrationModel: Time scale set to %.0f", timeScale);
    }

    // Optional direct input of material parameters
    le = 0;
    value = -1.;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_HydrationModel_hheat);
    if ( value >= 0 ) {
        le = value;
        OOFEM_LOG_INFO("HydrationModel: Latent heat of hydration set to %.0f", le);
    }

    value = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_HydrationModel_cv);
    if ( value >= 0 ) {
        cv = value;
        OOFEM_LOG_INFO("HydrationModel: Cement content set to %.0f kg/m3", cv);
        we = 0.23 * cv;
    }

    value = -1.;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_HydrationModel_water);
    if ( value >= 0 ) {
        we = value;
    }

    if ( cv || ( value >= 0 ) ) {
        OOFEM_LOG_INFO("HydrationModel: Water consumption for hydration set to %.0f kg/m3", we);
    }
}

void
HydrationModel :: setMixture(MixtureType mix)
{
    // === Initializes material parameters for given mixture ===
    mixture = mix;
    // Normalized chemical affinity regression function coefficients
    // (Lafarge mixture)
    if ( mix == mtLafarge ) {
        aa = 7.313;
        ba = 10.46;
        ca = 169.3;
        da = 4.370;

        e0 = 0.05; // ksi_0
    }
    // Huber mixture
    else if ( mix == mtHuber ) {
        aa = 15.93;
        ba = 7.33;
        ca = 0.855e8;
        da = 26.7;

        e0 = 0.10; // ksi_0
    }
    // Skanska C60/75 mixture
    else if ( mix == mtC60 || mix == mtC100 ) { // used for C100 too
        aa = 8.5;
        ba = 5.0;
        ca = 300;
        da = 10;

        e0 = 0.05; // ksi_0
    } else {
        OOFEM_ERROR("Unknown mixture type!");
    }

    ear = 4000; // activation term [K]
    if ( !le ) {
        le = 190000;   // latent heat [kJ/m3]
    }

    // rc = 2428; // heat capacity [kJ/m3K] - set in master material input
}


// === Material functions ===
double
HydrationModel :: affinity(double ksi) const
// Returns the normalized chemical affinity A~(ksi) [1/s]
{
    if ( ksi < e0 ) {
        ksi = e0;
    }

    return aa * ( 1 - exp(-ba * ksi) ) / ( 1 + ca * pow(ksi, da) );
}

double
HydrationModel :: dAdksi(double ksi) const
// Returns the derivation of chemical affinity dA~/dksi(ksi)
{
    if ( ksi < e0 ) {
        return 0;
    }

    double ksinad = pow(ksi, da);
    double enaksi = exp(-ba * ksi);

    return aa * ( ca * da * ksinad * ( enaksi - 1 ) / ksi + ba * enaksi * ( ca * ksinad + 1 ) ) / pow(1 + ca * ksinad, 2);
}

double
HydrationModel :: localResidual(double dks) const // G(dksi) 4.12
// !!! uses auxiliary ksi, dt, T in hydration model instance
{
    return dks - auxdt * affinity(auxksi + dks) * exp(-ear / auxT) * ( 1. + auxh * auxh ) / 2.;
}

double
HydrationModel :: dksidT(double ksi, double T, double h, double dt) const
/*
 * !!! should use state vector to enable adding state variables, but is OK on the level of dIntSource/dState
 * G = dksi - dt * aff(ksi) * exp(-ear/T) * (1+h^2)/2 = 0
 * dksi/dT = - dG/dT / dG/dksi, dksi->0                        dksi/dksi=1
 *       = - -dt*aff(ksi)*(--ear/T^2)*exp(-ear/T)*(1+h^2)/2  /  (1 - dt*daff/dksi * exp(-ear/T) * (1+h^2)/2)
 *
 * dh/dT, dh/dksi uz nejde, od toho je iterace ksi(h,T)
 */
{
    double aux = dt * exp(-ear / T) * ( 1 + h * h ) / 2;
    return ( aux * affinity(ksi) * ear / ( T * T ) ) / ( 1 - aux * dAdksi(ksi) );
}

double
HydrationModel :: dksidh(double ksi, double T, double h, double dt) const
/*
 * G = dks - dt * aff(ksi) * exp(-ear/T) * (1+h^2)/2 = 0
 * dksi/dh = - dG/dh / dG/dksi
 *       = - -dt*aff(ksi) *exp(-ear/T) * h  /  (1 - dt*daff/dksi * exp(-ear/T) * (1+h^2)/2)
 */
{
    double aux = dt * exp(-ear / T);
    return ( aux * affinity(ksi) * h ) / ( 1 - aux * dAdksi(ksi) * ( 1 + h * h ) / 2 );
}

MaterialStatus *
HydrationModel :: giveStatus(GaussPoint *gp) const
/*
 * Returns the hydration model status obtained from gp associated material status or from model associated status in case of isothermal analysis
 * Creates the hydration model status if necessary.
 */
{
    HydrationModelStatusInterface *hmi = static_cast< HydrationModelStatusInterface * >( this->giveStatus(gp)->giveInterface(HydrationModelStatusInterfaceType) );
    HydrationModelStatus *status = nullptr;
    if ( hmi ) {
        status = hmi->giveHydrationModelStatus();
        if ( !status ) {
            status = static_cast< HydrationModelStatus * >( this->CreateStatus(gp) );
            hmi->setHydrationModelStatus(status);
        }
    } else {
        OOFEM_ERROR("Master status undefined.");
    }

    return status;
}

void
HydrationModel :: computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const
/*
 * Returns the hydration heat generated in gp during time step tStep.
 * Can be overriden to return also water consumption as second component of the internal source vector.
 */
// !!! Works only for current TimeStep, no check done
// !!! Should suffice to return val(1) without moisture analysis.
{
    val.resize(2);
    val(0) = le * giveHydrationDegree(gp, tStep, mode); // heat SOURCE
    val(1) = -we *giveHydrationDegree(gp, tStep, mode); // water CONSUMPTION
}

double
HydrationModel :: _giveCharacteristicValue(double T, double h, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep) const
// Transport status needs to be obtained from master status, it's better to pass as a parameter
// to enable usage from structural material
{
    if ( rmode == IntSource || rmode == IntSource_hh || rmode == IntSource_ww || rmode == IntSource_wh || rmode == IntSource_hw ) {
        return computeIntSource(T, h, gp, tStep, rmode);
    } else {
        OOFEM_ERROR("wrong MatResponseMode.");
    }

    // return 0.;
}

double
HydrationModel :: computeHydrationDegreeIncrement(double ksi, double T, double h, double dt)
// Computes the hydration degree increment according to given hydration degree, temperature and time increment.
// !!! sets aux values in material
{
    double result = 0.0;

    if ( ksi < 1.0 ) {
        auxksi = ksi;
        auxT = T;
        auxh = h;
        auxdt = dt;
        switch ( useFindRoot ) {
        case frRegula:  result = regulafindroot();
            break;
        case frBinTree: result = bintreefindroot();
            break;
        case frMixed:   result = mixedfindroot();
            break;
#ifdef DEBUG
        default: {
            OOFEM_ERROR("unknown FindRootMethod %d", useFindRoot);
        }
#endif
        }

        if ( ksi + result > 1.0 ) {
#ifdef VERBOSE
            OOFEM_LOG_INFO("temp=%.12f, dksi %.15f -> %f\n", T, result, 1.0 - ksi);
#endif
            result = 1.0 - ksi;
        }
    } else {
        result = 0.;
    }

    return result;
}

double
HydrationModel :: computeIntSource(double T, double h, GaussPoint *gp, TimeStep *tStep, MatResponseMode rmode) const
/*
 * Computes and returns the derivatives of the material-generated Internal Source with respect to the tm state vector.
 * Used in IntSourceLHS matrix for quadratic convergency of global iteration.
 * State vector must contain relative humidity, not water content
 *
 * Called from giveCharacteristicValue - IntSource_hh, IntSource_ww, IntSource_wh, IntSource_hw.
 * IntSource_hh = -lksi * dksidT (dHeat / dT)
 * IntSource_ww = wksi * dksidh (dWater / dh)
 * IntSource_hw = -lksi * dksidh (dHeat / dh)
 * IntSource_wh = -wksi * dksidT (dWater / dT)
 */
{
    // prepare ksi, T, h, dt
    double ksi = giveHydrationDegree(gp, tStep, VM_Total);
    // !!! timeScale
    double dt = tStep->giveTimeIncrement() * timeScale;

    if ( ksi < 1.0 ) {
        switch ( rmode ) {
        case IntSource:
        case IntSource_hh: return -le *dksidT(ksi, T, h, dt);

        case IntSource_ww: return we *dksidh(ksi, T, h, dt);

        case IntSource_hw: return -le *dksidh(ksi, T, h, dt);

        case IntSource_wh: return -we *dksidT(ksi, T, h, dt);

        default: OOFEM_ERROR("Wrong MatResponseMode.");
        }
    }

    return 0;
}

// === public services ===

double
HydrationModel :: giveHydrationDegree(GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const
// returns the hydration degree in integration point gp
{
    HydrationModelStatus *status = static_cast< HydrationModelStatus * >( giveStatus(gp) );
    double ksi = status->giveTempHydrationDegree();
    if ( mode == VM_Incremental ) {
        ksi -= status->giveHydrationDegree();
    }

    return ksi;
}

void
HydrationModel :: updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *tStep)
/*
 * Recalculates the hydration degree according to the given new state vector and time increment, using equilib status
 *  State vector is supposed to contain [1]->temperature, [2]->relative humidity!
 *  caller should ensure that this is called only when state vector is changed
 */
{
    double dksi, T = 0., h = 1.;
    // get hydration model status associated with integration point
    HydrationModelStatus *status = static_cast< HydrationModelStatus * >( giveStatus(gp) );

    if ( vec.giveSize() ) {
        T = vec(0);
        if ( vec.giveSize() > 1 ) {
            h = vec(1);
        } else {
            h = 1; // assume saturated if undefined
        }
    } else {
        OOFEM_ERROR("undefined state vector.");
    }

    double ksi = status->giveHydrationDegree();
    if ( !ksi && initialHydrationDegree ) {
        ksi = initialHydrationDegree;
        status->setHydrationDegree(ksi);
    }

    // !!! timeScale
    if ( tStep->giveTimeIncrement() > 0. ) {
        dksi = computeHydrationDegreeIncrement(ksi, T, h, tStep->giveTimeIncrement() * timeScale);
    } else {
        dksi = 0.;
    }

    status->setTempHydrationDegree(ksi + dksi);
}


// === Auxiliary functions === (root finding)
double
HydrationModel :: regulafindroot() const
{
    double x0, y0, yl, xl = 0., xr = 1.;

    do {
        yl = localResidual(xl);
        x0 = yl * ( xl - xr ) / ( localResidual(xr) - yl ) + xl;
        y0 = localResidual(x0);
        if ( y0 < 0 ) {
            xl = x0;
        } else {
            xr = x0;
        }

#ifdef VERBOSEFINDROOT
        OOFEM_LOG_INFO("regulafindroot: x=%.15f, chyba %.15f \n", x0, y0);
#endif
    } while ( fabs(y0) > ROOT_PRECISION_DKSI );

    return x0;
}

double
HydrationModel :: bintreefindroot() const
{
    double xl = 0., xr = 1., x0, y0;

    do {
        x0 = ( xl + xr ) / 2;
        y0 = localResidual(x0);
        if ( y0 < 0 ) {
            xl = x0;
        } else {
            xr = x0;
        }

#ifdef VERBOSEFINDROOT
        OOFEM_LOG_INFO("bintreefindroot: x=%.15f, chyba %.15f \n", x0, y0);
#endif
    } while ( fabs(y0) > ROOT_PRECISION_DKSI );

    return x0;
}



double
HydrationModel :: mixedfindroot() const
{
    double x0 = 0., y0, xl = 0., xr = 1.;
    bool done = false;

    do {
        for ( int jcount = 0; ( jcount < BINARY_TREE_STEPS ); jcount++ ) {
            x0 = ( xl + xr ) / 2;
            y0 = localResidual(x0);
            if ( fabs(y0) < ROOT_PRECISION_DKSI ) {
                done = true;
                break;
            }

            if ( y0 < 0 ) {
                xl = x0;
            } else {
                xr = x0;
            }
        }

        if ( !done ) {
            double yl = localResidual(xl);
            x0 = yl * ( xl - xr ) / ( localResidual(xr) - yl ) + xl;
            y0 = localResidual(x0);

            if ( fabs(y0) < ROOT_PRECISION_DKSI ) {
                break;
            } else if ( y0 < 0 ) {
                xl = x0;
            } else {
                xr = x0;
            }

#ifdef VERBOSEFINDROOT
            OOFEM_LOG_INFO("mixedfindroot: x=%.15f, chyba %.15f \n", x0, y0);
#endif
        }
    } while ( !done );

    return x0;
}

MaterialStatus *
HydrationModel :: CreateStatus(GaussPoint *gp) const
{
    return new HydrationModelStatus(gp);
}

// ======= HydrationModelStatusInterface implementation =======
void
HydrationModelStatusInterface :: updateYourself(TimeStep *tStep)
{
    if ( hydrationModelStatus ) {
        hydrationModelStatus->updateYourself(tStep);
    }
}

void
HydrationModelStatusInterface :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    if ( hydrationModelStatus ) {
        hydrationModelStatus->printOutputAt(file, tStep);
    }
}

// ======= HydrationModelInterface implementation =======

void
HydrationModelInterface :: initializeFrom(InputRecord &ir)
{
    double value;

    // !!! should use separate field, e.g. hydramname #hydramnumber
    // Hydration>0  ->  Model starting at value, hydration<0 -> Constant at given value
    value = -2.;
    constantHydrationDegree = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_HydrationModelInterface_hydration);
    if ( value >= 0. ) {
        OOFEM_LOG_INFO("HydratingMaterial: creating HydrationModel.");
        hydrationModel = std::make_unique<HydrationModel>();
        if ( !hydrationModel ) {
            throw ValueInputException(ir, _IFT_HydrationModelInterface_hydration, "Could not create HydrationModel instance.");
        }

        hydrationModel->initializeFrom(ir);
    }
    // constant hydration degree
    else if ( value >= -1. ) {
        constantHydrationDegree = -value;
        OOFEM_LOG_INFO("HydratingMaterial: Hydration degree set to %.2f.", -value);
    } else {
        throw ValueInputException(ir, _IFT_HydrationModelInterface_hydration, "must be between 0 and 1");
    }

    // Material cast time - start of hydration
    // 11/3/2004 OK *unfinished in Hellmat, needs to be checked in hm_Interface->updateInternalState
    castAt = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, castAt, _IFT_HydrationModelInterface_castAt);
    if ( castAt >= 0. ) {
        OOFEM_LOG_INFO("HydratingMaterial: Hydration starts at time %.2g.", castAt);
    }
}

void
HydrationModelInterface :: updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *tStep)
{
    if ( hydrationModel ) {
        TimeStep hydraTime( ( const TimeStep ) *tStep );
        int notime = 0;
        if ( tStep->giveTargetTime() - tStep->giveTimeIncrement() < castAt ) {
            if ( tStep->giveTargetTime() >= castAt ) {
                hydraTime.setTimeIncrement(tStep->giveTargetTime() - castAt);
            } else {
                notime = 1;
            }
        }

        if ( !notime ) {
            hydrationModel->updateInternalState(vec, gp, &hydraTime);
        }
    }
}

double
HydrationModelInterface :: giveHydrationDegree(GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const
{
    if ( hydrationModel ) {
        return hydrationModel->giveHydrationDegree(gp, tStep, mode);
    } else {
        return constantHydrationDegree;
    }
}
} // end namespace oofem
