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

#ifndef AbaqusUserElement6d_h
#define AbaqusUserElement6d_h

#include "structuralelement.h"
#include "nlstructuralelement.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "timestep.h"

///@name Input fields for AbaqusUserElement6d
//@{
#define _IFT_AbaqusUserElement6d_Name "AbaqusUserElement6d"
#define _IFT_AbaqusUserElement6d_userElement "uel"
// #define _IFT_AbaqusUserElement6d_numcoords "coords"
//#define _IFT_AbaqusUserElement6d_dofs "dofs"
#define _IFT_AbaqusUserElement6d_refangle "refangle"
#define _IFT_AbaqusUserElement6d_numsvars "numsvars"
#define _IFT_AbaqusUserElement6d_properties "properties"
#define _IFT_AbaqusUserElement6d_iproperties "iproperties"
#define _IFT_AbaqusUserElement6d_type "type"
#define _IFT_AbaqusUserElement6d_name "name"
#define _IFT_AbaqusUserElement6d_macroElem "macroelem"
#define _IFT_AbaqusUserElement6d_orientation "dir"
#define _IFT_AbaqusUserElement6d_pos "pos"
//@}

namespace oofem {

/**
 * UEL interface from Abaqus user elements.
 * 
 * The function prototype for UEL is:
 * SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
 *      PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,
 *      DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
 *      PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,
 *      JPROPS,NJPROP,PERIOD)
 *
 * DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
 *      SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
 *      U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
 *      PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
 *      DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
 *      JPROPS(*)
 *
 * @author Giovanni
 */
class AbaqusUserElement6d : public NLStructuralElement
{
private:
    /// Dynamically loaded uel
    void *uelobj;

    /// Coord transf
    int nCoords;
    IntArray dofs;
    FloatMatrix coords;

    /// Number of status variables
    int numSvars;

    /// Element properties
    FloatArray props;

    /// Status variables
    FloatArray svars, tempSvars;

    /// Element type
    int jtype;

    /// Element rhs
    int nrhs = 2;
    FloatMatrix rhs, tempRHS;

    /// Element amatrx
    FloatMatrix amatrx, tempAmatrx;

    /// ndofel
    int ndofel = 0;

    /// mcrd
    int mcrd = 0;

    /// mlvarx
    int mlvarx = 1;

    /// Inputs to element routines. Velocity and Acceleration currently ignored.
    FloatArray U, V, A;
    FloatMatrix DU;

    /// LFlags
    IntArray lFlags;

    /// kinc, kstep
    int kinc = 1, kstep = 1;

    /// predef
    int npredef = 1;
    FloatArray predef;

    /// energy
    FloatArray energy;

    /// loads
    int ndLoad = 0;
    int mdLoad = 0;

    FloatMatrix adlmag, ddlmag;
    int *jdltype = NULL;                // Temporary init.

    /// params
    double params [ 3 ];

    /// jprops
    IntArray jprops;

    /// Keeps track of whether the tangent has been obtained already.
    bool hasTangentFlag;

    /// Pointer to the dynamically loaded uel-function (translated to C)
    void (*uel)(double *rhs, double *amatrx, double *svars, double energy [ 8 ], int *ndofel,            // 5
                int *nrhs, int *nsvars, double *props, int *nprops, double *coords, int *mcrd,           // 6
                int *nnode, double *u, double *du, double *v, double *a, int *jtype,                     // 6
                double time [ 2 ], double *dtime, int *kstep, int *kinc, int *jelem,                     // 5
                double params [ 3 ], int *ndload, int *jdltyp, double *adlmag, double *predef,           // 5
                int *npredef, int *lflags, int *mvarx, double *ddlmag, int *mdload,                      // 5
                double *pnewdt, int *jprops, int *njprop, double *period);                               // 4 - tot 36

    /// File containing the uel function
    std :: string filename;

	/// Pointer to the dynamically loaded uelInit function (translated to C)
	void(*uelInit)();

	// macro element number
	int macroElem;

public:
    /// Constructor
    AbaqusUserElement6d(int n, Domain *d);

    /// Destructor
    virtual ~AbaqusUserElement6d();
	virtual double computeLength();
    virtual void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL);
    //virtual void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, FloatArray &U, FloatMatrix &DU, int useUpdatedGpRecord);
    virtual int computeNumberOfDofs() { return 12; }
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual void computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
    { OOFEM_ERROR("Abaqus user element cannot support computation of local unknown vector\n"); }
    virtual void updateYourself(TimeStep *tStep);
    virtual void updateInternalState(TimeStep *tStep);

	virtual int computeNumberOfGlobalDofs() { return 12; }
	virtual void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep);
    bool hasTangent() const {
        return hasTangentFlag;
    }
    virtual void letTempTangentBe(FloatMatrix &src) {
        tempAmatrx = src;
        hasTangentFlag = true;
    }
    virtual void letTempRhsBe(FloatMatrix &src) {
		tempRHS = src;
    }
    virtual void letTempSvarsBe(FloatArray &src) {
		tempSvars = src;
    }
    virtual const FloatArray &giveStateVector() const {
        return svars;
    }
    virtual const FloatArray &giveTempStateVector() const {
        return tempSvars;
    }
    virtual const FloatMatrix &giveTempTangent() {
        return tempAmatrx;
    }

    virtual Interface *giveInterface(InterfaceType it);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
	virtual void postInitialize();

	virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition & identification
    virtual const char *giveClassName() const { return "AbaqusUserElement6d"; }
    virtual const char *giveInputRecordName() const { return _IFT_AbaqusUserElement6d_Name; }
    virtual integrationDomain giveIntegrationDomain() const {
        // return _Unknown_integrationDomain;
		return _UnknownIntegrationDomain;
    }
    virtual Element_Geometry_Type giveGeometryType() const {
        return EGT_line_1;
    }

protected:
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) {
        OOFEM_ERROR("function not defined for AbaqusUserElement6d and should never be called.");
    }

    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) {
        OOFEM_ERROR("function not defined for AbaqusUserElement6d and should never be called.");
    }

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) {
        OOFEM_ERROR("function not defined for AbaqusUserElement6d and should never be called.");
        return 0;
    }
	virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
	/**
	* Orientation vector. Defines orientation of spring element- for spring it defines the direction of spring,
	* for torsional spring it defines the axis of rotation.
	*/
	FloatArray dir;
	double referenceAngle = 0;
	double pos;
};

}// namespace oofem

#endif  // AbaqusUserElement6d_h
