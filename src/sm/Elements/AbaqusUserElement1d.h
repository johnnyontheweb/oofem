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

#ifndef AbaqusUserElement1d_h
#define AbaqusUserElement1d_h

#include "structuralelement.h"
// #include "nlstructuralelement.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "timestep.h"

///@name Input fields for AbaqusUserElement1d
//@{
#define _IFT_AbaqusUserElement1d_Name "AbaqusUserElement1d"
#define _IFT_AbaqusUserElement1d_userElement "uel"
// #define _IFT_AbaqusUserElement1d_numcoords "coords"
//#define _IFT_AbaqusUserElement1d_dofs "dofs"
#define _IFT_AbaqusUserElement1d_numsvars "numsvars"
#define _IFT_AbaqusUserElement1d_properties "properties"
#define _IFT_AbaqusUserElement1d_iproperties "iproperties"
#define _IFT_AbaqusUserElement1d_type "type"
#define _IFT_AbaqusUserElement1d_name "name"
#define _IFT_AbaqusUserElement1d_macroElem "macroelem"
#define _IFT_AbaqusUserElement1d_orientation "dir"
#define _IFT_AbaqusUserElement1d_mode "mode"
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
class AbaqusUserElement1d : public StructuralElement
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
    AbaqusUserElement1d(int n, Domain *d);

    /// Destructor
    virtual ~AbaqusUserElement1d();

    void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL) override;
    //virtual void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, FloatArray &U, FloatMatrix &DU, int useUpdatedGpRecord);
    int computeNumberOfDofs() override { return 2; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    void computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer) override
    { OOFEM_ERROR("Abaqus user element cannot support computation of local unknown vector\n"); }
    void updateYourself(TimeStep *tStep) override;
    void updateInternalState(TimeStep *tStep) override;
	void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep) override;
	int computeNumberOfGlobalDofs() override { return 6; }

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

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    void postInitialize() override;

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    // definition & identification
    const char *giveClassName() const override { return "AbaqusUserElement1d"; }
    const char *giveInputRecordName() const override { return _IFT_AbaqusUserElement1d_Name; }
    integrationDomain giveIntegrationDomain() const override {
        // return _Unknown_integrationDomain;
	return _UnknownIntegrationDomain;
    }
    Element_Geometry_Type giveGeometryType() const override {
        return EGT_line_1;
    }

protected:
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override {
        OOFEM_ERROR("function not defined for AbaqusUserElement1d and should never be called.");
    }

    void computeConstitutiveMatrixAt(FloatMatrix& answer, MatResponseMode rMode, GaussPoint* gp, TimeStep* tStep) override {
        OOFEM_ERROR("calling of this function is not allowed");
    }

    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override {
        OOFEM_ERROR("function not defined for AbaqusUserElement1d and should never be called.");
    }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override {
        OOFEM_ERROR("function not defined for AbaqusUserElement1d and should never be called.");
        return 0;
    }
	bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
	/**
	* Orientation vector. Defines orientation of spring element- for spring it defines the direction of spring,
	* for torsional spring it defines the axis of rotation.
	*/
	FloatArray dir;
	// trasl. 4, rot. 5
	int mode;
};

}// namespace oofem

#endif  // AbaqusUserElement1d_h
