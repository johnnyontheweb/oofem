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

#ifndef tr_shell01_h
#define tr_shell01_h

#include "../sm/Elements/structuralelement.h"
#include "zznodalrecoverymodel.h"
#include "../sm/ErrorEstimators/zzerrorestimator.h"
#include "../sm/Elements/Shells/cct3d.h"
#include "../sm/Elements/PlaneStress/trplanrot3d.h"
#include "spatiallocalizer.h"

#include <memory>

#define _IFT_TR_SHELL01_Name "tr_shell01"

///@name Additional Input fields for TR_SHELL01
//@{
#define _IFT_TR_SHELL01_FirstLocalAxis "lcs1"
//@}

namespace oofem {
/**
 * This class implements an triangular three-node shell finite element, composed of
 * cct3d and trplanrot3d elements.
 * Each node has 6 degrees of freedom.
 *
 * @author Ladislav Svoboda
 */
class TR_SHELL01 : public StructuralElement, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface, public ZZErrorEstimatorInterface, public SpatialLocalizerInterface
{
protected:
    /// Pointer to plate element.
    std :: unique_ptr< CCTPlate3d > plate;
    /// Pointer to membrane (plane stress) element.
    std :: unique_ptr< TrPlaneStrRot3d > membrane;
    /**
     * Element integraton rule (plate and membrane parts have their own integration rules)
     * this one used to integrate element error and perhaps can be (re)used for other putrposes.
     * Created on demand.
     */
    std :: unique_ptr< IntegrationRule > compositeIR;

    static IntArray loc_plate;
    static IntArray loc_membrane;

	// 1st local axis
	FloatArray la1;

public:
    /// Constructor
    TR_SHELL01(int n, Domain * d);
    /// Destructor
    virtual ~TR_SHELL01() {}

    virtual FEInterpolation *giveInterpolation() const { return plate->giveInterpolation(); }

    virtual int computeNumberOfDofs() { return 18; }
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const
    { plate->giveDofManDofIDMask(inode, answer); }
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_TR_SHELL01_Name; }
    virtual const char *giveClassName() const { return "TR_SHELL01"; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep);
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep);
    virtual double computeVolumeAround(GaussPoint *gp);
    virtual bool giveRotationMatrix(FloatMatrix &answer);
	int computeLoadGToLRotationMtrx(FloatMatrix &answer);

    virtual void updateYourself(TimeStep *tStep);
    virtual void updateInternalState(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual void postInitialize();
    void updateLocalNumbering(EntityRenumberingFunctor &f);
    void setCrossSection(int csIndx);
#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
#endif
    // the membrane and plate irules are same (chacked in initializeFrom)
    virtual int giveDefaultIntegrationRule() const { return plate->giveDefaultIntegrationRule(); }
    virtual IntegrationRule *giveDefaultIntegrationRulePtr() { return plate->giveDefaultIntegrationRulePtr(); }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }
    virtual integrationDomain giveIntegrationDomain() const { return _Triangle; }
    virtual MaterialMode giveMaterialMode() { return _Unknown; }

    virtual Interface *giveInterface(InterfaceType it);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep);

    virtual IntegrationRule *ZZErrorEstimatorI_giveIntegrationRule();
    virtual void ZZErrorEstimatorI_computeLocalStress(FloatArray &answer, FloatArray &sig);

    // SpatialLocalizerI
    virtual void SpatialLocalizerI_giveBBox(FloatArray &bb0, FloatArray &bb1);


    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) {
        return this->plate->computeGlobalCoordinates(answer, lcoords);
    }

    virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords) {
        return this->plate->computeLocalCoordinates(answer, gcoords);
    }

protected:
    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS)
    { OOFEM_ERROR("calling of this function is not allowed"); }
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &)
    { OOFEM_ERROR("calling of this function is not allowed"); }

    /// @todo In time delete
protected:
    virtual void computeGaussPoints()
    {
        this->membrane->computeGaussPoints();
        this->plate->computeGaussPoints();
    }
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
    { OOFEM_ERROR("calling of this function is not allowed"); }
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode);
	virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);

public:
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
    { OOFEM_ERROR("calling of this function is not allowed"); }
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { OOFEM_ERROR("calling of this function is not allowed"); }
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
    { OOFEM_ERROR("calling of this function is not allowed"); }

private:
	void giveNodeCoordinates(double &x1, double &x2, double &x3,
		double &y1, double &y2, double &y3,
		double &z1, double &z2, double &z3);
	void giveLocalCoordinates(FloatArray &answer, FloatArray &global);
	const FloatMatrix *computeGtoLRotationMatrix();
	virtual void computeSurfaceNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords);
	virtual void computeEdgeNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords);
	virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
	virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf);
};
} // end namespace oofem
#endif
