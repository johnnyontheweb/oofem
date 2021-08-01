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

#ifndef SpringElement3D_h
#define SpringElement3D_h

#include "sm/Elements/structuralelement.h"

///@name Input fields for spring element
//@{
#define _IFT_SpringElement3D_Name "spring3d"
// #define _IFT_SpringElement3D_mode "mode"
#define _IFT_SpringElement3D_orientation "orientation"
#define _IFT_SpringElement3D_springC1 "k1"
#define _IFT_SpringElement3D_springC2 "k2"
#define _IFT_SpringElement3D_springC3 "k3"
#define _IFT_SpringElement3D_springC4 "k4"
#define _IFT_SpringElement3D_springC5 "k5"
#define _IFT_SpringElement3D_springC6 "k6"
#define _IFT_SpringElement3D_mass "m"
#define _IFT_SpringElement3D_macroElem "macroelem"
#define _IFT_SpringElement3D_refangle "refangle"
#define _IFT_SpringElement3D_actAsRigidLink "rl"
//@}

namespace oofem {
/**
 * This class implements a simple spring element. Its purpose is to introduce
 * a longitudinal or torsional spring between two nodes. The spring element is defined
 * by its spring constant and orientation.
 */
class SpringElement3D : public StructuralElement
{
protected:
    /// The longitudinal spring constant [Force/Length], torsional spring constant [Force*Length/Radians].
    double springC1, springC2, springC3, springC4, springC5, springC6;
    /// total mass of the spring; to be distributed to nodes
    double mass;
    /**
     * Orientation vector. Defines orientation of spring element- for spring it defines the direction of spring,
     * for torsional spring it defines the axis of rotation.
     */
    FloatArray dir;
    // macro element number
    int macroElem;
    double referenceAngle = 0;
    // length of element for rigid link transport terms
    double d = 0;
public:
    SpringElement3D(int n, Domain * d);
    virtual ~SpringElement3D() { }
    virtual double computeLength();
    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { computeLumpedMassMatrix(answer, tStep); }
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    // { answer.clear(); }
    int giveLocalCoordinateSystem(FloatMatrix &answer) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;

    int computeNumberOfDofs() override { return 12; }
    int computeNumberOfGlobalDofs() override;

    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    void updateInternalState(TimeStep *tStep) override { }
    void updateYourself(TimeStep *tStep) override { }
    virtual int checkConsistency() { return 1; }
    void printOutputAt(FILE *file, TimeStep *tStep) override;
    virtual bool isCast(TimeStep *tStep) {return true;}
    
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_SpringElement3D_Name; }
    const char *giveClassName() const override { return "SpringElement3D"; }
    void initializeFrom(InputRecord &ir) override;
    void postInitialize() override;
    Element_Geometry_Type giveGeometryType() const override { return EGT_point; }

protected:
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override
    { answer.clear(); }

    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                                  int lowerIndx = 1, int upperIndx = ALL_STRAINS) override
    { answer.clear(); }
    void computeConstitutiveMatrixAt(FloatMatrix& answer, MatResponseMode rMode, GaussPoint* gp, TimeStep* tStep) override
    { OOFEM_ERROR("calling of this function is not allowed"); }
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override { answer.clear(); }
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
    FloatArray computeSpringInternalForce(TimeStep *tStep);
};
} // end namespace oofem
#endif // SpringElement3D_h
