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

#ifndef tutorialmaterial_h
#define tutorialmaterial_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"

#include <memory>

///@name Input fields for StructuralFE2Material
//@{
#define _IFT_StructuralFE2Material_Name "structfe2material"
#define _IFT_StructuralFE2Material_fileName "filename"
//@}

namespace oofem {
class EngngModel;
class PrescribedGradientHomogenization;

class StructuralFE2MaterialStatus : public StructuralMaterialStatus
{
protected:
    /// The RVE
    std :: unique_ptr< EngngModel > rve;
    /// Boundary condition in RVE that performs the computational homogenization.
    PrescribedGradientHomogenization *bc = nullptr;

    FloatMatrix tangent;
    bool oldTangent = true;

<<<<<<< HEAD
=======
    /// Interface normal direction
    FloatArray mNormalDir;

    std :: string mInputFile;

>>>>>>> bp2/master
public:
    StructuralFE2MaterialStatus(int rank, GaussPoint * g,  const std :: string & inputfile);

<<<<<<< HEAD
    EngngModel *giveRVE() { return this->rve.get(); }
    PrescribedGradientHomogenization *giveBC() { return this->bc; }
=======
    EngngModel *giveRVE() const { return this->rve.get(); }
    PrescribedGradientHomogenization *giveBC();// { return this->bc; }
>>>>>>> bp2/master

    void markOldTangent();
    void computeTangent(TimeStep *tStep);

    /// Creates/Initiates the RVE problem.
    bool createRVE(int n, const std :: string &inputfile, int rank);

    /// Copies time step data to RVE.
    void setTimeStep(TimeStep *tStep);

    FloatMatrix &giveTangent() { return tangent; }

    const char *giveClassName() const override { return "StructuralFE2MaterialStatus"; }

    void initTempStatus() override;

<<<<<<< HEAD
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
=======
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const FloatArray &giveNormal() const { return mNormalDir; }
    void letNormalBe(FloatArray iN) { mNormalDir = std :: move(iN); }

    double giveRveLength();

    /// Functions for MaterialStatusMapperInterface
    void copyStateVariables(const MaterialStatus &iStatus) override;
    void addStateVariables(const MaterialStatus &iStatus) override { OOFEM_ERROR("Not implemented."); }

    // For debugging only
    bool mNewlyInitialized = true;
>>>>>>> bp2/master
};


/**
 * Multiscale constitutive model for subscale structural problems.
 *
 * The material uses the PrescribedGradient boundary condition to perform computational homogenization.
 * The requirement for the supplied subscale problem is:
 * - It must have a PrescribedGradient boundary condition.
 * - It must be the first boundary condition
 *
 * @author Mikael Öhman 
 */
class StructuralFE2Material : public StructuralMaterial
{
protected:
    std :: string inputfile;
    static int n;
<<<<<<< HEAD

public:
    StructuralFE2Material(int n, Domain * d);
    virtual ~StructuralFE2Material();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual const char *giveInputRecordName() const { return _IFT_StructuralFE2Material_Name; }
    virtual const char *giveClassName() const { return "StructuralFE2Material"; }
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return true; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
    const void giveDeviatoricProjectionMatrix(FloatMatrix &answer);
    // stress computation methods
    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);
    
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
=======
    bool useNumTangent = false;

public:
    StructuralFE2Material(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    const char *giveInputRecordName() const override { return _IFT_StructuralFE2Material_Name; }
    const char *giveClassName() const override { return "StructuralFE2Material"; }
    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return true; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF< 3 > &strain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<4,4> givePlaneStrainStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;
>>>>>>> bp2/master
};

} // end namespace oofem
#endif // structuralfe2material
