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

#ifndef mitc4_h
#define mitc4_h

#include "sm/Elements/nlstructuralelement.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "spatiallocalizer.h"
#include "load.h"
<<<<<<< HEAD
=======
#include "floatmatrixf.h"
>>>>>>> bp2/master
//#include "eleminterpmapperinterface.h"//

#define _IFT_MITC4Shell_Name "mitc4shell"
#define _IFT_MITC4Shell_nipZ "nipz"
<<<<<<< HEAD
#define _IFT_MITC4Shell_FirstLocalAxis "lcs1"
#define _IFT_MITC4Shell_macroElem "macroelem"
=======
#define _IFT_MITC4Shell_directorType "directortype"
>>>>>>> bp2/master

namespace oofem {
	class FEI2dQuadLin;
#ifndef __CHARTENSOR
#define __CHARTENSOR
	enum CharTensor {
		LocalStrainTensor,
		GlobalStrainTensor,
		LocalForceTensor,
		GlobalForceTensor,
	};
#endif

<<<<<<< HEAD
	/**
	* This class implements an quad element based on Mixed Interpolation of Tensorial Components (MITC).
	* This element is a shell element suitable for both thin and thick shells.
	* The element has 24 DOFs (u,v,w-displacements and three rotations) in each node
	*
	* Tasks:
	* - calculating its B,D matrices and dV.
	*/
	class MITC4Shell : public NLStructuralElement, public ZZNodalRecoveryModelInterface,
		public SPRNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface,
		public SpatialLocalizerInterface
	{
	protected:
		/// Element geometry approximation
		static FEI2dQuadLin interp_lin;
		/// Ordering for the normal shell stiffness (everything but the out-of-plane rotations)
		static IntArray shellOrdering;
		/// Ordering for the drilling dofs (the out-of-plane rotations)
		static IntArray drillOrdering;

		/**
		* Transformation Matrix form GtoL(3,3) is stored
		* at the element level for computation efficiency.
		*/
		FloatMatrix GtoLRotationMatrix;
		int nPointsXY, nPointsZ;
		// macro element number
		int macroElem;
		double drillCoeff;
		double ang;
	private:
		//std::auto_ptr<FloatMatrix> NMatrix;
		std::vector<std::auto_ptr<FloatMatrix>> BMatrices;
		std::vector<std::auto_ptr<FloatMatrix>> LMatrices;

	public:

		MITC4Shell(int n, Domain *d);
		virtual ~MITC4Shell() { }

		FloatArray la1;
		
		virtual FEInterpolation *giveInterpolation() const;
		virtual FEInterpolation *giveInterpolation(DofIDItem id) const;
		virtual int testElementExtension(ElementExtension ext) { return (((ext == Element_EdgeLoadSupport) || (ext == Element_SurfaceLoadSupport)) ? 1 : 0); }

		virtual Interface *giveInterface(InterfaceType interface);

		virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
		virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
		virtual int SPRNodalRecoveryMI_giveNumberOfIP();
		virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();
		virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);

		virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
		void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override;

		// transformation
    bool computeGtoLRotationMatrix(FloatMatrix &answer);
    int computeLoadGToLRotationMtrx(FloatMatrix &answer);
		void computeLToDirectorRotationMatrix(FloatMatrix &answer1, FloatMatrix &answer2, FloatMatrix &answer3, FloatMatrix &answer4);
		virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);

		//virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
		virtual double computeArea();

	protected:

		virtual void computeGaussPoints();
		virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
		virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
		virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
		virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
		virtual void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep);


	private:
		void giveNodeCoordinates(double &x1, double &x2, double &x3, double &x4,
			double &y1, double &y2, double &y3, double &y4,
			double &z1, double &z2, double &z3, double &z4);
		void giveDirectorVectors(FloatArray &V1, FloatArray &V2, FloatArray &V3, FloatArray &V4);
		void giveLocalDirectorVectors(FloatArray &V1, FloatArray &V2, FloatArray &V3, FloatArray &V4);
		void giveThickness(double &a1, double &a2, double &a3, double &a4);
    void giveJacobian(FloatArray lcoords, FloatMatrix &jacobianMatrix);
		void giveLocalCoordinates(FloatArray &answer, FloatArray &global);
		const FloatMatrix *computeGtoLRotationMatrix();
		virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
		void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);
		virtual void printOutputAt(FILE *file, TimeStep *tStep);
		virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
		virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &coords);
		virtual double computeVolumeAround(GaussPoint *gp);
		void computeLocalBaseVectors(FloatArray &e1, FloatArray &e2, FloatArray &e3);
    void givedNdx(FloatArray &hkx, FloatArray &hky, FloatArray coords);

    void giveMidplaneIPValue(FloatArray &answer, int gpXY, InternalStateType type, TimeStep *tStep);

		// definition & identification
		virtual const char *giveClassName() const { return "MITC4Shell"; }
		virtual const char *giveInputRecordName() const { return _IFT_MITC4Shell_Name; }
		virtual IRResultType initializeFrom(InputRecord *ir);
		virtual int computeNumberOfDofs() { return 24; }
		virtual int computeNumberOfGlobalDofs() { return 24; }
    virtual integrationDomain giveIntegrationDomain() const { return _3dDegShell; }
		virtual MaterialMode giveMaterialMode() { return _3dDegeneratedShell; }


		// edge & body load
		virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
		virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
		virtual void giveDofManDofIDMask(int inode, IntArray &) const;
		virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf);
		virtual IntegrationRule *GetSurfaceIntegrationRule(int approxOrder);
		virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp);
		virtual void giveSurfaceDofMapping(IntArray &answer, int iSurf) const;

    virtual void computeSurfaceNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords);
    virtual void computeEdgeNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords);
	virtual void setupIRForMassMtrxIntegration(IntegrationRule &iRule);
 
	};
=======
/**
 * This class implements an quad element based on Mixed Interpolation of Tensorial Components (MITC).
 * This element is a shell element suitable for both thin and thick shells.
 * The element has 24 DOFs (u,v,w-displacements and three rotations) in each node
 *
 * Tasks:
 * - calculating its B,D matrices and dV.
 */

class MITC4Shell : public NLStructuralElement, public ZZNodalRecoveryModelInterface,
    public SPRNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface,
    public SpatialLocalizerInterface
{
protected:
    /// Element geometry approximation
    static FEI2dQuadLin interp_lin;
    /**
     * Transformation Matrix form GtoL(3,3) is stored
     * at the element level for computation efficiency.
     */
    FloatMatrixF<3,3> GtoLRotationMatrix;
    int nPointsXY = 0, nPointsZ = 0, directorType = 0;
    double drillCoeff = 0.;

public:
    MITC4Shell(int n, Domain *d);

    FEInterpolation *giveInterpolation() const override;
    FEInterpolation *giveInterpolation(DofIDItem id) const override;
    int testElementExtension(ElementExtension ext) override { return ( ( ( ext == Element_EdgeLoadSupport ) || ( ext == Element_SurfaceLoadSupport ) ) ? 1 : 0 ); }

    Interface *giveInterface(InterfaceType interface) override;
    // definition & identification
    const char *giveClassName() const override { return "MITC4Shell"; }
    const char *giveInputRecordName() const override { return _IFT_MITC4Shell_Name; }
    integrationDomain giveIntegrationDomain() const override { return _3dDegShell; }
    MaterialMode giveMaterialMode() override { return _3dDegeneratedShell; }
    void initializeFrom(InputRecord &ir) override;
    void postInitialize() override;
    int computeNumberOfDofs() override { return 24; }
    int computeNumberOfGlobalDofs() override { return 24; }

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override;
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override;
    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override;

    // transformation
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
    int computeLoadGToLRotationMtrx(FloatMatrix &answer) override;
    std::array<FloatMatrixF<3,3>, 4> computeLToDirectorRotationMatrix();
    int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp) override;

protected:
    void computeGaussPoints() override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;


private:
    std::array<FloatArrayF<3>, 4> giveNodeCoordinates();
    std::array<FloatArrayF<3>, 4> giveDirectorVectors();
    std::array<FloatArrayF<3>, 4> giveLocalDirectorVectors();
    std::array<double, 4> giveThickness();
    FloatMatrixF<3,3> giveJacobian(const FloatArrayF<3> &lcoords);
    FloatArrayF<3> giveLocalCoordinates(const FloatArrayF<3> &global);
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    FloatMatrix giveCharacteristicTensor(CharTensor type, GaussPoint *gp, TimeStep *tStep);
    void printOutputAt(FILE *file, TimeStep *tStep) override;
    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;
    bool computeLocalCoordinates(FloatArray &answer, const FloatArray &coords) override;
    double computeVolumeAround(GaussPoint *gp) override;
    std::array<FloatArrayF<3>, 3> computeLocalBaseVectors();
    std::array<FloatArrayF<4>, 2> givedNdx(const FloatArrayF<3> &coords);

    FloatArray giveMidplaneIPValue(int gpXY, InternalStateType type, TimeStep *tStep);

    // edge & body load
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) override;
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const override;
    void giveDofManDofIDMask(int inode, IntArray &) const override;
    double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf) override;
    void computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp);  ///@note This method isn't overloaded, it's never called for any class 
    void giveSurfaceDofMapping(IntArray &answer, int iSurf) const override;

    void computeSurfaceNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords) override;
    void computeEdgeNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords) override;
    void setupIRForMassMtrxIntegration(IntegrationRule &iRule) override;
};
>>>>>>> bp2/master
} // end namespace oofem
#endif // mitc4_h
