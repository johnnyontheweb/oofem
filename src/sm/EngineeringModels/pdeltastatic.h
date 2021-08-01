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

#ifndef PDeltaStatic_h
#define PDeltaStatic_h

#include "../sm/EngineeringModels/structengngmodel.h"
#include "sparselinsystemnm.h"
#include "sparsemtrxtype.h"

#define _IFT_PDeltaStatic_Name "PDeltaStatic"
#define _IFT_PDeltaStatic_rtolv "rtolv"

namespace oofem {
class SparseMtrx;

/**
 * This class implements linear static engineering problem.
 * Multiple loading works only if linear elastic material (such as isoLE)  is used.
 * (Other non-linear materials keep load history, so such multiple loading
 * will cause that next step will be assumed as new load increment,
 * not the total new load). Because they always compute real stresses according
 * to reached strain state, they are not able to respond to linear analysis.
 *
 * Solution of this problem is series of loading cases, maintained as sequence of
 * time-steps. This solution is in form of linear equation system Ax=b
 * Tasks:
 * - Creating Numerical method for solving @f$ K\cdot x=b @f$.
 * - Interfacing Numerical method to Elements.
 * - Managing time steps.
 */
class PDeltaStatic : public StructuralEngngModel
{
protected:
    std :: unique_ptr< SparseMtrx > stiffnessMatrix;
	std :: unique_ptr< SparseMtrx > initialStressMatrix;
    FloatArray loadVector;
    FloatArray displacementVector;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;
    /// Numerical method used to solve the problem.
    std :: unique_ptr< SparseLinearSystemNM > nMethod;
	double rtolv;
    int initFlag;
	void updateAfterStatic(TimeStep *tStep, Domain *domain);

public:
    PDeltaStatic(int i, EngngModel * _master = NULL);
    virtual ~PDeltaStatic();

    void solveYourself() override;
    void solveYourselfAt(TimeStep *tStep) override;

    double giveUnknownComponent(ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof) override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    virtual void updateDomainLinks();

    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);
	virtual void terminate(TimeStep *tStep);
    void initializeFrom(InputRecord &ir);

    // identification
    virtual const char *giveInputRecordName() const { return _IFT_PDeltaStatic_Name; }
    const char *giveClassName() const override { return "PDeltaStatic"; }
    virtual fMode giveFormulation() { return TL; }

    virtual int estimateMaxPackSize(IntArray &commMap, DataStream &buff, int packUnpackType);
};
} // end namespace oofem
#endif // PDeltaStatic_h
