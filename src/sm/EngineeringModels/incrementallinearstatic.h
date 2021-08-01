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

#ifndef incrementallinearstatic_h
#define incrementallinearstatic_h

#include "sm/EngineeringModels/structengngmodel.h"
#include "sparselinsystemnm.h"
#include "sparsemtrxtype.h"

#include <memory>

///@name Input fields for IncrementalLinearStatic
//@{
#define _IFT_IncrementalLinearStatic_Name "incrlinearstatic"
#define _IFT_IncrementalLinearStatic_endoftimeofinterest "endoftimeofinterest"
#define _IFT_IncrementalLinearStatic_prescribedtimes "prescribedtimes"
#define _IFT_IncrementalLinearStatic_deltat "deltat"
//@}

namespace oofem {
/**
 * This class implements Incremental LinearStatic Engineering problem.
 * Problem is solved as series of explicit linear solutions.
 * It will solve any nonlinear problem by writing the expressions explicitly
 * @f$ K({}^{n-1} r)\cdot \delta r = R({}^{n-1}r) - F @f$.
 * Explicit solver has requirements on the time step size to obtain a stable solution.
 *
 * This class can be used for solving linear creep problems.
 *
 * Supports the changes of static scheme (applying, removing and changing  boundary conditions)
 * during the analysis.
 *
 */
class IncrementalLinearStatic : public StructuralEngngModel
{
protected:
    std :: unique_ptr< SparseMtrx > stiffnessMatrix;
    FloatArray loadVector;
    FloatArray internalLoadVector;
    FloatArray incrementOfDisplacementVector;
    FloatArray discreteTimes;
    bool fixedSteps;
    double deltaT;
    double endOfTimeOfInterest;
    /// Numerical method used to solve the problem
    std :: unique_ptr< SparseLinearSystemNM > nMethod;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;

public:
    IncrementalLinearStatic(int i, EngngModel *master = nullptr);
    virtual ~IncrementalLinearStatic();

    void initializeFrom(InputRecord &ir) override;

    void solveYourself() override;
    void solveYourselfAt(TimeStep *tStep) override;

    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof) override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
    TimeStep *giveNextStep() override;

    /**
     * This function returns time valid for iStep time step, used in integration
     * of structure response.
     * This functions, when invoked for the first time, generates table of times.
     * Times in this table are generated according to:
     * - if there exists some load time function with abrupt change,
     *   then time just before and just at abrupt is included in table.
     * - between these steps under constant loads (or when load changes continuously)
     *   we use progressively increasing time step. They are best chosen so that time
     *   step be kept constant in the log (t-t') scale
     * @param iStep Time step number.
     */
    double giveDiscreteTime(int iStep);
    double giveEndOfTimeOfInterest() override { return endOfTimeOfInterest; }

    NumericalMethod *giveNumericalMethod(MetaStep *mStep) override;

    fMode giveFormulation() override { return TL; }

    const char *giveInputRecordName() const { return _IFT_IncrementalLinearStatic_Name; }
    const char *giveClassName() const override { return "IncrementalLinearStatic"; }

    int requiresUnknownsDictionaryUpdate() override { return true; }
    bool requiresEquationRenumbering(TimeStep *tStep) override { return true; }
    void updateDofUnknownsDictionary(DofManager *dman, TimeStep *tStep) override;
    // Here we store only total and incremental value; so hash is computed from mode value only
    int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // incrementallinearstatic_h
