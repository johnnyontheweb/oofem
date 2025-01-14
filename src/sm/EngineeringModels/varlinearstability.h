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

#ifndef linearstabilityvar_h
#define linearstabilityvar_h

#include "sm/EngineeringModels/structengngmodel.h"
#include "geneigvalsolvertype.h"
#include "sparsegeneigenvalsystemnm.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "nummet.h"

///@name Input fields for LinearStabilityVar
//@{
#define _IFT_VarLinearStability_Name "varlinearstability"
#define _IFT_VarLinearStability_nroot "nroot"
#define _IFT_VarLinearStability_rtolv "rtolv"
#define _IFT_VarLinearStability_stype "stype"
//@}

namespace oofem {

class EigenVectorPrimaryField;
// class DofDistributedPrimaryField;

/**
 * This class implements way for examining critical load of structure.
 *
 * Solution of this problem is base on equation in the form of: @f$ K\cdot y=w (K_\sigma)y @f$.
 * Currently eigenvalue problem is solved using subspace iteration.
 * The linear static solution, determining normal forces is done in time = 0
 * for constant loads, and time = 1 for the multiplicative load.
 *
 * Tasks:
 * - Assembling the governing equation in the form.
 * - Creating Numerical method for @f$ K\cdot y=w(K_\sigma)y @f$.
 * - Interfacing Numerical method to Elements.
 */
class VarLinearStability : public StructuralEngngModel
{
private:
    std ::unique_ptr<SparseMtrx> stiffnessMatrix;
    std ::unique_ptr<SparseMtrx> initialStressMatrix;
    std ::unique_ptr<EigenVectorPrimaryField> field;
    // std :: unique_ptr< DofDistributedPrimaryField > initialDisplacements;
    FloatArray eigVal;
    SparseMtrxType sparseMtrxType;

    int numberOfRequiredEigenValues;
    double rtolv;
    /// Numerical method used to solve the problem.
    GenEigvalSolverType solverType;
    std ::unique_ptr<SparseGeneralEigenValueSystemNM> nMethod;
    /// Numerical method used to solve the static problem.
    LinSystSolverType linStype;
    std ::unique_ptr<SparseLinearSystemNM> nMethodLS;

public:
    VarLinearStability( int i, EngngModel *master = nullptr );
    virtual ~VarLinearStability() {}

    void solveYourself() override;
    void solveYourselfAt( TimeStep *tStep ) override;

    void doStepOutput( TimeStep *tStep ) override;
    void printOutputAt( FILE *file, TimeStep *tStep ) override;
    // virtual void printOutputAt(FILE *file, TimeStep *tStep, const IntArray &nodeSets, const IntArray &elementSets);
    void terminateLinStatic( TimeStep *tStep );
    int requiresNewLsh() { return 0; }
    void updateYourself( TimeStep *tStep ) override;

    // the intrinsic time of time step defines active eigen value and corresponding vector,
    // for which values can be requested using
    // giveUnknownComponent method.
    // When DisplacementVector is requested, then if time==0 linear elastic solution displacement are returned,
    // otherwise corresponding eigen vector is considered as displacement vector
    double giveUnknownComponent( ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof ) override;
    int giveUnknownDictHashIndx( ValueModeType mode, TimeStep *tStep ) override;
    bool newDofHandling() override { return true; }
    void initializeFrom( InputRecord &ir ) override;
    void saveContext( DataStream &stream, ContextMode mode ) override;
    void restoreContext( DataStream &stream, ContextMode mode ) override;
    TimeStep *giveNextStep() override;

    double giveEigenValue( int eigNum ) override { return eigVal.at( eigNum ); }
    void setActiveVector( int i ) override;

    NumericalMethod *giveNumericalMethod( MetaStep *mStep ) override;
    SparseLinearSystemNM *giveNumericalMethodForLinStaticProblem( TimeStep *tStep );

    // identification
    const char *giveInputRecordName() const { return _IFT_VarLinearStability_Name; }
    const char *giveClassName() const override { return "VarLinearStability"; }
    fMode giveFormulation() override { return TL; }
};
} // end namespace oofem
#endif // linearstabilityvar_h
