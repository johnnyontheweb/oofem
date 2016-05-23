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

#ifndef responsespectrum_h
#define responsespectrum_h

#include "engngm.h"
#include "sparsegeneigenvalsystemnm.h"
#include "sparsemtrx.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "geneigvalsolvertype.h"
#include "linearstatic.h"
#include "integrationpointstatus.h"
#include "integrationrule.h"
#include "gausspoint.h"

///@name Input fields for ResponseSpectrum
//@{
#define _IFT_ResponseSpectrum_Name "responsespectrum"
#define _IFT_ResponseSpectrum_nroot "nroot"
#define _IFT_ResponseSpectrum_rtolv "rtolv"
#define _IFT_ResponseSpectrum_stype "stype"
#define _IFT_ResponseSpectrum_func "func"
#define _IFT_ResponseSpectrum_dir "dir"
//@}

namespace oofem {
/**
 * This class implements way for performing a response spectrum analysis.
 * Base for this class is EigenValueDynamic class, functionality
 * of which has been extended.
 *
 * @author Francesco Pontarin
 */
class ResponseSpectrum : public EngngModel
{
private:
    std :: unique_ptr< SparseMtrx > stiffnessMatrix;
    std :: unique_ptr< SparseMtrx > massMatrix;
    SparseMtrxType sparseMtrxType;
    FloatMatrix eigVec;
    FloatArray eigVal;
	FloatArray periods;
	FloatArray centroid;
	FloatArray totMass;
	FloatMatrix partFact;
	FloatMatrix massPart;
    int numberOfRequiredEigenValues;
    int activeVector;
    int restoreFlag;
    /// Relative tolerance.
    double rtolv;
    /// Numerical method used to solve the problem.
    std :: unique_ptr< SparseGeneralEigenValueSystemNM > nMethod;
    GenEigvalSolverType solverType;
	int func;
	int dir;
	FloatArray loadVector;
	FloatArray dummyDisps;
	FloatArray combReactions;
	FloatArray combDisps;
	IntArray dofManMap, dofidMap, eqnMap;


public:
	ResponseSpectrum(int i, EngngModel * _master = NULL) : EngngModel(i, _master)
    {
        numberOfSteps = 1;
        ndomains = 1;
    }
	virtual ~ResponseSpectrum() { }

    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void terminate(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);

    virtual double giveUnknownComponent(ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);
    virtual void setActiveVector(int i) { activeVector = i; }
    virtual int resolveCorrespondingEigenStepNumber(void *obj);

    virtual double giveEigenValue(int eigNum) { return eigVal.at(eigNum); }

	virtual void postInitialize();

	virtual void ResponseSpectrum::getGPOutputAt(GaussPoint *gp, TimeStep *tStep);
	virtual void ResponseSpectrum::getIntRuleOutputAt(IntegrationRule *iRule, TimeStep *tStep);
	virtual void ResponseSpectrum::getIntPointStatusOutputAt(IntegrationPointStatus *iStatus, TimeStep *tStep, MaterialMode materialMode);
	virtual double calcSpectrumOrdinate(double period);
	virtual void computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di);
	virtual void buildReactionTable(IntArray &restrDofMans, IntArray &restrDofs,
		IntArray &eqn, TimeStep *tStep, int di);
	virtual void computeReaction(FloatArray &answer, TimeStep *tStep, int di);
	virtual void updateInternalState(TimeStep *tStep);

    // identification
    virtual const char *giveClassName() const { return "ResponseSpectrum"; }
};
} // end namespace oofem
#endif // responsespectrum_h
