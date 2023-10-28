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
#include "eigenvectorprimaryfield.h"


///@name Input fields for ResponseSpectrum
//@{
#define _IFT_ResponseSpectrum_Name "responsespectrum"
#define _IFT_ResponseSpectrum_nroot "nroot"
#define _IFT_ResponseSpectrum_rtolv "rtolv"
#define _IFT_ResponseSpectrum_stype "stype"
#define _IFT_ResponseSpectrum_func "func"
#define _IFT_ResponseSpectrum_dir "dir"
#define _IFT_ResponseSpectrum_modalCombo "modalcombo"
#define _IFT_ResponseSpectrum_damp "damp"
//@}

namespace oofem {

enum RSpecComboType {
    RSC_CQC,
    RSC_SRSS
};

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
    std ::unique_ptr<EigenVectorPrimaryField> field;
    std ::unique_ptr<SparseMtrx> stiffnessMatrix;
    std ::unique_ptr<SparseMtrx> massMatrix;
    SparseMtrxType sparseMtrxType;
    FloatMatrix eigVec;
    FloatArray eigVal;
    FloatArray periods;
    FloatMatrix rhos;
    FloatArray centroid;
    FloatArray totMass;
    FloatMatrix partFact;
    FloatMatrix massPart;
    int dominantMode = 1;
    int numberOfRequiredEigenValues;
    int activeVector;
    int restoreFlag;
    /// Relative tolerance.
    double rtolv;
    /// Numerical method used to solve the problem.
    std ::unique_ptr<SparseGeneralEigenValueSystemNM> nMethod;
    GenEigvalSolverType solverType;
    /// Numerical method used to solve the static problem.
    LinSystSolverType linStype = ST_Direct;
    RSpecComboType modalCombo;
    double csi;
    int func;
    FloatArray dir;
    FloatArray loadVector;
    FloatArray dummyDisps;
    std::list<FloatArray> reactionsList;
    std::list<FloatArray> dispList;
    FloatArray combReactions;
    FloatArray combDisps;
    IntArray dofManMap, dofidMap, eqnMap;
    std::list<std::map<int, std::map<int, std::map<int, std::map<std::string, FloatArray> > > > > elemResponseList;
    std::list<std::map<int, std::map<std::string, FloatArray> > > beamResponseList;
    std::map<int, std::map<int, std::map<int, std::map<std::string, FloatArray> > > > combElemResponse;
    std::map<int, std::map<std::string, FloatArray> > combBeamResponse;


public:
    ResponseSpectrum( int i, EngngModel *_master = nullptr ) :
        EngngModel( i, _master )
    {
        numberOfSteps = 1;
        ndomains      = 1;
    }
    virtual ~ResponseSpectrum() {}

    void solveYourself() override;
    void terminate( TimeStep *tStep ) override;
    void updateYourself( TimeStep *tStep ) override;

    int giveUnknownDictHashIndx( ValueModeType mode, TimeStep *tStep ) override;
    double giveUnknownComponent( ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof ) override;
    void initializeFrom( InputRecord &ir ) override;
    void saveContext( DataStream &stream, ContextMode mode ) override;
    void restoreContext( DataStream &stream, ContextMode mode ) override;
    TimeStep *giveNextStep() override;
    NumericalMethod *giveNumericalMethod( MetaStep *mStep ) override;
    void setActiveVector( int i ) override { activeVector = i; }

    double giveEigenValue( int eigNum ) override { return eigVal.at( eigNum ); }

    void postInitialize() override;

    void getGPOutputAt( GaussPoint *gp, TimeStep *tStep, std::map<std::string, FloatArray> *&ips );
    void getIntRuleOutputAt( IntegrationRule *iRule, TimeStep *tStep, std::map<int, std::map<std::string, FloatArray> > *&ir );
    void getIntPointStatusOutputAt( IntegrationPointStatus *iStatus, TimeStep *tStep, MaterialMode materialMode, std::map<std::string, FloatArray> *&ir );
    double calcSpectrumOrdinate( double period );
    void computeExternalLoadReactionContribution( FloatArray &reactions, TimeStep *tStep, int di );
    void buildReactionTable( IntArray &restrDofMans, IntArray &restrDofs,
        IntArray &eqn, TimeStep *tStep, int di );
    void computeReaction( FloatArray &answer, TimeStep *tStep, int di );
    void updateInternalState( TimeStep *tStep );
    void SRSS();
    void CQC();
    void giveRhos( FloatMatrix &rhos );
    void giveDominantMode( int &mode );
    RSpecComboType giveComboType();

    // identification
    const char *giveClassName() const override { return "ResponseSpectrum"; }
};
} // end namespace oofem
#endif // responsespectrum_h
