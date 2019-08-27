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

#include "../sm/Elements/Bars/truss3d.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "fei3dlinelin.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"
#include "bctracker.h"

#include "bodyload.h"
#include "boundaryload.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(Truss3d);

FEI3dLineLin Truss3d :: interp;

Truss3d :: Truss3d(int n, Domain *aDomain) :
NLStructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this)
{
    numberOfDofMans = 2;
}


FEInterpolation *Truss3d :: giveInterpolation() const { return & interp; }


Interface *
Truss3d :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    }

    //OOFEM_LOG_INFO("Interface on Truss3d element not supported");
    return NULL;
}


void
Truss3d :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    answer.clear();
    OOFEM_WARNING("IP values will not be transferred to nodes. Use ZZNodalRecovery instead (parameter stype 1)");
}


void
Truss3d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    FloatMatrix dN;
    this->interp.evaldNdx( dN, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(1, 6);
    answer.at(1, 1) = dN.at(1, 1);
    answer.at(1, 2) = dN.at(1, 2);
    answer.at(1, 3) = dN.at(1, 3);
    answer.at(1, 4) = dN.at(2, 1);
    answer.at(1, 5) = dN.at(2, 2);
    answer.at(1, 6) = dN.at(2, 3);
}


void
Truss3d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], 1, this);
    }
}


double
Truss3d :: computeLength()
{
    return this->interp.giveLength( FEIElementGeometryWrapper(this) );
}


void
Truss3d :: computeInitialStressMatrix( FloatMatrix &answer, TimeStep *tStep )
{
	// computes initial stress matrix of receiver (or geometric stiffness matrix)

	FloatMatrix stiff;
	FloatArray endForces;

	double l = this->computeLength();
	double N;

	answer.resize(6, 6);
	answer.zero();

	answer.at(2, 2) = 1;
	answer.at(2, 5) = -1;

	answer.at(3, 3) = 1;
	answer.at(3, 6) = -1;

	answer.at(5, 2) = -1;
	answer.at(5, 5) = 1;

	answer.at(6, 3) = -1;
	answer.at(6, 6) = 1;

	GaussPoint* gp = integrationRulesArray[0]->getIntegrationPoint(0);
	double area = this->giveStructuralCrossSection()->give(CS_Area, gp);

	answer.symmetrized();
	// ask end forces in g.c.s
	this->giveEndForcesVector(endForces, tStep);

	N = (-endForces.at(1) + endForces.at(4)) / 2.;
	answer.times(N / l);
}


void
Truss3d :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    answer.resize(6, 6);
    answer.zero();
    if ( !this->isActivated(tStep) ) {
        return;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    double density = this->giveStructuralCrossSection()->give('d', gp);
    double halfMass = density * this->giveCrossSection()->give(CS_Area, gp) * this->computeLength() * 0.5;
    answer.at(1, 1) = halfMass;
    answer.at(2, 2) = halfMass;
    answer.at(3, 3) = halfMass;
    answer.at(4, 4) = halfMass;
    answer.at(5, 5) = halfMass;
    answer.at(6, 6) = halfMass;
}


void
Truss3d :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    FloatArray n;
    this->interp.evalN( n, iLocCoord, FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n, 3);
}


double
Truss3d :: computeVolumeAround(GaussPoint *gp)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double detJ = this->interp.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    double weight  = gp->giveWeight();
    return detJ *weight *this->giveCrossSection()->give(CS_Area, gp);
}


int
Truss3d :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    FloatArray lx, ly(3), lz;

    lx.beDifferenceOf( * this->giveNode(2)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );
    lx.normalize();

    ly(0) = lx(1);
    ly(1) = -lx(2);
    ly(2) = lx(0);

    // Construct orthogonal vector
    double npn = ly.dotProduct(lx);
    ly.add(-npn, lx);
    ly.normalize();
    lz.beVectorProductOf(ly, lx);

    answer.resize(3, 3);
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = lx.at(i);
        answer.at(2, i) = ly.at(i);
        answer.at(3, i) = lz.at(i);
    }

    return 1;
}


IRResultType
Truss3d :: initializeFrom(InputRecord *ir)
{
	IRResultType result;                    // Required by IR_GIVE_FIELD macro

	this->macroElem = 0;
	IR_GIVE_OPTIONAL_FIELD(ir, this->macroElem, _IFT_Truss3d_macroElem);

    return NLStructuralElement :: initializeFrom(ir);
}


void
Truss3d :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveRealStress_1d(answer, gp, strain, tStep);
}

void
Truss3d :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveStiffnessMatrix_1d(answer, rMode, gp, tStep);
}

void
Truss3d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, D_w};
}


void
Truss3d :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    if ( iEdge != 1 ) {
        OOFEM_ERROR("wrong edge number");
    }

    answer.resize(6);
    answer.at(1) = 1;
    answer.at(2) = 2;
    answer.at(3) = 3;
    answer.at(4) = 4;
    answer.at(5) = 5;
    answer.at(6) = 6;
}


double
Truss3d :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    if ( iEdge != 1 ) { // edge between nodes 1 2
        OOFEM_ERROR("wrong edge number");
    }

    double weight = gp->giveWeight();
    return this->interp.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) * weight;
}


int
Truss3d :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    FloatMatrix lcs;
    this->giveLocalCoordinateSystem(lcs);
    answer.beTranspositionOf(lcs);

    return 1;
}

void
Truss3d :: computeLocalForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// Why is this function taken separately ?
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further subtract part corresponding to non-nodal loading.
{
	FloatArray helpLoadVector(1);
	answer.clear();

	// loop over body load array first
	int nBodyLoads = this->giveBodyLoadArray()->giveSize();
	for (int i = 1; i <= nBodyLoads; i++) {
		int id = bodyLoadArray.at(i);
		Load *load = domain->giveLoad(id);
		bcGeomType ltype = load->giveBCGeoType();
		if ((ltype == BodyLoadBGT) && (load->giveBCValType() == ForceLoadBVT)) {
			this->computeBodyLoadVectorAt(helpLoadVector, load, tStep, mode);
			if (helpLoadVector.giveSize()) {
				answer.add(helpLoadVector);
			}
		}
		else {
			if (load->giveBCValType() != TemperatureBVT && load->giveBCValType() != EigenstrainBVT) {
				// temperature and eigenstrain is handled separately at computeLoadVectorAt subroutine
				OOFEM_ERROR("body load %d is of unsupported type (%d)", id, ltype);
			}
		}
	}

	// loop over boundary load array
	int nBoundaryLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
	for (int i = 1; i <= nBoundaryLoads; i++) {
		int n = boundaryLoadArray.at(1 + (i - 1) * 2);
		int id = boundaryLoadArray.at(i * 2);
		Load *load = domain->giveLoad(n);
		BoundaryLoad* bLoad;
		if ((bLoad = dynamic_cast<BoundaryLoad*> (load))) {
			bcGeomType ltype = load->giveBCGeoType();
			if (ltype == EdgeLoadBGT) {
				this->computeBoundaryEdgeLoadVector(helpLoadVector, bLoad, id, ExternalForcesVector, mode, tStep, false);
				if (helpLoadVector.giveSize()) {
					answer.add(helpLoadVector);
				}
			}
			else if (ltype == SurfaceLoadBGT) {
				this->computeBoundarySurfaceLoadVector(helpLoadVector, bLoad, id, ExternalForcesVector, mode, tStep, false);
				if (helpLoadVector.giveSize()) {
					answer.add(helpLoadVector);
				}
			}
			else if (ltype == PointLoadBGT) {
				// id not used
				this->computePointLoadVectorAt(helpLoadVector, load, tStep, mode, false);
				if (helpLoadVector.giveSize()) {
					answer.add(helpLoadVector);
				}
			}
			else {
				OOFEM_ERROR("boundary load %d is of unsupported type (%d)", id, ltype);
			}
		}
	}


	// add exact end forces due to nonnodal loading applied indirectly (via sets)
	BCTracker *bct = this->domain->giveBCTracker();
	BCTracker::entryListType bcList = bct->getElementRecords(this->number);
	FloatArray help;

	for (BCTracker::entryListType::iterator it = bcList.begin(); it != bcList.end(); ++it) {
		GeneralBoundaryCondition *bc = this->domain->giveBc((*it).bcNumber);
		BodyLoad *bodyLoad;
		BoundaryLoad *boundaryLoad;
		if (bc->isImposed(tStep)) {
			if ((bodyLoad = dynamic_cast<BodyLoad*>(bc))) { // body load
				this->computeBodyLoadVectorAt(help, bodyLoad, tStep, VM_Total); // this one is local
				answer.add(help);
			}
			else if ((boundaryLoad = dynamic_cast<BoundaryLoad*>(bc))) {
				// compute Boundary Edge load vector in GLOBAL CS !!!!!!!
				this->computeBoundaryEdgeLoadVector(help, boundaryLoad, (*it).boundaryId,
					ExternalForcesVector, VM_Total, tStep, false);
				// get it transformed back to local c.s.
				// this->computeGtoLRotationMatrix(t);
				// help.rotatedWith(t, 'n');
				answer.add(help);
			}
		}
	}
}

void
Truss3d :: giveEndForcesVector(FloatArray &answer, TimeStep *tStep)
{
	// computes exact global end-forces vector
	FloatArray loadEndForces;
	NLStructuralElement::giveInternalForcesVector(answer, tStep, false);

	// add exact end forces due to nonnodal loading
	this->computeLocalForceLoadVector(loadEndForces, tStep, VM_Total); // will compute only contribution of loads applied directly on receiver (not using sets)
	if (loadEndForces.giveSize()) {
		answer.subtract(loadEndForces);
	}
}

void
Truss3d::printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
	FloatArray Fl, strain;
	//FloatMatrix dir;
	//this->giveLocalCoordinateSystem(dir);
	this->computeStrainVector(strain, this->integrationRulesArray[0]->getIntegrationPoint(0), tStep);
	this->giveStructuralCrossSection()->giveRealStress_1d(Fl, this->integrationRulesArray[0]->getIntegrationPoint(0), strain, tStep);
	// fprintf(file, "truss3D %d (%8d) dir 3 %.4e %.4e %.4e macroelem %d : %.4e\n", this->giveLabel(), number, dir.at(1,1), dir.at(1,2), dir.at(1,3), this->macroElem, Fl.at(1)); // 1st component only
	fprintf(file, "truss3D %d (%8d) macroelem %d : %.4e\n", this->giveLabel(), number, this->macroElem, Fl.at(1)); // 1st component only
}

#ifdef __OOFEG
void Truss3d :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Truss3d :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;
    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}
#endif
} // end namespace oofem
