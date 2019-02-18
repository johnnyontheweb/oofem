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

#ifndef IntElSurfQuad_h
#define IntElSurfQuad_h

#include "../sm/Elements/Interfaces/structuralinterfaceelement.h"
#include "fei3dhexalin.h"

#define _IFT_IntElSurfQuad_Name "intelsurfquad"

namespace oofem {
/**
 * This class implements 3d triangular surface interface element with linear interpolation.
 */
class IntElSurfQuad : public StructuralInterfaceElement
{
protected:
    static FEI3dHexaLin interpolation;

public:
    IntElSurfQuad(int n, Domain *d);
    virtual ~IntElSurfQuad() { }

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    virtual void computeCovarBaseVectorsAt(IntegrationPoint *ip, FloatArray &G1, FloatArray &G2);
    virtual void computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer);

    virtual int computeNumberOfDofs() { return 24; }
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;

    virtual double computeAreaAround(IntegrationPoint *ip);
    
    virtual void giveEngTraction(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
    {
        this->giveInterfaceCrossSection()->giveEngTraction_3d(answer, gp, jump, tStep);
    }

    virtual void giveStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, IntegrationPoint *ip, TimeStep *tStep)
    {
        this->giveInterfaceCrossSection()->give3dStiffnessMatrix_Eng(answer, rMode, ip, tStep);
    }

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_IntElSurfQuad_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
	virtual Element_Geometry_Type giveGeometryType() const { return EGT_hexa_1; }

    #ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    #endif

protected:
    virtual void computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer);
    virtual void computeGaussPoints();

};
} // end namespace oofem
#endif // IntElSurfQuad_h
