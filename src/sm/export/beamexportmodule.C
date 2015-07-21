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

#include "beamexportmodule.h"
#include "timestep.h"
#include "element.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/Beams/Beam3d.h"
#include "gausspoint.h"
#include "engngm.h"
#include "material.h"
#include "classfactory.h"
#include "generalboundarycondition.h"
#include "constantedgeload.h"

namespace oofem {
REGISTER_ExportModule(BeamExportModule)

BeamExportModule :: BeamExportModule(int n, EngngModel *e) : ExportModule(n, e) { }

BeamExportModule :: ~BeamExportModule() { }

IRResultType
BeamExportModule :: initializeFrom(InputRecord *ir)
{
    //IRResultType result;                 // Required by IR_GIVE_FIELD macro
    return ExportModule :: initializeFrom(ir);
}

void
BeamExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }

    std :: vector< int >beamIDs;
    std :: map< int, std :: map< double, FloatArray > >BeamForces;
    IntArray temp;
    // loop through the beam elements
    Domain *d = emodel->giveDomain(1);
    for ( auto &elem : d->giveElements() ) {
        if ( strcmp(elem->giveClassName(), "Beam3d") == 0 || strcmp(elem->giveClassName(), "Beam2d") == 0 ) {   // check if elem is beam (LIbeam?)
            // store IDs of known beams
            beamIDs.push_back( elem->giveNumber() );

            StructuralElement *SElem;
            SElem = static_cast< StructuralElement * >( elem.get() );

            double ksi, l = elem->computeLength();
            FloatArray Fl, loadEndForces;

            SElem->giveInternalForcesVector(Fl, tStep);

            // add exact end forces due to nonnodal loading
            SElem->computeForceLoadVector(loadEndForces, tStep, VM_Total);
            if ( loadEndForces.giveSize() ) {
                Fl.subtract(loadEndForces);
            }

            //SElem->giveEndForcesVector(Fl, tStep);

            std :: map< double, FloatArray >Dict;

            //fprintf(this->stream, "\nBeam: %d \n", elem->giveNumber());
            //fprintf(this->stream, "length: %e \n", l);

            FloatArray I, E, Diff;

            //fprintf(this->stream, "End Forces: ");
            //for (auto &val : Fl) {
            //	fprintf(this->stream, " %.4e", val);
            //}

            I.resize(6);

            temp.resize(6);
            temp.at(1) = 1;
            temp.at(2) = 2;
            temp.at(3) = 3;
            temp.at(4) = 4;
            temp.at(5) = 5;
            temp.at(6) = 6;
            I.beSubArrayOf(Fl, temp);
            //I *= -1;
            Dict [ 0.0 ] = I;

            E.resize(6);
            for ( int i = 1; i < 7; i++ ) {
                temp.at(i) = temp.at(i) + 6;
            }
            E.beSubArrayOf(Fl, temp);

            Diff.beDifferenceOf(E, I);

            fprintf(this->stream, "\n");

            for ( GaussPoint *gp: *elem->giveDefaultIntegrationRulePtr() ) {
                //double dV = elem->computeVolumeAround(gp);
                FloatArray ipState;
                ksi = 0.5 + 0.5 * gp->giveNaturalCoordinate(1);

                // can't use this until beam is fixed?
                // elem->giveGlobalIPValue(ipState, gp, (InternalStateType)1, tStep); // IST_StressTensor
                ipState.zero();
                ipState.beScaled(ksi, Diff);
                ipState.add(I);
                Dict [ ksi * l ] = ipState;

                //fprintf(this->stream, "gp: %d pos: %e forces:" , gp->giveNumber(), ksi*l);
                //for (auto &val : ipState) {
                //	fprintf(this->stream, " %.4e", val);
                //}
            }

            Dict [ l ] = E;

            BeamForces [ elem->giveNumber() ] = Dict;

            //elem->giveBodyLoadArray
        }
    }

    //std::vector< std::unique_ptr< GeneralBoundaryCondition > > BCs = d->giveBcs();

    // tamper with stuff only if sets are defined.
    if ( d->giveNumberOfSets() ) {
        temp.resize(6);
        temp.at(1) = 1;
        temp.at(2) = 2;
        temp.at(3) = 3;
        temp.at(4) = 4;
        temp.at(5) = 5;
        temp.at(6) = 6;
        // loop through the loads
        for ( auto &bc : d->giveBcs() ) {
            // int bType = bc->giveBCValType(); // UNUSED: ConstantEdgeLoad is never == 2, they're all == 0 == unknown
            //if (bc->giveBCValType() == ForceLoadBVT) {
            if ( strcmp(bc->giveClassName(), "ConstantEdgeLoad") == 0 ) {
                ConstantEdgeLoad *CLoad = static_cast< ConstantEdgeLoad * >( bc.get() );

                // is it in a set?
                int nSet = CLoad->giveSetNumber();
                if ( nSet ) {
                    Set *mySet = d->giveSet(nSet);
                    // contains any of our beams?
                    const IntArray &EdgeList = mySet->giveEdgeList();
                    const int numEdges = EdgeList.giveSize() / 2;

                    int c = 1;
                    while ( c <= numEdges ) {
                        FloatArray compArr, tempArr;
                        FloatMatrix T;
                        int elNum = EdgeList.at(c), edgeNum = EdgeList.at(++c);
                        CLoad->computeValues(compArr, tStep, NULL, temp, VM_Total);
                        //d->giveElement(elNum)->computeBoundaryEdgeLoadVector(compArr, CLoad, edgeNum, ExternalForcesVector, VM_Total, tStep); // always vm_total???

                        // transform to local coordinates
                        d->giveElement(elNum)->computeGtoLRotationMatrix(T);
                        T.resizeWithData(6, 6);
                        compArr.rotatedWith(T, 'n');

                        // compute contribution to internal forces
                        std :: map< double, FloatArray >Dst = BeamForces [ elNum ];
                        for ( auto &PointVals : Dst ) {
                            const double &pos = PointVals.first;
                            FloatArray Vals = PointVals.second;

                            // tamper with values?
                            //Vals.add(pos, compArr);

                            // only these two parabolic. linear interpolation should be enough for the other directions
                            Vals.at(5) += compArr.at(3) * ( pos ) / 2 * ( d->giveElement(elNum)->computeLength() - pos );
                            Vals.at(6) -= compArr.at(2) * ( pos ) / 2 * ( d->giveElement(elNum)->computeLength() - pos );

                            // update in point-forces map
                            PointVals.second = Vals;
                        }
                        // update in beam-forces map
                        BeamForces [ elNum ] = Dst;
                    }
                }
            }
            //}
        }
    }

    //for (auto &set : d->giveSets()) {
    //	IntArray &ElEdges = set->giveEdgeList();
    //}


    //	d->giveSets or d->giveLoad ?

    double curTime = tStep->giveTargetTime();

    for ( auto &bForces : BeamForces ) {
        std :: map< double, FloatArray >pForces = bForces.second;
        int ID = bForces.first;
        for ( auto &vals : pForces ) {
            double pos = vals.first;
            FloatArray forces = vals.second;
            fprintf(this->stream, "%10.3e;%d;%10.3e;", curTime, ID, pos);
            for ( auto &val : forces ) {
                fprintf(this->stream, "%10.3e;", val);
            }
            fprintf(this->stream, "\n");
        }
    }

    // write file in the format:
    // elementNumber distanceFromIend N_x T_z T_y M_x M_y M_z
    // if 3 Gauss points are used, there would be 5 lines per beam (at distances 0, 0.1127*L, 0.5*L, 0.8873*L, L), ->>> to check

    //fprintf(this->stream, "%d ", avgState.giveSize());
    //for ( auto s: avgState ) {
    //    fprintf(this->stream, "%e ", s);
    //}
    //fprintf(this->stream, "    ");

    fprintf(this->stream, "\n");
    fflush(this->stream);
}

void
BeamExportModule :: initialize()
{
    std :: string fileName = emodel->giveOutputBaseFileName() + ".bem";
    if ( ( this->stream = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR( "failed to open file %s", fileName.c_str() );
    }
    // ";" as separator
    fprintf(this->stream, "#Time;BeamNo;DistanceFromI;N_x;T_y;T_z;M_x;M_y;M_z;");
    //for ( int var: this->ists ) {
    //    fprintf(this->stream, "%s    ", __InternalStateTypeToString( ( InternalStateType ) var) );
    //}
    fprintf(this->stream, "\n");
    fflush(this->stream);
}

void
BeamExportModule :: terminate()
{
    fclose(this->stream);
}
} // end namespace oofem
