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
#include "gausspoint.h"
#include "engngm.h"
#include "material.h"
#include "classfactory.h"

namespace oofem {
REGISTER_ExportModule(BeamExportModule)

BeamExportModule::BeamExportModule(int n, EngngModel *e) : ExportModule(n, e) { }

BeamExportModule :: ~BeamExportModule() { }

IRResultType
BeamExportModule::initializeFrom(InputRecord *ir)
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

	// loop through the beam elements
    FloatArray ipState;
    Domain *d = emodel->giveDomain(1);
    for ( auto &elem : d->giveElements() ) {
		if (strcmp(elem->giveClassName(), "beam3d") || strcmp(elem->giveClassName(), "beam2d")) { // check if elem is beam (LIbeam?)

		//	for ( GaussPoint *gp: *elem->giveDefaultIntegrationRulePtr() ) {
		//		double dV = elem->computeVolumeAround(gp);
		//		
		//		elem->giveGlobalIPValue(ipState, gp, (InternalStateType)1, tStep); // IST_StressTensor
		//		
		//	}

		}
    }

	// loop through the loads
	//	d->giveSets or d->giveLoad ?


	// write file in the format:
	// elementNumber distanceFromIend N_x T_z T_y M_x M_y M_z
	// if 3 Gauss points are used, there would be 5 lines per beam (at distances 0, 0.2254*L, 0.5*L, 0.7746*L, L), ->>> check

    //fprintf(this->stream, "%d ", avgState.giveSize());
    //for ( auto s: avgState ) {
    //    fprintf(this->stream, "%e ", s);
    //}
    fprintf(this->stream, "    ");
    
    fprintf(this->stream, "\n" );
    fflush(this->stream);
}

void
BeamExportModule :: initialize()
{
    std :: string fileName = emodel->giveOutputBaseFileName() + ".bem";
    if ( ( this->stream = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR( "failed to open file %s", fileName.c_str() );
    }

    fprintf(this->stream, "#Time      Volume       ");
    //for ( int var: this->ists ) {
    //    fprintf(this->stream, "%s    ", __InternalStateTypeToString( ( InternalStateType ) var) );
    //}
    fprintf(this->stream, "\n" );
    fflush(this->stream);
}

void
BeamExportModule :: terminate()
{
    fclose(this->stream);
}
} // end namespace oofem
