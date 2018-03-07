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

#ifndef BeamExportModule_h
#define BeamExportModule_h

#include "exportmodule.h"
#include "floatarray.h"
#include "../sm/EngineeringModels/responsespectrum.h"
#include <map>
#include <math.h>
#include <cstring>

///@name Input fields for Beam export module
//@{
#define _IFT_BeamExportModule_Name "bem"
#define _IFT_BeamExportModule_isrespspec "rspec"
//@}

namespace oofem {
/**
 * Represents beam (2D and 3D) export module. It gives the beam diagram values (N_x, T_z, T_y, M_x, M_y, M_z)
 * internal diplacements and rotations as well as soil pressure (if subsoil material is specified)
 * for all beam elements in the model, in local coordinate system.
 * Winkler formulation for timoshenko beam based upon:
 * Russo Spena, Ramaglia, Lignola, Prota, "Closed-Form solution for the Timoshenko beam on elastic soil", 2015
 *
 * @author Vit Smilauer
 * @author Mikael Öhman
 * @author Giovanni Rinaldin
 * @author Francesco Pontarin
 */
class OOFEM_EXPORT BeamExportModule : public ExportModule
{
protected:
    /// Stream for file.
    FILE *stream;
    /// Array for the beam diagrams
    FloatArray *res;
	bool isRespSpec=false; // tStep 0 does the magic
	FloatArray periods;
	RSpecComboType modalCombo;
	double csi;
	ResponseSpectrum *rs;

#ifdef MEMSTR
	bool usestream = true;
#endif

public:
    /// Constructor. Creates empty Output Manager.
    BeamExportModule(int n, EngngModel *e);
    /// Destructor.
    virtual ~BeamExportModule();
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void doOutput(TimeStep *tStep, bool forcedOutput = false);
    virtual void initialize();
    virtual void terminate();
    virtual const char *giveClassName() const { return "BeamExportModule"; }
    virtual const char *giveInputRecordName() const { return _IFT_BeamExportModule_Name; }

	static int checkValidType(const char* name) { return (strcmp(name, "Beam3d") == 0) || (strcmp(name, "Beam2d") == 0) || (strcmp(name, "beam3d") == 0) || (strcmp(name, "beam2d") == 0); };

private:
	
	std::map< int, std::map< double, FloatArray > >BeamForces;
	std::map< int, std::map< double, FloatArray > >BeamDisplacements;
	std::map< int, std::map< double, FloatArray > >BeamWinkler;

	std::list<std::map< int, std::map< double, FloatArray > > >BeamForcesList;
	std::list<std::map< int, std::map< double, FloatArray > > >BeamDisplacementsList;
	std::list<std::map< int, std::map< double, FloatArray > > >BeamWinklerList;

	std::map< int, std::map< double, FloatArray > >combBeamForces;
	std::map< int, std::map< double, FloatArray > >combBeamDisplacements;
	std::map< int, std::map< double, FloatArray > >combBeamWinkler;

	virtual void populateElResults(std::map<int, std::map<double, FloatArray>> &answer, std::map<int, std::map<double, FloatArray>> &src);
	virtual void addMultiply(std::map<int, std::map<double, FloatArray>> &answer, std::map<int, std::map<double, FloatArray>> &src, std::map<int, std::map<double, FloatArray>> &src2, double fact = 1.0);
	virtual void calcRoot(std::map<int, std::map<double, FloatArray>> &answer);
	virtual void SRSS();
	virtual void CQC();

};
} // end namespace oofem

#endif
