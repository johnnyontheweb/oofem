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

#ifndef NodalRecoveryModule_h
#define NodalRecoveryModule_h

#include "exportmodule.h"
#include "nodalrecoverymodel.h"
#include "floatarray.h"
#include "../sm/engineeringmodels/responseSpectrum.h"
#include <map>
#include "cltypes.h"

///@name Input fields for Beam export module
//@{
#define _IFT_NodalRecoveryModule_Name "nrm"
#define _IFT_NodalRecoveryModule_stype "stype"
#define _IFT_NodalRecoveryModule_rtypes "rtypes"
#define _IFT_NodalRecoveryModule_isrespspec "rspec"
//@}

namespace oofem {
/**
 * Represents generic export module. It gives the internal variables chosen
 * by the user for all elements in the model.
 * Parts adapted from VTKExportModule.
 *
 * @author Francesco Pontarin
 * @author Giovanni Rinaldin
 */
class OOFEM_EXPORT NodalRecoveryModule : public ExportModule
{
protected:
    /// Stream for file.
    FILE *stream;
    /// Array for the beam diagrams
    FloatArray *res;
	bool isRespSpec=false; // tStep 0 does the magic
	NodalRecoveryModel::NodalRecoveryModelType stype;
	RSpecComboType modalCombo;
	double csi;
	ResponseSpectrum *rs;
	NodalRecoveryModel* smoother;
	IntArray internalVarsToExport;
	Set* elemSet;

	std::list<std::string> valueTypesStr;
	int regionDofMans;

public:
    /// Constructor. Creates empty Output Manager.
    NodalRecoveryModule(int n, EngngModel *e);
    /// Destructor.
    virtual ~NodalRecoveryModule();
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void doOutput(TimeStep *tStep, bool forcedOutput = false);
    virtual void initialize();
    virtual void terminate();
    virtual const char *giveClassName() const { return "NodalRecoveryModule"; }
    virtual const char *giveInputRecordName() const { return _IFT_NodalRecoveryModule_Name; }

	static int checkValidType(const char* name) { return 0; };

	/**
	* Export internal variables.
	*/
	void exportIntVars(FILE *stream, TimeStep *tStep);
	/**
	* Exports single variable.
	*/
	void exportIntVarAs(InternalStateType valID, InternalStateValueType type, FILE *stream, TimeStep *tStep);
	/**
	* Assembles the region node map. Also computes the total number of nodes in region.
	* The region are numbered starting from offset+1.
	* if mode == 0 then regionNodalNumbers is array with mapping from global numbering to local region numbering.
	* The i-th value contains the corresponding local region number (or zero, if global numbar is not in region).
	* if mode == 1 then regionNodalNumbers is array with mapping from local to global numbering.
	* The i-th value contains the corresponding global node number.
	*/
	int initRegionNodeNumbering(IntArray &regionNodalNumbers, int &regionDofMans,
		int offset, Domain *domain, int reg, int mode);

private:
	
	std::map< int, std::map< double, FloatArray > >BeamForces;
	std::map< int, std::map< double, FloatArray > >BeamDisplacements;

	std::list<std::map< int, std::map< double, FloatArray > > >BeamForcesList;
	std::list<std::map< int, std::map< double, FloatArray > > >BeamDisplacementsList;

	std::map< int, std::map< double, FloatArray > >combBeamForces;
	std::map< int, std::map< double, FloatArray > >combBeamDisplacements;

	virtual void populateElResults(std::map<int, std::map<double, FloatArray>> &answer, std::map<int, std::map<double, FloatArray>> &src);
	virtual void addMultiply(std::map<int, std::map<double, FloatArray>> &answer, std::map<int, std::map<double, FloatArray>> &src, std::map<int, std::map<double, FloatArray>> &src2, double fact = 1.0);
	virtual void calcRoot(std::map<int, std::map<double, FloatArray>> &answer);
	virtual void SRSS();
	virtual void CQC();

};
} // end namespace oofem

#endif
