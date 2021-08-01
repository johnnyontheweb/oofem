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

#ifndef dummyenum_h
#define dummyenum_h

#include "structengngmodel.h"
#include "datastream.h"
#include "oofemtxtdatareader.h"
#include "elementgeometrytype.h"
#include <cstring>

///@name Input fields for DummyEnum
//@{
#define _IFT_DummyEnum_Name "dummyenum"
#define _IFT_DummyEnum_writeAll "writeall"
#define _IFT_DummyEnum_ist "ist"
#define _IFT_DummyEnum_unkType "unktype"
#define _IFT_DummyEnum_dofType "doftype"
#define _IFT_DummyEnum_domainType "domtype"
#define _IFT_DummyEnum_materialMode "matmode"
#define _IFT_DummyEnum_elementGeometry "elgeom"
#define _IFT_DummyEnum_valModeType "valmode"
#define _IFT_DummyEnum_dofIDItem "dofid"
#define _IFT_DummyEnum_charType "chartype"
#define _IFT_DummyEnum_materialResponseMode "matrespmode"
#define _IFT_DummyEnum_materialMappingAlgo "matmapalgo"
#define _IFT_DummyEnum_meshPackageType "meshpack"
#define _IFT_DummyEnum_xFEMStateType "xfemstate"
//@}

#define to_str(s) #s

namespace oofem {
/**
 * This class implements a dummy analysis that outputs all enum contents up to a requested id
 * This class is more a hack than a proper EngineeringModel overrides several EngineeringModel
 * just to allow the proper setup of the minimum number of objects required for oofem to work.
 * It may used in an automated debugging tool for a quick capability and compatibility check.
 * 
 * @author Francesco Pontarin
 */
class DummyEnum : public StructuralEngngModel
{
protected:
    int ist = 0;
    int unktype = 0;
    int doftype = 0;
    int domtype = 0;
    int matmode = 0;
    int elgeom = 0;
    int valmode = 0;
    int matrespmode = 0;
    int dofid = 0;
    int chartype = 0;
    int matmapalgo = 0;
    int meshpack = 0;
    int xfemstate = 0;

    int writeall = 1;

public:
    DummyEnum(int i, EngngModel * _master = NULL) : StructuralEngngModel(i, _master) { }
    virtual ~DummyEnum() { }

    void updateYourself(TimeStep *tStep) override;

    void terminate(TimeStep *tStep) override;

    int instanciateYourself(DataReader &dr, InputRecord &ir, const char* outFileName, const char* desc) override;

    void initializeFrom(InputRecord &ir) override;
    void postInitialize() override {} 
    virtual int checkProblemConsistency() { return 1; }
    //virtual void initStepIncrements() {};

    void init() override {}

    void solveYourself() override;

    // identification
    const char *giveClassName() const override { return "DummyEnum"; }
    virtual const char *giveInputRecordName() const { return _IFT_DummyEnum_Name; }

private:

    // the boolean is kept for eventual future writeall flags for each enum
    template <typename T, typename std::enable_if<std::is_enum<T>::value>::type* = nullptr>
    void doPrint(int num, bool all, const char* (*func)(T))
    {
	fprintf(outputStream, "\n\n");

	if (num){
	    fprintf(outputStream, "##############################");
	    fprintf(outputStream, "\n%s\n", typeid(T).name());
	    fprintf(outputStream, "%d - %d  ---  ", 0 , num);
	    if (all) {
		fprintf(outputStream, "complete\n", 0, num);
	    } else {
		fprintf(outputStream, "filtered\n", 0, num);
	    }

	    for (int c = 0; c <= num; c++){
		T enVal = (T)(c);
		const char *st = (*func)(enVal);
		// write everything only if requested, else filter out "Unknown" entries
		if (!all){
		    if (strcmp(st, "Unknown") == 0){ continue; }
		}
		// else write
		fprintf(outputStream, "  %4d : ", c);
		fprintf(outputStream, st);
		fprintf(outputStream, "\n");
	    }
	}
    }
}; // end of class

} // end namespace oofem
#endif // dummyenum_h
