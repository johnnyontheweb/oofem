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

#include "AbaqusUserElement6d.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "dofmanager.h"
#include "node.h"
#include "math.h"

#ifdef _WIN32 //_MSC_VER and __MINGW32__ included
 #include <Windows.h>
#else
 #include <dlfcn.h>
#endif
#include <cstring>

namespace oofem {
REGISTER_Element(AbaqusUserElement6d);

AbaqusUserElement6d :: AbaqusUserElement6d(int n, Domain *d) :
    NLStructuralElement(n, d), uelobj(NULL), hasTangentFlag(false), uel(NULL)
{}

AbaqusUserElement6d :: ~AbaqusUserElement6d()
{
#ifdef _WIN32
    if ( this->uelobj ) {
        FreeLibrary( ( HMODULE ) this->uelobj );
    }
#else
    if ( this->uelobj ) {
        dlclose(this->uelobj);
    }
#endif
}


IRResultType AbaqusUserElement6d :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                                        // Required by IR_GIVE_FIELD macro

    result = StructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    this->numberOfDofMans = dofManArray.giveSize();

    // necessary to prevent an array dimension error in Init
	this->nCoords = 1;
	// IR_GIVE_OPTIONAL_FIELD(ir, nCoords, _IFT_AbaqusUserElement6d_numcoords);

	//IR_GIVE_FIELD(ir, this->dofs, _IFT_AbaqusUserElement6d_dofs);
	this->dofs.resize(6); 
	this->dofs.at(1) = 1; 	this->dofs.at(2) = 2; 	this->dofs.at(3) = 3;
	this->dofs.at(4) = 4;   this->dofs.at(5) = 5; 	this->dofs.at(6) = 6;

	IR_GIVE_FIELD(ir, this->dir, _IFT_AbaqusUserElement6d_orientation);
	if ( this->dir.giveSize() ) {
		this->dir.normalize();
	}
	
	// mode
	// this->mode = 0;
	// IR_GIVE_OPTIONAL_FIELD(ir, mode, _IFT_AbaqusUserElement6d_mode);

    IR_GIVE_FIELD(ir, this->numSvars, _IFT_AbaqusUserElement6d_numsvars);
    if ( this->numSvars < 0 ) {
        OOFEM_ERROR("'numsvars' field has an invalid value");
    }
	IR_GIVE_FIELD(ir, this->props, _IFT_AbaqusUserElement6d_properties);

	IR_GIVE_OPTIONAL_FIELD(ir, this->jprops, _IFT_AbaqusUserElement6d_iproperties);

	this->jtype = 0;
	IR_GIVE_OPTIONAL_FIELD(ir, this->jtype, _IFT_AbaqusUserElement6d_type);
    if ( this->jtype < 0 ) {
        OOFEM_ERROR("'type' has an invalid value");
    }

	if (ir->hasField(_IFT_AbaqusUserElement6d_refangle)) {
		IR_GIVE_FIELD(ir, referenceAngle, _IFT_AbaqusUserElement6d_refangle);
	}

	this->macroElem = 0;
	IR_GIVE_OPTIONAL_FIELD(ir, this->macroElem, _IFT_AbaqusUserElement6d_macroElem);

    IR_GIVE_FIELD(ir, this->filename, _IFT_AbaqusUserElement6d_userElement);

#if 0
    uelname = "uel";
    IR_GIVE_OPTIONAL_FIELD(ir, uelname, _IFT_AbaqusUserElement6d_name);
#endif

#ifdef _WIN32
    this->uelobj = ( void * ) LoadLibrary( filename.c_str() );
    if ( !this->uelobj ) {
        OOFEM_ERROR( "couldn't load \"%s\",\ndlerror: %s", filename.c_str() );
    }

    * ( FARPROC * ) ( & this->uel ) = GetProcAddress( ( HMODULE ) this->uelobj, "uel3d_" );  //works for MinGW 32bit
    if ( !this->uel ) {
        // char *dlresult = GetLastError();
        DWORD dlresult = GetLastError();                 //works for MinGW 32bit
        OOFEM_ERROR("couldn't load symbol uel,\nerror: %s\n", dlresult);
    }
	// optional init call for initialization vectors
	* ( FARPROC * ) ( & this->uelInit ) = GetProcAddress( ( HMODULE ) this->uelobj, "init_" );  //works for MinGW 32bit
	if (this->uelInit) this->uelInit(); // call with no error management, this is optional.
#else
    this->uelobj = dlopen(filename.c_str(), RTLD_NOW);
    if ( !this->uelobj ) {
        OOFEM_ERROR( "couldn't load \"%s\",\ndlerror: %s", filename.c_str(), dlerror() );
    }

    * ( void ** ) ( & this->uel ) = dlsym(this->uelobj, "uel3d_");
    char *dlresult = dlerror();
    if ( dlresult ) {
        OOFEM_ERROR("couldn't load symbol uel,\ndlerror: %s\n", dlresult);
    }
#endif

    return IRRT_OK;
}


void AbaqusUserElement6d :: postInitialize()
{
    NLStructuralElement :: postInitialize();

    this->ndofel = 12; // 6d spring
    this->mlvarx = this->ndofel;
    this->nrhs = 2;
    this->rhs.resize(this->ndofel, this->nrhs);
    this->amatrx.resize(this->ndofel, this->ndofel);
    this->svars.resize(this->numSvars);
    this->lFlags.resize(5);
    this->predef.resize( this->npredef * this->numberOfDofMans * 2 );
    this->energy.resize(8);
    this->U.resize(this->ndofel);
    this->V.resize(this->ndofel);
    this->A.resize(this->ndofel);
    this->DU.resize(this->ndofel, this->nrhs);

	if (!this->coords.isNotEmpty()) {
		this->mcrd = 1;
		//for (auto j : this->dofs)
		//{
		//	switch ((DofIDItem)j)
		//	{
		//	case D_u:
		//	case D_v:
		//	case D_w:
		//		this->mcrd = (std::max)(this->mcrd, j);
		//	}
		//}

		//this->mcrd = (std::max)(this->mcrd, this->nCoords);

		this->coords.resize(this->mcrd, this->numberOfDofMans);
		// this->coords.resize(this->numberOfDofMans, this->mcrd);
        for ( int j = 1; j <= numberOfDofMans; j++ ) {
            Node *dm = this->giveNode(j);
            for ( int i = 1; i <= mcrd; i++ ) {
                this->coords.at(i, j) = dm->giveCoordinate(i);
            }
        }
    }
}


void AbaqusUserElement6d :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralElement :: giveInputRecord(input);

    // input.setField(this->coords, _IFT_AbaqusUserElement6d_numcoords);
    //input.setField(this->dofs, _IFT_AbaqusUserElement6d_dofs);
    input.setField(this->numSvars, _IFT_AbaqusUserElement6d_numsvars);
    input.setField(this->props, _IFT_AbaqusUserElement6d_properties);
    input.setField(this->jtype, _IFT_AbaqusUserElement6d_type);
    input.setField(this->filename, _IFT_AbaqusUserElement6d_userElement);
}

Interface *AbaqusUserElement6d :: giveInterface(InterfaceType it)
{
    return NULL;
}

void AbaqusUserElement6d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
	answer = { D_u, D_v, D_w, R_u, R_v, R_w };
}

bool
AbaqusUserElement6d::computeGtoLRotationMatrix(FloatMatrix &answer)
{
	/*
	* Spring is defined as 1D element along orientation axis (or around orientation axis for torsional springs)
	* The transformation from local (1d) to global system typically expand dimensions
	*/
	
	answer.resize(6, 6);
	FloatArray lx, ly, lz, help(3);

	lx.resize(3);
	lx.at(1) = dir.at(1); lx.at(2) = dir.at(2); lx.at(3) = dir.at(3);

	FloatMatrix rot(3, 3);
	double theta = referenceAngle * M_PI / 180.0;;

	rot.at(1, 1) = cos(theta) + pow(lx.at(1), 2) * (1 - cos(theta));
	rot.at(1, 2) = lx.at(1) * lx.at(2) * (1 - cos(theta)) - lx.at(3) * sin(theta);
	rot.at(1, 3) = lx.at(1) * lx.at(3) * (1 - cos(theta)) + lx.at(2) * sin(theta);

	rot.at(2, 1) = lx.at(2) * lx.at(1) * (1 - cos(theta)) + lx.at(3) * sin(theta);
	rot.at(2, 2) = cos(theta) + pow(lx.at(2), 2) * (1 - cos(theta));
	rot.at(2, 3) = lx.at(2) * lx.at(3) * (1 - cos(theta)) - lx.at(1) * sin(theta);

	rot.at(3, 1) = lx.at(3) * lx.at(1) * (1 - cos(theta)) - lx.at(2) * sin(theta);
	rot.at(3, 2) = lx.at(3) * lx.at(2) * (1 - cos(theta)) + lx.at(1) * sin(theta);
	rot.at(3, 3) = cos(theta) + pow(lx.at(3), 2) * (1 - cos(theta));

	help.at(3) = 1.0;         // up-vector
	// here is ly is used as a temp var
	// double prvect = acos(lx.dotProduct(help));
	// if (prvect < 0.001 || prvect > M_PI - 0.001) { // Check if it is vertical
	if (fabs(lx.dotProduct(help)) > 0.999) { // Check if it is vertical
		ly = { 0., 1., 0. };
	}
	else {
		ly.beVectorProductOf(lx, help);
	}
	lz.beProductOf(rot, ly);
	lz.normalize();

	ly.beVectorProductOf(lz, lx);
	ly.normalize();

	FloatMatrix lcs;
	lcs.resize(3, 3);
	lcs.zero();
	for (int i = 1; i <= 3; i++) {
		lcs.at(1, i) = lx.at(i);
		lcs.at(2, i) = ly.at(i);
		lcs.at(3, i) = lz.at(i);
	}

	answer.resize(computeNumberOfGlobalDofs(), computeNumberOfGlobalDofs()); // 12 x 12
	answer.zero();

	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			answer.at(i, j) = lcs.at(i, j);
			answer.at(i + 3, j + 3) = lcs.at(i, j);
			answer.at(i + 6, j + 6) = lcs.at(i, j);
			answer.at(i + 9, j + 9) = lcs.at(i, j);
		}
	}
	
	//answer.at(1, 1) = this->dir.at(1);
	//answer.at(2, 2) = this->dir.at(2);
	//answer.at(3, 3) = this->dir.at(3);
	//answer.at(4, 4) = this->dir.at(1);
	//answer.at(5, 5) = this->dir.at(2);
	//answer.at(6, 6) = this->dir.at(3);

	//answer.at(7, 7) = this->dir.at(1);
	//answer.at(8, 8) = this->dir.at(2);
	//answer.at(9, 9) = this->dir.at(3);
	//answer.at(10, 10) = this->dir.at(1);
	//answer.at(11, 11) = this->dir.at(2);
	//answer.at(12, 12) = this->dir.at(3);
	return true;
	
}

void AbaqusUserElement6d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    if ( !hasTangent() ) {
        // use uel to calculate the tangent
        FloatArray forces;
        giveInternalForcesVector(forces, tStep, U, DU, 0);
    }
    // give tangent
    answer = giveTempTangent();
    // add stuff to behave differently if mUseNumericalTangent is set?
}

void AbaqusUserElement6d :: updateYourself(TimeStep *tStep)
{
    StructuralElement :: updateYourself(tStep);
    svars = tempSvars;
    amatrx = tempAmatrx;
    rhs = tempRHS;
    hasTangentFlag = false;
}

void AbaqusUserElement6d :: updateInternalState(TimeStep *tStep)
{
    FloatArray tmp;
    this->giveInternalForcesVector(tmp, tStep, 0);
}

void AbaqusUserElement6d :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // init U vector
    //this->computeVectorOf(this->dofs, VM_Total, tStep, U);
	this->computeVectorOf(VM_Total, tStep, U, false);
    FloatArray tempIntVect;
    // init DU vector
    //this->computeVectorOf(this->dofs, VM_Incremental, tStep, tempIntVect);
	this->computeVectorOf(VM_Incremental, tStep, tempIntVect, false);
    //this->giveDomain()->giveClassName();
    DU.zero();
    DU.setColumn(tempIntVect, 1);
    //this->computeVectorOf(VM_Total, tStep, DU);
    this->giveInternalForcesVector(answer, tStep, U, DU, useUpdatedGpRecord);
}

void AbaqusUserElement6d :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, 
                                            FloatArray &U, FloatMatrix &DU, int useUpdatedGpRecord)
{
	answer.clear();
    if ( useUpdatedGpRecord ) {
        this->rhs.copyColumn(answer, 1);
    } else {
        this->lFlags.at(1) = 1;                 // 1 based access
        this->lFlags.at(3) = 1;                 // 1 based access
        this->lFlags.at(4) = 0;                 // 1 based access

        int nprops = props.giveSize();
        int njprops = jprops.giveSize();

        FloatMatrix loc_rhs(this->ndofel, this->nrhs);
        FloatMatrix loc_amatrx(this->ndofel, this->ndofel);
        FloatArray loc_svars = this->giveStateVector();

        //this->getSvars();
        double period = 0., pnewdt = 0.;
        double dtime = tStep->giveTimeIncrement();
		double ttime = tStep->giveTargetTime();
		double time[] = { ttime, ttime };

		// see description at http://abaqus.software.polimi.it/v6.12/books/sub/default.htm
        this->uel(
            loc_rhs.givePointer(),
            loc_amatrx.givePointer(),
            loc_svars.givePointer(),
            energy.givePointer(),
            & ndofel,
            & nrhs,
            & numSvars,
            props.givePointer(),
            & nprops,
            coords.givePointer(),
            & mcrd,
            & this->numberOfDofMans,
            U.givePointer(),
            DU.givePointer(),
            V.givePointer(),
            A.givePointer(),
            & jtype,
            time,
            & dtime,
            & kstep,
            & kinc,
            & ( this->number ),
            params,
            & ndLoad,
            jdltype,
            adlmag.givePointer(),
            predef.givePointer(),
            & npredef,
            lFlags.givePointer(),
            & mlvarx,
            ddlmag.givePointer(),
            & mdLoad,
            & pnewdt,
            jprops.givePointer(),
            & njprops,
            & period);

        loc_rhs.negated();                      //really needed???
        loc_rhs.copyColumn(answer, 1);
        letTempRhsBe(loc_rhs);
        letTempTangentBe(loc_amatrx);
        letTempSvarsBe(loc_svars);
    }
}


void
AbaqusUserElement6d :: computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity)
{
    answer.resize(ndofel, ndofel);
    answer.zero();
}


void
AbaqusUserElement6d::printOutputAt(FILE *File, TimeStep *tStep)
{
	FloatArray rl, Fl;
	// ask for global element displacement vector
	this->computeVectorOf(VM_Total, tStep, rl, false);
	// ask for global element end forces vector
	this->giveInternalForcesVector(Fl, tStep, 1);

	fprintf(File, "AbaqusUserElement6d %d (%8d) macroelem %d type %d dir 3 %.4e %.4e %.4e : %.4e %.4e %.4e %.4e %.4e %.4e \n", this->giveLabel(), this->giveNumber(), this->macroElem, this->jtype, this->dir.at(1), this->dir.at(2), this->dir.at(3), Fl.at(6 + 1), Fl.at(6 + 2), Fl.at(6 + 3), Fl.at(6 + 4), Fl.at(6 + 5), Fl.at(6 + 6));
	
	fprintf(File, "  local_displacements %d ", rl.giveSize());
	for (auto &val : rl) {
		fprintf(File, " %.4e", val);
	}
	
	fprintf(File, "\n  internal_forces %d ", Fl.giveSize());
	for (auto &val : Fl) {
		fprintf(File, " %.4e", val);
	}

	fprintf(File, "\n  element_svars %d ", this->numSvars);
	for (int i = 1; i <= this->numSvars; i++)
	{
		fprintf(File, " %.4e", this->svars.at(i));
	}
	fprintf(File, "\n");
}

}       // namespace oofem
