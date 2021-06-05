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

#include "AbaqusUserElement.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "dofmanager.h"
#include "dof.h"
#include "node.h"

#ifdef _WIN32 //_MSC_VER and __MINGW32__ included
 #include <Windows.h>
#else
 #include <dlfcn.h>
#endif
#include <cstring>

namespace oofem {
REGISTER_Element(AbaqusUserElement);


AbaqusUserElement :: AbaqusUserElement(int n, Domain *d) :
    NLStructuralElement(n, d), uelobj(NULL), hasTangentFlag(false), uel(NULL)
{}

AbaqusUserElement :: ~AbaqusUserElement()
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


void AbaqusUserElement :: initializeFrom(InputRecord &ir)
{
    StructuralElement :: initializeFrom(ir);

    this->numberOfDofMans = dofManArray.giveSize();

    // necessary to prevent an array dimension error in Init
	this->nCoords = 1;
	IR_GIVE_FIELD(ir, nCoords, _IFT_AbaqusUserElement_numcoords);

	IR_GIVE_FIELD(ir, this->dofs, _IFT_AbaqusUserElement_dofs);

    IR_GIVE_FIELD(ir, this->numSvars, _IFT_AbaqusUserElement_numsvars);
    if ( this->numSvars < 0 ) {
        OOFEM_ERROR("'numsvars' field has an invalid value");
    }
	IR_GIVE_FIELD(ir, this->props, _IFT_AbaqusUserElement_properties);

	IR_GIVE_OPTIONAL_FIELD(ir, this->jprops, _IFT_AbaqusUserElement_iproperties);

    IR_GIVE_FIELD(ir, this->jtype, _IFT_AbaqusUserElement_type);
    if ( this->jtype < 0 ) {
        OOFEM_ERROR("'type' has an invalid value");
    }

	this->macroElem = 0;
	IR_GIVE_OPTIONAL_FIELD(ir, this->macroElem, _IFT_AbaqusUserElement_macroElem);

    IR_GIVE_FIELD(ir, this->filename, _IFT_AbaqusUserElement_userElement);

#if 0
    uelname = "uel";
    IR_GIVE_OPTIONAL_FIELD(ir, uelname, _IFT_AbaqusUserElement_name);
#endif

#ifdef _WIN32
    this->uelobj = ( void * ) LoadLibrary( filename.c_str() );
    if ( !this->uelobj ) {
        OOFEM_ERROR( "couldn't load \"%s\",\ndlerror: %s", filename.c_str() );
    }

    * ( FARPROC * ) ( & this->uel ) = GetProcAddress( ( HMODULE ) this->uelobj, "uel_" );  //works for MinGW 32bit
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

    * ( void ** ) ( & this->uel ) = dlsym(this->uelobj, "uel_");
    char *dlresult = dlerror();
    if ( dlresult ) {
        OOFEM_ERROR("couldn't load symbol uel,\ndlerror: %s\n", dlresult);
    }
#endif
}


void AbaqusUserElement :: postInitialize()
{
    NLStructuralElement :: postInitialize();

    this->ndofel = this->numberOfDofMans * dofs.giveSize(); // this->nCoords;
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
		for (auto j : this->dofs)
		{
			switch ((DofIDItem)j)
			{
			case D_u:
			case D_v:
			case D_w:
				this->mcrd = (std::max)(this->mcrd, j);
			}
		}

		this->mcrd = (std::max)(this->mcrd, this->nCoords);

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


void AbaqusUserElement :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralElement :: giveInputRecord(input);

    input.setField(this->coords, _IFT_AbaqusUserElement_numcoords);
    input.setField(this->dofs, _IFT_AbaqusUserElement_dofs);
    input.setField(this->numSvars, _IFT_AbaqusUserElement_numsvars);
    input.setField(this->props, _IFT_AbaqusUserElement_properties);
    input.setField(this->jtype, _IFT_AbaqusUserElement_type);
    input.setField(this->filename, _IFT_AbaqusUserElement_userElement);
}

void AbaqusUserElement :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = this->dofs;
}

void AbaqusUserElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
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

void AbaqusUserElement :: updateYourself(TimeStep *tStep)
{
    StructuralElement :: updateYourself(tStep);
    svars = tempSvars;
    amatrx = tempAmatrx;
    rhs = tempRHS;
    hasTangentFlag = false;
}

void AbaqusUserElement :: updateInternalState(TimeStep *tStep)
{
    FloatArray tmp;
    this->giveInternalForcesVector(tmp, tStep, 0);
}

void AbaqusUserElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // init U vector
    // this->computeVectorOf(this->dofs, VM_Total, tStep, U);
	this->computeVectorOf(VM_Total, tStep, U, false);
	// get A and V
	this->computeVectorOf(VM_Velocity, tStep, V, false);
	this->computeVectorOf(VM_Acceleration, tStep, A, false);
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

void AbaqusUserElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, 
                                            FloatArray &U, FloatMatrix &DU, int useUpdatedGpRecord)
{
	answer.clear();
    if ( useUpdatedGpRecord ) {
        this->rhs.copyColumn(answer, 1);
    } else {
        this->lFlags.at(1) = 1;                 // 1 based access
        this->lFlags.at(3) = 1;                 // 1 based access
        this->lFlags.at(4) = 0;                 // 1 based access
		int label = this->giveLabel();
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
            & ( label ),
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
        //FloatArray vout;
        //vout.resize(12);
        //for (int i = 1; i <= 3; i++)
        //{
        //	vout.at(i) = rhs.at(i, 1);
        //	vout.at(i+6) = rhs.at(i+3, 1);
        //}
        //answer = vout;
        //this->rhs.copyColumn(answer, 1);
        //answer.negated();
        loc_rhs.negated();                      //really needed???
        loc_rhs.copyColumn(answer, 1);
        letTempRhsBe(loc_rhs);
        letTempTangentBe(loc_amatrx);
        letTempSvarsBe(loc_svars);
    }
}


void
AbaqusUserElement :: computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity)
{
    answer.resize(ndofel, ndofel);
    answer.zero();
}


void
AbaqusUserElement::printOutputAt(FILE *File, TimeStep *tStep)
{
	FloatArray rl, Fl;

	fprintf(File, "abaqususerelement %d (%8d) macroelem %d :\n", this->giveLabel(), this->giveNumber(), this->macroElem);

	// ask for global element displacement vector
	this->computeVectorOf(VM_Total, tStep, rl, false);
	// ask for global element end forces vector
	this->giveInternalForcesVector(Fl, tStep, 1);

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

void
AbaqusUserElement::computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
{
	// computes initial stress matrix of receiver (or geometric stiffness matrix)
	answer.resize(ndofel, ndofel);
	answer.zero();
}

}       // namespace oofem
