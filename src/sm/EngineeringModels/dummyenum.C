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

#include "sm/EngineeringModels/dummyenum.h"
#include "domain.h"
#include "classfactory.h"
#include "exportmodulemanager.h"
#include "xfem/xfemmanager.h"

namespace oofem {
REGISTER_EngngModel( DummyEnum );

void
DummyEnum::updateYourself( TimeStep *tStep ) { }

void
DummyEnum::terminate( TimeStep *tStep )
{
    int c = 0;

    fprintf( outputStream, "\n\n\tE N U M  O U T P U T:\n\t_______________________________\n" );

    this->doPrint<InternalStateType>( ist, (bool)( this->writeall ), &__InternalStateTypeToString );
    this->doPrint<UnknownType>( unktype, (bool)( this->writeall ), &__UnknownTypeToString );
    this->doPrint<dofType>( doftype, (bool)( this->writeall ), &__dofTypeToString );
    this->doPrint<domainType>( matmode, (bool)( this->writeall ), &__domainTypeToString );
    this->doPrint<MaterialMode>( unktype, (bool)( this->writeall ), &__MaterialModeToString );
    this->doPrint<Element_Geometry_Type>( elgeom, (bool)( this->writeall ), &__Element_Geometry_TypeToString );
    this->doPrint<ValueModeType>( valmode, (bool)( this->writeall ), &__ValueModeTypeToString );

    // damn
    fprintf( outputStream, "\n\n" );

    if ( dofid ) {
        int num = dofid;
        fprintf( outputStream, "##############################" );
        fprintf( outputStream, "\n%s\n", typeid( DofIDItem ).name() );
        fprintf( outputStream, "%d - %d  ---  ", 0, num );
        if ( writeall ) { fprintf( outputStream, "complete\n", 0, num ); } else { fprintf( outputStream, "filtered\n", 0, num ); }

        for ( int c = 0; c <= num; c++ ) {
            DofIDItem enVal = (DofIDItem)( c );
            std::string str = __DofIDItemToString( enVal ); // c_str does not live if std::string gets destroyed
            const char *st  = str.c_str();
            if ( !writeall ) { if ( strcmp( st, "Unknown" ) == 0 ) { continue; } }
            // else write
            fprintf( outputStream, "  %4d : ", c );
            fprintf( outputStream, st );
            fprintf( outputStream, "\n" );
        }
    }

    this->doPrint<CharType>( chartype, (bool)( this->writeall ), &__CharTypeToString );
    this->doPrint<MatResponseMode>( matrespmode, (bool)( this->writeall ), &__MatResponseModeToString );
    this->doPrint<MaterialMappingAlgorithmType>( matmapalgo, (bool)( this->writeall ), &__MaterialMappingAlgorithmTypeToString );
    this->doPrint<MeshPackageType>( meshpack, (bool)( this->writeall ), &__MeshPackageTypeToString );
    this->doPrint<XFEMStateType>( xfemstate, (bool)( this->writeall ), &__XFEMStateTypeToString );
}

void
DummyEnum::initializeFrom( InputRecord &ir )
// input from inputString
{
    writeall = 1;
    IR_GIVE_OPTIONAL_FIELD( ir, writeall, _IFT_DummyEnum_writeAll );

    ist = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, ist, _IFT_DummyEnum_ist );

    unktype = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, unktype, _IFT_DummyEnum_unkType );

    doftype = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, doftype, _IFT_DummyEnum_dofType );

    domtype = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, domtype, _IFT_DummyEnum_domainType );

    matmode = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, matmode, _IFT_DummyEnum_materialMode );

    elgeom = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, elgeom, _IFT_DummyEnum_elementGeometry );

    valmode = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, valmode, _IFT_DummyEnum_valModeType );

    dofid = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, dofid, _IFT_DummyEnum_dofIDItem );

    chartype = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, chartype, _IFT_DummyEnum_charType );

    matrespmode = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, matrespmode, _IFT_DummyEnum_materialResponseMode );

    matmapalgo = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, matmapalgo, _IFT_DummyEnum_materialMappingAlgo );

    meshpack = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, meshpack, _IFT_DummyEnum_meshPackageType );

    xfemstate = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, xfemstate, _IFT_DummyEnum_xFEMStateType );

    //StructuralEngngModel :: initializeFrom(ir);
}

int
DummyEnum::instanciateYourself( DataReader &dr, InputRecord &ir, const char* outFileName, const char* desc)
{
    referenceFileName = std::string( dr.giveReferenceName() ); 

    bool inputReaderFinish = true;

    this->coreOutputFileName = std::string( dataOutputFileName );
    this->dataOutputFileName = std::string( dataOutputFileName );

    if ( this->giveProblemMode() == _postProcessor ) {
        // modify output file name to prevent output to be lost
        this->dataOutputFileName.append( ".oofeg" );
    }

    if ( ( outputStream = fopen( this->dataOutputFileName.c_str(), "w" ) ) == NULL ) { OOFEM_ERROR( "Can't open output file %s", this->dataOutputFileName.c_str() ); }

    fprintf( outputStream, "%s", PRG_HEADER );
    this->startTime = time( NULL );
    fprintf( outputStream, "\nStarting analysis on: %s\n", ctime( &this->startTime ) );

    fprintf( outputStream, "%s\n", desc );

    // instanciate receiver
    this->initializeFrom( ir );

    exportModuleManager.initializeFrom( ir );
    // instanciate export module manager
    exportModuleManager.instanciateYourself( dr, ir );

    if ( inputReaderFinish ) { ir.finish(); }

    return 1;
}

void
DummyEnum::solveYourself()
{
    this->timer.startTimer( EngngModelTimer::EMTT_AnalysisTimer );
    this->terminate( NULL );
}
} // end namespace oofem
