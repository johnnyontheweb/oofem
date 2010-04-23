// New class created by mj
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2010   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

//
// class gpExportModule
//

#ifndef gpexportmodule_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "exportmodule.h"
#include "domain.h"
#include "engngm.h"

namespace oofem {

/**
 Represents GP (Gauss point) export module. 
 This module writes the coordinates of all Gauss points 
 along with the values of certain internal variables
 for further processing.
*/
class GPExportModule : public ExportModule
{
protected:

public:
 /// Constructor. Creates empty Output Manager. By default all components are selected.
 GPExportModule (EngngModel* e);

 /// Destructor
 ~GPExportModule();

 /// Initializes receiver acording to object description stored in input record
 virtual IRResultType initializeFrom (InputRecord* ir);
 /**
  Writes the output. Abstract service.
  @param tStep time step.
  */
 void              doOutput (TimeStep* tStep);
 /**
  Initializes receiver. 
  The init file messages should be printed.
  */
 void              initialize ();
 /**
  Terminates the receiver.
  The terminating messages should be printed.
  All the streams should be closed.
  */
 void              terminate ();
 /// Returns class name of the receiver.
 virtual const char* giveClassName () const { return "GPExportModuleClass" ;}

protected:
 /// returns the internal smoother
 NodalRecoveryModel* giveSmoother();

 /// returns the output stream for given solution step
 FILE* giveOutputStream (TimeStep*) ;
 /**
  Returns corresponding element cell_type.
  Some common element types are supported, others can be supported via interface concept.
  */
 int   giveCellType (Element*);
 /**
  Returns the number of elements vtk cells
  */
 int giveNumberOfElementCells (Element*);
 /**
  Returns number of nodes correpsonding to cell type
  */
 int giveNumberOfNodesPerCell (int cellType);
 /**
  Returns the element cell geometry.
  */
 void giveElementCell (IntArray& answer, Element* elem, int cell);
 /**
  export internal variables 
  */
 void exportIntVars (FILE* stream, TimeStep* tStep);
 /**
  export primary variables 
  */
 void exportPrimaryVars (FILE* stream, TimeStep* tStep);
 /** exports single variable */
 void exportIntVarAs (InternalStateType valID, InternalStateValueType type, FILE* stream, TimeStep* tStep);
 /** exports single variable */
 void exportPrimVarAs (UnknownType valID, FILE* stream, TimeStep* tStep);
 
};

#define gpexportmodule_h
#endif

} // namespace oofem

