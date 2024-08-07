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

#include "spoolessparsemtrx.h"
#include "engngm.h"
#include "sparsemtrxtype.h"
#include "classfactory.h"

namespace oofem {
REGISTER_SparseMtrx(SpoolesSparseMtrx, SMT_SpoolesMtrx);

void
SpoolesSparseMtrx :: times(const FloatArray &x, FloatArray &answer) const
{
    double alpha = 1.0, beta = 0.0;
    int result;
	double* temp = const_cast<double *>(x.givePointer());

    answer.resize( this->giveNumberOfColumns() );
    answer.zero();

    if ( sflag == SPOOLES_SYMMETRIC ) {
		result = InpMtx_sym_gmvm(this->mtrx, &beta, 1, answer.givePointer(), &alpha, 1, temp);
    } else if ( sflag == SPOOLES_NONSYMMETRIC ) {
        result = InpMtx_nonsym_gmvm( this->mtrx, & beta, 1, answer.givePointer(), & alpha, 1, temp);
    } else {
        OOFEM_ERROR("unsupported symmetry flag");
        exit(1);
    }
}

void
SpoolesSparseMtrx :: times(double x)
{
    OOFEM_ERROR("unsupported");
}

void
SpoolesSparseMtrx :: timesT(const FloatArray &x, FloatArray &answer) const
{
    double alpha = 1.0, beta = 0.0;
    int result;
	double* temp = const_cast<double *>(x.givePointer());

    answer.resize( this->giveNumberOfRows() );
    answer.zero();

    if ( sflag == SPOOLES_SYMMETRIC ) {
        result = InpMtx_sym_gmvm( this->mtrx, & beta, 1, answer.givePointer(), & alpha, 1, temp);
    } else if ( sflag == SPOOLES_NONSYMMETRIC ) {
        result = InpMtx_nonsym_gmvm_T( this->mtrx, & beta, 1, answer.givePointer(), & alpha, 1, temp);
    } else {
        OOFEM_ERROR("unsupported symmetry flag");
    }

    if ( result != 1 ) {
        OOFEM_ERROR("error code from InpMtx_(non)sym_gmvm %d", result);
    }
}

int
SpoolesSparseMtrx :: buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s)
{
    // Determine number of equations and estimate number of nonzero entries
    int neq = eModel->giveNumberOfDomainEquations(di, s);
    int nent = neq * 5;


    this->mtrx = InpMtx_new();
    InpMtx_init(this->mtrx, INPMTX_BY_ROWS, type, nent, neq);
    nRows = nColumns = neq;
    return true;
}

int
SpoolesSparseMtrx :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int ndofe = mat.giveNumberOfRows();

    for ( int i = 1; i <= ndofe; i++ ) {
        int ac1 = loc.at(i);
        if ( ac1 == 0 ) {
            continue;
        }

        for ( int j = 1; j <= ndofe; j++ ) {
            int ac2 = loc.at(j);
            if ( ac2 == 0 ) {
                continue;
            }

            if ( ac1 > ac2 ) {
                continue;
            }

            InpMtx_inputRealEntry( this->mtrx, ac1 - 1, ac2 - 1, mat.at(i, j) );
        }
    }

    this->version++;
    return 1;
}

int
SpoolesSparseMtrx :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
    int dim1 = mat.giveNumberOfRows();
    int dim2 = mat.giveNumberOfColumns();
    for ( int i = 1; i <= dim1; i++ ) {
        int ii = rloc.at(i);
        if ( ii ) {
            for ( int j = 1; j <= dim2; j++ ) {
                int jj = cloc.at(j);
                if ( jj ) {
                    InpMtx_inputRealEntry( this->mtrx, ii - 1, jj - 1, mat.at(i, j) );
                }
            }
        }
    }

    this->version++;

    return 1;
}

//void SpoolesSparseMtrx::add(double x, FloatMatrix &m)
//{
//	int dim1 = m.giveNumberOfRows();
//	int dim2 = m.giveNumberOfColumns();
//	int i, j;
//	for (i = 1; i <= dim1; i++) {
//	 for (j = 1; j <= dim2; j++)
//	 {
//		 //if (j > i) {
//			// continue; // symmetric lower triangular?
//		 //}
//		 InpMtx_inputRealEntry(this->mtrx, i - 1, j - 1, m.at(i, j));
//	 }
//	}
//	this->version++;
//}

void
SpoolesSparseMtrx :: zero()
{
    InpMtx_clearData(this->mtrx);
}

double &
SpoolesSparseMtrx :: at(int i, int j)
{
    OOFEM_ERROR("unsupported");
    abort();
}

double
SpoolesSparseMtrx :: at(int i, int j) const
{
    OOFEM_ERROR("unsupported");
    return 0.0;
}

void
SpoolesSparseMtrx :: printStatistics() const
{
    InpMtx_writeStats(this->mtrx, stdout);
}

void
SpoolesSparseMtrx :: printYourself() const
{
    InpMtx_writeForHumanEye(this->mtrx, stdout);
}
} // end namespace oofem
