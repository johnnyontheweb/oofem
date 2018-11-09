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

#include "eigensolvermatrix.h"
#include <Eigen/Sparse>
#include "error.h"
#include "engngm.h"
#include "domain.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "sparsemtrxtype.h"
#include "classfactory.h"
#include "activebc.h"
#include "eigensolver.h"
#include "unknownnumberingscheme.h"

#include <set>

namespace oofem {

REGISTER_SparseMtrx( EigenSolverMatrix, SMT_EigenSparse);

//typedef Eigen::SparseMatrix<double, 0, int> SparseMat;

EigenSolverMatrix::EigenSolverMatrix(int n) : SparseMtrx(n, n)
{
}

EigenSolverMatrix :: ~EigenSolverMatrix()
{
}

/*****************************/
/*  Copy constructor         */
/*****************************/

EigenSolverMatrix::EigenSolverMatrix(const EigenSolverMatrix &S) : SparseMtrx(S.nRows, S.nColumns)
{
	eigenMatrix.reset(new Eigen::SparseMatrix<double>(S.giveEigenMatrix()));
    OOFEM_ERROR("not implemented");
}

SparseMtrx *EigenSolverMatrix :: GiveCopy() const
{
    OOFEM_ERROR("not implemented");
    return NULL;
}

void EigenSolverMatrix :: times(const FloatArray &x, FloatArray &answer) const
{
  // Note: not really efficient. The sparse matrix is assembled directly into its block structure,
  // which is efficient for factorization, but unfortunately not efficient for implementing the multiplication,
  // as the blocks have to be identified (see implementation of ElementAt method) when traversing rows
  // Also note, that this method will yield correct results only before factorization, after that the blocks
  // contain factorized matrix.

  int i, j, dim;

  dim = this->nColumns;

  answer.resize(dim);
  answer.zero();

  for (i = 1; i <= dim; i++) {
	  for (j = 1; j <= dim; j++)
	  {
		  answer.at(i) += this->eigenMatrix->coeff(i - 1, j - 1) * x.at(j);
	  }
  }
}

void EigenSolverMatrix :: times(double x)
{
	(*this->eigenMatrix) *= x;
}

int EigenSolverMatrix :: buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s)
{
	int neq;
	if (s.isDefault()) {
		neq = eModel->giveNumberOfDomainEquations(di, s);
	}
	else {
		neq = s.giveRequiredNumberOfDomainEquation();
	}
	nRows = nColumns = neq;
	eigenMatrix.reset(new Eigen::SparseMatrix<double>(nRows, nColumns));
	return true;
}

int EigenSolverMatrix :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int i, j, ii, jj, dim;

 #  ifdef DEBUG
    dim = mat.giveNumberOfRows();
    if ( dim != loc.giveSize() ) {
        OOFEM_ERROR("dimension of 'k' and 'loc' mismatch");
    }

 #  endif

    dim = mat.giveNumberOfRows();

	for (j = 1; j <= dim; j++) {
		jj = loc.at(j);
		if (jj) {
			for (i = 1; i <= dim; i++) {
				ii = loc.at(i);
				if (ii) {
					this->eigenMatrix->coeffRef(ii - 1, jj - 1) += mat.at(i, j);
				}
			}
		}
	}

	// increment version
	this->version++;

	return 1;
}

int EigenSolverMatrix :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
    int dim1, dim2;

    // this->checkSizeTowards(rloc, cloc);

    dim1 = mat.giveNumberOfRows();
    dim2 = mat.giveNumberOfColumns();
	for (int i = 1; i <= dim1; i++) {
		int ii = rloc.at(i);
		if (ii) {
			for (int j = 1; j <= dim2; j++) {
				int jj = cloc.at(j);
				if (jj) {
					this->eigenMatrix->coeffRef(ii - 1, jj - 1) += mat.at(i, j);
				}
			}
		}
	}

    // increment version
    this->version++;

    return 1;
}

void EigenSolverMatrix :: zero()
{
	this->eigenMatrix->data().clear(); // don't know if this works as intended

    // increment version
    this->version++;
}

/*********************/
/*   Array access    */
/*********************/

double &EigenSolverMatrix :: at(int i, int j)
{
    // increment version
    this->version++;
	return this->eigenMatrix->coeffRef(i - 1, j - 1);
}


double EigenSolverMatrix :: at(int i, int j) const
{
	return this->eigenMatrix->coeff(i - 1, j - 1);
}

double EigenSolverMatrix :: operator() (int i, int j)  const
{
	return this->eigenMatrix->coeff(i, j);
}

double &EigenSolverMatrix :: operator() (int i, int j)
{
    // increment version
    this->version++;
	return this->eigenMatrix->coeffRef(i, j);
}
} // end namespace oofem

