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
#include <MatOp/SparseGenMatProd.h>

#include <set>

namespace oofem {

REGISTER_SparseMtrx( EigenSolverMatrix, SMT_EigenSparse);

//typedef Eigen::SparseMatrix<double, 0, int> SparseMat;

EigenSolverMatrix::EigenSolverMatrix(int n) : SparseMtrx(n, n)
{ }

EigenSolverMatrix :: ~EigenSolverMatrix()
{
}

/*****************************/
/*  Copy constructor         */
/*****************************/

EigenSolverMatrix::EigenSolverMatrix(const EigenSolverMatrix &S) : SparseMtrx(S.nRows, S.nColumns)
{
	eigenMatrix.reset(new Eigen::SparseMatrix<double>(S.giveEigenMatrix()));
	tripletList = S.tripletList;
//    OOFEM_ERROR("not implemented");
}

SparseMtrx *EigenSolverMatrix :: GiveCopy() const
{
	EigenSolverMatrix *answer;

	answer = new EigenSolverMatrix(*this);
	return answer;
}

void EigenSolverMatrix :: times(const FloatArray &x, FloatArray &answer) const
{
	//int i, j, dim;
	int dim = this->nColumns;

	answer.resize(dim);
	answer.zero();

	Spectra::SparseGenMatProd<double> op(*(eigenMatrix.get()));
	op.perform_op(x.givePointer(), answer.givePointer());

	//for (i = 1; i <= dim; i++) {
	// for (j = 1; j <= dim; j++)
	// {
	//  answer.at(i) += this->eigenMatrix->coeff(i - 1, j - 1) * x.at(j);
	// }
	//}
}

void EigenSolverMatrix :: times(double x)
{
	(*this->eigenMatrix) *= x;
}

void EigenSolverMatrix::add(double x, SparseMtrx &m)
{
	// for square matrices
	if (this->nColumns != m.giveNumberOfColumns()) {
		OOFEM_ERROR("dimension of 'k' and 'm' mismatch");
	}

	//for (i = 1; i <= dime; i++) {
	// for (j = 1; j <= dime; j++)
	// {
	//	 if (j > i) {
	//		 continue; // symmetric lower triangular
	//	 }
	//	 //this->eigenMatrix->coeffRef(i - 1, j - 1) += m.at(i, j) * x;
	//	 tripletList.push_back(Eigen::Triplet<double>(i - 1, j - 1, this->eigenMatrix->coeff(i - 1, j - 1) +  m.at(i, j) * x));
	// }
	//}
	//applyTriplets();

	this->version++;
	EigenSolverMatrix& temp = static_cast<EigenSolverMatrix&>(m);
	eigenMatrix->operator+=(temp.eigenMatrix->operator*(x));
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
	tripletList.reserve(neq / 20); //5% rough guess
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
					if (jj > ii) {
						continue; // symmetric
					}
					//this->eigenMatrix->coeffRef(ii - 1, jj - 1) += mat.at(i, j);
					tripletList.push_back(Eigen::Triplet<double>(ii - 1, jj - 1, mat.at(i, j)));
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
				if (jj && (jj <= ii) ) {
					//this->eigenMatrix->coeffRef(ii - 1, jj - 1) += mat.at(i, j);
					tripletList.push_back(Eigen::Triplet<double>(ii - 1, jj - 1, mat.at(i, j)));
				}
			}
		}
	}

    // increment version
    this->version++;

    return 1;
}

void EigenSolverMatrix :: applyTriplets()
{
	eigenMatrix->setFromTriplets(tripletList.begin(), tripletList.end());
	tripletList.clear();
}

int EigenSolverMatrix :: assembleEnd()
{
	applyTriplets();
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

void EigenSolverMatrix :: writeToFile(const char *fname) const
{
	FILE *file = fopen(fname, "w");
	for (int i = 1; i <= nRows; ++i) {
		for (int j = 1; j <= nColumns; ++j) {
			fprintf(file, "%10.3e  ", this->at(i, j));
		}
		fprintf(file, "\n");
	}
	fclose(file);
}


} // end namespace oofem

