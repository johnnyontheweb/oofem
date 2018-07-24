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

#ifndef eigensolvermatrix_h
#define eigensolvermatrix_h

#include <iostream>
#include <../unsupported/Eigen/ArpackSupport>
#include <../Eigen/SparseCore>
#include <../Eigen/SparseCholesky>
#include "sparsemtrx.h"
#include "intarray.h"
#include "floatarray.h"

#include <memory>

/* DSS module lives in global namespace, not in oofem namespace */
class EigenSolver;

namespace oofem {
/**
 * Interface to Direct Sparse Solver written by R.Vonracek.
 * This class represent the sparse matrix interface to DSS library. It allows to build internal structure,
 * assemble the DSS sparse matrix, and to factorize and back substitution operations.
 */
	class OOFEM_EXPORT EigenSolverMatrix : public SparseMtrx, public Eigen::SparseMatrix<double>
{
public:

protected:
	
public:
    /**
     * Constructor. Before any operation an internal profile must be built.
     * @param t Storage type
     * @param n Size of row and columns of square matrix
     * @see buildInternalStructure
     */
    EigenSolverMatrix(int n);
    /// Copy constructor
    EigenSolverMatrix(const EigenSolverMatrix &S);
    /// Destructor
    virtual ~EigenSolverMatrix();

    // Overloaded methods
    virtual SparseMtrx *GiveCopy() const;
    virtual void times(const FloatArray &x, FloatArray &answer) const;
    virtual void times(double x);
	virtual int buildInternalStructure(EngngModel *, int, const UnknownNumberingScheme & s) {};
    virtual int assemble(const IntArray &loc, const FloatMatrix &mat);
    virtual int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);
    virtual bool canBeFactorized() const { return false; }
    virtual void zero();
    virtual double &at(int i, int j);
    virtual double at(int i, int j) const;
	virtual SparseMtrxType giveType() const { return SMT_EigenSparse; }
    virtual bool isAsymmetric() const { return false; }

    virtual const char *giveClassName() const { return "EigenSolverMatrix"; }

	/// implements 0-based access
	double operator() (int i, int j) const;
	/// implements 0-based access
	double &operator() (int i, int j);
};

} // end namespace oofem

#endif // eigensolvermatrix_h

