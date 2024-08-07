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

// #include <../unsupported/Eigen/ArpackSupport>
#include <../Eigen/Core>
#include <../Eigen/SparseCore>
#include <../Eigen/SparseCholesky>
#include <SymGEigsSolver.h>
#include <MatOp/SparseGenMatProd.h>
#include <MatOp/SparseCholesky.h>
#include "eigensolvermatrix.h"
#include "eigensolver.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "sparsemtrx.h"
#include "mathfem.h"
#include "sparselinsystemnm.h"
#include "classfactory.h"
#include "domain.h"
#include "engngm.h"
#include <memory>

// using namespace Spectra;

namespace oofem {
REGISTER_GeneralizedEigenValueSolver( EigenSolver, GES_Eigen );


EigenSolver :: EigenSolver(Domain *d, EngngModel *m) :
    SparseGeneralEigenValueSystemNM(d, m)
{
    nitem = 100; // max number of iterations
}


EigenSolver :: ~EigenSolver() { }

ConvergedReason
EigenSolver :: solve(SparseMtrx &a, SparseMtrx &b, FloatArray &_eigv, FloatMatrix &_r, double rtol, int nroot)
{
	FILE *outStream;
	int size;

	outStream = domain->giveEngngModel()->giveOutputStream();

	// first check whether Lhs is defined

	if (a.giveNumberOfRows() != a.giveNumberOfColumns() ||
		b.giveNumberOfRows() != b.giveNumberOfRows() ||
		a.giveNumberOfColumns() != b.giveNumberOfColumns()) {
		OOFEM_ERROR("matrices size mismatch");
	}

	if (nroot <= 0) {
		OOFEM_ERROR("Invalid number of eigenvalues, %d", nroot);
	}

	// spectra supports n-1 eigenvalues
	if (nroot >= a.giveNumberOfRows()) {
		nroot = a.giveNumberOfRows() - 1;
		OOFEM_WARNING("The number of requested eigenvalues exceeds the problem size. The number is reduced down to %d", nroot);
	}

	//typedef Eigen::SparseMatrix<double,0,int> SparseMat;
	//typedef Eigen::SimplicialLDLT<SparseMat> SparseChol;
	//typedef Eigen::ArpackGeneralizedSelfAdjointEigenSolver <SparseMat, SparseChol> Arpack;
	//Arpack arpack;
	//// define sparse matrix A
	EigenSolverMatrix* A = dynamic_cast<EigenSolverMatrix*>(&a);
	EigenSolverMatrix* B = dynamic_cast<EigenSolverMatrix*>(&b);
	if (!A || !B)
		OOFEM_ERROR("Error casting matrices");
	////...
	//// calculate the two smallest eigenvalues
	//arpack.compute(*A, *B, nroot, "SM");

	Spectra::SparseSymMatProd<double> opB(B->giveEigenMatrix());
#if 1
	// cholesky
	Spectra::SparseCholesky<double> op(A->giveEigenMatrix());
    if ( op.info() != Spectra::SUCCESSFUL ) {
        fprintf( outStream, "Eigen-Spectra :: Stiffness matrix is singular.\n" );
        OOFEM_LOG_INFO( "Eigen-Spectra :: Stiffness matrix is singular.\n" );
        return CR_FAILED;
    }
	// Construct eigen solver object, requesting the largest three eigenvalues
	Spectra::SymGEigsSolver< double, Spectra::LARGEST_MAGN, Spectra::SparseSymMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEIGS_CHOLESKY  > eigs(&opB, &op, nroot, min(2 * nroot, a.giveNumberOfColumns()));
#else
	// with regular inverse
    Spectra::SparseRegularInverse<double> op( A->giveEigenMatrix() );
    // Construct eigen solver object, requesting the largest three eigenvalues
    Spectra::SymGEigsSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseSymMatProd<double>, Spectra::SparseRegularInverse<double>, Spectra::GEIGS_REGULAR_INVERSE> eigs( &opB, &op, nroot, min( 2 * nroot, a.giveNumberOfColumns() ) );
#endif
	// Initialize and compute
	eigs.init();
	int nconv = eigs.compute(1000, rtol);
	// Retrieve results
	Eigen::VectorXcd evalues;
	//if (eigs.info() == Spectra::SUCCESSFUL) // always copy
	evalues = eigs.eigenvalues();
	//_eigv.resize(evalues.size());
	_eigv.resize(nroot); _eigv.zero(); // return zero if not converged for all

	for (int i=0;i<evalues.size();++i)
		_eigv.at(i+1)=1/(evalues.coeff(i).real());

	Eigen::MatrixXcd evectors = eigs.eigenvectors();
	_r.resize(evectors.rows(), evectors.cols());
	for (int i = 0; i < evectors.rows(); ++i)
		for (int j = 0; j < evectors.cols(); ++j){
			_r.at(i + 1, j + 1) = evectors.coeff(i, j).real();
		}

#ifdef TIME_REPORT
	timer.stopTimer();
	OOFEM_LOG_INFO("EigenSolver info: user time consumed by solution: %.2fs\n", timer.getUtime());
#endif
	if (eigs.info() == Spectra::SUCCESSFUL) {
		fprintf(outStream, "Eigen-Spectra :: convergence reached\n");
	} else {
		fprintf(outStream, "Eigen-Spectra :: convergence not reached\n");
	}
	return CR_CONVERGED;
}
} // end namespace oofem