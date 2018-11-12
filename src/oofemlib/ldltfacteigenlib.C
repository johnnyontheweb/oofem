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

#include "ldltfacteigenlib.h"
#include "classfactory.h"
#include "eigensolvermatrix.h"

namespace oofem {
REGISTER_SparseLinSolver(LDLTFactEigenLib, ST_EigenLib)

LDLTFactEigenLib :: LDLTFactEigenLib(Domain *d, EngngModel *m) :
    SparseLinearSystemNM(d, m)
{
}

LDLTFactEigenLib :: ~LDLTFactEigenLib()
{
}

NM_Status
LDLTFactEigenLib :: solve(SparseMtrx &A, FloatArray &b, FloatArray &x)
{
    // check whether Lhs supports factorization
    if ( !A.canBeFactorized() ) {
        OOFEM_ERROR("Lhs not support factorization");
    }

	EigenSolverMatrix* opA = dynamic_cast<EigenSolverMatrix*>(&A);
	Eigen::VectorXd opb(b.giveSize());
	for (int i = 0; i < b.giveSize(); ++i)
		opb[i] = b[i];

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(opA->giveEigenMatrix());
	Eigen::VectorXd opx = sparseSolver.solve(opb);     // solving

	x.resize(opx.size());
	for (int i = 0; i < x.giveSize(); ++i)
		x[i] = opx[i];

	if (sparseSolver.info() != Eigen::Success)
		return NM_NoSuccess;

    return NM_Success;
}
} // end namespace oofem
