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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "inverseit.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "sparsemtrx.h"
#include "mathfem.h"
#include "sparselinsystemnm.h"
#include "classfactory.h"
#include "domain.h"
#include "engngm.h"
#include <memory>

namespace oofem {
REGISTER_GeneralizedEigenValueSolver(InverseIteration, GES_InverseIt);


InverseIteration :: InverseIteration(Domain *d, EngngModel *m) :
    SparseGeneralEigenValueSystemNM(d, m),
    nitem(100)
{
}


ConvergedReason
InverseIteration :: solve(SparseMtrx &a, SparseMtrx &b, FloatArray &eigv, FloatMatrix &r, double rtol, int nroot)
{
	FILE *outStream;
    if ( a.giveNumberOfColumns() != b.giveNumberOfColumns() ) {
        OOFEM_ERROR("matrices size mismatch");
    }

    auto solver = GiveClassFactory().createSparseLinSolver(ST_Direct, domain, engngModel);

    int nn = a.giveNumberOfColumns();

    if (nroot <= 0) {
        OOFEM_ERROR("Invalid number of eigenvalues, %d", nroot);
    }

    // limit to n eigenvalues
    if (nroot > nn) {
        nroot = nn;
        OOFEM_WARNING("The number of requested eigenvalues exceeds the problem size. The number is reduced down to %d", nroot);
    }

    int nc = min(2 * nroot, nroot + 8);
    nc = min(nc, nn);

    //// control of diagonal zeroes in mass matrix, to be avoided
    //int i;
    //for (i = 1; i <= nn; i++) {
    //	if (b.at(i, i) == 0) {
    //		b.at(i, i) = 1.0e-12;
    //	}
    //}

    FloatArray w(nc), ww(nc), t;
    std :: vector< FloatArray > z(nc, nn), zz(nc, nn), x(nc, nn);
    outStream = domain->giveEngngModel()->giveOutputStream();

    /*  initial setting  */
#if 0
    ww.add(1.0);
    for ( int j = 0; j < nc; j++ ) {
        z[j].add(1.0);
    }
#else
    {
        FloatArray ad(nn), bd(nn);
        for ( int i = 1; i <= nn; i++ ) {
            ad.at(i) = fabs(a.at(i, i));
            bd.at(i) = fabs(b.at(i, i));
        }
        IntArray order;
        order.enumerate(nn);
        std :: sort(order.begin(), order.end(), [&ad, &bd](int a, int b) { return bd.at(a) * ad.at(b) > bd.at(b) * ad.at(a); });
        for ( int i = 0; i < nc; i++ ) {
            x[i].at(order[i]) = 1.0;
            b.times(x[i], z[i]);
            ww.at(i + 1) = z[i].dotProduct(x[i]);
        }
    }
#endif

    int it;
    for ( it = 0; it < nitem; it++ ) {
        /*  copy zz=z  */
        for ( int j = 0; j < nc; j++ ) {
            zz[j] = z[j];
        }

        /*  solve matrix equation K.X = M.X  */
        for ( int j = 0; j < nc; j++ ) {
            solver->solve(a, z[j], x[j]);
        }

        /*  evaluation of Rayleigh quotients  */
        for ( int j = 0; j < nc; j++ ) {
            w.at(j + 1) = zz[j].dotProduct(x[j]);
        }

        for ( int j = 0; j < nc; j++ ) {
            b.times(x[j], z[j]);
        }

        for ( int j = 0; j < nc; j++ ) {
            w.at(j + 1) /= z[j].dotProduct(x[j]);
        }

        /*  check convergence  */
        int ac = 0;
        for ( int j = 1; j <= nc; j++ ) {
            if ( fabs( ww.at(j) - w.at(j) ) <= fabs( w.at(j) * rtol ) ) {
                ac++;
            }

            ww.at(j) = w.at(j);
        }

        //printf ("\n iteration  %d   %d",it,ac);
        //w.printYourself();

        /*  Gramm-Schmidt ortogonalization   */
        for ( int j = 0; j < nc; j++ ) {
            if ( j != 0 ) {
                b.times(x[j], t);
            }

            for ( int ii = 0; ii < j; ii++ ) {
                x[j].add( -x[ii].dotProduct(t), x[ii] );
            }

            b.times(x[j], t);
            x[j].times( 1.0 / sqrt( x[j].dotProduct(t) ) );
        }

        if ( ac > nroot ) {
            break;
        }

        /*  compute new approximation of Z  */
        for ( int j = 0; j < nc; j++ ) {
            b.times(x[j], z[j]);
        }
    }

    // copy results
    IntArray order;
    order.enumerate(w.giveSize());
    std :: sort(order.begin(), order.end(), [&w](int a, int b) { return w.at(a) < w.at(b); });

    eigv.resize(nroot);
    r.resize(nn, nroot);
    for ( int i = 1; i <= nroot; i++ ) {
        eigv.at(i) = w.at(order.at(i));
        r.setColumn(x[order.at(i) - 1], i);
    }

    if ( it < nitem ) {
	    fprintf(outStream, "InverseIteration :: convergence reached in %d iterations\n", it);
    } else {
        OOFEM_WARNING("convergence not reached after %d iterations\n", it);
        return CR_DIVERGED_ITS;
    }

    return CR_CONVERGED;
}
} // end namespace oofem