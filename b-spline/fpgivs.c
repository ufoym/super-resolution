/* .\parcur\fpgivs.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Subroutine */ int fpgivs_(real *piv, real *ww, real *cos__, real *sin__)
{
    /* System generated locals */
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real dd, one, store;

/*  subroutine fpgivs calculates the parameters of a givens */
/*  transformation . */
/*  .. */
/*  ..scalar arguments.. */
/*  ..local scalars.. */
/*  ..function references.. */
/*  .. */
    one = 1.f;
    store = dabs(*piv);
    if (store >= *ww) {
/* Computing 2nd power */
	r__1 = *ww / *piv;
	dd = store * sqrt(one + r__1 * r__1);
    }
    if (store < *ww) {
/* Computing 2nd power */
	r__1 = *piv / *ww;
	dd = *ww * sqrt(one + r__1 * r__1);
    }
    *cos__ = *ww / dd;
    *sin__ = *piv / dd;
    *ww = dd;
    return 0;
} /* fpgivs_ */

