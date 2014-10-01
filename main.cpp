#include <iostream>
#include <armadillo>
#include <math.h>
#include <unistd.h>
#include <limits>
#include <cfloat>

using namespace std;
using namespace arma;

// Port of the Fortran 77 LAPACK implementation of SSTERF
// Eigenvalues-only is assumed.
/*
colvec SSTERF (colvec D, colvec E, int MAXIT) {
*/

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// This is the tqlrat function from LAPACK ported to C++ Armadillo
// *** It is destructive to its arguments ***
//

int modified_tqlrat ( colvec diagonal, colvec offDiagonal ) {
    const int n = diagonal.n_elem; // n is the order of the matrix. at least this one deserves to be one letter...
    const int MAX_ITERATIONS = 30; // As it is in LAPACK, so it shall be in my code.

    // Forgive FORTRAN, for it has sinned
    // The same variable "name" is used multiple times in destinct parts of the code
    // to represent logically disconnected things! I am trying to fix this but it is
    // an uphill battle...
    int errorCode, indexOffset, j, eigenvalueIndex, eigenvalueIndexPlusOne, m, mMinusEigenvalueIndex;
    double b, c, f, g, h;
    double p, r, s, t;

    errorCode = 0;
    if ( n == 1 ) { // Check for degenerate matrix --> easy early abort
        return errorCode;
    }

    f = 0.0;
    t = 0.0;
    offDiagonal(n-1) = 0.0;
    for (eigenvalueIndex=0; eigenvalueIndex<n; eigenvalueIndex++) {
        iterationsRemaining = MAX_ITERATIONS;
        h = abs(diagonal(eigenvalueIndex)) + sqrt(offDiagonal(eigenvalueIndex));

        if ( t <= h ) {
            t = h;
            b = fabs (t) * DBL_EPSILON;
            c = b * b;
        }

        //  Look for small squared sub-diagonal element.
        for ( m = eigenvalueIndex; m < n; m++ ) {
            if ( offDiagonal(m) <= c ) {
                break; // m now holds the desired offdiagonal element
            }
        }

        if ( m != eigenvalueIndex ) { // There is simply nothing to be done if m's offdiagonal element is already less than 'c'
            // there are four places we can explicitly break out of the infinite loop below...
            for (;;) {
                if ( iterationsRemaining <= 0 ) { // Failure to converge :(
                    errorCode = eigenvalueIndex + 1;
                    return errorCode;
                }
                iterationsRemaining--;

                //  Form shift.
                eigenvalueIndexPlusOne = eigenvalueIndex + 1;
                s = sqrt ( offDiagonal(eigenvalueIndex) );
                g = diagonal(eigenvalueIndex);
                p = ( diagonal(eigenvalueIndexPlusOne) - g ) / ( 2.0 * s );
                r = hypot ( p, 1.0 );
                diagonal(eigenvalueIndex) = s / ( p + fabs ( r ) * sgn ( p ) );
                h = g - diagonal(eigenvalueIndex);
                for ( int i = eigenvalueIndexPlusOne; i < n; i++ ) {
                    diagonal(i) -= h; //diagonal(i) - h;
                }
                f += h;

                //  Rational QL transformation.
                g = diagonal(m);
                if ( g == 0.0 ) {
                    g = b;
                }

                h = g;
                s = 0.0;
                mMinusEigenvalueIndex = m - eigenvalueIndex;
                for ( indexOffset = 1; indexOffset <= mMinusEigenvalueIndex; indexOffset++ ) { // Sweeps 'i' from 'm-1' up to and including 'eigenvalueIndex'
                    int i = m - indexOffset;
                    p = g * h;
                    r = p + offDiagonal(i);
                    offDiagonal(i+1) = s * r;
                    s = offDiagonal(i) / r;
                    diagonal(i+1) = h + s * ( h + diagonal(i) );
                    g = diagonal(i) - offDiagonal(i) / g;
                    if (g == 0.0) {
                        g = b; // b is very tiny.
                    }
                    h = g * p / r;
                }
                offDiagonal(eigenvalueIndex) = s * g;
                diagonal(eigenvalueIndex) = h;
                //  Guard against underflow in convergence tests.
                if (h == 0.0) {
                    break;
                }
                //  (the aforementioned) convergence tests
                if ( abs(offDiagonal(eigenvalueIndex)) <= abs(c/h) ) {
                    break;
                }
                offDiagonal(eigenvalueIndex) *= h;
                if ( offDiagonal(eigenvalueIndex) == 0.0 ) {
                    break;
                }
            }
        }
        p = diagonal(eigenvalueIndex) + f;
    }

    return errorCode;
}








// Implementation of page 24 of
// people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf
// For unknown reasons, additional passes through this matrix decrease accuracy
// Interestingly, the offdiagonal is not converging to 0 but 1!
// The third element in the offdiagonal is large.

// How is that for a descriptive function name???
colvec symmetricTridiagonalFrancisAlgorithmWithWilkinsonShift(colvec a, colvec b) {
    float epsilon = 0.001f; // turned off currently
    int passes = 3;
    int pass = 0;
    unsigned long long int count = 0;

    int m = a.n_elem - 1;
    while (0<m) {
        ++pass;
        float s, c;
        float d = (a(m-1) - a(m)) * 0.5f;
        if (d == 0) {
            s = a(m) - abs(b(m));
        } else {
            s = a(m) - pow(b(m),2) / (d + sgn(d)*sqrt(d*d + b(m)*b(m)));
        }
        float x = a(0) - s;
        float y = b(1);
        for (int k=0; k<=m-1; ++k) {
            if (1 < m) {
                // givens rotation
                float r = hypot(y,x);
                c = x/r;
                s = y/r;
            } else {
                // Orthogonal diagonalization of a symmetric 2x2 matrix
                // http://scipp.ucsc.edu/~haber/ph116A/diag2x2_11.pdf
                float theta;
                theta = 0.5f * atan2(2*b(1),a(0)-a(1));
                s = sin(theta);
                c = cos(theta);
            }
            float w = c*x - s*y;
            d = a(k) - a(k+1);
            float z = (2*c*b(k+1) + d*s)*s;
            a(k)   -= z;
            a(k+1) += z;
            b(k+1) = d*c*s + (c*c - s*s)*b(k);
            x = b(k+1);
            if (0 < k) {
                b(k) = w;
            }
            if (k < m-1) {
                y = -s*b(k+2);
                b(k+2) *= c;
                if (abs(c) > 1) {
                    cout << c << endl;
                }
            }
            // Do not need to calculate Q. Do not need eigenvectors.
        }
        count += m;
//        cout << abs(b(m)) << "       " << epsilon*(abs(a(m-1)) + abs(a(m))) << endl;
        float smallThing = epsilon*(abs(a(m-1)) + abs(a(m)));
        if (isnan(smallThing)) {
            cout << "Not a number!" << endl;
        }
        if (abs(b(m)) < smallThing) { // check for convergence
//            --m;
//            temporarily turned off, to make iterations until m is decreased constant
        }
        if (pass == passes) {
            --m;
            pass = 0;
        }
    }
//    cout << "diagonal" << endl << a.t() << endl;
    cout << "off diagonal" << endl << b.t() << endl;
    cout << "inner loop iterations : " << count << endl; // total number of inner loop iterations performed
    return a;
}


// Finds eigenvalues of a tridiagonal matrix where the offdiagonals are identically 1
// no convergence check = unacceptable

colvec myEigenSolver(colvec diag) {
    cout << "myEigenSolver:" << endl;
    colvec l(diag.n_elem); // l := lambda, that is, eigenvalues
    float a;
    float theta;
    float s = 0.0f;
    float c = 1.0f;
    for (int i=0; i<diag.n_elem; ++i) {
        a = -s + c*diag(i);
        theta = atan(1/a);
        s = sin(theta);
        c = cos(theta);
        l(i) = 1.0f/s;
    }
    l(diag.n_elem-1) = a;
    l = sort(l);
    cout << l.t();
    return l;
}

// Fill diagonals

colvec fillDiagonals(int n, float h4) {
    colvec a(n);
    for (int i=0, j=1; i<n; ++i, ++j) {
        a(i) = (2.0f+(float)(j*j)*h4);
    }
    return a;
}

int main() {
    int n = 120;
    colvec d = fillDiagonals(n, 10);
    colvec o(n-1);
    mat A(n,n);
    A.diag() = d;
    A.diag(-1) = o.ones();
    A.diag( 1) = o.ones();
    mat eVectors(n,n);
    colvec eValues(n);
    eig_sym(eValues, eVectors, A);
    eValues = sort(eValues);
    cout << "Armadillo" << endl << eValues.t() << endl;
    o = o.ones(n);
    modified_tqlrat(d, o);
    colvec l = d;
//    colvec l = symmetricTridiagonalFrancisAlgorithmWithWilkinsonShift(d, o); //
//    colvec l = myEigenSolver(d);
    //l = sort(l);
    cout << "My algorithm" << endl << l.t() << endl << endl;
    colvec diff = abs(l - eValues);
    cout << "Absolute difference" << endl << diff.t() << endl;
    cout << "Relative difference" << endl << (diff.t() / eValues.t()) << endl;
    cout << endl << endl;

//    double*

    return 0;
}


