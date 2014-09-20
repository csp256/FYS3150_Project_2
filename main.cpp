#include <iostream>
#include <armadillo>
#include <math.h>
#include <unistd.h>

using namespace std;
using namespace arma;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// Implementation of page 24 of
// people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf
// For unknown reasons, additional passes through this matrix decrease accuracy

// How is that for a descriptive function name???
colvec symmetricTridiagonalFrancisAlgorithmWithWilkinsonShift(colvec a, colvec b) {
    float epsilon = 90.001f;
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
                // Diagonalizing a symmetric 2x2 matrix
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
    int n = 10;
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
    colvec l = symmetricTridiagonalFrancisAlgorithmWithWilkinsonShift(d, o); // myEigenSolver(d);
    l = sort(l);
    cout << "My algorithm" << endl << l.t() << endl << endl;
    colvec diff = abs(l - eValues);
    cout << "Absolute difference" << endl << diff.t() << endl;
    cout << "Relative difference" << endl << (diff.t() / eValues.t()) << endl;
    cout << endl << endl;

    return 0;
}


