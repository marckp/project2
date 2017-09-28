/*
  Solves the one-particle Schrodinger equation
  for a potential specified in function
  potential(). This example is for the harmonic oscillator in 3d
*/
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
//#define CATCH_CONFIG_MAIN
//#include <catch.hpp>

using namespace  std;
using namespace  arma;

// Forward Function declarations
double potential(double); // potential for interacting case
void output(double, double, int, vec& ); // Morten's code for output
void Jacobi_rotate ( mat& A, mat& R, int k, int l, int n); // actual Jacobi matrix rotation
double offdiag(mat& A, int *p, int *q, int n); // find the largest off diagonal matrix element
void initialize_non_interacting_case(mat& A, double rho_0, double rho_n, int n); // Initialize the matrix for non-interaction case
void initialize_interacting_case(mat& A, double rho_0, double rho_n, int n); // Initialize the matrix for interaction case
void diagonalize_matrix_jacobi(mat& A, mat& R, double tolerance, double iterations, double maxmondiag, int maxiter, int n );

// object for output files
ofstream ofile;

// Test cases for the following functions
// offdiag
// Jacobi_rotate
// Note: The main section must be commented out to perform test
//   The unit tests could be done in a separate file to make the process easier
//double test_offdiag_func(double testentry, int n, int i, int j)
//{

//    mat T(n,n);
//    T.zeros();
//    T(i,j)              = testentry;
//    T(j,i)              = testentry;
//    double outputval    = offdiag(T,&i,&j,n);
//    cout << "testentry: " << testentry << endl;
//    cout << "outputval: " << outputval << endl;
//    return outputval;
//}

//TEST_CASE("Test of off diagonal", "[test_offdiag_func]") {
//    REQUIRE(test_offdiag_func(100.0,10,2,1)==100);
//    REQUIRE(test_offdiag_func(4.50,10,2,1)==4.50);
//    REQUIRE(test_offdiag_func(56.0,5,2,1)==28.0);
//};

// Begin of main program

int main(int argc, char* argv[])
{
    //  we have defined a matrix A and a matrix R for the eigenvector, both of dim n x n
    //  The final matrix R has the eigenvectors in its row elements, it is set to one
    //  for the diagonal elements in the beginning, zero else.

    // size of the square matrix
    int n;

    // Output file stem
    string filename;


    // We read also the stem name for the output file and the dimension of the input matrix A we want
    if( argc <= 1 ){
          cout << "Bad Usage: " << argv[0] <<
              " read also file name on same line and max power 10^n" << endl;
    }
    else
    {
        filename    = argv[1]; // first command line argument after name of program
        n           = atoi(argv[2]);
        cout << "Output filename stem: " << filename << endl;
        cout << "size of matrix under consideration: " << n << endl;
    }

    // Initialize A for the non-interacting case
    mat A(n,n);
    A.fill(0.0);
    initialize_interacting_case(A,0,1,n);
    mat R(n,n);
    R.zeros();

    cout << "A= " << A << endl;
    cout << "R= " << R << endl;

    double  tolerance   = 1.0E-10;
    int     iterations  = 0;
    double  maxnondiag  = 1.0E10;
    int     maxiter     = 100;
    diagonalize_matrix_jacobi(A,R,tolerance,iterations,maxnondiag,maxiter,n);

//    while ( maxnondiag > tolerance && iterations <= maxiter)
//    {
//       int p, q;
//       maxnondiag   = offdiag(A, &p, &q, n);
//       // cout << "start state: " << p << " : " << q << endl;
//       // cout << "max non-diagonal element: " << maxnondiag << endl;
//       Jacobi_rotate(A, R, p, q, n);
//       iterations++;
//       // cout << "current iteration: " << iterations << endl;
//    }

    cout << "A= " << A << endl;
    cout << "R= " << R << endl;

}  //  end of main function

// Initialize A with tri-diag form modelling 1e no interaction
void initialize_non_interacting_case(mat& A, double rho_0, double rho_n, int n)
{
    // compute the step size
    double h = (rho_n - rho_0)/n;

    // initialization loop puts values in the diagonal and upper and lower diagonl elements
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i==j)
            {
                // diagonal elements
                A(i,j)  = 2/h*h;
            }
            else if (j==i+1)
            {
                // Upper diagonal elements
                A(i,j)  = -1/h*h;
            }
            else if (j==i-1)
            {
                // Lower diagonal elements
                A(i,j)  = -1/h*h;
            }
        }
    }
}

// Initialize A with tri-diag form modelling the interacting case (this may be wrong)
void initialize_interacting_case(mat& A, double rho_0, double rho_n, int n)
{
    // compute the step size
    // compute the position in the grid
    double h = (rho_n - rho_0)/n;
    double x;

    for (int i = 0; i < n; i++)
    {
        x = i * h;
        for (int j = 0; j < n; j++)
        {
            if (i==j)
            {
                // diagonal elements with potential as a function of x added
                A(i,j)  = 2/h*h + potential(x);
            }
            else if (j==i+1)
            {
                // Upper diagonal elements
                A(i,j)  = -1/h*h;
            }
            else if (j==i-1)
            {
                // Lower diagonal elements
                A(i,j)  = -1/h*h;
            }
        }
    }
}

/*
  The function potential()
  calculates and return the value of the
  potential for a given argument x.
  The potential here is for the hydrogen atom
*/
double potential(double x)
{
  return x*x;

} // End: function potential()


void output(double RMin , double RMax, int Dim, vec& d)
{
  int i;
  cout << "RESULTS:" << endl;
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout <<"Rmin = " << setw(15) << setprecision(8) << RMin << endl;
  cout <<"Rmax = " << setw(15) << setprecision(8) << RMax << endl;
  cout <<"Number of steps = " << setw(15) << Dim << endl;
  cout << "Five lowest eigenvalues:" << endl;
  for(i = 0; i < 5; i++)
  {
    cout << setw(15) << setprecision(8) << d[i] << endl;
  }
}  // end of function output

//  the offdiag function, using Armadillo find the largest off diagonal entry of the input matrix
double offdiag(mat& A, int *p, int *q, int n)
{
    // output max temporary processing matrix T
   double max;
   for (int i = 0; i < n; ++i)
   {
       for ( int j = i+1; j < n; ++j)
       {
           double aij = fabs(A(i,j));
           if ( aij > max)
           {
              cout << "The current max: " << aij << endl;
              max = aij;  *p = i; *q = j;
           }
       }
   }
   return max;
}

// diagonalize using the jacobi algorithm
void diagonalize_matrix_jacobi(mat& A, mat& R, double tolerance, double iterations, double maxnondiag, int maxiter, int n )
{
    while ( maxnondiag > tolerance && iterations <= maxiter)
    {
       int p, q;
       maxnondiag   = offdiag(A, &p, &q, n);
       // cout << "start state: " << p << " : " << q << endl;
       // cout << "max non-diagonal element: " << maxnondiag << endl;
       Jacobi_rotate(A, R, p, q, n);
       iterations++;
       // cout << "current iteration: " << iterations << endl;
    }
}

//Jacobi's method for eigenvalues
void Jacobi_rotate (mat& A, mat& R, int k, int l, int n )
{
    double c,s;
    if ( A(k,l) != 0.0 )
    {
        double
        t, tau;
        tau = (A(l,l) - A(k,k))/(2*A(k,l));
        if ( tau >= 0 )
        {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else
        {
            t = -1.0/(-tau +sqrt(1.0 + tau*tau));
        }
        c = 1/sqrt(1+t*t);
        s = c*t;
    }
    else
    {
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk    = A(k,k);
    a_ll    = A(l,l);
    A(k,k)  = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l)  = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l)  = 0.0;
    // hard-coding non-diagonal elements by hand
    A(l,k)  = 0.0;
    // same here
    for (int i = 0; i < n; i++ )
    {
        if ( i != k && i != l )
        {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }
        //  And finally the new eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
    return;
    }
    // end of function jacobi_rotate


//int       i, j, Dim, lOrbital;
//double    RMin, RMax, Step, DiagConst, NondiagConst, OrbitalFactor;
//// With spherical coordinates RMin = 0 always
//RMin = 0.0;

//RMax = 8.0;  lOrbital = 0;  Dim =2000;
//mat Hamiltonian = zeros<mat>(Dim,Dim);
//// Integration step length
//Step    = RMax/ Dim;
//DiagConst = 2.0 / (Step*Step);
//NondiagConst =  -1.0 / (Step*Step);
//OrbitalFactor = lOrbital * (lOrbital + 1.0);

//// local memory for r and the potential w[r]
//vec r(Dim); vec w(Dim);
//for(i = 0; i < Dim; i++) {
//  r(i) = RMin + (i+1) * Step;
//  w(i) = potential(r(i)) + OrbitalFactor/(r(i) * r(i));
//}


//// Setting up tridiagonal matrix and brute diagonalization using Armadillo
//Hamiltonian(0,0) = DiagConst + w(0);
//Hamiltonian(0,1) = NondiagConst;
//for(i = 1; i < Dim-1; i++) {
//  Hamiltonian(i,i-1)    = NondiagConst;
//  Hamiltonian(i,i)    = DiagConst + w(i);
//  Hamiltonian(i,i+1)    = NondiagConst;
//}
//Hamiltonian(Dim-1,Dim-2) = NondiagConst;
//Hamiltonian(Dim-1,Dim-1) = DiagConst + w(Dim-1);
//// diagonalize and obtain eigenvalues
//vec Eigval(Dim);
//eig_sym(Eigval, Hamiltonian);
//output(RMin , RMax, Dim, Eigval);

//return 0;
