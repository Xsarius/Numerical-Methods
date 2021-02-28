//
// Set of very basic, commonly used numerical methods.
// All of the methods are available in the web for free.
// Lib author: Kuba Jałoszyński
// Last update: Feb 2021
//
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define MAX_MATRIX_SIZE 100 // Arbitrarly chosen value - it's only necessary to code be compileable as C.

//////////////////////////
//  Struct declaration //
////////////////////////
//
// Data
//  Summary:
//     
//
struct Point
{
    double x, y;
};

  //////////////////////////
 // Function declaration //
//////////////////////////
//  
//  Equation solver //
//
// Muller method
//
//  Summary:
//      The actual root is aproximated by creating quadratic formula, based on 3 guesses (root 1,2,3).
//        If aproximation is to far away from real root, guess root's are adjusted and proccess iterates 
//        untill value in the tolerance range is obtained.
//  
//      Pros:
//       - Solves imaginary roots
//       - Doesn't require derivatives
//
//      Cons:
//       - Extraneous roots can be obtained
//       - Possible big error
//
//  Variables:
//      function - primary function
//
//      rootOne - first root guess 
//
//      rootTwo - second root guess
//
//      rootThree - third root guess
//  Returns:
//      Aproximation of the function root.   
//
double Muller(double (*function)(double), double rootOne, double rootTwo, double rootThree, double tolerance, int maxIterations);
// Bisec method
//  Summary:
//      Bisec aproximation of function root value in given a to b range.
//        Iterative method called interval halving. It's simple and robust, it'll always find the root aproximation,
//        but it's simultaneously slow.
//  Variables:
//      function - primary function
//
//      leftEnd - lower integral limit
//
//      rightEnd - higher integral limit
//
//      tolerance - root aproximation error
//
//      maxIterations - maximal number of function iteration (corresponds to accuracy of the aproximation)
//
//  Returns:
//      Aproximation of the function root.   
//  
double Bisec(double (*function)(double), double leftEnd, double rightEnd, double tolerance, int maxIterations);
//
// Falsi method
//  Summary:
//      Regula falsi aproximation of function root value in given a to b range.
//        It uses trial and error technique, equivalent to linear interpolation. 
//  Variables:
//      function - primary function
//
//      leftEnd - lower integral limit
//
//      rightEnd - higher integral limit
//
//      tolerance - root aproximation error
//
//      maxIterations - maximal number of function iteration (corresponds to accuracy of the aproximation)
//
//  Returns:
//      Aproximation of the function root.   
//  
double Falsi(double (*function)(double), double leftEnd, double rightEnd, double tolerance, int maxIterations);
//
// Newton method
//  Summary:
//      Newton aproximation of function root value nearest to the point of guess.
//        It uses tangents to aproximate true value of a root.
// 
//  Variables:
//      function - primary function
//
//      guess - potential root point (should be estimated as close as possible)
//
//      tolerance - root aproximation error
//
//      epsilon - 
//
//      maxIterations - maximal number of function iteration (corresponds to accuracy of the aproximation)
//
//  Returns:
//      Root of the function, nearest to the initial guess point.            
//  
double Newton(double (*function)(double), double guess, double tolerance, double epsilon, int maxIterations);
//
// Secant method
//  Summary:
//      Secant aproximation of function root value in given a to b range.
//        It uses recurrence method. It's simillar to newton method.
//
//  Variables:
//      function - primary function
//
//      leftEnd - lower integral limit
//
//      rightEnd - higher integral limit
//
//      tolerance - root aproximation error
//
//      maxIterations - maximal number of function iteration (corresponds to accuracy of the aproximation)
//
//  Returns:
//      Aproximation of the function root.     
//  
double Secant(double (*function)(double), double leftEnd, double rightEnd, double tolerance, int maxIterations);
//
//  Integration //
//
// Simpson method
//  Summary:
//      Simpsonal aproximation of a function integral value in given range a to b.
//        Integral is aproximated by quadratic function passing through calculated points.
//        
//  Variables:
//      function - primary function
//
//      leftEnd - lower integral limit
//
//      rightEnd - higher integral limit
//
//      tolerance - root aproximation error
//
//      substep - number of substeps of the integral aproximation
//
//  Returns:
//      Value of the function integral in given range a to b.         
//  
double Simpson(double (*function)(double), double leftEnd, double rightEnd, int substeps);
//
// Trapezoidal method
//  Summary:
//      Trapezosidal aproximation of a function integral value in given range a to b.
//        Namely, integral is split into sum of a trapezies, which folows the function curve in a given range.
//
//  Variables:
//      function - primary function
//
//      leftEnd - lower integral limit
//
//      rightEnd - higher integral limit
//
//      substep - number of substeps of the integral aproximation
//
//  Returns:
//      Value of the function integral in given range a to b.       
//  
double Trapeze(double (*function)(double), double leftEnd, double rightEnd, int substep);
//
//  Differentials //
//
// Derivative 
//  Summary:
//      Quick aproximation of the value of a function derivative in given point,
//          using derivative deffinition as a limit of a function with infinitesimally small step dx.
//
//  Variables:
//      function - primary function
//
//      x - point at which derivative value is wanted
//
//  Returns:
//      Value of the function derivative at given point x.    
// 
double Deriv(double (*function)(double), double);
//
// Euler method
//
//  Summary:
//      Method for solving first order ODEs with given initial conditions.
//        It have stict region of stability. If {h,k c C} |h*k| > 1 where h - step size, k - solution exponent coeff, 
//        then numerical solution will diverge. 
//  Variables:
//      function - primary function
//
//      x0 - initial argument
//
//      y0 - initial condition
//
//      maximalX - final argument
//
//      stepSize - resolution of function  
//
//  Returns:
//      Value of the ODE at given x (maximalX) point.
//
double Euler(double (*function)(double, double), double x0, double y0, double maximalX, double stepSize);
//
// Runge-Kutta 4th method
//
//  Summary:
//      Method for solving first order ODEs with given initial conditions.
//        It's much more accurate than euler method, but it's more resource and time consumming.      
//
//  Variables:
//      function - primary function
//
//      x0 - initial argument
//
//      y0 - initial condition
//
//      maximalX - final argument
//
//      stepSize - resolution of function 
//
//  Returns:
//      Value of the ODE at given x (maximalX) point.    
//
double RungeKutta(double (*function)(double, double), double x0, double y0, double maximalX, double stepSize);
//
//  Interpolation //
//
// Lagrange's Interpolation
//
//  Summary:
//      Method of finding new data points within the range of a known data points.
//
//  Variables:
//      point[] - array of the Point type, containing combined known arguments with their values.
//      
//      ptsNum - number of all known points.
//
//      x - point which is to be obtained in the interpolation.
//
//  Returns:
//      Value of a point at given argument x.
//
double Lagrange(Point point[], int x, int ptsNum);
//
// Linear Interpolation
//
//  Summary:
//      Method of finding new data points within the range of a known data points.
//
//  Variables:
//      point[] - array of the Point type, containing combined known arguments with their values 
//                  of two nearest points to wanted point.
//
//      x - wanted point.
//
//  Returns:
//      Value aproximation at point x.
//
double Linear(Point point[2], double x);
//
// GaussRBF Interpolation
//
//  Summary:
//      Method of finding new data points within the range of a known data points.
//  Variables:
//      points - array of known points combined, arguments with their values.
//
//      x - wanted point.
//
//      shapeParameter - arbitrarly chosen natural number, values << 1 are recommended if distance between points is big > 0.1
//                        and >> 1 if distance is small (1e-3 order).
//
//  Returns:
//      Value aproximation at point x.
//
double RBFIterpolation(Point points[], double x, double shapeParameter);
//
//  Functions //
//
// Gauss radial basis function (GaussRBF)
//
//  Summary:
//      Function whose value depends on the distance between the input and some fixed point.
//      
//  Variables:
//      r - distance between chosen point and reference point.
//
//      eps - shape factor, can be arbitrarly chosen natural number. 
//  Returns:
//
//
double GaussRBF(double r, double eps);
//
//  Matricies //
//
// Linear System Solver
//
//  Summary:
//      Function solves system of the linear equations AX = B.
//       It does not provide not quadratic matrix protection,
//       which means, the user should make sure to input quadratic matrix A
//       and X, B to be coresponding size.
//  
//  Variables:
//      matricies1DSize - size of arrays.
//
//      coefficientsMatrix - array with unknowns coefficients in each of systems equation.
//
//      valueMatrix - array with known values of the system. 
//
//      unknownsMatrix - array in which the results will be saved.
//
void SolveLinearSystem(int matricies1DSize,double coefficientsMatrix[][MAX_MATRIX_SIZE], double valueMatrix[], double unknownsMatrix[]);




