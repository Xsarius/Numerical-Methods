//
// Set of very basic, commonly used numerical methods.
// All of the methods are available in the web for free.
// Lib author: Kuba Jałoszyński
// Last update: Feb 2021
//
#include <math.h>

  //////////////////////////
 // Function declaration //
//////////////////////////
//  
//  Equation solver //
//
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
//      maxIter - maximal number of function iteration (corresponds to accuracy of the aproximation)
//
//  Returns:
//      Aproximation of the function root.   
//  
double Bisec(double (*function)(double), double leftEnd, double rightEnd, double tolerance, int maxIter);
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
//      maxIter - maximal number of function iteration (corresponds to accuracy of the aproximation)
//
//  Returns:
//      Aproximation of the function root.   
//  
double Falsi(double (*function)(double), double leftEnd, double rightEnd, double tolerance, int maxIter);
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
//      maxIter - maximal number of function iteration (corresponds to accuracy of the aproximation)
//
//  Returns:
//      Root of the function, nearest to the initial guess point.            
//  
double Newton(double (*function)(double), double guess, double tolerance, double epsilon, int maxIter);
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
//      maxIter - maximal number of function iteration (corresponds to accuracy of the aproximation)
//
//  Returns:
//      Aproximation of the function root.     
//  
double Secant(double (*function)(double), double leftEnd, double rightEnd, double tolerance, int maxIter);
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
//      First order method for solving ODEs with a given initial conditions.
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
//            
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



