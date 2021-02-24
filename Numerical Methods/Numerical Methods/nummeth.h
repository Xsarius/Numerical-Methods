//
// Commonly used numerical methods.
// All of the methods source code is available in the web for free.
// File author: Kuba Jałoszyński
// Last update: Feb 2021
//

#include <math.h>

  //////////////////////////
 // Function declaration //
//////////////////////////
// 
//  Equation solver //
//
//  Bisec method
//
double Bisec(double (*fun)(double), double leftEnd, double rightEnd, double tolerance, int maxIter);
//
//  Falsi method
//
double Falsi(double (*fun)(double), double leftEnd, double rightEnd, double tolerance, int maxIter);
//
//  Newton method
//
double Newton(double (*fun)(double), double guess, double tolerance, double epsilon, int maxIter);
//
//  Secant method
//
double Secant(double (*fun)(double), double leftEnd, double rightEnd, double tolerance, int maxIter);
//
//  Integration //
//
//  Simpson method
//
double Simpson();
//
//  Trapezoidal method
//
double Trapeze();
//
//  Numerical differentiation //
//
//  Derivative 
//
double Deriv(double (*fun)(double), double);
