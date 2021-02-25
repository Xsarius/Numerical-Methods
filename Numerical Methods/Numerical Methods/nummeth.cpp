#include "nummeth.h"

//
// Bisec root Method
//  
double Bisec(double (*f)(double), double a, double b, double tol, int n)
{

    for (int i = 0; i < n; i++)
    {
    double c = (a + b)/2.;

    if( f(c) == 0 || (b - a)/2. < tol)
    {
        return c;
    }

    else if(f(a) * f(c) > 0)
    {   
        a = c;
    }
    else if(f(b) * f(c) > 0)
    {
        b = c;
    }

    }
    
}
//
// Falsi root Method
//   
double Falsi(double (*f)(double), double a, double b, double tol, int n)
{
    double r = 0, fr;
    
    double fa = f(a);
    double fb = f(b);

    for(int n = 0, side = 0 ; n < n; n++)
    {
        r = (fa*b - fb*a) / (fa - fb);

        if(fabs(b - a) < tol*fabs(b+a))
        {
            break;
        }

        fr = f(r);

        if(fr*fb > 0)
        {
            b = r;
            fb = fr;

            if(side == -1)
            {
                fa /= 2.;
            }

            side = -1;
        }

        else if(fr*fa > 0)
        {
            a = r;
            fa = fr;

            if(side == 1)
            {
                fb /= 2.;
            }

            side = 1;
        }

        else
        { 
            break;
        }

    }

    return r;
}
//
// Newton root method
// 
double Newton(double (*f)(double), double x0, double tol, double eps, int n)
{
    double y, yp, x1;
    
    for (int i = 0; i < n; i++)
    {
        y = f(x0);
        yp = Deriv(f, x0);

        if(fabs(yp) < eps)
        {
            break;
        }

        x1 = x0 - y/yp;

        if(fabs(f(x1)) <= tol)
        {
            x0 = x1;
            break;
        }

        x0 = x1;
    }

    return x0;
}
//
// Secant root method
//    
double Secant(double (*f)(double), double a, double b, double tol, int n)
{
    double r = 0, fr;

    for (int i = 0; i < n; i++)
    {
        r = b - f(b) * (b - a) / (f(b) - f(a));
        fr = f(r);

        

        if (fabs(fr) <= tol)
        {
            return r;
        }

        a = b;
        b = r;
    }

}
//
// Derivative 
//  
double Deriv(double (*f)(double), double x)
{
    const double dx = 1e-6; 
    
    return (f(x + dx) - f(x))/dx;
}
//
// Simpson integration method
//
double Simpson(double (*f)(double), double a, double b, int n)
{
    double area = f(a) + f(b);
    double h = (b-a)/(n-1);

    for (int i = 0; i < n; i++)
    {
        if(i%2 == 0)
        {
            area += 2*f(a + i*h);
        }
        else
        {
            area += 4*f(a+i*h);
        }
    }
    
    area *= h/3.;

    return area;
}
//
// Trapeze integration method
//
double Trapeze(double (*f)(double), double a, double b, int n)
{
    double h = (b-a)/(n-1), area = 0, x0 = a, x1 = x0+h;

    for (int i = 0; i < n; i++)
    {
        area += 0.5*(f(x1) + f(x0))*(x1-x0);

        x0 = x1;
        x1 += h;
    }
    
    return area;
}
//
// Euler method
//
double Euler(double (*f)(double, double), double x0, double y0, double x, double h)
{
    
  while (x0 < x)
  {
      y0 += h*f(x0,y0);
      x0 += h;
  }
  
  return y0;
}
//
// Runge-Kutta method
//
double RungeKutta(double (*f)(double, double), double x0, double y0, double x, double h)
{
    int n = (int)((x - x0)/h);

    double k1,k2,k3,k4;

    double y = y0;

    for (int i = 0; i < n; i++)
    {
        k1 = h*f(x0, y);
        k2 = h*f(x0 + .5*h, y + .5*k1);
        k3 = h*f(x0 + .5*h, y + .5*k2);
        k4 = h*f(x0 + h, y + k3);

        y += (1./6.)*(k1 + 2*k2 + 2*k3 + k4);

        x0 = x0 +h;
    }
    
    return y;
}
