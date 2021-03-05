#include "nummeth.h"

//
// Muller method
//
double Muller(double (*f)(double), double p1, double p2, double p3, double tol, int n)
{
    double res = 0, a, b, c, delta, x1, x2;

    for (int i = 0; i < n; ++i)
    {
        a = (f(p3) - f(p1) - (p3-p1)*((f(p2)-f(p1))/(p2-p1)))/(pow(p3,2.)-pow(p1,2.)-(p3-p1)*((pow(p2,2.) - pow(p3,2.))/(p2-p1)));
        b = (f(p2)-f(p1)-a*(pow(p2,2.)-pow(p1,2.)))/(p2-p1);
        c = f(p1) - a*pow(p1,2.) - b * p1;
        
        delta = pow(b,2.)-4*a*c;

        x1 = (-b - sqrt(delta))/(2*a);
        x2 = (-b + sqrt(delta))/(2*a);

        if(fabs(f(x1)) < fabs(f(x2)))
        {
            res = x1;
        }
        else
        {
            res = x2;
        }

        if(fabs(f(res)) < tol)
        {
            return res;
        }

        p1 = res;
    }
}
//
// Bisec root method
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
// Falsi root method
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
    double h = (b-a)/double(n-1);

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
//
// PolynomialIter interpolation
//
double PolynomialIter(Point p[], int x, int n)
{
    double r = 0, term = 0;

    for (int i = 0; i < n; i++)
    {
        term = p[i].y;
        for (int j = 0; j < n; j++)
        {
            if(j != i)
            {
            term *= (x - p[j].x)/(double)(p[i].x - p[j].x);
            }   
        }
        
        r += term;
    }
    
    return r;
}
//
// Bilinear interpolation
//
double LinearIter(Point p[2], double x)
{
    double m = (p[1].y-p[0].y)/(p[1].x-p[0].x);
    double b = p[0].y - m*p[0].x;

    return b + x*m;
}
//
// Radial Basis Interpolation
//
double RBFIter(Point point[], double x0, double eps)
{
    int n = 3;

    double A[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];

    double *lambda = (double *)malloc(n * sizeof(double));
    double *b = (double *)malloc(n * sizeof(double));
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[j][i] = GaussRBF(fabs(point[j].x - point[i].x), eps);
        }
     
        b[i] = point[i].y;
    }

    SolveLinearSystem(n, A, b, lambda);

    double y0 = 0;

    for (int i = 0; i < n; i++)
    {
            y0 += lambda[i]*GaussRBF(fabs(x0 - point[i].x), eps);
    }
    

    free(lambda);
    free(b);

    return y0;
}
//
// Newton Forward Difference Interpolation 
// 
double NewtonFWDInter(Point p[], double x0)
{
    int n = 4;

    double **y = (double **)malloc(n * sizeof(double*));

    for (int i = 0; i < n; i++)
    {
        y[i] = (double *)malloc(n * sizeof(double));
        y[i][0] = p[i].y;
    }

    for(int i = 1; i < n; i++)
    {
        for(int j = 0; j < n-i; j++)
        {
            y[j][i] = y[j+1][i-1] - y[j][i-1];
        }
    }

    double sum = y[0][0];
    double u = (x0 - p[0].x) / (p[1].x - p[0].x);
    double temp;

    for (int i = 1; i < n; i++)
    {
        temp = u;
        for (int k = 1; k < i; k++)
        {
            temp *= (u - k);
        }
        
        sum += (temp * y[0][i]) / Factorial(i);
    }
    
    
    for (int i = 0; i < n; i++)
    {
        free(y[i]);
    }

    free(y);
    
    return sum;
}
//
// Newton Forward Difference Interpolation 
//
double NewtonRWDInter(Point p[], double x0)
{
    int n = 5;

    printf("%d\n", n);

    double **y = (double **)malloc(n * sizeof(double*));

    for (int i = 0; i < n; i++)
    {
        y[i] = (double *)malloc(n * sizeof(double));
        y[i][0] = p[i].y;
    }
    
    for (int i = 1; i < n; i++) 
    { 
        for (int j = n - 1; j >= i; j--)
        {
            y[j][i] = y[j][i - 1] - y[j - 1][i - 1]; 
        } 
    } 

    double sum = y[n-1][0];
    double temp;
    double u = (x0 - p[n-1].x)/(p[1].x - p[0].x);

    for (int i = 1; i < n; i++)
    {
        temp = u;
        for (int k = 1; k < i; k++)
        {
            temp *= (u + k);
        }

        sum += (temp*y[n-1][i]) / Factorial(i);
    }
    
    
    for (int i = 0; i < n; i++)
    {
        free(y[i]);
    }

    free(y);
    
    printf("NewtonRWD: %.4f\n", sum);

    return sum;
}
//
// Gaussian Radial Basis Function 
//
double GaussRBF(double r, double eps)
{
    return exp(-eps*pow(r, 2.));
}
//
// Factiorial function
//
int Factorial(int x)
{
    double temp = 1;
    
    for (int i = 2; i <= x; i++)
    {
        temp *= i;
    }

    return temp;
}
//
// Matrix equation solver
//
void SolveLinearSystem(int n ,double A[][MAX_MATRIX_SIZE], double B[], double X[])
{

    double **L = (double **)malloc(sizeof(double *) * n);
    double **U = (double **)malloc(sizeof(double *) * n);

    for (int i = 0; i < n; i++)
    {
        L[i] = (double *)malloc(n * sizeof(double));
        U[i] = (double *)malloc(n * sizeof(double));
    }

    double *Y = (double *)malloc(sizeof(double) * n);
    
    // Calculating lower L and upper U matricies
    for (int i = 0; i < n; i++)
    {
        for (int k = i; k < n; k++)
        {
            double sum = 0;

            for (int j = 0; j < i; j++)
            {
                sum += L[i][j] * U[j][k];
            }

            U[i][k] = A[i][k] - sum;
        }

        for (int k = 0; k < n; k++)
        {
            if( i == k)
            {
                L[i][i] = 1;
            }

            else
            {
                double sum =0;

                for (int j = 0; j < i; j++)
                {
                        sum += L[k][j] * U[j][i];
                }

                L[k][i] = (A[k][i] - sum) / U[i][i];
                
            }
        } 
    }

    // Calculating temporary Y matrix
    for(int i=0; i<n; i++)
    {
        Y[i]=B[i];
        for(int j=0; j<i; j++)
        {
            Y[i]-=L[i][j]*Y[j];
        }
    }
    
    // Calculating final value X matrix
    for(int i=n-1; i>=0; i--)
    {
        X[i]= Y[i];
        for(int j=i+1; j<n; j++)
        {
            X[i]-=U[i][j]*X[j];
        }
        X[i]/=U[i][i];
    }


    // Debugging and optimalisation purposes 
    //
    for (int i = 0; i < n; i++)
    {
        free(U[i]);
        free(L[i]);
    }
    
    free(L);
    free(U);
    free(Y);
    //
}
