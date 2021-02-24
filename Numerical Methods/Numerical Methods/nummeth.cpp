#include "nummeth.h"

//
// Bisec Method
//
double Bisec(double (*fun)(double), double leftEnd, double rightEnd, double tolerance, int maxIter)
{

    for (int i = 0; i < maxIter; i++)
    {
    double midPoint = (leftEnd + rightEnd)/2.;

    if( fun(midPoint) == 0 || (rightEnd - leftEnd)/2. < tolerance)
    {
        return midPoint;
    }

    else if(fun(leftEnd) * fun(midPoint) > 0)
    {   
        leftEnd = midPoint;
    }
    else if(fun(rightEnd) * fun(midPoint) > 0)
    {
        rightEnd = midPoint;
    }

    }
    
}
//
// Falsi Method
//
double Falsi(double (*fun)(double), double leftEnd, double rightEnd, double tolerance, int maxIter)
{
    double root = 0, fRoot;
    
    double fLeft = fun(leftEnd);
    double fRight = fun(rightEnd);

    for(int n = 0, side = 0 ; n < maxIter; n++)
    {
        root = (fLeft*rightEnd - fRight*leftEnd) / (fLeft - fRight);

        if(fabs(rightEnd - leftEnd) < tolerance*fabs(rightEnd+leftEnd))
        {
            break;
        }

        fRoot = fun(root);

        if(fRoot*fRight > 0)
        {
            rightEnd = root;
            fRight = fRoot;

            if(side == -1)
            {
                fLeft /= 2.;
            }

            side = -1;
        }

        else if(fRoot*fLeft > 0)
        {
            leftEnd = root;
            fLeft = fRoot;

            if(side == 1)
            {
                fRight /= 2.;
            }

            side = 1;
        }

        else
        { 
            break;
        }

    }

    return root;
}
//
// Newton method
//
double Newton(double (*fun)(double), double guess, double tolerance, double epsilon, int maxIter)
{
    double y, yPrime, x1;
    
    for (int i = 0; i < maxIter; i++)
    {
        y = fun(guess);
        yPrime = Deriv(fun, guess);

        if(fabs(yPrime) < epsilon)
        {
            break;
        }

        x1 = guess - y/yPrime;

        if(fabs(x1-guess) <= tolerance)
        {
            break;
        }

        guess = x1;
    }
    

    return guess;
}
//
// Secant method
//
double Secant(double (*fun)(double), double leftEnd, double rightEnd, double tolerance, int maxIter)
{
    double root = 0, fRoot;

    for (int i = 0; i < maxIter; i++)
    {
        root = rightEnd - fun(rightEnd) * (rightEnd - leftEnd) / (fun(rightEnd) - fun(leftEnd));
        fRoot = fun(root);

        if (fabs(fRoot) <= tolerance)
        {
            return root;
        }

        leftEnd = rightEnd;
        rightEnd = root;
    }

    return root; 
}
//
// Derivative 
//
double Deriv(double (*fun)(double), double x)
{
    const double delta = 1.0e-9; 
    double x1 = x - delta;
    double x2 = x + delta;
    double y1 = fun(x1);
    double y2 = fun(x2);
    
    return x* (y2 - y1) / (x2 - x1);
}


