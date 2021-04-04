using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GdGC
{
    public struct Point
    {
        public string name;
        public double x;
        public double y;
    }

    public class Calculator
    {
        public static double RBFInterpolation(Point[] point, double x0, double shapeFactor)
        {

            int matrixSize = point.Length;
    
           
            double y0 = 0;
            double[][] A = new double[matrixSize][];
            double[] b = new double[matrixSize];

            for (int i = 0; i < matrixSize; i++)
            {
                for (int j = 0; j < matrixSize; j++)
                {
                    A[j][i] = GaussRBF(Math.Abs(point[j].x - point[i].x), shapeFactor);
                }

                b[i] = point[i].y;
            }

            double[] lambda = SolveLinearSystem(A, b, matrixSize);

            for (int i = 0; i < matrixSize; i++)
            {
                y0 += lambda[i]*GaussRBF(Math.Abs(x0 - point[i].x), shapeFactor);
            }

            return y0;
        }

        public static double[] SolveLinearSystem(double[][] initMatrix, double[] valueMatrix, int matrix1Dsize)
        {
            int matrixSize = matrix1Dsize;

            double[][] L = new double[matrixSize][];
            double[][] U = new double[matrixSize][];
            double[] Y = new double[matrixSize];
            double[] sol = new double[matrixSize];

            for (int i = 0; i < matrixSize; i++)
            {
                for(int k = i; k < matrixSize; k++)
                {
                    double sum = 0;

                    for(int j = 0; j < i; j++)
                    {
                        sum += L[i][j] * U [j][k];
                    }

                    U[i][k] = initMatrix[i][k] - sum;
                }

                for (int k = 0; k < matrixSize; k++)
                {
                    if( i == k )
                    {
                        L[i][i] = 1;
                    }

                    else
                    {
                        double sum = 0;

                        for (int j = 0; j < i; j++)
                        {
                            sum += L[k][j] * U[j][i];
                        }

                        L[k][i] = (initMatrix[k][i] - sum)/ U[i][i];
                    }
                }
            }

            for (int i = 0; i < matrixSize; i++)
            {
                Y[i] = valueMatrix[i];
                for(int j = 0; j < i; j++)
                {
                    Y[i] -= L[i][j] * Y[j];
                }
            }

            for (int i = matrixSize - 1; i >= 0; i--)
            {
                sol[i] = Y[i];
                for (int j = i + 1; j < matrixSize; j++)
                {
                    sol[i] -= U[i][j] * sol[j];
                }
                sol[i] /= U[i][i];
            }

            return sol;
        }

        private static double GaussRBF(double r, double eps)
        {
            return Math.Exp(-eps*Math.Pow(r,2));
        }

    }

    
}
