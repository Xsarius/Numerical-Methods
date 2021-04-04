using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GdGC
{
    public class GD_Relation
    {
        public readonly string name;
        public readonly string xLabel;
        public readonly string yLabel;
        private readonly Point[] initialPoints;
        private readonly double resolution = 0.1;


        public GD_Relation(string name, string xLabel, string yLabel, Point[] point)
        {
            this.name = name;
            this.xLabel = xLabel;
            this.yLabel = yLabel;
            this.initialPoints = point;
        }

        public double YValue(double x0)
        {
            double y0 = Calculator.RBFInterpolation(initialPoints, x0, resolution);

            return y0;
        }

        public double XValue(double y0)
        {
            Point[] temp = Array.Empty<Point>();

            for (int i = 0; i < initialPoints.Length; i++)
            {
                temp[i].x = initialPoints[i].y;
                temp[i].y = initialPoints[i].x;
            }

            double x0 = 0;

            

            return x0;
        }


    }
}
