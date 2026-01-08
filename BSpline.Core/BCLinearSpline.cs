using System;
using System.Collections.Generic;

namespace BSpline.Core
{
    public sealed class BCLinearSpline : BCSpline
    {
        public BCLinearSpline()
        {
        }

        public BCLinearSpline(IList<double> points, IList<double> values, bool sorted = false)
            : base(points, values, sorted)
        {
        }

        public string ClassName => "BCLinearSpline";

        public double Evaluate(double argument)
        {
            var xValues = GetXValues();
            var yValues = GetYValues();
            if (xValues.Count == 0 || yValues.Count == 0)
            {
                BCException.Throw("BCLinearSpline.Evaluate: values not initialized");
            }

            var n = IndexFromValue(argument);
            if (n + 1 >= xValues.Count)
            {
                n = xValues.Count - 2;
            }

            if (n + 1 < xValues.Count)
            {
                var h = xValues[n + 1] - xValues[n];
                if (h > GetEpsilon())
                {
                    return yValues[n] + (yValues[n + 1] - yValues[n]) * (argument - xValues[n]) / h;
                }

                return (yValues[n] + yValues[n + 1]) / 2.0;
            }

            BCException.Throw("BCLinearSpline.Evaluate: argument out of range");
            return 0.0;
        }
    }
}
