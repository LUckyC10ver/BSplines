using System;
using System.Collections.Generic;

namespace BSpline.Core
{
    public sealed class BCCubicSpline : BCSpline
    {
        private readonly List<double> _d2y = new List<double>();

        public BCCubicSpline()
        {
        }

        public BCCubicSpline(IList<double> xValues, IList<double> yValues, bool sorted = false, double dyLeft = 1e30, double dyRight = 1e30)
        {
            Init(xValues, yValues, sorted, dyLeft, dyRight);
        }

        public string ClassName => "BCCubicSpline";

        public void Init(IList<double> xValues, IList<double> yValues, bool sorted = false, double dyLeft = 1e30, double dyRight = 1e30)
        {
            if (xValues == null || yValues == null)
            {
                BCException.Throw("BCCubicSpline.Init: values are null");
            }

            if (xValues.Count != yValues.Count)
            {
                BCException.Throw("BCCubicSpline.Init: vector sizes differ");
            }

            SetPoints(xValues, yValues, sorted);

            var n = GetXValues().Count;
            if (n < 2)
            {
                BCException.Throw("BCCubicSpline.Init: at least two points required");
            }

            _d2y.Clear();
            for (var i = 0; i < n; i++)
            {
                _d2y.Add(0.0);
            }

            var u = new double[n - 1];
            if (dyLeft > 0.99e30)
            {
                _d2y[0] = 0.0;
                u[0] = 0.0;
            }
            else
            {
                _d2y[0] = -0.5;
                u[0] = 3.0 / (GetXValues()[1] - GetXValues()[0]) *
                       ((GetYValues()[1] - GetYValues()[0]) / (GetXValues()[1] - GetXValues()[0]) - dyLeft);
            }

            for (var i = 1; i < n - 1; i++)
            {
                var sig = (GetXValues()[i] - GetXValues()[i - 1]) / (GetXValues()[i + 1] - GetXValues()[i - 1]);
                var p = sig * _d2y[i - 1] + 2.0;
                _d2y[i] = (sig - 1.0) / p;
                var uVal = (GetYValues()[i + 1] - GetYValues()[i]) / (GetXValues()[i + 1] - GetXValues()[i]) -
                           (GetYValues()[i] - GetYValues()[i - 1]) / (GetXValues()[i] - GetXValues()[i - 1]);
                uVal = (6.0 * uVal / (GetXValues()[i + 1] - GetXValues()[i - 1]) - sig * u[i - 1]) / p;
                u[i] = uVal;
            }

            double qn;
            double un;
            if (dyRight > 0.99e30)
            {
                qn = 0.0;
                un = 0.0;
            }
            else
            {
                qn = 0.5;
                un = 3.0 / (GetXValues()[n - 1] - GetXValues()[n - 2]) *
                     (dyRight - (GetYValues()[n - 1] - GetYValues()[n - 2]) / (GetXValues()[n - 1] - GetXValues()[n - 2]));
            }

            _d2y[n - 1] = (un - qn * u[n - 2]) / (qn * _d2y[n - 2] + 1.0);
            for (var k = n - 2; k >= 0; k--)
            {
                _d2y[k] = _d2y[k] * _d2y[k + 1] + u[k];
            }
        }

        public double Evaluate(double x)
        {
            var xValues = GetXValues();
            var yValues = GetYValues();
            var n = xValues.Count;
            if (n < 2)
            {
                BCException.Throw("BCCubicSpline.Evaluate: values not initialized");
            }

            var klo = 0;
            var khi = n - 1;
            while (khi - klo > 1)
            {
                var k = (khi + klo) / 2;
                if (xValues[k] > x)
                {
                    khi = k;
                }
                else
                {
                    klo = k;
                }
            }

            var h = xValues[khi] - xValues[klo];
            if (h > GetEpsilon())
            {
                var a = (xValues[khi] - x) / h;
                var b = (x - xValues[klo]) / h;
                var res = a * yValues[klo] + b * yValues[khi];
                res += ((a * a * a - a) * _d2y[klo] + (b * b * b - b) * _d2y[khi]) * h * h / 6.0;
                return res;
            }

            return (yValues[klo] + yValues[khi]) / 2.0;
        }

        public double GetDerivative(double x)
        {
            var xValues = GetXValues();
            var yValues = GetYValues();
            var n = xValues.Count;
            if (n < 2)
            {
                BCException.Throw("BCCubicSpline.GetDerivative: values not initialized");
            }

            var klo = 0;
            var khi = n - 1;
            while (khi - klo > 1)
            {
                var k = (khi + klo) / 2;
                if (xValues[k] > x)
                {
                    khi = k;
                }
                else
                {
                    klo = k;
                }
            }

            var h = xValues[khi] - xValues[klo];
            if (h > GetEpsilon())
            {
                var a = (xValues[khi] - x) / h;
                var b = (x - xValues[klo]) / h;
                var res = (yValues[khi] - yValues[klo]) / h;
                res -= (3 * a * a - 1) / 6.0 * h * _d2y[klo];
                res += (3 * b * b - 1) / 6.0 * h * _d2y[khi];
                return res;
            }

            return (yValues[klo] + yValues[khi]) / 2.0;
        }

        public double GetSecondDerivative(double x)
        {
            var xValues = GetXValues();
            var n = xValues.Count;
            if (n < 2)
            {
                BCException.Throw("BCCubicSpline.GetSecondDerivative: values not initialized");
            }

            var klo = 0;
            var khi = n - 1;
            while (khi - klo > 1)
            {
                var k = (khi + klo) / 2;
                if (xValues[k] > x)
                {
                    khi = k;
                }
                else
                {
                    klo = k;
                }
            }

            var h = xValues[khi] - xValues[klo];
            if (h > GetEpsilon())
            {
                var a = (xValues[khi] - x) / h;
                var b = (x - xValues[klo]) / h;
                return a * _d2y[klo] + b * _d2y[khi];
            }

            return (_d2y[klo] + _d2y[khi]) / 2.0;
        }

        public double GetThirdDerivative(double x)
        {
            var xValues = GetXValues();
            var n = xValues.Count;
            if (n < 2)
            {
                BCException.Throw("BCCubicSpline.GetThirdDerivative: values not initialized");
            }

            var klo = 0;
            var khi = n - 1;
            while (khi - klo > 1)
            {
                var k = (khi + klo) / 2;
                if (xValues[k] > x)
                {
                    khi = k;
                }
                else
                {
                    klo = k;
                }
            }

            var h = xValues[khi] - xValues[klo];
            if (h > GetEpsilon())
            {
                return (_d2y[khi] - _d2y[klo]) / h;
            }

            return 0.0;
        }
    }
}
