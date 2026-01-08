using System;
using System.Collections.Generic;

namespace BSpline.Core
{
    public sealed class BCGaussLobattoSpline
    {
        private const double Pi = Math.PI;

        private readonly List<double> _data = new List<double>();
        private double _xmin;
        private double _xmax;
        private int _mode;

        public BCGaussLobattoSpline(int mode = -1)
        {
            _xmin = 0.0;
            _xmax = 0.0;
            _mode = mode > -1 ? mode : 0;
        }

        public BCGaussLobattoSpline(IList<double> xValues, IList<double> yValues, bool sorted = false, int mode = -1)
        {
            _mode = mode > -1 ? mode : 0;
            Init(xValues, yValues, sorted, mode);
        }

        public string ClassName => "BCGaussLobattoSpline";

        public bool Empty()
        {
            return _data.Count == 0;
        }

        public IReadOnlyList<double> GetDataVector()
        {
            return _data.AsReadOnly();
        }

        public List<double> GetDataVectorMutable()
        {
            return _data;
        }

        public Tuple<double, double> GetInterval()
        {
            return Tuple.Create(_xmin, _xmax);
        }

        public int Size()
        {
            return _data.Count;
        }

        public double XMin()
        {
            return _xmin;
        }

        public double XMax()
        {
            return _xmax;
        }

        public int Mode()
        {
            return _mode;
        }

        public double Evaluate(double x)
        {
            if (_data.Count < 2)
            {
                BCException.Throw("BCGaussLobattoSpline.Evaluate: value size less than 2");
            }

            var index = (int)IndexOfX(x);
            if (index <= -1)
            {
                return _data[0];
            }

            if (index > _data.Count - 1)
            {
                return _data[_data.Count - 1];
            }

            var xLeft = XOfIndex(index);
            var xRight = XOfIndex(index + 1);
            var slope = (_data[index + 1] - _data[index]) / (xRight - xLeft);
            return _data[index] + slope * (x - xLeft);
        }

        public double GetIntegral(double a, double b)
        {
            if (a < _xmin || a > _xmax)
            {
                BCException.Throw($"BCGaussLobattoSpline.GetIntegral: input value a={a} out of range [{_xmin},{_xmax}]");
            }

            if (b < _xmin || b > _xmax)
            {
                BCException.Throw($"BCGaussLobattoSpline.GetIntegral: input value b={b} out of range [{_xmin},{_xmax}]");
            }

            if (_data.Count < 2)
            {
                BCException.Throw("BCGaussLobattoSpline.GetIntegral: value size less than 2");
            }

            var sign = 1.0;
            if (b < a)
            {
                var temp = a;
                a = b;
                b = temp;
                sign = -1.0;
            }

            var indexOfA = (int)IndexOfX(a);
            var indexOfB = (int)IndexOfX(b);
            EnsureInRange(indexOfA, 0, _data.Count - 1, "left boundary");
            EnsureInRange(indexOfB, 0, _data.Count - 1, "right boundary");

            double result;
            if (indexOfA == indexOfB)
            {
                var xLeft = XOfIndex(indexOfA);
                var xRight = XOfIndex(indexOfA + 1);
                var slope = (_data[indexOfA + 1] - _data[indexOfA]) / (xRight - xLeft);
                var ya = _data[indexOfA] + slope * (a - xLeft);
                var yb = _data[indexOfA] + slope * (b - xLeft);
                result = 0.5 * (ya + yb) * (b - a);
            }
            else
            {
                result = 0.0;
                result += 0.5 * (_data[indexOfA + 1] + Evaluate(a)) * (XOfIndex(indexOfA + 1) - a);
                result += 0.5 * (Evaluate(b) + _data[indexOfB]) * (b - XOfIndex(indexOfB));

                if (indexOfB - indexOfA > 1)
                {
                    for (var i = indexOfA + 1; i < indexOfB; i++)
                    {
                        result += 0.5 * (_data[i] + _data[i + 1]) * (XOfIndex(i + 1) - XOfIndex(i));
                    }
                }
            }

            return result * sign;
        }

        public void Init(int newSize, double initializer, double xmin, double xmax, int mode = -1)
        {
            if (xmin > xmax)
            {
                BCException.Throw($"BCGaussLobattoSpline.Init: invalid range [xmin,xmax]=[{xmin},{xmax}]");
            }

            if (newSize < 1)
            {
                BCException.Throw($"BCGaussLobattoSpline.Init: invalid size={newSize}");
            }

            if (mode > -1)
            {
                _mode = mode;
            }

            _xmin = xmin;
            _xmax = xmax;
            _data.Clear();
            for (var i = 0; i < newSize; i++)
            {
                _data.Add(initializer);
            }
        }

        public void Init(IList<double> xValues, IList<double> yValues, bool sorted = false, int mode = -1)
        {
            if (xValues == null || yValues == null)
            {
                BCException.Throw("BCGaussLobattoSpline.Init: values are null");
            }

            if (xValues.Count != yValues.Count)
            {
                BCException.Throw("BCGaussLobattoSpline.Init: vector sizes differ");
            }

            if (mode > -1)
            {
                _mode = mode;
            }

            var x = new List<double>(xValues);
            var y = new List<double>(yValues);
            if (!sorted)
            {
                SortPoints(x, y);
            }

            _xmin = x[0];
            _xmax = x[x.Count - 1];
            _data.Clear();
            _data.AddRange(y);

            var externalIndex = 0;
            var xLeft = x[externalIndex];
            var xRight = x[externalIndex + 1];
            var slope = (y[externalIndex + 1] - y[externalIndex]) / (xRight - xLeft);
            var result = new List<double>(_data.Count);
            for (var i = 0; i < _data.Count; i++)
            {
                var sampleX = XOfIndex(i);
                var intervalChanged = false;
                while (sampleX > x[externalIndex + 1] && externalIndex < _data.Count - 2)
                {
                    externalIndex++;
                    intervalChanged = true;
                }

                if (intervalChanged)
                {
                    xLeft = x[externalIndex];
                    xRight = x[externalIndex + 1];
                    slope = (y[externalIndex + 1] - y[externalIndex]) / (xRight - xLeft);
                }

                result.Add(y[externalIndex] + slope * (sampleX - x[externalIndex]));
            }

            _data.Clear();
            _data.AddRange(result);
        }

        public void Resize(int newSize)
        {
            if (newSize < 1)
            {
                BCException.Throw($"BCGaussLobattoSpline.Resize: invalid size={newSize}");
            }

            while (_data.Count < newSize)
            {
                _data.Add(0.0);
            }

            if (_data.Count > newSize)
            {
                _data.RemoveRange(newSize, _data.Count - newSize);
            }
        }

        public void SetInterval(double lowerBound, double upperBound)
        {
            if (lowerBound < upperBound && _data.Count > 0)
            {
                _xmin = lowerBound;
                _xmax = upperBound;
            }
        }

        public double XOfIndex(int index)
        {
            var xs = Pi * index / (_data.Count - 1.0);
            if (_mode > 0)
            {
                xs *= 0.5;
            }

            double fu;
            if (_mode == 2)
            {
                fu = Math.Sin(xs);
            }
            else
            {
                fu = 1.0 - Math.Cos(xs);
                if (_mode == 0)
                {
                    fu *= 0.5;
                }
            }

            return _xmin + (_xmax - _xmin) * fu;
        }

        public double IndexOfX(double x)
        {
            var fu = (x - _xmin) / (_xmax - _xmin);
            fu = Math.Min(1.0, Math.Max(0.0, fu));
            double xs;
            if (_mode == 2)
            {
                xs = Math.Asin(fu);
            }
            else
            {
                xs = _mode == 1 ? Math.Acos(1.0 - fu) : Math.Acos(1.0 - 2.0 * fu);
            }

            var result = xs / Pi * (_data.Count - 1.0);
            if (_mode > 0)
            {
                result *= 2.0;
            }

            return result;
        }

        private static void SortPoints(List<double> xValues, List<double> yValues)
        {
            var combined = new List<KeyValuePair<double, double>>(xValues.Count);
            for (var i = 0; i < xValues.Count; i++)
            {
                combined.Add(new KeyValuePair<double, double>(xValues[i], yValues[i]));
            }

            combined.Sort((a, b) => a.Key.CompareTo(b.Key));
            xValues.Clear();
            yValues.Clear();
            foreach (var pair in combined)
            {
                xValues.Add(pair.Key);
                yValues.Add(pair.Value);
            }
        }

        private static void EnsureInRange(int value, int min, int max, string label)
        {
            if (value < min || value > max)
            {
                BCException.Throw($"BCGaussLobattoSpline: {label} out of range");
            }
        }
    }
}
