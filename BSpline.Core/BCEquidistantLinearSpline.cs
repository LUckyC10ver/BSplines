using System;
using System.Collections.Generic;

namespace BSpline.Core
{
    public sealed class EquidistantLinearSpline
    {
        private const double EpsilonScale = 1e-11;

        private readonly List<double> _data = new List<double>();
        private double _xmin;
        private double _xmax;
        private double _dx;
        private double _epsilon;

        public EquidistantLinearSpline()
        {
            _xmin = 0.0;
            _xmax = 0.0;
            _dx = 1.0;
            _epsilon = EpsilonScale;
        }

        public EquidistantLinearSpline(IList<double> xValues, IList<double> yValues, bool sorted = false)
        {
            Init(xValues, yValues, sorted);
        }

        public string ClassName => "EquidistantLinearSpline";

        public bool Empty()
        {
            return _data.Count == 0;
        }

        public IReadOnlyList<double> GetDataVector()
        {
            return _data.AsReadOnly();
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

        public double Evaluate(double x)
        {
            if (_data.Count < 2)
            {
                throw new Exception("EquidistantLinearSpline.Evaluate: value size less than 2");
            }

            if (x < _xmin - _epsilon || x > _xmax + _epsilon)
            {
                throw new Exception($"EquidistantLinearSpline.Evaluate: input value {x} out of range [{_xmin},{_xmax}]");
            }

            var index = (int)((x - _xmin) / _dx);
            if (index <= -1)
            {
                return _data[0];
            }

            if (_data.Count - 1 <= index)
            {
                return _data[_data.Count - 1];
            }

            var xLeft = _xmin + index * _dx;
            var slope = (_data[index + 1] - _data[index]) / _dx;
            return _data[index] + slope * (x - xLeft);
        }

        public bool CheckEquidistance(IList<double> values)
        {
            if (values == null || values.Count < 2)
            {
                return false;
            }

            if (values.Count == 2)
            {
                return true;
            }

            var delta = values[1] - values[0];
            for (var i = 2; i < values.Count; i++)
            {
                var diff = values[i] - values[i - 1];
                if (Math.Abs(diff - delta) > 1e-10)
                {
                    return false;
                }
            }

            return true;
        }

        public double GetIntegral(double a, double b)
        {
            if (a < _xmin || a > _xmax)
            {
                throw new Exception($"EquidistantLinearSpline.GetIntegral: input value a={a} out of range [{_xmin},{_xmax}]");
            }

            if (b < _xmin || b > _xmax)
            {
                throw new Exception($"EquidistantLinearSpline.GetIntegral: input value b={b} out of range [{_xmin},{_xmax}]");
            }

            if (_data.Count < 2)
            {
                throw new Exception("EquidistantLinearSpline.GetIntegral: value size less than 2");
            }

            var sign = 1.0;
            if (b < a)
            {
                var temp = a;
                a = b;
                b = temp;
                sign = -1.0;
            }

            var minIndex = (a - _xmin - _dx / 2.0) / _dx;
            var maxIndex = (b - _xmin - _dx / 2.0) / _dx;
            var indexOfA = minIndex <= 0.0 ? 0 : (int)Math.Floor(minIndex);
            var indexOfB = maxIndex <= 0.0 ? 0 : (int)Math.Floor(maxIndex);

            double result;
            if (indexOfA == indexOfB)
            {
                var xLeft = _xmin + indexOfA * _dx;
                var slope = (_data[indexOfA + 1] - _data[indexOfA]) / _dx;
                var ya = _data[indexOfA] + slope * (a - xLeft);
                var yb = _data[indexOfA] + slope * (b - xLeft);
                result = 0.5 * (ya + yb) * (b - a);
            }
            else
            {
                result = 0.0;
                var xa = _xmin + indexOfA * _dx;
                var xb = _xmin + indexOfB * _dx;

                var slope = (_data[indexOfA + 1] - _data[indexOfA]) / _dx;
                var y0 = _data[indexOfA] + slope * (a - xa);
                result += 0.5 * (y0 + _data[indexOfA + 1]) * (xa + _dx - a);

                if (indexOfB < _data.Count - 1)
                {
                    slope = (_data[indexOfB + 1] - _data[indexOfB]) / _dx;
                }
                else
                {
                    slope = (_data[indexOfB] - _data[indexOfB - 1]) / _dx;
                }

                y0 = _data[indexOfB] + slope * (b - xb);
                result += 0.5 * (y0 + _data[indexOfB]) * (b - xb);

                if (indexOfB - indexOfA > 1)
                {
                    var sum = 0.5 * (_data[indexOfA + 1] + _data[indexOfB]);
                    for (var i = indexOfA + 2; i < indexOfB; i++)
                    {
                        sum += _data[i];
                    }

                    result += sum * _dx;
                }
            }

            return result * sign;
        }

        public void Init(IList<double> xValues, IList<double> yValues, bool sorted = false)
        {
            if (xValues == null || yValues == null)
            {
                throw new Exception("EquidistantLinearSpline.Init: values are null");
            }

            if (xValues.Count != yValues.Count)
            {
                throw new Exception("EquidistantLinearSpline.Init: vector sizes differ");
            }

            var x = new List<double>(xValues);
            var y = new List<double>(yValues);
            if (!sorted)
            {
                SortPoints(x, y);
            }

            if (!CheckEquidistance(x))
            {
                throw new Exception("EquidistantLinearSpline.Init: x values not equidistant");
            }

            _xmin = x[0];
            _xmax = x[x.Count - 1];
            _dx = (_xmax - _xmin) / (x.Count - 1);
            _epsilon = EpsilonScale * Math.Max(Math.Abs(_xmin), Math.Abs(_xmax));

            _data.Clear();
            _data.AddRange(y);
        }

        public void Init(int newSize, double initializer, double xmin, double xmax)
        {
            if (xmin > xmax)
            {
                throw new Exception($"EquidistantLinearSpline.Init: invalid range [xmin,xmax]=[{xmin},{xmax}]");
            }

            if (newSize < 1)
            {
                throw new Exception($"EquidistantLinearSpline.Init: invalid size={newSize}");
            }

            _xmin = xmin;
            _xmax = xmax;
            _dx = (_xmax - _xmin) / (newSize - 1);
            _epsilon = EpsilonScale * Math.Max(Math.Abs(_xmin), Math.Abs(_xmax));

            _data.Clear();
            for (var i = 0; i < newSize; i++)
            {
                _data.Add(initializer);
            }
        }

        public void Resize(int newSize)
        {
            if (newSize < 1)
            {
                throw new Exception($"EquidistantLinearSpline.Resize: invalid size={newSize}");
            }

            if (_data.Count == 0)
            {
                _data.Capacity = newSize;
            }
            else
            {
                while (_data.Count < newSize)
                {
                    _data.Add(0.0);
                }

                if (_data.Count > newSize)
                {
                    _data.RemoveRange(newSize, _data.Count - newSize);
                }
            }

            _dx = (_xmax - _xmin) / (newSize - 1);
        }

        public void SetInterval(double lowerBound, double upperBound)
        {
            if (lowerBound < upperBound && _data.Count > 0)
            {
                _xmin = lowerBound;
                _xmax = upperBound;
                _dx = (_xmax - _xmin) / (_data.Count - 1);
                _epsilon = EpsilonScale * Math.Max(Math.Abs(_xmin), Math.Abs(_xmax));
            }
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
    }
}
