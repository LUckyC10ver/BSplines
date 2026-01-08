using System;
using System.Collections.Generic;

namespace BSpline.Core
{
    public class BCSpline
    {
        private const double EpsilonScale = 1e-11;

        private bool _changed;
        private readonly List<double> _xValues = new List<double>();
        private readonly List<double> _yValues = new List<double>();
        private double _epsilon = EpsilonScale;

        public BCSpline()
        {
        }

        public BCSpline(IList<double> points, double value = 0.0, bool sorted = false)
        {
            Init(points, value, sorted);
        }

        public BCSpline(IList<double> points, IList<double> values, bool sorted = false)
        {
            Init(points, values, sorted);
        }

        public string ClassName => "BCSpline";

        public IReadOnlyList<double> GetXValues()
        {
            return _xValues.AsReadOnly();
        }

        public IReadOnlyList<double> GetYValues()
        {
            return _yValues.AsReadOnly();
        }

        public List<double> GetYValuesMutable()
        {
            return _yValues;
        }

        public int Size()
        {
            return _xValues.Count;
        }

        public bool IsChanged()
        {
            return _changed;
        }

        public void SetChanged()
        {
            _changed = true;
        }

        public void ResetChanged()
        {
            _changed = false;
        }

        public double XMin()
        {
            if (_xValues.Count == 0)
            {
                BCException.Throw("BCSpline.XMin(): Range not set");
            }

            return _xValues[0];
        }

        public double XMax()
        {
            if (_xValues.Count == 0)
            {
                BCException.Throw("BCSpline.XMax(): Range not set");
            }

            return _xValues[_xValues.Count - 1];
        }

        public int IndexFromValue(double x)
        {
            var n = _xValues.Count;
            if (n == 0)
            {
                return -1;
            }

            var klo = 0;
            var khi = n;
            while (khi - klo > 1)
            {
                var k = (khi + klo) / 2;
                if (_xValues[k] > x)
                {
                    khi = k;
                }
                else
                {
                    klo = k;
                }
            }

            return klo;
        }

        public void Init(IList<double> points, double value = 0.0, bool sorted = false)
        {
            if (points == null)
            {
                BCException.Throw("BCSpline.Init(points): points is null");
            }

            _xValues.Clear();
            _yValues.Clear();
            _xValues.AddRange(points);
            for (var i = 0; i < points.Count; i++)
            {
                _yValues.Add(value);
            }

            if (!sorted)
            {
                SortValues();
            }

            SetChanged();
            UpdateEpsilon();
        }

        public void Init(IList<double> points, IList<double> values, bool sorted = false)
        {
            if (points == null || values == null)
            {
                BCException.Throw("BCSpline.Init(points,values): points or values is null");
            }

            if (points.Count != values.Count)
            {
                BCException.Throw("BCSpline.Init(points,values): vector sizes differ");
            }

            _xValues.Clear();
            _yValues.Clear();
            _xValues.AddRange(points);
            _yValues.AddRange(values);

            if (!sorted)
            {
                SortValues();
            }

            SetChanged();
            UpdateEpsilon();
        }

        public void SetPoints(IList<double> points, double value = 0.0, bool sorted = false)
        {
            Init(points, value, sorted);
        }

        public void SetPoints(IList<double> points, IList<double> values, bool sorted = false)
        {
            Init(points, values, sorted);
        }

        public void SetValues(IList<double> values)
        {
            if (values == null)
            {
                BCException.Throw("BCSpline.SetValues: values is null");
            }

            if (_yValues.Count != values.Count)
            {
                BCException.Throw("BCSpline.SetValues: vector sizes differ");
            }

            _yValues.Clear();
            _yValues.AddRange(values);
            SetChanged();
        }

        public virtual void Update()
        {
            if (IsChanged())
            {
                ResetChanged();
            }
        }

        protected void SortValues()
        {
            var combined = new List<KeyValuePair<double, double>>(_xValues.Count);
            for (var i = 0; i < _xValues.Count; i++)
            {
                combined.Add(new KeyValuePair<double, double>(_xValues[i], _yValues[i]));
            }

            combined.Sort((a, b) => a.Key.CompareTo(b.Key));
            _xValues.Clear();
            _yValues.Clear();
            foreach (var pair in combined)
            {
                _xValues.Add(pair.Key);
                _yValues.Add(pair.Value);
            }
        }

        protected double GetEpsilon()
        {
            return _epsilon;
        }

        private void UpdateEpsilon()
        {
            if (_xValues.Count > 0)
            {
                var absMin = Math.Abs(_xValues[0]);
                var absMax = Math.Abs(_xValues[_xValues.Count - 1]);
                _epsilon = EpsilonScale * Math.Max(absMin, absMax);
            }
            else
            {
                _epsilon = EpsilonScale;
            }
        }
    }
}
