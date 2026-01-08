using System;
using System.IO;

namespace BSpline.Core
{
public sealed class XBSTools1Data<TInterval, TOrder>
        where TInterval : struct
        where TOrder : struct
    {
        private readonly BSToolsData _x = new BSToolsData();
        private int _lTpN;
        private int _uTpN;
        private int _numberOfW;
        private readonly double[] _w;

        public XBSTools1Data(int sizeOfW)
        {
            _w = new double[sizeOfW];
        }

        public double[] GetWeights() => _w;

        public int GetWeights(int maxElements, double[] elements)
        {
            if (_numberOfW <= maxElements)
            {
                return _x.Set(_numberOfW, _w, maxElements, elements);
            }

            return 0;
        }

        public int SetWeights(int numberElements, double[] elements)
        {
            if (numberElements <= _w.Length)
            {
                _numberOfW = numberElements;
                return _x.Set(numberElements, elements, _numberOfW, _w);
            }

            return 0;
        }

        public int GetOrderOfLowerTaylorPolynomial() => _lTpN;

        public void SetOrderOfLowerTaylorPolynomial(int order)
        {
            _lTpN = order;
        }

        public int GetOrderOfUpperTaylorPolynomial() => _uTpN;

        public void SetOrderOfUpperTaylorPolynomial(int order)
        {
            _uTpN = order;
        }

        public int GetCarrier(BSTools carrier)
        {
            return carrier.SetBinaryData(_x);
        }

        public int SetCarrier(BSTools carrier)
        {
            return carrier.GetBinaryData(_x);
        }

        public void Print(TextWriter writer, string path, string name)
        {
            if (writer == null)
            {
                return;
            }

            var parser = new SimpleParser("XBSTools1Data.Print", writer: writer);
            if (!parser.SetPath(path, name))
            {
                return;
            }

            parser.Write("a", _x.A);
            parser.Write("b", _x.B);
            parser.Write("n", _x.N);
            parser.Write("k", _x.K);
            parser.Write("lTp", _lTpN);
            parser.Write("uTp", _uTpN);
            parser.Write("w", _w, _numberOfW, "weights");
        }

        public bool IsValid(int maxInterval, int maxOrder)
        {
            if (BSToolsConstants.GetSizeTaylor(maxOrder) < _lTpN)
            {
                return false;
            }

            if (BSToolsConstants.GetSizeTaylor(maxOrder) < _uTpN)
            {
                return false;
            }

            return _x.IsValid(maxInterval, maxOrder);
        }
    }

    public abstract class BSTools1
    {
        protected readonly BSTools X;
        protected readonly int SizeOfW;
        protected readonly double[] W;
        protected readonly int[] WeightSelection;
        protected readonly int SizeOfWeightSelection;
        protected readonly int SizeOfLTp;
        protected readonly double[] LTp;
        protected readonly int SizeOfUTp;
        protected readonly double[] UTp;

        protected int LTpN;
        protected int UTpN;

        protected BSTools1(
            BSTools x,
            int sizeOfW,
            double[] w,
            int sizeOfLTp,
            double[] lTp,
            int sizeOfUTp,
            double[] uTp,
            int sizeOfWeightSelection,
            int[] weightSelection)
        {
            X = x;
            SizeOfW = sizeOfW;
            W = w;
            SizeOfLTp = sizeOfLTp;
            LTp = lTp;
            SizeOfUTp = sizeOfUTp;
            UTp = uTp;
            SizeOfWeightSelection = sizeOfWeightSelection;
            WeightSelection = weightSelection;
        }

        public void GetDomainSubWeightIndexes(double a, double b, out int ia, out int ib)
        {
            X.GetDomainSubWeightIndexes(a, b, out ia, out ib);
        }

        public int GetCarrierEdge(char edge, double[] buffer, int bufferLen)
        {
            return X.GetCarrierEdge(edge, buffer, bufferLen);
        }

        public void Init(BSTools carrier, double[] d = null, int numberOfWeights = 0, char lowerBoundary = 's', char upperBoundary = 's')
        {
            if (carrier == null)
            {
                return;
            }

            X.InitInternal(carrier.LowerBoundary, carrier.UpperBoundary, carrier.IntervalCount, carrier.Order);
            LTpN = X.GetFGradInternal(lowerBoundary);
            UTpN = X.GetFGradInternal(upperBoundary);

            if (numberOfWeights != 0)
            {
                var n = X.GetWLen();
                if (numberOfWeights == n && d != null)
                {
                    Array.Copy(d, W, n);
                }

                CompleteIt();
            }
        }

        public void CompleteIt()
        {
            X.BuildEquidistantCarrierInternal();

            for (var i = X.GetWLen(); i < SizeOfW; i++)
            {
                W[i] = 0.0;
            }

            X.TaylorPolyInternal(W, X.LowerBoundary, LTpN, LTp);
            var r = X.GetMainIndexInternal(X.GetIntervalIndexInternal(X.UpperBoundary));
            var slice = Slice(W, r, W.Length - r);
            X.TaylorPolyInternal(slice, X.UpperBoundary, UTpN, UTp);
        }

        public BSTools1 Multiply(double value)
        {
            var n = X.GetWLen();
            for (var i = 0; i < n; i++)
            {
                W[i] *= value;
            }

            CompleteIt();
            return this;
        }

        public double Evaluate(double x)
        {
            if (x < X.LowerBoundary)
            {
                return X.PolyValInternal(LTp, LTpN, x - X.LowerBoundary);
            }

            if (x <= X.UpperBoundary)
            {
                var index = X.GetIntervalIndexInternal(x);
                switch (X.Order)
                {
                    case 1:
                        return W[index];
                    case 2:
                        {
                            var r = X.GetMainIndexK2Internal(index);
                            var wSlice = Slice(W, r, W.Length - r);
                            var tSlice = Slice(X.Carrier, r, X.Carrier.Length - r);
                            return X.GetFctValueK2Internal(x, wSlice, tSlice);
                        }
                    case 3:
                        {
                            var r = X.GetMainIndexK3Internal(index);
                            var wSlice = Slice(W, r, W.Length - r);
                            var tSlice = Slice(X.Carrier, r, X.Carrier.Length - r);
                            return X.GetFctValueK3Internal(x, wSlice, tSlice);
                        }
                    default:
                        {
                            var r = X.GetMainIndexInternal(index);
                            var wSlice = Slice(W, r, W.Length - r);
                            var tSlice = Slice(X.Carrier, r, X.Carrier.Length - r);
                            return X.GetFctValueKnInternal(x, X.Order, wSlice, tSlice);
                        }
                }
            }

            return X.PolyValInternal(UTp, UTpN, x - X.UpperBoundary);
        }

        public double Derivative(double x)
        {
            var result = new double[1];
            Derivative(new[] { x }, 1, result);
            return result[0];
        }

        public double SecondDerivative(double x)
        {
            var result = new double[1];
            SecondDerivative(new[] { x }, 1, result);
            return result[0];
        }

        public double Inverse(double y)
        {
            var x1 = X.LowerBoundary;
            var x2 = X.UpperBoundary;
            var f1 = Evaluate(x1);
            var f2 = Evaluate(x2);

            if (f2 > f1)
            {
                if (y >= f2)
                {
                    return X.UpperBoundary;
                }

                if (y <= f1)
                {
                    return X.LowerBoundary;
                }
            }
            else
            {
                if (y >= f1)
                {
                    return X.LowerBoundary;
                }

                if (y <= f2)
                {
                    return X.UpperBoundary;
                }
            }

            f1 -= y;
            f2 -= y;

            const double machineEpsilon = 2.2204460492503131e-16;
            var epsilonX = machineEpsilon * Math.Abs(X.UpperBoundary - X.LowerBoundary);
            var epsilonF = BSToolsConstants.Delta * y + epsilonX;

            if (Math.Abs(f1) < epsilonF)
            {
                return x1;
            }

            if (Math.Abs(f2) < epsilonF)
            {
                return x2;
            }

            var z = 0.5 * (x1 + x2);
            var fz = Evaluate(z) - y;
            if (f2 * fz < 0.0)
            {
                x1 = x2;
                f1 = f2;
            }
            else
            {
                f1 *= 0.5;
            }

            x2 = z;
            f2 = fz;

            for (var i = 0; i < 30; i++)
            {
                if (Math.Abs(x2 - x1) < epsilonX || Math.Abs(f2) < epsilonF)
                {
                    break;
                }

                z = x1 - f1 * (x2 - x1) / (f2 - f1);
                fz = Evaluate(z) - y;

                if (fz * f2 < 0.0)
                {
                    x1 = x2;
                    f1 = f2;
                }
                else
                {
                    f1 *= 0.5;
                }

                x2 = z;
                f2 = fz;
            }

            var f2Prime = Derivative(x2);
            x2 -= f2 / f2Prime;

            if (epsilonF < Math.Abs(f2))
            {
                _ = Evaluate(x2) - y;
            }

            return x2;
        }

        public void Inverse(double[] x, int dim, double[] y)
        {
            for (var i = 0; i < dim; i++)
            {
                y[i] = Inverse(x[i]);
            }
        }

        public void Evaluate(double[] x, int dim, double[] y)
        {
            for (var i = 0; i < dim; i++)
            {
                y[i] = Evaluate(x[i]);
            }
        }

        public abstract void Derivative(double[] vx, int dim, double[] y);

        public abstract void SecondDerivative(double[] vx, int dim, double[] y);

        public BSTools1 Diff()
        {
            if (X.Order == 1)
            {
                for (var i = 0; i < SizeOfW; i++)
                {
                    W[i] = 0;
                }

                LTpN = 0;
                UTpN = 0;
            }
            else
            {
                X.DiffWeightsInternal(W, W);
                X.InitInternal(X.LowerBoundary, X.UpperBoundary, X.IntervalCount, X.Order - 1);

                for (var i = X.GetWLen(); i < SizeOfW; i++)
                {
                    W[i] = 0;
                }

                if (LTpN > 0)
                {
                    LTpN--;
                }

                if (UTpN > 0)
                {
                    UTpN--;
                }
            }

            CompleteIt();
            return this;
        }

        public BSTools1 Integral(double x = 0, double y = 0)
        {
            if (X.MaxOrder <= X.Order)
            {
                return this;
            }

            var n = X.GetWLen();
            X.InitInternal(X.LowerBoundary, X.UpperBoundary, X.IntervalCount, X.Order + 1);
            for (var i = n; i > 0; i--)
            {
                W[i] = W[i - 1];
            }

            var temp = new double[n];
            Array.Copy(W, 1, temp, 0, n);
            X.IntWeightsInternal(temp, W);
            CompleteIt();

            var c = y - Evaluate(x);
            if (c != 0)
            {
                for (var i = 0; i < n + 1; i++)
                {
                    W[i] += c;
                }
            }

            LTpN++;
            UTpN++;
            CompleteIt();
            return this;
        }

        public int GetNumberOfWeights()
        {
            return X.GetWLen();
        }

        public int GetNumberOfAllWeights()
        {
            return X.GetWLen();
        }

        public double GetWeightPosition(int index)
        {
            return X.GetWeightPosition(index);
        }

        public int GetWeightPosition(double[] buffer, int bufferLength)
        {
            var max = X.GetWLen();
            if (bufferLength < max)
            {
                max = bufferLength;
            }

            for (var i = 0; i < max; i++)
            {
                buffer[i] = X.GetWeightPosition(i);
            }

            return max;
        }

        public void SetWeight(int index, double value)
        {
            if (index < X.GetWLen())
            {
                W[index] = value;
            }
        }

        public void SetWeight(double[] values, int numberOfValues)
        {
            if (values == null)
            {
                return;
            }

            if (numberOfValues <= X.GetWLen())
            {
                Array.Copy(values, W, numberOfValues);
            }
        }

        public double GetWeight(int index)
        {
            if (index <= X.GetWLen())
            {
                return W[index];
            }

            return 0.0;
        }

        public double GetWeightFromAll(int index)
        {
            if (index <= X.GetWLen())
            {
                return W[index];
            }

            return 0.0;
        }

        public int GetWeight(double[] buffer, int lengthOfBuffer)
        {
            var max = Math.Min(lengthOfBuffer, X.GetWLen());
            Array.Copy(W, buffer, max);
            return max;
        }

        public void SetWeightMeanValue(double value)
        {
            var oldMean = GetWeightMeanValue();
            var delta = value - oldMean;
            var n = X.GetWLen();
            for (var i = 0; i < n; i++)
            {
                W[i] += delta;
            }
        }

        public double GetWeightMeanValue()
        {
            var mean = 0.0;
            var n = X.GetWLen();
            for (var i = 0; i < n; i++)
            {
                mean += W[i];
            }

            return mean / n;
        }

        public void SetWeightSelection(double[] values, int numberOfValues)
        {
            if (values == null || numberOfValues <= 0)
            {
                return;
            }

            var count = Math.Min(numberOfValues, SizeOfWeightSelection);
            for (var i = 0; i < count; i++)
            {
                WeightSelection[i] = (int)Math.Round(values[i]);
            }
        }

        public int Load(string filename)
        {
            if (string.IsNullOrEmpty(filename))
            {
                return 0;
            }

            using var reader = new StreamReader(filename);
            var parser = new SimpleParser("BSTools1.Load", reader: reader);
            return Load(parser) ? 1 : 0;
        }

        public bool Load(object parser)
        {
            if (parser is not SimpleParser simpleParser)
            {
                return false;
            }

            double? a = null;
            double? b = null;
            int? n = null;
            int? k = null;
            int? lTp = null;
            int? uTp = null;
            double[] weights = null;

            while (simpleParser.GetNewLine())
            {
                var line = simpleParser.GetData();
                var key = ExtractKey(line);
                if (string.IsNullOrEmpty(key))
                {
                    continue;
                }

                switch (key)
                {
                    case "a":
                        if (simpleParser.GetValue(out double aValue))
                        {
                            a = aValue;
                        }
                        break;
                    case "b":
                        if (simpleParser.GetValue(out double bValue))
                        {
                            b = bValue;
                        }
                        break;
                    case "n":
                        if (simpleParser.GetValue(out int nValue))
                        {
                            n = nValue;
                        }
                        break;
                    case "k":
                        if (simpleParser.GetValue(out int kValue))
                        {
                            k = kValue;
                        }
                        break;
                    case "lTp":
                        if (simpleParser.GetValue(out int lValue))
                        {
                            lTp = lValue;
                        }
                        break;
                    case "uTp":
                        if (simpleParser.GetValue(out int uValue))
                        {
                            uTp = uValue;
                        }
                        break;
                    case "w":
                        weights ??= new double[SizeOfW];
                        var count = simpleParser.GetArray(weights, weights.Length);
                        if (count == 0)
                        {
                            weights = null;
                        }
                        break;
                }
            }

            if (!a.HasValue || !b.HasValue || !n.HasValue || !k.HasValue)
            {
                return false;
            }

            X.InitInternal(a.Value, b.Value, n.Value, k.Value);
            LTpN = lTp ?? LTpN;
            UTpN = uTp ?? UTpN;
            if (weights != null)
            {
                Array.Copy(weights, W, Math.Min(weights.Length, W.Length));
            }

            CompleteIt();
            return true;
        }

        public int Save(string filename, string splineName = "data", string comment = "non")
        {
            if (string.IsNullOrEmpty(filename))
            {
                return 0;
            }

            using var writer = new StreamWriter(filename);
            return Save(writer, splineName, splineName) ? 1 : 0;
        }

        public bool Save(TextWriter writer, string path, string varName)
        {
            if (writer == null)
            {
                return false;
            }

            var parser = new SimpleParser("BSTools1.Save", writer: writer);
            if (!parser.SetPath(path, varName))
            {
                return false;
            }

            parser.Write("a", X.GetLowerBoundary());
            parser.Write("b", X.GetUpperBoundary());
            parser.Write("n", X.GetNumberOfIntervals());
            parser.Write("k", X.GetOrder());
            parser.Write("lTp", LTpN);
            parser.Write("uTp", UTpN);
            parser.Write("w", W, X.GetWLen(), "weights");
            return true;
        }

        public int SetCarrier(double a, double b, int n, int k, char lTp = 's', char uTp = 's')
        {
            X.InitInternal(a, b, n, k);
            Init(X, null, 0, lTp, uTp);
            return 0;
        }

        public int GetOrderOfLowerTaylorPolynomial() => LTpN;

        public void SetOrderOfLowerTaylorPolynomial(int order)
        {
            LTpN = order;
        }

        public int GetOrderOfUpperTaylorPolynomial() => UTpN;

        public void SetOrderOfUpperTaylorPolynomial(int order)
        {
            UTpN = order;
        }

        public void AddLine(double xFix, double xRef, double yRef)
        {
            if (xFix == xRef)
            {
                return;
            }

            var yOld = Evaluate(xRef);
            var m = (yRef - yOld) / (xRef - xFix);
            var b = ((yOld - yRef) * xFix) / (xRef - xFix);

            var n = GetNumberOfAllWeights();
            for (var i = 0; i < n; i++)
            {
                SetWeight(i, GetWeight(i) + m * GetWeightPosition(i) + b);
            }

            CompleteIt();
        }

        public bool IsApproximatedEqual(BSTools1 other, double maxDelta)
        {
            if (!X.IsApproximatedEqual(other.X, maxDelta)) return false;
            if (!BSTools.IsApproximatedEqual(W, X.GetWLen(), other.W, other.X.GetWLen(), maxDelta)) return false;
            if (maxDelta < Math.Abs(LTpN - other.LTpN)) return false;
            if (maxDelta < Math.Abs(UTpN - other.UTpN)) return false;
            return true;
        }

        public bool IsStaticDataEqual(BSTools1 other)
        {
            if (X.IsNotEqual(other.X)) return false;
            if (!BSTools.IsEqual(W, X.GetWLen(), other.W, other.X.GetWLen())) return false;
            if (!BSTools.IsEqual(LTp, LTpN, other.LTp, other.LTpN)) return false;
            if (!BSTools.IsEqual(UTp, UTpN, other.UTp, other.UTpN)) return false;
            return true;
        }

        public double LowerBound()
        {
            return X.LowerBoundary;
        }

        public double UpperBound()
        {
            return X.UpperBoundary;
        }

        public int GetNumberOfIntervals()
        {
            return X.GetNumberOfIntervals();
        }

        public int GetOrder()
        {
            return X.GetOrder();
        }

        public abstract double CalcWeightedNormOfBend(double weightsDirection = 1);

        public double CalcNormP(int p)
        {
            return CalcNormPPartial(p, 0, X.GetWLen());
        }

        protected double CalcNormPPartial(int p, int iMin, int iMax)
        {
            var ret = 0.0;
            for (var i = iMin; i <= iMax; i++)
            {
                var value = W[i];
                for (var j = 1; j < p; j++)
                {
                    value *= W[i];
                }

                ret += Math.Abs(value);
            }

            return ret;
        }

        private static double[] Slice(double[] source, int start, int length)
        {
            var maxLength = Math.Max(0, Math.Min(length, source.Length - start));
            var sliceLength = Math.Max(0, maxLength);
            var buffer = new double[sliceLength];
            if (sliceLength > 0)
            {
                Array.Copy(source, start, buffer, 0, sliceLength);
            }

            return buffer;
        }
    }

    public sealed class XBSTools1<TInterval, TOrder> : BSTools1
        where TInterval : struct
        where TOrder : struct
    {
        private readonly int _maxInterval;
        private readonly int _maxOrder;
        private readonly XBSTools<TInterval, TOrder> _xMem;
        private readonly double[] _wMem;
        private readonly double[] _lTpMem;
        private readonly double[] _uTpMem;

        public XBSTools1(int maxInterval, int maxOrder)
            : base(
                new XBSTools<TInterval, TOrder>(maxInterval, maxOrder),
                BSToolsConstants.GetSizeDeBoor(maxInterval, maxOrder),
                new double[BSToolsConstants.GetSizeDeBoor(maxInterval, maxOrder)],
                BSToolsConstants.GetSizeTaylor(maxOrder),
                new double[BSToolsConstants.GetSizeTaylor(maxOrder)],
                BSToolsConstants.GetSizeTaylor(maxOrder),
                new double[BSToolsConstants.GetSizeTaylor(maxOrder)],
                1,
                new int[1])
        {
            _maxInterval = maxInterval;
            _maxOrder = maxOrder;
            _xMem = (XBSTools<TInterval, TOrder>)X;
            _wMem = W;
            _lTpMem = LTp;
            _uTpMem = UTp;
        }

        public XBSTools1(int maxInterval, int maxOrder, BSTools input, double[] d = null, int numberOfWeights = 0, char lowerBoundary = 's', char upperBoundary = 's')
            : this(maxInterval, maxOrder)
        {
            Init(input, d, numberOfWeights, lowerBoundary, upperBoundary);
        }

        public override void Derivative(double[] vx, int dim, double[] y)
        {
            if (vx == null || y == null)
            {
                return;
            }

            var bs = new XBSTools1<TInterval, TOrder>(_maxInterval, _maxOrder, X, W, X.GetWLen());
            bs.Diff();
            bs.Evaluate(vx, dim, y);
        }

        public override void SecondDerivative(double[] vx, int dim, double[] y)
        {
            if (vx == null || y == null)
            {
                return;
            }

            var bs = new XBSTools1<TInterval, TOrder>(_maxInterval, _maxOrder, X, W, X.GetWLen());
            bs.Diff();
            bs.Diff();
            bs.Evaluate(vx, dim, y);
        }

        public int GetBinaryData(XBSTools1Data<TInterval, TOrder> data)
        {
            if (data == null)
            {
                return 0;
            }

            data.SetCarrier(X);
            data.SetOrderOfLowerTaylorPolynomial(LTpN);
            data.SetOrderOfUpperTaylorPolynomial(UTpN);
            data.SetWeights(X.GetWLen(), W);
            return 1;
        }

        public int SetBinaryData(XBSTools1Data<TInterval, TOrder> data)
        {
            if (data == null)
            {
                return 0;
            }

            data.GetCarrier(X);
            LTpN = data.GetOrderOfLowerTaylorPolynomial();
            UTpN = data.GetOrderOfUpperTaylorPolynomial();
            var weights = data.GetWeights();
            Array.Copy(weights, W, Math.Min(weights.Length, W.Length));
            CompleteIt();
            return 1;
        }

        public override double CalcWeightedNormOfBend(double weightsDirection = 1)
        {
            var n = X.GetWLen();
            if (n <= 0)
            {
                return 0.0;
            }

            var maxValue = 0.0;
            for (var i = 0; i < n; i++)
            {
                var value = Math.Abs(W[i] * weightsDirection);
                if (value > maxValue)
                {
                    maxValue = value;
                }
            }

            return maxValue;
        }

        private static string ExtractKey(string line)
        {
            if (string.IsNullOrWhiteSpace(line))
            {
                return string.Empty;
            }

            var equalsIndex = line.IndexOf('=');
            if (equalsIndex <= 0)
            {
                return string.Empty;
            }

            var left = line.Substring(0, equalsIndex).Trim();
            var dotIndex = left.LastIndexOf('.');
            return dotIndex >= 0 ? left.Substring(dotIndex + 1).Trim() : left;
        }
    }

    public abstract class BSTools2
    {
        protected sealed class ArrayBuffer
        {
            private readonly double[] _data;
            private BSTools2 _owner;

            public ArrayBuffer(int sizeOfD, double[] data)
            {
                _data = data;
            }

            public void SetBasicPointer(BSTools2 owner)
            {
                _owner = owner;
            }

            public double[] GetBuffer()
            {
                return _data;
            }

        public int GetLength()
        {
            return _data.Length;
        }

        public bool IsApproximatedEqual(ArrayBuffer other, double maxDelta)
        {
            return BSTools.IsApproximatedEqual(_data, GetLength(), other._data, other.GetLength(), maxDelta);
        }

        public bool IsNotEqual(ArrayBuffer other)
        {
            return !BSTools.IsEqual(_data, GetLength(), other._data, other.GetLength());
        }
        }

        protected readonly ArrayBuffer Weights;
        protected readonly BSTools X;
        protected readonly BSTools Y;
        protected readonly int SizeOfWeightSelection;
        protected readonly int[] WeightSelection;
        protected readonly int SizeOfChangedWeight;
        protected readonly int[] ChangedWeight;
        protected int LTpX;
        protected int UTpX;
        protected int LTpY;
        protected int UTpY;

        protected BSTools2(
            int sizeOfD,
            double[] d,
            BSTools x,
            BSTools y,
            int sizeOfWeightSelection,
            int[] weightSelection,
            int sizeOfChangedWeight,
            int[] changedWeight)
        {
            Weights = new ArrayBuffer(sizeOfD, d);
            X = x;
            Y = y;
            SizeOfWeightSelection = sizeOfWeightSelection;
            WeightSelection = weightSelection;
            SizeOfChangedWeight = sizeOfChangedWeight;
            ChangedWeight = changedWeight;
            Array.Clear(WeightSelection, 0, WeightSelection.Length);
            Array.Clear(ChangedWeight, 0, ChangedWeight.Length);
            Weights.SetBasicPointer(this);
        }

        public abstract BSTools2 DiffX();

        public abstract BSTools2 DiffY();

        public abstract BSTools2 IntX(double x = 0, double z = 0);

        public abstract BSTools2 IntY(double y = 0, double z = 0);

        public int Load(string filename)
        {
            if (string.IsNullOrEmpty(filename))
            {
                return 0;
            }

            using var reader = new StreamReader(filename);
            var parser = new SimpleParser("BSTools2.Load", reader: reader);
            return Load(parser) ? 1 : 0;
        }

        public bool Load(object parser)
        {
            if (parser is not SimpleParser simpleParser)
            {
                return false;
            }

            double? xA = null;
            double? xB = null;
            int? xN = null;
            int? xK = null;
            double? yA = null;
            double? yB = null;
            int? yN = null;
            int? yK = null;
            int? lTpX = null;
            int? uTpX = null;
            int? lTpY = null;
            int? uTpY = null;
            double[] weights = null;

            while (simpleParser.GetNewLine())
            {
                var line = simpleParser.GetData();
                var (prefix, key) = ExtractPrefixKey(line);
                if (string.IsNullOrEmpty(key))
                {
                    continue;
                }

                if (prefix == "x")
                {
                    switch (key)
                    {
                        case "a":
                            if (simpleParser.GetValue(out double valueA))
                            {
                                xA = valueA;
                            }
                            break;
                        case "b":
                            if (simpleParser.GetValue(out double valueB))
                            {
                                xB = valueB;
                            }
                            break;
                        case "n":
                            if (simpleParser.GetValue(out int valueN))
                            {
                                xN = valueN;
                            }
                            break;
                        case "k":
                            if (simpleParser.GetValue(out int valueK))
                            {
                                xK = valueK;
                            }
                            break;
                    }
                }
                else if (prefix == "y")
                {
                    switch (key)
                    {
                        case "a":
                            if (simpleParser.GetValue(out double valueA))
                            {
                                yA = valueA;
                            }
                            break;
                        case "b":
                            if (simpleParser.GetValue(out double valueB))
                            {
                                yB = valueB;
                            }
                            break;
                        case "n":
                            if (simpleParser.GetValue(out int valueN))
                            {
                                yN = valueN;
                            }
                            break;
                        case "k":
                            if (simpleParser.GetValue(out int valueK))
                            {
                                yK = valueK;
                            }
                            break;
                    }
                }

                switch (key)
                {
                    case "lTp_x":
                        if (simpleParser.GetValue(out int lValueX))
                        {
                            lTpX = lValueX;
                        }
                        break;
                    case "uTp_x":
                        if (simpleParser.GetValue(out int uValueX))
                        {
                            uTpX = uValueX;
                        }
                        break;
                    case "lTp_y":
                        if (simpleParser.GetValue(out int lValueY))
                        {
                            lTpY = lValueY;
                        }
                        break;
                    case "uTp_y":
                        if (simpleParser.GetValue(out int uValueY))
                        {
                            uTpY = uValueY;
                        }
                        break;
                    case "w":
                        weights ??= new double[Weights.GetLength()];
                        var count = simpleParser.GetArray(weights, weights.Length);
                        if (count == 0)
                        {
                            weights = null;
                        }
                        break;
                }
            }

            if (!xA.HasValue || !xB.HasValue || !xN.HasValue || !xK.HasValue ||
                !yA.HasValue || !yB.HasValue || !yN.HasValue || !yK.HasValue)
            {
                return false;
            }

            X.InitInternal(xA.Value, xB.Value, xN.Value, xK.Value);
            Y.InitInternal(yA.Value, yB.Value, yN.Value, yK.Value);
            LTpX = lTpX ?? LTpX;
            UTpX = uTpX ?? UTpX;
            LTpY = lTpY ?? LTpY;
            UTpY = uTpY ?? UTpY;
            if (weights != null)
            {
                Array.Copy(weights, Weights.GetBuffer(), Math.Min(weights.Length, Weights.GetLength()));
            }

            CompleteIt();
            return true;
        }

        public int Save(string filename, string splineName = "data", string comment = "non", string format4W = "%.16g")
        {
            if (string.IsNullOrEmpty(filename))
            {
                return 0;
            }

            using var writer = new StreamWriter(filename);
            return Save(writer, splineName, splineName, format4W) ? 1 : 0;
        }

        public bool Save(TextWriter writer, string path, string varName, string format4W = "%.16g")
        {
            if (writer == null)
            {
                return false;
            }

            var parser = new SimpleParser("BSTools2.Save", writer: writer);
            if (!parser.SetPath(path, varName))
            {
                return false;
            }

            X.Save(writer, parser.GetPath(), "x");
            Y.Save(writer, parser.GetPath(), "y");
            var weightFormat = NormalizeFormat(format4W);
            parser.Write("w", Weights.GetBuffer(), Weights.GetLength(), "weights", weightFormat);
            parser.Write("lTp_x", LTpX);
            parser.Write("uTp_x", UTpX);
            parser.Write("lTp_y", LTpY);
            parser.Write("uTp_y", UTpY);
            return true;
        }

        public double Evaluate(double x, double y)
        {
            var result = new double[1];
            Evaluate(new[] { x }, new[] { y }, 1, result);
            return result[0];
        }

        public double Dx(double x, double y)
        {
            var result = new double[1];
            Dx(new[] { x }, new[] { y }, 1, result);
            return result[0];
        }

        public double Dy(double x, double y)
        {
            var result = new double[1];
            Dy(new[] { x }, new[] { y }, 1, result);
            return result[0];
        }

        public void Evaluate(double[] vx, double[] vy, int dim, double[] result)
        {
            if (result == null || vx == null || vy == null || dim <= 0)
            {
                return;
            }

            var xOrder = X.Order;
            var yOrder = Y.Order;
            var wx = X.GetWLen();
            var weights = Weights.GetBuffer();
            var wy = new double[Math.Max(1, yOrder)];
            var taylor = new double[BSToolsConstants.BSTOOLS_TOTAL_SIZE_TAYLOR];

            for (var i = 0; i < dim; i++)
            {
                var ex = vx[i];
                var ey = vy[i];

                var posX = X.GetLocalInternal(ex);
                var posY = Y.GetLocalInternal(ey);

                var ry = posY switch
                {
                    BSTools.Position.Mitte => Y.GetMainIndexInternal(Y.GetIntervalIndexInternal(ey)),
                    BSTools.Position.Unter => Y.GetMainIndexInternal(Y.GetIntervalIndexInternal(Y.LowerBoundary)),
                    _ => Y.GetMainIndexInternal(Y.GetIntervalIndexInternal(Y.UpperBoundary))
                };

                if (posX == BSTools.Position.Mitte)
                {
                    var index = X.GetIntervalIndexInternal(ex);
                    var rx = xOrder switch
                    {
                        1 => index,
                        2 => X.GetMainIndexK2Internal(index),
                        3 => X.GetMainIndexK3Internal(index),
                        _ => X.GetMainIndexInternal(index)
                    };

                    for (var j = 0; j < yOrder; j++)
                    {
                        var rowOffset = (ry + j) * wx;
                        wy[j] = EvaluateRow(ex, rowOffset, rx, xOrder);
                    }
                }
                else
                {
                    var rx = X.GetMainIndexInternal(X.GetIntervalIndexInternal(posX == BSTools.Position.Unter ? X.LowerBoundary : X.UpperBoundary));
                    var tpN = posX == BSTools.Position.Unter ? LTpX : UTpX;
                    var boundary = posX == BSTools.Position.Unter ? X.LowerBoundary : X.UpperBoundary;

                    for (var j = 0; j < yOrder; j++)
                    {
                        var rowOffset = (ry + j) * wx;
                        var d = new double[xOrder];
                        Array.Copy(weights, rowOffset + rx, d, 0, xOrder);
                        X.TaylorPolyInternal(d, boundary, tpN, taylor);
                        wy[j] = X.PolyValInternal(taylor, tpN, ex - boundary);
                    }
                }

                result[i] = posY switch
                {
                    BSTools.Position.Mitte => EvaluateY(ey, ry, wy, yOrder),
                    BSTools.Position.Unter => EvaluateYTaylor(ey, Y.LowerBoundary, LTpY, wy),
                    _ => EvaluateYTaylor(ey, Y.UpperBoundary, UTpY, wy)
                };
            }
        }

        public void Dx(double[] vx, double[] vy, int dim, double[] result)
        {
            if (result == null || vx == null || vy == null || dim <= 0)
            {
                return;
            }

            var xOrder = X.Order;
            var yOrder = Y.Order;
            var wx = X.GetWLen();
            var weights = Weights.GetBuffer();
            var wy = new double[Math.Max(1, yOrder)];
            var taylor = new double[BSToolsConstants.BSTOOLS_TOTAL_SIZE_TAYLOR];

            for (var i = 0; i < dim; i++)
            {
                var ex = vx[i];
                var ey = vy[i];

                var posX = X.GetLocalInternal(ex);
                var posY = Y.GetLocalInternal(ey);

                var ry = posY switch
                {
                    BSTools.Position.Mitte => Y.GetMainIndexInternal(Y.GetIntervalIndexInternal(ey)),
                    BSTools.Position.Unter => Y.GetMainIndexInternal(Y.GetIntervalIndexInternal(Y.LowerBoundary)),
                    _ => Y.GetMainIndexInternal(Y.GetIntervalIndexInternal(Y.UpperBoundary))
                };

                if (posX == BSTools.Position.Mitte)
                {
                    var rx = X.GetMainIndexInternal(X.GetIntervalIndexInternal(ex));
                    for (var j = 0; j < yOrder; j++)
                    {
                        var rowOffset = (ry + j) * wx;
                        var d = new double[xOrder];
                        Array.Copy(weights, rowOffset + rx, d, 0, xOrder);
                        X.AblwertInternal(ex, xOrder, d, SliceCarrier(rx), out wy[j]);
                    }
                }
                else
                {
                    var rx = X.GetMainIndexInternal(X.GetIntervalIndexInternal(posX == BSTools.Position.Unter ? X.LowerBoundary : X.UpperBoundary));
                    var tpN = posX == BSTools.Position.Unter ? LTpX : UTpX;
                    var boundary = posX == BSTools.Position.Unter ? X.LowerBoundary : X.UpperBoundary;
                    for (var j = 0; j < yOrder; j++)
                    {
                        var rowOffset = (ry + j) * wx;
                        var d = new double[xOrder];
                        Array.Copy(weights, rowOffset + rx, d, 0, xOrder);
                        X.TaylorPolyInternal(d, boundary, tpN, taylor);
                        var tempN = tpN;
                        X.PolyDiffInternal(taylor, ref tempN);
                        wy[j] = X.PolyValInternal(taylor, tempN, ex - boundary);
                    }
                }

                result[i] = posY switch
                {
                    BSTools.Position.Mitte => EvaluateY(ey, ry, wy, yOrder),
                    BSTools.Position.Unter => EvaluateYTaylor(ey, Y.LowerBoundary, LTpY, wy),
                    _ => EvaluateYTaylor(ey, Y.UpperBoundary, UTpY, wy)
                };
            }
        }

        public void Dy(double[] vx, double[] vy, int dim, double[] result)
        {
            if (result == null || vx == null || vy == null || dim <= 0)
            {
                return;
            }

            var xOrder = X.Order;
            var yOrder = Y.Order;
            var wx = X.GetWLen();
            var weights = Weights.GetBuffer();
            var wy = new double[Math.Max(1, yOrder)];
            var taylor = new double[BSToolsConstants.BSTOOLS_TOTAL_SIZE_TAYLOR];

            for (var i = 0; i < dim; i++)
            {
                var ex = vx[i];
                var ey = vy[i];

                var posX = X.GetLocalInternal(ex);
                var posY = Y.GetLocalInternal(ey);

                var ry = posY switch
                {
                    BSTools.Position.Mitte => Y.GetMainIndexInternal(Y.GetIntervalIndexInternal(ey)),
                    BSTools.Position.Unter => Y.GetMainIndexInternal(Y.GetIntervalIndexInternal(Y.LowerBoundary)),
                    _ => Y.GetMainIndexInternal(Y.GetIntervalIndexInternal(Y.UpperBoundary))
                };

                var rx = posX == BSTools.Position.Mitte
                    ? X.GetMainIndexInternal(X.GetIntervalIndexInternal(ex))
                    : X.GetMainIndexInternal(X.GetIntervalIndexInternal(posX == BSTools.Position.Unter ? X.LowerBoundary : X.UpperBoundary));

                if (posX == BSTools.Position.Mitte)
                {
                    for (var j = 0; j < yOrder; j++)
                    {
                        var rowOffset = (ry + j) * wx;
                        wy[j] = EvaluateRow(ex, rowOffset, rx, xOrder);
                    }
                }
                else
                {
                    var tpN = posX == BSTools.Position.Unter ? LTpX : UTpX;
                    var boundary = posX == BSTools.Position.Unter ? X.LowerBoundary : X.UpperBoundary;
                    for (var j = 0; j < yOrder; j++)
                    {
                        var rowOffset = (ry + j) * wx;
                        var d = new double[xOrder];
                        Array.Copy(weights, rowOffset + rx, d, 0, xOrder);
                        X.TaylorPolyInternal(d, boundary, tpN, taylor);
                        wy[j] = X.PolyValInternal(taylor, tpN, ex - boundary);
                    }
                }

                result[i] = posY switch
                {
                    BSTools.Position.Mitte => EvaluateYDerivative(ey, ry, wy, yOrder),
                    BSTools.Position.Unter => EvaluateYDerivativeTaylor(ey, Y.LowerBoundary, LTpY, wy),
                    _ => EvaluateYDerivativeTaylor(ey, Y.UpperBoundary, UTpY, wy)
                };
            }
        }

        private double EvaluateRow(double ex, int rowOffset, int rx, int order)
        {
            var d = new double[order];
            Array.Copy(Weights.GetBuffer(), rowOffset + rx, d, 0, order);
            var t = SliceCarrier(rx);
            return order switch
            {
                1 => d[0],
                2 => X.GetFctValueK2Internal(ex, d, t),
                3 => X.GetFctValueK3Internal(ex, d, t),
                _ => X.GetFctValueKnInternal(ex, order, d, t)
            };
        }

        private double EvaluateY(double ey, int ry, double[] wy, int order)
        {
            var t = SliceCarrierY(ry);
            return order switch
            {
                1 => wy[0],
                2 => Y.GetFctValueK2Internal(ey, wy, t),
                3 => Y.GetFctValueK3Internal(ey, wy, t),
                _ => Y.GetFctValueInternal(ey, order, wy, t)
            };
        }

        private double EvaluateYTaylor(double ey, double boundary, int tpN, double[] wy)
        {
            var taylor = new double[BSToolsConstants.BSTOOLS_TOTAL_SIZE_TAYLOR];
            Y.TaylorPolyInternal(wy, boundary, tpN, taylor);
            return Y.PolyValInternal(taylor, tpN, ey - boundary);
        }

        private double EvaluateYDerivative(double ey, int ry, double[] wy, int order)
        {
            var t = SliceCarrierY(ry);
            Y.AblwertInternal(ey, order, wy, t, out var value);
            return value;
        }

        private double EvaluateYDerivativeTaylor(double ey, double boundary, int tpN, double[] wy)
        {
            var taylor = new double[BSToolsConstants.BSTOOLS_TOTAL_SIZE_TAYLOR];
            Y.TaylorPolyInternal(wy, boundary, tpN, taylor);
            var tempN = tpN;
            Y.PolyDiffInternal(taylor, ref tempN);
            return Y.PolyValInternal(taylor, tempN, ey - boundary);
        }

        private double[] SliceCarrier(int start)
        {
            var length = Math.Max(0, X.Carrier.Length - start);
            var buffer = new double[length];
            Array.Copy(X.Carrier, start, buffer, 0, length);
            return buffer;
        }

        private double[] SliceCarrierY(int start)
        {
            var length = Math.Max(0, Y.Carrier.Length - start);
            var buffer = new double[length];
            Array.Copy(Y.Carrier, start, buffer, 0, length);
            return buffer;
        }

        public double UpperBoundX()
        {
            return X.GetUpperBoundary();
        }

        public double LowerBoundX()
        {
            return X.GetLowerBoundary();
        }

        public double UpperBoundY()
        {
            return Y.GetUpperBoundary();
        }

        public double LowerBoundY()
        {
            return Y.GetLowerBoundary();
        }

        public int GetNumberOfWeights()
        {
            return Weights.GetLength();
        }

        public int GetNumberOfAllWeights()
        {
            return Weights.GetLength();
        }

        public int GetWeightPosition(int index, out double x, out double y)
        {
            if (Weights.GetLength() <= index)
            {
                x = 0;
                y = 0;
                return -1;
            }

            x = X.GetWeightPosition(index % X.GetWLen());
            y = Y.GetWeightPosition(index / X.GetWLen());
            return 0;
        }

        public void SetWeight(int index, double value)
        {
            if (index < Weights.GetLength())
            {
                Weights.GetBuffer()[index] = value;
            }
        }

        public void SetWeight(double[] values, int numberOfValues)
        {
            if (values == null)
            {
                return;
            }

            if (numberOfValues == Weights.GetLength())
            {
                Array.Copy(values, Weights.GetBuffer(), numberOfValues);
            }
        }

        public void ResetWeightSelection()
        {
            Array.Clear(WeightSelection, 0, WeightSelection.Length);
        }

        public int GetCarrierEdge(char axis, char edge, double[] buffer, int bufferLen)
        {
            return axis == 'x' ? X.GetCarrierEdge(edge, buffer, bufferLen) : Y.GetCarrierEdge(edge, buffer, bufferLen);
        }

        public void SetWeightSelection(double[] values, int numberOfValues)
        {
            if (values == null || numberOfValues <= 0)
            {
                return;
            }

            var count = Math.Min(numberOfValues, WeightSelection.Length);
            for (var i = 0; i < count; i++)
            {
                WeightSelection[i] = (int)Math.Round(values[i]);
            }
        }

        public int GetWeights(double[] buffer, int bufferLength)
        {
            if (bufferLength < Weights.GetLength())
            {
                return 0;
            }

            Array.Copy(Weights.GetBuffer(), buffer, Weights.GetLength());
            return Weights.GetLength();
        }

        public double[] GetWeightsBuffer()
        {
            return Weights.GetBuffer();
        }

        public int GetBinaryData<TIntervalX, TIntervalY, TOrderX, TOrderY>(XBSTools2Data<TIntervalX, TIntervalY, TOrderX, TOrderY> data)
            where TIntervalX : struct
            where TIntervalY : struct
            where TOrderX : struct
            where TOrderY : struct
        {
            if (data == null)
            {
                return 0;
            }

            data.SetCarrierX(X);
            data.SetCarrierY(Y);
            data.SetOrderOfLowerTaylorPolynomialX(LTpX);
            data.SetOrderOfUpperTaylorPolynomialX(UTpX);
            data.SetOrderOfLowerTaylorPolynomialY(LTpY);
            data.SetOrderOfUpperTaylorPolynomialY(UTpY);
            data.SetWeights(Weights.GetLength(), Weights.GetBuffer());
            return 1;
        }

        public int SetBinaryData<TIntervalX, TIntervalY, TOrderX, TOrderY>(XBSTools2Data<TIntervalX, TIntervalY, TOrderX, TOrderY> data)
            where TIntervalX : struct
            where TIntervalY : struct
            where TOrderX : struct
            where TOrderY : struct
        {
            if (data == null)
            {
                return 0;
            }

            data.GetCarrierX(X);
            data.GetCarrierY(Y);
            LTpX = data.GetOrderOfLowerTaylorPolynomialX();
            UTpX = data.GetOrderOfUpperTaylorPolynomialX();
            LTpY = data.GetOrderOfLowerTaylorPolynomialY();
            UTpY = data.GetOrderOfUpperTaylorPolynomialY();
            var weights = data.GetWeights();
            Array.Copy(weights, Weights.GetBuffer(), Math.Min(weights.Length, Weights.GetLength()));
            CompleteIt();
            return 1;
        }

        public BSTools GetCarrierY()
        {
            return Y;
        }

        public BSTools GetCarrierX()
        {
            return X;
        }

        public int GetOrderOfLowerTaylorPolynomialX()
        {
            return LTpX;
        }

        public int GetOrderOfLowerTaylorPolynomialY()
        {
            return LTpY;
        }

        public int GetOrderOfUpperTaylorPolynomialX()
        {
            return UTpX;
        }

        public int GetOrderOfUpperTaylorPolynomialY()
        {
            return UTpY;
        }

        public void SetOrderOfLowerTaylorPolynomialX(int order)
        {
            LTpX = order;
        }

        public void SetOrderOfUpperTaylorPolynomialX(int order)
        {
            UTpX = order;
        }

        public void SetOrderOfLowerTaylorPolynomialY(int order)
        {
            LTpY = order;
        }

        public void SetOrderOfUpperTaylorPolynomialY(int order)
        {
            UTpY = order;
        }

        public abstract double CalcWeightedNormOfBend(double weightsDirectionX = 1, double weightsDirectionY = 1, BSTools2 weightsPosition = null);

        public double CalcWeightedNormP(int p, double[] weights, int numberOfWeights)
        {
            var ret = 0.0;

            if (numberOfWeights == Weights.GetLength())
            {
                var d = Weights.GetBuffer();
                var sum = 0.0;

                for (var i = 0; i < numberOfWeights; i++)
                {
                    var value = d[i];
                    for (var j = 1; j < p; j++)
                    {
                        value *= d[i];
                    }

                    ret += value * weights[i];
                    sum += weights[i];
                }

                ret /= sum;
            }

            return ret;
        }

        public int SetCarrier(char nameOfCarrier, double a, double b, int n, int k, char lTp = 's', char uTp = 's')
        {
            if (nameOfCarrier == 'x')
            {
                X.InitInternal(a, b, n, k);
                LTpX = X.GetFGradInternal(lTp);
                UTpX = X.GetFGradInternal(uTp);
            }
            else if (nameOfCarrier == 'y')
            {
                Y.InitInternal(a, b, n, k);
                LTpY = Y.GetFGradInternal(lTp);
                UTpY = Y.GetFGradInternal(uTp);
            }

            return 0;
        }

        public int Init(BSTools x, BSTools y, double[] d = null, int numberOfWeights = 0, char lowerRestX = 's', char upperRestX = 's', char lowerRestY = 's', char upperRestY = 's')
        {
            X.InitInternal(x.LowerBoundary, x.UpperBoundary, x.IntervalCount, x.Order);
            Y.InitInternal(y.LowerBoundary, y.UpperBoundary, y.IntervalCount, y.Order);

            LTpX = X.GetFGradInternal(lowerRestX);
            UTpX = X.GetFGradInternal(upperRestX);
            LTpY = Y.GetFGradInternal(lowerRestY);
            UTpY = Y.GetFGradInternal(upperRestY);

            if (numberOfWeights != 0 && d != null)
            {
                var n = Y.GetWLen() * X.GetWLen();
                if (numberOfWeights == n)
                {
                    Array.Copy(d, Weights.GetBuffer(), n);
                }
            }

            return 0;
        }

        public void CompleteIt()
        {
            X.BuildEquidistantCarrierInternal();
            Y.BuildEquidistantCarrierInternal();
        }

        public bool IsStaticDataEqual(BSTools2 other)
        {
            if (X.IsNotEqual(other.X)) return false;
            if (Y.IsNotEqual(other.Y)) return false;
            if (Weights.IsNotEqual(other.Weights)) return false;
            if (LTpX != other.LTpX) return false;
            if (UTpX != other.UTpX) return false;
            if (LTpY != other.LTpY) return false;
            if (UTpY != other.UTpY) return false;
            return true;
        }

        public bool IsApproximatedEqual(BSTools2 other, double maxDelta)
        {
            if (!X.IsApproximatedEqual(other.X, maxDelta)) return false;
            if (!Y.IsApproximatedEqual(other.Y, maxDelta)) return false;
            if (!Weights.IsApproximatedEqual(other.Weights, maxDelta)) return false;
            if (maxDelta < Math.Abs(LTpX - other.LTpX)) return false;
            if (maxDelta < Math.Abs(UTpX - other.UTpX)) return false;
            if (maxDelta < Math.Abs(LTpY - other.LTpY)) return false;
            if (maxDelta < Math.Abs(UTpY - other.UTpY)) return false;
            return true;
        }

        private static (string Prefix, string Key) ExtractPrefixKey(string line)
        {
            if (string.IsNullOrWhiteSpace(line))
            {
                return (string.Empty, string.Empty);
            }

            var equalsIndex = line.IndexOf('=');
            if (equalsIndex <= 0)
            {
                return (string.Empty, string.Empty);
            }

            var left = line.Substring(0, equalsIndex).Trim();
            var parts = left.Split('.');
            if (parts.Length == 0)
            {
                return (string.Empty, string.Empty);
            }

            if (parts.Length == 1)
            {
                return (string.Empty, parts[0].Trim());
            }

            var prefix = parts[^2].Trim();
            var key = parts[^1].Trim();
            return (prefix, key);
        }

        private static string NormalizeFormat(string format)
        {
            if (string.IsNullOrWhiteSpace(format))
            {
                return "{0:0.################}";
            }

            if (format.Contains("{0", StringComparison.Ordinal))
            {
                return format;
            }

            if (!format.StartsWith("%", StringComparison.Ordinal))
            {
                return $"{{0:{format}}}";
            }

            var spec = format[^1];
            var dotIndex = format.IndexOf('.');
            var precision = -1;
            if (dotIndex >= 0 && dotIndex + 1 < format.Length - 1)
            {
                var span = format.AsSpan(dotIndex + 1, format.Length - dotIndex - 2);
                if (int.TryParse(span, out var parsed))
                {
                    precision = parsed;
                }
            }

            var netSpec = spec switch
            {
                'f' or 'F' => "F",
                'e' or 'E' => "E",
                'g' or 'G' => "G",
                _ => "G"
            };

            return precision >= 0 ? $"{{0:{netSpec}{precision}}}" : $"{{0:{netSpec}}}";
        }
    }
}
