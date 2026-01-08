using System;
using System.IO;

namespace BSpline.Core
{
public sealed class XBSTools3Data<TIntervalX, TIntervalY, TIntervalZ, TOrderX, TOrderY, TOrderZ>
        where TIntervalX : struct
        where TIntervalY : struct
        where TIntervalZ : struct
        where TOrderX : struct
        where TOrderY : struct
        where TOrderZ : struct
    {
        private readonly BSToolsData _x = new BSToolsData();
        private readonly BSToolsData _y = new BSToolsData();
        private readonly BSToolsData _z = new BSToolsData();
        private int _lTpX;
        private int _uTpX;
        private int _lTpY;
        private int _uTpY;
        private int _lTpZ;
        private int _uTpZ;
        private int _numberOfW;
        private readonly double[] _w;

        public XBSTools3Data(int sizeOfW)
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

        public int GetOrderOfLowerTaylorPolynomialX() => _lTpX;

        public void SetOrderOfLowerTaylorPolynomialX(int order)
        {
            _lTpX = order;
        }

        public int GetOrderOfUpperTaylorPolynomialX() => _uTpX;

        public void SetOrderOfUpperTaylorPolynomialX(int order)
        {
            _uTpX = order;
        }

        public int GetOrderOfLowerTaylorPolynomialY() => _lTpY;

        public void SetOrderOfLowerTaylorPolynomialY(int order)
        {
            _lTpY = order;
        }

        public int GetOrderOfUpperTaylorPolynomialY() => _uTpY;

        public void SetOrderOfUpperTaylorPolynomialY(int order)
        {
            _uTpY = order;
        }

        public int GetOrderOfLowerTaylorPolynomialZ() => _lTpZ;

        public void SetOrderOfLowerTaylorPolynomialZ(int order)
        {
            _lTpZ = order;
        }

        public int GetOrderOfUpperTaylorPolynomialZ() => _uTpZ;

        public void SetOrderOfUpperTaylorPolynomialZ(int order)
        {
            _uTpZ = order;
        }

        public int GetCarrierX(BSTools carrier)
        {
            return carrier.SetBinaryData(_x);
        }

        public int SetCarrierX(BSTools carrier)
        {
            return carrier.GetBinaryData(_x);
        }

        public int GetCarrierY(BSTools carrier)
        {
            return carrier.SetBinaryData(_y);
        }

        public int SetCarrierY(BSTools carrier)
        {
            return carrier.GetBinaryData(_y);
        }

        public int GetCarrierZ(BSTools carrier)
        {
            return carrier.SetBinaryData(_z);
        }

        public int SetCarrierZ(BSTools carrier)
        {
            return carrier.GetBinaryData(_z);
        }

        public void Print(TextWriter writer, string path, string name)
        {
            throw new NotImplementedException();
        }
    }

    public abstract class BSTools3Pre2
    {
        protected sealed class ArrayBuffer
        {
            private readonly double[] _data;
            private BSTools3Pre2 _owner;

            public ArrayBuffer(int sizeOfD, double[] data)
            {
                _data = data;
            }

            public void SetBasicPointer(BSTools3Pre2 owner)
            {
                _owner = owner;
            }

            public double GetValue(int ix, int iy, int iz)
            {
                return _data[GetFlatIndex(ix, iy, iz)];
            }

            public void SetValue(int ix, int iy, int iz, double value)
            {
                _data[GetFlatIndex(ix, iy, iz)] = value;
            }

            public int GetXIndex(int absIndex)
            {
                return (absIndex / (_owner.Y.GetWLen() * _owner.Z.GetWLen())) % _owner.X.GetWLen();
            }

            public int GetYIndex(int absIndex)
            {
                return (absIndex / _owner.Z.GetWLen()) % _owner.Y.GetWLen();
            }

            public int GetZIndex(int absIndex)
            {
                return absIndex % _owner.Z.GetWLen();
            }

            public bool IsApproximatedEqual(ArrayBuffer other, double maxDelta)
            {
                return BSTools.IsApproximatedEqual(_data, GetLength(), other._data, other.GetLength(), maxDelta);
            }

            public bool IsNotEqual(ArrayBuffer other)
            {
                return !BSTools.IsEqual(_data, GetLength(), other._data, other.GetLength());
            }

            public int GetLength()
            {
                return _data.Length;
            }

            public double[] GetBuffer()
            {
                return _data;
            }

            public void SetWeight(int index, double value)
            {
                if (index >= 0 && index < _data.Length)
                {
                    _data[index] = value;
                }
            }

            public double GetWeight(int index)
            {
                return _data[index];
            }

            private int GetFlatIndex(int ix, int iy, int iz)
            {
                return (iy + _owner.Y.GetWLen() * ix) * _owner.Z.GetWLen() + iz;
            }
        }

        protected readonly ArrayBuffer Weights;
        protected readonly BSTools X;
        protected readonly BSTools Y;
        protected readonly BSTools Z;
        protected readonly int SizeOfWx;
        protected readonly double[] Wx;
        protected readonly byte[] WxState;
        protected double Cy;
        protected double Cz;
        protected int LTpX;
        protected int UTpX;
        protected int LTpY;
        protected int UTpY;
        protected int LTpZ;
        protected int UTpZ;

        protected BSTools3Pre2(
            int sizeOfD,
            double[] d,
            BSTools x,
            BSTools y,
            BSTools z,
            int sizeOfWx,
            double[] wxMem,
            byte[] wxStateMem)
        {
            Weights = new ArrayBuffer(sizeOfD, d);
            X = x;
            Y = y;
            Z = z;
            SizeOfWx = sizeOfWx;
            Wx = wxMem;
            WxState = wxStateMem;
            Cy = 0;
            Cz = 0;
            ResetWx();
            Weights.SetBasicPointer(this);
        }

        private void BuildWX(int ix)
        {
            var ry = 0;
            var rz = 0;
            var wy = new double[BSToolsConstants.BSTOOLS_TOTAL_NUMBER_DEPENDENT];
            var tp = new double[BSToolsConstants.BSTOOLS_TOTAL_SIZE_TAYLOR];
            int tpN;

            var posY = Y.GetLocalInternal(Cy);
            var posZ = Z.GetLocalInternal(Cz);

            switch (posY)
            {
                case BSTools.Position.Mitte:
                    ry = Y.GetMainIndexInternal(Y.GetIntervalIndexInternal(Cy));
                    break;
                case BSTools.Position.Unter:
                    ry = Y.GetMainIndexInternal(Y.GetIntervalIndexInternal(Y.LowerBoundary));
                    break;
                case BSTools.Position.Ober:
                    ry = Y.GetMainIndexInternal(Y.GetIntervalIndexInternal(Y.UpperBoundary));
                    break;
            }

            switch (posZ)
            {
                case BSTools.Position.Mitte:
                    rz = Z.GetMainIndexInternal(Z.GetIntervalIndexInternal(Cz));
                    break;
                case BSTools.Position.Unter:
                    rz = Z.GetMainIndexInternal(Z.GetIntervalIndexInternal(Z.LowerBoundary));
                    break;
                case BSTools.Position.Ober:
                    rz = Z.GetMainIndexInternal(Z.GetIntervalIndexInternal(Z.UpperBoundary));
                    break;
            }

            for (var i = 0; i < X.Order; i++)
            {
                var rx = ix + i;
                if (X.GetWLen() == rx)
                {
                    break;
                }

                if (WxState[rx] == 0)
                {
                    switch (posZ)
                    {
                        case BSTools.Position.Mitte:
                            switch (Z.Order)
                            {
                                case 2:
                                    wy[0] = Z.GetFctValueK2Internal(Cz, GetWeightsSlice(ix + i, ry + 0, rz), Slice(Z.Carrier, rz, Z.Carrier.Length - rz));
                                    wy[1] = Z.GetFctValueK2Internal(Cz, GetWeightsSlice(ix + i, ry + 1, rz), Slice(Z.Carrier, rz, Z.Carrier.Length - rz));
                                    break;
                                case 3:
                                    wy[0] = Z.GetFctValueK3Internal(Cz, GetWeightsSlice(ix + i, ry + 0, rz), Slice(Z.Carrier, rz, Z.Carrier.Length - rz));
                                    wy[1] = Z.GetFctValueK3Internal(Cz, GetWeightsSlice(ix + i, ry + 1, rz), Slice(Z.Carrier, rz, Z.Carrier.Length - rz));
                                    wy[2] = Z.GetFctValueK3Internal(Cz, GetWeightsSlice(ix + i, ry + 2, rz), Slice(Z.Carrier, rz, Z.Carrier.Length - rz));
                                    break;
                                default:
                                    for (var j = ry; j < ry + Y.Order; j++)
                                    {
                                        wy[j - ry] = Z.GetFctValueKnInternal(Cz, Z.Order, GetWeightsSlice(ix + i, j, rz), Slice(Z.Carrier, rz, Z.Carrier.Length - rz));
                                    }
                                    break;
                            }
                            break;
                        case BSTools.Position.Unter:
                            for (var j = ry; j < ry + Y.Order; j++)
                            {
                                tpN = LTpZ;
                                Z.TaylorPolyInternal(GetWeightsSlice(ix + i, j, rz), Z.LowerBoundary, tpN, tp);
                                wy[j - ry] = Z.PolyValInternal(tp, tpN, Cz - Z.LowerBoundary);
                            }
                            break;
                        case BSTools.Position.Ober:
                            for (var j = ry; j < ry + Y.Order; j++)
                            {
                                tpN = UTpZ;
                                Z.TaylorPolyInternal(GetWeightsSlice(ix + i, j, rz), Z.UpperBoundary, tpN, tp);
                                wy[j - ry] = Z.PolyValInternal(tp, tpN, Cz - Z.UpperBoundary);
                            }
                            break;
                    }

                    switch (posY)
                    {
                        case BSTools.Position.Mitte:
                            Wx[ix + i] = Y.GetFctValueInternal(Cy, Y.Order, wy, Slice(Y.Carrier, ry, Y.Carrier.Length - ry));
                            break;
                        case BSTools.Position.Unter:
                            tpN = LTpY;
                            Y.TaylorPolyInternal(wy, Y.LowerBoundary, tpN, tp);
                            Wx[ix + i] = Y.PolyValInternal(tp, tpN, Cy - Y.LowerBoundary);
                            break;
                        case BSTools.Position.Ober:
                            tpN = UTpY;
                            Y.TaylorPolyInternal(wy, Y.UpperBoundary, tpN, tp);
                            Wx[ix + i] = Y.PolyValInternal(tp, tpN, Cy - Y.UpperBoundary);
                            break;
                    }
                }

                WxState[rx] = (byte)X.MaxInternal(X.Order - i, WxState[rx]);
            }
        }

        public abstract BSTools3Pre2 DiffX();

        protected BSTools3Pre2 AddSet(BSTools3Pre2 other, BSTools3Pre2 temp)
        {
            temp = this;
            ResetWx();

            SetMaxCarrier('x', this, other);
            SetMaxCarrier('y', this, other);
            SetMaxCarrier('z', this, other);

            var n = GetNumberOfAllWeights();
            for (var i = 0; i < n; i++)
            {
                GetWeightPosition(i, out var x, out var y, out var z);
                SetWeight(i, temp.Evaluate(x, y, z) + other.Evaluate(x, y, z));
            }

            return this;
        }

        protected BSTools3Pre2 MultWeights(BSTools3Pre2 other)
        {
            var n = GetNumberOfAllWeights();
            var dest = Weights.GetBuffer();
            var src = other.Weights.GetBuffer();
            for (var i = 0; i < n; i++)
            {
                dest[i] *= src[i];
            }

            ResetWx();
            return this;
        }

        protected BSTools3Pre2 MultFctValues(BSTools3Pre2 other, BSTools3Pre2 temp)
        {
            temp = this;

            SetMaxCarrier('x', this, other);
            SetMaxCarrier('y', this, other);
            SetMaxCarrier('z', this, other);
            CompleteIt();

            var n = GetNumberOfAllWeights();
            for (var i = 0; i < n; i++)
            {
                GetWeightPosition(i, out var x, out var y, out var z);
                SetWeight(i, temp.Evaluate(x, y, z) * other.Evaluate(x, y, z));
            }

            return this;
        }

        public int GetWeights(int bufferLen, double[] buffer)
        {
            var n = bufferLen;
            if (GetNumberOfAllWeights() < n)
            {
                n = GetNumberOfAllWeights();
            }

            var p = Weights.GetBuffer();
            Array.Copy(p, buffer, n);
            return n;
        }

        public void CompleteIt()
        {
            X.BuildEquidistantCarrierInternal();
            Y.BuildEquidistantCarrierInternal();
            Z.BuildEquidistantCarrierInternal();
            Cy = 0;
            Cz = 0;
            ResetWx();
        }

        public int Load(string filename)
        {
            throw new NotImplementedException();
        }

        public int Save(string filename, string splineName = "data", string comment = "non")
        {
            throw new NotImplementedException();
        }

        public void Preset(double y, double z)
        {
            if (Cy != y || Cz != z)
            {
                ResetWx();
                Cy = y;
                Cz = z;
            }
        }

        public double Evaluate(double x)
        {
            var result = new double[1];
            Evaluate(new[] { x }, 1, result);
            return result[0];
        }

        public double Evaluate(double x, double y, double z)
        {
            Preset(y, z);
            return Evaluate(x);
        }

        public void Evaluate(double[] x, int dim, double[] result)
        {
            switch (X.Order)
            {
                case 0:
                    return;
                case 2:
                    for (var i = 0; i < dim; i++)
                    {
                        var ex = x[i];
                        if (ex < X.LowerBoundary)
                        {
                            var tp = new double[BSToolsConstants.BSTOOLS_TOTAL_SIZE_TAYLOR];
                            var r = X.GetMainIndexK2Internal(X.GetIntervalIndexInternal(X.LowerBoundary));
                            X.TaylorPolyInternal(GetWX(r), X.LowerBoundary, LTpX, tp);
                            result[i] = X.PolyValInternal(tp, LTpX, ex - X.LowerBoundary);
                        }
                        else if (ex <= X.UpperBoundary)
                        {
                            var r = X.GetMainIndexK2Internal(X.GetIntervalIndexInternal(ex));
                            result[i] = X.GetFctValueK2Internal(ex, GetWX(r), Slice(X.Carrier, r, X.Carrier.Length - r));
                        }
                        else
                        {
                            var tp = new double[BSToolsConstants.BSTOOLS_TOTAL_SIZE_TAYLOR];
                            var r = X.GetMainIndexK2Internal(X.GetIntervalIndexInternal(X.UpperBoundary));
                            X.TaylorPolyInternal(GetWX(r), X.UpperBoundary, UTpX, tp);
                            result[i] = X.PolyValInternal(tp, UTpX, ex - X.UpperBoundary);
                        }
                    }
                    break;
                case 3:
                    for (var i = 0; i < dim; i++)
                    {
                        var ex = x[i];
                        if (ex < X.LowerBoundary)
                        {
                            var tp = new double[BSToolsConstants.BSTOOLS_TOTAL_SIZE_TAYLOR];
                            var r = X.GetMainIndexK3Internal(X.GetIntervalIndexInternal(X.LowerBoundary));
                            X.TaylorPolyInternal(GetWX(r), X.LowerBoundary, LTpX, tp);
                            result[i] = X.PolyValInternal(tp, LTpX, ex - X.LowerBoundary);
                        }
                        else if (ex <= X.UpperBoundary)
                        {
                            var r = X.GetMainIndexK3Internal(X.GetIntervalIndexInternal(ex));
                            result[i] = X.GetFctValueK3Internal(ex, GetWX(r), Slice(X.Carrier, r, X.Carrier.Length - r));
                        }
                        else
                        {
                            var tp = new double[BSToolsConstants.BSTOOLS_TOTAL_SIZE_TAYLOR];
                            var r = X.GetMainIndexK3Internal(X.GetIntervalIndexInternal(X.UpperBoundary));
                            X.TaylorPolyInternal(GetWX(r), X.UpperBoundary, UTpX, tp);
                            result[i] = X.PolyValInternal(tp, UTpX, ex - X.UpperBoundary);
                        }
                    }
                    break;
                default:
                    for (var i = 0; i < dim; i++)
                    {
                        var ex = x[i];
                        if (X.LowerBoundary <= ex && ex <= X.UpperBoundary)
                        {
                            var r = X.GetMainIndexInternal(X.GetIntervalIndexInternal(ex));
                            result[i] = X.GetFctValueKnInternal(ex, X.Order, GetWX(r), Slice(X.Carrier, r, X.Carrier.Length - r));
                        }
                        else if (ex < X.LowerBoundary)
                        {
                            var tp = new double[BSToolsConstants.BSTOOLS_TOTAL_SIZE_TAYLOR];
                            var r = X.GetMainIndexInternal(X.GetIntervalIndexInternal(X.LowerBoundary));
                            X.TaylorPolyInternal(GetWX(r), X.LowerBoundary, LTpX, tp);
                            result[i] = X.PolyValInternal(tp, LTpX, ex - X.LowerBoundary);
                        }
                        else
                        {
                            var tp = new double[BSToolsConstants.BSTOOLS_TOTAL_SIZE_TAYLOR];
                            var r = X.GetMainIndexInternal(X.GetIntervalIndexInternal(X.UpperBoundary));
                            X.TaylorPolyInternal(GetWX(r), X.UpperBoundary, UTpX, tp);
                            result[i] = X.PolyValInternal(tp, UTpX, ex - X.UpperBoundary);
                        }
                    }
                    break;
            }
        }

        public void Evaluate(double[] x, int dim, double y, double z, double[] result)
        {
            Preset(y, z);
            Evaluate(x, dim, result);
        }

        public double InvX(double r, double y, double z)
        {
            Preset(y, z);

            var x0 = X.LowerBoundary;
            var x1 = X.UpperBoundary;
            var f0 = Evaluate(x0) - r;
            var f1 = Evaluate(x1) - r;

            if (f0 * f1 < 0)
            {
                var rx = x0;
                var rf = f0;

                for (var i = 0; i < 100; i++)
                {
                    if (Math.Abs(f1) < BSToolsConstants.Delta)
                    {
                        break;
                    }

                    var df = f1 - f0;
                    if (Math.Abs(df) > 0)
                    {
                        var x = x1 - f1 * (x1 - x0) / df;
                        x0 = x1;
                        f0 = f1;
                        x1 = x;
                        f1 = Evaluate(x1) - r;

                        if (rf * f1 < 0)
                        {
                            rx = x0;
                            rf = f0;
                        }
                    }
                    else
                    {
                        var xn = x1 - f1 * (x1 - rx) / (f1 - rf);
                        var fn = Evaluate(xn) - r;

                        if (rf * fn < 0)
                        {
                            x1 = xn;
                            f1 = fn;
                        }
                        else
                        {
                            rx = xn;
                            rf = fn;
                            x0 = xn;
                            f0 = fn;
                        }
                    }
                }
            }

            return x1;
        }

        public BSTools3Pre2 Multiply(double value)
        {
            var n = GetNumberOfAllWeights();
            var weights = Weights.GetBuffer();
            for (var i = 0; i < n; i++)
            {
                weights[i] *= value;
            }

            return this;
        }

        public BSTools3Pre2 Multiply(BSTools3Pre2 other)
        {
            if (!X.IsNotEqual(other.X) && !Y.IsNotEqual(other.Y) && !Z.IsNotEqual(other.Z))
            {
                return MultWeights(other);
            }

            return MultFctValues(other, this);
        }

        public BSTools3Pre2 Add(double value)
        {
            var n = GetNumberOfAllWeights();
            var weights = Weights.GetBuffer();
            for (var i = 0; i < n; i++)
            {
                weights[i] += value;
            }

            return this;
        }

        public abstract BSTools3Pre2 Add(BSTools3Pre2 other);

        public bool IsApproximatedEqual(BSTools3Pre2 other, double maxDelta)
        {
            if (!X.IsApproximatedEqual(other.X, maxDelta)) return false;
            if (!Y.IsApproximatedEqual(other.Y, maxDelta)) return false;
            if (!Z.IsApproximatedEqual(other.Z, maxDelta)) return false;
            if (!Weights.IsApproximatedEqual(other.Weights, maxDelta)) return false;
            if (maxDelta < Math.Abs(LTpX - other.LTpX)) return false;
            if (maxDelta < Math.Abs(UTpX - other.UTpX)) return false;
            if (maxDelta < Math.Abs(LTpY - other.LTpY)) return false;
            if (maxDelta < Math.Abs(UTpY - other.UTpY)) return false;
            if (maxDelta < Math.Abs(LTpZ - other.LTpZ)) return false;
            if (maxDelta < Math.Abs(UTpZ - other.UTpZ)) return false;
            return true;
        }

        public bool IsStaticDataEqual(BSTools3Pre2 other)
        {
            if (X.IsNotEqual(other.X)) return false;
            if (Y.IsNotEqual(other.Y)) return false;
            if (Z.IsNotEqual(other.Z)) return false;
            if (Weights.IsNotEqual(other.Weights)) return false;
            if (LTpX != other.LTpX) return false;
            if (UTpX != other.UTpX) return false;
            if (LTpY != other.LTpY) return false;
            if (UTpY != other.UTpY) return false;
            if (LTpZ != other.LTpZ) return false;
            if (UTpZ != other.UTpZ) return false;
            return true;
        }

        public int GetOrderOfLowerTaylorPolynomialX()
        {
            return LTpX;
        }

        public int GetOrderOfLowerTaylorPolynomialY()
        {
            return LTpY;
        }

        public int GetOrderOfLowerTaylorPolynomialZ()
        {
            return LTpZ;
        }

        public int GetOrderOfUpperTaylorPolynomialX()
        {
            return UTpX;
        }

        public int GetOrderOfUpperTaylorPolynomialY()
        {
            return UTpY;
        }

        public int GetOrderOfUpperTaylorPolynomialZ()
        {
            return UTpZ;
        }

        public void SetOrderOfLowerTaylorPolynomialX(int order)
        {
            LTpX = order;
        }

        public void SetOrderOfLowerTaylorPolynomialY(int order)
        {
            LTpY = order;
        }

        public void SetOrderOfLowerTaylorPolynomialZ(int order)
        {
            LTpZ = order;
        }

        public void SetOrderOfUpperTaylorPolynomialX(int order)
        {
            UTpX = order;
        }

        public void SetOrderOfUpperTaylorPolynomialY(int order)
        {
            UTpY = order;
        }

        public void SetOrderOfUpperTaylorPolynomialZ(int order)
        {
            UTpZ = order;
        }

        public int GetNumberOfWeights()
        {
            return Weights.GetLength();
        }

        public int GetNumberOfAllWeights()
        {
            return Weights.GetLength();
        }

        public int GetWeightPosition(int index, out double x, out double y, out double z)
        {
            if (Weights.GetLength() <= index)
            {
                x = 0;
                y = 0;
                z = 0;
                return -1;
            }

            x = X.GetWeightPosition(Weights.GetXIndex(index));
            y = Y.GetWeightPosition(Weights.GetYIndex(index));
            z = Z.GetWeightPosition(Weights.GetZIndex(index));
            return 0;
        }

        public void SetWeight(int index, double value)
        {
            Weights.SetWeight(index, value);
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

            ResetWx();
        }

        public double GetWeight(int index)
        {
            return Weights.GetWeight(index);
        }

        private void ResetWx()
        {
            Array.Clear(Wx, 0, Wx.Length);
            Array.Clear(WxState, 0, WxState.Length);
        }

        public void SetCarrier(char nameOfCarrier, double a, double b, int n, int k, char lTp = 's', char uTp = 's')
        {
            switch (nameOfCarrier)
            {
                case 'x':
                    X.InitInternal(a, b, n, k);
                    LTpX = X.GetFGradInternal(lTp);
                    UTpX = X.GetFGradInternal(uTp);
                    break;
                case 'y':
                    Y.InitInternal(a, b, n, k);
                    LTpY = Y.GetFGradInternal(lTp);
                    UTpY = Y.GetFGradInternal(uTp);
                    break;
                case 'z':
                    Z.InitInternal(a, b, n, k);
                    LTpZ = Z.GetFGradInternal(lTp);
                    UTpZ = Z.GetFGradInternal(uTp);
                    break;
            }
        }

        public void SetMaxCarrier(char nameOfCarrier, BSTools3Pre2 bs1, BSTools3Pre2 bs2)
        {
            switch (nameOfCarrier)
            {
                case 'x':
                    X.InitInternal(bs1.X.LowerBoundary, bs1.X.UpperBoundary, bs1.X.IntervalCount, bs1.X.Order);
                    X.MaxInternal(bs2.X);
                    LTpX = X.MaxInternal(bs1.LTpX, bs2.LTpX);
                    UTpX = X.MaxInternal(bs1.UTpX, bs2.UTpX);
                    break;
                case 'y':
                    Y.InitInternal(bs1.Y.LowerBoundary, bs1.Y.UpperBoundary, bs1.Y.IntervalCount, bs1.Y.Order);
                    Y.MaxInternal(bs2.Y);
                    LTpY = Y.MaxInternal(bs1.LTpY, bs2.LTpY);
                    UTpY = Y.MaxInternal(bs1.UTpY, bs2.UTpY);
                    break;
                case 'z':
                    Z.InitInternal(bs1.Z.LowerBoundary, bs1.Z.UpperBoundary, bs1.Z.IntervalCount, bs1.Z.Order);
                    Z.MaxInternal(bs2.Z);
                    LTpZ = Z.MaxInternal(bs1.LTpZ, bs2.LTpZ);
                    UTpZ = Z.MaxInternal(bs1.UTpZ, bs2.UTpZ);
                    break;
            }
        }

        private double[] GetWX(int index)
        {
            if (WxState[index] != X.Order)
            {
                BuildWX(index);
            }

            return Slice(Wx, index, Wx.Length - index);
        }

        private double[] GetWeightsSlice(int ix, int iy, int iz)
        {
            var buffer = new double[Z.Order];
            for (var i = 0; i < buffer.Length; i++)
            {
                buffer[i] = Weights.GetValue(ix, iy, iz + i);
            }

            return buffer;
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

    public sealed class XBSTools3Pre2<TIntervalX, TIntervalY, TIntervalZ, TOrderX, TOrderY, TOrderZ> : BSTools3Pre2
        where TIntervalX : struct
        where TIntervalY : struct
        where TIntervalZ : struct
        where TOrderX : struct
        where TOrderY : struct
        where TOrderZ : struct
    {
        public XBSTools3Pre2(int maxIntervalX, int maxIntervalY, int maxIntervalZ, int maxOrderX, int maxOrderY, int maxOrderZ)
            : base(0, Array.Empty<double>(),
                new XBSTools<TIntervalX, TOrderX>(maxIntervalX, maxOrderX),
                new XBSTools<TIntervalY, TOrderY>(maxIntervalY, maxOrderY),
                new XBSTools<TIntervalZ, TOrderZ>(maxIntervalZ, maxOrderZ),
                0, Array.Empty<double>(), Array.Empty<byte>())
        {
        }

        public override BSTools3Pre2 DiffX()
        {
            throw new NotImplementedException();
        }

        public override BSTools3Pre2 Add(BSTools3Pre2 other)
        {
            throw new NotImplementedException();
        }
    }
}
}
