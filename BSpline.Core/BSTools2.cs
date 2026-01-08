using System;
using System.IO;

namespace BSpline.Core
{
public sealed class XBSTools2Data<TIntervalX, TIntervalY, TOrderX, TOrderY>
        where TIntervalX : struct
        where TIntervalY : struct
        where TOrderX : struct
        where TOrderY : struct
    {
        private readonly BSToolsData _x = new BSToolsData();
        private readonly BSToolsData _y = new BSToolsData();
        private int _lTpX;
        private int _uTpX;
        private int _lTpY;
        private int _uTpY;
        private int _numberOfW;
        private readonly double[] _w;

        public XBSTools2Data(int sizeOfW)
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

        public int GetCarrierX(BSTools carrier)
        {
            return carrier.SetBinaryData(_x);
        }

        public int GetCarrierY(BSTools carrier)
        {
            return carrier.SetBinaryData(_y);
        }

        public int SetCarrierX(BSTools carrier)
        {
            return carrier.GetBinaryData(_x);
        }

        public int SetCarrierY(BSTools carrier)
        {
            return carrier.GetBinaryData(_y);
        }

        public void Print(TextWriter writer, string path, string name)
        {
            throw new NotImplementedException();
        }

        public bool IsValid(int maxIntervalX, int maxOrderX, int maxIntervalY, int maxOrderY)
        {
            if (_numberOfW != _x.CalcWLen() * _y.CalcWLen())
            {
                return false;
            }

            for (var i = 0; i < _numberOfW; i++)
            {
                if (_w[i] < -1e50 || _w[i] > 1e50)
                {
                    return false;
                }
            }

            return _x.IsValid(maxIntervalX, maxOrderX) && _y.IsValid(maxIntervalY, maxOrderY);
        }
    }

    public sealed class XBSTools2<TIntervalX, TIntervalY, TOrderX, TOrderY> : BSTools2
        where TIntervalX : struct
        where TIntervalY : struct
        where TOrderX : struct
        where TOrderY : struct
    {
        public XBSTools2(
            int maxIntervalX,
            int maxIntervalY,
            int maxOrderX,
            int maxOrderY)
            : base(0, Array.Empty<double>(), new XBSTools<TIntervalX, TOrderX>(maxIntervalX, maxOrderX), new XBSTools<TIntervalY, TOrderY>(maxIntervalY, maxOrderY), 1, new[] { 0 }, 1, new[] { 0 })
        {
        }

        public override BSTools2 DiffX()
        {
            throw new NotImplementedException();
        }

        public override BSTools2 DiffY()
        {
            throw new NotImplementedException();
        }

        public override BSTools2 IntX(double x = 0, double z = 0)
        {
            throw new NotImplementedException();
        }

        public override BSTools2 IntY(double y = 0, double z = 0)
        {
            throw new NotImplementedException();
        }

        public override double CalcWeightedNormOfBend(double weightsDirectionX = 1, double weightsDirectionY = 1, BSTools2 weightsPosition = null)
        {
            throw new NotImplementedException();
        }
    }
}
