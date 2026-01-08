using System;

namespace BSpline.Core
{
    public sealed class BCEquidistantBSpline3
    {
        private readonly int _maxIntervalsX;
        private readonly int _maxIntervalsY;
        private readonly int _maxIntervalsZ;
        private readonly int _maxOrderX;
        private readonly int _maxOrderY;
        private readonly int _maxOrderZ;
        private readonly XBSTools3Pre2<int, int, int, int, int, int> _bspline;

        public BCEquidistantBSpline3(int maxIntervalsX, int maxIntervalsY, int maxIntervalsZ, int maxOrderX, int maxOrderY, int maxOrderZ)
        {
            _maxIntervalsX = maxIntervalsX;
            _maxIntervalsY = maxIntervalsY;
            _maxIntervalsZ = maxIntervalsZ;
            _maxOrderX = maxOrderX;
            _maxOrderY = maxOrderY;
            _maxOrderZ = maxOrderZ;
            _bspline = new XBSTools3Pre2<int, int, int, int, int, int>(maxIntervalsX, maxIntervalsY, maxIntervalsZ, maxOrderX, maxOrderY, maxOrderZ);
            BCEquidistantBSpline.Assert(maxOrderX <= BSToolsConstants.BSTOOLS_TOTAL_MAX_ORDER,
                "Maximal order of B-splines in x-direction too large.");
            BCEquidistantBSpline.Assert(maxOrderY <= BSToolsConstants.BSTOOLS_TOTAL_MAX_ORDER,
                "Maximal order of B-splines in y-direction too large.");
            BCEquidistantBSpline.Assert(maxOrderZ <= BSToolsConstants.BSTOOLS_TOTAL_MAX_ORDER,
                "Maximal order of B-splines in z-direction too large.");
        }

        public string ClassName => "BCEquidistantBSpline3";

        public int GetMaxIntervalsX() => _maxIntervalsX;

        public int GetMaxIntervalsY() => _maxIntervalsY;

        public int GetMaxIntervalsZ() => _maxIntervalsZ;

        public int GetMaxOrderX() => _maxOrderX;

        public int GetMaxOrderY() => _maxOrderY;

        public int GetMaxOrderZ() => _maxOrderZ;

        public double GetLowerBoundaryX() => _bspline.LowerBoundX();

        public double GetUpperBoundaryX() => _bspline.UpperBoundX();

        public double GetLowerBoundaryY() => _bspline.LowerBoundY();

        public double GetUpperBoundaryY() => _bspline.UpperBoundY();

        public double GetLowerBoundaryZ() => _bspline.LowerBoundZ();

        public double GetUpperBoundaryZ() => _bspline.UpperBoundZ();

        public int GetNumberOfAllWeights() => _bspline.GetNumberOfAllWeights();

        public XBSTools3Pre2<int, int, int, int, int, int> GetXBSTools() => _bspline;

        public int GetWeightPosition(int index, out double xValue, out double yValue, out double zValue)
        {
            return _bspline.GetWeightPosition(index, out xValue, out yValue, out zValue);
        }

        public int Save(string fileName, string splineName = "data", string client = "none")
        {
            return _bspline.Save(fileName, splineName, client);
        }

        public double Evaluate(double xValue, double yValue, double zValue)
        {
            return _bspline.Evaluate(xValue, yValue, zValue);
        }

        public double Evaluate(double xValue)
        {
            return _bspline.Evaluate(xValue);
        }

        public void Evaluate(double[] xVector, int number, double yValue, double zValue, double[] fVector)
        {
            _bspline.Evaluate(xVector, number, yValue, zValue, fVector);
        }

        public double InvX(double fValue, double yValue, double zValue)
        {
            return _bspline.InvX(fValue, yValue, zValue);
        }

        public void Evaluate(double[] xVector, int number, double[] fVector)
        {
            _bspline.Evaluate(xVector, number, fVector);
        }

        public void GetBinaryData(XBSTools3Data<int, int, int, int, int, int> data)
        {
            _bspline.GetBinaryData(data);
        }

        public bool IsStaticDataEqual(BCEquidistantBSpline3 other)
        {
            return _bspline.IsStaticDataEqual(other.GetXBSTools());
        }

        public void SetCarrier(char direction, double lower, double upper, int numSub, int order, char lowerPoly = 's', char upperPoly = 's')
        {
            _bspline.SetCarrier(direction, lower, upper, numSub, order, lowerPoly, upperPoly);
        }

        public void SetCarrierX(double lower, double upper, int numSub, int order, char lowerPoly = 's', char upperPoly = 's')
        {
            SetCarrier('x', lower, upper, numSub, order, lowerPoly, upperPoly);
        }

        public void SetCarrierY(double lower, double upper, int numSub, int order, char lowerPoly = 's', char upperPoly = 's')
        {
            SetCarrier('y', lower, upper, numSub, order, lowerPoly, upperPoly);
        }

        public void SetCarrierZ(double lower, double upper, int numSub, int order, char lowerPoly = 's', char upperPoly = 's')
        {
            SetCarrier('z', lower, upper, numSub, order, lowerPoly, upperPoly);
        }

        public void SetWeight(int index, double weight)
        {
            _bspline.SetWeight(index, weight);
        }

        public void CompleteIt()
        {
            _bspline.CompleteIt();
        }

        public int Load(string filename)
        {
            return _bspline.Load(filename);
        }

        public BCEquidistantBSpline3 DerivativeX()
        {
            _bspline.DiffX();
            return this;
        }

        public void SetBinaryData(XBSTools3Data<int, int, int, int, int, int> data)
        {
            _bspline.SetBinaryData(data);
        }
    }
}
