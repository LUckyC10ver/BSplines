using System;

namespace BSpline.Core
{
    public sealed class BCEquidistantBSpline2
    {
        private readonly int _maxIntervalsX;
        private readonly int _maxIntervalsY;
        private readonly int _maxOrderX;
        private readonly int _maxOrderY;
        private readonly XBSTools2<int, int, int, int> _bspline;

        public BCEquidistantBSpline2(int maxIntervalsX, int maxIntervalsY, int maxOrderX, int maxOrderY)
        {
            _maxIntervalsX = maxIntervalsX;
            _maxIntervalsY = maxIntervalsY;
            _maxOrderX = maxOrderX;
            _maxOrderY = maxOrderY;
            _bspline = new XBSTools2<int, int, int, int>(maxIntervalsX, maxIntervalsY, maxOrderX, maxOrderY);
            BCEquidistantBSpline.Assert(maxOrderX <= BSToolsConstants.BSTOOLS_TOTAL_MAX_ORDER,
                "Maximal order of B-splines in x-direction too large.");
            BCEquidistantBSpline.Assert(maxOrderY <= BSToolsConstants.BSTOOLS_TOTAL_MAX_ORDER,
                "Maximal order of B-splines in y-direction too large.");
        }

        public string ClassName => "BCEquidistantBSpline2";

        public int GetMaxIntervalsX() => _maxIntervalsX;

        public int GetMaxIntervalsY() => _maxIntervalsY;

        public int GetMaxOrderX() => _maxOrderX;

        public int GetMaxOrderY() => _maxOrderY;

        public double GetLowerBoundaryX() => _bspline.LowerBoundX();

        public double GetUpperBoundaryX() => _bspline.UpperBoundX();

        public double GetLowerBoundaryY() => _bspline.LowerBoundY();

        public double GetUpperBoundaryY() => _bspline.UpperBoundY();

        public int GetNumberOfAllWeights() => _bspline.GetNumberOfAllWeights();

        public XBSTools2<int, int, int, int> GetXBSTools() => _bspline;

        public int GetWeightPosition(int index, out double xValue, out double yValue)
        {
            return _bspline.GetWeightPosition(index, out xValue, out yValue);
        }

        public int Save(string fileName, string splineName = "data", string client = "none")
        {
            return _bspline.Save(fileName, splineName, client);
        }

        public double Evaluate(double xValue, double yValue)
        {
            return _bspline.Evaluate(xValue, yValue);
        }

        public void Evaluate(double[] xVector, double[] yVector, int number, double[] fVector)
        {
            _bspline.Evaluate(xVector, yVector, number, fVector);
        }

        public void GetBinaryData(XBSTools2Data<int, int, int, int> data)
        {
            _bspline.GetBinaryData(data);
        }

        public bool IsStaticDataEqual(BCEquidistantBSpline2 other)
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

        public BCEquidistantBSpline2 DerivativeX()
        {
            _bspline.DiffX();
            return this;
        }

        public BCEquidistantBSpline2 DerivativeY()
        {
            _bspline.DiffY();
            return this;
        }

        public void SetBinaryData(XBSTools2Data<int, int, int, int> data)
        {
            _bspline.SetBinaryData(data);
        }
    }
}
