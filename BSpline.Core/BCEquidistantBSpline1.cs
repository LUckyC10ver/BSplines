using System;

namespace BSpline.Core
{
    public sealed class BCEquidistantBSpline1
    {
        private readonly int _maxIntervals;
        private readonly int _maxOrder;
        private readonly XBSTools1<int, int> _bspline;

        public BCEquidistantBSpline1(int maxIntervals, int maxOrder)
        {
            _maxIntervals = maxIntervals;
            _maxOrder = maxOrder;
            _bspline = new XBSTools1<int, int>(maxIntervals, maxOrder);
            BCEquidistantBSpline.Assert(maxOrder <= BSToolsConstants.BSTOOLS_TOTAL_MAX_ORDER,
                "Maximal order of B-splines is too large.");
        }

        public string ClassName => "BCEquidistantBSpline1";

        public int GetMaxIntervals() => _maxIntervals;

        public int GetMaxOrder() => _maxOrder;

        public double GetLowerBoundary() => _bspline.LowerBound();

        public double GetUpperBoundary() => _bspline.UpperBound();

        public int GetNumberOfAllWeights() => _bspline.GetNumberOfAllWeights();

        public XBSTools1<int, int> GetXBSTools() => _bspline;

        public double GetWeightPosition(int index) => _bspline.GetWeightPosition(index);

        public int Save(string fileName, string splineName = "data", string client = "none")
        {
            return _bspline.Save(fileName, splineName, client);
        }

        public double Evaluate(double xValue)
        {
            return _bspline.Evaluate(xValue);
        }

        public void Evaluate(double[] xVector, int number, double[] fVector)
        {
            _bspline.Evaluate(xVector, number, fVector);
        }

        public void Derivative(double[] xVector, int number, double[] fVector)
        {
            _bspline.Derivative(xVector, number, fVector);
        }

        public void SecondDerivative(double[] xVector, int number, double[] fVector)
        {
            _bspline.SecondDerivative(xVector, number, fVector);
        }

        public double Inverse(double fValue)
        {
            return _bspline.Inverse(fValue);
        }

        public void GetBinaryData(XBSTools1Data<int, int> data)
        {
            _bspline.GetBinaryData(data);
        }

        public bool IsStaticDataEqual(BCEquidistantBSpline1 other)
        {
            return _bspline.IsStaticDataEqual(other.GetXBSTools());
        }

        public void SetCarrier(double lower, double upper, int numSub, int order, char lowerPoly = 's', char upperPoly = 's')
        {
            BCEquidistantBSpline.Assert(lower < upper, "Lower boundary is not less than upper one.");
            BCEquidistantBSpline.Assert(numSub > 0, "Number of subintervals is not greater than 0.");
            BCEquidistantBSpline.Assert(numSub <= _maxIntervals, "Number of subintervals of B-splines is oversized.");
            BCEquidistantBSpline.Assert(order > 0, "Order of B-splines is not greater than 0.");
            BCEquidistantBSpline.Assert(order <= _maxOrder, "Order of B-splines is oversized.");
            _bspline.SetCarrier(lower, upper, numSub, order, lowerPoly, upperPoly);
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

        public BCEquidistantBSpline1 Derivative()
        {
            _bspline.Diff();
            return this;
        }

        public BCEquidistantBSpline1 AntiDerivative(double x0 = 0, double f0 = 0)
        {
            _bspline.Integral(x0, f0);
            return this;
        }

        public void SetBinaryData(XBSTools1Data<int, int> data)
        {
            _bspline.SetBinaryData(data);
        }
    }
}
