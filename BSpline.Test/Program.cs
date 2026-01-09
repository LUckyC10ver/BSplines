using System;
using BSpline.Core;

namespace BSpline.Test
{
    internal static class Program
    {
        private const int MaxIntervalsX = 40;
        private const int MaxOrderX = 6;
        private const int MaxIntervalsY = 20;
        private const int MaxOrderY = 5;
        private const int MaxIntervalsZ = 10;
        private const int MaxOrderZ = 4;
        private const double Epsilon = 1e-9;

        private static int Main(string[] args)
        {
            try
            {
                Console.WriteLine("--- Create B-spline functions.");
                var bcbspline1 = new EquidistantBSpline1(MaxIntervalsX, MaxOrderX);
                var bspline1 = new XBSTools1<int, int>(MaxIntervalsX, MaxOrderX);
                var bcbspline3 = new EquidistantBSpline3(MaxIntervalsZ, MaxIntervalsZ, MaxIntervalsZ, MaxOrderX, MaxOrderY, MaxOrderZ);
                var bspline3 = new XBSTools3Pre2<int, int, int, int, int, int>(MaxIntervalsZ, MaxIntervalsZ, MaxIntervalsZ, MaxOrderX, MaxOrderY, MaxOrderZ);

                Console.WriteLine("\n------- Fill B-spline functions.");
                BuildBSpline1(bcbspline1, bspline1);
                BuildBSpline3(bcbspline3, bspline3);

                Console.WriteLine("\n------- Examples of evaluation of B-spline functions.");
                EvalBSpline1(bcbspline1, bspline1);
                EvalBSpline2();
                EvalBSpline3(bcbspline3, bspline3);

                Console.WriteLine("\n------- Examples for transformations of B-spline functions.");
                TrafoBSpline1(bcbspline1, bspline1);
                TrafoBSpline2();
                TrafoBSpline3();

                Console.WriteLine("\n------- Basic 1-dimensional B-splines for editing.");
                BasicBSplines();

                Console.WriteLine("\n------- Test time for old and new interfaces.");
                Console.WriteLine("Time tests are not implemented in the managed port.");

                return 0;
            }
            catch (Exception ex)
            {
                Console.Error.WriteLine($"Error: {ex.Message}");
                return 1;
            }
        }

        private static double F1(double x)
        {
            return Math.Sin(x * 3.14);
        }

        private static double F2(double x, double y)
        {
            var sum = x + y;
            return sum * sum;
        }

        private static void BuildBSpline1(EquidistantBSpline1 bcbspline, XBSTools1<int, int> bspline)
        {
            Console.WriteLine("--- Initialize 1-dimensional B-spline functions:");
            Console.WriteLine("    x-working area = [{0}, {1}], {2} subintervals, order = {3}.", -1.0, 2.0, 30, 4);
            bcbspline.SetCarrier(-1.0, 2.0, 30, 4);
            bspline.SetCarrier(-1.0, 2.0, 30, 4);

            var bcCount = bcbspline.GetNumberOfAllWeights();
            var count = bspline.GetNumberOfAllWeights();
            AssertEqual(bcCount, count, "Different number of weights.");
            Console.WriteLine("--- Fill 1-dimensional B-spline functions with {0} weights.", count);
            for (var i = 0; i < count; i++)
            {
                var x = bcbspline.GetWeightPosition(i);
                var r = bspline.GetWeightPosition(i);
                AssertClose(x, r, "Different weight positions.");
                var w = F1(x);
                if (i < 10 || count - 10 <= i)
                {
                    Console.WriteLine("    Set weight {0} ({1}) to {2}.", i, x, w);
                }

                bcbspline.SetWeight(i, w);
                bspline.SetWeight(i, w);
            }

            Console.WriteLine("--- Complete initialization of 1-dimensional B-spline.");
            bcbspline.CompleteIt();
            bspline.CompleteIt();
            Console.WriteLine("--- Save is not implemented in the managed port; skipping.");
        }

        private static void BuildBSpline3(EquidistantBSpline3 bcbspline, XBSTools3Pre2<int, int, int, int, int, int> bspline)
        {
            Console.WriteLine("--- Initialize 3-dimensional B-spline functions:");

            Console.WriteLine("    {0}-working area = [{1}, {2}], {3} subintervals, order = {4}.", 'x', -1.0, 2.0, 6, 4);
            bcbspline.SetCarrierX(-1.0, 2.0, 6, 4);
            bspline.SetCarrier('x', -1.0, 2.0, 6, 4);

            Console.WriteLine("    {0}-working area = [{1}, {2}], {3} subintervals, order = {4}.", 'y', -1.5, 1.0, 8, 4);
            bcbspline.SetCarrierY(-1.5, 1.0, 8, 4);
            bspline.SetCarrier('y', -1.5, 1.0, 8, 4);

            Console.WriteLine("    {0}-working area = [{1}, {2}], {3} subintervals, order = {4}.", 'z', 0.0, 3.0, 10, 4);
            bcbspline.SetCarrierZ(0.0, 3.0, 10, 4);
            bspline.SetCarrier('z', 0.0, 3.0, 10, 4);

            var bcCount = bcbspline.GetNumberOfAllWeights();
            var count = bspline.GetNumberOfAllWeights();
            AssertEqual(bcCount, count, "Different number of weights.");
            Console.WriteLine("--- Fill 3-dimensional B-spline functions with {0} weights.", count);
            for (var i = 0; i < count; i++)
            {
                bcbspline.GetWeightPosition(i, out var x, out var y, out var z);
                bspline.GetWeightPosition(i, out var r, out var s, out var t);
                AssertClose(x, r, "Different weight positions.");
                AssertClose(y, s, "Different weight positions.");
                AssertClose(z, t, "Different weight positions.");
                var w = F2(x, y);
                if (i < 10 || count - 10 <= i)
                {
                    Console.WriteLine("    Set weight {0} ({1}, {2}, {3}) to {4}.", i, x, y, z, w);
                }

                bcbspline.SetWeight(i, w);
                bspline.SetWeight(i, w);
            }

            Console.WriteLine("--- Complete initialization of 3-dimensional B-spline.");
            bcbspline.CompleteIt();
            bspline.CompleteIt();
            Console.WriteLine("--- Save is not implemented in the managed port; skipping.");
        }

        private static void EvalBSpline1(EquidistantBSpline1 bcbspline, XBSTools1<int, int> bspline)
        {
            const int n = 10;
            var x = new[] { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
            var f = new double[n];
            var g = new double[n];

            Console.WriteLine("--- Scalar evaluation of 1-dimensional B-spline function.");
            var fs = bcbspline.Evaluate(2.0);
            var gs = bspline.Evaluate(2.0);
            AssertClose(fs, gs, "Different B-spline values.");
            Console.WriteLine("    {0} = B-spline(2.0)", fs);

            Console.WriteLine("--- Vectorial evaluation of 1-dimensional B-spline function.");
            bcbspline.Evaluate(x, n, f);
            bspline.Evaluate(x, n, g);
            for (var i = 0; i < n; i++)
            {
                AssertClose(f[i], g[i], "Different B-spline values.");
                Console.WriteLine("    {0} = B-spline({1})", f[i], x[i]);
            }

            Console.WriteLine("--- Scalar evaluation of inverse function.");
            var invX = bcbspline.Inverse(0.5);
            var invR = bspline.Inverse(0.5);
            AssertClose(invX, invR, "Different values for inverse of B-spline function.");
            var fv = bcbspline.Evaluate(invX);
            var gv = bspline.Evaluate(invR);
            AssertClose(fv, gv, "Different B-spline values.");
            Console.WriteLine("    {0} = inv(0.5), {1} = B-spline({0})", invX, fv);

            Console.WriteLine("--- Konstruktion values:");
            Console.WriteLine("    Maximal number of subintervals [in x-direction]: {0}", bcbspline.GetMaxIntervals());
            Console.WriteLine("    Maximal order of B-splines     [in x-direction]: {0}", bcbspline.GetMaxOrder());
            Console.WriteLine("    Lower boundary                 [in x-direction]: {0}", bcbspline.GetLowerBoundary());
            Console.WriteLine("    Upper boundary                 [in x-direction]: {0}", bcbspline.GetUpperBoundary());
        }

        private static void EvalBSpline2()
        {
            Console.WriteLine("--- 2-dimensional evaluation skipped: BSTools2 is not implemented in this port.");
        }

        private static void EvalBSpline3(EquidistantBSpline3 bcbspline, XBSTools3Pre2<int, int, int, int, int, int> bspline)
        {
            Console.WriteLine("--- Scalar evaluations of 3-dimensional B-spline function.");
            var f = bcbspline.Evaluate(2.0, 1.0, 0.5);
            var g = bspline.Evaluate(2.0, 1.0, 0.5);
            AssertClose(f, g, "Different B-spline values.");
            Console.WriteLine("    {0} = B-spline(2.0, 1.0, 0.5)", f);

            f = bcbspline.Evaluate(1.5);
            g = bspline.Evaluate(1.5);
            AssertClose(f, g, "Different B-spline values.");
            Console.WriteLine("    {0} = B-spline(1.5)", f);

            const int n = 10;
            var x = new[] { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
            var fv = new double[n];
            var gv = new double[n];

            Console.WriteLine("--- Vectorial evaluations of 3-dimensional B-spline function.");
            bcbspline.Evaluate(x, n, 1.0, 0.5, fv);
            bspline.Evaluate(x, n, 1.0, 0.5, gv);
            for (var i = 0; i < n; i++)
            {
                AssertClose(fv[i], gv[i], "Different B-spline values.");
                Console.WriteLine("    {0} = B-spline({1}, 1.0, 0.5)", fv[i], x[i]);
            }

            Console.WriteLine("--- Vectorial evaluations with predefined y- and z-value.");
            bcbspline.Evaluate(x, n, fv);
            bspline.Evaluate(x, n, gv);
            for (var i = 0; i < n; i++)
            {
                AssertClose(fv[i], gv[i], "Different B-spline values.");
                Console.WriteLine("    {0} = B-spline({1})", fv[i], x[i]);
            }

            Console.WriteLine("--- Scalar evaluation of x-inverse function.");
            var invX = bcbspline.InvX(5.0, 1.5, 0.5);
            var invR = bspline.InvX(5.0, 1.5, 0.5);
            AssertClose(invX, invR, "Different values for inverse of B-spline function.");
            var f1 = bcbspline.Evaluate(invX, 1.5, 0.5);
            var g1 = bspline.Evaluate(invR, 1.5, 0.5);
            AssertClose(f1, g1, "Different B-spline values.");
            Console.WriteLine("    {0} = invX(5.0,1.5,0.5), {1} = B-spline({0}, 1.5, 0.5)", invX, f1);

            Console.WriteLine("--- Konstruktion values:");
            Console.WriteLine("    Maximal number of subintervals in x-direction: {0}", bcbspline.GetMaxIntervalsX());
            Console.WriteLine("    Maximal order of B-splines     in x-direction: {0}", bcbspline.GetMaxOrderX());
            Console.WriteLine("    Lower boundary                 in x-direction: {0}", bcbspline.GetLowerBoundaryX());
            Console.WriteLine("    Upper boundary                 in x-direction: {0}", bcbspline.GetUpperBoundaryX());
            Console.WriteLine("    Maximal number of subintervals in y-direction: {0}", bcbspline.GetMaxIntervalsY());
            Console.WriteLine("    Maximal order of B-splines     in y-direction: {0}", bcbspline.GetMaxOrderY());
            Console.WriteLine("    Lower boundary                 in y-direction: {0}", bcbspline.GetLowerBoundaryY());
            Console.WriteLine("    Upper boundary                 in y-direction: {0}", bcbspline.GetUpperBoundaryY());
            Console.WriteLine("    Maximal number of subintervals in z-direction: {0}", bcbspline.GetMaxIntervalsZ());
            Console.WriteLine("    Maximal order of B-splines     in z-direction: {0}", bcbspline.GetMaxOrderZ());
            Console.WriteLine("    Lower boundary                 in z-direction: {0}", bcbspline.GetLowerBoundaryZ());
            Console.WriteLine("    Upper boundary                 in z-direction: {0}", bcbspline.GetUpperBoundaryZ());
        }

        private static void TrafoBSpline1(EquidistantBSpline1 bcbspline, XBSTools1<int, int> bspline)
        {
            Console.WriteLine("--- Replace 1-dim. B-spline function with its derivative.");
            var dbcbspline = new EquidistantBSpline1(MaxIntervalsX, MaxOrderX);
            var dbspline = new XBSTools1<int, int>(MaxIntervalsX, MaxOrderX);
            CopySpline1(bspline, dbcbspline, dbspline);
            dbcbspline.Derivative();
            dbspline.Diff();

            Console.WriteLine("--- Replace 1-dim. B-spline function with its antiderivative.");
            var ibcbspline = new EquidistantBSpline1(MaxIntervalsX, MaxOrderX);
            var ibspline = new XBSTools1<int, int>(MaxIntervalsX, MaxOrderX);
            CopySpline1(bspline, ibcbspline, ibspline);
            ibcbspline.AntiDerivative(0.0, 1.0);
            ibspline.Integral(0.0, 1.0);

            Console.WriteLine("--- Binary data transfer is not implemented in the managed port; skipping.");
        }

        private static void TrafoBSpline2()
        {
            Console.WriteLine("--- 2-dimensional transformations skipped: BSTools2 is not implemented in this port.");
        }

        private static void TrafoBSpline3()
        {
            Console.WriteLine("--- 3-dimensional transformations skipped: not implemented in this port.");
        }

        private static void BasicBSplines()
        {
            Console.WriteLine("--- Basic B-splines output skipped: Save is not implemented in this port.");
        }

        private static void CopySpline1(
            XBSTools1<int, int> source,
            EquidistantBSpline1 destWrapper,
            XBSTools1<int, int> dest)
        {
            destWrapper.SetCarrier(source.LowerBound(), source.UpperBound(), source.GetNumberOfIntervals(), source.GetOrder());
            dest.SetCarrier(source.LowerBound(), source.UpperBound(), source.GetNumberOfIntervals(), source.GetOrder());
            var count = source.GetNumberOfAllWeights();
            for (var i = 0; i < count; i++)
            {
                var weight = source.GetWeight(i);
                destWrapper.SetWeight(i, weight);
                dest.SetWeight(i, weight);
            }

            destWrapper.CompleteIt();
            dest.CompleteIt();
        }

        private static void AssertEqual(int expected, int actual, string message)
        {
            if (expected != actual)
            {
                throw new InvalidOperationException(message);
            }
        }

        private static void AssertClose(double expected, double actual, string message)
        {
            if (Math.Abs(expected - actual) > Epsilon)
            {
                throw new InvalidOperationException(message);
            }
        }
    }
}
