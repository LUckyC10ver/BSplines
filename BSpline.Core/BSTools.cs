using System;
using System.IO;

namespace BSpline.Core
{
    public static class BSToolsConstants
    {
        public const int BSTOOLS_TOTAL_MAX_ORDER = 6;
        public const int BSTOOLS_TOTAL_NUMBER_DEPENDENT = BSTOOLS_TOTAL_MAX_ORDER;
        public const int BSTOOLS_TOTAL_SIZE_TAYLOR = BSTOOLS_TOTAL_MAX_ORDER;
        public const double Delta = 1.0e-8;

        public static int GetSizeTaylor(int maxOrder) => maxOrder;

        public static int GetNumberDependent(int maxOrder) => maxOrder;

        public static int GetSizeDeBoor(int maxInterval, int maxOrder) => maxInterval + maxOrder - 1;

        public static int GetSizeCarrier(int maxInterval, int maxOrder) => maxInterval + 2 * maxOrder - 3;
    }

    public sealed class BSToolsData
    {
        public double A { get; private set; }
        public double B { get; private set; }
        public int K { get; private set; }
        public int N { get; private set; }

        public int Set(int nSrc, double[] src, int nDest, double[] dest)
        {
            if (src == null || dest == null)
            {
                return 0;
            }

            if (nSrc == nDest)
            {
                Array.Copy(src, dest, nSrc);
                return nSrc;
            }

            return 0;
        }

        public int Get(out double a, out double b, out int n, out int k)
        {
            a = A;
            b = B;
            n = N;
            k = K;
            return 1;
        }

        public int Set(double a, double b, int n, int k)
        {
            A = a;
            B = b;
            N = n;
            K = k;
            return 1;
        }

        public void OutTaylor(char axis, TextWriter writer, int lowerOrder, int upperOrder)
        {
            writer.Write('T');
            writer.Write(axis);
            writer.Write(' ');
            writer.Write(lowerOrder);
            writer.Write(' ');
            writer.Write(upperOrder);
            writer.Write(' ');
        }

        public void InTaylor(char axis, TextReader reader, out int lowerOrder, out int upperOrder)
        {
            var id = ReadToken(reader);
            if (id == "T")
            {
                var axisToken = ReadToken(reader);
                if (!string.IsNullOrEmpty(axisToken) && axisToken[0] == axis)
                {
                    lowerOrder = ReadInt(reader);
                    upperOrder = ReadInt(reader);
                    return;
                }
            }

            lowerOrder = 0;
            upperOrder = 0;
        }

        public void OutWeight(TextWriter writer, double[] weights, int length)
        {
            writer.Write('W');
            writer.Write(' ');
            writer.Write(length);
            writer.Write(' ');
            for (var i = 0; i < length; i++)
            {
                writer.Write(weights[i]);
                writer.Write(' ');
            }
        }

        public void InWeight(TextReader reader, int maxLength, double[] weights, out int length)
        {
            var id = ReadToken(reader);
            if (id == "W")
            {
                length = ReadInt(reader);
                if (length <= maxLength)
                {
                    for (var i = 0; i < maxLength; i++)
                    {
                        weights[i] = ReadDouble(reader);
                    }

                    return;
                }
            }

            length = 0;
            if (weights.Length > 0)
            {
                weights[0] = 0;
            }
        }

        public void OutCarrier(char axis, TextWriter writer)
        {
            writer.Write('C');
            writer.Write(axis);
            writer.Write(' ');
            writer.Write(A);
            writer.Write(' ');
            writer.Write(B);
            writer.Write(' ');
            writer.Write(N);
            writer.Write(' ');
            writer.Write(K);
            writer.Write(' ');
        }

        public void InCarrier(char axis, TextReader reader)
        {
            var id = ReadToken(reader);
            if (id == "C")
            {
                var axisToken = ReadToken(reader);
                if (!string.IsNullOrEmpty(axisToken) && axisToken[0] == axis)
                {
                    A = ReadDouble(reader);
                    B = ReadDouble(reader);
                    N = ReadInt(reader);
                    K = ReadInt(reader);
                    return;
                }
            }

            A = 0;
            B = 1;
            N = 1;
            K = 1;
        }

        public void Print(TextWriter writer, string path, string name)
        {
            writer.WriteLine($"{path}/{name}: a={A} b={B} k={K} n={N}");
        }

        public bool IsValid(int maxInterval, int maxOrder)
        {
            if (B <= A)
            {
                return false;
            }

            if (N == 0 || maxInterval < N)
            {
                return false;
            }

            if (K == 0 || maxOrder < K)
            {
                return false;
            }

            return true;
        }

        public int CalcWLen()
        {
            return BSTools.CalcWLen(N, K);
        }

        private static string ReadToken(TextReader reader)
        {
            int ch;
            do
            {
                ch = reader.Read();
                if (ch == -1)
                {
                    return string.Empty;
                }
            } while (char.IsWhiteSpace((char)ch));

            var buffer = new System.Text.StringBuilder();
            while (ch != -1 && !char.IsWhiteSpace((char)ch))
            {
                buffer.Append((char)ch);
                ch = reader.Read();
            }

            return buffer.ToString();
        }

        private static int ReadInt(TextReader reader)
        {
            var token = ReadToken(reader);
            return int.TryParse(token, out var value) ? value : 0;
        }

        private static double ReadDouble(TextReader reader)
        {
            var token = ReadToken(reader);
            return double.TryParse(token, out var value) ? value : 0.0;
        }
    }

    public class BSTools
    {
        protected readonly int MaxK;
        protected readonly int MaxN;
        protected readonly int SizeOfT;
        protected readonly int SizeOfDeltaTInv;
        protected readonly double[] T;
        protected readonly double[] DeltaTInv;

        protected double A;
        protected double B;
        protected int K;
        protected int N;
        protected int WLen;
        protected double H;
        protected int MaxIntervalIndex;

        protected BSTools(int sizeOfT, double[] t, int sizeOfDeltaTInv, double[] deltaTInv, int maxInterval, int maxOrder)
        {
            SizeOfT = sizeOfT;
            T = t;
            SizeOfDeltaTInv = sizeOfDeltaTInv;
            DeltaTInv = deltaTInv;
            MaxN = maxInterval;
            MaxK = maxOrder;
            A = 0;
            B = 1;
            K = 1;
            N = 1;
        }

        public static int CalcWLen(int n, int k)
        {
            return BSToolsConstants.GetSizeDeBoor(n, k);
        }

        public int GetTLen()
        {
            if (K == 1)
            {
                return N + 1;
            }

            return N + 2 * K - 3;
        }

        public int GetWLen() => WLen;

        public void SetWLen()
        {
            WLen = CalcWLen(N, K);
        }

        public bool Load(object simpleParser)
        {
            if (simpleParser is not SimpleParser parser)
            {
                return false;
            }

            double? a = null;
            double? b = null;
            int? n = null;
            int? k = null;

            while (parser.GetNewLine())
            {
                var line = parser.GetData();
                var key = ExtractKey(line);
                if (string.IsNullOrEmpty(key))
                {
                    continue;
                }

                switch (key)
                {
                    case "a":
                        if (parser.GetValue(out double aValue))
                        {
                            a = aValue;
                        }
                        break;
                    case "b":
                        if (parser.GetValue(out double bValue))
                        {
                            b = bValue;
                        }
                        break;
                    case "n":
                        if (parser.GetValue(out int nValue))
                        {
                            n = nValue;
                        }
                        break;
                    case "k":
                        if (parser.GetValue(out int kValue))
                        {
                            k = kValue;
                        }
                        break;
                }
            }

            if (!a.HasValue || !b.HasValue || !n.HasValue || !k.HasValue)
            {
                return false;
            }

            InitInternal(a.Value, b.Value, n.Value, k.Value);
            return true;
        }

        public bool Save(TextWriter writer, string path, string varName)
        {
            if (writer == null)
            {
                return false;
            }

            var parser = new SimpleParser("BSTools.Save", writer: writer);
            if (!parser.SetPath(path, varName))
            {
                return false;
            }

            parser.Write("a", A);
            parser.Write("b", B);
            parser.Write("n", N);
            parser.Write("k", K);
            return true;
        }

        public double GetLowerBoundary() => A;

        public double GetUpperBoundary() => B;

        public int GetOrder() => K;

        public int GetNumberOfIntervals() => N;

        public int GetBinaryData(BSToolsData data)
        {
            if (data == null)
            {
                return 0;
            }

            return data.Set(A, B, N, K);
        }

        public int SetBinaryData(BSToolsData data)
        {
            if (data == null)
            {
                return 0;
            }

            return data.Get(out A, out B, out N, out K);
        }

        public static bool IsEqual(double[] x, int lenX, double[] y, int lenY)
        {
            if (lenX != lenY || x == null || y == null)
            {
                return false;
            }

            for (var i = 0; i < lenX; i++)
            {
                if (x[i] != y[i])
                {
                    return false;
                }
            }

            return true;
        }

        public static bool IsApproximatedEqual(double[] x, int lenX, double[] y, int lenY, double maxDelta)
        {
            if (lenX != lenY || x == null || y == null)
            {
                return false;
            }

            for (var i = 0; i < lenX; i++)
            {
                if (Math.Abs(x[i] - y[i]) > maxDelta)
                {
                    return false;
                }
            }

            return true;
        }

        protected int Max(int x1, int x2) => Math.Max(x1, x2);

        protected BSTools Max(BSTools bs)
        {
            var oldA = A;
            var oldB = B;
            A = A < bs.A ? A : bs.A;
            B = B < bs.B ? bs.B : B;
            K = Max(K, bs.K);

            if (N < bs.N)
            {
                N = (int)(bs.N * (B - A) / (bs.B - bs.A));
            }
            else
            {
                N = (int)(N * (B - A) / (oldB - oldA));
            }

            BuildEquidistantCarrier();
            return this;
        }

        internal int MaxInternal(int x1, int x2)
        {
            return Max(x1, x2);
        }

        internal BSTools MaxInternal(BSTools bs)
        {
            return Max(bs);
        }

        protected void BuildEquidistantCarrier()
        {
            SetWLen();
            var tLen = GetTLen();
            var dim = SizeOfT < tLen ? SizeOfT : tLen;
            var outOf = 0.0;

            if (K > 2)
            {
                outOf = (B - A) / (GetWLen() + 1 - K) * (K - 2);
            }

            var length = B - A + 2 * outOf;
            for (var i = 0; i < dim; i++)
            {
                T[i] = A - outOf + length * i / (dim - 1.0);
            }

            CompleteIt();
            InitEquidistant();
        }

        internal void BuildEquidistantCarrierInternal()
        {
            BuildEquidistantCarrier();
        }

        internal double LowerBoundary => A;

        internal double UpperBoundary => B;

        internal int Order => K;

        internal int IntervalCount => N;

        internal double StepSize => H;

        internal int MaxOrder => MaxK;

        internal double[] Carrier => T;

        protected void Init(double a, double b, int n, int k)
        {
            if (SizeOfT < BSToolsConstants.GetSizeCarrier(n, k) ||
                SizeOfDeltaTInv < BSToolsConstants.GetNumberDependent(k))
            {
                return;
            }

            A = a;
            B = b;
            N = n;
            K = k;
            BuildEquidistantCarrier();
        }

        internal void InitInternal(double a, double b, int n, int k)
        {
            Init(a, b, n, k);
        }

        protected void DiffCarrier()
        {
            if (K > 1)
            {
                K--;
                Array.Copy(T, 1, T, 0, GetTLen());
            }

            CompleteIt();
            InitEquidistant();
        }

        internal void DiffCarrierInternal()
        {
            DiffCarrier();
        }

        protected void DiffWeights(double[] sourceWeights, double[] destWeights)
        {
            if (K == 1)
            {
                var n = GetWLen();
                for (var i = 0; i < n; i++)
                {
                    destWeights[i] = 0.0;
                }
            }
            else
            {
                var odh = 1 / H;
                var n = GetWLen();
                for (var i = 1; i < n; i++)
                {
                    destWeights[i - 1] = (sourceWeights[i] - sourceWeights[i - 1]) * odh;
                }
            }
        }

        internal void DiffWeightsInternal(double[] sourceWeights, double[] destWeights)
        {
            DiffWeights(sourceWeights, destWeights);
        }

        protected void IntWeights(double[] sourceWeights, double[] destWeights)
        {
            destWeights[0] = 0.0;

            for (var i = 1; i < GetWLen(); i++)
            {
                destWeights[i] = H * sourceWeights[i - 1] + destWeights[i - 1];
            }
        }

        internal void IntWeightsInternal(double[] sourceWeights, double[] destWeights)
        {
            IntWeights(sourceWeights, destWeights);
        }

        protected void CompleteIt()
        {
            SetWLen();

            for (var i = GetTLen(); i < SizeOfT; i++)
            {
                T[i] = 0.0;
            }
        }

        internal double GetWeightPosition(int index)
        {
            if (K % 2 == 1)
            {
                if (K == 1)
                {
                    return T[index];
                }

                var i = index + K / 2 - 1;
                return (T[i] + T[i + 1]) / 2;
            }

            return T[index + K / 2 - 1];
        }

        internal enum Position
        {
            Ober,
            Mitte,
            Unter
        }

        protected Position GetLocal(double x)
        {
            if (x < A)
            {
                return Position.Unter;
            }

            if (x > B)
            {
                return Position.Ober;
            }

            return Position.Mitte;
        }

        internal Position GetLocalInternal(double x)
        {
            return GetLocal(x);
        }

        protected int GetFGrad(char f)
        {
            switch (f)
            {
                case 'a':
                    return K;
                case 's':
                    return K == 0 ? 0 : K - 1;
                default:
                    var n = f - '0';
                    return K >= n ? n : 0;
            }
        }

        internal int GetFGradInternal(char f)
        {
            return GetFGrad(f);
        }

        protected void ErrorMsg(string file, string function, string message, double errno)
        {
            var fullMessage = $"{file}:{function}: {message}: {errno}";
            Console.Error.WriteLine(fullMessage);
        }

        protected int Check(string proc)
        {
            if (K == 0 || K > MaxK)
            {
                return 1;
            }

            if (N == 0 || N > MaxN)
            {
                return 1;
            }

            if (A > B)
            {
                return 1;
            }

            return 0;
        }

        protected int GetIntervalIndex(double x)
        {
            var i = (int)((x - T[0]) * DeltaTInv[1]);
            return i <= MaxIntervalIndex ? i : i - 1;
        }

        internal int GetIntervalIndexInternal(double x)
        {
            return GetIntervalIndex(x);
        }

        protected void DeBoor(double x, int k, double[] d, double[] t, out double y)
        {
            for (var i = 0; i < k - 1; i++)
            {
                for (var j = k - 2; j + 1 > i; j--)
                {
                    var alpha = (x - t[j]) / (t[j + k - i - 1] - t[j]);
                    d[j + 1] = (1 - alpha) * d[j] + alpha * d[j + 1];
                }
            }

            y = d[k - 1];
        }

        protected double DeBoorEqu(double x, int k, double[] d, double[] t)
        {
            for (var i = 0; i < k - 1; i++)
            {
                for (var j = 0; j < k - i - 1; j++)
                {
                    var alpha = (x - t[j + i]) * DeltaTInv[k - 1 - i];
                    d[j] = (1 - alpha) * d[j] + alpha * d[j + 1];
                }
            }

            return d[0];
        }

        protected double GetFctValue(double x, int k, double[] d, double[] t)
        {
            switch (k)
            {
                case 1:
                    return d[0];
                case 2:
                    return GetFctValueK2(x, d, t);
                case 3:
                    return GetFctValueK3(x, d, t);
                default:
                    return GetFctValueKn(x, k, d, t);
            }
        }

        protected double GetFctValueKn(double x, int k, double[] d, double[] t)
        {
            var dd = new double[BSToolsConstants.BSTOOLS_TOTAL_NUMBER_DEPENDENT];
            Array.Copy(d, dd, k);
            return DeBoorEqu(x, k, dd, t);
        }

        protected double GetFctValueK2(double x, double[] d, double[] t)
        {
            var alpha = (x - t[0]) * DeltaTInv[1];
            return (1 - alpha) * d[0] + alpha * d[1];
        }

        protected double GetFctValueK3(double x, double[] d, double[] t)
        {
            var a0 = (x - t[0]) * DeltaTInv[2];
            var xT1 = x - t[1];
            var a1 = xT1 * DeltaTInv[2];
            var a01 = xT1 * DeltaTInv[1];
            return (1 - a01) * ((1 - a0) * d[0] + a0 * d[1]) + a01 * ((1 - a1) * d[1] + a1 * d[2]);
        }

        protected void Ablwert(double x, int k, double[] d, double[] t, out double y)
        {
            switch (k)
            {
                case 1:
                    y = 0;
                    break;
                case 2:
                    y = (d[1] - d[0]) / (t[1] - t[0]);
                    break;
                default:
                    var dd = new double[BSToolsConstants.BSTOOLS_TOTAL_NUMBER_DEPENDENT];
                    for (var i = 1; i < k; i++)
                    {
                        dd[i - 1] = (k - 1.0) * (d[i] - d[i - 1]) / (t[i + k - 2] - t[i - 1]);
                    }

                    var tSlice = Slice(t, 1, t.Length - 1);
                    DeBoor(x, k - 1, dd, tSlice, out y);
                    break;
            }
        }

        protected int LowIndex()
        {
            return K == 1 ? 0 : K - 2;
        }

        protected int HighIndex()
        {
            return LowIndex() + N - 1;
        }

        protected int GetMainIndex(int i)
        {
            switch (K)
            {
                case 1:
                    return i;
                case 2:
                    return GetMainIndexK2(i);
                case 3:
                    return GetMainIndexK3(i);
                default:
                    var ip2 = i + 2;
                    return ip2 < K ? 0 : ip2 - K;
            }
        }

        protected int GetMainIndexK2(int i)
        {
            return i;
        }

        protected int GetMainIndexK3(int i)
        {
            return i == 0 ? 0 : i - 1;
        }

        protected void TaylorPoly(double[] d, double x, int n, double[] poly)
        {
            if (K < n)
            {
                return;
            }

            if (n > 0)
            {
                if (K == 1)
                {
                    poly[0] = d[0];
                }
                else
                {
                    var r = GetMainIndex(GetIntervalIndex(x));
                    var dd = new double[BSToolsConstants.BSTOOLS_TOTAL_NUMBER_DEPENDENT];
                    var kk = K;
                    var nn = n - 1;
                    var f = 1.0;

                    Array.Copy(d, dd, K);
                    var tSlice = Slice(T, r, T.Length - r);
                    poly[nn] = GetFctValue(x, K, dd, tSlice);

                    for (var i = 1; i < n; i++)
                    {
                        f *= i;

                        for (var j = 1; j < kk; j++)
                        {
                            dd[j - 1] = (kk - 1.0) * (dd[j] - dd[j - 1]) / (T[r + j + kk - 3 + i] - T[r + j + i - 2]);
                        }

                        kk--;
                        nn--;
                        var tSliceStep = Slice(T, r + i, T.Length - (r + i));
                        poly[nn] = GetFctValue(x, kk, dd, tSliceStep) / f;
                    }
                }
            }
            else
            {
                poly[0] = 0.0;
            }
        }

        protected double PolyVal(double[] p, int n, double x)
        {
            var y = 0.0;
            for (var i = 0; i < n; i++)
            {
                y = y * x + p[i];
            }

            return y;
        }

        protected void PolyDiff(double[] p, ref int n)
        {
            var g = 0;
            if (n > 0)
            {
                n--;
            }

            for (var i = n; i > 0; i--)
            {
                p[g] = i * p[g];
                g++;
            }
        }

        internal int GetCarrierEdge(char edge, double[] buffer, int bufferLen)
        {
            var n = GetWLen();
            if (bufferLen < n)
            {
                n = bufferLen;
            }

            for (var i = 0; i < n; i++)
            {
                if (K == 1)
                {
                    if (edge == 'u')
                    {
                        buffer[i] = i + 1 < GetTLen() - 1 ? T[i + 1] : 1e38;
                    }
                    else
                    {
                        buffer[i] = i > 0 ? T[i] : -1e38;
                    }
                }
                else
                {
                    if (edge == 'u')
                    {
                        if (i + K - 1 < GetTLen())
                        {
                            buffer[i] = T[i + K - 1];
                            if (GetLocal(buffer[i]) == Position.Ober)
                            {
                                buffer[i] = 1e38;
                            }
                        }
                        else
                        {
                            buffer[i] = 1e38;
                        }
                    }
                    else
                    {
                        if (i > 0)
                        {
                            buffer[i] = T[i - 1];
                            if (GetLocal(buffer[i]) == Position.Unter)
                            {
                                buffer[i] = -1e38;
                            }
                        }
                        else
                        {
                            buffer[i] = -1e38;
                        }
                    }
                }
            }

            return n;
        }

        protected double GetLowerCarrierEdge(int index)
        {
            if (K == 1)
            {
                return index > 0 ? T[index] : -1e38;
            }

            if (index > 0)
            {
                var value = T[index - 1];
                return GetLocal(value) == Position.Unter ? -1e38 : value;
            }

            return -1e38;
        }

        protected double GetUpperCarrierEdge(int index)
        {
            if (K == 1)
            {
                return index + 1 < GetTLen() - 1 ? T[index + 1] : 1e38;
            }

            if (index + K - 1 < GetTLen())
            {
                var value = T[index + K - 1];
                return GetLocal(value) == Position.Ober ? 1e38 : value;
            }

            return 1e38;
        }

        internal void GetDomainSubWeightIndexes(double a, double b, out int ia, out int ib)
        {
            var n = GetWLen();
            var i = 0;

            for (i = 0; i < n; i++)
            {
                if (a <= GetLowerCarrierEdge(i))
                {
                    break;
                }
            }

            ia = i < n ? i : n - 1;

            for (i = n; i > 0; i--)
            {
                if (GetUpperCarrierEdge(i - 1) <= b)
                {
                    break;
                }
            }

            ib = i > 0 ? i - 1 : 0;
        }

        protected void InitEquidistant()
        {
            SetWLen();

            H = T[1] - T[0];
            DeltaTInv[0] = 1;

            if (K < 2)
            {
                DeltaTInv[1] = 1 / H;
            }
            else
            {
                for (var i = 1; i < K; i++)
                {
                    DeltaTInv[i] = 1 / (i * H);
                }
            }

            MaxIntervalIndex = (int)((B - 0.5 * H - T[0]) * DeltaTInv[1]);
        }

        internal int GetMainIndexInternal(int i) => GetMainIndex(i);

        internal int GetMainIndexK2Internal(int i) => GetMainIndexK2(i);

        internal int GetMainIndexK3Internal(int i) => GetMainIndexK3(i);

        internal double GetFctValueKnInternal(double x, int k, double[] d, double[] t) => GetFctValueKn(x, k, d, t);

        internal double GetFctValueK2Internal(double x, double[] d, double[] t) => GetFctValueK2(x, d, t);

        internal double GetFctValueK3Internal(double x, double[] d, double[] t) => GetFctValueK3(x, d, t);

        internal double GetFctValueInternal(double x, int k, double[] d, double[] t) => GetFctValue(x, k, d, t);

        internal double PolyValInternal(double[] p, int n, double x) => PolyVal(p, n, x);

        internal void TaylorPolyInternal(double[] d, double x, int n, double[] poly) => TaylorPoly(d, x, n, poly);

        internal void PolyDiffInternal(double[] p, ref int n) => PolyDiff(p, ref n);

        internal void AblwertInternal(double x, int k, double[] d, double[] t, out double y) => Ablwert(x, k, d, t, out y);

        public bool IsNotEqual(BSTools other)
        {
            if (other == null)
            {
                return true;
            }

            if (A != other.A || B != other.B || N != other.N || K != other.K)
            {
                return true;
            }

            return false;
        }

        public bool IsApproximatedEqual(BSTools other, double maxDelta)
        {
            if (other == null)
            {
                return false;
            }

            if (maxDelta < Math.Abs(A - other.A)) return false;
            if (maxDelta < Math.Abs(B - other.B)) return false;
            if (maxDelta < Math.Abs(N - other.N)) return false;
            if (maxDelta < Math.Abs(K - other.K)) return false;

            return IsApproximatedEqual(T, GetTLen(), other.T, other.GetTLen(), maxDelta);
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

    public sealed class XBSTools<TInterval, TOrder> : BSTools
        where TInterval : struct
        where TOrder : struct
    {
        private readonly double[] _tMem;
        private readonly double[] _deltaTInvMem;

        public XBSTools(int maxInterval, int maxOrder)
            : base(
                BSToolsConstants.GetSizeCarrier(maxInterval, maxOrder),
                new double[BSToolsConstants.GetSizeCarrier(maxInterval, maxOrder)],
                BSToolsConstants.GetNumberDependent(maxOrder),
                new double[BSToolsConstants.GetNumberDependent(maxOrder)],
                maxInterval,
                maxOrder)
        {
            _tMem = T;
            _deltaTInvMem = DeltaTInv;
            Init(0, 1, 1, 1);
        }

        public XBSTools(int maxInterval, int maxOrder, double a, double b, int n, int k)
            : this(maxInterval, maxOrder)
        {
            Init(a, b, n, k);
        }

        public XBSTools(int maxInterval, int maxOrder, BSTools input)
            : this(maxInterval, maxOrder)
        {
            if (input == null)
            {
                return;
            }

            Init(input.GetLowerBoundary(), input.GetUpperBoundary(), input.GetNumberOfIntervals(), input.GetOrder());
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
}
