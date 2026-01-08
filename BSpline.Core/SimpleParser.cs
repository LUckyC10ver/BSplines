using System;
using System.Globalization;
using System.IO;
using System.Text;
namespace BSpline.Core
{
    public sealed class SimpleParser
    {
        private const int PathLength = 256;
        private const int DataBufferLength = 256;

        private readonly string _appealFrom;
        private string _appealFromSubroutine;
        private TextReader _reader;
        private TextWriter _writer;
        private string _path = string.Empty;
        private string _data = string.Empty;

        public SimpleParser(string appealFrom, TextReader reader = null, TextWriter writer = null)
        {
            _appealFrom = appealFrom;
            _reader = reader;
            _writer = writer;
        }

        public void ResetCurrentSubroutine()
        {
            _appealFromSubroutine = null;
        }

        public void SetCurrentSubroutine(string functionName)
        {
            _appealFromSubroutine = functionName;
        }

        public void SetData(string data)
        {
            _data = data ?? string.Empty;
        }

        public string GetData()
        {
            return _data;
        }

        public string GetPath()
        {
            return _path;
        }

        public bool SetPath(string prePath, string varName = "")
        {
            var combined = string.IsNullOrEmpty(varName) ? prePath : $"{prePath}.{varName}";
            if (combined.Length + 1 >= PathLength)
            {
                return false;
            }

            _path = combined;
            return true;
        }

        public bool WriteFctHeader(string filename, string structName = "")
        {
            if (_writer == null)
            {
                return false;
            }

            var name = Path.GetFileNameWithoutExtension(filename);
            _writer.WriteLine($"function {structName} = {name}");
            return true;
        }

        public bool Write(string varName, double[] array, int arrayLength, string comment = "", string format = "{0:0.################}")
        {
            if (_writer == null)
            {
                return false;
            }

            _writer.Write($"{_path}.{varName} = [");
            for (var i = 0; i < arrayLength; i++)
            {
                _writer.Write(string.Format(format, array[i]));
                if (i < arrayLength - 1)
                {
                    _writer.Write(",");
                }
            }
            _writer.WriteLine($"]; % {comment}");
            return true;
        }

        public bool Write(string varName, int[] array, int arrayLength, string comment = "")
        {
            if (_writer == null)
            {
                return false;
            }

            _writer.Write($"{_path}.{varName} = [");
            for (var i = 0; i < arrayLength; i++)
            {
                _writer.Write(array[i]);
                if (i < arrayLength - 1)
                {
                    _writer.Write(",");
                }
            }
            _writer.WriteLine($"]; % {comment}");
            return true;
        }

        public bool Write(string varName, double value, string comment = "", string format = "{0:0.################}")
        {
            if (_writer == null)
            {
                return false;
            }

            _writer.WriteLine($"{_path}.{varName} = {string.Format(format, value)}; % {comment}");
            return true;
        }

        public bool Write(string varName, int value, string comment = "")
        {
            if (_writer == null)
            {
                return false;
            }

            _writer.WriteLine($"{_path}.{varName} = {value}; % {comment}");
            return true;
        }

        public bool Write(string varName, string value, string comment = "")
        {
            if (_writer == null)
            {
                return false;
            }

            _writer.WriteLine($"{_path}.{varName} = '{value}'; % {comment}");
            return true;
        }

        public bool GetValue(out double value)
        {
            value = 0;
            var token = ExtractValueToken();
            if (string.IsNullOrWhiteSpace(token) || token.IndexOf('[') >= 0)
            {
                return false;
            }

            return double.TryParse(token, NumberStyles.Float, CultureInfo.InvariantCulture, out value);
        }

        public bool GetValue(out int value)
        {
            value = 0;
            var token = ExtractValueToken();
            if (string.IsNullOrWhiteSpace(token) || token.IndexOf('[') >= 0)
            {
                return false;
            }

            if (int.TryParse(token, NumberStyles.Integer, CultureInfo.InvariantCulture, out value))
            {
                return true;
            }

            if (double.TryParse(token, NumberStyles.Float, CultureInfo.InvariantCulture, out var doubleValue))
            {
                value = (int)Math.Round(doubleValue);
                return true;
            }

            return false;
        }

        public int GetArray(double[] dest, int maxLength)
        {
            if (dest == null || maxLength <= 0)
            {
                return 0;
            }

            var content = ExtractArrayToken();
            if (string.IsNullOrWhiteSpace(content))
            {
                return 0;
            }

            var count = 0;
            var parts = SplitArray(content);
            foreach (var part in parts)
            {
                if (count >= maxLength)
                {
                    break;
                }

                if (double.TryParse(part, NumberStyles.Float, CultureInfo.InvariantCulture, out var value))
                {
                    dest[count++] = value;
                }
            }

            return count;
        }

        public int GetArray(int[] dest, int maxLength)
        {
            if (dest == null || maxLength <= 0)
            {
                return 0;
            }

            var content = ExtractArrayToken();
            if (string.IsNullOrWhiteSpace(content))
            {
                return 0;
            }

            var count = 0;
            var parts = SplitArray(content);
            foreach (var part in parts)
            {
                if (count >= maxLength)
                {
                    break;
                }

                if (int.TryParse(part, NumberStyles.Integer, CultureInfo.InvariantCulture, out var value))
                {
                    dest[count++] = value;
                }
                else if (double.TryParse(part, NumberStyles.Float, CultureInfo.InvariantCulture, out var doubleValue))
                {
                    dest[count++] = (int)Math.Round(doubleValue);
                }
            }

            return count;
        }

        public bool GetNewLine()
        {
            if (_reader == null)
            {
                return false;
            }

            _data = _reader.ReadLine();
            return _data != null;
        }

        public TextReader OpenForLoading(string filename)
        {
            _reader = new StreamReader(filename);
            return _reader;
        }

        public void CantRecognise(string filename)
        {
            _writer?.WriteLine($"{_appealFrom}: Can't recognise line from {filename}.");
        }

        public void TooManyErrors(string filename)
        {
            _writer?.WriteLine($"{_appealFrom}: too many errors! \"{filename}\"");
        }

        private string ExtractValueToken()
        {
            if (string.IsNullOrWhiteSpace(_data))
            {
                return string.Empty;
            }

            var line = StripComment(_data);
            var index = line.IndexOf('=');
            if (index < 0 || index + 1 >= line.Length)
            {
                return string.Empty;
            }

            var valuePart = line.Substring(index + 1).Trim();
            if (valuePart.EndsWith(";"))
            {
                valuePart = valuePart.Substring(0, valuePart.Length - 1).Trim();
            }

            if (valuePart.Length >= 2 && valuePart.StartsWith("'") && valuePart.EndsWith("'"))
            {
                valuePart = valuePart.Substring(1, valuePart.Length - 2);
            }

            return valuePart;
        }

        private string ExtractArrayToken()
        {
            var token = ExtractValueToken();
            if (string.IsNullOrWhiteSpace(token))
            {
                return string.Empty;
            }

            var start = token.IndexOf('[');
            var end = token.LastIndexOf(']');
            if (start < 0 || end <= start)
            {
                return string.Empty;
            }

            return token.Substring(start + 1, end - start - 1);
        }

        private static string StripComment(string line)
        {
            var index = line.IndexOf('%');
            return index >= 0 ? line.Substring(0, index).Trim() : line.Trim();
        }

        private static string[] SplitArray(string content)
        {
            if (string.IsNullOrWhiteSpace(content))
            {
                return Array.Empty<string>();
            }

            var builder = new StringBuilder();
            foreach (var ch in content)
            {
                builder.Append(ch == ',' ? ' ' : ch);
            }

            return builder.ToString()
                .Split(new[] { ' ', '\t', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
        }
    }
}
