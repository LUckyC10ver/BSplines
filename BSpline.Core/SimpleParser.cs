using System;
using System.IO;
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
            return false;
        }

        public bool GetValue(out int value)
        {
            value = 0;
            return false;
        }

        public int GetArray(double[] dest, int maxLength)
        {
            return 0;
        }

        public int GetArray(int[] dest, int maxLength)
        {
            return 0;
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
    }
}
