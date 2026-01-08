using System;

namespace BSpline.Core
{
    public sealed class BCException : Exception
    {
        public BCException(string message)
            : base(message)
        {
        }

        public static void Throw(string message)
        {
            throw new BCException(message);
        }
    }
}
