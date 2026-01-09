using System;

namespace BSpline.Core
{
    public static class EquidistantBSpline
    {
        public static void Assert(bool condition, string message)
        {
            if (!condition)
            {
                throw new Exception(message);
            }
        }
    }
}
