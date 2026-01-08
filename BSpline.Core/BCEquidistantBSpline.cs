namespace BSpline.Core
{
    public static class BCEquidistantBSpline
    {
        public static void Assert(bool condition, string message)
        {
            if (!condition)
            {
                BCException.Throw(message);
            }
        }
    }
}
