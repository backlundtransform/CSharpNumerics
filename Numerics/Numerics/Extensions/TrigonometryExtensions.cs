namespace System
{
    public static class TrigonometryExtensions
    {
        public static double DegreeToRadians(this double angle)
        {
            return (angle * Math.PI) / 180;
        }
    }
}
