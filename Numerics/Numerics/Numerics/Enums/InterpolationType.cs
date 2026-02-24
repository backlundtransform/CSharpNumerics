namespace CSharpNumerics.Numerics.Enums;

    public enum InterpolationType
    {
        Linear,
        Logarithmic, 
        LinLog,
        LogLin,

        /// <summary>Lagrange polynomial interpolation through all data points.</summary>
        Polynomial,

        /// <summary>Natural cubic spline interpolation (smooth, C² continuous).</summary>
        CubicSpline,

        /// <summary>Trigonometric interpolation (best for periodic data).</summary>
        Trigonometric
    }

