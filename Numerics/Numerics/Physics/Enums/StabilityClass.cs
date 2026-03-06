namespace CSharpNumerics.Physics.Enums
{
    /// <summary>
    /// Pasquill–Gifford atmospheric stability classes used by the Gaussian plume model.
    /// </summary>
    public enum StabilityClass
    {
        /// <summary>Very unstable — strong insolation, light wind.</summary>
        A,
        /// <summary>Unstable — moderate insolation.</summary>
        B,
        /// <summary>Slightly unstable — weak insolation.</summary>
        C,
        /// <summary>Neutral — overcast day or night, moderate wind.</summary>
        D,
        /// <summary>Slightly stable — night, light wind, partial cloud.</summary>
        E,
        /// <summary>Stable — night, light wind, clear sky.</summary>
        F
    }
}
