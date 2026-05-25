using System;

namespace CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

/// <summary>
/// Dynamic pressure computation and Max-Q detection for rocket ascent.
/// q = ½ρv² — the key structural load parameter during atmospheric flight.
/// Max-Q (peak dynamic pressure) occurs when the product ρ·v² is maximized,
/// typically between 11–14 km altitude for orbital launch vehicles.
/// </summary>
public static class DynamicPressure
{
    /// <summary>
    /// Computes dynamic pressure q = ½ρv².
    /// </summary>
    /// <param name="density">Air density in kg/m³.</param>
    /// <param name="velocity">True airspeed in m/s.</param>
    /// <returns>Dynamic pressure in Pascals (Pa).</returns>
    public static double Compute(double density, double velocity)
    {
        return 0.5 * density * velocity * velocity;
    }

    /// <summary>
    /// Computes dynamic pressure from altitude and velocity using ISA atmosphere.
    /// </summary>
    /// <param name="altitude">Altitude in metres.</param>
    /// <param name="velocity">True airspeed in m/s.</param>
    /// <returns>Dynamic pressure in Pascals.</returns>
    public static double FromAltitudeAndVelocity(double altitude, double velocity)
    {
        double rho = AtmosphereModel.Density(altitude);
        return 0.5 * rho * velocity * velocity;
    }

    /// <summary>
    /// Estimates whether Max-Q has been passed by checking the sign of dq/dt.
    /// Max-Q occurs when dq/dt crosses zero from positive to negative.
    /// </summary>
    /// <param name="qPrevious">Dynamic pressure at previous time step (Pa).</param>
    /// <param name="qCurrent">Dynamic pressure at current time step (Pa).</param>
    /// <returns>True if Max-Q was just passed (q is now decreasing).</returns>
    public static bool IsMaxQPassed(double qPrevious, double qCurrent)
    {
        return qCurrent < qPrevious && qPrevious > 0;
    }
}
