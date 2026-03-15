using CSharpNumerics.Physics.Materials.Nuclear.Isotopes;
using System;

namespace CSharpNumerics.Physics.Materials.Nuclear.Decay
{
    /// <summary>
    /// Static helpers for radioactive decay calculations:
    /// activity, remaining mass, and decay constant.
    /// </summary>
    public static class Decay
    {
        /// <summary>
        /// Returns the activity (Bq) of a sample of the given mass.
        /// A = m × specificActivity.
        /// </summary>
        /// <param name="massKg">Mass of the radioactive material (kg).</param>
        /// <param name="isotope">The isotope.</param>
        public static double Activity(double massKg, Isotope isotope)
        {
            if (isotope.IsStable) return 0;
            return massKg * isotope.SpecificActivity;
        }

        /// <summary>
        /// Returns the activity (Bq) at a given time after an initial activity A₀.
        /// A(t) = A₀ × exp(-λt).
        /// </summary>
        public static double ActivityAtTime(double initialActivityBq, double timeSeconds, Isotope isotope)
        {
            if (isotope.IsStable) return 0;
            return initialActivityBq * Math.Exp(-isotope.Lambda * timeSeconds);
        }

        /// <summary>
        /// Returns the remaining mass (kg) after radioactive decay.
        /// m(t) = m₀ × exp(-λt).
        /// </summary>
        /// <param name="initialMassKg">Initial mass (kg).</param>
        /// <param name="timeSeconds">Elapsed time (seconds).</param>
        /// <param name="isotope">The isotope.</param>
        public static double RemainingMass(double initialMassKg, double timeSeconds, Isotope isotope)
        {
            if (isotope.IsStable) return initialMassKg;
            return initialMassKg * Math.Exp(-isotope.Lambda * timeSeconds);
        }

        /// <summary>
        /// Decay constant λ = ln(2) / t½  (s⁻¹).
        /// </summary>
        public static double DecayConstant(Isotope isotope) => isotope.Lambda;

        /// <summary>
        /// Converts mass concentration (kg/m³) to activity concentration (Bq/m³).
        /// </summary>
        public static double ConcentrationToActivity(double concentrationKgM3, Isotope isotope)
        {
            if (isotope.IsStable) return 0;
            return concentrationKgM3 * isotope.SpecificActivity;
        }

        /// <summary>
        /// Converts activity concentration (Bq/m³) to mass concentration (kg/m³).
        /// </summary>
        public static double ActivityToConcentration(double activityBqM3, Isotope isotope)
        {
            if (isotope.IsStable || isotope.SpecificActivity <= 0) return 0;
            return activityBqM3 / isotope.SpecificActivity;
        }

        /// <summary>
        /// Number of half-lives elapsed in the given time.
        /// </summary>
        public static double HalfLivesElapsed(double timeSeconds, Isotope isotope)
        {
            if (isotope.IsStable) return 0;
            return timeSeconds / isotope.HalfLife;
        }
    }
}
