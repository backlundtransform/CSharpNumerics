using CSharpNumerics.Physics.Materials.Nuclear.Isotopes;
using System;

namespace CSharpNumerics.Physics.Materials.Nuclear.RadiationDose
{
    /// <summary>
    /// Radiation dose calculations: external dose rate, ground-shine,
    /// and inhalation committed effective dose.
    /// Uses ICRP-72 dose coefficients stored on the <see cref="Isotope"/> struct.
    /// </summary>
    public static class RadiationDose
    {
        /// <summary>Default adult breathing rate (m³/s). ~1.2 m³/h = 3.33e-4 m³/s.</summary>
        public const double DefaultBreathingRate = 3.33e-4;

        /// <summary>
        /// External dose rate from a point source at a given distance.
        /// Simplified inverse-square model using gamma energy.
        /// Returns Sv/h.
        /// </summary>
        /// <param name="activityBq">Total source activity (Bq).</param>
        /// <param name="distanceM">Distance from source (metres).</param>
        /// <param name="isotope">The isotope (for gamma energy).</param>
        /// <returns>Dose rate in Sv/h. Zero for pure beta/alpha emitters.</returns>
        public static double DoseRate(double activityBq, double distanceM, Isotope isotope)
        {
            if (isotope.GammaEnergy <= 0) return 0;
            if (distanceM <= 0) throw new ArgumentException("Distance must be positive.", nameof(distanceM));

            // Simplified: Ḋ = Γ × A / d²
            // Gamma dose constant Γ ≈ 5.76e-14 × E_γ (Sv·m²·s / Bq·h) for a point source
            double gamma = 5.76e-14 * isotope.GammaEnergy; // Sv·m²/(Bq·h)
            return gamma * activityBq / (distanceM * distanceM);
        }

        /// <summary>
        /// Ground-shine committed effective dose from deposited radioactive material.
        /// D = depositionBq/m² × doseCoeff × exposureTime.
        /// </summary>
        /// <param name="depositionBqM2">Surface deposition (Bq/m²).</param>
        /// <param name="isotope">The isotope.</param>
        /// <param name="exposureSeconds">Duration of exposure (seconds).</param>
        /// <returns>Committed effective dose in Sv.</returns>
        public static double GroundShineDose(double depositionBqM2, Isotope isotope, double exposureSeconds)
        {
            if (isotope.GroundShineDoseCoeff <= 0) return 0;
            // Dose = deposition × coefficient × time
            // coefficient is in Sv·m²/(Bq·s)
            return depositionBqM2 * isotope.GroundShineDoseCoeff * exposureSeconds;
        }

        /// <summary>
        /// Inhalation committed effective dose.
        /// D = concentration (Bq/m³) × breathing rate (m³/s) × exposure time (s) × dose coefficient (Sv/Bq).
        /// </summary>
        /// <param name="concentrationBqM3">Airborne activity concentration (Bq/m³).</param>
        /// <param name="breathingRateM3s">Breathing rate (m³/s). Use <see cref="DefaultBreathingRate"/>.</param>
        /// <param name="exposureSeconds">Duration of exposure (seconds).</param>
        /// <param name="isotope">The isotope.</param>
        /// <returns>Committed effective dose in Sv.</returns>
        public static double InhalationDose(
            double concentrationBqM3,
            double breathingRateM3s,
            double exposureSeconds,
            Isotope isotope)
        {
            if (isotope.InhalationDoseCoeff <= 0) return 0;
            // Intake = concentration × breathing rate × time  [Bq]
            // Dose = Intake × dose coefficient  [Sv]
            double intake = concentrationBqM3 * breathingRateM3s * exposureSeconds;
            return intake * isotope.InhalationDoseCoeff;
        }

        /// <summary>
        /// Inhalation dose using the default adult breathing rate.
        /// </summary>
        public static double InhalationDose(
            double concentrationBqM3,
            double exposureSeconds,
            Isotope isotope)
            => InhalationDose(concentrationBqM3, DefaultBreathingRate, exposureSeconds, isotope);
    }
}
