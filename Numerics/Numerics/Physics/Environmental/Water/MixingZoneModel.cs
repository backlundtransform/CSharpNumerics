using System;

namespace CSharpNumerics.Physics.Environmental.Water;

/// <summary>
/// Mixing-zone calculations for tributary confluences and lateral mixing
/// in open channels.
/// <para>
/// All methods are pure functions. No grid or simulation state.
/// </para>
/// </summary>
public static class MixingZoneModel
{
    // ═══════════════════════════════════════════════════════════════
    //  Tributary mixing (mass balance)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Downstream concentration after a tributary joins the main channel,
    /// assuming instantaneous complete mixing (mass balance):
    /// <para>C_ds = (C_main · Q_main + C_trib · Q_trib) / (Q_main + Q_trib)</para>
    /// Returns 0 if total discharge is zero.
    /// </summary>
    /// <param name="mainConcentration">Main stem concentration C_main (mg/L or any consistent unit).</param>
    /// <param name="mainDischarge">Main stem volumetric discharge Q_main (m³/s).</param>
    /// <param name="tributaryConcentration">Tributary concentration C_trib.</param>
    /// <param name="tributaryDischarge">Tributary discharge Q_trib (m³/s).</param>
    public static double TributaryMixing(
        double mainConcentration,
        double mainDischarge,
        double tributaryConcentration,
        double tributaryDischarge)
    {
        double totalQ = mainDischarge + tributaryDischarge;
        if (totalQ <= 0) return 0;

        return (mainConcentration * mainDischarge + tributaryConcentration * tributaryDischarge) / totalQ;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Mixing length
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Mixing length L_mix (m) — the downstream distance required for
    /// complete lateral mixing after a side discharge:
    /// <para>L_mix ≈ 0.4 · u · W² / E_t</para>
    /// where u is mean velocity, W is channel width, and E_t is
    /// transverse dispersion coefficient. The factor 0.4 corresponds
    /// to ≈ 95 % mixing (Fischer et al. 1979).
    /// Returns 0 if any input is non-positive.
    /// </summary>
    /// <param name="velocity">Mean channel velocity u (m/s).</param>
    /// <param name="width">Channel width W (m).</param>
    /// <param name="transverseDispersion">Transverse dispersion coefficient E_t (m²/s).
    /// Typical: E_t ≈ 0.15–0.60 · H · u* for straight channels.</param>
    public static double MixingLength(double velocity, double width, double transverseDispersion)
    {
        if (velocity <= 0 || width <= 0 || transverseDispersion <= 0)
            return 0;

        return 0.4 * velocity * width * width / transverseDispersion;
    }
}
