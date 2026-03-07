using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics;
using System;

namespace CSharpNumerics.Physics
{
    /// <summary>
    /// Provides extension methods for heat transfer calculations:
    /// Fourier conduction, Newton's law of cooling, Stefan–Boltzmann radiation,
    /// the heat equation, thermal resistance, dimensionless numbers,
    /// and lumped-capacitance transient analysis.
    /// <para>
    /// Bridges <see cref="ScalarField"/> (∇, ∇²) and
    /// <see cref="VectorFieldExtensions"/> to thermal physics.
    /// </para>
    /// </summary>
    public static class HeatExtensions
    {
        // ═══════════════════════════════════════════════════════════════
        //  Conduction — Fourier's law
        // ═══════════════════════════════════════════════════════════════

        #region Conduction

        /// <summary>
        /// Fourier's law (scalar, 1-D): computes the conductive heat flux
        /// through a slab: q = k · ΔT / L.
        /// </summary>
        /// <param name="thermalConductivity">Thermal conductivity k in W/(m·K).</param>
        /// <param name="temperatureDifference">Temperature difference ΔT = T_hot − T_cold in K.</param>
        /// <param name="thickness">Thickness L of the slab in metres.</param>
        /// <returns>Heat flux in W/m².</returns>
        public static double ConductiveHeatFlux(
            this double thermalConductivity,
            double temperatureDifference,
            double thickness)
        {
            if (thickness <= 0) throw new ArgumentException("Thickness must be greater than zero.");
            return thermalConductivity * temperatureDifference / thickness;
        }

        /// <summary>
        /// Fourier's law (scalar, 1-D): computes the conductive heat transfer rate
        /// through a slab of area A: Q̇ = k·A·ΔT / L.
        /// </summary>
        /// <param name="thermalConductivity">Thermal conductivity k in W/(m·K).</param>
        /// <param name="area">Cross-sectional area A in m².</param>
        /// <param name="temperatureDifference">Temperature difference ΔT in K.</param>
        /// <param name="thickness">Thickness L in metres.</param>
        /// <returns>Heat transfer rate in watts.</returns>
        public static double ConductiveHeatRate(
            this double thermalConductivity,
            double area,
            double temperatureDifference,
            double thickness)
        {
            if (thickness <= 0) throw new ArgumentException("Thickness must be greater than zero.");
            if (area <= 0) throw new ArgumentException("Area must be greater than zero.");
            return thermalConductivity * area * temperatureDifference / thickness;
        }

        /// <summary>
        /// Fourier's law (vector, 3-D): computes the heat flux vector field
        /// q = −k ∇T from a temperature <see cref="ScalarField"/>.
        /// </summary>
        /// <param name="temperature">Temperature field T(r) as a ScalarField.</param>
        /// <param name="thermalConductivity">Thermal conductivity k in W/(m·K).</param>
        /// <returns>Heat flux vector field q(r) in W/m².</returns>
        public static VectorField ConductiveHeatFlux(
            this ScalarField temperature,
            double thermalConductivity)
        {
            var grad = temperature.GradientField();
            return new VectorField(
                r => -thermalConductivity * grad.fx(r),
                r => -thermalConductivity * grad.fy(r),
                r => -thermalConductivity * grad.fz(r));
        }

        /// <summary>
        /// Conductive thermal resistance for a flat slab: R = L / (k·A).
        /// </summary>
        /// <param name="thickness">Slab thickness L in metres.</param>
        /// <param name="thermalConductivity">Thermal conductivity k in W/(m·K).</param>
        /// <param name="area">Cross-sectional area A in m².</param>
        /// <returns>Thermal resistance in K/W.</returns>
        public static double SlabThermalResistance(
            this double thickness,
            double thermalConductivity,
            double area)
        {
            if (thermalConductivity <= 0) throw new ArgumentException("Thermal conductivity must be greater than zero.");
            if (area <= 0) throw new ArgumentException("Area must be greater than zero.");
            return thickness / (thermalConductivity * area);
        }

        /// <summary>
        /// Conductive thermal resistance for a cylindrical shell:
        /// R = ln(r₂/r₁) / (2π·k·L).
        /// </summary>
        /// <param name="innerRadius">Inner radius r₁ in metres.</param>
        /// <param name="outerRadius">Outer radius r₂ in metres.</param>
        /// <param name="thermalConductivity">Thermal conductivity k in W/(m·K).</param>
        /// <param name="length">Cylinder length L in metres.</param>
        /// <returns>Thermal resistance in K/W.</returns>
        public static double CylindricalThermalResistance(
            this double innerRadius,
            double outerRadius,
            double thermalConductivity,
            double length)
        {
            if (innerRadius <= 0) throw new ArgumentException("Inner radius must be greater than zero.");
            if (outerRadius <= innerRadius) throw new ArgumentException("Outer radius must be greater than inner radius.");
            if (thermalConductivity <= 0) throw new ArgumentException("Thermal conductivity must be greater than zero.");
            if (length <= 0) throw new ArgumentException("Length must be greater than zero.");
            return Math.Log(outerRadius / innerRadius) / (2 * Math.PI * thermalConductivity * length);
        }

        /// <summary>
        /// Conductive thermal resistance for a spherical shell:
        /// R = (1/r₁ − 1/r₂) / (4π·k).
        /// </summary>
        /// <param name="innerRadius">Inner radius r₁ in metres.</param>
        /// <param name="outerRadius">Outer radius r₂ in metres.</param>
        /// <param name="thermalConductivity">Thermal conductivity k in W/(m·K).</param>
        /// <returns>Thermal resistance in K/W.</returns>
        public static double SphericalThermalResistance(
            this double innerRadius,
            double outerRadius,
            double thermalConductivity)
        {
            if (innerRadius <= 0) throw new ArgumentException("Inner radius must be greater than zero.");
            if (outerRadius <= innerRadius) throw new ArgumentException("Outer radius must be greater than inner radius.");
            if (thermalConductivity <= 0) throw new ArgumentException("Thermal conductivity must be greater than zero.");
            return (1.0 / innerRadius - 1.0 / outerRadius) / (4 * Math.PI * thermalConductivity);
        }

        /// <summary>
        /// Total thermal resistance for resistances in series:
        /// R_total = R₁ + R₂ + … + Rₙ.
        /// </summary>
        /// <param name="resistances">Array of thermal resistances in K/W.</param>
        /// <returns>Total series thermal resistance in K/W.</returns>
        public static double SeriesThermalResistance(this double[] resistances)
        {
            double sum = 0;
            for (int i = 0; i < resistances.Length; i++)
                sum += resistances[i];
            return sum;
        }

        /// <summary>
        /// Total thermal resistance for resistances in parallel:
        /// 1/R_total = 1/R₁ + 1/R₂ + … + 1/Rₙ.
        /// </summary>
        /// <param name="resistances">Array of thermal resistances in K/W.</param>
        /// <returns>Total parallel thermal resistance in K/W.</returns>
        public static double ParallelThermalResistance(this double[] resistances)
        {
            double sum = 0;
            for (int i = 0; i < resistances.Length; i++)
            {
                if (resistances[i] <= 0) throw new ArgumentException("Each resistance must be greater than zero.");
                sum += 1.0 / resistances[i];
            }
            return 1.0 / sum;
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Convection — Newton's law of cooling
        // ═══════════════════════════════════════════════════════════════

        #region Convection

        /// <summary>
        /// Newton's law of cooling: heat flux q = h·(T_s − T_∞).
        /// </summary>
        /// <param name="heatTransferCoefficient">Convective heat transfer coefficient h in W/(m²·K).</param>
        /// <param name="surfaceTemperature">Surface temperature T_s in K.</param>
        /// <param name="fluidTemperature">Bulk fluid temperature T_∞ in K.</param>
        /// <returns>Heat flux in W/m².</returns>
        public static double ConvectiveHeatFlux(
            this double heatTransferCoefficient,
            double surfaceTemperature,
            double fluidTemperature)
        {
            return heatTransferCoefficient * (surfaceTemperature - fluidTemperature);
        }

        /// <summary>
        /// Convective heat transfer rate: Q̇ = h·A·(T_s − T_∞).
        /// </summary>
        /// <param name="heatTransferCoefficient">Convective heat transfer coefficient h in W/(m²·K).</param>
        /// <param name="area">Surface area A in m².</param>
        /// <param name="surfaceTemperature">Surface temperature T_s in K.</param>
        /// <param name="fluidTemperature">Bulk fluid temperature T_∞ in K.</param>
        /// <returns>Heat transfer rate in watts.</returns>
        public static double ConvectiveHeatRate(
            this double heatTransferCoefficient,
            double area,
            double surfaceTemperature,
            double fluidTemperature)
        {
            if (area <= 0) throw new ArgumentException("Area must be greater than zero.");
            return heatTransferCoefficient * area * (surfaceTemperature - fluidTemperature);
        }

        /// <summary>
        /// Convective thermal resistance: R = 1 / (h·A).
        /// </summary>
        /// <param name="heatTransferCoefficient">Convective heat transfer coefficient h in W/(m²·K).</param>
        /// <param name="area">Surface area A in m².</param>
        /// <returns>Thermal resistance in K/W.</returns>
        public static double ConvectiveThermalResistance(
            this double heatTransferCoefficient,
            double area)
        {
            if (heatTransferCoefficient <= 0) throw new ArgumentException("Heat transfer coefficient must be greater than zero.");
            if (area <= 0) throw new ArgumentException("Area must be greater than zero.");
            return 1.0 / (heatTransferCoefficient * area);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Radiation — Stefan–Boltzmann law
        // ═══════════════════════════════════════════════════════════════

        #region Radiation

        /// <summary>
        /// Stefan–Boltzmann law: emissive power from a surface.
        /// q = ε·σ·T⁴.
        /// </summary>
        /// <param name="surfaceTemperature">Absolute surface temperature T in K.</param>
        /// <param name="emissivity">Surface emissivity ε (0–1, default 1 = blackbody).</param>
        /// <returns>Emissive power in W/m².</returns>
        public static double RadiativeHeatFlux(
            this double surfaceTemperature,
            double emissivity = 1.0)
        {
            if (surfaceTemperature < 0) throw new ArgumentException("Absolute temperature must be non-negative.");
            if (emissivity < 0 || emissivity > 1) throw new ArgumentException("Emissivity must be between 0 and 1.");
            double T4 = surfaceTemperature * surfaceTemperature * surfaceTemperature * surfaceTemperature;
            return emissivity * PhysicsConstants.StefanBoltzmannConstant * T4;
        }

        /// <summary>
        /// Net radiative heat flux between a surface and its surroundings:
        /// q_net = ε·σ·(T_s⁴ − T_sur⁴).
        /// </summary>
        /// <param name="surfaceTemperature">Surface temperature T_s in K.</param>
        /// <param name="surroundingTemperature">Surrounding temperature T_sur in K.</param>
        /// <param name="emissivity">Surface emissivity ε (0–1, default 1).</param>
        /// <returns>Net radiative heat flux in W/m².</returns>
        public static double NetRadiativeHeatFlux(
            this double surfaceTemperature,
            double surroundingTemperature,
            double emissivity = 1.0)
        {
            if (surfaceTemperature < 0 || surroundingTemperature < 0)
                throw new ArgumentException("Absolute temperatures must be non-negative.");
            if (emissivity < 0 || emissivity > 1) throw new ArgumentException("Emissivity must be between 0 and 1.");
            double Ts4 = surfaceTemperature * surfaceTemperature * surfaceTemperature * surfaceTemperature;
            double Tsur4 = surroundingTemperature * surroundingTemperature * surroundingTemperature * surroundingTemperature;
            return emissivity * PhysicsConstants.StefanBoltzmannConstant * (Ts4 - Tsur4);
        }

        /// <summary>
        /// Radiative heat transfer rate between a surface and its surroundings:
        /// Q̇ = ε·σ·A·(T_s⁴ − T_sur⁴).
        /// </summary>
        /// <param name="surfaceTemperature">Surface temperature T_s in K.</param>
        /// <param name="surroundingTemperature">Surrounding temperature T_sur in K.</param>
        /// <param name="area">Surface area A in m².</param>
        /// <param name="emissivity">Surface emissivity ε (0–1, default 1).</param>
        /// <returns>Radiative heat transfer rate in watts.</returns>
        public static double RadiativeHeatRate(
            this double surfaceTemperature,
            double surroundingTemperature,
            double area,
            double emissivity = 1.0)
        {
            if (area <= 0) throw new ArgumentException("Area must be greater than zero.");
            return area * surfaceTemperature.NetRadiativeHeatFlux(surroundingTemperature, emissivity);
        }

        /// <summary>
        /// Radiative thermal resistance (linearised form):
        /// R_rad = 1 / (h_rad · A), where h_rad = ε·σ·(T_s² + T_sur²)(T_s + T_sur).
        /// Valid when ΔT is moderate relative to the mean temperature.
        /// </summary>
        /// <param name="surfaceTemperature">Surface temperature T_s in K.</param>
        /// <param name="surroundingTemperature">Surrounding temperature T_sur in K.</param>
        /// <param name="area">Surface area A in m².</param>
        /// <param name="emissivity">Surface emissivity ε (0–1, default 1).</param>
        /// <returns>Thermal resistance in K/W.</returns>
        public static double RadiativeThermalResistance(
            this double surfaceTemperature,
            double surroundingTemperature,
            double area,
            double emissivity = 1.0)
        {
            if (area <= 0) throw new ArgumentException("Area must be greater than zero.");
            double Ts = surfaceTemperature;
            double Tsur = surroundingTemperature;
            double hRad = emissivity * PhysicsConstants.StefanBoltzmannConstant
                        * (Ts * Ts + Tsur * Tsur) * (Ts + Tsur);
            if (hRad <= 0) throw new ArgumentException("Linearised radiative coefficient is non-positive; temperatures may be invalid.");
            return 1.0 / (hRad * area);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Heat equation — ∂T/∂t = α ∇²T
        // ═══════════════════════════════════════════════════════════════

        #region Heat Equation

        /// <summary>
        /// Computes the thermal diffusivity: α = k / (ρ·c_p).
        /// </summary>
        /// <param name="thermalConductivity">Thermal conductivity k in W/(m·K).</param>
        /// <param name="density">Density ρ in kg/m³.</param>
        /// <param name="specificHeat">Specific heat capacity c_p in J/(kg·K).</param>
        /// <returns>Thermal diffusivity in m²/s.</returns>
        public static double ThermalDiffusivity(
            this double thermalConductivity,
            double density,
            double specificHeat)
        {
            if (density <= 0) throw new ArgumentException("Density must be greater than zero.");
            if (specificHeat <= 0) throw new ArgumentException("Specific heat must be greater than zero.");
            return thermalConductivity / (density * specificHeat);
        }

        /// <summary>
        /// Evaluates the heat equation: ∂T/∂t = α ∇²T.
        /// Returns a <see cref="ScalarField"/> representing the local rate of
        /// temperature change due to conduction.
        /// </summary>
        /// <param name="temperature">Temperature field T(r) as a ScalarField.</param>
        /// <param name="thermalDiffusivity">Thermal diffusivity α in m²/s.</param>
        /// <returns>∂T/∂t as a ScalarField.</returns>
        public static ScalarField HeatEquationRate(
            this ScalarField temperature,
            double thermalDiffusivity)
        {
            var func = temperature.f;
            return new ScalarField(r =>
                thermalDiffusivity * func.Laplacian((r.x, r.y, r.z)));
        }

        /// <summary>
        /// Evaluates the heat equation with a volumetric source term:
        /// ∂T/∂t = α ∇²T + q̇/(ρ·c_p).
        /// </summary>
        /// <param name="temperature">Temperature field T(r) as a ScalarField.</param>
        /// <param name="thermalDiffusivity">Thermal diffusivity α in m²/s.</param>
        /// <param name="volumetricSource">Volumetric heat source q̇(r) in W/m³.</param>
        /// <param name="density">Density ρ in kg/m³.</param>
        /// <param name="specificHeat">Specific heat capacity c_p in J/(kg·K).</param>
        /// <returns>∂T/∂t as a ScalarField.</returns>
        public static ScalarField HeatEquationRate(
            this ScalarField temperature,
            double thermalDiffusivity,
            ScalarField volumetricSource,
            double density,
            double specificHeat)
        {
            if (density <= 0) throw new ArgumentException("Density must be greater than zero.");
            if (specificHeat <= 0) throw new ArgumentException("Specific heat must be greater than zero.");
            var func = temperature.f;
            var srcFunc = volumetricSource.f;
            double rhoCp = density * specificHeat;
            return new ScalarField(r =>
                thermalDiffusivity * func.Laplacian((r.x, r.y, r.z)) + srcFunc(r) / rhoCp);
        }

        /// <summary>
        /// Analytical solution for 1-D heat conduction in a semi-infinite solid
        /// with surface held at T_s:
        /// T(x,t) = T_i + (T_s − T_i) · erfc(x / (2√(αt))).
        /// </summary>
        /// <param name="initialTemperature">Initial uniform temperature T_i in K.</param>
        /// <param name="surfaceTemperature">Surface temperature T_s in K (applied at x = 0).</param>
        /// <param name="thermalDiffusivity">Thermal diffusivity α in m²/s.</param>
        /// <param name="depth">Depth x into the solid in metres.</param>
        /// <param name="time">Time t since the surface temperature was applied in seconds.</param>
        /// <returns>Temperature at depth x and time t in K.</returns>
        public static double SemiInfiniteTemperature(
            this double initialTemperature,
            double surfaceTemperature,
            double thermalDiffusivity,
            double depth,
            double time)
        {
            if (time <= 0) throw new ArgumentException("Time must be greater than zero.");
            if (thermalDiffusivity <= 0) throw new ArgumentException("Thermal diffusivity must be greater than zero.");
            if (depth < 0) throw new ArgumentException("Depth must be non-negative.");
            double eta = depth / (2 * Math.Sqrt(thermalDiffusivity * time));
            return initialTemperature + (surfaceTemperature - initialTemperature) * Erfc(eta);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Dimensionless numbers
        // ═══════════════════════════════════════════════════════════════

        #region Dimensionless Numbers

        /// <summary>
        /// Biot number: Bi = h·L_c / k.
        /// <para>
        /// Bi ≪ 1 → lumped-capacitance model is valid (uniform internal temperature).
        /// Bi ≫ 1 → significant internal temperature gradients.
        /// </para>
        /// </summary>
        /// <param name="heatTransferCoefficient">Convective coefficient h in W/(m²·K).</param>
        /// <param name="characteristicLength">Characteristic length L_c in metres (volume / surface area).</param>
        /// <param name="thermalConductivity">Solid thermal conductivity k in W/(m·K).</param>
        public static double BiotNumber(
            this double heatTransferCoefficient,
            double characteristicLength,
            double thermalConductivity)
        {
            if (thermalConductivity <= 0) throw new ArgumentException("Thermal conductivity must be greater than zero.");
            return heatTransferCoefficient * characteristicLength / thermalConductivity;
        }

        /// <summary>
        /// Nusselt number: Nu = h·L / k_fluid.
        /// Represents the ratio of convective to conductive heat transfer
        /// across the boundary layer.
        /// </summary>
        /// <param name="heatTransferCoefficient">Convective coefficient h in W/(m²·K).</param>
        /// <param name="characteristicLength">Characteristic length L in metres.</param>
        /// <param name="fluidConductivity">Thermal conductivity of the fluid k_f in W/(m·K).</param>
        public static double NusseltNumber(
            this double heatTransferCoefficient,
            double characteristicLength,
            double fluidConductivity)
        {
            if (fluidConductivity <= 0) throw new ArgumentException("Fluid conductivity must be greater than zero.");
            return heatTransferCoefficient * characteristicLength / fluidConductivity;
        }

        /// <summary>
        /// Fourier number: Fo = α·t / L².
        /// Represents the ratio of conductive heat transfer rate to
        /// the rate of thermal energy storage.
        /// </summary>
        /// <param name="thermalDiffusivity">Thermal diffusivity α in m²/s.</param>
        /// <param name="time">Time t in seconds.</param>
        /// <param name="characteristicLength">Characteristic length L in metres.</param>
        public static double FourierNumber(
            this double thermalDiffusivity,
            double time,
            double characteristicLength)
        {
            if (characteristicLength <= 0) throw new ArgumentException("Characteristic length must be greater than zero.");
            return thermalDiffusivity * time / (characteristicLength * characteristicLength);
        }

        /// <summary>
        /// Prandtl number: Pr = ν / α = c_p·μ / k.
        /// Represents the ratio of momentum diffusivity to thermal diffusivity.
        /// </summary>
        /// <param name="specificHeat">Specific heat capacity c_p in J/(kg·K).</param>
        /// <param name="dynamicViscosity">Dynamic viscosity μ in Pa·s.</param>
        /// <param name="thermalConductivity">Thermal conductivity k in W/(m·K).</param>
        public static double PrandtlNumber(
            this double specificHeat,
            double dynamicViscosity,
            double thermalConductivity)
        {
            if (thermalConductivity <= 0) throw new ArgumentException("Thermal conductivity must be greater than zero.");
            return specificHeat * dynamicViscosity / thermalConductivity;
        }

        /// <summary>
        /// Grashof number: Gr = g·β·ΔT·L³ / ν².
        /// Represents the ratio of buoyancy to viscous forces in natural convection.
        /// </summary>
        /// <param name="temperatureDifference">Temperature difference ΔT = T_s − T_∞ in K.</param>
        /// <param name="characteristicLength">Characteristic length L in metres.</param>
        /// <param name="thermalExpansionCoefficient">Volumetric thermal expansion coefficient β in 1/K.</param>
        /// <param name="kinematicViscosity">Kinematic viscosity ν in m²/s.</param>
        public static double GrashofNumber(
            this double temperatureDifference,
            double characteristicLength,
            double thermalExpansionCoefficient,
            double kinematicViscosity)
        {
            if (kinematicViscosity <= 0) throw new ArgumentException("Kinematic viscosity must be greater than zero.");
            double L3 = characteristicLength * characteristicLength * characteristicLength;
            return PhysicsConstants.GravitationalAcceleration * thermalExpansionCoefficient
                 * Math.Abs(temperatureDifference) * L3 / (kinematicViscosity * kinematicViscosity);
        }

        /// <summary>
        /// Rayleigh number: Ra = Gr · Pr = g·β·ΔT·L³ / (ν·α).
        /// <para>
        /// Ra &lt; 10⁹ → laminar natural convection.
        /// Ra &gt; 10⁹ → turbulent natural convection.
        /// </para>
        /// </summary>
        /// <param name="temperatureDifference">Temperature difference ΔT in K.</param>
        /// <param name="characteristicLength">Characteristic length L in metres.</param>
        /// <param name="thermalExpansionCoefficient">Volumetric thermal expansion coefficient β in 1/K.</param>
        /// <param name="kinematicViscosity">Kinematic viscosity ν in m²/s.</param>
        /// <param name="thermalDiffusivity">Thermal diffusivity α in m²/s.</param>
        public static double RayleighNumber(
            this double temperatureDifference,
            double characteristicLength,
            double thermalExpansionCoefficient,
            double kinematicViscosity,
            double thermalDiffusivity)
        {
            if (kinematicViscosity <= 0) throw new ArgumentException("Kinematic viscosity must be greater than zero.");
            if (thermalDiffusivity <= 0) throw new ArgumentException("Thermal diffusivity must be greater than zero.");
            double L3 = characteristicLength * characteristicLength * characteristicLength;
            return PhysicsConstants.GravitationalAcceleration * thermalExpansionCoefficient
                 * Math.Abs(temperatureDifference) * L3 / (kinematicViscosity * thermalDiffusivity);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Lumped-capacitance transient analysis
        // ═══════════════════════════════════════════════════════════════

        #region Lumped Capacitance

        /// <summary>
        /// Lumped-capacitance model: temperature of a body at time t.
        /// T(t) = T_∞ + (T₀ − T_∞) · exp(−t / τ),
        /// where τ = ρ·V·c_p / (h·A).
        /// <para>
        /// Valid when Bi ≪ 1 (typically Bi &lt; 0.1).
        /// </para>
        /// </summary>
        /// <param name="initialTemperature">Initial body temperature T₀ in K.</param>
        /// <param name="ambientTemperature">Ambient temperature T_∞ in K.</param>
        /// <param name="time">Elapsed time t in seconds.</param>
        /// <param name="timeConstant">Thermal time constant τ = ρVc_p / (hA) in seconds.</param>
        /// <returns>Body temperature at time t in K.</returns>
        public static double LumpedCapacitanceTemperature(
            this double initialTemperature,
            double ambientTemperature,
            double time,
            double timeConstant)
        {
            if (timeConstant <= 0) throw new ArgumentException("Time constant must be greater than zero.");
            return ambientTemperature + (initialTemperature - ambientTemperature) * Math.Exp(-time / timeConstant);
        }

        /// <summary>
        /// Computes the thermal time constant for the lumped-capacitance model:
        /// τ = ρ·V·c_p / (h·A).
        /// </summary>
        /// <param name="density">Density ρ in kg/m³.</param>
        /// <param name="volume">Volume V in m³.</param>
        /// <param name="specificHeat">Specific heat capacity c_p in J/(kg·K).</param>
        /// <param name="heatTransferCoefficient">Convective coefficient h in W/(m²·K).</param>
        /// <param name="surfaceArea">Surface area A in m².</param>
        /// <returns>Time constant τ in seconds.</returns>
        public static double ThermalTimeConstant(
            this double density,
            double volume,
            double specificHeat,
            double heatTransferCoefficient,
            double surfaceArea)
        {
            if (heatTransferCoefficient <= 0) throw new ArgumentException("Heat transfer coefficient must be greater than zero.");
            if (surfaceArea <= 0) throw new ArgumentException("Surface area must be greater than zero.");
            return density * volume * specificHeat / (heatTransferCoefficient * surfaceArea);
        }

        /// <summary>
        /// Computes the time required for a body to reach a target temperature
        /// under the lumped-capacitance model:
        /// t = −τ · ln((T_target − T_∞) / (T₀ − T_∞)).
        /// </summary>
        /// <param name="initialTemperature">Initial body temperature T₀ in K.</param>
        /// <param name="targetTemperature">Target temperature T_target in K.</param>
        /// <param name="ambientTemperature">Ambient temperature T_∞ in K.</param>
        /// <param name="timeConstant">Thermal time constant τ in seconds.</param>
        /// <returns>Time in seconds to reach the target temperature.</returns>
        public static double LumpedCapacitanceTime(
            this double initialTemperature,
            double targetTemperature,
            double ambientTemperature,
            double timeConstant)
        {
            if (timeConstant <= 0) throw new ArgumentException("Time constant must be greater than zero.");
            double ratio = (targetTemperature - ambientTemperature) / (initialTemperature - ambientTemperature);
            if (ratio <= 0) throw new ArgumentException("Target temperature is not between initial and ambient temperatures.");
            return -timeConstant * Math.Log(ratio);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Fin heat transfer
        // ═══════════════════════════════════════════════════════════════

        #region Fins

        /// <summary>
        /// Computes the fin parameter: m = √(h·P / (k·A_c)),
        /// where P is the fin perimeter and A_c is the cross-sectional area.
        /// </summary>
        /// <param name="heatTransferCoefficient">Convective coefficient h in W/(m²·K).</param>
        /// <param name="perimeter">Fin perimeter P in metres.</param>
        /// <param name="thermalConductivity">Fin thermal conductivity k in W/(m·K).</param>
        /// <param name="crossSectionalArea">Fin cross-sectional area A_c in m².</param>
        /// <returns>Fin parameter m in 1/m.</returns>
        public static double FinParameter(
            this double heatTransferCoefficient,
            double perimeter,
            double thermalConductivity,
            double crossSectionalArea)
        {
            if (thermalConductivity <= 0) throw new ArgumentException("Thermal conductivity must be greater than zero.");
            if (crossSectionalArea <= 0) throw new ArgumentException("Cross-sectional area must be greater than zero.");
            return Math.Sqrt(heatTransferCoefficient * perimeter / (thermalConductivity * crossSectionalArea));
        }

        /// <summary>
        /// Heat transfer from an infinitely long fin:
        /// Q̇ = √(h·P·k·A_c) · (T_b − T_∞).
        /// </summary>
        /// <param name="baseTemperature">Fin base temperature T_b in K.</param>
        /// <param name="ambientTemperature">Ambient temperature T_∞ in K.</param>
        /// <param name="heatTransferCoefficient">Convective coefficient h in W/(m²·K).</param>
        /// <param name="perimeter">Fin perimeter P in metres.</param>
        /// <param name="thermalConductivity">Fin thermal conductivity k in W/(m·K).</param>
        /// <param name="crossSectionalArea">Fin cross-sectional area A_c in m².</param>
        /// <returns>Heat transfer rate from the fin in watts.</returns>
        public static double InfiniteFinHeatRate(
            this double baseTemperature,
            double ambientTemperature,
            double heatTransferCoefficient,
            double perimeter,
            double thermalConductivity,
            double crossSectionalArea)
        {
            double m = heatTransferCoefficient.FinParameter(perimeter, thermalConductivity, crossSectionalArea);
            return m * thermalConductivity * crossSectionalArea * (baseTemperature - ambientTemperature);
        }

        /// <summary>
        /// Heat transfer from a finite-length fin with an insulated tip:
        /// Q̇ = √(h·P·k·A_c) · (T_b − T_∞) · tanh(m·L).
        /// </summary>
        /// <param name="baseTemperature">Fin base temperature T_b in K.</param>
        /// <param name="ambientTemperature">Ambient temperature T_∞ in K.</param>
        /// <param name="heatTransferCoefficient">Convective coefficient h in W/(m²·K).</param>
        /// <param name="perimeter">Fin perimeter P in metres.</param>
        /// <param name="thermalConductivity">Fin thermal conductivity k in W/(m·K).</param>
        /// <param name="crossSectionalArea">Fin cross-sectional area A_c in m².</param>
        /// <param name="finLength">Fin length L in metres.</param>
        /// <returns>Heat transfer rate from the fin in watts.</returns>
        public static double InsulatedTipFinHeatRate(
            this double baseTemperature,
            double ambientTemperature,
            double heatTransferCoefficient,
            double perimeter,
            double thermalConductivity,
            double crossSectionalArea,
            double finLength)
        {
            double m = heatTransferCoefficient.FinParameter(perimeter, thermalConductivity, crossSectionalArea);
            return m * thermalConductivity * crossSectionalArea
                 * (baseTemperature - ambientTemperature) * Math.Tanh(m * finLength);
        }

        /// <summary>
        /// Fin efficiency for a fin with an insulated tip:
        /// η = tanh(m·L) / (m·L).
        /// </summary>
        /// <param name="finParameter">Fin parameter m in 1/m.</param>
        /// <param name="finLength">Fin length L in metres.</param>
        /// <returns>Fin efficiency (0–1).</returns>
        public static double InsulatedTipFinEfficiency(
            this double finParameter,
            double finLength)
        {
            double mL = finParameter * finLength;
            if (mL == 0) return 1.0;
            return Math.Tanh(mL) / mL;
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Internal helpers
        // ═══════════════════════════════════════════════════════════════

        #region Helpers

        /// <summary>
        /// Complementary error function erfc(x) = 1 − erf(x).
        /// Uses the Abramowitz and Stegun rational approximation (max error ≈ 1.5 × 10⁻⁷).
        /// </summary>
        private static double Erfc(double x)
        {
            return 1.0 - Erf(x);
        }

        /// <summary>
        /// Error function erf(x). Approximation from Abramowitz and Stegun (7.1.26).
        /// </summary>
        private static double Erf(double x)
        {
            bool negative = x < 0;
            x = Math.Abs(x);

            const double a1 = 0.254829592;
            const double a2 = -0.284496736;
            const double a3 = 1.421413741;
            const double a4 = -1.453152027;
            const double a5 = 1.061405429;
            const double p = 0.3275911;

            double t = 1.0 / (1.0 + p * x);
            double poly = t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))));
            double result = 1.0 - poly * Math.Exp(-x * x);

            return negative ? -result : result;
        }

        #endregion
    }
}
