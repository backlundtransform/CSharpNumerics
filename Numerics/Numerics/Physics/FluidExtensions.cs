using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Enums;
using CSharpNumerics.Numerics;
using System;

namespace CSharpNumerics.Physics
{
    /// <summary>
    /// Provides extension methods for fluid dynamics calculations:
    /// Navier–Stokes equations, Bernoulli's principle, continuity equation,
    /// vorticity, stream functions, drag/lift forces, and viscous flows.
    /// <para>
    /// Bridges <see cref="VectorField"/> (∇·, ∇×) and
    /// <see cref="ScalarField"/> (∇, ∇²) to classical fluid dynamics.
    /// </para>
    /// </summary>
    public static class FluidExtensions
    {
        // ═══════════════════════════════════════════════════════════════
        //  Navier–Stokes equations (incompressible)
        //
        //  ρ(∂v/∂t + (v·∇)v) = −∇p + μ∇²v + f
        //  ∇·v = 0  (incompressibility constraint)
        // ═══════════════════════════════════════════════════════════════

        #region Navier–Stokes

        /// <summary>
        /// Computes the convective acceleration term (v·∇)v of the Navier–Stokes
        /// equations at a point. This is the nonlinear advection term.
        /// </summary>
        /// <param name="velocity">The velocity vector field v(r).</param>
        /// <param name="point">The point at which to evaluate.</param>
        /// <param name="h">Step size for numerical differentiation (default 1×10⁻⁶).</param>
        public static Vector ConvectiveAcceleration(
            this VectorField velocity, (double, double, double) point, double h = 1e-6)
        {
            var p = new Vector(point);
            var vx = velocity.fx(p);
            var vy = velocity.fy(p);
            var vz = velocity.fz(p);

            // ∂vᵢ/∂xⱼ via central differences
            double Partial(Func<Vector, double> comp, int axis)
            {
                var pPlus = new Vector(
                    point.Item1 + (axis == 0 ? h : 0),
                    point.Item2 + (axis == 1 ? h : 0),
                    point.Item3 + (axis == 2 ? h : 0));
                var pMinus = new Vector(
                    point.Item1 - (axis == 0 ? h : 0),
                    point.Item2 - (axis == 1 ? h : 0),
                    point.Item3 - (axis == 2 ? h : 0));
                return (comp(pPlus) - comp(pMinus)) / (2.0 * h);
            }

            // (v·∇)vᵢ = vx ∂vᵢ/∂x + vy ∂vᵢ/∂y + vz ∂vᵢ/∂z
            double ax = vx * Partial(velocity.fx, 0) + vy * Partial(velocity.fx, 1) + vz * Partial(velocity.fx, 2);
            double ay = vx * Partial(velocity.fy, 0) + vy * Partial(velocity.fy, 1) + vz * Partial(velocity.fy, 2);
            double az = vx * Partial(velocity.fz, 0) + vy * Partial(velocity.fz, 1) + vz * Partial(velocity.fz, 2);

            return new Vector(ax, ay, az);
        }

        /// <summary>
        /// Computes the viscous diffusion term μ∇²v of the Navier–Stokes equations
        /// at a point, where ∇²v is the vector Laplacian applied component-wise.
        /// </summary>
        /// <param name="velocity">The velocity vector field v(r).</param>
        /// <param name="dynamicViscosity">Dynamic viscosity μ in Pa·s.</param>
        /// <param name="point">The point at which to evaluate.</param>
        /// <param name="h">Step size for numerical differentiation (default 1×10⁻⁶).</param>
        public static Vector ViscousTerm(
            this VectorField velocity, double dynamicViscosity,
            (double, double, double) point, double h = 1e-6)
        {
            double LaplacianComponent(Func<Vector, double> comp)
            {
                return new ScalarField(comp).Laplacian(point);
            }

            return dynamicViscosity * new Vector(
                LaplacianComponent(velocity.fx),
                LaplacianComponent(velocity.fy),
                LaplacianComponent(velocity.fz));
        }

        /// <summary>
        /// Computes the pressure gradient term −∇p of the Navier–Stokes equations.
        /// </summary>
        /// <param name="pressure">The pressure scalar field p(r).</param>
        /// <param name="point">The point at which to evaluate.</param>
        public static Vector PressureGradientForce(
            this ScalarField pressure, (double, double, double) point)
        {
            var grad = pressure.Gradient(point);
            return new Vector(-grad.x, -grad.y, -grad.z);
        }

        /// <summary>
        /// Computes the full Navier–Stokes acceleration (right-hand side / ρ)
        /// for an incompressible Newtonian fluid at a steady state:
        /// a = (−∇p + μ∇²v + f) / ρ − (v·∇)v.
        /// <para>
        /// This returns what ∂v/∂t should equal if the flow is not in steady state.
        /// For steady-state flow, the return value should be approximately zero.
        /// </para>
        /// </summary>
        /// <param name="velocity">The velocity vector field v(r).</param>
        /// <param name="pressure">The pressure scalar field p(r).</param>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        /// <param name="dynamicViscosity">Dynamic viscosity μ in Pa·s.</param>
        /// <param name="point">The point at which to evaluate.</param>
        /// <param name="bodyForce">External body force per unit volume f (e.g. gravity). Pass zero vector if none.</param>
        public static Vector NavierStokesResidual(
            this VectorField velocity, ScalarField pressure,
            double density, double dynamicViscosity,
            (double, double, double) point, Vector bodyForce)
        {
            var pressureGrad = pressure.PressureGradientForce(point);
            var viscous = velocity.ViscousTerm(dynamicViscosity, point);
            var convective = velocity.ConvectiveAcceleration(point);

            // ∂v/∂t = (−∇p + μ∇²v + f) / ρ − (v·∇)v
            return (1.0 / density) * (pressureGrad + viscous + bodyForce) - convective;
        }

        /// <summary>
        /// Computes the full Navier–Stokes acceleration with no external body force.
        /// Overload of <see cref="NavierStokesResidual(VectorField, ScalarField, double, double, ValueTuple{double,double,double}, Vector)"/>.
        /// </summary>
        public static Vector NavierStokesResidual(
            this VectorField velocity, ScalarField pressure,
            double density, double dynamicViscosity,
            (double, double, double) point)
        {
            return velocity.NavierStokesResidual(pressure, density, dynamicViscosity, point, new Vector(0, 0, 0));
        }

        /// <summary>
        /// Checks the incompressibility constraint ∇·v = 0.
        /// Returns the divergence of the velocity field, which should be
        /// zero for an incompressible flow.
        /// </summary>
        /// <param name="velocity">The velocity vector field v(r).</param>
        /// <param name="point">The point at which to evaluate.</param>
        public static double IncompressibilityResidual(
            this VectorField velocity, (double, double, double) point)
        {
            return velocity.Divergence(point);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Euler equations (inviscid flow)
        //
        //  ρ(∂v/∂t + (v·∇)v) = −∇p + f
        // ═══════════════════════════════════════════════════════════════

        #region Euler Equations

        /// <summary>
        /// Computes the Euler equation residual for inviscid flow:
        /// ∂v/∂t = (−∇p + f) / ρ − (v·∇)v.
        /// Equivalent to Navier–Stokes with μ = 0.
        /// </summary>
        /// <param name="velocity">The velocity vector field v(r).</param>
        /// <param name="pressure">The pressure scalar field p(r).</param>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        /// <param name="point">The point at which to evaluate.</param>
        /// <param name="bodyForce">External body force per unit volume f. Pass zero vector if none.</param>
        public static Vector EulerEquationResidual(
            this VectorField velocity, ScalarField pressure,
            double density, (double, double, double) point, Vector bodyForce)
        {
            var pressureGrad = pressure.PressureGradientForce(point);
            var convective = velocity.ConvectiveAcceleration(point);

            return (1.0 / density) * (pressureGrad + bodyForce) - convective;
        }

        /// <summary>
        /// Computes the Euler equation residual with no external body force.
        /// Overload of <see cref="EulerEquationResidual(VectorField, ScalarField, double, ValueTuple{double,double,double}, Vector)"/>.
        /// </summary>
        public static Vector EulerEquationResidual(
            this VectorField velocity, ScalarField pressure,
            double density, (double, double, double) point)
        {
            return velocity.EulerEquationResidual(pressure, density, point, new Vector(0, 0, 0));
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Bernoulli's principle
        //
        //  p + ½ρv² + ρgh = const   (along a streamline, steady, inviscid)
        // ═══════════════════════════════════════════════════════════════

        #region Bernoulli

        /// <summary>
        /// Computes Bernoulli's constant along a streamline:
        /// B = p + ½ρv² + ρgh.
        /// For steady, inviscid, incompressible flow B is constant along streamlines.
        /// </summary>
        /// <param name="pressure">Static pressure at the point in Pa.</param>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        /// <param name="velocity">Flow speed at the point in m/s.</param>
        /// <param name="height">Height above reference datum in m.</param>
        public static double BernoulliConstant(
            this double pressure, double density, double velocity, double height = 0)
        {
            return pressure
                + 0.5 * density * velocity * velocity
                + density * PhysicsConstants.GravitationalAcceleration * height;
        }

        /// <summary>
        /// Computes the pressure at a second point along a streamline using
        /// Bernoulli's equation, given conditions at a first point.
        /// p₂ = p₁ + ½ρ(v₁² − v₂²) + ρg(h₁ − h₂).
        /// </summary>
        /// <param name="p1">Pressure at point 1 in Pa.</param>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        /// <param name="v1">Speed at point 1 in m/s.</param>
        /// <param name="v2">Speed at point 2 in m/s.</param>
        /// <param name="h1">Height at point 1 in m.</param>
        /// <param name="h2">Height at point 2 in m.</param>
        public static double BernoulliPressure(
            this double p1, double density, double v1, double v2,
            double h1 = 0, double h2 = 0)
        {
            return p1
                + 0.5 * density * (v1 * v1 - v2 * v2)
                + density * PhysicsConstants.GravitationalAcceleration * (h1 - h2);
        }

        /// <summary>
        /// Computes the dynamic pressure: q = ½ρv².
        /// </summary>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        /// <param name="speed">Flow speed in m/s.</param>
        public static double DynamicPressure(this double density, double speed)
        {
            return 0.5 * density * speed * speed;
        }

        /// <summary>
        /// Computes the stagnation (total) pressure: p₀ = p + ½ρv².
        /// </summary>
        /// <param name="staticPressure">Static pressure in Pa.</param>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        /// <param name="speed">Flow speed in m/s.</param>
        public static double StagnationPressure(this double staticPressure, double density, double speed)
        {
            return staticPressure + 0.5 * density * speed * speed;
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Continuity equation
        //
        //  ∂ρ/∂t + ∇·(ρv) = 0
        //  For incompressible flow: ∇·v = 0
        // ═══════════════════════════════════════════════════════════════

        #region Continuity

        /// <summary>
        /// Computes the mass flux ρv as a vector field from a constant-density
        /// velocity field.
        /// </summary>
        /// <param name="velocity">The velocity vector field v(r).</param>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        public static VectorField MassFlux(this VectorField velocity, double density)
        {
            return new VectorField(
                r => density * velocity.fx(r),
                r => density * velocity.fy(r),
                r => density * velocity.fz(r));
        }

        /// <summary>
        /// Computes the mass flux ρv as a vector field with spatially varying density.
        /// </summary>
        /// <param name="velocity">The velocity vector field v(r).</param>
        /// <param name="density">The density scalar field ρ(r).</param>
        public static VectorField MassFlux(this VectorField velocity, ScalarField density)
        {
            return new VectorField(
                r => density.f(r) * velocity.fx(r),
                r => density.f(r) * velocity.fy(r),
                r => density.f(r) * velocity.fz(r));
        }

        /// <summary>
        /// Evaluates the continuity equation residual ∇·(ρv) for constant density.
        /// For incompressible flow this equals ρ∇·v, which should be zero.
        /// </summary>
        /// <param name="velocity">The velocity vector field v(r).</param>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        /// <param name="point">The point at which to evaluate.</param>
        public static double ContinuityResidual(
            this VectorField velocity, double density, (double, double, double) point)
        {
            return velocity.MassFlux(density).Divergence(point);
        }

        /// <summary>
        /// Computes the volume flow rate Q = A · v through a cross-section
        /// of known area, given a uniform speed perpendicular to the section.
        /// </summary>
        /// <param name="area">Cross-sectional area in m².</param>
        /// <param name="speed">Average flow speed normal to the section in m/s.</param>
        public static double VolumeFlowRate(this double area, double speed)
        {
            return area * speed;
        }

        /// <summary>
        /// Computes the speed at a second cross-section from the continuity
        /// equation for incompressible flow: A₁v₁ = A₂v₂ → v₂ = A₁v₁ / A₂.
        /// </summary>
        /// <param name="area1">Cross-sectional area at section 1 in m².</param>
        /// <param name="speed1">Speed at section 1 in m/s.</param>
        /// <param name="area2">Cross-sectional area at section 2 in m².</param>
        public static double ContinuitySpeed(this double area1, double speed1, double area2)
        {
            if (area2 <= 0) throw new ArgumentException("Cross-sectional area must be greater than zero.");
            return area1 * speed1 / area2;
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Vorticity & circulation
        //
        //  ω = ∇ × v
        // ═══════════════════════════════════════════════════════════════

        #region Vorticity

        /// <summary>
        /// Computes the vorticity ω = ∇ × v at a point.
        /// </summary>
        /// <param name="velocity">The velocity vector field v(r).</param>
        /// <param name="point">The point at which to evaluate.</param>
        public static Vector Vorticity(
            this VectorField velocity, (double, double, double) point)
        {
            return velocity.Curl(point);
        }

        /// <summary>
        /// Computes the enstrophy density ½|ω|² at a point, a measure of
        /// rotational kinetic energy per unit volume (divided by ρ).
        /// </summary>
        /// <param name="velocity">The velocity vector field v(r).</param>
        /// <param name="point">The point at which to evaluate.</param>
        public static double Enstrophy(
            this VectorField velocity, (double, double, double) point)
        {
            var omega = velocity.Curl(point);
            return 0.5 * omega.Dot(omega);
        }

        /// <summary>
        /// Computes the helicity density v·ω at a point, quantifying the
        /// degree of linkage between velocity and vorticity.
        /// </summary>
        /// <param name="velocity">The velocity vector field v(r).</param>
        /// <param name="point">The point at which to evaluate.</param>
        public static double HelicityDensity(
            this VectorField velocity, (double, double, double) point)
        {
            var v = new Vector(velocity.fx(new Vector(point)),
                               velocity.fy(new Vector(point)),
                               velocity.fz(new Vector(point)));
            var omega = velocity.Curl(point);
            return v.Dot(omega);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Stream function (2-D incompressible)
        //
        //  vx = ∂ψ/∂y,  vy = −∂ψ/∂x
        // ═══════════════════════════════════════════════════════════════

        #region Stream Function

        /// <summary>
        /// Creates a 2-D velocity vector field from a stream function ψ(x, y):
        /// vx = ∂ψ/∂y, vy = −∂ψ/∂x, vz = 0.
        /// The resulting field is automatically divergence-free.
        /// </summary>
        /// <param name="streamFunction">The stream function ψ(r), depending on x and y.</param>
        public static VectorField VelocityFromStreamFunction(this ScalarField streamFunction)
        {
            // vx = ∂ψ/∂y
            Func<Vector, double> vx = r =>
            {
                var sf = streamFunction.f;
                return (sf(new Vector(r.x, r.y + 1e-6, r.z)) - sf(new Vector(r.x, r.y - 1e-6, r.z))) / 2e-6;
            };

            // vy = −∂ψ/∂x
            Func<Vector, double> vy = r =>
            {
                var sf = streamFunction.f;
                return -(sf(new Vector(r.x + 1e-6, r.y, r.z)) - sf(new Vector(r.x - 1e-6, r.y, r.z))) / 2e-6;
            };

            return new VectorField(vx, vy, _ => 0.0);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Drag & lift
        // ═══════════════════════════════════════════════════════════════

        #region Drag & Lift

        /// <summary>
        /// Computes the drag force on a body: F_D = ½ρv²C_D A.
        /// </summary>
        /// <param name="dragCoefficient">Drag coefficient C_D (dimensionless).</param>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        /// <param name="speed">Flow speed relative to the body in m/s.</param>
        /// <param name="referenceArea">Reference (frontal) area A in m².</param>
        public static double DragForce(
            this double dragCoefficient, double density, double speed, double referenceArea)
        {
            return 0.5 * density * speed * speed * dragCoefficient * referenceArea;
        }

        /// <summary>
        /// Computes the lift force on a body: F_L = ½ρv²C_L A.
        /// </summary>
        /// <param name="liftCoefficient">Lift coefficient C_L (dimensionless).</param>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        /// <param name="speed">Flow speed relative to the body in m/s.</param>
        /// <param name="referenceArea">Reference (wing) area A in m².</param>
        public static double LiftForce(
            this double liftCoefficient, double density, double speed, double referenceArea)
        {
            return 0.5 * density * speed * speed * liftCoefficient * referenceArea;
        }

        /// <summary>
        /// Computes the terminal velocity of a falling object:
        /// v_t = √(2mg / (ρC_D A)).
        /// </summary>
        /// <param name="mass">Mass of the object in kg.</param>
        /// <param name="dragCoefficient">Drag coefficient C_D.</param>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        /// <param name="referenceArea">Reference (frontal) area A in m².</param>
        public static double TerminalVelocity(
            this double mass, double dragCoefficient, double density, double referenceArea)
        {
            if (dragCoefficient <= 0 || density <= 0 || referenceArea <= 0)
                throw new ArgumentException("Drag coefficient, density, and area must be positive.");
            return Math.Sqrt(2.0 * mass * PhysicsConstants.GravitationalAcceleration
                / (density * dragCoefficient * referenceArea));
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Dimensionless numbers
        // ═══════════════════════════════════════════════════════════════

        #region Dimensionless Numbers

        /// <summary>
        /// Computes the Reynolds number: Re = ρvL / μ = vL / ν.
        /// Characterises the ratio of inertial to viscous forces.
        /// </summary>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        /// <param name="speed">Characteristic flow speed in m/s.</param>
        /// <param name="characteristicLength">Characteristic length L in m.</param>
        /// <param name="dynamicViscosity">Dynamic viscosity μ in Pa·s.</param>
        public static double ReynoldsNumber(
            this double density, double speed, double characteristicLength, double dynamicViscosity)
        {
            if (dynamicViscosity <= 0) throw new ArgumentException("Dynamic viscosity must be positive.");
            return density * speed * characteristicLength / dynamicViscosity;
        }

        /// <summary>
        /// Computes the Mach number: Ma = v / c.
        /// </summary>
        /// <param name="speed">Flow speed in m/s.</param>
        /// <param name="speedOfSound">Speed of sound in the medium in m/s
        /// (defaults to <see cref="PhysicsConstants.SpeedOfSoundAir"/>).</param>
        public static double MachNumber(this double speed, double speedOfSound = 0)
        {
            if (speedOfSound <= 0) speedOfSound = PhysicsConstants.SpeedOfSoundAir;
            return speed / speedOfSound;
        }

        /// <summary>
        /// Computes the Froude number: Fr = v / √(gL).
        /// Characterises the ratio of flow inertia to gravity.
        /// </summary>
        /// <param name="speed">Flow speed in m/s.</param>
        /// <param name="characteristicLength">Characteristic length L in m (e.g. water depth).</param>
        public static double FroudeNumber(this double speed, double characteristicLength)
        {
            if (characteristicLength <= 0) throw new ArgumentException("Characteristic length must be positive.");
            return speed / Math.Sqrt(PhysicsConstants.GravitationalAcceleration * characteristicLength);
        }

        /// <summary>
        /// Computes the Strouhal number: St = fL / v.
        /// Characterises oscillating flow mechanisms (e.g. vortex shedding).
        /// </summary>
        /// <param name="frequency">Oscillation frequency f in Hz.</param>
        /// <param name="characteristicLength">Characteristic length L in m.</param>
        /// <param name="speed">Flow speed v in m/s.</param>
        public static double StrouhalNumber(this double frequency, double characteristicLength, double speed)
        {
            if (speed <= 0) throw new ArgumentException("Speed must be positive.");
            return frequency * characteristicLength / speed;
        }

        /// <summary>
        /// Computes the Weber number: We = ρv²L / σ.
        /// Characterises the ratio of inertia to surface tension.
        /// </summary>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        /// <param name="speed">Flow speed in m/s.</param>
        /// <param name="characteristicLength">Characteristic length L in m.</param>
        /// <param name="surfaceTension">Surface tension σ in N/m.</param>
        public static double WeberNumber(
            this double density, double speed, double characteristicLength, double surfaceTension)
        {
            if (surfaceTension <= 0) throw new ArgumentException("Surface tension must be positive.");
            return density * speed * speed * characteristicLength / surfaceTension;
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Viscous / pipe flow
        // ═══════════════════════════════════════════════════════════════

        #region Viscous Flow

        /// <summary>
        /// Computes the Hagen–Poiseuille volumetric flow rate through a
        /// circular pipe of radius R and length L under pressure drop ΔP:
        /// Q = πR⁴ΔP / (8μL).
        /// </summary>
        /// <param name="radius">Pipe inner radius R in m.</param>
        /// <param name="pressureDrop">Pressure drop ΔP across the pipe in Pa.</param>
        /// <param name="dynamicViscosity">Dynamic viscosity μ in Pa·s.</param>
        /// <param name="length">Pipe length L in m.</param>
        public static double PoiseuilleFlowRate(
            this double radius, double pressureDrop, double dynamicViscosity, double length)
        {
            if (radius <= 0 || dynamicViscosity <= 0 || length <= 0)
                throw new ArgumentException("Radius, viscosity, and length must be positive.");
            return Math.PI * Math.Pow(radius, 4) * pressureDrop / (8.0 * dynamicViscosity * length);
        }

        /// <summary>
        /// Computes the Poiseuille velocity profile in a circular pipe:
        /// v(r) = (ΔP / 4μL)(R² − r²), where r is the radial distance from centre.
        /// </summary>
        /// <param name="radius">Pipe inner radius R in m.</param>
        /// <param name="pressureDrop">Pressure drop ΔP in Pa.</param>
        /// <param name="dynamicViscosity">Dynamic viscosity μ in Pa·s.</param>
        /// <param name="length">Pipe length L in m.</param>
        /// <param name="radialDistance">Radial distance r from the pipe centre in m.</param>
        public static double PoiseuilleVelocity(
            this double radius, double pressureDrop, double dynamicViscosity,
            double length, double radialDistance)
        {
            if (radialDistance < 0 || radialDistance > radius)
                throw new ArgumentException("Radial distance must be between 0 and the pipe radius.");
            return pressureDrop / (4.0 * dynamicViscosity * length) * (radius * radius - radialDistance * radialDistance);
        }

        /// <summary>
        /// Computes the Stokes drag force on a sphere in creeping flow:
        /// F = 6πμRv.
        /// </summary>
        /// <param name="dynamicViscosity">Dynamic viscosity μ in Pa·s.</param>
        /// <param name="radius">Sphere radius R in m.</param>
        /// <param name="speed">Relative speed v in m/s.</param>
        public static double StokesDrag(this double dynamicViscosity, double radius, double speed)
        {
            return 6.0 * Math.PI * dynamicViscosity * radius * speed;
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Hydrostatics
        // ═══════════════════════════════════════════════════════════════

        #region Hydrostatics

        /// <summary>
        /// Computes the hydrostatic pressure at a depth below the surface:
        /// p = p₀ + ρgh.
        /// </summary>
        /// <param name="depth">Depth h below the surface in m.</param>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        /// <param name="surfacePressure">Pressure at the surface p₀ in Pa
        /// (defaults to <see cref="PhysicsConstants.StandardAtmosphericPressure"/>).</param>
        public static double HydrostaticPressure(
            this double depth, double density,
            double surfacePressure = 0)
        {
            if (surfacePressure <= 0) surfacePressure = PhysicsConstants.StandardAtmosphericPressure;
            return surfacePressure + density * PhysicsConstants.GravitationalAcceleration * depth;
        }

        /// <summary>
        /// Computes the buoyant force (Archimedes' principle):
        /// F_b = ρ_fluid · V_displaced · g.
        /// </summary>
        /// <param name="fluidDensity">Density of the surrounding fluid in kg/m³.</param>
        /// <param name="displacedVolume">Volume of fluid displaced in m³.</param>
        public static double BuoyantForce(this double fluidDensity, double displacedVolume)
        {
            return fluidDensity * displacedVolume * PhysicsConstants.GravitationalAcceleration;
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Kinetic energy & momentum
        // ═══════════════════════════════════════════════════════════════

        #region Fluid Energy

        /// <summary>
        /// Computes the kinetic energy density of a flow: e_k = ½ρ|v|².
        /// </summary>
        /// <param name="velocity">Velocity vector at the point.</param>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        public static double KineticEnergyDensity(this Vector velocity, double density)
        {
            return 0.5 * density * velocity.Dot(velocity);
        }

        /// <summary>
        /// Computes the momentum density of a flow: ρv.
        /// </summary>
        /// <param name="velocity">Velocity vector at the point.</param>
        /// <param name="density">Fluid density ρ in kg/m³.</param>
        public static Vector MomentumDensity(this Vector velocity, double density)
        {
            return density * velocity;
        }

        #endregion
    }
}
