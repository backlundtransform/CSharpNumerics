namespace CSharpNumerics.Engines.Game.Fluids;

/// <summary>
/// Configuration for a game fluid solver (2D or 3D).
/// Controls grid resolution, physical parameters, solver quality, and boundary conditions.
/// </summary>
public class FluidConfig
{
    /// <summary>Grid cells along X axis.</summary>
    public int GridX { get; set; } = 64;

    /// <summary>Grid cells along Y axis.</summary>
    public int GridY { get; set; } = 64;

    /// <summary>Grid cells along Z axis (ignored for 2D solvers).</summary>
    public int GridZ { get; set; } = 64;

    /// <summary>Physical size of a single cell in metres.</summary>
    public double CellSize { get; set; } = 1.0;

    /// <summary>Kinematic viscosity ν in m²/s. Higher = more diffusive / smoother. Default 0.0001.</summary>
    public double Viscosity { get; set; } = 0.0001;

    /// <summary>Diffusion rate for the density (smoke/dye) field. Default 0.0001.</summary>
    public double DensityDiffusion { get; set; } = 0.0001;

    /// <summary>Number of Gauss-Seidel iterations for the pressure Poisson solve. More = more accurate but slower.</summary>
    public int PoissonIterations { get; set; } = 20;

    /// <summary>Number of Gauss-Seidel iterations for diffusion. Default 20.</summary>
    public int DiffusionIterations { get; set; } = 20;

    /// <summary>Boundary mode for the domain edges.</summary>
    public FluidBoundaryMode BoundaryMode { get; set; } = FluidBoundaryMode.Closed;

    /// <summary>Vorticity confinement strength ε. 0 = disabled, typical 0.1–5.0. Default 0.</summary>
    public double VorticityConfinementStrength { get; set; } = 0;

    /// <summary>Enable buoyancy force for temperature/density-driven flows. Default false.</summary>
    public bool EnableBuoyancy { get; set; } = false;

    /// <summary>Ambient temperature for buoyancy calculation. Default 300 K.</summary>
    public double AmbientTemperature { get; set; } = 300.0;

    /// <summary>Buoyancy thermal expansion coefficient β (1/K). Default 1/300.</summary>
    public double BuoyancyBeta { get; set; } = 1.0 / 300.0;

    /// <summary>
    /// Quality preset — adjusts Poisson iterations and diffusion automatically.
    /// </summary>
    public FluidQuality Quality
    {
        set
        {
            switch (value)
            {
                case FluidQuality.Low:
                    PoissonIterations = 10;
                    DiffusionIterations = 10;
                    break;
                case FluidQuality.Medium:
                    PoissonIterations = 20;
                    DiffusionIterations = 20;
                    break;
                case FluidQuality.High:
                    PoissonIterations = 40;
                    DiffusionIterations = 40;
                    break;
            }
        }
    }
}

/// <summary>Boundary condition mode for domain edges.</summary>
public enum FluidBoundaryMode
{
    /// <summary>Solid walls on all sides (no-slip or free-slip).</summary>
    Closed,
    /// <summary>Fluid wraps around edges (periodic).</summary>
    Periodic,
    /// <summary>Fluid flows out freely (open / Neumann).</summary>
    Open
}

/// <summary>Quality preset for solver accuracy vs. performance.</summary>
public enum FluidQuality
{
    /// <summary>Fewer solver iterations — fastest, most diffusive.</summary>
    Low,
    /// <summary>Balanced quality and performance.</summary>
    Medium,
    /// <summary>More solver iterations — most accurate, slowest.</summary>
    High
}
