using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

namespace CSharpNumerics.Engines.Game.Flight;

/// <summary>
/// Configuration data for an aircraft: geometry, mass properties, aerodynamic models,
/// and propulsion. Immutable after construction — swap configs to change aircraft type.
/// </summary>
public class AircraftConfig
{
    /// <summary>Aircraft name / identifier.</summary>
    public string Name { get; }

    // ─── Mass properties ─────────────────────────────────────────
    /// <summary>Total mass in kg.</summary>
    public double Mass { get; }

    /// <summary>Moment of inertia about body X axis (roll) in kg·m².</summary>
    public double Ixx { get; }

    /// <summary>Moment of inertia about body Y axis (pitch) in kg·m².</summary>
    public double Iyy { get; }

    /// <summary>Moment of inertia about body Z axis (yaw) in kg·m².</summary>
    public double Izz { get; }

    /// <summary>Product of inertia Ixz in kg·m² (typically nonzero for symmetric aircraft). Default 0.</summary>
    public double Ixz { get; }

    // ─── Wing geometry ───────────────────────────────────────────
    /// <summary>Wing reference area (m²).</summary>
    public double WingArea { get; }

    /// <summary>Wingspan (m).</summary>
    public double WingSpan { get; }

    /// <summary>Mean aerodynamic chord (m).</summary>
    public double MeanChord { get; }

    // ─── Aerodynamic models ──────────────────────────────────────
    /// <summary>Airfoil lift/drag model.</summary>
    public AirfoilModel Airfoil { get; }

    /// <summary>Elevator control surface.</summary>
    public ControlSurface Elevator { get; }

    /// <summary>Aileron control surface.</summary>
    public ControlSurface Aileron { get; }

    /// <summary>Rudder control surface.</summary>
    public ControlSurface Rudder { get; }

    // ─── Propulsion ──────────────────────────────────────────────
    /// <summary>Engine / propulsion model.</summary>
    public PropulsionModel Propulsion { get; }

    // ─── Misc ────────────────────────────────────────────────────
    /// <summary>Additional parasitic drag area (m²) for gear, flaps, etc. Default 0.</summary>
    public double GearDragArea { get; }

    /// <summary>Additional Cl increment per unit flap setting. Default 0.4.</summary>
    public double FlapClIncrement { get; }

    /// <summary>Additional Cd increment per unit flap setting. Default 0.02.</summary>
    public double FlapCdIncrement { get; }

    /// <summary>
    /// Roll damping derivative Clp (dimensionless). Typical: −0.4 to −0.6.
    /// Produces roll-rate damping moment Cl = Clp · (p · b) / (2V).
    /// </summary>
    public double Clp { get; }

    /// <summary>
    /// Pitch damping derivative Cmq (dimensionless). Typical: −10 to −20.
    /// Produces pitch-rate damping moment Cm = Cmq · (q · c̄) / (2V).
    /// </summary>
    public double Cmq { get; }

    /// <summary>
    /// Yaw damping derivative Cnr (dimensionless). Typical: −0.1 to −0.3.
    /// Produces yaw-rate damping moment Cn = Cnr · (r · b) / (2V).
    /// </summary>
    public double Cnr { get; }

    /// <summary>
    /// Creates an aircraft configuration.
    /// </summary>
    public AircraftConfig(
        string name,
        double mass,
        double ixx, double iyy, double izz, double ixz,
        double wingArea, double wingSpan, double meanChord,
        AirfoilModel airfoil,
        ControlSurface elevator,
        ControlSurface aileron,
        ControlSurface rudder,
        PropulsionModel propulsion,
        double gearDragArea = 0,
        double flapClIncrement = 0.4,
        double flapCdIncrement = 0.02,
        double clp = -0.5,
        double cmq = -15.0,
        double cnr = -0.15)
    {
        Name = name;
        Mass = mass;
        Ixx = ixx; Iyy = iyy; Izz = izz; Ixz = ixz;
        WingArea = wingArea;
        WingSpan = wingSpan;
        MeanChord = meanChord;
        Airfoil = airfoil;
        Elevator = elevator;
        Aileron = aileron;
        Rudder = rudder;
        Propulsion = propulsion;
        GearDragArea = gearDragArea;
        FlapClIncrement = flapClIncrement;
        FlapCdIncrement = flapCdIncrement;
        Clp = clp;
        Cmq = cmq;
        Cnr = cnr;
    }

    /// <summary>Aspect ratio: AR = b² / S.</summary>
    public double AspectRatio => WingSpan * WingSpan / WingArea;

    // ═══════════════════════════════════════════════════════════════
    //  Presets
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// A generic single-engine light aircraft resembling a Cessna 172.
    /// Useful for testing and prototyping.
    /// </summary>
    public static AircraftConfig GenericLightAircraft()
    {
        return new AircraftConfig(
            name: "Generic Light Aircraft",
            mass: 1043,
            ixx: 1285, iyy: 1824, izz: 2666, ixz: 0,
            wingArea: 16.17,
            wingSpan: 10.92,
            meanChord: 1.49,
            airfoil: AirfoilModel.NACASymmetric(),
            elevator: ControlSurface.Elevator(),
            aileron: ControlSurface.Aileron(),
            rudder: ControlSurface.Rudder(),
            propulsion: PropulsionModel.Propeller(
                maxPower: 119_000,       // 160 HP
                propellerEfficiency: 0.80,
                staticThrust: 3800),
            gearDragArea: 0.3,
            flapClIncrement: 0.5,
            flapCdIncrement: 0.025);
    }

    /// <summary>
    /// A generic twin-engine jet resembling a business jet.
    /// </summary>
    public static AircraftConfig GenericJet()
    {
        return new AircraftConfig(
            name: "Generic Jet",
            mass: 9000,
            ixx: 12000, iyy: 26000, izz: 36000, ixz: 1400,
            wingArea: 30.0,
            wingSpan: 15.0,
            meanChord: 2.0,
            airfoil: AirfoilModel.NACASymmetric(
                clAlpha: 5.7, alphaStall: 0.27, clMax: 1.6, cd0: 0.015, k: 0.045),
            elevator: ControlSurface.Elevator(clDelta: 4.5),
            aileron: ControlSurface.Aileron(clDelta: 3.8),
            rudder: ControlSurface.Rudder(clDelta: 3.2),
            propulsion: PropulsionModel.Jet(
                maxThrustSeaLevel: 22000,
                engineCount: 2),
            gearDragArea: 0.6);
    }
}
