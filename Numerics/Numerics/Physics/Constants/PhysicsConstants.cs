namespace CSharpNumerics.Physics.Constants;

public static class PhysicsConstants
{
    // Mechanics / Gravitation
    public const double GravitationalAcceleration = 9.80665; // m/s²
    public const double GravitationalConstant = 6.67430e-11; // m³/(kg·s²)
    public const double EscapeVelocityEarth = 11186; // m/s
    public const double EarthMass = 5.97237e24; // kg
    public const double MoonMass = 7.342e22; // kg
    public const double SolarMass = 1.98847e30; // kg

    // Electromagnetism
    public const double SpeedOfLight = 299792458; // m/s
    public const double VacuumPermittivity = 8.854187817e-12; // F/m
    public const double VacuumPermeability = 1.25663706212e-6; // H/m
    public const double ElementaryCharge = 1.602176634e-19; // C
    public const double FaradayConstant = 96485.33212; // C/mol
    public const double MagneticFluxQuantum = 2.067833848e-15; // Wb
    public const double QuantumOfConductance = 7.748091729e-5; // S

    // Quantum / Atomic / Particle
    public const double PlancksConstant = 6.62607015e-34; // J·s
    public const double FineStructureConstant = 7.2973525693e-3; // dimensionless
    public const double RydbergConstant = 10973731.568160; // 1/m
    public const double RydbergEnergy = 2.1798723611035e-18; // J
    public const double ElectronVolt = 1.602176634e-19; // J
    public const double BohrRadius = 5.29177210903e-11; // meters
    public const double ComptonWavelength = 2.42631023867e-12; // meters
    public const double ClassicalElectronRadius = 2.8179403262e-15; // meters
    public const double ThomsonCrossSection = 6.6524587321e-29; // m²
    public const double ElectronMass = 9.1093837015e-31; // kg
    public const double ProtonMass = 1.67262192369e-27; // kg
    public const double NeutronMass = 1.67492749804e-27; // kg

    // Thermodynamics / Statistical mechanics
    public const double BoltzmannConstant = 1.380649e-23; // J/K
    public const double AvogadrosNumber = 6.02214076e23; // 1/mol
    public const double GasConstant = 8.314462618; // J/(mol·K)
    public const double StefanBoltzmannConstant = 5.670374419e-8; // W/(m²·K⁴)

    // Atmosphere / Fluids / Materials
    public const double StandardAtmosphericPressure = 101325; // Pa
    public const double SpeedOfSoundAir = 343; // m/s at 20°C
    public const double MolarMassOfAir = 0.02897; // kg/mol
    public const double SpecificHeatCapacityOfWater = 4184; // J/(kg·K)
    public const double LatentHeatOfFusionWater = 334000; // J/kg
    public const double LatentHeatOfVaporizationWater = 2260000; // J/kg
    public const double StandardTemperature = 273.15; // Kelvin (0°C)
    public const double TriplePointOfWater = 273.16; // Kelvin
    public const double CriticalTemperatureOfWater = 647.096; // Kelvin

    // Astronomy / Cosmology / Astrophysics
    public const double AstronomicalUnit = 1.495978707e11; // meters
    public const double LightYear = 9.4607e15; // meters
    public const double Parsec = 3.0857e16; // meters
    public const double EarthRadius = 6.371e6; // meters
    public const double MoonRadius = 1.7371e6; // meters
    public const double SolarLuminosity = 3.828e26; // Watts
    public const double HubbleConstant = 2.2e-18; // 1/s
}
