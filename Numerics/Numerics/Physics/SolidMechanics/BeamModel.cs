using CSharpNumerics.Physics.SolidMechanics.Interfaces;
using System;

namespace CSharpNumerics.Physics.SolidMechanics;

/// <summary>
/// Default Euler–Bernoulli beam model using analytical deflection
/// formulas from <see cref="BeamExtensions"/>.
/// </summary>
public class BeamModel : IBeamModel
{
    public (double deflection, double moment, double shear) CantileverPointLoad(
        double P, double a, double L, double EI, double x)
    {
        double deflection, moment, shear;
        if (x <= a)
        {
            deflection = P * x * x * (3.0 * a - x) / (6.0 * EI);
            moment = P * (a - x);
            shear = P;
        }
        else
        {
            deflection = P * a * a * (3.0 * x - a) / (6.0 * EI);
            moment = 0;
            shear = 0;
        }
        return (deflection, moment, shear);
    }

    public (double deflection, double moment, double shear) CantileverUniformLoad(
        double q, double L, double EI, double x)
    {
        double deflection = q.CantileverUniformLoadDeflection(L, EI, x);
        double moment = q * (L - x) * (L - x) / 2.0;
        double shear = q * (L - x);
        return (deflection, moment, shear);
    }

    public (double deflection, double moment, double shear) SimplySupportedPointLoad(
        double P, double a, double L, double EI, double x)
    {
        double b = L - a;
        double deflection, moment, shear;
        if (x <= a)
        {
            deflection = P * b * x * (L * L - b * b - x * x) / (6.0 * L * EI);
            moment = P * b * x / L;
            shear = P * b / L;
        }
        else
        {
            deflection = P * a * (L - x) * (2.0 * L * x - a * a - x * x) / (6.0 * L * EI);
            moment = P * a * (L - x) / L;
            shear = -P * a / L;
        }
        return (deflection, moment, shear);
    }

    public (double deflection, double moment, double shear) SimplySupportedUniformLoad(
        double q, double L, double EI, double x)
    {
        double deflection = q.SimplySupportedUniformLoadDeflection(L, EI, x);
        double moment = q * x * (L - x) / 2.0;
        double shear = q * (L / 2.0 - x);
        return (deflection, moment, shear);
    }

    public (double deflection, double moment, double shear) FixedFixedPointLoad(
        double P, double a, double L, double EI, double x)
    {
        double b = L - a;
        double deflection;
        if (x <= a)
        {
            deflection = P * b * b * x * x * (3.0 * a * L - x * (3.0 * a + b))
                         / (6.0 * EI * L * L * L);
        }
        else
        {
            deflection = P * a * a * (L - x) * (L - x) * (3.0 * b * L - (L - x) * (3.0 * b + a))
                         / (6.0 * EI * L * L * L);
        }

        double Ma = -P * a * b * b / (L * L);
        double Ra = P * b * b * (3.0 * a + b) / (L * L * L);
        double moment, shear;
        if (x <= a)
        {
            shear = Ra;
            moment = Ma + Ra * x;
        }
        else
        {
            shear = Ra - P;
            moment = Ma + Ra * x - P * (x - a);
        }
        return (deflection, moment, shear);
    }

    public (double deflection, double moment, double shear) FixedFixedUniformLoad(
        double q, double L, double EI, double x)
    {
        double Lx = L - x;
        double deflection = q * x * x * Lx * Lx / (24.0 * EI);
        double moment = q * (L * x - x * x) / 2.0 - q * L * L / 12.0;
        double shear = q * (L / 2.0 - x);
        return (deflection, moment, shear);
    }

    public double BendingStress(double moment, double halfHeight, double secondMoment)
    {
        if (secondMoment <= 0 || halfHeight <= 0) return 0;
        return Math.Abs(moment) * halfHeight / secondMoment;
    }
}
