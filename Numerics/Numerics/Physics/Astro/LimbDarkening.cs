using System;

namespace CSharpNumerics.Physics.Astro;

public static class LimbDarkening
{
    /// <summary>
    /// Linear limb darkening law: I(μ) = 1 − u1·(1 − μ).
    /// </summary>
    /// <param name="mu">Cosine of angle between line-of-sight and surface normal (0–1).</param>
    /// <param name="u1">Linear coefficient.</param>
    public static double Linear(double mu, double u1)
    {
        return 1.0 - u1 * (1.0 - mu);
    }

    /// <summary>
    /// Quadratic limb darkening law: I(μ) = 1 − u1·(1−μ) − u2·(1−μ)².
    /// </summary>
    public static double Quadratic(double mu, double u1, double u2)
    {
        double oneMinusMu = 1.0 - mu;
        return 1.0 - u1 * oneMinusMu - u2 * oneMinusMu * oneMinusMu;
    }

    /// <summary>
    /// Nonlinear four-parameter limb darkening law (Claret 2000):
    /// I(μ) = 1 − c1·(1−μ^(1/2)) − c2·(1−μ) − c3·(1−μ^(3/2)) − c4·(1−μ²).
    /// </summary>
    public static double NonlinearFourParam(double mu, double c1, double c2, double c3, double c4)
    {
        double sqrtMu = Math.Sqrt(mu);
        return 1.0
            - c1 * (1.0 - sqrtMu)
            - c2 * (1.0 - mu)
            - c3 * (1.0 - sqrtMu * mu)
            - c4 * (1.0 - mu * mu);
    }

    /// <summary>
    /// Computes the intensity profile for an array of μ values.
    /// </summary>
    /// <param name="model">Limb darkening model type.</param>
    /// <param name="coefficients">Model coefficients (1 for Linear, 2 for Quadratic, 4 for NonlinearFourParam; 0 for Uniform).</param>
    /// <param name="muArray">Array of μ values (0–1).</param>
    /// <returns>Intensity at each μ value.</returns>
    public static double[] IntensityProfile(LimbDarkeningModel model, double[] coefficients, double[] muArray)
    {
        if (muArray == null) throw new ArgumentNullException(nameof(muArray));
        if (coefficients == null) throw new ArgumentNullException(nameof(coefficients));

        double[] result = new double[muArray.Length];

        for (int i = 0; i < muArray.Length; i++)
        {
            double mu = muArray[i];
            switch (model)
            {
                case LimbDarkeningModel.Uniform:
                    result[i] = 1.0;
                    break;
                case LimbDarkeningModel.Linear:
                    result[i] = Linear(mu, coefficients[0]);
                    break;
                case LimbDarkeningModel.Quadratic:
                    result[i] = Quadratic(mu, coefficients[0], coefficients[1]);
                    break;
                case LimbDarkeningModel.NonlinearFourParam:
                    result[i] = NonlinearFourParam(mu, coefficients[0], coefficients[1], coefficients[2], coefficients[3]);
                    break;
                default:
                    result[i] = 1.0;
                    break;
            }
        }

        return result;
    }
}
