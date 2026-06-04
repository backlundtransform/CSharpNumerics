using System;

namespace CSharpNumerics.Statistics.TimeSeriesAnalysis;

/// <summary>
/// Seasonal component model for <see cref="HoltWintersSmoothing"/>.
/// </summary>
public enum SeasonalType
{
    /// <summary>Season is added to the level + trend (constant-amplitude seasonality).</summary>
    Additive,

    /// <summary>Season multiplies the level + trend (amplitude grows with the level).</summary>
    Multiplicative
}

/// <summary>
/// Holt–Winters triple exponential smoothing: models a time series as the sum (or product)
/// of a slowly varying <em>level</em>, a linear <em>trend</em>, and a repeating <em>seasonal</em>
/// component of fixed period. A lightweight forecaster for series with strong periodicity
/// (e.g. a 24-hour diurnal wastewater baseline).
/// <code>
///   Additive:        ŷ_{t+1} = level_t + trend_t + season_{t+1-m}
///   Multiplicative:  ŷ_{t+1} = (level_t + trend_t) · season_{t+1-m}
/// </code>
/// Initial level, trend, and seasonal indices are estimated automatically from the first
/// two full seasons of data.
/// </summary>
public class HoltWintersSmoothing
{
    private readonly int _seasonLength;
    private readonly double _alpha;
    private readonly double _beta;
    private readonly double _gamma;
    private readonly SeasonalType _seasonalType;

    private double[] _seasonal;
    private double[] _fitted;
    private bool _isFitted;

    /// <summary>
    /// Creates a Holt–Winters smoother.
    /// </summary>
    /// <param name="seasonLength">Number of samples per season (e.g. 24 for hourly data with a daily cycle). Must be ≥ 2.</param>
    /// <param name="alpha">Level smoothing factor ∈ [0, 1].</param>
    /// <param name="beta">Trend smoothing factor ∈ [0, 1].</param>
    /// <param name="gamma">Seasonal smoothing factor ∈ [0, 1].</param>
    /// <param name="seasonalType">Additive or multiplicative seasonality.</param>
    public HoltWintersSmoothing(int seasonLength, double alpha, double beta, double gamma,
        SeasonalType seasonalType = SeasonalType.Additive)
    {
        if (seasonLength < 2)
            throw new ArgumentException("Season length must be at least 2.", nameof(seasonLength));
        ValidateFactor(alpha, nameof(alpha));
        ValidateFactor(beta, nameof(beta));
        ValidateFactor(gamma, nameof(gamma));

        _seasonLength = seasonLength;
        _alpha = alpha;
        _beta = beta;
        _gamma = gamma;
        _seasonalType = seasonalType;
    }

    /// <summary>Final level component after the last observation.</summary>
    public double Level { get; private set; }

    /// <summary>Final trend (per-step slope) component after the last observation.</summary>
    public double Trend { get; private set; }

    /// <summary>The most recent <c>seasonLength</c> seasonal indices.</summary>
    public double[] Seasonal => _isFitted ? (double[])_seasonal.Clone() : throw new InvalidOperationException("Call Fit first.");

    /// <summary>In-sample one-step-ahead fitted values, aligned with the input series.</summary>
    public double[] Fitted => _isFitted ? (double[])_fitted.Clone() : throw new InvalidOperationException("Call Fit first.");

    /// <summary>
    /// Fits the model to the data, estimating the level, trend, and seasonal components.
    /// Requires at least two full seasons of observations.
    /// </summary>
    /// <param name="data">The observed time series.</param>
    public void Fit(double[] data)
    {
        if (data == null) throw new ArgumentNullException(nameof(data));
        if (data.Length < 2 * _seasonLength)
            throw new ArgumentException($"At least two full seasons ({2 * _seasonLength} samples) are required.", nameof(data));

        int m = _seasonLength;

        // --- Initialisation from the first two seasons ---
        double firstSeasonMean = Mean(data, 0, m);
        double secondSeasonMean = Mean(data, m, m);

        double level = firstSeasonMean;
        double trend = (secondSeasonMean - firstSeasonMean) / m;

        var seasonal = new double[m];
        for (int i = 0; i < m; i++)
        {
            seasonal[i] = _seasonalType == SeasonalType.Additive
                ? data[i] - firstSeasonMean
                : (firstSeasonMean != 0.0 ? data[i] / firstSeasonMean : 1.0);
        }

        _fitted = new double[data.Length];

        // First season has no prior one-step forecast — use the seasonal initial guess.
        for (int t = 0; t < m; t++)
        {
            _fitted[t] = _seasonalType == SeasonalType.Additive
                ? level + trend + seasonal[t]
                : (level + trend) * seasonal[t];
        }

        // --- Recursive smoothing ---
        for (int t = m; t < data.Length; t++)
        {
            int sIdx = t % m;
            double prevSeasonal = seasonal[sIdx];
            double prevLevel = level;

            // One-step-ahead forecast made before seeing data[t]
            _fitted[t] = _seasonalType == SeasonalType.Additive
                ? prevLevel + trend + prevSeasonal
                : (prevLevel + trend) * prevSeasonal;

            if (_seasonalType == SeasonalType.Additive)
            {
                level = _alpha * (data[t] - prevSeasonal) + (1 - _alpha) * (prevLevel + trend);
                trend = _beta * (level - prevLevel) + (1 - _beta) * trend;
                seasonal[sIdx] = _gamma * (data[t] - level) + (1 - _gamma) * prevSeasonal;
            }
            else
            {
                double safeSeasonal = prevSeasonal != 0.0 ? prevSeasonal : 1e-12;
                level = _alpha * (data[t] / safeSeasonal) + (1 - _alpha) * (prevLevel + trend);
                trend = _beta * (level - prevLevel) + (1 - _beta) * trend;
                double safeLevel = level != 0.0 ? level : 1e-12;
                seasonal[sIdx] = _gamma * (data[t] / safeLevel) + (1 - _gamma) * prevSeasonal;
            }
        }

        Level = level;
        Trend = trend;
        _seasonal = seasonal;
        _isFitted = true;
    }

    /// <summary>
    /// Forecasts <paramref name="horizon"/> steps beyond the end of the fitted series.
    /// </summary>
    /// <param name="horizon">Number of steps to forecast (≥ 1).</param>
    public double[] Forecast(int horizon)
    {
        if (!_isFitted) throw new InvalidOperationException("Call Fit before Forecast.");
        if (horizon < 1) throw new ArgumentException("Horizon must be at least 1.", nameof(horizon));

        int m = _seasonLength;
        var forecast = new double[horizon];

        for (int h = 1; h <= horizon; h++)
        {
            // Seasonal index cycles through the last estimated season.
            double seasonal = _seasonal[(h - 1) % m];
            forecast[h - 1] = _seasonalType == SeasonalType.Additive
                ? Level + h * Trend + seasonal
                : (Level + h * Trend) * seasonal;
        }

        return forecast;
    }

    private static double Mean(double[] data, int start, int count)
    {
        double sum = 0.0;
        for (int i = 0; i < count; i++) sum += data[start + i];
        return sum / count;
    }

    private static void ValidateFactor(double value, string name)
    {
        if (value < 0.0 || value > 1.0)
            throw new ArgumentException($"{name} must be in [0, 1].", name);
    }
}
