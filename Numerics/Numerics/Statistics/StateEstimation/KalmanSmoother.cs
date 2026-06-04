using System;
using System.Collections.Generic;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Statistics.StateEstimation;

/// <summary>
/// Result of a Rauch–Tung–Striebel (RTS) smoothing pass: the smoothed state
/// estimates and covariances at every time step, together with the forward
/// (filter-only) estimates for comparison.
/// </summary>
public class KalmanSmootherResult
{
    /// <summary>Smoothed state estimates <c>x̂ₖ|ₙ</c> for k = 0 … n-1.</summary>
    public VectorN[] SmoothedStates { get; }

    /// <summary>Smoothed state covariances <c>Pₖ|ₙ</c> for k = 0 … n-1.</summary>
    public Matrix[] SmoothedCovariances { get; }

    /// <summary>Forward (filter-only) state estimates <c>x̂ₖ|ₖ</c>.</summary>
    public VectorN[] FilteredStates { get; }

    /// <summary>Forward (filter-only) state covariances <c>Pₖ|ₖ</c>.</summary>
    public Matrix[] FilteredCovariances { get; }

    public KalmanSmootherResult(
        VectorN[] smoothedStates, Matrix[] smoothedCovariances,
        VectorN[] filteredStates, Matrix[] filteredCovariances)
    {
        SmoothedStates = smoothedStates;
        SmoothedCovariances = smoothedCovariances;
        FilteredStates = filteredStates;
        FilteredCovariances = filteredCovariances;
    }
}

/// <summary>
/// Rauch–Tung–Striebel (RTS) fixed-interval Kalman smoother for offline analysis.
/// Runs a forward Kalman filter over the full measurement sequence and then a
/// backward recursion that refines each estimate using all later measurements.
/// Smoothed covariances are never larger than the forward-only covariances, so
/// the smoothed track is less noisy — at the cost of needing the whole batch up front.
/// Assumes time-invariant model matrices F, Q, H, R.
/// </summary>
public class KalmanSmoother
{
    /// <summary>
    /// Smooths a sequence of measurements under a linear-Gaussian model.
    /// </summary>
    /// <param name="measurements">Measurement vectors z₀ … z_{n-1} (each length m).</param>
    /// <param name="initialState">Initial state estimate <c>x̂₀</c>.</param>
    /// <param name="initialCovariance">Initial state covariance <c>P₀</c> (n×n).</param>
    /// <param name="F">State transition matrix (n×n).</param>
    /// <param name="Q">Process noise covariance (n×n).</param>
    /// <param name="H">Measurement matrix (m×n).</param>
    /// <param name="R">Measurement noise covariance (m×m).</param>
    public KalmanSmootherResult Smooth(
        IReadOnlyList<VectorN> measurements,
        VectorN initialState, Matrix initialCovariance,
        Matrix F, Matrix Q, Matrix H, Matrix R)
    {
        if (measurements == null) throw new ArgumentNullException(nameof(measurements));
        if (measurements.Count == 0) throw new ArgumentException("At least one measurement is required.", nameof(measurements));

        int n = measurements.Count;
        int dim = initialState.Length;

        var filteredStates = new VectorN[n];
        var filteredCovariances = new Matrix[n];
        var predictedStates = new VectorN[n];      // x̂ₖ|ₖ₋₁
        var predictedCovariances = new Matrix[n];  // Pₖ|ₖ₋₁

        Matrix Ft = F.Transpose();
        Matrix Ht = H.Transpose();
        Matrix identity = KalmanFilter.Identity(dim);

        // ---- Forward pass (standard Kalman filter, storing the prediction priors) ----
        VectorN prevState = initialState;
        Matrix prevCovariance = initialCovariance;

        for (int k = 0; k < n; k++)
        {
            // Predict
            VectorN xPred = F * prevState;
            Matrix pPred = F * prevCovariance * Ft + Q;

            predictedStates[k] = xPred;
            predictedCovariances[k] = pPred;

            // Update
            VectorN innovation = measurements[k] - H * xPred;
            Matrix S = H * pPred * Ht + R;
            Matrix K = pPred * Ht * S.Inverse();

            VectorN xFilt = xPred + K * innovation;
            Matrix pFilt = (identity - K * H) * pPred;

            filteredStates[k] = xFilt;
            filteredCovariances[k] = pFilt;

            prevState = xFilt;
            prevCovariance = pFilt;
        }

        // ---- Backward RTS pass ----
        var smoothedStates = new VectorN[n];
        var smoothedCovariances = new Matrix[n];

        smoothedStates[n - 1] = filteredStates[n - 1];
        smoothedCovariances[n - 1] = filteredCovariances[n - 1];

        for (int k = n - 2; k >= 0; k--)
        {
            // Smoother gain: C = P_f[k] Fᵀ (P_pred[k+1])⁻¹
            Matrix C = filteredCovariances[k] * Ft * predictedCovariances[k + 1].Inverse();

            smoothedStates[k] = filteredStates[k]
                + C * (smoothedStates[k + 1] - predictedStates[k + 1]);

            smoothedCovariances[k] = filteredCovariances[k]
                + C * (smoothedCovariances[k + 1] - predictedCovariances[k + 1]) * C.Transpose();
        }

        return new KalmanSmootherResult(
            smoothedStates, smoothedCovariances,
            filteredStates, filteredCovariances);
    }
}
