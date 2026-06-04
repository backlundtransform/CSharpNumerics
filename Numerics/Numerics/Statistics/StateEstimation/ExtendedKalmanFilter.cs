using System;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Statistics.StateEstimation;

/// <summary>
/// Extended Kalman filter (EKF) for state estimation in a non-linear state-space model:
/// <code>
///   x_k = f(x_{k-1}) + w,   w ~ N(0, Q)   (process model)
///   z_k = h(x_k)     + v,   v ~ N(0, R)   (measurement model)
/// </code>
/// The covariance is propagated through the Jacobians of <c>f</c> and <c>h</c>, which the
/// caller supplies as functions evaluated at the current state estimate.
/// </summary>
public class ExtendedKalmanFilter
{
    private VectorN _state;
    private Matrix _covariance;

    /// <summary>
    /// Creates a filter with an initial state estimate and covariance.
    /// </summary>
    /// <param name="initialState">Initial state estimate <c>x̂₀</c>.</param>
    /// <param name="initialCovariance">Initial state covariance <c>P₀</c> (n×n).</param>
    public ExtendedKalmanFilter(VectorN initialState, Matrix initialCovariance)
    {
        if (initialCovariance.rowLength != initialState.Length ||
            initialCovariance.columnLength != initialState.Length)
            throw new ArgumentException("Covariance must be n×n where n is the state dimension.");

        _state = initialState;
        _covariance = initialCovariance;
    }

    /// <summary>Current state estimate <c>x̂</c>.</summary>
    public VectorN State => _state;

    /// <summary>Current state covariance <c>P</c>.</summary>
    public Matrix Covariance => _covariance;

    /// <summary>State dimension.</summary>
    public int Dimension => _state.Length;

    /// <summary>
    /// Prediction step through the non-linear process model:
    /// <code>x̂ = f(x̂),   P = F P Fᵀ + Q</code>
    /// where <c>F = ∂f/∂x</c> is evaluated at the prior state.
    /// </summary>
    /// <param name="f">Non-linear state transition function <c>f(x)</c>.</param>
    /// <param name="jacobianF">Jacobian <c>F(x) = ∂f/∂x</c> (n×n) evaluated at the prior state.</param>
    /// <param name="Q">Process noise covariance (n×n).</param>
    public void Predict(Func<VectorN, VectorN> f, Func<VectorN, Matrix> jacobianF, Matrix Q)
    {
        if (f == null) throw new ArgumentNullException(nameof(f));
        if (jacobianF == null) throw new ArgumentNullException(nameof(jacobianF));

        Matrix F = jacobianF(_state);            // linearise about the prior state
        _state = f(_state);
        _covariance = F * _covariance * F.Transpose() + Q;
    }

    /// <summary>
    /// Measurement update through the non-linear measurement model:
    /// <code>
    ///   y = z - h(x̂)
    ///   H = ∂h/∂x |_{x̂}
    ///   S = H P Hᵀ + R
    ///   K = P Hᵀ S⁻¹
    ///   x̂ = x̂ + K y
    ///   P  = (I - K H) P
    /// </code>
    /// </summary>
    /// <param name="h">Non-linear measurement function <c>h(x)</c>.</param>
    /// <param name="jacobianH">Jacobian <c>H(x) = ∂h/∂x</c> (m×n) evaluated at the current state.</param>
    /// <param name="R">Measurement noise covariance (m×m).</param>
    /// <param name="z">Measurement vector (length m).</param>
    public void Update(Func<VectorN, VectorN> h, Func<VectorN, Matrix> jacobianH, Matrix R, VectorN z)
    {
        if (h == null) throw new ArgumentNullException(nameof(h));
        if (jacobianH == null) throw new ArgumentNullException(nameof(jacobianH));

        Matrix H = jacobianH(_state);
        if (H.columnLength != Dimension)
            throw new ArgumentException("Jacobian H must have as many columns as the state dimension.", nameof(jacobianH));

        Matrix Ht = H.Transpose();

        VectorN innovation = z - h(_state);
        Matrix S = H * _covariance * Ht + R;
        Matrix K = _covariance * Ht * S.Inverse();

        _state = _state + K * innovation;
        _covariance = (KalmanFilter.Identity(Dimension) - K * H) * _covariance;
    }
}
