using System;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Statistics.StateEstimation;

/// <summary>
/// Linear Kalman filter for recursive state estimation in a linear-Gaussian
/// state-space model:
/// <code>
///   x_k = F x_{k-1} + w,   w ~ N(0, Q)   (process model)
///   z_k = H x_k     + v,   v ~ N(0, R)   (measurement model)
/// </code>
/// Use <see cref="Predict"/> to propagate the state through the process model and
/// <see cref="Update"/> to fold in a new measurement. The current state estimate
/// <c>x̂</c> and its covariance <c>P</c> are exposed via <see cref="State"/> and
/// <see cref="Covariance"/>.
/// </summary>
public class KalmanFilter
{
    private VectorN _state;
    private Matrix _covariance;

    /// <summary>
    /// Creates a filter with an initial state estimate and covariance.
    /// </summary>
    /// <param name="initialState">Initial state estimate <c>x̂₀</c>.</param>
    /// <param name="initialCovariance">Initial state covariance <c>P₀</c> (n×n).</param>
    public KalmanFilter(VectorN initialState, Matrix initialCovariance)
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
    /// Prediction step: propagates the state and covariance through the process model.
    /// <code>x̂ = F x̂,   P = F P Fᵀ + Q</code>
    /// </summary>
    /// <param name="F">State transition matrix (n×n).</param>
    /// <param name="Q">Process noise covariance (n×n).</param>
    public void Predict(Matrix F, Matrix Q)
    {
        ValidateSquare(F, nameof(F));
        ValidateSquare(Q, nameof(Q));

        _state = F * _state;
        _covariance = F * _covariance * F.Transpose() + Q;
    }

    /// <summary>
    /// Prediction step with a control input <c>u</c>:
    /// <code>x̂ = F x̂ + B u,   P = F P Fᵀ + Q</code>
    /// </summary>
    /// <param name="F">State transition matrix (n×n).</param>
    /// <param name="Q">Process noise covariance (n×n).</param>
    /// <param name="B">Control matrix (n×c).</param>
    /// <param name="u">Control vector (length c).</param>
    public void Predict(Matrix F, Matrix Q, Matrix B, VectorN u)
    {
        ValidateSquare(F, nameof(F));
        ValidateSquare(Q, nameof(Q));

        _state = F * _state + B * u;
        _covariance = F * _covariance * F.Transpose() + Q;
    }

    /// <summary>
    /// Measurement update (correction) step. Folds measurement <c>z</c> into the estimate:
    /// <code>
    ///   y = z - H x̂                  (innovation)
    ///   S = H P Hᵀ + R               (innovation covariance)
    ///   K = P Hᵀ S⁻¹                 (Kalman gain)
    ///   x̂ = x̂ + K y
    ///   P  = (I - K H) P
    /// </code>
    /// </summary>
    /// <param name="H">Measurement matrix (m×n).</param>
    /// <param name="R">Measurement noise covariance (m×m).</param>
    /// <param name="z">Measurement vector (length m).</param>
    public void Update(Matrix H, Matrix R, VectorN z)
    {
        if (H.columnLength != Dimension)
            throw new ArgumentException("H must have as many columns as the state dimension.", nameof(H));
        if (z.Length != H.rowLength)
            throw new ArgumentException("Measurement length must match the number of rows in H.", nameof(z));
        ValidateSquare(R, nameof(R));

        Matrix Ht = H.Transpose();

        VectorN innovation = z - H * _state;          // y
        Matrix S = H * _covariance * Ht + R;          // innovation covariance
        Matrix K = _covariance * Ht * S.Inverse();    // Kalman gain (n×m)

        _state = _state + K * innovation;

        Matrix identity = Identity(Dimension);
        _covariance = (identity - K * H) * _covariance;
    }

    private static void ValidateSquare(Matrix m, string name)
    {
        if (m.rowLength != m.columnLength)
            throw new ArgumentException($"{name} must be square.", name);
    }

    internal static Matrix Identity(int n)
    {
        var values = new double[n, n];
        for (int i = 0; i < n; i++)
            values[i, i] = 1.0;
        return new Matrix(values);
    }
}
