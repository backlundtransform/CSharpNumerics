using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Engines.GIS.Grid
{
    /// <summary>
    /// A lightweight value representing a single spatial cell in a <see cref="GeoGrid"/>
    /// with an associated scalar value (e.g. concentration) and time-step index.
    /// </summary>
    public readonly struct GeoCell
    {
        /// <summary>Centre position of the cell in world coordinates.</summary>
        public Vector Position { get; }

        /// <summary>Scalar value at this cell (e.g. concentration in kg/m³).</summary>
        public double Value { get; }

        /// <summary>Zero-based time-step index this cell belongs to.</summary>
        public int TimeIndex { get; }

        /// <summary>Flat grid index (row-major iz → iy → ix).</summary>
        public int GridIndex { get; }

        public GeoCell(Vector position, double value, int timeIndex, int gridIndex)
        {
            Position = position;
            Value = value;
            TimeIndex = timeIndex;
            GridIndex = gridIndex;
        }

        public override string ToString() =>
            $"({Position.x:F1}, {Position.y:F1}, {Position.z:F1}) = {Value:E3} [t={TimeIndex}]";
    }
}
