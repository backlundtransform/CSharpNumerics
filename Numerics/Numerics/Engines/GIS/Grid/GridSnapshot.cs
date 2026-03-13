using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Grid
{
    /// <summary>
    /// An immutable snapshot of scalar values over every cell of a <see cref="GeoGrid"/>
    /// at a single point in time. Values are stored in a flat array aligned with
    /// <see cref="GeoGrid.Index(int, int, int)"/>.
    /// </summary>
    public class GridSnapshot
    {
        private readonly double[] _values;

        /// <summary>The grid this snapshot is defined over.</summary>
        public GeoGrid Grid { get; }

        /// <summary>Simulation time in seconds for this snapshot.</summary>
        public double Time { get; }

        /// <summary>Zero-based time-step index.</summary>
        public int TimeIndex { get; }

        /// <summary>Number of cells (same as <see cref="GeoGrid.CellCount"/>).</summary>
        public int Count => _values.Length;

        /// <summary>
        /// Creates a snapshot from pre-computed values.
        /// The array length must equal <see cref="GeoGrid.CellCount"/>.
        /// </summary>
        public GridSnapshot(GeoGrid grid, double[] values, double time, int timeIndex)
        {
            if (grid == null) throw new ArgumentNullException(nameof(grid));
            if (values == null) throw new ArgumentNullException(nameof(values));
            if (values.Length != grid.CellCount)
                throw new ArgumentException(
                    $"Values length ({values.Length}) must equal grid cell count ({grid.CellCount}).");

            Grid = grid;
            _values = values;
            Time = time;
            TimeIndex = timeIndex;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Value access
        // ═══════════════════════════════════════════════════════════════

        /// <summary>Gets the value at a flat index.</summary>
        public double this[int flatIndex] => _values[flatIndex];

        /// <summary>Gets the value at 3-D indices.</summary>
        public double this[int ix, int iy, int iz] => _values[Grid.Index(ix, iy, iz)];

        /// <summary>Returns the raw values array (copy).</summary>
        public double[] GetValues() => (double[])_values.Clone();

        /// <summary>
        /// Returns the value at the cell nearest to <paramref name="position"/>.
        /// </summary>
        public double ValueAt(Vector position) => _values[Grid.NearestFlatIndex(position)];

        // ═══════════════════════════════════════════════════════════════
        //  Queries
        // ═══════════════════════════════════════════════════════════════

        /// <summary>Maximum value in the snapshot.</summary>
        public double Max()
        {
            double max = double.MinValue;
            for (int i = 0; i < _values.Length; i++)
                if (_values[i] > max) max = _values[i];
            return max;
        }

        /// <summary>Minimum value in the snapshot.</summary>
        public double Min()
        {
            double min = double.MaxValue;
            for (int i = 0; i < _values.Length; i++)
                if (_values[i] < min) min = _values[i];
            return min;
        }

        /// <summary>
        /// Returns all cells whose value exceeds <paramref name="threshold"/>.
        /// </summary>
        public List<GeoCell> CellsAbove(double threshold)
        {
            var result = new List<GeoCell>();
            for (int i = 0; i < _values.Length; i++)
            {
                if (_values[i] > threshold)
                {
                    result.Add(new GeoCell(
                        Grid.CellCentre(i), _values[i], TimeIndex, i));
                }
            }
            return result;
        }

        /// <summary>
        /// Enumerates all cells as <see cref="GeoCell"/> instances.
        /// </summary>
        public IEnumerable<GeoCell> AllCells()
        {
            for (int i = 0; i < _values.Length; i++)
            {
                yield return new GeoCell(
                    Grid.CellCentre(i), _values[i], TimeIndex, i);
            }
        }
    }
}
