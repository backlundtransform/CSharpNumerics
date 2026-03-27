using CSharpNumerics.Numerics.Objects;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Analysis
{
    /// <summary>
    /// The type of exposure metric used to generate a polygon.
    /// </summary>
    public enum ExposureType
    {
        /// <summary>Peak instantaneous value across all time steps.</summary>
        Peak,

        /// <summary>Time-integrated (cumulative) exposure over the simulation duration.</summary>
        Integrated
    }

    /// <summary>
    /// A closed polygon boundary delineating the area where an exposure metric
    /// (peak or time-integrated) exceeds a specified threshold.
    /// The polygon is computed on the ground-level (z = zMin) slice of the grid.
    /// </summary>
    public class ExposurePolygon
    {
        /// <summary>Ordered ring of boundary vertices (closed: first == last).</summary>
        public List<Vector> Boundary { get; }

        /// <summary>The threshold value that defines this polygon.</summary>
        public double Threshold { get; }

        /// <summary>Name of the data layer used (e.g. "concentration", "activity", "dose").</summary>
        public string LayerName { get; }

        /// <summary>Whether this is a peak or integrated exposure polygon.</summary>
        public ExposureType Type { get; }

        /// <summary>Number of ground-level cells that exceed the threshold.</summary>
        public int ExceedanceCellCount { get; }

        /// <summary>
        /// Approximate area of the exceedance zone in square metres
        /// (cell count × cell area).
        /// </summary>
        public double AreaSquareMetres { get; }

        internal ExposurePolygon(
            List<Vector> boundary,
            double threshold,
            string layerName,
            ExposureType type,
            int exceedanceCellCount,
            double areaSquareMetres)
        {
            Boundary = boundary;
            Threshold = threshold;
            LayerName = layerName;
            Type = type;
            ExceedanceCellCount = exceedanceCellCount;
            AreaSquareMetres = areaSquareMetres;
        }
    }
}
