using CSharpNumerics.ML.DimensionalityReduction.Algorithms;
using CSharpNumerics.Numerics.Objects;

namespace NumericTest;

/// <summary>
/// Tests for Principal Component Analysis.
/// </summary>
/// <remarks>
/// Written before the eigensolver behind PCA was swapped from power iteration with deflation to
/// the shared <c>EigenDecomposition</c>. They assert mathematical properties of the
/// decomposition — orthonormal components, variance ordering, alignment with known structure —
/// so they hold for any correct eigensolver rather than pinning one implementation's output.
/// </remarks>
[TestClass]
public class PcaTests
{
    /// <summary>
    /// Points along the line y = 2x, so the leading principal direction is (1, 2) normalised.
    /// </summary>
    private static Matrix CollinearData()
    {
        var values = new double[7, 2];
        for (var i = 0; i < 7; i++)
        {
            var t = i - 3.0;
            values[i, 0] = t;
            values[i, 1] = 2.0 * t;
        }

        return new Matrix(values);
    }

    /// <summary>Axis-aligned data whose variance is much larger along the first column.</summary>
    private static Matrix AxisAlignedData()
    {
        var random = new Random(17);
        var values = new double[60, 3];

        for (var i = 0; i < 60; i++)
        {
            values[i, 0] = random.NextDouble() * 100.0;
            values[i, 1] = random.NextDouble() * 10.0;
            values[i, 2] = random.NextDouble();
        }

        return new Matrix(values);
    }

    [TestMethod]
    public void Fit_CentersOnTheColumnMeans()
    {
        var data = new Matrix(new double[,] { { 1, 10 }, { 3, 20 }, { 5, 30 } });

        var pca = new PCA { NComponents = 1, Seed = 1 };
        pca.Fit(data);

        Assert.AreEqual(3.0, pca.Mean[0], 1e-12);
        Assert.AreEqual(20.0, pca.Mean[1], 1e-12);
    }

    [TestMethod]
    public void FirstComponent_AlignsWithTheDirectionOfGreatestVariance()
    {
        var pca = new PCA { NComponents = 1, Seed = 1 };
        pca.Fit(CollinearData());

        // Expected direction (1, 2)/sqrt(5). Sign is arbitrary for an eigenvector, so compare
        // the absolute value of the dot product with the unit expectation.
        var expectedX = 1.0 / Math.Sqrt(5.0);
        var expectedY = 2.0 / Math.Sqrt(5.0);

        var alignment = Math.Abs(pca.Components[0, 0] * expectedX + pca.Components[0, 1] * expectedY);

        Assert.AreEqual(1.0, alignment, 1e-6);
    }

    [TestMethod]
    public void FirstComponent_PicksTheHighestVarianceAxis()
    {
        var pca = new PCA { NComponents = 1, Seed = 1 };
        pca.Fit(AxisAlignedData());

        // Variance is dominated by column 0, so the leading component must point essentially
        // along it.
        Assert.AreEqual(1.0, Math.Abs(pca.Components[0, 0]), 1e-3);
    }

    [TestMethod]
    public void Components_AreUnitLength()
    {
        var pca = new PCA { NComponents = 3, Seed = 1 };
        pca.Fit(AxisAlignedData());

        for (var k = 0; k < pca.NComponents; k++)
        {
            var norm = 0.0;
            for (var j = 0; j < 3; j++)
            {
                norm += pca.Components[k, j] * pca.Components[k, j];
            }

            Assert.AreEqual(1.0, Math.Sqrt(norm), 1e-8, $"Component {k}");
        }
    }

    [TestMethod]
    public void Components_AreMutuallyOrthogonal()
    {
        var pca = new PCA { NComponents = 3, Seed = 1 };
        pca.Fit(AxisAlignedData());

        for (var a = 0; a < pca.NComponents; a++)
        {
            for (var b = a + 1; b < pca.NComponents; b++)
            {
                var dot = 0.0;
                for (var j = 0; j < 3; j++)
                {
                    dot += pca.Components[a, j] * pca.Components[b, j];
                }

                Assert.AreEqual(0.0, dot, 1e-6, $"Components {a} and {b}");
            }
        }
    }

    [TestMethod]
    public void ExplainedVariance_IsDescending()
    {
        var pca = new PCA { NComponents = 3, Seed = 1 };
        pca.Fit(AxisAlignedData());

        for (var k = 1; k < pca.NComponents; k++)
        {
            Assert.IsTrue(
                pca.ExplainedVariance[k] <= pca.ExplainedVariance[k - 1] + 1e-9,
                $"Component {k} explains more variance than {k - 1}: " +
                $"{pca.ExplainedVariance[k]} > {pca.ExplainedVariance[k - 1]}");
        }
    }

    [TestMethod]
    public void ExplainedVarianceRatio_SumsToOneWhenAllComponentsAreKept()
    {
        var pca = new PCA { NComponents = 3, Seed = 1 };
        pca.Fit(AxisAlignedData());

        var total = 0.0;
        for (var k = 0; k < pca.NComponents; k++)
        {
            total += pca.ExplainedVarianceRatio[k];
        }

        Assert.AreEqual(1.0, total, 1e-6);
    }

    [TestMethod]
    public void RankDeficientData_PutsAllVarianceInTheFirstComponent()
    {
        var pca = new PCA { NComponents = 2, Seed = 1 };
        pca.Fit(CollinearData());

        // The data lies on a line, so the second component explains nothing.
        Assert.AreEqual(1.0, pca.ExplainedVarianceRatio[0], 1e-6);
        Assert.AreEqual(0.0, pca.ExplainedVarianceRatio[1], 1e-6);
    }

    [TestMethod]
    public void Transform_ProjectsOntoTheRequestedNumberOfComponents()
    {
        var data = AxisAlignedData();

        var pca = new PCA { NComponents = 2, Seed = 1 };
        var reduced = pca.FitTransform(data);

        Assert.AreEqual(data.rowLength, reduced.rowLength);
        Assert.AreEqual(2, reduced.columnLength);
    }

    [TestMethod]
    public void Transform_ProducesUncorrelatedZeroMeanScores()
    {
        var data = AxisAlignedData();

        var pca = new PCA { NComponents = 3, Seed = 1 };
        var scores = pca.FitTransform(data);

        // Projections onto principal directions are centred and uncorrelated by construction.
        for (var k = 0; k < 3; k++)
        {
            var mean = 0.0;
            for (var i = 0; i < scores.rowLength; i++)
            {
                mean += scores.values[i, k];
            }

            Assert.AreEqual(0.0, mean / scores.rowLength, 1e-8, $"Score column {k} mean");
        }

        for (var a = 0; a < 3; a++)
        {
            for (var b = a + 1; b < 3; b++)
            {
                var covariance = 0.0;
                for (var i = 0; i < scores.rowLength; i++)
                {
                    covariance += scores.values[i, a] * scores.values[i, b];
                }

                Assert.AreEqual(0.0, covariance / scores.rowLength, 1e-6,
                    $"Score columns {a} and {b}");
            }
        }
    }

    [TestMethod]
    public void Transform_ReconstructsFullRankDataWhenAllComponentsAreKept()
    {
        var data = AxisAlignedData();

        var pca = new PCA { NComponents = 3, Seed = 1 };
        var scores = pca.FitTransform(data);

        // X ≈ mean + scores · components
        for (var i = 0; i < data.rowLength; i++)
        {
            for (var j = 0; j < 3; j++)
            {
                var reconstructed = pca.Mean[j];
                for (var k = 0; k < 3; k++)
                {
                    reconstructed += scores.values[i, k] * pca.Components[k, j];
                }

                Assert.AreEqual(data.values[i, j], reconstructed, 1e-8, $"Element ({i}, {j})");
            }
        }
    }

    [TestMethod]
    public void WideData_UsesTheDualPathAndStillProducesUnitComponents()
    {
        // Fewer samples than features drives the Gram-matrix branch.
        var random = new Random(5);
        var values = new double[4, 12];
        for (var i = 0; i < 4; i++)
        {
            for (var j = 0; j < 12; j++)
            {
                values[i, j] = random.NextDouble();
            }
        }

        var pca = new PCA { NComponents = 2, Seed = 1 };
        pca.Fit(new Matrix(values));

        Assert.AreEqual(2, pca.NComponents);

        for (var k = 0; k < 2; k++)
        {
            var norm = 0.0;
            for (var j = 0; j < 12; j++)
            {
                norm += pca.Components[k, j] * pca.Components[k, j];
            }

            Assert.AreEqual(1.0, Math.Sqrt(norm), 1e-8, $"Component {k}");
        }
    }

    [TestMethod]
    public void NComponents_IsClampedToTheAvailableRank()
    {
        var data = new Matrix(new double[,] { { 1, 2 }, { 3, 5 } });

        var pca = new PCA { NComponents = 10, Seed = 1 };
        pca.Fit(data);

        Assert.AreEqual(2, pca.NComponents);
    }

    [TestMethod]
    public void TransformBeforeFit_Throws()
    {
        var pca = new PCA { NComponents = 1 };

        Assert.ThrowsException<InvalidOperationException>(
            () => pca.Transform(AxisAlignedData()));
    }

    [TestMethod]
    public void ZeroComponents_Throws()
    {
        var pca = new PCA { NComponents = 0 };

        Assert.ThrowsException<ArgumentException>(() => pca.Fit(AxisAlignedData()));
    }

    [TestMethod]
    public void Clone_CarriesTheFittedState()
    {
        var pca = new PCA { NComponents = 2, Seed = 1 };
        pca.Fit(AxisAlignedData());

        var clone = (PCA)pca.Clone();

        Assert.AreEqual(pca.NComponents, clone.NComponents);
        Assert.AreEqual(pca.Components[0, 0], clone.Components[0, 0], 1e-15);
        Assert.AreEqual(pca.Mean[0], clone.Mean[0], 1e-15);
    }
}
