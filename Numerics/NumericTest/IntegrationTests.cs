using CSharpNumerics.Numerics;


namespace NumericTest
{
    [TestClass]
    public class IntegrationTests
    {
        #region Simpson's 1/3 Rule

        [TestMethod]
        public void IntegrateSimpson_Polynomial_ExactForCubic()
        {
            // Simpson's rule is exact for polynomials up to degree 3
            // ∫₀¹ x³ dx = 1/4
            Func<double, double> func = x => x * x * x;
            double result = func.IntegrateSimpson(0, 1, 2);
            Assert.AreEqual(0.25, result, 1e-14);
        }

        [TestMethod]
        public void IntegrateSimpson_Exponential_HighAccuracy()
        {
            // ∫₀⁵ eˣ dx = e⁵ - 1
            Func<double, double> func = x => Math.Exp(x);
            double expected = Math.Exp(5) - 1;
            double result = func.IntegrateSimpson(0, 5, 1000);
            Assert.AreEqual(expected, result, 1e-8);
        }

        [TestMethod]
        public void IntegrateSimpson_Sin_ZeroToPI()
        {
            // ∫₀π sin(x) dx = 2
            Func<double, double> func = x => Math.Sin(x);
            double result = func.IntegrateSimpson(0, Math.PI, 200);
            Assert.AreEqual(2.0, result, 1e-8);
        }

        #endregion

        #region Simpson's 3/8 Rule

        [TestMethod]
        public void IntegrateSimpson38_Polynomial()
        {
            // ∫₀¹ x³ dx = 1/4 — also exact for cubics
            Func<double, double> func = x => x * x * x;
            double result = func.IntegrateSimpson38(0, 1, 3);
            Assert.AreEqual(0.25, result, 1e-14);
        }

        [TestMethod]
        public void IntegrateSimpson38_Exponential()
        {
            Func<double, double> func = x => Math.Exp(x);
            double expected = Math.Exp(5) - 1;
            double result = func.IntegrateSimpson38(0, 5, 999);
            Assert.AreEqual(expected, result, 1e-7);
        }

        #endregion

        #region Gauss-Legendre Quadrature

        [TestMethod]
        public void IntegrateGaussLegendre_Polynomial_NearExact()
        {
            // 5-point GL is exact for polynomials up to degree 9
            // ∫₀¹ x⁹ dx = 1/10
            Func<double, double> func = x => Math.Pow(x, 9);
            double result = func.IntegrateGaussLegendre(0, 1, 1);
            Assert.AreEqual(0.1, result, 1e-14);
        }

        [TestMethod]
        public void IntegrateGaussLegendre_Exponential_HighAccuracy()
        {
            Func<double, double> func = x => Math.Exp(x);
            double expected = Math.Exp(5) - 1;
            double result = func.IntegrateGaussLegendre(0, 5, 10);
            Assert.AreEqual(expected, result, 1e-12);
        }

        [TestMethod]
        public void IntegrateGaussLegendre_Sin_HighAccuracy()
        {
            // ∫₀π sin(x) dx = 2
            Func<double, double> func = x => Math.Sin(x);
            double result = func.IntegrateGaussLegendre(0, Math.PI, 5);
            Assert.AreEqual(2.0, result, 1e-12);
        }

        #endregion

        #region Romberg Integration

        [TestMethod]
        public void IntegrateRomberg_Exponential_HighAccuracy()
        {
            Func<double, double> func = x => Math.Exp(x);
            double expected = Math.Exp(5) - 1;
            double result = func.IntegrateRomberg(0, 5, 10);
            Assert.AreEqual(expected, result, 1e-10);
        }

        [TestMethod]
        public void IntegrateRomberg_Polynomial_Exact()
        {
            // ∫₀¹ (3x² + 2x + 1) dx = 1 + 1 + 1 = 3
            Func<double, double> func = x => 3 * x * x + 2 * x + 1;
            double result = func.IntegrateRomberg(0, 1, 5);
            Assert.AreEqual(3.0, result, 1e-12);
        }

        [TestMethod]
        public void IntegrateRomberg_Sin()
        {
            // ∫₀π sin(x) dx = 2
            Func<double, double> func = x => Math.Sin(x);
            double result = func.IntegrateRomberg(0, Math.PI, 8);
            Assert.AreEqual(2.0, result, 1e-12);
        }

        #endregion

        #region Adaptive Simpson

        [TestMethod]
        public void IntegrateAdaptive_Exponential()
        {
            Func<double, double> func = x => Math.Exp(x);
            double expected = Math.Exp(5) - 1;
            double result = func.IntegrateAdaptive(0, 5, 1e-12);
            Assert.AreEqual(expected, result, 1e-10);
        }

        [TestMethod]
        public void IntegrateAdaptive_Sqrt_NonSmooth()
        {
            // ∫₀¹ √x dx = 2/3 — tests adaptive subdivision on a function with infinite derivative at 0
            Func<double, double> func = x => Math.Sqrt(x);
            double result = func.IntegrateAdaptive(0, 1, 1e-8);
            Assert.AreEqual(2.0 / 3.0, result, 1e-7);
        }

        [TestMethod]
        public void IntegrateAdaptive_Sin()
        {
            Func<double, double> func = x => Math.Sin(x);
            double result = func.IntegrateAdaptive(0, Math.PI, 1e-12);
            Assert.AreEqual(2.0, result, 1e-10);
        }

        #endregion

        #region Method Comparison

        [TestMethod]
        public void AllMethods_AgreeOn_GaussianIntegral()
        {
            // ∫₋₃³ e^(-x²) dx = √π * erf(3) ≈ 1.772414697...
            Func<double, double> func = x => Math.Exp(-x * x);
            double expected = 1.7724146966; // √π * erf(3)

            double simpson = func.IntegrateSimpson(-3, 3, 10000);
            double gauss = func.IntegrateGaussLegendre(-3, 3, 20);
            double romberg = func.IntegrateRomberg(-3, 3, 10);
            double adaptive = func.IntegrateAdaptive(-3, 3, 1e-12);

            Assert.AreEqual(expected, simpson, 1e-6, $"Simpson: {simpson}");
            Assert.AreEqual(expected, gauss, 1e-10, $"Gauss-Legendre: {gauss}");
            Assert.AreEqual(expected, romberg, 1e-10, $"Romberg: {romberg}");
            Assert.AreEqual(expected, adaptive, 1e-10, $"Adaptive: {adaptive}");
        }

        #endregion
    }
}
