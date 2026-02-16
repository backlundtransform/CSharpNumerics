using CSharpNumerics.Objects;
using CSharpNumerics.Physics.Constants;
using Numerics.Objects;
using Xunit.Sdk;

namespace NumericsTests
{
    [TestClass]
    public class DifferentialEquationTest
    {
        [TestMethod]
        public void TestRungeKutta()
        {
            Func<(double y, double t), double> func = ((double t, double y) v) => Math.Tan(v.y) + 1;

            var result = func.RungeKutta(1, 1.1, 0.025, 1);
            Assert.IsTrue(Math.Round(result, 2) == 1.34);
        }

        [TestMethod]
        public void TestRungeKuttaMatrix()
        {
            Func<(double y, double t), double> func = ((double t, double y) v) => Math.Tan(v.y) + 1;

            var result = func.RungeKutta(1, 1.1, 0.025, 1, new Matrix(new double[,] { { 0, 0 }, { 2.0/3.0, 0 }}),new double[] { 1.0/4.0, 3.0/4.0 }, new double[] { 0.0, 2.0/3.0 });
            Assert.IsTrue(Math.Round(result, 2) == 1.34);
        }

        [TestMethod]
        public void TestTrapets()
        {
            Func<(double y, double t), double> func = ((double t, double y) v) => Math.Tan(v.y) + 1;

            var result = func.TrapezoidalRule(1, 1.1, 0.00025, 1);
            Assert.IsTrue(Math.Round(result,2) == 1.34);
        }

        [TestMethod]
        public void TestLinearSystem()
        {
            var matrix = new Matrix(new double[,] { { 1, -2, 3 }, { -1, 1, -2 }, { 2, -1, -1 } });

            var vector = new Vector(7, -5, 4);

            var result = matrix.LinearSystemSolver(vector);
            Assert.IsTrue(result.x == 2);
            Assert.IsTrue(result.y == -1);
            Assert.IsTrue(result.z== 1);
        }


        [TestMethod]
        public void TestLinearSystemN()
        {
            var matrix = new Matrix(new double[,] { { 1, -2, 3 }, { -1, 1, -2 }, { 2, -1, -1 } });

            var vector = new VectorN(new double[] { 7, -5, 4 });

            var result = matrix.LinearSystemSolver(vector);
            Assert.IsTrue(result[0] == 2);
            Assert.IsTrue(result[1] == -1);
            Assert.IsTrue(result[2] == 1);
        }

        [TestMethod]
        public void TestLinearSystemRange()
        {
            var matrix = new Matrix(new double[,] { { 1, -2, 3 }, { -1, 1, -2 }, { 2, -1, -1 } });

            var vector = new List<double>() {  7 ,-5 ,  4  };

            var result = matrix.LinearSystemSolver(vector);
            Assert.IsTrue(result[0] == 2);
            Assert.IsTrue(result[1] == -1);
            Assert.IsTrue(result[2] == 1);
        }
        [TestMethod]
        public void TestGaussElimination()
        {
            var matrix = new Matrix(new double[,] { { 1, -2, 3 }, { -1, 1, -2 }, { 2, -1, -1 } });

            var vector = new Vector(7, -5, 4);

            var result = matrix.GaussElimination(vector);
            Assert.IsTrue(result.x == 2);
            Assert.IsTrue(result.y == -1);
            Assert.IsTrue(result.z == 1);
        }
   

        [TestMethod]
        public void TestGaussEliminationRange()
        {
            var matrix = new Matrix(new double[,] { { 1, -2, 3 }, { -1, 1, -2 }, { 2, -1, -1 } });

            var vector = new List<double>() { { 7 }, { -5 }, { 4 } };

            var result = matrix.GaussElimination(vector);
            Assert.IsTrue(result[0] == 2);
            Assert.IsTrue(result[1] == -1);
            Assert.IsTrue(result[2] == 1);
        }


        [TestMethod]
        public void TestEigenVector()
        {
            var matrix = new Matrix(new double[,] { { 3, -4 }, {4, -7 } });
            var result = matrix.EigenVector(1);
            Assert.IsTrue(Math.Round(result[0], 1) ==2);
            Assert.IsTrue(Math.Round(result[1], 1) == 1);
            var result2 = matrix.EigenVector(-5);
            Assert.IsTrue(Math.Round(result2[0], 1) == 1);
            Assert.IsTrue(Math.Round(result2[1], 1) == 2);

        }


        [TestMethod]
        public void TestDominantEigenVector()
        {
           var matrix = new Matrix(new double[,] { { 3, -4 }, { 4, -7 } });
           var result = matrix.DominantEigenVector();

            Assert.IsTrue(Math.Round(result[0], 1) == 1);
            Assert.IsTrue(Math.Round(result[1], 1) == 2);

            result =  matrix.Inverse().DominantEigenVector();
            Assert.IsTrue(Math.Round(result[0], 1) == 2);
            Assert.IsTrue(Math.Round(result[1], 1) == 1);

        }

        [TestMethod]
        public void TestEigenvalues()
        {

          
            var matrix = new Matrix(new double[,] { { 3, -1 }, { 4, -2 } });

            var result = matrix.EigenValues();
            Assert.IsTrue(result.First() == 2);
            Assert.IsTrue(result.Last() == -1);
            matrix = new Matrix(new double[,] { { 2, 1 }, { 1, 2 } });
            result = matrix.EigenValues();
            Assert.IsTrue(result.First() == 3);
            Assert.IsTrue(result.Last() == 1);


           

        }

      



        [TestMethod]
        public void TestSolveOde()
        {
            var matrix = new Matrix(new double[,] { { 3, -4 }, { 4, -7 } });
            var result = matrix.OdeSolver(1);

            Assert.IsTrue(Math.Round(result[0](2),5) == Math.Round(2.0 /3.0* Math.Exp(2)+ 1.0 / 3.0 * Math.Exp(-10),5));
            Assert.IsTrue(Math.Round(result[1](2),5) == Math.Round(1.0 / 3.0 * Math.Exp(2) + 2.0 / 3.0 * Math.Exp(-10), 5));
                    
     
        }

        #region Vector ODE Solvers

        [TestMethod]
        public void RungeKutta_Vector_ExponentialDecay()
        {
            // dy/dt = -y, y(0) = (1,2,3) → y(t) = (e^-t, 2e^-t, 3e^-t)
            Func<(double t, Vector y), Vector> func = v => -1.0 * v.y;
            var y0 = new Vector(1, 2, 3);

            var result = func.RungeKutta(0, 1, 0.001, y0);

            Assert.AreEqual(Math.Exp(-1), result.x, 1e-6);
            Assert.AreEqual(2 * Math.Exp(-1), result.y, 1e-6);
            Assert.AreEqual(3 * Math.Exp(-1), result.z, 1e-6);
        }

        [TestMethod]
        public void EulerMethod_Vector_ConstantDerivative()
        {
            // dy/dt = (1, 0, -9.8), y(0) = (0,0,100) → y(t) = (t, 0, 100 - 9.8t)
            Func<(double t, Vector y), Vector> func = _ => new Vector(1, 0, -9.8);
            var y0 = new Vector(0, 0, 100);

            var result = func.EulerMethod(0, 5, 0.001, y0);

            Assert.AreEqual(5, result.x, 1e-3);
            Assert.AreEqual(0, result.y, 1e-10);
            Assert.AreEqual(100 - 9.8 * 5, result.z, 1e-3);
        }

        [TestMethod]
        public void RungeKuttaTrajectory_Vector_ReturnsCorrectLength()
        {
            Func<(double t, Vector y), Vector> func = v => -1.0 * v.y;
            var y0 = new Vector(1, 0, 0);

            var traj = func.RungeKuttaTrajectory(0, 1, 0.1, y0);

            Assert.IsTrue(traj.Count >= 10);
            Assert.AreEqual(0, traj[0].t, 1e-10);
            Assert.AreEqual(1, traj[0].y.x, 1e-10);
        }

        [TestMethod]
        public void RungeKutta_Vector_CircularMotion()
        {
            // dy/dt = (-y.y, y.x, 0) → circular motion with ω=1
            // y(0) = (1, 0, 0) → y(t) = (cos t, sin t, 0)
            Func<(double t, Vector y), Vector> func = v => new Vector(-v.y.y, v.y.x, 0);
            var y0 = new Vector(1, 0, 0);

            double t = 2 * Math.PI; // full circle
            var result = func.RungeKutta(0, t, 0.001, y0);

            Assert.AreEqual(1, result.x, 1e-4);
            Assert.AreEqual(0, result.y, 1e-4);
        }

        #endregion

        #region System ODE Solvers (double[])

        [TestMethod]
        public void RungeKutta_System_SimpleHarmonic()
        {
            // x'' = -x → system: y = [x, v], y' = [v, -x]
            // y(0) = [1, 0] → x(t) = cos(t)
            Func<(double t, double[] y), double[]> func = v =>
                [v.y[1], -v.y[0]];
            var y0 = new double[] { 1, 0 };

            var result = func.RungeKutta(0, Math.PI, 0.001, y0);

            Assert.AreEqual(Math.Cos(Math.PI), result[0], 1e-5); // x(π) = -1
            Assert.AreEqual(Math.Sin(Math.PI), -result[1], 1e-5); // v(π) ≈ 0
        }

        [TestMethod]
        public void RungeKutta_System_FreeFallDynamics()
        {
            // State: [x, y, z, vx, vy, vz]
            // Derivative: [vx, vy, vz, 0, 0, -g]
            double g = PhysicsConstants.GravitationalAcceleration;
            Func<(double t, double[] y), double[]> dynamics = v =>
                [v.y[3], v.y[4], v.y[5], 0, 0, -g];

            // Launch at 45° with 20 m/s
            double speed = 20;
            double angle = Math.PI / 4;
            var y0 = new double[] { 0, 0, 0, speed * Math.Cos(angle), 0, speed * Math.Sin(angle) };

            double T = 2 * speed * Math.Sin(angle) / g; // time of flight
            var result = dynamics.RungeKutta(0, T, 0.001, y0);

            // Should return to z ≈ 0
            Assert.AreEqual(0, result[2], 0.1);
            // Range should match analytical: R = v²sin(2θ)/g
            double expectedRange = speed * speed * Math.Sin(2 * angle) / g;
            Assert.AreEqual(expectedRange, result[0], 0.1);
        }

        [TestMethod]
        public void EulerMethod_System_LinearODE()
        {
            // y' = y, y(0) = [1, 2] → y(t) = [e^t, 2e^t]
            Func<(double t, double[] y), double[]> func = v =>
                [v.y[0], v.y[1]];
            var y0 = new double[] { 1, 2 };

            var result = func.EulerMethod(0, 1, 0.0001, y0);

            Assert.AreEqual(Math.E, result[0], 0.01);
            Assert.AreEqual(2 * Math.E, result[1], 0.02);
        }

        [TestMethod]
        public void RungeKuttaTrajectory_System_ConservesEnergy()
        {
            // Simple harmonic oscillator: E = 0.5*(v² + x²) should be constant
            Func<(double t, double[] y), double[]> func = v =>
                [v.y[1], -v.y[0]];
            var y0 = new double[] { 1, 0 }; // E = 0.5

            var traj = func.RungeKuttaTrajectory(0, 10, 0.01, y0);

            double E0 = 0.5 * (y0[0] * y0[0] + y0[1] * y0[1]);
            foreach (var (t, y) in traj)
            {
                double E = 0.5 * (y[0] * y[0] + y[1] * y[1]);
                Assert.AreEqual(E0, E, 1e-6);
            }
        }

        [TestMethod]
        public void RungeKuttaTrajectory_System_ReturnsTimeSteps()
        {
            Func<(double t, double[] y), double[]> func = v =>
                [v.y[0]];
            var y0 = new double[] { 1 };

            var traj = func.RungeKuttaTrajectory(0, 1, 0.1, y0);

            Assert.IsTrue(traj.Count >= 10);
            Assert.AreEqual(0, traj[0].t, 1e-10);
            Assert.AreEqual(1, traj.Last().t, 1e-10);
        }

        [TestMethod]
        public void RungeKutta_System_OrbitalMotion()
        {
            // Circular orbit: a = -GM/r³ * r
            // State: [x, y, vx, vy] (2D orbit)
            double GM = PhysicsConstants.GravitationalConstant * PhysicsConstants.EarthMass;
            double R = 7e6; // 7000 km
            double vOrb = Math.Sqrt(GM / R);

            Func<(double t, double[] y), double[]> dynamics = v =>
            {
                double x = v.y[0], y2 = v.y[1], vx = v.y[2], vy = v.y[3];
                double r3 = Math.Pow(x * x + y2 * y2, 1.5);
                return [vx, vy, -GM * x / r3, -GM * y2 / r3];
            };

            var y0 = new double[] { R, 0, 0, vOrb };
            double T = 2 * Math.PI * R / vOrb; // orbital period

            var result = dynamics.RungeKutta(0, T, 1, y0);

            // Should return close to start after one full orbit
            Assert.AreEqual(R, result[0], R * 1e-4);
            Assert.AreEqual(0, result[1], R * 1e-4);
        }

        #endregion

        #region VectorN ODE Solvers

        [TestMethod]
        public void RungeKutta_VectorN_SimpleHarmonic()
        {
            // x'' = -x → system: y = [x, v], y' = [v, -x]
            Func<(double t, VectorN y), VectorN> func = v =>
                new VectorN([v.y[1], -v.y[0]]);
            var y0 = new VectorN([1, 0]);

            VectorN result = func.RungeKutta(0, Math.PI, 0.001, y0);

            Assert.AreEqual(-1, result[0], 1e-5); // cos(π) = -1
            Assert.AreEqual(0, result[1], 1e-4);  // sin(π) ≈ 0
        }

        [TestMethod]
        public void RungeKutta_VectorN_FreeFall()
        {
            double g = PhysicsConstants.GravitationalAcceleration;
            Func<(double t, VectorN y), VectorN> dynamics = v =>
                new VectorN([v.y[3], v.y[4], v.y[5], 0, 0, -g]);

            double speed = 20;
            double angle = Math.PI / 4;
            var y0 = new VectorN([0, 0, 0, speed * Math.Cos(angle), 0, speed * Math.Sin(angle)]);

            double T = 2 * speed * Math.Sin(angle) / g;
            VectorN result = dynamics.RungeKutta(0, T, 0.001, y0);

            Assert.AreEqual(0, result[2], 0.1);
            double expectedRange = speed * speed * Math.Sin(2 * angle) / g;
            Assert.AreEqual(expectedRange, result[0], 0.1);
        }

        [TestMethod]
        public void EulerMethod_VectorN_LinearODE()
        {
            Func<(double t, VectorN y), VectorN> func = v =>
                new VectorN([v.y[0], v.y[1]]);
            var y0 = new VectorN([1, 2]);

            VectorN result = func.EulerMethod(0, 1, 0.0001, y0);

            Assert.AreEqual(Math.E, result[0], 0.01);
            Assert.AreEqual(2 * Math.E, result[1], 0.02);
        }

        [TestMethod]
        public void RungeKuttaTrajectory_VectorN_ConservesEnergy()
        {
            Func<(double t, VectorN y), VectorN> func = v =>
                new VectorN([v.y[1], -v.y[0]]);
            var y0 = new VectorN([1, 0]);

            var traj = func.RungeKuttaTrajectory(0, 10, 0.01, y0);

            double E0 = 0.5 * y0.Dot(y0);
            foreach (var (t, y) in traj)
            {
                double E = 0.5 * y.Dot(y);
                Assert.AreEqual(E0, E, 1e-6);
            }
        }

        [TestMethod]
        public void RungeKutta_VectorN_OrbitalMotion()
        {
            double GM = PhysicsConstants.GravitationalConstant * PhysicsConstants.EarthMass;
            double R = 7e6;
            double vOrb = Math.Sqrt(GM / R);

            Func<(double t, VectorN y), VectorN> dynamics = v =>
            {
                double x = v.y[0], y2 = v.y[1], vx = v.y[2], vy = v.y[3];
                double r3 = Math.Pow(x * x + y2 * y2, 1.5);
                return new VectorN([vx, vy, -GM * x / r3, -GM * y2 / r3]);
            };

            var y0 = new VectorN([R, 0, 0, vOrb]);
            double T = 2 * Math.PI * R / vOrb;

            VectorN result = dynamics.RungeKutta(0, T, 1, y0);

            Assert.AreEqual(R, result[0], R * 1e-4);
            Assert.AreEqual(0, result[1], R * 1e-4);
        }

        [TestMethod]
        public void RungeKutta_VectorN_UsesVectorNOperations()
        {
            // Verify that Norm, indexer, etc. work naturally in the ODE context
            Func<(double t, VectorN y), VectorN> func = v =>
            {
                double norm = v.y.Norm();
                // Damped: dy/dt = -y/|y| (normalizes then negates)
                return norm > 0 ? (-1.0 / norm) * v.y : new VectorN(v.y.Length);
            };

            var y0 = new VectorN([3, 4]); // |y0| = 5
            VectorN result = func.RungeKutta(0, 2, 0.001, y0);

            // After t=2 with unit speed decay, |y| should be |y0| - t = 5 - 2 = 3
            Assert.AreEqual(3.0, result.Norm(), 0.01);
        }

        #endregion

        #region Semi-Implicit Euler

        [TestMethod]
        public void SemiImplicitEuler_SimpleHarmonic_Stable()
        {
            // x'' = -x → state = [x, v], derivative = [v, -x]
            // Semi-implicit Euler should be stable (symplectic) even over long time
            Func<(double t, VectorN y), VectorN> func = v =>
                new VectorN([v.y[1], -v.y[0]]);
            var y0 = new VectorN([1, 0]);

            var result = func.SemiImplicitEuler(0, 100, 0.01, y0);

            // Energy should be approximately conserved: E = 0.5*(x² + v²) ≈ 0.5
            double E = 0.5 * (result[0] * result[0] + result[1] * result[1]);
            Assert.AreEqual(0.5, E, 0.05);
        }

        [TestMethod]
        public void SemiImplicitEuler_FreeFall()
        {
            double g = PhysicsConstants.GravitationalAcceleration;
            Func<(double t, VectorN y), VectorN> func = v =>
                new VectorN([v.y[1], -g]);
            var y0 = new VectorN([100, 0]); // height=100, v=0

            var result = func.SemiImplicitEuler(0, 2, 0.001, y0);

            double expectedH = 100 - 0.5 * g * 4;
            Assert.AreEqual(expectedH, result[0], 0.1);
        }

        [TestMethod]
        public void SemiImplicitEuler_DoubleArray_Convenience()
        {
            Func<(double t, double[] y), double[]> func = v =>
                [v.y[1], -v.y[0]];
            var y0 = new double[] { 1, 0 };

            var result = func.SemiImplicitEuler(0, Math.PI, 0.001, y0);

            Assert.AreEqual(-1, result[0], 0.05);
        }

        [TestMethod]
        public void SemiImplicitEulerTrajectory_ReturnsSteps()
        {
            Func<(double t, VectorN y), VectorN> func = v =>
                new VectorN([v.y[1], -v.y[0]]);
            var y0 = new VectorN([1, 0]);

            var traj = func.SemiImplicitEulerTrajectory(0, 1, 0.1, y0);

            Assert.IsTrue(traj.Count >= 10);
            Assert.AreEqual(0, traj[0].t, 1e-10);
            Assert.AreEqual(1, traj[0].y[0], 1e-10);
        }

        [TestMethod]
        public void SemiImplicitEuler_BetterThanExplicitEuler_Energy()
        {
            // After 100 time units, semi-implicit should conserve energy much better
            Func<(double t, VectorN y), VectorN> func = v =>
                new VectorN([v.y[1], -v.y[0]]);
            var y0 = new VectorN([1, 0]);

            var semiImplicit = func.SemiImplicitEuler(0, 50, 0.01, y0);
            var euler = func.EulerMethod(0, 50, 0.01, y0);

            double E0 = 0.5;
            double eSemi = 0.5 * (semiImplicit[0] * semiImplicit[0] + semiImplicit[1] * semiImplicit[1]);
            double eEuler = 0.5 * (euler[0] * euler[0] + euler[1] * euler[1]);

            Assert.IsTrue(Math.Abs(eSemi - E0) < Math.Abs(eEuler - E0));
        }

        #endregion

        #region Velocity Verlet

        [TestMethod]
        public void VelocityVerlet_SimpleHarmonic_ExcellentEnergy()
        {
            // x'' = -x, acceleration function returns just [-x]
            Func<(double t, VectorN y), VectorN> accel = v =>
                new VectorN([-v.y[0]]);
            var y0 = new VectorN([1, 0]);

            var result = accel.VelocityVerlet(0, 100, 0.01, y0);

            double E = 0.5 * (result[0] * result[0] + result[1] * result[1]);
            Assert.AreEqual(0.5, E, 1e-4);
        }

        [TestMethod]
        public void VelocityVerlet_FreeFall_SecondOrderAccurate()
        {
            double g = PhysicsConstants.GravitationalAcceleration;
            // State: [height, velocity], acceleration = [-g]
            Func<(double t, VectorN y), VectorN> accel = _ =>
                new VectorN([-g]);
            var y0 = new VectorN([100, 0]);

            var result = accel.VelocityVerlet(0, 2, 0.1, y0);

            double expectedH = 100 - 0.5 * g * 4;
            double expectedV = -g * 2;
            Assert.AreEqual(expectedH, result[0], 1e-10); // exact for constant accel
            Assert.AreEqual(expectedV, result[1], 1e-10);
        }

        [TestMethod]
        public void VelocityVerlet_Orbital_ReturnsToStart()
        {
            double GM = PhysicsConstants.GravitationalConstant * PhysicsConstants.EarthMass;
            double R = 7e6;
            double vOrb = Math.Sqrt(GM / R);

            Func<(double t, VectorN y), VectorN> accel = v =>
            {
                double x = v.y[0], y2 = v.y[1];
                double r3 = Math.Pow(x * x + y2 * y2, 1.5);
                return new VectorN([-GM * x / r3, -GM * y2 / r3]);
            };

            var y0 = new VectorN([R, 0, 0, vOrb]);
            double T = 2 * Math.PI * R / vOrb;

            var result = accel.VelocityVerlet(0, T, 0.5, y0);

            Assert.AreEqual(R, result[0], R * 1e-4);
            Assert.AreEqual(0, result[1], R * 1e-4);
        }

        [TestMethod]
        public void VelocityVerlet_DoubleArray_Convenience()
        {
            Func<(double t, double[] y), double[]> accel = v =>
                [-v.y[0]];
            var y0 = new double[] { 1, 0 };

            var result = accel.VelocityVerlet(0, Math.PI, 0.01, y0);

            Assert.AreEqual(-1, result[0], 0.01);
        }

        [TestMethod]
        public void VelocityVerletTrajectory_ConservesEnergy()
        {
            Func<(double t, VectorN y), VectorN> accel = v =>
                new VectorN([-v.y[0]]);
            var y0 = new VectorN([1, 0]);

            var traj = accel.VelocityVerletTrajectory(0, 20, 0.01, y0);

            double E0 = 0.5;
            foreach (var (t, y) in traj)
            {
                double E = 0.5 * (y[0] * y[0] + y[1] * y[1]);
                Assert.AreEqual(E0, E, 1e-4);
            }
        }

        [TestMethod]
        public void VelocityVerlet_SymplecticEnergyConservation()
        {
            // Verlet is symplectic: energy oscillates around true value instead of drifting
            // Over 10000 time units, RK4 energy drifts monotonically while Verlet stays bounded
            Func<(double t, VectorN y), VectorN> accel = v =>
                new VectorN([-v.y[0]]);
            Func<(double t, VectorN y), VectorN> rk4Func = v =>
                new VectorN([v.y[1], -v.y[0]]);
            var y0 = new VectorN([1, 0]);

            var verlet = accel.VelocityVerlet(0, 10000, 0.5, y0);
            var rk4 = rk4Func.RungeKutta(0, 10000, 0.5, y0);

            double E0 = 0.5;
            double eVerlet = 0.5 * (verlet[0] * verlet[0] + verlet[1] * verlet[1]);
            double eRk4 = 0.5 * (rk4[0] * rk4[0] + rk4[1] * rk4[1]);

            Assert.IsTrue(Math.Abs(eVerlet - E0) < Math.Abs(eRk4 - E0));
        }

        #endregion
    }
}
