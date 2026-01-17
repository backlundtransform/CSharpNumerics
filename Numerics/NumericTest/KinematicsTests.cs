using CSharpNumerics.Physics;
using Numerics.Objects;


namespace NumericTest
{
    [TestClass]
    public class KinematicsTests
    {
        [TestMethod]
        public void FreeFallVelocity_From10Meters_IsCorrect()
        {
            var v = 10.0.FreeFallVelocity();
            Assert.AreEqual(Math.Sqrt(2 * 9.80665 * 10), v, 6);
        }

        [TestMethod]
        public void PositionFromConstantAcceleration_Scalar()
        {
            var s = 2.0.PositionFromConstantAcceleration(
                time: 3,
                initialVelocity: 1,
                initialPosition: 0);

            // s = 0 + 1*3 + 0.5*2*9 = 12
            Assert.AreEqual(12.0, s, 6);
        }

        [TestMethod]
        public void PositionFromConstantAcceleration_Vector()
        {
            var a = new Vector(1, 0, 0);
            var v0 = new Vector(0, 0, 0);
            var s0 = new Vector(0, 0, 0);

            var s = a.PositionFromConstantAcceleration(2, v0, s0);

            // s = 0 + 0 + 0.5 * 1 * 4 = (2,0,0)
            Assert.AreEqual(2, s.x, 6);
        }

        [TestMethod]
        public void VelocityFromConstantAcceleration_Vector()
        {
            var a = new Vector(1, 0, 0);
            var v0 = new Vector(0, 0, 0);
            var v = a.VelocityFromConstantAcceleration(3, v0);
            // v = 0 + 1 * 3 = (3,0,0)
            Assert.AreEqual(3, v.x, 6);
        }
        [TestMethod]
        public void VelocityFromConstantAcceleration_Scalar()
        {
            var v = 2.0.VelocityFromConstantAcceleration(
                time: 4,
                initialVelocity: 1);
            // v = 1 + 2*4 = 9
            Assert.AreEqual(9.0, v, 6);
        }
        [TestMethod]
        public void PositionFromConstantVelocity_Scalar()
        {
            var s = 3.0.PositionFromConstantVelocity(
                time: 5,
                initialPosition: 2);
            // s = 2 + 3*5 = 17
            Assert.AreEqual(17.0, s, 6);
        }
        [TestMethod]
        public void PositionFromConstantVelocity_Vector()
        {
            var v = new Vector(2, 0, 0);
            var s0 = new Vector(1, 0, 0);
            var s = v.PositionFromConstantVelocity(4, s0);
            // s = (1,0,0) + (2,0,0)*4 = (9,0,0)
            Assert.AreEqual(9, s.x, 6);
        }
    }
}
