using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
namespace NumericsTests
{
    [TestClass]
    public class TrigonometryTests
    {
    
        [TestMethod]
        public void TestVelocityToAcceleration()
        {
            var result = 180d.DegreeToRadians();
            Assert.IsTrue(result == Math.PI);
        }
    }
}
