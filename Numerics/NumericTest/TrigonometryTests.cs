using Xunit.Sdk;

namespace NumericsTests
{
    [TestClass]
    public class TrigonometryTests
    {
    
        [TestMethod]
        public void TestDegreeToRadians()
        {
            var result = 180d.DegreeToRadians();
            Assert.IsTrue(result == Math.PI);
        }
    }
}
