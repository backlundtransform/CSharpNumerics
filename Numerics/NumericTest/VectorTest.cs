using Xunit.Sdk;
using CSharpNumerics.Numerics.Objects;


namespace NumericsTests
{
    [TestClass]
    public class VectorTest
    {
        [TestMethod]
        public void TestVectorMultipliedByScalar()
        {
            var v = new Vector(2,6,1);
            v = 2 * v;
            Assert.IsTrue(v.x == 4);
            Assert.IsTrue(v.y == 12);
            Assert.IsTrue(v.z == 2);
          
        }


        [TestMethod]
        public void TestMagnitude()
        {
            var v = new Vector(2, 2, 1);
            var length = v.GetMagnitude();
            Assert.IsTrue(length == 3);
        }

        [TestMethod]
        public void TestUnitVector()
        {
            var v = new Vector(2, 2, 1);
            var v2 =v.GetUnitVector();
            Assert.IsTrue(v2.x == (double)2 /3);
            Assert.IsTrue(v2.y == (double)2 /3);
            Assert.IsTrue(v2.z == (double)1 /3);
            Assert.IsTrue(v2.GetMagnitude() == 1);

        }


        [TestMethod]
        public void TestVectorAddition()
        {
            var v = new Vector(2, 3, 0);
            var v2 = new Vector(3, 1, 0);
            var v3 = v+v2;
            Assert.IsTrue(v3.x ==5);
            Assert.IsTrue(v3.y ==4);
            Assert.IsTrue(v3.z ==0);
        }

        [TestMethod]
        public void TestVectorSubstraction()
        {
            var v = new Vector(2, 3, 0);
            var v2 = new Vector(3, 1, 0);
            var v3 = v-v2;
            Assert.IsTrue(v3.x ==-1);
            Assert.IsTrue(v3.y == 2);
            Assert.IsTrue(v3.z == 0);
        }

        [TestMethod]
        public void TestDotProduct()
        {
            var v = new Vector(5, 3, 0);
            var v2 = new Vector(2, 6, 0);
            var skalar = v.Dot(v2);
            Assert.IsTrue(skalar == 28);
       
        }

        [TestMethod]
        public void TestVectorProduct()
        {
            var v = new Vector(1, 3, 2);
            var v2 = new Vector(-2, 1, 0);
            var v3 = v.Cross(v2);
            Assert.IsTrue(v3.x == -2);
            Assert.IsTrue(v3.y == -4);
            Assert.IsTrue(v3.z == 7);
        }

        [TestMethod]
        public void TestPoints()
        {
            var v = new Vector((5,4,0),(3,7,0)); 
            Assert.IsTrue(v.x == -2);
            Assert.IsTrue(v.y == 3);
            Assert.IsTrue(v.z == 0);
        }

        [TestMethod]
        public void TestAngle()
        {
            var v = new Vector(2, 2, 0);
            var v2 = new Vector(0, 3, 0);
            var angle =v.GetAngle(v2);
            Assert.IsTrue(Math.Round(angle,2) == Math.Round(Math.PI/4,2));
        }


        [TestMethod]
        public void TestProjection()
        {
            var a= new Vector(-7, 4, 1);
            var b = new Vector(3, 2, -1);
            var v3 = a.Projection(b);
            var v4 = (-7.0 / 33.0)*a;
            Assert.IsTrue(Math.Round(v3.x,4) == Math.Round(v4.x,4));
            Assert.IsTrue(Math.Round(v3.y, 4) == Math.Round(v4.y,4));
            Assert.IsTrue(Math.Round(v3.z,4) == Math.Round(v4.z,4));
        }

        [TestMethod]
        public void TestReflection()
        {
            var v= new Vector(2, 4, 2);
            var v2= new Vector(-6, 2, 4);
            var v3 =v2.Reflection(v);
            var v4 = (-1.0 / 7.0) * new Vector(20, 26, 10);
            Assert.IsTrue(Math.Round(v3.x, 4) == Math.Round(v4.x, 4));
            Assert.IsTrue(Math.Round(v3.y, 4) == Math.Round(v4.y, 4));
            Assert.IsTrue(Math.Round(v3.z, 4) == Math.Round(v4.z, 4));
        }

        [TestMethod]
        public void TestSphericalCoordinates()
        {
            var v = new Vector(2, 4, 2);
            var v2 = v.ToSphericalCoordinates();
            var v3 = Vector.FromSphericalCoordinates(v2.x, v2.y, v2.z);
          
            Assert.IsTrue(Math.Round(v3.x) == Math.Round(v.x));
            Assert.IsTrue(Math.Round(v3.y) == Math.Round(v.y));
            Assert.IsTrue(Math.Round(v3.z) == Math.Round(v.z));

        }

    }
}
