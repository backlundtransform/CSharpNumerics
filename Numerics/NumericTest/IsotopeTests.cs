using CSharpNumerics.Physics.Materials.Nuclear.Isotopes;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Linq;

namespace NumericTest
{
    [TestClass]
    public class IsotopeTests
    {
        [TestMethod]
        public void Cs137_HasCorrectProperties()
        {
            var cs = Isotope.Cs137;
            Assert.AreEqual("Cs137", cs.Name);
            Assert.AreEqual(55, cs.AtomicNumber);
            Assert.AreEqual(137, cs.MassNumber);
            Assert.AreEqual(DecayMode.BetaMinus, cs.DecayMode);
            Assert.IsTrue(cs.HalfLife > 9e8 && cs.HalfLife < 1e9, $"HalfLife: {cs.HalfLife}");
            Assert.IsFalse(cs.IsStable);
        }

        [TestMethod]
        public void Xe131_IsStable()
        {
            Assert.IsTrue(Isotope.Xe131.IsStable);
            Assert.AreEqual(0.0, Isotope.Xe131.Lambda);
        }

        [TestMethod]
        public void Lambda_MatchesHalfLife()
        {
            var cs = Isotope.Cs137;
            double expected = Math.Log(2) / cs.HalfLife;
            Assert.AreEqual(expected, cs.Lambda, 1e-20);
        }

        [TestMethod]
        public void Equality_SameIsotope()
        {
            var a = Isotope.Cs137;
            var b = Isotope.Cs137;
            Assert.AreEqual(a, b);
            Assert.IsTrue(a == b);
        }

        [TestMethod]
        public void Equality_DifferentIsotope()
        {
            Assert.AreNotEqual(Isotope.Cs137, Isotope.I131);
            Assert.IsTrue(Isotope.Cs137 != Isotope.I131);
        }

        [TestMethod]
        public void ToString_ContainsName()
        {
            string s = Isotope.Cs137.ToString();
            Assert.IsTrue(s.Contains("Cs137"), s);
        }

        [TestMethod]
        public void Library_GetByName()
        {
            var cs = IsotopeLibrary.Get("Cs137");
            Assert.AreEqual(Isotope.Cs137, cs);
        }

        [TestMethod]
        public void Library_CaseInsensitive()
        {
            var cs = IsotopeLibrary.Get("cs137");
            Assert.AreEqual(55, cs.AtomicNumber);
        }

        [TestMethod]
        [ExpectedException(typeof(System.Collections.Generic.KeyNotFoundException))]
        public void Library_UnknownThrows()
        {
            IsotopeLibrary.Get("Unobtainium99");
        }

        [TestMethod]
        public void Library_TryGet_ReturnsFalseForUnknown()
        {
            bool found = IsotopeLibrary.TryGet("Nope", out _);
            Assert.IsFalse(found);
        }

        [TestMethod]
        public void Library_AllContainsKnownIsotopes()
        {
            var all = IsotopeLibrary.All().ToList();
            Assert.IsTrue(all.Any(i => i.Name == "Cs137"));
            Assert.IsTrue(all.Any(i => i.Name == "I131"));
            Assert.IsTrue(all.Any(i => i.Name == "Sr90"));
            Assert.IsTrue(all.Count >= 8);
        }

        [TestMethod]
        public void Library_ByElement_FindsCaesium()
        {
            var caesiums = IsotopeLibrary.ByElement(55).ToList();
            Assert.IsTrue(caesiums.Any(i => i.Name == "Cs137"));
        }

        [TestMethod]
        public void Library_RegisterCustomIsotope()
        {
            var custom = new Isotope("TestIso99", 99, 200, 3600, 1e10, DecayMode.Alpha);
            IsotopeLibrary.Register(custom);
            var retrieved = IsotopeLibrary.Get("TestIso99");
            Assert.AreEqual(99, retrieved.AtomicNumber);
        }

        [TestMethod]
        public void SpecificActivity_Cs137_IsReasonable()
        {
            // Cs-137 specific activity ~3.2e12 Bq/kg
            Assert.IsTrue(Isotope.Cs137.SpecificActivity > 3e12);
            Assert.IsTrue(Isotope.Cs137.SpecificActivity < 4e12);
        }

        [TestMethod]
        public void GammaEnergy_Cs137()
        {
            Assert.AreEqual(0.662, Isotope.Cs137.GammaEnergy, 0.001);
        }

        [TestMethod]
        public void InhalationDoseCoeff_Cs137()
        {
            Assert.AreEqual(3.9e-8, Isotope.Cs137.InhalationDoseCoeff, 1e-10);
        }
    }
}
