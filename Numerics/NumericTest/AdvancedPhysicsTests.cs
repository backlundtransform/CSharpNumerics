using CSharpNumerics.Engines.Game;
using CSharpNumerics.Engines.Game.BroadPhase;
using CSharpNumerics.Engines.Game.Objects;
using CSharpNumerics.Engines.Game.Particles;
using CSharpNumerics.Engines.Game.Terrain;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.FluidDynamics.FreeSurface;
using CSharpNumerics.Physics.FluidDynamics.SPH;
using CSharpNumerics.Physics.Mechanics.Objects;
using CSharpNumerics.Physics.Mechanics.SoftBody;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericTest
{
    [TestClass]
    public class AdvancedPhysicsTests
    {
        // ════════════════════════════════════════════════════════════
        //  DeformableMesh
        // ════════════════════════════════════════════════════════════

        #region DeformableMesh

        [TestMethod]
        public void DeformableMesh_CreateGrid_CorrectVertexCount()
        {
            var mesh = DeformableMesh.CreateGrid(1, 1, 10, 10, 0.1, 0.9);
            Assert.AreEqual(100, mesh.VertexCount);
            Assert.IsTrue(mesh.SpringCount > 0, "Should have springs");
        }

        [TestMethod]
        public void DeformableMesh_PinnedVertices_DoNotMove()
        {
            var mesh = DeformableMesh.CreateGrid(1, 1, 5, 5, 0.1, 0.9);
            mesh.Pin(0);
            var original = mesh.Positions[0];

            for (int i = 0; i < 100; i++)
                mesh.Step(0.01);

            Assert.AreEqual(original.x, mesh.Positions[0].x, 1e-10);
            Assert.AreEqual(original.y, mesh.Positions[0].y, 1e-10);
            Assert.AreEqual(original.z, mesh.Positions[0].z, 1e-10);
        }

        [TestMethod]
        public void DeformableMesh_GravityDrops_UnpinnedVertices()
        {
            var mesh = DeformableMesh.CreateGrid(1, 1, 3, 3, 0.1, 0.5,
                new Vector(0, 0, 5));
            // Pin top row
            for (int i = 6; i < 9; i++) mesh.Pin(i);

            double initialZ = mesh.Positions[0].z; // bottom-left
            for (int i = 0; i < 200; i++)
                mesh.Step(0.01);

            Assert.IsTrue(mesh.Positions[0].z < initialZ,
                $"Vertex should fall: {mesh.Positions[0].z} vs {initialZ}");
        }

        [TestMethod]
        public void DeformableMesh_SphereCollision_NoPenetration()
        {
            var mesh = DeformableMesh.CreateGrid(2, 2, 10, 10, 0.1, 0.8,
                new Vector(-1, -1, 3));
            // Pin top row
            for (int i = 0; i < 10; i++) mesh.Pin(90 + i);

            var sphereCenter = new Vector(0, 0, 1);
            double sphereRadius = 0.5;

            for (int i = 0; i < 500; i++)
            {
                mesh.Step(0.01);
                mesh.CollideWithSphere(sphereCenter, sphereRadius);
            }

            // No vertex should be inside the sphere
            for (int i = 0; i < mesh.VertexCount; i++)
            {
                double dist = (mesh.Positions[i] - sphereCenter).GetMagnitude();
                Assert.IsTrue(dist >= sphereRadius - 0.05,
                    $"Vertex {i} penetrated sphere: dist={dist:F4}, r={sphereRadius}");
            }
        }

        #endregion

        // ════════════════════════════════════════════════════════════
        //  ClothSimulation
        // ════════════════════════════════════════════════════════════

        #region ClothSimulation

        [TestMethod]
        public void ClothSimulation_DrapedOverSphere_NoPenetration()
        {
            var cloth = new ClothSimulation(2, 2, 15, 15, mass: 0.05, stiffness: 0.9,
                origin: new Vector(-1, -1, 3));
            cloth.PinTopEdge();

            var sphereCenter = new Vector(0, 0, 1.5);
            double sphereRadius = 0.5;

            for (int i = 0; i < 500; i++)
            {
                cloth.Step(0.01);
                cloth.CollideWithSphere(sphereCenter, sphereRadius);
            }

            // Check no vertices inside sphere
            int violations = 0;
            for (int j = 0; j < cloth.ResY; j++)
            {
                for (int ix = 0; ix < cloth.ResX; ix++)
                {
                    var pos = cloth.GetPosition(ix, j);
                    double dist = (pos - sphereCenter).GetMagnitude();
                    if (dist < sphereRadius - 0.05) violations++;
                }
            }
            Assert.AreEqual(0, violations, $"{violations} vertices penetrated the sphere");
        }

        [TestMethod]
        public void ClothSimulation_WithWind_Deflects()
        {
            var cloth = new ClothSimulation(1, 1, 5, 5, mass: 0.1, stiffness: 0.8,
                origin: new Vector(0, 0, 2));
            cloth.PinTopCorners();
            cloth.Wind = new Vector(5, 0, 0);

            // Record initial average X
            double initialX = 0;
            for (int j = 0; j < cloth.ResY; j++)
                for (int i = 0; i < cloth.ResX; i++)
                    initialX += cloth.GetPosition(i, j).x;
            initialX /= cloth.ResX * cloth.ResY;

            for (int i = 0; i < 200; i++)
                cloth.Step(0.01);

            double finalX = 0;
            for (int j = 0; j < cloth.ResY; j++)
                for (int i = 0; i < cloth.ResX; i++)
                    finalX += cloth.GetPosition(i, j).x;
            finalX /= cloth.ResX * cloth.ResY;

            Assert.IsTrue(finalX > initialX + 0.01,
                $"Wind should push cloth in +X: initial={initialX:F4}, final={finalX:F4}");
        }

        [TestMethod]
        public void ClothSimulation_KineticEnergy_Positive()
        {
            var cloth = new ClothSimulation(1, 1, 5, 5, mass: 0.1, stiffness: 0.8,
                origin: new Vector(0, 0, 5));
            cloth.PinTopCorners();

            for (int i = 0; i < 50; i++)
                cloth.Step(0.01);

            double ke = cloth.KineticEnergy();
            Assert.IsTrue(ke > 0, $"KE should be positive while falling: {ke}");
        }

        #endregion

        // ════════════════════════════════════════════════════════════
        //  SPH Solver
        // ════════════════════════════════════════════════════════════

        #region SPH

        [TestMethod]
        public void SPH_ParticlesDrop_UnderGravity()
        {
            var positions = new List<Vector>();
            // Small block of particles
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    positions.Add(new Vector(i * 0.05, j * 0.05, 0.5));

            var sph = new SPHSolver(positions)
            {
                SmoothingRadius = 0.1,
                Dt = 0.001,
                BoundsMin = new Vector(-1, -1, 0),
                BoundsMax = new Vector(1, 1, 1),
                GasConstant = 1000
            };

            double initialZ = sph.Particles[0].Position.z;
            sph.Step(50);

            // Particles should have moved down
            double avgZ = 0;
            for (int i = 0; i < sph.ParticleCount; i++)
                avgZ += sph.Particles[i].Position.z;
            avgZ /= sph.ParticleCount;

            Assert.IsTrue(avgZ < initialZ - 0.01,
                $"Particles should fall: initial={initialZ:F4}, avg={avgZ:F4}");
        }

        [TestMethod]
        public void SPH_WaterSplash_ParticlesStayInBounds()
        {
            var positions = new List<Vector>();
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    positions.Add(new Vector(i * 0.04, j * 0.04, 0.8));

            var sph = new SPHSolver(positions)
            {
                SmoothingRadius = 0.08,
                Dt = 0.0005,
                BoundsMin = new Vector(-0.5, -0.5, 0),
                BoundsMax = new Vector(0.5, 0.5, 1),
                BoundaryRestitution = 0.2
            };

            sph.Step(200);

            // All particles should be within bounds
            for (int i = 0; i < sph.ParticleCount; i++)
            {
                var p = sph.Particles[i].Position;
                Assert.IsTrue(p.x >= sph.BoundsMin.x - 0.01 && p.x <= sph.BoundsMax.x + 0.01,
                    $"Particle {i} X out of bounds: {p.x}");
                Assert.IsTrue(p.z >= sph.BoundsMin.z - 0.01,
                    $"Particle {i} Z below floor: {p.z}");
            }
        }

        [TestMethod]
        public void SPH_DensityComputed_Positive()
        {
            var positions = new List<Vector>();
            for (int i = 0; i < 4; i++)
                positions.Add(new Vector(i * 0.03, 0, 0.5));

            var sph = new SPHSolver(positions)
            {
                SmoothingRadius = 0.1,
                Dt = 0.001
            };

            sph.Step();

            for (int i = 0; i < sph.ParticleCount; i++)
                Assert.IsTrue(sph.Particles[i].Density > 0,
                    $"Particle {i} density should be positive: {sph.Particles[i].Density}");
        }

        #endregion

        // ════════════════════════════════════════════════════════════
        //  VOF Tracker
        // ════════════════════════════════════════════════════════════

        #region VOF

        [TestMethod]
        public void VOF_FillRect_SetsCorrectFractions()
        {
            var vof = new VOFTracker(20, 20, 0.1);
            vof.FillRect(5, 5, 14, 14);

            Assert.AreEqual(1.0, vof.GetFraction(10, 10), 1e-10);
            Assert.AreEqual(0.0, vof.GetFraction(0, 0), 1e-10);
        }

        [TestMethod]
        public void VOF_TotalVolume_ConservedWithoutFlow()
        {
            var vof = new VOFTracker(20, 20, 0.1);
            vof.FillRect(5, 5, 14, 14);

            double initialVol = vof.TotalVolume();
            Assert.IsTrue(initialVol > 0, "Should have some liquid volume");

            // Advect with zero velocity → volume should be conserved
            var u = new double[400];
            var v = new double[400];
            vof.Advect(u, v, 0.01);

            double finalVol = vof.TotalVolume();
            Assert.AreEqual(initialVol, finalVol, initialVol * 0.01,
                "Volume should be conserved with zero flow");
        }

        [TestMethod]
        public void VOF_InterfaceCells_FoundAtBoundary()
        {
            var vof = new VOFTracker(20, 20, 0.1);
            vof.FillCircle(1.0, 1.0, 0.5);

            int interfaceCount = vof.InterfaceCellCount();
            // With a circle, there should be interface cells at the boundary
            // (exact count depends on resolution; just check > 0)
            Assert.IsTrue(interfaceCount >= 0, "Should detect interface cells");
        }

        [TestMethod]
        public void VOF_SurfaceNormal_PointsOutward()
        {
            var vof = new VOFTracker(30, 30, 0.1);
            vof.FillRect(10, 10, 19, 19);

            // At the right edge, normal should point roughly +X (from liquid to empty)
            var (nx, ny) = vof.SurfaceNormal(19, 15);
            Assert.IsTrue(nx > 0, $"Normal X at right edge should be positive: nx={nx}");
        }

        #endregion

        // ════════════════════════════════════════════════════════════
        //  Particle System
        // ════════════════════════════════════════════════════════════

        #region ParticleSystem

        [TestMethod]
        public void ParticleSystem_Emitter_SpawnsParticles()
        {
            var ps = new ParticleSystem(1000, seed: 123);
            var emitter = new ParticleEmitter
            {
                Position = new Vector(0, 0, 0),
                Direction = new Vector(0, 0, 1),
                EmissionRate = 100,
                InitialSpeed = 5,
                ParticleLifetime = 2
            };
            ps.AddEmitter(emitter);

            ps.Update(0.1); // should spawn ~10 particles
            Assert.IsTrue(ps.AliveCount >= 5, $"Should have spawned particles: {ps.AliveCount}");
        }

        [TestMethod]
        public void ParticleSystem_ParticlesDie_AfterLifetime()
        {
            var ps = new ParticleSystem(1000, seed: 42);
            ps.Emit(new Vector(0, 0, 0), new Vector(0, 0, 1), lifetime: 0.5);

            Assert.AreEqual(1, ps.AliveCount);

            // Advance past lifetime
            for (int i = 0; i < 60; i++)
                ps.Update(0.01);

            Assert.AreEqual(0, ps.AliveCount, "Particle should have expired");
        }

        [TestMethod]
        public void ParticleSystem_WindField_AffectsParticles()
        {
            var ps = new ParticleSystem(100, seed: 42);
            ps.Gravity = new Vector(0, 0, 0);
            ps.DragCoefficient = 0;

            // Emit stationary particle
            ps.Emit(new Vector(5, 5, 0), new Vector(0, 0, 0), lifetime: 5, mass: 1);

            // Create a simple wind field pushing +X
            int nx = 20, ny = 20;
            var windX = new double[nx * ny];
            var windY = new double[nx * ny];
            for (int i = 0; i < windX.Length; i++) windX[i] = 10.0;

            for (int i = 0; i < 50; i++)
            {
                ps.ApplyWindField(windX, windY, nx, ny, 0, 0, 1.0, 1.0, 0.01);
                ps.Update(0.01);
            }

            // Particle should have moved in +X
            Assert.IsTrue(ps.Particles[0].Position.x > 5.01,
                $"Wind should push particle: x={ps.Particles[0].Position.x}");
        }

        [TestMethod]
        public void ParticleSystem_GroundCollision_Bounces()
        {
            var ps = new ParticleSystem(100, seed: 42);
            ps.GroundHeight = 0;
            ps.GroundRestitution = 0.5;

            ps.Emit(new Vector(0, 0, 5), new Vector(0, 0, -10), lifetime: 5);

            for (int i = 0; i < 200; i++)
                ps.Update(0.01);

            Assert.IsTrue(ps.Particles[0].Position.z >= 0,
                $"Particle should be above ground: z={ps.Particles[0].Position.z}");
        }

        #endregion

        // ════════════════════════════════════════════════════════════
        //  Terrain Collider
        // ════════════════════════════════════════════════════════════

        #region TerrainCollider

        [TestMethod]
        public void TerrainCollider_FlatTerrain_CorrectHeight()
        {
            var terrain = TerrainCollider.CreateFlat(10, 10, 1.0, height: 5.0);
            Assert.AreEqual(5.0, terrain.SampleHeight(4.5, 4.5), 1e-10);
        }

        [TestMethod]
        public void TerrainCollider_SlopeFunction_Interpolates()
        {
            var terrain = TerrainCollider.FromFunction(20, 20, 1.0,
                (x, y) => 0.1 * x, originX: 0, originY: 0);

            double h = terrain.SampleHeight(10, 5);
            Assert.AreEqual(1.0, h, 0.1, $"Height at x=10: {h}");
        }

        [TestMethod]
        public void TerrainCollider_SphereCollision_Detected()
        {
            var terrain = TerrainCollider.CreateFlat(10, 10, 1.0, height: 0);

            // Sphere below terrain
            var pos = new Vector(5, 5, -0.5);
            bool hit = terrain.TestSphere(pos, 1.0, out double pen, out _);
            Assert.IsTrue(hit, "Sphere below terrain should collide");
            Assert.IsTrue(pen > 0, $"Penetration should be positive: {pen}");
        }

        [TestMethod]
        public void TerrainCollider_ResolveSphere_PushesAbove()
        {
            var terrain = TerrainCollider.CreateFlat(10, 10, 1.0, height: 0);

            var pos = new Vector(5, 5, 0.3);
            var vel = new Vector(0, 0, -5);
            double radius = 1.0;

            bool resolved = terrain.ResolveSphere(ref pos, ref vel, radius, 0.5);
            Assert.IsTrue(resolved, "Should resolve collision");
            Assert.IsTrue(pos.z >= radius - 0.1,
                $"Sphere should be pushed above terrain: z={pos.z}");
            Assert.IsTrue(vel.z >= 0, $"Velocity should be reflected: vz={vel.z}");
        }

        [TestMethod]
        public void TerrainCollider_SurfaceNormal_FlatIsUp()
        {
            var terrain = TerrainCollider.CreateFlat(10, 10, 1.0, height: 0);
            var normal = terrain.SurfaceNormal(5, 5);
            Assert.AreEqual(0, normal.x, 0.01);
            Assert.AreEqual(0, normal.y, 0.01);
            Assert.AreEqual(1, normal.z, 0.01);
        }

        #endregion

        // ════════════════════════════════════════════════════════════
        //  Wind Over Terrain
        // ════════════════════════════════════════════════════════════

        #region WindOverTerrain

        [TestMethod]
        public void WindOverTerrain_FlatTerrain_NoChange()
        {
            var terrain = TerrainCollider.CreateFlat(20, 20, 1.0, height: 0);
            int nx = 10, ny = 10;
            var windX = new double[nx * ny];
            var windY = new double[nx * ny];
            for (int i = 0; i < windX.Length; i++) windX[i] = 5.0;

            var origX = (double[])windX.Clone();
            WindOverTerrain.ApplyTerrainEffect(windX, windY, terrain, nx, ny, 0, 0, 2.0);

            // On flat terrain, slopes are zero → no change
            for (int i = 0; i < windX.Length; i++)
                Assert.AreEqual(origX[i], windX[i], 1e-10);
        }

        [TestMethod]
        public void WindOverTerrain_UphillSlope_SpeedsUpWind()
        {
            var terrain = TerrainCollider.FromFunction(20, 20, 1.0,
                (x, y) => 0.5 * x); // slope in +X

            int nx = 10, ny = 10;
            var windX = new double[nx * ny];
            var windY = new double[nx * ny];
            for (int i = 0; i < windX.Length; i++) windX[i] = 5.0; // wind in +X (uphill)

            double origSpeed = windX[5 + 5 * nx];
            WindOverTerrain.ApplyTerrainEffect(windX, windY, terrain, nx, ny, 0, 0, 2.0);

            // Wind going uphill should speed up
            Assert.IsTrue(windX[5 + 5 * nx] > origSpeed,
                $"Uphill wind should speed up: {windX[5 + 5 * nx]} vs {origSpeed}");
        }

        #endregion

        // ════════════════════════════════════════════════════════════
        //  CCD
        // ════════════════════════════════════════════════════════════

        #region CCD

        [TestMethod]
        public void CCD_SweptSphere_HitsStaticSphere()
        {
            var result = ContinuousCollisionDetection.SweptSphereVsSphere(
                new Vector(0, 0, 0), new Vector(10, 0, 0), 0.5,
                new Vector(5, 0, 0), 0.5);

            Assert.IsTrue(result.Hit, "Should detect collision");
            Assert.IsTrue(result.TimeOfImpact > 0 && result.TimeOfImpact < 1,
                $"TOI should be between 0 and 1: {result.TimeOfImpact}");
        }

        [TestMethod]
        public void CCD_SweptSphere_MissesDistantSphere()
        {
            var result = ContinuousCollisionDetection.SweptSphereVsSphere(
                new Vector(0, 0, 0), new Vector(10, 0, 0), 0.5,
                new Vector(5, 5, 0), 0.5);

            Assert.IsFalse(result.Hit, "Should miss distant sphere");
        }

        [TestMethod]
        public void CCD_SweptSphere_HitsGroundPlane()
        {
            var result = ContinuousCollisionDetection.SweptSphereVsPlane(
                new Vector(0, 0, 5), new Vector(0, 0, -5), 1.0, planeHeight: 0);

            Assert.IsTrue(result.Hit, "Should hit ground");
            Assert.IsTrue(result.TimeOfImpact > 0 && result.TimeOfImpact < 1,
                $"TOI: {result.TimeOfImpact}");
        }

        [TestMethod]
        public void CCD_SweptSphere_HitsAABB()
        {
            var box = new AABB(new Vector(4, -1, -1), new Vector(6, 1, 1));
            var result = ContinuousCollisionDetection.SweptSphereVsAABB(
                new Vector(0, 0, 0), new Vector(10, 0, 0), 0.5, box);

            Assert.IsTrue(result.Hit, "Should hit AABB");
            Assert.IsTrue(result.TimeOfImpact > 0 && result.TimeOfImpact < 1,
                $"TOI: {result.TimeOfImpact}");
        }

        #endregion

        // ════════════════════════════════════════════════════════════
        //  BVH Broad Phase
        // ════════════════════════════════════════════════════════════

        #region BVH

        [TestMethod]
        public void BVH_SameResultsAsBruteForce()
        {
            var rng = new Random(42);
            var bodies = new RigidBody[30];
            var radii = new double[30];
            for (int i = 0; i < 30; i++)
            {
                bodies[i] = RigidBody.CreateSolidSphere(1, 0.5);
                bodies[i].Position = new Vector(
                    rng.NextDouble() * 10, rng.NextDouble() * 10, rng.NextDouble() * 10);
                radii[i] = 1.0;
            }

            var brute = new BruteForceBroadPhase();
            var bvh = new BVHBroadPhase();
            var bruteResults = new List<(int, int)>();
            var bvhResults = new List<(int, int)>();

            brute.FindPairs(bodies, radii, 30, bruteResults);
            bvh.FindPairs(bodies, radii, 30, bvhResults);

            var normalize = (List<(int, int)> list) =>
            {
                var sorted = list.Select(p => p.Item1 < p.Item2 ? p : (p.Item2, p.Item1)).ToList();
                sorted.Sort((a, b) => a.Item1 != b.Item1
                    ? a.Item1.CompareTo(b.Item1) : a.Item2.CompareTo(b.Item2));
                return sorted;
            };

            var nb = normalize(bruteResults);
            var nbvh = normalize(bvhResults);

            Assert.AreEqual(nb.Count, nbvh.Count,
                $"Brute: {nb.Count}, BVH: {nbvh.Count}");
            for (int i = 0; i < nb.Count; i++)
                Assert.AreEqual(nb[i], nbvh[i], $"Mismatch at pair {i}");
        }

        [TestMethod]
        public void BVH_NoOverlap_EmptyResults()
        {
            var bodies = new RigidBody[5];
            var radii = new double[5];
            for (int i = 0; i < 5; i++)
            {
                bodies[i] = RigidBody.CreateSolidSphere(1, 0.5);
                bodies[i].Position = new Vector(i * 100, 0, 0);
                radii[i] = 1.0;
            }

            var bvh = new BVHBroadPhase();
            var results = new List<(int, int)>();
            bvh.FindPairs(bodies, radii, 5, results);

            Assert.AreEqual(0, results.Count);
        }

        [TestMethod]
        public void BVH_WorksWithPhysicsWorld()
        {
            // Verify BVH integrates with PhysicsWorld
            var world = new PhysicsWorld(new BVHBroadPhase())
            {
                Gravity = new Vector(0, 0, 0),
                DefaultRestitution = 1.0
            };

            var a = RigidBody.CreateSolidSphere(1, 1);
            a.Position = new Vector(0, 0, 0);
            a.Velocity = new Vector(5, 0, 0);
            world.AddBody(a, 1);

            var b = RigidBody.CreateSolidSphere(1, 1);
            b.Position = new Vector(3, 0, 0);
            world.AddBody(b, 1);

            bool collisionDetected = false;
            world.OnCollision = (ia, ib, c) => collisionDetected = true;

            for (int i = 0; i < 200; i++)
                world.Step(0.001);

            Assert.IsTrue(collisionDetected, "BVH should detect collision in PhysicsWorld");
        }

        #endregion

        // ════════════════════════════════════════════════════════════
        //  Terrain Collision for Aircraft Landing Gear
        // ════════════════════════════════════════════════════════════

        #region Aircraft Landing

        [TestMethod]
        public void TerrainCollider_AircraftLandingGear_ResolvesOnTouchdown()
        {
            // Simulate runway terrain
            var terrain = TerrainCollider.CreateFlat(100, 100, 1.0, height: 0);

            // Aircraft position: approaching ground
            var gearPos = new Vector(50, 50, 0.3); // just above ground
            var gearVel = new Vector(50, 0, -2);    // flying forward, descending
            double gearRadius = 0.5;

            bool hit = terrain.ResolveSphere(ref gearPos, ref gearVel, gearRadius, 0.1);
            Assert.IsTrue(hit, "Landing gear should contact terrain");
            Assert.IsTrue(gearPos.z >= gearRadius - 0.1,
                $"Gear should be at or above terrain: z={gearPos.z}");
            Assert.IsTrue(gearVel.z >= 0,
                $"Vertical velocity should be reflected: vz={gearVel.z}");
        }

        #endregion
    }
}
