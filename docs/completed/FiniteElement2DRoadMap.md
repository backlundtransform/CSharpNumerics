# 2D Finite Element Analysis Roadmap

## Overview

Extend the existing 1D FEA framework (`FiniteElement/`) with **2D plane-stress/plane-strain elements** and a sparse linear algebra backend. The existing `PlaneStressSolver` in `Engines/Multiphysics/Solvers/` solves plane stress via iterative FD relaxation — this roadmap replaces its internals with a proper FEM assembly/solve pipeline, reusing the same `SimulationBuilder` API and `SimulationResult` output fields (`Ux`, `Uy`, `StressXX`, `StressYY`, `StressXY`, von Mises).

### What exists today

| Component | Location | Status |
|-----------|----------|--------|
| 1D FEA elements (`BarElement`, `BeamElement`) | `Numerics/FiniteElement/` | Done |
| `IElement1D`, `Mesh1D`, `Assembler1D` | `Numerics/FiniteElement/` | Done |
| `PlaneStressSolver` (FD-based, Gauss-Seidel SOR) | `Engines/Multiphysics/Solvers/` | Done — to be refactored |
| `SimulationResult` with `Ux`, `Uy`, `StressXX/YY/XY` | `Engines/Multiphysics/` | Done — no changes needed |
| `MultiphysicsType.PlaneStress` | `Engines/Multiphysics/Enums/` | Done — no changes needed |
| Dense `Matrix`, `VectorN` | `Numerics/Objects/` | Done |
| Sparse matrix / iterative solvers | — | **Missing** |
| 2D elements (tri, quad) | — | **Missing** |
| 2D mesh | — | **Missing** |
| 2D assembler | — | **Missing** |

### Architecture

```
Numerics/FiniteElement/           ← math primitives (section: Numerics)
  Interfaces/
    IElement1D.cs                 ← existing
    IElement2D.cs                 ← NEW: 2D element contract
  BarElement.cs                   ← existing
  BeamElement.cs                  ← existing
  TriElement.cs                   ← NEW: 3-node CST
  QuadElement.cs                  ← NEW: 4-node bilinear quad (Q4)
  Mesh1D.cs                      ← existing
  Mesh2D.cs                      ← NEW: structured triangular/quad mesh
  Assembler1D.cs                  ← existing
  Assembler2D.cs                  ← NEW: sparse assembly + CG solve

Numerics/Objects/
  SparseMatrix.cs                 ← NEW: CSR storage + SpMV

Engines/Multiphysics/Solvers/
  PlaneStressSolver.cs            ← MODIFY: delegate to Assembler2D
```

### Design decisions

1. **Sparse matrix in `Numerics/Objects/`** — it's a math primitive like `Matrix` and `VectorN`, usable across all sections.
2. **CG solver inside `Assembler2D`** — Conjugate Gradient with diagonal preconditioning, sufficient for symmetric positive-definite stiffness systems. No separate solver class needed.
3. **`PlaneStressSolver` becomes a thin adapter** — builds a `Mesh2D`, creates elements with material properties from `SimulationBuilder`, calls `Assembler2D.Solve()`, and maps results back to `SimulationResult` fields. Same API, better accuracy.
4. **Plane strain support** — same elements, different constitutive matrix (D). A `PlaneType` enum selects between plane stress and plane strain. The `SimulationBuilder` already carries `YoungsModulus` and `PoissonsRatio` from `EngineeringMaterial`.

### Estimated size

| Component | LOC |
|-----------|-----|
| `SparseMatrix` (CSR) | ~200 |
| `IElement2D` | ~30 |
| `TriElement` (CST) | ~120 |
| `QuadElement` (Q4) | ~180 |
| `Mesh2D` | ~150 |
| `Assembler2D` + CG solver | ~250 |
| `PlaneStressSolver` refactor | ~80 (net reduction) |
| Tests | ~250 |
| **Total** | **~1 260** |

---

## Phase 1 — Sparse Matrix Infrastructure

- [x] Create `SparseMatrix` in `Numerics/Objects/SparseMatrix.cs`
  - [x] CSR (Compressed Sparse Row) storage: `double[] Values`, `int[] ColIndices`, `int[] RowPointers`
  - [x] Build from triplet list: `SparseMatrix.FromTriplets(int rows, int cols, List<(int row, int col, double val)>)` — accumulates duplicates
  - [x] Sparse matrix–vector product: `Multiply(VectorN x) → VectorN`
  - [x] Diagonal extraction: `Diagonal() → VectorN` (for Jacobi preconditioning)
  - [x] Row/column elimination for Dirichlet BCs: `ApplyDirichlet(Dictionary<int, double> fixedDofs, VectorN rhs)`
  - [x] `Rows`, `Cols`, `NonZeroCount` properties
- [x] Unit tests for `SparseMatrix`
  - [x] Round-trip: dense → triplets → sparse → multiply matches dense multiply
  - [x] Dirichlet application preserves symmetry
  - [x] Empty-row handling

## Phase 2 — 2D Element Interface & CST Triangle

- [x] Create `IElement2D` in `Numerics/FiniteElement/Interfaces/IElement2D.cs`
  - [x] `int NodesPerElement { get; }`
  - [x] `int DofsPerNode => 2;` (ux, uy)
  - [x] `int TotalDofs { get; }`
  - [x] `Matrix LocalStiffness(double[,] nodeCoords, double thickness, double E, double nu, bool planeStress)`
  - [x] `VectorN LocalLoad(double[,] nodeCoords, double thickness, double qx, double qy)` — consistent body-force load vector
  - [x] `double[] ShapeFunctions(double xi, double eta)`
- [x] Create `TriElement` (CST) in `Numerics/FiniteElement/TriElement.cs`
  - [x] 3-node Constant Strain Triangle, 2 DOFs/node → 6×6 stiffness
  - [x] Area from cross product of edge vectors
  - [x] B-matrix (strain-displacement) from shape function derivatives: B = (1/2A)[∂N/∂x, ∂N/∂y]
  - [x] Constitutive matrix D: plane stress `E/(1-ν²) × [1,ν,0; ν,1,0; 0,0,(1-ν)/2]` or plane strain variant
  - [x] Ke = t · A · Bᵀ D B (closed-form, no numerical integration needed for CST)
  - [x] Consistent load vector: fe = t · A/3 · [qx, qy, qx, qy, qx, qy]
- [x] Unit tests for `TriElement`
  - [x] Stiffness symmetry
  - [x] Rigid body modes: Ke · [rigid translation] = 0
  - [x] Known analytical stiffness for unit right triangle

## Phase 3 — Q4 Quadrilateral Element

- [x] Create `QuadElement` (Q4) in `Numerics/FiniteElement/QuadElement.cs`
  - [x] 4-node bilinear quadrilateral, 2 DOFs/node → 8×8 stiffness
  - [x] Isoparametric mapping: shape functions N_i(ξ,η) = ¼(1±ξ)(1±η)
  - [x] Jacobian matrix J from ∂N/∂ξ, ∂N/∂η and physical node coordinates
  - [x] B-matrix via J⁻¹: ∂N/∂x = J⁻¹ · ∂N/∂ξ
  - [x] 2×2 Gauss quadrature for stiffness integration: Ke = Σ_gp (t · Bᵀ D B · det(J) · w_gp)
  - [x] Consistent load via same quadrature
- [x] Unit tests for `QuadElement`
  - [x] Stiffness symmetry
  - [x] Rigid body modes
  - [x] Rectangular element matches analytical (degenerate to simple formula)
  - [x] Patch test: uniform stress state recovered exactly

## Phase 4 — 2D Mesh & Assembler

- [x] Create `Mesh2D` in `Numerics/FiniteElement/Mesh2D.cs`
  - [x] Structured rectangular domain: `Mesh2D(double width, double height, int nx, int ny, ElementType type)` where `ElementType` is `Tri` or `Quad`
  - [x] Node coordinates: `double[,] Nodes` (nNodes × 2)
  - [x] Element connectivity: `int[,] Elements` (nElements × nodesPerElement)
  - [x] For Tri: each rectangle split into 2 triangles (diagonal split)
  - [x] For Quad: one element per rectangle
  - [x] `NodeCount`, `ElementCount` properties
  - [x] `GetElementNodes(int elemIndex) → double[,]` — extract node coordinates for an element
- [x] Create `Assembler2D` in `Numerics/FiniteElement/Assembler2D.cs`
  - [x] Constructor: `Assembler2D(Mesh2D mesh, IElement2D elementTemplate, double E, double nu, double thickness, bool planeStress)`
  - [x] `Assemble()` — build global stiffness in triplet form, convert to `SparseMatrix`
  - [x] `ApplyNodalLoad(int nodeIndex, int direction, double value)` — direction: 0=x, 1=y
  - [x] `ApplyBodyForce(double qx, double qy)` — distribute to all elements
  - [x] `Solve(Dictionary<int, double> fixedDofs) → VectorN` — apply Dirichlet BCs, solve via Preconditioned Conjugate Gradient
  - [x] Internal PCG: diagonal (Jacobi) preconditioning, tolerance 1e-10, max 10 000 iterations
- [x] Unit tests for `Mesh2D` and `Assembler2D`
  - [x] Mesh node count: (nx+1)×(ny+1) nodes, correct element count
  - [x] Global stiffness dimensions: 2×nodeCount
  - [x] Patch test (uniform strain under linear displacement BCs)

## Phase 5 — PlaneStressSolver Refactor

- [x] Refactor `PlaneStressSolver` to use FEM pipeline
  - [x] Build `Mesh2D` from `cfg.Nx`, `cfg.Ny`, `cfg.GeomWidth`, `cfg.GeomHeight`
  - [x] Create `Assembler2D` with material from `cfg.Material` (E, ν), thickness = 1.0 (plane stress default)
  - [x] Map `cfg.Sources2D` point loads to `ApplyNodalLoad()` — find nearest node to (ix, iy) grid coordinate
  - [x] Map `cfg.UniformLoad` to `ApplyBodyForce(0, uniformLoad)`
  - [x] Map boundary conditions: left fixed (all nodes at x=0 → ux=uy=0), bottom roller (uy=0), matching current FD solver
  - [x] Solve → extract displacement, compute stresses per element, project to nodes
  - [x] Populate `SimulationResult` fields: `Ux`, `Uy`, `StressXX`, `StressYY`, `StressXY`, von Mises `Field`, `MaxValue`, `MinValue`
  - [x] Remove old FD-based Gauss-Seidel iteration loop
- [x] Verify backward compatibility: existing `PlaneStress` tests still pass
- [x] Add `PlaneType` option to `SimulationBuilder` (`.WithPlaneType(PlaneType.PlaneStress | PlaneStrain)`) — default: PlaneStress

## Phase 6 — Tests & Documentation

- [x] Integration tests in `NumericTest/FiniteElement2DTests.cs` (extend existing file)
  - [x] Cantilever beam (2D FEM) vs analytical tip deflection: w = PL³/(3EI) — verify convergence with mesh refinement
  - [x] Plate with hole: stress concentration factor K ≈ 3.0 at hole edge (qualitative)
  - [x] Cook's membrane (skewed quad mesh): standard FEM benchmark for Q4 elements
  - [x] Plane strain thick cylinder under internal pressure: σ_rr, σ_θθ vs Lamé solution
- [x] PlaneStressSolver regression tests
  - [x] Existing cantilever-like setup produces same-order results
  - [x] von Mises stress field matches expectations
- [x] Update `Numerics/README.md` with 2D FEA section
  - [x] Element catalogue (Bar, Beam, CST, Q4)
  - [x] Usage examples for 2D FEM
  - [x] Sparse matrix API summary

---

## Key Files

### New files
| File | Section | LOC |
|------|---------|-----|
| `Numerics/Objects/SparseMatrix.cs` | Numerics | ~200 |
| `Numerics/FiniteElement/Interfaces/IElement2D.cs` | Numerics | ~30 |
| `Numerics/FiniteElement/TriElement.cs` | Numerics | ~120 |
| `Numerics/FiniteElement/QuadElement.cs` | Numerics | ~180 |
| `Numerics/FiniteElement/Mesh2D.cs` | Numerics | ~150 |
| `Numerics/FiniteElement/Assembler2D.cs` | Numerics | ~250 |
| `Numerics/FiniteElement/Enums/ElementType.cs` | Numerics | ~10 |
| `Numerics/FiniteElement/Enums/PlaneType.cs` | Numerics | ~10 |

### Existing — to modify
| File | Change |
|------|--------|
| `Engines/Multiphysics/Solvers/PlaneStressSolver.cs` | Replace FD internals with `Assembler2D` pipeline |
| `Engines/Multiphysics/SimulationBuilder.cs` | Add `.WithPlaneType()`, `.WithThickness()` (optional) |
| `NumericTest/FiniteElementTests.cs` | Add 2D FEM tests |
| `Numerics/README.md` | Add 2D FEA documentation |

### Existing — to reuse (no changes)
| File | Usage |
|------|-------|
| `Numerics/Objects/Matrix.cs` | Local element stiffness matrices (small dense) |
| `Numerics/Objects/VectorN.cs` | Load vectors, solution vectors |
| `Numerics/FiniteElement/Assembler1D.cs` | Design reference for assembly pattern |
| `Engines/Multiphysics/SimulationResult.cs` | Already has `Ux`, `Uy`, `StressXX/YY/XY` — no changes needed |
| `Engines/Multiphysics/Enums/MultiphysicsType.cs` | `PlaneStress` already exists |

---

## Usage Example

### Direct FEM API (Numerics layer)

```csharp
using CSharpNumerics.Numerics.FiniteElement;
using CSharpNumerics.Numerics.FiniteElement.Enums;

// Rectangular plate: 2m × 0.5m, 40×10 quad elements
var mesh = new Mesh2D(width: 2.0, height: 0.5, nx: 40, ny: 10, ElementType.Quad);

// Assemble: steel, plane stress, 1mm thickness
var asm = new Assembler2D(mesh, new QuadElement(), E: 200e9, nu: 0.3, thickness: 0.001, planeStress: true);
asm.Assemble();

// Tip load: 1 kN downward at top-right corner node
int tipNode = mesh.NodeCount - 1;
asm.ApplyNodalLoad(tipNode, direction: 1, value: -1000.0);

// Fixed left edge
var bcs = new Dictionary<int, double>();
for (int i = 0; i <= 10; i++)  // nodes along x=0
{
    bcs[i * 2] = 0.0;     // ux = 0
    bcs[i * 2 + 1] = 0.0; // uy = 0
}

VectorN u = asm.Solve(bcs);
```

### Via Multiphysics SimulationBuilder (Engine layer)

```csharp
using CSharpNumerics.Engines.Multiphysics;
using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Physics.Materials.Engineering;

var result = SimulationType.Create(MultiphysicsType.PlaneStress)
    .WithMaterial(EngineeringLibrary.Steel)
    .WithGeometry(width: 2.0, height: 0.5, nx: 40, ny: 10)
  .AddSource(ix: 40, iy: 10, value: -1000.0)
    .Solve();

double[,] vonMises = result.Field;       // von Mises stress field
double[,] ux = result.Ux;                // x-displacement
double[,] uy = result.Uy;                // y-displacement
double maxStress = result.MaxValue;       // peak von Mises (Pa)
```
