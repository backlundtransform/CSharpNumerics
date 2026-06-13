using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteElement;
using CSharpNumerics.Numerics.FiniteElement.Enums;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves 2D linear elasticity using a finite-element pipeline on a structured triangular mesh.
/// Public output stays aligned with the legacy Nx×Ny field layout for compatibility.
/// </summary>
internal class PlaneStressSolver : IMultiphysicsSolver
{
    public SimulationResult Solve(SimulationBuilder cfg)
    {
        var mat = cfg.Material.Value;
        double E = mat.YoungsModulus;
        double nu = mat.PoissonsRatio;

        int nx = cfg.Nx;
        int ny = cfg.Ny;
        bool planeStress = cfg.PlaneType == PlaneType.PlaneStress;

        var mesh = new Mesh2D(cfg.GeomWidth, cfg.GeomHeight, nx, ny, ElementType.Tri);
        var assembler = new Assembler2D(mesh, new TriElement(), E, nu, cfg.Thickness, planeStress);
        assembler.Assemble();

        foreach (var (ix, iy, value) in cfg.Sources2D)
        {
            int nodeX = Math.Clamp(ix, 0, nx);
            int nodeY = Math.Clamp(iy, 0, ny);
            int nodeIndex = mesh.GetNodeIndex(nodeX, nodeY);
            assembler.ApplyNodalLoad(nodeIndex, 0, value);
        }

        if (cfg.UniformLoad != 0.0)
            assembler.ApplyBodyForce(0.0, cfg.UniformLoad);

        var fixedDofs = BuildBoundaryConditions(mesh);
        var displacement = assembler.Solve(fixedDofs);
        assembler.ComputeElementStresses(displacement, out var elementSxx, out var elementSyy, out var elementSxy);

        BuildFieldOutputs(mesh, displacement, elementSxx, elementSyy, elementSxy,
            out var ux, out var uy, out var sxx, out var syy, out var sxy, out var vonMises,
            out double minV, out double maxV);

        return new SimulationResult
        {
            Type = MultiphysicsType.PlaneStress,
            Field = vonMises,
            Ux = ux,
            Uy = uy,
            StressXX = sxx,
            StressYY = syy,
            StressXY = sxy,
            MaxValue = maxV,
            MinValue = minV,
            Iterations = 1,
            FinalTime = 0.0
        };
    }

    private static Dictionary<int, double> BuildBoundaryConditions(Mesh2D mesh)
    {
        var fixedDofs = new Dictionary<int, double>();

        for (int iy = 0; iy <= mesh.Ny; iy++)
        {
            int node = mesh.GetNodeIndex(0, iy);
            fixedDofs[node * 2] = 0.0;
            fixedDofs[node * 2 + 1] = 0.0;
        }

        for (int ix = 0; ix <= mesh.Nx; ix++)
        {
            int node = mesh.GetNodeIndex(ix, 0);
            fixedDofs[node * 2 + 1] = 0.0;
        }

        return fixedDofs;
    }

    private static void BuildFieldOutputs(
        Mesh2D mesh,
        VectorN displacement,
        double[] elementSxx,
        double[] elementSyy,
        double[] elementSxy,
        out double[,] ux,
        out double[,] uy,
        out double[,] sxx,
        out double[,] syy,
        out double[,] sxy,
        out double[,] vonMises,
        out double minV,
        out double maxV)
    {
        int nx = mesh.Nx;
        int ny = mesh.Ny;
        int nodesPerRow = nx + 1;

        ux = new double[nx, ny];
        uy = new double[nx, ny];
        sxx = new double[nx, ny];
        syy = new double[nx, ny];
        sxy = new double[nx, ny];
        vonMises = new double[nx, ny];

        var nodeSxx = new double[nodesPerRow * (ny + 1)];
        var nodeSyy = new double[nodesPerRow * (ny + 1)];
        var nodeSxy = new double[nodesPerRow * (ny + 1)];
        var nodeCounts = new int[nodesPerRow * (ny + 1)];

        for (int e = 0; e < mesh.ElementCount; e++)
        {
            for (int local = 0; local < mesh.NodesPerElement; local++)
            {
                int node = mesh.Elements[e, local];
                nodeSxx[node] += elementSxx[e];
                nodeSyy[node] += elementSyy[e];
                nodeSxy[node] += elementSxy[e];
                nodeCounts[node]++;
            }
        }

        for (int node = 0; node < nodeCounts.Length; node++)
        {
            if (nodeCounts[node] == 0)
                continue;

            nodeSxx[node] /= nodeCounts[node];
            nodeSyy[node] /= nodeCounts[node];
            nodeSxy[node] /= nodeCounts[node];
        }

        minV = double.MaxValue;
        maxV = double.MinValue;

        for (int iy = 0; iy < ny; iy++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                int node = mesh.GetNodeIndex(ix, iy);
                ux[ix, iy] = displacement[node * 2];
                uy[ix, iy] = displacement[node * 2 + 1];
                sxx[ix, iy] = nodeSxx[node];
                syy[ix, iy] = nodeSyy[node];
                sxy[ix, iy] = nodeSxy[node];

                double vm = Math.Sqrt(sxx[ix, iy] * sxx[ix, iy]
                                    - sxx[ix, iy] * syy[ix, iy]
                                    + syy[ix, iy] * syy[ix, iy]
                                    + 3.0 * sxy[ix, iy] * sxy[ix, iy]);
                vonMises[ix, iy] = vm;
                if (vm < minV) minV = vm;
                if (vm > maxV) maxV = vm;
            }
        }
    }
}
