using System;
using System.Collections.Generic;
using System.Linq;
using CSharpNumerics.Physics.Materials.Fire.Enums;

namespace CSharpNumerics.Physics.Materials.Fire;

/// <summary>
/// Static registry of known <see cref="FuelModel"/> instances.
/// Pre-registers all 13 Anderson standard fuel models.
/// Custom models can be added at runtime via <see cref="Register"/>.
/// </summary>
public static class FuelLibrary
{
    private static readonly Dictionary<FuelModelType, FuelModel> _registry =
        new Dictionary<FuelModelType, FuelModel>();

    static FuelLibrary()
    {
        Register(FuelModel.ShortGrass);
        Register(FuelModel.TimberGrassUnderstory);
        Register(FuelModel.TallGrass);
        Register(FuelModel.Chaparral);
        Register(FuelModel.Brush);
        Register(FuelModel.DormantBrush);
        Register(FuelModel.SouthernRough);
        Register(FuelModel.ClosedTimberLitter);
        Register(FuelModel.HardwoodLitter);
        Register(FuelModel.TimberLitterUnderstory);
        Register(FuelModel.LightLoggingSlash);
        Register(FuelModel.MediumLoggingSlash);
        Register(FuelModel.HeavyLoggingSlash);
    }

    /// <summary>
    /// Retrieves a fuel model by type. Throws <see cref="KeyNotFoundException"/> if not found.
    /// </summary>
    public static FuelModel Get(FuelModelType type)
    {
        if (_registry.TryGetValue(type, out var model))
            return model;
        throw new KeyNotFoundException($"Fuel model '{type}' not found in library.");
    }

    /// <summary>
    /// Tries to retrieve a fuel model by type. Returns true if found.
    /// </summary>
    public static bool TryGet(FuelModelType type, out FuelModel model) =>
        _registry.TryGetValue(type, out model);

    /// <summary>Returns all registered fuel models.</summary>
    public static IReadOnlyCollection<FuelModel> All() =>
        _registry.Values.ToList().AsReadOnly();

    /// <summary>
    /// Registers a custom fuel model. Overwrites any existing entry with the same type.
    /// </summary>
    public static void Register(FuelModel model)
    {
        _registry[model.Type] = model;
    }
}
