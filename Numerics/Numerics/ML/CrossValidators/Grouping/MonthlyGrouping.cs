using System;
using System.Linq;

namespace CSharpNumerics.ML.CrossValidators.Grouping;

/// <summary>
/// Groups timestamps by month (year + month).
/// </summary>
public class MonthlyGrouping : ITimeGrouping
{
    /// <summary>
    /// Returns an integer group id per timestamp where all timestamps in the same month share the same id.
    /// </summary>
    /// <param name="time">Input timestamps.</param>
    /// <returns>Group ids aligned with <paramref name="time"/>.</returns>
    public int[] GetGroups(DateTime[] time)
    {
        return [.. time.Select(t => HashCode.Combine(t.Year, t.Month))];
    }
}
