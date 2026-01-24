using System;

using System.Linq;


namespace CSharpNumerics.ML.CrossValidators.Grouping;

public class DailyGrouping : ITimeGrouping
{
    public int[] GetGroups(DateTime[] time)
    {
        return [.. time.Select(t => t.Date.GetHashCode())];
    }
}
