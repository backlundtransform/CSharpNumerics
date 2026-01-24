using System;
using System.Globalization;
using System.Linq;

namespace CSharpNumerics.ML.CrossValidators.Grouping;

/// <summary>
/// Groups timestamps by ISO-8601 week.
/// ISO weeks start on Monday and week 1 is the week with the year's first Thursday.
/// </summary>
public class WeeklyGrouping : ITimeGrouping
{
    /// <summary>
    /// Returns an integer group id per timestamp where all timestamps in the same ISO week share the same id.
    /// </summary>
    /// <param name="time">Input timestamps.</param>
    /// <returns>Group ids aligned with <paramref name="time"/>.</returns>
    public int[] GetGroups(DateTime[] time)
    {
        return [.. time.Select(GetIsoWeekKey)];
    }

    private static int GetIsoWeekKey(DateTime t)
    {
        var date = t.Date;
        var week = ISOWeek.GetWeekOfYear(date);
        var year = ISOWeek.GetYear(date);

        // Stable and collision-resistant integer key for (isoYear, isoWeek)
        return HashCode.Combine(year, week);
    }
}
