using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.ML.CrossValidators.Grouping;

public interface ITimeGrouping
{
    int[] GetGroups(DateTime[] time);

}
