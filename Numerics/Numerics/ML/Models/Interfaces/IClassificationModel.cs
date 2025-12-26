

namespace CSharpNumerics.ML.Models.Interfaces;

public interface IClassificationModel: IModel
{
    int NumClasses { get; }

}
