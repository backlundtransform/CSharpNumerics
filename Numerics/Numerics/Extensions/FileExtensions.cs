using Numerics.Objects;
using System.Collections.Generic;
using System.Text;

namespace System.IO
{
    public static class FileExtensions
    {
        public static void Save(this IDictionary<Vector,Vector> data, string path) => data.Save(path, Encoding.Default);

        public static void Save(this IDictionary<Vector, Vector> data, string path, Encoding encoding)
        {
            var csv = new StringBuilder();

            var newLine = $"X,Y,Z,x,y,z";
            csv.AppendLine(newLine);

            foreach (var item in data)
            {
                newLine = $"{item.Value.x},{item.Value.y},{item.Value.z},{item.Key.x},{item.Key.y},{item.Key.z}";
                csv.AppendLine(newLine);

            }
            File.WriteAllText(path, csv.ToString(), encoding);

        }
    }
}
