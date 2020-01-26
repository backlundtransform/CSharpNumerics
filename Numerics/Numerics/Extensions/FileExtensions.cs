using Numerics.Objects;
using System.Collections.Generic;
using System.Text;

namespace System.IO
{
    public static class FileExtensions
    {
        public static void Save(this IEnumerable<Vector> data, string path) => data.Save(path, Encoding.Default);

        public static void Save(this IEnumerable<Vector> data, string path, Encoding encoding)
        {
            var csv = new StringBuilder();

            var newLine = $"X,Y,Z";
            csv.AppendLine(newLine);

            foreach (var item in data)
            {
                newLine = $"{item.x},{item.y},{item.z}";
                csv.AppendLine(newLine);

            }
            File.WriteAllText(path, csv.ToString(), encoding);

        }
    }
}
