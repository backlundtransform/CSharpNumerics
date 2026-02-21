using CSharpNumerics.Numerics.Objects;
using System.Collections.Generic;
using System.Globalization;
using System.Text;

namespace System.IO
{
    public static class FileExtensions
    {
        public static void Save(this IDictionary<Vector,Vector> data, string path) => data.Save(path, Encoding.Default);
        public static void Save<T>(this IEnumerable<T> data, string path) => data.Save(path, Encoding.Default);

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


        public static void Save(this IDictionary<double, Vector> data, string path, Encoding encoding)
        {
            var csv = new StringBuilder();

            var newLine = $"X,Y,Z,t";
            csv.AppendLine(newLine);

            foreach (var item in data)
            {
                newLine = $"{item.Value.x},{item.Value.y},{item.Value.z},{item.Key}";
                csv.AppendLine(newLine);

            }
            File.WriteAllText(path, csv.ToString(), encoding);

        }

        public static void Save<T>(this IEnumerable<T> data, string path, Encoding encoding)
        {
            var csv = new StringBuilder();

    
            foreach (var item in data)
            {

                var newLine = "";
                var properties = typeof(T).GetProperties();
                for (var i = 0; i < properties.Length; i++)
                {
                    newLine += Convert.ToString(properties[i].GetValue(item),CultureInfo.InvariantCulture);
                    if (i != properties.Length - 1)
                    {
                        newLine += ",";
                    }
                }
           
                csv.AppendLine(newLine );

            }
            File.WriteAllText(path, csv.ToString(), encoding);

        }

    }
}
