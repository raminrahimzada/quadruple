using System;

namespace Quadruple.Example
{
    class Program
    {
        static void Main(string[] args)
        { 
            Console.WriteLine(Quad.Sqrt(2).ToString(QuadrupleStringFormat.ScientificExact));
            Console.WriteLine(Math.Sqrt(2));
            //returns me
            //1.4142135623730951272e+0
            //1.4142135623731
            Console.WriteLine();
            Console.WriteLine(Quad.Cos(2).ToString(QuadrupleStringFormat.ScientificExact));
            Console.WriteLine(Math.Cos(2));
            //returns me
            //   -4.161468365471424375175e-1
            //  -0.416146836547142

            Console.Read();
        }
    }
}
