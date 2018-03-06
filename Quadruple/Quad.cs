/*
    Copyright (c) 2011 Jeff Pasternack.  All rights reserved.
 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

using System;
using System.Text;

namespace Quadruple
{
    /// <inheritdoc>
    ///     <cref></cref>
    /// </inheritdoc>
    /// <summary>
    /// Quad is a signed 128-bit floating point number, stored internally as a 64-bit significand (with the most significant bit as the sign bit) and
    /// a 64-bit signed exponent, with a value == significand * 2^exponent.  Quads have both a higher precision (64 vs. 53 effective significand bits)
    /// and a much higher range (64 vs. 11 exponent bits) than doubles, but also support NaN and PositiveInfinity/NegativeInfinity values and can be generally
    /// used as a drop-in replacement for doubles, much like double is a drop-in replacement for float.  Operations are checked and become +/- infinity in the
    /// event of overflow (values larger than ~8E+2776511644261678592) and 0 in the event of underflow (values less than ~4E-2776511644261678592).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Exponents &gt;= long.MaxValue - 64 and exponents &lt;= long.MinValue + 64 are reserved
    /// and constitute overflow and underflow, respectively.  Zero, PositiveInfinity, NegativeInfinity and NaN are
    /// defined by significand bits == 0 and an exponent of long.MinValue + 0, + 1, + 2, and + 3, respectively.
    /// </para>
    /// <para>
    /// Quad multiplication and division operators are slightly imprecise for the sake of efficiency; specifically,
    /// they may assign the wrong least significant bit, such that the precision is effectively only 63 bits rather than 64.
    /// </para>    
    /// <para>
    /// For speed, consider using instance methods (like Multiply and Add) rather
    /// than the operators (like * and +) when possible, as the former are significantly faster (by as much as 50%).
    /// </para>    
    /// </remarks>
    [System.Diagnostics.DebuggerDisplay("{" + nameof(ToString) + "(),nq}")] //this attributes makes the debugger display the value without braces or quotes
    public struct Quad : IComparable<Quad>
    {
        #region public math constants

        static Quad()
        {
            PI = Parse("3.14159265358979323846264338327950288419716939937510");
            PIx2 = Parse("6.28318530717958647692528676655900576839433879875021");
            E = Parse("2.7182818284590452353602874713526624977572470936999595749");
            PIdiv2 = Parse("1.570796326794896619231321691639751442098584699687552910487");
            PIdiv4 = Parse("0.785398163397448309615660845819875721049292349843776455243");
            Einv = Parse("0.3678794411714423215955237701614608674458111310317678");
            LOG2 = Parse("0.693147180559945309417232121458176568075500134360255254120");
            Log10Inv = Parse("0.434294481903251827651128918916605082294397005803666566114");
            OneDiv2 = Parse("0.5");
        }

        /// <summary>
        /// represents PI
        /// </summary>
        public static readonly Quad PI;
        /// <summary>
        /// represents 2*PI
        /// </summary>
        public static readonly Quad PIx2;

        /// <summary>
        /// represents E
        /// </summary>
        public static readonly Quad E; 
        /// <summary>
        /// represents PI/2
        /// </summary>
        public static readonly Quad PIdiv2;
        /// <summary>
        /// represents PI/4
        /// </summary>
        public static readonly Quad PIdiv4  ;

        /// <summary>
        /// represents 1.0/E
        /// </summary>
        public static readonly Quad Einv;
        /// <summary>
        /// represents Logarithm 2 from base E
        /// </summary>
        public static readonly Quad LOG2 ;
        /// <summary>
        /// log(10,E) factor
        /// </summary>
        public static readonly Quad Log10Inv  ;

        /// <summary>
        /// just 0.5
        /// </summary>
        public static readonly Quad OneDiv2  ;
        /// <summary>
        /// Max iterations count in Taylor series
        /// </summary>
        public static int MaxIteration = 100;
        #endregion

        #region new math functions

        /// <summary>
        /// Analogy of Math.Exp method
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Quad Exp(Quad x)
        {
            int count = 0;
            while (x > One)
            {
                x--;
                count++;
            }
            while (x < Zero)
            {
                x++;
                count--;
            }
            int iteration = 1;
            Quad result = One;
            Quad fatorial = One;
            Quad cachedResult;
            do
            {
                cachedResult = result;
                fatorial *= x / iteration++;
                result += fatorial;
            } while (cachedResult != result);
            if (count != 0) result = result * PowerN(E, count);
            return result;
        }
        /// <summary>
        /// Analogy of Math.Pow method
        /// </summary>
        /// <param name="value"></param>
        /// <param name="pow"></param>
        /// <returns></returns>
        public static Quad Power(Quad value, Quad pow)
        {
            return Exp(pow * Log(value));
        }
        /// <summary>
        /// Power to the integer value
        /// </summary>
        /// <param name="value"></param>
        /// <param name="power"></param>
        /// <returns></returns>
        public static Quad PowerN(Quad value, int power)
        {
            if (power == 0) return One;
            if (power < 0) return PowerN(One / value, -power);

            var q = power;
            var prod = One;
            var current = value;
            while (q > 0)
            {
                if (q % 2 == 1)
                {
                    // detects the 1s in the binary expression of power
                    prod = current * prod; // picks up the relevant power
                    q--;
                }
                current = current * current; // value^i -> value^(2*i)
                q = q / 2;
            }

            return prod;
        }

        /// <summary>
        /// Analogy of Math.Log10
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Quad Log10(Quad x)
        {
            return Log(x) * Log10Inv;
        }
        /// <summary>
        /// Analogy of Math.Log
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Quad Log(Quad x)
        {
            if (x <= Zero)
            {
                throw new ArgumentException("x must be greater than zero");
            }
            int count = 0;
            while (x >= 1)
            {
                x *= Einv;
                count++;
            }
            while (x <= Einv)
            {
                x *= E;
                count--;
            }
            x--;
            if (x == 0) return count;
            Quad result = Zero;
            int iteration = 0;
            Quad y =One;
            Quad cacheResult = result - One;
            while (cacheResult != result && iteration < MaxIteration)
            {
                iteration++;
                cacheResult = result;
                y *= -x;
                result += y / iteration;
            }
            return count - result;
        }
        /// <summary>
        /// Analogy of Math.Cos
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Quad Cos(Quad x)
        {
            while (x > PIx2)
            {
                x -= PIx2;
            }
            while (x < -PIx2)
            {
                x += PIx2;
            }
            // now x in (-2pi,2pi)
            if (x >= PI && x <= PIx2)
            {
                return -Cos(x - PI);
            }
            if (x >= -PIx2 && x <= -PI)
            {
                return -Cos(x + PI);
            }
            x = x * x;
            //y=1-x/2!+x^2/4!-x^3/6!...
            Quad xx = -x * OneDiv2;
            Quad y = One + xx;
            Quad cachedY = y - One;//init cache  with different value
            for (int i = 1; cachedY != y && i < MaxIteration; i++)
            {
                cachedY = y;
                Quad factor = i * (i + i + 3) + 1; //2i^2+2i+i+1=2i^2+3i+1
                factor = -OneDiv2 / factor;
                xx *= x * factor;
                y += xx;
            }
            return y;
        }
        /// <summary>
        /// Analogy of Math.Tan
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Quad Tan(Quad x)
        {
            return Sin(x) / Cos(x);
        }
        /// <summary>
        /// Analogy of Math.Sin
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Quad Sin(Quad x)
        {
            var cos = Cos(x);
            var real = Math.Sin((double)x);
            return Sqrt(One - cos * cos) * Math.Sign(real);
        }

        public static Quad Sqrt(Quad x )
        {
            return Sqrt(x, Epsilon);
        }
        /// <summary>
        /// Analogy of Math.Sqrt
        /// </summary>
        /// <param name="x"></param>
        /// <param name="epsilon">lasts iteration while error less than this epsilon</param>
        /// <returns></returns>
        public static Quad Sqrt(Quad x, Quad epsilon)
        {
            if (x < 0) throw new OverflowException("Cannot calculate square root from a negative number");
            //initial approximation
            Quad current = Math.Sqrt((double)x), previous;
            do
            {
                previous = current;
                if (previous == Zero) return Zero;
                current = (previous + x / previous) * OneDiv2;
            } while (Abs(previous - current) > epsilon);
            return current;
        }
        /// <summary>
        /// Analogy of Math.Sinh
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Quad Sinh(Quad x)
        {
            var y = Exp(x);
            var yy = One / y;
            return (y - yy) * OneDiv2;
        }
        /// <summary>
        /// Analogy of Math.Cosh
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Quad Cosh(Quad x)
        {
            var y = Exp(x);
            var yy = One / y;
            return (y + yy) * OneDiv2;
        }
        /// <summary>
        /// Analogy of Math.Tanh
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Quad Tanh(Quad x)
        {
            var y = Exp(x);
            var yy = One / y;
            return (y - yy) / (y + yy);
        }

        /// <summary>
        /// Analogy of Math.Asin
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Quad Asin(Quad x)
        {
            if (x > One || x < -One)
            {
                throw new ArgumentException("x must be in [-1,1]");
            }
            //known values
            if (x == Zero) return 0;
            if (x == One) return PIdiv2;
            //asin function is odd function
            if (x < Zero) return -Asin(-x);

            //my optimize trick here

            // used a math formula to speed up :
            // asin(x)=0.5*(pi/2-asin(1-2*x*x)) 
            // if x>=0 is true

            var newX = 1 - 2 * x * x;

            //for calculating new value near to zero than current
            //because we gain more speed with values near to zero
            if (Abs(x) > Abs(newX))
            {
                var t = Asin(newX);
                return OneDiv2 * (PIdiv2 - t);
            }
            Quad y = Zero;
            Quad result = x;
            Quad cachedResult;
            int i = 1;
            y += result;
            var xx = x * x;
            do
            {
                cachedResult = result;
                result *= xx * (One - One / (i + i));
                y += result / (2 * i + 1);
                i++;
            } while (cachedResult != result);
            return y;
        }
        /// <summary>
        /// Analogy of Math.Atan
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Quad ATan(Quad x)
        {
            if (x == Zero) return Zero;
            if (x == One) return PIdiv4;
            return Asin(x / Sqrt(1 + x * x));
        }
        /// <summary>
        /// Analogy of Math.Acos
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Quad Acos(Quad x)
        {
            if (x == Zero) return PIdiv2;
            if (x == One) return Zero;
            if (x < Zero) return PI - Acos(-x);
            return PIdiv2 - Asin(x);
        }

        /// <summary>
        /// Analogy of Math.Atan2
        /// for more see this
        /// <seealso cref="http://i.imgur.com/TRLjs8R.png"/>
        /// </summary>
        /// <param name="y"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Quad Atan2(Quad y, Quad x)
        {
            if (x > Zero)
            {
                return ATan(y / x);
            }
            if (x < Zero && y >= Zero)
            {
                return ATan(y / x) + PI;
            }
            if (x < Zero && y < Zero)
            {
                return ATan(y / x) - PI;
            }
            if (x == Zero && y > Zero)
            {
                return PIdiv2;
            }
            if (x == Zero && y < Zero)
            {
                return -PIdiv2;
            }
            throw new ArgumentException("invalid atan2 arguments");
        }
        #endregion

        #region Public constants
        /// <summary>
        /// 0.  Equivalent to (Quad)0.
        /// </summary>
        public static readonly Quad Zero = new Quad(0UL, long.MinValue); //there is only one zero; all other significands with exponent long.MinValue are invalid.

        /// <summary>
        /// 1.  Equivalent to (Quad)1.
        /// </summary>
        public static readonly Quad One = (Quad)1UL; //used for increment/decrement operators

        /// <summary>
        /// Positive infinity.  Equivalent to (Quad)double.PositiveInfinity.
        /// </summary>
        public static readonly Quad PositiveInfinity = new Quad(0UL, InfinityExponent);

        /// <summary>
        /// Negative infinity.  Equivalent to (Quad)double.NegativeInfinity.
        /// </summary>
        public static readonly Quad NegativeInfinity = new Quad(0UL, NegativeInfinityExponent);

        /// <summary>
        /// The Not-A-Number value.  Equivalent to (Quad)double.NaN.
        /// </summary>
        public static readonly Quad NaN = new Quad(0UL, NotANumberExponent);

        /// <summary>
        /// The maximum value representable by a Quad, (2 - 1/(2^63)) * 2^(long.MaxValue-65)
        /// </summary>
        public static readonly Quad MaxValue = new Quad(~HighestBit, ExponentUpperBound);

        /// <summary>
        /// The minimum value representable by a Quad, -(2 - 1/(2^63)) * 2^(long.MaxValue-65)
        /// </summary>
        public static readonly Quad MinValue = new Quad(ulong.MaxValue, ExponentUpperBound);

        /// <summary>
        /// The smallest positive value greater than zero representable by a Quad, 2^(long.MinValue+65)
        /// </summary>
        public static readonly Quad Epsilon = new Quad(0UL, ExponentLowerBound);
        #endregion

        #region Public fields
        /// <summary>
        /// The first (most significant) bit of the significand is the sign bit; 0 for positive values, 1 for negative.
        /// The remainder of the bits represent the fractional part (after the binary point) of the significant; there is always an implicit "1"
        /// preceding the binary point, just as in IEEE's double specification.  For "special" values 0, PositiveInfinity, NegativeInfinity, and NaN,
        /// SignificantBits == 0.
        /// </summary>
        public ulong SignificandBits;

        /// <summary>
        /// The value of the Quad == (-1)^[first bit of significant] * 1.[last 63 bits of significand] * 2^exponent.
        /// Exponents >= long.MaxValue - 64 and exponents &lt;= long.MinValue + 64 are reserved.
        /// Exponents of long.MinValue + 0, + 1, + 2 and + 3 are used to represent 0, PositiveInfinity, NegativeInfinity, and NaN, respectively.
        /// </summary>
        public long Exponent;
        #endregion

        #region Constructors
        /// <summary>
        /// Creates a new Quad with the given significand bits and exponent.  The significand has a first (most significant) bit
        /// corresponding to the quad's sign (1 for positive, 0 for negative), and the rest of the bits correspond to the fractional
        /// part of the significand value (immediately after the binary point).  A "1" before the binary point is always implied.
        /// </summary>
        /// <param name="significandBits"></param>
        /// <param name="exponent"></param>
        public Quad(ulong significandBits, long exponent)
        {
            SignificandBits = significandBits;
            Exponent = exponent;
        }

        /// <summary>
        /// Creates a new Quad with the given significand value and exponent.
        /// </summary>
        /// <param name="significand"></param>
        /// <param name="exponent"></param>
        public Quad(long significand, long exponent)
        {
            if (significand == 0) //handle 0
            {
                SignificandBits = 0;
                Exponent = long.MinValue;
                return;
            }

            if (significand < 0)
            {
                if (significand == long.MinValue) //corner case
                {
                    SignificandBits = HighestBit;
                    Exponent = 0;
                    return;
                }

                significand = -significand;
                SignificandBits = HighestBit;
            }
            else
                SignificandBits = 0;

            int shift = Nlz((ulong)significand); //we must normalize the value such that the most significant bit is 1
            SignificandBits |= ~HighestBit & (((ulong)significand) << shift); //mask out the highest bit--it's implicit
            Exponent = exponent - shift;
        }
        #endregion

        #region Helper functions and constants
        #region "Special" arithmetic tables for zeros, infinities, and NaN's
        //first index = first argument to the operation; second index = second argument
        //One's are used as placeholders when dividing a finite by a finite; these will not be used as the actual result of division, of course.
        //arguments are in the order: 0, positive infinity, negative infinity, NaN, positive finite, negative finite
        private static readonly Quad[,] SpecialDivisionTable = new Quad[,]{
            { NaN, Zero, Zero, NaN, Zero, Zero }, // 0 divided by something
            { PositiveInfinity, NaN, NaN, NaN, PositiveInfinity, NegativeInfinity }, // +inf divided by something
            { NegativeInfinity, NaN, NaN, NaN, NegativeInfinity, PositiveInfinity }, // -inf divided by something
            { NaN, NaN, NaN, NaN, NaN, NaN }, // NaN divided by something
            { PositiveInfinity, Zero, Zero, NaN, One, One }, //positive finite divided by something
            { NegativeInfinity, Zero, Zero, NaN, One, One } //negative finite divided by something
        };

        private static readonly Quad[,] SpecialMultiplicationTable = new Quad[,]{
            { Zero, NaN, NaN, NaN, Zero, Zero }, // 0 * something
            { NaN, PositiveInfinity, NegativeInfinity, NaN, PositiveInfinity, NegativeInfinity }, // +inf * something
            { NaN, NegativeInfinity, PositiveInfinity, NaN, NegativeInfinity, PositiveInfinity }, // -inf * something
            { NaN, NaN, NaN, NaN, NaN, NaN }, // NaN * something
            { Zero, PositiveInfinity, NegativeInfinity, NaN, One, One }, //positive finite * something
            { Zero, NegativeInfinity, PositiveInfinity, NaN, One, One } //negative finite * something
        };

        private static readonly bool[,] SpecialGreaterThanTable = new bool[,]{
            { false, false, true, false, false, true }, // 0 > something
            { true, false, true, false, true, true }, // +inf > something
            { false, false, false, false, false, false }, // -inf > something
            { false, false, false, false, false, false }, // NaN > something
            { true, false, true, false, false, true }, //positive finite > something
            { false, false, true, false, false, false } //negative finite > something
        };

        private static readonly bool[,] SpecialGreaterEqualThanTable = new bool[,]{
            { true, false, true, false, false, true }, // 0 >= something
            { true, true, true, false, true, true }, // +inf >= something
            { false, false, true, false, false, false }, // -inf >= something
            { false, false, false, false, false, false }, // NaN >= something
            { true, false, true, false, false, true }, //positive finite >= something
            { false, false, true, false, false, false } //negative finite >= something
        };

        private static readonly bool[,] SpecialLessThanTable = new bool[,]{
            { false, true, false, false, true, false }, // 0 < something
            { false, false, false, false, false, false }, // +inf < something
            { true, true, false, false, true, true }, // -inf < something
            { false, false, false, false, false, false }, // NaN < something
            { false, true, false, false, false, false }, //positive finite < something
            { true, true, false, false, true, false } //negative finite < something
        };

        private static readonly bool[,] SpecialLessEqualThanTable = new bool[,]{
            { true, true, false, false, true, false }, // 0 < something
            { false, true, false, false, false, false }, // +inf < something
            { true, true, true, false, true, true }, // -inf < something
            { false, false, false, false, false, false }, // NaN < something
            { false, true, false, false, false, false }, //positive finite < something
            { true, true, false, false, true, false } //negative finite < something
        };

        private static readonly Quad[,] SpecialSubtractionTable = new Quad[,]{
            {Zero, NegativeInfinity, PositiveInfinity, NaN, One, One}, //0 - something
            {PositiveInfinity, NaN, PositiveInfinity, NaN, PositiveInfinity, PositiveInfinity}, //+Infinity - something
            {NegativeInfinity, NegativeInfinity, NaN, NaN,NegativeInfinity,NegativeInfinity}, //-Infinity - something
            { NaN, NaN, NaN, NaN, NaN, NaN }, //NaN - something
            { One, NegativeInfinity, PositiveInfinity, NaN, One, One }, //+finite - something
            { One, NegativeInfinity, PositiveInfinity, NaN, One, One } //-finite - something
        };

        private static readonly Quad[,] SpecialAdditionTable = new Quad[,]{
            {Zero, PositiveInfinity, NegativeInfinity, NaN, One, One}, //0 + something
            {PositiveInfinity, PositiveInfinity, NaN, NaN, PositiveInfinity, PositiveInfinity}, //+Infinity + something
            {NegativeInfinity, NaN, NegativeInfinity, NaN,NegativeInfinity,NegativeInfinity}, //-Infinity + something
            { NaN, NaN, NaN, NaN, NaN, NaN }, //NaN + something
            { One, PositiveInfinity, NegativeInfinity, NaN, One, One }, //+finite + something
            { One, PositiveInfinity, NegativeInfinity, NaN, One, One } //-finite + something
        };

        private static readonly double[] SpecialDoubleLogTable = new double[] { double.NegativeInfinity, double.PositiveInfinity, double.NaN, double.NaN };

        private static readonly string[] SpecialStringTable = new string[] { "0", "Infinity", "-Infinity", "NaN" };
        #endregion

        private const long ZeroExponent = long.MinValue;
        private const long InfinityExponent = long.MinValue + 1;
        private const long NegativeInfinityExponent = long.MinValue + 2;
        private const long NotANumberExponent = long.MinValue + 3;

        private const long ExponentUpperBound = long.MaxValue - 65; //no exponent should be higher than this
        private const long ExponentLowerBound = long.MinValue + 65; //no exponent should be lower than this

        private const double Base2To10Multiplier = 0.30102999566398119521373889472449; //Math.Log(2) / Math.Log(10);
        private const ulong HighestBit = 1UL << 63;
        private const ulong SecondHighestBit = 1UL << 62;
        private const ulong LowWordMask = 0xffffffff; //lower 32 bits
        //private const ulong HighWordMask = 0xffffffff00000000; //upper 32 bits

        private const ulong B = 4294967296; // Number base (32 bits).

        private static readonly Quad E19 = (Quad)10000000000000000000UL;
        private static readonly Quad E10 = (Quad)10000000000UL;
        private static readonly Quad E5 = (Quad)100000UL;
        private static readonly Quad E3 = (Quad)1000UL;
        private static readonly Quad E1 = (Quad)10UL;

        //private static readonly Quad En19 = One / E19;
        //private static readonly Quad En10 = One / E10;
        //private static readonly Quad En5 = One / E5;
        //private static readonly Quad En3 = One / E3;
        //private static readonly Quad En1 = One / E1;

        private static readonly Quad En18 = One / (Quad)1000000000000000000UL;
        private static readonly Quad En9 = One / (Quad)1000000000UL;
        private static readonly Quad En4 = One / (Quad)10000UL;
        private static readonly Quad En2 = One / (Quad)100UL;


        /// <summary>
        /// Returns the position of the highest set bit, counting from the most significant bit position (position 0).
        /// Returns 64 if no bit is set.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private static int Nlz(ulong x)
        {
            //Future work: might be faster with a huge, explicit nested if tree, or use of an 256-element per-byte array.            

            int n;

            if (x == 0) return (64);
            n = 0;
            if (x <= 0x00000000FFFFFFFF) { n = n + 32; x = x << 32; }
            if (x <= 0x0000FFFFFFFFFFFF) { n = n + 16; x = x << 16; }
            if (x <= 0x00FFFFFFFFFFFFFF) { n = n + 8; x = x << 8; }
            if (x <= 0x0FFFFFFFFFFFFFFF) { n = n + 4; x = x << 4; }
            if (x <= 0x3FFFFFFFFFFFFFFF) { n = n + 2; x = x << 2; }
            if (x <= 0x7FFFFFFFFFFFFFFF) { n = n + 1; }
            return n;
        }

        #endregion

        #region Struct-modifying instance arithmetic functions
        public unsafe void Multiply(double multiplierDouble)
        {
            Quad multiplier;
            #region Parse the double
            // Implementation note: the use of goto is generally discouraged,
            // but here the idea is to copy-paste the casting call for double -> Quad
            // to avoid the expense of an additional function call
            // and the use of a single "return" goto target keeps things simple

            // Translate the double into sign, exponent and mantissa.
            //long bits = BitConverter.DoubleToInt64Bits(value); // doing an unsafe pointer-conversion to get the bits is faster
            ulong bits = *((ulong*)&multiplierDouble);

            // Note that the shift is sign-extended, hence the test against -1 not 1                
            long exponent = (((long)bits >> 52) & 0x7ffL);
            ulong mantissa = (bits) & 0xfffffffffffffUL;

            if (exponent == 0x7ffL)
            {
                if (mantissa == 0)
                {
                    multiplier = bits >= HighestBit ? NegativeInfinity : PositiveInfinity;
                    goto Parsed;
                }
                multiplier = NaN;
                goto Parsed;
            }

            // Subnormal numbers; exponent is effectively one higher,
            // but there's no extra normalisation bit in the mantissa
            if (exponent == 0)
            {
                if (mantissa == 0)
                {
                    multiplier = Zero;
                    goto Parsed;
                }
                exponent++;

                int firstSetPosition = Nlz(mantissa);
                mantissa <<= firstSetPosition;
                exponent -= firstSetPosition;
            }
            else
            {
                mantissa = mantissa << 11;
                exponent -= 11;
            }

            exponent -= 1075;

            multiplier.SignificandBits = (HighestBit & bits) | mantissa;
            multiplier.Exponent = exponent;

            Parsed:
            #endregion

            #region Multiply
            if (Exponent <= NotANumberExponent) //zero/infinity/NaN * something            
            {
                Quad result = SpecialMultiplicationTable[(int)(Exponent - ZeroExponent), multiplier.Exponent > NotANumberExponent ? (int)(4 + (multiplier.SignificandBits >> 63)) : (int)(multiplier.Exponent - ZeroExponent)];
                SignificandBits = result.SignificandBits;
                Exponent = result.Exponent;
                return;
            }
            if (multiplier.Exponent <= NotANumberExponent) //finite * zero/infinity/NaN            
            {
                Quad result = SpecialMultiplicationTable[(int)(4 + (SignificandBits >> 63)), (int)(multiplier.Exponent - ZeroExponent)];
                SignificandBits = result.SignificandBits;
                Exponent = result.Exponent;
                return;
            }

            ulong high1 = (SignificandBits | HighestBit) >> 32; //de-implicitize the 1
            ulong high2 = (multiplier.SignificandBits | HighestBit) >> 32;

            //because the MSB of both significands is 1, the MSB of the result will also be 1, and the product of low bits on both significands is dropped (and thus we can skip its calculation)
            ulong significandBits = high1 * high2 + (((SignificandBits & LowWordMask) * high2) >> 32) + ((high1 * (multiplier.SignificandBits & LowWordMask)) >> 32);

            long qd2Exponent;
            long qd1Exponent = Exponent;
            if (significandBits < (1UL << 63))
            {
                SignificandBits = ((SignificandBits ^ multiplier.SignificandBits) & HighestBit) | ((significandBits << 1) & ~HighestBit);
                qd2Exponent = multiplier.Exponent - 1 + 64;
                Exponent = Exponent + qd2Exponent;
            }
            else
            {
                SignificandBits = ((SignificandBits ^ multiplier.SignificandBits) & HighestBit) | (significandBits & ~HighestBit);
                qd2Exponent = multiplier.Exponent + 64;
                Exponent = Exponent + qd2Exponent;
            }

            if (qd2Exponent < 0 && Exponent > qd1Exponent) //did the exponent get larger after adding something negative?
            {
                SignificandBits = 0;
                Exponent = ZeroExponent;
            }
            else if (qd2Exponent > 0 && Exponent < qd1Exponent) //did the exponent get smaller when it should have gotten larger?
            {
                SignificandBits = 0;
                Exponent = SignificandBits >= HighestBit ? NegativeInfinityExponent : InfinityExponent; //overflow
            }
            else if (Exponent < ExponentLowerBound) //check for underflow
            {
                SignificandBits = 0;
                Exponent = ZeroExponent;
            }
            else if (Exponent > ExponentUpperBound) //overflow
            {
                SignificandBits = 0;
                Exponent = SignificandBits >= HighestBit ? NegativeInfinityExponent : InfinityExponent; //overflow
            }
            #endregion
        }

        public void Multiply(Quad multiplier)
        {
            if (Exponent <= NotANumberExponent) //zero/infinity/NaN * something            
            {
                Quad result = SpecialMultiplicationTable[(int)(Exponent - ZeroExponent), multiplier.Exponent > NotANumberExponent ? (int)(4 + (multiplier.SignificandBits >> 63)) : (int)(multiplier.Exponent - ZeroExponent)];
                SignificandBits = result.SignificandBits;
                Exponent = result.Exponent;
                return;
            }
            if (multiplier.Exponent <= NotANumberExponent) //finite * zero/infinity/NaN            
            {
                Quad result = SpecialMultiplicationTable[(int)(4 + (SignificandBits >> 63)), (int)(multiplier.Exponent - ZeroExponent)];
                SignificandBits = result.SignificandBits;
                Exponent = result.Exponent;
                return;
            }

            ulong high1 = (SignificandBits | HighestBit) >> 32; //de-implicitize the 1
            ulong high2 = (multiplier.SignificandBits | HighestBit) >> 32;

            //because the MSB of both significands is 1, the MSB of the result will also be 1, and the product of low bits on both significands is dropped (and thus we can skip its calculation)
            ulong significandBits = high1 * high2 + (((SignificandBits & LowWordMask) * high2) >> 32) + ((high1 * (multiplier.SignificandBits & LowWordMask)) >> 32);

            long qd2Exponent;
            long qd1Exponent = Exponent;
            if (significandBits < (1UL << 63))
            {
                SignificandBits = ((SignificandBits ^ multiplier.SignificandBits) & HighestBit) | ((significandBits << 1) & ~HighestBit);
                qd2Exponent = multiplier.Exponent - 1 + 64;
            }
            else
            {
                SignificandBits = ((SignificandBits ^ multiplier.SignificandBits) & HighestBit) | (significandBits & ~HighestBit);
                qd2Exponent = multiplier.Exponent + 64;
            }

            Exponent = Exponent + qd2Exponent;

            if (qd2Exponent < 0 && Exponent > qd1Exponent) //did the exponent get larger after adding something negative?
            {
                SignificandBits = 0;
                Exponent = ZeroExponent;
            }
            else if (qd2Exponent > 0 && Exponent < qd1Exponent) //did the exponent get smaller when it should have gotten larger?
            {
                SignificandBits = 0;
                Exponent = SignificandBits >= HighestBit ? NegativeInfinityExponent : InfinityExponent; //overflow
            }
            else if (Exponent < ExponentLowerBound) //check for underflow
            {
                SignificandBits = 0;
                Exponent = ZeroExponent;
            }
            else if (Exponent > ExponentUpperBound) //overflow
            {
                SignificandBits = 0;
                Exponent = SignificandBits >= HighestBit ? NegativeInfinityExponent : InfinityExponent; //overflow
            }

            #region Multiply with reduced branching (slightly faster?)
            //zeros
            ////if (this.Exponent == long.MinValue)// || multiplier.Exponent == long.MinValue)
            ////{
            ////    this.Exponent = long.MinValue;
            ////    this.Significand = 0;
            ////    return;
            ////}  

            //ulong high1 = (this.Significand | highestBit ) >> 32; //de-implicitize the 1
            //ulong high2 = (multiplier.Significand | highestBit) >> 32;

            ////because the MSB of both significands is 1, the MSB of the result will also be 1, and the product of low bits on both significands is dropped (and thus we can skip its calculation)
            //ulong significandBits = high1 * high2 + (((this.Significand & lowWordMask) * high2) >> 32) + ((high1 * (multiplier.Significand & lowWordMask)) >> 32);

            //if (significandBits < (1UL << 63)) //first bit clear?
            //{
            //    long zeroMask = ((this.Exponent ^ -this.Exponent) & (multiplier.Exponent ^ -multiplier.Exponent)) >> 63;                    
            //    this.Significand = (ulong)zeroMask & ((this.Significand ^ multiplier.Significand) & highestBit) | ((significandBits << 1) & ~highestBit);
            //    this.Exponent = (zeroMask & (this.Exponent + multiplier.Exponent - 1 + 64)) | (~zeroMask & long.MinValue);
            //}
            //else
            //{
            //    this.Significand = ((this.Significand ^ multiplier.Significand) & highestBit) | (significandBits & ~highestBit);
            //    this.Exponent = this.Exponent + multiplier.Exponent + 64;
            //}

            ////long zeroMask = ((isZeroBit1 >> 63) & (isZeroBit2 >> 63));
            ////this.Significand = (ulong)zeroMask & ((this.Significand ^ multiplier.Significand) & highestBit) | ((significandBits << (int)(1 ^ (significandBits >> 63))) & ~highestBit);
            ////this.Exponent = (zeroMask & (this.Exponent + multiplier.Exponent - 1 + 64 + (long)(significandBits >> 63))) | (~zeroMask & long.MinValue);

            #endregion
        }

        /// <summary>
        /// Multiplies this Quad by a given multiplier, but does not check for underflow or overflow in the result.
        /// This is substantially (~20%) faster than the standard Multiply() method.
        /// </summary>
        /// <param name="multiplier"></param>
        public void MultiplyUnchecked(Quad multiplier)
        {
            if (Exponent <= NotANumberExponent) //zero/infinity/NaN * something            
            {
                Quad result = SpecialMultiplicationTable[(int)(Exponent - ZeroExponent), multiplier.Exponent > NotANumberExponent ? (int)(4 + (multiplier.SignificandBits >> 63)) : (int)(multiplier.Exponent - ZeroExponent)];
                SignificandBits = result.SignificandBits;
                Exponent = result.Exponent;
                return;
            }
            if (multiplier.Exponent <= NotANumberExponent) //finite * zero/infinity/NaN            
            {
                Quad result = SpecialMultiplicationTable[(int)(4 + (SignificandBits >> 63)), (int)(multiplier.Exponent - ZeroExponent)];
                SignificandBits = result.SignificandBits;
                Exponent = result.Exponent;
                return;
            }

            ulong high1 = (SignificandBits | HighestBit) >> 32; //de-implicitize the 1
            ulong high2 = (multiplier.SignificandBits | HighestBit) >> 32;

            //because the MSB of both significands is 1, the MSB of the result will also be 1, and the product of low bits on both significands is dropped (and thus we can skip its calculation)
            ulong significandBits = high1 * high2 + (((SignificandBits & LowWordMask) * high2) >> 32) + ((high1 * (multiplier.SignificandBits & LowWordMask)) >> 32);

            long qd2Exponent;
            if (significandBits < (1UL << 63))
            {
                SignificandBits = ((SignificandBits ^ multiplier.SignificandBits) & HighestBit) | ((significandBits << 1) & ~HighestBit);
                qd2Exponent = multiplier.Exponent - 1 + 64;
                Exponent = Exponent + qd2Exponent;
            }
            else
            {
                SignificandBits = ((SignificandBits ^ multiplier.SignificandBits) & HighestBit) | (significandBits & ~HighestBit);
                qd2Exponent = multiplier.Exponent + 64;
                Exponent = Exponent + qd2Exponent;
            }
        }

        public unsafe void Add(double valueDouble)
        {
            #region Parse the double
            // Implementation note: the use of goto is generally discouraged,
            // but here the idea is to copy-paste the casting call for double -> Quad
            // to avoid the expense of an additional function call
            // and the use of a single "return" goto target keeps things simple

            Quad value;
            {
                // Translate the double into sign, exponent and mantissa.
                //long bits = BitConverter.DoubleToInt64Bits(value); // doing an unsafe pointer-conversion to get the bits is faster
                ulong bits = *((ulong*)&valueDouble);

                // Note that the shift is sign-extended, hence the test against -1 not 1                
                long exponent = (((long)bits >> 52) & 0x7ffL);
                ulong mantissa = (bits) & 0xfffffffffffffUL;

                if (exponent == 0x7ffL)
                {
                    if (mantissa == 0)
                    {
                        value = bits >= HighestBit ? NegativeInfinity : PositiveInfinity;
                        goto Parsed;
                    }
                    value = NaN;
                    goto Parsed;
                }

                // Subnormal numbers; exponent is effectively one higher,
                // but there's no extra normalisation bit in the mantissa
                if (exponent == 0)
                {
                    if (mantissa == 0)
                    {
                        value = Zero;
                        goto Parsed;
                    }
                    exponent++;

                    int firstSetPosition = Nlz(mantissa);
                    mantissa <<= firstSetPosition;
                    exponent -= firstSetPosition;
                }
                else
                {
                    mantissa = mantissa << 11;
                    exponent -= 11;
                }

                exponent -= 1075;

                value.SignificandBits = (HighestBit & bits) | mantissa;
                value.Exponent = exponent;
            }
            Parsed:
            #endregion
            #region Addition
            {
                if (Exponent <= NotANumberExponent) //zero or infinity or NaN + something
                {
                    if (Exponent == ZeroExponent)
                    {
                        SignificandBits = value.SignificandBits;
                        Exponent = value.Exponent;
                    }
                    else
                    {
                        Quad result = SpecialAdditionTable[(int)(Exponent - ZeroExponent), value.Exponent > NotANumberExponent ? (int)(4 + (value.SignificandBits >> 63)) : (int)(value.Exponent - ZeroExponent)];
                        SignificandBits = result.SignificandBits;
                        Exponent = result.Exponent;
                    }

                    return;
                }
                if (value.Exponent <= NotANumberExponent) //finite + (infinity or NaN)
                {
                    if (value.Exponent != ZeroExponent)
                    {
                        Quad result = SpecialAdditionTable[(int)(4 + (SignificandBits >> 63)), (int)(value.Exponent - ZeroExponent)];
                        SignificandBits = result.SignificandBits;
                        Exponent = result.Exponent;
                    }
                    return; //if value == 0, no need to change
                }

                if ((SignificandBits ^ value.SignificandBits) >= HighestBit) //this and value have different signs--use subtraction instead
                {
                    Subtract(new Quad(value.SignificandBits ^ HighestBit, value.Exponent));
                    return;
                }

                if (Exponent > value.Exponent)
                {
                    if (Exponent >= value.Exponent + 64)
                        return; //value too small to make a difference
                    ulong bits = (SignificandBits | HighestBit) + ((value.SignificandBits | HighestBit) >> (int)(Exponent - value.Exponent));

                    if (bits < HighestBit) //this can only happen in an overflow  
                    {
                        SignificandBits = (SignificandBits & HighestBit) | (bits >> 1);
                        Exponent = Exponent + 1;
                    }
                    else
                    {
                        SignificandBits = (SignificandBits & HighestBit) | (bits & ~HighestBit);
                        //this.Exponent = this.Exponent; //exponent stays the same
                    }
                }
                else if (Exponent < value.Exponent)
                {
                    if (value.Exponent >= Exponent + 64)
                    {
                        SignificandBits = value.SignificandBits; //too small to matter
                        Exponent = value.Exponent;
                    }
                    else
                    {
                        ulong bits = (value.SignificandBits | HighestBit) + ((SignificandBits | HighestBit) >> (int)(value.Exponent - Exponent));

                        if (bits < HighestBit) //this can only happen in an overflow                    
                        {
                            SignificandBits = (value.SignificandBits & HighestBit) | (bits >> 1);
                            Exponent = value.Exponent + 1;
                        }
                        else
                        {
                            SignificandBits = (value.SignificandBits & HighestBit) | (bits & ~HighestBit);
                            Exponent = value.Exponent;
                        }
                    }
                }
                else //expDiff == 0
                {
                    //the MSB must have the same sign, so the MSB will become 0, and logical overflow is guaranteed in this situation (so we can shift right and increment the exponent).
                    SignificandBits = ((SignificandBits + value.SignificandBits) >> 1) | (SignificandBits & HighestBit);
                    Exponent = Exponent + 1;
                }
            }
            #endregion
        }

        public void Add(Quad value)
        {
            #region Addition

            if (Exponent <= NotANumberExponent) //zero or infinity or NaN + something
            {
                if (Exponent == ZeroExponent)
                {
                    SignificandBits = value.SignificandBits;
                    Exponent = value.Exponent;
                }
                else
                {
                    Quad result = SpecialAdditionTable[(int)(Exponent - ZeroExponent), value.Exponent > NotANumberExponent ? (int)(4 + (value.SignificandBits >> 63)) : (int)(value.Exponent - ZeroExponent)];
                    SignificandBits = result.SignificandBits;
                    Exponent = result.Exponent;
                }

                return;
            }
            if (value.Exponent <= NotANumberExponent) //finite + (infinity or NaN)
            {
                if (value.Exponent != ZeroExponent)
                {
                    Quad result = SpecialAdditionTable[(int)(4 + (SignificandBits >> 63)), (int)(value.Exponent - ZeroExponent)];
                    SignificandBits = result.SignificandBits;
                    Exponent = result.Exponent;
                }
                return; //if value == 0, no need to change
            }

            if ((SignificandBits ^ value.SignificandBits) >= HighestBit) //this and value have different signs--use subtraction instead
            {
                Subtract(new Quad(value.SignificandBits ^ HighestBit, value.Exponent));
                return;
            }

            if (Exponent > value.Exponent)
            {
                if (Exponent >= value.Exponent + 64)
                    return; //value too small to make a difference
                ulong bits = (SignificandBits | HighestBit) + ((value.SignificandBits | HighestBit) >> (int)(Exponent - value.Exponent));

                if (bits < HighestBit) //this can only happen in an overflow  
                {
                    SignificandBits = (SignificandBits & HighestBit) | (bits >> 1);
                    Exponent = Exponent + 1;
                }
                else
                {
                    SignificandBits = (SignificandBits & HighestBit) | (bits & ~HighestBit);
                    //this.Exponent = this.Exponent; //exponent stays the same
                }
            }
            else if (Exponent < value.Exponent)
            {
                if (value.Exponent >= Exponent + 64)
                {
                    SignificandBits = value.SignificandBits; //too small to matter
                    Exponent = value.Exponent;
                }
                else
                {
                    ulong bits = (value.SignificandBits | HighestBit) + ((SignificandBits | HighestBit) >> (int)(value.Exponent - Exponent));

                    if (bits < HighestBit) //this can only happen in an overflow                    
                    {
                        SignificandBits = (value.SignificandBits & HighestBit) | (bits >> 1);
                        Exponent = value.Exponent + 1;
                    }
                    else
                    {
                        SignificandBits = (value.SignificandBits & HighestBit) | (bits & ~HighestBit);
                        Exponent = value.Exponent;
                    }
                }
            }
            else //expDiff == 0
            {
                //the MSB must have the same sign, so the MSB will become 0, and logical overflow is guaranteed in this situation (so we can shift right and increment the exponent).
                SignificandBits = ((SignificandBits + value.SignificandBits) >> 1) | (SignificandBits & HighestBit);
                Exponent = Exponent + 1;
            }

            #endregion
        }

        public unsafe void Subtract(double valueDouble)
        {
            #region Parse the double
            // Implementation note: the use of goto is generally discouraged,
            // but here the idea is to copy-paste the casting call for double -> Quad
            // to avoid the expense of an additional function call
            // and the use of a single "return" goto target keeps things simple

            Quad value;
            {
                // Translate the double into sign, exponent and mantissa.
                //long bits = BitConverter.DoubleToInt64Bits(value); // doing an unsafe pointer-conversion to get the bits is faster
                ulong bits = *((ulong*)&valueDouble);

                // Note that the shift is sign-extended, hence the test against -1 not 1                
                long exponent = (((long)bits >> 52) & 0x7ffL);
                ulong mantissa = (bits) & 0xfffffffffffffUL;

                if (exponent == 0x7ffL)
                {
                    if (mantissa == 0)
                    {
                        value = bits >= HighestBit ? NegativeInfinity : PositiveInfinity;
                        goto Parsed;
                    }
                    value = NaN;
                    goto Parsed;
                }

                // Subnormal numbers; exponent is effectively one higher,
                // but there's no extra normalisation bit in the mantissa
                if (exponent == 0)
                {
                    if (mantissa == 0)
                    {
                        value = Zero;
                        goto Parsed;
                    }
                    exponent++;

                    int firstSetPosition = Nlz(mantissa);
                    mantissa <<= firstSetPosition;
                    exponent -= firstSetPosition;
                }
                else
                {
                    mantissa = mantissa << 11;
                    exponent -= 11;
                }

                exponent -= 1075;

                value.SignificandBits = (HighestBit & bits) | mantissa;
                value.Exponent = exponent;
            }
            Parsed:
            #endregion

            #region Subtraction
            if (Exponent <= NotANumberExponent) //infinity or NaN - something
            {
                if (Exponent == ZeroExponent)
                {
                    SignificandBits = value.SignificandBits ^ HighestBit; //negate value
                    Exponent = value.Exponent;
                }
                else
                {
                    Quad result = SpecialSubtractionTable[(int)(Exponent - ZeroExponent), value.Exponent > NotANumberExponent ? (int)(4 + (value.SignificandBits >> 63)) : (int)(value.Exponent - ZeroExponent)];
                    SignificandBits = result.SignificandBits;
                    Exponent = result.Exponent;
                }

                return;
            }
            if (value.Exponent <= NotANumberExponent) //finite - (infinity or NaN)
            {
                if (value.Exponent != ZeroExponent)
                {
                    Quad result = SpecialSubtractionTable[(int)(4 + (SignificandBits >> 63)), (int)(value.Exponent - ZeroExponent)];
                    SignificandBits = result.SignificandBits;
                    Exponent = result.Exponent;
                }

                return;
            }

            if ((SignificandBits ^ value.SignificandBits) >= HighestBit) //this and value have different signs--use addition instead            
            {
                Add(new Quad(value.SignificandBits ^ HighestBit, value.Exponent));
                return;
            }

            if (Exponent > value.Exponent)
            {
                if (Exponent >= value.Exponent + 64)
                    return; //value too small to make a difference
                ulong bits = (SignificandBits | HighestBit) - ((value.SignificandBits | HighestBit) >> (int)(Exponent - value.Exponent));

                //make sure MSB is 1                       
                int highestBitPos = Nlz(bits);
                SignificandBits = ((bits << highestBitPos) & ~HighestBit) | (SignificandBits & HighestBit);
                Exponent = Exponent - highestBitPos;
            }
            else if (Exponent < value.Exponent) //must subtract our significand from value, and switch the sign
            {
                if (value.Exponent >= Exponent + 64)
                {
                    SignificandBits = value.SignificandBits ^ HighestBit;
                    Exponent = value.Exponent;
                    return;
                }
                ulong bits = (value.SignificandBits | HighestBit) - ((SignificandBits | HighestBit) >> (int)(value.Exponent - Exponent));

                //make sure MSB is 1                       
                int highestBitPos = Nlz(bits);
                SignificandBits = ((bits << highestBitPos) & ~HighestBit) | (~value.SignificandBits & HighestBit);
                Exponent = value.Exponent - highestBitPos;
            }
            else // (this.Exponent == value.Exponent)
            {
                if (value.SignificandBits > SignificandBits) //must switch sign
                {
                    ulong bits = value.SignificandBits - SignificandBits; //notice that we don't worry about de-implicitizing the MSB--it'd be eliminated by subtraction anyway
                    int highestBitPos = Nlz(bits);
                    SignificandBits = ((bits << highestBitPos) & ~HighestBit) | (~value.SignificandBits & HighestBit);
                    Exponent = value.Exponent - highestBitPos;
                }
                else if (value.SignificandBits < SignificandBits) //sign remains the same
                {
                    ulong bits = SignificandBits - value.SignificandBits; //notice that we don't worry about de-implicitizing the MSB--it'd be eliminated by subtraction anyway
                    int highestBitPos = Nlz(bits);
                    SignificandBits = ((bits << highestBitPos) & ~HighestBit) | (SignificandBits & HighestBit);
                    Exponent = Exponent - highestBitPos;
                }
                else //this == value
                {
                    //result is 0
                    SignificandBits = 0;
                    Exponent = ZeroExponent;
                    return;
                }
            }

            if (Exponent < ExponentLowerBound) //catch underflow
            {
                SignificandBits = 0;
                Exponent = ZeroExponent;
            }

            #endregion
        }

        public void Subtract(Quad value)
        {
            #region Subtraction
            if (Exponent <= NotANumberExponent) //infinity or NaN - something
            {
                if (Exponent == ZeroExponent)
                {
                    SignificandBits = value.SignificandBits ^ HighestBit; //negate value
                    Exponent = value.Exponent;
                }
                else
                {
                    Quad result = SpecialSubtractionTable[(int)(Exponent - ZeroExponent), value.Exponent > NotANumberExponent ? (int)(4 + (value.SignificandBits >> 63)) : (int)(value.Exponent - ZeroExponent)];
                    SignificandBits = result.SignificandBits;
                    Exponent = result.Exponent;
                }

                return;
            }
            if (value.Exponent <= NotANumberExponent) //finite - (infinity or NaN)
            {
                if (value.Exponent != ZeroExponent)
                {
                    Quad result = SpecialSubtractionTable[(int)(4 + (SignificandBits >> 63)), (int)(value.Exponent - ZeroExponent)];
                    SignificandBits = result.SignificandBits;
                    Exponent = result.Exponent;
                }

                return;
            }

            if ((SignificandBits ^ value.SignificandBits) >= HighestBit) //this and value have different signs--use addition instead            
            {
                Add(new Quad(value.SignificandBits ^ HighestBit, value.Exponent));
                return;
            }

            if (Exponent > value.Exponent)
            {
                if (Exponent >= value.Exponent + 64)
                    return; //value too small to make a difference
                ulong bits = (SignificandBits | HighestBit) - ((value.SignificandBits | HighestBit) >> (int)(Exponent - value.Exponent));

                //make sure MSB is 1                       
                int highestBitPos = Nlz(bits);
                SignificandBits = ((bits << highestBitPos) & ~HighestBit) | (SignificandBits & HighestBit);
                Exponent = Exponent - highestBitPos;
            }
            else if (Exponent < value.Exponent) //must subtract our significand from value, and switch the sign
            {
                if (value.Exponent >= Exponent + 64)
                {
                    SignificandBits = value.SignificandBits ^ HighestBit;
                    Exponent = value.Exponent;
                    return;
                }
                ulong bits = (value.SignificandBits | HighestBit) - ((SignificandBits | HighestBit) >> (int)(value.Exponent - Exponent));

                //make sure MSB is 1                       
                int highestBitPos = Nlz(bits);
                SignificandBits = ((bits << highestBitPos) & ~HighestBit) | (~value.SignificandBits & HighestBit);
                Exponent = value.Exponent - highestBitPos;
            }
            else // (this.Exponent == value.Exponent)
            {
                if (value.SignificandBits > SignificandBits) //must switch sign
                {
                    ulong bits = value.SignificandBits - SignificandBits; //notice that we don't worry about de-implicitizing the MSB--it'd be eliminated by subtraction anyway
                    int highestBitPos = Nlz(bits);
                    SignificandBits = ((bits << highestBitPos) & ~HighestBit) | (~value.SignificandBits & HighestBit);
                    Exponent = value.Exponent - highestBitPos;
                }
                else if (value.SignificandBits < SignificandBits) //sign remains the same
                {
                    ulong bits = SignificandBits - value.SignificandBits; //notice that we don't worry about de-implicitizing the MSB--it'd be eliminated by subtraction anyway
                    int highestBitPos = Nlz(bits);
                    SignificandBits = ((bits << highestBitPos) & ~HighestBit) | (SignificandBits & HighestBit);
                    Exponent = Exponent - highestBitPos;
                }
                else //this == value
                {
                    //result is 0
                    SignificandBits = 0;
                    Exponent = ZeroExponent;
                    return;
                }
            }

            if (Exponent < ExponentLowerBound) //catch underflow
            {
                SignificandBits = 0;
                Exponent = ZeroExponent;
            }

            #endregion
        }

        public unsafe void Divide(double divisorDouble)
        {
            #region Parse the double
            // Implementation note: the use of goto is generally discouraged,
            // but here the idea is to copy-paste the casting call for double -> Quad
            // to avoid the expense of an additional function call
            // and the use of a single "return" goto target keeps things simple

            Quad divisor;
            {
                // Translate the double into sign, exponent and mantissa.
                //long bits = BitConverter.DoubleToInt64Bits(divisor); // doing an unsafe pointer-conversion to get the bits is faster
                ulong bits = *((ulong*)&divisorDouble);

                // Note that the shift is sign-extended, hence the test against -1 not 1                
                long exponent = (((long)bits >> 52) & 0x7ffL);
                ulong mantissa = (bits) & 0xfffffffffffffUL;

                if (exponent == 0x7ffL)
                {
                    if (mantissa == 0)
                    {
                        divisor = bits >= HighestBit ? NegativeInfinity : PositiveInfinity;
                        goto Parsed;
                    }
                    divisor = NaN;
                    goto Parsed;
                }

                // Subnormal numbers; exponent is effectively one higher,
                // but there's no extra normalisation bit in the mantissa
                if (exponent == 0)
                {
                    if (mantissa == 0)
                    {
                        divisor = Zero;
                        goto Parsed;
                    }
                    exponent++;

                    int firstSetPosition = Nlz(mantissa);
                    mantissa <<= firstSetPosition;
                    exponent -= firstSetPosition;
                }
                else
                {
                    mantissa = mantissa << 11;
                    exponent -= 11;
                }

                exponent -= 1075;

                divisor.SignificandBits = (HighestBit & bits) | mantissa;
                divisor.Exponent = exponent;
            }
            Parsed:
            #endregion

            #region Division
            if (Exponent <= NotANumberExponent) //zero/infinity/NaN divided by something
            {
                Quad result = SpecialDivisionTable[(int)(Exponent - ZeroExponent), divisor.Exponent > NotANumberExponent ? (int)(4 + (divisor.SignificandBits >> 63)) : (int)(divisor.Exponent - ZeroExponent)];
                SignificandBits = result.SignificandBits;
                Exponent = result.Exponent;
                return;
            }
            if (divisor.Exponent <= NotANumberExponent) //finite divided by zero/infinity/NaN
            {
                Quad result = SpecialDivisionTable[(int)(4 + (SignificandBits >> 63)), (int)(divisor.Exponent - ZeroExponent)];
                SignificandBits = result.SignificandBits;
                Exponent = result.Exponent;
                return;
            }

            ulong un1 = 0,     // Norm. dividend LSD's.
                     vn1, vn0,        // Norm. divisor digits.
                     q1, q0,          // Quotient digits.
                     un21,// Dividend digit pairs.
                     rhat;            // A remainder.            

            //result.Significand = highestBit & (this.Significand ^ divisor.Significand); //determine the sign bit

            //this.Significand |= highestBit; //de-implicitize the 1 before the binary point
            //divisor.Significand |= highestBit;

            long adjExponent = 0;
            ulong thisAdjSignificand = SignificandBits | HighestBit;
            ulong divisorAdjSignificand = divisor.SignificandBits | HighestBit;

            if (thisAdjSignificand >= divisorAdjSignificand)
            {
                //need to make this's significand smaller than divisor's
                adjExponent = 1;
                un1 = (SignificandBits & 1) << 31;
                thisAdjSignificand = thisAdjSignificand >> 1;
            }

            vn1 = divisorAdjSignificand >> 32;            // Break divisor up into
            vn0 = divisor.SignificandBits & 0xFFFFFFFF;         // two 32-bit digits.            

            q1 = thisAdjSignificand / vn1;            // Compute the first
            rhat = thisAdjSignificand - q1 * vn1;     // quotient digit, q1.
            again1:
            if (q1 >= B || q1 * vn0 > B * rhat + un1)
            {
                q1 = q1 - 1;
                rhat = rhat + vn1;
                if (rhat < B) goto again1;
            }

            un21 = thisAdjSignificand * B + un1 - q1 * divisorAdjSignificand;  // Multiply and subtract.

            q0 = un21 / vn1;            // Compute the second
            rhat = un21 - q0 * vn1;     // quotient digit, q0.
            again2:
            if (q0 >= B || q0 * vn0 > B * rhat)
            {
                q0 = q0 - 1;
                rhat = rhat + vn1;
                if (rhat < B) goto again2;
            }

            thisAdjSignificand = q1 * B + q0; //convenient place to store intermediate result

            //if (this.Significand == 0) //the final significand should never be 0
            //    result.Exponent = 0;
            //else

            long originalExponent;
            long divisorExponent;

            if (thisAdjSignificand < (1UL << 63))
            {
                SignificandBits = (~HighestBit & (thisAdjSignificand << 1)) | ((SignificandBits ^ divisor.SignificandBits) & HighestBit);

                originalExponent = Exponent - 1 + adjExponent;
                divisorExponent = divisor.Exponent + 64;
            }
            else
            {
                SignificandBits = (~HighestBit & thisAdjSignificand) | ((SignificandBits ^ divisor.SignificandBits) & HighestBit);

                originalExponent = Exponent + adjExponent;
                divisorExponent = divisor.Exponent + 64;
            }

            Exponent = originalExponent - divisorExponent;

            //now check for underflow or overflow
            if (divisorExponent > 0 && Exponent > originalExponent) //underflow
            {
                SignificandBits = 0;
                Exponent = ZeroExponent; //new value is 0
            }
            else if (divisorExponent < 0 && Exponent < originalExponent) //overflow
            {
                SignificandBits = 0;// (this.SignificandBits & highestBit);
                Exponent = SignificandBits >= HighestBit ? NegativeInfinityExponent : InfinityExponent;
            }
            else if (Exponent < ExponentLowerBound)
            {
                SignificandBits = 0;
                Exponent = ZeroExponent; //new value is 0
            }
            else if (Exponent > ExponentUpperBound)
            {
                SignificandBits = 0;// (this.SignificandBits & highestBit);
                Exponent = SignificandBits >= HighestBit ? NegativeInfinityExponent : InfinityExponent;
            }

            #endregion
        }

        public void Divide(Quad divisor)
        {
            #region Division
            if (Exponent <= NotANumberExponent) //zero/infinity/NaN divided by something
            {
                Quad result = SpecialDivisionTable[(int)(Exponent - ZeroExponent), divisor.Exponent > NotANumberExponent ? (int)(4 + (divisor.SignificandBits >> 63)) : (int)(divisor.Exponent - ZeroExponent)];
                SignificandBits = result.SignificandBits;
                Exponent = result.Exponent;
                return;
            }
            if (divisor.Exponent <= NotANumberExponent) //finite divided by zero/infinity/NaN
            {
                Quad result = SpecialDivisionTable[(int)(4 + (SignificandBits >> 63)), (int)(divisor.Exponent - ZeroExponent)];
                SignificandBits = result.SignificandBits;
                Exponent = result.Exponent;
                return;
            }

            ulong un1 = 0,     // Norm. dividend LSD's.
                     vn1, vn0,        // Norm. divisor digits.
                     q1, q0,          // Quotient digits.
                     un21,// Dividend digit pairs.
                     rhat;            // A remainder.            

            //result.Significand = highestBit & (this.Significand ^ divisor.Significand); //determine the sign bit

            //this.Significand |= highestBit; //de-implicitize the 1 before the binary point
            //divisor.Significand |= highestBit;

            long adjExponent = 0;
            ulong thisAdjSignificand = SignificandBits | HighestBit;
            ulong divisorAdjSignificand = divisor.SignificandBits | HighestBit;

            if (thisAdjSignificand >= divisorAdjSignificand)
            {
                //need to make this's significand smaller than divisor's
                adjExponent = 1;
                un1 = (SignificandBits & 1) << 31;
                thisAdjSignificand = thisAdjSignificand >> 1;
            }

            vn1 = divisorAdjSignificand >> 32;            // Break divisor up into
            vn0 = divisor.SignificandBits & 0xFFFFFFFF;         // two 32-bit digits.            

            q1 = thisAdjSignificand / vn1;            // Compute the first
            rhat = thisAdjSignificand - q1 * vn1;     // quotient digit, q1.
            again1:
            if (q1 >= B || q1 * vn0 > B * rhat + un1)
            {
                q1 = q1 - 1;
                rhat = rhat + vn1;
                if (rhat < B) goto again1;
            }

            un21 = thisAdjSignificand * B + un1 - q1 * divisorAdjSignificand;  // Multiply and subtract.

            q0 = un21 / vn1;            // Compute the second
            rhat = un21 - q0 * vn1;     // quotient digit, q0.
            again2:
            if (q0 >= B || q0 * vn0 > B * rhat)
            {
                q0 = q0 - 1;
                rhat = rhat + vn1;
                if (rhat < B) goto again2;
            }

            thisAdjSignificand = q1 * B + q0; //convenient place to store intermediate result

            //if (this.Significand == 0) //the final significand should never be 0
            //    result.Exponent = 0;
            //else

            long originalExponent;
            long divisorExponent;

            if (thisAdjSignificand < (1UL << 63))
            {
                SignificandBits = (~HighestBit & (thisAdjSignificand << 1)) | ((SignificandBits ^ divisor.SignificandBits) & HighestBit);

                originalExponent = Exponent - 1 + adjExponent;
                divisorExponent = divisor.Exponent + 64;
            }
            else
            {
                SignificandBits = (~HighestBit & thisAdjSignificand) | ((SignificandBits ^ divisor.SignificandBits) & HighestBit);

                originalExponent = Exponent + adjExponent;
                divisorExponent = divisor.Exponent + 64;
            }

            Exponent = originalExponent - divisorExponent;

            //now check for underflow or overflow
            if (divisorExponent > 0 && Exponent > originalExponent) //underflow
            {
                SignificandBits = 0;
                Exponent = ZeroExponent; //new value is 0
            }
            else if (divisorExponent < 0 && Exponent < originalExponent) //overflow
            {
                SignificandBits = 0;// (this.SignificandBits & highestBit);
                Exponent = SignificandBits >= HighestBit ? NegativeInfinityExponent : InfinityExponent;
            }
            else if (Exponent < ExponentLowerBound)
            {
                SignificandBits = 0;
                Exponent = ZeroExponent; //new value is 0
            }
            else if (Exponent > ExponentUpperBound)
            {
                SignificandBits = 0;// (this.SignificandBits & highestBit);
                Exponent = SignificandBits >= HighestBit ? NegativeInfinityExponent : InfinityExponent;
            }

            #endregion
        }
        #endregion

        #region Operators
        /// <summary>
        /// Efficiently multiplies the Quad by 2^shift.
        /// </summary>
        /// <param name="qd"></param>
        /// <param name="shift"></param>
        /// <returns></returns>        
        public static Quad operator <<(Quad qd, int shift)
        {
            if (qd.Exponent <= NotANumberExponent)
                return qd; //finite * infinity == infinity, finite * NaN == NaN, finite * 0 == 0
            return new Quad(qd.SignificandBits, qd.Exponent + shift);
        }

        /// <summary>
        /// Efficiently divides the Quad by 2^shift.
        /// </summary>
        /// <param name="qd"></param>
        /// <param name="shift"></param>
        /// <returns></returns>
        public static Quad operator >>(Quad qd, int shift)
        {
            if (qd.Exponent <= NotANumberExponent)
                return qd; //infinity / finite == infinity, NaN / finite == NaN, 0 / finite == 0
            return new Quad(qd.SignificandBits, qd.Exponent - shift);
        }

        /// <summary>
        /// Efficiently multiplies the Quad by 2^shift.
        /// </summary>
        /// <param name="qd"></param>
        /// <param name="shift"></param>
        /// <returns></returns>        
        public static Quad LeftShift(Quad qd, long shift)
        {
            if (qd.Exponent <= NotANumberExponent)
                return qd; //finite * infinity == infinity, finite * NaN == NaN, finite * 0 == 0
            return new Quad(qd.SignificandBits, qd.Exponent + shift);
        }

        /// <summary>
        /// Efficiently divides the Quad by 2^shift.
        /// </summary>
        /// <param name="qd"></param>
        /// <param name="shift"></param>
        /// <returns></returns>
        public static Quad RightShift(Quad qd, long shift)
        {
            if (qd.Exponent <= NotANumberExponent)
                return qd; //infinity / finite == infinity, NaN / finite == NaN, 0 / finite == 0
            return new Quad(qd.SignificandBits, qd.Exponent - shift);
        }

        /// <summary>
        /// Divides one Quad by another and returns the result
        /// </summary>
        /// <param name="qd1"></param>
        /// <param name="qd2"></param>
        /// <returns></returns>
        /// <remarks>
        /// This code is a heavily modified derivation of a division routine given by http://www.hackersdelight.org/HDcode/divlu.c.txt ,
        /// which has a very liberal (public domain-like) license attached: http://www.hackersdelight.org/permissions.htm
        /// </remarks>
        public static Quad operator /(Quad qd1, Quad qd2)
        {
            if (qd1.Exponent <= NotANumberExponent) //zero/infinity/NaN divided by something            
                return SpecialDivisionTable[(int)(qd1.Exponent - ZeroExponent), qd2.Exponent > NotANumberExponent ? (int)(4 + (qd2.SignificandBits >> 63)) : (int)(qd2.Exponent - ZeroExponent)];
            if (qd2.Exponent <= NotANumberExponent) //finite divided by zero/infinity/NaN            
                return SpecialDivisionTable[(int)(4 + (qd1.SignificandBits >> 63)), (int)(qd2.Exponent - ZeroExponent)];

            if (qd2.Exponent == long.MinValue)
                throw new DivideByZeroException();
            if (qd1.Exponent == long.MinValue)
                return Zero;

            ulong un1 = 0,     // Norm. dividend LSD's.
                     vn1, vn0,        // Norm. divisor digits.
                     q1, q0,          // Quotient digits.
                     un21,// Dividend digit pairs.
                     rhat;            // A remainder.                        

            long adjExponent = 0;
            ulong qd1AdjSignificand = qd1.SignificandBits | HighestBit;  //de-implicitize the 1 before the binary point
            ulong qd2AdjSignificand = qd2.SignificandBits | HighestBit;  //de-implicitize the 1 before the binary point

            if (qd1AdjSignificand >= qd2AdjSignificand)
            {
                // need to make qd1's significand smaller than qd2's
                // If we were faithful to the original code this method derives from,
                // we would branch on qd1AdjSignificand > qd2AdjSignificand instead.
                // However, this results in undesirable results like (in binary) 11/11 = 0.11111...,
                // where the result should be 1.0.  Thus, we branch on >=, which prevents this problem.
                adjExponent = 1;
                un1 = (qd1.SignificandBits & 1) << 31;
                qd1AdjSignificand = qd1AdjSignificand >> 1;
            }

            vn1 = qd2AdjSignificand >> 32;            // Break divisor up into
            vn0 = qd2.SignificandBits & 0xFFFFFFFF;         // two 32-bit digits.            

            q1 = qd1AdjSignificand / vn1;            // Compute the first
            rhat = qd1AdjSignificand - q1 * vn1;     // quotient digit, q1.
            again1:
            if (q1 >= B || q1 * vn0 > B * rhat + un1)
            {
                q1 = q1 - 1;
                rhat = rhat + vn1;
                if (rhat < B) goto again1;
            }

            un21 = qd1AdjSignificand * B + un1 - q1 * qd2AdjSignificand;  // Multiply and subtract.

            q0 = un21 / vn1;            // Compute the second
            rhat = un21 - q0 * vn1;     // quotient digit, q0.
            again2:
            if (q0 >= B || q0 * vn0 > B * rhat)
            {
                q0 = q0 - 1;
                rhat = rhat + vn1;
                if (rhat < B) goto again2;
            }

            qd1AdjSignificand = q1 * B + q0; //convenient place to store intermediate result

            //if (qd1.Significand == 0) //the final significand should never be 0
            //    result.Exponent = 0;
            //else

            long originalExponent;
            long divisorExponent;
            Quad result;

            if (qd1AdjSignificand < (1UL << 63))
            {
                result.SignificandBits = (~HighestBit & (qd1AdjSignificand << 1)) | ((qd1.SignificandBits ^ qd2.SignificandBits) & HighestBit);

                originalExponent = qd1.Exponent - 1 + adjExponent;
                divisorExponent = qd2.Exponent + 64;
            }
            else
            {
                result.SignificandBits = (~HighestBit & qd1AdjSignificand) | ((qd1.SignificandBits ^ qd2.SignificandBits) & HighestBit);

                originalExponent = qd1.Exponent + adjExponent;
                divisorExponent = qd2.Exponent + 64;
            }

            result.Exponent = originalExponent - divisorExponent;

            //now check for underflow or overflow
            if (divisorExponent > 0 && result.Exponent > originalExponent) //underflow
                return Zero;
            if (divisorExponent < 0 && result.Exponent < originalExponent) //overflow            
                return result.SignificandBits >= HighestBit ? NegativeInfinity : PositiveInfinity;
            if (result.Exponent < ExponentLowerBound)
                return Zero;
            if (result.Exponent > ExponentUpperBound)
                return result.SignificandBits >= HighestBit ? NegativeInfinity : PositiveInfinity;
            return result;
        }

        /// <summary>
        /// Divides two numbers and gets the remainder.
        /// This is equivalent to qd1 - (qd2 * Truncate(qd1 / qd2)).
        /// </summary>
        /// <param name="qd1"></param>
        /// <param name="qd2"></param>
        /// <returns></returns>
        public static Quad operator %(Quad qd1, Quad qd2)
        {
            if (qd2.Exponent == InfinityExponent || qd2.Exponent == NegativeInfinityExponent)
            {
                if (qd1.Exponent == InfinityExponent || qd1.Exponent == NegativeInfinityExponent)
                    return NaN;
                return qd1;
            }

            return qd1 - (qd2 * Truncate(qd1 / qd2));
        }

        public static Quad operator -(Quad qd)
        {
            if (qd.Exponent <= NotANumberExponent)
            {
                if (qd.Exponent == InfinityExponent) return NegativeInfinity;
                if (qd.Exponent == NegativeInfinityExponent) return PositiveInfinity;
                return qd;
            }
            return new Quad(qd.SignificandBits ^ HighestBit, qd.Exponent); //just swap the sign bit                        
        }

        public static Quad operator +(Quad qd1, Quad qd2)
        {
            if (qd1.Exponent <= NotANumberExponent) //zero or infinity or NaN + something
            {
                if (qd1.Exponent == ZeroExponent) return qd2;
                return SpecialAdditionTable[(int)(qd1.Exponent - ZeroExponent), qd2.Exponent > NotANumberExponent ? (int)(4 + (qd2.SignificandBits >> 63)) : (int)(qd2.Exponent - ZeroExponent)];
            }
            if (qd2.Exponent <= NotANumberExponent) //finite + (infinity or NaN)
            {
                if (qd2.Exponent == ZeroExponent) return qd1;
                return SpecialAdditionTable[(int)(4 + (qd1.SignificandBits >> 63)), (int)(qd2.Exponent - ZeroExponent)];
            }

            if ((qd1.SignificandBits ^ qd2.SignificandBits) >= HighestBit) //qd1 and qd2 have different signs--use subtraction instead
            {
                return qd1 - new Quad(qd2.SignificandBits ^ HighestBit, qd2.Exponent);
            }

            Quad result;
            if (qd1.Exponent > qd2.Exponent)
            {
                if (qd1.Exponent >= qd2.Exponent + 64)
                    return qd1; //qd2 too small to make a difference
                ulong bits = (qd1.SignificandBits | HighestBit) + ((qd2.SignificandBits | HighestBit) >> (int)(qd1.Exponent - qd2.Exponent));

                if (bits < HighestBit) //this can only happen in an overflow                    
                    result = new Quad((qd1.SignificandBits & HighestBit) | (bits >> 1), qd1.Exponent + 1);
                else
                    return new Quad((qd1.SignificandBits & HighestBit) | (bits & ~HighestBit), qd1.Exponent);
            }
            else if (qd1.Exponent < qd2.Exponent)
            {
                if (qd2.Exponent >= qd1.Exponent + 64)
                    return qd2; //qd1 too small to matter
                ulong bits = (qd2.SignificandBits | HighestBit) + ((qd1.SignificandBits | HighestBit) >> (int)(qd2.Exponent - qd1.Exponent));

                if (bits < HighestBit) //this can only happen in an overflow                    
                    result = new Quad((qd2.SignificandBits & HighestBit) | (bits >> 1), qd2.Exponent + 1);
                else
                    return new Quad((qd2.SignificandBits & HighestBit) | (bits & ~HighestBit), qd2.Exponent);
            }
            else //expDiff == 0
            {
                //the MSB must have the same sign, so the MSB will become 0, and logical overflow is guaranteed in this situation (so we can shift right and increment the exponent).
                result = new Quad(((qd1.SignificandBits + qd2.SignificandBits) >> 1) | (qd1.SignificandBits & HighestBit), qd1.Exponent + 1);
            }

            if (result.Exponent > ExponentUpperBound) //overflow check
                return result.SignificandBits >= HighestBit ? NegativeInfinity : PositiveInfinity;
            return result;
        }

        public static Quad operator -(Quad qd1, Quad qd2)
        {
            if (qd1.Exponent <= NotANumberExponent) //infinity or NaN - something
            {
                if (qd1.Exponent == ZeroExponent) return -qd2;
                return SpecialSubtractionTable[(int)(qd1.Exponent - ZeroExponent), qd2.Exponent > NotANumberExponent ? (int)(4 + (qd2.SignificandBits >> 63)) : (int)(qd2.Exponent - ZeroExponent)];
            }
            if (qd2.Exponent <= NotANumberExponent) //finite - (infinity or NaN)            
            {
                if (qd2.Exponent == ZeroExponent) return qd1;
                return SpecialSubtractionTable[(int)(4 + (qd1.SignificandBits >> 63)), (int)(qd2.Exponent - ZeroExponent)];
            }

            if ((qd1.SignificandBits ^ qd2.SignificandBits) >= HighestBit) //qd1 and qd2 have different signs--use addition instead
            {
                return qd1 + new Quad(qd2.SignificandBits ^ HighestBit, qd2.Exponent);
            }

            Quad result;
            if (qd1.Exponent > qd2.Exponent)
            {
                if (qd1.Exponent >= qd2.Exponent + 64)
                    return qd1; //qd2 too small to make a difference
                ulong bits = (qd1.SignificandBits | HighestBit) - ((qd2.SignificandBits | HighestBit) >> (int)(qd1.Exponent - qd2.Exponent));

                //make sure MSB is 1                       
                int highestBitPos = Nlz(bits);
                result = new Quad(((bits << highestBitPos) & ~HighestBit) | (qd1.SignificandBits & HighestBit), qd1.Exponent - highestBitPos);
            }
            else if (qd1.Exponent < qd2.Exponent) //must subtract qd1's significand from qd2, and switch the sign
            {
                if (qd2.Exponent >= qd1.Exponent + 64)
                    return new Quad(qd2.SignificandBits ^ HighestBit, qd2.Exponent); //qd1 too small to matter, switch sign of qd2 and return

                ulong bits = (qd2.SignificandBits | HighestBit) - ((qd1.SignificandBits | HighestBit) >> (int)(qd2.Exponent - qd1.Exponent));

                //make sure MSB is 1                       
                int highestBitPos = Nlz(bits);
                result = new Quad(((bits << highestBitPos) & ~HighestBit) | (~qd2.SignificandBits & HighestBit), qd2.Exponent - highestBitPos);
            }
            else // (qd1.Exponent == qd2.Exponent)
            {
                if (qd2.SignificandBits > qd1.SignificandBits) //must switch sign
                {
                    ulong bits = qd2.SignificandBits - qd1.SignificandBits; //notice that we don't worry about de-implicitizing the MSB--it'd be eliminated by subtraction anyway
                    int highestBitPos = Nlz(bits);
                    result = new Quad(((bits << highestBitPos) & ~HighestBit) | (~qd2.SignificandBits & HighestBit), qd2.Exponent - highestBitPos);
                }
                else if (qd2.SignificandBits < qd1.SignificandBits) //sign remains the same
                {
                    ulong bits = qd1.SignificandBits - qd2.SignificandBits; //notice that we don't worry about de-implicitizing the MSB--it'd be eliminated by subtraction anyway
                    int highestBitPos = Nlz(bits);
                    result = new Quad(((bits << highestBitPos) & ~HighestBit) | (qd1.SignificandBits & HighestBit), qd1.Exponent - highestBitPos);
                }
                else //qd1 == qd2
                    return Zero;
            }

            if (result.Exponent < ExponentLowerBound) //handle underflow
                return Zero;
            return result;
        }

        public static Quad operator *(Quad qd1, Quad qd2)
        {
            if (qd1.Exponent <= NotANumberExponent) //zero/infinity/NaN * something            
                return SpecialMultiplicationTable[(int)(qd1.Exponent - ZeroExponent), qd2.Exponent > NotANumberExponent ? (int)(4 + (qd2.SignificandBits >> 63)) : (int)(qd2.Exponent - ZeroExponent)];
            if (qd2.Exponent <= NotANumberExponent) //finite * zero/infinity/NaN            
                return SpecialMultiplicationTable[(int)(4 + (qd1.SignificandBits >> 63)), (int)(qd2.Exponent - ZeroExponent)];

            ulong high1 = (qd1.SignificandBits | HighestBit) >> 32; //de-implicitize the 1
            ulong high2 = (qd2.SignificandBits | HighestBit) >> 32;

            //because the MSB of both significands is 1, the MSB of the result will also be 1, and the product of low bits on both significands is dropped (and thus we can skip its calculation)
            ulong significandBits = high1 * high2 + (((qd1.SignificandBits & LowWordMask) * high2) >> 32) + ((high1 * (qd2.SignificandBits & LowWordMask)) >> 32);

            long qd2Exponent;
            Quad result;
            if (significandBits < (1UL << 63))
            {
                qd2Exponent = qd2.Exponent - 1 + 64;
                result = new Quad(((qd1.SignificandBits ^ qd2.SignificandBits) & HighestBit) | ((significandBits << 1) & ~HighestBit), qd1.Exponent + qd2Exponent);
            }
            else
            {
                qd2Exponent = qd2.Exponent + 64;
                result = new Quad(((qd1.SignificandBits ^ qd2.SignificandBits) & HighestBit) | (significandBits & ~HighestBit), qd1.Exponent + qd2Exponent);
            }

            if (qd2Exponent < 0 && result.Exponent > qd1.Exponent) //did the exponent get larger after adding something negative?
                return Zero; //underflow
            if (qd2Exponent > 0 && result.Exponent < qd1.Exponent) //did the exponent get smaller when it should have gotten larger?
                return result.SignificandBits >= HighestBit ? NegativeInfinity : PositiveInfinity; //overflow
            if (result.Exponent < ExponentLowerBound) //check for underflow
                return Zero;
            if (result.Exponent > ExponentUpperBound) //overflow
                return result.SignificandBits >= HighestBit ? NegativeInfinity : PositiveInfinity; //overflow
            return result;
        }

        public static Quad operator ++(Quad qd)
        {
            return qd + One;
        }

        public static Quad operator --(Quad qd)
        {
            return qd - One;
        }
        #endregion

        #region Comparison
        /// <summary>
        /// Determines if qd1 is the same value as qd2.
        /// The same rules for doubles are used, e.g. PositiveInfinity == PositiveInfinity, but NaN != NaN.
        /// </summary>        
        public static bool operator ==(Quad qd1, Quad qd2)
        {
            return (qd1.SignificandBits == qd2.SignificandBits && qd1.Exponent == qd2.Exponent && qd1.Exponent != NotANumberExponent);// || (qd1.Exponent == long.MinValue && qd2.Exponent == long.MinValue);
        }

        /// <summary>
        /// Determines if qd1 is different from qd2.
        /// Always true if qd1 or qd2 is NaN.  False if both qd1 and qd2 are infinity with the same polarity (e.g. PositiveInfinities).
        /// </summary>        
        public static bool operator !=(Quad qd1, Quad qd2)
        {
            return (qd1.SignificandBits != qd2.SignificandBits || qd1.Exponent != qd2.Exponent || qd1.Exponent == NotANumberExponent);// && (qd1.Exponent != long.MinValue || qd2.Exponent != long.MinValue);
        }

        public static bool operator >(Quad qd1, Quad qd2)
        {
            if (qd1.Exponent <= NotANumberExponent) //zero/infinity/NaN * something            
                return SpecialGreaterThanTable[(int)(qd1.Exponent - ZeroExponent), qd2.Exponent > NotANumberExponent ? (int)(4 + (qd2.SignificandBits >> 63)) : (int)(qd2.Exponent - ZeroExponent)];
            if (qd2.Exponent <= NotANumberExponent) //finite * zero/infinity/NaN            
                return SpecialGreaterThanTable[(int)(4 + (qd1.SignificandBits >> 63)), (int)(qd2.Exponent - ZeroExponent)];

            //There is probably a faster way to accomplish this by cleverly exploiting signed longs
            switch ((qd1.SignificandBits & HighestBit) | ((qd2.SignificandBits & HighestBit) >> 1))
            {
                case HighestBit: //qd1 is negative, qd2 positive
                    return false;
                case SecondHighestBit: //qd1 positive, qd2 negative
                    return true;
                case HighestBit | SecondHighestBit: //both negative
                    return qd1.Exponent < qd2.Exponent || (qd1.Exponent == qd2.Exponent && qd1.SignificandBits < qd2.SignificandBits);
                default: //both positive
                    return qd1.Exponent > qd2.Exponent || (qd1.Exponent == qd2.Exponent && qd1.SignificandBits > qd2.SignificandBits);
            }
        }

        public static bool operator <(Quad qd1, Quad qd2)
        {
            if (qd1.Exponent <= NotANumberExponent) //zero/infinity/NaN * something            
                return SpecialLessThanTable[(int)(qd1.Exponent - ZeroExponent), qd2.Exponent > NotANumberExponent ? (int)(4 + (qd2.SignificandBits >> 63)) : (int)(qd2.Exponent - ZeroExponent)];
            if (qd2.Exponent <= NotANumberExponent) //finite * zero/infinity/NaN            
                return SpecialLessThanTable[(int)(4 + (qd1.SignificandBits >> 63)), (int)(qd2.Exponent - ZeroExponent)];

            switch ((qd1.SignificandBits & HighestBit) | ((qd2.SignificandBits & HighestBit) >> 1))
            {
                case HighestBit: //qd1 is negative, qd2 positive
                    return true;
                case SecondHighestBit: //qd1 positive, qd2 negative
                    return false;
                case HighestBit | SecondHighestBit: //both negative
                    return qd1.Exponent > qd2.Exponent || (qd1.Exponent == qd2.Exponent && qd1.SignificandBits > qd2.SignificandBits);
                default: //both positive
                    return qd1.Exponent < qd2.Exponent || (qd1.Exponent == qd2.Exponent && qd1.SignificandBits < qd2.SignificandBits);
            }

        }

        public static bool operator >=(Quad qd1, Quad qd2)
        {
            if (qd1.Exponent <= NotANumberExponent) //zero/infinity/NaN * something            
                return SpecialGreaterEqualThanTable[(int)(qd1.Exponent - ZeroExponent), qd2.Exponent > NotANumberExponent ? (int)(4 + (qd2.SignificandBits >> 63)) : (int)(qd2.Exponent - ZeroExponent)];
            if (qd2.Exponent <= NotANumberExponent) //finite * zero/infinity/NaN            
                return SpecialGreaterEqualThanTable[(int)(4 + (qd1.SignificandBits >> 63)), (int)(qd2.Exponent - ZeroExponent)];

            switch ((qd1.SignificandBits & HighestBit) | ((qd2.SignificandBits & HighestBit) >> 1))
            {
                case HighestBit: //qd1 is negative, qd2 positive
                    return false;
                case SecondHighestBit: //qd1 positive, qd2 negative
                    return true;
                case HighestBit | SecondHighestBit: //both negative
                    return qd1.Exponent < qd2.Exponent || (qd1.Exponent == qd2.Exponent && qd1.SignificandBits <= qd2.SignificandBits);
                default: //both positive
                    return qd1.Exponent > qd2.Exponent || (qd1.Exponent == qd2.Exponent && qd1.SignificandBits >= qd2.SignificandBits);
            }
        }

        public static bool operator <=(Quad qd1, Quad qd2)
        {
            if (qd1.Exponent <= NotANumberExponent) //zero/infinity/NaN * something            
                return SpecialLessEqualThanTable[(int)(qd1.Exponent - ZeroExponent), qd2.Exponent > NotANumberExponent ? (int)(4 + (qd2.SignificandBits >> 63)) : (int)(qd2.Exponent - ZeroExponent)];
            if (qd2.Exponent <= NotANumberExponent) //finite * zero/infinity/NaN            
                return SpecialLessEqualThanTable[(int)(4 + (qd1.SignificandBits >> 63)), (int)(qd2.Exponent - ZeroExponent)];

            switch ((qd1.SignificandBits & HighestBit) | ((qd2.SignificandBits & HighestBit) >> 1))
            {
                case HighestBit: //qd1 is negative, qd2 positive
                    return true;
                case SecondHighestBit: //qd1 positive, qd2 negative
                    return false;
                case HighestBit | SecondHighestBit: //both negative
                    return qd1.Exponent > qd2.Exponent || (qd1.Exponent == qd2.Exponent && qd1.SignificandBits >= qd2.SignificandBits);
                default: //both positive
                    return qd1.Exponent < qd2.Exponent || (qd1.Exponent == qd2.Exponent && qd1.SignificandBits <= qd2.SignificandBits);
            }
        }
        #endregion

        #region String parsing
        /// <summary>
        /// Parses decimal number strings in the form of "1234.5678".  Does not presently handle exponential/scientific notation.
        /// </summary>
        /// <param name="number"></param>
        /// <returns></returns>
        public static Quad Parse(string number)
        {
            if (number.Equals(SpecialStringTable[1], StringComparison.OrdinalIgnoreCase))
                return PositiveInfinity;
            if (number.Equals(SpecialStringTable[2], StringComparison.OrdinalIgnoreCase))
                return NegativeInfinity;
            if (number.Equals(SpecialStringTable[3], StringComparison.OrdinalIgnoreCase))
                return NaN;

            //Can piggyback on BigInteger's parser for this, but this is inefficient.
            //Smarter way is to break the numeric string into chunks and parse each of them using long's parse method, then combine.

            bool negative = number.StartsWith("-");
            if (negative) number = number.Substring(1);

            string left = number, right = null;
            int decimalPoint = number.IndexOf('.');
            if (decimalPoint >= 0)
            {
                left = number.Substring(0, decimalPoint);
                right = number.Substring(decimalPoint + 1);
            }

            System.Numerics.BigInteger leftInt = System.Numerics.BigInteger.Parse(left);

            Quad result = (Quad)leftInt;
            if (right != null)
            {
                System.Numerics.BigInteger rightInt = System.Numerics.BigInteger.Parse(right);
                Quad fractional = (Quad)rightInt;

                // we implicitly multiplied the stuff right of the decimal point by 10^(right.length) to get an integer;
                // now we must reverse that and add this quantity to our results.
                result += fractional * (Pow(new Quad(10L, 0), -right.Length));
            }

            return negative ? -result : result;
        }
        #endregion

        #region Math functions
        /// <summary>
        /// Removes any fractional part of the provided value (rounding down for positive numbers, and rounding up for negative numbers)            
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static Quad Truncate(Quad value)
        {
            if (value.Exponent <= NotANumberExponent) return value;

            if (value.Exponent <= -64) return Zero;
            if (value.Exponent >= 0) return value;
            //clear least significant "-value.exponent" bits that come after the binary point by shifting
            return new Quad((value.SignificandBits >> (int)(-value.Exponent)) << (int)(-value.Exponent), value.Exponent);
        }

        /// <summary>
        /// Returns only the fractional part of the provided value.  Equivalent to value % 1.
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static Quad Fraction(Quad value)
        {
            if (value.Exponent >= 0) return Zero; //no fraction
            if (value.Exponent <= -64)
            {
                if (value.Exponent == InfinityExponent || value.Exponent == NegativeInfinityExponent)
                    return NaN;
                return value; //all fraction (or zero or NaN)
            }
            //clear most significant 64+value.exponent bits before the binary point
            ulong bits = (value.SignificandBits << (int)(64 + value.Exponent)) >> (int)(64 + value.Exponent);
            if (bits == 0) return Zero; //value is an integer

            int shift = Nlz(bits); //renormalize                

            return new Quad((~HighestBit & (bits << shift)) | (HighestBit & value.SignificandBits), value.Exponent - shift);
        }

        /// <summary>
        /// Calculates the log (base 2) of a Quad.            
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static double Log2(Quad value)
        {
            if (value.SignificandBits >= HighestBit) return double.NaN;
            if (value.Exponent <= NotANumberExponent) return SpecialDoubleLogTable[(int)(value.Exponent - ZeroExponent)];

            return Math.Log(value.SignificandBits | HighestBit, 2) + value.Exponent;
        }

        ///// <summary>
        ///// Calculates the natural log (base e) of a Quad.
        ///// </summary>
        ///// <param name="value"></param>
        ///// <returns></returns>
        //public static double Log(Quad value)
        //{
        //    if (value.SignificandBits >= highestBit) return double.NaN;
        //    if (value.Exponent <= notANumberExponent) return specialDoubleLogTable[(int)(value.Exponent - zeroExponent)];

        //    return Math.Log(value.SignificandBits | highestBit) + value.Exponent * 0.69314718055994530941723212145818;
        //}

        /// <summary>
        /// Raise a Quad to a given exponent.  Pow returns 1 for x^0 for all x >= 0.  An exception is thrown
        /// if 0 is raised to a negative exponent (implying division by 0), or if a negative value is raised
        /// by a non-integer exponent (yielding an imaginary number).
        /// </summary>
        /// <param name="value"></param>
        /// <param name="exponent"></param>
        /// <returns></returns>
        /// <remarks>Internally, Pow uses Math.Pow.  This effectively limits the precision of the output to a double's 53 bits.</remarks>
        public static Quad Pow(Quad value, double exponent)
        {
            if (value.Exponent <= NotANumberExponent)
            {
                //check NaN
                if (value.Exponent == NotANumberExponent || double.IsNaN(exponent)) return NaN;

                //anything ^ 0 == 1
                if (Math.Abs(exponent) < double.Epsilon) return One;

                //0 ^ y
                if (value.Exponent == ZeroExponent)
                    return exponent < 0 ? PositiveInfinity : Zero;

                //PositiveInfinity ^ y
                if (value.Exponent == InfinityExponent)
                    return exponent < 0 ? Zero : PositiveInfinity;

                if (value.Exponent == NegativeInfinityExponent)
                    return Math.Pow(double.NegativeInfinity, exponent); //lots of weird special cases
            }

            if (double.IsNaN(exponent)) return NaN;
            if (double.IsInfinity(exponent))
            {
                if (value < -2)
                    return Math.Pow(-2, exponent);
                if (value > 2)
                    return Math.Pow(2, exponent);
                return Math.Pow((double)value, exponent);
            }

            if (Math.Abs(exponent) < double.Epsilon) return One;

            if (value.SignificandBits >= HighestBit && Math.Abs(exponent % 1) > double.Epsilon)
                return NaN; //result is an imaginary number--negative value raised to non-integer exponent

            double resultSignificand = Math.Pow((double)new Quad(value.SignificandBits, -63), exponent);
            double resultExponent = (value.Exponent + 63) * exponent; //exponents multiply

            resultSignificand *= Math.Pow(2, resultExponent % 1); //push the fractional exponent into the significand

            Quad result = resultSignificand;
            result.Exponent += (long)Math.Truncate(resultExponent);

            return result;
        }

        public static Quad Max(Quad qd1, Quad qd2)
        {
            if (qd1.Exponent == NotANumberExponent) return NaN;
            return qd1 > qd2 ? qd1 : qd2;
        }

        public static Quad Min(Quad qd1, Quad qd2)
        {
            if (qd1.Exponent == NotANumberExponent) return NaN;
            return qd1 < qd2 ? qd1 : qd2;
        }

        public static Quad Abs(Quad qd)
        {
            if (qd.Exponent == NegativeInfinityExponent) return PositiveInfinity;
            return new Quad(qd.SignificandBits & ~HighestBit, qd.Exponent); //clear the sign bit
        }

        public static int Sign(Quad qd)
        {
            if (qd.Exponent >= ExponentLowerBound) //regular number
            {
                if ((qd.SignificandBits & HighestBit) == 0) //positive?
                    return 1;
                return -1;
            }
            if (qd.Exponent == ZeroExponent) return 0;
            if (qd.Exponent == InfinityExponent) return 1;
            if (qd.Exponent == NegativeInfinityExponent) return -1;
            throw new ArithmeticException("Cannot find the Sign of a Quad that is NaN");
        }
        #endregion

        #region IsInfinity/IsNaN static test methods
        public static bool IsNaN(Quad quad)
        {
            return quad.Exponent == NotANumberExponent;
        }

        public static bool IsInfinity(Quad quad)
        {
            return quad.Exponent == InfinityExponent || quad.Exponent == NegativeInfinityExponent;
        }

        public static bool IsPositiveInfinity(Quad quad)
        {
            return quad.Exponent == InfinityExponent;
        }

        public static bool IsNegativeInfinity(Quad quad)
        {
            return quad.Exponent == NegativeInfinityExponent;
        }
        #endregion

        #region Casts
        public static explicit operator Quad(System.Numerics.BigInteger value)
        {
            bool positive = value.Sign >= 0;
            if (!positive)
                value = -value; //don't want 2's complement!

            if (value == System.Numerics.BigInteger.Zero)
                return Zero;

            if (value <= ulong.MaxValue) //easy
            {
                ulong bits = (ulong)value;
                int shift = Nlz(bits);
                return new Quad((bits << shift) & ~HighestBit | (positive ? 0 : HighestBit), -shift);
            }
            byte[] bytes = value.ToByteArray(); //least significant byte is first

            if (bytes[bytes.Length - 1] == 0) //appended when the MSB is set to differentiate from negative values
                return new Quad((positive ? 0 : HighestBit) | (~HighestBit & ((ulong)bytes[bytes.Length - 2] << 56 | (ulong)bytes[bytes.Length - 3] << 48 | (ulong)bytes[bytes.Length - 4] << 40 | (ulong)bytes[bytes.Length - 5] << 32 | (ulong)bytes[bytes.Length - 6] << 24 | (ulong)bytes[bytes.Length - 7] << 16 | (ulong)bytes[bytes.Length - 8] << 8 | bytes[bytes.Length - 9])), (bytes.Length - 9) * 8);
            {
                ulong bits = (ulong)bytes[bytes.Length - 1] << 56 | (ulong)bytes[bytes.Length - 2] << 48 | (ulong)bytes[bytes.Length - 3] << 40 | (ulong)bytes[bytes.Length - 4] << 32 | (ulong)bytes[bytes.Length - 5] << 24 | (ulong)bytes[bytes.Length - 6] << 16 | (ulong)bytes[bytes.Length - 7] << 8 | bytes[bytes.Length - 8];
                int shift = Nlz(bits);
                bits = (bits << shift) | (((ulong)bytes[bytes.Length - 9]) >> (8 - shift));
                return new Quad((positive ? 0 : HighestBit) | (~HighestBit & bits), (bytes.Length - 8) * 8 - shift);
            }
        }

        public static explicit operator System.Numerics.BigInteger(Quad value)
        {
            if (value.Exponent == NegativeInfinityExponent
                || value.Exponent == InfinityExponent)
                throw new InvalidCastException("Cannot cast infinity to BigInteger");
            if (value.Exponent == NotANumberExponent)
                throw new InvalidCastException("Cannot cast NaN to BigInteger");

            if (value.Exponent <= -64) //fractional or zero
                return System.Numerics.BigInteger.Zero;

            if (value.Exponent < 0)
            {
                if ((value.SignificandBits & HighestBit) == HighestBit)
                    return -new System.Numerics.BigInteger((value.SignificandBits) >> ((int)-value.Exponent));
                return new System.Numerics.BigInteger((value.SignificandBits | HighestBit) >> ((int)-value.Exponent));
            }

            if (value.Exponent > int.MaxValue) //you can presumably get a BigInteger bigger than 2^int.MaxValue bits, but you probably don't want to (it'd be several hundred MB).
                throw new InvalidCastException("BigIntegers do not permit left-shifts by more than int.MaxValue bits.  Since the exponent of the quad is more than this, the conversion cannot be performed.");

            if ((value.SignificandBits & HighestBit) == HighestBit) //negative number?
                return -(new System.Numerics.BigInteger(value.SignificandBits) << (int)value.Exponent);
            return (new System.Numerics.BigInteger(value.SignificandBits | HighestBit) << (int)value.Exponent);
        }

        public static explicit operator ulong(Quad value)
        {
            if (value.Exponent == NegativeInfinityExponent
                || value.Exponent == InfinityExponent)
                throw new InvalidCastException("Cannot cast infinity to 64-bit unsigned integer");
            if (value.Exponent == NotANumberExponent)
                throw new InvalidCastException("Cannot cast NaN to 64-bit unsigned integer");

            if (value.SignificandBits >= HighestBit) throw new ArgumentOutOfRangeException(nameof(value));

            if (value.Exponent > 0)
                throw new InvalidCastException("Value too large to fit in 64-bit unsigned integer");

            if (value.Exponent <= -64) return 0;

            return (HighestBit | value.SignificandBits) >> (int)(-value.Exponent);
        }

        public static explicit operator long(Quad value)
        {
            if (value.Exponent == NegativeInfinityExponent
                || value.Exponent == InfinityExponent)
                throw new InvalidCastException("Cannot cast infinity to 64-bit signed integer");
            if (value.Exponent == NotANumberExponent)
                throw new InvalidCastException("Cannot cast NaN to 64-bit signed integer");

            if (value.SignificandBits == HighestBit && value.Exponent == 0) //corner case
                return long.MinValue;

            if (value.Exponent >= 0)
                throw new InvalidCastException("Value too large to fit in 64-bit signed integer");

            if (value.Exponent <= -64) return 0;

            if (value.SignificandBits >= HighestBit) //negative
                return -(long)(value.SignificandBits >> (int)(-value.Exponent));
            return (long)((value.SignificandBits | HighestBit) >> (int)(-value.Exponent));
        }

        public static unsafe explicit operator double(Quad value)
        {
            switch (value.Exponent)
            {
                case ZeroExponent: return 0;
                case InfinityExponent: return double.PositiveInfinity;
                case NegativeInfinityExponent: return double.NegativeInfinity;
                case NotANumberExponent: return double.NaN;
            }

            if (value.Exponent <= -1086)
            {
                if (value.Exponent > -1086 - 52) //can create subnormal double value
                {
                    ulong bits = (value.SignificandBits & HighestBit) | ((value.SignificandBits | HighestBit) >> (int)(-value.Exponent - 1086 + 12));
                    return *((double*)&bits);
                }
                return 0;
            }
            {

                ulong bits = (ulong)(value.Exponent + 1086);
                if (bits >= 0x7ffUL) return value.SignificandBits >= HighestBit ? double.NegativeInfinity : double.PositiveInfinity; //too large

                bits = (value.SignificandBits & HighestBit) | (bits << 52) | (value.SignificandBits & (~HighestBit)) >> 11;

                return *((double*)&bits);
            }
        }

        /// <summary>
        /// Converts a 64-bit unsigned integer into a Quad.  No data can be lost, nor will any exception be thrown, by this cast;
        /// however, it is marked explicit in order to avoid ambiguity with the implicit long-to-Quad cast operator.
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static explicit operator Quad(ulong value)
        {
            if (value == 0) return Zero;
            int firstSetPosition = Nlz(value);
            return new Quad((value << firstSetPosition) & ~HighestBit, -firstSetPosition);
        }

        public static implicit operator Quad(long value)
        {
            return new Quad(value, 0);
        }

        public static unsafe implicit operator Quad(double value)
        {
            // Translate the double into sign, exponent and mantissa.
            //long bits = BitConverter.DoubleToInt64Bits(value); // doing an unsafe pointer-conversion to get the bits is faster
            ulong bits = *((ulong*)&value);

            // Note that the shift is sign-extended, hence the test against -1 not 1                
            long exponent = (((long)bits >> 52) & 0x7ffL);
            ulong mantissa = (bits) & 0xfffffffffffffUL;

            if (exponent == 0x7ffL)
            {
                if (mantissa == 0)
                {
                    if (bits >= HighestBit) //sign bit set?
                        return NegativeInfinity;
                    return PositiveInfinity;
                }
                return NaN;
            }

            // Subnormal numbers; exponent is effectively one higher,
            // but there's no extra normalisation bit in the mantissa
            if (exponent == 0)
            {
                if (mantissa == 0) return Zero;
                exponent++;

                int firstSetPosition = Nlz(mantissa);
                mantissa <<= firstSetPosition;
                exponent -= firstSetPosition;
            }
            else
            {
                mantissa = mantissa << 11;
                exponent -= 11;
            }

            exponent -= 1075;

            return new Quad((HighestBit & bits) | mantissa, exponent);
        }
        #endregion

        #region ToString
        /// <summary>
        /// Returns this number as a decimal, or in scientific notation where a decimal would be excessively long.
        /// Equivalent to ToString(QuadrupleFormat.ScientificApproximate).
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return ToString(QuadrupleStringFormat.ScientificExact);
        }

        /// <summary>
        /// Obtains a string representation for this Quad according to the specified format.
        /// </summary>
        /// <param name="format"></param>
        /// <returns></returns>
        /// <remarks>
        /// ScientificExact returns the value in scientific notation as accurately as possible, but is still subject to imprecision due to the conversion from 
        /// binary to decimal and the divisions or multiplications used in the conversion.  It does not use rounding, which can lead to odd-looking outputs
        /// that would otherwise be rounded by double.ToString() or the ScientificApproximate format (which uses double.ToString()).  For example, 0.1 will be rendered
        /// as the string "9.9999999999999999981e-2".
        /// </remarks>
        public string ToString(QuadrupleStringFormat format)
        {
            if (Exponent <= NotANumberExponent) return SpecialStringTable[(int)(Exponent - ZeroExponent)];

            switch (format)
            {
                case QuadrupleStringFormat.HexExponential:
                    if (SignificandBits >= HighestBit)
                        return "-" + SignificandBits.ToString("x") + "*2^" + (Exponent >= 0 ? Exponent.ToString("x") : "-" + (-Exponent).ToString("x"));
                    else
                        return (SignificandBits | HighestBit).ToString("x") + "*2^" + (Exponent >= 0 ? Exponent.ToString("x") : "-" + (-Exponent).ToString("x"));

                case QuadrupleStringFormat.DecimalExponential:
                    if (SignificandBits >= HighestBit)
                        return "-" + SignificandBits.ToString() + "*2^" + Exponent.ToString();
                    else
                        return (SignificandBits | HighestBit).ToString() + "*2^" + Exponent.ToString();

                case QuadrupleStringFormat.ScientificApproximate:
                    if (Exponent >= -1022 && Exponent <= 1023) //can be represented as double (albeit with a precision loss)
                        return ((double)this).ToString(System.Globalization.CultureInfo.InvariantCulture);

                    double dVal = (double)new Quad(SignificandBits, -61);
                    double dExp = Base2To10Multiplier * (Exponent + 61);

                    string sign = "";
                    if (dVal < 0)
                    {
                        sign = "-";
                        dVal = -dVal;
                    }

                    if (dExp >= 0)
                        dVal *= Math.Pow(10, (dExp % 1));
                    else
                        dVal *= Math.Pow(10, -((-dExp) % 1));

                    long iExp = (long)Math.Truncate(dExp);

                    while (dVal >= 10) { iExp++; dVal /= 10; }
                    while (dVal < 1) { iExp--; dVal *= 10; }

                    if (iExp >= -10 && iExp < 0)
                    {
                        string dValString = dVal.ToString(System.Globalization.CultureInfo.InvariantCulture);
                        if (dValString[1] != '.')
                            goto returnScientific; //unexpected formatting; use default behavior.
                        return sign + "0." + new string('0', (int)((-iExp) - 1)) + dVal.ToString(System.Globalization.CultureInfo.InvariantCulture).Remove(1, 1);
                    }
                    else if (iExp >= 0 && iExp <= 10)
                    {
                        string dValString = dVal.ToString(System.Globalization.CultureInfo.InvariantCulture);
                        if (dValString[1] != '.')
                            goto returnScientific; //unexpected formating; use default behavior.
                        dValString = dValString.Remove(1, 1);
                        if (iExp < dValString.Length - 1)
                            return sign + dValString.Substring(0, 1 + (int)iExp) + "." + dValString.Substring(1 + (int)iExp);
                        return sign + dValString + new string('0', (int)iExp - (dValString.Length - 1)) + ".0";
                    }

                    returnScientific:
                    return sign + dVal.ToString(System.Globalization.CultureInfo.InvariantCulture) + "E" + (iExp >= 0 ? "+" + iExp : iExp.ToString());

                case QuadrupleStringFormat.ScientificExact:
                    if (this == Zero) return "0";
                    if (Fraction(this) == Zero && Exponent <= 0) //integer value that we can output directly
                        return (SignificandBits >= HighestBit ? "-" : "") + ((SignificandBits | HighestBit) >> (int)(-Exponent)).ToString();

                    Quad absValue = Abs(this);

                    long e = 0;
                    if (absValue < One)
                    {
                        while (true)
                        {
                            if (absValue < En18)
                            {
                                absValue.Multiply(E19);
                                e -= 19;
                            }
                            else if (absValue < En9)
                            {
                                absValue.Multiply(E10);
                                e -= 10;
                            }
                            else if (absValue < En4)
                            {
                                absValue.Multiply(E5);
                                e -= 5;
                            }
                            else if (absValue < En2)
                            {
                                absValue.Multiply(E3);
                                e -= 3;
                            }
                            else if (absValue < One)
                            {
                                absValue.Multiply(E1);
                                e -= 1;
                            }
                            else
                                break;
                        }
                    }
                    else
                    {
                        while (true)
                        {
                            if (absValue >= E19)
                            {
                                absValue.Divide(E19);
                                e += 19;
                            }
                            else if (absValue >= E10)
                            {
                                absValue.Divide(E10);
                                e += 10;
                            }
                            else if (absValue >= E5)
                            {
                                absValue.Divide(E5);
                                e += 5;
                            }
                            else if (absValue >= E3)
                            {
                                absValue.Divide(E3);
                                e += 3;
                            }
                            else if (absValue >= E1)
                            {
                                absValue.Divide(E1);
                                e += 1;
                            }
                            else
                                break;
                        }
                    }

                    //absValue is now in the interval [1,10)
                    StringBuilder result = new StringBuilder();

                    result.Append(IntegerString(absValue, 1) + ".");

                    while ((absValue = Fraction(absValue)) > Zero)
                    {
                        absValue.Multiply(E19);
                        result.Append(IntegerString(absValue, 19));
                    }

                    string resultString = result.ToString().TrimEnd('0'); //trim excess 0's at the end
                    if (resultString[resultString.Length - 1] == '.') resultString += "0"; //e.g. 1.0 instead of 1.

                    return (SignificandBits >= HighestBit ? "-" : "") + resultString + "e" + (e >= 0 ? "+" : "") + e;

                default:
                    throw new ArgumentException("Unknown format requested");
            }
        }

        /// <summary>
        /// Retrieves the integer portion of the quad as a string,
        /// assuming that the quad's value is less than long.MaxValue.
        /// No sign ("-") is prepended to the result in the case of negative values.
        /// </summary>
        /// <returns></returns>
        private static string IntegerString(Quad quad, int digits)
        {
            if (quad.Exponent > 0) throw new ArgumentOutOfRangeException(nameof(quad));
            if (quad.Exponent <= -64) return "0";

            ulong significand = quad.SignificandBits | HighestBit; //make explicit the implicit bit
            return (significand >> (int)(-quad.Exponent)).ToString(new string('0', digits));
        }
        #endregion

        #region GetHashCode and Equals
        public override int GetHashCode()
        {
            // ReSharper disable once NonReadonlyMemberInGetHashCode
            int expHash = Exponent.GetHashCode();
            // ReSharper disable once NonReadonlyMemberInGetHashCode
            return SignificandBits.GetHashCode() ^ (expHash << 16 | expHash >> 16); //rotate expHash's bits 16 places
        }

        public override bool Equals(object obj)
        {
            if (obj == null) return false;

            try
            {
                return this == (Quad)obj;
            }
            catch
            {
                return false;
            }
        }
        #endregion

        #region IComparable<Quad>
        /// <inheritdoc />
        /// <summary>
        /// Returns 1 if this Quad is greater than the argument, or the argument is NaN; 0 if they are both equal or both NaN/PositiveInfinity/NegativeInfinity;
        /// and -1 if this Quad is less than the argument, or this Quad is NaN.
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public int CompareTo(Quad other)
        {
            if (Exponent == NotANumberExponent) //special value
            {
                return other.Exponent == NotANumberExponent ? 0 : -1; //If both NaN, return 0; otherwise, this NaN is "less than" everything else
            }
            if (other.Exponent == NotANumberExponent)
                return 1; //this non-NaN "greater" than other (NaN)

            if (this == other) return 0;
            if (this > other) return 1;
            return -1; //this < other
        }
        #endregion
    }
}
