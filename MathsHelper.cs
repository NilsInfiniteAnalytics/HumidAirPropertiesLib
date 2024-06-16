namespace HumidAirPropertiesLib
{
    internal class MathsHelper
    {
        public static double CalculateSummationOfPowerTerms(double[] coefficients, double baseValue, Func<double, double> coefficientPowerFunc)
        {
            return coefficients
                .Sum(t => t * Math.Pow(baseValue, coefficientPowerFunc(t)));
        }

        public static double CalculateSummationOfPowerTerms(double[] coefficients, double baseValue, double power)
        {
            return coefficients
                .Sum(coefficient => coefficient * Math.Pow(baseValue, power));
        }

        public static double CalculateSummationOfPowerTerms(double[] coefficients, double baseValueOne, double[] powerOne, double baseValueTwo, double[] powerTwo)
        {
            // Ensure the length of all arrays are the same
            if (coefficients.Length != powerOne.Length || coefficients.Length != powerTwo.Length)
            {
                throw new ArgumentException("The length of the coefficients array must be the same as the length of the power arrays.");
            }
            // Calculate the sum of the power terms
            return coefficients
                .Select((coefficient, index) => coefficient * Math.Pow(baseValueOne, powerOne[index]) * Math.Pow(baseValueTwo, powerTwo[index]))
                .Sum();

        }
    }
}
