namespace HumidAirPropertiesLib
{
    public class HumidAirPropertiesCalculator
    {
        /// <summary>
        /// From RP-1485, molar mass of dry air for use in psychrometric calculations.
        /// </summary>
        public static readonly double MolarMassDryAir = 0.02896546; // kg / mol

        public static readonly double MolarMassWater = 0.018015268; // kg / mol

        public static readonly double UniversalGasConstantLemmon = 8.31451; // J / mol K

        public static readonly double CelsiusZeroTemperature = 273.15; // K

        public static readonly double NormalPressure = 101325.0; // Pa

        public static readonly double MaxCondenThermTemperature = 132.6312; // K

        public static readonly double MaxCondenThermDensity = 10447.7; // mol / m^3

        private static readonly Dictionary<int, double> DimensionlessDryAirIdealGasHelmholtzCoefficients = new()
        {
            { 1,   0.605719400e-7 },
            { 2,  -0.210274769e-4 },
            { 3,  -0.158860716e-3 },
            { 4,  -13.841928076 },
            { 5,   17.275266575 },
            { 6,  -0.195363420e-3 },
            { 7,   2.490888032 },
            { 8,   0.791309509 },
            { 9,   0.212236768 },
            { 10, -0.197938904 },
            { 11,  25.36365 },
            { 12,  16.90741 },
            { 13,  87.31279 }
        };

        private static readonly List<(double Nk, int ik, double jk, int lk)> DimensionlessResidualHelmholtzCoefficients =
        [
            (0.118160747229, 1, 0.00, 0),
            (0.713116392079, 1, 0.33, 0),
            (-0.161824192067e1, 1, 1.01, 0),
            (0.714140178971e-1, 2, 0.00, 0),
            (-0.865421396646e-1, 3, 0.00, 0),
            (0.134211176704, 3, 0.15, 0),
            (0.112626704218e-1, 4, 0.00, 0),
            (-0.420533228842e-1, 4, 0.20, 0),
            (0.349008431982e-1, 4, 0.35, 0),
            (0.164957183186e-3, 6, 1.35, 0),
            (-0.101365037912, 1, 1.60, 1),
            (-0.173813690970, 3, 0.80, 1),
            (-0.472103183731e-1, 5, 0.95, 1),
            (-0.122523554253e-1, 6, 1.25, 1),
            (-0.146629609713, 1, 3.60, 2),
            (-0.316055879821e-1, 3, 6.00, 2),
            (0.233594806142e-3, 11, 3.25, 2),
            (0.148287891978e-1, 1, 3.50, 3),
            (-0.938782884667e-2, 3, 15.0, 3)
        ];

        /// <summary>
        /// Calculates the ideal gas contribution to the Helmholtz energy equation from Lemmon [2000] for dry air.
        /// </summary>
        /// <param name="reciprocalReducedTemperature">Dimensionless temperature, scaled by MaxCondenThermTemperature</param>
        /// <param name="reducedDensity">Dimensionless density, scaled by MaxCondenThermDensity </param>
        /// <returns>Ideal gas dimensionless Helmholtz energy (A0 / RT)</returns>
        public static double CalculateDimensionlessIdealGasHelmholtzEnergy(double reciprocalReducedTemperature, double reducedDensity)
        {
            var termOne = Math.Log(reducedDensity);

            var termTwo = 0.0;
            for (var i = 1; i < 6; i++)
            {
                termTwo += DimensionlessDryAirIdealGasHelmholtzCoefficients[i] * Math.Pow(reciprocalReducedTemperature, i - 4.0);
            }

            var termThree = DimensionlessDryAirIdealGasHelmholtzCoefficients[6] * Math.Pow(reciprocalReducedTemperature, 1.5);

            var termFour = DimensionlessDryAirIdealGasHelmholtzCoefficients[7] * Math.Log(reciprocalReducedTemperature);

            var termFive = DimensionlessDryAirIdealGasHelmholtzCoefficients[8] * Math.Log(1.0 - Math.Exp(-DimensionlessDryAirIdealGasHelmholtzCoefficients[11] * reciprocalReducedTemperature));

            var termSix = DimensionlessDryAirIdealGasHelmholtzCoefficients[9] * Math.Log(1.0 - Math.Exp(-DimensionlessDryAirIdealGasHelmholtzCoefficients[12] * reciprocalReducedTemperature));

            const double termSevenLogTermOne = 2.0 / 3.0;
            var termSevenLogTermTwo = DimensionlessDryAirIdealGasHelmholtzCoefficients[13] * reciprocalReducedTemperature;
            var expTermSevenLogTermTwo = Math.Exp(termSevenLogTermTwo);
            var termSeven = DimensionlessDryAirIdealGasHelmholtzCoefficients[10] * Math.Log(termSevenLogTermOne + expTermSevenLogTermTwo);

            var termSummation = termOne + termTwo + termThree + termFour + termFive + termSix + termSeven;
            return termSummation;
        }

        /// <summary>
        /// Calculates the residual contribution to the Helmholtz energy equation from Lemmon [2000] for dry air.
        /// </summary>
        public static double CalculateDimensionlessResidualHelmholtzEnergy(double reciprocalReducedTemperature, double reducedDensity)
        {
            var termOneSum = 0.0;
            for (var i = 0; i < 10; i++)
            {
                var coefficient = DimensionlessResidualHelmholtzCoefficients[i];
                termOneSum += coefficient.Nk * Math.Pow(reducedDensity, coefficient.ik) * Math.Pow(reciprocalReducedTemperature, coefficient.jk);
            }
            var termTwoSum = 0.0;
            for (var i = 10; i < DimensionlessResidualHelmholtzCoefficients.Count; i++)
            {
                var coefficient = DimensionlessResidualHelmholtzCoefficients[i];
                termTwoSum += coefficient.Nk * Math.Pow(reducedDensity, coefficient.ik) * Math.Pow(reciprocalReducedTemperature, coefficient.jk) * Math.Exp(-Math.Pow(reducedDensity, coefficient.lk));
            }
            return termOneSum + termTwoSum;
        }

        public static double CalculateDimensionlessHelmholtzEnergy(double reciprocalReducedTemperature, double reducedDensity)
        {
            var idealGasHelmholtzEnergy = CalculateDimensionlessIdealGasHelmholtzEnergy(reciprocalReducedTemperature, reducedDensity);
            var residualHelmholtzEnergy = CalculateDimensionlessResidualHelmholtzEnergy(reciprocalReducedTemperature, reducedDensity);
            return idealGasHelmholtzEnergy + residualHelmholtzEnergy;
        }

        public static double CalculateDryAirSpecificHelmholtzEnergy(double temperature, double density)
        {
            var reciprocalReducedTemperature = MaxCondenThermTemperature / temperature;
            var molarDensity = density / MolarMassDryAir;
            var reducedDensity = molarDensity / MaxCondenThermDensity;
            var dimensionlessHelmholtzEnergy = CalculateDimensionlessHelmholtzEnergy(reciprocalReducedTemperature, reducedDensity);
            var scaleConstant = UniversalGasConstantLemmon * temperature / MolarMassDryAir;
            return scaleConstant * dimensionlessHelmholtzEnergy;
        }
    }
}
