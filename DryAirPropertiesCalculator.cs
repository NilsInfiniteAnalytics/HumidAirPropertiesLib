﻿namespace HumidAirPropertiesLib;
/// <summary>
/// Class for calculating properties of dry air using the Lemmon [2000] equation of state.
/// </summary>
public class DryAirPropertiesCalculator
{
    public static readonly double MolarMassDryAir = 0.02896546; // kg / mol

    public static readonly double MolarMassWater = 0.018015268; // kg / mol

    public static readonly double UniversalGasConstantLemmon = 8.314510; // J / mol K

    public static readonly double CelsiusZeroTemperature = 273.15; // K

    public static readonly double NormalPressure = 101325.0; // Pa

    public static readonly double DryAirMaxCondenThermTemperature = 132.6312; // K

    public static readonly double DryAirMaxCondenThermDensity = 10447.7; // mol / m^3

    public static readonly double DryAirMaxCondenThermPressure = 3.78502e6; // Pa

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

    private static readonly List<(double Nk, int ik, double jk, int lk)> DimensionlessDryAirResidualHelmholtzCoefficients =
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
    /// Calculates the ideal gas contribution to the Reduced Helmholtz energy equation from Lemmon [2000] for dry air.
    /// </summary>
    /// <param name="reciprocalReducedTemperature">Dimensionless temperature, scaled by DryAirMaxCondenThermTemperature</param>
    /// <param name="reducedDensity">Dimensionless density, scaled by DryAirMaxCondenThermDensity </param>
    /// <returns>Ideal gas dimensionless Helmholtz energy (A0 / RT)</returns>
    public static double CalculateDryAirDimensionlessIdealGasHelmholtzEnergy(double reciprocalReducedTemperature, double reducedDensity)
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
    /// Calculates the residual contribution to the Reduced Helmholtz energy equation from Lemmon [2000] for dry air.
    /// </summary>
    public static double CalculateDryAirDimensionlessResidualHelmholtzEnergy(double reciprocalReducedTemperature, double reducedDensity)
    {
        var termOneSum = 0.0;
        for (var i = 0; i < 10; i++)
        {
            var coefficient = DimensionlessDryAirResidualHelmholtzCoefficients[i];
            termOneSum += coefficient.Nk * Math.Pow(reducedDensity, coefficient.ik) * Math.Pow(reciprocalReducedTemperature, coefficient.jk);
        }
        var termTwoSum = 0.0;
        for (var i = 10; i < DimensionlessDryAirResidualHelmholtzCoefficients.Count; i++)
        {
            var coefficient = DimensionlessDryAirResidualHelmholtzCoefficients[i];
            termTwoSum += coefficient.Nk * Math.Pow(reducedDensity, coefficient.ik) * Math.Pow(reciprocalReducedTemperature, coefficient.jk) * Math.Exp(-Math.Pow(reducedDensity, coefficient.lk));
        }
        return termOneSum + termTwoSum;
    }

    /// <summary>
    /// Calculates the first derivative of the residual contribution to the Reduced Helmholtz energy equation from Lemmon [2000] for dry air.
    /// </summary>
    public static double CalculateDryAirDimensionlessResidualHelmholtzEnergyFirstDerivativeWithRespectToReducedDensity(double reciprocalReducedTemperature, double reducedDensity)
    {
        var termOneSum = 0.0;
        for (var i = 0; i < 10; i++)
        {
            var coefficient = DimensionlessDryAirResidualHelmholtzCoefficients[i];
            termOneSum += coefficient.Nk * coefficient.ik * Math.Pow(reducedDensity, (coefficient.ik - 1)) * Math.Pow(reciprocalReducedTemperature, coefficient.jk);
        }
        var termTwoSum = 0.0;
        for (var i = 10; i < DimensionlessDryAirResidualHelmholtzCoefficients.Count; i++)
        {
            var coefficient = DimensionlessDryAirResidualHelmholtzCoefficients[i];
            termTwoSum += coefficient.Nk * (coefficient.ik - coefficient.lk * Math.Pow(reducedDensity, coefficient.lk)) * Math.Pow(reducedDensity, coefficient.ik - 1) * Math.Exp(-Math.Pow(reducedDensity, coefficient.lk))
                          * Math.Pow(reciprocalReducedTemperature, coefficient.jk);
        }
        return termOneSum + termTwoSum;
    }

    /// <summary>
    /// Calculates the second derivative of the residual contribution to the Reduced Helmholtz energy equation from Lemmon [2000] for dry air.
    /// </summary>
    public static double CalculateDryAirDimensionlessResidualHelmholtzEnergySecondDerivativeWithRespectToReducedDensity(double reciprocalReducedTemperature, double reducedDensity)
    {
        var termOneSum = 0.0;
        for (var i = 0; i < 10; i++)
        {
            var coefficient = DimensionlessDryAirResidualHelmholtzCoefficients[i];
            termOneSum += coefficient.Nk * coefficient.ik * (coefficient.ik - 1) * Math.Pow(reducedDensity, (coefficient.ik - 2)) * Math.Pow(reciprocalReducedTemperature, coefficient.jk);
        }
        var termTwoSum = 0.0;
        for (var i = 10; i < DimensionlessDryAirResidualHelmholtzCoefficients.Count; i++)
        {
            var coefficient = DimensionlessDryAirResidualHelmholtzCoefficients[i];
            termTwoSum += coefficient.Nk * Math.Pow(reducedDensity, coefficient.ik - 2) * Math.Pow(reciprocalReducedTemperature, coefficient.jk)
                          * ((coefficient.ik - coefficient.lk * Math.Pow(reducedDensity, coefficient.lk)) * (coefficient.ik - 1.0 - coefficient.lk * Math.Pow(reducedDensity, coefficient.lk)) - Math.Pow(coefficient.lk, 2) * Math.Pow(reducedDensity, coefficient.lk))
                          * Math.Exp(-Math.Pow(reducedDensity, coefficient.lk));
        }
        return termOneSum + termTwoSum;
    }

    public static double CalculateDryAirDimensionlessIdealGasHelmholtzEnergyFirstDerivativeWithRespectToReciprocalReducedTemperature(double reciprocalReducedTemperature, double reducedDensity)
    {
        var termTwo = 0.0;
        for (var i = 1; i < 6; i++)
        {
            termTwo += DimensionlessDryAirIdealGasHelmholtzCoefficients[i] * (i - 4.0) * Math.Pow(reciprocalReducedTemperature, i - 5.0);
        }

        var termThree = 1.5 * DimensionlessDryAirIdealGasHelmholtzCoefficients[6] * Math.Sqrt(reciprocalReducedTemperature);

        var termFour = DimensionlessDryAirIdealGasHelmholtzCoefficients[7] / reciprocalReducedTemperature;

        var termFive = (DimensionlessDryAirIdealGasHelmholtzCoefficients[8] * DimensionlessDryAirIdealGasHelmholtzCoefficients[11]) / (Math.Exp(DimensionlessDryAirIdealGasHelmholtzCoefficients[11] * reciprocalReducedTemperature) - 1.0);

        var termSix = (DimensionlessDryAirIdealGasHelmholtzCoefficients[9] * DimensionlessDryAirIdealGasHelmholtzCoefficients[12]) / (Math.Exp(DimensionlessDryAirIdealGasHelmholtzCoefficients[12] * reciprocalReducedTemperature) - 1.0);

        const double twoThirds = 2.0 / 3.0;
        var termSevenLogTermTwo = DimensionlessDryAirIdealGasHelmholtzCoefficients[13] * reciprocalReducedTemperature;
        var expTermSevenLogTermTwo = Math.Exp(-termSevenLogTermTwo);
        var termSeven = (DimensionlessDryAirIdealGasHelmholtzCoefficients[10] * DimensionlessDryAirIdealGasHelmholtzCoefficients[13]) / (twoThirds * expTermSevenLogTermTwo + 1);

        var termSummation = termTwo + termThree + termFour + termFive + termSix + termSeven;
        return termSummation;
    }

    public static double CalculateDryAirDimensionlessIdealGasHelmholtzEnergySecondDerivativeWithRespectToReciprocalReducedTemperature(double reciprocalReducedTemperature, double reducedDensity)
    {
        var termTwo = 0.0;
        for (var i = 1; i < 6; i++)
        {
            termTwo += DimensionlessDryAirIdealGasHelmholtzCoefficients[i] * (i - 4.0) * (i - 5.0) * Math.Pow(reciprocalReducedTemperature, i - 6.0);
        }

        var termThree = 0.75 * DimensionlessDryAirIdealGasHelmholtzCoefficients[6] * Math.Pow(reciprocalReducedTemperature, -0.5);

        var termFour = -DimensionlessDryAirIdealGasHelmholtzCoefficients[7] / Math.Pow(reciprocalReducedTemperature, 2.0);

        var termFive = -(DimensionlessDryAirIdealGasHelmholtzCoefficients[8] * Math.Pow(DimensionlessDryAirIdealGasHelmholtzCoefficients[11], 2) * Math.Exp(DimensionlessDryAirIdealGasHelmholtzCoefficients[11] * reciprocalReducedTemperature))
                       / Math.Pow(Math.Exp(DimensionlessDryAirIdealGasHelmholtzCoefficients[11] * reciprocalReducedTemperature) - 1.0, 2.0);

        var termSix = -(DimensionlessDryAirIdealGasHelmholtzCoefficients[9] * Math.Pow(DimensionlessDryAirIdealGasHelmholtzCoefficients[12], 2) * Math.Exp(DimensionlessDryAirIdealGasHelmholtzCoefficients[12] * reciprocalReducedTemperature))
            / Math.Pow(Math.Exp(DimensionlessDryAirIdealGasHelmholtzCoefficients[12] * reciprocalReducedTemperature) - 1.0, 2.0);

        const double twoThirds = 2.0 / 3.0;
        var termSeven = (twoThirds * DimensionlessDryAirIdealGasHelmholtzCoefficients[10] * Math.Pow(DimensionlessDryAirIdealGasHelmholtzCoefficients[13], 2.0) * Math.Exp(-DimensionlessDryAirIdealGasHelmholtzCoefficients[13] * reciprocalReducedTemperature))
                        / Math.Pow(twoThirds * Math.Exp(-DimensionlessDryAirIdealGasHelmholtzCoefficients[13] * reciprocalReducedTemperature) + 1.0, 2.0);
        return termTwo + termThree + termFour + termFive + termSix + termSeven;
    }

    public static double CalculateDryAirDimensionlessResidualHelmholtzEnergyFirstDerivativeWithRespectToReciprocalReducedTemperature(double reciprocalReducedTemperature, double reducedDensity)
    {
        var termOneSum = 0.0;
        for (var i = 0; i < 10; i++)
        {
            var coefficient = DimensionlessDryAirResidualHelmholtzCoefficients[i];
            termOneSum += coefficient.Nk * Math.Pow(reducedDensity, coefficient.ik) * coefficient.jk * Math.Pow(reciprocalReducedTemperature, coefficient.jk - 1);
        }
        var termTwoSum = 0.0;
        for (var i = 10; i < DimensionlessDryAirResidualHelmholtzCoefficients.Count; i++)
        {
            var coefficient = DimensionlessDryAirResidualHelmholtzCoefficients[i];
            termTwoSum += coefficient.Nk * Math.Pow(reducedDensity, coefficient.ik) * coefficient.jk * Math.Pow(reciprocalReducedTemperature, coefficient.jk - 1) * Math.Exp(-Math.Pow(reducedDensity, coefficient.lk));
        }
        return termOneSum + termTwoSum;
    }

    public static double CalculateDryAirDimensionlessResidualHelmholtzEnergySecondDerivativeWithRespectToReciprocalReducedTemperature(double reciprocalReducedTemperature, double reducedDensity)
    {
        var termOneSum = 0.0;
        for (var i = 0; i < 10; i++)
        {
            var coefficient = DimensionlessDryAirResidualHelmholtzCoefficients[i];
            termOneSum += coefficient.Nk * Math.Pow(reducedDensity, coefficient.ik) * coefficient.jk * (coefficient.jk - 1) * Math.Pow(reciprocalReducedTemperature, coefficient.jk - 2);
        }
        var termTwoSum = 0.0;
        for (var i = 10; i < DimensionlessDryAirResidualHelmholtzCoefficients.Count; i++)
        {
            var coefficient = DimensionlessDryAirResidualHelmholtzCoefficients[i];
            termTwoSum += coefficient.Nk * Math.Pow(reducedDensity, coefficient.ik) * coefficient.jk * (coefficient.jk - 1) * Math.Pow(reciprocalReducedTemperature, coefficient.jk - 2) * Math.Exp(-Math.Pow(reducedDensity, coefficient.lk));
        }
        return termOneSum + termTwoSum;
    }

    public static double CalculateDryAirDimensionlessResidualHelmholtzEnergySecondMixedDerivativeWithRespectToReciprocalReducedTemperatureThenReducedDensity(double reciprocalReducedTemperature, double reducedDensity)
    {
        var termOneSum = 0.0;
        for (var i = 0; i < 10; i++)
        {
            var coefficient = DimensionlessDryAirResidualHelmholtzCoefficients[i];
            termOneSum += coefficient.Nk * coefficient.ik * Math.Pow(reducedDensity, coefficient.ik - 1) * coefficient.jk * Math.Pow(reciprocalReducedTemperature, coefficient.jk - 1);
        }
        var termTwoSum = 0.0;
        for (var i = 10; i < DimensionlessDryAirResidualHelmholtzCoefficients.Count; i++)
        {
            var coefficient = DimensionlessDryAirResidualHelmholtzCoefficients[i];
            termTwoSum += coefficient.Nk * (coefficient.ik - coefficient.lk * Math.Pow(reducedDensity, coefficient.lk)) * Math.Pow(reducedDensity, coefficient.ik - 1) * Math.Exp(-Math.Pow(reducedDensity, coefficient.lk))
                                         * coefficient.jk * Math.Pow(reciprocalReducedTemperature, coefficient.jk - 1);
        }
        return termOneSum + termTwoSum;
    }

    public static double CalculateDryAirDimensionlessHelmholtzEnergy(double reciprocalReducedTemperature, double reducedDensity)
    {
        var idealGasHelmholtzEnergy = CalculateDryAirDimensionlessIdealGasHelmholtzEnergy(reciprocalReducedTemperature, reducedDensity);
        var residualHelmholtzEnergy = CalculateDryAirDimensionlessResidualHelmholtzEnergy(reciprocalReducedTemperature, reducedDensity);
        return idealGasHelmholtzEnergy + residualHelmholtzEnergy;
    }

    public static (double, double) CalculateDryAirDensityError(double temperature, double molarDensity, double specifiedPressure)
    {
        var (reciprocalReducedTemperature, reducedDensity) = CalculateDryAirReducedTerms(temperature, molarDensity);

        var calculatedPressure = CalculateDryAirPressure(temperature, molarDensity);
        var error = calculatedPressure - specifiedPressure;

        var dAdDelta = CalculateDryAirDimensionlessResidualHelmholtzEnergyFirstDerivativeWithRespectToReducedDensity(reciprocalReducedTemperature, reducedDensity);
        var dA2dDelta2 = CalculateDryAirDimensionlessResidualHelmholtzEnergySecondDerivativeWithRespectToReducedDensity(reciprocalReducedTemperature, reducedDensity);
        var RT = UniversalGasConstantLemmon * temperature;
        // Calculate the derivative of the pressure with respect to density
        var dPdrho = RT * (1.0 + 2 * reducedDensity * dAdDelta + reducedDensity * reducedDensity * dA2dDelta2);
        return (error, dPdrho);
    }

    /// <summary>
    /// Calculates the molar density of dry air using the Lemmon [2000] equation of state.
    /// </summary>
    /// <param name="temperature">Temperature in Kelvin</param>
    /// <param name="pressure">Pressure in Pascal</param>
    /// <returns>Molar density of dry air (mol m-3)</returns>
    /// <exception cref="Exception">Root finding convergence failure exception</exception>
    public static double CalculateDryAirMolarDensity(double temperature, double pressure)
    {
        // Set up coefficients for the iterative solution. Use a midpoint root-finding algorithm to solve for the density.
        var a = 1e-20; // Value where error is always negative
        var b = 5000000.0; // Value where error is always positive

        var molarDensity = pressure / (UniversalGasConstantLemmon * temperature); // Initial guess using ideal gas law
        var densityOld = 100000.0;

        var iterations = 0;
        while (iterations < 2 || (Math.Abs(densityOld - molarDensity) > Math.Abs(1e-13 * molarDensity) && iterations < 100))
        {
            var (error, derror) = CalculateDryAirDensityError(temperature, molarDensity, pressure);
            if (error < 0)
            {
                a = molarDensity;
            }
            else
            {
                b = molarDensity;
            }
            var dDensity = -error / derror;
            densityOld = molarDensity;
            molarDensity += dDensity;
            if (molarDensity < a || molarDensity > b)
            {
                molarDensity = 0.5 * (a + b);
            }
            iterations++;
        }
        if (iterations >= 99) throw new Exception("Failed to converge");
        return molarDensity;
    }

    /// <summary>
    /// Calculates the dry air pressure using the Lemmon [2000] equation of state.
    /// </summary>
    /// <param name="temperature">Temperature in Kelvin</param>
    /// <param name="molarDensity">Molar density in mol m-3</param>
    /// <returns>The calculated pressure in Pascal (Pa)</returns>
    public static double CalculateDryAirPressure(double temperature, double molarDensity)
    {
        var (reciprocalReducedTemperature, reducedDensity) = CalculateDryAirReducedTerms(temperature, molarDensity);
        // Use calculus notation for the sake of not having to type out the full variable name
        var dAdDelta = CalculateDryAirDimensionlessResidualHelmholtzEnergyFirstDerivativeWithRespectToReducedDensity(reciprocalReducedTemperature, reducedDensity);
        return (1.0 + dAdDelta * reducedDensity) * molarDensity * (UniversalGasConstantLemmon * temperature);
    }

    /// <summary>
    /// Calculates the internal energy of dry air using the Lemmon [2000] equation of state.
    /// </summary>
    /// <param name="temperature">Required temperature in Kelvin</param>
    /// <param name="molarDensity">Optional if providing pressure, molar density in mol m-3</param>
    /// <param name="pressure">Optional if providing molar density, pressure in Pascals</param>
    /// <returns>The calculated molar internal energy in J mol-1</returns>
    /// <exception cref="ArgumentException">Throws if invalid or empty arguments for molar density or pressure are provided</exception>
    public static double CalculateDryAirMolarInternalEnergy(double temperature, double? molarDensity = null, double? pressure = null)
    {
        if (molarDensity is null && pressure is null)
        {
            throw new ArgumentException("Either molar density or pressure must be specified");
        }
        if (molarDensity is not null && pressure is not null)
        {
            throw new ArgumentException("Only one of molar density or pressure may be specified");
        }
        if (pressure is not null && molarDensity is null)
        {
            molarDensity = CalculateDryAirMolarDensity(temperature, pressure.Value);
        }
        if (molarDensity is null)
        {
            throw new ArgumentException("Molar density could not be calculated");
        }
        var (reciprocalReducedTemperature, reducedDensity) = CalculateDryAirReducedTerms(temperature, (double)molarDensity);
        var dAdTau = CalculateDryAirDimensionlessIdealGasHelmholtzEnergyFirstDerivativeWithRespectToReciprocalReducedTemperature(reciprocalReducedTemperature, reducedDensity);
        var daRdtau = CalculateDryAirDimensionlessResidualHelmholtzEnergyFirstDerivativeWithRespectToReciprocalReducedTemperature(reciprocalReducedTemperature, reducedDensity);
        var RT = UniversalGasConstantLemmon * temperature;
        return RT * (reciprocalReducedTemperature * (dAdTau + daRdtau));
    }

    /// <summary>
    /// Calculates the molar enthalpy of dry air using the Lemmon [2000] equation of state.
    /// </summary>
    /// <param name="temperature">Required temperature in Kelvin</param>
    /// <param name="molarDensity">Optional if providing pressure, molar density in mol m-3</param>
    /// <param name="pressure">Optional if providing molar density, pressure in Pascals</param>
    /// <returns>The calculated molar enthalpy in J mol-1</returns>
    /// <exception cref="ArgumentException">Throws if invalid or empty arguments for molar density or pressure are provided</exception>
    public static double CalculateDryAirMolarEnthalpy(double temperature, double? molarDensity = null, double? pressure = null)
    {
        if (molarDensity is null && pressure is null)
        {
            throw new ArgumentException("Either molar density or pressure must be specified");
        }
        if (molarDensity is not null && pressure is not null)
        {
            throw new ArgumentException("Only one of molar density or pressure may be specified");
        }
        if (pressure is not null && molarDensity is null)
        {
            molarDensity = CalculateDryAirMolarDensity(temperature, pressure.Value);
        }
        if (molarDensity is null)
        {
            throw new ArgumentException("Molar density could not be calculated");
        }
        var (reciprocalReducedTemperature, reducedDensity) = CalculateDryAirReducedTerms(temperature, (double)molarDensity);
        var RT = UniversalGasConstantLemmon * temperature;
        var internalEnergy = CalculateDryAirMolarInternalEnergy(temperature, molarDensity);
        var dArdDelta = CalculateDryAirDimensionlessResidualHelmholtzEnergyFirstDerivativeWithRespectToReducedDensity(reciprocalReducedTemperature, reducedDensity);
        var pressureEnergy = RT * (reducedDensity * dArdDelta + 1.0);
        return internalEnergy + pressureEnergy;
    }

    /// <summary>
    /// Calculates the molar entropy of dry air using the Lemmon [2000] equation of state.
    /// </summary>
    /// <param name="temperature">Required temperature in Kelvin</param>
    /// <param name="molarDensity">Optional if providing pressure, molar density in mol m-3</param>
    /// <param name="pressure">Optional if providing molar density, pressure in Pascals</param>
    /// <returns>The calculated molar entropy in J mol-1 K-1</returns>
    /// <exception cref="ArgumentException">Throws if invalid or empty arguments for molar density or pressure are provided</exception>
    public static double CalculateDryAirMolarEntropy(double temperature, double? molarDensity = null, double? pressure = null)
    {
        if (molarDensity is null && pressure is null)
        {
            throw new ArgumentException("Either molar density or pressure must be specified");
        }
        if (molarDensity is not null && pressure is not null)
        {
            throw new ArgumentException("Only one of molar density or pressure may be specified");
        }
        if (pressure is not null && molarDensity is null)
        {
            molarDensity = CalculateDryAirMolarDensity(temperature, pressure.Value);
        }
        if (molarDensity is null)
        {
            throw new ArgumentException("Molar density could not be calculated");
        }
        var (reciprocalReducedTemperature, reducedDensity) = CalculateDryAirReducedTerms(temperature, (double)molarDensity);
        var dAdTau = CalculateDryAirDimensionlessIdealGasHelmholtzEnergyFirstDerivativeWithRespectToReciprocalReducedTemperature(reciprocalReducedTemperature, reducedDensity);
        var daRdtau = CalculateDryAirDimensionlessResidualHelmholtzEnergyFirstDerivativeWithRespectToReciprocalReducedTemperature(reciprocalReducedTemperature, reducedDensity);
        var a0 = CalculateDryAirDimensionlessIdealGasHelmholtzEnergy(reciprocalReducedTemperature, reducedDensity);
        var aR = CalculateDryAirDimensionlessResidualHelmholtzEnergy(reciprocalReducedTemperature, reducedDensity);
        return UniversalGasConstantLemmon * (reciprocalReducedTemperature * (dAdTau + daRdtau) - a0 - aR);
    }

    /// <summary>
    /// Calculates the molar isochoric specific heat of dry air using the Lemmon [2000] equation of state.
    /// </summary>
    /// <param name="temperature">Required temperature in Kelvin</param>
    /// <param name="molarDensity">Optional if providing pressure, molar density in mol m-3</param>
    /// <param name="pressure">Optional if providing molar density, pressure in Pascals</param>
    /// <returns>The calculated molar isochoric specific heat in  J mol-1 K-1</returns>
    /// <exception cref="ArgumentException"></exception>
    public static double CalculateDryAirIsochoricSpecificHeat(double temperature, double? molarDensity = null, double? pressure = null)
    {
        if (molarDensity is null && pressure is null)
        {
            throw new ArgumentException("Either molar density or pressure must be specified");
        }
        if (molarDensity is not null && pressure is not null)
        {
            throw new ArgumentException("Only one of molar density or pressure may be specified");
        }
        if (pressure is not null && molarDensity is null)
        {
            molarDensity = CalculateDryAirMolarDensity(temperature, pressure.Value);
        }
        if (molarDensity is null)
        {
            throw new ArgumentException("Molar density could not be calculated");
        }
        var (reciprocalReducedTemperature, reducedDensity) = CalculateDryAirReducedTerms(temperature, (double)molarDensity);
        var d2AdTau2 = CalculateDryAirDimensionlessIdealGasHelmholtzEnergySecondDerivativeWithRespectToReciprocalReducedTemperature(reciprocalReducedTemperature, reducedDensity);
        var d2aRdtau2 = CalculateDryAirDimensionlessResidualHelmholtzEnergySecondDerivativeWithRespectToReciprocalReducedTemperature(reciprocalReducedTemperature, reducedDensity);
        return -UniversalGasConstantLemmon * Math.Pow(reciprocalReducedTemperature, 2.0) * (d2AdTau2 + d2aRdtau2);
    }

    public static double CalculateDryAirIsobaricSpecificHeat(double temperature, double? molarDensity = null, double? pressure = null)
    {
        if (molarDensity is null && pressure is null)
        {
            throw new ArgumentException("Either molar density or pressure must be specified");
        }
        if (molarDensity is not null && pressure is not null)
        {
            throw new ArgumentException("Only one of molar density or pressure may be specified");
        }
        if (pressure is not null && molarDensity is null)
        {
            molarDensity = CalculateDryAirMolarDensity(temperature, pressure.Value);
        }
        if (molarDensity is null)
        {
            throw new ArgumentException("Molar density could not be calculated");
        }
        var (reciprocalReducedTemperature, reducedDensity) = CalculateDryAirReducedTerms(temperature, (double)molarDensity);
        var isochoricSpecificHeat = CalculateDryAirIsochoricSpecificHeat(temperature, molarDensity);
        var dAdDelta = CalculateDryAirDimensionlessResidualHelmholtzEnergyFirstDerivativeWithRespectToReducedDensity(reciprocalReducedTemperature, reducedDensity);
        var d2AdDelta2 = CalculateDryAirDimensionlessResidualHelmholtzEnergySecondDerivativeWithRespectToReducedDensity(reciprocalReducedTemperature, reducedDensity);
        var d2ARdDeltadTau = CalculateDryAirDimensionlessResidualHelmholtzEnergySecondMixedDerivativeWithRespectToReciprocalReducedTemperatureThenReducedDensity(reciprocalReducedTemperature, reducedDensity);
        var numerator = 1 + reducedDensity * dAdDelta - reducedDensity * reciprocalReducedTemperature * d2ARdDeltadTau;
        var denominator = 1 + 2 * reducedDensity * dAdDelta + reducedDensity * reciprocalReducedTemperature * d2AdDelta2;
        var termTwo = Math.Pow(numerator, 2.0) / denominator;
        return isochoricSpecificHeat + UniversalGasConstantLemmon * termTwo;
    }


    /// <summary>
    /// Calculates the reduced terms for use in the dimensionless Helmholtz energy equations.
    /// </summary>
    /// <param name="temperature">Temperature in Kelvin</param>
    /// <param name="molarDensity">Molar density in mol m-3</param>
    /// <returns>A tuple containing (Reciprocal Reduced Temperature, Reduced Density)</returns>
    public static (double, double) CalculateDryAirReducedTerms(double temperature, double molarDensity)
    {
        return (DryAirMaxCondenThermTemperature / temperature, molarDensity / DryAirMaxCondenThermDensity);
    }
}