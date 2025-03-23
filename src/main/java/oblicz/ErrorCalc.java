package oblicz;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;

public class ErrorCalc {

    private static void calc(double a, double deltaA,
                             double b, double deltaB,
                             double c, double deltaC) {
        DerivativeStructure varA = new DerivativeStructure(3, 1, 0, a);
        DerivativeStructure varB = new DerivativeStructure(3, 1, 1, b);
        DerivativeStructure varC = new DerivativeStructure(3, 1, 2, c);

        DerivativeStructure x = (varA.pow(2).multiply(varB)).divide(varC);
        double value = x.getValue();

        double dXda = x.getPartialDerivative(1, 0, 0);
        double dXdb = x.getPartialDerivative(0, 1, 0);
        double dXdc = x.getPartialDerivative(0, 0, 1);

        double deltaX = Math.abs(dXda) * deltaA
                + Math.abs(dXdb) * deltaB
                + Math.abs(dXdc) * deltaC;
        double errorX = deltaX / Math.abs(value);

        System.out.printf("x = %.5f (+-)%.5f%n", value, deltaX);
        System.out.printf("Błąd względny x = %.5f (%.2f%%)%n", errorX, errorX * 100);
    }

    public static void main(String[] args) {
        double a = 11.7, deltaA = 0.04;
        double b = 2.745, deltaB = 0.001;
        double c = 10.536, deltaC = 0.002;

        calc(a, deltaA, b, deltaB, c, deltaC);
    }
}
