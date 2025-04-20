package oblicz;

import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

public class NumericalIntegration {

    static boolean userInput = false;

    static final double a = 1.4, b = 3.2;
    static double f(double x) {
        return (Math.sin(1.1 * x - 0.3)) / (2.5 + Math.pow(x, 2)); // 14
    }

    static int n;

    private static void input() {
        if (userInput) {
            Scanner sc = new Scanner(System.in);
            System.out.print("n = ");
            n = sc.nextInt();
        } else {
            n = 10;
        }
        System.out.println("-------------------------------------------------------------------------------------------------------------------------------------------");
        System.out.printf("a = %.1f, b = %.1f, n = %d%n", a, b, n);
        System.out.println("-------------------------------------------------------------------------------------------------------------------------------------------");
    }

    private static double trapez() {
        double h = (b - a) / n;
        double sum = 0.5 * (f(a) + f(b));

        for (int i = 1; i < n; i++) {
            sum += f(a + i * h);
        }

        return sum * h;
    }

    private static double simpson() {
        double h = (b - a) / n;
        double result = 0.0;

        for (int i = 0; i < n; i++) {
            double xi = a + i * h;
            double xi1 = xi + h;
            double ti = (xi + xi1) / 2;

            result += (h / 6.0) * (f(xi) + 4 * f(ti) + f(xi1));
        }

        return result;
    }

    private static double gaussLegendre() {
        if (!gaussData.containsKey(n)) {
            throw new IllegalArgumentException("Proszę wybrać n spośród {2, 4, 6, 8, 10, 16} żeby użyć metodę Gaussa-Legendre'a");
        }
        double[][] data = gaussData.get(n);

        double[] t = data[0];
        double[] w = data[1];
        final double base = (b - a) / 2;

        double result = 0;
        for (int i = 0; i < t.length; i++) {
            double xi = base * t[i] + ((a + b) / 2);
            result += w[i] * f(xi);
        }

        return result * base;
    }

    private static void getResult(String method, double result) {
        System.out.printf("%s:\tSdx = %.16f%n", method, result);
    }

    public static void main(String[] args) {
        input();
        getResult("Trapezow", trapez());
        getResult("Simpsona", simpson());
        getResult("Gaussa-Legendre'a", gaussLegendre());
    }

    private static final Map<Integer, double[][]> gaussData = new HashMap<>();
    static {
        gaussData.put(2, new double[][]{
                {-0.5773502691896257, 0.5773502691896257},
                { 1.0000000000000000, 1}
        });
        gaussData.put(4, new double[][]{
                {-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526},
                { 0.6521451548625461, 0.6521451548625461,  0.3478548451374538, 0.3478548451374538}
        });
        gaussData.put(6, new double[][]{
                {-0.6612093864662645, 0.6612093864662645, -0.2386191860831969, 0.2386191860831969, -0.9324695142031521, 0.9324695142031521},
                { 0.3607615730481386, 0.3607615730481386,  0.4679139345726910, 0.4679139345726910,  0.1713244923791704, 0.1713244923791704}
        });
        gaussData.put(8, new double[][]{
                {-0.1834346424956498, 0.1834346424956498, -0.5255324099163290, 0.5255324099163290, -0.7966664774136267, 0.7966664774136267, -0.9602898564975363, 0.9602898564975363},
                { 0.3626837833783620, 0.3626837833783620,  0.3137066458778873, 0.3137066458778873,  0.2223810344533745, 0.2223810344533745,  0.1012285362903763, 0.1012285362903763}
        });
        gaussData.put(10, new double[][]{
                {-0.1488743389816312, 0.1488743389816312, -0.4333953941292472, 0.4333953941292472, -0.6794095682990244, 0.6794095682990244, -0.8650633666889845, 0.8650633666889845, -0.9739065285171717, 0.9739065285171717},
                { 0.2955242247147529, 0.2955242247147529,  0.2692667193099963, 0.2692667193099963,  0.2190863625159820, 0.2190863625159820,  0.1494513491505806, 0.1494513491505806,  0.0666713443086881, 0.0666713443086881}
        });
        gaussData.put(16, new double[][]{
                {-0.0950125098376374, 0.0950125098376374, -0.2816035507792589, 0.2816035507792589, -0.4580167776572274, 0.4580167776572274, -0.6178762444026438, 0.6178762444026438, -0.7554044083550030, 0.7554044083550030, -0.8656312023878318, 0.8656312023878318, -0.9445750230732326, 0.9445750230732326, -0.9894009349916499, 0.9894009349916499},
                { 0.1894506104550685, 0.1894506104550685,  0.1826034150449236, 0.1826034150449236,  0.1691565193950025, 0.1691565193950025,  0.1495959888165767, 0.1495959888165767,  0.1246289712555339, 0.1246289712555339,  0.0951585116824928, 0.0951585116824928,  0.0622535239386479, 0.0622535239386479,  0.0271524594117541, 0.0271524594117541}
        });
    }
}
