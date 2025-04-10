package oblicz;

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

    private static void getResult(String method, double result) {
        System.out.printf("%s:\tSdx = %f%n", method, result);
    }

    public static void main(String[] args) {
        input();
        getResult("Trapezow", trapez());
        getResult("Simpsona", simpson());
    }
}
