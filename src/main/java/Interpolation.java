import java.util.*;

public class Interpolation {

    static boolean userInput = false;

    static Map<Double, Double> pairs = new HashMap<>();
    static double x;

    private static void input() {
        if (userInput) {
            Scanner sc = new Scanner(System.in);

            for (int i = 0; i < 10000; i++) {
                System.out.printf("x[%d] = ", i);
                double x = sc.nextDouble();
                System.out.printf("f(%.1f) = ", x);
                double y = sc.nextDouble();
                pairs.put(x, y);

                System.out.print("Continue? Y/n: ");
                String s = sc.next();
                if (!s.toLowerCase().startsWith("y")) break;
            }
            System.out.print("x = ");
            x = sc.nextDouble();
        } else {
            pairs.put(-4.0, 1127.0);
            pairs.put(-2.0, 81.0);
            pairs.put(0.0, 3.0);
            pairs.put(2.0, 77.0);
            pairs.put(4.0, 1023.0);
            x = 3;
        }
        System.out.println("----------------------------------------------------------------------------------");
        for (Double xn : pairs.keySet()) {
            System.out.printf("f(%.1f) = %.1f, ", xn, pairs.get(xn));
        }
        System.out.println("x = " + x);
        System.out.println("----------------------------------------------------------------------------------");
    }

    private static double lagrunge() {
        double w = 0.0;
        for (Double xn : pairs.keySet()) {
            double fracU = 1.0;
            double fracD = 1.0;
            for (Double xf : pairs.keySet()) {
                if (!xf.equals(xn)) {
                    fracU *= x - xf;
                    fracD *= xn - xf;
                }
            }
            w += pairs.get(xn) * (fracU / fracD);
        }
        return w;
    }

    private static double newtonIR() {
        int n = pairs.size();
        List<Double> xVals = new ArrayList<>(pairs.keySet());
        double[][] diffQuot = new double[n][n];

        for (int i = 0; i < n; i++) {
            diffQuot[i][0] = pairs.get(xVals.get(i));
        }
        for (int j = 1; j < n; j++) {
            for (int i = 0; i < n - j; i++) {
                diffQuot[i][j] = (diffQuot[i + 1][j - 1] - diffQuot[i][j - 1]) / (xVals.get(i + j) - xVals.get(i));
            }
        }

        double w = diffQuot[0][0];
        for (int i = 1; i < n; i++) {
            double product = 1;
            for (int j = 0; j < i; j++) {
                product *= (x - xVals.get(j));
            }
            w += diffQuot[0][i] * product;
        }
        return w;
    }

    private static void getResult(String method, double result) {
        System.out.printf("%s: W(x) = %.1f%n", method, result);
    }

    public static void main(String[] args) {
        input();
        getResult("LaGrunge", lagrunge());
        getResult("NewtonIR", newtonIR());
    }
}