package oblicz;

import java.util.*;

public class Interpolation {

    static boolean userInput = false;

    static Map<Double, Double> pairs = new HashMap<>();
    static List<Double> dfs = new ArrayList<>();
    static double x;

    static List<Double> xList;
    static int n;

    private static void input() {
        if (userInput) {
            Scanner sc = new Scanner(System.in);

            for (int i = 0; ; i++) {
                System.out.printf("x[%d] = ", i);
                double x = sc.nextDouble();
                System.out.printf("f(%.1f) = ", x);
                double y = sc.nextDouble();
                pairs.put(x, y);

                if (i == 0) {
                    System.out.printf("f'(%.1f) = ", x);
                    dfs.add(sc.nextDouble());
                }

                System.out.print("Continue? Y/n: ");
                String s = sc.next();
                if (!s.toLowerCase().startsWith("y")) {
                    System.out.printf("f'(%.1f) = ", x);
                    dfs.add(sc.nextDouble());
                    break;
                }
            }
            System.out.print("x = ");
            x = sc.nextDouble();
        } else {
            pairs.put(-4.0, 1127.0);
            pairs.put(-2.0, 81.0);
            pairs.put(0.0, 3.0);
            pairs.put(2.0, 77.0);
            pairs.put(4.0, 1023.0);
            dfs.add(-1093.0);
            dfs.add(1003.0);
            x = 3;
        }
        xList = new ArrayList<>(pairs.keySet());
        Collections.sort(xList);
        n = pairs.size();

        System.out.println("-------------------------------------------------------------------------------------------------------------------------------------------");
        for (Double xn : xList) {
            System.out.printf("f(%.1f) = %.1f, ", xn, pairs.get(xn));
        }
        System.out.printf("f'(%.1f) = %.1f, f'(%.1f) = %.1f, ", xList.get(0), dfs.get(0), xList.get(xList.size() - 1), dfs.get(1));
        System.out.println("x = " + x);
        System.out.println("-------------------------------------------------------------------------------------------------------------------------------------------");
    }

    private static double lagrunge() {
        double w = 0.0;
        for (Double xn : xList) {
            double fracU = 1.0;
            double fracD = 1.0;
            for (Double xf : xList) {
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
        double[][] diffQuot = new double[n][n];

        for (int i = 0; i < n; i++) {
            diffQuot[i][0] = pairs.get(xList.get(i));
        }
        for (int j = 1; j < n; j++) {
            for (int i = 0; i < n - j; i++) {
                diffQuot[i][j] = (diffQuot[i + 1][j - 1] - diffQuot[i][j - 1]) / (xList.get(i + j) - xList.get(i));
            }
        }

        double w = diffQuot[0][0];
        for (int i = 1; i < n; i++) {
            double product = 1;
            for (int j = 0; j < i; j++) {
                product *= (x - xList.get(j));
            }
            w += diffQuot[0][i] * product;
        }
        return w;
    }

    private static double newtonRP() {
        if (!Helper.isValidH(xList)) {
            System.out.println("Punkty nie są równoodległe. Metoda Newtona z różnicami progresywnymi nie może być użyta.");
            return Double.NaN;
        }

        double h = xList.get(1) - xList.get(0);
        double[][] diff = new double[n][n];
        for (int i = 0; i < n; i++) {
            diff[i][0] = pairs.get(xList.get(i));
        }
        for (int j = 1; j < n; j++) {
            for (int i = 0; i < n - j; i++) {
                diff[i][j] = diff[i + 1][j - 1] - diff[i][j - 1];
            }
        }

        double result = diff[0][0];
        double product = 1;
        for (int k = 1; k < n; k++) {
            product *= (x - xList.get(k - 1));
            result += (diff[0][k] / (Helper.factorial(k) * Math.pow(h, k))) * product;
        }
        return result;
    }

    private static double[] splineCoeffs;
    private static double spline() {
        // 4 (a_0...a_3) + (n-2) (alpha_1...alpha_n-2) = n + 2
        double[][] A = new double[n + 2][n + 2];
        double[] b = new double[n + 2];

        // x_0: S3(x_0) = W3(x_0)
        A[0][0] = 1;
        A[0][1] = xList.get(0);
        A[0][2] = Math.pow(xList.get(0), 2);
        A[0][3] = Math.pow(xList.get(0), 3);
        b[0] = pairs.get(xList.get(0));

        // x_1: S3(x_1) = W3(x_1)
        A[1][0] = 1;
        A[1][1] = xList.get(1);
        A[1][2] = Math.pow(xList.get(1), 2);
        A[1][3] = Math.pow(xList.get(1), 3);
        b[1] = pairs.get(xList.get(1));

        // S3(x_i) = W3(x_i) + sum{j=1...i–1} α_j*(x_i – x_j)^3
        for (int i = 2; i < n; i++) {
            A[i][0] = 1;
            A[i][1] = xList.get(i);
            A[i][2] = Math.pow(xList.get(i), 2);
            A[i][3] = Math.pow(xList.get(i), 3);
            for (int j = 1; j <= Math.min(i - 1, n - 2); j++) {
                A[i][3 + j] = Math.pow(xList.get(i) - xList.get(j), 3);
            }
            b[i] = pairs.get(xList.get(i));
        }

        // x_0: S3'(x_0) = W3'(x_0)
        A[n][1] = 1;
        A[n][2] = 2 * xList.get(0);
        A[n][3] = 3 * Math.pow(xList.get(0), 2);
        b[n] = dfs.get(0);

        // x_n-1: S3'(x_n-1) = W3'(x_n-1) + sum{j=1...n–2} 3α_j*(x_n-1 – x_j)^2
        A[n + 1][1] = 1;
        A[n + 1][2] = 2 * xList.get(n - 1);
        A[n + 1][3] = 3 * Math.pow(xList.get(n - 1), 2);
        for (int j = 1; j <= n - 2; j++) {
            A[n + 1][3 + j] = 3 * Math.pow(xList.get(n - 1) - xList.get(j), 2);
        }
        b[n + 1] = dfs.get(1);

        double[] coeffs = Helper.gaussElimination(A, b);
        splineCoeffs = coeffs.clone();

        double w3 = coeffs[0] + coeffs[1] * x + coeffs[2] * x * x + coeffs[3] * x * x * x;
        int interval = 0;
        for (int i = 0; i < n - 1; i++) {
            if (x >= xList.get(i) && x < xList.get(i + 1)) {
                interval = i;
                break;
            }
            if (i == n - 2 && x >= xList.get(i + 1)) {
                interval = n - 1;
            }
        }
        if (interval >= 1) {
            for (int j = 1; j <= interval && j <= (n - 2); j++) {
                w3 += coeffs[3 + j] * Math.pow(x - xList.get(j), 3);
            }
        }
        return w3;
    }

    private static void getResult(String method, double result) {
        System.out.printf("%s:\tW(x) = %.1f%n", method, result);
        if (method.equals("Sklejane")) {
            System.out.printf("a0 = %.6f%n", splineCoeffs[0]);
            System.out.printf("a1 = %.6f%n", splineCoeffs[1]);
            System.out.printf("a2 = %.6f%n", splineCoeffs[2]);
            System.out.printf("a3 = %.6f%n", splineCoeffs[3]);
            for (int j = 1; j <= pairs.size() - 2; j++) {
                System.out.printf("alpha%d = %.6f%n", j, splineCoeffs[3 + j]);
            }
        }
    }

    public static void main(String[] args) {
        input();
        getResult("LaGrunge", lagrunge());
        getResult("Newton IR", newtonIR());
        getResult("Newton RP", newtonRP());
        getResult("Sklejane", spline());
    }
}