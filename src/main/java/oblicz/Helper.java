package oblicz;

import java.util.List;

public class Helper {

    public static double[] gaussElimination(double[][] A, double[] b) {
        int n = b.length;
        for (int i = 0; i < n; i++) {
            int max = i;
            for (int j = i + 1; j < n; j++) {
                if (Math.abs(A[j][i]) > Math.abs(A[max][i])) {
                    max = j;
                }
            }
            double[] temp = A[i];
            A[i] = A[max];
            A[max] = temp;
            double t = b[i];
            b[i] = b[max];
            b[max] = t;

            for (int j = i + 1; j < n; j++) {
                double factor = A[j][i] / A[i][i];
                b[j] -= factor * b[i];
                for (int k = i; k < n; k++) {
                    A[j][k] -= factor * A[i][k];
                }
            }
        }

        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }

    public static boolean isValidH(List<Double> xList) {
        double h = xList.get(1) - xList.get(0);
        for (int i = 1; i < xList.size() - 1; i++) {
            if (Math.abs(xList.get(i + 1) - xList.get(i)) != h) {
                return false;
            }
        }
        return true;
    }

    public static int factorial(int n) {
        int result = 1;
        for (int i = 2; i <= n; i++) result *= i;
        return result;
    }
}
