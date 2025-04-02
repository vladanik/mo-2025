package optym;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;

import java.util.ArrayList;
import java.util.List;

public class Optimum {
    public static final Boolean IS_TEST = true;

    public static final double epsilonZero = 0;
    public static final double[] epsilons_S = { 0.1, 0.01, 0.001, 0.00001, 0.000000001, 0.0000000000001 };
    public static final double[] epsilons_T = { 0.001 };
    public static final double[] epsilons = IS_TEST ? epsilons_T : epsilons_S;
    public static final int maxIter = 100000;
    public static final int[] iterations_S = { 3, 5, 7, 10, 15, 30, 50 };
    public static final int[] iterations_T = { };
    public static final int[] iterations = IS_TEST ? iterations_T : iterations_S;

    public static final double X = (-41.0 - Math.sqrt(5281)) / 20.0;
    public static final Boolean FIND_MIN = false;

    public static double a;
    public static double b;
    
    public static void setDefaultRegion() {
        a = -10;
        b = -3;
    }

    public static double f(double x) {
        return (1.0/3.0) * Math.pow(x, 3) + 2.05 * Math.pow(x, 2) - 9 * x + 15; // 1/3x^3 + 2.05x^2 - 9x + 15
    }

    public static double[] fp(double x) {
        DerivativeStructure xDS = new DerivativeStructure(1, 3, 0, x);
        DerivativeStructure term1 = xDS.pow(3).multiply(1.0/3.0);
        DerivativeStructure term2 = xDS.pow(2).multiply(2.05);
        DerivativeStructure term3 = xDS.multiply(-9);
        DerivativeStructure result = term1.add(term2).add(term3).add(15);
        return new double[]{
                result.getValue(),
                result.getPartialDerivative(1),
                result.getPartialDerivative(2),
                result.getPartialDerivative(3)
        };
    }

    public static class OptimizationResult {
        public double optimum;
        public double optimumValue;
        public double errorPercentage;
        public int iterations;
        public String details;

        public OptimizationResult(double optimum, double optimumValue, int iterations, String details) {
            this.optimum = optimum;
            this.optimumValue = optimumValue;
            this.iterations = iterations;
            this.details = details;
        }

        @Override
        public String toString() {
            return "x: " + optimum + ", f(x): " + optimumValue + ", iterations: " + iterations
                    + ", error: " + this.getErrorPercentage() + (IS_TEST && details != null ? ("\n" + details) : "");
        }

        private String getErrorPercentage() {
            double errorP = Math.abs((X - this.optimum) / X) * 100.0;
            return String.format("%.6f%%", errorP);
        }
    }

    public static long fib(int n) {
        if(n <= 1) return 1;
        long f0 = 1, f1 = 1, fn = 1;
        for (int i = 2; i <= n; i++) {
            fn = f0 + f1;
            f0 = f1;
            f1 = fn;
        }
        return fn;
    }

    public static long[] computeFibSequence(int n) {
        long[] fibs = new long[n + 1];
        fibs[0] = 1;
        if(n >= 1) fibs[1] = 1;
        for (int i = 2; i <= n; i++) {
            fibs[i] = fibs[i - 1] + fibs[i - 2];
        }
        return fibs;
    }

    public static OptimizationResult dwudzielnaSearch(double epsilon, int iterationLimit) {
        double x_s = (a + b) / 2.0;
        double L = b - a;
        int iteration = 0;
        StringBuilder details = new StringBuilder();
        details.append("iter: x, f(x), a, b, L\n");
        details.append(iteration + ": " + x_s + ", " + f(x_s) + ", " + a + ", " + b + ", " + L + "\n");
        List<double[]> points = new ArrayList<>();

        while (L > epsilon && iteration < iterationLimit) {
            iteration++;
            L = b - a;
            double x1 = a + L / 4.0;
            double x2 = b - L / 4.0;
            double f_xs = f(x_s);
            double f_x1 = f(x1);
            double f_x2 = f(x2);

            if (FIND_MIN) {
                if (f_x1 < f_xs) {
                    b = x_s;
                    x_s = x1;
                } else if (f_x2 < f_xs) {
                    a = x_s;
                    x_s = x2;
                } else {
                    a = x1;
                    b = x2;
                }
            } else {
                if (f_x1 > f_xs) {
                    b = x_s;
                    x_s = x1;
                } else if (f_x2 > f_xs) {
                    a = x_s;
                    x_s = x2;
                } else {
                    a = x1;
                    b = x2;
                }
            }
            details.append(iteration + ": " + x_s + ", " + f(x_s) + ", " + a + ", " + b + ", " + L + "\n");
            points.add(new double[]{x_s, f(x_s)});
        }
        showChart("Dwudzielna", points, x_s);
        return new OptimizationResult(x_s, f(x_s), iteration, details.toString());
    }

    public static OptimizationResult fibonacciSearch(double epsilon, int iterationLimit) {
        double L = b - a;
        int n;
        if (epsilon > 0) {
            n = 1;
            while ( (L / (double)fib(n)) >= 2 * epsilon ) {
                n++;
            }
            n--;
        } else {
            n = iterationLimit;
        }
        long[] fibSeq = computeFibSequence(n);

        double x1 = b - ((double)fibSeq[n - 1] / fibSeq[n]) * (b - a);
        double x2 = a + ((double)fibSeq[n - 1] / fibSeq[n]) * (b - a);

        int iter = 0;
        StringBuilder details = new StringBuilder();
        details.append("iter: x1, f(x1), x2, f(x2)\n");
        details.append(iter + ": " + x1 + ", " + f(x1) + ", " + x2 + ", " + f(x2) + "\n");
        List<double[]> points = new ArrayList<>();
        points.add(new double[]{x1, f(x1)});
        points.add(new double[]{x2, f(x2)});

        while (n > 1) {
            iter++;
            if (epsilon > 0 && Math.abs(x2 - x1) < epsilon) {
                break;
            }
            if (FIND_MIN) {
                if (f(x1) < f(x2)) {
                    b = x2;
                    x2 = x1;
                    n = n - 1;
                    x1 = b - ((double)fibSeq[n - 1] / fibSeq[n]) * (b - a);
                } else {
                    a = x1;
                    x1 = x2;
                    n = n - 1;
                    x2 = a + ((double)fibSeq[n - 1] / fibSeq[n]) * (b - a);
                }
            } else {
                if (f(x1) > f(x2)) {
                    b = x2;
                    x2 = x1;
                    n = n - 1;
                    x1 = b - ((double)fibSeq[n - 1] / fibSeq[n]) * (b - a);
                } else {
                    a = x1;
                    x1 = x2;
                    n = n - 1;
                    x2 = a + ((double)fibSeq[n - 1] / fibSeq[n]) * (b - a);
                }
            }
            details.append(iter + ": " + x1 + ", " + f(x1) + ", " + x2 + ", " + f(x2) + "\n");
            points.add(new double[]{x1, f(x1)});
            points.add(new double[]{x2, f(x2)});
            if (epsilon == 0 && iter >= iterationLimit) {
                break;
            }
        }
        double x = (a + b) / 2.0;
        showChart("Fibonacciego", points, x);
        return new OptimizationResult(x, f(x), iter, details.toString());
    }

    public static OptimizationResult bisectionSearch(double epsilon, int iterationLimit) {
        if (fp(a)[1] * fp(b)[1] >= 0) {
            return new OptimizationResult(Double.NaN, Double.NaN, 0, null);
        }
        int iter = 0;
        double x_mid;
        StringBuilder details = new StringBuilder();
        details.append("iter: x, f(x), a, b\n");
        List<double[]> points = new ArrayList<>();

        while (true) {
            iter++;
            x_mid = (a + b) / 2.0;
            double fpm = fp(x_mid)[1];
            if (epsilon > 0) {
                if (Math.abs(fpm) < epsilon) {
                    break;
                }
            } else {
                if (iter >= iterationLimit || Math.abs(fpm) <= epsilon) {
                    break;
                }
            }
            if (fp(a)[1] * fpm < 0) {
                b = x_mid;
            } else {
                a = x_mid;
            }
            details.append(iter + ": " + x_mid + ", " + f(x_mid) + ", " + a + ", " + b + "\n");
            points.add(new double[]{x_mid, f(x_mid)});
        }
        showChart("Bisekcji", points, x_mid);
        return new OptimizationResult(x_mid, f(x_mid), iter, details.toString());
    }

    public static OptimizationResult newtonSearch(double epsilon, int iterationLimit) {
        if (fp(a)[2] * fp(b)[2] < 0 || fp(a)[3] * fp(b)[3] < 0) {
            return new OptimizationResult(Double.NaN, Double.NaN, 0, null);
        }
        double x0 = (fp(a)[3] * fp(a)[1] > 0) ? a : b;
        int iter = 0;
        StringBuilder details = new StringBuilder();
        details.append("iter: x0, f(x0)\n");
        details.append(iter + ": " + x0 + ", " + f(x0) + "\n");
        List<double[]> points = new ArrayList<>();
        points.add(new double[]{x0, f(x0)});

        while (true) {
            iter++;
            double x1 = x0 - fp(x0)[1] / fp(x0)[2];
            if (epsilon > 0) {
                if (Math.abs(fp(x1)[1]) < epsilon || Math.abs(x1 - x0) < epsilon) {
                    x0 = x1;
                    break;
                }
            } else {
                if (iter >= iterationLimit) {
                    x0 = x1;
                    break;
                }
            }
            x0 = x1;
            details.append(iter + ": " + x0 + ", " + f(x0) + "\n");
            points.add(new double[]{x0, f(x0)});
        }
        showChart("Stycznych (Newtona)", points, x0);
        return new OptimizationResult(x0, f(x0), iter, details.toString());
    }

    public static OptimizationResult secantSearch(double epsilon, int iterationLimit) {
        double x0 = fp(a)[1] > 0 ? b : a;
        int iter = 0;
        StringBuilder details = new StringBuilder();
        details.append("iter: x0, f(x0)\n");
        details.append(iter + ": " + x0 + ", " + f(x0) + "\n");
        List<double[]> points = new ArrayList<>();
        points.add(new double[]{x0, f(x0)});

        while (true) {
            iter++;
            double x1;
            if (fp(a)[1] > 0) {
                x1 = x0 - fp(x0)[1] / (fp(x0)[1] - fp(a)[1]) * (x0 - a);
            } else {
                x1 = x0 - fp(x0)[1] / (fp(b)[1] - fp(x0)[1]) * (b - x0);
            }
            if (epsilon > 0) {
                if (Math.abs(fp(x1)[1]) < epsilon || Math.abs(x1 - x0) < epsilon) {
                    x0 = x1;
                    break;
                }
            } else {
                if (iter >= iterationLimit) {
                    x0 = x1;
                    break;
                }
            }
            x0 = x1;
            details.append(iter + ": " + x0 + ", " + f(x0) + "\n");
            points.add(new double[]{x0, f(x0)});
        }
        showChart("Siecznych", points, x0);
        return new OptimizationResult(x0, f(x0), iter, details.toString());
    }

    public static OptimizationResult goldenSearch(double epsilon, int iterationLimit) {
        final double k = (Math.sqrt(5) - 1) / 2;
        double x1 = b - k * (b - a);
        double x2 = a + k * (b - a);
        int iter = 0;
        StringBuilder details = new StringBuilder();
        details.append("iter: x1, f(x1), x2, f(x2)\n");
        details.append(iter + ": " + x1 + ", " + f(x1) + ", " + x2 + ", " + f(x2) + "\n");
        List<double[]> points = new ArrayList<>();
        points.add(new double[]{x1, f(x1)});
        points.add(new double[]{x2, f(x2)});

        while (Math.abs(x2 - x1) >= epsilon && iter < iterationLimit) {
            iter++;
            if (FIND_MIN) {
                if (f(x1) < f(x2)) {
                    b = x2;
                    x2 = x1;
                    x1 = b - k * (b - a);
                } else {
                    a = x1;
                    x1 = x2;
                    x2 = a + k * (b - a);
                }
            } else {
                if (f(x1) > f(x2)) {
                    b = x2;
                    x2 = x1;
                    x1 = b - k * (b - a);
                } else {
                    a = x1;
                    x1 = x2;
                    x2 = a + k * (b - a);
                }
            }
            details.append(iter + ": " + x1 + ", " + f(x1) + ", " + x2 + ", " + f(x2) + "\n");
            points.add(new double[]{x1, f(x1)});
            points.add(new double[]{x2, f(x2)});
        }
        double x = (a + b) / 2.0;
        showChart("Zlotego Podzialu", points, x);
        return new OptimizationResult(x, f(x), iter, details.toString());
    }

    public static void main(String[] args) {
        run("Dwudzielna");
        run("Fibonacciego");
        run("Bisekcji");
        run("Stycznych");
        run("Siecznych");
        run("Zlotego Podzialu");
    }

    private static void run(String name) {
        System.out.println("\n======= Metoda " + name + " =======");
        setDefaultRegion();
        for (double e : epsilons) {
            System.out.println("\n=== Epsilon: " + e + " ===");
            OptimizationResult result;
            switch (name) {
                case "Dwudzielna":
                    result = dwudzielnaSearch(e, maxIter);
                    break;
                case "Fibonacciego":
                    result = fibonacciSearch(e, maxIter);
                    break;
                case "Bisekcji":
                    result = bisectionSearch(e, maxIter);
                    break;
                case "Stycznych":
                    result = newtonSearch(e, maxIter);
                    break;
                case "Siecznych":
                    result = secantSearch(e, maxIter);
                    break;
                case "Zlotego Podzialu":
                    result = goldenSearch(e, maxIter);
                    break;
                default:
                    throw new IllegalArgumentException("Unknown method: " + name);
            }
            System.out.println(result);
            setDefaultRegion();
        }
        for (int i : iterations) {
            System.out.println("\n=== Max iterations: " + i + " ===");
            OptimizationResult result;
            switch (name) {
                case "Dwudzielna":
                    result = dwudzielnaSearch(epsilonZero, i);
                    break;
                case "Fibonacciego":
                    result = fibonacciSearch(epsilonZero, i);
                    break;
                case "Bisekcji":
                    result = bisectionSearch(epsilonZero, i);
                    break;
                case "Stycznych":
                    result = newtonSearch(epsilonZero, i);
                    break;
                case "Siecznych":
                    result = secantSearch(epsilonZero, i);
                    break;
                case "Zlotego Podzialu":
                    result = goldenSearch(epsilonZero, i);
                    break;
                default:
                    throw new IllegalArgumentException("Unknown method: " + name);
            }
            System.out.println(result);
            setDefaultRegion();
        }
    }

    private static List<double[]> getFunctionPoints() {
        List<double[]> points = new ArrayList<>();
        setDefaultRegion();
        for (double i = a - 1; i <= b + 1; i += 0.1) {
            points.add(new double[]{i, f(i)});
        }
        return points;
    }

    private static void showChart(String name, List<double[]> points, double x) {
        if (IS_TEST) {
            Chart chart = new Chart(name, getFunctionPoints(), points, new double[]{x, f(x)});
            chart.show();
        }
    }
}
