package optym;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;

public class Optimum {

    public static double epsilonZero = 0;
    public static double[] epsilons = { 1, 0.1, 0.001, 0.00001, 0.000000001, 0.0000000000001, 0.000000000000001 };
    public static int maxIter = 100000;
    public static int[] iterations = { 3, 4, 5, 7, 10, 15, 20, 30, 40, 50, 70, 100 };

    public static double a = -10, b = -3;
    public static double f(double x) {
        return (1.0/3.0) * Math.pow(x, 3) + 2.05 * Math.pow(x, 2) - 9 * x + 15;
    }

    public static double[] fp(double x) {
        DerivativeStructure xDS = new DerivativeStructure(1, 3, 0, x);
        DerivativeStructure term1 = xDS.pow(3).multiply(1.0/3.0);
        DerivativeStructure term2 = xDS.pow(2).multiply(2.05);
        DerivativeStructure term3 = xDS.multiply(-9);
        DerivativeStructure result = term1.add(term2).add(term3).add(15);
        double value = result.getValue();
        double derivative1 = result.getPartialDerivative(1);
        double derivative2 = result.getPartialDerivative(2);
        double derivative3 = result.getPartialDerivative(3);
        return new double[]{value, derivative1, derivative2, derivative3};
    }

    public static class OptimizationResult {
        public double optimum;
        public double optimumValue;
        public int iterations;
        public double a, b;

        public OptimizationResult(double optimum, double optimumValue, int iterations, double a, double b) {
            this.optimum = optimum;
            this.optimumValue = optimumValue;
            this.iterations = iterations;
            this.a = a;
            this.b = b;
        }

        @Override
        public String toString() {
            return "x: " + optimum + ", f(x): " + optimumValue + ", iterations: " + iterations;
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

    public static OptimizationResult dwudzielnaSearch(double a, double b, double epsilon, int iterationLimit, boolean findMin) {
        double x_s = (a + b) / 2.0;
        double L = b - a;
        int iteration = 0;

        while (L > epsilon && iteration < iterationLimit) {
            iteration++;
            L = b - a;
            double x1 = a + L / 4.0;
            double x2 = b - L / 4.0;
            double f_xs = f(x_s);
            double f_x1 = f(x1);
            double f_x2 = f(x2);

            if (findMin) {
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
        }
        return new OptimizationResult(x_s, f(x_s), iteration, a, b);
    }

    public static OptimizationResult fibonacciSearch(double a, double b, double epsilon, int iterationLimit, boolean findMin) {
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
        while (n > 1) {
            iter++;
            if (epsilon > 0 && Math.abs(x2 - x1) < epsilon) {
                break;
            }
            if (findMin) {
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
            if (epsilon == 0 && iter >= iterationLimit) {
                break;
            }
        }
        double x = (a + b) / 2.0;
        return new OptimizationResult(x, f(x), iter, a, b);
    }

    public static OptimizationResult bisectionMethod(double a, double b, double epsilon, int iterationLimit) {
        if (fp(a)[1] * fp(b)[1] >= 0) {
            return new OptimizationResult(Double.NaN, Double.NaN, 0, a, b);
        }
        int iter = 0;
        double x_mid;
        while (true) {
            iter++;
            x_mid = (a + b) / 2.0;
            double fpm = fp(x_mid)[1];
            if (epsilon > 0) {
                if (Math.abs(fpm) < epsilon) {
                    break;
                }
            } else {
                if (iter >= iterationLimit) {
                    break;
                }
            }
            if (fp(a)[1] * fpm < 0) {
                b = x_mid;
            } else {
                a = x_mid;
            }
        }
        return new OptimizationResult(x_mid, f(x_mid), iter, a, b);
    }

    public static OptimizationResult newtonMethod(double a, double b, double epsilon, int iterationLimit) {
        double x0 = (fp(a)[3] * fp(a)[1] > 0) ? a : b;
        int iter = 0;
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
        }
        return new OptimizationResult(x0, f(x0), iter, a, b);
    }

    public static OptimizationResult secantMethod(double a, double b, double epsilon, int iterationLimit) {
        double x0 = fp(a)[1] > 0 ? b : a;
        int iter = 0;
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
        }
        return new OptimizationResult(x0, f(x0), iter, a, b);
    }

    public static void main(String[] args) {
        runDwudzielna();
        runFibonacciego();
        runBisekcji();
        runStycznych();
        runSiecznych();
    }

    private static void runDwudzielna() {
        System.out.println("\n======= Metoda Dwudzielna =======");
        for (double e : epsilons) {
            System.out.println("\n=== Epsilon: " + e + " ===");

            System.out.print("Maximum: ");
            OptimizationResult resultMax = dwudzielnaSearch(a, b, e, maxIter, false);
            System.out.println(resultMax);

            System.out.print("Minimum: ");
            OptimizationResult resultMin = dwudzielnaSearch(a, b, e, maxIter, true);
            System.out.println(resultMin);
        }
        for (int i : iterations) {
            System.out.println("\n=== Max iterations: " + i + " ===");

            System.out.print("Maximum: ");
            OptimizationResult resultMax = dwudzielnaSearch(a, b, epsilonZero, i, false);
            System.out.println(resultMax);

            System.out.print("Minimum: ");
            OptimizationResult resultMin = dwudzielnaSearch(a, b, epsilonZero, i, true);
            System.out.println(resultMin);
        }
    }

    private static void runFibonacciego() {
        System.out.println("\n======= Metoda Fibonacciego =======");
        for (double e : epsilons) {
            System.out.println("\n=== Epsilon: " + e + " ===");

            System.out.print("Maximum: ");
            OptimizationResult resultMax = fibonacciSearch(a, b, e, maxIter, false);
            System.out.println(resultMax);

            System.out.print("Minimum: ");
            OptimizationResult resultMin = fibonacciSearch(a, b, e, maxIter, true);
            System.out.println(resultMin);
        }
        for (int i : iterations) {
            System.out.println("\n=== Max iterations: " + i + " ===");

            System.out.print("Maximum: ");
            OptimizationResult resultMax = fibonacciSearch(a, b, epsilonZero, i, false);
            System.out.println(resultMax);

            System.out.print("Minimum: ");
            OptimizationResult resultMin = fibonacciSearch(a, b, epsilonZero, i, true);
            System.out.println(resultMin);
        }
    }

    private static void runBisekcji() {
        System.out.println("\n======= Metoda Bisekcji =======");
        for (double e : epsilons) {
            System.out.println("\n=== Epsilon: " + e + " ===");

            OptimizationResult result = bisectionMethod(a, b, e, maxIter);
            System.out.println(result);
        }
        for (int i : iterations) {
            System.out.println("\n=== Max iterations: " + i + " ===");

            OptimizationResult result = bisectionMethod(a, b, epsilonZero, i);
            System.out.println(result);
        }
    }

    private static void runStycznych() {
        System.out.println("\n======= Metoda Stycznych =======");
        for (double e : epsilons) {
            System.out.println("\n=== Epsilon: " + e + " ===");

            OptimizationResult result = newtonMethod(a, b, e, maxIter);
            System.out.println(result);
        }
        for (int i : iterations) {
            System.out.println("\n=== Max iterations: " + i + " ===");

            OptimizationResult result = newtonMethod(a, b, epsilonZero, i);
            System.out.println(result);
        }
    }

    private static void runSiecznych() {
        System.out.println("\n======= Metoda Siecznych =======");
        for (double e : epsilons) {
            System.out.println("\n=== Epsilon: " + e + " ===");

            OptimizationResult result = secantMethod(a, b, e, maxIter);
            System.out.println(result);
        }
        for (int i : iterations) {
            System.out.println("\n=== Max iterations: " + i + " ===");

            OptimizationResult result = secantMethod(a, b, epsilonZero, i);
            System.out.println(result);
        }
    }
}
