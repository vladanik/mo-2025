package optym;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;

public class Optimum {

    public static final Boolean example1 = false; // (100 - x)^2
    public static final Boolean example2 = false; // 1/3x^3 + 1/2x^2 - 5x + 2
    public static final Boolean example3 = true; // 1/3x^3 + 2.05x^2 - 9x + 15

    public static final double epsilonZero = 0;
    public static final double[] epsilons = { 0.1, 0.01, 0.001, 0.00001, 0.000000001, 0.0000000000001 };
    public static final int maxIter = 100000;
    public static final int[] iterations = { 3, 5, 7, 10, 15, 30, 50 };

    public static final double X = example1 ? 100.0 : (example2 ? (-1.0 + Math.sqrt(21)) / 2.0 : (-41.0 - Math.sqrt(5281)) / 20.0);
    public static final Boolean FIND_MIN = example1 || example2;
    public static final Boolean INCLUDE_DETAILS = false;

    public static double a;
    public static double b;
    
    public static void setDefaultRegion() {
        a = example1 ? 60 : (example2 ? 1 : -10);
        b = example1 ? 150 : (example2 ? 2 : -3);
    }

    public static double f(double x) {
        if (example1) {
            return Math.pow(100 - x, 2);
        }
        if (example2) {
            return (1.0/3.0) * Math.pow(x, 3) + (1.0/2.0) * Math.pow(x, 2) - 5 * x + 2;
        }
        if (example3) {
            return (1.0/3.0) * Math.pow(x, 3) + 2.05 * Math.pow(x, 2) - 9 * x + 15;
        }
        return 0;
    }

    public static double[] fp(double x) {
        if (example1) {
            DerivativeStructure xDS = new DerivativeStructure(1, 2, 0, x);
            DerivativeStructure result = xDS.multiply(-1).add(100).pow(2);

            double value = result.getValue();
            double derivative1 = result.getPartialDerivative(1);
            double derivative2 = result.getPartialDerivative(2);
            return new double[]{value, derivative1, derivative2};
        }
        if (example2) {
            DerivativeStructure xDS = new DerivativeStructure(1, 3, 0, x);
            DerivativeStructure term1 = xDS.pow(3).multiply(1.0/3.0);
            DerivativeStructure term2 = xDS.pow(2).multiply(1.0/2.0);
            DerivativeStructure term3 = xDS.multiply(-5);
            DerivativeStructure result = term1.add(term2).add(term3).add(2);

            double value = result.getValue();
            double derivative1 = result.getPartialDerivative(1);
            double derivative2 = result.getPartialDerivative(2);
            double derivative3 = result.getPartialDerivative(3);
            return new double[]{value, derivative1, derivative2, derivative3};
        }
        if (example3) {
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
        return null;
    }

    public static class OptimizationResult {
        public double optimum;
        public double optimumValue;
        public int iterations;
        public String details;

        public OptimizationResult(double optimum, double optimumValue, int iterations) {
            this.optimum = optimum;
            this.optimumValue = optimumValue;
            this.iterations = iterations;
        }

        public OptimizationResult(double optimum, double optimumValue, int iterations, String details) {
            this.optimum = optimum;
            this.optimumValue = optimumValue;
            this.iterations = iterations;
            this.details = details;
        }

        @Override
        public String toString() {
            return "x: " + optimum + ", f(x): " + optimumValue + ", iterations: " + iterations + (INCLUDE_DETAILS ? ("\n" + details) : "");
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
        }
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
            if (epsilon == 0 && iter >= iterationLimit) {
                break;
            }
        }
        double x = (a + b) / 2.0;
        return new OptimizationResult(x, f(x), iter, details.toString());
    }

    public static OptimizationResult bisectionMethod(double epsilon, int iterationLimit) {
        if (fp(a)[1] * fp(b)[1] >= 0) {
            return new OptimizationResult(Double.NaN, Double.NaN, 0);
        }
        int iter = 0;
        double x_mid;
        StringBuilder details = new StringBuilder();
        details.append("iter: x, f(x), a, b\n");
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
            details.append(iter + ": " + x_mid + ", " + f(x_mid) + ", " + a + ", " + b + "\n");
        }
        return new OptimizationResult(x_mid, f(x_mid), iter, details.toString());
    }

    public static OptimizationResult newtonMethod(double epsilon, int iterationLimit) {
        if (fp(a)[2] * fp(b)[2] < 0 || fp(a)[3] * fp(b)[3] < 0) {
            return new OptimizationResult(Double.NaN, Double.NaN, 0);
        }
        double x0 = (fp(a)[3] * fp(a)[1] > 0) ? a : b;
        int iter = 0;
        StringBuilder details = new StringBuilder();
        details.append("iter: x0, f(x0)\n");
        details.append(iter + ": " + x0 + ", " + f(x0) + "\n");
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
        }
        return new OptimizationResult(x0, f(x0), iter, details.toString());
    }

    public static OptimizationResult secantMethod(double epsilon, int iterationLimit) {
        double x0 = fp(a)[1] > 0 ? b : a;
        int iter = 0;
        StringBuilder details = new StringBuilder();
        details.append("iter: x0, f(x0)\n");
        details.append(iter + ": " + x0 + ", " + f(x0) + "\n");
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
        }
        return new OptimizationResult(x0, f(x0), iter, details.toString());
    }

    public static void main(String[] args) {
        if (example1 || example3) {
            runDwudzielna();
            runFibonacciego();
        }
        if (example2 || example3) {
            runBisekcji();
            runStycznych();
            runSiecznych();
        }
    }

    private static void runDwudzielna() {
        System.out.println("\n======= Metoda Dwudzielna =======");
        setDefaultRegion();
        for (double e : epsilons) {
            System.out.println("\n=== Epsilon: " + e + " ===");

            OptimizationResult result = dwudzielnaSearch(e, maxIter);
            System.out.println(result);
            setDefaultRegion();
        }
        for (int i : iterations) {
            System.out.println("\n=== Max iterations: " + i + " ===");

            OptimizationResult result = dwudzielnaSearch(epsilonZero, i);
            System.out.println(result);
            setDefaultRegion();
        }
    }

    private static void runFibonacciego() {
        System.out.println("\n======= Metoda Fibonacciego =======");
        setDefaultRegion();
        for (double e : epsilons) {
            System.out.println("\n=== Epsilon: " + e + " ===");

            OptimizationResult result = fibonacciSearch(e, maxIter);
            System.out.println(result);
            setDefaultRegion();
        }
        for (int i : iterations) {
            System.out.println("\n=== Max iterations: " + i + " ===");

            OptimizationResult result = fibonacciSearch(epsilonZero, i);
            System.out.println(result);
            setDefaultRegion();
        }
    }

    private static void runBisekcji() {
        System.out.println("\n======= Metoda Bisekcji =======");
        setDefaultRegion();
        for (double e : epsilons) {
            System.out.println("\n=== Epsilon: " + e + " ===");

            OptimizationResult result = bisectionMethod(e, maxIter);
            System.out.println(result);
            setDefaultRegion();
        }
        for (int i : iterations) {
            System.out.println("\n=== Max iterations: " + i + " ===");

            OptimizationResult result = bisectionMethod(epsilonZero, i);
            System.out.println(result);
            setDefaultRegion();
        }
    }

    private static void runStycznych() {
        System.out.println("\n======= Metoda Stycznych =======");
        setDefaultRegion();
        for (double e : epsilons) {
            System.out.println("\n=== Epsilon: " + e + " ===");

            OptimizationResult result = newtonMethod(e, maxIter);
            System.out.println(result);
            setDefaultRegion();
        }
        for (int i : iterations) {
            System.out.println("\n=== Max iterations: " + i + " ===");

            OptimizationResult result = newtonMethod(epsilonZero, i);
            System.out.println(result);
            setDefaultRegion();
        }
    }

    private static void runSiecznych() {
        System.out.println("\n======= Metoda Siecznych =======");
        setDefaultRegion();
        for (double e : epsilons) {
            System.out.println("\n=== Epsilon: " + e + " ===");

            OptimizationResult result = secantMethod(e, maxIter);
            System.out.println(result);
            setDefaultRegion();
        }
        for (int i : iterations) {
            System.out.println("\n=== Max iterations: " + i + " ===");

            OptimizationResult result = secantMethod(epsilonZero, i);
            System.out.println(result);
            setDefaultRegion();
        }
    }
}
