import java.util.*;

public class Interpolation {

    static boolean userInput = false;

    static Map<Double, Double> pairs = new HashMap<>();
    static double x;
    static double W = 0;

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
    }

    private static void lagrunge() {
        for (Double xn : pairs.keySet()) {
            double fracU = 1.0;
            double fracD = 1.0;
            for (Double xf : pairs.keySet()) {
                if (!xf.equals(xn)) {
                    fracU *= x - xf;
                    fracD *= xn - xf;
                }
            }
            W += pairs.get(xn) * (fracU / fracD);
        }
    }

    private static void getResult() {
        System.out.println("----------------------------------------------------------------------------------");
        for (Double xn : pairs.keySet()) {
            System.out.print("f(" + xn + ") = " + pairs.get(xn) + ", ");
        }
        System.out.println("x = " + x);
        System.out.println();
        System.out.println("W(x) = " + W);
        System.out.println("----------------------------------------------------------------------------------");
    }

    public static void main(String[] args) {
        input();
        lagrunge();
        getResult();
    }
}