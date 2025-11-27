import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.lang.Math;

public class Ex1 {
    /** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
    public static final double EPS = 0.001;
    /** The zero polynomial function is represented as an array with a single (0) entry. */
    public static final double[] ZERO = {0};
    public enum AreaMode { SIGNED, SIGNED_ABS, INTEGRAL_ABS }


    public static double f(double[] poly, double x) {
        double ans = 0;
        for(int i=0;i<poly.length;i++) {
            double c = Math.pow(x, i);
            ans += c*poly[i];
        }
        return ans;
    }

    public static double root_rec(double[] p, double x1, double x2, double eps) {
        double f1 = f(p,x1);

        if (x2 - x1 < eps) {
            return (x1 + x2) / 2.0;
        }

        double x12 = (x1+x2)/2;
        double f12 = f(p,x12);

        if (Math.abs(f12) < eps) {
            return x12;
        }

        if (f1 * f12 <= 0) {
            return root_rec(p, x1, x12, eps);
        } else {
            return root_rec(p, x12, x2, eps);
        }
    }

    public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        double [] ans = null;
        int lx = xx.length;
        int ly = yy.length;

        if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
            if (lx == 2) {
                // קו ישר: y = ax + b
                double x1 = xx[0], y1 = yy[0];
                double x2 = xx[1], y2 = yy[1];
                if (Math.abs(x1 - x2) < EPS) return null;

                double a = (y2 - y1) / (x2 - x1);
                double b = y1 - a * x1;
                ans = new double[] {b, a};

            } else if (lx == 3) {
                // פרבולה: y = ax^2 + bx + c
                double x1 = xx[0], y1 = yy[0];
                double x2 = xx[1], y2 = yy[1];
                double x3 = xx[2], y3 = yy[2];

                double x1_sq = x1 * x1;
                double x2_sq = x2 * x2;
                double x3_sq = x3 * x3;

                // הנוסחאות מפתרון מטריצה 3x3 - שיטת דטרמיננטות
                double det = x1_sq*(x2 - x3) + x2_sq*(x3 - x1) + x3_sq*(x1 - x2);

                if (Math.abs(det) < EPS) return null;

                double a = (y1*(x2 - x3) + y2*(x3 - x1) + y3*(x1 - x2)) / det;
                double b = (y1*(x3_sq - x2_sq) + y2*(x1_sq - x3_sq) + y3*(x2_sq - x1_sq)) / det;
                double c = (y1*x2*x3*(x3 - x2) + y2*x1*x3*(x2 - x1) + y3*x1*x2*(x1 - x3)) / det;

                ans = new double[] {c, b, a};
            }
        }
        return ans;
    }

    public static boolean equals(double[] p1, double[] p2) {
        int n = Math.max(p1.length, p2.length) - 1;
        for (int i = 0; i <= n + 1; i++) {
            double x = i;
            if (Math.abs(f(p1, x) - f(p2, x)) > EPS) return false;
        }
        return true;
    }

    /**
     * Computes a String representing the polynomial function.
     */
    public static String poly(double[] poly) {
        if (poly == null || poly.length == 0) return "0";

        StringBuilder sb = new StringBuilder();
        boolean first = true;

        for (int i = poly.length - 1; i >= 0; i--) {
            double coef = poly[i];
            if (Math.abs(coef) < 1e-12) continue; // דלג על 0

            // סימן
            if (!first) {
                sb.append(coef >= 0 ? " +" : " -");
            } else {
                if (coef < 0) sb.append("-");
                first = false;
            }

            double abs = Math.abs(coef);

            // דרגה 0
            if (i == 0) {
                sb.append(String.format("%.1f", abs));
            }
            // דרגה 1
            else if (i == 1) {
                if (Math.abs(abs - 1) > 1e-12)
                    sb.append(String.format("%.1f", abs));
                sb.append("x");
            }
            // דרגה ≥ 2
            else {
                if (Math.abs(abs - 1) > 1e-12)
                    sb.append(String.format("%.1f", abs));
                sb.append("x^").append(i);
            }
        }

        return first ? "0" : sb.toString();
    }
    // פונקציה להערכת פולינום בנקודה x (שווה ל-f)
    public static double evaluate(double[] p, double x) {
        return f(p, x);
    }
    private static double root_rec2(double[] p1, double[] p2,
                                   double left, double right, double eps) {

        double fLeft = f(p1, left) - f(p2, left);
        double fRight = f(p1, right) - f(p2, right);

        // אם אין שינוי סימן, לא נכנסים לפה
        if (Math.abs(right - left) < eps)
            return (left + right) / 2.0;

        double mid = (left + right) / 2.0;
        double fMid = f(p1, mid) - f(p2, mid);

        if (fLeft * fMid <= 0) {
            return root_rec2(p1, p2, left, mid, eps);
        } else {
            return root_rec2(p1, p2, mid, right, eps);
        }
    }

    public static double sameValue(double[] p1, double[] p2,
                                   double x1, double x2, double eps) {

        double prevX = x1;
        double prevVal = f(p1, prevX) - f(p2, prevX);

        int samples = 500000;         // דגימה מאוד צפופה
        double step = (x2 - x1) / samples;

        for (int i = 1; i <= samples; i++) {
            double x = x1 + i * step;
            double currVal = f(p1, x) - f(p2, x);

            if (prevVal * currVal <= 0) {
                // מצאנו קטע עם שינוי סימן -> עושים חיפוש בינארי על הקטע הזה
                return rootBinary(p1, p2, prevX, x, eps);
            }

            prevX = x;
            prevVal = currVal;
        }

        return Double.NaN;
    }



    public static double length(double[] p, double x1, double x2, int numberOfSegments) {
        if (numberOfSegments <= 0) throw new IllegalArgumentException("numberOfSegments must be positive");

        double dx = (x2 - x1) / numberOfSegments;
        double length = 0.0;

        double prevX = x1;
        double prevY = f(p, prevX);

        for (int i = 1; i <= numberOfSegments; i++) {
            double currX = x1 + i * dx;
            double currY = f(p, currX);

            double segment = Math.sqrt(Math.pow(dx, 2) + Math.pow(currY - prevY, 2));
            length += segment;

            prevX = currX;
            prevY = currY;
        }

        return length;
    }

    public static double area(double[] p1, double[] p2, double x1, double x2, int n) {
        if (n <= 0) throw new IllegalArgumentException("n must be positive");
        if (x1 > x2) { double t = x1; x1 = x2; x2 = t; }

        double dx = (x2 - x1) / n;
        double totalArea = 0.0;

        for (int i = 0; i < n; i++) {
            double left = x1 + i * dx;
            double right = left + dx;

            double fLeft = f(p1, left) - f(p2, left);
            double fRight = f(p1, right) - f(p2, right);

            if (fLeft * fRight < 0) {
                // יש חיתוך בתוך הטרפז
                double root = rootBinary(p1, p2, left, right, EPS);
                // חישוב שטח בנפרד לכל חלק
                double area1 = (Math.abs(f(p1, left) - f(p2, left)) + Math.abs(f(p1, root) - f(p2, root))) / 2.0 * (root - left);
                double area2 = (Math.abs(f(p1, root) - f(p2, root)) + Math.abs(f(p1, right) - f(p2, right))) / 2.0 * (right - root);
                totalArea += area1 + area2;
            } else {
                // טרפז רגיל
                totalArea += (Math.abs(fLeft) + Math.abs(fRight)) / 2.0 * dx;
            }
        }

        return totalArea;
    }


    // אינטגרציה על קטע קטן עם ערך מוחלט
    private static double areaSubsegment(double[] p1, double[] p2, double x1, double x2, int segments) {
        double dx = (x2 - x1) / segments;
        double area = 0.0;
        double prevX = x1;
        double prevY = Math.abs(f(p1, prevX) - f(p2, prevX));

        for (int i = 1; i <= segments; i++) {
            double currX = x1 + i * dx;
            double currY = Math.abs(f(p1, currX) - f(p2, currX));
            area += (prevY + currY) / 2.0 * dx;
            prevY = currY;
        }

        return area;
    }



    // rootBinary מדויק
    private static double rootBinary(double[] p1, double[] p2, double left, double right, double eps) {
        while (right - left > eps) {
            double mid = (left + right) / 2.0;
            double fLeft = f(p1, left) - f(p2, left);
            double fMid = f(p1, mid) - f(p2, mid);

            if (fLeft * fMid <= 0) right = mid;
            else left = mid;
        }
        return (left + right) / 2.0;
    }


    private static double simpson(java.util.function.DoubleUnaryOperator f,
                                  double a, double b) {
        double c = (a + b) / 2.0;
        return (b - a) / 6.0 * (f.applyAsDouble(a)
                + 4.0 * f.applyAsDouble(c)
                + f.applyAsDouble(b));
    }

    private static double adaptiveSimpson(java.util.function.DoubleUnaryOperator f,
                                          double a, double b, double eps,
                                          int maxRecDepth) {

        double c = (a + b) / 2.0;
        double S = simpson(f, a, b);
        double Sleft = simpson(f, a, c);
        double Sright = simpson(f, c, b);

        return adaptiveSimpsonAux(f, a, b, eps, S, Sleft, Sright, 0, maxRecDepth);
    }

    private static double adaptiveSimpsonAux(
            java.util.function.DoubleUnaryOperator f,
            double a, double b, double eps,
            double S, double Sleft, double Sright,
            int depth, int maxDepth) {

        double c = (a + b) / 2.0;
        double S2 = Sleft + Sright;

        if (depth >= maxDepth || Math.abs(S2 - S) <= 15.0 * eps) {
            return S2 + (S2 - S) / 15.0;
        }

        double leftMid = (a + c) / 2.0;
        double rightMid = (c + b) / 2.0;

        double Sll = simpson(f, a, c);
        double Slr = simpson(f, c, b);

        double left = adaptiveSimpsonAux(f, a, c, eps / 2.0,
                Sleft, simpson(f, a, leftMid), simpson(f, leftMid, c),
                depth + 1, maxDepth);

        double right = adaptiveSimpsonAux(f, c, b, eps / 2.0,
                Sright, simpson(f, c, rightMid), simpson(f, rightMid, b),
                depth + 1, maxDepth);

        return left + right;
    }



    // הפונקציה sub נשארת זהה
    private static double[] sub(double[] p1, double[] p2) {
        int len = Math.max(p1.length, p2.length);
        double[] r = new double[len];
        for (int i = 0; i < len; i++) {
            double a = (i < p1.length ? p1[i] : 0);
            double b = (i < p2.length ? p2[i] : 0);
            r[i] = a - b;
        }
        return r;
    }


    public static double[] getPolynomFromString(String p) {
        if (p == null || p.isEmpty()) return new double[]{0};

        String cleaned = p.replaceAll("\\s+", "");
        if (cleaned.isEmpty()) return new double[]{0};

        // החלפת - ב-+-, אך רק אם הוא לא בתחילת המחרוזת
        if (cleaned.startsWith("-")) {
            cleaned = "0" + cleaned;
        }
        cleaned = cleaned.replace("-", "+-");
        String[] monoms = cleaned.split("\\+");

        int maxDegree = 0;
        List<double[]> parsedMonoms = new ArrayList<>();

        for (String m : monoms) {
            if (m.isEmpty() || m.equals("0")) continue;

            int degree = 0;
            double coeff = 0.0;
            String cleanMonom = m;

            boolean isNegative = m.startsWith("-");
            if (isNegative) cleanMonom = m.substring(1);

            if (cleanMonom.contains("x^")) {
                int indexX = cleanMonom.indexOf("x");
                degree = Integer.parseInt(cleanMonom.substring(cleanMonom.indexOf("^") + 1));
                String coefStr = cleanMonom.substring(0, indexX);

                if (coefStr.isEmpty()) coeff = 1.0;
                else coeff = Double.parseDouble(coefStr);
            }
            else if (cleanMonom.contains("x")) {
                degree = 1;
                int indexX = cleanMonom.indexOf("x");
                String coefStr = cleanMonom.substring(0, indexX);

                if (coefStr.isEmpty()) coeff = 1.0;
                else coeff = Double.parseDouble(coefStr);
            }
            else {
                degree = 0;
                coeff = Double.parseDouble(cleanMonom);
            }

            if (isNegative) coeff = -coeff;

            if (degree > maxDegree) maxDegree = degree;
            parsedMonoms.add(new double[]{coeff, degree});
        }

        double[] ans = new double[maxDegree + 1];
        for (double[] mono : parsedMonoms) {
            ans[(int) mono[1]] += mono[0];
        }

        // צמצום אפסים מובילים
        int effectiveLen = ans.length;
        while (effectiveLen > 1 && Math.abs(ans[effectiveLen - 1]) < EPS) {
            effectiveLen--;
        }

        if (effectiveLen < ans.length) {
            double[] trimmedAns = new double[effectiveLen];
            System.arraycopy(ans, 0, trimmedAns, 0, effectiveLen);
            return trimmedAns;
        }

        return ans;
    }

    public static double[] add(double[] p1, double[] p2) {
        int maxLen = Math.max(p1.length, p2.length);
        double [] ans = new double[maxLen];

        for (int i = 0; i < maxLen; i++) {
            double c1 = (i < p1.length) ? p1[i] : 0;
            double c2 = (i < p2.length) ? p2[i] : 0;
            ans[i] = c1 + c2;
        }

        // צמצום אפסים מובילים
        int effectiveLen = ans.length;
        while (effectiveLen > 1 && Math.abs(ans[effectiveLen - 1]) < EPS) {
            effectiveLen--;
        }

        // אם הצמצום הפך את המערך לפולינום האפס (אורך 1, ערך 0)
        if (effectiveLen == 1 && Math.abs(ans[0]) < EPS) {
            return ZERO;
        }

        if (effectiveLen < ans.length) {
            double[] trimmedAns = new double[effectiveLen];
            System.arraycopy(ans, 0, trimmedAns, 0, effectiveLen);
            return trimmedAns;
        }

        return ans;
    }

    public static double[] mul(double[] p1, double[] p2) {
        if (p1 == null || p2 == null || p1.length == 0 || p2.length == 0) return ZERO;

        double [] ans = new double[p1.length + p2.length - 1];

        for (int i = 0; i < p1.length; i++) {
            for (int j = 0; j < p2.length; j++) {
                ans[i+j] += p1[i] * p2[j];
            }
        }

        // צמצום אפסים מובילים (כפל אפס בפולינום צריך להחזיר רק {0})
        int effectiveLen = ans.length;
        while (effectiveLen > 1 && Math.abs(ans[effectiveLen - 1]) < EPS) {
            effectiveLen--;
        }

        if (effectiveLen < ans.length) {
            double[] trimmedAns = new double[effectiveLen];
            System.arraycopy(ans, 0, trimmedAns, 0, effectiveLen);
            return trimmedAns;
        }

        return ans;
    }

    public static double[] derivative(double[] po) {
        if (po.length <= 1) return new double[]{0};

        double[] ans = new double[po.length - 1];

        for (int i = 1; i < po.length; i++) {
            ans[i - 1] = po[i] * i;
        }

        return ans;
    }
}