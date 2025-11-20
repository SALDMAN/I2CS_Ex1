
/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class
Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans += c*poly[i];
		}
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
		/** add you code below

		/////////////////// */
		}
		return ans;
	}
	/** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynomial function
	 * @param p2 second polynomial function
	 * @return true iff p1 represents the same polynomial function as p2.
	 */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true;
        /** add you code below

         /////////////////// */
		return ans;
	}

	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
	 */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {
            /** add you code below

             /////////////////// */
		}
		return ans;
	}
    // פונקציה להערכת פולינום בנקודה x
    public static double evaluate(double[] p, double x) {
        double sum = 0;
        for (int i = 0; i < p.length; i++) {
            sum += p[i] * Math.pow(x, i);
        }
        return sum;
    }

    // הפונקציה הראשית למציאת x שבו |p1(x)-p2(x)| < eps
    public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        double f1 = evaluate(p1, x1) - evaluate(p2, x1);
        double f2 = evaluate(p1, x2) - evaluate(p2, x2);

        if (f1 == 0) return x1;
        if (f2 == 0) return x2;

        // ביסקשן
        while ((x2 - x1) > eps) {
            double xm = (x1 + x2) / 2;
            double fm = evaluate(p1, xm) - evaluate(p2, xm);

            if (Math.abs(fm) < eps) {
                return xm; // מצאנו פתרון בתוך הטולרנס
            }

            if (f1 * fm <= 0) {
                x2 = xm;
                f2 = fm;
            } else {
                x1 = xm;
                f1 = fm;
            }
        }

        // מחזירים את האמצע אם לא מצאנו בדיוק
        return (x1 + x2) / 2;
    }

    /**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = x1;
        /** add you code below

         /////////////////// */
		return ans;
	}
	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
	 */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
		double ans = 0;
        /** add you code below

         /////////////////// */
		return ans;
	}
	/**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynomial function.
	 * @return
	 */
    private static String cleanPolynomialString(String p) {
        p = p.replace("-", "+-"); // סימנים
        return p.replaceAll("\\s+", ""); // מסיר רווחים
    }

    private static String[] splitMonoms(String p) {
        return p.split("\\+");
    }

    private static double parseCoefficient(String coefStr) {
        if (coefStr.isEmpty() || coefStr.equals("+")) return 1.0;
        if (coefStr.equals("-")) return -1.0;
        return Double.parseDouble(coefStr);
    }

    private static double[] parseMonom(String m) {
        if (m.isEmpty()) return new double[]{0, 0};

        int degree;
        double coeff;

        if (m.contains("x^")) {
            degree = Integer.parseInt(m.substring(m.indexOf("^") + 1));
            String coefStr = m.substring(0, m.indexOf("x"));
            coeff = parseCoefficient(coefStr);
        } else if (m.contains("x")) {
            degree = 1;
            String coefStr = m.substring(0, m.indexOf("x"));
            coeff = parseCoefficient(coefStr);
        } else {
            degree = 0;
            coeff = Double.parseDouble(m);
        }

        return new double[]{coeff, degree};
    }

    public static double[] getPolynomFromString(String p) {
        if (p == null || p.isEmpty()) return new double[]{0};
        String cleaned = cleanPolynomialString(p);
        String[] monoms = splitMonoms(cleaned);

        // חיפוש דרגה מקסימלית
        int maxDegree = 0;
        for (String m : monoms) {
            double[] parsed = parseMonom(m);
            int deg = (int) parsed[1];
            if (deg > maxDegree) maxDegree = deg;
        }

        double[] ans = new double[maxDegree + 1];

        for (String m : monoms) {
            double[] parsed = parseMonom(m);
            ans[(int) parsed[1]] += parsed[0];
        }

        return ans;
    }

    /**
         * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
         * @param p1
         * @param p2
         * @return
         */
        public static double[] add(double[] p1, double[] p2) {
            double [] ans;
            int l=0;
            if (p1.length>p2.length) {ans=new double[p1.length];}
            else if (p2.length>p1.length) {ans=new double[p2.length];}
            else {ans=new double[p1.length];}
            for (int i = 0; i < p1.length&&i<p2.length; i++) {
                ans[i]=p1[i]+p2[i];
                l=i;
            }
            if (p1.length>p2.length) {
                for (int i =l; i<p1.length;i++) {ans[i]=p1[i];
                }
            }
            else{
                for (int i=l; i<p2.length;i++) {ans[i]=p2[i];}
            }
            return ans;
        }
	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] mul(double[] p1, double[] p2) {
        double [] ans = new double[p1.length*p2.length];
        for (int i = 0; i < p1.length; i++) {
            for (int j = 0; j < p2.length; j++) {
                ans[i+j]+=p1[i]*p2[j];
            }
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
