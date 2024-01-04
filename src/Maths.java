import java.util.Random;
/**
 * {@link Maths} is a static class containing low level methods for mathematical operations that are not included in the
 * basic java library
 *
 * @author      Antonio Mucherino
 * @author      Simon Hengeveld
 * @author      JAMA
 * @since       2 November 2021
 * @version     December, 2023
 * @see         Pair
 * @package     JavaDGP
 */
public class Maths {
   /**
    * Division of two complex numbers x and y
    * @param xr             the real part of x
    * @param xi             the imaginary part of x
    * @param yr             the real part of y
    * @param yi             the imaginary part of y
    * @return               a Pair containing the real (._1()) and imaginary (._2()) of the solution
    */
   public static Pair<Double,Double> cdiv(double xr,double xi,double yr,double yi)
   {
      double r,d;
      double cdivr,cdivi;
      if (Math.abs(yr) > Math.abs(yi))
      {
         r = yi/yr;
         d = yr + r*yi;
         cdivr = (xr + r*xi)/d;
         cdivi = (xi - r*xr)/d;
      }
      else
      {
         r = yr/yi;
         d = yi + r*yr;
         cdivr = (r*xr + xi)/d;
         cdivi = (r*xi - xr)/d;
      }
      return new Pair<> (cdivr,cdivi);
   }

   /**
    * Multiplication of two complex numbers x and y
    * @param xr             the real part of x
    * @param xi             the imaginary part of x
    * @param yr             the real part of y
    * @param yi             the imaginary part of y
    * @return               a Pair containing the real (._1()) and imaginary (._2()) of the solution
    */
   public static Pair<Double,Double> cmul(double xr,double xi,double yr,double yi)
   {
      double real = xr * yr - xi * yi; //the real part
      double imag = xr * yi + xi * yr; //the imaginary part
      return new Pair<> (real,imag);
   }

   /**
    * Computes the combination number (nCr) (this could be further optimized by DP, but not necessary for our uses)
    * @param n              the size of the set
    * @param r              the size of the subsets
    * @return               the number of possible (unordered) subsets we can make from the larger number
    */
   public static int combinationNumber(int n, int r){
      if(r > n - r) r = n - r; // because C(n, r) == C(n, n - r)
      int ans = 1;
      int i;
      for(i = 1; i <= r; i++) {
         ans *= n - r + i;
         ans /= i;
      }
      return ans;
   }

   /**
    * Computes the square root of the sum of squares of two numbers without under or overflow
    * @param a       the first number
    * @param b       the second number
    * @return        the square root of the sum of squares
    */
   public static double hypot(double a, double b) {
      double r;
      if (Math.abs(a) > Math.abs(b)) {
         r = b/a;
         r = Math.abs(a)*Math.sqrt(1+r*r);
      } else if (b != 0) {
         r = a/b;
         r = Math.abs(b)*Math.sqrt(1+r*r);
      } else {
         r = 0.0;
      }
      return r;
   }

   /**
    * Checks whether two doubles are equal given a tolerance epsilon
    * @param a       the first double
    * @param b       the second double
    * @param eps     the tolerance epsilon
    * @return        true when a == b given the tolerance, false otherwise
    */
   public static boolean epsEqual(double a, double b, double eps){
      return Math.abs(a - b) <= eps;
   }
}
