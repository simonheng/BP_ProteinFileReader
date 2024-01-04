
/**
 * {@link Distance} is a class which represents the distance between two Coordinates objects.
 * The distances between these objects may be exact, interval or expected distances.
 *
 * @author   Antonio Mucherino
 * @author   Simon Hengeveld
 * @since    November 2nd, 2021
 * @version  December, 2023
 * @see      Cartesian
 * @package  ProteinFileReader
 */

import java.util.Objects;
import java.util.Random;

public class Distance<T extends Cartesian>
{
   /**
    * The first object associated with the distance
    */
   private T a;

   /**
    * The second object associated with the distance
    */
   private T b;

   /**
    * The computed distance value
    */
   private Double value;

   /**
    * If the distance is an interval distance, represents the lowerbound of the interval.
    * If an upperbound {@link Distance#ub} is not given, represents the expected value of the distance
    */
   private Double lb;

   /**
    * If the distance is an interval distance, represents the upperbound of the interval, null otherwise
    */
   private Double ub;

   /**
    * An optional priority value assigned to the distance
    */
   private Double priority;

   /**
    * An optional label assigned to the distance
    */
   private String label = "";


   /**
    * Constructs the {@link Distance} object given the two involved objects {@link Distance#a} and {@link Distance#b}.
    * The {@link Distance#value} and the expected distance {@link Distance#lb} is set to the computed value based on the coordinates
    * @param a          the first object of the distance
    * @param b          the second object of the distance
    */
   public Distance(T a,T b)
   {
      this.a = a;
      this.b = b;
      this.check();  // cannot create the Distance object if a and b are null because no other info is given
      this.value = null;
      if (a.areCoordsDefined() && b.areCoordsDefined())  this.compute();  // sets up this.value
      this.lb = this.value;
      this.ub = null;
      this.priority = null;
   }

   /**
    * Constructs the {@link Distance} object given the two objects {@link Distance#a} and {@link Distance#b} and a {@link Distance#priority} value.
    * The {@link Distance#value} and the expected distance {@link Distance#lb} is set to the computed value based on the coordinates
    * @param b          the second object of the distance
    * @param priority   desired priority value to be paired with this distance
    */
   public Distance(T a,T b,double priority)
   { 
      this(a,b);
      this.priority = priority;
   }

   /**
    * Constructs the {@link Distance} object given the two objects {@link Distance#a} and {@link Distance#b}.
    * The {@link Distance#value} is set to the computed value based on the coordinates if they are defined.
    * The expected distance is set to the given value (saved in {@link Distance#lb})
    *
    * @param expected   the expected value on the distance
    * @param a          the first object of the distance
    * @param b          the second object of the distance
    */
   public Distance(Double expected,T a,T b)
   {
      this.a = a;
      this.b = b;
      this.check();
      this.lb = expected;
      if (a.areCoordsDefined() && b.areCoordsDefined())  this.compute();  // sets up this.value
   }

   /**
    * Constructs the {@link Distance} object given the two objects {@link Distance#a} and {@link Distance#b} and a {@link Distance#priority} value.
    * The {@link Distance#value} is set to the computed value based on the coordinates.
    * The expected distance is set to the given value (saved in {@link Distance#lb})
    *
    * @param expected   the expected value on the distance
    * @param a          the first object of the distance
    * @param b          the second object of the distance
    * @param priority   desired priority value to be paired with this distance
    */
   public Distance(Double expected,T a,T b,double priority)
   {
      this(expected,a,b);
      this.priority = priority;
   }

   /**
    * Constructs the {@link Distance} object given the two objects {@link Distance#a} and {@link Distance#b}.
    * The {@link Distance#value} is set to the computed value based on the coordinates.
    * The bounds of the interval are set to the given values (saved in {@link Distance#lb} and {@link Distance#ub})
    *
    * @param a          the first object of the distance
    * @param b          the second object of the distance
    * @param lb         the lower bound of the interval distance
    * @param ub         the upper bound of the interval distance
    */
   public Distance(T a,T b,double lb,double ub)
   {
      try
      {
         this.a = a;
         this.b = b;
         this.check();
         this.lb = lb;
         if (this.lb > ub)
            throw new Exception("Lower bound is strictly larger than upper bound");
         else if (this.lb < ub)
            this.ub = ub;
         if (a.areCoordsDefined() && b.areCoordsDefined())  this.compute();  // sets up this.value
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
   }

   /**
    * Constructs the {@link Distance} object given the two objects {@link Distance#a} and {@link Distance#b} and a {@link Distance#priority} value.
    * The {@link Distance#value} is set to the computed value based on the coordinates.
    * The bounds of the interval are set to the given values (saved in {@link Distance#lb} and {@link Distance#ub})
    *
    * @param a          the first object of the distance
    * @param b          the second object of the distance
    * @param lb         the lower bound of the interval distance
    * @param ub         the upper bound of the interval distance
    * @param priority   desired priority value to be paired with this distance
    */
   public Distance(T a,T b,double lb,double ub,double priority)
   {
      this(a,b,lb,ub);
      this.priority = priority;
   }

   /**
    * Accessor for {@link Distance#a}
    * @return The reference to {@link Distance#a}
    */
   public T getA(){
      return a;
   }

   /**
    * Accessor for {@link Distance#b}
    * @return The reference to {@link Distance#b}
    */
   public T getB(){
      return b;
   }

   /**
    * Checks whether the distance is an interval distance
    * @return true if it has a defined interval, false otherwise
    */
   public boolean hasBounds()
   {
      return this.ub != null;
   }

   /**
    * Accessor for the lower bound {@link Distance#lb} (throws exception if the distance is not an interval distance
    * @return the lower bound {@link Distance#lb}
    */
   public double getLowerBound()
   {
      try
      {
         if (!this.hasBounds()) throw new Exception("Bounds on distance value are unknown");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return this.lb;
   }

   /**
    * Accessor for the upper bound {@link Distance#ub} (throws exception if the distance is not an interval distance
    * @return the upper bound {@link Distance#ub}
    */
   public double getUpperBound()
   {
      try
      {
         if (!this.hasBounds()) throw new Exception("Bounds on distance value are unknown");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return this.ub;
   }

   /**
    * Checks whether or not the distance has an associated expected value
    * @return if the expected value is defined, false otherwise
    */
   public boolean hasExpectedValue()
   {
      return this.lb != null && this.ub == null;
   }

   /**
    * Accessor for the expected value (throws exception if the distance is an interval distance, or there is no expected value defined)
    * @return the expected value {@link Distance#lb}
    */
   public double getExpectedValue()
   {
      try
      {
         if (!this.hasExpectedValue()) throw new Exception("Expected value of distance is unknown");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return this.lb;
   }

   /**
    * Checks whether or not the distance has an associated {@link Distance#priority} value
    * @return if the {@link Distance#priority} value is defined, false otherwise
    */
   public boolean hasPriority()
   {
      return this.priority != null;
   }

   /**
    * Accessor for the {@link Distance#priority} value (throws exception if the distance has no associated priority value)
    * @return the associated {@link Distance#priority} value
    */
   public double getPriority()
   {
      try
      {
         if (!this.hasPriority()) throw new Exception("Distance priority was not set");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return this.priority;
   }

   /**
    * Mutator for the {@link Distance#label} variable
    * @param label      the new label to update to
    */
   public void setLabel(String label)
   {
      this.label = label;
   }

   /**
    * Accessor for the {@link Distance#label} value
    * @return           the associated {@link Distance#label}
    */
   public String getLabel()
   {
      return label;
   }

   /**
    * Introducing noise on the distance values
    * @param eps the noise amplitude (must be nonnegative)
    * @exception IllegalArgumentException The given eps value is negative
    */
   public void noisify(double eps)
   {
      try
      {
         if (eps < 0.0) throw new IllegalArgumentException("The specified noise amplitude is negative");
         if (eps == 0.0)  return;
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      double r = new Random().nextDouble();
      if (this.hasExpectedValue())
      {
         // the distance was originally exact
         this.ub = this.lb + r*eps;
         this.lb = this.lb - (1.0 - r)*eps;
      }
      else
      {
         // the distance was already an interval
         this.lb = this.lb - r*eps;
         this.ub = this.ub + (1.0 - r)*eps;
      }
   }

   /**
    * Checks whether the current distance {@link Distance#value} satisfies the bounds or the expected value, given a certain tolerance.
    * NOTE: the bounds are inclusive
    * @param eps        the given tolerance
    * @return           true if the distance is satisfied, false otherwise
    */
   public boolean isSatisfied(double eps)
   {
      // We slightly increase the tolerance to allow for a rounding error
      double roundingError = 1e-12;
      eps += roundingError;
      this.compute();

      try
      {
         if (this.value == null) throw new Exception("Distance cannot be verified: object coordinates not defined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      if (this.hasBounds() && this.value >= lb - eps && this.value <= ub + eps)  return true;
      return this.hasExpectedValue() && Maths.epsEqual(this.value, this.lb, eps);
   }

   /**
    * Checks whether the current distance value satisfies the bounds or the expected value, given a tolerance of 0
    * @return true if the distance is satisfied, false otherwise.
    */
   public boolean isSatisfied()
   {
      return this.isSatisfied(0.0);
   }

   /**
    * Gives an estimation of the current error on the distance (0 if it is satisfied)
    *  @return the estimated error
    */
   public double absError()
   {
      this.compute();
      if (this.hasExpectedValue())  return Math.abs(this.getExpectedValue() - this.value);
      double v1 = 0.0;
      if (this.value < this.lb)  v1 = this.lb - this.value;
      double v2 = 0.0;
      if (this.value > this.ub)  v2 = this.value - this.ub;
      return Math.max(v1,v2);
   }

   /**
    * Gives an estimation of the current relative error on the distance (0 if it is satisfied)
    *  @return the estimated relative error
    */
   public double relativError()
   {
      double err = this.absError();
      if (this.hasExpectedValue())
      {
         if (this.getExpectedValue() != 0.0)  err = err/this.getExpectedValue();
      }
      else
      {
         double div = 0.5*(this.getLowerBound() + this.getUpperBound());
         if (div != 0.0)  err = err/div;
      }
      return err;
   }

   /**
    * Gives an estimation of the current "squared" error on the distance (0 if it is satisfied)
     * @return the estimated squared error
     */
   public double squaredError()
   {
      try
      {
         if (this.value == null) throw new Exception("Distance cannot be verified: object coordinates not defined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      // updating the distance value (if necessary)
      this.compute();

      // expected value case
      if (this.hasExpectedValue())
      {
         double a = this.getExpectedValue();
         double b = this.value;
         double c = b - a;
         return c*c;
      }

      // interval case
      double v1 = 0.0; 
      if (this.value < this.lb)
      {
         v1 = this.lb - this.value;
         v1 = v1*v1;
      }
      double v2 = 0.0; 
      if (this.value > this.ub)
      {
         v2 = this.value - this.ub;
         v2 = v2*v2;
      }
      return Math.max(v1,v2);
   }

   /**
    * Computes the absolute value difference between the current distance value and the input argument
    */
   public double absDiff(double d)  // [NOT TESTED IN THE MAIN]
   {
      this.compute();
      return Math.abs(this.value - d);
   }

   /**
    * Computes and updates the distance between a and b
    * This distance is saved in {@link Distance#value}
    */
   public void compute()
   {
      this.check();
      this.value = a.distanceTo(b);
   }

   /**
    * Projects the given distance on the real distance
    * (bounds or expected distance)
    * exception on eps not checked!
    * @param d the distance to project
    * @param eps small error on projection
    * @return the projected distance value
    */
   public double project(double d,double eps)
   {
      double projection = 0.0;
      if (this.hasExpectedValue())
      {
         projection = this.getExpectedValue();
      }
      else
      {
         double lb = this.getLowerBound();
         double ub = this.getUpperBound();
         if (d < lb)
            projection = lb - eps;
         else if (d > ub)
            projection = ub + eps;
         else
            projection = d; 
      }
      return projection;
   }

   /**
    * Converts the {@link Distance} object to a String
    * @return a String value representing Distance object
    */
   public String toString()
   {
      String print = "[";
      if (this.hasExpectedValue())
         print = print + this.lb + " -> " + this.value;
      else
         print = print + this.lb + " (" + this.value + ") " + this.ub;
      print = print + "]";
      if (this.hasPriority())
         print = print + " {" + this.priority + "}";
      return print;
   }

   /**
    * Checks whether at least one of the dirty bits of {@link Distance#a} and {@link Distance#b} are set to true.
    * It throws an exception when {@link Distance#a} or {@link Distance#b} is undefined.
    */
   private void check()
   {
      try
      {
         if (this.a == null) throw new Exception("First object involved in distance computation is null");
         if (this.b == null) throw new Exception("Second object involved in distance computation is null"); 
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
   }

   /**
    * We override the hashcode so that it only depends on a and b, and not on other parameters
    * @return the hashcode
    */
   @Override
   public int hashCode()
   {
      return Objects.hash(this.a,this.b);
   }
}

