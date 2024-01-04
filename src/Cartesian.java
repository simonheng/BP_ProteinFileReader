
/**
 * {@link Cartesian} Represents an object with Cartesian coordinates
 *
 * @author      Antonio Mucherino
 * @author      Simon Hengeveld
 * @since       December 9th, 2021
 * @version     December, 2023
 * @package     ProteinFileReader
 */

import java.util.Arrays;

public class Cartesian
{
   /**
    * A static variable which contains a seed for the generation of unique hashcodes for every {@link Cartesian} object.
    * Gets incremented every time the nextHashCode method is invoked.
    */
   private static int seed = 1;

   /**
    * The hashcode variable
    */
   protected int h;


   /**
    * The dimension (number of coordinates)
    */
   protected int K;

   /**
    * The Cartesian coordinates (available in all implementations)
    */
   protected Double[] c;

   /**
    * A String label for this instance of {@link Cartesian} object.
    * It can be set to the hashcode by default; but implementations may allow it to be changed after initialization.
    */
   protected String name;

   /**
    * Constructs a {@link Cartesian} instance from another instance.
    * The hashcode can also be identical (optional).
    * @param C The input Cartesian instance.
    * @param sameHash indicates whether the hashcode needs to be copied as well (true) or not (false).
    * @exception IllegalArgumentException The input Cartesian instance is null.
    */
   public Cartesian(Cartesian C,boolean sameHash)
   {
      try
      {
         if (C == null) throw new IllegalArgumentException("The input Cartesian instance is null");
         this.K = C.K;
         this.c = new Double[this.K];
         for (int k = 0; k < this.K; k++)  this.c[k] = C.c[k];
         if (sameHash)
            this.h = C.h;
         else
            this.h = this.nextHashCode();
         this.name = String.valueOf(this.h);
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
   }

   /**
    * Constructs a {@link Cartesian} instance from another instance, but with different hashcode.
    * @param C The input Cartesian instance.
    * @exception IllegalArgumentException The input Cartesian instance is null.
    */
   public Cartesian(Cartesian C)
   {
      this(C,false);
   }

   /**
    * Constructs a {@link Cartesian} instance with known dimension.
    * @param K The dimension of the Cartesian instance.
    * @exception IllegalArgumentException Dimension cannot be nonpositive.
    */
   public Cartesian(int K)
   {
      try
      {
         if (K <= 0) throw new IllegalArgumentException("Dimension of Cartesian coordinates cannot be nonpositive");
         this.K = K;
         this.c = new Double[K];
         this.name = String.valueOf(this.h);
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
   }

   /**
    * Constructs a {@link Cartesian} instance from a given array of double.
    * @param coords The array of double with the coordinates.
    * @exception IllegalArgumentException The array of double is null.
    */
   public Cartesian(double[] coords)
   {
      try
      {
         if (coords == null) throw new IllegalArgumentException("The array of double is null");
         this.K = coords.length;
         this.c = new Double[coords.length];
         for (int k = 0; k < this.K; k++)  this.c[k] = coords[k];
         this.h = this.nextHashCode();
         this.name = String.valueOf(this.h);
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
   }

   /**
    * Gives a unique hashcode for Coordinates instances.
    * @return the hashcode
    */
   protected int nextHashCode()
   {
      int h = Cartesian.seed;
      Cartesian.seed++;
      return h;
   }

   /**
    * Verifies whether the internal Cartesian coordinates are defined or not.
    * @return True if the coordinates are defined; False otherwise (boolean).
    */
   public boolean areCoordsDefined()
   {
      for (int k = 0; k < this.K; k++)  if (this.c[k] == null)  return false;
      return true;
   }

   /**
    * Sets the variables of the internal representation.
    * @param coords An array of double containing the new coordinates.
    * @exception IllegalArgumentException The array of double is null.
    * @exception IllegalArgumentException The length of the array of double does not match with this Cartesian instance dimension.
    */
   public void setCoords(double[] coords)
   {
      try
      {
         if (coords == null) throw new IllegalArgumentException("The array of double is null");
         if (this.K != coords.length)
            throw new IllegalArgumentException("The length of the array of double does not match with this Cartesian instance dimension");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      // setting the coordinates
      for (int k = 0; k < this.K; k++)  this.c[k] = coords[k];
   }


   /**
    * Resets the internal variables
    */
   public void reset()
   {
      for (int k = 0; k < this.K; k++)  this.c[k] = null;
   }

   /**
    * Gives a copy of the internal coordinates.
    * @exception IllegalStateException Internal coordinates are not defined.
    * @return an array of double containing a copy of the internal coordinates.
    */
   public double[] getCoords()
   {
      try
      {
         if (!this.areCoordsDefined()) throw new IllegalStateException("Internal coordinates are not defined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      // copying
      double[] coords = new double[this.K];
      for (int k = 0; k < this.K; k++)  coords[k] = this.c[k];
      return coords;
   }

   /**
    * The distance from one {@link Cartesian} object to another, by default: Euclidean distance.
    * @param other The other {@link Cartesian} object.
    * @return the Euclidean distance from (this) to another {@link Cartesian} object.
    * @exception IllegalStateException The Cartesian coordinates vectors do not have the same dimensions.
    */
   public double distanceTo(Cartesian other)
   {
      return Geometry.euclideanDistance(this.getCoords(),other.getCoords());
   }

   /**
    * Creates a clone of this Cartesian instance.
    * @return The clone (Object).
    */
   @Override
   public Object clone()
   {
      return new Cartesian(this,true);
   }

   /**
    * Gives a standard String representation of this Cartesian instance.
    * @return The String representation (String).
    */
   @Override
   public String toString()
   {
      if (!this.areCoordsDefined())  return "[Internal coordinates not defined]";
      return Arrays.toString(this.c);
   }

   /**
    * Sets this instance name (String).
    * @exception IllegalArgumentException The input name is null.
    */
   public void setName(String name)
   {
      try
      {
         if (name == null) throw new IllegalArgumentException("The input name is null");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
      this.name = name;
   }

   /**
    * Gets a copy of this instance name (String).
    */
   public String getName()
   {
      return this.name;
   }

}

