/**
 * {@link Pair} is a class which represents a Tuple
 *
 * @author   Antonio Mucherino
 * @author   Simon Hengeveld
 * @since    2021
 * @version  December, 2023
 * @package  ProteinFileReader
 */

public class Pair<Tx,Ty>
{
   /**
    * The first object
    */
   private final Tx x;

   /**
    * The second object
    */
   private final Ty y;

   /**
    * Constructs the pair
    * @param x the first object
    * @param y the second object
    */
   public Pair(Tx x,Ty y)
   {
      this.x = x;
      this.y = y;
   }


   /**
    * Getter for first object
    * @return the first object
    */
   public Tx _1()
   {
      return this.x;
   }

   /**
    * Getter for first object class
    * @return the first object class
    */
   public Class _class_1()
   {
      if (this.x == null)  return null;
      return this.x.getClass();
   }

   /**
    * Getter for second object
    * @return the second object
    */
   public Ty _2()
   {
      return this.y;
   }

   /**
    * Getter for second object class
    * @return the second object class
    */
   public Class _class_2()
   {
      if (this.y == null)  return null;
      return this.y.getClass();
   }

   /**
    * Checks if this is equal to another object o
    * @param o the other object
    * @return true if they are equal, false o
    */
   @Override
   public boolean equals(Object o)
   {
      boolean E1 = false;
      boolean E2 = false;

      if (o instanceof Pair)
      {
         Pair p = (Pair) o;
         if (this.x == null || p.x == null)
         {
            if (this.x == null && p.x == null)  E1 = true;
         }
         else
         {
            if (this.x.getClass() == p.x.getClass())  E1 = this.x.equals(p.x);
         }

         if (this.y == null || p.y == null)
         {
            if (this.y == null && p.y == null)  E2 = true;
         }
         else
         {
            if (this.y.getClass() == p.y.getClass())  E2 = this.y.equals(p.y);
         }
      }
      return E1 && E2;
   }

   /**
    * Computes a hashcode based on the two objects in the pair. Note: (x,y) != (y,x).
    * @return the generated hashcode
    */
   @Override
   public int hashCode()
   {
      return this.x.hashCode() + 2*this.y.hashCode();
   }

   /**
    * Prints the pair
    * @return a description of the two objects
    */
   public String toString()
   {
      return "(" + this.x + "," + this.y + ")";
   }
}

