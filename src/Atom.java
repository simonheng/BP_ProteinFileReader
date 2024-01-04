/**
 * {@link Atom} is a class extending {@link Cartesian} adding behaviour for the atomic case.
 * @author  Antonio Mucherino
 * @author  Simon Hengeveld
 * @since   December 9th, 2021
 * @version December 2023
 * @see     Cartesian
 * @see     Coordinates
 * @package ProteinFileReader
 */

import java.util.*;

public class Atom extends Cartesian
{
   /**
    * An enumerated type offering the predefined constant atomic types that we will allow
    */
   public enum AtomType {C, H, N, O, P, S}

   /**
    * A static variable which maps an {@link AtomType} to a van der Waals radius
    * The van der Waals radii are taken from:
    * A.Bondi, "van der Waals Volumes and Radii", 1964.
    */
   public static final Map<Atom.AtomType, Double> vanDerWaalsRadius = new HashMap<Atom.AtomType, Double>() {{
      // Atom type     // Radius
      put(AtomType.C,  1.7);
      put(AtomType.H,  1.2);
      put(AtomType.N,  1.55);
      put(AtomType.O,  1.52);
      put(AtomType.P,  1.9);
      put(AtomType.S,  1.80);
   }};

   /**
    * The (immutable) {@link AtomType} of the atom
    */
   private AtomType type;


   /**
    * Optional: the name of the amino acid it is in
    */
   private String residue;

   /**
    * Optional: the ID of the amino acid in the chain of a protein
    */
   private Integer residueID;

   /**
    * Optional: chain code (PDB information)
    */
   private Character chainCode;

   /**
    * Optional: PDB id (PDB information)
    */
   private Integer PDBid;

   /**
    * Optional: occupancy (PDB information)
    */
   private Double occupancy;

   /**
    * Optional: temperature (PDB information)
    */
   private Double temperature;

   /**
    * Constructs an {@link Atom} given another {@link Atom} (copy constructor).
    * @param a The {@link Atom} to copy from.
    * @param sameHash True if the copy instance needs to share the same hashcode of the original instance (boolean).
    * @exception IllegalArgumentException The input Cartesian instance is null.
    */
   public Atom(Atom a,boolean sameHash)
   {
      super(a,sameHash);
      this.type = a.type;
      this.name = a.name;
   }

   /**
    * Constructs an {@link Atom} given another {@link Atom} (copy constructor assigning a new hashcode).
    * @param a The {@link Atom} to copy from.
    * @exception IllegalArgumentException The input Cartesian instance is null.
    */
   public Atom(Atom a)
   {
      this(a,false);
   }

   /**
    * Constructs an {@link Atom} with the minimal required information
    * @param type      a {@link String} matching an {@link AtomType}
    * @param name      a {@link String} denoting the name of the {@link Atom}
    */
   public Atom(String type,String name) {
      super(3);
      try {
         if (type == null) throw new Exception("Type String is null");
         if (type.length() == 0) throw new Exception("Type String is empty");
         this.type = AtomType.valueOf(type);
         if (name == null) throw new Exception("Specified Atom name is null");
         if (name.length() == 0) throw new Exception("Specified Atom name contains no characters");
         this.name = name;
      } catch (Exception e) {
         e.printStackTrace();
         System.exit(1);
      }
   }

   /**
    * Constructs an {@link Atom} with coordinates and a biased hashcode
    * @param type      a {@link String} matching a {@link AtomType}
    * @param name      a {@link String} denoting the name of the {@link Atom}
    * @param x,y,z     doubles describing the three coordinates of the position of the {@link Atom}
    */
   public Atom(String type,String name,double x,double y,double z)
   {
      super(new double[]{x,y,z});
      try
      {
         if (type == null) throw new Exception("Type String is null");
         if (type.length() == 0) throw new Exception("Type String is empty");
         this.type = AtomType.valueOf(type);
         if (name == null) throw new Exception("Specified Atom name is null");
         if (name.length() == 0) throw new Exception("Specified Atom name contains no characters");
         this.name = name;
      } catch (Exception e) {
         e.printStackTrace();
         System.exit(1);
      }
   }

   /**
    * Constructs an {@link Atom} with coordinates and a biased hashcode
    * @param type      a {@link String} matching a {@link AtomType}
    * @param name      a {@link String} denoting the name of the {@link Atom}
    * @param x,y,z     doubles describing the three coordinates of the position of the {@link Atom}
    * @param residue   the three letter code of the residue
    * @param residueID the residue ID in the amino acid chain
    * */
   public Atom(String type,String name,double x,double y,double z, String residue, int residueID)
   {
      // we use the residue name + id as hashcode!
      this(type, name, x,y,z);
      this.residue = residue;
      this.residueID = residueID;
   }

   /**
    * Constructs an {@link Atom} with residue information
    * @param type      a {@link String} matching a {@link AtomType}
    * @param name      a {@link String} denoting the name of the {@link Atom}
    * @param residue   the three letter code of the residue
    * @param residueID the residue ID in the amino acid chain
    * */
   public Atom(String type,String name, String residue, int residueID)
   {
      // we use the residue name + id as hashcode!
      this(type, name);
      this.residue = residue;
      this.residueID = residueID;
   }

   /**
    * Constructor from PDB, extra information about the {@link Atom}
    * @param type      a {@link String} matching a {@link AtomType}
    * @param name      a {@link String} denoting the name of the {@link Atom}
    * @param x,y,z     doubles describing the three coordinates of the position of the {@link Atom}
    * @param residue   the three letter code of the residue
    * @param residueID the residue ID in the amino acid chain
    * @param PDBid     the PDB ID of the atom
    * @param chainCode     the PDB identifier for the chain
    */
   public Atom(String type,String name,double x,double y,double z, String residue, int residueID, int PDBid, char chainCode)
   {
      // we use the residue name + id as hashcode!
      this(type, name, x,y,z, residue, residueID);
      this.PDBid = PDBid;
      this.chainCode = chainCode;
   }

   /**
    * The accessor of the {@link AtomType} {@link Atom#type} variable.
    * @return the {@link Atom#type} associated with the {@link Atom} object.
    */
   public AtomType getType()
   {
      return this.type;
   }

   /**
    * The accessor of the {@link String} version of the {@link AtomType} {@link Atom#type} variable.
    * @return the {@link String} value of the {@link Atom#type} associated with the {@link Atom} object.
    */
   public String getTypeName()
   {
      return this.type.name();
   }

   /**
    * Sets the {@link Atom#residue} variable.
    * @param name The residue name (String).
    * @exception IllegalArgumentException The input String is null.
    */
   public void setResidueName(String name)
   {
      try
      {
         if (name == null) throw new IllegalArgumentException("The input String is null");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      this.residue = name;
   }

   /**
    * The accessor of the {@link Atom#residue} variable.
    * @return The {@link Atom#residue} associated with the {@link Atom} object.
    * @exception IllegalStateException The residue name is undefined.
    */
   public String getResidueName()
   {
      try
      {
         if (this.residue == null) throw new IllegalStateException("The residue name is undefined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return new String(this.residue);
   }

   /**
    * The accessor of the {@link Atom#residueID} variable.
    * @return The {@link Atom#residueID} associated with the {@link Atom} object.
    * @exception IllegalStateException The residue id is undefined.
    */
   public int getResidueID()
   {
      try
      {
         if (this.residueID == null) throw new IllegalStateException("The residue id is undefined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return (int) this.residueID;
   }

   /**
    * The accessor of the {@link Atom#chainCode} variable.
    * @return The {@link Atom#chainCode} associated with the {@link Atom} object.
    * @exception IllegalStateException The chain code is undefined.
    */
   public char getChainCode()
   {
      try
      {
         if (this.chainCode == null) throw new IllegalStateException("The chain code is undefined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return (char) this.chainCode;
   }

   /**
    * The accessor of the {@link Atom#PDBid} variable.
    * @return The {@link Atom#PDBid} associated with the {@link Atom} object.
    * @exception IllegalStateException The PDB id is undefined.
    */
   public int getPDBid()
   {
      try
      {
         if (this.PDBid == null) throw new IllegalStateException("The PDB id is undefined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return (int) this.PDBid;
   }

   /**
    * Sets the occupancy variable.
    * @param occupancy The new occupancy value (double).
    */
   public void setOccupancy(double occupancy)
   {
      this.occupancy = occupancy;
   }

   /**
    * The accessor of the {@link Atom#occupancy} variable.
    * @return The {@link Atom#occupancy} associated with the {@link Atom} object.
    * @exception IllegalStateException The occupancy is undefined.
    */
   public double getOccupancy()
   {
      try
      {
         if (this.occupancy == null) throw new IllegalStateException("The occupancy is undefined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return (double) this.occupancy;
   }

   /**
    * Sets the temperature variable.
    * @param temperature The new temperature value (double).
    */
   public void setTemperature(double temperature)
   {
      this.temperature = temperature;
   }

   /**
    * The accessor of the {@link Atom#temperature} variable.
    * @return The {@link Atom#temperature} associated with the {@link Atom} object.
    * @exception IllegalStateException The temperature is undefined.
    */
   public double getTemperature()
   {
      try
      {
         if (this.temperature == null) throw new IllegalStateException("The temperature is undefined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return (double) this.temperature;
   }

   /**
    * The accessor of the x-coordinate of the {@link Atom} object.
    * @return The x-coordinate associated with the {@link Atom} object.
    * @exception IllegalStateException The Atom's Cartesian coordinates are undefined.
    */
   public double getX()
   {
      try
      {
         if (!this.areCoordsDefined()) throw new IllegalStateException("The Atom's Cartesian coordinates are undefined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return this.c[0];
   }

   /**
    * The accessor of the y-coordinate of the {@link Atom} object.
    * @return The y-coordinate associated with the {@link Atom} object.
    * @exception IllegalStateException The Atom's Cartesian coordinates are undefined.
    */
   public double getY()
   {
      try
      {
         if (!this.areCoordsDefined()) throw new IllegalStateException("The Atom's Cartesian coordinates are undefined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return this.c[1];
   }

   /**
    * The accessor of the z-coordinate of the {@link Atom} object.
    * @return The z-coordinate associated with the {@link Atom} object.
    * @exception IllegalStateException The Atom's Cartesian coordinates are undefined.
    */
   public double getZ()
   {
      try
      {
         if (!this.areCoordsDefined()) throw new IllegalStateException("The Atom's Cartesian coordinates are undefined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return this.c[2];
   }

   /**
    * Creates a clone of this Atom instance.
    * @return The clone (Object).
    */
   @Override
   public Object clone()
   {
      return new Atom(this,true);
   }

   /**
    * Provides a {@link String} description of the {@link Atom} object.
    * @return {@link String} describing the {@link Atom}.
    */
   public String toString()
   {
      String print = "(" + this.getName();
      if (this.residue != null)  print += " " + this.residue;
      if (this.residueID != null)  print += " " + this.residueID;
      return print + ")";
   }
}

