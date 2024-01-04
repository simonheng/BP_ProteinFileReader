
/* javaDGP project
 *
 * Simon Hengeveld and Antonio Mucherino
 *
 * last update: July 24th, 2022.
 */

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Class that will read and write PDB files
 */

public class PDB
{
   /**
    * The original file name
    */
   private String filename;

   /**
    * The target chain
    */
   private char chain;

   /**
    * The header (first line of the PDB file)
    */
   private String header;

   /**
    * The ID of the loaded model
    */
   private int modelID;

   /**
    * The list of Atoms
    */
   private LinkedList<Atom> atoms;

   /**
    * The master line (last line of the PDB file)
    */
   private String master;

   /**
    * This list contains the lines that describe the secondary structures, and of which amino acids they consist
    */
   public ArrayList<String> structures;

   /**
    * This is the maximal inter-atomic distance present in the PDB file
    */
   private Double maxDistance;

   /**
    * The first residue ID encountered in the original PDB file given the chain provided (not necessarily 1!)
    */
   public Integer firstResidueID = null;

   /**
    * The last stored residue ID. Does not use the original file numbering, but the numbering starting at 1!
    */
   public int lastResidueID = -1;

   /**
    * Constructs the PDB object by reading from a file
    * No chain provided (default = 'n') and no folder provided
    * @param filename   the PDB file location
    * @param modelID    the model to load (0 if any)
    */
   public PDB(String filename,int modelID){
      this(filename, modelID, 'n');
   }

   /**
    * Constructs the PDB object by reading from a file
    * No folder provided (default = 'PDB_files/')
    * @param filename   the PDB file location
    * @param modelID    the model to load (0 if any)
    * @param chain      the target chain code
    */
   public PDB(String filename,int modelID, char chain)
   {
      this(filename, modelID, chain, "PDB_files/");
   }

   /**
    * Constructs the PDB object by reading from a file
    * Note on residue IDs: when stored in the instance, they are always from 1 to this.lastResidueID
    * Instead, when we save
    * @param filename   the PDB file location
    * @param modelID    the model to load (0 if any)
    * @param chain      the target chain code, 'n' in case we want all chains considered!
    * @param folder     the name of the folder containing the PDB file...
    */
   public PDB(String filename,int modelID, char chain, String folder)
   {
      try
      {
         if (filename == null) throw new IllegalArgumentException("The String supposing to contain the file name is null");
         if (filename.length() == 0) throw new IllegalArgumentException("The String supposing to contain the file name is empty");
         this.filename = filename;
         if (modelID < 0) throw new IllegalArgumentException("The model ID cannot be a negative integer");
         this.modelID = modelID;
         this.chain = chain;

         // initializing the lists
         this.atoms = new LinkedList<Atom> ();
         this.structures = new ArrayList<String> ();

         // reading the file
         File input = new File(folder + filename);
         Scanner scan = new Scanner(input);
         boolean skipCurrentModel = false;
         while (scan.hasNext())
         {
            // the current line
            String line = scan.nextLine();

            // split on multiple whitespaces
            String[] parts = line.split("\\s+");
            if (parts[0].equals("MODEL"))
            {

               // we may not be interested in the current model
               if (modelID != 0 && parts[1].length() > 0)
               {
                  int ID = Integer.parseInt(parts[1]);

                  //We went too far!
                  if(ID > modelID)
                     break;
                  skipCurrentModel = (ID != modelID);
               }
            }
            else if (parts[0].equals("ENDMDL"))
            {
               // we make sure to load only one model
               if (modelID == 0)  skipCurrentModel = true;
            }
            else if (parts[0].equals("HEADER"))
            {
               this.header = line;
            }
            else if (parts[0].equals("MASTER"))
            {
               this.master = line;
            }
            else if (parts[0].equals("SEQRES") || parts[0].equals("HELIX") || parts[0].equals("SHEET"))
            {
               // For now we just copy the structure strings.
               // Perhaps add more here later, so that we can use some of this information?
               this.structures.add(line);
            }
            else if (!skipCurrentModel)
            {
               // loading the atom info
               if(parts[0].equals("ANISOU")){
                  //System.out.println("Check....");
                  //throw new IllegalStateException("Oops");
               }
               if (parts[0].equals("ATOM"))
               {
                  if (parts.length < 9 && parts[2].length() > 4)
                  {
                     String[] temp = new String[parts.length+1];
                     temp[0] = parts[0];
                     temp[1] = parts[1];
                     temp[2] = parts[2].substring(0, 4);
                     temp[3] = parts[2].substring(4);
                     for (int i = 4; i < temp.length; i++)  temp[i] = parts[i - 1];
                     parts = temp;
                  }
                  if (parts.length < 9 && parts[6].length() > 9)
                  {
                     String[] temp = new String[9];
                     for (int i = 0; i < 6; i++)  temp[i] = parts[i];
                     temp[6] = parts[6].substring(0, parts[6].length() - 8);
                     temp[7] = parts[6].substring(parts[6].length() - 8);
                     temp[8] = parts[7];
                     parts = temp;
                  }

                  String type;
                  if(parts.length == 13){
                     type = parts[12];
                  }
                  else if (parts.length > 11)
                     type = parts[11];
                  else
                     type = Character.toString(parts[2].charAt(0));

                  if(type.length() > 1){
                     type = Character.toString(type.charAt(0));
                  }

                  if(parts[2].length() > 6){
                     String[] oldparts = parts;
                     parts = new String[parts.length + 1];

                     parts[0] = oldparts[0];
                     parts[1] = oldparts[1];

                     if(oldparts[2].length() == 7) {
                        parts[2] = oldparts[2].substring(0, 3);
                        parts[3] = oldparts[2].substring(3);
                     }
                     else{
                        parts[2] = oldparts[2].substring(0, 4);
                        parts[3] = oldparts[2].substring(4);
                     }
                     for(int i = 4; i < oldparts.length;i++)
                        parts[i] = oldparts[i-1];
                  }

                  int index = 4;

                  char resChain = 'n';
                  if (Character.isLetter(parts[index].charAt(0))) //sometimes the chain info is missing...
                     resChain = parts[index++].charAt(0);


                  if(resChain != this.chain && this.chain != 'n')
                     continue;

                  int resID;

                  if(type.equals("A") || type.equals("D"))
                     continue;

                  int chainAttached = 0;
                  if(parts[index-1].length() > 3){
                     String resIDstr = parts[index-1].substring(1);
                     resID = Integer.parseInt(resIDstr);
                     chainAttached = -1;
                  }
                  else{
                     String resIDstr = parts[index++];

                     if(!resIDstr.contains("-"))
                        resIDstr = resIDstr.replaceAll("[^\\d.]", "");

                     resID = Integer.parseInt(resIDstr);
                  }

                  if(filename.contains("5SVY") && resID == 0)
                     continue;


                  if(firstResidueID == null) {
                     firstResidueID = resID;
                  }

                  if(firstResidueID != null)
                     resID = resID - firstResidueID + 1;

                  if(resID > lastResidueID)
                     lastResidueID = resID;

                  if(parts[7].length() > 7){
                     String[] oldParts = parts;
                     parts = new String[parts.length + 1];

                     for(int i = 0; i < 7; i++)
                        parts[i] = oldParts[i];
                     parts[7] = oldParts[7].substring(0,7);
                     parts[8] = oldParts[7].substring(7);

                     for(int i = 9; i < oldParts.length; i++)
                        parts[i] = oldParts[i-1];
                  }

                  Atom a = new Atom(
                           type,
                          parts[2],
                           Double.parseDouble(parts[6 + chainAttached]),
                           Double.parseDouble(parts[7 + chainAttached]),
                           Double.parseDouble(parts[8 + chainAttached]),
                           parts[3],
                           resID,
                           Integer.parseInt(parts[1]),
                           resChain
                   );
                  index = index + 3;
                  if (parts.length > 9) // not all pdb files have these two values!
                  {
                     a.setOccupancy(Double.parseDouble(parts[index++]));
                     a.setTemperature(Double.parseDouble(parts[index]));
                  }

                  // adding atom in the list
                  this.atoms.add(a);
               }
            }
         }

         // model loaded?
         if (skipCurrentModel && this.atoms.size() == 0) throw new IllegalArgumentException("The requested model was not found in the PDB file");
         if(this.atoms.size() == 0) throw new IllegalStateException("No atoms were found in the PDB file");

      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
   }

   /**
    * Counts the number of atoms in the PDB instance
    * @return the number of atoms
    */
   public int numberOfAtoms()
   {
      return this.atoms.size();
   }

   /**
    * Accesses the original PDB file name
    * @return the file name
    */
   public String getOriginalFileName()
   {
      return this.filename;
   }

   /**
    * Computes the phi dihedral angle for a given residue
    * @param residueID the residue ID (from 1 to this.lastResidueID)
    * @return the phi angle in radians
    */
   public Double computePhi(int residueID){
      Atom C = this.getAtom("C", residueID);
      Atom Np = this.getAtom("N", residueID+1);
      Atom CAp = this.getAtom("CA", residueID+1);
      Atom Cp = this.getAtom("C", residueID+1);

      if(C == null || Np == null || CAp == null | Cp == null) return null;
      return Geometry.signedDihedralAngle(C.getCoords(),Np.getCoords(), CAp.getCoords(), Cp.getCoords());
   }

   /**
    * Computes the psi dihedral angle for a given residue
    * @param residueID the residue ID (from 1 to this.lastResidueID)
    * @return the psi angle in radians
    */   public Double computePsi(int residueID){
      Atom N = this.getAtom("N", residueID);
      Atom CA = this.getAtom("CA", residueID);
      Atom C = this.getAtom("C", residueID);
      Atom Np = this.getAtom("N", residueID+1);

      if(N == null || CA == null || C == null | Np == null) return null;
      return Geometry.signedDihedralAngle(N.getCoords(),CA.getCoords(), C.getCoords(), Np.getCoords());
   }

   /**
    * Computes the omega dihedral angle for a given residue
    * @param residueID the residue ID (from 1 to this.lastResidueID)
    * @return the omega angle in radians
    */
   public Double computeOmega(int residueID){
      Atom CA = this.getAtom("CA", residueID);
      Atom C = this.getAtom("C", residueID);
      Atom Np = this.getAtom("N", residueID+1);
      Atom CAp = this.getAtom("CA", residueID+1);

      if(CA == null || C == null || Np == null | CAp == null) return null;
      return Geometry.signedDihedralAngle(CA.getCoords(),C.getCoords(), Np.getCoords(), CAp.getCoords());
   }

   /**
    * Get header String in PDB file
    * @return the header String
    */
   public String getHeader()
   {
      return this.header;
   }

   /**
    * Get model ID
    * @return the model ID
    */
   public int getOriginalModelID()
   {
      return this.modelID;
   }

   /**
    * Get master String in PDB file
    * @return the master String
    */
   public String getMaster()
   {
      return this.master;
   }

   /**
    * Get structures String in PDB file
    * @return the structures String
    */
   public List<String> getStructures()
   {
      return this.structures;
   }

   /**
    * Get list of all {@link Atom}s in model
    * @return the list of {@link Atom}s
    */
   public List<Atom> getAtoms()
   {
      return this.atoms;
   }

   /**
    * Returns the total number of residues in the PDB that matched the model and chain code
    * @return the total number of atoms
    */
   public int numberOfResidues(){
      return  this.lastResidueID;
   }

   /**
    * Access an {@link Atom} given its name and the residue ID that it is in
    * @param atomName the name of the {@link Atom}
    * @param residueID the residue ID
    * @return the {@link Atom} in the PDB model
    */
   public Atom getAtom(String atomName, int residueID){
      for(Atom a : getAtoms()) if(a.getName().equals(atomName) && a.getResidueID() == residueID && (a.getChainCode() == this.chain || this.chain == 'n')) return a;
      return null;
   }

   /**
    * Constructs a {@link Matrix} of (cartesian) coordinates from the atoms inside the PDB file
    * @return the Coordinate matrix
    */
   public Matrix coordinateMatrix(List<Atom> atoms)
   {
      double[][] coordinates = new double[atoms.size()][3];
      for (int i = 0; i < atoms.size(); i++)
         coordinates[i] = atoms.get(i).getCoords();
      return new Matrix(coordinates);
   }

   /**
    * Constructs a {@link Matrix} of (cartesian) coordinates from the atoms inside the PDB file
    * @return the Coordinate matrix
    */
   public Matrix coordinateMatrix()
   {
      return this.coordinateMatrix(this.atoms);
   }


   /**
    * Computes the RMSD score between this model and another, up to a certain number of residues and only counting certain atoms
    * @param other the other PDB model
    * @param atomNames the atom names that we want to count
    * @param numberOfRes up to this residue
    * @return the computed RMSD score
    */
   public double RMSD(PDB other, List<String> atomNames, int numberOfRes){
      List<Atom> filteredA = new ArrayList<>();
      List<Atom> filteredB = new ArrayList<>();

      Matrix A = null;
      Matrix B = null;

      try
      {
         for(Atom a : this.atoms) {
            if((a.getChainCode() == this.chain || this.chain == 'n') && atomNames.contains(a.getName()) && a.getResidueID() <= numberOfRes){
               //corresponding resID

               //find corresponding atom in other PDB
               Atom b = other.getAtom(a.getName(), a.getResidueID());

               if(b != null) {
                  filteredA.add(a);
                  filteredB.add(b);
               }
            }
         }
         A = this.coordinateMatrix(filteredA);
         B = other.coordinateMatrix(filteredB);
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
      Geometry.kabschAlignment(A,B);
      return Geometry.RMSD3D(B,A);
   }

   /**
    * Computes the RMSD score between this model and another only counting certain atoms
    * @param other the other PDB model
    * @param atomNames the atom names that we want to count
    * @return the computed RMSD score
    */
   public double RMSD(PDB other, List<String> atomNames)
   {
     int numberOfRes = this.numberOfResidues();
     return this.RMSD(other, atomNames, numberOfRes);
   }

   /**
    * Gets the max distance
    * If it is null, computes it first.
    */
   public double getMaxDistance()
   {
      if (this.maxDistance == null)
      {
         // compute the max distance
         this.maxDistance = 0.0;
         for (int i = 0; i < this.atoms.size(); i++)
            for (int j = i + 1; j < this.atoms.size(); j++)
               this.maxDistance = Math.max(this.maxDistance, atoms.get(i).distanceTo(atoms.get(j)));
      }
      return this.maxDistance;
   }

   /**
    * Saves the PDB object to a PDB file or a MDjeep file
    * @param format the output format (pdb or mdjeep)
    * @param filename the output filename
    * @throws IOException
    * @throws IllegalArgumentException Not a valid format
    */
   public void save(String format, String filename)
   {
      Locale.setDefault(Locale.US);
      try {
         if(!format.equalsIgnoreCase("pdb") && !format.equalsIgnoreCase("mdjeep"))
            throw new IllegalArgumentException("Not a valid format");

         File outputFile = new File(filename);
         outputFile.createNewFile();
         FileWriter myWriter = new FileWriter(outputFile);

         if(format.equalsIgnoreCase("pdb")) {
            // write headers and structures
            if (this.header != null) myWriter.write(this.header + "\n");
            if (this.structures.size() > 0) myWriter.write(String.join("\n", this.structures) + "\n");
            myWriter.write(String.format("%-6s%5d\n", "MODEL", 1));

            for (Atom a : this.atoms) {
               // write all atom information, with correct format
               myWriter.write(String.format("%-6s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%-2s", "ATOM", a.getPDBid(), a.getName(), "", a.getResidueName(),
                       a.getChainCode(), a.getResidueID(), "", a.getX(), a.getY(), a.getZ(), 1.0, 0.0, a.getChainCode(), a.getTypeName(), ""));
               myWriter.write("\n");
            }

            // end model
            myWriter.write(String.format("%-6s%5d  ", "TER", this.atoms.get(this.atoms.size() - 1).getPDBid() + 1));
            myWriter.write(String.format("%-4s", ""));
            myWriter.write(String.format("%-3s\n", this.atoms.get(this.atoms.size() - 1).getResidueName()));

            // write master and end lines
            if (this.master != null) myWriter.write(this.master + "\n");
            myWriter.write("ENDMDL");
         }
         else{
            //mdjeep
            // write distance data
            int n = this.atoms.size();
            for (int i = 0; i < n; i++)
            {
               for (int j = i + 1; j < n; j++)
               {
                  Atom A = this.atoms.get(i);
                  Atom B = this.atoms.get(j);

                  myWriter.write(String.format("%5s",i + 1)); //1-indexed
                  myWriter.write(String.format("%5s",j + 1)); //1-indexed
                  myWriter.write(" ");
                  myWriter.write(String.format("%5s",A.getResidueID()));
                  myWriter.write(String.format("%5s",B.getResidueID()));
                  myWriter.write(" ");

                  double distance = A.distanceTo(B);
                  myWriter.write(String.format("%21.16f",distance));
                  myWriter.write(String.format("%21.16f",distance));
                  myWriter.write("  ");
                  myWriter.write(String.format("%-5s",A.getName()));
                  myWriter.write(String.format("%-5s",B.getName()));
                  myWriter.write(" ");
                  myWriter.write(String.format("%-4s",A.getResidueName()));
                  myWriter.write(" ");
                  myWriter.write(String.format("%-4s",B.getResidueName()));
                  myWriter.write("\n");
               }
            }

         }
         // ending
         myWriter.close();
      }
      catch (Exception e) {
         e.printStackTrace();
         System.exit(1);
      }
   }
   /**
    * Gets the number of PDB models in a PDB file
    * @param filename      the PDB file name
    */
   public static int numberOfModels(String filename, String folder)
   {
      int models = 0;
      try
      {
         if (filename == null) throw new IllegalArgumentException("The String supposing to contain the file name is null");
         if (filename.length() == 0) throw new IllegalArgumentException("The String supposing to contain the file name is empty");

         // reading the file
         File input = new File(folder + filename);
         Scanner scan = new Scanner(input);
         while (scan.hasNext()) {
            // the current line
            String line = scan.nextLine();

            // split on multiple whitespaces
            String[] parts = line.split("\\s+");

            if (parts[0].equals("MODEL")) {
               models++;
            }
         }

         //When there is only one model, the MODEL tag never shows up in the file.
         if(models == 0)
            return 1;
      }
      catch (Exception e)
      {
         System.out.println("Could not read file at " + "PDB_files/" + filename);
         e.printStackTrace();
         System.exit(1);
      }
      return models;
   }
}

