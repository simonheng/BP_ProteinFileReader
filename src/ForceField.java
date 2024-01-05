/**
 * {@link ForceField} is a class that reads in force field data and provide functions that output data
 * NOTE: this class will also read files to compute bond angles based on Ramachandran plots!
 * @author  Antonio Mucherino
 * @author  Simon Hengeveld
 * @since   November 25th, 2022
 * @version December, 2023
 * @package ProteinFileReader
 */

import java.io.File;
import java.util.*;

public class ForceField {
    /**
     * The force field we are reading from
     */
    String source;

    /**
     * This maps a pair <AtomName, AminoAcidName> to a forcefield name
     */
    Map<Pair<String, String>, String> ffNames;

    /**
     * This maps a pair of ffnames to an expected bond length
     */
    Map<Pair<String, String>, Double> bonds;

    /**
     * This maps a triplet of ffnames to an expected bond angle
     */
    Map<Pair<String, Pair<String, String>>, Double> bondAngles;

    /**
     * This maps a quadruplet of ffnames to an expected improper angle
     */
    Map<Pair<Pair<String, String>, Pair<String, String>>, Double> impropers;

    private Map<Pair<Integer, Integer>, String> RamachandranIdentifiers;

    /**
     * A map that contains bond angle values based on Ramachandran region.
     * Key = Pair<a,b>, where a is ramachandran identifier (e.g. p) and b is angle identifier (e.g. CACO, for the angle CA-C-O)
     */
    private Map<Pair<String,String>, Double> RamachandranAngles;

    /**
     * Create a {@link ForceField} object from a force field data file
     * @param source A string indicating which force field is desired. Supported force-fields: amber and charmm
     * @throws IllegalArgumentException The String source must correspond to a supported force-field name
     */
    public ForceField(String source){
        try {
            this.source = source;
            if (!source.equalsIgnoreCase("amber") && !source.equalsIgnoreCase("charmm"))
                throw new IllegalArgumentException("The force-field name: " + source + " is not supported (we only support amber and charmm)");

            //Read the information necessary to translate from atom name + residue code to force-field name
            //We expect the following format for the .ffnames file:
            //RESI <residueCode>
            //ATOM <AtomName> <ffName>
            //ATOM <AtomName> <ffName>
            //RESI <residueCode>
            //ATOM <AtomName> <ffName>
            //ATOM <AtomName> <ffName>
            //etc...
            File input = new File("forcefields/" + this.source + ".ffnames");
            Scanner scan = new Scanner(input);
            String residueCode = null;
            this.ffNames = new HashMap<>();
            List<String> residues = new ArrayList<>();
            while (scan.hasNext()) {
                // the current line
                String line = scan.nextLine();
                // split on multiple whitespaces
                String[] parts = line.split("\\s+");
                if (parts[0].equals("RESI")) {
                    residueCode = parts[1].equals("HSD") ? "HIS" : parts[1]; //CHARM calls HISTIDINE as HSP, HSD and HSE, the general code is HIS
                    residues.add(residueCode);
                }
                else if (parts[0].equals("ATOM")) {
                    this.ffNames.put(new Pair<>(parts[1], residueCode), parts[2]);
                }
            }

            input = new File("forcefields/" + this.source + ".dat");
            scan = new Scanner(input);
            this.bonds = new HashMap<>();
            this.bondAngles = new HashMap<>();
            boolean doingBonds = false;
            boolean doingAngles = false;
            //Now that we know the force-field names we can read the parameter files!
            while (scan.hasNext()) {
                String line = scan.nextLine().replaceAll("\\*", "star");
                String[] parts = line.split("\\s*(-|\\s)\\s*");
                if (parts.length > 0 && !parts[0].isEmpty()) {
                    if (parts[0].contains("BOND")) {
                        doingBonds = true;
                    } else if (parts[0].contains("ANGLE")) {
                        doingAngles = true;
                        doingBonds = false;
                    } else if (parts[0].equals("EXTRA") || parts[0].equals("DIHEDRALS")){
                        doingAngles = false;
                    } else if (doingBonds && this.ffNames.containsValue(parts[0])) {
                        this.bonds.put(new Pair<>(parts[0], parts[1]), Double.parseDouble(parts[3]));
                        //parts[2] includes "The harmonic force constants for the angle"
                    } else if (doingAngles && this.ffNames.containsValue(parts[0])) {
                        this.bondAngles.put(new Pair<>(parts[0], new Pair<>(parts[1], parts[2])), Double.parseDouble(parts[4]));
                        //parts[5] includes "The harmonic force constants for the angle"
                    }
                }
            }

            //Missing mappings specific to amino acids...
            if(this.source.equals("charmm")) {
                for(String residue : residues){
                    this.ffNames.put(new Pair<>("H", residue), "H");
                    this.ffNames.put(new Pair<>("H1", residue), "H");
                    this.ffNames.put(new Pair<>("H2", residue), "H");
                    this.ffNames.put(new Pair<>("H3", residue), "H");
                    this.ffNames.put(new Pair<>("H3", residue), "H");
                    this.ffNames.put(new Pair<>("OXT", residue), "OC");

                    if(!this.ffNames.containsKey(new Pair<>("HB3", residue)))
                        this.ffNames.put(new Pair<>("HB3", residue), "HA2"); //named HB1 in charmm
                    this.ffNames.put(new Pair<>("HG3", residue), "HA2"); //named HG1 in charmm
                    this.ffNames.put(new Pair<>("HD3", residue), "HA2"); //named HD1 in charmm

                    if(!this.ffNames.containsKey(new Pair<>("HE3", residue)))
                        this.ffNames.put(new Pair<>("HE3", residue), "HA2"); //named HE1 in charmm
                    this.ffNames.put(new Pair<>("N", residue), "NH2"); //we rename all NH1 to NH2 for consistency
                }
                this.ffNames.put(new Pair<>("HG", "SER"), "H"); //named HG1 instead of HG in charmm
                this.ffNames.put(new Pair<>("CD1", "ILE"), "CT3"); //named CD instead of CD1 in charmm
                this.ffNames.put(new Pair<>("HG13", "ILE"), "HA2"); //missing fom ILE in charmm
                this.ffNames.put(new Pair<>("HD13", "ILE"), "HA3"); //missing fom ILE in charmm
                this.ffNames.put(new Pair<>("HD11", "ILE"), "HA3"); //missing fom ILE in charmm
                this.ffNames.put(new Pair<>("HD12", "ILE"), "HA3"); //missing fom ILE in charmm
                this.ffNames.put(new Pair<>("HA3", "GLY"), "HB2"); //named HA1 in charmm
                this.ffNames.put(new Pair<>("HA", "GLY"), "HB2"); //named HA1 in charmm

                //NH1 and NH2 are the same in charmm
                List<Pair<String,String>> bondKeys = new ArrayList<>(this.bonds.keySet());
                for(Pair<String, String> ffnames : bondKeys){
                    if(ffnames._1().equals("NH1")) this.bonds.put(new Pair<>("NH2", ffnames._2()), this.bonds.get(ffnames));
                    else if(ffnames._2().equals("NH1")) this.bonds.put(new Pair<>(ffnames._1(),"NH2"), this.bonds.get(ffnames));

                    if(ffnames._1().equals("N")) this.bonds.put(new Pair<>("NH2", ffnames._2()), this.bonds.get(ffnames));
                    else if(ffnames._2().equals("N")) this.bonds.put(new Pair<>(ffnames._1(),"NH2"), this.bonds.get(ffnames));
                }
                this.bonds.put(new Pair<>("OC", "C"), this.bonds.get(new Pair<>("OC", "CC")));
                this.bonds.put(new Pair<>("CP1", "NH2"), this.bonds.get(new Pair<>("N", "CP1")));
                this.bonds.put(new Pair<>("CP3", "NH2"), this.bonds.get(new Pair<>("N", "CP3")));

                List<Pair<String, Pair<String,String>>> bondAngleKeys = new ArrayList<>(this.bondAngles.keySet());
                for(Pair<String, Pair<String,String>> ffnames : bondAngleKeys){
                    if(ffnames._1().equals("NH1")) this.bondAngles.put(new Pair<>("NH2", ffnames._2()), this.bondAngles.get(ffnames));
                    else if(ffnames._2()._2().equals("NH1")) this.bondAngles.put(new Pair<>(ffnames._1(),new Pair<>(ffnames._2()._1(), "NH2")), this.bondAngles.get(ffnames));
                    else if(ffnames._2()._1().equals("NH1")) this.bondAngles.put(new Pair<>(ffnames._1(),new Pair<>("NH2", ffnames._2()._2())), this.bondAngles.get(ffnames));

                    if(ffnames._1().equals("N")) this.bondAngles.put(new Pair<>("NH2", ffnames._2()), this.bondAngles.get(ffnames));
                    else if(ffnames._2()._2().equals("N")) this.bondAngles.put(new Pair<>(ffnames._1(),new Pair<>(ffnames._2()._1(), "NH2")), this.bondAngles.get(ffnames));
                    else if(ffnames._2()._1().equals("N")) this.bondAngles.put(new Pair<>(ffnames._1(),new Pair<>("NH2", ffnames._2()._2())), this.bondAngles.get(ffnames));

                    if(ffnames._1().equals("OC") && ffnames._2()._1().equals("CC"))
                        this.bondAngles.put(new Pair<>("OC", new Pair<>("C", ffnames._2()._2())),this.bondAngles.get(ffnames));
                }
                this.bondAngles.put(new Pair<>("H", new Pair<>("NH2","CP1")), this.bondAngles.get(new Pair<>("HC", new Pair<>("NP","CP1"))));
                this.bondAngles.put(new Pair<>("H", new Pair<>("NH2","CP3")), this.bondAngles.get(new Pair<>("HC", new Pair<>("NP","CP3"))));
                this.bondAngles.put(new Pair<>("OC", new Pair<>("C","CT1")), this.bondAngles.get(new Pair<>("OC", new Pair<>("CC","CT1"))));
                this.bondAngles.put(new Pair<>("OC", new Pair<>("C","O")), this.bondAngles.get(new Pair<>("OC", new Pair<>("CC","OC"))));
            }
            else if(this.source.equals("amber")) {
                for (String residue : residues) {
                    this.ffNames.put(new Pair<>("H", residue), "H");
                    this.ffNames.put(new Pair<>("HA", residue), "H1");
                }
            }

            //Read Ramachandran identifiers
            RamachandranIdentifiers = new HashMap<>();
            input = new File("ramachandran/Ramachandran_matrix.txt");
            scan = new Scanner(input);

            int psi = 180;
            while (scan.hasNext()) {
                String line = scan.nextLine();
                String[] parts = line.split(" ");
                int phi = -180;
                for(String part : parts){
                    if(part.equals("-"))
                        part = "null";
                    RamachandranIdentifiers.put(new Pair<>(phi,psi), part);
                    phi = phi + 5;
                }
                psi = psi - 5;
            }

            //Read Ramachandran values
            RamachandranAngles = new HashMap<>();
            String[] angles = new String[]{"CACNp","CACO","CmNCA","NCAC","OCNp"};
            for(String angle : angles){
                input = new File("ramachandran/averaged_"+ angle + ".dat");
                scan = new Scanner(input);
                while (scan.hasNext()) {
                    String line = scan.nextLine();
                    String[] parts = line.split(" ");
                    RamachandranAngles.put(new Pair<>(parts[0],angle), Double.parseDouble(parts[2]));
                }
            }
        }
        catch (Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Finds the forcefield name of an atom
     * @param a the atom
     * @return the name of the atom within the current force-field
     * @throws IllegalStateException Cannot find a force-field name for this atom
     */
    public String ffName(Atom a) {
        //We remove numbers from the amino acid
        Pair<String, String> nameResPair = new Pair<String,String>(a.getName(), a.getResidueName().replaceAll("\\d",""));
        try {

            if(!this.ffNames.containsKey(nameResPair))
                throw new IllegalStateException("Cannot find a force-field name for this atom: " + nameResPair);
            else
                return this.ffNames.get(nameResPair);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return null;
    }

    /**
     * Looks up the expected bond distance between two atoms
     * @param a the first atom
     * @param b the second atom
     * @return The expected bond distance between a and b
     * @throws IllegalStateException We do not have information about a bond between a and b
     */
    public double expectedBondDistance(Atom a, Atom b) {
        Pair<String, String> ffPair = new Pair<>(this.ffName(a), this.ffName(b));
        Pair<String, String> ffPair2 = new Pair<>(this.ffName(b), this.ffName(a));
        try {
            if (this.bonds.containsKey(ffPair))
                return this.bonds.get(ffPair);
            else if(this.bonds.containsKey(ffPair2))
                return this.bonds.get(ffPair2);
            else
                throw new IllegalStateException("We do not have information about a bond between: " + a + "(" + ffPair._1() + ") and " + b + "(" + ffPair._2() + ")");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return -1; //will never reach here
    }

    /**
     * Looks up the expected bond angle distance between three atoms
     * @param a the first atom
     * @param b the second atom
     * @param c the third atom
     * @return The expected bond angle derived distance between a, b and c, in degrees!
     * @throws IllegalStateException We do not have information about a bond between a and b
     */
    public double getBondAngle(Atom a,Atom b,Atom c){
        String ffA = ffName(a);
        String ffB = ffName(b);
        String ffC = ffName(c);
        Pair<String, Pair<String,String>> ffPair = new Pair<>(ffA, new Pair<>(this.ffName(b), this.ffName(c)));
        Pair<String, Pair<String,String>> ffPair2 = new Pair<>(ffC, new Pair<>(this.ffName(b), this.ffName(a)));
        try {
            if (this.bondAngles.containsKey(ffPair))
                return this.bondAngles.get(ffPair);
            else if(this.bondAngles.containsKey(ffPair2))
                return this.bondAngles.get(ffPair2);
            else
                throw new Exception("We do not have information about an angle between: " + a + "(" + ffA + "), " + b + "(" + ffB + ") and " + c + "(" + ffC + ")");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return -1; //will never reach here
    }

    /**
     * Looks up the expected bond angle distance between three atoms
     * @param a the first atom
     * @param b the second atom
     * @param c the third atom
     * @return The expected bond angle derived distance between a, b and c
     */
    public double expectedBondAngleDistance(Atom a,Atom b,Atom c)
    {
        double ab = this.expectedBondDistance(a, b);
        double bc = this.expectedBondDistance(b, c);
        double abc = this.getBondAngle(a, b, c);
        return Geometry.solveSAS(ab, bc, Math.toRadians(abc));
    }

    /**
     * Looks up the expected (proper) torsion angle distance between four atoms, given a torsion angle
     * @param a the first atom
     * @param b the second atom
     * @param c the third atom
     * @return The expected bond angle derived distance between a, b and c
     */
    public double expectedTorsionAngleDistance(Atom a,Atom b,Atom c, Atom d, double abcdDeg){
        //compute the necessary angles and distances
        double abcd = Math.toRadians(abcdDeg);
        double ab = this.expectedBondDistance(a, b);
        double bc = this.expectedBondDistance(b, c);
        double cd = this.expectedBondDistance(c, d);
        double abc = Math.toRadians(this.getBondAngle(a,b,c));
        double bcd = Math.toRadians(this.getBondAngle(b,c,d));
        return Geometry.dihedralDistance(ab, bc, cd, abc, bcd, abcd);
    }

    public double vanDerWaalsDistance(Atom a, Atom b) {
        return Atom.vanDerWaalsRadius.get(a.getType()) + Atom.vanDerWaalsRadius.get(b.getType());
    }

    /**
     * Gives the RamachandranAngleIdentifier, as well as the which residueID is the starter for the phi/psi pair
     * @param a The first atom
     * @param b The second atom
     * @param c The third atom
     * @return null if not found, the angle in radians otherwise
     */
    public Pair<String, Integer> ramachandranAngleIdentifer(Atom a, Atom b, Atom c){
        //thet starting resID for the phi/psi pair
        int startingRes = b.getResidueID();
        String id = null;
        //Check if we are {"CACNp","CACO","CmNCA","NCAC","OCNp"};
        if(a.getResidueID() == b.getResidueID() && b.getResidueID() == c.getResidueID()){
            //all in same residue: can only be NCAC or CACO
            if(b.getName().equals("C")){
                if(a.getName().equals("CA") && c.getName().equals("O"))
                    id = "CACO";
                if(a.getName().equals("O") && c.getName().equals("CA"))
                    id = "CACO";
            }
            else if(b.getName().equals("CA")){
                if(a.getName().equals("N") && c.getName().equals("C"))
                    id = "NCAC";
                if(a.getName().equals("C") && c.getName().equals("N"))
                    id = "NCAC";
            }
            else return null;
        }
        else if(a.getResidueID() == b.getResidueID() && c.getResidueID() == b.getResidueID() + 1) {
            if(a.getName().equals("O") && b.getName().equals("C") & c.getName().equals("N"))
                id = "OCNp";
            if(a.getName().equals("CA") && b.getName().equals("C") & c.getName().equals("N"))
                id = "CACNp";
        }
        else if(c.getResidueID() == b.getResidueID() && a.getResidueID() == b.getResidueID() + 1) {
            if(c.getName().equals("O") && b.getName().equals("C") & a.getName().equals("N"))
                id = "OCNp";
            if(c.getName().equals("CA") && b.getName().equals("C") & a.getName().equals("N"))
                id = "CACNp";
        }
        else if(c.getResidueID() == b.getResidueID() && a.getResidueID() == b.getResidueID() - 1) {
            //the phi/psi pair is this + next
            startingRes = a.getResidueID();
            if(a.getName().equals("C") && b.getName().equals("N") && c.getName().equals("CA"))
                id = "CmNCA";
        }
        else if(a.getResidueID() == b.getResidueID() && c.getResidueID() == b.getResidueID() - 1) {
            startingRes = c.getResidueID();
            if(a.getName().equals("CA") && b.getName().equals("N") && c.getName().equals("C"))
                id = "CmNCA";
        }

        if (id == null) return null;
        return new Pair<>(id, startingRes);
    }

    /**
     * Returns the angle according to the Ramachadran plot in radians (!)
     * @param angleIdentifier the identfier of the angle
     * @param phi the phi angle
     * @param psi the psi angle
     * @return the bond angle
     */
    public double ramachandranAngle(String angleIdentifier, double phi, double psi){
        double phiDeg = Math.toDegrees(phi);
        double psiDeg = Math.toDegrees(psi);
        int phiR = -180;
        int psiR = -180;

        while(phiR <= 180){
            if(Math.abs(phiDeg - phiR) <= 2.5)
                break;
            phiR += 5;
        }
        while(psiR <= 180){
            if(Math.abs(psiDeg - psiR) <= 2.5)
                break;
            psiR += 5;
        }

        String ramachandranID = this.RamachandranIdentifiers.get(new Pair<>(phiR,psiR));
        return Math.toRadians(this.RamachandranAngles.get(new Pair<>(ramachandranID, angleIdentifier)));
    }
}
