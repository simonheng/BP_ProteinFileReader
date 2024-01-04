/**
 * {@link NMR} is a class that can read NMR files, storing all the relevant distance and angle information.
 * @author  Antonio Mucherino
 * @author  Simon Hengeveld
 * @since   June 7, 2021
 * @version December 2023
 * @package ProteinFileReader
 */

import java.io.File;
import java.util.*;

/**
 * Subclass used to group the NMR restraint information!
 */
class NMRDistance{
    /**
     * The name of the first atom
     */
    public String atomName1;

    /**
     * The name of the second atom
     */
    public String atomName2;

    /**
     * The residue ID of the first atom in the NMR file
     */
    public int seqID1;

    /**
     * The residue ID of the second atom in the NMR file
     */
    public int seqID2;

    /**
     * The lowerbound of the restraint
     */
    public double minDistance;

    /**
     * The upperbound of the restraint
     */
    public double maxDistance;

    /**
     * A third "distance" value in case the precise value is known.
     */
    public double distanceVal;

    /**
     * Check if one NMR distance equals another (used to attempt to recognize uncertain restraints!)
     * @param other the other NMR distance
     * @return true if they are the same restraint, false otherwise
     */
    public boolean equals(NMRDistance other){
        //check atomNames
        boolean uncertainAtomName = false;
        String[] uncertainty = new String[]{"#", "HN", "M", "Q"};
        for(String tag : uncertainty){
            if(atomName1.contains(tag) || atomName2.contains(tag) || other.atomName1.contains(tag) || other.atomName2.contains(tag))
            {
                uncertainAtomName = true;
                break;
            }
        }

        //If we have an uncertain atom name it is impossible to find a matching pair...
        if(!uncertainAtomName) {
            if(!atomName1.equals(other.atomName1) || !atomName2.equals(other.atomName2)){
                //maybe they are switched up?
                if(!atomName1.equals(other.atomName2) && !atomName2.equals(other.atomName1))
                    return false;
                else {
                    if (seqID1 != other.seqID2 || seqID2 != other.seqID1) return false;
                }
            }
            else
            {
                if(seqID1 != other.seqID1 || seqID2 != other.seqID2)
                {
                    //it could be that atomName1 == atomName2
                    if(atomName1.equals(atomName2)){
                        if (seqID1 != other.seqID2 || seqID2 != other.seqID1) return false;
                    }
                    else
                        return false;
                }
            }
        }

        double epsilon = 0.0001;
        //make sure to use an epsilon for floating imprecision


        if(!Maths.epsEqual(minDistance, other.minDistance, epsilon) ) return false;
        if(!Maths.epsEqual(maxDistance, other.maxDistance, epsilon)) return false;
        if(!Maths.epsEqual(distanceVal, other.distanceVal, epsilon) ) return false;
        return true;
    }
}

/**
 * Subclass used to group the NMR torsion angle information (usually predicted using TALOS+)
 */
class NMRAngle{
    /**
     * The name of the first atom of the quadruplet
     */
    public String atomName1;

    /**
     * The name of the second atom of the quadruplet
     */
    public String atomName2;

    /**
     * The name of the third atom of the quadruplet
     */
    public String atomName3;

    /**
     * The name of the fourth atom of the quadruplet
     */
    public String atomName4;

    /**
     * The residueID of the first atom of the quadruplet
     */
    public int seqID1;

    /**
     * The residueID of the second atom of the quadruplet
     */
    public int seqID2;

    /**
     * The residueID of the third atom of the quadruplet
     */
    public int seqID3;

    /**
     * The residueID of the fourth atom of the quadruplet
     */
    public int seqID4;

    /**
     * The lowerbound of the torsion angle restraint (in degrees!)
     */
    public double minAngle;

    /**
     * The upperbound of the torsion angle restraint (in degrees!)
     */
    public double maxAngle;

    /**
     * Checks if this torsion angle restraints equals another
     * @param other the other torsion angle restraint
     * @return true if they are equal, false otherwise
     */
    public boolean equals(NMRAngle other){
        if(!atomName1.equals(other.atomName1)) return false;
        if(!atomName2.equals(other.atomName2)) return false;
        if(!atomName3.equals(other.atomName3)) return false;
        if(!atomName4.equals(other.atomName4)) return false;
        if(seqID1 != other.seqID1) return false;
        if(seqID2 != other.seqID2) return false;
        if(seqID3 != other.seqID3) return false;
        if(seqID4 != other.seqID4) return false;

        double epsilon = 0.0001;
        //make sure to use an epsilon for floating imprecision
        if(Math.abs(minAngle - other.minAngle) > epsilon ) return false;
        if(Math.abs(maxAngle - other.maxAngle) > epsilon ) return false;
        return true;
    }
}


/**
 * Class that will read files that hold NMR data.
 */
public class NMR {
    /**
     * The chain of amino acids, as read from the NMR file.
     */
    private List<String> residueChain;

    /**
     * The list of NMR restraints
     */
    private List<NMRDistance> NMRdistances;

    /**
     * The list of NMR (TALOS+) torsion angle restraints
     */
    private List<NMRAngle> NMRangles;

    /**
     * The residue chain read from the NMR file
     * @return the residue chain
     */
    public List<String> getResidueChain() {
        return residueChain;
    }

    /**
     * The accessor for the NMR distance list
     * @return the list of NMR restraints
     */
    public List<NMRDistance> getNMRdistances() {
        return NMRdistances;
    }

    /**
     * The accessor for the NMR torsion angle list
     * @return the list of NMR (TALOS+) torsion angle restraints
     */
    public List<NMRAngle> getNMRangles() {
        return NMRangles;
    }

    /**
     * Constructs the NMR object by reading from an NMR file (STAR or .mr)
     * @param filename      the NMR file location
     */
    public NMR(String filename) {
        String fileExt = filename.substring(filename.lastIndexOf(".") + 1);
        boolean isStarFile = fileExt.equals("str");

        NMRdistances = new ArrayList<>();
        NMRangles = new ArrayList<>();
        residueChain = new ArrayList<>();

        try {
            File input = new File("NMR_files/" + filename);
            Scanner scan = new Scanner(input);

            if (isStarFile) {
                /*
                 *  A STAR file contains different "saves"\
                 *  For every save, we are given a list of identifiers which map to a table, such as torsion angles or distance constraints
                 *  Such a table then contains an array of maps, with string identifiers. If there are n observed distances, the length of this array is n
                 */
                ArrayList<HashMap<String, String>> chainEntries = new ArrayList<>();
                ArrayList<HashMap<String, String>> distanceEntries = new ArrayList<>();
                ArrayList<HashMap<String, String>> angleEntries = new ArrayList<>();
                boolean loopMode = false;
                boolean tableFound = false;
                ArrayList<HashMap<String, String>> currentTable = distanceEntries;
                List<String> identifiers = new ArrayList<>();

                while (scan.hasNext()){
                    //split on multiple whitespaces
                    String line = scan.nextLine().trim();
                    if (line.isEmpty())
                        continue;
                    String[] parts = line.split("\\s+");

                    if (loopMode) {
                        if (parts.length == 1) {
                            //We are adding identifiers
                            if (!identifiers.contains(parts[0]))
                                identifiers.add(parts[0]);
                        } else {
                            //We only care about torsion angles, distances and the amino acid chain
                            if (!identifiers.get(0).equals("_Entity_comp_index.ID")
                                && !identifiers.get(0).equals("_Gen_dist_constraint.ID")
                                && !identifiers.get(0).equals("_Torsion_angle_constraint.ID"))
                            {
                                loopMode = false;
                                continue;
                            }

                            //We are adding an entry
                            HashMap<String, String> entry = new HashMap<>();
                            for (int j = 0; j < parts.length; j++)
                                entry.put(identifiers.get(j), parts[j]);
                            currentTable.add(entry);
                        }
                    }

                    if (parts[0].equals("_Entity.Sf_category")) {
                        currentTable = chainEntries;
                        tableFound = true;
                    } else if (parts[0].contains("distance_constraints")) {
                        currentTable = distanceEntries;
                        tableFound = true;
                    } else if (parts[0].contains("dihedral")) {
                        currentTable = angleEntries;
                        tableFound = true;
                    } else if (parts[0].equals("loop_")) {
                        if (tableFound) {
                            loopMode = true;
                            identifiers.clear();
                        }
                    } else if (parts[0].equals("stop_"))
                        loopMode = false;
                }

                if(distanceEntries.isEmpty())
                    throw new IllegalStateException("No distance restraints found!");
                if(angleEntries.isEmpty())
                    System.out.println("We have no torsion angle information in the STAR FILE!");
                    //throw new IllegalStateException("No torsion angle restraints found!");

                //Get the Distance info from the hashtables
                for (HashMap<String, String> entry : distanceEntries){
                    NMRDistance dist = new NMRDistance();
                    dist.atomName1 = entry.get("_Gen_dist_constraint.Atom_ID_1");
                    dist.atomName2 = entry.get("_Gen_dist_constraint.Atom_ID_2");
                    dist.seqID1 = Integer.parseInt(entry.get("_Gen_dist_constraint.Seq_ID_1"));
                    dist.seqID2 = Integer.parseInt(entry.get("_Gen_dist_constraint.Seq_ID_2"));
                    dist.minDistance = Double.parseDouble(entry.get("_Gen_dist_constraint.Distance_lower_bound_val"));
                    dist.maxDistance = Double.parseDouble(entry.get("_Gen_dist_constraint.Distance_upper_bound_val"));
                    dist.distanceVal = Double.parseDouble(entry.get("_Gen_dist_constraint.Distance_val"));
                    if(dist.atomName1.isEmpty() || dist.atomName2.isEmpty())
                        throw new IllegalStateException("Two atoms for the NMR distance should be defined, please check the lines pertaining to the distance restraint information!");
                    NMRdistances.add(dist);
                }

                //Get the Angle  info from the hashtables
                for (HashMap<String, String> entry : angleEntries) {
                    NMRAngle angle = new NMRAngle();
                    angle.atomName1 = entry.get("_Torsion_angle_constraint.Atom_ID_1");
                    angle.atomName2 = entry.get("_Torsion_angle_constraint.Atom_ID_2");
                    angle.atomName3 = entry.get("_Torsion_angle_constraint.Atom_ID_3");
                    angle.atomName4 = entry.get("_Torsion_angle_constraint.Atom_ID_4");
                    angle.seqID1 = Integer.parseInt(entry.get("_Torsion_angle_constraint.Seq_ID_1"));
                    angle.seqID2 = Integer.parseInt(entry.get("_Torsion_angle_constraint.Seq_ID_2"));
                    angle.seqID3 = Integer.parseInt(entry.get("_Torsion_angle_constraint.Seq_ID_3"));
                    angle.seqID4 = Integer.parseInt(entry.get("_Torsion_angle_constraint.Seq_ID_4"));
                    angle.minAngle = Double.parseDouble(entry.get("_Torsion_angle_constraint.Angle_lower_bound_val"));
                    angle.maxAngle = Double.parseDouble(entry.get("_Torsion_angle_constraint.Angle_upper_bound_val"));

                    if(angle.atomName1.isEmpty() || angle.atomName2.isEmpty() || angle.atomName3.isEmpty() || angle.atomName4.isEmpty())
                        throw new IllegalStateException("All of the four atoms for the dihedral angles should be defined, please check the lines pertaining to the NMR torsion angle information!");
                    NMRangles.add(angle);


                }

                //Get the amino acid chain from the read file
                for (HashMap<String, String> chainEntry : chainEntries) {
                    String letterCode = chainEntry.get("_Entity_comp_index.Comp_ID"); //the three letter residue code
                    residueChain.add(letterCode);
                }
            } else if (fileExt.equals("res")) {
                while (scan.hasNext()){
                    //split on multiple whitespaces
                    String line = scan.nextLine().trim();
                    if (line.isEmpty())
                        continue;

                    line = line.replaceAll("(\\(|\\))", ""); //remove parenthesis

                    String[] parts = line.split("\\s+");

                    if (parts.length == 25) {
                        NMRAngle currentAngle = new NMRAngle();
                        currentAngle.atomName1 = parts[5];
                        currentAngle.atomName2 = parts[10];
                        currentAngle.atomName3 = parts[15];
                        currentAngle.atomName4 = parts[20];

                        currentAngle.seqID1 = Integer.parseInt(parts[2]);
                        currentAngle.seqID2 = Integer.parseInt(parts[7]);
                        currentAngle.seqID3 = Integer.parseInt(parts[12]);
                        currentAngle.seqID4 = Integer.parseInt(parts[17]);

                        double target = Double.parseDouble(parts[22]);
                        double delta = Double.parseDouble(parts[23]);
                        currentAngle.minAngle = target - delta;
                        currentAngle.maxAngle = target + delta;

                        NMRangles.add(currentAngle);
                    }
                }
            } else {
                //Parse MR file
                boolean doingAngleRestraint = false;
                NMRAngle currentAngle = new NMRAngle();
                while (scan.hasNext()){
                    //split on multiple whitespaces
                    String line = scan.nextLine().trim();
                    line = line.replaceAll("(\\(|\\))", ""); //remove parenthesis
                    if (line.isEmpty())
                        continue;
                    String[] parts = line.split("\\s+");

                    //We are doing a restraint
                    if (parts[0].equals("assign")) {
                        if (parts.length == 11) {
                            //torsion angle restraint
                            doingAngleRestraint = true;
                            currentAngle = new NMRAngle();
                            currentAngle.seqID1 = Integer.parseInt(parts[2]);
                            currentAngle.atomName1 = parts[5];
                            currentAngle.seqID2 = Integer.parseInt(parts[7]);
                            currentAngle.atomName2 = parts[10];
                        }
                        if (parts.length == 14 || parts.length == 16) {
                            //distance restraint
                            NMRDistance dist = new NMRDistance();
                            dist.seqID1 = Integer.parseInt(parts[2]);
                            dist.atomName1 = parts[5];
                            dist.seqID2 = Integer.parseInt(parts[7]);
                            dist.atomName2 = parts[10];

                            dist.distanceVal = Double.parseDouble(parts[11]);
                            double minus = Double.parseDouble(parts[12]);
                            double plus = Double.parseDouble(parts[13]);
                            dist.minDistance = dist.distanceVal - minus;
                            dist.maxDistance = dist.distanceVal + plus;
                            NMRdistances.add(dist);
                        }
                    } else if (doingAngleRestraint) {
                        if(parts.length < 12)
                            throw new IllegalStateException("The second half of the NMR distance constraint should at least have 12 parts!");
                        //second half of angle restraint
                        currentAngle.seqID3 = Integer.parseInt(parts[1]);
                        currentAngle.atomName3 = parts[4];
                        currentAngle.seqID4 = Integer.parseInt(parts[6]);
                        currentAngle.atomName4 = parts[9];

                        double angleMiddle = Double.parseDouble(parts[11]);
                        double angleDev = Double.parseDouble(parts[12]);

                        currentAngle.minAngle = angleMiddle - angleDev;
                        currentAngle.maxAngle = angleMiddle + angleDev;
                        NMRangles.add(currentAngle);
                        doingAngleRestraint = false;
                    }
                }
            }

            if(NMRangles.isEmpty() && NMRdistances.isEmpty())
                throw new IllegalStateException("No NMR information could be read, please check the NMR file!");
        }
        catch (Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
}

