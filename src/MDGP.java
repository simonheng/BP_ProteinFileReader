/**
 * {@link MDGP} is the class that constructs the molecular DGP instances, using several constructors which will use various data sources
 * @author  Antonio Mucherino
 * @author  Simon Hengeveld
 * @since   November 25th, 2022
 * @version January, 2024
 * @see     DGP
 * @package ProteinFileReader
 */

import java.io.File;
import java.io.FileWriter;
import java.util.*;

public class MDGP extends DGP<Atom> {
    /**
     * TEMPORARY: data-structure for dihedral angles
     */
    public Map<Atom, Pair<List<Atom>, double[]>> dihedralAngles;

    /**
     * Label for distances derived from bonds
     */
    final static String bondDistance = "bond";

    /**
     * Label for distances derived from bond angles
     */
    final static String angleDistance = "angle";

    /**
     * Label for distances derived from a proper torsion angle
     */
    final static String torsionDistance = "angle";

    /**
     * Label for distances derived from NMR distance restraints
     */
    final static String NMRDistance = "NMR_dist_entry";

    /**
     * Label for distances derived from NMR torsion angle restraints
     */
    final static String NMRTorsionDistance = "NMR_torsion_entry";

    /**
     * Label for distances derived from van der Waals radii
     */
    final static String waalsDistance = "waals";

    /**
     * The {@link ForceField} object we use to obtain force-field data!
     */
    private ForceField ff;

    /**
     * Do we include sidechains?
     */
    private boolean sidechain;

    /**
     * Do we include hydrogens?
     */
    private boolean hydrogens;

    /**
     * Do we want only bond distances?
     */
    private boolean onlybonds;

    /**
     * The associated PDB model: can be null
     */
    public PDB model;

    /**
     * The chain char within the PDB model (default = 'n')
     */
    public char chain;

    /**
     * Creates a {@link MDGP} based on a list of amino acids (in order)
     *
     * @param aminoAcids a list of Strings naming the amino acids in the {@link MDGP} in order
     * @param sidechain  whether we want to include sidechains
     * @param hydrogens  whether we want to include hydrogens
     * @param nTerm      denoting whether we want an N-terminus at the beginning of the {@link MDGP}
     * @param cTerm      denoting whether we want a C-terminus at the end of the {@link MDGP}
     * @param onlybonds  denoting whether we want only bond distances or also angle distances
     */
    public MDGP(List<String> aminoAcids, ForceField ff, boolean sidechain, boolean hydrogens, boolean nTerm, boolean cTerm, boolean onlybonds) {
        //initialize dgp at fixed dimension 3
        super(3);
        this.ff = ff;
        this.sidechain = sidechain;
        this.hydrogens = hydrogens;
        this.onlybonds = onlybonds;

        //NOTE: i starts at 1, because PDB  and NMR follow standard that residue IDs start at 1!
        for (int residueID = 1; residueID <= aminoAcids.size(); residueID++) {
            String aminoCode = aminoAcids.get(residueID - 1);

            //add the amino acid to the MDGP with sidechain if desired
            if (aminoCode.equals("NH2"))
                continue;

            this.addBackbone(aminoCode, residueID, nTerm && residueID == 1, cTerm && residueID == aminoAcids.size());
            if (sidechain)
                this.addSideChain(aminoCode, residueID);

            //Create peptide bonds to the previous amino acid
            if (residueID > 1) {
                //Add N to the C atom of the previous amino acid
                Atom N = this.getAtom("N", residueID);
                Atom previousC = this.getAtom("C", residueID - 1);

                this.addAtom(N, previousC);
                this.addAtom(previousC, N);

                for (Atom a : this.starOf(N))
                    if (this.getDistance(a, N).getLabel().equals(MDGP.bondDistance))
                        addAtom(a, N);
                for (Atom a : this.starOf(previousC))
                    if (this.getDistance(a, previousC).getLabel().equals(MDGP.bondDistance))
                        addAtom(a, previousC);
            }

        }

        //Add dihedral angles from force field parameters!
        this.dihedralAngles = new HashMap<>();

        System.out.println("Constructed MDGP object with " + aminoAcids.size() + " amino acids and " + this.numberOfVertices() + " atoms!");
    }

    /**
     * Creates a {@link MDGP} based on an NMR object
     *
     * @param nmr       the NMR object to get the amino acid chain and the distances from (needs sidechain!)
     * @param sidechain whether we want to include sidechains
     * @param hydrogens whether we want to include hydrogens
     * @param nTerm     denoting whether we want an N-terminus at the beginning of the {@link MDGP}
     * @param cTerm     denoting whether we want a C-terminus at the end of the {@link MDGP}
     * @param pdb       Exact distance found in dgp matching a torsion angle distance: this might be an improper dihedral
     */
    public MDGP(NMR nmr, ForceField ff, boolean sidechain, boolean hydrogens, boolean nTerm, boolean cTerm, PDB pdb, char chain) {
        //Construct amino acid chain
        this(nmr.getResidueChain(), ff, sidechain, hydrogens, nTerm, cTerm, false);

        //save chain and PDB model in instance
        this.model = pdb;
        this.chain = chain;

        double totalRanges = 0;
        int count = 0;
        int i = 0;
        //Add NMR information
        for (NMRDistance nmrDist : nmr.getNMRdistances()) {
            i++;
            Atom a1 = this.getAtom(nmrDist.atomName1, nmrDist.seqID1);
            Atom a2 = this.getAtom(nmrDist.atomName2, nmrDist.seqID2);

            if (a1 == null || a2 == null) {
                //System.out.println("Atom with name " + atomName1 + " or " + atomName2 + " not found!");
                continue;
            }

            boolean weHaveBounds = true;

            double expected = -1;

            if (nmrDist.minDistance == nmrDist.maxDistance) {
                expected = nmrDist.minDistance;
                weHaveBounds = false;
            } else {
                double range = nmrDist.maxDistance - nmrDist.minDistance;
                count++;
                totalRanges += range;
            }


            //System.out.println(nmrDist.minDistance + " " + nmrDist.maxDistance);
            Distance<Atom> dist;

            if (!weHaveBounds) {
                dist = new Distance<>(expected, a1, a2);
                if (this.contains(a1, a2))
                    this.removeEdge(a1, a2);
            } else {
                //Create the new distance
                dist = new Distance<>(a1, a2, nmrDist.minDistance, nmrDist.maxDistance);

                //Add the new edge to the dgp (if it does not exist)
                if (this.contains(a1, a2)) {
                    Distance<Atom> existing = this.getDistance(a1, a2);
                    if (existing.hasBounds()) {
                        //intersect with bounds in file
                        //Get existing bounds
                        double elb = existing.getLowerBound();
                        double eub = existing.getUpperBound();

                        if (elb > nmrDist.maxDistance || nmrDist.minDistance > eub) {
                            System.out.println("No overlap :(");

                            //This is problematic, because there is no intersection
                            //Throw away the existing information?
                        } else {
                            //Compute intersection
                            double nlb = Math.max(nmrDist.minDistance, elb);
                            double nub = Math.min(nmrDist.maxDistance, eub);
                            dist = new Distance<>(a1, a2, nlb, nub);
                        }
                    } else if (existing.hasExpectedValue()) {
                        //we do not use the NMR information...
                        dist = existing;
                        //throw new IllegalStateException("Exact distance found in dgp matching a torsion angle distance: this might be an improper dihedral.");
                    }
                    this.removeEdge(a1, a2);
                }
            }
            dist.setLabel(MDGP.NMRDistance);
            this.addEdge(dist);
        }

        System.out.println(count + " NMR restraints added, average restraint range: " + (totalRanges / count));

        totalRanges = 0;
        count = 0;

        //Save the NMR torsion angle based distances in the dgp
        for (NMRAngle nmrAngle : nmr.getNMRangles()) {
            Atom a1 = this.getAtom(nmrAngle.atomName1, nmrAngle.seqID1);
            Atom a2 = this.getAtom(nmrAngle.atomName2, nmrAngle.seqID2);
            Atom a3 = this.getAtom(nmrAngle.atomName3, nmrAngle.seqID3);
            Atom a4 = this.getAtom(nmrAngle.atomName4, nmrAngle.seqID4);

            if (a1 == null || a2 == null || a3 == null || a4 == null)
                continue;

            double angle1 = Math.abs(nmrAngle.minAngle);

            double angle2 = Math.abs(nmrAngle.maxAngle);

            if (angle1 > 180)
                angle1 = 180 - angle1;
            if (angle2 > 180)
                angle2 = 180 - angle2;

            double minAngle = Math.min(angle1, angle2);
            double maxAngle = Math.max(angle1, angle2);

            double range = maxAngle - minAngle;
            totalRanges += range;

            double lb = this.ff.expectedTorsionAngleDistance(a1, a2, a3, a4, minAngle);
            double ub = this.ff.expectedTorsionAngleDistance(a1, a2, a3, a4, maxAngle);

            if (lb > ub) {
                double tmp = lb;
                lb = ub;
                ub = tmp;
            }


            //Because angles can be negative, the angle lower bound does not always correspond to the dihedral distance lowerbound
            Distance<Atom> dist = new Distance<>(a1, a4, lb, ub);

            //Check the pdb file to check if this NMR distance is valid!
            if (pdb != null) {
                Atom pdbA1 = pdb.getAtom(a1.getName(), a1.getResidueID());
                Atom pdbA4 = pdb.getAtom(a4.getName(), a4.getResidueID());

                double d = pdbA1.distanceTo(pdbA4);
                if (lb > d || d > ub)
                    continue; //this is an invalid NMR distance!
            }

            if (lb == ub) {
                //we actually only have an expected value...
                if (this.contains(a1, a4))
                    this.removeEdge(a1, a4);
            } else {
                //Add the new edge to the dgp (if it does not exist)
                if (this.contains(a1, a4)) {
                    Distance<Atom> existing = this.getDistance(a1, a4);
                    if (existing.hasBounds()) {
                        //intersect with bounds in file
                        //Get existing bounds
                        double elb = existing.getLowerBound();
                        double eub = existing.getUpperBound();

                        if (elb > dist.getUpperBound() || dist.getLowerBound() > eub) {
                            System.out.println("No overlap :(");
                            //This is problematic, because there is no intersection
                            //Throw away the existing information?
                        } else {
                            //Compute intersection
                            double nlb = Math.max(dist.getLowerBound(), elb);
                            double nub = Math.min(dist.getUpperBound(), eub);
                            dist = new Distance<>(a1, a4, nlb, nub);
                        }
                    } else if (existing.hasExpectedValue()) {
                        double elb = existing.getExpectedValue();
                        if (elb > dist.getUpperBound()) {
                            System.out.println("No overlap :(");
                            //This is problematic, because there is no intersection
                            //Throw away the existing information?
                        } else {
                            double nlb = Math.max(dist.getLowerBound(), elb);
                            dist = new Distance<>(a1, a4, nlb, dist.getUpperBound());
                        }
                    }
                    this.removeEdge(a1, a4);
                }
            }
            dist.setLabel(MDGP.NMRTorsionDistance);
            this.addEdge(dist);

            this.dihedralAngles.put(a4, new Pair<>(List.of(a3, a2, a1), new double[]{Math.toRadians(nmrAngle.minAngle), Math.toRadians(nmrAngle.maxAngle)}));
            count++;
        }

        System.out.println(count + " NMR torsion angle restraints added, average range: " + (totalRanges / count) + " degrees");

        this.addOmegas(0);
        this.addImpropers();
    }

    /**
     * Creates a {@link MDGP} based on a PDB file (for x-ray crystallography experiments)
     *
     * @param pdb               Exact distance found in dgp matching a torsion angle distance: this might be an improper dihedral
     * @param chain             The chain code of the used model
     * @param ff                The forcefield object to get the parameters from
     * @param bonds             true --> from pdb, false from ff
     * @param angles            2 --> from pdb, 1 --> from ramachandran, 0 --> from ff
     * @param dihedral_interval interval size for phi/psi in radians
     * @param seed              random seed used for random placement of dihedral angles
     */
    public MDGP(PDB pdb, char chain, ForceField ff, boolean bonds, int angles, double dihedral_interval, long seed) {
        super(3);
        this.model = pdb;
        this.chain = chain;
        this.ff = ff;
        this.dihedralAngles = new HashMap<>();

        Atom previousC = null;

        //we start building the residues one by one
        for (int resID = 1; resID <= this.model.lastResidueID; resID++) {
            //the PDB is from X-ray, so there are no hydrogens!
            Atom pdbN = this.model.getAtom("N", resID);

            String residue;
            if (pdbN != null)
                residue = pdbN.getResidueName();
            else
                residue = this.model.getAtom("C", resID).getResidueName();

            if (residue.length() > 3)
                residue = residue.substring(residue.length() - 3);

            //Create the atoms of the backbone
            Atom N = new Atom("N", "N", residue, resID);
            Atom CA = new Atom("C", "CA", residue, resID);
            Atom C = new Atom("C", "C", residue, resID);
            Atom O = new Atom("O", "O", residue, resID);

            if (resID > 1) {
                //do not forget peptide bond
                this.addAtomPDB(N, previousC, bonds, angles);
                this.addAtomPDB(previousC, N, bonds, angles);
            } else {
                //We start with N
                this.addVertex(N);
            }
            //In case the last amino acid is just an "N" atom in the PDB
            if (resID == this.model.lastResidueID && this.model.getAtom("O", resID) == null)
                continue;

            this.addAtomPDB(CA, N, bonds, angles);
            this.addAtomPDB(C, CA, bonds, angles);
            this.addAtomPDB(O, C, bonds, angles);


            //C terminus OXT
            if (resID == this.model.lastResidueID) {
                Atom pdbOXT = this.model.getAtom("OXT", resID);
                if (pdbOXT == null)
                    pdbOXT = this.model.getAtom("O2", resID);

                if (pdbOXT != null) {
                    Atom OXT = new Atom("O", pdbOXT.getName(), residue, resID); //OXT is only present at C-term
                    this.addAtomPDB(OXT, C, bonds, angles);
                }
                this.addEdge(new Distance<>(pdbN.distanceTo(this.model.getAtom("O", resID)), O, N));
                double imp = Geometry.signedDihedralAngle(pdbN.getCoords(),
                        this.model.getAtom("CA", resID).getCoords(),
                        this.model.getAtom("C", resID).getCoords(),
                        this.model.getAtom("O", resID).getCoords());
                this.dihedralAngles.put(O, new Pair<>(List.of(C, CA, N), new double[]{imp, imp}));
            }

            // hydrogens
            Atom H1;
            if (resID == 1)
                H1 = new Atom("H", "H1", residue, resID);  // bonded to N
            else
                H1 = new Atom("H", "H", residue, resID);  // bonded to N
            Atom H2 = new Atom("H", "H2", residue, resID);  // bonded to N (only present in N-terminus)
            Atom HA = new Atom("H", "HA", residue, resID);  // bonded to CA

            this.addAtomPDB(H1, N, false, 0);
            if (resID == 1) {
                this.addAtomPDB(H2, N, false, 0);
            }
            this.addAtomPDB(HA, CA, false, 0);

            //we have to add a (-180,180) distance from H1 to HA (H1, N, CA, HA dihedral angle)
            double ab = this.getExpectedValue(H1, N);
            double bc = this.getExpectedValue(N, CA);
            double ac = this.getExpectedValue(H1, CA);
            double cd = this.getExpectedValue(CA, HA);
            double bd = this.getExpectedValue(N, HA);

            double abc = Geometry.solveSSS(ab, bc, ac);
            double bcd = Geometry.solveSSS(bc, cd, bd);
            double min = Geometry.dihedralDistance(ab, bc, cd, abc, bcd, Math.toRadians(5));
            double max = Geometry.dihedralDistance(ab, bc, cd, abc, bcd, Math.toRadians(179.99));
            Distance<Atom> H1HA = new Distance<>(H1, HA, min, max);
            H1HA.setLabel(MDGP.torsionDistance);
            this.addEdge(H1HA);

            //update previousC
            previousC = C;
        }

        //relax H1-O distances
        for (Distance<Atom> dist : this.distanceSet()) {
            if (dist.hasExpectedValue()) {
                Atom a = dist.getA();
                Atom b = dist.getB();
                if (a.getName().equals("O") && b.getName().equals("H1")) {
                    //relax this distance
                    Distance<Atom> relaxed = new Distance<>(a, b, dist.getExpectedValue() - 0.5, dist.getExpectedValue() + 0.5);
                    dist = relaxed;
                    this.removeEdge(a, b);
                    this.addEdge(dist);
                }
            }
        }
        Random r = new Random(seed);
        System.out.println("MDGP constructed with " + this.numberOfVertices() + " atoms");

        //add omega angle exactly
        this.simulateDihedrals(true, false, false, 0, r);

        //add other dihedral angles in interval
        this.simulateDihedrals(false, true, true, dihedral_interval, r);
        System.out.println(this.dihedralAngles.size() + " dihedral angles simulated");
        this.addImpropers();
    }

    public MDGP(String filename, String format, ForceField ff) {
        this(filename, format, ff, "Id1 Id2 groupId1 groupId2 lb ub Name1 Name2 groupName1 groupName2");
    }

    /**
     * Creates a MDGP by reading a file (for instance, an MDjeep file)
     *
     * @param filename   the filename to read
     * @param format     the format of the file, supported: "mdjeep", "dmdgp", "pdb"
     * @param lineFormat optional parameter: mdjeep line format:
     *                   possible values in format, where * is required: Id1*, Id2* (IDs of the vertices), groupId1, groupId2 (group information for vertices, e.g., amino acids),
     *                   Name1, Name2 (name of vertices), groupName1, groupName2, (name of groups), lb*, ub* (lower and upperbound on the distance between the vertices).
     *                   Example line format = "Id1 Id2 groupId1 groupId2 lb ub Name1 Name2 groupName1 groupName2"
     * @throws IllegalArgumentException No parser for this format implemented
     */
    public MDGP(String filename, String format, ForceField ff, String lineFormat) {
        super(3);
        this.ff = ff;
        try {
            if (format.equalsIgnoreCase("mdjeep")) {
                //parse format String, make sure required elements are there and that we do not have superfluous elements
                if (!lineFormat.contains("Id1") || !lineFormat.contains("Id2") || !lineFormat.contains("lb") || !lineFormat.contains("ub"))
                    throw new IllegalArgumentException("Impossible to create DDGP object: format is missing required elements (Id1, Id2, lb, ub)");
                String[] formatParts = lineFormat.split(" ");
                Set<String> allowedParts = Set.of("Id1", "Id2", "Name1", "Name2", "groupId1", "groupId2", "groupName1", "groupName2", "ub", "lb");
                Map<String, Integer> formatIndex = new HashMap<>();
                for (int i = 0; i < formatParts.length; i++) {
                    String part = formatParts[i];
                    if (!allowedParts.contains(part))
                        throw new IllegalArgumentException("Impossible to create DDGP object: illegal elements in format: " + part);
                    formatIndex.put(part, i);
                }

                File input = new File(filename);
                Scanner scan = new Scanner(input);

                Map<String, Integer> residueID = new HashMap<>();
                while (scan.hasNext()) {
                    String line = scan.nextLine();
                    String[] parts = line.trim().replaceAll("( )+", " ").split("\\s+");
                    int id1 = Integer.parseInt(parts[formatIndex.get("Id1")]);
                    int id2 = Integer.parseInt(parts[formatIndex.get("Id2")]);

                    double lb = Double.parseDouble(parts[formatIndex.get("lb")]);
                    double ub = Double.parseDouble(parts[formatIndex.get("ub")]);

                    String name1 = parts[formatIndex.get("Name1")];
                    String name2 = parts[formatIndex.get("Name2")];

                    Atom a;
                    Atom b;

                    int residueID1 = -1;
                    int residueID2 = -1;

                    if (formatIndex.containsKey("groupId1") && formatIndex.containsKey("groupId2")) {
                        //read residue ID from file
                        residueID1 = Integer.parseInt(parts[formatIndex.get("groupId1")]);
                        residueID2 = Integer.parseInt(parts[formatIndex.get("groupId2")]);
                    } else {
                        //make the atom unique by adding id to name...
                        name1 = name1 + "-" + id1;
                        name2 = name2 + "-" + id2;
                    }


                    String residue1 = parts[formatIndex.get("groupName1")];
                    String residue2 = parts[formatIndex.get("groupName2")];

                    if (!this.containsAtom(name1, residueID1)) {
                        a = new Atom(name1.substring(0, 1), name1, residue1, residueID1);
                        //add to dgp
                        this.addVertex(a);
                    } else a = this.getAtom(name1, residueID1);

                    if (!this.containsAtom(name2, residueID2)) {
                        b = new Atom(name2.substring(0, 1), name2, residue2, residueID2);
                        //add to dgp
                        this.addVertex(b);
                    } else b = this.getAtom(name2, residueID2);

                    //add edge to dgp
                    Distance<Atom> dist;
                    if (lb == ub) dist = new Distance<>(lb, a, b);
                    else dist = new Distance<>(a, b, lb, ub);
                    this.addEdge(dist);
                }
            } else if (format.equalsIgnoreCase("dmdgp")) {
                File input = new File(filename);
                Scanner scan = new Scanner(input);
                List<Atom> atoms = new ArrayList<Atom>();

                //Get PDB code
                filename = filename.replace("DMDGP_files/", "");
                this.chain = filename.charAt(4);


                this.dihedralAngles = new HashMap<>();

                int pseudo = 0;
                int dihedrals = 0;

                int mode = -1; //0 = vertex, 1 = edges, 2 = dihedral angles
                while (scan.hasNext()) {
                    String line = scan.nextLine();
                    if (line.startsWith("end ")) {
                        mode = -1;
                    } else if (mode > -1) {
                        String[] parts = line.split("\\s+");

                        //handle the different modes
                        if (mode == 0) {
                            if (parts.length != 8)
                                throw new IllegalStateException("Error while reading an atom: the line is too short or too long! (" + parts.length + ") but should be 8");
                            String name = parts[6];
                            String residue = parts[5];
                            String residueID = residue.replaceAll("\\D+", "");
                            String residueName = residue.replace(residueID, "");
                            Atom a = new Atom(name.substring(0, 1), name, residueName, Integer.parseInt(residueID));
                            this.addVertex(a);
                            atoms.add(a);
                        } else if (mode == 1) {
                            if (parts.length < 10 || parts.length > 11)
                                throw new IllegalStateException("Error while reading a distance: the line is too short or too long! (" + parts.length + ") but should be 10 or 11 parts");

                            int i = 2;
                            String type = parts[i++];
                            boolean exact = type.equals("D");
                            double lb = Double.parseDouble(parts[i++]);
                            double ub = lb;

                            if (!exact) {
                                ub = Double.parseDouble(parts[i++]);
                            }

                            i++; //skip the #
                            String residue1 = parts[i++];
                            String name1 = parts[i++];
                            String residueID1 = residue1.replaceAll("\\D+", "");

                            i++; //skip the --
                            String residue2 = parts[i++];
                            String name2 = parts[i++];
                            String residueID2 = residue2.replaceAll("\\D+", "");

                            Distance<Atom> dist;
                            Atom a = this.getAtom(name1, Integer.parseInt(residueID1));
                            Atom b = this.getAtom(name2, Integer.parseInt(residueID2));

                            if (lb == ub) {
                                dist = new Distance<>(lb, a, b);
                                if (a.getName().equals("O") && b.getName().equals("H1")) {
                                    //relax this distance
                                    Distance<Atom> relaxed = new Distance<>(a, b, dist.getExpectedValue() - 0.5, dist.getExpectedValue() + 0.5);
                                    dist = relaxed;
                                }
                            } else
                                dist = new Distance<>(a, b, lb, ub);
                            this.addEdge(dist);
                        } else if (mode == 2) {
                            if (parts.length < 6 || parts.length > 7)
                                throw new IllegalStateException("Error while reading a dihedral angle: the line is too short or too long! (" + parts.length + ") but should be 6 or 7 parts");

                            int i = 0;
                            //-1 because our array is 0-indexed
                            Atom a = atoms.get(Integer.parseInt(parts[i++]) - 1);
                            Atom b = atoms.get(Integer.parseInt(parts[i++]) - 1);
                            Atom c = atoms.get(Integer.parseInt(parts[i++]) - 1);
                            Atom d = atoms.get(Integer.parseInt(parts[i++]) - 1);
                            boolean exact = parts[i++].equals("D");
                            double minAngle = Double.parseDouble(parts[i++]);
                            double maxAngle = minAngle;
                            Distance<Atom> existing = this.getDistance(a, d);

                            if (!exact)
                                maxAngle = Double.parseDouble(parts[i++]);

                            //skip the useless information
                            if (minAngle == -180 && maxAngle == 180)
                                continue;
                            else {
                                minAngle = Math.toRadians(minAngle);
                                maxAngle = Math.toRadians(maxAngle);
                            }

                            if (existing.hasExpectedValue()) {
                                dihedralAngles.put(d, new Pair<>(List.of(c, b, a), new double[]{minAngle, maxAngle}));
                                //N, CA, HA, C
                                pseudo++;
                            } else if (!exact) {
                                dihedralAngles.put(d, new Pair<>(List.of(c, b, a), new double[]{minAngle, maxAngle}));

                                dihedrals++;
                                double ab = this.getExpectedValue(a, b);
                                double ac = this.getExpectedValue(a, c);
                                double bc = this.getExpectedValue(b, c);
                                double bd = this.getExpectedValue(b, d);
                                double cd = this.getExpectedValue(c, d);

                                double abc = Geometry.solveSSS(ab, bc, ac);
                                double bcd = Geometry.solveSSS(bc, cd, bd);

                                //compute torsion distance
                                double lb = Geometry.dihedralDistance(ab, bc, cd, abc, bcd, minAngle);
                                double ub = Geometry.dihedralDistance(ab, bc, cd, abc, bcd, maxAngle);

                                Distance<Atom> dist;
                                if (lb == ub) {
                                    continue;
                                } else if (lb < ub) dist = new Distance<>(a, d, lb, ub);
                                else dist = new Distance<>(a, d, ub, lb);

                                if (this.contains(a, d)) {
                                    if (existing.hasBounds()) {
                                        double elb = existing.getLowerBound();
                                        double eub = existing.getUpperBound();
                                        if (eub < dist.getLowerBound() || dist.getUpperBound() < elb) {
                                            double nlb = Math.min(lb, elb);
                                            double nub = Math.max(ub, eub);
                                            dist = new Distance<>(a, d, nlb, nub);
                                            //throw new IllegalStateException("There is no intersection between a distance calculated using the angles " + minAngle + " and " + maxAngle + " and the distance " + lb + " and " + ub);
                                            //This is problematic, because there is no intersection
                                            //Throw away the existing information?
                                        } else {
                                            //Compute intersection
                                            double nlb = Math.max(dist.getLowerBound(), elb);
                                            double nub = Math.min(dist.getUpperBound(), eub);
                                            dist = new Distance<>(a, d, nlb, nub);
                                        }
                                    }
                                    this.removeEdge(a, d);
                                    this.addEdge(dist);
                                }
                            }
                        }
                    }

                    //recognize the mode
                    if (line.equalsIgnoreCase("begin vertices"))
                        mode = 0;
                    else if (line.equalsIgnoreCase("begin edges")) {
                        if (this.numberOfVertices() == 0)
                            throw new IllegalStateException("We could not find any atoms: please make sure the lines describing the vertices are above the distances in the DMDGP file!");
                        mode = 1;
                    } else if (line.equalsIgnoreCase("begin dihedral_angles")) {
                        if (this.numberOfVertices() == 0)
                            throw new IllegalStateException("We could not find any distances: please make sure the lines describing the distances are above the torsion angles in the DMDGP file!");
                        mode = 2;
                    }
                }
                //System.out.println("Chain code: " + chain);
                //System.out.println("Number of edges: " + this.numberOfEdges());
                //System.out.println("Number of impropers: " + pseudo);
                //System.out.println("Number of real dihedrals: " + dihedrals);
            } else if (format.equalsIgnoreCase("pdb")) {
                //read PDB file
                PDB pdb = new PDB(filename, 1);

                //create a complete instance from the PDB file
                for (Atom a : pdb.getAtoms()) {
                    this.addVertex(a);
                    for (Atom b : this.vertexList()) {
                        if (a != b)
                            this.addEdge(new Distance<>(a, b));
                    }
                }
            } else throw new IllegalArgumentException("No parser for this format implemented");
        } catch (Exception e) {
            System.out.println("Something went wrong reading the input file...");
            e.printStackTrace();
            System.exit(1);
        }

        System.out.println("MDGP contructed with " + this.numberOfVertices() + " atoms and " + this.numberOfEdges() + " distances!");
    }


    /**
     * Adds an atom to the {@link MDGP}, connected to a previously added atom (NMR instances)
     *
     * @param next The atom to add to the dgp
     * @param prev The atom already in the dgp that next is bonded to
     */
    private void addAtom(Atom next, Atom prev) {
        //Add the new atom to the amino acid as a vertex
        this.addVertex(next);

        //Using the previous atom, identify all the atoms (already in the dgp) that we need to create an edge to
        Set<Atom> prevBondNeighbours = new HashSet<>();
        Map<Atom, List<Atom>> prevAngleNeighbours = new HashMap<>();

        //identify direct neighbours and ADD BOND distances
        for (Atom a : this.starOf(prev)) {
            if (a.equals(next))
                continue;

            Distance<Atom> prevDist = this.getDistance(prev, a);
            if (prevDist.getLabel().equals(MDGP.bondDistance)) {
                prevBondNeighbours.add(a);
            } else if (prevDist.getLabel().equals(MDGP.angleDistance)) {
                //find the corresponding bond neighbour of prev that created this angle
                for (Atom b : this.starOf(prev)) {
                    if (this.contains(b, a)) {
                        Distance<Atom> prevToB = this.getDistance(prev, b);
                        Distance<Atom> bToA = this.getDistance(b, a);

                        if (prevToB.getLabel().equals(MDGP.bondDistance) && bToA.getLabel().equals(MDGP.bondDistance)) {
                            if (!prevAngleNeighbours.containsKey(b))
                                prevAngleNeighbours.put(b, new ArrayList<>());
                            prevAngleNeighbours.get(b).add(a);
                        }
                    }
                }
            }
        }

        //Add the connection between the new atom and the previous
        Distance<Atom> bond = new Distance<>(this.ff.expectedBondDistance(next, prev), next, prev);
        bond.setLabel(MDGP.bondDistance);

        //A bond edge carries more importance than angle or torsion angle edges
        //This is important for amino acids with cycles!
        if (this.contains(next, prev))
            this.removeEdge(next, prev);
        this.addEdge(bond);

        //Only if we want angles too!
        if(!onlybonds) {
            //ADD bond angles
            for (Atom a : prevBondNeighbours) {
                //double delta = 0.25;
                double angleDist = this.ff.expectedBondAngleDistance(next, prev, a);
                Distance<Atom> angleDistance = new Distance<>(angleDist, next, a);
                angleDistance.setLabel(MDGP.angleDistance);

                //We save this distance if next,a are not connected yet, or if connected by a torsion angle (we much prefer to use the ambertype angles)
                if (!this.contains(next, a))
                    this.addEdge(angleDistance);
                else if (this.contains(next, a) && this.getDistance(next, a).getLabel().equals(MDGP.torsionDistance)) {
                    this.removeEdge(next, a);
                    this.addEdge(angleDistance);
                }
            }

            //ADD torsion distances
            for (Map.Entry<Atom, List<Atom>> entry : prevAngleNeighbours.entrySet()) {
                Atom C = entry.getKey();
                for (Atom D : entry.getValue()) {
                    if (!this.contains(next, D)) {
                        double lb = this.ff.expectedTorsionAngleDistance(next, prev, C, D, 10);
                        double ub = this.ff.expectedTorsionAngleDistance(next, prev, C, D, 170);
                        Distance<Atom> torsionDist = new Distance<>(next, D, lb, ub);
                        torsionDist.setLabel(MDGP.torsionDistance);
                        this.addEdge(torsionDist);
                    }
                }
            }
        }
    }

    /**
     * Adds an atom to the {@link MDGP}, connected to a previously added atom (X-RAY instances)
     *
     * @param next The atom to add to the dgp
     * @param prev The atom already in the dgp that next is bonded to
     * @param bonds If true, we take the bond from the PDB, false from ff
     * @param angles If 2, we take the bond from the PDB, 1 from Ramachandran, 0 from ff
     */
    private void addAtomPDB(Atom next, Atom prev, boolean bonds, int angles) {
        Atom PDBnext = this.model.getAtom(next.getName(), next.getResidueID());
        Atom PDBprev = this.model.getAtom(prev.getName(), prev.getResidueID());

        //Add the new atom to the amino acid as a vertex
        this.addVertex(next);

        //Using the previous atom, identify all the atoms (already in the dgp) that we need to create an edge to
        Set<Atom> prevBondNeighbours = new HashSet<>();
        Map<Atom, List<Atom>> prevAngleNeighbours = new HashMap<>();

        //identify direct neighbours and ADD BOND distances
        for (Atom a : this.starOf(prev)) {
            if (a.equals(next))
                continue;
            Distance<Atom> prevDist = this.getDistance(prev, a);
            if (prevDist.getLabel().equals(MDGP.bondDistance)) {
                prevBondNeighbours.add(a);
            }
        }

        //Add the connection between the new atom and the previous, based on force field or on PDB
        Distance<Atom> bond;
        if (!bonds || PDBnext == null || PDBprev == null)
            bond = new Distance<>(this.ff.expectedBondDistance(next, prev), next, prev);
        else
            bond = new Distance<>(PDBnext.distanceTo(PDBprev), next, prev);

        bond.setLabel(MDGP.bondDistance);

        //A bond edge carries more importance than angle or torsion angle edges
        //This is important for amino acids with cycles!
        if (this.contains(next, prev))
            this.removeEdge(next, prev);
        this.addEdge(bond);

        //ADD bond angles
        for (Atom a : prevBondNeighbours) {
            Atom PDBa = this.model.getAtom(a.getName(),a.getResidueID());
            Double angleDist = null;
            if (angles == 0 || PDBnext == null || PDBa == null || PDBprev == null)
                angleDist = this.ff.expectedBondAngleDistance(next, prev, a);
            else if(angles == 1){
                //Derive from Ramachandran plot
                Pair<String, Integer> angleRes = this.ff.ramachandranAngleIdentifer(a, prev, next);
                if(angleRes != null){
                    Double phi = this.model.computePhi(angleRes._2());
                    Double psi = this.model.computePhi(angleRes._2());
                    if(phi != null && psi != null){
                        double angle = this.ff.ramachandranAngle(angleRes._1(), phi, psi);
                        double abG = this.getExpectedValue(a, prev);
                        double bcG = this.getExpectedValue(prev,next);
                        angleDist = Geometry.solveSAS(abG, bcG, angle);
                    }
                }
                //we failed to use the Ramachandran value...
                if(angleDist == null)
                    angleDist = this.ff.expectedBondAngleDistance(next, prev, a);
            }
            else {
                //First: calculate angle from PDB
                double ab = PDBa.distanceTo(PDBprev);
                double bc = PDBprev.distanceTo(PDBnext);
                double ac = PDBa.distanceTo(PDBnext);
                double angle = Geometry.solveSSS(ab, bc, ac);

                //Next: calculate ac using two distances in the dgp + abc
                double abG = this.getExpectedValue(a, prev);
                double bcG = this.getExpectedValue(prev,next);
                angleDist = Geometry.solveSAS(abG, bcG, angle);

                //System.out.println(angleDist + " VS " + ac);
            }

            Distance<Atom> angleDistance = new Distance<>(angleDist, next, a);
            angleDistance.setLabel(MDGP.angleDistance);

            //We save this distance if next,a are not connected yet, or if connected by a torsion angle (we much prefer to use the ambertype angles)
            this.addEdge(angleDistance);
        }
    }

    /**
     * Retrieves an atom by name and residueID in the {@link MDGP}
     *
     * @param name      the name of the atom
     * @param residueID the ID of the residue
     * @return The atom if it exists, null otherwise
     * @throws IllegalStateException The atom was not found in the {@link MDGP}
     */
    public Atom getAtom(String name, int residueID) {
        for (Atom a : this.vertexList())
            if (a.getName().equals(name) && a.getResidueID() == residueID) return a;
        return null;
    }

    /**
     * Checks if the {@link MDGP}contains an atom by name and residueID
     *
     * @param name      the name of the atom
     * @param residueID the ID of the residue
     * @return True if the atom exists, false otherwise
     * @throws IllegalStateException The atom was not found in the {@link MDGP}
     */
    private boolean containsAtom(String name, int residueID) {
        return this.getAtom(name, residueID) != null;
    }

    /**
     * Counts the number of atoms by type in the {@link MDGP}
     *
     * @param type the type to look for
     * @return The number of atoms of type in the {@link MDGP}
     * @throws IllegalStateException The atom was not found in the {@link MDGP}
     */
    public int countType(String type) {
        int count = 0;
        for (Atom a : this.vertexList())
            if (a.getTypeName().equals(type)) count++;
        return count;
    }

    /**
     * Adds the van der Waals lower bounds to the instance
     * @param multiplier the multiplier of the sum of the radii
     * @param atomNames the atoms from which a van der Waals distance is allowed to originate
     */
    public void addVanDerWaals(boolean replace, double multiplier, Set<String> atomNames){
        //add van der Waals
        double maxDistance = 100;

        int count = 0;

        //remove unnecessary VDW constraints
        for(int i = 0; i < this.vertexList().size(); i++) {
            Atom a = this.vertexList().get(i);

            for (int j = i + 1; j < this.vertexList().size(); j++) {
                Atom b = this.vertexList().get(j);

                if (!this.contains(a, b)) {
                    if(atomNames == null || (atomNames.contains(a.getName()) && atomNames.contains(b.getName()))) {
                        double vanDerWaals = this.ff.vanDerWaalsDistance(a, b);
                        //Use 80% of the van der Waals distance as lowerbound and 120% of the largest expected distance as upperbound
                        Distance<Atom> dist = new Distance<>(a, b, vanDerWaals * multiplier, maxDistance);
                        dist.setLabel(MDGP.waalsDistance);
                        this.addEdge(dist);
                        count++;
                    }
                }
                else if(replace) {
                    Distance<Atom> oldDist = this.getDistance(a,b);
                    if(oldDist.hasBounds()) {
                        double range = oldDist.getUpperBound() - oldDist.getLowerBound();
                        if(range > 1.2){
                            if(Math.abs(a.getResidueID() - b.getResidueID()) < 3 ||(atomNames != null && !atomNames.contains(a.getName())) || ((atomNames != null && !atomNames.contains(b.getName()))))
                                this.removeEdge(a,b);
                            else {
                                //make stricter
                                Distance<Atom> newDist = new Distance<>(a, b, this.ff.vanDerWaalsDistance(a, b) * 0.8, oldDist.getUpperBound() + 0.5);
                                newDist.setLabel(MDGP.waalsDistance);
                                this.removeEdge(a, b);
                                this.addEdge(newDist);
                                count++;
                            }
                        }
                    }
                }
            }
        }

        System.out.println(count + " van der Waals lowerbounds added!");
    }

    /**
     * Add omega angle information to the input
     * @param interval the size of the created interval in radians [pi/2 - interval, pi/2]
     */
    public void addOmegas(double interval) {
        int noResidues = this.numberOfResidues();
        for(int i = 1; i < noResidues;i++){
            Atom CA = this.getAtom("CA", i);
            Atom C = this.getAtom("C", i);
            Atom N = this.getAtom("N", i + 1);
            Atom CA_n = this.getAtom("CA", i + 1);

            double ab = this.getExpectedValue(CA, C);
            double bc = this.getExpectedValue(C, N);
            double ac = this.getExpectedValue(CA,N);
            double abc = Geometry.solveSSS(ab,bc,ac);

            double bd = this.getExpectedValue(C, CA_n);
            double cd = this.getExpectedValue(N, CA_n);
            double bcd = Geometry.solveSSS(bc,cd,bd);

            double lower = Math.toRadians(180) - interval;
            double upper = Math.toRadians(180);

            double lower_dist = Geometry.dihedralDistance(ab, bc, cd, abc, bcd, lower);
            double upper_dist = Geometry.dihedralDistance(ab, bc, cd, abc, bcd, upper);

            Distance<Atom> omega_dist;
            if(interval == 0)
                omega_dist = new Distance<>(lower_dist, CA,CA_n);
            else
                omega_dist = new Distance<>(CA,CA_n, lower_dist, upper_dist);

            this.removeEdge(CA,CA_n);
            this.addEdge(omega_dist);

            this.dihedralAngles.put(CA_n, new Pair<>(List.of(N,C,CA), new double[]{lower,upper}));
        }

    }

    /**
     * Add improper angle information (from force field)
     * NOTE: these angles are hardcoded because there are only three relevant improper angles and they are always the same value for each of the force-fields.
     */
    public void addImpropers(){
        //number of residues:
        int numberOfRes = this.numberOfResidues();

        //add improper for HA
        for(int i = 1; i <= numberOfRes; i++) {
            Atom N = this.getAtom("N", i);
            Atom CA = this.getAtom("CA", i);
            Atom C = this.getAtom("C", i);
            Atom HA = this.getAtom("HA", i);

            if(HA == null)
                continue;

            //this angle is always around 119 degrees
            double improper = Math.toRadians(119);
            this.dihedralAngles.put(HA, new Pair<>(List.of(C,CA,N), new double[]{improper,improper}));
        }

        //add improper for H
        for(int i = 2; i <= numberOfRes; i++){
            Atom C = this.getAtom("C", i-1);
            Atom N = this.getAtom("N", i);
            Atom CA = this.getAtom("CA", i);
            Atom H = this.getAtom("H", i);

            if(H == null)
                continue;

            //this angle is always around 180 degrees
            double improper = Math.toRadians(180);
            this.dihedralAngles.put(H, new Pair<>(List.of(CA,N,C), new double[]{improper,improper}));
        }

        //add improper for O
        for(int i = 1; i < numberOfRes; i++) {
            Atom C = this.getAtom("C", i);
            Atom CA = this.getAtom("CA", i);
            Atom N = this.getAtom("N", i + 1);
            Atom O = this.getAtom("O", i);

            //this angle is around 0 degrees
            //double imp = Geometry.signedDihedralAngle(pdbC.getCartesians(), pdbCA.getCartesians(), pdbN.getCartesians(), pdbO.getCartesians());
            //double imp = Math.toRadians(0);
            this.dihedralAngles.put(O, new Pair<>(List.of(N, CA, C), new double[]{Math.toRadians(0), Math.toRadians(0)}));
        }
    }

    /**
     * Add improper angle information to the input (from PDB model)
     * NOTE: this only works without sidechain!
     */
    public void simulateImpropers(){
        //number of residues:
        int numberOfRes = this.model.numberOfResidues();

        System.out.println("Next HA:");

        //add improper for HA
        for(int i = 1; i < numberOfRes; i++){
            Atom N = this.getAtom("N", i);
            Atom CA = this.getAtom("CA", i);
            Atom C = this.getAtom("C", i);
            Atom HA = this.getAtom("HA", i);
            if(HA == null)
                HA = this.getAtom("HA2",i);


            //glycine
            int pdbResID = i;
            Atom pdbN = this.model.getAtom("N", pdbResID);
            Atom pdbCA = this.model.getAtom("CA", pdbResID);
            Atom pdbC = this.model.getAtom("C", pdbResID);
            Atom pdbHA = this.model.getAtom("HA", pdbResID);

            if(HA.getResidueName().equals("GLY"))
                pdbHA = this.model.getAtom("HA2", pdbResID);

            double imp = Geometry.signedDihedralAngle(pdbN.getCoords(), pdbCA.getCoords(), pdbC.getCoords(), pdbHA.getCoords());
            this.dihedralAngles.put(HA, new Pair<>(List.of(C,CA,N), new double[]{imp,imp}));

            System.out.println(Math.toDegrees(imp));

        }

        System.out.println("Next H:");

        //add improper for H
        for(int i = 2; i < numberOfRes; i++){
            Atom C = this.getAtom("C", i-1);
            Atom N = this.getAtom("N", i);
            Atom CA = this.getAtom("CA", i);
            Atom H = this.getAtom("H", i);

            int pdbResID = i;
            Atom pdbC = this.model.getAtom("C", pdbResID-1);
            Atom pdbN = this.model.getAtom("N", pdbResID);
            Atom pdbCA = this.model.getAtom("CA", pdbResID);
            Atom pdbH = this.model.getAtom("H", pdbResID);

            if(pdbH == null)
                continue;

            double imp = Geometry.signedDihedralAngle(pdbC.getCoords(), pdbN.getCoords(), pdbCA.getCoords(), pdbH.getCoords());
            this.dihedralAngles.put(H, new Pair<>(List.of(CA,N,C), new double[]{imp,imp}));
            System.out.println(Math.toDegrees(imp));

        }

        System.out.println("Next O:");
        //add improper for O
        for(int i = 1; i < numberOfRes; i++){
            Atom C = this.getAtom("C", i);
            Atom CA = this.getAtom("CA", i);
            Atom N = this.getAtom("N", i+1);
            Atom O = this.getAtom("O", i);

            int pdbResID = i;
            Atom pdbC = this.model.getAtom("C", pdbResID);
            Atom pdbCA = this.model.getAtom("CA", pdbResID);
            Atom pdbN = this.model.getAtom("N", pdbResID+1);
            Atom pdbO = this.model.getAtom("O", pdbResID);

            double imp = Geometry.signedDihedralAngle(pdbC.getCoords(), pdbCA.getCoords(), pdbN.getCoords(), pdbO.getCoords());
            this.dihedralAngles.put(O, new Pair<>(List.of(N,CA,C), new double[]{imp,imp}));

            System.out.println(Math.toDegrees(imp));
        }

    }

    /**
     * Simulate dihedral angles from PDB, in a given range
     * @param omega whether we simulate omega
     * @param psi whether we simulate psi
     * @param phi whether we simulate phi
     * @param range the range of the intervals (for phi and psi) in radians
     */
    public void simulateDihedrals(boolean omega, boolean psi, boolean phi, double range, Random r){
        try {
            if (this.model == null)
                throw new IllegalStateException("We can only simulate dihedral angles when we have access to the PDB model!");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        int noResidues = this.numberOfResidues();
        for(int i = 1; i < noResidues;i++) {
            Atom N = this.getAtom("N", i);
            Atom CA = this.getAtom("CA", i);
            Atom C = this.getAtom("C", i);
            Atom N_n = this.getAtom("N", i + 1);
            Atom CA_n = this.getAtom("CA", i + 1);
            Atom C_n = this.getAtom("C", i+1);

            Atom pdbCA = this.model.getAtom("CA", i);
            Atom pdbCA_n = this.model.getAtom("CA", i + 1);
            Atom pdbN = this.model.getAtom("N", i);
            Atom pdbN_n = this.model.getAtom("N", i + 1);
            Atom pdbC = this.model.getAtom("C", i);
            Atom pdbC_n = this.model.getAtom("C", i + 1);


            if(omega && CA_n != null) {
                Double dihedral = this.model.computeOmega(i);

                //we do not have PDB info available, simply use VDW for distance then...
                if (dihedral == null)
                    this.addEdge(new Distance<>(CA, CA_n, ff.vanDerWaalsDistance(CA, CA_n), 5.0));
                else {
                    this.dihedralAngles.put(CA_n, new Pair<>(List.of(N_n, C, CA), new double[]{dihedral, dihedral}));
                    Distance<Atom> dist = new Distance<>(pdbCA.distanceTo(pdbCA_n), CA, CA_n);
                    this.removeEdge(CA, CA_n);
                    this.addEdge(dist);
                }
            }
            if(phi && C_n != null) {
                Double dihedral = this.model.computePhi(i);
                double below = r.nextDouble() * range;
                double above = range - below;
                double lower = Math.max(-Math.PI, dihedral - below);
                double upper = Math.min(Math.PI, dihedral + above);

                //compute lower and upperbound on distances based on angle range
                Atom a = pdbC, b = pdbN_n, c = pdbCA_n, d = pdbC_n;
                double ab = a.distanceTo(b);
                double bc = b.distanceTo(c);
                double ac = a.distanceTo(c);
                double cd = c.distanceTo(d);
                double bd = b.distanceTo(d);
                double abc = Geometry.solveSSS(ab, bc, ac);
                double bcd = Geometry.solveSSS(bc, cd, bd);
                double one = Geometry.dihedralDistance(ab, bc, cd, abc, bcd, lower);
                double two = Geometry.dihedralDistance(ab, bc, cd, abc, bcd, upper);
                double min = Math.min(one,two);
                double max = Math.max(one,two);



                this.dihedralAngles.put(C_n, new Pair<>(List.of(CA_n, N_n, C), new double[]{lower, upper}));
                Distance<Atom> dist = new Distance<>(C, C_n, min, max);
                this.removeEdge(C, C_n);
                this.addEdge(dist);
            }
            if(psi) {
                //we do not have PDB info available, simply use VDW for distance then...
                if ((pdbN == null || pdbCA == null) && !this.contains(N, N_n)) {
                    this.addEdge(new Distance<>(N, N_n, ff.vanDerWaalsDistance(N, N_n), 5.0));
                }
                else {
                    Double dihedral = this.model.computePsi(i);
                    double below = r.nextDouble() * range;
                    double above = range - below;
                    double lower = Math.max(-Math.PI, dihedral - below);
                    double upper = Math.min(Math.PI, dihedral + above);
                    this.dihedralAngles.put(N_n, new Pair<>(List.of(C, CA, N), new double[]{lower, upper}));

                    //compute lower and upperbound on distances based on angle range
                    Atom a = pdbN, b = pdbCA, c = pdbC, d = pdbN_n;
                    double ab = a.distanceTo(b);
                    double bc = b.distanceTo(c);
                    double ac = a.distanceTo(c);
                    double cd = c.distanceTo(d);
                    double bd = b.distanceTo(d);
                    double abc = Geometry.solveSSS(ab, bc, ac);
                    double bcd = Geometry.solveSSS(bc, cd, bd);
                    double one = Geometry.dihedralDistance(ab, bc, cd, abc, bcd, lower);
                    double two = Geometry.dihedralDistance(ab, bc, cd, abc, bcd, upper);
                    double min = Math.min(one,two);
                    double max = Math.max(one,two);

                    Distance<Atom> dist = new Distance<>(N, N_n, min, max);
                    this.removeEdge(N, N_n);
                    this.addEdge(dist);
                }
            }
        }
    }

    /**
     * Adds the long range Cbeta-Cbeta distances to the instance by simulating them.
     * Distances between neighbouring and next-neighbouring residues are not included.
     * @param maxDistance only distances below maxDistance are included in the instance
     * @param interval the size of the created interval
     */
    public void simulateCbetas(double maxDistance, double interval, boolean illegal, String CBCB_file) {
        try {
            if (this.model == null)
                throw new IllegalStateException("We can only simulate cBeta distances when we have access to the PDB model!");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        List<Atom> pdbCbetas = new ArrayList<>();
        List<Atom> cbetas = new ArrayList<>();
        int maxResidues = numberOfResidues();
        //Add C-beta atoms to instance
        int count = 0;

        for(int i = 1; i <= maxResidues; i++) {
            Atom CA = this.getAtom("CA", i);
            Atom N = this.getAtom("N", i);
            Atom C = this.getAtom("C", i);
            Atom CB = this.getAtom("CB", i);

            Atom pdbCB = this.model.getAtom("CB", i);
            if (pdbCB == null)
                continue;


            Atom pdbCA = this.model.getAtom("CA", i);
            Atom pdbN = this.model.getAtom("N", i);
            Atom pdbC = this.model.getAtom("C", i);

            //if we do not have a Cbeta atom...
            if (CB == null) {
                //...we create it!
                CB = new Atom("C", "CB", CA.getResidueName(), i);
                this.addVertex(CB);

                //and we add necessary distances!
                Distance<Atom> CaCb = new Distance<>(pdbCB.distanceTo(pdbCA), CA, CB);
                this.addEdge(CaCb);

                this.addEdge(new Distance<>(pdbCB.distanceTo(pdbN), N, CB));
                this.addEdge(new Distance<>(pdbCB.distanceTo(pdbC), C, CB));
            }

            //add the dihedral angle too!
            double torsion = Geometry.signedDihedralAngle(pdbN.getCoords(), pdbCA.getCoords(), pdbC.getCoords(), pdbCB.getCoords());
            this.dihedralAngles.put(CB, new Pair<>(List.of(C, CA, N), new double[]{torsion, torsion}));

            count++;
            cbetas.add(CB);
            pdbCbetas.add(pdbCB);
        }

        for(int i = 1; i <= maxResidues; i++) {
            Atom CB = this.getAtom("CB", i);
            Atom PDB_CB = this.model.getAtom("CB", i);

            if (PDB_CB == null)
                continue;

            for (int j = i + 2; j <= maxResidues && j <= i + 5; j++) {
                Atom CB_next = this.getAtom("CB", j);
                Atom PDB_CB_next = this.model.getAtom("CB", j);

                if (PDB_CB_next != null) {
                    double true_distance = PDB_CB.distanceTo(PDB_CB_next);
                    this.addEdge(new Distance<>(CB, CB_next, Math.max(true_distance - interval / 2, this.ff.vanDerWaalsDistance(CB, CB_next) * 0.8), true_distance + interval / 2));
                }
            }
        }



//        //read the distances from the CBCB file
//        try {
//            File input = new File(CBCB_file);
//            Scanner scan = new Scanner(input);
//            while (scan.hasNext()) {
//                String line = scan.nextLine();
//                if(!line.contains("cb1"))
//                    continue;
//
//                String[] parts = line.trim().replaceAll("( )+", " ").split("\\s+");
//
//                double true_distance = Double.parseDouble(parts[5]);
//                int resID1 = Math.abs(Integer.parseInt(parts[1]) - this.model.firstResidueID + 1);
//                int resID2 = Math.abs(Integer.parseInt(parts[3]) - this.model.firstResidueID + 1);
//
//                Atom CB1 = this.getAtom("CB",resID1);
//                Atom CB2 = this.getAtom("CB",resID2);
//
//                if(!this.contains(CB1,CB2))
//                    this.addDistance(new Distance<>(CB1, CB2, Math.max(true_distance - interval/2,this.ff.vanDerWaalsDistance(CB1, CB2) * 0.8), true_distance + interval/2));
//            }
//
//        } catch (Exception e) {
//            e.printStackTrace();
//            System.exit(1);
//        }
    }

    /**
     * Creates a hardcoded atom order for the experiments with the x-ray crystallography proteins. This is the best (discretizable) order to exploit dihedral angles!
     * @return the order in terms of a list of atoms
     */
    public List<Atom> dmdgpOrder(){
        List<Atom> newOrder = new ArrayList<>();

        int numberOfResidues = this.numberOfResidues();

        //first amino acid
        newOrder.add(this.getAtom("N", 1));

        //Add H1, H2 and H3 if they are not null
        Atom H1 = this.getAtom("H1", 1);
        Atom H2 = this.getAtom("H2", 1);
        Atom H3 = this.getAtom("H3", 1);

        if(H1 == null) H1 = this.getAtom("H",1);
        if(H1 != null) newOrder.add(H1);
        if(H2 != null) newOrder.add(H2);
        if(H3 != null) newOrder.add(H3);

        newOrder.add(this.getAtom("CA", 1));
        newOrder.add(this.getAtom("HA", 1));
        newOrder.add(this.getAtom("C", 1));

        //add CB if it exists
        Atom CB = this.getAtom("CB", 1);
        if(CB != null) newOrder.add(CB);

        //from 2 onwards
        int count = 0;
        for(int i = 2; i <= numberOfResidues; i++){
            newOrder.add(this.getAtom("N", i));
            newOrder.add(this.getAtom("O", i-1));
            newOrder.add(this.getAtom("CA", i));
            newOrder.add(this.getAtom("C", i));

            //H1 does not always exist and sometimes it is called "H"
            H1 = this.getAtom("H1", i);
            if(H1 == null) H1 = this.getAtom("H", i);
            if(H1 != null) newOrder.add(H1);

            newOrder.add(this.getAtom("HA", i));
            //add CB if it exists
            CB = this.getAtom("CB", i);
            if(CB != null){
                newOrder.add(CB);
                count++;
            }
        }

        //add the last 0 and 02
        newOrder.add(this.getAtom("O", numberOfResidues));

        //try to add the last oxygen...
        Atom OXT = this.getAtom("OXT", numberOfResidues);
        if(OXT == null)
            OXT = this.getAtom("O2", numberOfResidues);
        if(OXT != null) newOrder.add(OXT);

        //remove atoms that did not exist in this specific model
        while (newOrder.remove(null));
        return newOrder;
    }

    /**
     * Adds the atoms of a residue backbone to the {@link MDGP} instance
     * @param residue   the residue three letter code
     * @param residueID the position of the residue in the {@link MDGP}
     * @param atNTerm   are we at an N-terminus?
     * @param atCTerm   are we at a C-terminus?
     */
    private void addBackbone(String residue, int residueID, boolean atNTerm, boolean atCTerm) {
        Atom N;
//        if (atNTerm)
//            N = new Atom("N", "N", "NXT", residue, residueID);
//        else
            N = new Atom("N", "N",residue, residueID);
        Atom CA = new Atom("C", "CA", residue, residueID);
        Atom C = new Atom("C", "C",residue, residueID);
        Atom O = new Atom("O", "O", residue, residueID); //O is the one that remains when the amino acids form a peptide bond

        // hydrogens
        Atom H1;
        if (atNTerm)
            H1 = new Atom("H", "H1", residue, residueID);  // bonded to N
        else
            H1 = new Atom("H", "H", residue, residueID);  // bonded to N
        Atom H2 = new Atom("H", "H2", residue, residueID);  // bonded to N (only present in N-terminus)
        //Atom H3 = new Atom("H", "H3", residue, residueID);  // bonded to N (only present in N-terminus)
        Atom HA = new Atom("H", "HA", residue, residueID);  // bonded to CA

        //We start with N
        this.addVertex(N);
        this.addAtom(CA, N);
        this.addAtom(C, CA);
        this.addAtom(O, C);

        if (atCTerm) {
            Atom OXT = new Atom("O", "OXT", residue, residueID); //OXT is only present at C-term
            this.addAtom(OXT, C);
        }

        if (hydrogens) {
            this.addAtom(H1, N);
            if (atNTerm) {
                this.addAtom(H2, N);
                //this.addAtom(H3, N);
            }
            this.addAtom(HA, CA);
        }
    }

    /**
     * Adds the side chain of an amino acid to the MDGP instance
     * Note: we construct the full amino acid - the N-atom has three H-atoms and the C-atom has two connected O-atoms
     *
     * @param residue   the residue three letter code
     * @param residueID the position of the residue in the {@link MDGP}
     */
    private void addSideChain(String residue, int residueID) {
        //extract existing atoms from backbone
        Atom CA = this.getAtom("CA", residueID);
        Atom N = this.getAtom("N", residueID);

        if (residue.equals("ALA")) { //alanine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);

            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB1 = new Atom("H", "HB1", residue, residueID);
                Atom HB2 = new Atom("H", "HB2", residue, residueID);
                Atom HB3 = new Atom("H", "HB3", residue, residueID);

                //Connect hydrogens
                addAtom(HB1, CB);
                addAtom(HB2, CB);
                addAtom(HB3, CB);
            }
        } else if (residue.equals("ARG")) { //arganine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG = new Atom("C", "CG", residue, residueID);
            Atom CD = new Atom("C", "CD", residue, residueID);
            Atom NE = new Atom("N", "NE", residue, residueID);
            Atom CZ = new Atom("C", "CZ", residue, residueID);
            Atom NH1 = new Atom("N", "NH1", residue, residueID);
            Atom NH2 = new Atom("N", "NH2", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG, CB);
            addAtom(CD, CG);
            addAtom(NE, CD);
            addAtom(CZ, NE);
            addAtom(NH1, CZ);
            addAtom(NH2, CZ);

            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB2 = new Atom("H", "HB2", residue, residueID);
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                Atom HG2 = new Atom("H", "HG2", residue, residueID);
                Atom HG3 = new Atom("H", "HG3", residue, residueID);
                Atom HD2 = new Atom("H", "HD2", residue, residueID);
                Atom HD3 = new Atom("H", "HD3", residue, residueID);
                Atom HE = new Atom("H", "HE", residue, residueID);
                Atom HH11 = new Atom("H", "HH11", residue, residueID);
                Atom HH12 = new Atom("H", "HH12", residue, residueID);
                Atom HH21 = new Atom("H", "HH21", residue, residueID);
                Atom HH22 = new Atom("H", "HH22", residue, residueID);

                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                addAtom(HG2, CG);
                addAtom(HG3, CG);
                addAtom(HD2, CD);
                addAtom(HD3, CD);
                addAtom(HE, NE);
                addAtom(HH11, NH1);
                addAtom(HH12, NH1);
                addAtom(HH21, NH2);
                addAtom(HH22, NH2);
            }
        } else if (residue.equals("ASN")) { //asparagine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG = new Atom("C", "CG", residue, residueID);
            Atom OD1 = new Atom("O", "OD1", residue, residueID);
            Atom ND2 = new Atom("N", "ND2", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG, CB);
            addAtom(OD1, CG);
            addAtom(ND2, CG);

            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB2 = new Atom("H", "HB2", residue, residueID);
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                Atom HD21 = new Atom("H", "HD21", residue, residueID);
                Atom HD22 = new Atom("H", "HD22", residue, residueID);

                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                addAtom(HD21, ND2);
                addAtom(HD22, ND2);
            }
        } else if (residue.equals("ASP")) { //aspartic
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG = new Atom("C", "CG", residue, residueID);
            Atom OD1 = new Atom("O", "OD1", residue, residueID);
            Atom OD2 = new Atom("O", "OD2", residue, residueID); //amber type becomes OH if it has a hydrogen attached

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG, CB);
            addAtom(OD1, CG);
            addAtom(OD2, CG);

            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB2 = new Atom("H", "HB2", residue, residueID);
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                //Atom HD2 = new Atom("H","HD2", residue, residueID);

                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                //addAtom(HD2, OD2, CG);
            }
        } else if (residue.equals("CYS")) { //cysteine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom SG = new Atom("S", "SG", residue, residueID);

            //Connect side-chain
            addAtom(CB, CA);
            addAtom(SG, CB);

            if (hydrogens) {

                //The hydrogens in the side chain
                Atom HB2 = new Atom("H", "HB2", residue, residueID); //for some reason this starts at 2?
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                //Atom HG = new Atom("H","HG", "HS", residue, residueID);

                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                //addAtom(HG, SG); //no HG in cysteine?
            }
        } else if (residue.equals("GLU")) { //glutamic acid
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG = new Atom("C", "CG", residue, residueID);
            Atom CD = new Atom("C", "CD", residue, residueID);
            Atom OE1 = new Atom("O", "OE1", residue, residueID);
            Atom OE2 = new Atom("O", "OE2", residue, residueID); //amber type becomes OH if it has a hyrdogen attached

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG, CB);
            addAtom(CD, CG);
            addAtom(OE1, CD);
            addAtom(OE2, CD);

            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB2 = new Atom("H", "HB2", residue, residueID);
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                Atom HG2 = new Atom("H", "HG2", residue, residueID);
                Atom HG3 = new Atom("H", "HG3", residue, residueID);
                //Atom HE2 = new Atom("H","HD2", residue, residueID);
                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                addAtom(HG2, CG);
                addAtom(HG3, CG);
                //addAtom(HE2, OE2, CD);
            }
        } else if (residue.equals("GLN")) { //glutamine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG = new Atom("C", "CG", residue, residueID);
            Atom CD = new Atom("C", "CD", residue, residueID);
            Atom OE1 = new Atom("O", "OE1", residue, residueID);
            Atom NE2 = new Atom("N", "NE2", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG, CB);
            addAtom(CD, CG);
            addAtom(OE1, CD);
            addAtom(NE2, CD);

            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB2 = new Atom("H", "HB2", residue, residueID);
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                Atom HG2 = new Atom("H", "HG2", residue, residueID);
                Atom HG3 = new Atom("H", "HG3", residue, residueID);
                Atom HE21 = new Atom("H", "HE21", residue, residueID);
                Atom HE22 = new Atom("H", "HE22", residue, residueID);

                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                addAtom(HG2, CG);
                addAtom(HG3, CG);
                addAtom(HE21, NE2);
                addAtom(HE22, NE2);
            }
        } else if (residue.equals("GLY")) { //glycine
            //we must rename HA to HA2
            Atom HA = this.getAtom("HA", residueID);
            HA.setName("HA2");

            //Add new atoms
            Atom HA3 = new Atom("H", "HA3", residue, residueID);

            if (hydrogens)
                addAtom(HA3, CA);
        } else if (residue.equals("HIS")) { //histidine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG = new Atom("C", "CG",  residue, residueID);
            Atom ND1 = new Atom("N", "ND1", residue, residueID);
            Atom CE1 = new Atom("C", "CE1", residue, residueID);
            Atom CD2 = new Atom("C", "CD2", residue, residueID);
            Atom NE2 = new Atom("N", "NE2", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG, CB);

            //Cycle forwards
            addAtom(ND1, CG);
            addAtom(CE1, ND1);
            addAtom(NE2, CE1);
            addAtom(CD2, NE2);
            addAtom(CG, CD2);

            //Cycle backwards
            addAtom(CD2, CG);
            addAtom(NE2, CD2);
            addAtom(CE1, NE2);
            addAtom(ND1, CE1);
            addAtom(CG, ND1);

            if (hydrogens) {
                //Hydrogen atoms
                Atom HB2 = new Atom("H", "HB2", residue, residueID);
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                Atom HD1 = new Atom("H", "HD1", residue, residueID);
                Atom HD2 = new Atom("H", "HD2", residue, residueID);
                Atom HE1 = new Atom("H", "HE1", residue, residueID);
                //Atom HE2 = new Atom("H", "HE2", residue, residueID);

                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                addAtom(HD1, ND1);
                addAtom(HD2, CD2);
                addAtom(HE1, CE1); //because of the cycle, HE1 is connected to both nitrogens!
                //addAtom(HE2, NE2, CE1, CD2); //because of the cycle, HE2 is connected to both CE1 and CD2!
            }
        } else if (residue.equals("ILE")) { //isoleucine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG1 = new Atom("C", "CG1", residue, residueID);
            Atom CG2 = new Atom("C", "CG2", residue, residueID);
            Atom CD1 = new Atom("C", "CD1", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG1, CB);
            addAtom(CG2, CB);
            addAtom(CD1, CG1);

            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB = new Atom("H", "HB", residue, residueID);
                Atom HG12 = new Atom("H", "HG12", residue, residueID);
                Atom HG13 = new Atom("H", "HG13", residue, residueID);
                Atom HG21 = new Atom("H", "HG21", residue, residueID);
                Atom HG22 = new Atom("H", "HG22", residue, residueID);
                Atom HG23 = new Atom("H", "HG23", residue, residueID);
                Atom HD11 = new Atom("H", "HD11", residue, residueID);
                Atom HD12 = new Atom("H", "HD12", residue, residueID);
                Atom HD13 = new Atom("H", "HD13", residue, residueID);

                //Connect hydrogens
                addAtom(HB, CB);
                addAtom(HG12, CG1);
                addAtom(HG13, CG1);
                addAtom(HG21, CG2);
                addAtom(HG22, CG2);
                addAtom(HG23, CG2);
                addAtom(HD11, CD1);
                addAtom(HD12, CD1);
                addAtom(HD13, CD1);
            }
        } else if (residue.equals("LEU")) { //leucine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG = new Atom("C", "CG", residue, residueID);
            Atom CD1 = new Atom("C", "CD1", residue, residueID);
            Atom CD2 = new Atom("C", "CD2", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG, CB);
            addAtom(CD1, CG);
            addAtom(CD2, CG);

            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB2 = new Atom("H", "HB2", residue, residueID); //for some reason this starts at 2?
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                Atom HG = new Atom("H", "HG", residue, residueID);
                Atom HD11 = new Atom("H", "HD11", residue, residueID);
                Atom HD12 = new Atom("H", "HD12", residue, residueID);
                Atom HD13 = new Atom("H", "HD13", residue, residueID);
                Atom HD21 = new Atom("H", "HD21", residue, residueID);
                Atom HD22 = new Atom("H", "HD22", residue, residueID);
                Atom HD23 = new Atom("H", "HD23", residue, residueID);

                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                addAtom(HG, CG);
                addAtom(HD11, CD1);
                addAtom(HD12, CD1);
                addAtom(HD13, CD1);
                addAtom(HD21, CD2);
                addAtom(HD22, CD2);
                addAtom(HD23, CD2);
            }
        } else if (residue.equals("LYS")) { //lysine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG = new Atom("C", "CG", residue, residueID);
            Atom CD = new Atom("C", "CD", residue, residueID);
            Atom CE = new Atom("C", "CE", residue, residueID);
            Atom NZ = new Atom("N", "NZ", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG, CB);
            addAtom(CD, CG);
            addAtom(CE, CD);
            addAtom(NZ, CE);

            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB2 = new Atom("H", "HB2", residue, residueID);
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                Atom HG2 = new Atom("H", "HG2", residue, residueID);
                Atom HG3 = new Atom("H", "HG3", residue, residueID);
                Atom HD2 = new Atom("H", "HD2", residue, residueID);
                Atom HD3 = new Atom("H", "HD3", residue, residueID);
                Atom HE2 = new Atom("H", "HE2", residue, residueID);
                Atom HE3 = new Atom("H", "HE3", residue, residueID);
                Atom HZ1 = new Atom("H", "HZ1", residue, residueID);
                Atom HZ2 = new Atom("H", "HZ2", residue, residueID);
                Atom HZ3 = new Atom("H", "HZ3", residue, residueID);

                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                addAtom(HG2, CG);
                addAtom(HG3, CG);
                addAtom(HD2, CD);
                addAtom(HD3, CD);
                addAtom(HE2, CE);
                addAtom(HE3, CE);
                addAtom(HZ1, NZ);
                addAtom(HZ2, NZ);
                addAtom(HZ3, NZ);
            }
        } else if (residue.equals("MET")) { //methionine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG = new Atom("C", "CG", residue, residueID);
            Atom SD = new Atom("S", "SD", residue, residueID);
            Atom CE = new Atom("C", "CE", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG, CB);
            addAtom(SD, CG);
            addAtom(CE, SD);

            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB2 = new Atom("H", "HB2", residue, residueID);
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                Atom HG2 = new Atom("H", "HG2", residue, residueID);
                Atom HG3 = new Atom("H", "HG3", residue, residueID);
                Atom HE1 = new Atom("H", "HE1", residue, residueID);
                Atom HE2 = new Atom("H", "HE2", residue, residueID);
                Atom HE3 = new Atom("H", "HE3", residue, residueID);

                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                addAtom(HG2, CG);
                addAtom(HG3, CG);
                addAtom(HE1, CE);
                addAtom(HE2, CE);
                addAtom(HE3, CE);
            }
        } else if (residue.equals("PHE")) { //phenylalanine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG = new Atom("C", "CG", residue, residueID);
            Atom CD1 = new Atom("C", "CD1", residue, residueID);
            Atom CD2 = new Atom("C", "CD2", residue, residueID);
            Atom CE1 = new Atom("C", "CE1", residue, residueID);
            Atom CE2 = new Atom("C", "CE2", residue, residueID);
            Atom CZ = new Atom("C", "CZ", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG, CB);
            addAtom(CD1, CG);
            addAtom(CE1, CD1);
            addAtom(CZ, CE1);
            addAtom(CE2, CZ);
            addAtom(CD2, CE2);
            addAtom(CG, CD2);

            //now do the cycle backwards
            addAtom(CD2, CG);
            addAtom(CE2, CD2);
            addAtom(CZ, CE2);
            addAtom(CE1, CZ);
            addAtom(CD1, CE1);
            addAtom(CG, CD1);

            if (hydrogens) {
                //Hydrogen atoms
                Atom HB2 = new Atom("H", "HB2", residue, residueID);
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                Atom HD1 = new Atom("H", "HD1", residue, residueID);
                Atom HD2 = new Atom("H", "HD2", residue, residueID);
                Atom HE1 = new Atom("H", "HE1", residue, residueID);
                Atom HE2 = new Atom("H", "HE2", residue, residueID);
                Atom HZ = new Atom("H", "HZ", residue, residueID); //attached to CZ

                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                addAtom(HD1, CD1);
                addAtom(HD2, CD2);
                addAtom(HE1, CE1);
                addAtom(HE2, CE2);
                addAtom(HZ, CZ); //connected to CE1 and CE2 because cycle
            }
        } else if (residue.equals("PRO")) { //proline

            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG = new Atom("C", "CG", residue, residueID);
            Atom CD = new Atom("C", "CD", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG, CB);
            addAtom(CD, CG);

            //Complete the cycle
            addAtom(CD, N);
            addAtom(CG, CD);
            addAtom(CB, CG);
            addAtom(CA, CB);
            if (hydrogens) {
                //REMOVE H3 from the amino acid (CD is connected to N, so one hydrogen disappears)
                Atom H2 = this.getAtom("H2", residueID);
                Atom H3 = this.getAtom("H3", residueID);

                if (H3 == null) {
                    //NOTE: when proline is not the N-terminus, the N is not bonded to any H-atoms
                    Atom H = this.getAtom("H", residueID);
                    this.removeVertex(H);
                } else {
                    // proline.removeVertex(H2);
                    this.removeVertex(H3);
                }


                //The hydrogens in the side chain
                Atom HB2 = new Atom("H", "HB2", residue, residueID);
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                Atom HG2 = new Atom("H", "HG2", residue, residueID);
                Atom HG3 = new Atom("H", "HG3", residue, residueID);
                Atom HD2 = new Atom("H", "HD2", residue, residueID);
                Atom HD3 = new Atom("H", "HD3", residue, residueID);

                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                addAtom(HG2, CG);
                addAtom(HG3, CG);
                addAtom(HD2, CD);
                addAtom(HD3, CD);
            }
        } else if (residue.equals("SER")) { //serine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom OG = new Atom("O", "OG", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(OG, CB);

            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB2 = new Atom("H", "HB2", residue, residueID); //for some reason this starts at 2?
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                Atom HG = new Atom("H", "HG", residue, residueID);


                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                addAtom(HG, OG);
            }
        } else if (residue.equals("THR")) { //threonine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom OG1 = new Atom("O", "OG1", residue, residueID);
            Atom CG2 = new Atom("C", "CG2", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(OG1, CB);
            addAtom(CG2, CB);

            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB = new Atom("H", "HB", residue, residueID);
                Atom HG1 = new Atom("H", "HG1", residue, residueID);
                Atom HG21 = new Atom("H", "HG21", residue, residueID);
                Atom HG22 = new Atom("H", "HG22", residue, residueID);
                Atom HG23 = new Atom("H", "HG23", residue, residueID);

                //Connect hydrogens
                addAtom(HB, CB);
                addAtom(HG1, OG1);
                addAtom(HG21, CG2);
                addAtom(HG22, CG2);
                addAtom(HG23, CG2);
            }
        } else if (residue.equals("TRP")) { //tryptophan
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG = new Atom("C", "CG", residue, residueID);
            Atom CD1 = new Atom("C", "CD1", residue, residueID);
            Atom CD2 = new Atom("C", "CD2", residue, residueID);
            Atom NE1 = new Atom("N", "NE1", residue, residueID);
            Atom CE2 = new Atom("C", "CE2", residue, residueID);
            Atom CE3 = new Atom("C", "CE3", residue, residueID);
            Atom CZ2 = new Atom("C", "CZ2", residue, residueID);
            Atom CZ3 = new Atom("C", "CZ3", residue, residueID);
            Atom CH2 = new Atom("C", "CH2", residue, residueID);

            //Connect side chain
            //First cycle forwards
            addAtom(CB, CA);
            addAtom(CG, CB);
            addAtom(CD1, CG);
            addAtom(NE1, CD1);
            addAtom(CE2, NE1);
            addAtom(CD2, CE2);
            addAtom(CG, CD2);

            //Now first cycle backwards
            addAtom(CD2, CG);
            addAtom(CE2, CD2);
            addAtom(NE1, CE2);
            addAtom(CD1, NE1);
            addAtom(CG, CD1);

            //second cycle forwards
            addAtom(CZ2, CE2);
            addAtom(CH2, CZ2);
            addAtom(CZ3, CH2);
            addAtom(CE3, CZ3);
            addAtom(CD2, CE3);

            //second cycle backwards
            addAtom(CE3, CD2);
            addAtom(CZ3, CE3);
            addAtom(CH2, CZ3);
            addAtom(CZ2, CH2);
            addAtom(CE2, CZ2);
            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB2 = new Atom("H", "HB2", residue, residueID);
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                Atom HD1 = new Atom("H", "HD1", residue, residueID);
                Atom HE1 = new Atom("H", "HE1", residue, residueID);
                Atom HE3 = new Atom("H", "HE3", residue, residueID);
                Atom HZ2 = new Atom("H", "HZ2", residue, residueID);
                Atom HZ3 = new Atom("H", "HZ3", residue, residueID);
                Atom HH2 = new Atom("H", "HH2", residue, residueID);

                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                addAtom(HD1, CD1);
                addAtom(HE1, NE1);
                addAtom(HE3, CE3);
                addAtom(HZ2, CZ2);
                addAtom(HZ3, CZ3);
                addAtom(HH2, CH2);
            }
        } else if (residue.equals("TYR")) { //tyrosine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG = new Atom("C", "CG", residue, residueID);
            Atom CD1 = new Atom("C", "CD1", residue, residueID);
            Atom CD2 = new Atom("C", "CD2", residue, residueID);
            Atom CE1 = new Atom("C", "CE1", residue, residueID);
            Atom CE2 = new Atom("C", "CE2", residue, residueID);
            Atom CZ = new Atom("C", "CZ", residue, residueID);
            Atom OH = new Atom("O", "OH", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG, CB);
            addAtom(CD1, CG);
            addAtom(CE1, CD1);
            addAtom(CZ, CE1);
            addAtom(CE2, CZ);
            addAtom(CD2, CE2);
            addAtom(CG, CD2);

            //now do the cycle backwards
            addAtom(CD2, CG);
            addAtom(CE2, CD2);
            addAtom(CZ, CE2);
            addAtom(CE1, CZ);
            addAtom(CD1, CE1);
            addAtom(CG, CD1);

            //Connect OH
            addAtom(OH, CZ);

            if (hydrogens) {
                //Hydrogen atoms
                Atom HB2 = new Atom("H", "HB2", residue, residueID);
                Atom HB3 = new Atom("H", "HB3", residue, residueID);
                Atom HD1 = new Atom("H", "HD1", residue, residueID);
                Atom HD2 = new Atom("H", "HD2", residue, residueID);
                Atom HE1 = new Atom("H", "HE1", residue, residueID);
                Atom HE2 = new Atom("H", "HE2", residue, residueID);
                Atom HH = new Atom("H", "HH", residue, residueID); //attached to O

                //Connect hydrogens
                addAtom(HB2, CB);
                addAtom(HB3, CB);
                addAtom(HD1, CD1);
                addAtom(HD2, CD2);
                addAtom(HE1, CE1);
                addAtom(HE2, CE2);
                addAtom(HH, OH);
            }
        } else if (residue.equals("VAL")) { //valine
            //Add new atoms
            Atom CB = new Atom("C", "CB", residue, residueID);
            Atom CG1 = new Atom("C", "CG1", residue, residueID);
            Atom CG2 = new Atom("C", "CG2", residue, residueID);

            //Connect side chain
            addAtom(CB, CA);
            addAtom(CG1, CB);
            addAtom(CG2, CB);

            if (hydrogens) {
                //The hydrogens in the side chain
                Atom HB = new Atom("H", "HB", residue, residueID);
                Atom HG11 = new Atom("H", "HG11", residue, residueID);
                Atom HG12 = new Atom("H", "HG12", residue, residueID);
                Atom HG13 = new Atom("H", "HG13", residue, residueID);
                Atom HG21 = new Atom("H", "HG21", residue, residueID);
                Atom HG22 = new Atom("H", "HG22", residue, residueID);
                Atom HG23 = new Atom("H", "HG23", residue, residueID);

                //Connect hydrogens
                addAtom(HB, CB);
                addAtom(HG11, CG1);
                addAtom(HG12, CG1);
                addAtom(HG13, CG1);
                addAtom(HG21, CG2);
                addAtom(HG22, CG2);
                addAtom(HG23, CG2);
            }
        }
    }

    private int numberOfResidues(){
        int maxRes = 0;
        for(Atom a : this.vertexList())
            maxRes = Math.max(maxRes, a.getResidueID());
        return maxRes;
    }

    /**
     * Saves the MDGP to a file, suitable for either MDjeep or ibp-ng
     * @param format    the desired format of the file (only "mdjeep" and "dmdgp") supported
     * @param filename  the name of the file
     * @param iterator  the ordering of the Atoms
     * @return  true if saving the file was a success, false otherwise
     */
    public boolean save(String format, String filename, DGP<Atom>.GIterator iterator, int precision) {
        Locale.setDefault(Locale.US);
        try {
            File outputFile = new File(filename);
            outputFile.createNewFile();
            FileWriter myWriter = new FileWriter(outputFile);

            if (format.equalsIgnoreCase("mdjeep")) {
                List<Atom> vertices = iterator.getOrder(); //vertices in the right vertex order!
                for (int i = 0; i < vertices.size(); i++) {
                    for (int j = i + 1; j < vertices.size(); j++) {
                        if (contains(vertices.get(i), vertices.get(j))) {
                            Atom A = vertices.get(i);
                            Atom B = vertices.get(j);

                            myWriter.write(String.format("%5s", i + 1)); //1-indexed
                            myWriter.write(String.format("%5s", j + 1)); //1-indexed
                            myWriter.write(" ");
                            myWriter.write(String.format("%5s", A.getResidueID()));
                            myWriter.write(String.format("%5s", B.getResidueID()));
                            myWriter.write(" ");

                            Distance<Atom> d = this.getDistance(A, B);

                            //variable decimal precision:
                            int dLen = 6 + precision;
                            String dFormat = "%" + dLen + "." + precision + "f";

                            if (d.hasBounds()) {
                                myWriter.write(String.format(dFormat, d.getLowerBound()));
                                myWriter.write(String.format(dFormat, d.getUpperBound()));
                            } else if (d.hasExpectedValue()) {
                                myWriter.write(String.format(dFormat, d.getExpectedValue()));
                                myWriter.write(String.format(dFormat, d.getExpectedValue()));
                            } else {
                                throw new Exception("Edge between " + A.getName() + " and " + B.getName() + " does not have bounds, this should not happen!");
                            }

                            myWriter.write("  ");
                            myWriter.write(String.format("%-5s", A.getName()));
                            myWriter.write(String.format("%-5s", B.getName()));
                            myWriter.write(" ");
                            myWriter.write(String.format("%-4s", A.getResidueName()));
                            myWriter.write(" ");
                            myWriter.write(String.format("%-4s", B.getResidueName()));
                            myWriter.write("\n");
                        }
                    }
                }
                myWriter.close();
                return true;
            }
            if (format.equalsIgnoreCase("dmdgp")) {
                myWriter.write("# " + filename + "\n"); //1-indexed
                myWriter.write("# generated by BP_ProteinFileReader\n"); //1-indexed
                myWriter.write("\n");
                myWriter.write("# sequence\n");

                //writing sequence data
                int res = 1;
                int numberOfRes = this.numberOfResidues();
                while(res <= numberOfRes){
                    myWriter.write("#");

                    //15 residues per line...
                    int limit = res + 15;
                    while(res < limit && res <= numberOfRes){
                        //all residues should have N
                        Atom N = this.getAtom("N",res);

                        String resName = N.getResidueName();
                        myWriter.write(" "  + resName);
                        res++;
                    }
                    myWriter.write("\n");
                }

                //writing vertices
                myWriter.write("# vertices: " + this.numberOfVertices() + "\n");
                myWriter.write("begin vertices\n");

                List<Atom> vertices = iterator.getOrder(); //vertices in the right vertex order!
                int i = 1;
                for(Atom a : vertices){
                    myWriter.write(String.format("%-4d  *   *   *   # %s%-4d %-4s (%s)\n", i, a.getResidueName(),a.getResidueID(), a.getName(), a.getTypeName()));
                    i++;
                }
                //myWriter.write("\n");

                int exact = 0;
                int interval = 0;

                for(Distance<Atom> dist : this.distanceSet())
                    if(dist.hasExpectedValue())
                        exact++;
                    else
                        interval++;
                myWriter.write("end vertices\n\n");

                //writing edges
                myWriter.write(String.format("%-18s%-10d\n", "# exact edges: ",exact));
                myWriter.write(String.format("%-18s%-10d\n", "# interval edges: ",interval));

                myWriter.write("begin edges\n");
                for (i = 0; i < vertices.size(); i++) {
                    for (int j = i + 1; j < vertices.size(); j++) {
                        if (contains(vertices.get(i), vertices.get(j))) {
                            Atom a = vertices.get(i);
                            Atom b = vertices.get(j);

                            Distance<Atom> d = this.getDistance(a, b);

                            if(d.hasExpectedValue())
                                myWriter.write(String.format("%-4d%-4dD %11.6f             # %s%-4d %-4s -- %s%-4d %-4s\n",
                                        i+1,
                                        j+1,
                                        d.getExpectedValue(),
                                        a.getResidueName(),
                                        a.getResidueID(),
                                        a.getName(),
                                        b.getResidueName(),
                                        b.getResidueID(),
                                        b.getName()));
                            else
                                myWriter.write(String.format("%-4d%-4dI %11.6f %11.6f # %s%-4d %-4s -- %s%-4d %-4s\n",
                                        i+1,
                                        j+1,
                                        d.getLowerBound(),
                                        d.getUpperBound(),
                                        a.getResidueName(),
                                        a.getResidueID(),
                                        a.getName(),
                                        b.getResidueName(),
                                        b.getResidueID(),
                                        b.getName()));
                        }
                    }
                }
                myWriter.write("end edges\n\n");

                myWriter.write("# atoms: " + this.numberOfVertices() + "\n");
                myWriter.write("begin atom_names\n");
                Set<String> uniqueAtomNames = new HashSet<>();
                for(Atom a : vertices){
                    boolean isNew = uniqueAtomNames.add(a.getName());
                    if(isNew){
                        myWriter.write(String.format("%-5s",a.getName()));
                        for(i = 0; i < vertices.size(); i++)
                            if(vertices.get(i).getName().equals(a.getName()))
                                myWriter.write(String.format("%-4d",i));
                        myWriter.write("\n");
                    }
                }
                myWriter.write("end atom_names\n\n");

                myWriter.write("# residues: " + numberOfRes + "\n");
                myWriter.write("begin residues\n");
                Set<String> uniqueResNames = new HashSet<>();
                for(res = 1; res <= numberOfRes; res++) {
                    Atom N = this.getAtom("N",res);
                    String resName = N.getResidueName();
                    boolean isNew = uniqueResNames.add(resName);
                    if(isNew){
                        myWriter.write(String.format("%-5s",resName));
                        for(i = 1; i <= numberOfRes; i++){
                            Atom N_2 = this.getAtom("N",i);
                            if(N_2.getResidueName().equals(resName))
                                myWriter.write(String.format("%-4d",i));
                        }
                        myWriter.write("\n");
                    }

                }
                myWriter.write("end residues\n\n");

                int impropers = 0;
                int dihdedrals = 0;
                List<Pair<List<Atom>, double[]>> exactDihed = new ArrayList<>();
                List<Pair<List<Atom>, double[]>> intervalDihed = new ArrayList<>();
                Set<String> dihedralAtomNames = Set.of("CA","C","N");
                for(Atom a : this.dihedralAngles.keySet()){
                    Pair<List<Atom>, double[]> dihed = this.dihedralAngles.get(a);
                    double[] angles = dihed._2();
                    List<Atom> atoms = new ArrayList<>(dihed._1());
                    atoms.add(0, a);

                    Pair<List<Atom>, double[]> toAdd = new Pair<>(atoms, angles);
                    if(angles[0] == angles[1])
                        exactDihed.add(toAdd);
                    else
                        intervalDihed.add(toAdd);

                    if(dihedralAtomNames.contains(a.getName()) && a.getName().equals(dihed._1().get(2).getName()))
                        dihdedrals++;
                    else
                        impropers++;
                }
                myWriter.write("# impropers: " + impropers + "\n");
                myWriter.write("# dihedrals: " + dihdedrals + "\n");
                myWriter.write("begin dihedral_angles\n");
                for(Pair<List<Atom>, double[]> dihed : exactDihed) {
                    List<Atom> atoms = dihed._1();
                    Collections.reverse(atoms);
                    myWriter.write(String.format("%-4d%-4d%-4d%-4dD %11.6f\n",
                            vertices.indexOf(atoms.get(0)) + 1,
                            vertices.indexOf(atoms.get(1)) + 1,
                            vertices.indexOf(atoms.get(2)) + 1,
                            vertices.indexOf(atoms.get(3)) + 1,
                            Math.toDegrees(dihed._2()[0])));
                }

                for(Pair<List<Atom>, double[]> dihed : intervalDihed) {
                    List<Atom> atoms = dihed._1();
                    Collections.reverse(atoms);
                    myWriter.write(String.format("%-4d%-4d%-4d%-4dI %11.6f %11.6f\n",
                            vertices.indexOf(atoms.get(0)) + 1,
                            vertices.indexOf(atoms.get(1)) + 1,
                            vertices.indexOf(atoms.get(2)) + 1,
                            vertices.indexOf(atoms.get(3)) + 1,
                            Math.toDegrees(dihed._2()[0]),
                            Math.toDegrees(dihed._2()[1])));
                }
                myWriter.write("end dihedral_angles\n\n");
                myWriter.close();
                return true;
            }
            else{
                throw new IllegalArgumentException("The output format must equal dmdgp or mdjeep!");
            }
        } catch (Exception e) {
            e.printStackTrace();
            return false;
        }
    }
}

