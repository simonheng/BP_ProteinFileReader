/**
 * The main class of the tool: input/output handling
 *
 * @author   Antonio Mucherino
 * @author   Simon Hengeveld
 * @since    November 2nd, 2021
 * @version  January, 2024
 * @see      Distance
 * @see      Cartesian
 * @package  ProteinFileReader
 */

import java.util.Random;

public class Main {
    public static void main(String[] args)
    {
        //args[0] decides the source: ("nmr", "pdb", "dmdgp", "list")
        String source = args[0];

        String inputFile = null;
        String outputFile = null;
        String ffName = "charmm";

        //try to locate input file
        for (int i = 1; i < args.length; i++) {
            if ("--ff".equals(args[i])) {
                if (i + 1 < args.length) {
                    ffName = args[++i];
                } else {
                    System.err.println("Error: missing forcefield name");
                    System.exit(1);
                }
            }
            if ("--input".equals(args[i])) {
                if (i + 1 < args.length) {
                    inputFile = args[++i];
                } else {
                    System.err.println("Error: missing input file");
                    System.exit(1);
                }
            }
            else if ("--output".equals(args[i])) {
                if (i + 1 < args.length) {
                    outputFile = args[++i];
                } else {
                    System.err.println("Error: missing output file name");
                    System.exit(1);
                }
            }
        }
        if(inputFile == null){
            System.err.println("Error: no input file specified!");
            System.exit(1);
        }
        if(outputFile == null){
            System.err.println("Error: no output file specified!");
            System.exit(1);
        }
        ForceField ff = new ForceField(ffName);

        //initialize variables
        MDGP mdgp = null;
        Revorder<Atom> order = null;
        boolean waals = false;

        if(source.equalsIgnoreCase("nmr")){
            //Note: input files must be located in NMR_files/ folder!
            System.out.println("Parsing " + inputFile  + " as NMR file...");
            NMR nmr = new NMR(inputFile);

            boolean sidechain = false;
            boolean hydrogens = true;
            waals = true;

            for (int i = 1; i < args.length; i++) {
                if ("--input".equals(args[i]) || "--output".equals(args[i]) || "--ff".equals(args[i])) {
                    i++; //skip
                }
                else  if ("--hydro".equals(args[i])) {
                    if (i + 1 < args.length) {
                        hydrogens = args[++i].equalsIgnoreCase("true");
                    } else {
                        System.err.println("Error: missing hydro option value");
                        System.exit(1);
                    }
                }
                else  if ("--waals".equals(args[i])) {
                    if (i + 1 < args.length) {
                        waals = args[++i].equalsIgnoreCase("true");
                    } else {
                        System.err.println("Error: missing waals option value");
                        System.exit(1);
                    }
                }
                else  if ("--sidechain".equals(args[i])) {
                    if (i + 1 < args.length) {
                        sidechain = args[++i].equalsIgnoreCase("true");
                    } else {
                        System.err.println("Error: missing sidechain option value");
                        System.exit(1);
                    }
                }
                else {
                    System.err.println("Error: unrecognized argument: " + args[i]);
                    System.exit(1);
                }
            }
            System.out.println("Constructing MDGP object with parameters: " + "--sidechain: " + sidechain + " --hydrogens: " + hydrogens + " --forcefield: " + ffName + " --waals: " + waals);
            mdgp = new MDGP(nmr, ff, sidechain, hydrogens, false, true, null,'n');

            order = new Revorder<>(mdgp);
        }
        else if(source.equalsIgnoreCase("pdb")){
            int model = 1;
            char chain = 'n';
            boolean bonds = false;
            int angles = 0;
            int seed = new Random().nextInt();
            double range = 10.0;
            waals = true;

            for (int i = 1; i < args.length; i++) {
                if ("--chain".equalsIgnoreCase(args[i])) {
                    if (i + 1 < args.length) {
                        chain = args[++i].charAt(0);
                    } else {
                        System.err.println("Error: missing chain value");
                        System.exit(1);
                    }
                }
                else if ("--model".equalsIgnoreCase(args[i])) {
                    if (i + 1 < args.length) {
                        model = Integer.parseInt(args[++i]);
                    } else {
                        System.err.println("Error: missing moel number");
                        System.exit(1);
                    }
                }
                else if ("--bonds".equalsIgnoreCase(args[i])) {
                    if (i + 1 < args.length) {
                        String val = args[++i];
                        if(val.equalsIgnoreCase("pdb"))
                            bonds = true;
                        else
                            bonds = false;
                    } else {
                        System.err.println("Error: no value for --bonds option");
                        System.exit(1);
                    }
                }
                else if ("--angles".equalsIgnoreCase(args[i])) {
                    if (i + 1 < args.length) {
                        String val = args[++i];
                        if(val.equalsIgnoreCase("pdb"))
                            angles = 2;
                        else if(val.contains("rama"))
                            angles = 1;
                    } else {
                        System.err.println("Error: no value for --bonds option");
                        System.exit(1);
                    }
                }
                else if ("--seed".equalsIgnoreCase(args[i])) {
                    if (i + 1 < args.length) {
                        seed = Integer.parseInt(args[++i]);
                    } else {
                        System.err.println("Error: no value for --seed option");
                        System.exit(1);
                    }
                }
                else if ("--torsrange".equalsIgnoreCase(args[i])) {
                    if (i + 1 < args.length) {
                        range = Double.parseDouble(args[++i]);
                    } else {
                        System.err.println("Error: no value for --seed option");
                        System.exit(1);
                    }
                }
                else if ("--input".equalsIgnoreCase(args[i]) || "--output".equalsIgnoreCase(args[i]) || "--ff".equalsIgnoreCase(args[i])) {
                    i++; //skip
                }
                else {
                    System.err.println("Error: unrecognized argument: " + args[i]);
                    System.exit(1);
                }
            }

            //read the model...
            System.out.println("Reading PDB file: " + inputFile);
            PDB pdb = new PDB(inputFile,model, chain);

            String angStr = "";
            if(angles == 0) angStr = "from FF";
            else if(angles == 1) angStr = "from Ramachandran";
            else angStr = "from PDB";

            System.out.println("Constructing MDGP object with parameters: " + "--bonds: " + (bonds ? "from PDB" : "from FF") + " --angles: " +angStr + " --forcefield: " + ffName + " --torsrange: " + range + " degrees");
            System.out.println("Using random seed: " + seed);
            mdgp = new MDGP(pdb, chain, ff, bonds,angles,true, Math.toRadians(range), seed);
            order = new Revorder<>(mdgp,mdgp.dmdgpOrder());
        }
        else if(source.equalsIgnoreCase("dmdgp")){
            System.out.println("Constructing MDGP object from .dmdgp file: " + inputFile);
            mdgp = new MDGP("DMDGP_files/" + inputFile, "dmdgp", ff);
            order = new Revorder<>(mdgp, mdgp.dmdgpOrder());
        }
        else{
            System.out.println("Source not recognised, must be: nmr, pdb or dmdgp!");
        }

        //optional things...
        if(waals)
            mdgp.addVanDerWaals(true, 0.8, null);

        boolean mdjeep = !outputFile.endsWith(".dmdgp");

        if(mdjeep) {
            System.out.println("Saving " + outputFile + " as MDjeep file!");
            mdgp.save("mdjeep", outputFile,order, 4);
        }
        else{
            System.out.println("Saving " + outputFile + " as .dmdgp (ibp-ng) file!");
            mdgp.save("dmdgp", outputFile,order, 4);
        }
        System.out.println("Done!");
    }
}
