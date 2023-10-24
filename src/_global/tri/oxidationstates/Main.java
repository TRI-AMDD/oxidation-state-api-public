package _global.tri.oxidationstates;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import _global.tri.oxidationstates.calculator.BondValenceCalculator;
import _global.tri.oxidationstates.calculator.OxidationStateCalculator;
import _global.tri.oxidationstates.calculator.OxidationStateSet;
import _global.tri.oxidationstates.calculator.frequency.FrequencyCalculator;
import _global.tri.oxidationstates.calculator.likelihood.LikelihoodCalculator;
import _global.tri.oxidationstates.fitting.OxidationStateData;
import _global.tri.oxidationstates.fitting.ParamOptimizer;
import _global.tri.oxidationstates.fitting.OxidationStateData.Entry;
import _global.tri.oxidationstates.ion.CompositionIonFinder;
import _global.tri.oxidationstates.ion.IIonFinder;
import _global.tri.oxidationstates.ion.IonAssigner;
import _global.tri.oxidationstates.ion.IonFactory;
import _global.tri.oxidationstates.ion.KnownIonFinder;
import _global.tri.oxidationstates.ion.OxideIonFinder;
import _global.tri.oxidationstates.ion.ZintlIonFinder;
import _global.tri.oxidationstates.ion.IonFactory.Ion;
import _global.tri.oxidationstates.ion.IonFactory.IonType;
import _global.tri.oxidationstates.util.Composition;
import _global.tri.oxidationstates.util.CompositionParseException;
import _global.tri.oxidationstates.util.IonTools;
import _global.tri.oxidationstates.util.VaspFileFilter;
import _global.tri.oxidationstates.webapi.PageData;
import _global.tri.oxidationstates.webapi.PotentialMapper;
import _global.tri.oxidationstates.webapi.WebOxidationAnalyzer;
import matsci.Element;
import matsci.Species;
import matsci.engine.minimizer.CGMinimizer;
import matsci.io.app.log.Status;
import matsci.io.clusterexpansion.PRIM;
import matsci.io.structure.CIF;
import matsci.io.structure.XYZFile;
import matsci.io.vasp.POSCAR;
import matsci.structure.PartiallyOccupiedStructure;
import matsci.structure.Structure;
import matsci.util.arrays.ArrayUtils;

/**
 * This is the main class containing the starting points for various oxidaiton analyzer routines.  
 * 
 * @author timmueller
 *
 */
public class Main {

	/**
	 * All input and output files should be contained under this directory
	 */
	public static String ROOT_DIR = "/home/timmueller/Projects/Oxidation state analysis/Public/";
	
	/**
	 *  This contains all of the structures from the broad ICSD data set
	 */
	public static String STRUCT_DIR = ROOT_DIR + "/Structures/";
	
	/**
	 * This is where we store the training data sets
	 */
	public static String TRAINING_DATA_DIR = ROOT_DIR + "/Training_Data/";
	
	/**
	 *  This is where we store the model parameters (oxidation state boundaries)
	 */
	public static String PARAMETER_DIR = ROOT_DIR + "/Parameters/";
	
	/**
	 * This subdirectory contains structures with assigned oxidation states, as assigned
	 * by different methods.  It also should contain files for oxidation state assignments output
	 * from pymatgen and BERTOS.
	 */
	public static String GII_DIR = ROOT_DIR + "/GII_Structures/";
	
	/**
	 * This directory contains the data set from "Novel inorganic crystal 
	 * structures predicted using autonomous simulation agents" (https://doi.org/10.1038/s41597-022-01438-8)
	 */
	public static String CAMD_DIR = ROOT_DIR + "/CAMD_data/";
	
	/**
	 * These represent various version of the data set.  "ChargeBalanced" represents the subset of 
	 * ICSD data with charge-balanced structures.  The reset are increasingly smaller subsets of this 
	 * initial set.  
	 * 
	 * "KnownIons" means only strcutures that contain known elements (e.g. no deuterium)
	 * 
	 * "Integer" means all ICSD oxidation states are integers
	 * 
	 * "NonZero" means no elemental oxidation states are zero; "EAH" means that energy above 
	 * the hull information, from the Materials Project, has been added when it's available; 
	 * 
	 * "WebIons" means that structures that contain polyatomic ions are in the data set twice:
	 *  once with the composition written in terms of elements and once with composiiton written
	 *  in terms of polyatomic ions; 
	 *  
	 *  "NoRareIons" means ions that all entries containing rare ions have been removed, where an
	 *  ion is considered to be rare if it occurs in fewer than 25 entries in the data set.
	 *  
	 *  "NoGIIDecrease" means that all structures for which the model found a set of oxidation states
	 *  that yields a lower global instability index than the ICSD oxidation states were remvoed.  The
	 *  numbers after these entries indicate successive iterations of this cleaning process, where the 
	 *  model was retrained after each cleaning.
	 */
	enum DataType {
		ChargeBalanced,
		ChargeBalanced_KnownIons_Integer_NonZero_EAH,
		ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons,
		ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons_NoRareIons,
		ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons_NoRareIons_NoGIIDecrease,
		ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons_NoRareIons_NoGIIDecrease2,
		ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons_NoRareIons_NoGIIDecrease3,
		ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons_NoRareIons_NoGIIDecrease4,
		ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons_NoRareIons_NoGIIDecrease5,
		ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons_NoRareIons_NoGIIDecrease6,
		ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons_NoRareIons_NoGIIDecrease7,
	}

	/**
	 * This is the set of ions that were removed by the data cleaning process for the data
	 * set used in the paper.  It was manually constructed from the cleaning logs.
	 */
	static String[] REMOVED_IONS = new String[] {
		"(AsO4)5-"," (C2O4)12-"," (C2O4)4-"," (CN)2-"," (CN)5-"," (CN)5+"," (CO)1-"," (CO)1+"," (CO)4-"," (CO3)1-"," (CrO4)4-"," (MoO4)3-"," (MoO4)4-"," (NO2)1+"," (NO2)7-"," (NO3)3-"," (NO3)7-"," (NO3)8-"," (NO3)9-"," (PO3)2-"," (PO4)1-"," (PO4)2-"," (PO4)4-"," (PS4)4-"," (SeO3)3-"," (SeO3)4-"," (SO2)1-"," (SO2)1+"," (SO2)2-"," (SO3)1-"," (SO4)1-"," Ac3+"," Ag3+"," Al1+"," Al2+"," Al3-"," Am2+"," Am4+"," Am5+"," Am6+"," As1+"," As5-"," Au1-"," Au2+"," B1-"," B1+"," B2+"," Bi1-"," Bi1+"," Bi2+"," Bi4+"," Bk2+"," Bk3+"," Bk4+"," Br1+"," Br7+"," Cd1+"," Cf2+"," Cf3+"," Cf4+"," Cl1+"," Cl2+"," Cl4-"," Cl4+"," Cm2+"," Cm4+"," Co1-"," Co1+"," Cr1+"," Cu4+"," Dy2+"," Dy4+"," Er2+"," Er4+"," Es3+"," Eu4+"," Eu6+"," F2+"," Fe1-"," Fe1+"," Fe5+"," Fe6+"," Ga1+"," Ga3-"," Ga4-"," Gd1+"," Gd2+"," Gd4+"," Ge1-"," Ge1+"," Ge3-"," Ge5+"," Hf1+"," Hf2+"," Hf3+"," Ho2+"," Ho4+"," I1+"," I2-"," In1-"," In2-"," In3-"," Ir1+"," Ir2+"," Ir6+"," Kr2+"," La1+"," La4+"," Lu2+"," Lu4+"," Mn1-"," Mn1+"," Mn6+"," Mo1+"," N4-"," N4+"," N5-"," Nb1+"," Nb6+"," Nd4+"," Np2+"," Np7+"," O2+"," Os1+"," Os3+"," Os7+"," Os8+"," P1+"," P2+"," P4-"," P6+"," Pa2+"," Pa3+"," Pa4+"," Pa5+"," Pb1-"," Pb1+"," Pb4-"," Pd1+"," Pd3+"," Pm3+"," Po2-"," Po2+"," Po4+"," Pr2+"," Pt1+"," Pt2-"," Pt3+"," Pt5+"," Pt6+"," Pu2+"," Pu5+"," Pu6+"," Pu7+"," Ra2+"," Re1-"," Re1+"," Re2+"," Rh1-"," Rh1+"," Rh2+"," Rh5+"," Rh6+"," Ru1+"," Ru3-"," Ru7+"," Ru8+"," S1+"," S3-"," Sb1+"," Sb2+"," Sb4+"," Sc1+"," Sc2+"," Sc6+"," Se1+"," Se3+"," Se5+"," Si1+"," Si2-"," Si2+"," Si3-"," Si3+"," Sm4+"," Sn1-"," Sn1+"," Sn2-"," Sn3-"," Sn3+"," Sn4-"," Sr2-"," Ta1+"," Tb1+"," Tb2+"," Tc1+"," Tc2+"," Tc3+"," Tc5+"," Tc6+"," Te1+"," Te2+"," Te3-"," Te3+"," Te5+"," Th2+"," Ti1+"," Ti6+"," Tl1-"," Tl2-"," Tl2+"," Tl5+"," Tm2+"," Tm4+"," U2+"," U3-"," V1+"," W3+"," Xe3+"," Xe4+"," Xe8+"," Y1+"," Y2+"," Y4+"," Yb4+"," Zn1+"," Zn3+"," Zn4+"," Zr1+","Au5+","B2-","Ce2+","Nd2+","Ta2+",
	};


	/**
	 * This is the main entry point for the program.  Other routines are called from here.
	 * 
	 * @param args Command line arguments
	 */
	public static void main(String[] args) {

		/*
		 * This loads up the data set of polyatomic ions used for the paper (and web site).  In
		 * general you should call it before calling any of the below methods.  There's no harm in
		 * calling it twice.
		 */
		//loadWebIons();

		/**
		 * This reads through the structures extracted from the ICSD to create an initial data file
		 */
		//prepareDataFromICSD(DataType.ChargeBalanced);
		
		/**
		 * This searches through the ICSD data to find potential polyatomic ions.
		 */
		//findIons(DataType.ChargeBalanced);
		
		/**
		 * This places the most commonions in groups with structureally similar ions and 
		 * calculates representative ("mean") ions for each group.
		 */
		//groupIons(DataType.ChargeBalanced, 100, 100);
		
		/**
		 * This groups the most common ions from the previous step by oxidation state.
		 */
		//groupIonsByOxidationState(0.003, true, 100);
		
		/**
		 * This extracts the subset of training data that has integer oxidation states with no
		 * zero oxidation states and adds energies above the convex hull as extracted from teh Materials
		 * Project
		 */
		//getInitialTrainingData();
		
		/**
		 * This identifies all structures in the training data that contain common polyatomic ions
		 * and adds entries to the training data set for each of those structures with compositions
		 * written in terms of polyatomic ions.
		 */
		//addWebIonsToData(DataType.ChargeBalanced_KnownIons_Integer_NonZero_EAH, DataType.ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons);

		/**
		 * This removes all structures from the data set that contain rare ions, where rare ions are 
		 * defined as those that appear fewer than 25 times in the data set.
		 */
		//removeRareIons(DataType.ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons, DataType.ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons_NoRareIons);

		/**
		 * This method can be used to train the model on a HPC cluster.  You can uncomment 
		 * the below line and export a runnable jar to create an app that can be called from the 
		 * command line.
		 */
		
		// Use this for testing / debugging
		//DataType dataType = DataType.ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons_NoRareIons;
		//trainModel(new String[] {getDataFileName(dataType), getParamFileName(dataType), getWebIonDirName()});

		// Use this for command line interface, and pass the arguments in for debugging
		trainModel(args); 
		
		/**
		 * This assigns oxidation states to all entries in the training set using 10-fold
		 * cross validation.
		 */
		//getValidationAssignments(DataType.ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons, true);
		
		/**
		 * This assigns oxidation states to all entries in the training set using the model trained 
		 * on that set.
		 */
		//getTrainingAssigments(DataType.ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons, true);
		
		/**
		 * This assigns oxidation states to all entries in the given data set based on the oxidation 
		 * states from the ICSD.
		 */
		//getICSDAssignments(DataType.ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons, true);

		/**
		 * This removes all entries from the first data type for which the GII is lower than the GII calculated
		 * using the ICSD oxidation states.  The remaining entries are stored in the data set for the second
		 * data type.
		 */
		//cleanDataGII(DataType.ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons_NoRareIons, DataType.ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons_NoRareIons_NoGIIDecrease);

		/**
		 * Assigns oxidation states to structures based on the BERTOS assignments
		 */
		//assignBertosStates();
		
		/**
		 * Calculates the global instability indices for all structures assigned oxidation states by pymatgen
		 */
		//calcGIIForPymatgenStructs();

		/**
		 * The rest of the methods perform analysis on the final model, usually to generate data for the paper.
		 * This is the data type for the final model.
		 */
		DataType finalDataType = DataType.ChargeBalanced_KnownIons_Integer_NonZero_EAH_WebIons_NoRareIons_NoGIIDecrease6;

		/**
		 * Compares two different sets of assignments for the given data type and prints summary statistics
		 * of the comparison.
		 */
		//compareAssignments(finalDataType);

		/**
		 * Assigns oxidation states for structures in the CAMD data set using the final model
		 */
		//getCAMDAssignments(finalDataType, true);
		
		/**
		 * Generates data for a histogram for the fraction of materials on the convex hull as a function
		 * of the likelihood score for the CAMD dataset.
		 */
		//getCAMDHullHistogram(finalDataType);
		
		/**
		 * Generates data from the CAMD data set showing the cumulative number of new materials on
		 * the convex hull as compositions are evaluated from highest to lowest likelihood scores.
		 */
		//getCAMDDiscoveryCurve(finalDataType);

		/**
		 * Writes xzy files containing the coordinates of the polyatomic ions used in the paper
		 * and the web app
		 */
		//writeXYZFiles();
		
		/**
		 * Writes a JSON file containing the mean boundary positions for the final model
		 */
	    //writeBoundaryJSON(finalDataType);
	    
	    /**
	     * Writes out all of the included and removed oxidation states in the final model
	     * for each ion type 
	     */
		//printAllKnownStates(finalDataType);

		/**
		 * Prints the full electrochemical series of boundary values for the final model
		 */
		//printElectrochemicalSeries(finalDataType);

		/**
		 * Writes the files that can be fed into a python script for generating those wavy bar plots.
		 * (or the straight versions of them)
		 */
		//writeDataFilesForVisualization(finalDataType, false);
		
		/**
		 * Generate a table like the ones generated by the web API
		 */
		//testWebAPI();
	}
	
	/**
	 * This method is designed for use in a command-line interface to fit the model.
	 * 
	 * @param args The arguments that can be passed to this method.  They will be parsed as follows:
	 * 	args[0]:  The name of the file containing the training data.
	 *	args[1]:  The name of the file to be written containing parameters of the trained model.  The file name should start with "parameters".  A corresponding file containing boundary values will also be written, where the name of the boundaries file will be name of the parameters file with "parameters" replaced with "boundaries".
	 *  args[2]: The name of the directory containing the atomic structures of allowed polyatomic ions.
	 *  args[3] (optional):  This many optimization steps will be run before updating the files "boundaries.txt" and "parameters.txt" and restarting the optimization from the written values.  The default value is 1000.
	 *  args[4] (optional):  The parameter to be multiplied by the sum of the spreads in the minimum and maximum boundary values for regularization.  The default value is 0.
	 *  args[5] (optional):  The number of threads to run simultaneously for the optimization.  The default value is 1.
	 */
	public static void trainModel(String[] args) {
		
		/**
		 * Print a help statement if someone doesn't seem to know what they're doing
		 */
		if (args.length < 3) {
			System.out.println("Allowed command line arguments: ");
			System.out.println("dataFileName:  The name of the file containing the training data.");
			System.out.println("paramFileName:  The name of the file to be written containing parameters of the trained model.  The file name should start with \"parameters\".  A corresponding file containing boundary values will also be written, where the name of the boundaries file will be name of the parameters file with \"parameters\" replaced with \"boundaries\".");
			System.out.println("ionStructureDirName: The name of the directory containing the atomic structures of allowed polyatomic ions.");
			System.out.println("numStepsPerPass (optional):  This many optimization steps will be run before updating the files \"boundaries.txt\" and \"parameters.txt\" and restarting the optimization from the written values.  The default value is 1000.");
			System.out.println("regParam (optional):  The parameter to be multiplied by the sum of the spreads in the minimum and maximum boundary values for regularization.  The default value is 0.");
			System.out.println("numThreads (optional):  The number of threads to run simultaneously for the optimization.  The default value is 1.");
			System.out.println();
			System.out.println("To run this program, provide values for the above variables in order.  Do not include the variable names.");
			System.exit(-1);
		}

		/**
		 * The arguments as described above.
		 */
		String dataFileName = args[0];
		String paramFileName = args[1];
		String ionStructureDirName = args[2];
		int numStepsPerPass = args.length < 4 ? 1000 : Integer.parseInt(args[3]);
		double regParam = args.length < 5 ? 0 : Double.parseDouble(args[4]);
		int numThreads = args.length < 6 ? 1 : Integer.parseInt(args[5]);
		
		// This is necessary because we replaced "parameters" with "boundaries" to create the boundaries file.
		if (! (new File(paramFileName)).getName().startsWith("parameters")) {
			Status.error("Data file name must begin with the word \"parameters\".");
		}

		// Set the level of logging
		Status.includeTime(true);
		Status.setLogLevelDetail();

		// Load the known polyatomic ions
		Status.flow("Reading ion information from " + ionStructureDirName);
		IonFactory.loadPolyatomicIons(ionStructureDirName);

		// Read the training data
		Status.flow("Loading oxidation state data from " + dataFileName);
		OxidationStateData data = new OxidationStateData(dataFileName, STRUCT_DIR);

		// Train the model
		fitParameters(data, numStepsPerPass, regParam, paramFileName, numThreads);
	}
	
	/**
	 * Fits the model parameters for all ions contained in the training data, randomly initializing the parameters to values between zero and 1
	 * 
	 * @param data The training data
	 * @param maxIterationsPerPass After this many optimization steps have been performed, the code will write the current parameters to a file and restart optimization from those parameters.
	 * @param regParameter The regularization parameter to determine how much weight to give the regularization term.
	 * @param paramFileName The name of the output file containing the parameters.  The file name should start with "parameters".  A separate file containing the boundary values will also be 
	 * written, where the name of the boundary file will be the same as the name of the parameters file, with the word "parameters" replaced by the word "boundaries".
	 * @param numThreads The number of threads to use when training the model using parallel processing.
	 */
	public static void fitParameters(OxidationStateData data, int maxIterationsPerPass, double regParameter, String paramFileName, int numThreads) {
		
		// Create a new calculator based on all the oxidation states contained in the training data
		HashMap<String, int[]> oxidationStates = data.getKnownOxidationStates();
		LikelihoodCalculator calculator = new LikelihoodCalculator(oxidationStates);

		// Initialize the parameters randomly
		Random random = new Random();
		double[] initParams = random.doubles(calculator.numParameters()).toArray();
		calculator = calculator.setParameters(initParams);	
		
		// Fit the parameters
		fitParameters(data, calculator, maxIterationsPerPass, regParameter, paramFileName, numThreads);
		
	}
	
	/**
	 * Fits the model parameters for all ions contained in the training data, randomly initializing the parameters to values between zero and 1
	 * 
	 * @param data The training data
	 * @param initCalculator The Likelihood calculator containing an initial version of the model
	 * @param maxIterationsPerPass After this many optimization steps have been performed, the code will write the current parameters to a file and restart optimization from those parameters.
	 * @param regParameter The regularization parameter to determine how much weight to give the regularization term.
	 * @param paramFileName The name of the output file containing the parameters.  The file name should start with "parameters".  A separate file containing the boundary values will also be 
	 * written, where the name of the boundary file will be the same as the name of the parameters file, with the word "parameters" replaced by the word "boundaries".
	 * @param numThreads The number of threads to use when training the model using parallel processing.
	 */
	public static void fitParameters(OxidationStateData data, LikelihoodCalculator initCalculator, int maxIterationsPerPass, double regParameter, String paramFileName, int numThreads) {

		Status.flow("Running on " + numThreads + " threads.");
		
		// Initialize the conjugate gradient optimizer
		CGMinimizer optimizer = new CGMinimizer();
		optimizer.setMaxAllowedIterations(maxIterationsPerPass);
		optimizer.setMinAllowedStep(1E-6); // We want to hit gradient convergence, not step size

		// Initialize the parameter optimizer
		ParamOptimizer paramOptimizer = new ParamOptimizer(data, numThreads);
		paramOptimizer.setRegularizationParameter(regParameter);
		ParamOptimizer.ParameterState state = paramOptimizer.getParameterState(initCalculator);

		// At some point the number of passes should probably become user-definable, but 20 seems to work well.
		for (int passNum = 0; passNum < 20; passNum++) {
			
			// Report the initial state for this pass
			double score = state.getValue();
			double regularizer = state.getRegularizer();
			Status.flow("Pass number " + passNum + ".  Initial score: " + score + ", Initial score with no regularizer: " + (score - regularizer));
			
			// Optimize the parameters
			optimizer.minimize(state);
			
			// Get the optimized parameters
			state = (ParamOptimizer.ParameterState) optimizer.getMinimumState();
			
			// Write the boundaries to the console
			state.getCalculator().writeBoundaries(System.out);
			
			// Write the boundaries and parameters to files
			String boundaryFileName = paramFileName.replaceFirst("/parameters", "/boundaries");
			Status.flow("Writing parameters to " + paramFileName);
			Status.flow("Writing boundaries to " + boundaryFileName);
			state.getCalculator().writeParameters(paramFileName);
			state.getCalculator().writeBoundaries(boundaryFileName);			
		}
		
		// Report the final state
		double finalScore = state.getValue();
		double finalRegularizer = state.getRegularizer();
		Status.basic("Final score: " + finalScore + ", Final score with no regularizer: " + (finalScore - finalRegularizer));
		
	}

	/**
	 * Reads a directory of CIF files exported from the ICSD and builds an initial training data
	 * file for all ordered, charge-balanced structures.  The CIF files should be in ICSD_CIFS/All CIFs/.
	 * POSCAR-formatted structures will be written to the Structures directory.
	 * 
	 * @param outDataType The dataType for the generated training data
	 */
	public static void prepareDataFromICSD(DataType outDataType) {

		// Where we read the CIF files from 
		File cifDir = new File(ROOT_DIR + "/ICSD_CIFS/All CIFs/");
		
		// Get the name of the training data file to generate
		String outputFileName = getDataFileName(outDataType);
		
		// Log everything to file
		try {
			Status.useOutForErr();
			Status.setOutfile(ROOT_DIR + "/TrainingDataGeneration_" + outDataType + ".log");
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		
		// How far the site occupancy can deviate from an integer before we consider this to be disordered
		double occupancyTolerance = 0.02;
		
		// How far the sum of oxidation states can deviate from 0 before we consider it not to be charge balanced
		double chargeBalanceTolerance = 1E-6;

		// Some values we keep track of for reporting purposes
		int numSuccess = 0;
		int numNotChargeBalanced = 0;
		int numNotOrdered = 0;
		int numNotAllIntegers = 0;
		int numAllZeros = 0;
		int numParseErrors = 0;

		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName));

			for (File cifFile : cifDir.listFiles()) {
				
				// Only process CIF files
				if (!cifFile.isFile() || !cifFile.getName().endsWith(".cif")) {continue;}

				CIF cif = null;
				Structure structure = null;
				try {					
					// Check to make sure we can parse the CIF and this is an ordered structure
					cif = new CIF(cifFile.getAbsolutePath());
					if (cif.numDefiningSites() == 0) {continue;}
					if (!cif.isOrdered(occupancyTolerance)) {
						numNotOrdered++;
						continue;
					}
					structure = new Structure(cif);
				} catch (Exception e) {
					Status.warning("Failed to parse CIF: " + cifFile.getName());
					numParseErrors++;
					continue;
				}
				
				// Strip the suffix to get a unique ID for this entry
				String id = cifFile.getName().replace(".cif", "");

				// Some information for the training data file
				String compositionString = structure.getCompositionString();
				String origins = "ICSD";

				// Gather the oxidation states some reporting information on them
				boolean hasNonZeroOxidationStates = false;
				boolean notAllIntegers = false;
				Species[] allSpecies = structure.getDistinctSpecies();
				String speciesString = "";
				for (Species species : allSpecies) {
					
					// All of the species in the structure, including their oxidation states
					speciesString += " " + species.getSymbol();
					
					// Some reporting on oxidation states.  We screen out these out in another method.
					double oxidationState = species.getOxidationState();
					if (oxidationState != 0) {
						hasNonZeroOxidationStates = true;
					}
					
					// We screen these out in another method.
					if (!(oxidationState == Math.round(oxidationState))) {
						notAllIntegers = true;
					}
				}
				
				// Collect data for reporting.
				if (!hasNonZeroOxidationStates) {
					numAllZeros++;
				}

				if (notAllIntegers) {
					numNotAllIntegers++;
				}

				// Require the structure to be charge balanced
				if (!structure.isChargeBalanced(chargeBalanceTolerance)) {
					numNotChargeBalanced++;
					continue;
				}

				// Write a POSCAR-formatted file of the ordered structure
				PRIM outfile = new PRIM(structure);
				outfile.writeSpeciesWithDescription(false);
				outfile.writeElementSymbolsOnly(false);
				outfile.writeFile(STRUCT_DIR + id + ".vasp");

				// Calculate the global instability index for later use
				double gii = new BondValenceCalculator(structure).getGlobalInstabilityIndex();

				// Write the information to the training data file
				String outputString = numSuccess + ": " + id + "\t" + compositionString + "\t" + speciesString + "\t" + origins + "\t" + "NaN" + "\t" + gii;
				Status.basic(outputString);
				writer.write(outputString + "\n");
				numSuccess++;
			}

			writer.flush();
			writer.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		// Some reporting.
		Status.basic(numNotOrdered + " disordered structures removed.");
		Status.basic(numNotChargeBalanced + " non-charge-balanced structures removed.");
		Status.basic(numAllZeros + " structures with zero oxidation states.");
		Status.basic(numNotAllIntegers + " structures with non-integer oxidation states.");
		Status.basic(numParseErrors + " structures had parse errors.");
		Status.basic(numSuccess + " structures successfully added.");
	}

	/**
	 * Looks through the training set for the given data type to extract polyatomic ions, as 
	 * determined by the OxideIonFinder, ZintlIonFinder, or CompositionIonFinder
	 * 
	 * @param dataType The data set to search through.
	 */
	public static void findIons(DataType dataType) {

		// Load the data set
		OxidationStateData data = getData(dataType);

		// This is where we will store the ions we find.
		String ionDirName = getIonStructureDirName(dataType);
		new File(ionDirName).mkdirs();

		
		for (int entryNum = 0; entryNum < data.numEntries(); entryNum++) {

			// Load the structure in most compact form
			Entry entry = data.getEntry(entryNum);
			String fileName = STRUCT_DIR + entry.getID() + ".vasp";

			POSCAR infile = new POSCAR(fileName);
			Structure structure = new Structure(infile);
			structure = structure.findPrimStructure().getCompactStructure();

			/**
			 *  These are the types of ions we'll look for.  Extra information 
			 *  about these can be found in their class files.
			 */
			IIonFinder[] finders = new IIonFinder[] {
					new ZintlIonFinder(structure),
					new OxideIonFinder(structure),
					new CompositionIonFinder(structure),
			};

			for (IIonFinder ionFinder : finders) {
				if (ionFinder.numFoundPolyatomicIons() > 0) { // If we found any ions...
					
					// Prepare a message for the user
					String message = fileName + ", " + structure.getCompositionString();
					
					for (int ionNum = 0; ionNum < ionFinder.numFoundPolyatomicIons(); ionNum++) {
						Structure ion = ionFinder.getFoundPolyatomicIon(ionNum);
						
						/**
						 * Write the ion to a structure file, organized by the number of periodic 
						 * dimensions and composition
						 */
						int numDimensions = ion.numPeriodicDimensions();
						String composition = ion.getCompositionString("_", false);
						String outDirName = ionDirName + "/" + numDimensions + "_periodic_dimensions/" + composition + "/";
						File outDir = new File(outDirName);
						outDir.mkdirs();
						String ionFileName = outDirName + "/" + outDir.listFiles().length + ".vasp";
						IonTools.writeIonStructure(ion, ionFileName, new File(fileName).getName());
						
						// Add more information to the message
						message += ", " + ionFinder.getFoundPolyatomicIon(ionNum).getCompositionString("_", false);
					}
					
					// Report the message to the user
					Status.basic(message);
				}
			}
		}
	}

	/**
	 * Places ions of the same composition into groups of structurally similar ions, where the 
	 * atomic oxidation states of all ions in the set need to match.
	 * 
	 * @param dataType The data set for which we are grouping the ions
	 * @param minOccurrences The method will ignore any compositions for which the number of found ions is less than this value.
	 * This is useful for screening out rare, very large ions that take a long time to group.
	 * @param numSamples The number of structures to sample when calculating the mean structure for the group
	 */
	public static void groupIons(DataType dataType, int minOccurrences, int numSamples) {

		// The output directory
		String ionTypeDirName = getIonStructureDirName(dataType);
		File ionTypeDir = new File(ionTypeDirName);

		// Loop through all of the folders of discovered ions
		File[] dimensionDirs = ionTypeDir.listFiles();
		for (File dimensionDir : dimensionDirs) {
			if (!dimensionDir.isDirectory()) {continue;}
			File[] compositionDirs = dimensionDir.listFiles();
			for (File compositionDir : compositionDirs) {
				if (!compositionDir.isDirectory()) {continue;}
				
				// This selects everything that ends with ".vasp"
				File[] molecules = compositionDir.listFiles(new VaspFileFilter());
				
				// If there are only a few ions of this type, don't bother
				if (molecules.length < minOccurrences) {
					Status.flow("Skipping " + compositionDir + " because it only has " + molecules.length + " occurrences.");
					continue;
				}

				Status.flow("Grouping ions for " + compositionDir);
				HashSet<HashSet<Structure>> knownMolecules= new HashSet<>();

				for (File moleculeFile : molecules) {
					
					// Read thie structure for each found ion
					if (!moleculeFile.getName().endsWith(".vasp")) {continue;}
					Structure molecule = IonTools.loadIonStructureFromFile(moleculeFile.getAbsolutePath());

					/**
					 * Loop through all sets of grouped ions to see if the ion we just found is a structural match
					 * in any of the sets.  The "matchedMoleculeSet" is the master set containing all ions that are
					 * similar to this one.  If multiple sets of ions contain ions that are simlar to this one, they
					 * are all merged into the matchedMoleculeSet.
					 */
					HashSet<Structure> matchedMoleculeSet = null;
					for (HashSet<Structure> moleculeSet : knownMolecules) {
						for (Structure knownMolecule : moleculeSet) {
							if (IonTools.compareStructures(knownMolecule, molecule, false)) {
								
								// We don't have a master set, so create one
								if (matchedMoleculeSet == null) {
									moleculeSet.add(molecule);
									matchedMoleculeSet = moleculeSet;
								} else { // Merge this set into the master set.
									matchedMoleculeSet.addAll(moleculeSet);
									moleculeSet.clear();
								}
								break;
							}
						}
					}
					
					if (matchedMoleculeSet == null) { // We found no matches to other ions, so create a new set for this ion.
						HashSet<Structure> moleculeSet = new HashSet<>();
						moleculeSet.add(molecule);
						knownMolecules.add(moleculeSet);
					}
				}

				// Write the sets to subdirectories and create a mean structure in each subdirectory 
				int setNum = 0;
				for (HashSet<Structure> moleculeSet : knownMolecules) {
					if (moleculeSet.size()== 0) {continue;}
					String dirName = compositionDir.getAbsolutePath() + "/set_" + (setNum++) + "/";
					new File(dirName).mkdirs();
					int moleculeNum = 0;
					for (Structure molecule : moleculeSet) {
						POSCAR outfile = new POSCAR(molecule, true);
						outfile.writeSpeciesWithDescription(false);
						outfile.writeElementSymbolsOnly(false);
						String outfileName = dirName + (moleculeNum++) + ".vasp";
						outfile.writeFile(outfileName);
					}
					
					// Create a mean structure by randomly sampling numSamples structures from this set.
					Structure meanStructure = IonTools.makeRepresentativeStructure(dirName, numSamples);
					String description = "Mean structure";
					String meanFileName = dirName + "mean.vasp";
					IonTools.writeIonStructure(meanStructure, meanFileName, description);
				}
			}
		}
	}

	/**
	 * Place the mean structures found by {@link groupIons(DataType, int, int)} in groups by total oxidation state, calculated
	 * by adding the oxidation states of all atoms in the ion.  The atomic oxidation states do not need to match
	 * each other for ions to be placed in the same group.  For each group, selects a representative structure that
	 * is structurally similar to all other ions in the group.
	 * 
	 * @param minAllowedFraction A representative structure will only be generated if the ratio of the number of structurally similar ions 
	 * (ignoring total oxidation states) to the total number of ions with the same composition is at least this.
	 * @param removeZeroes True if any ions with zero oxidation state should be removed
	 * @param minAllowedOccurrences A representative structure will only be generated if the ratio of the number of structurally similar ions 
	 * (ignoring total oxidation states) to the total number of ions with the same composition is greater than this.
	 */
	public static void groupIonsByOxidationState(double minAllowedFraction, boolean removeZeroes, int minAllowedOccurrences) {

		// The data set for which we are grouping the ions
		DataType dataType = DataType.ChargeBalanced;
		String ionStructureDirName = getIonStructureDirName(dataType);
		
		// The output directory
		String outDirName = getPolyatomicIonDirName(dataType);
		File ionTypeDir = new File(ionStructureDirName);
		
		/**
		 * Group mean ions by structural similarity, ignoring atomic oxidation states, and count
		 * the number we have in each group
		 */
		for (File periodicityDir : ionTypeDir.listFiles()) {
			for (File compositionDir : periodicityDir.listFiles()) {
					
				ArrayList<HashMap<Structure, Integer>> knownMolecules= new ArrayList<>();
				for (File setDir : compositionDir.listFiles()) {
					
					// Read the mean ion structure
					File meanFile = new File(setDir + "/mean.vasp");
					if (!meanFile.exists()) {continue;}
					Structure ion = IonTools.loadIonStructureFromFile(setDir + "/mean.vasp");
						
					// Count the number of ions in this group
					int count = setDir.listFiles(new VaspFileFilter()).length - 1;
					ion.setDescription("" + count);
					
					// This is the master set that will keep track of all matching mean structures
					HashMap<Structure, Integer> matchedMoleculeMap = null;
					
					// Check to see if any previously-found sets match the mean structure of this one
					for (HashMap<Structure, Integer> moleculeMap : knownMolecules) {
						for (Structure knownIon : moleculeMap.keySet()) {
							if (IonTools.compareStructures(knownIon, ion, true)) {
								moleculeMap.put(ion, count);
								if (matchedMoleculeMap == null) { // This is the master set of matching ions
									matchedMoleculeMap = moleculeMap;
								} else { // Ad the previously-discovered set of matching ions to the master set for this ion.
									matchedMoleculeMap.putAll(moleculeMap);
									moleculeMap.clear();
								}
								break;
							}
						}
					}
					
					// We didn't find any existing sets of structures to merge, so create a new set.
					if (matchedMoleculeMap == null) {
						HashMap<Structure, Integer> moleculeMap = new HashMap<>();
						moleculeMap.put(ion, count);
						knownMolecules.add(moleculeMap);
					}
				}
				
				// Trim the rare ions and those with zero oxidation states
				int numKeepers = 0;
				for (HashMap<Structure, Integer> moleculeMap : knownMolecules) {
					
					// The total count is the number of found ions with the same structure, ignoring atomic oxidation states
					int totalCount = 0;
					for (Integer count : moleculeMap.values()) {
						totalCount += count;
					}
					
					// Calculate the minimum number of found ions a particular structure must have to be considered.
					int minAllowedCount = (int) Math.ceil(totalCount * minAllowedFraction);
					
					// Keep track of which sets of molecules we should remove.
					HashSet<Structure> removedMolecules = new HashSet<>();
					for (Structure molecule : moleculeMap.keySet()) {
						
						// First see if we have any ions with zero oxidation state
						boolean hasZero = (molecule.getOxidationPerUnitCell() == 0);
						for (int siteNum = 0; siteNum < molecule.numDefiningSites(); siteNum++) {
							hasZero |= (molecule.getSiteSpecies(siteNum).getOxidationState() == 0);
						}
						if (removeZeroes && hasZero) {
							removedMolecules.add(molecule);
						// Next check if the number of found ions is sufficiently large
						} else if (moleculeMap.get(molecule) < minAllowedCount) {
							removedMolecules.add(molecule);
						}
					}
					
					/**
					 *  Remove everything we flagged above.  We do it in two steps because we don't
					 *  want to be removing things from the moleculeMap while we're iterating over it.
					 */
					for (Structure molecule : removedMolecules) {
						moleculeMap.remove(molecule);
					}
					numKeepers += moleculeMap.size();
				}
				
				// There's nothing here, so on to the next composition
				if (numKeepers == 0) {continue;}

				// Create the output directory
				String outCompositionDirName = outDirName + "/" + periodicityDir.getName() + "/" + compositionDir.getName() + "/";
				File outCompositionDir = new File(outCompositionDirName);
				
				// Write out all groups that are large enough.
				for (HashMap<Structure, Integer> moleculeMap : knownMolecules) {
					
					// Calculate the size of the group.
					int totalCount = 0;
					for (Integer count : moleculeMap.values()) {
						totalCount += count;
					}

					// If the group is not large enough, skip it
					if (totalCount < minAllowedOccurrences) {continue;}
					
					/**
					 *  Write a representative file for each group, organized by oxidation state
					 */
					outCompositionDir.mkdirs();
					int typeNum  = outCompositionDir.listFiles().length;
					String setDirName = outCompositionDirName + "/type_" + typeNum + "/";
					for (Structure molecule : moleculeMap.keySet()) {
						double oxidationState = molecule.getOxidationPerUnitCell();
						String oxidationStateDirName = setDirName + "/" + oxidationState + "/";
						File oxidationStateDir = new File(oxidationStateDirName);
						oxidationStateDir.mkdirs();
						int fileNum = oxidationStateDir.listFiles().length;
						String outFileName = oxidationStateDirName + "/" + fileNum + ".vasp";
						IonTools.writeIonStructure(molecule, outFileName, molecule.getDescription());
					}
				}
			}
		}
	}

	/**
	 * Starting with a directory of CIF files extracted from the ICSD, select all charge calanced
	 * structures and use them to construct the initial training data file
	 */
	public static void getInitialTrainingData() {

		DataType inDataType = DataType.ChargeBalanced;
		DataType outDataType = DataType.ChargeBalanced_KnownIons_Integer_NonZero_EAH;
		OxidationStateData data = getData(inDataType);
		data.removeEntriesWithZeroOxidationStates();
		data.removeEntriesWithNonIntegerStates();
		data.setEnergiesAboveHull(getEnergiesAboveHull());
		data.writeFile(getDataFileName(outDataType));

	}

	/**
	 * Reads energies above the convex hull from a file extracted from the Materials Project
	 * and creates a Map (dictionary, for you Python folks) keyed by the icsd_id and with the 
	 * energy above the hull as the value
	 * 
	 * @return A Map keyed by the icsd_id and with the energy above the hull as the value.
	 */
	public static Map<String, Double> getEnergiesAboveHull() {

		// We'll read the energies from this file
		String fileName = TRAINING_DATA_DIR + "MP_EnergyAboveHull_ICSD.txt";
		HashMap<String, Double> returnMap = new HashMap<>();

		try {
			LineNumberReader reader = new LineNumberReader(new FileReader(fileName));
			String line = reader.readLine();
			while (line != null && line.trim().length() > 0) {
				String[] fields = line.split(",");
				
				// Get the icsd id from the file name
				String id = fields[0].replace(".vasp", "").split(":")[1].trim();
				
				// Get the energy above hull
				double energy = Double.parseDouble(fields[2]);
				returnMap.put(id, energy);
				line = reader.readLine();
			}
			reader.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		return returnMap;
	}

	/**
	 * The "web ions" are the polyatomic ions used in the manuscript and web site.  This method 
	 * identifies all entries in the "inDataType" data set that contain these ions, creates a new version 
	 * of that entry that contains the composition written in terms of polyatomic ions, and adds the
	 * new entry to the data set.  The new data set is written to the file corresponding to "outDataType".
	 * 
	 * TODO this whole method should probably be added to the OxidationStateData class.
	 * 
	 * @param inDataType The data set to which entries with polyatomic ions should be added
	 * @param outDataType The name of the new data set 
	 */
	public static void addWebIonsToData(DataType inDataType, DataType outDataType) {

		// Read the polyatomic ions into memory
		loadWebIons();

		// This is the new file we will write
		String outDataFileName = getDataFileName(outDataType);

		// This is the data set to which we will add entries
		OxidationStateData data = getData(inDataType);

		/**
		 * Write the new data file line-by-line.  TODO it would be better to create new entries 
		 * and then just write the file using the appropriate method on OxidationStateData.
		 */	
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outDataFileName));

			int outEntryNum = 0;
			for (int entryNum = 0; entryNum < data.numEntries(); entryNum++) {
				Entry entry = data.getEntry(entryNum);
				String id = entry.getID();

				String outString = outEntryNum + ": ";
				outString += id + "\t";
				outString += entry.getGivenCompositionString() + "\t";
				for (Ion ion : entry.getAllIons()) {
					outString += " " + ion.getSymbol();
				}
				outString += "\tICSD";
				outString += "\t" + entry.getEnergyAboveHull();
				//outString += "\tfalse";
				outString += "\t" + entry.getGlobalInstabilityIndex();
				Status.basic(outString);
				writer.write(outString);
				writer.newLine();
				
				// Now we look for polyatomic ions and add a new entry if we find any.

				/**
				 * This reads the corresponding structure from a POSCAR-formatted file 
				 * and identifies all polyatomic ions in the structure.
				 */
				String structFileName = STRUCT_DIR + id + ".vasp";
				Structure structure = new Structure(new POSCAR(structFileName));
				KnownIonFinder finder = new KnownIonFinder(structure, true);

				if (finder.numFoundPolyatomicIons() > 0) { // Figure out the new composition

					outEntryNum++;
					outString = outEntryNum + ": ";
					outString += id + "\t";

					// This gets the composition from the finder
					Composition ionTypeComposition = finder.getIonTypeComposition();
					ionTypeComposition = ionTypeComposition.getReducedComposition();
					outString += ionTypeComposition.getStandardizedCompositionString();

					outString += "\t";

					// We also write out the individual ions with oxidaiton states
					Composition ionComposition = finder.getIonComposition();
					for (String ionSymbol : ionComposition.getSymbols()) {
						double count = ionComposition.getCount(ionSymbol);
						if (count < 1E-7) {continue;}
						outString += " " + ionSymbol;
					}
					
					// Finish writing the line for the new entry.
					outString += "\tICSD";
					outString += "\t" + entry.getEnergyAboveHull();
					outString += "\t" + entry.getGlobalInstabilityIndex();
					Status.basic(outString);
					writer.write(outString);
					writer.newLine();
				}
				outEntryNum++;
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		OxidationStateData outData = getData(outDataType);
		outData.removeEntriesWithZeroOxidationStates();
		outData.writeFile(outDataFileName);
	}

	/**
	 * Removes all entries with rare ions from a data set, where an ion is considered "rare" if it 
	 * appears in fewer than 25 entries.  The ion removal process is run repeatedly until no more 
	 * ions are removed.
	 * 
	 * @param inDataType The data set from which ions should be removed
	 * @param outDataType The new data set to be written.
	 */
	public static void removeRareIons(DataType inDataType, DataType outDataType) {

		loadWebIons();

		OxidationStateData data = getData(inDataType);
		
		// Iteratively repeat the removal process until no more entries are removed.
		int numEntries = data.numEntries();
		do {
			numEntries = data.numEntries();
			data.removeUncommonIonsByCount(25);
		} while (data.numEntries() < numEntries);

		data.writeFile(getDataFileName(outDataType));
	}

	/**
	 * Assigns oxidation states calculated using 10-fold cross validation to all of the structures
	 * in the given dataset and writes out corresponding CIF files to a directory.  The description
	 * of each CIF file will be the GII for that assignment, and oxidation states will be assigned 
	 * to sites in a way that minimizes the GII.
	 * 
	 * @param dataType The data set to be used.  Fitted parameters will be read for this data set when calculated the Likelihood Score.
	 * @param frequency True if the Frequency Score should be used to assign oxidation states, false if the Likelihood score should be used.
	 */
	public static void getValidationAssignments(DataType dataType, boolean frequency) {

		// Load the polyatomic ions into memory
		loadWebIons();
				
		Status.flow("Generating validation assignments for dataType = " + dataType + ", frequency = " + frequency);

		// Where everythign will be written
		String outBaseDirName = GII_DIR + dataType + "/";
		String outDirName = outBaseDirName + (frequency ? "Frequency_validation" : "Likelihood_validation");
		String outFileName = outDirName + ".txt";
		outDirName += "/";

		new File(outDirName).mkdirs();

		String paramDirName = PARAMETER_DIR + "/splits/" + dataType + "/5E-6/";
		String dataDirName = TRAINING_DATA_DIR + dataType + "/";

		// We will write out a new data set with GII values and oxidation states assigned
		OxidationStateData newData = new OxidationStateData(new ArrayList<Entry>(), STRUCT_DIR);

		/**
		 * Cycle through all of the test-training splits, and get the results of the model trained
		 * on the training split but evaluated on the test split
		 */
		for (int splitNum = 0; splitNum < 10; splitNum++) {
			String paramFileName = paramDirName + "/split_" + splitNum + "/parameters.txt";
			String validationDataName = dataDirName + "/test_" + splitNum + ".txt";
			String trainingDataName = dataDirName + "/training_" + splitNum + ".txt";
			OxidationStateCalculator calculator = frequency ? new FrequencyCalculator(new OxidationStateData(trainingDataName, STRUCT_DIR)) : new LikelihoodCalculator(paramFileName);
			OxidationStateData validationData = new OxidationStateData(validationDataName, STRUCT_DIR);
			
			for (int entryNum = 0; entryNum < validationData.numEntries(); entryNum++) {
				Entry entry = validationData.getEntry(entryNum);
				Structure structure  = entry.getStructure();
				IonAssigner assigner = null;
				
				/**
				 * Check to see if this entry is the polyatomic representation of the ICSD entry.  If it is,
				 * then we will assign oxidation states to atoms and calculate the GII based on the polyatomic
				 * oxidation states.  Otherwise just use elemental oxidation states.
				 */
				boolean noPoly = !entry.getGivenCompositionString().contains("(");
				OxidationStateSet stateSet = null;
				if (noPoly) {
					stateSet = calculator.getLikelyOxidationStates(structure);
					assigner = new IonAssigner(structure, stateSet);
				} else {
					KnownIonFinder finder = new KnownIonFinder(structure, true) ;//polyatomicIons);
					stateSet = calculator.getLikelyOxidationStates(finder.getIonTypeComposition());
					assigner = new IonAssigner(finder, stateSet);
				}
				
				/**
				 * Don't calculate the GII if some atoms have zero bond valence sums.  Sometimes this
				 * indicates something is wrong with the structure.
				 */
				boolean zeroSums = assigner.getBondValenceCalculator().hasZeroSums();
				double gii = zeroSums ? Double.NaN : assigner.getGlobalInstabilityIndex();
				
				// Save and write out the new entry
				newData.addEntry(entry.getID(), entry.getGivenCompositionString(), stateSet.getIons(), new String[] {dataType.toString()}, entry.getEnergyAboveHull(), gii);
				String outfileName = outDirName + entry.getID() + (noPoly ? "" : "_poly") + ".cif";
				writeAssignedFile(assigner.getStructure(), entry, outfileName, gii, zeroSums);
			}
		}
		newData.writeFile(outFileName);
	}

	/**
	 * Assigns oxidation states to all of the structures in the given dataset using the model trained
	 * on that dataset and writes out corresponding CIF files to a directory.  The description
	 * of each CIF file will be the GII for that assignment, and oxidation states will be assigned 
	 * to sites in a way that minimizes the GII.
	 * 
	 * @param dataType The data set to be used.  Fitted parameters will be read for this data set when calculated the Likelihood Score.
	 * @param frequency True if the Frequency Score should be used to assign oxidation states, false if the Likelihood score should be used.
	 */
	public static void getTrainingAssigments(DataType dataType, boolean frequency) {

		// Load the polyatomic ions into memory
		loadWebIons();
		Status.flow("Generating training assignments for dataType = " + dataType + ", frequency = " + frequency);

		// Set up the input and output files
		String outBaseDirName = GII_DIR + dataType + "/";
		String outDirName = outBaseDirName + (frequency ? "Frequency" : "Likelihood");
		String outFileName = outDirName + ".txt";
		outDirName += "/";
		new File(outDirName).mkdirs();

		// We will write out a new data set with GII values and oxidation states assigned
		OxidationStateData newData = new OxidationStateData(new ArrayList<Entry>(), STRUCT_DIR);

		// Load the data set and oxidation state calculator
		OxidationStateData trainingData = getData(dataType);
		OxidationStateCalculator calculator = frequency ? new FrequencyCalculator(trainingData) : getLikelihoodCalculator(dataType);
		
		/**
		 * Cycle through all entries in the data set and calculate oxidation states
		 * and GII for each of them.
		 */
		for (int entryNum = 0; entryNum < trainingData.numEntries(); entryNum++) {
			Entry entry = trainingData.getEntry(entryNum);
			Structure structure  = entry.getStructure();
			IonAssigner assigner = null;
			
			/**
			 * If the composition is expressed in terms of polyatomic ions, then assign oxidation
			 * states and calculate the GII using the polyatomic ions.
			 */
			boolean noPoly = !entry.getGivenCompositionString().contains("(");
			OxidationStateSet stateSet = null;
			if (noPoly) {
				stateSet = calculator.getLikelyOxidationStates(structure);
				assigner = new IonAssigner(structure, stateSet);
			} else {
				KnownIonFinder finder = new KnownIonFinder(structure, true);//polyatomicIons);
				stateSet = calculator.getLikelyOxidationStates(finder.getIonTypeComposition());
				assigner = new IonAssigner(finder, stateSet);
			}
			
			/**
			 * We don't assign GII if any of the atoms have bond valence sums of zero, as 
			 * this sometimes indicates something is wrong with the structure.
			 */
			boolean zeroSums = assigner.getBondValenceCalculator().hasZeroSums();
			double gii = zeroSums ? Double.NaN : assigner.getGlobalInstabilityIndex();
			
			// Save and write out the new entries with oxidation states and GII assigned.
			newData.addEntry(entry.getID(), entry.getGivenCompositionString(), stateSet.getIons(), new String[] {dataType.toString()}, entry.getEnergyAboveHull(), gii);
			String outfileName = outDirName + entry.getID() + (noPoly ? "" : "_poly") + ".cif";
			writeAssignedFile(assigner.getStructure(), entry, outfileName, gii, zeroSums);
		}
		newData.writeFile(outFileName);
	}

	/**
	 * Assigns oxidation states to all of the structures in the given dataset using oxidation states
	 * provided in the ICSD. The description of each CIF file will be the GII for that assignment.
	 * 
	 * @param dataType The data set to be used.  Fitted parameters will be read for this data set when calculated the Likelihood Score.
	 * @param reassign True if the oxidation states should be reassigned to sites in a way the minimizes the GII,
	 * which is usually (but not always) how they are assigned in the ICSD. False if we should just use the ICSD assignments.
	 */
	public static void getICSDAssignments(DataType dataType, boolean reassign) {

		// Load the polyatomic ions
		loadWebIons();

		Status.flow("Generating icsd assignments for dataType = " + dataType + ", reassign = " + reassign);

		// Set up the input and output files
		String outBaseDirName = GII_DIR + dataType + "/";
		String outDirName = outBaseDirName + "ICSD" + (reassign ? ".reassigned" : "");
		String outFileName = outDirName + ".txt";
		outDirName += "/";
		new File(outDirName).mkdirs();

		// Read the training data from the ICSD for this data set
		OxidationStateData icsdData = getData(dataType);

		// This file will contain the training data with oxidation states reassigned and GII calculated
		OxidationStateData newData = new OxidationStateData(new ArrayList<Entry>(), STRUCT_DIR);

		// Loop through the entries and reassign oxidation states and calculate GII
		for (int entryNum = 0; entryNum < icsdData.numEntries(); entryNum++) {
			Entry entry = icsdData.getEntry(entryNum);
			Structure structure  = entry.getStructure();
			
			/**
			 * If this entry has a composition written in terms of polyatomic ions, find the 
			 * oxidation states using polyatomic ions.  Otherwise use monatomic ions.
			 */
			boolean noPoly = !entry.getGivenCompositionString().contains("(");

			KnownIonFinder finder = new KnownIonFinder(structure, !noPoly); 
			Composition ionComposition = finder.getIonComposition();
			OxidationStateSet stateSet = new OxidationStateSet(ionComposition);

			IonAssigner assigner = null;
			double gii = Double.NaN;
			boolean zeroSums = false;
			
			/**
			 * IF we need to reassign the oxidation states, then reassign them here and use the reassigned
			 * states to calculate the GII.
			 */
			if (reassign) {
				assigner = noPoly ? new IonAssigner(structure, stateSet) : new IonAssigner(finder, stateSet);
				zeroSums = assigner.getBondValenceCalculator().hasZeroSums();
				gii = zeroSums ? Double.NaN : assigner.getGlobalInstabilityIndex();
			} else { // Otherwise just calculate the GII for the ICSD states.
				BondValenceCalculator bvCalculator = new BondValenceCalculator(structure);
				zeroSums = bvCalculator.hasZeroSums();
				gii = zeroSums ? Double.NaN : bvCalculator.getGlobalInstabilityIndex();
			}

			// Save and write out the new entry with GII and (if requestsed) reassigned states.
			newData.addEntry(entry.getID(), entry.getGivenCompositionString(), stateSet.getIons(), new String[] {dataType.toString()}, entry.getEnergyAboveHull(), gii);

			String cifFileName = outDirName + entry.getID() + (noPoly ? "" : "_poly") + ".cif";
			PartiallyOccupiedStructure outStructure = reassign ? assigner.getStructure() : new PartiallyOccupiedStructure(structure);
			writeAssignedFile(outStructure, entry, cifFileName, gii, zeroSums);
		}

		newData.writeFile(outFileName);
	}

	/**
	 * Reads the states calculated by BERTOS based on composition and assigns them to sites in a way
	 * that minimizes the GII.  The assigned states are used to calculate the GII
	 */
	public static void assignBertosStates() {

		// Input and output files
		String inputFileName = GII_DIR + "/BERTOS_CN_states.csv";
		String outDirName = GII_DIR + "/BERTOS_CN/";
		String outFileName = GII_DIR + "/BERTOS_CN.txt";

		// We'll write out a new data file with the BERTOS-assigned states and corresponding GII
		OxidationStateData newData = new OxidationStateData(new ArrayList<Entry>(), STRUCT_DIR);
		try {
			
			// Read a file containing the BERTOS assignments
			LineNumberReader reader = new LineNumberReader(new FileReader(inputFileName));
			String line = reader.readLine();
			while (line != null && line.trim().length() > 0) {
				
				// Parse the BERTOS assignment
				String[] tokens = line.split(",");
				String id = tokens[0];
				String composition = tokens[1];
				String bertosAssignment = tokens[2];
				String structFileName = STRUCT_DIR + id + ".vasp";
				OxidationStateSet stateSet = stateSetFromBertosString(bertosAssignment);
				if (stateSet == null) {
					Status.warning("Failed to parse assignment for structure " + id + ": " + bertosAssignment);
					line = reader.readLine();
					continue;
				}
				
				// Assign the BERTOS-calculated states to sites
				Structure structure = new Structure(new POSCAR(structFileName));
				IonAssigner assigner = new IonAssigner(structure, stateSet);
				
				/**
				 *  If some sites have a bond valence sum of zero, then don't calculate the GII
				 *  as this sometimes indicates a problem with the structure.
				 */
				boolean zeroSums = assigner.getBondValenceCalculator().hasZeroSums();
				double gii = zeroSums ? Double.NaN : assigner.getGlobalInstabilityIndex();
				
				// Save and write out the new entries.
				Entry entry = newData.addEntry(id, composition, stateSet.getIons(), new String[] {"BERTOS"}, Double.NaN, gii);
				String cifFileName = outDirName + id + ".cif";
				writeAssignedFile(assigner.getStructure(), entry, cifFileName, gii, zeroSums);
				line = reader.readLine();
			}
			newData.writeFile(outFileName);
			reader.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Reads an oxidation state prediction as output by BERTOS and converts it into an 
	 * OxidationStateSet for use with this code.
	 * 
	 * @param bertosString A set of oxidation states for a given composition as calculated
	 * using BERTOS
	 * @return An OxidationStateSet representing the assignment in the bertosString.
	 */
	public static OxidationStateSet stateSetFromBertosString(String bertosString) {
		HashMap<Ion, Integer> counts = new HashMap<>();

		String[] tokens = bertosString.trim().split(" ");
		for (String token : tokens) {
			
			// Get the number of atoms for this ion
			String compositionString = token.substring(0, token.indexOf("("));
			Composition composition = new Composition(compositionString);
			try {
				if (composition.getElementalCompositionMap().size() == 0) {return null;} // Probably an unknown element, like D
			} catch (CompositionParseException e) {
				return null;
			}
			String elementSymbol = composition.getSymbols().iterator().next();
			int count = (int) Math.round(composition.getCount(elementSymbol));

			// Read the corresponding oxidation state
			String stateString = token.substring(token.indexOf("(") + 1, token.indexOf(":"));
			Integer oxidationState = Integer.parseInt(stateString);

			// Increment the total number of ions of this type
			Ion ion = IonFactory.get(Element.getElement(elementSymbol), oxidationState);
			if (ion == null) {return null;} // Probably an unknown element, like D
			Integer knownCount = counts.get(ion);
			if (knownCount == null) {
				knownCount = 0;
			}
			knownCount += count;
			counts.put(ion, knownCount);
		}

		// Generate the OxidationStateSet to return
		double[] weights = new double[counts.size()];
		Ion[] ions = counts.keySet().toArray(new Ion[0]);
		for (int ionNum = 0; ionNum < weights.length; ionNum++) {
			weights[ionNum] = counts.get(ions[ionNum]);
		}
		return new OxidationStateSet(ions, weights);

	}

	/**
	 * Reads the states calculated by PyMatGen and assigns them to sites in a way
	 * that minimizes the GII.  The assigned states are used to calculate the GII
	 */
	public static void calcGIIForPymatgenStructs() {

		// These files are output by simple Python scripts and contain the PyMatGen-calculated states 

		// The PyMatGen BVAnalyzer tended to time out, so we increated the default timeout to see how much it helped.
		//String inDirName = GII_DIR + "pymatgen_bv_fromPOSCAR.noGII/";
		//String outDirName = GII_DIR + "pymatgen_bv_fromPOSCAR/";

		// This is the version used for the paper.
		String inDirName = GII_DIR + "pymatgen_bv_fromPOSCAR_1milPerm.noGII/";
		String outDirName = GII_DIR + "pymatgen_bv_fromPOSCAR_1milPerm/";

		// We'll write a new data file with PyMatGen-calculated states
		String outFileName = GII_DIR + "pymatgen_bv_fromPOSCAR.txt";
		OxidationStateData newData = new OxidationStateData(new ArrayList<Entry>(), STRUCT_DIR);

		/**
		 * Read through the structures with PyMatGen-assigned states
		 */
		File inDir = new File(inDirName);
		for (File infile : inDir.listFiles()) {
			if (!infile.getName().endsWith(".cif")) {continue;}
			
			// Read the structure
			CIF cif = new CIF(infile.getAbsolutePath());
			Structure structure = new Structure(cif);
			
			// Get the atomic and ionic compositions of the structure, sorted alphabetically
			Species[] species = structure.getDistinctSpecies();
			Ion[] allIons = new Ion[species.length];
			String[] symbols = new String[species.length];
			for (int specNum= 0; specNum < species.length; specNum++) {
				symbols[specNum] = IonFactory.get(species[specNum]).getSymbol();
			}
			Arrays.sort(symbols);
			for (int specNum= 0; specNum < species.length; specNum++) {
				allIons[specNum] = IonFactory.get(symbols[specNum]);
			}
			
			// Calculate the GII
			BondValenceCalculator bvCalculator = new BondValenceCalculator(structure);
			boolean disordered = !cif.isOrdered(0);
			
			// We only caluclate GII for ordered structures for now
			if (disordered) {
				Status.warning("Disordered structure: " + infile.getName());
				continue;
			}
			
			/** 
			 * We don't calculate the GII if any of the sites has a zero bond valence sum, as that
			 * sometimes indicates a problem with the structure.
			 */
			boolean zeroSums = bvCalculator.hasZeroSums();
			double gii = (zeroSums || disordered) ? Double.NaN : bvCalculator.getGlobalInstabilityIndex();
			
			// Save and output the results.
			String cifFileName = outDirName + infile.getName();
			Entry entry = newData.addEntry(infile.getName().replace(".cif",""), structure.getCompositionString(), allIons, new String[] {"PyMatGen"}, Double.NaN, gii);

			writeAssignedFile(new PartiallyOccupiedStructure(cif), entry, cifFileName, gii, zeroSums);
		}
		newData.writeFile(outFileName);

	}

	/**
	 * Looks through all entries of the data set given by inType and determines whether the GII is lower
	 * than the GII for the ICSD (calculated by assinging ICSD states to sites in a way that minimizes
	 * the GII).  A new data set, outType, is created with such entries removed.  Rare ions are also 
	 * removed.  If an 
	 * 
	 * To facilitate calculations and analysis, additional files are written listing the "remainder" entries
	 * (those that were removed) based on GII and rare ions.  Test-training splits for leave-10-out cross-validation
	 * are also generated.
	 * 
	 * @param inType The data set that we are cleaning
	 * @param outType The cleaned data set
	 */
	public static void cleanDataGII(DataType inType, DataType outType) {

		Status.includeTime(true);
		Status.setLogLevelDetail();
		
		// Load the polyatomic ions
		loadWebIons();

		// Input and output files
		String dataFileName = TRAINING_DATA_DIR + inType + ".txt";
		String outDataFileName = TRAINING_DATA_DIR + outType + ".txt";
		String giiRemainderDataFileName = TRAINING_DATA_DIR + outType + "_remainder_gii.txt";
		String ionRemainderDataFileName = TRAINING_DATA_DIR + outType + "_remainder_ions.txt";

		// This contains the reference GII -- typically the ICSD states, reassigned to minimize GII
		// TODO just use the data files instead of the directories of structure files
		String refDirectoryName = GII_DIR+ inType + "/ICSD.reassigned/";
		
		// If these GII are lower than the reference GII, the corresponding entry will be removed.
		String directoryName = GII_DIR + inType + "/Likelihood/";

		// Load the data to be cleaned
		Status.flow("Loading oxidation state data from " + dataFileName);
		OxidationStateData data = new OxidationStateData(dataFileName, STRUCT_DIR);

		// The calculator is used to track changes in likelihood while logging
		LikelihoodCalculator calculator = getLikelihoodCalculator(inType);

		// The GII remainder stores the entries that were removed because of the GII screen
		OxidationStateData gii_remainder = data.copy();
		
		// The ion remainder stores the entries that were removed because they contained rare ions
		OxidationStateData ion_remainder = data.copy();

		// Remove all entries where the GII in directoryName is lower than the GII in refDirectoryName
		// TODO Better to just load data files and use the GII field on entries.  (This was written before that field existed.)
		data.removeGIIDecrease(directoryName, refDirectoryName, calculator);

		// Calculate the remainder after structures were removed due to the global instability index
		gii_remainder.removeEntries(data);

		// An ion is considered rare if it appears in fewer than this many entries
		int minNumSamples = 25;
		
		// Iteratively repeat the removal process until no more entries are removed.
		int numEntries = data.numEntries();
		do {
			numEntries = data.numEntries();
			data.removeUncommonIonsByCount(minNumSamples);
		} while (data.numEntries() < numEntries);

		// Calculate the remainder after removing entries because they contained rare ions
		ion_remainder.removeEntries(data);
		ion_remainder.removeEntries(gii_remainder);

		// Write the output
		data.writeFile(outDataFileName);
		gii_remainder.writeFile(giiRemainderDataFileName);
		ion_remainder.writeFile(ionRemainderDataFileName);
		
		// Generate the leave-10-out cross-validation splits
		String splitDir = outDataFileName.replace(".txt", "/");
		splitData(data, 10, splitDir);
	}

	/**
	 * Randomly splits the data for leave-k-out cross validation in a way that ensures that no 
	 * composition appears in both the test and training set for any split.
	 * 
	 * @param data The data set to be split
	 * @param numSplits The number of splits to generate (What "k" is in leave-k-out cross validation).
	 * @param outDirName Where to write the splits being generated
	 */
	public static void splitData(OxidationStateData data, int numSplits, String outDirName) {

		// Load the polyatomic ions
		loadWebIons();
		
		// Create the output directory
		new File(outDirName).mkdirs();
		Status.flow("Writing file for unsplit data with " + data.numEntries() + " entries.");
		
		// Keep a record of the data set we split, before it was split
		data.writeFile(outDirName + "/unsplit_data.txt");

		// OxidationStateData.splitData does the heavy lifting for generating the splits.
		OxidationStateData[] splits = data.splitData(numSplits);
		
		// Now we just write data files for each of the splits.
		int splitNum = 0;
		for (OxidationStateData testSet : splits) {
			Status.flow("Writing files for split " + splitNum + " with " + testSet.numEntries() + " entries.");
			OxidationStateData trainingSet = data.copy();
			trainingSet.removeEntries(testSet);
			trainingSet.writeFile(outDirName + "/training_" + splitNum + ".txt");
			testSet.writeFile(outDirName + "/test_" + splitNum + ".txt");
			splitNum++;
		}
	}

	/**
	 * This method prints out summary statistics comparing sets of oxidation states assignments (e.g.
	 * from different methods for predicting oxidation states) for a given data set.
	 * 
	 * @param dataType The data set for which we are comparing assignments.
	 */
	public static void compareAssignments(DataType dataType) {

		// Load the polyatomic ions.
		loadWebIons();

		/**
		 * The energy range for which we are comparing assignments, where "energy" is the energy
		 * above the hull in eV as calculated by the Materials Project.
		 */
		double minEnergy = 0;
		double increment = 1000.025;
		double maxEnergy = minEnergy + increment;
		
		// Whether we should compare entries for which no energy was avaiable in the Materials Project.
		boolean allowMissingEnergies = false;
		
		// True if we should use the compositions written in terms of polyatomic ions; false otherwise
		boolean polyIon = false;
		
		// Whether we should skip entries that contain ions missing from the PyMatGen chemical space
		boolean skipPMGMissing = true;
		
		/**
		 * Whether we should only compare entries for which oxidation states could be successfully assigned
		 * by all methods in the paper
		 */
		boolean onlyCommonSetEntries = true;

		// Below are data files containing information about the assigned oxidation states.
		
		String dataFileName1 = GII_DIR + dataType + "/ICSD.reassigned.txt";
		//String dataFileName2 = GII_DIR + dataType + "/Frequency.txt";
		//String dataFileName2 = GII_DIR + dataType + "/Frequency_validation.txt";
		//String dataFileName2 = GII_DIR + dataType + "/Likelihood_validation.txt";
		//String dataFileName2 = GII_DIR + dataType + "/Likelihood.txt";
		//String dataFileName2 = GII_DIR + BERTOS_CN.txt";
		String dataFileName2 = GII_DIR  + "pymatgen_bv_fromPOSCAR_1milPerm.txt";
		//String dataFileName2 = GII_DIR + "pymatgen_bv_fromPOSCAR.txt";

		// This file contains information about which ions exist in our model but not in PyMatGen's chemical space
		String pymatgenMissingIonsFile = GII_DIR + "pymatgen_missing_ions.txt";

		// This is the set of entries that could be assigned oxidation states by BERTOS.  It will be used to define the common set.
		String bertosFileName = GII_DIR + "BERTOS_CN.txt";
		
		// This is the set of entries that could be assigned oxidation states by pymatgen.  It will be used to define the common set.
		String pmgFileName = GII_DIR + "pymatgen_bv_fromPOSCAR_1milPerm.txt";

		// These are ions that don't exist in pymatgen, so they need special treatment
		HashSet<Ion> pymatgenMissingIons = new HashSet<>();
		try {
			LineNumberReader reader = new LineNumberReader(new FileReader(pymatgenMissingIonsFile));
			String line = reader.readLine();
			while (line != null && line.trim().length() > 0) {
				Ion ion = IonFactory.get(line.trim());
				pymatgenMissingIons.add(ion);
				line = reader.readLine();
			}
			reader.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		// Load the two data sets to compare.
		OxidationStateData data1 = new OxidationStateData(dataFileName1, STRUCT_DIR);
		OxidationStateData data2 = new OxidationStateData(dataFileName2, STRUCT_DIR);

		// We use these data sets to create the common set
		OxidationStateData bertosData = new OxidationStateData(bertosFileName, STRUCT_DIR);
		OxidationStateData pmgData = new OxidationStateData(pmgFileName, STRUCT_DIR);

		/**
		 * Strip out any duplicate ICSD entries (some might be in the data set once with atomic 
		 * composition and once with polyatomic composition).
		 */
		HashMap<String, Entry> map1 = data1.getUniqueEntries(polyIon);
		HashMap<String, Entry> map2 = data2.getUniqueEntries(polyIon);
		HashMap<String, Entry> bertosMap = bertosData.getUniqueEntries(polyIon);
		HashMap<String, Entry> pmgMap = pmgData.getUniqueEntries(polyIon);

		Status.basic(map1.size() + " remaining entries in map 1.");
		Status.basic(map2.size() + " remaining entries in map 2.");
		Status.basic(bertosMap.size() + " remaining entries in bertos map.");
		Status.basic(pmgMap.size() + " remaining entries in pmg map.");

		// How many times the two data sets agree
		double numMatches = 0;
		
		// How many times an entry exists in data set 1 but not in data set 2.  These are generally failures to assign oxidation states.
		double numFailed = 0;
		
		// How many entries get screen out because they aren't in the pymatgen chemical space.
		double numPMGMissing = 0;
		
		// The number of entries with lower GII in data set 2 than in data set 1.
		double numGIILower = 0;
		
		// The total number of entries compared.
		double numTotal = 0;
		
		// The number of entries that are not in the common set
		double numNonCommonEntries = 0;
		
		// The number of entries that were not in the BERTOS set.
		double numNotInBertos = 0;
		
		// The number of entries that were not in the pymatgen set.
		double numNotInPMG = 0;
		
		// The number of entries for which the oxidation states in set 1 and set 2 did not match.
		double numNotMatched = 0;
		
		// Cycle through all of the entires in set 1.  This is usually the ICSD, so this is usually the more comprehensive set.
		for (String id : map1.keySet()) {
			
			Entry entry1 = map1.get(id);
			
			// Perform the enregy screen
			double eah = entry1.getEnergyAboveHull();
			if ((eah< minEnergy) || (eah > maxEnergy) || (!allowMissingEnergies && Double.isNaN(eah))) {continue;}

			// Perform the common set screen and count the number not in BERTOS and/or pymatgen.
			Entry bertosEntry = bertosMap.get(id);
			if (bertosEntry == null) {
				numNotInBertos++;
			}

			Entry pmgEntry = pmgMap.get(id);
			if (pmgEntry == null) {
				numNotInPMG++;
			}

			boolean nonCommonSet = false;

			if (bertosEntry == null || pmgEntry == null) {
				numNonCommonEntries++;
				nonCommonSet = true;
			}

			if (onlyCommonSetEntries && nonCommonSet) {
				continue;
			}

			// Peform the screen for ions in the pymatgen chemical space
			Ion[] ions1 = entry1.getAllIons();
			boolean pmgMissing = false;
			for (Ion ion : ions1) {
				if (pymatgenMissingIons.contains(ion)) {
					pmgMissing = true;
					Status.basic("PMG Missing: " + entry1.getGivenCompositionString());
					numPMGMissing++;
					break;
				}
			}
			if (skipPMGMissing && pmgMissing) {continue;}

			// This is the total number of entries that passed all screens and will be compared.
			numTotal++;

			// Get matching entry 2.  If it doesn't exist, mark this as a failed assignment and go to the next one.
			Entry entry2 = map2.get(id);
			if (entry2 == null) {
				numFailed++;
				continue;
			}
			
			// Compare the assignments
			Ion[] ions2 = entry2.getAllIons();
			if (Arrays.equals(ions1, ions2)) {numMatches++;} else {numNotMatched++;}
			
			// Compare the GII
			// TODO something like this would probably be a more efficient way to clean the data.
			if (entry2.getGlobalInstabilityIndex() < entry1.getGlobalInstabilityIndex()- 1E-6) {
				Status.basic(entry2.getID() + "\t" + entry2.getGivenCompositionString());
				numGIILower++;
			}
		}

		// Report the results.
		Status.basic("");
		Status.basic(numTotal + " entries in energy range.");
		Status.basic("Num not matched: " + (numTotal - numMatches - numFailed));
		Status.basic("Num failed: " + numFailed);
		Status.basic("Num PyMatGen missing: " + numPMGMissing);
		Status.basic("Num not in Bertos: " + numNotInBertos);
		Status.basic("Num not in PMG: " + numNotInPMG);
		Status.basic("Num not universal: " + numNonCommonEntries);
		Status.basic("");
		Status.basic("Match rate: " + numMatches / numTotal);
		Status.basic("Failure rate: " + numFailed / numTotal);
		Status.basic("PMGMissing rate: " + numPMGMissing / numTotal);
		Status.basic("Failure - PMGMissing rate: " + (numFailed - numPMGMissing) / numTotal);
		Status.basic("Matches / succeeded: " + numMatches / (numTotal - numFailed));
		Status.basic("No matches / succeeded: " + (numTotal - numFailed - numMatches) / (numTotal - numFailed));
		Status.basic("Matches + failed: " + (numMatches + numFailed) / numTotal);
		Status.basic("Fraction GII Lower: " + numGIILower / numTotal);

		Status.basic("");
		Status.basic("Total entries: " + numTotal);
		Status.basic("Fraction matched: " + numMatches / numTotal);
		Status.basic("Fraction failed: " + numFailed / numTotal);
		Status.basic("Fraction unmatched: " + numNotMatched / numTotal);
		Status.basic("Fraction sum: " + (numMatches + numFailed + numNotMatched) / numTotal);
	}

	/**
	 * Assign oxidation states to the data set used in the CAMD search for new materials.
	 * 
	 * @param dataType The data set that we will use to generate the assignments (we will use parameters fit to this data).
	 * @param noPoly True if we should write the CAMD compositions in terms of polyatomic ions, false otherwise.
	 */
	public static void getCAMDAssignments(DataType dataType, boolean noPoly) {

		Status.flow("Generating CAMD assignments for dataType = " + dataType);

		// The input and output files
		String outBaseDirName = GII_DIR + dataType + "/";
		String outDirName = outBaseDirName + "CAMD" + (noPoly ? "" : "_poly");
		String outFileName = outDirName + ".txt";
		outDirName += "/";
		new File(outDirName).mkdirs();

		// This is where we store the CAMD data
		String camdDirName = CAMD_DIR + "published/camd2022/files/";
		String camdFileName = camdDirName + "data_camd.csv";

		// We will write the oxidation state assignments to this file
		OxidationStateData newData = new OxidationStateData(new ArrayList<Entry>(), STRUCT_DIR);

		// Load an oxidation state calculator with parameters fit to the given data set.
		OxidationStateCalculator calculator = getLikelihoodCalculator(dataType);

		// Loop through the CAMD data and assign oxidation states where we can
		try {
			LineNumberReader reader = new LineNumberReader(new FileReader(camdFileName));
			for (String line = reader.readLine(); line != null && line.trim().length() != 0; line = reader.readLine()) {

				// Parse the CAMD data
				String[] fields = line.split(",");
				String compositionString = fields[0];
				if (new Composition(compositionString).getDistinctElements().length < 2) {continue;}
				double energy = Double.parseDouble(fields[1]);
				String fileName = fields[2].trim();

				// Now we assign the oxidation states
				IonAssigner assigner = null;
				OxidationStateSet stateSet = null;

				// Load the structure to which we will assign oxidation states
				String structFileName = camdDirName + fileName;
				String id = new File(structFileName).getName().replace(".cif", "");
				CIF infile = new CIF(structFileName);
				Structure structure = new Structure(infile);

				// Assign oxidation states to the structure
				if (noPoly) { // Use monatomic ions
					stateSet = calculator.getLikelyOxidationStates(structure);
					if (stateSet == null) {continue;}
					assigner = new IonAssigner(structure, stateSet);
				} else { // Use polyatomic ions
					KnownIonFinder finder = new KnownIonFinder(structure, true);
					stateSet = calculator.getLikelyOxidationStates(finder.getIonTypeComposition());
					if (stateSet == null) {continue;}
					assigner = new IonAssigner(finder, stateSet);
					
					// We get a new composition string in terms of the polyatomic ions
					compositionString = finder.getIonTypeComposition().getReducedComposition().getStandardizedCompositionString();
				}

				/**
				 * Don't calculate the global instability index if some bond valence sums are zero, as 
				 * this could indicate a problem with the structure
				 */
				boolean zeroSums = assigner.getBondValenceCalculator().hasZeroSums();
				double gii = zeroSums ? Double.NaN : assigner.getGlobalInstabilityIndex();

				// Save and write out new entry and assigned states
				Entry entry = newData.addEntry(id, compositionString, stateSet.getIons(), new String[] {dataType.toString()}, energy, gii);

				String cifFileName = outDirName + id + ".cif";
				PartiallyOccupiedStructure outStructure = assigner.getStructure();
				writeAssignedFile(outStructure, entry, cifFileName, gii, zeroSums);
			}
			reader.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		// Write out the new data file.
		newData.writeFile(outFileName);
	}

	/**
	 * This prints the data showing the percent of compositions with stable structures found vs the 
	 * percent of total compositions evaluated using the CAMD data sorted by Likelihood score.
	 * 
	 * @param dataType The data set for which analysis should be performed.  (The oxidation analyzer
	 * was trained using this data set.)
	 */
	public static void getCAMDDiscoveryCurve(DataType dataType) {

		// Load the polyatomic ions
		loadWebIons();

		// How many likelihood score bins should be in the resulting plot.
		int numBins = 100;
		
		// True if we should use the frequency score, false if we should use the likelihood score
		boolean useFrequency = false;
		
		// True if we should write compositions in terms of polyatomic ions, false otherwise
		boolean polyIon = false;
		
		// Only evaluate compounds with at least this many distinct elements.  Should be no fewer than 2.
		int minNumElements = 2;
		
		// Only evaluate compounds with no more than this many distinct elements.
		int maxNumElements = 12;

		// We read the CAMD data from here.
		String camdDirName = CAMD_DIR + "published/camd2022/files/";
		String camdFileName = camdDirName + "data_camd.csv";

		// Build the calculator
		OxidationStateCalculator evaluator = useFrequency ? new FrequencyCalculator(getData(dataType)) : getLikelihoodCalculator(dataType);

		// We store information for the lowest-energy structure at each composition
		HashMap<String, double[]> dataByComposition = new HashMap<>();

		try {
			
			// Read the CAMD data and assign oxidation states
			LineNumberReader reader = new LineNumberReader(new FileReader(camdFileName));
			String line = reader.readLine();
			int structNum;
			int numFailed = 0;
			for (structNum = 1; line != null && line.trim().length() != 0; structNum++) {

				// Parse the CAM data
				String[] fields = line.split(",");
				line = reader.readLine();
				String compositionString = fields[0];
				Composition composition = new Composition(compositionString).getReducedComposition();
				if ((composition.getDistinctElements().length < minNumElements) || (composition.getDistinctElements().length > maxNumElements)) {continue;}

				// Get the composition in terms of polyatomic ions, if polyIon is true.
				if (polyIon) {
					String structFileName = camdDirName + fields[2].trim();
					CIF infile = new CIF(structFileName);
					KnownIonFinder finder = new KnownIonFinder(new Structure(infile), true); //polyatomicIons);
					compositionString = finder.getIonTypeComposition().getReducedComposition().getStandardizedCompositionString();
				} else {
					compositionString = composition.getStandardizedCompositionString();
				}

				double energy = Double.parseDouble(fields[1]);

				double[] knownData = dataByComposition.get(compositionString);
				if (knownData != null) { // We found a lower-energy structure with the same composition, so just update the appropriate fields.
					double prevEnergy = knownData[0];
					if (energy < prevEnergy) {
						knownData[0] = energy;
						knownData[2] = structNum;
					}
					continue;
				}

				// New composition, so calculate the oxidation states
				OxidationStateSet stateSet = evaluator.getLikelyOxidationStates(compositionString);
				if (stateSet == null) {
					numFailed++;
					continue;
				}
				
				// Get the score for the oxidation states
				Ion[] likelyStates = stateSet.getIons();
				double score = useFrequency ? ((FrequencyCalculator) evaluator).getFrequencyScore(likelyStates) : ((LikelihoodCalculator) evaluator).optimizeLikelihood(likelyStates).getMaxLikelihood();

				// Store the information about the oxidation states at this composition
				double[] newData = new double[] {energy, score, structNum, composition.getDistinctElements().length};
				dataByComposition.put(compositionString, newData);
			}
			reader.close();

			// Sort the dataset by likelihood in descending order and count the total number of compositions on the convex hull
			double[][] dataArray = dataByComposition.values().toArray(new double[0][]);
			double[] negLikelihoods = new double[dataArray.length];
			int numOnHull = 0;
			for (int setNum = 0; setNum < dataArray.length; setNum++) {
				double[] dataSet = dataArray[setNum];
				negLikelihoods[setNum] = -dataSet[1];
				if (dataSet[0] == 0) {numOnHull++;}
			}
			int[] map = ArrayUtils.getSortPermutation(negLikelihoods);

			// Keep track of the number on the hull for each sorted bin
			double cumulativeNumOnHull = 0;
			int entryNum = 0;
			for (double binNum = 0; binNum < numBins; binNum++) {

				while(entryNum < map.length * (binNum + 1) / numBins) {
					if (dataArray[map[entryNum]][0] == 0) {cumulativeNumOnHull++;}
					entryNum++;
				}
				Status.basic("" + (cumulativeNumOnHull / numOnHull));
			}
			
			/// Report results
			Status.basic(numFailed + " failed.");
			Status.basic(structNum + " structures evaluated.");
			Status.basic(map.length + " total entries.");
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Get the histogram of what percentage of structures are on the hull as a function of likelihood score
	 * for the CAMD data.
	 * 
	 * @param dataType THe data set used to parameterize the oxidation predictions.
	 */
	public static void getCAMDHullHistogram(DataType dataType) {

		// Read the polyatomic ions.
		loadWebIons();

		// How many likelihood score bins will the histogram have.
		int numBins = 5;
		
		// True if the frequency score should be used, but false if the likelihood score should be used.
		boolean useFrequency = false;
		
		// The range of the number of distinct elements for structures to be analyzed.  minNumElements should be at least 2.
		int minNumElements = 2;
		int maxNumElements = 12;

		// Where to find the CAMD data
		String camdFileName = CAMD_DIR + "published/camd2022/files/data_camd.csv";

		// Build the calculator
		OxidationStateCalculator evaluator = useFrequency ? new FrequencyCalculator(getData(dataType)) : getLikelihoodCalculator(dataType); // new LikelihoodCalculator(paramFileName);

		// Store information for the lowest-energy structure at each composition
		HashMap<String, Double[]> dataByComposition = new HashMap<>();
		HashMap<String, Integer> countsByComposition = new HashMap<>();

		/**
		 * Loop through the CAMD data set and count how many structures are on the hull, as well as the 
		 * total number of structures, by bin.
		 */
		int numFailed = 0;
		int numSucceeded = 0;
		try {
			LineNumberReader reader = new LineNumberReader(new FileReader(camdFileName));
			String line = reader.readLine();
			for (int structNum = 0; line != null && line.trim().length() != 0; structNum++) {

				String[] fields = line.split(",");
				line = reader.readLine();

				// Get the composition
				String compositionString = fields[0];
				Composition composition = new Composition(compositionString).getReducedComposition();
				if ((composition.getDistinctElements().length < minNumElements) || (composition.getDistinctElements().length > maxNumElements)) {continue;}
				compositionString = composition.getStandardizedCompositionString();

				// Energy above the hull
				double energy = Double.parseDouble(fields[1]);

				// How many structures we have by composition
				Integer count = countsByComposition.get(compositionString);
				if (count == null) {count = 1;}
				countsByComposition.put(compositionString, count + 1);

				/**
				 * If we've seen something with this composition, update the relevant fields
				 * if the energy above the hull is lower.
				 */
				Double[] knownData = dataByComposition.get(compositionString);
				if (knownData != null) {
					double prevEnergy = knownData[0];
					if (energy < prevEnergy) {
						knownData[0] = energy;
						knownData[2] = structNum * 1.0;
					}
					continue;
				}

				// Evaluate the oxidation states
				OxidationStateSet stateSet = evaluator.getLikelyOxidationStates(fields[0]);
				if (stateSet == null) {
					numFailed++;
					continue;
				}
				numSucceeded++;
				
				// Get the score
				Ion[] likelyStates = stateSet.getIons();
				double score = useFrequency ? ((FrequencyCalculator) evaluator).getFrequencyScore(likelyStates) : ((LikelihoodCalculator) evaluator).optimizeLikelihood(likelyStates).getMaxLikelihood();

				// Add the data for this composition
				Double[] newData = new Double[] {energy, score, structNum * 1.0};
				dataByComposition.put(compositionString, newData);
			}
			reader.close();

			// Now create the histogram
			double[] onHullBins = new double[numBins];
			double[] offHullBins = new double[numBins];

			for (String composition : dataByComposition.keySet()) {
				Double[] values = dataByComposition.get(composition);
				double energy = values[0];
				double likelihood = values[1];

				// Figure out which bin we are in
				double[] bins = (energy == 0) ? onHullBins : offHullBins;
				int binNum = (int) Math.floor(likelihood * numBins);
				binNum = Math.min(binNum, numBins - 1);
				
				// Increment the total number of structures for this bin
				bins[binNum]++;
			}

			// Print the bin data
			Status.basic(ArrayUtils.toString(onHullBins));
			Status.basic(ArrayUtils.toString(offHullBins));

			// Get the fraction of structures on the hull by bin
			for (int binNum = 0; binNum < numBins; binNum++) {
				double total = onHullBins[binNum] + offHullBins[binNum];
				if (total == 0) {continue;}
				onHullBins[binNum] /= total;
				offHullBins[binNum] /= total;
			}

			// Report the results.
			Status.basic(ArrayUtils.toString(onHullBins));
			Status.basic(ArrayUtils.toString(offHullBins));

			Status.basic(numFailed + " compositions failed.");
			Status.basic(numSucceeded + " compositions succeeded.");

		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Write the JSON file containing the oxidation state boundaries.
	 * 
	 * @param dataType The data set for which we will write the boundaries.
	 */
	public static void writeBoundaryJSON(DataType dataType) {

		// Load the polyatomic ions.
		loadWebIons();

		// Write the boundaries
		String outFileName = "input_files/oxidation_boundaries.json";

		LikelihoodCalculator calculator = getLikelihoodCalculator(dataType);
		calculator.writeBoundaryJSON(outFileName);

	}

	/**
	 * Write files in xyz format for each of the polyatomic ions used for the web site nad paper
	 */
	public static void writeXYZFiles() {

		// Load the polyatomic ions
		loadWebIons();

		// This is used to create the file names
		DecimalFormat formatter = new DecimalFormat("+#;-#");
		for (IonType ionType : IonFactory.getKnownIonTypes().values()) {
			
			// Get the representative structure for the ion type
			Map<Integer, Structure> ionStructures = IonFactory.getKnownRepresentativeStructures(ionType.getSymbol()); 
			
			// Loop through all known oxidation states
			for (Integer oxidationState : ionStructures.keySet()) {
				Structure ionStructure = ionStructures.get(oxidationState);
				
				// Create the filename 
				String symbol = ionType.getSymbol(); 
				String name = symbol + formatter.format(oxidationState) + ".xyz";
				
				// Write the file
				XYZFile outfile = new XYZFile(ionStructure);
				outfile.writeFile(ROOT_DIR+ "../Paper/polyatomic_ion_xyz_files/" + name);
			}
		}
	}

	/**
	 * Loads  the polyatomic ions used for the web site and paper to the IonFactory.  This should usually be
	 * called first, as the code frequently checks IonFactory to get the list of known polyatomic ions.
	 * It's OK to call this multiple times on the same directory of ions.
	 */
	public static void loadWebIons() {

		String ionStructureDirName = getWebIonDirName();

		// Where the ions are located
		Status.flow("Reading ion information from " + ionStructureDirName);
		
		// Most of the work is done here.
		IonFactory.loadPolyatomicIons(ionStructureDirName);
	}
	
	/**
	 * Returns the directory name for the polyatomic ions used for the paper and web site
	 * 
	 * @return the directory name for the polyatomic ions used for the paper and web site
	 */
	public static String getWebIonDirName() {
		return  ROOT_DIR + "/ion_types_web/";		
	}

	/**
	 * Writes a CIF file containing the calculated GII in the description field and ssigned oxidation states.
	 * 
	 * @param structure The structure to be written
	 * @param entry Then corresponding entry in the data set.
	 * @param fileName The name of the file to be written.
	 * @param gii The calculated global instability index for the structure.
	 * @param zeroSums True if the structure contains some sites with bond valence sums of zero, false
	 * otherwise.  If true, the structure file will not be written.
	 */
	private static void writeAssignedFile(PartiallyOccupiedStructure structure, Entry entry, String fileName, double gii, boolean zeroSums) {

		// Don't write the file if some sites have bond valence sums of zero.
		if (zeroSums) {
			Status.warning("Found zero bond valence sums for " + entry.getID() + ", " + entry.getGivenCompositionString());
			return;
		}
		
		// Don't write the structure if the GII couldn't be calculated.
		if (Double.isNaN(gii)) {
			Status.warning("Couldn't calculate gii for " + entry.getID() + ", " + entry.getGivenCompositionString());
			return;
		}
		
		// Write the file
		CIF outfile = new CIF(structure);
		outfile.setDescription("" + gii);
		outfile.writeFile(fileName);
	}

	/**
	 * Write all of the distinct oxidation states in a given data set
	 * 
	 * @param dataType The data set for which we should print out the oxidation states
	 */
	public static void printAllKnownStates(DataType dataType) {

		// Load the polyatomic ions.
		loadWebIons();

		// A simple way to get all of the oxidation states for the data set.
		LikelihoodCalculator calculator = getLikelihoodCalculator(dataType);
		Ion[] allIons = calculator.getKnownOxidationStates();

		/**
		 *  Sort the states by ion type symbol (alphabetically), then by oxidation state.
		 *  The TreeMap and TreeSet do the sorting.
		 */
		TreeMap<String, TreeSet<Integer>> knownStates = new TreeMap<>();
		for (Ion ion : allIons) {
			int oxidationState = (int) ion.getOxidationState();
			String structSymbol = ion.getIonType().getSymbol(); //.getStructureSymbol();
			TreeSet<Integer> statesForIon = knownStates.get(structSymbol);
			if (statesForIon == null) {statesForIon = new TreeSet<>();}
			statesForIon.add(oxidationState);
			knownStates.put(structSymbol, statesForIon);
		}

		// Generate a similar sorted list for the ions that have been remvoed.
		TreeMap<String, TreeSet<Integer>> removedStates = new TreeMap<>();
		for (String ionSymbol : REMOVED_IONS) {
			Ion ion = IonFactory.get(ionSymbol);
			int oxidationState = (int) ion.getOxidationState();
			String structSymbol = ion.getIonType().getSymbol(); //.getStructureSymbol();
			TreeSet<Integer> statesForIon = removedStates.get(structSymbol);
			if (statesForIon == null) {statesForIon = new TreeSet<>();}
			statesForIon.add(oxidationState);
			removedStates.put(structSymbol, statesForIon);
		}

		// Print out the included states
		// Regex from https://stackoverflow.com/questions/5243316/format-a-number-with-leading-sign
		DecimalFormat formatter = new DecimalFormat("+#;-#");
		for (String structSymbol : knownStates.keySet()) {
			TreeSet<Integer> states = knownStates.get(structSymbol);
			System.out.print(structSymbol + "\t");
			for (int state: states) {
				System.out.print(formatter.format(state) + ", ");
			}

			states = removedStates.get(structSymbol);
			if (states != null) {
				System.out.print("\t");
				for (int state: states) {
					System.out.print(formatter.format(state) + ", ");
				}
			}
			System.out.println();
		}

		// Print out the removed states
		Status.basic("");
		Status.basic("All states removed:");
		for (String symbol : removedStates.keySet()) {
			if (!knownStates.containsKey(symbol)) {
				TreeSet<Integer> states = removedStates.get(symbol);
				System.out.print(symbol + "\t");
				for (int state: states) {
					System.out.print(formatter.format(state) + ", ");
				}
				System.out.println();
			}
		}

	}

	/**
	 * Generate a table for the electrochemical series, consisting of boundaries for all redox
	 * pairs in the data set.  In general it will not be sorted, so you may want to sort the output.
	 * 
	 * @param dataType The data set for which we will print out the electrochemical series.
	 */
	public static void printElectrochemicalSeries(DataType dataType) {
		
		// Load the polyatomic ions.
		loadWebIons();
		
		// A quick way to get all of the ions and boundaries
		LikelihoodCalculator calculator = getLikelihoodCalculator(dataType);

		// We print out everything in terms of the mapped potential (ICSD derived electrochemical potential)
		PotentialMapper potentialMapper = new PotentialMapper();

		// Print out all of the boundaries
		for (int ionTypeIndex = 0; ionTypeIndex < calculator.numIonTypes(); ionTypeIndex++) {
			double[] boundaries = calculator.getBoundaries(ionTypeIndex);
			String ionSymbol = calculator.getIonTypeSymbol(ionTypeIndex);
			int[] oxidationStates = calculator.getOxidationStates(ionSymbol);
			for (int boundNum = 0; boundNum < boundaries.length; boundNum++) {

				String leftState = (boundNum == 0) ? "N/A" : IonFactory.get(ionSymbol, oxidationStates[boundNum - 1]).toString();
				String rightState = (boundNum == boundaries.length - 1) ? "N/A" : IonFactory.get(ionSymbol, oxidationStates[boundNum]).toString();
				Status.basic(leftState + "\t" + rightState + "\t" + potentialMapper.toMappedPotential(boundaries[boundNum]));

			}
		}

	}

	/**
	 * Geneate the files we use to generate the wavy bar plots (or equivalent straight line plots).
	 * 
	 * @param dataType The data set for which we are generating the plots.
	 * @param smoothCutoff True if the plots should show the logistic function cutoff, and false if 
	 * they should just show the mean boundary values.
	 */
	public static void writeDataFilesForVisualization(DataType dataType, boolean smoothCutoff) {

		// Load the polyatomic ions.
		loadWebIons();

		// Read the parameters (boundareis) form this
		LikelihoodCalculator calculator = getLikelihoodCalculator(dataType);
		
		// Where we write everything
		String outDirName = ROOT_DIR + "/Visualization/" + dataType + "/DataFiles/";
		outDirName += smoothCutoff ? "SmoothCutoff/" : "SharpCutoff/";
		new File(outDirName).mkdirs();
		
		// Write the files.
		calculator.writeDataFiles(outDirName, smoothCutoff);

	}

	/**
	 * Return the data set for the given data type.
	 * 
	 * @param dataType The type of data set we want
	 * @return The data set containing entries with compositions and oxidation states.
	 */
	public static OxidationStateData getData(DataType dataType) {

		String dataFileName = getDataFileName(dataType);
		return new OxidationStateData(dataFileName, STRUCT_DIR);
	}

	/**
	 * Return the file name for the given data set
	 * 
	 * @param dataType The data set we are interested in
	 * @return The file name
	 */
	public static String getDataFileName(DataType dataType) {
		return TRAINING_DATA_DIR + dataType + ".txt";
	}

	/**
	 * Return the polyatomic ion structure directory name for the given data set
	 * 
	 * @param dataType The data set we are interested in
	 * @return The directory name
	 */
	public static String getIonStructureDirName(DataType dataType) {
		return ROOT_DIR + "/Ion_Structures/" + dataType + "/";
	}

	/**
	 * Return the name of directory of found polyatomic ions for the given data set
	 * 
	 * @param dataType The data set we are interested in
	 * @return The directory name
	 */
	public static String getPolyatomicIonDirName(DataType dataType) {
		return ROOT_DIR + "/Polyatomic_Ions/" + dataType + "/";
	}

	/**
	 * Return a likelihood calculator for the given data type.
	 * 
	 * @param dataType The data set we are interested in.
	 * @return A likelihood calculator using parameters trained on the given data set.
	 */
	public static LikelihoodCalculator getLikelihoodCalculator(DataType dataType) {

		String paramFileName = getParamFileName(dataType);
		return new LikelihoodCalculator(paramFileName);

	}
	
	/**
	 * Returns the name of the fitted parameter file for the given data type.  This method assumes
	 * a regularization parameter of 5E-6.
	 * 
	 * @param dataType The given data type
	 * @return the name of the fitted parameter file for the given data type
	 */
	public static String getParamFileName(DataType dataType) {

		return PARAMETER_DIR + dataType + "/5E-6/parameters.txt";
	}

	/**
	 * This method provides an example of how to call the API to replicate the table-generating functionality of
	 * the web app.
	 */
	public static void testWebAPI() {
		
		String paramFileName = "input_files/oxidation_boundaries.json";
		String polyIonDir = "input_files/polyatomic_ions_web";
		WebOxidationAnalyzer analyzer = new WebOxidationAnalyzer(paramFileName, polyIonDir);
		
		// Analyzing oxidation states from composition
		PageData pageData = analyzer.getPageDataFromComposition("LiMn2O4");
		
		// Alternative, for reading structure info from a file
		/*String fileString = null;
		Path path = Paths.get("PATH_HERE");
		try {
			fileString = Files.readString(path);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		PageData pageData = analyzer.getPageDataFromStructure(fileString);*/

		System.out.println(pageData);
	}
}
