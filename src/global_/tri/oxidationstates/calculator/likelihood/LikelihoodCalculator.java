package global_.tri.oxidationstates.calculator.likelihood;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.TreeSet;

import com.fasterxml.jackson.databind.ObjectMapper;

import global_.tri.oxidationstates.calculator.MixedValenceChargeBalanceFilter;
import global_.tri.oxidationstates.calculator.OxidationStateCalculator;
import global_.tri.oxidationstates.fitting.OxidationStateData.Entry;
import global_.tri.oxidationstates.ion.IonFactory;
import global_.tri.oxidationstates.ion.IonFactory.Ion;
import global_.tri.oxidationstates.ion.IonFactory.IonType;
import global_.tri.oxidationstates.util.Composition;
import global_.tri.oxidationstates.webapi.PotentialMapper;
import global_.tri.oxidationstates.webapi.RangeData;
import matsci.Element;
import matsci.Species;
import matsci.io.app.log.Status;
import matsci.util.arrays.ArrayIndexer;
import matsci.util.arrays.ArrayUtils;

/**
 * This is the main class for calculating likelihood scores
 * 
 * @author timmueller
 *
 */
public class LikelihoodCalculator extends OxidationStateCalculator {

	// Which oxidation states are allowed for each ion type, in order of the ion
	// tyep indices
	private int[][] m_OxidationStatesByIonType = new int[Element.numKnownElements() + 1][0];

	// These are the allowed ion types, in order of the ion type indices
	private String[] m_IonTypes;

	// The ion type indices by ion type id
	private HashMap<String, Integer> m_IonTypeIndices = new HashMap<String, Integer>();

	/**
	 * The main index here is the parameter index. Each entry contains a 2-element
	 * array, where the first is the ion type index and the second is the oxidation
	 * state index for that type.
	 */
	private int[][] m_StateIndexByParam;

	/**
	 * The parameter index for each oxidation state. The main index is the ion type
	 * index, and the seconary index is the oxidation state index for that ion type.
	 */
	private int[][] m_ParamIndexByState;

	// How small the change in electronic chemical potential should get before we
	// stop trying to optimzie it
	private double m_PotentialResolution = 1E-6;

	// The current set of parameters used to calculate oxidation states
	private double[] m_Parameters;

	/**
	 * Create a new likelihood calculator, cloning an old one and updating the
	 * paramters
	 * 
	 * @param source        The old calculator
	 * @param newParameters The new parameters to use. They should be in the same
	 *                      order as the parameters for the old calculator.
	 */
	private LikelihoodCalculator(LikelihoodCalculator source, double[] newParameters) {
		m_OxidationStatesByIonType = source.m_OxidationStatesByIonType;
		m_IonTypes = source.m_IonTypes;
		m_IonTypeIndices = source.m_IonTypeIndices;
		m_StateIndexByParam = source.m_StateIndexByParam;
		m_ParamIndexByState = source.m_ParamIndexByState;
		m_Parameters = newParameters.clone();
	}

	/**
	 * Initialize a likelihood calculator from a parameters file. This is an older
	 * constructor; to initialize from a JSON file, use the
	 * {@link #LikelihoodCalculator(String, boolean) LikelihoodCalculator(String,
	 * boolean)} constructor.
	 * 
	 * @param paramFileName The name of the parameter file to read
	 */
	public LikelihoodCalculator(String paramFileName) {
		this(paramFileName, false);
	}

	/**
	 * Initialize a likelihood calculator from a parameters file.
	 * 
	 * @param inputFileName      The name of the parameter file
	 * @param isJSONBoundaryFile True if the parameter file is a JSON-formated
	 *                           boundaries file, false if it's a parameters file
	 */
	public LikelihoodCalculator(String inputFileName, boolean isJSONBoundaryFile) {

		m_Parameters = isJSONBoundaryFile ? readParametersFromBoundaryJSON(inputFileName)
				: readParameters(inputFileName);

	}

	/**
	 * Initialized the likelihood calculator from a map of molecule types to allowed
	 * oxidation states
	 * 
	 * @param oxidationStateMap The map, in which keys are molecule type ids and the
	 *                          values are the allowed oxidation states
	 */
	public LikelihoodCalculator(HashMap<String, int[]> oxidationStateMap) {

		ArrayList<String> moleculeIDs = new ArrayList<String>(oxidationStateMap.keySet());
		AtomicMassComparator comparator = new AtomicMassComparator();
		moleculeIDs.sort(comparator);

		ArrayList<int[]> oxidationStates = new ArrayList<int[]>();
		for (String moleculeID : moleculeIDs) {
			oxidationStates.add(oxidationStateMap.get(moleculeID));
		}

		this.initializeOxidationStates(oxidationStates, moleculeIDs);
	}

	/**
	 * Construct the internal data structures for oxidation states
	 * 
	 * @param oxidationStates A list of allowed oxidation states, ordered in the
	 *                        same way as the list of ion types
	 * @param ionTypesA       list of allowed ion types, ordered in the same way as
	 *                        the list of oxidation states
	 */
	private void initializeOxidationStates(List<int[]> oxidationStates, List<String> ionTypes) {

		// Create the ion type index maps
		int numParams = 0;
		m_OxidationStatesByIonType = new int[oxidationStates.size()][];
		m_IonTypes = new String[oxidationStates.size()];
		for (int ionTypeIndex = 0; ionTypeIndex < ionTypes.size(); ionTypeIndex++) {
			String ionType = ionTypes.get(ionTypeIndex);
			m_IonTypeIndices.put(ionType, ionTypeIndex);
			m_IonTypes[ionTypeIndex] = ionType;
			m_OxidationStatesByIonType[ionTypeIndex] = oxidationStates.get(ionTypeIndex);
			Arrays.sort(m_OxidationStatesByIonType[ionTypeIndex]);
			int numStates = m_OxidationStatesByIonType[ionTypeIndex].length;
			numParams += (numStates == 0) ? 0 : numStates + 1;
		}

		m_Parameters = new double[numParams];
		m_StateIndexByParam = new int[numParams][];
		m_ParamIndexByState = new int[m_OxidationStatesByIonType.length][];

		int paramNum = 0;
		for (int ionTypeIndex = 0; ionTypeIndex < m_OxidationStatesByIonType.length; ionTypeIndex++) {
			int numStates = m_OxidationStatesByIonType[ionTypeIndex].length;
			int numParamsForIonType = (numStates == 0) ? 0 : numStates + 1;
			m_ParamIndexByState[ionTypeIndex] = new int[numParamsForIonType];
			for (int paramNumForElement = 0; paramNumForElement < numParamsForIonType; paramNumForElement++) {
				m_StateIndexByParam[paramNum] = new int[] { ionTypeIndex, paramNumForElement };
				m_ParamIndexByState[ionTypeIndex][paramNumForElement] = paramNum;
				paramNum++;
			}
		}

	}

	/**
	 * Whether this calculator has parameters for the given ion
	 * 
	 * @param ion The given ion
	 * @return True if this calculator has parameters for the ion, and false
	 *         otherwise.
	 */
	public boolean hasParamsForIon(Ion ion) {
		if (!m_IonTypeIndices.containsKey(ion.getIonType().getSymbol())) {
			return false;
		}
		int[] oxidationStates = this.getOxidationStates(ion.getIonType().getSymbol());
		for (int oxidationState : oxidationStates) {
			if (Math.abs(oxidationState - ion.getOxidationState()) < 1E-6) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Whether this calculator has parameters for the given ion type
	 * 
	 * @param symbol The symbol for the given ion type
	 * @return True if this calculator has parameters for the given ion type, and
	 *         false otherwise
	 */
	public boolean hasParamsForIonType(String symbol) {
		return m_IonTypeIndices.containsKey(symbol);
	}

	/**
	 * Whether this calculator has parameters for the given element
	 * 
	 * @param element The given element
	 * @return True if the calculator has parameters for the given element, false
	 *         otherwise
	 */
	public boolean hasParamsForElement(Element element) {
		return m_IonTypeIndices.containsKey(element.getSymbol());
	}

	/**
	 * Returns the total number of parameters used by this calculator
	 * 
	 * @return The total number of parameters used by this calculator
	 */
	public int numParameters() {
		return m_Parameters.length;
	}

	/**
	 * Returns the parameters for this calculator
	 * 
	 * @param template If this array is provided, the parameters will be copied into
	 *                 it. Otherwise a new array is created.
	 * @return The parameters for this calculator.
	 */
	public double[] getParameters(double[] template) {

		double[] returnArray = (template != null) ? template : new double[m_Parameters.length];

		System.arraycopy(m_Parameters, 0, returnArray, 0, m_Parameters.length);
		return returnArray;
	}

	/**
	 * Returns a new calculator copied with this one, with parameters given by the
	 * provided parameters.
	 * 
	 * @param parameters The parameters for the new calculator
	 * @return A copy of this calculator with the parameters changed to the given
	 *         parameters.
	 */
	public LikelihoodCalculator setParameters(double[] parameters) {

		return new LikelihoodCalculator(this, parameters);

	}

	/**
	 * Read the parameters from a parameters file (the old format, not the JSON
	 * format).
	 * 
	 * @param fileName The name of the file to be read
	 * @return An array of parameters
	 */
	public double[] readParameters(String fileName) {

		// Initializes lists for the read data
		ArrayList<Double> parameters = new ArrayList<Double>();
		ArrayList<int[]> oxidationStates = new ArrayList<int[]>();
		ArrayList<String> moleculeIDs = new ArrayList<String>();

		// Loop through the lines in the file and read the parameters.
		try {
			LineNumberReader reader = new LineNumberReader(new FileReader(fileName));
			String line = reader.readLine();
			while (line != null && line.trim().length() != 0) {
				String[] fields = line.split(" ");

				double firstParam = Double.parseDouble(fields[0]);
				parameters.add(firstParam);
				int[] statesForMolecule = new int[(fields.length - 1) / 2];
				String subStructureID = IonFactory.getIonTypeSymbol(fields[1]);

				for (int stateIndex = 2; stateIndex < fields.length; stateIndex += 2) {
					Ion ion = IonFactory.get(fields[stateIndex - 1]);
					statesForMolecule[stateIndex / 2 - 1] = (int) Math.round(ion.getOxidationState());
					double parameter = Double.parseDouble(fields[stateIndex]);
					parameters.add(parameter);
				}

				oxidationStates.add(statesForMolecule);
				moleculeIDs.add(subStructureID);
				line = reader.readLine();
			}
			reader.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		this.initializeOxidationStates(oxidationStates, moleculeIDs);
		return parameters.stream().mapToDouble(i -> i).toArray();
	}

	/**
	 * Read parameters from a JSON-formatted boundaries file
	 * 
	 * @param fileName The name of the JSON-formatted file
	 * @return And array of the read parameters
	 */
	public double[] readParametersFromBoundaryJSON(String fileName) {
		ObjectMapper mapper = new ObjectMapper();
		PotentialMapper potentialMapper = new PotentialMapper();

		try {
			RangeData[] rangeData = mapper.readValue(new File(fileName), RangeData[].class);
			ArrayList<Double> parameters = new ArrayList<Double>();
			ArrayList<int[]> oxidationStates = new ArrayList<int[]>();
			ArrayList<String> ionTypeSymbols = new ArrayList<String>();

			for (RangeData data : rangeData) {
				ionTypeSymbols.add(data.getIonTypeSymbol());
				oxidationStates.add(data.getOxidationStates());

				// The JSON file will have parameters in units of the mapped potential, so they
				// need to be converted
				double[] paramArray = boundariesToParameters(
						potentialMapper.fromMappedPotential(data.getRangeBoundaries()));
				for (double param : paramArray) {
					parameters.add(param);
				}
			}

			this.initializeOxidationStates(oxidationStates, ionTypeSymbols);
			return parameters.stream().mapToDouble(i -> i).toArray();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}

	/**
	 * Returns the minimum and maximum boundary values across all ion types
	 * 
	 * @return the minimum and maximum boundary values across all ion types
	 */
	public double[] getMinAndMaxBoundary() {
		double minBoundary = Double.POSITIVE_INFINITY;
		double maxBoundary = Double.NEGATIVE_INFINITY;
		for (int ionTypeIndex = 0; ionTypeIndex < m_IonTypes.length; ionTypeIndex++) {
			double[] boundaries = this.getBoundaries(ionTypeIndex);
			minBoundary = Math.min(minBoundary, boundaries[0]);
			maxBoundary = Math.max(maxBoundary, boundaries[boundaries.length - 1]);
		}
		return new double[] { minBoundary, maxBoundary };
	}

	/**
	 * Get the midpoint boundaries (lower and upper) for the given ion
	 * 
	 * @param ion The ion for which we want boundaries
	 * @return the midpoint boundaries (lower and upper) for the given ion
	 */
	public double[] getBoundaries(Ion ion) {
		int ionTypeIndex = this.getIonTypeIndex(ion.getIonType().getSymbol());
		double[] boundaries = this.getBoundaries(ionTypeIndex);
		for (int stateNum = 0; stateNum < boundaries.length - 1; stateNum++) {
			double oxidationState = this.getOxidationState(ionTypeIndex, stateNum);
			if (Math.abs(oxidationState - ion.getOxidationState()) < 1E-6) {
				return new double[] { boundaries[stateNum], boundaries[stateNum + 1] };
			}
		}
		return null;
	}

	/**
	 * Return all of the midpoint boundary values for the given ion type
	 * 
	 * @param ionTypeIndex The index of the ion type
	 * @return all of the midpoint boundary values for the given ion type
	 */
	public double[] getBoundaries(int ionTypeIndex) {

		int numOxidationStates = this.numOxidationStates(ionTypeIndex);
		if (numOxidationStates == 0) {
			return new double[0];
		}

		return parametersToBoundaries(this.getParameters(ionTypeIndex));
	}

	/**
	 * Converts a set of parameters for an ion type into oxidation state boundaries
	 * 
	 * @param parameters The parameters for the model
	 * @return The oxidation state boundaries for that ion type, ordered from lowest
	 *         to highest
	 */
	public static double[] parametersToBoundaries(double[] parameters) {

		double[] returnArray = new double[parameters.length];

		double center = parameters[0];
		double width = 0;
		for (int paramNum = 1; paramNum < parameters.length; paramNum++) {
			double param = parameters[paramNum];
			width += param * param;
		}
		returnArray[0] = center - (width / 2);

		for (int boundNum = 1; boundNum < returnArray.length; boundNum++) {
			double param = parameters[boundNum];
			returnArray[boundNum] = returnArray[boundNum - 1] + param * param;
		}

		return returnArray;

	}

	/**
	 * Converts a set of oxidation state boundaries for an ion type into parameters
	 * 
	 * @param boundaries The oxidation state boundaries
	 * @return The parameters representing those boundaries, suitable for use in the
	 *         conjugate gradient optimizer
	 */
	public static double[] boundariesToParameters(double[] boundaries) {

		double[] parameters = new double[boundaries.length];
		double width = boundaries[parameters.length - 1] - boundaries[0];
		double center = boundaries[0] + width / 2;

		parameters[0] = center;
		for (int paramNum = 1; paramNum < parameters.length; paramNum++) {
			double delta = boundaries[paramNum] - boundaries[paramNum - 1];
			parameters[paramNum] = Math.sqrt(delta);
		}

		return parameters;
	}

	/**
	 * Write the boundaries to a simple text file. This will not be the
	 * JSON-formatted file.
	 * 
	 * @param fileName The name of the file to be written.
	 */
	public void writeBoundaries(String fileName) {

		try {
			OutputStream stream = new FileOutputStream(fileName);
			this.writeBoundaries(stream);
			stream.flush();
			stream.close();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}

	}

	/**
	 * Write the boundaries to an output stream in a simple (non-JSON) format
	 * 
	 * @param stream The output stream to which boundaries should be written.
	 */
	public void writeBoundaries(OutputStream stream) {
		PrintWriter writer = new PrintWriter(stream);
		for (int moleculeIndex = 0; moleculeIndex < m_IonTypes.length; moleculeIndex++) {
			int numOxidationStates = this.numOxidationStates(moleculeIndex);
			if (numOxidationStates == 0) {
				continue;
			}

			double[] bounds = getBoundaries(moleculeIndex);
			writer.print(bounds[0] + " ");

			for (int stateNum = 0; stateNum < numOxidationStates; stateNum++) {
				int oxidationState = this.getOxidationState(moleculeIndex, stateNum);
				Ion ion = IonFactory.get(m_IonTypes[moleculeIndex], oxidationState);
				writer.print(ion.getSymbol() + " ");

				double bound = bounds[stateNum + 1];
				writer.print(bound + " ");
			}

			writer.println();
		}
		writer.flush();
	}

	/**
	 * Write the boundaries to a JSON-formatted file.
	 * 
	 * @param fileName The name of the file to be written.
	 */
	public void writeBoundaryJSON(String fileName) {

		try {
			OutputStream stream = new FileOutputStream(fileName);
			this.writeBoundaryJSON(stream);
			stream.flush();
			stream.close();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Write the boundaries in JSON format to an output stream.
	 * 
	 * @param stream The output stream to which the bouandaries should be written.
	 */
	public void writeBoundaryJSON(OutputStream stream) {

		PrintWriter writer = new PrintWriter(stream);
		Ion[] knownIons = this.getKnownOxidationStates();
		TreeSet<String> ionTypeSymbols = new TreeSet<String>();
		for (Ion ion : knownIons) {
			ionTypeSymbols.add(ion.getIonType().getSymbol());
		}

		RangeData[] rangeData = new RangeData[ionTypeSymbols.size()];
		int symbolNum = 0;
		for (String symbol : ionTypeSymbols) {
			rangeData[symbolNum++] = new RangeData(symbol, this);
		}

		ObjectMapper mapper = new ObjectMapper();
		try {
			mapper.writeValue(writer, rangeData);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Write files that can be used to visualize the boundaries
	 * 
	 * @param directoryName The name of the directory in which the files should be
	 *                      written.
	 * @param smoothCutoff  True if the plots should have smooth boundaries,
	 *                      reflecting the logistic functions, false if the plots
	 *                      should have sharp boundaries.
	 */
	public void writeDataFiles(String directoryName, boolean smoothCutoff) {

		File outDir = new File(directoryName);
		if (!outDir.exists()) {
			outDir.mkdirs();
		}

		for (int moleculeNumber = 0; moleculeNumber < m_IonTypes.length; moleculeNumber++) {

			int[] oxidationStates = m_OxidationStatesByIonType[moleculeNumber];
			if (oxidationStates.length < 1) {
				continue;
			}
			String moleculeID = m_IonTypes[moleculeNumber];
			String fileName = directoryName + "/" + moleculeID + ".csv";
			this.writeDataFile(moleculeID, fileName, smoothCutoff);
		}
	}

	/**
	 * Write files that can be used to visualize the boundaries
	 * 
	 * @param directoryName The name of the directory in which the files should be
	 *                      written.
	 * @param lowerBound    The lower electronic chemical potential bound for the
	 *                      plots
	 * @param upperBound    The upper electronic chemical potential bound for the
	 *                      plots
	 * @param increment     The width of the electronic chemical potential bins used
	 *                      to generate the plots
	 * @param smoothCutoff  True if the plots should have smooth boundaries,
	 *                      reflecting the logistic functions, false if the plots
	 *                      should have sharp boundaries.
	 */
	public void writeDataFiles(String directoryName, double lowerBound, double upperBound, double increment,
			boolean smoothCutoff) {

		File outDir = new File(directoryName);
		if (!outDir.exists()) {
			outDir.mkdirs();
		}

		for (int moleculeIndex = 0; moleculeIndex < m_IonTypes.length; moleculeIndex++) {

			int[] oxidationStates = m_OxidationStatesByIonType[moleculeIndex];
			if (oxidationStates.length <= 1) {
				continue;
			}
			String moleculeID = m_IonTypes[moleculeIndex];
			String fileName = directoryName + "/" + moleculeID + ".csv";
			this.writeDataFile(moleculeID, fileName, lowerBound, upperBound, increment, smoothCutoff);

		}

	}

	/**
	 * Write a file that can be used to visualize the boundaries for a particular
	 * ion type
	 * 
	 * @param ionType      The ion type
	 * @param fileName     The name of the file to be written
	 * @param smoothCutoff True if the plots should have smooth boundaries,
	 *                     reflecting the logistic functions, false if the plots
	 *                     should have sharp boundaries.
	 */
	public void writeDataFile(String ionType, String fileName, boolean smoothCutoff) {

		double minLowerBound = Double.POSITIVE_INFINITY;
		double maxUpperBound = Double.NEGATIVE_INFINITY;
		double nextMinLowerBound = Double.POSITIVE_INFINITY;
		double nextMaxUpperBound = Double.NEGATIVE_INFINITY;

		for (int moleculeIndex = 0; moleculeIndex < m_ParamIndexByState.length; moleculeIndex++) {
			double[] bounds = getBoundaries(moleculeIndex);
			if (bounds == null || bounds.length == 0) {
				continue;
			}
			double lowerBound = bounds[0];
			if (lowerBound < minLowerBound) {
				minLowerBound = lowerBound;
			} else if (lowerBound < nextMinLowerBound) {
				nextMinLowerBound = lowerBound;
			}

			double upperBound = bounds[bounds.length - 1];
			if (upperBound > maxUpperBound) {
				maxUpperBound = upperBound;
			} else if (upperBound > nextMaxUpperBound) {
				nextMaxUpperBound = upperBound;
			}
		}

		// So we can visualize the extreme values
		minLowerBound -= 4;
		maxUpperBound += 4;
		double increment = (maxUpperBound - minLowerBound) / 1000;
		this.writeDataFile(ionType, fileName, minLowerBound, maxUpperBound, increment, smoothCutoff);
	}

	/**
	 * Write a file that can be used to visualize the boundaries for a particular
	 * ion type
	 * 
	 * @param ionType      The ion type
	 * @param fileName     The name of the file to be written
	 * @param lowerBound   The lower electronic chemical potential bound for the
	 *                     plots
	 * @param upperBound   The upper electronic chemical potential bound for the
	 *                     plots
	 * @param increment    The width of the electronic chemical potential bins used
	 *                     to generate the plots
	 * @param smoothCutoff True if the plots should have smooth boundaries,
	 *                     reflecting the logistic functions, false if the plots
	 *                     should have sharp boundaries.
	 */
	public void writeDataFile(String ionType, String fileName, double lowerBound, double upperBound, double increment,
			boolean smoothCutoff) {

		int moleculeIndex = m_IonTypeIndices.get(ionType);
		int[] oxidationStates = m_OxidationStatesByIonType[moleculeIndex];

		// Long form
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
			String header = "potential,oxidationState,likelihood";
			writer.write(header);
			writer.newLine();
			for (double potential = lowerBound; potential <= upperBound + 1E-7; potential += increment) {

				double[] bounds = getBoundaries(moleculeIndex);
				double cutoff = bounds[0];
				double likelihood = smoothCutoff ? 1 / (1 + Math.exp((potential - cutoff)))
						: (potential < cutoff ? 1 : 0);
				writer.write(potential + ",0," + likelihood); // This is left of the lower bound
				writer.newLine();

				for (int stateNum = 0; stateNum < oxidationStates.length; stateNum++) {
					int oxidationState = oxidationStates[stateNum];
					cutoff = bounds[stateNum + 1];
					likelihood = smoothCutoff ? 1 / (1 + Math.exp((potential - cutoff))) : (potential < cutoff ? 1 : 0);
					writer.write(potential + "," + oxidationState + "," + likelihood);
					writer.newLine();
				}
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}

	/**
	 * Write the parameters to a text file (not JSON formatted).
	 * 
	 * @param fileName The name of the file to be written.
	 */
	public void writeParameters(String fileName) {

		try {
			OutputStream stream = new FileOutputStream(fileName);
			this.writeParameters(stream);
			stream.flush();
			stream.close();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}

	}

	/**
	 * Write the parameters to an output stream (not JSON formatted).
	 * 
	 * @param stream The output stream to which the parameters should be written.
	 */
	public void writeParameters(OutputStream stream) {
		PrintWriter writer = new PrintWriter(stream);
		for (int moleculeIndex = 0; moleculeIndex < m_IonTypes.length; moleculeIndex++) {
			String moleculeID = m_IonTypes[moleculeIndex];
			int numOxidationStates = this.numOxidationStates(moleculeIndex);
			if (numOxidationStates == 0) {
				continue;
			}
			double param = this.getParameter(moleculeIndex, 0);
			writer.print(param + " ");
			for (int stateNum = 0; stateNum < numOxidationStates; stateNum++) {
				int oxidationState = this.getOxidationState(moleculeIndex, stateNum);
				Ion ion = IonFactory.get(moleculeID, oxidationState);
				writer.print(ion + " ");
				param = this.getParameter(moleculeIndex, stateNum + 1);
				writer.print(param + " ");
			}
			writer.println();
		}
		writer.flush();
	}

	/**
	 * Sum of the maximum boundaries for all ion types, subtract the sum of the
	 * minimum boundaries, and report the difference.
	 * 
	 * @return The difference between the sum of the maximum boundaries and the sum
	 *         of the minimum boundaries.
	 */
	public double getSumOfSpreads() {
		double returnValue = 0;
		for (int ionIndex = 0; ionIndex < m_ParamIndexByState.length; ionIndex++) {
			double[] bounds = getBoundaries(ionIndex);
			returnValue += bounds[bounds.length - 1] - bounds[0];
		}
		return returnValue;
	}

	/**
	 * Return the likelihood score for a given ion at a give electronic chemical
	 * potential (fermiLevel)
	 * 
	 * @param fermiLevel The electronic chemical potential at which to calculate the
	 *                   likelihood score
	 * @param ion        The ion for which the likelihood score should be calculated
	 * @return the likelihood score for a given ion at a give electronic chemical
	 *         potential (fermiLevel)
	 */
	public double getLikelihood(double fermiLevel, Ion ion) {
		return getLikelihood(fermiLevel, new Ion[] { ion });
	}

	/**
	 * Finds the electronic chemical potential that maximizes the likelihood score
	 * for the given set of species
	 * 
	 * @param allSpecies The distinct species for which we want to optimize the
	 *                   likelihood score
	 * @return the electronic chemical potential that maximizes the likelihood score
	 *         for the given set of species
	 */
	public PotentialOptimizer optimizeLikelihood(Species[] allSpecies) {

		Ion[] ions = new Ion[0];
		for (int specNum = 0; specNum < allSpecies.length; specNum++) {
			Species species = allSpecies[specNum];
			if (species == Species.vacancy) {
				continue;
			}
			if (!m_IonTypeIndices.containsKey(species.getElementSymbol())) {
				return null;
			}
			Ion ion = IonFactory.get(species.getElementSymbol(), species.getOxidationState());
			ions = (Ion[]) ArrayUtils.appendElement(ions, ion);
		}
		return new PotentialOptimizer(ions);
	}

	/**
	 * Finds the electronic chemical potential that maximizes the likelihood score
	 * for the given set of ions
	 * 
	 * @param ions The distinct ions for which we want to optimize the likelihood
	 *             score
	 * @return the PotentialOptimizer containing the electronic chemical potential
	 *         that maximizes the likelihood score for the given set of ions
	 */
	public PotentialOptimizer optimizeLikelihood(Ion[] ions) {
		return new PotentialOptimizer(ions);
	}

	/**
	 * Finds the electronic chemical potential that maximizes the likelihood score
	 * for the given set of ions
	 * 
	 * @param ions              The distinct ions for which we want to optimize the
	 *                          likelihood score
	 * @param initialFermiLevel The chemical potential at which we should start the
	 *                          optimization process
	 * @return the electronic chemical potential that maximizes the likelihood score
	 *         for the given set of ions
	 */
	public PotentialOptimizer optimizeLikelihood(Ion[] ions, double initialFermiLevel) {
		return new PotentialOptimizer(ions, initialFermiLevel);
	}

	/**
	 * Returns the ion type index (used by this class) for ion types with the given
	 * symbol.
	 * 
	 * @param ionTypeSymbol The symbol for the ion type
	 * @return the ion type index (used by this class) for ion types with the given
	 *         symbol.
	 */
	public int getIonTypeIndex(String ionTypeSymbol) {
		return m_IonTypeIndices.get(ionTypeSymbol);
	}

	/**
	 * Calculates the likelihood score for a given set of ions at a given electronic
	 * chemical potential
	 * 
	 * @param fermiLevel The given electronic chemical potential
	 * @param allIons    The ions for which we are calculating the likelihood score
	 * @return the likelihood score for a given set of ions at a given electronic
	 *         chemical potential
	 */
	public double getLikelihood(double fermiLevel, Ion[] allIons) {

		double returnValue = 1;
		for (Ion ion : allIons) {
			int oxidationState = (int) Math.round(ion.getOxidationState());
			String ionStructureSymbol = ion.getIonType().getSymbol();
			int ionTypeIndex = m_IonTypeIndices.get(ionStructureSymbol);
			int[] oxidationStates = m_OxidationStatesByIonType[ionTypeIndex];
			int stateNum = Arrays.binarySearch(oxidationStates, oxidationState);

			double[] bounds = getBoundaries(ionTypeIndex);
			double lowerBound = bounds[stateNum];
			double upperBound = bounds[stateNum + 1];

			double lowerDelta = fermiLevel - lowerBound;
			double lowerValue = 1 / (1 + Math.exp(-lowerDelta));

			double upperDelta = fermiLevel - upperBound;
			double upperValue = 1 / (1 + Math.exp(upperDelta));

			returnValue *= Math.min(upperValue, lowerValue);
		}

		return returnValue;

	}

	/**
	 * Determine the electronic chemical potential that maximizes the likelihood
	 * score for the given entry
	 * 
	 * @param entry The entry for which we will maximize the likelihood score
	 * @return the electronic chemical potential that maximizes the likelihood score
	 *         for the given entry
	 */
	public PotentialOptimizer optimizeLikelihood(Entry entry) {
		return optimizeLikelihood(entry.getAllIons());
	}

	/**
	 * Determine the electronic chemical potential that maximizes the likelihood
	 * score for the given ion type IDS, where weights indicate the composition of
	 * each ion type.
	 * 
	 * @param ionTypeIDs The ion type IDs for which we are maximizing the likelihood
	 *                   score.
	 * @param weights    The compositions of each of the ion type ids, in the same
	 *                   order.
	 * @return the electronic chemical potential that maximizes the likelihood score
	 *         for the given ion type IDS, where weights indicate the composition of
	 *         each ion type.
	 */
	public PotentialOptimizer optimizeLikelihood(String[] ionTypeIDs, double[] weights) {
		return new PotentialOptimizer(this.getLikelyOxidationStates(ionTypeIDs, weights).getIons());
	}

	/**
	 * Get all possible oxidation state assignments for the given composition,
	 * subject to the given constraints.
	 * 
	 * @param composition          The composition for which we want oxidation
	 *                             states
	 * @param minAllowedLikelihood The minimum allowed likelihood score for the
	 *                             returned states. Set this to something non-zero
	 *                             to improve performance.
	 * @param minNumberToReturn    The minimum number of states to return, even if
	 *                             some of those states have a likelihood score
	 *                             below minAllowedLikelihood
	 * @return all possible oxidation state assignments for the given composition,
	 *         subject to the given constraints.
	 */
	public LikelihoodStateSet[] getAllOxidationStates(String composition, double minAllowedLikelihood,
			int minNumberToReturn) {
		Composition parser = new Composition(composition);
		return this.getAllOxidationStates(parser, minAllowedLikelihood, minNumberToReturn);
	}

	/**
	 * Returns an array possible oxidation state assignments for the given
	 * composition, subject to the given constraints.
	 * 
	 * @param composition          The composition for which we want oxidation
	 *                             states
	 * @param minAllowedLikelihood The minimum allowed likelihood score for the
	 *                             returned states. Set this to something non-zero
	 *                             to improve performance.
	 * @param minNumberToReturn    The minimum number of states to return, even if
	 *                             some of those states have a likelihood score
	 *                             below minAllowedLikelihood
	 * @return an array possible oxidation state assignments for the given
	 *         composition, subject to the given constraints.
	 */
	public LikelihoodStateSet[] getAllOxidationStates(Composition composition, double minAllowedLikelihood,
			int minNumberToReturn) {

		String[] moleculeIDs = composition.getSymbols().toArray(new String[0]);
		double[] counts = new double[moleculeIDs.length];
		for (int idNum = 0; idNum < moleculeIDs.length; idNum++) {
			counts[idNum] = composition.getCount(moleculeIDs[idNum]);
		}
		ArrayList<LikelihoodStateSet> aboveThreshholdList = new ArrayList<LikelihoodStateSet>();
		ArrayList<LikelihoodStateSet> belowThreshholdList = new ArrayList<LikelihoodStateSet>();

		int[] numStates = new int[moleculeIDs.length];
		Ion[][] allowedIons = new Ion[moleculeIDs.length][];
		for (int moleculeNum = 0; moleculeNum < moleculeIDs.length; moleculeNum++) {
			String moleculeID = moleculeIDs[moleculeNum];
			if (!m_IonTypeIndices.containsKey(moleculeID)) {
				Status.warning("No parameters available for molecule: " + moleculeID);
				return null;
			}
			int moleculeIndex = m_IonTypeIndices.get(moleculeID);
			numStates[moleculeNum] = this.numOxidationStates(moleculeIndex);
			allowedIons[moleculeNum] = new Ion[numStates[moleculeNum]];
			for (int stateNum = 0; stateNum < numStates[moleculeNum]; stateNum++) {
				double oxidationState = this.getOxidationState(moleculeIndex, stateNum);
				allowedIons[moleculeNum][stateNum] = IonFactory.get(moleculeID, oxidationState);
			}
		}
		MinAllowedLikelihoodFilter likelihoodThreshholdFilter = new MinAllowedLikelihoodFilter(this, allowedIons, 0);
		MixedValenceChargeBalanceFilter chargeBalanceFilter = new MixedValenceChargeBalanceFilter(allowedIons, counts);
		ArrayIndexer.Filter[] filters = new ArrayIndexer.Filter[] { likelihoodThreshholdFilter, chargeBalanceFilter };
		ArrayIndexer indexer = new ArrayIndexer(numStates);

		int[] stateIndices = indexer.getInitialState(filters);
		if (stateIndices == null) {
			return null;
		} // No match possible
		Ion[] baseIons = new Ion[moleculeIDs.length];
		LikelihoodStateSet result = null;
		do {

			double likelihood = likelihoodThreshholdFilter.getLastValidLikelihood();
			double fermiLevel = likelihoodThreshholdFilter.getLastValidFermiLevel();

			// Set the oxidation states and calculate the net charge
			double netCharge = 0;
			// double likelihood = 1;
			for (int moleculeNum = 0; moleculeNum < stateIndices.length; moleculeNum++) {
				int moleculeIndex = m_IonTypeIndices.get(moleculeIDs[moleculeNum]);
				int stateNum = stateIndices[moleculeNum];
				int oxidationState = this.getOxidationState(moleculeIndex, stateNum);
				baseIons[moleculeNum] = allowedIons[moleculeNum][stateNum];
				netCharge += oxidationState * counts[moleculeNum];
			}

			if (Math.abs(netCharge) < 1E-2) {
				result = new LikelihoodStateSet(baseIons, counts, likelihood, fermiLevel);
				double newMinLikelihood = addNewResult(result, aboveThreshholdList, belowThreshholdList,
						minNumberToReturn, minAllowedLikelihood);
				likelihoodThreshholdFilter.setMinAllowedLikelihood(newMinLikelihood);
			} else { // Consider mixed valence. We only need to consider lowering valence to cover
						// all cases
				for (int elementNum = 0; elementNum < moleculeIDs.length; elementNum++) {
					int moleculeIndex = m_IonTypeIndices.get(moleculeIDs[elementNum]);
					int stateIndex = stateIndices[elementNum];
					int oxidationState = this.getOxidationState(moleculeIndex, stateIndex);
					for (int decrement = 1; decrement <= stateIndex; decrement++) {
						int newStateIndex = stateIndex - decrement;
						int newOxidationState = this.getOxidationState(moleculeIndex, newStateIndex);
						double delta = newOxidationState - oxidationState;
						double maxChange = delta * counts[elementNum];
						if (netCharge + maxChange > 0) {
							continue;
						} // Can't reach zero
						Ion newIon = allowedIons[elementNum][newStateIndex];
						Ion[] allIons = (Ion[]) ArrayUtils.appendElement(baseIons, newIon);

						// Balance the charge
						double[] newCounts = new double[counts.length + 1];
						System.arraycopy(counts, 0, newCounts, 0, counts.length);
						newCounts[newCounts.length - 1] = -netCharge / delta;
						newCounts[elementNum] -= -netCharge / delta;
						if (Math.abs(newCounts[elementNum]) < 1E-7) { // We've essentially completely removed it
							continue;
						}
						// Recalculate the likelihood
						PotentialOptimizer optimizer = this.optimizeLikelihood(allIons);
						likelihood = optimizer.getMaxLikelihood();
						fermiLevel = optimizer.getOptimalFermiLevel();

						result = new LikelihoodStateSet(allIons, newCounts, likelihood, fermiLevel);
						double newMinLikelihood = addNewResult(result, aboveThreshholdList, belowThreshholdList,
								minNumberToReturn, minAllowedLikelihood);
						likelihoodThreshholdFilter.setMinAllowedLikelihood(newMinLikelihood);
					}
				}
			}

		} while (indexer.increment(stateIndices, filters));

		aboveThreshholdList.addAll(belowThreshholdList);
		return aboveThreshholdList.toArray(new LikelihoodStateSet[0]);
	}

	/**
	 * Used to keep track of which results should be included in the list of
	 * oxidation state assignments
	 * 
	 * @param result              The set of oxidation states to consider adding
	 * @param aboveThreshholdList THe list of known states above or equal to the
	 *                            minimum likelihood score
	 * @param belowThreshholdList The list of known states below the minimum
	 *                            likelihood score
	 * @param minCount            The minimum number of states that should be
	 *                            included
	 * @param minThreshhold       The minimum likelihood score to be considered for
	 *                            adding states to the aboveThreshholdList
	 * @return The likelihood score below which we no longer need to consider
	 *         results
	 */
	private double addNewResult(LikelihoodStateSet result, ArrayList<LikelihoodStateSet> aboveThreshholdList,
			ArrayList<LikelihoodStateSet> belowThreshholdList, int minCount, double minThreshhold) {
		double likelihood = result.getMaxLikelihood();

		// Keep track of whether this state is above or below the minimum threshhold
		if (likelihood >= minThreshhold) {
			aboveThreshholdList.add(result);
		} else {
			belowThreshholdList.add(result);
		}

		// If we don't have enough states yet, just keep adding more
		if (aboveThreshholdList.size() + belowThreshholdList.size() < minCount) {
			return 0;
		}

		// We have too many states, and there are some on the belowThreshooldList, so we
		// remove the least likely state from that list
		if ((belowThreshholdList.size() > 0) && (aboveThreshholdList.size() + belowThreshholdList.size() > minCount)) {
			int minIndex = getMinLikelihoodIndex(belowThreshholdList);
			belowThreshholdList.remove(minIndex);
		}

		// We have all the states we need, so we can get rid of everything in the future
		// below minThreshhold
		if (aboveThreshholdList.size() >= minCount) {
			return minThreshhold;
		}

		int minRemainingIndex = getMinLikelihoodIndex(belowThreshholdList);

		// We have enough states, so we don't need to accept anything lower than this.
		return belowThreshholdList.get(minRemainingIndex).getMaxLikelihood();
	}

	/**
	 * Return the index of the state set with the minimum likelihood score from the
	 * given list
	 * 
	 * @param list The given list
	 * @return the index of the state set with the minimum likelihood score from the
	 *         given list
	 */
	public int getMinLikelihoodIndex(List<LikelihoodStateSet> list) {

		int returnValue = -1;
		double minLikelihood = Double.POSITIVE_INFINITY;
		for (int setNum = 0; setNum < list.size(); setNum++) {
			LikelihoodStateSet set = list.get(setNum);
			if (set.getMaxLikelihood() < minLikelihood) {
				minLikelihood = set.getMaxLikelihood();
				returnValue = setNum;
			}
		}
		return returnValue;
	}

	@Override
	public LikelihoodStateSet getLikelyOxidationStates(String[] ionTypeIDs, double[] weights) {

		// Set up a branch and bound search
		int[] numStates = new int[ionTypeIDs.length];
		Ion[][] allowedIons = new Ion[ionTypeIDs.length][];
		for (int ionTypeNum = 0; ionTypeNum < ionTypeIDs.length; ionTypeNum++) {
			String ionTypeID = ionTypeIDs[ionTypeNum];
			if (!m_IonTypeIndices.containsKey(ionTypeID)) {
				Status.warning("No parameters available for molecule: " + ionTypeID);
				return null;
			}
			int ionTypeIndex = m_IonTypeIndices.get(ionTypeID);
			numStates[ionTypeNum] = this.numOxidationStates(ionTypeIndex);
			allowedIons[ionTypeNum] = new Ion[numStates[ionTypeNum]];
			for (int stateNum = 0; stateNum < numStates[ionTypeNum]; stateNum++) {
				double oxidationState = this.getOxidationState(ionTypeIndex, stateNum);
				allowedIons[ionTypeNum][stateNum] = IonFactory.get(ionTypeID, oxidationState);
			}
		}
		MinAllowedLikelihoodFilter minAllowedLikelihoodFilter = new MinAllowedLikelihoodFilter(this, allowedIons, 0);
		MixedValenceChargeBalanceFilter chargeBalanceFilter = new MixedValenceChargeBalanceFilter(allowedIons, weights);
		ArrayIndexer.Filter[] filters = new ArrayIndexer.Filter[] { minAllowedLikelihoodFilter, chargeBalanceFilter };
		ArrayIndexer indexer = new ArrayIndexer(numStates);

		int[] stateIndices = indexer.getInitialState(filters);
		if (stateIndices == null) {
			return null;
		} // No match possible
		Ion[] baseIons = new Ion[ionTypeIDs.length];
		LikelihoodStateSet result = null;
		do {

			double likelihood = minAllowedLikelihoodFilter.getLastValidLikelihood();
			double fermiLevel = minAllowedLikelihoodFilter.getLastValidFermiLevel();

			// Set the oxidation states and calculate the net charge
			double netCharge = 0;
			for (int ionTypeNum = 0; ionTypeNum < stateIndices.length; ionTypeNum++) {
				int ionTypeIndex = m_IonTypeIndices.get(ionTypeIDs[ionTypeNum]);
				int stateNum = stateIndices[ionTypeNum];
				int oxidationState = this.getOxidationState(ionTypeIndex, stateNum);
				// TODO why not just get the oxidation state from the ion?
				baseIons[ionTypeNum] = allowedIons[ionTypeNum][stateNum];
				netCharge += oxidationState * weights[ionTypeNum];
			}

			if (Math.abs(netCharge) < 1E-2) {
				result = new LikelihoodStateSet(baseIons, weights, likelihood, fermiLevel);
				minAllowedLikelihoodFilter.setMinAllowedLikelihood(likelihood);
			} else { // Consider mixed valence. We only need to consider lowering valence to cover
						// all cases
				for (int elementNum = 0; elementNum < ionTypeIDs.length; elementNum++) {
					int ionTypeIndex = m_IonTypeIndices.get(ionTypeIDs[elementNum]);
					int stateIndex = stateIndices[elementNum];
					int oxidationState = this.getOxidationState(ionTypeIndex, stateIndex);
					for (int decrement = 1; decrement <= stateIndex; decrement++) {
						int newStateIndex = stateIndex - decrement;
						int newOxidationState = this.getOxidationState(ionTypeIndex, newStateIndex);
						double delta = newOxidationState - oxidationState;
						double maxChange = delta * weights[elementNum];
						if (netCharge + maxChange > 0) {
							continue;
						} // Can't reach zero
						Ion newIon = allowedIons[elementNum][newStateIndex];
						Ion[] allIons = (Ion[]) ArrayUtils.appendElement(baseIons, newIon);

						// Balance the charge
						double[] newWeights = new double[weights.length + 1];
						System.arraycopy(weights, 0, newWeights, 0, weights.length);
						newWeights[newWeights.length - 1] = -netCharge / delta;
						newWeights[elementNum] -= -netCharge / delta;

						// Recalculate the likelihood
						PotentialOptimizer optimizer = this.optimizeLikelihood(allIons);
						likelihood = optimizer.getMaxLikelihood();
						fermiLevel = optimizer.getOptimalFermiLevel();

						if (likelihood > minAllowedLikelihoodFilter.getMinAllowedLikelihood()) {
							minAllowedLikelihoodFilter.setMinAllowedLikelihood(likelihood);
							result = new LikelihoodStateSet(allIons, newWeights, likelihood, fermiLevel);
						}
					}
				}
			}

		} while (indexer.increment(stateIndices, filters));

		return result;
	}

	/**
	 * Get a an array of all ions known to this calculator
	 * 
	 * @return a an array of all ions known to this calculator
	 */
	public Ion[] getKnownOxidationStates() {

		int numOxidationStates = 0;
		for (int ionTypeIndex = 0; ionTypeIndex < m_OxidationStatesByIonType.length; ionTypeIndex++) {
			numOxidationStates += m_OxidationStatesByIonType[ionTypeIndex].length;
		}

		Ion[] returnArray = new Ion[numOxidationStates];
		int returnIndex = 0;
		for (int ionTypeIndex = 0; ionTypeIndex < m_OxidationStatesByIonType.length; ionTypeIndex++) {
			String symbol = m_IonTypes[ionTypeIndex];
			for (int stateNum = 0; stateNum < m_OxidationStatesByIonType[ionTypeIndex].length; stateNum++) {
				int oxidationState = m_OxidationStatesByIonType[ionTypeIndex][stateNum];
				returnArray[returnIndex++] = IonFactory.get(symbol, oxidationState);
			}
		}
		return returnArray;

	}

	/**
	 * Get the ion type symbol corresponding to the given index
	 * 
	 * @param ionTypeIndex The given ion type index
	 * @return the ion type symbol corresponding to the given index
	 */
	public String getIonTypeSymbol(int ionTypeIndex) {
		return m_IonTypes[ionTypeIndex];
	}

	/**
	 * Returns the number of ion types known to this calculator
	 * 
	 * @return the number of ion types known to this calculator
	 */
	public int numIonTypes() {
		return m_OxidationStatesByIonType.length;
	}

	/**
	 * Returns the number of oxidation states known to this calculator for the ion
	 * type with the given index
	 * 
	 * @param ionTypeIndex The given ion type index
	 * @return the number of oxidation states known to this calculator for the ion
	 *         type with the given index
	 */
	public int numOxidationStates(int ionTypeIndex) {

		return m_OxidationStatesByIonType[ionTypeIndex].length;

	}

	/**
	 * Returns the oxidation states known to this calculator for the ion type with
	 * the given symbol
	 * 
	 * @param ionTypeSymbol The given symbol
	 * @return the oxidation states known to this calculator for the ion type with
	 *         the given symbol
	 */
	public int[] getOxidationStates(String ionTypeSymbol) {
		int ionTypeIndex = getIonTypeIndex(ionTypeSymbol);
		return m_OxidationStatesByIonType[ionTypeIndex].clone();
	}

	/**
	 * Returns the stateNum'th oxidation state for the ion type with the given index
	 * 
	 * @param ionTypeIndex The given ion type index
	 * @param stateNum     The number of the oxidation state to return
	 * @return the stateNum'th oxidation state for the ion type with the given index
	 */
	public int getOxidationState(int ionTypeIndex, int stateNum) {
		return m_OxidationStatesByIonType[ionTypeIndex][stateNum];
	}

	/**
	 * Return the paramNumForMolecule'th parameter for the ion type with the given
	 * index
	 * 
	 * @param ionTypeIndex       The given ion type index
	 * @param paramNumForIonType The parameter number for this ion type
	 * @return the paramNumForMolecule'th parameter for the ion type with the given
	 *         index
	 */
	public double getParameter(int ionTypeIndex, int paramNumForIonType) {
		int paramNum = m_ParamIndexByState[ionTypeIndex][paramNumForIonType];
		return m_Parameters[paramNum];
	}

	/**
	 * Returns the parameters for the ion type with the given index
	 * 
	 * @param ionTypeIndex The given ion type index
	 * @return the parameters for the ion type with the given index
	 */
	public double[] getParameters(int ionTypeIndex) {
		double[] returnArray = new double[m_ParamIndexByState[ionTypeIndex].length];
		for (int paramNum = 0; paramNum < returnArray.length; paramNum++) {
			returnArray[paramNum] = this.getParameter(ionTypeIndex, paramNum);
		}
		return returnArray;
	}

	/**
	 * This class is used to idnetify the electronic chemical potential that
	 * maximizes the likelihood score
	 * 
	 * @author timmueller
	 *
	 */
	public class PotentialOptimizer {

		// The set of ions for which we are maximizing the likelihood score
		private Ion[] m_AllIons;

		// Start searching from this electronic chemical potential
		private double m_InitialFermiLevel;

		// How large the initial step is in the search
		private double m_InitialStepSize = 0.5;

		// The electronic chemical potential that maximizes the likelihood score
		private double m_FinalFermiLevel = Double.NaN;

		// The likelihood score at the optimal electronic chemical potential (i.e. the
		// maximum score)
		private double m_FinalLikelihood = Double.NaN;

		/**
		 * Constructs a potential optimizer for the given set of ions
		 * 
		 * @param allIons The set of ions for which we are optimizing the electronic
		 *                chemical potential
		 */
		private PotentialOptimizer(Ion[] allIons) {

			m_AllIons = allIons.clone();

			m_InitialFermiLevel = findMiddleFermiLevel();
			if (Double.isNaN(m_FinalFermiLevel)) {
				this.optimize();
			}

		}

		/**
		 * Constructs a potential optimizer for the given set of ions with the
		 * electronic chemical potentially initially set to the given Fermi level.
		 * 
		 * @param allIons           The set of ions for which we are optimizing the
		 *                          electronic chemical potential
		 * @param initialFermiLevel The initial chemical potential for the optimization
		 *                          procedure
		 */
		private PotentialOptimizer(Ion[] allIons, double initialFermiLevel) {

			m_AllIons = allIons.clone();

			if (!Double.isInfinite(initialFermiLevel)) {
				m_InitialFermiLevel = initialFermiLevel;
			}

			double middleFermiLevel = findMiddleFermiLevel();
			if (Double.isNaN(m_FinalFermiLevel)) {
				m_InitialStepSize = Math.abs(m_InitialFermiLevel - middleFermiLevel);
				this.optimize();
			}

		}

		/**
		 * Returns the electronic chemical potential that maximizes the likelihood score
		 * 
		 * @return the electronic chemical potential that maximizes the likelihood score
		 */
		public double getOptimalFermiLevel() {
			return m_FinalFermiLevel;
		}

		/**
		 * Returns the likelihood score at the optimal electronic chemical potential
		 * (i.e. the maximum likelihood score)
		 * 
		 * @return the likelihood score at the optimal electronic chemical potential
		 *         (i.e. the maximum likelihood score)
		 */
		public double getMaxLikelihood() {
			return m_FinalLikelihood;
		}

		/**
		 * This method finds an initial guess for the optimal electronic chemical
		 * potential by finding a point that is at the mid point of the minimum upper
		 * bound and the maximum lower bound for all ions in this set.
		 * 
		 * @return an initial guess for the optimal electronic chemical potential that
		 *         is at the mid point of the minimum upper bound and the maximum lower
		 *         bound for all ions in this set.
		 */
		private double findMiddleFermiLevel() {
			// First find a good guess
			double maxLowerBound = Double.NEGATIVE_INFINITY;
			double minUpperBound = Double.POSITIVE_INFINITY;
			for (Ion ion : m_AllIons) {
				int oxidationState = (int) Math.round(ion.getOxidationState());
				String moleculeID = ion.getIonType().getSymbol();
				int moleculeIndex = m_IonTypeIndices.get(moleculeID);
				int[] oxidationStates = m_OxidationStatesByIonType[moleculeIndex];
				int stateNum = Arrays.binarySearch(oxidationStates, oxidationState);
				if (stateNum < 0) {
					return Double.NaN;
				}

				double[] bounds = getBoundaries(moleculeIndex);
				double lowerBound = bounds[stateNum];
				double upperBound = bounds[stateNum + 1];

				maxLowerBound = Math.max(maxLowerBound, lowerBound);
				minUpperBound = Math.min(upperBound, minUpperBound);
			}

			return (maxLowerBound + minUpperBound) / 2;

		}

		/**
		 * Finds the electronic chemical potential that maximizes the likelihood score
		 */
		private void optimize() {

			// Binary search
			double fermiLevel = m_InitialFermiLevel;
			double stepSize = m_InitialStepSize;
			Ion[] allIons = m_AllIons;

			// TODO use a line minimize here (e.g. GRLinearMinimizer2)
			double centerValue = getLikelihood(fermiLevel, allIons);
			while (stepSize > m_PotentialResolution) {

				double leftValue = getLikelihood(fermiLevel - stepSize, allIons);
				double rightValue = getLikelihood(fermiLevel + stepSize, allIons);

				if ((centerValue >= leftValue) && (centerValue >= rightValue)) {
					stepSize /= 2;
					continue;
				}

				if (leftValue > centerValue) {
					fermiLevel -= stepSize;
					centerValue = leftValue;
				} else {
					fermiLevel += stepSize;
					centerValue = rightValue;
				}

			}

			m_FinalFermiLevel = fermiLevel;
			m_FinalLikelihood = centerValue;
		}

	}

	/**
	 * Compares the atomic mass of two ion types, which is mostly useful for
	 * reporting purposes (it allows the ion types to be sorted by atomic mass). The
	 * comparison is done so that when sorted, the masses will be increasing.
	 * 
	 * @author timmueller
	 *
	 */
	private class AtomicMassComparator implements Comparator<String> {

		@Override
		public int compare(String ionType1ID, String ionType2ID) {

			IonType ionType1 = IonFactory.getKnownIonType(ionType1ID);
			IonType ionType2 = IonFactory.getKnownIonType(ionType2ID);

			if (ionType1 == null) {
				throw new RuntimeException("No substructure known for ID " + ionType1ID);
			}

			if (ionType2 == null) {
				throw new RuntimeException("No substructure known for ID " + ionType2ID);
			}

			return (int) Math.signum(ionType1.getAtomicWeight() - ionType2.getAtomicWeight());
		}
	}

}
