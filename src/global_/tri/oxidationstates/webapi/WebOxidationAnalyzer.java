package global_.tri.oxidationstates.webapi;

import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import global_.tri.oxidationstates.calculator.BondValenceCalculator;
import global_.tri.oxidationstates.calculator.likelihood.LikelihoodCalculator;
import global_.tri.oxidationstates.calculator.likelihood.LikelihoodStateSet;
import global_.tri.oxidationstates.ion.IonAssigner;
import global_.tri.oxidationstates.ion.IonFactory;
import global_.tri.oxidationstates.ion.KnownIonFinder;
import global_.tri.oxidationstates.ion.IonFactory.IonType;
import global_.tri.oxidationstates.ion.filters.PolyatomicSubstitutionFilter;
import global_.tri.oxidationstates.util.Composition;
import global_.tri.oxidationstates.util.CompositionParseException;
import matsci.Element;
import matsci.io.structure.CIF;
import matsci.io.vasp.POSCAR;
import matsci.structure.Structure;
import matsci.util.arrays.ArrayIndexer;
import matsci.util.arrays.ArrayUtils;

/**
 * This is the primary class for the Web API.  The main methods to call are:
 * 
 * getTableData():  This will return a JSON-formatted string containing all of the data required 
 * to build the table on the web app.  Two methods are provided:  one that can take a structure file, 
 * and one that can take a composition string
 * 
 * assignOxidationStates():  This will assign the selected oxidation states to the structure
 * 
 * @author timmueller
 *
 */
public class WebOxidationAnalyzer {
	
	private LikelihoodCalculator m_Calculator;
	private double m_MinAllowedLikelihood = 5E-3;
	private int m_MinListLength = 1;
	
	/**
	 * Initializes the analyzer.  Ideally this would only be done once, to avoid initialization overhead with each request.
	 * 
	 * @param jsonBoundaryFileLocation The location of the file containing the oxidation state analysis boundaries in JSON format.
	 * @param polyIonLocation The location of the directory containing the information on known polyatomic ions.
	 */
	public WebOxidationAnalyzer(String jsonBoundaryFileLocation, String polyIonLocation) {		
		IonFactory.loadPolyatomicIons(polyIonLocation); // Important to to this first to load the polyatomic ions
		m_Calculator = new LikelihoodCalculator(jsonBoundaryFileLocation, true);
	}
	
	/**
	 * Sets the minimum allowed likelihood.  No combinations of oxidation states with optimal likelihood
	 * below this will be returned.
	 * 
	 * @param value The minimum allowed likelihood.
	 */
	public void setMinAllowedLikelihood(double value) {
		m_MinAllowedLikelihood = value;
	}
	
	/**
	 * Returns the minimum allowed likelihood.  No combinations of oxidation states with optimal likelihood
	 * below this will be returned.
	 * 
	 * @return the minimum allowed likelihood.  The default value is 0.01.
	 */
	public double getMinAllowedLikelihood() {
		return m_MinAllowedLikelihood;
	}

	/**
	 * Set the minimum length of the list. At least this many rows will be returned, even if some 
	 * of them have likelihoods below {@link #getMinAllowedLikelihood()}. Fewer rows may be returned 
	 * if there are fewer total possible combinations of oxidation states.
	 * 
	 * @param value The minimum length of the list. At least this many rows will be returned, even if some 
	 * of them have likelihoods below {@link #getMinAllowedLikelihood()}. Fewer rows may be returned 
	 * if there are fewer total possible combinations of oxidation states.
	 */
	public void setMinListLength(int value) {
		m_MinListLength = value;
	}
	
	/**
	 * Return the minimum length of the list. At least this many rows will be returned, even if some 
	 * of them have likelihoods below {@link #getMinAllowedLikelihood()}. Fewer rows may be returned 
	 * if there are fewer total possible combinations of oxidation states.
	 * 
	 * @return the minimum length of the list. At least this many rows will be returned, even if some 
	 * of them have likelihoods below {@link #getMinAllowedLikelihood()}. Fewer rows may be returned 
	 * if there are fewer total possible combinations of oxidation states.
	 */
	public int minListLength() {
		return m_MinListLength;
	}
	
	/**
	 * Returns all combinations of oxidation states with optimal likelihood above the specified lower bound.
	 * 
	 * @param elementalComposition The composition of the structure. Specifying polyatomic ions in parenthesis is allowed,
	 * If any polyatomic ions are recognized, the oxidation state ranges for those ions will be used instead of those
	 * for the constitutive elements.
	 * @return A JSON-formatted string with all of the possible oxidation states for the table in the Web UI.
	 */
	public String getJSONFromComposition(String elementalComposition) {
		return this.getPageDataFromComposition(elementalComposition).toJSON();
	}	

	/**
	 * Returns all combinations of oxidation states with optimal likelihood above the specified lower bound.
	 * 
	 * @param compositionString The composition of the structure. Specifying polyatomic ions in parenthesis is allowed,
	 * If any polyatomic ions are recognized, the oxidation state ranges for those ions will be used instead of those
	 * for the constitutive elements.
	 * @return A TableData object with all of the possible oxidation states for the table in the Web UI.
	 */
	public PageData getPageDataFromComposition(String compositionString) {
		return this.getPageDataFromComposition(compositionString, false);
	}
	
	/**
	 * Returns all combinations of oxidation states with optimal likelihood above the provided lower bound.
	 * 
	 * @param compositionString The composition of the structure. Specifying polyatomic ions in parenthesis is allowed,
	 * If any polyatomic ions are recognized, the oxidation state ranges for those ions will be used instead of those
	 * for the constitutive elements.
	 * @param enumeratePossiblePolyanions Set this to "true" if you want the table to evaluate combinations of integer numbers of polyatomic 
	 * ions that could exist in a the given formula unit.  For example, if the given composition string is Li2VPO5, oxidation states
	 * will be evaluated for "Li O (PO4)", "Li 2O (PO3)", and "Li 5O P".  WARNING:  This can result in a table with a large number
	 * (>10,000) of rows.
	 * @return A TableData object with all of the possible oxidation states for the table in the Web UI.
	 */
	public PageData getPageDataFromComposition(String compositionString, boolean enumeratePossiblePolyanions) {
		
		Composition composition = null;
		Map<String, Double> knownIonComposition = null;
		try {
			// First we parse the composition and determine how many of each element we have.
			composition = new Composition(compositionString);
			knownIonComposition = composition.getKnownIonCompositionMap();
		} catch (CompositionParseException e) {
			return error(compositionString, "Failed to parse composition: " + e.getMessage());
		} catch (Exception e) {
			return error(compositionString, "Failed to parse composition.");
		}
		
		/**
		 * Handle the case where there are unknown elements 
		 */
		for (String symbol : knownIonComposition.keySet()) {
			if (!m_Calculator.hasParamsForIonType(symbol)) {
				return error(compositionString, "Species " + subscriptNumbers(symbol) + " is not supported by the oxidation analyzer.");
			}
		}
		
		/**
		 * Check to see if there are any polyatomic ions
		 */
		ArrayList<String> polySymbols = new ArrayList<String>();
		for (String symbol : knownIonComposition.keySet()) {
			IonType ionType = IonFactory.getKnownIonType(symbol);
			//if (m_PolyatomicIons.containsIonType(symbol)) {
			if (ionType != null && !ionType.isElement()) {
				polySymbols.add(symbol);
			}
		}	
		
		/**
		 * Check to see if there could be any polyatomic ions
		 */
		Map<String, Double> elementalCounts = composition.getElementalCompositionMap();
		ArrayList<String> possiblePolySymbols = new ArrayList<String>();
		//for (int ionNum = 0; ionNum < m_PolyatomicIons.numIons(); ionNum++) {
		for (IonType ionType : IonFactory.getKnownIonTypes().values()) {
			if (ionType.isElement()) {continue;}
			String symbol = ionType.getSymbol(); 
			if (polySymbols.contains(symbol)) {continue;} 
			Element[] elements = ionType.getElements();
			boolean matchesComposition = true;
			for (Element element : elements) {
				matchesComposition &= elementalCounts.containsKey(element.getSymbol());
			}
			if (matchesComposition) {possiblePolySymbols.add(symbol);}
		}	
		
		/**
		 * This is much faster than enumerating the ions
		 */
		if (!enumeratePossiblePolyanions) {
			TableRow[] rows = this.getTableRows(new Composition(knownIonComposition));
			rows = this.sortAndTrimTableRows(rows);
			PageData returnData = new PageData(compositionString, rows, m_Calculator);
			
			// Explain to the user why the table is empty.
			if (rows.length == 0) {
				returnData.addMessage("No assignment of common, non-zero oxidation states maintains charge neutrality.", false);
			}
			
			// Add a message explaining polyatomic ions
			if (polySymbols.size() > 0) {
				String message = "The following species are recognized as polyatomic ions and will be assigned their own oxidation state ranges: ";
				for (String symbol : polySymbols) {
					message = message + "(" + subscriptNumbers(symbol) + "), ";
				}
				message = message.substring(0, message.length() - 2);
				returnData.addMessage(message, false);
			}
			
			// Add a message explaining potential polyatomic ions
			if (possiblePolySymbols.size() > 0) {
				String message = "Some of the following polyatomic ions may be present in this material: ";
				for (String symbol : possiblePolySymbols) {
					message = message + "(" + subscriptNumbers(symbol) + "), ";
				}
				message = message.substring(0, message.length() - 2);
				message += ". If so, then expressing your composition in terms of the polyatomic ions will usually produce more accurate results.";
				returnData.addMessage(message, false);
			}
			return returnData;
		}
				
		/**
		 * We need to figure out how many of each polyatomic ion we could put in a 
		 * single formula unit of the given composition.
		 */
		IonType[] polyatomicTypes = IonFactory.getPolyatomicTypes();
		int[] maxAllowedIons = new int[polyatomicTypes.length]; 
		for (int ionNum = 0; ionNum < maxAllowedIons.length; ionNum++) {
			
			/**
			 * We get the elements in this ion, as well as the number of atoms
			 * of each element.
			 */
			Element[] ionElements = polyatomicTypes[ionNum].getElements();
			int[] ionCounts = polyatomicTypes[ionNum].getCounts();
			
			/**
			 * This value will tell us how many of this ion we can include
			 * in the unit cell before we use up all of the atoms of at least one
			 * element.
			 */
			double minRatio = Double.POSITIVE_INFINITY;
			
			/**
			 * Calculate the minimum ratio among all of the elements to figure out the maximum number of allowed ions.
			 */
			for (int elementNum = 0; elementNum < ionElements.length; elementNum++) {
				Element element = ionElements[elementNum];
				double structureCount = elementalCounts.containsKey(element.getSymbol()) ? elementalCounts.get(element.getSymbol()) : 0;
				int ionCount = ionCounts[elementNum];
				minRatio = Math.min(minRatio, structureCount / ionCount);
			}
			
			/**
			 *  We add an extra 1 here to account for the situation in which we have
			 *  zero polyatomic ions of this type.
			 */
			maxAllowedIons[ionNum] = (int) Math.floor(minRatio) + 1;  
		}
		
		// We'll keep adding our results to this collection
		ArrayList<TableRow> returnRows = new ArrayList<TableRow>();
		
		/**
		 * The filter is used to speed up the combinatorial search.  It will
		 * only allow combinations of polyatomic ions that don't use more atoms 
		 * of any given element than there is available in the formula unit.
		 */
		double[] countsByAtomicNumber = new double[Element.numKnownElements() + 1];
		for (String symbol : elementalCounts.keySet()) {
			int atomicNumber = Element.getElement(symbol).getAtomicNumber();
			countsByAtomicNumber[atomicNumber] = elementalCounts.get(symbol);
		}
		PolyatomicSubstitutionFilter filter = new PolyatomicSubstitutionFilter(countsByAtomicNumber, polyatomicTypes);
		ArrayIndexer.Filter[] filters = new ArrayIndexer.Filter[] {filter};
		
		/**
		 * This is a class I use for rapidly doing combinatorial searches.
		 * It's very useful.
		 */
		ArrayIndexer indexer = new ArrayIndexer(maxAllowedIons);
		int[] currState = indexer.getInitialState(filters);
		
		/**
		 * Each iteration through this loop analyzes a new composition and adds all
		 * generated sets of oxidation states to the collection of table rows.
		 */
		do {
			HashMap<String, Double> newCompositionMap = new HashMap<String, Double>();
			
			// Add the polyatomic ions to the composition
			for (int ionNum = 0; ionNum < currState.length; ionNum++) {
				int count = currState[ionNum];
				if (count == 0) {continue;}
				String ionSymbol = polyatomicTypes[ionNum].getSymbol();
				newCompositionMap.put(ionSymbol, 1.0 * count);
			}
			
			// Add the elements
			filter.addRemainingElementsForLastSeenState(newCompositionMap);
			
			// Now convert the result to table rows and add them to the collection
			Composition newComposition = new Composition(newCompositionMap);
			TableRow[] rows = this.getTableRows(newComposition);
			for (TableRow row : rows) {
				returnRows.add(row);
			}
		} while (indexer.increment(currState, filters)); // This jumps to the next combination that passes all filters.
		
		/**
		 * Convert the collection of rows to a table
		 */
		TableRow[] rowArray = returnRows.toArray(new TableRow[0]);
		
		// Sort the results and remove any unnecessary results
		rowArray = this.sortAndTrimTableRows(rowArray);

		PageData tableData = new PageData(compositionString, rowArray, m_Calculator);
		return tableData;
	}
	
	/**
	 * Returns an object containing the data needed to generate the web site pagefor the given structure.
	 * 
	 * @param structureFileString A string representation of a text file that contains an atomic structure.  
	 * Currently POSCAR and CIF are supported, and the format will be auto-detected.  
	 * @return The data needed to generate the oxidation state table for the given structure.
	 */
	public PageData getPageDataFromStructure(String structureFileString) {
		
		/**
		 * Convert the given file into a structure object
		 */
		Structure structure = null;
		try {
			structure = getStructure(structureFileString);
		} catch (DisorderedStructureException e) {
			return error("", "All sites must be fully occupied in the input structure.");
		} catch (Exception e) {
			return error("", "Unable to parse structure file.  Currently supported file types are POSCAR (VASP) and CIF.");
		}
		
		/**
		 * If oxidation states were provided in the structure, remove them.
		 */
		structure = structure.removeOxidationStates();
		
		/**
		 * Find the polyatomic ions in this structure.
		 */
		KnownIonFinder ionFinder = new KnownIonFinder(structure, true); //m_PolyatomicIons);
		String compositionString = ionFinder.getIonTypeComposition().getStandardizedCompositionString("");

		/**
		 * Check for unknown elements
		 */
		Element[] elements = structure.getDistinctElements();
		Element[] bvUnsupportedElements = new Element[0];
		for (Element element : elements) {
			String symbol = element.getSymbol();
			if (!m_Calculator.hasParamsForIonType(symbol)) {
				return error(compositionString, "Element " + symbol + " is not supported by the oxidation analyzer.");
			}
			if (!BondValenceCalculator.hasParametersForElement(element)) {
				bvUnsupportedElements = (Element[]) ArrayUtils.appendElement(bvUnsupportedElements, element);
			}
		}

		/**
		 * Given the composition, including polyatomic ions, calculate the most likely sets
		 * of oxidation states for the web table.
		 */
		TableRow[] tableRows = getTableRows(ionFinder);
		
		/**
		 * This will put the most likely oxidation states first and trim any unnecessary results.
		 */
		tableRows = this.sortAndTrimTableRows(tableRows);
				
		// Return the JSON string
		PageData returnData = new PageData(compositionString, tableRows, m_Calculator);
		
		// Explain to the user why the table is empty.
		if (tableRows.length == 0) {
			returnData.addMessage("No assignment of common, non-zero oxidation states maintains charge neutrality.", false);
		} else if (bvUnsupportedElements.length > 0) {

			/**
			 * Inform the user if there were some elements for which we did not have bond valence paraemters.
			 */
			String unsupportedElementMessage = "Could not calculate the global instability index or assign oxidation states to structure because no bond valence parameters are available for the following elements: ";
			int numRemainingElements = bvUnsupportedElements.length;
			for (Element element : bvUnsupportedElements) {
				unsupportedElementMessage += element.getSymbol();
				numRemainingElements--;
				if (numRemainingElements > 0) {
					unsupportedElementMessage += ", ";
				}
			}
			returnData.addMessage(unsupportedElementMessage, false);
		} else if (new BondValenceCalculator(structure).hasZeroSums()) {
			String message = "Some sites have a bond valence sum of zero. This might be caused by a problem with the input structure (e.g. a lattice parameter that is too large) and could affect GII calculations and the assignment of oxidation states to sites.";
			returnData.addMessage(message, false);
		}
		
		return returnData;		
		
	}
	
	private TableRow[] sortAndTrimTableRows(TableRow[] rows) {
		
		Arrays.sort(rows);
		
		int newLength = rows.length; 
		while ((newLength > m_MinListLength) && (rows[newLength - 1].getOptimalLikelihood() < m_MinAllowedLikelihood)) {
			newLength--;
		}
		
		return (TableRow[]) ArrayUtils.truncateArray(rows, newLength);
		
	}
	
	/**
	 * Returns a JSON-formatted string containing the data needed to generate the oxidation state
	 * table for the given structure.
	 * 
	 * @param structureFileString A string representation of a text file that contains an atomic structure.  
	 * Currently POSCAR and CIF are supported, and the format will be auto-detected.  
	 * @return A string containing the data needed to generate the oxidation state table for the given structure.
	 */
	public String getJSONFromStructure(String structureFileString) {
		return this.getPageDataFromStructure(structureFileString).toJSON();
	}
			
	/**
	 * Generate table rows for the structure covered by the ion finder, including polyatomci ions.
	 * 
	 * @param ionFinder An ion finder identifying any polyatomic ions in the structure
	 * @return  table rows for the structure covered by the ion finder, including polyatomci ions.
	 */
	private TableRow[] getTableRows(KnownIonFinder ionFinder) {
		
		/**
		 * First get the composition for the structure
		 */
		String composition = ionFinder.getIonTypeComposition().getStandardizedCompositionString();
		
		/**
		 * Now caclulate the corresponding oxidation states
		 */
		LikelihoodStateSet[] sets = m_Calculator.getAllOxidationStates(composition, m_MinAllowedLikelihood, m_MinListLength);
		
		/**
		 * If we couldn't find any oxidation states
		 */
		if (sets == null) {return new TableRow[0];} 
		
		/**
		 * For each set of oxidation states, generate the corresponding row in the table.
		 */
		TableRow[] returnRows = new TableRow[sets.length];
		for (int setNum = 0; setNum < sets.length; setNum++) {
			
			LikelihoodStateSet set = sets[setNum];
			
			// We need to assign oxidation states to sites to calculate the global instability index.
			IonAssigner assigner = new IonAssigner(ionFinder, set);
			double gii = assigner.getGlobalInstabilityIndex();
			
			// Create the new table row
			returnRows[setNum] = new TableRow(m_Calculator, set, gii);
			
			// We couldn't calculate the GII, so we probably couldn't assign oxidation states.
			if (Double.isNaN(gii)) {continue;}

			// We pre-generate all of the CIF strings, since we have the ion assignment handy
			StringWriter writer = new StringWriter();
			String cifString = null;
			try {
				new CIF(assigner.getStructure()).write(writer);
				cifString = writer.toString();
			} catch (IOException e) {}
			returnRows[setNum].setCIFString(cifString.trim());
		}
		
		// Return the table rows
		return returnRows;
	}
	
	/**
	 * Returns the table rows for the given composition
	 * 
	 * @param composition The given composition
	 * @return the table rows for the given composition
	 */
	private TableRow[] getTableRows(Composition composition) {
		
		LikelihoodStateSet[] sets = m_Calculator.getAllOxidationStates(composition, m_MinAllowedLikelihood, m_MinListLength);
		if (sets == null) {return new TableRow[0];}
		TableRow[] returnRows = new TableRow[sets.length];
		for (int setNum = 0; setNum < sets.length; setNum++) {
			returnRows[setNum] = new TableRow(m_Calculator, sets[setNum], Double.NaN);
		}
		return returnRows;
	}
		
	/**
	 * Returns a structure from parsing a CIF- or POSCAR-formatted input file
	 * 
	 * @param fileText The text of the input file
	 * @return The structure object
	 */
	private Structure getStructure(String fileText) {
		
		StringReader reader = new StringReader(fileText); 
		
		try {
			POSCAR poscar = new POSCAR(reader);
			reader.close();
			return new Structure(poscar);
		} catch (Exception e) {
			// Failed to parse, so try another type
		}
		
		CIF cif = null;
		try {
			reader.reset();
			cif = new CIF(reader);
			reader.close();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}		
		
		if (!cif.isOrdered(0)) {
			throw new DisorderedStructureException();
		}
		return new Structure(cif);
	}
	
	/**
	 * Add an error message to the page data and return it
	 * 
	 * @param composition The composition for which we are generating data
	 * @param message The error message
	 * @return The page data with the error message
	 */
	private PageData error(String composition, String message) {
		PageData returnData = new PageData(composition, m_Calculator);
		returnData.addMessage(message, true);
		return returnData;
	}
	
	/**
	 * Use this to format messages for display, so they don't have to be parsed client side.
	 * 
	 * @param symbol
	 * @return The given symbols with numbers in subscript
	 */
	private static String subscriptNumbers(String symbol) {
		
		symbol = symbol.replace("0", "₀");
		symbol = symbol.replace("1", "₁");
		symbol = symbol.replace("2", "₂");
		symbol = symbol.replace("3", "₃");
		symbol = symbol.replace("4", "₄");
		symbol = symbol.replace("5", "₅");
		symbol = symbol.replace("6", "₆");
		symbol = symbol.replace("7", "₇");
		symbol = symbol.replace("8", "₈");
		symbol = symbol.replace("9", "₉");
		
		return symbol;
	}
	
	/**
	 * An exception indicating that the given structure was disordered.
	 * 
	 * @author timmueller
	 *
	 */
	private class DisorderedStructureException extends RuntimeException {};

}
