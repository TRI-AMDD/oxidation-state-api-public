package global_.tri.oxidationstates.ion;

import java.io.File;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;

import global_.tri.oxidationstates.util.CompositionParseException;
import global_.tri.oxidationstates.util.IonTools;
import global_.tri.oxidationstates.util.VaspFileFilter;
import matsci.Element;
import matsci.Species;
import matsci.io.app.log.Status;
import matsci.structure.Structure;

/**
 * This class keeps track of available polyatomic and monatomci ions.
 * 
 * @author timmueller
 */
public class IonFactory {

	private static HashMap<String, IonType> KNOWN_ION_TYPES = new HashMap<String, IonType>();
	private static double OXIDATION_TOLERANCE = 0.001;

	private static HashSet<String> LOADED_POLYATOMIC_ION_DIRS = new HashSet<String>();

	// Initialize with known elements
	static {
		for (int atomicNumber = 1; atomicNumber <= Element.numKnownElements(); atomicNumber++) {
			Element element = Element.getElement(atomicNumber);
			KNOWN_ION_TYPES.put(element.getSymbol(), new IonType(element));
		}
	}

	/**
	 * Returns the ion matching the given symbol, or creates a new ion if no such
	 * ion exists. If the ion type for the provided ion is not known, return null.
	 * There is only a single ion object for each symbol.
	 * 
	 * @param symbol The symbol for the ion to return, including both a structure
	 *               type symbol and an oxidation state.
	 * @return the ion matching the given symbol, or creates a new ion if no such
	 *         ion exists. If the ion type for the provided ion is not known, return
	 *         null.
	 */
	public static Ion get(String symbol) {

		symbol = symbol.trim();

		double oxidationState = getOxidationState(symbol);
		if (!symbol.contains("(") && (oxidationState == 0) && !symbol.endsWith("0")) {
			symbol = "(" + symbol + ")";
		}

		String ionTypeSymbol = getIonTypeSymbol(symbol);
		return get(ionTypeSymbol, oxidationState);
	}

	/**
	 * Returns an ion corresponding to the given species. There is only a single ion
	 * object for each species.
	 * 
	 * @param species The species for which we want the matching ion.
	 * @return the ion corresponding to the given species.
	 */
	public static Ion get(Species species) {
		return get(species.getElementSymbol(), species.getOxidationState());
	}

	/**
	 * Returns an ion corresponding to the given element and oxidation state. There
	 * is only a single ion object for each combination of element and oxidation
	 * state.
	 * 
	 * @param element        The element for which we want the ion.
	 * @param oxidationState The oxidation state of the ion.
	 * @return an ion corresponding to the given element and oxidation state. There
	 *         is only a single ion object for each combination of element and
	 *         oxidation state.
	 */
	public static Ion get(Element element, double oxidationState) {
		return get(element.getSymbol(), oxidationState);
	}

	/**
	 * Returns an ion corresponding to the given ion type and oxidation state. There
	 * is only a single ion object for each combination of element and oxidation
	 * state. If no ionType matches the provided symbol, null is returned.
	 * 
	 * @param ionTypeSymbol  The symbol of the ion type for which we want the ion.
	 * @param oxidationState The oxidation state of the ion.
	 * @return an ion corresponding to the given ion type and oxidation state. There
	 *         is only a single ion object for each combination of element and
	 *         oxidation state. If no ionType matches the provided symbol, null is
	 *         returned.
	 */
	public static Ion get(String ionTypeSymbol, double oxidationState) {

		IonType ionType = KNOWN_ION_TYPES.get(ionTypeSymbol);
		if (ionType == null) {
			return null;
		}

		return ionType.getIon(oxidationState);

	}

	/**
	 * Creates a standardized symbol for an ion from an ion type symbol and and
	 * oxidation state.
	 * 
	 * @param ionTypeSymbol  The symbol for the ion type.
	 * @param oxidationState The oxidation state.
	 * @return a standardized symbol for an ion from an ion type symbol and and
	 *         oxidation state.
	 */
	public static String getIonSymbol(String ionTypeSymbol, double oxidationState) {
		ionTypeSymbol = ionTypeSymbol.trim();
		IonType ionType = KNOWN_ION_TYPES.get(ionTypeSymbol);
		if (ionType == null) {
			return null;
		}

		if (!ionType.isElement()) {
			return "(" + ionTypeSymbol + ")" + getOxidationString(oxidationState);
		}
		return ionTypeSymbol + getOxidationString(oxidationState);
	}

	/**
	 * Returns true if two oxidation states are within a given tolerance, which by
	 * default is 0.001. Returns false otherwise.
	 * 
	 * @param oxidation1 The first oxidation state to compare
	 * @param oxidation2 The second oxidation state to compare
	 * @return true if two oxidation states are within a given tolerance, which by
	 *         default is 0.001. Returns false otherwise.
	 */
	protected static boolean oxidationCloseEnough(double oxidation1, double oxidation2) {
		double diff = Math.abs(oxidation1 - oxidation2);
		return diff < OXIDATION_TOLERANCE;
	}

	/**
	 * Returns the ion type symbol corresponding to the given ion symbol.
	 * 
	 * @param ionSymbol A symbol for an ion.
	 * @return the ion type symbol corresponding to the given ion symbol.
	 */
	public static String getIonTypeSymbol(String ionSymbol) {
		if (ionSymbol.contains("(")) {
			int start = ionSymbol.indexOf("(");
			int end = ionSymbol.indexOf(")");
			try {
				return ionSymbol.substring(start + 1, end).trim();
			} catch (Exception e) {
				throw new CompositionParseException("Could not process symbol : " + ionSymbol);
			}
		}
		int start = 0;
		int end = 0;
		while (end < ionSymbol.length() && Character.isLetter(ionSymbol.charAt(end))) {
			end++;
		}
		return ionSymbol.substring(start, end).trim();
	}

	/**
	 * Returns the oxidation state of the ion with the given ion symbol.
	 * 
	 * @param ionSymbol A symbol for an ion.
	 * @return the oxidation state of the ion with the given ion symbol.
	 */
	public static double getOxidationState(String ionSymbol) {
		String endString = ionSymbol.substring(ionSymbol.length() - 1);
		double sign = (endString.equals("+") ? 1 : (endString.equals("-") ? -1 : 0));
		if (sign == 0) {
			return 0;
		}
		;

		int start = ionSymbol.indexOf(")") + 1;
		if (start == 0) {
			while (Character.isLetter(ionSymbol.charAt(start))) {
				start++;
			}
		}
		String stateString = ionSymbol.substring(start, ionSymbol.length() - 1).trim();
		double magnitude = stateString.length() == 0 ? 1 : Double.parseDouble(stateString);
		return sign * magnitude;
	}

	/**
	 * Returns a standardized string representation of the given oxidation state
	 * 
	 * @param oxidationState An oxidation state
	 * @return a standardized string representation of the given oxidation state
	 */
	public static String getOxidationString(double oxidationState) {
		if (oxidationState == 0) {
			return "";
		}
		double absValue = Math.abs(oxidationState);
		boolean isPositive = (oxidationState > 0);
		int intState = (int) Math.round(absValue);
		if (oxidationCloseEnough(intState, absValue)) {
			return isPositive ? intState + "+" : intState + "-";
		}
		return isPositive ? absValue + "+" : absValue + "-";
	}

	/**
	 * Loads into memory the polyatomic ions whose structure files are provided in
	 * the given top-level directory. The top-level directory should contain
	 * subdirectories, each of which has a name that matches an ion type. Each of
	 * these subdirectories should contain POSCAR-formatted files with a name
	 * matching an oxidation state with the extension of ".vasp". The contents of
	 * each of these files should be a representation of the structure of the ion
	 * with the matching ion type and oxidation state.
	 * 
	 * This method should usually be called near the beginning of the program,
	 * before any methods that use polyatomic ions are called, so that those methods
	 * have access to the structures of known polyatomic ions.
	 * 
	 * If this method is called more than once on the same top-level directory, then
	 * the ions will only be loaded into memory the first time it is called.
	 * 
	 * @param ionTypesDirName The top-level directory containing information about
	 *                        polyatomic ion structures. The top-level directory
	 *                        should contain subdirectories, each of which has a
	 *                        name that matches an ion type. Each of these
	 *                        subdirectories should contain POSCAR-formatted files
	 *                        with a name matching an oxidation state with the
	 *                        extension of ".vasp". The contents of each of these
	 *                        files should be a representation of the structure of
	 *                        the ion with the matching ion type and oxidation
	 *                        state.
	 */
	public static void loadPolyatomicIons(String ionTypesDirName) {

		ionTypesDirName = new File(ionTypesDirName).getAbsolutePath();
		if (LOADED_POLYATOMIC_ION_DIRS.contains(ionTypesDirName)) {
			Status.detail("Already loaded polyatomic ions from " + ionTypesDirName
					+ ", so they will not be loaded again from the same location.");
			return;
		}
		LOADED_POLYATOMIC_ION_DIRS.add(ionTypesDirName);

		File[] ionTypeDirs = new File(ionTypesDirName).listFiles();
		for (File ionTypeDir : ionTypeDirs) {
			String ionTypeSymbol = ionTypeDir.getName();
			if (KNOWN_ION_TYPES.containsKey(ionTypeSymbol)) {
				throw new RuntimeException("Loading ion type " + ionTypeSymbol + " from " + ionTypeDir.getAbsolutePath()
						+ " but name is already in use.");
			}

			IonType ionType = new IonType(ionTypeSymbol);
			KNOWN_ION_TYPES.put(ionTypeSymbol, ionType); // Need to do this first so that this is recognized as an ion
															// type when building ions.
			ionType.loadIonStructures(ionTypesDirName);
		}
	}

	/**
	 * Returns an array of known ion types for polyatomic ions.
	 * 
	 * @return an array of known ion types for polyatomic ions.
	 */
	public static IonType[] getPolyatomicTypes() {

		TreeMap<String, IonType> treeMap = new TreeMap<String, IonType>();
		for (IonType ionType : KNOWN_ION_TYPES.values()) {
			if (ionType.isElement()) {
				continue;
			}
			treeMap.put(ionType.getSymbol(), ionType);
		}

		return (IonType[]) treeMap.values().toArray(new IonType[0]);

	}

	/**
	 * Returns a map in which the keys are the ion type symbols and the values are
	 * the corresponding ion types for all known ion types, monatomic and
	 * polyatomic.
	 * 
	 * @return a map in which the keys are the ion type symbols and the values are
	 *         the corresponding ion types for all known ion types, monatomic and
	 *         polyatomic.
	 * 
	 */
	public static Map<String, IonType> getKnownIonTypes() {
		return Collections.unmodifiableMap(KNOWN_ION_TYPES);
	}

	/**
	 * Returns the ion type for the given symbol, or null if no such ion type
	 * exists.
	 * 
	 * @param ionTypeSymbol The symbol for the ion type
	 * @return the ion type for the given symbol, or null if no such ion type
	 *         exists.
	 */
	public static IonType getKnownIonType(String ionTypeSymbol) {
		return KNOWN_ION_TYPES.get(ionTypeSymbol);
	}

	/**
	 * Returns a map in which the keys are integer oxidation states and the values
	 * are the representative structures for ions of the given ion type with the
	 * corresponding oxidation state. Only integer oxidation states for which the
	 * representative structure is known will be included in the returned map.
	 * 
	 * @param ionTypeSymbol The symbol for the the ion type.
	 * @return a map in which the keys are integer oxidation states and the values
	 *         are the representative structures for ions of the given ion type with
	 *         the corresponding oxidation state.
	 */
	public static Map<Integer, Structure> getKnownRepresentativeStructures(String ionTypeSymbol) {

		HashMap<Integer, Structure> returnMap = new HashMap<Integer, Structure>();

		IonType ionType = KNOWN_ION_TYPES.get(ionTypeSymbol);
		if (ionType == null) {
			return returnMap;
		}

		Map<String, Ion> knownIons = ionType.getKnownIons();
		for (Ion ion : knownIons.values()) {
			if (!ion.hasIntegerOxidationState()) {
				continue;
			} // For this method we only return integer oxidation states
			if (ion.getRepresentativeStructure() == null) {
				continue;
			} // We don't know the structure
			returnMap.put((int) Math.round(ion.getOxidationState()), ion.getRepresentativeStructure());
		}

		return returnMap;
	}

	/**
	 * Returns a map in which the keys are known ion types and the values are the
	 * corresponding representative structures. This method only returns structures
	 * for ions with integer oxidation states.
	 * 
	 * @param polyatomic True if only polyatomic ions should be included in the
	 *                   returned map, false if all known ions should be included in
	 *                   the returned map.
	 * @return a map in which the keys are known ion types and the values are the
	 *         corresponding representative structures. This method only returns
	 *         structures for ions with integer oxidation states.
	 */
	public static HashMap<String, Structure> getAllKnownRepresentativeStructures(boolean polyatomic) {

		HashMap<String, Structure> returnMap = new HashMap<String, Structure>();
		for (IonType ionType : KNOWN_ION_TYPES.values()) {
			if (polyatomic && ionType.isElement()) {
				continue;
			}
			String typeSymbol = ionType.getSymbol();
			Map<Integer, Structure> representativeStructures = getKnownRepresentativeStructures(typeSymbol);
			for (Integer oxidationState : representativeStructures.keySet()) {
				String ionSymbol = IonFactory.getIonSymbol(typeSymbol, oxidationState);
				Structure structure = representativeStructures.get(oxidationState);
				returnMap.put(ionSymbol, structure);
			}
		}
		return returnMap;
	}

	/**
	 * This class represents an ion type, which is an atom (or cluster of atoms)
	 * that can be assigned an oxidation state.
	 * 
	 * @author timmueller
	 *
	 */
	public static class IonType {

		private String m_Symbol;
		private int m_NumSites;
		private Element[] m_Elements;
		private int[] m_Counts;
		private int m_NumPeriodicDimensions = -1;
		private boolean m_IsElement;
		private HashMap<String, Ion> m_Ions = new HashMap<String, Ion>();

		/**
		 * Creates an ion type for an atom of the given element.
		 * 
		 * @param element THe given element
		 */
		private IonType(Element element) {
			m_Symbol = element.getSymbol();
			m_Elements = new Element[] { element };
			m_Counts = new int[] { 1 };
			m_IsElement = true;
			m_NumPeriodicDimensions = 0;
		}

		/**
		 * Creates an ion type with the given symbol
		 * 
		 * @param typeSymbol The ion type symbol
		 */
		private IonType(String typeSymbol) {
			m_Symbol = typeSymbol;
			m_IsElement = Element.isKnownElement(typeSymbol);
		}

		/**
		 * Reads all representative structures for different oxidation states of a
		 * polyatomic ion.
		 * 
		 * @param baseDirName The base directory containing the representative
		 *                    structures. The top-level directory should contain
		 *                    subdirectories, each of which has a name that matches an
		 *                    ion type. Each of these subdirectories should contain
		 *                    POSCAR-formatted files with a name matching an oxidation
		 *                    state with the extension of ".vasp". The contents of each
		 *                    of these files should be a representation of the structure
		 *                    of the ion with the matching ion type and oxidation state.
		 * 
		 */
		private void loadIonStructures(String baseDirName) {
			m_IsElement = false;

			File dir = new File(baseDirName + "/" + m_Symbol);
			for (File structFile : dir.listFiles(new VaspFileFilter())) {
				String fileName = structFile.getName();
				if (fileName.equals("mean.vasp")) {
					continue;
				}
				String oxidationString = fileName.replace(".vasp", "");
				Integer oxidationState = Integer.parseInt(oxidationString);

				Structure structure = IonTools.loadIonStructureFromFile(structFile.getAbsolutePath());

				if (m_Elements == null) {
					m_Elements = structure.getDistinctElements();
					m_Counts = new int[m_Elements.length];
					for (int elementNum = 0; elementNum < m_Elements.length; elementNum++) {
						Element element = m_Elements[elementNum];
						m_Counts[elementNum] = structure.numDefiningSitesWithElement(element);
						m_NumPeriodicDimensions = structure.numPeriodicDimensions();

						/**
						 * For now we assume that the total number of atoms needed to represent each
						 * oxidation state is the same. This might not be true for periodic structures.
						 */
						m_NumSites += m_Counts[elementNum];
					}
				} else {
					if (m_NumSites != structure.numDefiningSites()) {
						throw new RuntimeException(
								"Found structures with differing numbers of defining sites for ion type " + m_Symbol
										+ ".");
					}
					if (m_NumPeriodicDimensions != structure.numPeriodicDimensions()) {
						throw new RuntimeException(
								"Found structures with differing numbers of periodic dimensions for ion type "
										+ m_Symbol + ".");
					}
					for (int elementNum = 0; elementNum < m_Elements.length; elementNum++) {
						Element element = m_Elements[elementNum];
						int count = m_Counts[elementNum];
						if (structure.numDefiningSitesWithElement(element) != count) {
							throw new RuntimeException(
									"Structures representing ion type " + m_Symbol + " have different compositions.");
						}
					}
				}

				Ion ion = new Ion(this, oxidationState, structure);

				m_Ions.put(getOxidationString(oxidationState), ion);
			}
		}

		/**
		 * Return a map in which the keys are ion symbols and the values are the known
		 * ions of this type.
		 * 
		 * @return a map in which the keys are ion symbols and the values are the known
		 * ions of this type.
		 */
		public Map<String, Ion> getKnownIons() {
			return Collections.unmodifiableMap(m_Ions);
		}

		/**
		 * Returns the distinct elements in this ion type
		 * 
		 * @return the distinct elements in this ion type
		 */
		public Element[] getElements() {
			return m_Elements.clone();
		}

		/**
		 *
		 * Returns the counts of each element in this ion type, in the same order as the
		 * array returned by {@link getElements()}.
		 *
		 * @return the counts of each element in this ion type, in the same order as the
		 *         array returned by {@link getElements()}.
		 */
		public int[] getCounts() {
			return m_Counts.clone();
		}

		/**
		 * Returns true if this ion type is an element, and false if otherwise (e.g. for
		 * polyatomic ion types).
		 * 
		 * @return true if this ion type is an element, and false if otherwise (e.g. for
		 *         polyatomic ion types).
		 */
		public boolean isElement() {
			return m_IsElement;
		}

		/**
		 * Return the symbol for this ion type. It should be unique.
		 * 
		 * @return the symbol for this ion type. It should be unique.
		 */
		public String getSymbol() {
			return m_Symbol;
		}

		/**
		 * Some ions are crystalline. This method returns the number of periodic vectors
		 * in the Bravais lattice for this ion type.
		 * 
		 * @return the number of periodic vectors in the Bravais lattice for this ion
		 *         type.
		 */
		public int numPeriodicDimensions() {
			return m_NumPeriodicDimensions;
		}

		/**
		 * Returns the atomic weight for this ion type. All returned values are per unit
		 * cell for periodic ions
		 * 
		 * @return the atomic weight for this ion type. All returned values are per unit
		 *         cell for periodic ions
		 */
		public double getAtomicWeight() {

			double atomicWeight = 0;
			for (int elementNum = 0; elementNum < m_Elements.length; elementNum++) {
				atomicWeight += m_Elements[elementNum].getAtomicWeight() * m_Counts[elementNum];
			}
			return atomicWeight;

		}

		/**
		 * Returns the ion of this type with the given oxidation state. Returns an
		 * existing object if one exists, or creates (and stores) a new ion object if
		 * one does not. This ensures that there is one global object for each ion.
		 * 
		 * @param oxidationState The given oxidation state.
		 * @return the ion of this type with the given oxidation state. Returns an
		 *         existing object if one exists, or creates (and stores) a new ion
		 *         object if one does not. This ensures that there is one global object
		 *         for each ion.
		 */
		public Ion getIon(double oxidationState) {
			String oxidationString = getOxidationString(oxidationState);
			Ion returnIon = m_Ions.get(oxidationString);
			if (returnIon != null) {
				return returnIon;
			}

			returnIon = new Ion(this, oxidationState);
			m_Ions.put(oxidationString, returnIon); // TODO move this into ion constructor?
			return returnIon;
		}

	}

	/**
	 * An Ion is a combination of an ion type (e.g. element or cluster of elements)
	 * with an oxidation state.
	 * 
	 * @author timmueller
	 *
	 */
	public static class Ion {

		private IonType m_IonType;
		private Structure m_RepresentativeStructure;
		private String m_Symbol;
		private double m_OxidationState;

		/**
		 * Creates an ion with the given ion type and oxidation state.
		 * 
		 * @param ionType        The given on type
		 * @param oxidationState The given oxidation state
		 */
		private Ion(IonType ionType, double oxidationState) {
			m_IonType = ionType;
			m_Symbol = IonFactory.getIonSymbol(ionType.getSymbol(), oxidationState);
			m_OxidationState = oxidationState;
			m_RepresentativeStructure = ionType.isElement() ? ionType.getElements()[0].getAtomStructure() : null;
		}

		/**
		 * Creates an ion with the given ion type, oxidation state, and representative
		 * structure.
		 * 
		 * @param ionType        The given on type
		 * @param oxidationState The given oxidation state
		 * @param structure      A representative structure for this ion.
		 */
		private Ion(IonType ionType, double oxidationState, Structure structure) {
			m_IonType = ionType;
			m_RepresentativeStructure = structure;
			m_Symbol = IonFactory.getIonSymbol(ionType.getSymbol(), oxidationState);
			m_OxidationState = oxidationState;
		}

		/**
		 * Returns the representative structure for this ion.
		 * 
		 * @return the representative structure for this ion
		 */
		public Structure getRepresentativeStructure() {
			return m_RepresentativeStructure;
		}

		/**
		 * Returns the oxidation state for this ion.
		 * 
		 * @return the oxidation state for this ion.
		 */
		public double getOxidationState() {
			return m_OxidationState;
		}

		/**
		 * Returns the ion type for this ion.
		 * 
		 * @return the ion type for this ion.
		 */
		public IonType getIonType() {
			return m_IonType;
		}

		/**
		 * Returns the symbol for this ion.
		 * 
		 * @return the symbol for this ion.
		 */
		public String getSymbol() {
			return m_Symbol;
		}

		@Override
		public String toString() {
			return this.getSymbol();
		}

		/**
		 * Returns true if this is a polyatomic ion (as determined by the representative
		 * structure), and false otherwise.
		 * 
		 * @return true if this is a polyatomic ion (as determined by the representative
		 *         structure), and false otherwise.
		 */
		public boolean isPolyatomic() {
			return (m_RepresentativeStructure.numDefiningSites() > 1);
		}

		/**
		 * Returns true if the oxidation state is within the oxidation tolerance (0.001
		 * by default) of an integer, and false otherwise.
		 * 
		 * @return true if the oxidation state is within the oxidation tolerance (0.001
		 *         by default) of an integer, and false otherwise.
		 */
		public boolean hasIntegerOxidationState() {
			return oxidationCloseEnough(m_OxidationState, Math.round(m_OxidationState));
		}
	}

}
