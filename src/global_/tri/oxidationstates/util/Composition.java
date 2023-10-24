package global_.tri.oxidationstates.util;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import global_.tri.oxidationstates.ion.IonFactory;
import global_.tri.oxidationstates.ion.IonFactory.Ion;
import matsci.Element;
import matsci.util.MSMath;

/**
 * A utility class to represent a composition
 * 
 * @author timmueller
 *
 */
public class Composition {

	private final TreeMap<String, Double> m_CompositionMap = new TreeMap<String, Double>();

	/**
	 * Create a composition from a map in which the keys are ion type symbols and
	 * the values are amounts.
	 * 
	 * @param map a map in which the keys are symbols and the values are amounts.
	 */
	public Composition(Map map) {

		this.addSymbols(map, 1);

	}

	/**
	 * Create a composition from the given symbols and counts.
	 * 
	 * @param symbols The symbols for tne ion types in this composition (or similar
	 *                molecular units)
	 * @param counts  The amounds for each symbol. The indices of this array
	 *                correspond to those for the symbols array.
	 */
	public Composition(String[] symbols, int[] counts) {
		for (int countNum = 0; countNum < counts.length; countNum++) {
			this.addSymbol(symbols[countNum], counts[countNum]);
		}
	}

	/**
	 * Create a composition from the given symbols and counts.
	 * 
	 * @param symbols The symbols for tne ion types in this composition (or similar
	 *                molecular units)
	 * @param counts  The amounds for each symbol. The indices of this array
	 *                correspond to those for the symbols array.
	 */
	public Composition(String[] symbols, double[] counts) {
		for (int countNum = 0; countNum < counts.length; countNum++) {
			this.addSymbol(symbols[countNum], counts[countNum]);
		}
	}

	/**
	 * Create a composition from the given composition string
	 * 
	 * @param composition A string representation of the composition
	 */
	public Composition(String composition) {

		this.parseComposition(composition);

	}

	/**
	 * Adds the composition in the given map to this one, where all quantities will
	 * be multiplied by the group count.
	 * 
	 * @param map        a map in which the keys are symbols and the values are
	 *                   amounts.
	 * @param groupCount The composition in the given map will be multiplied by this
	 *                   factor before being added to this composition.
	 */
	private void addSymbols(Map map, double groupCount) {

		for (Object symbol : map.keySet()) {
			double symbolCount = Double.parseDouble(map.get(symbol).toString());
			this.addSymbol((String) symbol, groupCount * symbolCount);
		}

	}

	/**
	 * Adds the given symbol to this composition.
	 * 
	 * @param symbol The symbol to be added.
	 * @param count  The quantity of the given symbol to be added.
	 */
	private void addSymbol(String symbol, double count) {

		if (Math.abs(count) < 1E-7) {
			return;
		}
		if (m_CompositionMap.containsKey(symbol)) {
			count += m_CompositionMap.get(symbol);
		}
		m_CompositionMap.put(symbol.trim(), count);

	}

	/**
	 * Parses a composition string to create the composition. Nested parenthesis and
	 * polyatomic ions are allowed. All polyatomic ions should be indicated with the
	 * appropriate symbol enclosed in parenthesis.
	 * 
	 * @param compositionString
	 */
	private void parseComposition(String compositionString) {

		// compositionString = compositionString.trim();
		compositionString = compositionString.replaceAll("\\s", ""); // Remove all whitespace. TODO penalize species
																		// separated by whitespace?

		Set<String> knownSymbols = IonFactory.getKnownIonTypes().keySet();

		int charIndex = 0;
		Character currChar = nextCharacter(compositionString, charIndex++);
		while (currChar != null) {

			Composition subComposition = null;
			if (currChar == '(') {
				int parenthesisCount = 1;
				String subCompositionString = "";
				currChar = nextCharacter(compositionString, charIndex++);
				while (currChar != null) {
					if (currChar == '(') {
						parenthesisCount++;
					}
					if (currChar == ')') {
						parenthesisCount--;
					}
					if (parenthesisCount == 0) {
						break;
					}
					subCompositionString += currChar;
					currChar = nextCharacter(compositionString, charIndex++);
				}
				if (parenthesisCount > 0) {
					throw new CompositionParseException(
							"Found unclosed parenthesis in composition " + compositionString + ".");
				}
				currChar = nextCharacter(compositionString, charIndex++); // Get past closing parenthesis
				subCompositionString = subCompositionString.trim();
				if (knownSymbols.contains(subCompositionString)) {
					subComposition = new Composition(new String[] { subCompositionString }, new int[] { 1 });
				} else {
					subComposition = new Composition(subCompositionString);
				}
			} else {
				String moleculeID = "" + currChar;
				currChar = nextCharacter(compositionString, charIndex++);
				while (currChar != null && Character.isLowerCase(currChar)) {
					moleculeID += currChar;
					currChar = nextCharacter(compositionString, charIndex++);
				}
				if (!knownSymbols.contains(moleculeID)) {
					throw new CompositionParseException("Found unknown symbol in composition: " + moleculeID + ".");
				}
				subComposition = new Composition(new String[] { moleculeID }, new int[] { 1 });
			}

			String countString = "";
			while (currChar != null && (Character.isDigit(currChar) || currChar == '.')) {
				countString += currChar;
				currChar = nextCharacter(compositionString, charIndex++);
			}
			double count = countString.length() == 0 ? 1 : Double.parseDouble(countString);

			this.addSymbols(subComposition.m_CompositionMap, count);
		}
	}

	/**
	 * Returns the character in the given string starting at the given index, or
	 * null if the index is past the end of the string.
	 * 
	 * @param string The given string
	 * @param index  The given index
	 * @return the character in the given string starting at the given index, or
	 *         null if the index is past the end of the string.
	 */
	private Character nextCharacter(String string, int index) {
		return (index < string.length()) ? string.charAt(index) : null;
	}

	/**
	 * For compositions with integer (or near-integer) counts for all symbols,
	 * returns the composition with all counts divided by their greater common
	 * factor.
	 * 
	 * @return the composition with all counts divided by their greater common
	 *         factor for compositions with integer (or near-integer) counts for all
	 *         symbols
	 */
	public Composition getReducedComposition() {

		int[] counts = new int[m_CompositionMap.size()];
		String[] symbols = new String[counts.length];

		int symbolNum = 0;
		for (String symbol : m_CompositionMap.keySet()) {
			symbols[symbolNum] = symbol;
			double count = m_CompositionMap.get(symbol);
			int intCount = (int) Math.round(count);
			if (Math.abs(count - intCount) > 1E-7) {
				return this;
			} // Not reducible
			counts[symbolNum] = intCount;
			symbolNum++;
		}

		int gcf = MSMath.GCF(counts);
		if (gcf == 1) {
			return this;
		}

		counts = MSMath.arrayDivide(counts, gcf);
		return new Composition(symbols, counts);
	}

	/**
	 * Returns the symbols for this composition
	 * 
	 * @return the symbols for this composition
	 */
	public Set<String> getSymbols() {
		return m_CompositionMap.keySet();
	}

	/**
	 * Returns the amount of the given symbol in this composition.
	 * 
	 * @param symbol The given symbol (typically an element symbol or ion type).
	 * @return the amount of the given symbol in this composition.
	 */
	public double getCount(String symbol) {
		if (m_CompositionMap.containsKey(symbol)) {
			return m_CompositionMap.get(symbol);
		}
		return 0;
	}

	/**
	 * Returns a map representation of this composition in which the keys are
	 * symbols and the values are amounts. The map is sorted by symbol.
	 * 
	 * @return a map representation of this composition in which the keys are
	 *         symbols and the values are amounts. The map is sorted by symbol.
	 */
	public Map<String, Double> getCompositionMap() {
		return new TreeMap<String, Double>(m_CompositionMap);
	}

	/**
	 * Return the composition in terms of elements. E.g. all symbols for polyatomic
	 * ion types will be replaced by the symbols for the corresponding elements.
	 * 
	 * @return the composition in terms of elements. E.g. all symbols for polyatomic
	 *         ion types will be replaced by the symbols for the corresponding
	 *         elements.
	 */
	public Composition getElementalComposition() {
		return new Composition(this.getElementalCompositionMap());
	}

	/**
	 * Returns the distinct elements in the compound represented by this
	 * composition.
	 * 
	 * @return the distinct elements in the compound represented by this
	 *         composition.
	 */
	public Element[] getDistinctElements() {
		String[] symbols = this.getElementalCompositionMap().keySet().toArray(new String[0]);
		Element[] returnArray = new Element[symbols.length];
		for (int symbolNum = 0; symbolNum < symbols.length; symbolNum++) {
			returnArray[symbolNum] = Element.getElement(symbols[symbolNum]);
		}
		return returnArray;
	}

	/**
	 * Returns a map in which the keys are the distinct elements in the compound
	 * represented by this composition, and the values are the corresponding
	 * amounts.
	 * 
	 * @return a map in which the keys are the distinct elements in the compound
	 *         represented by this composition, and the values are the corresponding
	 *         amounts.
	 */
	public Map<String, Double> getElementalCompositionMap() {

		// This approach does not rely on the symbol for a polyatomic ion matching its
		// composition, but it does rely on the Ion Factory knowing the structures.

		// Simplified. TODO simplify more by requiring only known ions?
		TreeMap<String, Double> returnMap = new TreeMap<String, Double>();
		for (String symbol : m_CompositionMap.keySet()) {
			Ion ion = IonFactory.get(symbol);
			Map<String, Double> ionCounts = new HashMap<String, Double>();
			if (ion != null) {
				Element[] elements = ion.getIonType().getElements();
				int[] counts = ion.getIonType().getCounts();
				for (int elementNum = 0; elementNum < elements.length; elementNum++) {
					Element element = elements[elementNum];
					ionCounts.put(element.getSymbol(), 1.0 * counts[elementNum]);
				}
			} else {
				Composition composition = new Composition(symbol);
				if (composition.getSymbols().size() == 1) {
					throw new CompositionParseException("Species " + symbol + " is not recognized.");
				}
				ionCounts = composition.getElementalCompositionMap();
			}
			for (String elementSymbol : ionCounts.keySet()) {
				double elementCountInIon = ionCounts.get(elementSymbol);
				double count = returnMap.containsKey(elementSymbol) ? returnMap.get(elementSymbol) : 0;
				count += m_CompositionMap.get(symbol) * elementCountInIon;
				returnMap.put(elementSymbol, count);
			}
		}

		return returnMap;
	}

	/**
	 * Returns this composition in terms of known ion types, including polyatomic ions.
	 * 
	 * @return this composition in terms of known ion types, including polyatomic ions.
	 */
	public Composition getKnownIonComposition() {
		return new Composition(this.getKnownIonCompositionMap());
	}

	/**
	 * Returns a composition in terms of known ion types, including polyatomic ions.
	 * 
	 * @return a composition in terms of known ion types, including polyatomic ions.
	 */
	public Map<String, Double> getKnownIonCompositionMap() {

		// This approach does not rely on the symbol for a polyatomic ion matching its
		// composition, but it does rely on the Ion Factory knowing the structures.
		TreeMap<String, Double> returnMap = new TreeMap<String, Double>();
		for (String symbol : m_CompositionMap.keySet()) {
			Ion ion = IonFactory.get(symbol);
			Map<String, Double> ionCounts = new HashMap<String, Double>();
			if (ion != null) {
				ionCounts.put(symbol, 1.0);
			} else {
				Composition composition = new Composition(symbol);
				if (composition.getSymbols().size() == 1) {
					throw new CompositionParseException("Unknown species: " + symbol);
				}
				ionCounts = composition.getElementalCompositionMap();
			}
			for (String subSymbol : ionCounts.keySet()) {
				double numSitesInIon = ionCounts.get(subSymbol);
				double count = returnMap.containsKey(subSymbol) ? returnMap.get(subSymbol) : 0;
				count += m_CompositionMap.get(symbol) * numSitesInIon;
				returnMap.put(subSymbol, count);
			}
		}

		return returnMap;

	}

	/**
	 * Returns a standardized string representation of this composition.
	 * 
	 * @return a standardized string representation of this composition.
	 */
	public String getStandardizedCompositionString() {
		return getStandardizedCompositionString("");
	}

	/**
	 * Returns a standardized string representation of this composition, using the
	 * given delimiter to separate terms.
	 * 
	 * @param delimiter The delimiter used to separate terms.
	 * @return a standardized string representation of this composition, using the
	 *         given delimiter to separate terms.
	 */
	public String getStandardizedCompositionString(String delimiter) {
		String returnString = "";

		String[] symbols = m_CompositionMap.keySet().toArray(new String[0]);
		for (int symbolNum = 0; symbolNum < symbols.length; symbolNum++) {
			String symbol = symbols[symbolNum];
			if (!Element.isKnownElement(symbol)) {
				symbols[symbolNum] = "~" + symbol + "~";
			}
		}
		Arrays.sort(symbols);

		for (String symbol : symbols) {
			String strippedSymbol = symbol.replaceAll("~", "");
			symbol = symbol.replaceFirst("~", "(");
			symbol = symbol.replaceFirst("~", ")");
			double count = m_CompositionMap.get(strippedSymbol);
			int intCount = (int) Math.round(count);
			returnString += symbol;
			if (Math.abs(count - intCount) < 1E-7) {
				returnString += intCount;
			} else {
				returnString += count;
			}
			returnString += delimiter;
		}

		// Remove the last delimiter
		return returnString.substring(0, returnString.length() - delimiter.length());

	}

	/**
	 * Returns true if the given composition matches this one, within the given
	 * tolerance for amounts.
	 * 
	 * @param composition The given composition.
	 * @param tolerance   The given tolerance for amounts.
	 * @return true if the given composition matches this one, within the given
	 *         tolerance for amounts.
	 */
	public boolean matches(Composition composition, double tolerance) {

		for (String symbol : m_CompositionMap.keySet()) {
			double thisCount = m_CompositionMap.get(symbol);
			double thatCount = composition.getCount(symbol);
			if (Math.abs(thisCount - thatCount) > tolerance) {
				return false;
			}
		}
		return true;

	}

	/**
	 * Returns true if this composition is charge balanced within the given charge
	 * sum tolerance, or false if it is not.
	 * 
	 * @param tolerance The given charge sum tolerance
	 * @return true if this composition is charge balanced within the given charge
	 *         sum tolerance, or false if it is not.
	 */
	public boolean isChargeBalanced(double tolerance) {

		double sum = 0;
		for (String symbol : m_CompositionMap.keySet()) {
			Ion ion = IonFactory.get(symbol);
			double weight = m_CompositionMap.get(symbol);
			sum += weight * ion.getOxidationState();
		}
		return (Math.abs(sum) < tolerance);

	}

}
