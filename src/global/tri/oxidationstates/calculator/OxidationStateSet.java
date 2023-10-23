package global.tri.oxidationstates.calculator;

import java.util.HashSet;
import java.util.Map;

import global.tri.oxidationstates.ion.IonFactory;
import global.tri.oxidationstates.ion.IonFactory.Ion;
import global.tri.oxidationstates.util.Composition;
import matsci.util.arrays.ArrayUtils;

/**
 * Contains a set of oxidation states assigned to ion types (i.e. a set of ions)
 * 
 * @author timmueller
 *
 */
public class OxidationStateSet {

	// The ions in this set
	private Ion[] m_Ions;

	// The amounts of those ions, in the same order as m_Ions
	private double[] m_Weights;

	/**
	 * Initialize the OxidationStateSet from a given composition containing the ions
	 * and weights
	 * 
	 * @param composition The given composition
	 */
	public OxidationStateSet(Composition composition) {
		Map<String, Double> map = composition.getCompositionMap();
		m_Ions = new Ion[map.size()];
		m_Weights = new double[map.size()];
		int ionNum = 0;
		for (String symbol : map.keySet()) {
			m_Ions[ionNum] = IonFactory.get(symbol);
			m_Weights[ionNum] = map.get(symbol);
			ionNum++;
		}
	}

	/**
	 * Create an OxidationStateSet from the given ions and weights
	 * 
	 * @param ions    The given ions
	 * @param weights The amount of each of those ions, in the same order as the
	 *                ions array.
	 */
	public OxidationStateSet(Ion[] ions, double[] weights) {
		m_Ions = ions.clone();
		m_Weights = weights.clone();
	}

	/**
	 * Returns the set of ions
	 * 
	 * @return the set of ions
	 */
	public Ion[] getIons() {
		return m_Ions.clone();
	}

	/**
	 * Returns the amounts of each of the ions, in the same order as the array
	 * returned by getIons()
	 * 
	 * @return the amounts of each of the ions, in the same order as the array
	 *         returned by {@link getIons()}
	 */
	public double[] getWeights() {
		return m_Weights.clone();
	}

	/**
	 * Return the amount for the given ion
	 * 
	 * @param ion THe given ion
	 * @return the amount for the given ion
	 */
	public double getWeight(Ion ion) {
		int ionIndex = ArrayUtils.findIndex(m_Ions, ion);
		return m_Weights[ionIndex];
	}

	/**
	 * Returns true if the same ion type has more than one oxidation state in this
	 * set, false otherwise.
	 * 
	 * @return true if the same ion type has more than one oxidation state in this
	 *         set, false otherwise.
	 */
	public boolean isMixedValence() {

		HashSet<String> seenSymbols = new HashSet();
		for (int ionNum = 0; ionNum < m_Ions.length; ionNum++) {
			if (m_Weights[ionNum] < 1E-6) {
				continue;
			}
			String ionTypeSymbol = m_Ions[ionNum].getIonType().getSymbol();
			if (seenSymbols.contains(ionTypeSymbol)) {
				return true;
			}
			seenSymbols.add(ionTypeSymbol);
		}
		return false;
	}

	@Override
	public String toString() {
		String returnString = "";
		String[] symbols = new String[m_Ions.length];
		for (int ionNum = 0; ionNum < m_Ions.length; ionNum++) {
			symbols[ionNum] = m_Ions[ionNum].getSymbol();
		}
		int[] map = ArrayUtils.getSortPermutation(symbols);

		for (int ionNum = 0; ionNum < m_Ions.length; ionNum++) {
			Ion ion = m_Ions[map[ionNum]];
			double weight = m_Weights[map[ionNum]];
			returnString += "[" + ion.getSymbol() + "]" + weight + " ";
		}
		return returnString.trim();
	}
}