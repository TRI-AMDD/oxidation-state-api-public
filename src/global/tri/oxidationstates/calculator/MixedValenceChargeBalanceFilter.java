package global.tri.oxidationstates.calculator;

import global.tri.oxidationstates.ion.IonFactory.Ion;
import matsci.util.arrays.ArrayIndexer.Filter;
import matsci.util.arrays.ArrayUtils;

/**
 * Filters out all combinations of ions that aren't charge balanced, or could
 * not be made charge balanced by allowing one ion to have mixed valence. The
 * mixed-valence ion would have one oxidation state corresponding to the given
 * set of states and one oxidation state that is more negative.
 * 
 * @author timmueller
 *
 */
public class MixedValenceChargeBalanceFilter implements Filter {

	private Ion[][] m_AllowedIons;
	private double[] m_Weights;

	// If the total charge (or possible charge, allowing for mixed valence),
	// starting from the end, is outside these bounds then it cannot be charge
	// balanced
	private double[] m_MinAllowedCharges;
	private double[] m_MaxAllowedCharges;

	private double m_NonSigmaCharge;

	/**
	 * Initializes the filter. It's important that allAllowedIons are sorted from
	 * most negative to least negative oxidation state.
	 * 
	 * @param allAllowedIons Each element of this array is the set of all allowed
	 *                       ions for a particular ion type. IMPORTANT: These arrays
	 *                       should be sorted from most negative to least negative
	 *                       oxidation state.
	 * @param weights        The relative compositions of each ion type, in order
	 *                       corresponding to the order of allAllowedIons.
	 */
	public MixedValenceChargeBalanceFilter(Ion[][] allAllowedIons, double[] weights) {

		m_AllowedIons = (Ion[][]) ArrayUtils.copyArray(allAllowedIons);
		m_Weights = weights.clone();

		m_MinAllowedCharges = new double[allAllowedIons.length];
		m_MaxAllowedCharges = new double[allAllowedIons.length];

		double tolerance = 1E-2; // TOOD consider making this flexible

		double maxAllowedCharge = tolerance;
		double minAllowedCharge = -tolerance;
		for (int moleculeNum = 0; moleculeNum < m_MinAllowedCharges.length; moleculeNum++) {
			m_MinAllowedCharges[moleculeNum] = minAllowedCharge;
			m_MaxAllowedCharges[moleculeNum] = maxAllowedCharge;
			double minForSite = Double.POSITIVE_INFINITY;
			double maxForSite = Double.NEGATIVE_INFINITY;
			Ion[] allowedIons = m_AllowedIons[moleculeNum];
			for (int specNum = 0; specNum < allowedIons.length; specNum++) {
				double oxidationState = allowedIons[specNum].getOxidationState();
				minForSite = Math.min(minForSite, oxidationState * weights[moleculeNum]);
				maxForSite = Math.max(maxForSite, oxidationState * weights[moleculeNum]);
			}
			maxAllowedCharge -= minForSite;
			minAllowedCharge -= maxForSite;
		}
	}

	@Override
	public int getBranchIndex(int[] currentState) {

		double totalCharge = m_NonSigmaCharge;
		double minDeltaCharge = 0; // To account for mixed valence
		for (int ionTypeIndex = currentState.length - 1; ionTypeIndex >= 0; ionTypeIndex--) {
			int stateIndex = currentState[ionTypeIndex];
			double charge = m_AllowedIons[ionTypeIndex][stateIndex].getOxidationState();
			totalCharge += charge * m_Weights[ionTypeIndex];

			/**
			 * This figures out how much the total charge could shift if we allowed this ion
			 * type to have mixed valence, where the additional oxidation state is more
			 * negative than this one.
			 */
			if (stateIndex > 0) {
				double prevCharge = m_AllowedIons[ionTypeIndex][0].getOxidationState();
				double deltaCharge = (prevCharge - charge) * m_Weights[ionTypeIndex];
				minDeltaCharge = Math.min(minDeltaCharge, deltaCharge);
			}

			// Check to see if the total charge is within bounds, accounting for possible
			// mixed valence.
			if (totalCharge < m_MinAllowedCharges[ionTypeIndex]) {
				return ionTypeIndex;
			}
			if (totalCharge + minDeltaCharge > m_MaxAllowedCharges[ionTypeIndex]) {
				return ionTypeIndex;
			}
		}
		return -1;
	}

}
