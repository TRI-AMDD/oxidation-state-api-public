package global.tri.oxidationstates.calculator.frequency;

import matsci.util.arrays.ArrayIndexer.Filter;
import global.tri.oxidationstates.ion.IonFactory.Ion;
import matsci.util.arrays.ArrayUtils;

/**
 * This is used to screen out families of oxidation state assignments that can't
 * possibly be better than the best one found so far. For example, if the best
 * frequency score so far is 0.9, then all set of ions containing an ion with
 * frequency less than 0.9 should be filtered out.
 * 
 * @author timmueller
 *
 */
public class MaxFrequencyFilter implements Filter {

	private FrequencyCalculator m_Evaluator;
	private Ion[][] m_AllowedIons;
	private double m_MaxKnownFrequencyScore = 0;
	private double m_LastValidFrequencyScore = 0;

	/**
	 * Initialize the filter.
	 * 
	 * @param calculator  The calculator that will be used to evaluate the
	 *                    frequencies of ions
	 * @param allowedIons The allowed ions in this search. The indexing of this
	 *                    array should correspond to the currentState array, so that
	 *                    allowedIons[i][currenState[i]] is the ion that would be
	 *                    assigned in that state.
	 */
	public MaxFrequencyFilter(FrequencyCalculator calculator, Ion[][] allowedIons) {
		m_Evaluator = calculator;
		m_AllowedIons = (Ion[][]) ArrayUtils.copyArray(allowedIons);
	}

	/**
	 * The maximum valid frequency score found so far. We set this externally in
	 * case a state that passes this filter is screened out by some other filter.
	 * Only frequency scores for states that pass all filters should be set here.
	 * 
	 * @param value The maximum frequency score for a valid state encountered so
	 *              far.
	 */
	public void setMaxKnownFrequencyScore(double value) {
		m_MaxKnownFrequencyScore = value;
	}

	/**
	 * Returns the maximum frequency score for a valid state encountered so far, as
	 * set by the user.
	 * 
	 * @return the maximum frequency score for a valid state encountered so far, as
	 *         set by the user.
	 */
	public double getMaxKnownFrequencyScore() {
		return m_MaxKnownFrequencyScore;
	}

	/**
	 * Get the last frequency score that passed this filter.
	 * 
	 * @return the last frequency score that passed this filter.
	 */
	public double getLastValidFrequencyScore() {
		return m_LastValidFrequencyScore;
	}

	@Override
	public int getBranchIndex(int[] currentState) {

		double frequencyScore = 1;
		for (int stateIndex = currentState.length - 1; stateIndex >= 0; stateIndex--) {
			int speciesIndex = currentState[stateIndex];
			Ion ion = m_AllowedIons[stateIndex][speciesIndex];
			frequencyScore *= m_Evaluator.getFrequency(ion);

			if (frequencyScore < m_MaxKnownFrequencyScore) { // Since all frequencies are <= 1, we can't possibly beat
																// the best score by multiplying more frequencies
				return stateIndex;
			}
		}

		m_LastValidFrequencyScore = frequencyScore;
		return -1;
	}
}