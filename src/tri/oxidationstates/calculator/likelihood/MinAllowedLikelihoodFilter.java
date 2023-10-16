package tri.oxidationstates.calculator.likelihood;

import matsci.util.arrays.ArrayIndexer.Filter;
import matsci.util.arrays.ArrayUtils;
import tri.oxidationstates.calculator.likelihood.LikelihoodCalculator.PotentialOptimizer;
import tri.oxidationstates.ion.IonFactory.Ion;

/**
 * This filter only allows combinations of oxidation states that have a
 * likelihood score above a given threshold
 * 
 * @author timmueller
 */
public class MinAllowedLikelihoodFilter implements Filter {

	private LikelihoodCalculator m_Evaluator;
	private Ion[][] m_AllowedIons;
	private double m_MinAllowedLikelihood;
	private double m_LastValidLikelihood = 0;
	private double m_LastValidFermiLevel = Double.NaN;

	/**
	 * Initialize the filter
	 * 
	 * @param calculator           The calculator used to calculate the likelihood
	 *                             score
	 * @param allowedIons          Each element in this array is an array of ions of
	 *                             the same type with different oxidation states.
	 *                             For each state, the set of ions is given by
	 *                             ions[i] = allowedIons[i][state[i]]
	 * @param minAllowedLikelihood This filter will screen out all combinations of
	 *                             ions with likelihood score below this value
	 */
	public MinAllowedLikelihoodFilter(LikelihoodCalculator calculator, Ion[][] allowedIons,
			double minAllowedLikelihood) {
		m_Evaluator = calculator;
		m_AllowedIons = (Ion[][]) ArrayUtils.copyArray(allowedIons);
		m_MinAllowedLikelihood = minAllowedLikelihood;
	}

	/**
	 * Set the minimum allowed likelihood score. Any set of ions with a likelihood
	 * score below this value will be screened out.
	 * 
	 * @param value the minimum allowed likelihood score. Any set of ions with a
	 *              likelihood score below this value will be screened out.
	 */
	public void setMinAllowedLikelihood(double value) {
		m_MinAllowedLikelihood = value;
	}

	/**
	 * Returns the minimum allowed likelihood score. Any set of ions with a
	 * likelihood score below this value will be screened out.
	 * 
	 * @return the minimum allowed likelihood score. Any set of ions with a
	 *         likelihood score below this value will be screened out.
	 */
	public double getMinAllowedLikelihood() {
		return m_MinAllowedLikelihood;
	}

	/**
	 * Returns the likelihood score of the last state to pass this filter. Note that
	 * this state may not have passed other filters.
	 * 
	 * @return the likelihood score of the last state to pass this filter. Note that
	 *         this state may not have passed other filters.
	 */
	public double getLastValidLikelihood() {
		return m_LastValidLikelihood;
	}

	/**
	 * Returns the optimal electronic chemical potential for the last state to pass
	 * this filter. Note that this state may not have passed other filters.
	 * 
	 * @return the optimal electronic chemical potential for the last state to pass
	 *         this filter. Note that this state may not have passed other filters.
	 */
	public double getLastValidFermiLevel() {
		return m_LastValidFermiLevel;
	}

	@Override
	public int getBranchIndex(int[] currentState) {

		PotentialOptimizer optimizer = null;
		
		/**
		 * Could really start from 2 if we rule-out single-element compounds, but this leaves open the 
		 * possibility of mixed-valence single-element compounds.
		 */
		for (int numIons = 1; numIons <= currentState.length; numIons++) { 
			Ion[] ions = new Ion[numIons];
			for (int ionNum = 0; ionNum < ions.length; ionNum++) {
				int stateIndex = currentState.length - ionNum - 1;
				int speciesIndex = currentState[stateIndex];
				ions[ionNum] = m_AllowedIons[stateIndex][speciesIndex];
			}
			optimizer = m_Evaluator.optimizeLikelihood(ions);

			if (optimizer.getMaxLikelihood() < m_MinAllowedLikelihood) {
				return currentState.length - numIons;
			}
		}

		m_LastValidLikelihood = optimizer.getMaxLikelihood();
		m_LastValidFermiLevel = optimizer.getOptimalFermiLevel();
		return -1;
	}
}
