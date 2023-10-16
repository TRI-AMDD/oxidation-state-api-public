package tri.oxidationstates.calculator.likelihood;

import tri.oxidationstates.calculator.OxidationStateSet;
import tri.oxidationstates.ion.IonFactory.Ion;

/**
 * An OxidationStateSet that keeps track of properties specific to the
 * likelihood score, such as the likelihood score and electronic chemical
 * potential.
 * 
 * @author timmueller
 *
 */
public class LikelihoodStateSet extends OxidationStateSet implements Comparable<LikelihoodStateSet> {

	private double m_MaxLikelihood;
	private double m_FermiLevel;

	/**
	 * Initialize the LikelihoodStateSet
	 * 
	 * @param ions              The ions in this set
	 * @param weights           The composition of each ion, in the same order as
	 *                          the ions array
	 * @param maxLikelihood     The likelihood score at the optimal electronic
	 *                          chemical potential
	 * @param optimalFermiLevel The optimal electronic chemical potential (i.e. the
	 *                          one that maximizes the likelihood score)
	 */
	public LikelihoodStateSet(Ion[] ions, double[] weights, double maxLikelihood, double optimalFermiLevel) {
		super(ions, weights);
		m_MaxLikelihood = maxLikelihood;
		m_FermiLevel = optimalFermiLevel;
	}

	/**
	 * Returns the likelihood score at the optimal electronic chemical potential
	 * 
	 * @return the likelihood score at the optimal electronic chemical potential
	 */
	public double getMaxLikelihood() {
		return m_MaxLikelihood;
	}

	/**
	 * Returns the optimal electronic chemical potential (i.e. the one that
	 * maximizes the likelihood score)
	 * 
	 * @return the optimal electronic chemical potential (i.e. the one that
	 *         maximizes the likelihood score)
	 */
	public double getOptimalFermiLevel() {
		return m_FermiLevel;
	}

	@Override
	public int compareTo(LikelihoodStateSet o) {
		return Double.compare(o.getMaxLikelihood(), m_MaxLikelihood);
	}
}
