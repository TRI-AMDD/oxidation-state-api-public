package _global.tri.oxidationstates.fitting;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.Callable;

import _global.tri.oxidationstates.calculator.likelihood.LikelihoodCalculator;
import _global.tri.oxidationstates.calculator.likelihood.LikelihoodCalculator.PotentialOptimizer;
import _global.tri.oxidationstates.fitting.OxidationStateData.Entry;

/**
 * This class is used for multi-threaded model fitting
 * 
 * @author timmueller
 *
 */
public class CallableCalculator implements Callable<Double> {

	private LikelihoodCalculator m_Calculator;
	private Entry[] m_Entries;

	/**
	 * Create a CallableCalculator that will work on the given subset of entries.
	 * Each thread should have its own subset.
	 * 
	 * @param entries The given subset of entries.
	 */
	public CallableCalculator(Entry[] entries) {
		m_Entries = entries.clone();
	}

	/**
	 * Set the calculator to be used to calculate the likelihood score.
	 * 
	 * @param calculator the likelihood calculator to be used to calculate the
	 *                   likelihood score.
	 */
	public void setCalculator(LikelihoodCalculator calculator) {
		m_Calculator = calculator;
	}

	@Override
	public Double call() throws Exception {
		double returnValue = 0;
		for (Entry entry : m_Entries) {
			returnValue += Math.log(m_Calculator.optimizeLikelihood(entry).getMaxLikelihood());
		}
		return returnValue;
	}

	/**
	 * Useful for debugging, appends information about the current fit to a file
	 * with the given name
	 * 
	 * @param fileName The given filename
	 */
	public void appendFile(String fileName) {

		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(fileName, true));

			for (Entry entry : m_Entries) {
				PotentialOptimizer optimizer = m_Calculator.optimizeLikelihood(entry);
				double potential = optimizer.getOptimalFermiLevel();
				double likelihood = optimizer.getMaxLikelihood();
				writer.write(entry.getID() + "\t" + entry.getGivenCompositionString() + "\t" + likelihood + "\t"
						+ potential + "\n");
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}
}
