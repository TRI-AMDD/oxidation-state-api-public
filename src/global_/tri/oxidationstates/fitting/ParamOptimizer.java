package global_.tri.oxidationstates.fitting;

import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import global_.tri.oxidationstates.calculator.likelihood.LikelihoodCalculator;
import global_.tri.oxidationstates.fitting.OxidationStateData.Entry;
import matsci.engine.IContinuousFunctionState;

/**
 * This class is used for optimizing the model parameters
 * 
 * @author timmueller
 *
 */
public class ParamOptimizer {

	private OxidationStateData m_Data;
	private double m_RegularizationParameter = 1E-5;

	// The below fields are for multithreading
	private int m_NumThreads;
	private ExecutorService m_ExecutorService;
	private CallableCalculator[] m_Calculators; // These don't actually save the state. Just create them once for speed.

	/**
	 * Construct a parameter optimizer for a given set of training data
	 * 
	 * @param data The set of training data
	 */
	public ParamOptimizer(OxidationStateData data) {
		this(data, 1);
	}

	/**
	 * Construct a parameter optimizer for a given set of training data and number
	 * of threads
	 * 
	 * @param data       The set of training data
	 * @param numThreads The number of threads to run on (shared memory).
	 */
	public ParamOptimizer(OxidationStateData data, int numThreads) {
		m_Data = data;
		m_NumThreads = numThreads;
		m_ExecutorService = Executors.newFixedThreadPool(m_NumThreads);
		Entry[][] batches = batchEntries(m_NumThreads);
		m_Calculators = new CallableCalculator[m_NumThreads];
		for (int threadNum = 0; threadNum < m_Calculators.length; threadNum++) {
			m_Calculators[threadNum] = new CallableCalculator(batches[threadNum]);
		}

	}

	/**
	 * Divide the entries evenly into batches, each to be handled on its own thread
	 * 
	 * @param numBatches The number of batches (should match the number of threads)
	 * @return The entries in the training set divided evenly into batches
	 */
	public Entry[][] batchEntries(int numBatches) {

		Entry[][] returnArray = new Entry[numBatches][];

		int entryNum = 0;
		for (int batchNum = 0; batchNum < returnArray.length; batchNum++) {
			int maxEntryNum = (int) Math.round(m_Data.numEntries() * 1.0 * (batchNum + 1) / numBatches);
			int numEntries = maxEntryNum - entryNum;
			returnArray[batchNum] = new Entry[numEntries];
			while (entryNum < maxEntryNum) {
				int returnIndex = numEntries - (maxEntryNum - entryNum);
				returnArray[batchNum][returnIndex] = m_Data.getEntry(entryNum++);
			}
		}

		return returnArray;

	}

	/**
	 * Sets the regularization parameter
	 * 
	 * @param value The regularization parameter
	 */
	public void setRegularizationParameter(double value) {
		m_RegularizationParameter = value;
	}

	/**
	 * Returns the regularization parameter
	 * 
	 * @return the regularization parameter
	 */
	public double getRegularizationParamaeter() {
		return m_RegularizationParameter;
	}

	/**
	 * Returns the training data
	 * 
	 * @return the training data
	 */
	public OxidationStateData getData() {
		return m_Data;
	}

	/**
	 * Returns a state representing the parameters contained in the provided
	 * calculator.
	 * 
	 * @param calculator The calculator containing the parameters and used to
	 *                   evaluate them.
	 * @return a state representing a current snapshot of the parameters
	 */
	public ParameterState getParameterState(LikelihoodCalculator calculator) {
		return new ParameterState(calculator);
	}

	/**
	 * Returns a state representing a set of parameters read from a file
	 * 
	 * @param paramFileName The name of the file to read
	 * @return a state representing a set of parameters read from a file
	 */
	public ParameterState getParameterState(String paramFileName) {
		LikelihoodCalculator calculator = new LikelihoodCalculator(paramFileName);
		return this.getParameterState(calculator);
	}

	/**
	 * This represents a snapshot of the parameters at one point in the optimization
	 * algorithm
	 * 
	 * @author timmueller
	 *
	 */
	public class ParameterState implements IContinuousFunctionState {

		private LikelihoodCalculator m_Calculator;

		/**
		 * Initialize the state with a calculator used to evaluate the parameters. The
		 * parameters are actually stored in the calculator.
		 * 
		 * @param calculator a calculator used to evaluate the parameters. The
		 *                   parameters are actually stored in the calculator.
		 */
		private ParameterState(LikelihoodCalculator calculator) {
			this.setCalculator(calculator);
		}

		/**
		 * Creates a new parameter state by updating the parameters of an old state
		 * 
		 * @param oldState      The state to be updated
		 * @param newParameters The new parameters for the state
		 */
		private ParameterState(ParameterState oldState, double[] newParameters) {
			this(oldState.m_Calculator.setParameters(newParameters));
		}

		@Override
		public int numParameters() {
			return m_Calculator.numParameters();
		}

		@Override
		public double[] getUnboundedParameters(double[] template) {

			return m_Calculator.getParameters(template);
		}

		/**
		 * Sets the calculator (and parameters) to the provided calculator. Used
		 * internally for calculating gradients with finite differences as well as
		 * constructing a parallel array of calculators.
		 * 
		 * @param calculator
		 */
		private void setCalculator(LikelihoodCalculator calculator) {

			m_Calculator = calculator;
			for (CallableCalculator callableCalculator : m_Calculators) {
				callableCalculator.setCalculator(calculator);
			}
		}

		@Override
		public double[] getGradient(double[] template) {

			double[] parameters = m_Calculator.getParameters(null);
			double[] returnArray = (template != null) ? template : new double[parameters.length];

			LikelihoodCalculator origCalculator = m_Calculator;

			// Use finite differences
			double delta = 1E-5;
			for (int paramNum = 0; paramNum < returnArray.length; paramNum++) {

				double oldParam = parameters[paramNum];

				parameters[paramNum] = oldParam + delta;
				this.setCalculator(m_Calculator.setParameters(parameters));
				double plusValue = this.getValue();

				parameters[paramNum] = oldParam - delta;
				this.setCalculator(m_Calculator.setParameters(parameters));
				double minusValue = this.getValue();

				parameters[paramNum] = oldParam;
				this.setCalculator(origCalculator);

				returnArray[paramNum] = (plusValue - minusValue) / (2 * delta);
			}

			return returnArray;
		}

		@Override
		public IContinuousFunctionState setUnboundedParameters(double[] parameters) {
			return new ParameterState(this, parameters);
		}

		/**
		 * Write the current parameters to a text file
		 * 
		 * @param fileName The name of the file to write
		 */
		public void writeFile(String fileName) {

			for (CallableCalculator calc : m_Calculators) {
				calc.appendFile(fileName);
			}
		}

		@Override
		public double getValue() {

			double score = 0;
			ArrayList<Future<Double>> futures = new ArrayList<Future<Double>>();

			// Submit the parallel tasks
			for (CallableCalculator calculator : m_Calculators) {
				futures.add(m_ExecutorService.submit(calculator));
			}

			// Collect the results
			for (Future<Double> result : futures) {
				try {
					score += result.get();
				} catch (Exception e) {
					throw new RuntimeException(e);
				}
			}

			double regularizer = this.getRegularizer();
			return -score / m_Data.numEntries() + regularizer;
		}

		/**
		 * Gets the regularization term, which is the regularization parameter times the
		 * sum of spreads
		 * 
		 * @return the regularization term, which is the regularization parameter times
		 *         the sum of spreads
		 */
		public double getRegularizer() {
			return m_Calculator.getSumOfSpreads() * m_RegularizationParameter;
		}

		/**
		 * Returns the likelihood calculator containing the current parameters
		 * 
		 * @return the likelihood calculator containing the current parameters
		 */
		public LikelihoodCalculator getCalculator() {
			return m_Calculator;
		}

	}
}
