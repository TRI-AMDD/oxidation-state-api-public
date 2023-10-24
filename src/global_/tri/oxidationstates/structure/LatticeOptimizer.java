package global_.tri.oxidationstates.structure;

import matsci.engine.IContinuousFunctionState;
import matsci.location.Vector;
import matsci.location.basis.CartesianBasis;
import matsci.structure.BravaisLattice;
import global.tri.structure.mapper.GeneralStructureMapper;

/**
 * This class is used to find an "average" lattice given a set of target
 * lattices. The average lattice is found by minimizing the mean squared skew
 * factor between the target lattices and the other lattices, where the skew
 * factor is defined in the GeneralStructureMapper. The parameters are the
 * coordinate of the lattice vectors.
 * 
 * @author timmueller
 *
 */
public class LatticeOptimizer implements IContinuousFunctionState {

	private BravaisLattice[] m_TargetLattices;
	private int m_NumPeriodicDimensions;
	private double[] m_Parameters;

	/**
	 * Initialize a Lattice Optimzier by copying an old optimizer and setting a new
	 * set of parameters.
	 * 
	 * @param oldOptimizer  The old lattice optimizer
	 * @param newParameters The parameters for the new optimizer
	 */
	private LatticeOptimizer(LatticeOptimizer oldOptimizer, double[] newParameters) {
		m_TargetLattices = oldOptimizer.m_TargetLattices;
		m_NumPeriodicDimensions = oldOptimizer.m_NumPeriodicDimensions;
		m_Parameters = newParameters;
	}

	/**
	 * Initialize the lattice optimizer for a search starting from the
	 * initialLattice. The optimizer will try to find the lattice closest to the
	 * target lattices.
	 * 
	 * @param initialLattice The lattice from which we should start the search.
	 * @param targetLattices We are trying to find the lattice that is the best
	 *                       representation of these lattices.
	 */
	public LatticeOptimizer(BravaisLattice initialLattice, BravaisLattice[] targetLattices) {
		m_TargetLattices = targetLattices.clone();
		m_NumPeriodicDimensions = targetLattices[0].numPeriodicVectors();
		m_Parameters = new double[m_NumPeriodicDimensions * 3];

		// Initial guess of the parameters
		BravaisLattice lattice = initialLattice;
		Vector[] periodicVectors = lattice.getPeriodicVectors();
		int paramNum = 0;
		for (Vector vector : periodicVectors) {
			double[] coordArray = vector.getDirectionArray(CartesianBasis.getInstance());
			for (double value : coordArray) {
				m_Parameters[paramNum++] = value; // Make the parameters unbounded
			}
		}
	}

	@Override
	public int numParameters() {
		return m_Parameters.length;
	}

	@Override
	public double[] getUnboundedParameters(double[] template) {
		return m_Parameters.clone();
	}

	@Override
	public double[] getGradient(double[] template) {
		double[] returnArray = (template == null) ? new double[this.numParameters()] : template;
		double increment = 1E-5;
		double[] parameters = m_Parameters.clone();
		double baseValue = this.getValue();
		for (int paramNum = 0; paramNum < template.length; paramNum++) {
			double param = parameters[paramNum];
			parameters[paramNum] = param + increment;
			double plusValue = this.setUnboundedParameters(parameters).getValue();
			parameters[paramNum] = param - increment;
			double minusValue = this.setUnboundedParameters(parameters).getValue();
			returnArray[paramNum] = (plusValue - minusValue) / (2 * increment);
			parameters[paramNum] = param;
		}

		return returnArray;
	}

	@Override
	public IContinuousFunctionState setUnboundedParameters(double[] parameters) {
		for (double value : parameters) {
			if (Double.isNaN(value)) {
				System.currentTimeMillis();
			}
		}
		return new LatticeOptimizer(this, parameters);
	}

	@Override
	public double getValue() {
		BravaisLattice lattice = this.getLattice();
		double totalSkewFactor = 0;
		for (BravaisLattice knownLattice : m_TargetLattices) {
			double skewFactor = GeneralStructureMapper.getSkewFactor(knownLattice, lattice);
			totalSkewFactor += skewFactor * skewFactor; // It's important to square this to hit the middle; otherwise
														// multiple values are equally likely
		}
		return totalSkewFactor / m_TargetLattices.length;
	}

	/**
	 * Returns the optimized lattice
	 * 
	 * @return the optimized lattice
	 */
	public BravaisLattice getLattice() {

		Vector[] cellVectors = new Vector[m_Parameters.length / 3];
		int paramNum = 0;
		for (int vecNum = 0; vecNum < cellVectors.length; vecNum++) {
			double[] direction = new double[3];
			for (int dimNum = 0; dimNum < direction.length; dimNum++) {
				direction[dimNum] = m_Parameters[paramNum++];
			}
			cellVectors[vecNum] = new Vector(direction, CartesianBasis.getInstance());
		}
		return new BravaisLattice(cellVectors);

	}

}
