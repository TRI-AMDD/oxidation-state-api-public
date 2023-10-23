package global.tri.oxidationstates.webapi;

/**
 * This class performs linear transformations between the electronic chemical potential and the mapped potential
 * 
 * @author timmueller
 *
 */
public class PotentialMapper {
	
	private double m_Slope = 0.197921181759699;
	private double m_Intercept= 0.681805717542049;
	
	/**
	 * Create a Potential Mapper with the default slope and intercept
	 */
	public PotentialMapper() {}
	
	/**
	 * Create a PotentialMapper with the given slope and intercept
	 * 
	 * @param slope The slope
	 * @param intercept The intercept
	 */
	public PotentialMapper(double slope, double intercept) {
		m_Slope = slope;
		m_Intercept = intercept;
	}
	
	/**
	 * Converts the mapped potential to the electronic chemical potential used internally by the oxidation analyzer.
	 * 
	 * @param value An electronic chemical potential value
	 * @return The corresponding mapped potential value
	 */
	public double toMappedPotential(double value) {
		return value * m_Slope + m_Intercept;		
	}
	
	/**
	 * Converts the mapped potential to the electronic chemical potential used internally by the oxidation analyzer.
	 * 
	 * @param values An array of electronic chemical potential values
	 * @return An array of the corresponding mapped potential value
	 */
	public double[] toMappedPotential(double[] values) {
		double[] returnArray = new double[values.length];
		for (int valNum = 0; valNum < returnArray.length; valNum++) {
			returnArray[valNum] = toMappedPotential(values[valNum]);
		}
		return returnArray;
	}
	
	/**
	 * Converts the mapped potential to the electronic chemical potential used internally by the oxidation analyzer.
	 * 
	 * @param mappedPotentials An array of mapped potential values
	 * @return An array of the corresponding electronic chemical potential values
	 */
	public double[] fromMappedPotential(double[] mappedPotentials) {
		double[] returnArray = new double[mappedPotentials.length];
		for (int valNum = 0; valNum < returnArray.length; valNum++) {
			returnArray[valNum] = fromMappedPotential(mappedPotentials[valNum]);
		}
		return returnArray;
	}
	
	/**
	 * Converts the mapped potential to the electronic chemical potential used internally by the oxidation analyzer.
	 * 
	 * @param mappedPotential A mapped potential value
	 * @return The corresponding electronic chemical potential value
	 */
	public double fromMappedPotential(double mappedPotential) {
		return (mappedPotential - m_Intercept) / m_Slope;
	}
	
	/**
	 * Returns the "slope" value in the equation mapped_potential = electronic_chemical_potential * slope + intercept,
	 * where the electronic chemical potential is the potential used internally by the oxidation analyzer.
	 * 
	 * @return the "slope" value in the equation mapped_potential = electronic_chemical_potential * slope + intercept,
	 * where the electronic chemical potential is the potential used internally by the oxidation analyzer.
	 * 
	 */
	public double getSlope() {
		return m_Slope;
	}
	
	/**
	 * Returns the "intercept" value in the equation mapped_potential = electronic_chemical_potential * slope + intercept,
	 * where the electronic chemical potential is the potential used internally by the oxidation analyzer.
	 * 
	 * @return the "intercept" value in the equation mapped_potential = electronic_chemical_potential * slope + intercept,
	 * where the electronic chemical potential is the potential used internally by the oxidation analyzer.
	 * 
	 */
	public double getIntercept() {
		return m_Intercept;
	}

}
