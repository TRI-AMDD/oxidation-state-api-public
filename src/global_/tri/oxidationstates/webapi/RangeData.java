package global_.tri.oxidationstates.webapi;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;

import global_.tri.oxidationstates.calculator.likelihood.LikelihoodCalculator;

/**
 * This class contains the data needed to generate the oxidation state range plots
 * 
 * @author timmueller
 *
 */
public class RangeData {
	
	private PotentialMapper m_PotentialMapper = new PotentialMapper();
	
	private String m_IonTypeSymbol;
	private double[] m_Boundaries = new double[0];
	private int[] m_OxidationStates = new int[0];
	
	/**
	 * Creates the data necessary to create the plot for a given ion type.
	 * 
	 * @param ionTypeSymbol  The symbol for the ion type (e.g. "Li", "SO4").  For polyatomic ions, the symbol may be 
	 * enclosed in parenthesis, but it doesn't have to be.
	 * @param calculator  The LikelihoodCalculator used to calculate the oxidation states.  The range boundary information
	 * will be read from this calculator.
	 */
	public RangeData(String ionTypeSymbol, LikelihoodCalculator calculator) {
		
		m_IonTypeSymbol = ionTypeSymbol.replace("(", "").replace(")", "");
		int ionTypeIndex = calculator.getIonTypeIndex(m_IonTypeSymbol);
		m_Boundaries = calculator.getBoundaries(ionTypeIndex);
		m_OxidationStates = new int[m_Boundaries.length - 1];
		for (int stateNum = 0; stateNum < m_OxidationStates.length; stateNum++) {
			m_OxidationStates[stateNum] = calculator.getOxidationState(ionTypeIndex, stateNum);
		}
		
	}
	
	/**
	 * This constructor is used by the Jackson library to deserialize a JSON file containing the boundaries.
	 * 
	 * @param ionTypeSymbol  The symbol for the ion type (e.g. "Li", "SO4").  For polyatomic ions, the symbol may be 
	 * enclosed in parenthesis, but it doesn't have to be.
	 * @param oxidationStates The allowed oxidation states, in ascending order
	 * @param rangeBoundaries The boundaries between the oxidation states in units of the mapped potential, in the same order as the oxidation states.  There should be one more boundary than oxidation state.
	 */
	@JsonCreator
	public RangeData(@JsonProperty("ionTypeSymbol")String ionTypeSymbol, @JsonProperty("oxidationStates")int[] oxidationStates, @JsonProperty("rangeBoundaries")double[] rangeBoundaries) {
		m_IonTypeSymbol = ionTypeSymbol.replace("(", "").replace(")", "");
		m_OxidationStates = oxidationStates.clone(); 
		m_Boundaries = new double[rangeBoundaries.length];
		for (int boundaryNum = 0; boundaryNum < m_Boundaries.length; boundaryNum++) {
			m_Boundaries[boundaryNum] = m_PotentialMapper.fromMappedPotential(rangeBoundaries[boundaryNum]);
		}
	}
	
	/**
	 * Returns the given symbol for the ion type (e.g. "Li", "SO4") for which the plot will be generated.
	 * 
	 * @return the given symbol for the ion type (e.g. "Li", "SO4") for which the plot will be generated.
	 */
	public String getIonTypeSymbol() {
		return m_IonTypeSymbol;
	}
	
	/**
	 * Returns the boundaries of the oxidation states in terms of the mapped potential, starting from the left-most and ending with the right-most
	 * 
	 * @return  the boundaries of the oxidation states in terms of the mapped potential, starting from the left-most and ending with the right-most
	 */
	public double[] getRangeBoundaries() {
		double[] mappedBoundaries = new double[m_Boundaries.length];
		for (int boundaryNum = 0; boundaryNum < m_Boundaries.length; boundaryNum++) {
			mappedBoundaries[boundaryNum] = m_PotentialMapper.toMappedPotential(m_Boundaries[boundaryNum]);
		}
		return mappedBoundaries;
	}
	
	/**
	 * Returns the oxidation states to be plotted.  If there are N oxidation states, then there 
	 * will be N-1 oxidation states.  The nth oxidation state is bounded by the nth boundary on 
	 * the left, and the (n+1)th boundary on the right.
	 * 
	 * @return the oxidation states to be plotted.  If there are N oxidation states, then there 
	 * will be N-1 oxidation states.  The nth oxidation state is bounded by the nth boundary on 
	 * the left, and the (n+1)th boundary on the right.
	 */
	public int[] getOxidationStates() {
		return m_OxidationStates.clone();
	}
	
	@Override
	public String toString() {
		String returnString = m_IonTypeSymbol;
		returnString += " " + m_Boundaries[0];
		for (int stateNum = 0; stateNum < m_OxidationStates.length; stateNum++) {
			returnString += " | " + m_OxidationStates[stateNum];
			returnString += " | " + m_Boundaries[stateNum + 1];			
		}
		return returnString;
	}
}
