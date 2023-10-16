package tri.oxidationstates.webapi;

import java.text.DecimalFormat;
import java.util.HashSet;

import matsci.util.arrays.ArrayUtils;
import tri.oxidationstates.calculator.OxidationStateSet;
import tri.oxidationstates.calculator.likelihood.LikelihoodCalculator;
import tri.oxidationstates.calculator.likelihood.LikelihoodStateSet;
import tri.oxidationstates.ion.IonFactory;
import tri.oxidationstates.ion.IonFactory.Ion;

/**
 * This class contains the data in a single row of the table on the web app
 * 
 * @author timmueller
 *
 */
public class TableRow implements Comparable<TableRow> {
	
	private PotentialMapper m_PotentialMapper = new PotentialMapper();
	
	private double[] m_Counts;
	private String[] m_Symbols;
	private int[] m_OxidationStates;
	private double m_Optimalikelihood;
	private double m_OptimalChemicalPotential;
	private double m_GlobalInstabilityIndex;
	private String m_CIFString;
	private double[][] m_BoundaryPairs;
	
	/**
	 * Constructs a row based on an OxidationStateSet (calculated by the oxidation state tool
	 * we've developed) and the global instability index (gii).
	 * 
	 * @param calculator The Likelihood calculator used to generate the oxidation states.  The oxidation state range boundaries will be read from this calculator.
	 * @param stateSet A class containing information about a particular set of oxidation states.
	 * @param gii The global instability index.
	 */
	public TableRow(LikelihoodCalculator calculator, LikelihoodStateSet stateSet, double gii) {
		Ion[] ions = stateSet.getIons();
		
		String[] symbols = new String[ions.length];
		double[][] boundaryPairs = new double[ions.length][2];
		int[] oxidationStates = new int[ions.length];
		for (int ionNum = 0; ionNum < ions.length; ionNum++) {
			Ion ion = ions[ionNum];
			symbols[ionNum] = ion.getIonType().getSymbol();
			if (!ion.getIonType().isElement()) { // ion.getStructure().numDefiningSites() > 1) {
				symbols[ionNum] = "~" + symbols[ionNum] + "~"; // We do this so polyatomic ions end up at the end of the list
			}
			oxidationStates[ionNum] = (int) Math.round(ion.getOxidationState()); // TODO warning if it's not an int
			boundaryPairs[ionNum] = calculator.getBoundaries(ion);
		}
		double[] counts = stateSet.getWeights();
		
		int[] map = ArrayUtils.getSortPermutation(symbols);
		m_Symbols = new String[map.length];
		m_Counts = new double[map.length];
		m_OxidationStates = new int[map.length];
		m_BoundaryPairs = new double[map.length][];
		for (int ionNum = 0; ionNum < map.length; ionNum++) {
			m_Symbols[ionNum] = symbols[map[ionNum]];
			m_Symbols[ionNum] = m_Symbols[ionNum].replaceFirst("~", "("); // We do this so polyatomic ions end up at the end of the list
			m_Symbols[ionNum] = m_Symbols[ionNum].replaceFirst("~", ")");
			m_Counts[ionNum] = counts[map[ionNum]];
			m_OxidationStates[ionNum] = oxidationStates[map[ionNum]];
			m_BoundaryPairs[ionNum] = boundaryPairs[map[ionNum]];
		}
		
		m_Optimalikelihood = stateSet.getMaxLikelihood();
		m_OptimalChemicalPotential = stateSet.getOptimalFermiLevel();
		m_GlobalInstabilityIndex = gii;
	}

	/**
	 * Returns an array of values indicating the count of each ion in a formula unit for the 
	 * given material.  Each element in the array corresponds to an ion for which we have
	 * calculated oxidation states.  The indices in this array match the indices in the arrays returned by the
	 * {@link #getOxidationStates()} and {@link #getSymbols()} methods.
	 * 
	 * @return  An array of values indicating the count of each ion in a formula unit for the 
	 * given material.  Each element in the array corresponds to an ion for which we have
	 * calculated oxidation states.  The indices in this array match the indices in the arrays returned by the
	 * {@link #getOxidationStates()} and {@link #getSymbols()} methods.
	 */
	public double[] getCounts() {
		return m_Counts.clone();
	}
	
	/**
	 * Returns the symbols for each ion in the material.  Oxidation states are not included in 
	 * these symbols.  The indices in this array match the indices in the arrays returned by the
	 * {@link #getCounts()} and {@link #getOxidationStates()} methods.
	 * 
	 * @return the symbols for each ion in the material.  Oxidation states are not included in 
	 * these symbols.  The indices in this array match the indices in the arrays returned by the
	 * {@link #getCounts()} and {@link #getOxidationStates()} methods.
	 */
	public String[] getSymbols() {
		return m_Symbols.clone();
	}
	
	/**
	 * Returns the oxidation states for each ion in the material.  The indices in this array match the indices in the arrays returned by the
	 * {@link #getCounts()} and {@link #getSymbols()} methods.
	 * 
	 * @return the oxidation states for each ion in the material.  The indices in this array match the indices in the arrays returned by the
	 * {@link #getCounts()} and {@link #getSymbols()} methods.
	 */
	public int[] getOxidationStates() {
		return m_OxidationStates.clone();
	}

	/**
	 * Returns the likelihood at the optimal electronic chemical potential for the 
	 * given set of oxidation states.
	 * 
	 * @return the likelihood at the optimal electronic chemcial potential for the 
	 * given set of oxidation states.
	 */
	public double getOptimalLikelihood() {
		return m_Optimalikelihood;
	}
	
	/**
	 * Returns the value of the mapped potential that maximizes the likelihood
	 * of the set of oxidation states coexisting.
	 * 
	 * @return  the value of the mapped potential that maximizes the likelihood
	 * of the set of oxidation states coexisting.
	 */
	public double getOptimalMappedPotential() {
		return m_PotentialMapper.toMappedPotential(m_OptimalChemicalPotential);
	}
	
	/**
	 * Returns the global instability index (gii) for the given set of oxidation states applied to 
	 * a particular structure.
	 * 
	 * @return the global instability index (gii) for the given set of oxidation states applied to 
	 * a particular structure.
	 */
	public double getGlobalInstabilityIndex() {
		return m_GlobalInstabilityIndex;
	}
	
	/**
	 * Returns the minimum and maximum boundaries for each ion in this row.  
	 * 
	 * @return the minimum and maximum boundaries for each ion in this row.  The pairs are arranged as an array pairs of boundary values,
	 * in the same order as the list of ions returned by {@link #getSymbols()}.  The first element in the pair is the lower bound, and the second element
	 * in the pair is the upper bound.
	 */
	public double[][] getBoundaryPairs() {
		return m_BoundaryPairs.clone();
	}
	
	/**
	 * Sets an optional string representation of a CIF file with all oxidation states assigned.
	 * 
	 * @param cifString A string containing contents of the CIF file.
	 */
	protected void setCIFString(String cifString) {
		m_CIFString = cifString;
	}
	
	/**
	 * Returns a string representation of a CIF file with all oxidation states assigned, or null if no such file could be generated.
	 * 
	 * @return a string representation of a CIF file with all oxidation states assigned, or null if no such file could be generated.
	 */
	public String getCIFString() {
		return m_CIFString;
	}
	
	/**
	 * Returns true if the generated oxidation states are "mixed valence", meaning the same element
	 * (or polyatomic ion type) is associated with more than one oxidation state.
	 * 
	 * @return true if the generated oxidation states are "mixed valence", meaning the same element
	 * (or polyatomic ion type) is associated with more than one oxidation state.
	 */
	public boolean isMixedValence() {

		HashSet<String> seenSymbols = new HashSet();
		for (int symbolNum = 0; symbolNum < m_Symbols.length; symbolNum++) {
			if (m_Counts[symbolNum] < 1E-3) {continue;}
			String symbol = m_Symbols[symbolNum];
			String ionTypeSymbol = IonFactory.getIonTypeSymbol(symbol);
			if (seenSymbols.contains(ionTypeSymbol)) {return true;}
			seenSymbols.add(ionTypeSymbol);
		}
		return false;
	}
	
	@Override
	public int compareTo(TableRow row) {
		return Double.compare(row.getOptimalLikelihood(), m_Optimalikelihood);
	}
	
	@Override
	public String toString() {
		
		String returnString = "";
		DecimalFormat countFormatter = new DecimalFormat("#.##");
		for (int symbolNum = 0; symbolNum < m_OxidationStates.length; symbolNum++) {
			double count = m_Counts[symbolNum];
			String symbol = m_Symbols[symbolNum];
			int oxidationState = m_OxidationStates[symbolNum];
			String countString = countFormatter.format(count);
			returnString += countString;
			returnString += symbol;
			returnString += Math.abs(oxidationState);
			returnString += (oxidationState < 0) ? "-" : "+";
			returnString += " ";
		}
		
		returnString += "\t";
		returnString += this.getOptimalLikelihood() + "\t";
		returnString += this.getOptimalMappedPotential() + "\t";
		returnString += this.getGlobalInstabilityIndex() + "\t";
		returnString += this.isMixedValence();
		return returnString;
	}
	
}
