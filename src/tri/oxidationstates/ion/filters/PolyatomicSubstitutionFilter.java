package tri.oxidationstates.ion.filters;

import java.util.HashMap;

import matsci.Element;
import matsci.util.arrays.ArrayIndexer.Filter;
import tri.oxidationstates.ion.IonFactory.IonType;

/**
 * This filter is used to consider possible combinations of polyatomic ions that can be created from
 * an atomic composition.  It screens out any combinations of polyatomic ions that can't be formed given
 * the atomic composition (assuming an integer number of polyatomic ions per formula unit).
 * 
 * @author timmueller
 *
 */
public class PolyatomicSubstitutionFilter implements Filter {
	
	private double[] m_MaxCountsByAtomicNumber;
	private double[] m_LastCountsByAtomicNumber;
	private IonType[] m_PolyatomicIonTypes;

	/**
	 * Initialize the filter
	 * 
	 * @param maxCountsByAtomicNumber The indices of this array correspond to atomic numbers, and the values should
	 * correspond to the number of atoms of that type per formula unit.  The formula unit does not need to be fully reduced.
	 * @param polyatomicIonTypes The various types of polyatomic ion types that could be formed using the atoms in this system.
	 * The order of the polyatomic ions should correspond to the state being evaluated, so state[i] is the number of polyatomic
	 * ions of type polyatomicIonTypes[i].
	 */
	public PolyatomicSubstitutionFilter(double[] maxCountsByAtomicNumber, IonType[] polyatomicIonTypes) {
		m_MaxCountsByAtomicNumber = maxCountsByAtomicNumber.clone();
		m_LastCountsByAtomicNumber = maxCountsByAtomicNumber.clone();
		m_PolyatomicIonTypes = polyatomicIonTypes.clone();
	}
	
	@Override
	public int getBranchIndex(int[] currentState) {
		
		System.arraycopy(m_MaxCountsByAtomicNumber, 0, m_LastCountsByAtomicNumber, 0, m_MaxCountsByAtomicNumber.length);
		for (int ionNum = currentState.length - 1; ionNum >= 0; ionNum--) {
			if (currentState[ionNum] == 0) {continue;}
			Element[] elements = m_PolyatomicIonTypes[ionNum].getElements(); 
			int[] counts = m_PolyatomicIonTypes[ionNum].getCounts();
			for (int elementNum = 0; elementNum < elements.length; elementNum++) {
				Element element = elements[elementNum];
				int atomicNumber = element.getAtomicNumber();
				int count = counts[elementNum];
				m_LastCountsByAtomicNumber[atomicNumber] -= count * currentState[ionNum];
				if (m_LastCountsByAtomicNumber[atomicNumber] < 0) { // We've used up more elements in the polyatomic ions than there are available elements
					return ionNum;
				}
			}
		}
		return -1;
	}
	
	/**
	 * Given a composition map representing the polyatomic ion composition, adds the composition of monatomic
	 * ions according to the last state seen by this filter.
	 * 
	 * @param polyatomicIonComposition A map in which the keys are the ion type IDs for polyatomic ions
	 * and the values are amounts.  The elemental amounts will be added to this map.
	 */
	public void addRemainingElementsForLastSeenState(HashMap<String, Double> polyatomicIonComposition) {
		
		for (int atomicNumber = 0; atomicNumber < m_LastCountsByAtomicNumber.length; atomicNumber++) {
			double count = m_LastCountsByAtomicNumber[atomicNumber];
			if (count < 1E-7) {continue;}
			String symbol = Element.getElement(atomicNumber).getSymbol();
			polyatomicIonComposition.put(symbol, count);
		}	
	}
}
