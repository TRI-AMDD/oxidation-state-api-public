package global.tri.oxidationstates.ion;

import matsci.structure.Structure;

/**
 * This is a general interface for all ion finders
 * 
 * @author timmueller
 *
 */
public interface IIonFinder {

	/**
	 * Returns the total number of found polyatomic ions. If more than one ion of
	 * the same type is found in a unit cell in different locations, they are all
	 * counted as separate found ions.
	 * 
	 * @return the total number of found polyatomic ions. If more than one ion of
	 *         the same type is found in a unit cell in different locations, they
	 *         are all counted as separate found ions.
	 */
	public int numFoundPolyatomicIons();

	/**
	 * Return the ionNum'th discovered polyatomic ion
	 * 
	 * @param ionNum The number of the polyatomic ion to return
	 * @return the ionNum'th discovered ion
	 */
	public Structure getFoundPolyatomicIon(int ionNum);

}
