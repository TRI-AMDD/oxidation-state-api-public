package global.tri.oxidationstates.calculator;

import global.tri.oxidationstates.util.Composition;
import matsci.Element;
import matsci.structure.Structure;
import matsci.util.arrays.ArrayUtils;

/**
 * A general class for calculating the most likely sets of oxidation states for
 * a given set of ion types.
 * 
 * @author timmueller
 *
 */
public abstract class OxidationStateCalculator {

	/**
	 * Calculate the most likely oxidation states for the elements in a given
	 * structure. This method will not identify polyatomic ions. If you would like
	 * to use polyatomic ions, calculate the likely oxidation states by composition
	 * instead.
	 * 
	 * @param structure The given structure
	 * @return the most likely oxidation states for the elements in a given
	 *         structure.
	 */
	public OxidationStateSet getLikelyOxidationStates(Structure structure) {

		// TODO smarter identification of ions
		Element[] elements = structure.getDistinctElements();
		String[] ionTypeIDs = new String[0];
		double[] counts = new double[0];
		for (Element element : elements) {
			if (element == Element.vacancy) {
				continue;
			}
			ionTypeIDs = (String[]) ArrayUtils.appendElement(ionTypeIDs, element.getSymbol());
			double count = structure.numDefiningSitesWithElement(element);
			counts = (double[]) ArrayUtils.appendElement(counts, count);
		}

		return this.getLikelyOxidationStates(ionTypeIDs, counts);
	}

	/**
	 * Calculate the most likely oxidation states for the given composition
	 * 
	 * @param compositionString The given composition
	 * @return the most likely oxidation states for the given composition
	 */
	public OxidationStateSet getLikelyOxidationStates(String compositionString) {
		Composition parser = new Composition(compositionString);
		return this.getLikelyOxidationStates(parser);
	}

	/**
	 * Calculate the most likely oxidation states for the given composition
	 * 
	 * @param composition The given composition
	 * @return the most likely oxidation states for the given composition
	 */
	public OxidationStateSet getLikelyOxidationStates(Composition composition) {

		String[] moleculeIDs = composition.getSymbols().toArray(new String[0]);
		double[] weights = new double[moleculeIDs.length];
		for (int idNum = 0; idNum < moleculeIDs.length; idNum++) {
			weights[idNum] = composition.getCount(moleculeIDs[idNum]);
		}
		return this.getLikelyOxidationStates(moleculeIDs, weights);
	}

	/**
	 * Calculate the most likely oxidation states for the given ion types with the
	 * given composition.
	 * 
	 * @param ionTypeIDs The IDS of the ion types.
	 * @param counts     The composition of each ion type, ordered in the way that
	 *                   corresponds to the ordering of ionTypeIDs.
	 * @return the most likely oxidation states for the given ion types with the
	 *         given composition.
	 */
	public abstract OxidationStateSet getLikelyOxidationStates(String[] ionTypeIDs, double[] counts);

}
