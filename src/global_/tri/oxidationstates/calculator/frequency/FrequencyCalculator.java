package global_.tri.oxidationstates.calculator.frequency;

import java.util.Collection;
import java.util.HashMap;

import global_.tri.oxidationstates.calculator.MixedValenceChargeBalanceFilter;
import global_.tri.oxidationstates.calculator.OxidationStateCalculator;
import global_.tri.oxidationstates.calculator.OxidationStateSet;
import global_.tri.oxidationstates.fitting.OxidationStateData;
import global_.tri.oxidationstates.fitting.OxidationStateData.Entry;
import global_.tri.oxidationstates.ion.IonFactory;
import global_.tri.oxidationstates.ion.IonFactory.Ion;
import matsci.io.app.log.Status;
import matsci.util.MSMath;
import matsci.util.arrays.ArrayIndexer;
import matsci.util.arrays.ArrayUtils;

/**
 * Calculates the frequency score for possible sets of oxidation states
 * 
 * @author timmueller
 *
 */
public class FrequencyCalculator extends OxidationStateCalculator {

	// Allowed oxidation states by ion type ID
	private HashMap<String, int[]> m_OxidationStatesByID;

	// Relative frequency by ion. All frequencies for a particular ion type should
	// add up to 1.
	private HashMap<Ion, Double> m_FrequenciesByIon;

	/**
	 * Initialize the frequency score calculator, generating frequencies per ion
	 * from the user-provided data
	 * 
	 * @param data The frequency of each ion will be determined based on the number
	 *             of entries in this data set that contain the ion.
	 * 
	 */
	public FrequencyCalculator(OxidationStateData data) {

		// Get the data written in terms of polyatomic ions (including entries that
		// contain no polyatomic ions).
		Collection<Entry> polyIonSet = data.getUniqueEntries(true).values();
		OxidationStateData polyIonData = new OxidationStateData(polyIonSet, data.getStructDir());

		// Get the data written in terms of monatomic ions
		Collection<Entry> noPolyIonSet = data.getUniqueEntries(false).values();
		OxidationStateData noPolyIonData = new OxidationStateData(noPolyIonSet, data.getStructDir());

		// Get the allowed oxidation states
		m_OxidationStatesByID = noPolyIonData.getKnownOxidationStates();

		// Add the polyatomic ions
		HashMap<String, int[]> polyIonOxidationStatesByID = polyIonData.getKnownOxidationStates();
		for (String id : polyIonOxidationStatesByID.keySet()) {
			if (m_OxidationStatesByID.containsKey(id)) {
				continue;
			} // We already added this from the monatomic data set
			m_OxidationStatesByID.put(id, polyIonOxidationStatesByID.get(id));
		}

		HashMap<Ion, Integer> countsByIon = noPolyIonData.getCountsByIon();
		HashMap<Ion, Integer> countsByPolyIon = polyIonData.getCountsByIon();

		// Add the counts for the polyatomic ions from the polyatomic data set
		for (Ion ion : countsByPolyIon.keySet()) {
			if (countsByIon.containsKey(ion)) {
				continue;
			}
			countsByIon.put(ion, countsByPolyIon.get(ion));
		}

		// Normalize how many times each ion occurs, so the sum for each type is 1.
		m_FrequenciesByIon = new HashMap<Ion, Double>();
		for (String ionTypeID : m_OxidationStatesByID.keySet()) {
			int[] oxidationStates = m_OxidationStatesByID.get(ionTypeID);

			// Don't allow the zero oxidation state
			int zeroIndex = ArrayUtils.findIndex(oxidationStates, 0);
			if (zeroIndex >= 0) {
				oxidationStates = ArrayUtils.removeElement(oxidationStates, zeroIndex);
				m_OxidationStatesByID.put(ionTypeID, oxidationStates);
			}
			double[] frequencies = new double[oxidationStates.length];

			for (int stateNum = 0; stateNum < oxidationStates.length; stateNum++) {
				int oxidationState = oxidationStates[stateNum];
				Ion ion = IonFactory.get(ionTypeID, oxidationState);
				frequencies[stateNum] = countsByIon.get(ion);
			}

			// Normalization occurs here
			frequencies = MSMath.normalizeToUnitSum(frequencies);
			for (int stateNum = 0; stateNum < oxidationStates.length; stateNum++) {
				int oxidationState = oxidationStates[stateNum];
				Ion ion = IonFactory.get(ionTypeID, oxidationState);
				m_FrequenciesByIon.put(ion, frequencies[stateNum]);
			}
		}

	}

	@Override
	public OxidationStateSet getLikelyOxidationStates(String[] ionTypeIDs, double[] weights) {

		// Set up a branch and bound search
		int[] numStates = new int[ionTypeIDs.length];
		Ion[][] allowedIons = new Ion[ionTypeIDs.length][];
		for (int ionTypeNum = 0; ionTypeNum < ionTypeIDs.length; ionTypeNum++) {
			String ionTypeID = ionTypeIDs[ionTypeNum];
			int[] oxidationStates = m_OxidationStatesByID.get(ionTypeID);
			if (oxidationStates == null) {
				Status.warning("Oxidation states not available for " + ionTypeID);
				return null;
			}

			numStates[ionTypeNum] = oxidationStates.length;
			allowedIons[ionTypeNum] = new Ion[numStates[ionTypeNum]];
			for (int stateNum = 0; stateNum < numStates[ionTypeNum]; stateNum++) {
				double oxidationState = oxidationStates[stateNum];
				allowedIons[ionTypeNum][stateNum] = IonFactory.get(ionTypeID, oxidationState);
			}
		}

		// This is used to avoid combinations that have no hope of reaching the maximum
		// frequency
		MaxFrequencyFilter maxFrequencyFilter = new MaxFrequencyFilter(this, allowedIons);

		// Ensure charge balance while allowing for mixed valence
		MixedValenceChargeBalanceFilter chargeBalanceFilter = new MixedValenceChargeBalanceFilter(allowedIons, weights);
		ArrayIndexer.Filter[] filters = new ArrayIndexer.Filter[] { maxFrequencyFilter, chargeBalanceFilter };
		ArrayIndexer indexer = new ArrayIndexer(numStates);

		int[] stateIndices = indexer.getInitialState(filters);
		if (stateIndices == null) {
			return null;
		} // No match possible
		Ion[] baseIons = new Ion[ionTypeIDs.length];
		OxidationStateSet result = null;
		do {

			double frequencyScore = maxFrequencyFilter.getLastValidFrequencyScore();

			// Set the oxidation states and calculate the net charge
			double netCharge = 0;
			for (int moleculeNum = 0; moleculeNum < stateIndices.length; moleculeNum++) {
				int stateNum = stateIndices[moleculeNum];
				baseIons[moleculeNum] = allowedIons[moleculeNum][stateNum];
				double oxidationState = allowedIons[moleculeNum][stateNum].getOxidationState();
				netCharge += oxidationState * weights[moleculeNum];
			}

			if (Math.abs(netCharge) < 1E-2) { // It's charge balanced
				result = new OxidationStateSet(baseIons, weights);
				maxFrequencyFilter.setMaxKnownFrequencyScore(frequencyScore);
			} else { // Consider mixed valence. We only need to consider lowering valence to cover
						// all cases
				for (int elementNum = 0; elementNum < ionTypeIDs.length; elementNum++) {
					int[] oxidationStates = m_OxidationStatesByID.get(ionTypeIDs[elementNum]);
					int stateIndex = stateIndices[elementNum];
					int oxidationState = oxidationStates[stateIndex];
					for (int decrement = 1; decrement <= stateIndex; decrement++) {
						int newStateIndex = stateIndex - decrement;
						int newOxidationState = oxidationStates[newStateIndex];
						double delta = newOxidationState - oxidationState;
						double maxChange = delta * weights[elementNum];
						if (netCharge + maxChange > 0) {
							continue;
						} // Can't reach zero
						Ion newIon = allowedIons[elementNum][newStateIndex];
						Ion[] allIons = (Ion[]) ArrayUtils.appendElement(baseIons, newIon);

						// Balance the charge
						double[] newWeights = new double[weights.length + 1];
						System.arraycopy(weights, 0, newWeights, 0, weights.length);
						newWeights[newWeights.length - 1] = -netCharge / delta;
						newWeights[elementNum] -= -netCharge / delta;

						// Recalculate the freuqency score
						frequencyScore = 1;
						for (Ion ion : allIons) {
							frequencyScore *= this.getFrequency(ion);
						}

						if (frequencyScore > maxFrequencyFilter.getMaxKnownFrequencyScore()) {
							maxFrequencyFilter.setMaxKnownFrequencyScore(frequencyScore);
							result = new OxidationStateSet(allIons, newWeights);
						}
					}
				}
			}

		} while (indexer.increment(stateIndices, filters));

		return result;
	}

	/**
	 * Returns the frequency for a particular ion, where all frequencies for an ion
	 * type sum to 1.
	 * 
	 * @param ion The ion for which we want the frequency
	 * @return The frequency for the ion, where all frequencies for an ion type sum
	 *         to 1.
	 */
	public double getFrequency(Ion ion) {
		Double returnValue = m_FrequenciesByIon.get(ion);
		return (returnValue == null) ? 0 : returnValue;
	}

	/**
	 * Calculates the frequency for a set of ions by multiplying the frequencies for
	 * all ionsin the set
	 * 
	 * @param ions The set of ions for which we want the frequency
	 * @return The frequency score for the set of ions
	 */
	public double getFrequencyScore(Ion[] ions) {
		double returnValue = 1;
		for (Ion ion : ions) {
			returnValue *= this.getFrequency(ion);
		}
		return returnValue;
	}

}
