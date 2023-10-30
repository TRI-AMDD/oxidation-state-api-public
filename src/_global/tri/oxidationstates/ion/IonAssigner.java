package _global.tri.oxidationstates.ion;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;

import _global.tri.oxidationstates.calculator.BondValenceCalculator;
import _global.tri.oxidationstates.calculator.OxidationStateSet;
import _global.tri.oxidationstates.ion.IonFactory.Ion;
import _global.tri.structure.mapper.GeneralStructureMapper;
import matsci.Element;
import matsci.Species;
import matsci.structure.PartiallyOccupiedStructure;
import matsci.structure.Structure;
import matsci.structure.Structure.Site;
import matsci.util.MSMath;
import matsci.util.arrays.ArrayUtils;

/**
 * This class assigns ions to sites in a way that minimizes the global
 * instability index
 * 
 * @author timmueller
 *
 */
public class IonAssigner {

	private static double BOND_VALENCE_TOLERANCE = 1E-6; // Two sites need to have nearly identical bv sums to be
															// considered equivalent
	private static double OXIDATION_STATE_TOLERANCE = 0.01;

	private Structure m_GivenStructure;
	private PartiallyOccupiedStructure m_DisorderedStructure;
	private OxidationStateSet m_StateSet;

	private BondValenceCalculator m_BondValenceCalculator;

	private boolean[] m_IsInPolyatomicIon;
	private HashMap<String, ArrayList<int[]>> m_SiteIndicesByIonType = new HashMap<String, ArrayList<int[]>>();

	private double m_GlobalInstabilityIndex = Double.NaN;

	/**
	 * Given an ionFinder (which identifies polyatomic ions) and a set of oxidation
	 * states, assign the oxidation states to atoms in a way that minimizes the GII.
	 * The atomic oxidation states for the polyatomic ions are read from the
	 * polyatomic ion structure files read by {@link IonFactory}. If multiple
	 * different assignments yield nearly identical global instability indices, this
	 * method will generate a partially-occupied structure with occupancies given by
	 * the averages of all of those assignments.
	 * 
	 * @param ionFinder An ion finder that has identified polyatomic ions in the
	 *                  structure to which ions should be assigned.
	 * @param stateSet  The oxidation states to be assigned to the structure.
	 */
	public IonAssigner(KnownIonFinder ionFinder, OxidationStateSet stateSet) {

		m_StateSet = stateSet;
		m_GivenStructure = ionFinder.getStructure();

		m_IsInPolyatomicIon = new boolean[m_GivenStructure.numDefiningSites()];
		for (int foundIonNum = 0; foundIonNum < ionFinder.numFoundPolyatomicIons(); foundIonNum++) {

			String ionType = ionFinder.getFoundIonType(foundIonNum);
			ArrayList<int[]> siteIndices = m_SiteIndicesByIonType.get(ionType);
			if (siteIndices == null) {
				siteIndices = new ArrayList<int[]>();
				m_SiteIndicesByIonType.put(ionType, siteIndices);
			}
			GeneralStructureMapper.Map map = ionFinder.getFoundIonMap(foundIonNum);
			int[] mappedIndices = getMappedIndices(map);
			siteIndices.add(mappedIndices);

			for (int hostIndex : mappedIndices) {
				m_IsInPolyatomicIon[hostIndex] = true;
			}
		}

		m_BondValenceCalculator = new BondValenceCalculator(m_GivenStructure);
		if (m_BondValenceCalculator.getBondValenceSums() != null) {
			m_DisorderedStructure = this.makePartiallyOccupiedStructure(m_GivenStructure);
			double bvDeltaSqSum = this.assignPolyatomicOxidationStates();
			bvDeltaSqSum += this.assignElementalOxidationStates();

			// Important to do this only after all occupancies are assigned.
			this.calculateAverageOccupancies(m_BondValenceCalculator.getBondValenceSums());
			m_GlobalInstabilityIndex = Math.sqrt(bvDeltaSqSum / m_GivenStructure.numDefiningSites());
		}
	}

	/**
	 * Assign the given oxidation states to atom in the structure in a way that
	 * minimizes the GII. If multiple different assignments yield nearly identical
	 * global instability indices, this method will generate a partially-occupied
	 * structure with occupancies given by the averages of all of those assignments.
	 * 
	 * This constructor does not attempt to identify polyatomic ions, so the
	 * provided oxidation states should be for single-atom ions only.
	 * 
	 * @param structure The structure to which oxidation states should be assigned.
	 * @param stateSet  The oxidation states to be assigned to the structure. The
	 *                  state set should not contain any polyatomic ion types when
	 *                  using this construction.
	 */
	public IonAssigner(Structure structure, OxidationStateSet stateSet) {

		m_StateSet = stateSet;
		m_GivenStructure = structure;
		m_IsInPolyatomicIon = new boolean[structure.numDefiningSites()];

		m_BondValenceCalculator = new BondValenceCalculator(structure);
		if (m_BondValenceCalculator.getBondValenceSums() != null) {
			m_DisorderedStructure = this.makePartiallyOccupiedStructure(structure);
			double bvDeltaSqSum = this.assignElementalOxidationStates();
			this.calculateAverageOccupancies(m_BondValenceCalculator.getBondValenceSums());
			m_GlobalInstabilityIndex = Math.sqrt(bvDeltaSqSum / structure.numDefiningSites());
		}
	}

	/**
	 * Makes a structure in which multiple ions can be assigned to each site. This
	 * will be the structure to which oxidation states will be assigned. It is
	 * possible for the ionic fraction at some sites to be between zero and one, if
	 * multiple different ways of assigning oxidation states to the sites yield
	 * nearly identical global instability indices.
	 * 
	 * @param structure The fully occupied structure to which we are assigning
	 *                  oxidation states.
	 * @return a structure in which multiple ions can be assigned to each site. All
	 *         of the sites will be fully occupied, but the returned structure is
	 *         used internally to keep track of which ions can be assigned to which
	 *         sites.
	 */
	private PartiallyOccupiedStructure makePartiallyOccupiedStructure(Structure structure) {

		HashMap<Element, HashSet<Species>> oxidationStatesByElement = new HashMap<Element, HashSet<Species>>();
		Element[] allElements = m_GivenStructure.getDistinctElements();
		for (Element element : allElements) {
			oxidationStatesByElement.put(element, new HashSet<Species>());
		}
		oxidationStatesByElement.put(Element.vacancy, new HashSet<Species>());

		// First deal with the polyatomic ions
		for (String ionType : m_SiteIndicesByIonType.keySet()) {
			for (Ion ion : m_StateSet.getIons()) {
				if (!ion.getIonType().getSymbol().equals(ionType)) {
					continue;
				}
				int intOxidationState = (int) Math.round(ion.getOxidationState());
				if (Math.abs(ion.getOxidationState() - intOxidationState) > OXIDATION_STATE_TOLERANCE) {
					throw new RuntimeException("Non-integer oxidation state found: " + ion.getSymbol());
				}
				Structure ionStructure = IonFactory.getKnownRepresentativeStructures(ionType).get(intOxidationState); // m_PolyatomicIons.getKnownRepresentativeStructures(ionType).get(intOxidationState);
				for (int siteNum = 0; siteNum < ionStructure.numDefiningSites(); siteNum++) {
					Species siteSpecies = ionStructure.getSiteSpecies(siteNum);
					oxidationStatesByElement.get(siteSpecies.getElement()).add(siteSpecies);
				}
			}
		}

		// Now the elemental ions
		for (int siteNum = 0; siteNum < structure.numDefiningSites(); siteNum++) {

			Element element = structure.getSiteSpecies(siteNum).getElement();
			
			// Handle vacancies
			if (element == Element.vacancy) {
				oxidationStatesByElement.get(element).add(Species.vacancy);
			}
			
			// Find the allowed oxidation states
			for (Ion ion : m_StateSet.getIons()) {
				if (!ion.getIonType().isElement()) {
					continue;
				}
				Species neutralSpecies = Species.get(ion.getIonType().getElements()[0]);
				if (neutralSpecies.getElement() != element) {
					continue;
				}
				Species ionSpecies = neutralSpecies.setOxidationState(ion.getOxidationState());
				oxidationStatesByElement.get(element).add(ionSpecies);
			}
		}

		HashMap<Element, Species[]> allowedSpeciesByElement = new HashMap<Element, Species[]>();
		for (Element element : allElements) {
			Species[] allowedSpecies = oxidationStatesByElement.get(element).toArray(new Species[0]);
			allowedSpeciesByElement.put(element, allowedSpecies);
		}

		Species[][] allowedSpeciesBySite = new Species[structure.numDefiningSites()][];
		double[][] occupancies = new double[structure.numDefiningSites()][];
		for (int siteNum = 0; siteNum < allowedSpeciesBySite.length; siteNum++) {
			Element element = m_GivenStructure.getSiteSpecies(siteNum).getElement();
			allowedSpeciesBySite[siteNum] = allowedSpeciesByElement.get(element);
			occupancies[siteNum] = new double[allowedSpeciesBySite[siteNum].length];
			occupancies[siteNum][0] = 1;
		}

		return new PartiallyOccupiedStructure(structure, allowedSpeciesBySite, occupancies);
	}

	/**
	 * Returns an array of site indices correponding to the sites to which a
	 * polyatomic ion has been mapped in the given map. In other words, the site
	 * indices of the atoms in one of the polyatomic ions in the structure. The
	 * indices are returned in the same order as the sites in GeneralStructureMapper.Map.getMatchingHostSites().
	 * 
	 * @param map The map containing the locations of the atoms in the polyatomic
	 *            ions.
	 * @return an array of site indices correponding to the sites to which a
	 *         polyatomic ion has been mapped in the given map. In other words, the
	 *         site indices of the atoms in one of the polyatomic ions in the
	 *         structure. The indices are returned in the same order as the sites in
	 *         GeneralStructureMapper.Map.getMatchingHostSites().
	 */
	public int[] getMappedIndices(GeneralStructureMapper.Map map) {

		Site[] matchingHostSites = map.getMatchingHostSites();
		int[] returnArray = new int[matchingHostSites.length];
		for (int siteNum = 0; siteNum < returnArray.length; siteNum++) {
			Site mappedSite = matchingHostSites[siteNum];
			Site origSite = m_GivenStructure.getDefiningSite(mappedSite.getCoords());
			if (origSite.getSpecies().getElement() != mappedSite.getSpecies().getElement()) {
				throw new RuntimeException("Polyatomic mapping failure: element " + origSite.getSpecies().getElement()
						+ " found in structure where " + mappedSite.getSpecies().getElement() + " was expected.");
			}
			returnArray[siteNum] = origSite.getIndex();
		}

		return returnArray;
	}

	/**
	 * This method is used to simplify the process of assigning oxidation states in
	 * a way that minimizes the GII for polyatomic ions. In the polyatomic ions
	 * considered for the manuscript and web site, when the ion changes oxidation
	 * state only one atom in the ion alsoc hanges oxidation states. All of the
	 * other atoms in the polyatomic ions always have the same oxidation state,
	 * regardless of what the oxidation state of the polyatomic ion is. The atom
	 * that changes oxidation states is called the "active atom", and here we return
	 * the index of that atom in the reptresentative structure of the given
	 * polyatomic ion.
	 * 
	 * Note that another approach will need to be used if we extend the set of
	 * polyatomic ions to include ions that don't have a single "active atom".
	 * 
	 * @param polyatomicSymbol The ion type symbol of the polyatomic ion for which
	 *                         we are trying to find the active atom.
	 * @return The index of the active atom in the representative structure for the given ion.
	 */
	public int getActiveAtomIndex(String polyatomicSymbol) {

		Map<Integer, Structure> structures = IonFactory.getKnownRepresentativeStructures(polyatomicSymbol); // m_PolyatomicIons.getKnownRepresentativeStructures(polyatomicSymbol);

		Iterator<Structure> iterator = structures.values().iterator();
		Structure baseStructure = iterator.next();
		while (iterator.hasNext()) {
			Structure structure = iterator.next();
			for (int siteNum = 0; siteNum < structure.numDefiningSites(); siteNum++) {
				double baseOxidationState = baseStructure.getSiteSpecies(siteNum).getOxidationState();
				double oxidationState = structure.getSiteSpecies(siteNum).getOxidationState();
				if (Math.abs(baseOxidationState - oxidationState) > OXIDATION_STATE_TOLERANCE) {
					return siteNum;
				}
			}
		}

		return -1;
	}

	/**
	 * Assigns oxidation states for polyatomic ions to the given structure in a way
	 * that minimizes the global instability index.
	 * 
	 * @return The sum of the squared difference between the assigned oxidation
	 *         states and the bond valence sums on each of the atoms assigned an
	 *         oxidation state.
	 */
	private double assignPolyatomicOxidationStates() {

		/**
		 * We assign states by first calculating bond valence sums, and then assigning
		 * the states from smallest to largest in order of increasing bond valence sums
		 */
		double[] bvSums = m_BondValenceCalculator.getBondValenceSums();
		if (bvSums == null) {
			return Double.NaN;
		}

		/**
		 * Used to calculate the Global Instability Index (gii), a measure of how well
		 * the assigned states match up with the bond valence sums.
		 */
		double bvDeltaSqSum = 0;

		Ion[] oxidationIons = m_StateSet.getIons();

		for (String ionType : m_SiteIndicesByIonType.keySet()) {
			ArrayList<int[]> maps = m_SiteIndicesByIonType.get(ionType);
			ArrayList<Ion> matchingIons = new ArrayList<Ion>();
			for (Ion ion : oxidationIons) {
				if (ion.getIonType().getSymbol().equals(ionType)) {
					matchingIons.add(ion);
				}
			}

			if (matchingIons.size() == 0) {
				continue;
			}
			if (matchingIons.size() == 1) {
				Ion ion = matchingIons.get(0);
				double oxidationState = ion.getOxidationState();
				int intOxidationState = (int) Math.round(oxidationState);
				if (Math.abs(oxidationState - intOxidationState) > OXIDATION_STATE_TOLERANCE) {
					throw new RuntimeException("Non-integer oxidation state found: " + ion.getSymbol());
				}
				for (int[] siteIndices : maps) {
					bvDeltaSqSum += this.assignPolyatomicIon(ionType, new int[] { intOxidationState }, siteIndices,
							bvSums, new double[] { 1 });
				}
				continue;
			}

			/**
			 * Here we use the fact that if only one atom changes its oxidation state every
			 * time the ion changes its oxidation state, then we can just sort the bv sums
			 * for the polyatomic ions like we do for single-atom ions.
			 */
			int activeIndex = this.getActiveAtomIndex(ionType);

			// For now we work with integer weights
			double[] weights = new double[matchingIons.size()];
			for (int ionNum = 0; ionNum < matchingIons.size(); ionNum++) {
				Ion ion = matchingIons.get(ionNum);
				weights[ionNum] = m_StateSet.getWeight(ion);
			}
			weights = normalizeWeights(weights, maps.size());

			/**
			 * Initialize the set of bond valence values for each polyatomic ion in the
			 * structure
			 */
			double[] activeBVByMolecule = new double[maps.size()];
			for (int mapNum = 0; mapNum < maps.size(); mapNum++) {
				int[] siteIndices = maps.get(mapNum);
				activeBVByMolecule[mapNum] = bvSums[siteIndices[activeIndex]]; // The bond valence sum should be the
																				// same for all oxidation states
			}
			int[] bvSortMap = ArrayUtils.getSortPermutation(activeBVByMolecule);

			/**
			 * Now get the oxidation states of each site for each possible oxidation state
			 * of the ion
			 */
			int[] oxidationStates = new int[matchingIons.size()];
			for (int ionNum = 0; ionNum < matchingIons.size(); ionNum++) {
				Ion ion = matchingIons.get(ionNum);
				oxidationStates[ionNum] = (int) Math.round(ion.getOxidationState());
			}
			int[] stateMap = ArrayUtils.getSortPermutation(oxidationStates);

			int ionNum = 0;
			for (int mapNum = 0; mapNum < bvSortMap.length; mapNum++) {
				int[] siteIndices = maps.get(bvSortMap[mapNum]);
				double[] occupancies = new double[oxidationStates.length];

				/**
				 * We figure out the occupancies for this site based on ascending BV sums
				 */
				double remainingWeight = 1; // This is how much total weight is available for the current oxidation
											// state
				while (remainingWeight > 0) {
					double assignedWeight = Math.min(weights[stateMap[ionNum]], remainingWeight);
					occupancies[stateMap[ionNum]] += assignedWeight;
					weights[stateMap[ionNum]] -= assignedWeight;
					if (weights[stateMap[ionNum]] == 0) {
						ionNum++;
					}
					remainingWeight -= assignedWeight;
				}

				bvDeltaSqSum += this.assignPolyatomicIon(ionType, oxidationStates, siteIndices, bvSums, occupancies);
			}
		}

		return bvDeltaSqSum;

	}

	/**
	 * Assigns oxidation states to all atoms in the given polyatomic ion
	 * 
	 * @param ionType         The type of polyatomic ion for which we are assigning
	 *                        oxidation states
	 * @param oxidationStates The allowed oxidation states for the polyatomic ion
	 * @param siteIndices     The indices of the sites in the structure to which we
	 *                        are assigning oxidation states, corresponding to the
	 *                        respective sites in the representative structure for
	 *                        the polyatomic ion.
	 * @param bvSums          The bond valence sums for all sites in a unit cell of
	 *                        the structure to which we are assigning oxidation
	 *                        states, sorted by the site indices of the sites.
	 * @param occupancies     An array corresponding to the oxidation states array,
	 *                        where every element gives the fractional occupancy of
	 *                        the corresponding oxidation state.
	 *
	 * @return The sum of the squared difference between the assigned oxidation
	 *         states and the bond valence sums on each of the atoms assigned an
	 *         oxidation state.
	 */
	private double assignPolyatomicIon(String ionType, int[] oxidationStates, int[] siteIndices, double[] bvSums,
			double[] occupancies) {

		double bvDeltaSqSum = 0;

		double[][] siteOccupancies = new double[siteIndices.length][];
		for (int siteNum = 0; siteNum < siteOccupancies.length; siteNum++) {
			siteOccupancies[siteNum] = new double[m_DisorderedStructure.numAllowedSpecies(siteIndices[siteNum])];
		}

		for (int stateNum = 0; stateNum < oxidationStates.length; stateNum++) {
			int oxidationState = oxidationStates[stateNum];
			double occupancy = occupancies[stateNum];
			Structure representativeStructure = IonFactory.getKnownRepresentativeStructures(ionType)
					.get(oxidationState);

			for (int siteNum = 0; siteNum < siteIndices.length; siteNum++) {
				int siteIndex = siteIndices[siteNum];
				Species species = representativeStructure.getSiteSpecies(siteNum);
				Site sourceHostSite = m_GivenStructure.getDefiningSite(siteIndex);
				if (species.getElement() != sourceHostSite.getSpecies().getElement()) {
					throw new RuntimeException("Cannot assign species " + species + " to site occupied by element "
							+ sourceHostSite.getSpecies().getElementSymbol() + " when assigning oxidation states.");
				}

				int specIndex = m_DisorderedStructure.getIndexForAllowedSpecie(siteIndex, species);
				siteOccupancies[siteNum][specIndex] += occupancy;
				double bvValue = bvSums[sourceHostSite.getIndex()];
				double delta = bvValue - species.getOxidationState();
				bvDeltaSqSum += delta * delta * occupancy;
			}
		}

		for (int siteNum = 0; siteNum < siteIndices.length; siteNum++) {
			int hostSiteIndex = siteIndices[siteNum];
			m_DisorderedStructure.setOccupancies(hostSiteIndex, siteOccupancies[siteNum]);
		}

		return bvDeltaSqSum;
	}

	/**
	 * Assigns oxidation states to the monatomic ion types in the given structure
	 * 
	 * @return The sum of the squared difference between the assigned oxidation
	 *         states and the bond valence sums on each of the atoms assigned an
	 *         oxidation state.
	 */
	private double assignElementalOxidationStates() {

		/**
		 * We assign states by first calculating bond valence sums, and then assigning
		 * the states from smallest to largest in order of increasing bond valence sums
		 */
		double[] bvSums = m_BondValenceCalculator.getBondValenceSums();
		if (bvSums == null) {
			return Double.NaN;
		}

		/**
		 * The is is used to calculate the Global Instability Index (gii), a measure of
		 * how well the assigned states match up with the bond valence sums.
		 */
		double bvDeltaSqSum = 0;

		// This sorts the indices from smallest BV sum to largest
		int[] bvMap = ArrayUtils.getSortPermutation(bvSums);

		Element[] elements = m_GivenStructure.getDistinctElements();
		Ion[] ions = m_StateSet.getIons();

		/**
		 * Find the ions corresponding to each element in the structure. We need to
		 * update this to work with polyatomic ions.
		 */
		for (Element element : elements) {
			ArrayList<Ion> matchingIons = new ArrayList<Ion>();
			for (Ion ion : ions) {
				if (ion.getIonType().getSymbol().equals(element.getSymbol())) {
					matchingIons.add(ion);
				}
			}

			if (matchingIons.size() == 0) { // Can't assign the oxidation state. Probably in a polyatomic ion
				continue;
			}

			// There is only one possibility, so we make the assignment.
			if (matchingIons.size() == 1) {

				// This is the only possible oxidation state for this ion.
				double oxidationState = matchingIons.get(0).getOxidationState();

				// We get the corresponding species (element + oxidation state);
				Species species = Species.get(element).setOxidationState(oxidationState);

				/**
				 * Cycle through all of the sites in the structure and find all that match this
				 * element.
				 */
				for (int siteIndex = 0; siteIndex < m_GivenStructure.numDefiningSites(); siteIndex++) {

					// Don't assign the site if it's in a polyatomic ion
					if (m_IsInPolyatomicIon[siteIndex]) {
						continue;
					}

					// We will use the bond valence sum to update the global instability index
					double bvSum = bvSums[siteIndex];

					// If site has the matching element, assign the oxidation state and update the
					// global instability index
					if (m_GivenStructure.getSiteSpecies(siteIndex).getElement() == element) {
						double delta = bvSum - oxidationState;
						bvDeltaSqSum += delta * delta;
						double[] occupancies = new double[m_DisorderedStructure.numAllowedSpecies(siteIndex)];
						int stateNum = m_DisorderedStructure.getIndexForAllowedSpecie(siteIndex, species);
						occupancies[stateNum] = 1;
						m_DisorderedStructure.setOccupancies(siteIndex, occupancies);
					}
				}

				// Go on to the next element
				continue;
			}

			/**
			 * There are multiple oxidation states corresponding to this element. So we need
			 * to assign the smallest oxidation states to the sites with the smallest bond
			 * valence sums, and so on up the line.
			 */

			// The weights tell us the relative fraction of ions with each oxidation state
			double[] weights = new double[matchingIons.size()];

			// These should correspond to the weights
			double[] oxidationStates = new double[matchingIons.size()];

			/**
			 * We initialize the weights and oxidation states arrays, and keep track of the
			 * total weight for normalization later.
			 */
			for (int ionNum = 0; ionNum < matchingIons.size(); ionNum++) {
				Ion ion = matchingIons.get(ionNum);
				weights[ionNum] = m_StateSet.getWeight(ion);
				oxidationStates[ionNum] = ion.getOxidationState();
			}
			int numElementalSites = 0;
			for (int siteNum = 0; siteNum < m_GivenStructure.numDefiningSites(); siteNum++) {
				if (m_IsInPolyatomicIon[siteNum]) {
					continue;
				}
				if (m_GivenStructure.getSiteSpecies(siteNum).getElement() == element) {
					numElementalSites++;
				}
			}
			weights = normalizeWeights(weights, numElementalSites);

			// Gets the order of oxidation states from smallest to largest.
			int[] stateMap = ArrayUtils.getSortPermutation(oxidationStates);

			/**
			 * Now we go through the sites from smallest bond valence sum to largest and
			 * assign the ions accordingly.
			 */
			int ionNum = 0;
			for (int siteNum = 0; siteNum < bvMap.length; siteNum++) {

				// This orders the sites from smallest BV sum to largest.
				int siteIndex = bvMap[siteNum];

				// Don't assign the site if it's in a polyatomic ion
				if (m_IsInPolyatomicIon[siteIndex]) {
					continue;
				}

				// Ignore sites that don't match the current element.
				Species baseSpecies = m_GivenStructure.getSiteSpecies(siteIndex);
				if (baseSpecies.getElement() != element) {
					continue;
				}

				/**
				 * We figure out the occupancies for this site based on ascending BV sums
				 */
				double[] occupancies = new double[m_DisorderedStructure.numAllowedSpecies(siteIndex)];
				double remainingWeight = 1; // This is how much total weight is available for the current oxidation
											// state
				while (remainingWeight > 1E-10) { // For numerical noise
					double assignedWeight = Math.min(weights[stateMap[ionNum]], remainingWeight);
					weights[stateMap[ionNum]] -= assignedWeight;
					Species species = baseSpecies.setOxidationState(oxidationStates[stateMap[ionNum]]);
					int specNum = m_DisorderedStructure.getIndexForAllowedSpecie(siteIndex, species);
					occupancies[specNum] += assignedWeight;
					double delta = oxidationStates[stateMap[ionNum]] - bvSums[siteIndex];
					bvDeltaSqSum += delta * delta * assignedWeight;
					if (weights[stateMap[ionNum]] == 0) {
						ionNum++;
					}
					remainingWeight -= assignedWeight;
				}
				// TODO round off occupancies to get rid of numerical issues?
				// Assign the oxidation state
				m_DisorderedStructure.setOccupancies(siteIndex, occupancies);
			}

		}

		return bvDeltaSqSum;

	}

	/**
	 * Return the structure with assigned oxidation states
	 * 
	 * @return the structure with assigned oxidation states
	 */
	public PartiallyOccupiedStructure getStructure() {
		return m_DisorderedStructure;
	}

	/**
	 * Returns the structure provided in the constructor. Oxidation states are not
	 * assigned directly to this structure. If you want a structure with assigned
	 * oxidation states, call the {@link getStructure()} method.
	 * 
	 * @return the structure provided in the constructor. Oxidation states are not
	 * assigned directly to this structure. If you want a structure with assigned
	 * oxidation states, call the {@link getStructure()} method.
	 */
	public Structure getGivenStructure() {
		return m_GivenStructure;
	}

	/**
	 * Returns the calculated global instability index for the assigned oxidation
	 * states.
	 * 
	 * @return the calculated global instability index for the assigned oxidation
	 *         states.
	 */
	public double getGlobalInstabilityIndex() {
		return m_GlobalInstabilityIndex;
	}

	/**
	 * Returns the bond valence calculator used to calculate bond valence sums.
	 * 
	 * @return the bond valence calculator used to calculate bond valence sums.
	 */
	public BondValenceCalculator getBondValenceCalculator() {
		return m_BondValenceCalculator;
	}

	/**
	 * If there are multiple different ways to assign bond valence sums that yield
	 * nearly identical global instability indices, then rather than choose just one
	 * we assign partial occupancies to the sites of the structure by taking the
	 * averages over all of the different ionic arrangements that have nearly the
	 * same global instability index.
	 * 
	 * @param bvSums The bond valence sums for the structure to which we are
	 *               assigning oxidation states, sorted by the site index.
	 */
	private void calculateAverageOccupancies(double[] bvSums) {

		int[] map = ArrayUtils.getSortPermutation(bvSums);
		for (Element element : m_GivenStructure.getDistinctElements()) {
			ArrayList<Integer> siteList = new ArrayList<Integer>();
			double bvAvg = 0;
			for (int siteNum = 0; siteNum < bvSums.length; siteNum++) {
				int siteIndex = map[siteNum];
				if (m_GivenStructure.getSiteSpecies(siteIndex).getElement() != element) {
					continue;
				}
				double bvSum = bvSums[siteIndex];
				// Add it to the current set of sites
				if (siteList.size() == 0 || Math.abs(bvSum - bvAvg) < BOND_VALENCE_TOLERANCE) {
					bvAvg = ((bvAvg * siteList.size()) + bvSum) / (siteList.size() + 1);
					siteList.add(siteIndex);
					continue;
				}

				// We've found a new list, so record the old one.
				double[] avgOccupancies = new double[m_DisorderedStructure.numAllowedSpecies(siteList.get(0))];
				for (Integer seenSiteIndex : siteList) {
					double[] occupancies = m_DisorderedStructure.getSiteOccupancies(seenSiteIndex);
					avgOccupancies = MSMath.fastArrayAdd(avgOccupancies, occupancies);
				}
				MSMath.arrayDivideInPlace(avgOccupancies, siteList.size());
				for (Integer seenSiteIndex : siteList) {
					m_DisorderedStructure.setOccupancies(seenSiteIndex, avgOccupancies);
				}

				// Start the new site list
				siteList.clear();
				siteList.add(siteIndex);
				bvAvg = bvSum;
			}
		}

	}

	/**
	 * Normalizes the given weights to that the sum of all weights is equal to the
	 * given ionTypeCount.
	 * 
	 * @param weights      Weights corresponding to the relative number of atoms in
	 *                     the material with each oxidation state.
	 * @param ionTypeCount The sum of the returned weights will be this value.
	 *                     Typically this would be the total number of ions of this
	 *                     type per unit cell, so the elements of the returned array
	 *                     given the total number of atoms per unit cell with the
	 *                     corresponding oxidation state.
	 * @return An array in which the weights are normalized so that they sum to the
	 *         ionTypeCount.
	 */
	public static double[] normalizeWeights(double[] weights, int ionTypeCount) {

		double totalWeight = MSMath.arraySum(weights);
		double[] returnArray = new double[weights.length];

		double weightPerAtom = totalWeight / ionTypeCount;
		for (int ionNum = 0; ionNum < weights.length; ionNum++) {
			returnArray[ionNum] = weights[ionNum] / weightPerAtom;
		}

		return returnArray;

	}
}
