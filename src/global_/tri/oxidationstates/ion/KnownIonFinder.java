package global_.tri.oxidationstates.ion;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import global.tri.structure.mapper.GeneralStructureMapper;
import global_.tri.oxidationstates.ion.IonFactory.Ion;
import global_.tri.oxidationstates.ion.IonFactory.IonType;
import global_.tri.oxidationstates.util.Composition;
import global_.tri.oxidationstates.util.IonTools;
import matsci.Element;
import matsci.Species;
import matsci.location.Coordinates;
import matsci.location.Vector;
import matsci.location.basis.CartesianBasis;
import matsci.structure.BravaisLattice;
import matsci.structure.IStructureData;
import matsci.structure.Structure;
import matsci.structure.Structure.Site;
import matsci.util.arrays.ArrayUtils;

/**
 * This class identifies known polyatomic ions in a given structure. If the
 * known ion exists in a larger cluster containing the same atoms, it is not
 * added to the list of found ions. If two known ions share share a common site
 * (or sites), the larger is added to the list of found ions. If the two known
 * ions are the same size, neither is added to the list of known ions.
 * 
 * @author timmueller
 *
 */
public class KnownIonFinder implements IIonFinder {

	private final Structure m_Structure;

	private ArrayList<FoundIonData> m_FoundIonData = new ArrayList<FoundIonData>();
	private FoundIonData[] m_IonDataBySite;
	private Composition m_IonTypeComposition;
	private Composition m_IonComposition;

	/**
	 * Creates an ion finder and finds the ions.
	 * 
	 * @param structure          The structure in which we are searching for ions.
	 * @param findPolyatomicIons True if this finder should search for polyatomic
	 *                           ions, and false otherwise.
	 */
	public KnownIonFinder(Structure structure, boolean findPolyatomicIons) {

		m_Structure = structure;

		m_IonDataBySite = new FoundIonData[structure.numDefiningSites()];

		if (findPolyatomicIons) {
			Collection<IonType> ionTypes = IonFactory.getKnownIonTypes().values();
			for (IonType ionType : ionTypes) {
				if (ionType.isElement()) {
					continue;
				}
				this.findPolyatomicIons(ionType);
			}
		}
	}

	/**
	 * Find polyatomic ions of the given type. There may be multiple such ions per
	 * unit cell. If two found ions share the same site(s), the larger one is kept.
	 * If they are the same size, neither is kept.
	 * 
	 * @param ionType The type of ion to be found.
	 */
	private void findPolyatomicIons(IonType ionType) {

		Element[] elements = ionType.getElements();
		int[] targetCounts = ionType.getCounts();
		String ionTypeSymbol = ionType.getSymbol();
		Map<Integer, Structure> ionStructures = IonFactory.getKnownRepresentativeStructures(ionType.getSymbol());

		// Only look for clusters of atoms with sites containing the elements in this
		// ion
		boolean[] allowedSites = new boolean[m_Structure.numDefiningSites()];
		for (int siteNum = 0; siteNum < allowedSites.length; siteNum++) {
			Element element = m_Structure.getSiteSpecies(siteNum).getElement();
			allowedSites[siteNum] = ArrayUtils.arrayContains(elements, element);
		}

		for (int siteNum = 0; siteNum < allowedSites.length; siteNum++) {
			if (!allowedSites[siteNum]) {
				continue;
			}
			Structure.Site site = m_Structure.getDefiningSite(siteNum);
			IonBuilder builder = new IonBuilder(site, allowedSites);
			if (builder.numDefiningSites() <= 1 && builder.isAperiodic()) {
				continue;
			} // Don't find monatomic ions
			if (!builder.matchesCounts(elements, targetCounts)) {
				continue;
			}
			if (builder.checkForOverlappingIons()) {
				Structure foundMolecule = new Structure(builder);
				double bestScore = Double.POSITIVE_INFINITY;
				GeneralStructureMapper.Map bestMap = null;
				for (Structure ionStructure : ionStructures.values()) {
					GeneralStructureMapper mapper = new GeneralStructureMapper(foundMolecule.removeOxidationStates(),
							ionStructure.removeOxidationStates());
					IonTools.setMappingParams(mapper);
					GeneralStructureMapper.Map map = mapper.getBestMap();
					if (map == null) {
						continue;
					}
					if (map.scoreMap() < bestScore) {
						bestMap = map;
						bestScore = map.scoreMap();
					}
				}
				if (bestMap == null) {
					continue;
				}
				FoundIonData ionData = new FoundIonData(ionTypeSymbol, foundMolecule, bestMap);

				for (int builderSiteIndex = 0; builderSiteIndex < builder.numDefiningSites(); builderSiteIndex++) {
					m_IonDataBySite[builder.getSite(builderSiteIndex).getIndex()] = ionData;
				}

				m_FoundIonData.add(ionData);
			}
		}
	}

	/**
	 * Returns true if the atom at the given site is contained within a polyatomic
	 * ion, and false otherwise.
	 * 
	 * @param siteIndex THe index of the given site
	 * @return true if the atom at the given site is contained within a polyatomic
	 *         ion, and false otherwise.
	 */
	public boolean isInPolyatomicIon(int siteIndex) {
		return (m_IonDataBySite[siteIndex] == null);
	}

	@Override
	public int numFoundPolyatomicIons() {
		return m_FoundIonData.size();
	}

	@Override
	public Structure getFoundPolyatomicIon(int index) {
		return m_FoundIonData.get(index).getFoundStructure();
	}

	/**
	 * Returns the ion type of the index'th found polyatomic ion.
	 * 
	 * @param index The index of the found ion
	 * @return the ion type of the index'th found polyatomic ion.
	 */
	public String getFoundIonType(int index) {
		return m_FoundIonData.get(index).getIonStructureSymbol();
	}

	/**
	 * Returns a map that connecting the sites in the representative structure for
	 * the ion to sites in the given structure for the index'th found polyatomic
	 * ion.
	 * 
	 * @param index The index of the polyatomic ion.
	 * @return a map that connecting the sites in the representative structure for
	 *         the ion to sites in the given structure for the index'th found
	 *         polyatomic ion.
	 */
	public GeneralStructureMapper.Map getFoundIonMap(int index) {
		return m_FoundIonData.get(index).getMap();
	}

	/**
	 * Return the structure in which we are searching for ions.
	 * 
	 * @return the structure in which we are searching for ions.
	 */
	public Structure getStructure() {
		return m_Structure;
	}

	/**
	 * Returns the composition in terms of found ions.
	 * 
	 * @return the composition in terms of found ions.
	 */
	public Composition getIonComposition() {
		if (m_IonComposition == null) {
			this.findCompositions();
		}
		return m_IonComposition;
	}

	/**
	 * Returns the composition in terms of found ion types (i.e. without including
	 * oxidation states).
	 * 
	 * @return the composition in terms of found ion types (i.e. without including
	 *         oxidation states).
	 */
	public Composition getIonTypeComposition() {
		if (m_IonTypeComposition == null) {
			this.findCompositions();
		}
		return m_IonTypeComposition;
	}

	/**
	 * Determines both the composition in terms of found polyatomic ions and the
	 * composition in terms of found polyatomic ion types (i.e. without oxidation
	 * states).
	 */
	private void findCompositions() {

		Map<String, Structure> knownIons = IonFactory.getAllKnownRepresentativeStructures(true); // m_PolyatomicIons.getAllKnownRepresentativeStructures();
		TreeMap<String, Integer> ionCounts = new TreeMap<String, Integer>();
		TreeMap<String, Integer> ionTypeCounts = new TreeMap<String, Integer>();
		for (int siteNum = 0; siteNum < m_Structure.numDefiningSites(); siteNum++) {
			Species species = m_Structure.getSiteSpecies(siteNum);
			String ionSymbol = species.getSymbol();
			int ionCount = ionCounts.containsKey(ionSymbol) ? ionCounts.get(ionSymbol) : 0;
			ionCounts.put(ionSymbol, ionCount + 1);

			String ionTypeSymbol = IonFactory.getIonTypeSymbol(ionSymbol); // ion.getStructureSymbol();
			int ionTypeCount = ionTypeCounts.containsKey(ionTypeSymbol) ? ionTypeCounts.get(ionTypeSymbol) : 0;
			ionTypeCounts.put(ionTypeSymbol, ionTypeCount + 1);
		}

		for (int ionNum = 0; ionNum < this.numFoundPolyatomicIons(); ionNum++) {
			Structure foundIon = this.getFoundPolyatomicIon(ionNum);
			for (String knownIonID : knownIons.keySet()) {
				Structure knownIon = knownIons.get(knownIonID);
				if (!IonTools.compareStructures(foundIon, knownIon, true)) {
					continue;
				}
				double oxidationState = 0;
				String ionTypeID = IonFactory.getIonTypeSymbol(knownIonID);
				for (int siteNum = 0; siteNum < foundIon.numDefiningSites(); siteNum++) {
					Species species = foundIon.getSiteSpecies(siteNum);
					Ion ion = IonFactory.get(species);
					int ionCount = ionCounts.get(ion.getSymbol());
					ionCounts.put(ion.getSymbol(), ionCount - 1);
					int ionTypeCount = ionTypeCounts.get(ion.getIonType().getSymbol());
					ionTypeCounts.put(ion.getIonType().getSymbol(), ionTypeCount - 1);
					oxidationState += ion.getOxidationState();
				}
				Ion ion = IonFactory.get(ionTypeID, oxidationState);
				int count = ionCounts.containsKey(ion.getSymbol()) ? ionCounts.get(ion.getSymbol()) : 0;
				ionCounts.put(ion.getSymbol(), count + 1);

				String ionTypeSymbol = ion.getIonType().getSymbol();
				count = ionTypeCounts.containsKey(ionTypeSymbol) ? ionTypeCounts.get(ionTypeSymbol) : 0;
				ionTypeCounts.put(ionTypeSymbol, count + 1);

				break; // We only find one match, then stop.
			}
		}

		m_IonComposition = new Composition(ionCounts);
		m_IonTypeComposition = new Composition(ionTypeCounts);
	}

	/**
	 * This class starts with the given site and recursively searches for any known
	 * polyatomic ion containing the site.
	 * 
	 * @author timmueller
	 *
	 */
	private class IonBuilder implements IStructureData {

		// Atoms will be considered neigbhors if there distance is less than the sum of
		// their covalent radii times this factor.
		private static double RADIUS_FACTOR = 1.1;

		private Site[] m_Sites = new Site[m_Structure.numDefiningSites()];
		private int m_NumSites = 0;
		private BravaisLattice m_Lattice = new BravaisLattice(new Vector[0]);

		/**
		 * Constructs an ion builder to search for ions from the given initial site.
		 * 
		 * @param initialSite  The site from which the search should be initiated.
		 * @param allowedSites The sites available to be added to this ion. The index of
		 *                     this array corresponds to the site indices, and the value
		 *                     is "true" if the site is available for this ion and
		 *                     "false" otherwise.
		 */
		private IonBuilder(Site initialSite, boolean[] allowedSites) {

			double maxCovalentRadius = 0;
			for (int siteNum = 0; siteNum < allowedSites.length; siteNum++) {
				if (!allowedSites[siteNum]) {
					continue;
				}
				Element element = m_Structure.getDefiningSite(siteNum).getSpecies().getElement();
				maxCovalentRadius = Math.max(maxCovalentRadius, element.getCovalentRadius());
			}
			this.addSites(initialSite, allowedSites, maxCovalentRadius);
		}

		/**
		 * Checks to see if this ion overlaps with any others. If so, the larger ion is
		 * kept. If they are the same size, neither is kept.
		 * 
		 * @return True if this ion should be kept, false otherwise.
		 */
		public boolean checkForOverlappingIons() {
			boolean keeper = true;
			for (Site site : m_Sites) {
				if (site == null) {
					continue;
				}
				FoundIonData prevIonData = m_IonDataBySite[site.getIndex()];
				if (prevIonData != null) {
					if (prevIonData.getFoundStructure().numDefiningSites() > this.m_NumSites) {
						return false;
					} else {
						for (int siteNum = 0; siteNum < m_IonDataBySite.length; siteNum++) {
							if (m_IonDataBySite[siteNum] == prevIonData) {
								m_IonDataBySite[siteNum] = null;
							}
						}
						m_FoundIonData.remove(prevIonData);
						// Only keep the new one if it's larger than the old one
						keeper |= (prevIonData.getFoundStructure().numDefiningSites() < this.m_NumSites);
					}
				}
			}
			return keeper;
		}

		/**
		 * Checks to see whether the composition of the found cluster of sites matches
		 * the composition given by elements and targetCounts.
		 * 
		 * @param elements     The elements we are comparing.
		 * @param targetCounts The number of atoms of each element. The indices of this
		 *                     array correspond to the indices of the elements array.
		 * @return True if the found ion matches the given composition, and false
		 *         otherwise.
		 */
		private boolean matchesCounts(Element[] elements, int[] targetCounts) {

			for (int elementIndex = 0; elementIndex < elements.length; elementIndex++) {
				Element element = elements[elementIndex];
				int targetCount = targetCounts[elementIndex];
				int count = 0;
				for (Site site : m_Sites) {
					if (site == null) {
						continue;
					}
					if (site.getSpecies().getElement() != element) {
						continue;
					}
					count++;
				}
				if (count != targetCount) {
					return false;
				}
			}
			return true;
		}

		/**
		 * Tries to add this site to the ion. If successful, recursively tries to add
		 * neighboring sites (as determined by covalent radii).
		 * 
		 * @param currentSite       The site we are trying to add.
		 * @param allowedSites      THe indices of this array are site indices. The
		 *                          values are "true" if the corresponding site is
		 *                          available to be added to tihs ion, and "false"
		 *                          otherwise.
		 * @param maxCovalentRadius The maximum covalent radius of any of the allowed
		 *                          sites, used for searching for neighbors.
		 * @return true of this site was added, false otherwise
		 */
		private boolean addSites(Site currentSite, boolean[] allowedSites, double maxCovalentRadius) {

			Site knownSite = m_Sites[currentSite.getIndex()];
			if (knownSite == null) {
				if (!allowedSites[currentSite.getIndex()]) {
					return false;
				}
				m_Sites[currentSite.getIndex()] = currentSite;
				m_NumSites++;
				allowedSites[currentSite.getIndex()] = false;
				double currRadius = currentSite.getSpecies().getElement().getCovalentRadius();
				double searchDistance = RADIUS_FACTOR * (currRadius + maxCovalentRadius);
				Structure.Site[] neighbors = m_Structure.getNearbySites(currentSite.getCoords(), searchDistance, false);
				for (Site neighbor : neighbors) {
					double neighborRadius = neighbor.getSpecies().getElement().getCovalentRadius();
					double maxAllowedDistance = RADIUS_FACTOR * (currRadius + neighborRadius);
					double distance = neighbor.distanceFrom(currentSite);
					if (distance <= maxAllowedDistance) {
						addSites(neighbor, allowedSites, maxCovalentRadius);
					}
				}
				return true;
			}

			Vector translation = new Vector(currentSite.getCoords(), knownSite.getCoords());
			Vector remainder = m_Lattice.removeLattice(translation);
			if (remainder.length() > CartesianBasis.getPrecision()) { // A new dimension to the lattice
				Vector[] periodicVectors = m_Lattice.getPeriodicVectors();
				periodicVectors = (Vector[]) ArrayUtils.appendElement(periodicVectors, translation);
				m_Lattice = new BravaisLattice(periodicVectors);
			}
			return false;
		}

		@Override
		public String getDescription() {
			return "";
		}

		@Override
		public Vector[] getCellVectors() {
			return m_Lattice.getCellVectors();
		}

		@Override
		public boolean[] getVectorPeriodicity() {
			return m_Lattice.getDimensionPeriodicity();
		}

		@Override
		public int numDefiningSites() {
			return m_NumSites;
		}

		@Override
		public Coordinates getSiteCoords(int index) {
			return this.getSite(index).getCoords();
		}

		@Override
		public Species getSiteSpecies(int index) {
			return this.getSite(index).getSpecies();
		}

		/**
		 * REturns true if this ion does not repeat periodically in any direction, and
		 * false otherwise.
		 * 
		 * @return true if this ion does not repeat periodically in any direction, and
		 *         false otherwise.
		 */
		public boolean isAperiodic() {
			return (!ArrayUtils.arrayContains(getVectorPeriodicity(), true));
		}

		/**
		 * Returns the index'th site in the found ion.
		 * 
		 * @param index The index of the site to find.
		 * @return the index'th site in the found ion.
		 */
		private Site getSite(int index) {

			int currIndex = 0;
			for (Site site : m_Sites) {
				if (site == null) {
					continue;
				}
				if (currIndex == index) {
					return site;
				}
				currIndex++;
			}
			return null;
		}
	}

	/**
	 * Contains information about a found ion.
	 * 
	 * @author timmueller
	 *
	 */
	private class FoundIonData {

		private String m_Symbol;
		private Structure m_Structure;
		private GeneralStructureMapper.Map m_Map;

		/**
		 * Initialize an object with the given data.
		 * 
		 * @param ionTypeSymbol  The symbol of the ion type.
		 * @param foundStructure The sub-structure containing the sites in the found
		 *                       ion.
		 * @param map            A map between the sub-structure containing the sites of
		 *                       the found ion and the representative structure for this
		 *                       ion.
		 */
		private FoundIonData(String ionTypeSymbol, Structure foundStructure, GeneralStructureMapper.Map map) {
			m_Symbol = ionTypeSymbol;
			m_Structure = foundStructure;
			m_Map = map;
		}

		/**
		 * Returns the symbol of the ion type.
		 * 
		 * @return the symbol of the ion type.
		 */
		public String getIonStructureSymbol() {
			return m_Symbol;
		}

		/**
		 * Returns the sub-structure containing the sites in the found ion.
		 * 
		 * @return the sub-structure containing the sites in the found ion.
		 */
		public Structure getFoundStructure() {
			return m_Structure;
		}

		/**
		 * Returns a map between the sub-structure containing the sites of the found ion
		 * and the representative structure for this ion.
		 * 
		 * @returns a map between the sub-structure containing the sites of the found
		 *          ion and the representative structure for this ion.
		 */
		public GeneralStructureMapper.Map getMap() {
			return m_Map;
		}
	}
}
