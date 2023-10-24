package _global.tri.oxidationstates.ion;

import java.util.ArrayList;
import java.util.Set;

import _global.tri.oxidationstates.util.Composition;
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
 * Find networks of atoms that are polyatomic ions, where the composition of the
 * network(per unit cell) can be found in a list of compositions. If a network
 * can be extended to include additional atoms of a type already in the network,
 * it will not be returned.
 * 
 * @author timmueller
 *
 */
public class CompositionIonFinder implements IIonFinder {

	// Manually update this list with any compositions you want to include
	private static String[] COMPOSITIONS = new String[] { "O2", "CN", "S2O3", "S2O4", "S2O6", "C2O4", "OCN" };

	// For each composition, break it down into the distinct elements and how many
	// of each of them are in the composition
	private static Element[][] ELEMENTS;
	private static int[][] COUNTS;

	static {
		ELEMENTS = new Element[COMPOSITIONS.length][];
		COUNTS = new int[COMPOSITIONS.length][];
		for (int compositionNum = 0; compositionNum < COMPOSITIONS.length; compositionNum++) {
			Composition parser = new Composition(COMPOSITIONS[compositionNum]);
			Set<String> elementSymbols = parser.getSymbols();
			ELEMENTS[compositionNum] = new Element[elementSymbols.size()];
			COUNTS[compositionNum] = new int[elementSymbols.size()];
			int elementNum = 0;
			for (String symbol : elementSymbols) {
				ELEMENTS[compositionNum][elementNum] = Element.getElement(symbol);
				COUNTS[compositionNum][elementNum] = (int) Math.round(parser.getCount(symbol));
				elementNum++;
			}
		}
	}

	private Structure m_Structure;

	// All of the networks we have found so far with the required composition
	private ArrayList<Structure> m_FoundIons = new ArrayList<Structure>();

	/**
	 * Initialize the composition finder for the given structure and find the ions
	 * 
	 * @param structure The structure in which we will find the polyatomic ions
	 */
	public CompositionIonFinder(Structure structure) {

		m_Structure = structure;
		this.findIons();

	}

	/**
	 * Find the polyatomic ions for all compositions
	 */
	private void findIons() {
		for (int compositionNum = 0; compositionNum < COMPOSITIONS.length; compositionNum++) {
			this.findIons(compositionNum);
		}
	}

	/**
	 * Find the ions for a given composition
	 * 
	 * @param compositionNum The index of the given composition
	 */
	private void findIons(int compositionNum) {

		Element[] elements = ELEMENTS[compositionNum];
		int[] targetCounts = COUNTS[compositionNum];

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
			if (builder.numDefiningSites() <= 1) {
				continue;
			}
			if (!builder.matchesCounts(elements, targetCounts)) {
				continue;
			}
			m_FoundIons.add(new Structure(builder));
		}
	}

	@Override
	public int numFoundPolyatomicIons() {
		return m_FoundIons.size();
	}

	@Override
	public Structure getFoundPolyatomicIon(int index) {
		return m_FoundIons.get(index);
	}

	/**
	 * Finds networks of atoms with a given composition and extracts them as
	 * polyatomic ions
	 * 
	 * @author timmueller
	 *
	 */
	private class IonBuilder implements IStructureData {

		private static double RADIUS_FACTOR = 1.1;

		private Site[] m_Sites = new Site[m_Structure.numDefiningSites()];
		private int m_NumSites = 0;
		private BravaisLattice m_Lattice = new BravaisLattice(new Vector[0]);

		/**
		 * Find any ion that contains the given site
		 * 
		 * @param initialSite  The site from which we start to build the network of
		 *                     atoms
		 * @param allowedSites The indices of this array correspond to site indices in
		 *                     the structure. The value is true if the site is available
		 *                     for the network, and false otherwise.
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
		 * Return true if the ion found by this builder matches the composition given by
		 * elements and targetCounts
		 * 
		 * @param elements     The elements to be matched, in order
		 * @param targetCounts The target number of atoms for each element. The indices
		 *                     of this array correspond to the indices in the elements
		 *                     array.
		 * @return True if the found ion matches the given composition, false otherwise
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
		 * Try to add the current site to the network of sites found so far. This is a
		 * recursive method that will attempt to continue to grow the network by
		 * searching for neighbors (based on covalent radii) if the site is added.
		 * 
		 * @param currentSite       A site to be added to the network
		 * @param allowedSites      The indices of this array correspond to site indices
		 *                          in the structure. The value is true if the site is
		 *                          available for the network, and false otherwise.
		 * @param maxCovalentRadius The maximum possible covalent radius for a
		 *                          neighboring site that could be added to this
		 *                          network.
		 * @return true if this site is added, false if it isn't
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
		 * Returns the site in this ion with the given index
		 * 
		 * @param index The index of the site in this ion
		 * @return the site in this ion with the given index
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
}
