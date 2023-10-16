package tri.oxidationstates.ion;

import java.util.ArrayList;

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
 * This class finds polyatomic ions that are networks of atoms in consisting of
 * oxygen atoms boudn to one of the following elements: boron, carbon, nitrogen,
 * aluminum, silicon, germanium, tin, lead, phosphorus, sulfur, arsenic,
 * selenium, chlorine, bromine, chromium, molybdenum, tungsten,
 * 
 * @author timmueller
 *
 */
public class OxideIonFinder implements IIonFinder {

	// The list of cations allowed in the oxide ions.
	private static Element[] CATION_ELEMENTS = new Element[] { Element.boron, Element.carbon, Element.nitrogen,
			Element.aluminum, Element.silicon, Element.germanium, Element.tin, Element.lead, Element.phosphorus,
			Element.sulfur, Element.arsenic, Element.selenium, Element.chlorine, Element.bromine, Element.chromium,
			Element.molybdenum, Element.tungsten, };

	private static double MAX_COVALENT_RADIUS = 0;

	static {
		for (Element element : CATION_ELEMENTS) {
			MAX_COVALENT_RADIUS = Math.max(MAX_COVALENT_RADIUS, element.getCovalentRadius());
		}
	}

	private Structure m_Structure;
	private ArrayList<Structure> m_FoundIons = new ArrayList<Structure>();

	/**
	 * Initialize the finder and find the oxide ions in the given strucutre.
	 * 
	 * @param structure The given structure.
	 */
	public OxideIonFinder(Structure structure) {

		m_Structure = structure;
		this.findIons();

	}

	/**
	 * Searches for the oxide ions in the given structure.
	 */
	private void findIons() {

		Element[] knownElements = m_Structure.getDistinctElements();
		if (!ArrayUtils.arrayContains(knownElements, Element.oxygen)) {
			return;
		}

		boolean[] allowedSites = new boolean[m_Structure.numDefiningSites()];
		for (int siteNum = 0; siteNum < allowedSites.length; siteNum++) {
			Element element = m_Structure.getSiteSpecies(siteNum).getElement();
			allowedSites[siteNum] = ArrayUtils.arrayContains(CATION_ELEMENTS, element) || (element == Element.oxygen);
		}

		for (int siteNum = 0; siteNum < allowedSites.length; siteNum++) {
			if (!allowedSites[siteNum]) {
				continue;
			}
			Structure.Site site = m_Structure.getDefiningSite(siteNum);
			IonBuilder builder = new IonBuilder(site, allowedSites);
			if (builder.numDefiningSites() <= 1) {
				continue;
			} // Don't allow monatomic ions
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
	 * This class starts from a given site and tries to build a polyatomic oxide ion
	 * that includes that site.
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
			this.addSites(initialSite, allowedSites);
		}

		/**
		 * Tries to add this site to the ion. If successful, recursively tries to add
		 * neighboring sites. Sites are considered to be neighbors if their bond
		 * distance is less than the sum of their covalent radii times a factor (1.1 by
		 * default).
		 * 
		 * @param currentSite  The site we are trying to add.
		 * @param allowedSites THe indices of this array are site indices. The values
		 *                     are "true" if the corresponding site is available to be
		 *                     added to tihs ion, and "false" otherwise.
		 * @return true of this site was added, false otherwise
		 */
		private boolean addSites(Site currentSite, boolean[] allowedSites) {

			Site knownSite = m_Sites[currentSite.getIndex()];
			boolean isOxygen = (currentSite.getSpecies().getElement() == Element.oxygen);
			if (knownSite == null) {
				if (!allowedSites[currentSite.getIndex()]) {
					return false;
				}
				m_Sites[currentSite.getIndex()] = currentSite;
				m_NumSites++;
				allowedSites[currentSite.getIndex()] = false;
				double currRadius = currentSite.getSpecies().getElement().getCovalentRadius();
				double maxNeighborRadius = isOxygen ? MAX_COVALENT_RADIUS : Element.oxygen.getCovalentRadius();
				double searchDistance = RADIUS_FACTOR * (currRadius + maxNeighborRadius);
				Structure.Site[] neighbors = m_Structure.getNearbySites(currentSite.getCoords(), searchDistance, false);
				for (Site neighbor : neighbors) {
					Element neighborElement = neighbor.getSpecies().getElement();
					if ((isOxygen) == (neighborElement == Element.oxygen)) {
						continue;
					}
					double neighborRadius = neighbor.getSpecies().getElement().getCovalentRadius();
					double maxAllowedDistance = RADIUS_FACTOR * (currRadius + neighborRadius);
					double distance = neighbor.distanceFrom(currentSite);
					if (distance <= maxAllowedDistance) {
						addSites(neighbor, allowedSites);
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
}
