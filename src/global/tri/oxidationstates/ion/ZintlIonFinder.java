package global.tri.oxidationstates.ion;

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
 * Finds possible Zintl ions in the given structure, where a possible Zintl ion
 * is defined as a covalently bound network of ions containing one or more of
 * boron, carbon, silicon, phosphorus, sulfur, gallium, germanium, arsenic,
 * selenium, indium, tin, antimony, tellurium, lead, bismuth, polonium.
 * 
 * @author timmueller
 *
 */
public class ZintlIonFinder implements IIonFinder {

	// The elements allowed in possible Zintl ions
	private static Element[] ZINTL_ELEMENTS = new Element[] { Element.boron, Element.carbon, Element.silicon,
			Element.phosphorus, Element.sulfur, Element.gallium, Element.germanium, Element.arsenic, Element.selenium,
			Element.indium, Element.tin, Element.antimony, Element.tellurium, Element.lead, Element.bismuth,
			Element.polonium };

	private static double MAX_COVALENT_RADIUS = 0;
	static {
		for (Element element : ZINTL_ELEMENTS) {
			MAX_COVALENT_RADIUS = Math.max(MAX_COVALENT_RADIUS, element.getCovalentRadius());
		}
	}

	private Structure m_Structure;
	private ArrayList<Structure> m_FoundIons = new ArrayList<Structure>();

	/**
	 * Initializes a finder to search for ions in the given structure.
	 * 
	 * @param structure The given structure.
	 */
	public ZintlIonFinder(Structure structure) {

		m_Structure = structure;
		this.findIons();

	}

	/**
	 * Find candidate Zintl ions in this structure
	 */
	private void findIons() {

		boolean[] allowedSites = new boolean[m_Structure.numDefiningSites()];
		for (int siteNum = 0; siteNum < allowedSites.length; siteNum++) {
			Element element = m_Structure.getSiteSpecies(siteNum).getElement();
			allowedSites[siteNum] = ArrayUtils.arrayContains(ZINTL_ELEMENTS, element);
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
	 * This class starts from a given site and tries to build a polyatomic potential
	 * Zintl ion that includes that site.
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
			if (knownSite == null) {
				if (!allowedSites[currentSite.getIndex()]) {
					return false;
				}
				m_Sites[currentSite.getIndex()] = currentSite;
				m_NumSites++;
				allowedSites[currentSite.getIndex()] = false;
				double currRadius = currentSite.getSpecies().getElement().getCovalentRadius();
				double searchDistance = RADIUS_FACTOR * (currRadius + MAX_COVALENT_RADIUS);
				Structure.Site[] neighbors = m_Structure.getNearbySites(currentSite.getCoords(), searchDistance, false);
				for (Site neighbor : neighbors) {
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
