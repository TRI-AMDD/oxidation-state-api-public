package _global.tri.oxidationstates.calculator;

import java.util.Arrays;
import java.util.HashSet;

import matsci.Element;
import matsci.structure.Structure;
import matsci.structure.Structure.Site;
import matsci.util.arrays.ArrayUtils;

/**
 * This is a simplified bond valence calculator that should be easier to use for
 * this app than the one in the matsci oxidation analyzer. The basic code and
 * parameters are borrowed from that one.
 * 
 * @author timmueller
 *
 */
public class BondValenceCalculator {

	private static double MAX_BOND_LENGTH_PAD = 1.25;

	private Structure m_Structure;
	private double[] m_BVSums;
	private String m_FailureReason;

	/**
	 * Create a bond valence calculator for the given structure
	 * 
	 * @param structure The given structure
	 */
	public BondValenceCalculator(Structure structure) {
		m_Structure = structure;
		m_BVSums = calculateBVSums();
	}

	/**
	 * Returns true if some the bond valence sum for at least one site is zero,
	 * false otherwise
	 * 
	 * @return true if some the bond valence sum for at least one site is zero,
	 *         false otherwise
	 */
	public boolean hasZeroSums() {
		if (m_BVSums == null) {
			return false;
		}
		return ArrayUtils.arrayContains(m_BVSums, 0);
	}

	/**
	 * Returns the global instability index, based on the oxidations states of the
	 * provided structure and the calculated bond valence sums.
	 * 
	 * @return the global instability index, based on the oxidations states of the
	 *         provided structure and the calculated bond valence sums.
	 */
	public double getGlobalInstabilityIndex() {
		if (m_BVSums == null) {
			return Double.NaN;
		}
		double returnValue = 0;
		for (int siteNum = 0; siteNum < m_BVSums.length; siteNum++) {
			double bvSum = m_BVSums[siteNum];
			double oxidationState = m_Structure.getSiteSpecies(siteNum).getOxidationState();
			double delta = bvSum - oxidationState;
			returnValue += delta * delta;
		}
		return Math.sqrt(returnValue / m_BVSums.length);
	}

	/**
	 * Returns the calculated bond valence sums, with the array indices
	 * corresponding to the site indices in the provided structure. If bond valence
	 * sums could not be calculated, returns null.
	 * 
	 * @return the calculated bond valence sums, with the array indices
	 *         corresponding to the site indices in the provided structure. If bond
	 *         valence sums could not be calculated, returns null.
	 */
	public double[] getBondValenceSums() {
		return (m_BVSums == null) ? null : m_BVSums.clone();
	}

	/**
	 * Returns the reason sums couldn't be calculated, if calculating sums failed.
	 * 
	 * @return the reason sums couldn't be calculated, if calculating sums failed.
	 */
	public String getFailureReason() {
		return m_FailureReason;
	}

	/**
	 * Calculates bond valence sums. If bond valence sums could not be calculated,
	 * returns null.
	 * 
	 * @return the calculated bond valence sums, with the array indices
	 *         corresponding to the site indices in the provided structure. If bond
	 *         valence sums could not be calculated, returns null.
	 */
	private double[] calculateBVSums() {

		Element[] elements = m_Structure.getDistinctElements();
		if (elements.length < 2) {
			m_FailureReason = "Only one element in structure.";
			return null;
		}
		double[] rValues = new double[elements.length];
		for (int elementNum = 0; elementNum < elements.length; elementNum++) {
			Element element = elements[elementNum];
			int atomicNumber = element.getAtomicNumber();
			double rValue = (atomicNumber < 0) ? Double.NaN : R_PARAMETERS[atomicNumber];
			if (Double.isNaN(rValue)) {
				m_FailureReason = "Unknown element: " + element.getSymbol();
				return null;
			}
			rValues[elementNum] = rValue;
		}

		Arrays.sort(rValues);
		double maxSearchRadius = rValues[rValues.length - 1] + rValues[rValues.length - 2];
		maxSearchRadius += MAX_BOND_LENGTH_PAD; // For possible errors

		double[] bvSums = new double[m_Structure.numDefiningSites()];
		for (int siteNum = 0; siteNum < m_Structure.numDefiningSites(); siteNum++) {

			Site site = m_Structure.getDefiningSite(siteNum);
			Element element1 = site.getSpecies().getElement();
			double eNeg1 = element1.getElectronegativityAllredRochow();
			int atomicNumber1 = element1.getAtomicNumber();
			double r1 = R_PARAMETERS[atomicNumber1];
			double c1 = C_PARAMETERS[atomicNumber1];
			Site[] nearbySites = m_Structure.getNearbySites(site.getCoords(), maxSearchRadius, false);

			for (Site neighbor : nearbySites) {

				Element element2 = neighbor.getSpecies().getElement();
				int atomicNumber2 = element2.getAtomicNumber();
				if (atomicNumber1 == atomicNumber2) {
					continue;
				} // bond valence is zero

				double r2 = R_PARAMETERS[atomicNumber2];
				double c2 = C_PARAMETERS[atomicNumber2];

				// Check to see if the two sites are close enough
				double maxDistance = r1 + r2;
				maxDistance += MAX_BOND_LENGTH_PAD;
				double distance = neighbor.distanceFrom(site);

				if (distance > maxDistance) {
					continue;
				}

				boolean anion1 = BV_ELECTRONEGATIVE_ELEMENTS.contains(element1);
				boolean anion2 = BV_ELECTRONEGATIVE_ELEMENTS.contains(element2);
				if (!anion1 && !anion2) {
					continue;
				}

				Element anion = anion1 ? element1 : element2;

				if (anion1 && anion2) {
					double eNeg2 = element2.getElectronegativityAllredRochow();
					if (Double.isNaN(eNeg1) || Double.isNaN(eNeg2) || (eNeg1 == eNeg2)) {
						double paulingENeg1 = element1.getElectronegativityPaulingCRC();
						double paulingEneg2 = element2.getElectronegativityPaulingCRC();
						anion = (paulingENeg1 < paulingEneg2) ? element2 : element1;
					} else {
						anion = (eNeg1 < eNeg2) ? element2 : element1;
					}
				}

				double cTerm = Math.sqrt(c1) - Math.sqrt(c2);
				cTerm *= cTerm;
				double bigR = r1 + r2 - (r1 * r2 * cTerm) / ((c1 * r1) + (c2 * r2));
				double bondValence = Math.exp((bigR - distance) / 0.37);

				double sign = (anion == element2) ? 1 : -1;

				bvSums[site.getIndex()] += bondValence * sign / 2; // Prevent double-counting
				bvSums[neighbor.getIndex()] -= bondValence * sign / 2; // Prevent double-counting
			}
		}
		return bvSums;
	}

	/**
	 * Return true if this calculator has bond valence parameters or the given
	 * element, false otherwise.
	 * 
	 * @param element The given element
	 * @return true if this calculator has bond valence parameters or the given
	 *         element, false otherwise.
	 */
	public static boolean hasParametersForElement(Element element) {
		return !Double.isNaN(C_PARAMETERS[element.getAtomicNumber()]);
	}

	// From J. Am. Chem. Soc. 1991, 113, 9, 3226–3229
	private static double[] C_PARAMETERS = new double[] { 0, 0.89, Double.NaN, 0.97, 1.47, 1.6, 2.0, 2.61, 3.15, 3.98,
			Double.NaN, 1.01, 1.23, 1.47, 1.58, 1.96, 2.35, 2.74, Double.NaN, 0.91, 1.04, 1.2, 1.32, 1.45, 1.56, 1.6,
			1.64, 1.7, 1.75, 1.75, 1.66, 1.82, 1.51, 2.23, 2.51, 2.58, Double.NaN, 0.89, 0.99, 1.11, 1.22, 1.23, 1.3,
			Double.NaN, 1.42, 1.54, 1.35, 1.42, 1.46, 1.49, 1.72, 1.72, 2.72, 2.38, Double.NaN, 0.86, 0.97, 1.08, 1.08,
			1.07, 1.07, Double.NaN, 1.07, 1.01, 1.11, 1.1, 1.1, 1.1, 1.11, 1.11, 1.06, 1.14, 1.23, 1.33, 1.4, 1.46,
			Double.NaN, 1.55, Double.NaN, Double.NaN, 1.44, 1.44, 1.55, 1.67, Double.NaN, Double.NaN, Double.NaN,
			Double.NaN, Double.NaN, Double.NaN, 1.11, Double.NaN, 1.22, };

	private static double[] R_PARAMETERS = new double[] { 0, 0.38, Double.NaN, 1.0, 0.81, 0.79, 0.78, 0.72, 0.63, 0.58,
			Double.NaN, 1.36, 1.21, 1.13, 1.12, 1.09, 1.03, 0.99, Double.NaN, 1.73, 1.5, 1.34, 1.27, 1.21, 1.16, 1.17,
			1.16, 1.09, 1.04, 0.87, 1.07, 1.14, 1.21, 1.21, 1.18, 1.13, Double.NaN, 1.84, 1.66, 1.52, 1.43, 1.4, 1.37,
			Double.NaN, 1.21, 1.18, 1.11, 1.12, 1.28, 1.34, 1.37, 1.41, 1.4, 1.33, Double.NaN, 2.05, 1.88, 1.71, 1.68,
			1.66, 1.64, Double.NaN, 1.61, 1.62, 1.58, 1.56, 1.54, 1.53, 1.51, 1.5, 1.49, 1.47, 1.42, 1.39, 1.38, 1.37,
			Double.NaN, 1.37, Double.NaN, Double.NaN, 1.32, 1.62, 1.53, 1.54, Double.NaN, Double.NaN, Double.NaN,
			Double.NaN, Double.NaN, Double.NaN, 1.7, Double.NaN, 1.59, };

	private static HashSet<Element> BV_ELECTRONEGATIVE_ELEMENTS = new HashSet<Element>();

	static {
		// Make sure we have R and C values for each of the known elements
		int oldCLength = C_PARAMETERS.length;
		C_PARAMETERS = ArrayUtils.growArray(C_PARAMETERS, Element.numKnownElements() + 1 - oldCLength);
		for (int paramNum = oldCLength; paramNum < C_PARAMETERS.length; paramNum++) {
			C_PARAMETERS[paramNum] = Double.NaN;
		}

		int oldRLength = R_PARAMETERS.length;
		R_PARAMETERS = ArrayUtils.growArray(R_PARAMETERS, Element.numKnownElements() + 1 - oldRLength);
		for (int paramNum = oldRLength; paramNum < R_PARAMETERS.length; paramNum++) {
			R_PARAMETERS[paramNum] = Double.NaN;
		}

		// Create the set of known electronegative elements
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("H"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("B"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("C"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("Si"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("N"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("P"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("As"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("Sb"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("O"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("S"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("Se"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("Te"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("F"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("Cl"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("Br"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("I"));
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("Ge")); // This is in table 1 but not listed in the main
																	// text. Main text is likely a typo, as they say
																	// there should be 17 elements.

		// These were not included in the original paper, but are potential anions in
		// our data set
		BV_ELECTRONEGATIVE_ELEMENTS.add(Element.getElement("Bi"));
	}
}
