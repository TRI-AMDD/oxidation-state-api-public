package _global.tri.oxidationstates.fitting;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.Reader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import _global.tri.oxidationstates.calculator.likelihood.LikelihoodCalculator;
import _global.tri.oxidationstates.ion.IonFactory;
import _global.tri.oxidationstates.ion.ZintlIonFinder;
import _global.tri.oxidationstates.ion.IonFactory.Ion;
import _global.tri.oxidationstates.util.Composition;
import matsci.Species;
import matsci.io.app.log.Status;
import matsci.io.clusterexpansion.PRIM;
import matsci.io.structure.CIF;
import matsci.io.vasp.POSCAR;
import matsci.structure.PartiallyOccupiedStructure;
import matsci.structure.Structure;
import matsci.util.MSMath;
import matsci.util.arrays.ArrayUtils;

/**
 * This is the main class for data sets (e.g. testing, training data).
 * 
 * @author timmueller
 *
 */
public class OxidationStateData {

	private ArrayList<Entry> m_Entries = new ArrayList<Entry>();
	private String m_StructDir;

	/**
	 * Create a data set with the provided entries
	 * 
	 * @param entries   The entries to included in this data set
	 * @param structDir A directory that contains structure files for each in the
	 *                  entries, in VASP POSCAR format
	 */
	public OxidationStateData(Collection<Entry> entries, String structDir) {
		m_Entries = new ArrayList<Entry>(entries);
		m_StructDir = structDir;
		Status.basic(this.numEntries() + " entries in initial data set.");
	}

	/**
	 * Read a date set from a given file
	 * 
	 * @param fileName  The name of the given file
	 * @param structDir A directory that contains structure files for each in the
	 *                  entries, in VASP POSCAR format
	 */
	public OxidationStateData(String fileName, String structDir) {
		m_StructDir = structDir;
		readFile(fileName);
		Status.basic(this.numEntries() + " entries in initial data set.");
	}

	/**
	 * Read a data set from a given file and remove entries according to the given
	 * options
	 * 
	 * @param fileName            The name of the given file
	 * @param removeNonInteger    Remove all entries that contain oxidation states
	 *                            with non-integer values
	 * @param hullCutoff          An energy in eV / atom. All entries with energies
	 *                            above the convex hull above this value will be
	 *                            removed.
	 * @param removeZeroOxidation Remove entries with oxidaiton states of zero
	 * @param removeZintl         Remove entries that contain Zintl ions, as
	 *                            determined by the {@link ZintlIonFinder}
	 * @param structDir           A directory that contains structure files for each
	 *                            in the entries, in VASP POSCAR format
	 */
	public OxidationStateData(String fileName, boolean removeNonInteger, double hullCutoff, boolean removeZeroOxidation,
			boolean removeZintl, String structDir) {

		this(fileName, structDir);

		if (removeNonInteger) {
			this.removeEntriesWithNonIntegerStates();
		}

		if (hullCutoff >= 0) {
			this.removeUnstableEntries(hullCutoff);
		}

		if (removeZeroOxidation) {
			this.removeEntriesWithZeroOxidationStates();
		}

		if (removeZintl) {
			this.removeEntriesWithZintlIons();
		}
	}

	/**
	 * Read data from a file with the given name
	 * 
	 * @param fileName The name of the data file
	 */
	private void readFile(String fileName) {

		try {
			FileReader reader = new FileReader(fileName);
			readFile(reader);
			reader.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Read data from the given Reader
	 * 
	 * @param reader The given reader
	 * @throws IOException if an IO error is encountered while reading the file.
	 */
	private void readFile(Reader reader) throws IOException {

		LineNumberReader lineReader = new LineNumberReader(reader);
		String line = lineReader.readLine();
		while (line != null && line.trim().length() > 0) {
			String[] fields = line.split("\t");
			String[] idFields = fields[0].split(" ");
			String id = (idFields.length == 1) ? idFields[0] : idFields[1];
			String composition = fields[1];
			String[] ionSymbols = fields[2].trim().split(" ");
			String[] sources = fields[3].trim().split(" ");
			double energyAboveHull = (fields.length > 4) ? Double.parseDouble(fields[4]) : Double.NaN;
			double gii = (fields.length > 5) ? Double.parseDouble(fields[5]) : Double.NaN;

			Ion[] ions = new Ion[ionSymbols.length];

			// TODO consider whether this is really necessary
			boolean knownIons = true;
			for (int ionNum = 0; ionNum < ions.length; ionNum++) {
				ions[ionNum] = IonFactory.get(ionSymbols[ionNum]);
				knownIons &= (ions[ionNum] != null);
			}

			if (knownIons) {
				Entry entry = new Entry(id, composition, ions, sources, energyAboveHull, gii);
				m_Entries.add(entry);
			} else {
				String output = "Unknown ion types in entry " + id + ": ";
				for (String ionSymbol : ionSymbols) {
					output += ionSymbol + ", ";
				}
				Status.warning(output);
			}

			line = lineReader.readLine();
		}
	}

	/**
	 * Writes a file containing this data set
	 * 
	 * @param fileName The name of the file to be written
	 */
	public void writeFile(String fileName) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
			this.writeFile(writer);
			writer.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Writes a file containing this data set
	 * 
	 * @param writer The file will be written to this writer
	 * @throws IOException if there is an I/O error
	 */
	public void writeFile(Writer writer) throws IOException {

		for (Entry entry : m_Entries) {
			String output = "";
			output = output + entry.getID() + "\t";
			output = output + entry.getGivenCompositionString() + "\t";
			Ion[] ions = entry.getAllIons();
			String[] symbols = new String[ions.length];
			for (int ionNum = 0; ionNum < ions.length; ionNum++) {
				symbols[ionNum] = ions[ionNum].getSymbol();
			}
			int[] map = ArrayUtils.getSortPermutation(symbols);

			for (int specNum = 0; specNum < ions.length; specNum++) {
				output = output + ions[map[specNum]] + " ";
			}
			output = output + "\t";
			String[] sources = entry.getAllSources();
			Arrays.sort(sources);
			for (String source : sources) {
				output = output + source + " ";
			}
			output = output + "\t";
			output = output + entry.getEnergyAboveHull() + "\t";
			output = output + entry.getGlobalInstabilityIndex();
			writer.write(output + "\n");
		}
		writer.flush();

	}

	/**
	 * Removes a random subset of this data set
	 * 
	 * @param percentToRemove The percent of entries to remove (rounded off).
	 */
	public void removeRandomEntries(double percentToRemove) {
		ArrayList<Entry> newEntryList = new ArrayList<Entry>();
		int numToKeep = (int) Math.round((1.0 - percentToRemove) * m_Entries.size());
		int[] randomIndexOrder = MSMath.getRandomShuffle(m_Entries.size());
		for (int entryNum = 0; entryNum < numToKeep; entryNum++) {
			newEntryList.add(m_Entries.get(randomIndexOrder[entryNum]));
		}

		m_Entries = newEntryList;
		Status.basic(this.numEntries() + " remaining entries in data set after " + (percentToRemove * 100)
				+ " percent of entries were randomly removed.");

	}

	/**
	 * Returns a copy of this data set
	 * 
	 * @return a copy of this data set
	 */
	public OxidationStateData copy() {
		return new OxidationStateData(m_Entries, m_StructDir);
	}

	/**
	 * Removes the entries in the given set. Note that the entry objects need to be
	 * exactly the same; i.e. both this data set and the "entriesToRemove" data set
	 * should be derived from the some data set.
	 * 
	 * @param entriesToRemove A data set containing the entries to be removed. Note that the entry objects need to be
	 * exactly the same; i.e. both this data set and the "entriesToRemove" data set
	 * should be derived from the some data set.
	 */
	public void removeEntries(OxidationStateData entriesToRemove) {
		HashSet<Entry> entrySet = new HashSet<Entry>(m_Entries);
		for (int entryNum = 0; entryNum < entriesToRemove.numEntries(); entryNum++) {
			Entry entry = entriesToRemove.getEntry(entryNum);
			entrySet.remove(entry);
		}
		m_Entries = new ArrayList<Entry>(entrySet);
		Status.basic(
				this.numEntries() + " remaining entries in data set after entries in provided data set were removed.");

	}

	/**
	 * Only keep entries for which all ions are in the given set
	 * 
	 * @param allowedIons Entries will only be kept if all ions in the entry are in
	 *                    this set.
	 */
	public void dataKeepOnlyIons(Set<Ion> allowedIons) {

		HashSet<Entry> entrySet = new HashSet<Entry>();
		for (Entry entry : m_Entries) {
			Ion[] ions = entry.getAllIons();
			boolean allowed = true;
			for (Ion ion : ions) {
				allowed &= allowedIons.contains(ion);
			}
			if (allowed) {
				entrySet.add(entry);
			}
		}
		m_Entries = new ArrayList<Entry>(entrySet);
		Status.basic(this.numEntries()
				+ " remaining entries in data set after keeping only entries that had all ions in the following set: "
				+ allowedIons);
	}

	/**
	 * Randomly split the data into numSplits test sets. The union of all of the
	 * tests sets will be this complete data set, and all test sets will be
	 * approximately the same size. The split is done so that no composition will
	 * appear in more than one test set, so there is never the same composition in a
	 * test and training set.
	 * 
	 * @param numSplits The number of test sets to generate
	 * @return An array of generated test sets
	 */
	public OxidationStateData[] splitData(int numSplits) {

		HashSet<String>[] compositionSplits = new HashSet[numSplits];
		ArrayList<Entry>[] splits = new ArrayList[numSplits];
		for (int splitNum = 0; splitNum < compositionSplits.length; splitNum++) {
			compositionSplits[splitNum] = new HashSet<String>();
			splits[splitNum] = new ArrayList<Entry>();
		}

		int[] permutation = MSMath.getRandomShuffle(m_Entries.size());
		int splitIndex = 0;
		for (int entryNum : permutation) {
			Entry entry = m_Entries.get(entryNum);
			String compositionString = entry.getComposition().getElementalComposition().getReducedComposition()
					.getStandardizedCompositionString();
			boolean foundMatch = false;
			for (int splitNum = 0; splitNum < splits.length; splitNum++) {
				if (!compositionSplits[splitNum].contains(compositionString)) {
					continue;
				}
				splits[splitNum].add(entry);
				foundMatch = true;
				break;
			}
			if (foundMatch) {
				continue;
			}
			compositionSplits[splitIndex].add(compositionString);
			splits[splitIndex++].add(entry);
			splitIndex %= numSplits;
		}

		OxidationStateData[] returnArray = new OxidationStateData[numSplits];
		for (int splitNum = 0; splitNum < returnArray.length; splitNum++) {
			returnArray[splitNum] = new OxidationStateData(splits[splitNum], m_StructDir);
		}

		return returnArray;

	}

	/**
	 * Returns the directory with atomic structure files for the entries
	 * 
	 * @return the directory with atomic structure files for the entries
	 */
	public String getStructDir() {
		return m_StructDir;
	}

	/**
	 * Returns a map in which the keys are the ions contained in this data set and
	 * the values are the number of entries that contain the corresponding ion.
	 * 
	 * @return a map in which the keys are the ions contained in this data set and
	 *         the values are the number of entries that contain the corresponding
	 *         ion.
	 */
	public HashMap<Ion, Integer> getCountsByIon() {

		HashMap<Ion, Integer> counts = new HashMap<Ion, Integer>();

		for (Entry entry : m_Entries) {
			for (int ionNum = 0; ionNum < entry.numIons(); ionNum++) {
				Ion ion = entry.getIon(ionNum);
				int count = counts.containsKey(ion) ? counts.get(ion) : 0;
				counts.put(ion, count + 1);
			}
		}
		return counts;
	}

	/**
	 * Removes all entries that contain a rate ions, where "rare" ions are those
	 * that appear in fewer than minAllowedCount entries
	 * 
	 * @param minAllowedCount The minimum number of entries an ion must appear in to
	 *                        not be considered rare.
	 */
	public void removeUncommonIonsByCount(int minAllowedCount) {

		ArrayList<Entry> newEntryList = new ArrayList<Entry>();
		HashMap<Ion, Integer> counts = getCountsByIon();

		HashSet<Ion> forbiddenIons = new HashSet<Ion>();
		for (Ion symbol : counts.keySet()) {
			int count = counts.get(symbol);
			if (count < minAllowedCount) {
				forbiddenIons.add(symbol);
			}
		}

		for (Entry entry : m_Entries) {
			boolean isAllowed = true;
			for (int ionNum = 0; ionNum < entry.numIons(); ionNum++) {
				Ion ion = entry.getIon(ionNum);
				if (forbiddenIons.contains(ion)) {// .getSymbol())) {
					isAllowed = false;
					break;
				}
			}
			if (isAllowed) {
				newEntryList.add(entry);
			}
		}

		m_Entries = newEntryList;
		Status.basic(this.numEntries()
				+ " remaining entries in data set after removing entries with ions with counts that were less than "
				+ minAllowedCount + ": " + forbiddenIons);
	}

	/**
	 * Prints to standard output the number of entries containing each ion in this
	 * data set.
	 */
	public void printNumEntriesByIon() {

		TreeMap<String, Integer> counts = new TreeMap<String, Integer>();

		for (Entry entry : m_Entries) {
			for (int ionNum = 0; ionNum < entry.numIons(); ionNum++) {
				Ion ion = entry.getIon(ionNum);
				String symbol = ion.getSymbol();
				int count = counts.containsKey(symbol) ? counts.get(symbol) : 0;
				counts.put(symbol, count + 1);
			}
		}

		for (String symbol : counts.keySet()) {
			int count = counts.get(symbol);
			System.out.println(symbol + "\t" + count);
		}

	}

	/**
	 * Removes entries containing rare ions, where an ion is rare if the fraction of
	 * entries it appears in for its ion type is less than minAllowedFraction
	 * 
	 * @param minAllowedFraction An ion will be considered rare if the fraction of
	 *                           entries it appears in for its ion type is less than
	 *                           this value. For example, if A2+ appears in 10
	 *                           entries and A3+ appears in 90 entries, then all
	 *                           entries containing A2+ will be removed if
	 *                           minAllowedFraction is less than 0.1.
	 */
	public void removeUncommonOxidationStates(double minAllowedFraction) {

		ArrayList<Entry> newEntryList = new ArrayList<Entry>();
		HashMap<String, Integer> ionCounts = new HashMap<String, Integer>();
		HashMap<String, Integer> ionTypeCounts = new HashMap<String, Integer>();

		for (Entry entry : m_Entries) {
			for (int ionNum = 0; ionNum < entry.numIons(); ionNum++) {
				Ion ion = entry.getIon(ionNum);
				String symbol = ion.getSymbol();
				int count = ionCounts.containsKey(symbol) ? ionCounts.get(symbol) : 0;
				ionCounts.put(symbol, count + 1);

				String typeSymbol = ion.getIonType().getSymbol();
				int typeCount = ionTypeCounts.containsKey(typeSymbol) ? ionTypeCounts.get(typeSymbol) : 0;
				ionTypeCounts.put(typeSymbol, typeCount + 1);
			}
		}

		HashSet<String> forbiddenIons = new HashSet<String>();
		for (String typeSymbol : ionTypeCounts.keySet()) {
			int numTotalIons = ionTypeCounts.get(typeSymbol);
			int minAllowedCount = (int) Math.ceil(numTotalIons * minAllowedFraction);
			for (String symbol : ionCounts.keySet()) {
				if (!IonFactory.getIonTypeSymbol(symbol).equals(typeSymbol)) {
					continue;
				}
				int count = ionCounts.get(symbol);
				if (count < minAllowedCount) {
					forbiddenIons.add(symbol);
				}
			}
		}

		for (Entry entry : m_Entries) {
			boolean isAllowed = true;
			for (int ionNum = 0; ionNum < entry.numIons(); ionNum++) {
				Ion ion = entry.getIon(ionNum);
				if (forbiddenIons.contains(ion.getSymbol())) {
					isAllowed = false;
					break;
				}
			}
			if (isAllowed) {
				newEntryList.add(entry);
			}
		}

		m_Entries = newEntryList;
		Status.basic(this.numEntries()
				+ " remaining entries in data set after removing entries with ions with counts that were less than "
				+ (minAllowedFraction * 100) + "% of ions of each type: " + forbiddenIons);
	}

	/**
	 * When entries with compositions written in terms of polyatomic ions are added
	 * to the data set, there will be two entries with the same ID: one with a
	 * composition written in terms of monatomic ions, and one with composition
	 * written in terms of polyatomic ions. This method removes of the two entries
	 * with the same ID. This method does not change the data set, but returns a map
	 * of the remaining entries keyed by entry ID.
	 * 
	 * @param keepPolyIons If true, remove the entries with duplicate ID that have
	 *                     monatomic ions. If false, remove the entries with
	 *                     duplicate ID that have polyatomic ions.
	 * @return A map of the remaining entries keyed by entry ID.
	 */
	public HashMap<String, Entry> getUniqueEntries(boolean keepPolyIons) {
		HashMap<String, Entry> keepers = new HashMap<String, Entry>();
		HashSet<String> idsToRemove = new HashSet<String>();
		for (Entry entry : m_Entries) {
			if (entry.getGivenCompositionString().contains(")")) { // Possible polyIon
				idsToRemove.add(entry.getID()); // We need to remove one of the two entries with this ID
			}
		}

		for (Entry entry : m_Entries) {

			// Figure out if this is an entry we should remove.
			if (idsToRemove.contains(entry.getID())
					&& (keepPolyIons == !entry.getGivenCompositionString().contains(")"))) {
				continue; // It is.
			}
			if (keepers.containsKey(entry.getID())) {
				throw new RuntimeException("Cannot determine which entry with ID " + entry.getID() + " to keep.");
			}
			keepers.put(entry.getID(), entry);
		}
		return keepers;
	}

	/**
	 * Removes all entries for which at least one of the ions has an oxidation state
	 * of zero.
	 */
	public void removeEntriesWithZeroOxidationStates() {
		ArrayList<Entry> newEntryList = new ArrayList<Entry>();
		for (Entry entry : m_Entries) {
			int numSpecies = entry.numIons();
			boolean allNonZero = true;
			for (int ionNum = 0; ionNum < numSpecies; ionNum++) {
				double oxidationState = entry.getIon(ionNum).getOxidationState();
				if (Math.abs(oxidationState) < 1E-2) {
					allNonZero = false;
					break;
				}
			}
			if (allNonZero) {
				newEntryList.add(entry);
			}
		}

		m_Entries = newEntryList;
		Status.basic(this.numEntries() + " remaining entries in data set with nonzero oxidation state.");

	}

	/**
	 * Removes all entries for which the GII in structDirectory is less than the GII
	 * in refStructDirectory. A tolerance of 1E-6 is used when comparing GIIs.
	 * 
	 * TODO re-write this method so that just reads the GII from the entry (that
	 * field wasn't there when this was written).
	 * 
	 * @param structDirectory    A directory containing structures, where the
	 *                           description field gives the GII.
	 * @param refStructDirectory A directory containing structures, where the
	 *                           description field gives the GII.
	 * @param calculator         A likelihood calculator used for logging purposes
	 *                           (tracking the likelihood score of the removed
	 *                           entries).
	 */
	public void removeGIIDecrease(String structDirectory, String refStructDirectory, LikelihoodCalculator calculator) {

		HashSet<String> idsToRemove = new HashSet<String>();

		for (Entry entry : m_Entries) {

			String fileName = structDirectory + "/" + entry.getID() + ".cif";
			String refFileName = refStructDirectory + "/" + entry.getID() + ".cif";

			String polyFileName = fileName.replace(".cif", "_poly.cif");
			String polyRefFileName = refFileName.replace(".cif", "_poly.cif");

			if (!entry.getGivenCompositionString().contains(")")) {

				// If the poly model is available, use that as it's more accurate.
				if (new File(polyFileName).exists()) {
					continue;
				}
				boolean remove = this.cleanGIIEntry(entry, calculator, fileName, refFileName);
				if (remove) {
					idsToRemove.add(entry.getID());
				}
			} else {

				boolean remove = this.cleanGIIEntry(entry, calculator, polyFileName, polyRefFileName);
				if (remove) {
					idsToRemove.add(entry.getID());
				}
			}
		}

		ArrayList<Entry> newEntryList = new ArrayList<Entry>();
		for (Entry entry : m_Entries) {
			if (idsToRemove.contains(entry.getID())) {
				continue;
			}
			newEntryList.add(entry);
		}

		m_Entries = newEntryList;
		Status.basic("Removed " + idsToRemove.size() + " unique IDs.");
		Status.basic(this.numEntries() + " remaining entries in data set that do not have lower GII in "
				+ structDirectory + " than in " + refStructDirectory);

	}

	/**
	 * Determine whether an entry should be removed based on the GII criterion. Note
	 * that the logging in this method can be a little misleading, as it calculated
	 * the likelihood score (but not the reference likelihood score) using only
	 * monatomic ions.
	 * 
	 * @param entry       The entry to be evaluated
	 * @param calculator  A likelihood calculator to be used for loggin purposes.
	 * @param fileName    The name of a file containing the atomic structure. The
	 *                    description field of this structure should have the file
	 *                    name.
	 * @param refFileName The name of a file containing the atomic structure. The
	 *                    description field of this structure should have the file
	 *                    name.
	 * @return True if the GII for the structure represented by fileName is less
	 *         than the GII for the structure represented by refFileName. A
	 *         tolerance of 1E-6 is used when comparing GIIs.
	 */
	private boolean cleanGIIEntry(Entry entry, LikelihoodCalculator calculator, String fileName, String refFileName) {
		if (!new File(fileName).exists()) {
			return false;
		}
		if (!new File(refFileName).exists()) {
			return false;
		}

		CIF infile = new CIF(fileName);
		CIF refInfile = new CIF(refFileName);

		double gii = Double.parseDouble(infile.getDescription());
		double refGII = Double.parseDouble(refInfile.getDescription());

		double delta = gii - refGII;

		if (delta > -1E-6) {
			return false;
		}

		PartiallyOccupiedStructure refStructure = new PartiallyOccupiedStructure(refInfile);
		PartiallyOccupiedStructure structure = new PartiallyOccupiedStructure(infile);

		String message = "Removing id " + entry.getID() + " with fileName " + new File(fileName).getName()
				+ " and composition " + structure.getCompositionString(" ") + " with delta " + delta
				+ " | Reference species: ";
		Species[] refSpecies = refStructure.getAllAllowedSpecies();
		for (Species species : refSpecies) {
			message += species + ". ";
		}
		message += "| Species: ";
		Species[] allSpecies = structure.getAllAllowedSpecies();
		for (Species species : allSpecies) {
			message += species + ". ";
		}

		double refLikelihood = calculator.optimizeLikelihood(entry).getMaxLikelihood();
		double likelihood = calculator.optimizeLikelihood(allSpecies).getMaxLikelihood();

		message += "| Ref likelihood: " + refLikelihood + " | likelihood: " + likelihood;
		message += " | EAH: " + entry.getEnergyAboveHull();
		Status.detail(message);
		return true;

	}

	/**
	 * Removes all entries that do not have charge neutral structures, defined as
	 * structures for which all of the oxidation states of the atoms in each unit
	 * cell add up to zero.
	 * 
	 * @param structDirName The name of the directory containing the structure files
	 *                      in VASP POSCAR format.
	 */
	public void removeNonChargeBalancedStructures(String structDirName) {

		ArrayList<Entry> newEntryList = new ArrayList<Entry>();
		for (Entry entry : m_Entries) {
			String fileName = structDirName + entry.getID() + ".vasp";

			PRIM infile = new PRIM(fileName);
			Structure structure = new Structure(infile);

			for (int siteNum = 0; siteNum < structure.numDefiningSites(); siteNum++) {
				structure.getDefiningSite(siteNum).setSpecies(infile.getAllowedSpecies()[siteNum][0]);
			}
			if (structure.isChargeBalanced()) {
				newEntryList.add(entry);
			}
		}
		m_Entries = newEntryList;
		Status.basic(this.numEntries() + " remaining entries in data set that are charge balanced.");

	}

	/**
	 * Removes all entries with ZintlIons, as determined by the
	 * {@link ZintlIonFinder}.
	 */
	public void removeEntriesWithZintlIons() {

		ArrayList<Entry> newEntryList = new ArrayList<Entry>();
		for (Entry entry : m_Entries) {

			Structure structure = entry.getStructure();
			ZintlIonFinder ionFinder = new ZintlIonFinder(structure);
			if (ionFinder.numFoundPolyatomicIons() == 0) {
				newEntryList.add(entry);
			}
		}
		m_Entries = newEntryList;
		Status.basic(this.numEntries() + " remaining entries in data set with no Zintl ions.");
	}

	/**
	 * Returns a map of oxidation states for each ion type in this data set. The map
	 * is keyed by the ion type ID and the values are the oxidation states, in
	 * ascending order.
	 * 
	 * @return a map of oxidation states for each ion type in this data set. The map
	 *         is keyed by the ion type ID and the values are the oxidation states,
	 *         in ascending order.
	 */
	public HashMap<String, int[]> getKnownOxidationStates() {

		// Get the set of all oxidation states.
		HashMap<String, TreeSet<Integer>> oxidationStates = new HashMap();

		int numEntries = this.numEntries();
		for (int entryNum = 0; entryNum < numEntries; entryNum++) {
			Entry entry = this.getEntry(entryNum);
			for (int ionNum = 0; ionNum < entry.numIons(); ionNum++) {
				Ion ion = entry.getIon(ionNum);
				double oxidationState = ion.getOxidationState();
				int intOxidationState = (int) Math.round(oxidationState);
				if (Math.abs(intOxidationState - oxidationState) > 1E-2) {
					continue;
				}

				// The tree set sorts the oxidation states
				TreeSet<Integer> knownOxidationStates = oxidationStates.get(ion.getIonType().getSymbol());
				if (knownOxidationStates == null) {
					knownOxidationStates = new TreeSet<Integer>();
					oxidationStates.put(ion.getIonType().getSymbol(), knownOxidationStates);
				}
				if (!knownOxidationStates.contains(intOxidationState)) {
					knownOxidationStates.add(intOxidationState);
				}
			}
		}

		HashMap<String, int[]> returnMap = new HashMap();
		for (String moleculeID : oxidationStates.keySet()) {
			int[] intArray = oxidationStates.get(moleculeID).stream().mapToInt(i -> i).toArray();
			returnMap.put(moleculeID, intArray);
		}
		return returnMap;
	}

	/**
	 * Removes all entries that contain oxidation states that are not within 0.01 of
	 * an integer.
	 */
	public void removeEntriesWithNonIntegerStates() {
		this.removeEntriesWithNonIntegerStates(1E-2);
	}

	/**
	 * Removes all entries that contain oxidation states that are not within
	 * "tolerance" of an integer.
	 * 
	 * @param tolerance The maximum allowed difference between the oxidation state
	 *                  and an integer to be considered an integer oxidation state.
	 */
	public void removeEntriesWithNonIntegerStates(double tolerance) {
		ArrayList<Entry> newEntryList = new ArrayList<Entry>();
		for (Entry entry : m_Entries) {
			int numIons = entry.numIons();
			boolean allInts = true;
			for (int ionNum = 0; ionNum < numIons; ionNum++) {
				double oxidationState = entry.getIon(ionNum).getOxidationState();
				int intValue = (int) Math.round(oxidationState);

				if (Math.abs(oxidationState - intValue) > tolerance) {
					allInts = false;
					break;
				}
			}
			if (allInts) {
				newEntryList.add(entry);
			}
		}

		m_Entries = newEntryList;
		Status.basic(this.numEntries() + " remaining entries in data set with all-integer oxidation states.");

	}

	/**
	 * Removes all structures with energy above hull greater than the provided
	 * value. If the energy above the hull is not defined, the entry is removed.
	 * 
	 * @param energyAboveHull The minimum allowed energy above the hull, in eV /
	 *                        atom.
	 */
	public void removeStructuresNotNearHull(double energyAboveHull) {

		ArrayList<Entry> newEntryList = new ArrayList<Entry>();
		for (Entry entry : m_Entries) {

			// If entry.getEnergyAboveHull() is Double.NaN, it will be removed.
			if (entry.getEnergyAboveHull() <= energyAboveHull) {
				newEntryList.add(entry);
			}
		}
		m_Entries = newEntryList;
		Status.basic(this.numEntries() + " remaining entries in data set that are within " + energyAboveHull
				+ " eV of the convex hull.");

	}

	/**
	 * Removes all structures with energy above hull greater than the provided
	 * value. If the energy above the hull is not defined, the entry is not removed.
	 * 
	 * @param maxEnergyAboveHull The minimum allowed energy above the hull, in eV /
	 *                           atom.
	 */
	public void removeUnstableEntries(double maxEnergyAboveHull) {
		ArrayList<Entry> newEntryList = new ArrayList<Entry>();
		for (Entry entry : m_Entries) {
			double energyAboveHull = entry.getEnergyAboveHull();
			if (energyAboveHull > maxEnergyAboveHull) {
				continue;
			}
			newEntryList.add(entry);
		}

		m_Entries = newEntryList;
		Status.basic(this.numEntries() + " remaining entries in data set within " + maxEnergyAboveHull
				+ " eV/atom of convex hull.");
	}

	/**
	 * Gets the lowest integer oxidation state in this data set, where any oxidation
	 * state within 0.01 of an integer is rounded to that integer.
	 * 
	 * @return the lowest integer oxidation state in this data set, where any
	 *         oxidation state within 0.01 of an integer is rounded to that integer.
	 */
	public int getMinIntegerOxidationState() {
		int returnValue = Integer.MAX_VALUE;
		for (Entry entry : m_Entries) {
			int numIons = entry.numIons();
			for (int ionNum = 0; ionNum < numIons; ionNum++) {
				double oxidationState = entry.getIon(ionNum).getOxidationState();
				int intValue = (int) Math.round(oxidationState);
				if (Math.abs(oxidationState - intValue) > 1E-2) {
					continue;
				}
				returnValue = Math.min(intValue, returnValue);
			}
		}
		return returnValue;
	}

	/**
	 * Gets the highest integer oxidation state in this data set, where any
	 * oxidation state within 0.01 of an integer is rounded to that integer.
	 * 
	 * @return the highest integer oxidation state in this data set, where any
	 *         oxidation state within 0.01 of an integer is rounded to that integer.
	 */
	public int getMaxIntegerOxidationState() {
		int returnValue = Integer.MIN_VALUE;
		for (Entry entry : m_Entries) {
			int numSpecies = entry.numIons();
			for (int ionNum = 0; ionNum < numSpecies; ionNum++) {
				double oxidationState = entry.getIon(ionNum).getOxidationState();
				int intValue = (int) Math.round(oxidationState);
				if (Math.abs(oxidationState - intValue) > 1E-2) {
					continue;
				}
				returnValue = Math.max(intValue, returnValue);
			}
		}
		return returnValue;
	}

	/**
	 * The total number of entries in this data set.
	 * 
	 * @return the total number of entries in this data set.
	 */
	public int numEntries() {
		return m_Entries.size();
	}

	/**
	 * Returns the "entryNum"'th entry in this data set.
	 * 
	 * @param entryNum The index of the entry to be returned.
	 * @return the "entryNum"'th entry in this data set.
	 */
	public Entry getEntry(int entryNum) {
		return m_Entries.get(entryNum);
	}

	/**
	 * Add an entry to this data set
	 * 
	 * @param structureID     The ID for this entry. Entries do not need to have
	 *                        unique IDs in the case of monatomic / polyatomic
	 *                        compositions for the same structure, but if non-unique
	 *                        IDs are used in other contexts some functionality
	 *                        might not work as expected.
	 * @param composition     The composition for this entry.
	 * @param ions            The ions (including oxidation states) in this entry.
	 * @param sources         Where this entry came from. Multiple sources are
	 *                        allowed.
	 * @param energyAboveHull The energy above the convex hull, in eV / atom.
	 *                        Double.NaN if unknown.
	 * @param gii             The global instability index for this entry.
	 *                        Double.NaN if unknown.
	 * @return the entry that was added.
	 */
	public Entry addEntry(String structureID, String composition, Ion[] ions, String[] sources, double energyAboveHull,
			double gii) {
		Entry newEntry = new Entry(structureID, composition, ions, sources, energyAboveHull, gii);
		m_Entries.add(newEntry);
		return newEntry;
	}

	/**
	 * Sets the energies above the hull for entries in the given map. The energies
	 * are set for all entries with the given ID, even if multiple entries share the
	 * same ID.
	 * 
	 * @param energiesByID A map in which the key is an entry ID, the value is the
	 *                     energy above the hull in eV / atom, and the key is the
	 *                     entry ID. The energies are set for all entries with the
	 *                     given ID, even if multiple entries share the same ID.
	 */
	public void setEnergiesAboveHull(Map<String, Double> energiesByID) {
		for (Entry entry : m_Entries) {
			String id = entry.getID();
			if (!energiesByID.containsKey(id)) {
				continue;
			}
			entry.setEnergyAboveHull(energiesByID.get(id));
		}
	}

	/**
	 * Writes the calculated likelihood scores, along with information about the
	 * composition and ions, for all entries in this data set to the given file.
	 * 
	 * @param calculator The calculator used to calculate the likelihood score.
	 * @param fileName   The name of the file to be written.
	 */
	public void writeLikelihoods(LikelihoodCalculator calculator, String fileName) {

		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
			for (int entryNum = 0; entryNum < this.numEntries(); entryNum++) {
				Entry entry = this.getEntry(entryNum);
				double value = calculator.optimizeLikelihood(entry.getAllIons()).getMaxLikelihood();
				String outputString = entry.getID() + "\t" + entry.getGivenCompositionString() + "\t" + value + "\t";
				for (int specNum = 0; specNum < entry.numIons(); specNum++) {
					outputString = outputString + entry.getIon(specNum) + ", ";
				}
				outputString = outputString + "\n";
				writer.write(outputString);
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}

	/**
	 * Represents a single data point in the data set
	 * 
	 * @author timmueller
	 *
	 */
	public class Entry {

		private String m_Structure_ID;
		private String m_GivenCompositionString;
		private Composition m_Composition;
		private Ion[] m_Ions;
		private String[] m_Sources;
		private double m_EnergyAboveHull = Double.NaN;
		private double m_GII = Double.NaN;
		private Structure m_Structure;

		/**
		 * Create an entry
		 * 
		 * @param structureID       The ID for this entry's structure. These IDs do not
		 *                          need to be unique in the data set; e.g. two entries
		 *                          can have the same ID if they have the same structure
		 *                          but one has a composition written in terms of only
		 *                          atoms and the other has a composition written in
		 *                          terms of polyatomic ions.
		 * @param compositionString The composition of this entry.
		 * @param ions              A set of ions (including oxidation states) found in
		 *                          this structure
		 * @param sources           Where this information came from. More than one
		 *                          source is allowed
		 * @param energyAboveHull   The energy above the thermodynamic convex hull, in
		 *                          eV / atom. Double.NaN if it is not available.
		 * @param gii               The global instability index. Double.NaN if it is
		 *                          not available.
		 */
		private Entry(String structureID, String compositionString, Ion[] ions, String[] sources,
				double energyAboveHull, double gii) {

			m_Structure_ID = structureID;
			m_GivenCompositionString = compositionString;
			m_Ions = ions;
			m_Sources = sources;
			m_EnergyAboveHull = energyAboveHull;
			m_GII = gii;
		}

		// Adapted from
		// https://stackoverflow.com/questions/8058768/superscript-in-java-string
		private static String removeSubrscript(String str) {
			str = str.replaceAll("₀", "0");
			str = str.replaceAll("₁", "1");
			str = str.replaceAll("₂", "2");
			str = str.replaceAll("₃", "3");
			str = str.replaceAll("₄", "4");
			str = str.replaceAll("₅", "5");
			str = str.replaceAll("₆", "6");
			str = str.replaceAll("₇", "7");
			str = str.replaceAll("₈", "8");
			str = str.replaceAll("₉", "9");
			return str;
		}

		/**
		 * Sets the energy above the thermodynamic hull in eV / atom.
		 * 
		 * @param value The energy above the thermodynamic hull in eV / atom
		 */
		private void setEnergyAboveHull(double value) {
			m_EnergyAboveHull = value;
		}

		/**
		 * Returns the energy above the thermodynamic hull, in eV / atom, or Double.NaN
		 * if it is not available
		 * 
		 * @return the energy above the thermodynamic hull, in eV / atom, or Double.NaN
		 *         if it is not available
		 */
		public double getEnergyAboveHull() {
			return m_EnergyAboveHull;
		}

		/**
		 * Returns the global instability index, or Double.NaN if it is not available
		 * 
		 * @return the global instability index, or Double.NaN if it is not available
		 */
		public double getGlobalInstabilityIndex() {
			return m_GII;
		}

		/**
		 * The structure ID for this entry. Each entry does not need to have a unique
		 * structure ID; entries representing the same structure may have the same ID.
		 * 
		 * @return the structure ID for this entry. Each entry does not need to have a
		 *         unique structure ID; entries representing the same structure may have
		 *         the same ID.
		 */
		public String getID() {
			return m_Structure_ID;
		}

		/**
		 * The atomic structure of this entry.
		 * 
		 * @return the atomic structure of this entry.
		 */
		public Structure getStructure() {
			if (m_Structure == null) {
				String fileName = m_StructDir + this.getID() + ".vasp";
				POSCAR infile = new POSCAR(fileName);
				m_Structure = new Structure(infile);
			}
			return m_Structure;
		}

		/**
		 * The composition of this entry
		 * 
		 * @return the composition of this entry
		 */
		public Composition getComposition() {
			if (m_Composition == null) {
				m_Composition = new Composition(removeSubrscript(m_GivenCompositionString));
			}
			return m_Composition;
		}

		/**
		 * The composition of this entry
		 * 
		 * @return the composition of this entry
		 */
		public String getGivenCompositionString() {
			return m_GivenCompositionString;
		}

		/**
		 * THe number of distinct ions in this entry
		 * 
		 * @return the number of distinct ions in this entry
		 */
		public int numIons() {
			return m_Ions.length;
		}

		/**
		 * Return the specNum'th ion in this entry
		 * 
		 * @param specNum The index of the ion to return
		 * @return the specNum'th ion in this entry
		 */
		public Ion getIon(int specNum) {
			return m_Ions[specNum];
		}

		/**
		 * The set of distinct ions in this entry
		 * 
		 * @return the set of distinct ions in this entry
		 */
		public Ion[] getAllIons() {
			return m_Ions.clone();
		}

		/**
		 * The sources of data for this entry
		 * 
		 * @return the sources of data for this entry
		 */
		public String[] getAllSources() {
			return m_Sources.clone();
		}
	}
}
