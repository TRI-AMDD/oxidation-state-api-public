package _global.tri.oxidationstates.util;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import _global.tri.oxidationstates.structure.LatticeOptimizer;
import _global.tri.oxidationstates.structure.StructureOptimizer;
import _global.tri.structure.mapper.GeneralStructureMapper;
import matsci.Species;
import matsci.engine.minimizer.CGMinimizer;
import matsci.io.vasp.POSCAR;
import matsci.location.Coordinates;
import matsci.structure.BravaisLattice;
import matsci.structure.Structure;
import matsci.structure.StructureBuilder;
import matsci.util.MSMath;

/**
 * This class contains some utility methods for identifying and working with
 * polyatomic ions.
 * 
 * @author timmueller
 *
 */
public class IonTools {

	/**
	 * Randomly samples POSCAR-formatted files from the given directory and creates
	 * an "average" structure from those files.
	 * 
	 * @param dirName    The directory containing the files to be sampled.
	 * @param numSamples The number of files to use to create the average structure.
	 *                   Files are sampled without replacement.
	 * @return the average structure
	 */
	public static Structure makeRepresentativeStructure(String dirName, int numSamples) {

		File dir = new File(dirName);
		File[] structureFiles = dir.listFiles();
		if (structureFiles.length == 0) {
			return null;
		}
		ArrayList<Structure> structures = new ArrayList<Structure>();
		for (File structureFile : structureFiles) {
			if (!structureFile.getName().endsWith(".vasp")) {
				continue;
			}
			if (structureFile.getName().endsWith("mean.vasp")) {
				continue;
			}
			Structure structure = loadIonStructureFromFile(structureFile.getAbsolutePath());
			structure.setDescription(structure.getDescription() + " " + structureFile.getName());
			structures.add(structure);
		}
		if (structures.size() == 1) {
			return structures.get(0);
		}

		Structure template = findRepresentativeStructure(structures, numSamples, false);

		/**
		 * We do this because the mapper will do it anyway. It ensures that the mapped
		 * lattices correspond to the template lattice.
		 */
		template = template.getCompactStructure();
		BravaisLattice templateLattice = template.normalizeNonPeriodicVectors().getDefiningLattice();
		BravaisLattice bestLattice = templateLattice;
		CGMinimizer cgEngine = new CGMinimizer();

		// Find the representative lattice
		if (template.numPeriodicDimensions() > 0) {
			GeneralStructureMapper.Map[] bestMaps = new GeneralStructureMapper.Map[structures.size()];
			for (int structNum = 0; structNum < structures.size(); structNum++) {
				Structure structure = structures.get(structNum);
				GeneralStructureMapper mapper = new GeneralStructureMapper(template, structure);
				mapper.setDistanceTolerance(0.07 * 5); // TODO come up with a better way to determine these
				mapper.setSkewTolerance(0.08 * 4);
				mapper.setVolumeTolerance(0.15 * 5);
				GeneralStructureMapper.Map[] maps = mapper.getMaps();

				double minSkewFactor = Double.POSITIVE_INFINITY;
				for (GeneralStructureMapper.Map map : maps) {
					double skewFactor = map.getSkewFactor();
					if (skewFactor < minSkewFactor) {
						minSkewFactor = skewFactor;
						bestMaps[structNum] = map;
					}
				}
				if (bestMaps[structNum] == null) {
					mapper = new GeneralStructureMapper(template, structure);
					mapper.setDistanceTolerance(0.07 * 5);
					mapper.setSkewTolerance(0.08 * 4);
					mapper.setVolumeTolerance(0.15 * 5);
					maps = mapper.getMaps();
				}
			}

			ArrayList<BravaisLattice> lattices = new ArrayList<BravaisLattice>();
			for (int mapNum = 0; mapNum < bestMaps.length; mapNum++) {
				lattices.add(bestMaps[mapNum].getSubLattice());
			}

			LatticeOptimizer optimizer = new LatticeOptimizer(template.getDefiningLattice(),
					lattices.toArray(new BravaisLattice[0]));
			cgEngine.minimize(optimizer);

			LatticeOptimizer bestState = (LatticeOptimizer) cgEngine.getMinimumState();
			bestLattice = bestState.getLattice();
		}

		// TODO assumes that the template was also used to find best lattice -- enforce
		// this.
		StructureBuilder builder = new StructureBuilder(template);
		builder.setCellVectors(bestLattice.getCellVectors());
		builder.setVectorPeriodicity(bestLattice.getDimensionPeriodicity());
		for (int siteNum = 0; siteNum < builder.numDefiningSites(); siteNum++) {
			double[] coordArray = builder.getSiteCoords(siteNum).getCoordArray(templateLattice.getLatticeBasis());
			Coordinates coords = new Coordinates(coordArray, bestLattice.getLatticeBasis());
			builder.setSiteCoordinates(siteNum, coords);
		}
		template = new Structure(builder);

		GeneralStructureMapper.Map[] bestMaps = new GeneralStructureMapper.Map[structures.size()];
		for (int structNum = 0; structNum < structures.size(); structNum++) {
			Structure structure = structures.get(structNum);
			GeneralStructureMapper mapper = new GeneralStructureMapper(template, structure);
			mapper.setDistanceTolerance(0.07 * 5);
			mapper.setSkewTolerance(0.08 * 5);
			mapper.setVolumeTolerancePerDimension(0.07 * 5);
			GeneralStructureMapper.Map[] maps = mapper.getMaps();
			double minDistanceError = Double.POSITIVE_INFINITY;
			for (GeneralStructureMapper.Map map : maps) {
				double distanceError = map.getDistanceFactor();
				if (distanceError < minDistanceError) {
					minDistanceError = distanceError;
					bestMaps[structNum] = map;
				}
			}
		}

		ArrayList<Structure> bestStructures = new ArrayList<Structure>();
		for (int mapNum = 0; mapNum < bestMaps.length; mapNum++) {
			bestStructures.add(bestMaps[mapNum].getTransformedSubStructure(true));
		}

		StructureOptimizer structOptimizer = new StructureOptimizer(bestStructures.toArray(new Structure[0]));
		cgEngine.minimize(structOptimizer);

		StructureOptimizer bestState2 = (StructureOptimizer) cgEngine.getMinimumState();
		Structure bestStructure = bestState2.getStructure();
		return bestStructure;
	}

	/**
	 * Returns true if two structures are equivalent based on the standard mapping
	 * parameters set in {@link setMappingParams}, and false otherwise.
	 * 
	 * @param structure1            One of the structures to be compared.
	 * @param structure2            One of the structures to be compared.
	 * @param ignoreOxidationStates True if oxidation states should be ignored when
	 *                              making the comparison (Fe2+ and Fe3+ are the
	 *                              same), and false otherwise (Fe2+ and Fe3+ are
	 *                              different).
	 * @return true if two structures are equivalent based on the standard mapping
	 *         parameters set in {@link setMappingParams}, and false otherwise.
	 */
	public static boolean compareStructures(Structure structure1, Structure structure2, boolean ignoreOxidationStates) {

		if (ignoreOxidationStates) {
			structure1 = structure1.removeOxidationStates();
			structure2 = structure2.removeOxidationStates();
		}

		double scaleFactor = Math.pow(structure1.getDefiningVolume() / structure2.getDefiningVolume(),
				1.0 / structure1.numPeriodicDimensions());

		if (structure1.numPeriodicDimensions() != structure2.numPeriodicDimensions()) {
			return false;
		}

		if (structure1.numDefiningSites() != structure2.numDefiningSites()) {
			return false;
		}

		Species[] allSpecies = structure1.getDistinctSpecies();
		for (Species species : allSpecies) {
			int numAtoms1 = structure1.numDefiningSitesWithSpecies(species);
			int numAtoms2 = structure2.numDefiningSitesWithSpecies(species);
			if (numAtoms1 != numAtoms2) {
				return false;
			}
		}

		GeneralStructureMapper mapper = new GeneralStructureMapper(structure1, structure2);
		setMappingParams(mapper);
		return (mapper.numMaps() > 0);
	}

	/**
	 * From a list of given structures, finds the one that is closest to the mean of
	 * random subset of the others based on the sum of the scores of the maps
	 * between the structures. If not such structure can be found (usually because
	 * there are no reasonable maps between any one structure and all other
	 * structures), return null.
	 * 
	 * @param structures            The set of structures for which we will
	 *                              calculate the average.
	 * @param numSamples            The number of files to use to create the average
	 *                              structure. Structures are sampled without
	 *                              replacement.
	 * @param ignoreOxidationStates True if oxidation states should be ignored when
	 *                              making the comparison (Fe2+ and Fe3+ are the
	 *                              same), and false otherwise (Fe2+ and Fe3+ are
	 *                              different).
	 * @return The average structure
	 */
	public static Structure findRepresentativeStructure(List<Structure> structures, int numSamples,
			boolean ignoreOxidationStates) {

		structures = new ArrayList<Structure>(structures);
		if (ignoreOxidationStates) {
			for (int structNum = 0; structNum < structures.size(); structNum++) {
				Structure structure = structures.get(structNum);
				structure = structure.removeOxidationStates();
				structures.set(structNum, structure);
			}
		}

		Structure bestStructure = null;
		double bestScore = Double.POSITIVE_INFINITY;
		int[] structOrder = MSMath.getRandomShuffle(structures.size());
		int numStructs = Math.min(structOrder.length, numSamples); // We sample up to this many structures
		for (int structNum = 0; structNum < numStructs; structNum++) {
			Structure structure1 = structures.get(structOrder[structNum]);
			double scoreSum = 0;
			for (Structure structure2 : structures) {
				if (structure1 == structure2) {
					continue;
				}
				GeneralStructureMapper mapper = new GeneralStructureMapper(structure1, structure2);
				mapper.setDistanceTolerance(0.07 * 5); // TODO come up with a better way to determine these
				mapper.setSkewTolerance(0.08 * 4);
				mapper.setVolumeTolerance(0.15 * 5);
				GeneralStructureMapper.Map bestMap = mapper.getBestMap();
				if (bestMap == null) {
					scoreSum = Double.POSITIVE_INFINITY;
					break;
				}
				scoreSum += bestMap.scoreMap();
				if (scoreSum > bestScore) {
					break;
				}
			}
			if (scoreSum == Double.POSITIVE_INFINITY) { // This isn't as accurate but it's a lot faster.
				continue;
			}
			if (scoreSum < bestScore) {
				bestScore = scoreSum;
				bestStructure = structure1;
			}
		}

		return bestStructure;

	}

	/**
	 * Sets standard structure comparison parameters for the structure mapper
	 * 
	 * @param mapper The mapper for which we are setting the standard parameters.
	 */
	public static void setMappingParams(GeneralStructureMapper mapper) {

		// Used for original ion identification
		mapper.setDistanceTolerance(0.07);
		mapper.setSkewTolerance(0.08);
		mapper.setVolumeTolerance(0.15);
		mapper.setOxidationStateTolerance(0.01);

		// This gets some sulfate ions the other one misses, but might also confuse
		// tetrahedral with trigonal planar
		/*
		 * mapper.setDistanceTolerance(0.12); mapper.setSkewTolerance(0.08);
		 * mapper.setVolumeTolerance(0.15); mapper.setOxidationStateTolerance(0.01);
		 */
	}

	/**
	 * Writes an ion structure to a special type of POSCAR-formatted file that
	 * contains information about the lattice vector periodicity in the first line
	 * of the file.
	 * 
	 * @param ion         The structure file will be written for this ion.
	 * @param ionFileName The name of the file to write.
	 * @param description A description string to be added to the file, along with
	 *                    the information about lattice vector periodicity.
	 */
	public static void writeIonStructure(Structure ion, String ionFileName, String description) {
		boolean[] periodicity = ion.getVectorPeriodicity();
		String periodicityString = "";
		periodicityString += periodicity[0] ? "T" : "F";
		periodicityString += periodicity[1] ? "_T" : "_F";
		periodicityString += periodicity[2] ? "_T" : "_F";
		ion.setDescription(periodicityString + " " + description);
		ion = ion.padNonPeriodicDimensions(10);
		ion = ion.centerNonPeriodicDimensions();
		POSCAR outfile = new POSCAR(ion, true);
		outfile.writeSpeciesWithDescription(false);
		outfile.writeElementSymbolsOnly(false);
		outfile.writeFile(ionFileName);
	}

	/**
	 * Reads an ion structure from a file written in the format defined by {2link
	 * writeIonStructure()}. Sets the lattice vector periodicity appropriately, as
	 * defined in the first line of the file.
	 * 
	 * @param path The path to the file to be read.
	 * @return The structure read from the file, with lattice vector periodicity set
	 *         as specified in the file.
	 */
	public static Structure loadIonStructureFromFile(String path) {
		POSCAR infile = new POSCAR(path);
		String description = infile.getDescription();
		String periodicityString = description.split(" ")[0];
		String[] periodicity = periodicityString.split("_");
		for (int dimNum = 0; dimNum < periodicity.length; dimNum++) {
			infile.setVectorPeriodicity(dimNum, periodicity[dimNum].toLowerCase().equals("t"));
		}
		return new Structure(infile);
	}

}
