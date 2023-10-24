package global_.tri.oxidationstates.util;

import java.io.File;
import java.io.FilenameFilter;

/**
 * Returns only files whose name ends with ".vasp". Used for quickly filtering
 * non-structure files from a directory.
 * 
 * @author timmueller
 *
 */
public class VaspFileFilter implements FilenameFilter {

	@Override
	public boolean accept(File dir, String name) {
		return (name.endsWith(".vasp"));
	}

}
