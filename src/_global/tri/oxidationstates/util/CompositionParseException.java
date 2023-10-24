package _global.tri.oxidationstates.util;

/**
 * This exception is thrown if there is a problem parsing a composition
 * 
 * @author timmueller
 *
 */
public class CompositionParseException extends RuntimeException {

	/**
	 * Creates a composition parse exception with the given message.
	 * 
	 * @param message The given message.
	 */
	public CompositionParseException(String message) {
		super(message);
	}
}