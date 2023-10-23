package global.tri.oxidationstates.webapi;

/**
 * This class represents an message to be returned to the user.
 * 
 * @author timmueller
 *
 */
public class Message {
	
	private String m_Message;
	private boolean m_IsErrorMessage;
	
	/**
	 * 
	 * @param message The message content
	 * @param isErrorMessage True if this is an error, false otherwise
	 */
	public Message(String message, boolean isErrorMessage) {
		m_Message = message;
		m_IsErrorMessage = isErrorMessage;
	}
	
	/**
	 * Returns the message content.
	 * 
	 * @return the message content.
	 */
	public String getMessageString() {
		return m_Message;
	}
	
	/**
	 * Returns true if this is an error message, false otherwise
	 * 
	 * @return true if this is an error message, false otherwise
	 */
	public boolean getIsErrorMessage() {
		return m_IsErrorMessage;
	}

}
