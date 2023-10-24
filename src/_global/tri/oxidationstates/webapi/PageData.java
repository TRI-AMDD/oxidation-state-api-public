package _global.tri.oxidationstates.webapi;

import java.util.Set;
import java.util.TreeSet;

import com.fasterxml.jackson.databind.ObjectMapper;

import _global.tri.oxidationstates.calculator.likelihood.LikelihoodCalculator;
import _global.tri.oxidationstates.util.Composition;
import matsci.util.arrays.ArrayUtils;

/**
 * This class contains all of the data that is needed to render a response to a request
 * on the oxidation state app.
 * 
 * @author timmueller
 *
 */
public class PageData {
	
	private TableData m_TableData;
	private String m_Composition;

	private RangeData[] m_OxidationStateRangeData = new RangeData[0];
	private double m_MinBoundaryValue;
	private double m_MaxBoundaryValue;
	
	private Message[] m_Messages = new Message[0];
	private boolean[] m_IsErrorMessage = new boolean[0];
	
	private PotentialMapper m_PotentialMapper = new PotentialMapper();
	
	/**
	 * Construct object containing data for this page
	 * 
	 * @param composition The composition for which this page was generated.
	 * @param rows The rows to be included in the table on this page.
	 * @param calculator The calculator used to generate this page.  Oxidation state ranges will be
	 * read from this calculator.
	 */
	public PageData(String composition, TableRow[] rows, LikelihoodCalculator calculator) {
		this(composition, calculator);
		this.setTableData(rows);
	}
	
	/**
	 * Construct object containing data for this page
	 * 
	 * @param composition The composition for which this page was generated.
	 * @param calculator The calculator used to generate this page.  Oxidation state ranges will be
	 * read from this calculator.
	 */
	public PageData(String composition, LikelihoodCalculator calculator) {
		m_Composition = composition;
		
		try {
			Set<String> symbols = new Composition(composition).getKnownIonComposition().getCompositionMap().keySet();
			TreeSet<String> allSymbols = new TreeSet<String>();
			for (String symbol : symbols) {
				symbol = symbol.replace("(", "~"); // Hack to get the polyatomic ions to sort to last place
				allSymbols.add(symbol);
			}
			m_OxidationStateRangeData = new RangeData[allSymbols.size()];
			int symbolNum = 0;
			for (String symbol : allSymbols) {
				symbol = symbol.replace("~", "(");
				m_OxidationStateRangeData[symbolNum++] = new RangeData(symbol, calculator);
			}
		} catch (Exception e) {
			// Do nothing, as the error will be set by the user.  TODO:  Handle these internally?
		}
		
		double[] rangeMinAndMax = calculator.getMinAndMaxBoundary();
		m_MinBoundaryValue = rangeMinAndMax[0]; // For visualization
		m_MaxBoundaryValue = rangeMinAndMax[1]; // For visualization
		
		this.setTableData(new TableRow[0]);
	}
	
	/**
	 * Returns the composition of the material for which the oxidation states were generated.
	 * 
	 * @return the composition of the material for which the oxidation states were generated.
	 */
	public String getComposition() {
		return m_Composition;
	}
	
	/**
	 * Returns the information needed to generate the oxidation state ranges for the ions used in this table.
	 * 
	 * @return the information needed to generate the oxidation state ranges for the ions used in this table.
	 * 
	 */
	public RangeData[] getOxidationStateRangeData() {
		return m_OxidationStateRangeData.clone();
	}
	
	/**
	 * Returns a list of messages (e.g. error) to be passed to the user
	 * 
	 * @return a list of messages to be passed to the user
	 */
	public Message[] getMessages() {
		return m_Messages.clone();
	}
	
	/**
	 * Returns an array of booleans, set to true if the corresponding messages returned by {@link #getMessages()} contains an error that should
	 * be highlighted for the user, false otherwise.
	 * 
	 * @return an array of booleans, set to true if the corresponding messages returned by {@link #getMessages()} contains an error that should
	 * be highlighted for the user, false otherwise.
	 */
	public boolean[] isErrorMessage() {
		return m_IsErrorMessage.clone();
	}
	
	/**
	 * Returns the minimum of boundary values for all possible ion types in terms of the mapped potential
	 * 
	 * @return the minimum of boundary values for all possible ion types in terms of the mapped potential
	 */
	public double getMinBoundaryValue() {
		return m_PotentialMapper.toMappedPotential(m_MinBoundaryValue);
	}

	/**
	 * Returns the maximum of boundary values for all possible ion types in terms of the mapped potential
	 * 
	 * @return the maximum of boundary values for all possible ion types in terms of the mapped potential
	 */
	public double getMaxBoundaryValue() {
		return m_PotentialMapper.toMappedPotential(m_MaxBoundaryValue);
	}
	
	/**
	 * This method can be used to set a message to be returned to the user, such as an error or warning encountered when
	 * generating the table.
	 * 
	 * @param messageString A message (e.g. an error message) to be returned to the user.  Returnns null if there is no message.
	 * @param isErrorMessage True if the message contains an error that should be highlighted for the user.
	 */
	public void addMessage(String messageString, boolean isErrorMessage) {
		Message message = new Message(messageString, isErrorMessage);
		m_Messages = (Message[]) ArrayUtils.appendElement(m_Messages, message); 
	}
	
	/**
	 * Sets the table data to a new object containing the given rows
	 * 
	 * @param rows The rows to be included in the table
	 */
	public void setTableData(TableRow[] rows) {
		m_TableData = new TableData(rows);
	}
	
	/**
	 * Returns the table data for this page.
	 * 
	 * @return the table data for this page.
	 */
	public TableData getTableData() {
		return m_TableData;
	}
	
	/**
	 * Returns a class that can transform the potentials used internally by the oxidation analyzer to and from mapped potentials
	 * 
	 * @return A class that can transform the potentials used internally by the oxidation analyzer to and from mapped potentials
	 */
	public PotentialMapper getPotentialMapper() {
		return m_PotentialMapper;
	}
	
	/**
	 * Returns a JSON-formatted string containing all of the information on this page.
	 * 
	 * @return a JSON-formatted string containing all of the information on this page.
	 */
	public String toJSON() {
		ObjectMapper mapper = new ObjectMapper();
		
		try {
			return mapper.writeValueAsString(this);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
	
	@Override
	public String toString() {
		String returnString = "Table data for " + m_Composition + "\n";
		for (int messageNum = 0; messageNum < m_Messages.length; messageNum++) {
			Message message = m_Messages[messageNum];
			if (message.getIsErrorMessage()) {
				returnString += "Error: " + message.getMessageString() + "\n";
			} else {
				returnString += "Message to the user: " + message.getMessageString() + "\n";
			}
		}
		returnString += m_TableData.toString();
		return returnString;
	}
	
}

