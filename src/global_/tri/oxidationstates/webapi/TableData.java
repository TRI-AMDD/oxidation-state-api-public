package global_.tri.oxidationstates.webapi;

/**
 * This class contains all of the data required for the table in the web app.
 * 
 * @author timmueller
 *
 */
public class TableData {

	private TableRow[] m_TableRows = new TableRow[0];
	
	/**
	 * Construct the table data from an array of table rows.
	 * 
	 * @param rows The rows to be included in this table.
	 */
	public TableData(TableRow[] rows) {
		m_TableRows = rows.clone();
	}
	
	/**
	 * Returns the array of rows in this table.
	 * 
	 * @return the array of rows in this table.
	 */
	public TableRow[] getTableRows() {
		return m_TableRows.clone();
	}
	
	@Override
	public String toString() {
		String returnString = "Oxidation States\tOptimal Likelihood\tOptimal electronic chemical potential\tGII\tMixed Valence\n";
		for (TableRow row : m_TableRows) {
			returnString += row + "\n";
		}
		return returnString;
	}
	
}
