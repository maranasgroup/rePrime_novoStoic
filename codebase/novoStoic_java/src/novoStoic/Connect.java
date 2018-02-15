package paper6;

import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Properties;

/*import chemaxon.util.ConnectionHandler;*/

/*import chemaxon.util.ConnectionHandler;*/



public class Connect
{

    /**
     * @param args
     */
    public static void main(String[] args)
    {
	// TODO Auto-generated method stub

    }

    public static Connection getConnection()
    {
	Connection conn;
	try
	{
	    Properties properties = new Properties();
	    try
	    {
		properties.load(new FileInputStream("DB_local.properties"));
	    } catch (IOException e)
	    {
		System.out.println("Unable to read the properties file");
		traceFile(e, "propertyFile");
		System.exit(0);
	    }

	    String getP = properties.getProperty("driver").trim();
		Class.forName(getP)
		    .newInstance();
	    String url = properties.getProperty("url").trim();
	    conn = DriverManager.getConnection(url, properties.getProperty(
		    "user").trim(), properties.getProperty("password").trim());
	    
	    /*if (conn.isValid(0))
	    {*/
		System.out.println("Connection open");
		return conn;
	    /*}*/

	} catch (ClassNotFoundException ex)
	{
	    System.err.println(ex.getMessage());
	} catch (IllegalAccessException ex)
	{
	    System.err.println(ex.getMessage());
	} catch (InstantiationException ex)
	{
	    System.err.println(ex.getMessage());
	} catch (SQLException ex)
	{
	    System.err.println(ex.getMessage());
	}
	System.out.println("Unable to connect");
	return null;

    }

    public static void traceFile(Exception err, String pointOfOrigin)
    {
	String dateTime = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss")
		.format(new Date());
	try
	{

	    PrintWriter outputStream = new PrintWriter(new FileWriter(
		    "TraceFile.txt", true)); // change this in the property
	    // file.
	    outputStream.println(pointOfOrigin
		    + "+---------- Another Error encountered ----------+"
		    + dateTime);
	    outputStream.println(err.toString());
	    outputStream.close();

	} catch (IOException e)
	{

	    System.out.println("IOException:");
	    e.printStackTrace();

	}
    }

    public static int[][] split(int[] val, int splitSize)
    {

	int newBlockSize = (int) Math.ceil(val.length / splitSize + 1);
	System.out.println(newBlockSize);
	int[][] K = new int[splitSize][newBlockSize];
	int blocksCovered = 0;
	int currentval = 0;
	for (int i = 0; i < splitSize; i++)
	{
	    System.out.println(i);
	    for (int l = 0; l < newBlockSize; l++)
	    {
		currentval = l + blocksCovered;
		System.out.println("->" + currentval);
		try
		{
		    K[i][l] = val[l + blocksCovered];
		} catch (ArrayIndexOutOfBoundsException e)
		{
		    // TODO: handle exception
		}
	    }
	    blocksCovered = blocksCovered + newBlockSize;
	    if (blocksCovered > val.length)
	    {
		blocksCovered = val.length;
	    }

	}

	return K;

    }

/*	public static ConnectionHandler connectToJbase() throws SQLException,
			ClassNotFoundException {
		ConnectionHandler conh = new ConnectionHandler();
		conh.setDriver("org.gjt.mm.mysql.Driver");
		conh.setUrl("jdbc:mysql://130.203.234.86/Jchem");
		conh.setLoginName("root");
		conh.setPassword("login@123");
		conh.connectToDatabase();
		return conh;
	}*/

} 