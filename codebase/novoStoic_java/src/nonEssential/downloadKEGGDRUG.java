package nonEssential;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.HttpURLConnection;
import java.net.ProtocolException;
import java.net.URL;
import java.sql.Connection;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.io.IOUtils;

import chemaxon.descriptors.FPRetriever;

import com.opencsv.CSVReader;

public class downloadKEGGDRUG {
	static Connection conn;
	static Statement st;
	private static Pattern p;
	private static Matcher m;
	static String cpd;
	private static HttpURLConnection httpConn;
	private static PrintWriter outputStream;
	
	
	static String fileExt = ".mol";

	public static void main(String[] args) {
		try {
			
			
				CSVReader reader = new CSVReader(new FileReader("DRUG.txt"), '\t', '\'');
				String[] nextLine;
				String drugId;
				String urlString;
				URL url;
				String smilesFileName = new File(System.getProperty("user.dir")).getName()+".smiles"; 
				File dir = new File(System.getProperty("user.dir"));
				File[] listFiles = dir.listFiles();
				HashSet<String> completedFile = new HashSet<String>();
				for (int i = 0; i < listFiles.length; i++) {
					if (listFiles[i].getName().contains(".mol")) {
						completedFile.add(listFiles[i].getName().replace(".mol", ""));
					}
				}
				while ((nextLine = reader.readNext()) != null) {
					drugId = nextLine[0];
					try {
						drugId = drugId.replace("dr:", "");
						urlString = "http://rest.kegg.jp/get/"+drugId+"/mol";
						if (completedFile.contains(drugId))
						{
							System.err.println(drugId + " Already downloaded");
							continue;
						}						
						
						url = new URL(urlString);					
						callCreateFile(url,drugId);
					} catch (Exception e) {
						System.err.println(drugId);						
					}
				}
				reader.close();		

		} 
		catch (FileNotFoundException e1) {			
			e1.printStackTrace();
		}
		catch (IOException e) {
			
			e.printStackTrace();
		} 
	}

	private static void callCreateFile(URL url, String abbreviation) throws SQLException,
			IOException, ProtocolException {
		/*conn = Connect.getConnection();
		st = conn.createStatement();*/
		InputStream resolvedIpStream = resolve(url);
		String output;
		
		StringWriter writer = new StringWriter();
		IOUtils.copy(resolvedIpStream, writer);
		String molFileContent = writer.toString();
		
		
		createFile(molFileContent, abbreviation, fileExt);			
		
		httpConn.disconnect();
	}

	private static InputStream resolve(URL url) throws IOException, ProtocolException {
		BufferedReader br;
		httpConn = (HttpURLConnection) url.openConnection();
		httpConn.setRequestMethod("GET");
		httpConn.setRequestProperty("Accept", "application/json");

		if (httpConn.getResponseCode() != 200) {
			throw new RuntimeException("Failed : HTTP error code : "
					+ httpConn.getResponseCode());
		}
			
		return httpConn.getInputStream();
	}

	private static void createFile(String FileContent,String fileName, String fileExt) {
		
		try {			
			outputStream = new PrintWriter(new FileWriter(
					fileName.concat(fileExt)));
			outputStream
					.print(FileContent);
			outputStream.close();

		} catch (IOException e) {

			System.out.println("IOException:");
			e.printStackTrace();

		}
	}

	private static String extractRegex(String stringToFindPatternIn,
			String regex2) {
		p = Pattern.compile(regex2);
		m = p.matcher(stringToFindPatternIn);
		while (m.find()) {
			return m.group();
		}
		return "";

	}

}
