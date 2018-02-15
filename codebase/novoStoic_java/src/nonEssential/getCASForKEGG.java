package nonEssential;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.ProtocolException;
import java.net.URL;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.StringUtils;

import com.opencsv.CSVWriter;

public class getCASForKEGG {	
	private static HttpURLConnection httpConn;
	private static PrintWriter outputStream;
	private static Matcher m;

	public static void main(String[] args) {
		processForCAS("compound");
		processForCAS("drug");
	}
	private static void processForCAS(String db) {
		try {
			String urlString = "http://rest.kegg.jp/list/"+db;
			URL url = new URL(urlString);
			InputStream resolvedIpStream = resolve(url);
			HashSet<HashSet<String>> store = new HashSet<HashSet<String>>();
			Pattern entryPattern = Pattern.compile("(ENTRY.+)");
			Pattern idPattern = Pattern.compile("(C[0-9]+|D[0-9]+)");
			Pattern casPattern = Pattern.compile("[C][A][S].+");
			try {
				CSVWriter writer = new CSVWriter(new BufferedWriter(new FileWriter(db+".cas")), '\t', '"');
				String[] toPrint = new String[2];
				toPrint[0] = "id";
				toPrint[1] = "CAS";				
				
				BufferedReader inputStream = new BufferedReader(new InputStreamReader(resolvedIpStream));
				String inLine;
				int j = 0;
				HashSet<String> ToStore = new HashSet<String>();
				while ((inLine = inputStream.readLine()) != null) {
					if (j == 0) {
						ToStore = new HashSet<String>();
					}
					j++;
					if (j > 9) {
						j = 0;
						store.add(ToStore);
					}
					ToStore.add(inLine.split("\t")[0]);					
				}
				inputStream.close();
				
				for (HashSet<String> hashSet : store) {
					urlString = "http://rest.kegg.jp/get/"+StringUtils.join(hashSet.toArray(),"+");
					url = new URL(urlString);
					resolvedIpStream = resolve(url);
					StringWriter resolvedWriter = new StringWriter();
					IOUtils.copy(resolvedIpStream, resolvedWriter);
					String listOf = resolvedWriter.toString();
					String[] split = listOf.split("[/][/][/]");
					for (int i = 0; i < split.length; i++) {
						
						if (null != regexSearch(split[i], casPattern)) {							
							toPrint[0] = regexSearch(split[i], casPattern).replace("CAS:", "").trim();
							toPrint[1] = regexSearch(split[i], idPattern);
							System.out.println(toPrint[1]);
							writer.writeNext(toPrint);
						}
						
					}
				}
				writer.close();
				
			}
			catch (IOException e) {
				System.err.println("IOException:");
				e.printStackTrace();
			}
			
			
			
		}
		catch (MalformedURLException e) {			
			e.printStackTrace();
		}
		catch (ProtocolException e) {			
			e.printStackTrace();
		}
		catch (IOException e) {		
			e.printStackTrace();
		}
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
	private static String regexSearch(String input,Pattern r) {		
		m = r.matcher(input);
		while (m.find()) {
			return m.group();
		}
		return null;
	}
}
