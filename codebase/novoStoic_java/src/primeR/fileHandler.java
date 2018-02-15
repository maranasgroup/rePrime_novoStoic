package primeR;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.apache.tika.utils.CharsetUtils;

import com.opencsv.CSVReader;

public class fileHandler {	
	public static CSVReader readTSV(String fileName,int skipLines, char quoteChar) throws IOException {
		CSVReader reader = null;
		try {			
			reader = new CSVReader(new FileReader(fileName), '\t',quoteChar,skipLines);			
			return reader;			
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		} 
		return null;
	}
}
