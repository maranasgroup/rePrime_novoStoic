package nonEssential;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import com.google.common.collect.HashMultimap;
import com.opencsv.CSVReader;

public class getMoleculesWithRGroup {

	

	public static void main(String[] args) {
		String nameFile = "names.tab";
		HashMultimap<String, String> names = getNames(nameFile);
		
	
		try {
			BufferedReader inputStream = new BufferedReader(new FileReader("ALL.smiles"));
			PrintWriter outputStream = new PrintWriter(new FileWriter("RMolecules.smilesAndName"));
			String inLine;
			String[] split;
			String id;
			String source;
			String smiles;
			String type;
			while ((inLine = inputStream.readLine()) != null) {
				split = inLine.split("\t");
				id = split[2];
				source = split[4];
				smiles = split[0];
				
				if (smiles.contains("[H-")) {
					System.out.println(smiles);
				}
				
				if (smiles.contains("*")) {					
					outputStream.println(smiles+"\t"+id+"\t"+source+"\t"+names.get(id+""+source).toString());
				/*	if (source.equalsIgnoreCase("kegg")) {
						System.out.println(smiles+"\t"+id+"\t"+source+"\t"+names.get(id+""+source).toString());
					}*/
				}
				
			}
			inputStream.close();
			outputStream.close();
		}
		catch (IOException e) {
			System.err.println("IOException:");
			e.printStackTrace();
		}

	}

	private static HashMultimap<String, String> getNames(String nameFile) {
		HashMultimap<String, String> returnVal = HashMultimap.create();
		try {
			CSVReader reader = new CSVReader(new FileReader(nameFile), '\t', '\'', 3);
			String[] nextLine;
			while ((nextLine = reader.readNext()) != null) {
				returnVal.put(nextLine[0].replace("cpd:", "")+""+nextLine[2], nextLine[1]);			
			}
			reader.close();
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return returnVal;
	}

}
