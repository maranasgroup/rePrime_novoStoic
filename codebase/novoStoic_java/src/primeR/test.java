package primeR;

import java.awt.Component;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.calculations.pKaPlugin;
import chemaxon.marvin.calculations.pKaPluginDisplay;
import chemaxon.marvin.plugin.CalculatorPluginDisplay;
import chemaxon.marvin.plugin.PluginException;
import chemaxon.struc.Molecule;

public class test {

	public static void main(String[] args) {
	
		BufferedReader br = null;
		String filePath = "";
		try {

			String ln;

			br = new BufferedReader(new FileReader("primes/primes1.txt"));
			while ((ln = br.readLine()) != null) {
				System.out.println(ln);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)
					br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
	
	}
	
}
