package nonEssential;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.concurrent.Semaphore;

import org.apache.commons.io.FilenameUtils;

import chemaxon.calculations.clean.Clean2D;
import chemaxon.calculations.clean.Cleaner;
import chemaxon.calculations.hydrogenize.Hydrogenize;
import chemaxon.fixers.NeutralizeChargeFixer;
import chemaxon.formats.MolConverter;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.calculations.MajorMicrospeciesPlugin;
import chemaxon.marvin.plugin.PluginException;
import chemaxon.standardizer.Standardizer;
import chemaxon.standardizer.actions.NeutralizeAction;
import chemaxon.struc.Molecule;
import chemaxon.util.standardizer.StandardizerUtil;

import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;

public class combineMolsToSmilesFile {

	private static final int THREADCOUNT = 12;
	private static CSVWriter writer;
	

	public static void main(String[] args) {
		
		try {
			String smilesFileName = new File(System.getProperty("user.dir")).getName()+".smiles"; 
			File dir = new File(System.getProperty("user.dir"));
			File[] listFiles = dir.listFiles();
			writer = new CSVWriter(new BufferedWriter(new FileWriter(smilesFileName)), '\t','\0');
			
			final Semaphore s = new Semaphore(THREADCOUNT);
			
			
			
			for (int i = 0; i < listFiles.length; i++) {
				s.acquireUninterruptibly();
				final File strucFile = listFiles[i];
				new Thread() {
					public void run() {
						try {
							try {
								final String[] toPrint = new String[2];
								if (strucFile.getName().contains(".mol"))
								{toPrint[0] = strucFile.getName().replace(".mol", "");
								toPrint[1] = new Scanner(strucFile).useDelimiter("\\A").next();				
								process(toPrint);
								}
							}
							catch (Exception e) {					
								e.printStackTrace();
							}
						} finally {
							s.release();
						}
					}
				}.start();
						
			}
			
			while (s.availablePermits() < THREADCOUNT) {				
			}
			
			writer.close();
		}
		catch (Exception e) {			
			e.printStackTrace();
		}

	}

	private static void process(String[] nextLine) {
		try {
			if (nextLine.length < 2) return;
			String id = nextLine[0].trim();
			Molecule mol = new Molecule();
			String structure = nextLine[1];
			String[] split;			
			int lineCount = 0;
			try {
				mol = MolImporter.importMol(structure);
			}
			catch (Exception e) {
				structure = "\n" + nextLine[1];
				
				try {
					mol = MolImporter.importMol(structure);
				}
				catch (Exception e1) {
					split = structure.split("\n");

					for (String str : split) {
						lineCount++;
						if (str.toLowerCase().contains("v2000")) {
							break;
						}
					}
					for (int k = lineCount; k < 4; k++) {
						structure = "\n" + structure;
					}
					try {
						mol = MolImporter.importMol(structure, "mol");
					}
					catch (Exception e2) {
						System.err.println(id);
					}
				}
			}
			
			mol = getPH_0(mol);
			Hydrogenize.convertExplicitHToImplicit(mol);		
			Cleaner.clean(mol, 2);  
			
			String smiles = MolExporter.exportToFormat(mol, "smiles:u0");
			/*make changes here*/
			int smilesLength = 0;
			if (!smiles.trim().equalsIgnoreCase("")) {				
				if (!smiles.contains("*")) {
					split = new String[2];
					split[1] = id;
					String[] split2 = smiles.split("[.]");
					for (int i = 0; i < split2.length; i++) {
						if (split2[i].length() > smilesLength) {
							smilesLength = split2[i].length();
							split[0] = split2[i];
						}
					}
					writer.writeNext(split);
				}
			}			
		}
		catch (IOException e) {		
			e.printStackTrace();
		}
	}

	private static Molecule getPH_0(Molecule mol) {
		try {
			MajorMicrospeciesPlugin majorMicrospeciesPlugin = new MajorMicrospeciesPlugin();
			majorMicrospeciesPlugin.setMolecule(mol);
			majorMicrospeciesPlugin.setpH(0.00);
			majorMicrospeciesPlugin.run();
			return majorMicrospeciesPlugin.getMajorMicrospecies();
		}
		catch (Exception e) {			
			e.printStackTrace();
		}
		return mol;
	}
}
