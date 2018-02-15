package primeR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.Semaphore;

import chemaxon.calculations.clean.Cleaner;
import chemaxon.formats.MolExporter;
import chemaxon.struc.Molecule;

import com.opencsv.CSVReader;
/*alias = "e"
ETYPE = "e-h", "e-lp", "e-c", "2-e","3-e"

alias = "en"
INDEX has the property that needs to be checked, "-1,0,1"

alias = "rad", it is its own 

alias = "chrl" has multiple, need to list out these
INDEX has the property that needs to be checked

alias = "hbz" has multiple, need to list out these
INDEX has the property that needs to be checked

alias = "rng", it is its own
alias = "erng", it is its own

Then we have the atoms and their atomic numbers.

The we have the bonds, which will always get the prime 2*/

public class console {

	
	private static final String BEARCLAW = "|||";
	private static final int THREADCOUNT = 12;
	
	private static BufferedWriter outputStream;
	private static BufferedWriter processedStream;
	private static boolean processOnlyConnected = true;
	private static int level;
	private static Set<Molecule> processedMolecules = Collections.synchronizedSet(new HashSet<Molecule>());;

	public static void main(String[] args) {
		/* workflows will be designed here, most functions will be written in other classes */
		/**/
		
		String metaboliteSmilesFile = "All.smiles"; /*we wont need all since we need to run this only for connected metabolites*/
		String connectedMetabolitesFile = "connectedMetabolites.tab";
		String processedFile = "processedMetabolites.tab";
		String[] nextLine;
		
		HashSet<String> processedMetabolites = getProcessedMetabolitesSet(processedFile);
		HashSet<String> connectedMetabolites = getConnectedMetabolitesSet(connectedMetabolitesFile,processedMetabolites);
		
		HashMap<String,Integer> productPrimeLookup = new HashMap<String, Integer>();
		 
		
		int record = 0;		
		
		try {
			outputStream = new BufferedWriter(new FileWriter("mrv.db",true));
			processedStream = new BufferedWriter(new FileWriter(processedFile,true));
			
			CSVReader readTSV = fileHandler.readTSV(metaboliteSmilesFile, 0, '\0');
			final Semaphore s = new Semaphore(THREADCOUNT);
			
			while ((nextLine = readTSV.readNext()) != null) {
				
								
				if (processOnlyConnected) {
					if (!connectedMetabolites.contains(nextLine[2]+BEARCLAW+nextLine[4])) continue;
				}
				record++;
				
				if (record > 2) {
					return;
				}
				
				final String[] inLine = nextLine;
				outputStream.flush();				
				s.acquireUninterruptibly();
				new Thread(String.valueOf(record)) {
					public void run() {
						try {
							process_zero(inLine);
							/*storage will happen here*/
						} finally {
							s.release();							
						}
					}
				}.start();
				
				if (record > 1000) {					
					/* wait for the threads to complete, so that all the products are known*/
					while (s.availablePermits() < THREADCOUNT) {				
					}
					record = 0;
					level = 0;
					
					
					
					for (int i = 1; i <= molProcess.MAXLEVEL; i++) {
						while (s.availablePermits() < THREADCOUNT) {				
						}
						/*read existing products*/
						productPrimeLookup = molProcess.updatedLookup();						
						final int level = i;
						for (final Molecule p_mol : processedMolecules) {							
							s.acquireUninterruptibly();
							new Thread() {
								public void run() {
									try {
										molProcess.updatePrimes(p_mol,level);
										molProcess.calcProduct(p_mol, level);
									} finally {
										s.release();
									}
								}
							}.start();									
						}
						/*pause for all threads to close*/
						while (s.availablePermits() < THREADCOUNT) {				
						}					
						
						/*update new products with primes*/
						/*assign primes to molecules*/
						/*calculate new products*/
					}
					
					
					/* each product will now get a prime number, assign them into a lookup table, print them as an update into a local file*/
					/**/
				}
				
				
			}
			
			while (s.availablePermits() < THREADCOUNT) {				
			}
			
			
			outputStream.close();
			
			PrintWriter outputStream2 = new PrintWriter(new FileWriter(PRODUCT_LOOKUP));
			for (String prop : molProcess.sortedProducts) {
				outputStream2.println(prop);
			}
			outputStream2.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}

	}



	private static HashSet<String> getConnectedMetabolitesSet(String connectedMetabolitesFile, HashSet<String> processedMetabolites) {
		HashSet<String> returnVal = new HashSet<String>();
		try {
			CSVReader reader = new CSVReader(new FileReader(connectedMetabolitesFile), '\t', '\0', 1);
			String[] nextLine;
			while ((nextLine = reader.readNext()) != null) {
				// nextLine[] is an array of values from the line
				if (!processedMetabolites.contains(nextLine[0]+BEARCLAW+nextLine[1])) {
					returnVal.add(nextLine[0]+BEARCLAW+nextLine[1]);
				}
				
			}
			reader.close();
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		catch (IOException e) {			
			e.printStackTrace();
		}
		return returnVal;
	}	

	private static HashSet<String> getProcessedMetabolitesSet(String processedMetabolitesFile) {
		HashSet<String> returnVal = new HashSet<String>();
		try {
			BufferedReader inputStream = new BufferedReader(new FileReader(processedMetabolitesFile));
			String inLine;
			while ((inLine = inputStream.readLine()) != null) {
				returnVal.add(inLine);
			}
			inputStream.close();
		}
		catch (IOException e) {
			System.err.println("IOException:");
			e.printStackTrace();
		}
		return returnVal;
	}

	private static void process_zero(String[] nextLine) {
		String smiles = nextLine[0].replace(" tab", "").trim();				
		String id = nextLine[2];
		String type = nextLine[3];
		String source = nextLine[4];		
		
		Molecule mol = molProcess.getMolObject(smiles);
		
		/*setting global properties*/				
		mol.setProperty("id", id);
		mol.setProperty("source", source);
		mol.setProperty("structure format", type);
		
		mol = molProcess.standardizeMol(mol);
		mol = molProcess.annotatedMolecule(mol, true);
		
		HashMap<Molecule, Double> mmDistribution = molProcess.getMMDistribution(mol);
		Set<Molecule> molSet = mmDistribution.keySet();		
		
		for (Molecule mol2 : molSet) {
			/*TODO: combine all the information into a single cml file,for e.g, the protons, changes radical etc which are not common for all atoms will be given phSpecific information*/
				mol2 = molProcess.annotatedMolecule(mol2, false);
				molProcess.addPseudo(mol2, mmDistribution.get(mol2));
				molProcess.addPrimeZero(mol2);
				molProcess.calcProduct(mol2,0);
				
				/*store it somewhere*/
				processedMolecules.add(mol2);
				
/*				new Cleaner();
				Cleaner.clean(mol2, 2);
				String out = MolExporter.exportToFormat(mol2, "mrv:S");
				outputStream.flush();
				outputStream.newLine();
				outputStream.write(out);
				outputStream.newLine();
				outputStream.write("****record seperator****");
				outputStream.newLine();
				outputStream.flush();*/
		}
		try {
			processedStream.flush();
			processedStream.write(id+BEARCLAW+source);
			processedStream.newLine();
			processedStream.flush();
		}
		catch (IOException e) {			
			e.printStackTrace();
		}
		
	}

	

}
