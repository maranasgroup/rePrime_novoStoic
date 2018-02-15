package nonEssential;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Set;
import java.util.Vector;

import com.google.common.collect.HashMultimap;
import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolImporter;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;
import chemaxon.struc.RxnMolecule;

public class SCJToSmiles {

	private static HashMap<String, String> naomiSmiles;
	private static HashMap<String, Molecule> cpdSmiles;
	private static HashMap<String, RxnMolecule> rxnMols;

	public static void main(String[] args) {

		/* read Naomi */

		String folderPath = "/home/azk172/metrxn0200/input/metabolites/structure/";
		String cpdfile = "/home/azk172/fromWindows/Akhil/D/paper 6/alchemist/sets/cpd.index";
		String scjFile = "/home/azk172/fromWindows/Akhil/D/paper 6/alchemist/parameters/rule.scj";
		String rxnSmiFile = "/home/azk172/fromWindows/Akhil/D/paper 6/alchemist/mapping/rxn.smi";

		naomiSmiles = getNaomiSmiles(folderPath);

		cpdSmiles = getCpdSmiles(cpdfile);

		rxnMols = readSCJ(scjFile);
		
		print(rxnSmiFile);

		/* if a smiles is not read, return fail and store as error */
		/*
		 * number the atoms, also store it in the smiles and put the number
		 * where the product starts.
		 */

	}



	private static void print(String rxnSmiFile) {
		try {
			CSVWriter writer = new CSVWriter(new FileWriter(rxnSmiFile), '\t',CSVWriter.NO_QUOTE_CHARACTER);
			Set<String> keySet = rxnMols.keySet();
			String[] entries;
			for (String key : keySet) {
				entries = new String[3];
				entries[0] = MolExporter.exportToFormat(rxnMols.get(key), "smiles");
				entries[1] = key.split("\\t")[0];
				entries[2] = key.split("\\t")[1];
				writer.writeNext(entries);
			}			
			writer.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}



	private static HashMap<String, RxnMolecule> readSCJ(String scjFile) {
		BufferedReader br = null;
		String filePath = scjFile;
		HashMultimap<String, Molecule> reactants = HashMultimap.create();
		HashMultimap<String, Molecule> products = HashMultimap.create();
		
		HashMap<String, RxnMolecule> returnVal = new HashMap<String, RxnMolecule>();
		
		try {

			String ln;
			String rid;
			String smiles;
			String mid;
			String direction;
			br = new BufferedReader(new FileReader(filePath));
			String[] split;
			while ((ln = br.readLine()) != null) {
				split = ln.split("[.]|[ ]|[']");
				mid = split[1];
				rid = split[4];
				direction = split[6];
				/*TODO: when the stoichiometry is not 1*/
				if (direction.contains("-")) {
					reactants.put(rid, cpdSmiles.get(mid).clone());
				} else {
					products.put(rid, cpdSmiles.get(mid).clone());
				}				
			}
			Set<String> keySet = reactants.keySet();
			for (String key : keySet) {
				RxnMolecule rxnMol = new RxnMolecule();
				int atomNum = 0;
				Set<Molecule> set = reactants.get(key);
				for (Molecule molecule : set) {
					MolAtom[] atomArray = molecule.getAtomArray();
					for (MolAtom molAtom : atomArray) {
						atomNum++;
						molAtom.setAtomMap(atomNum);
					}
					rxnMol.addComponent(molecule, RxnMolecule.REACTANTS);
				}
				
				
				int productAtomIndex = atomNum + 1;
				set = products.get(key);
				for (Molecule molecule : set) {
					MolAtom[] atomArray = molecule.getAtomArray();
					for (MolAtom molAtom : atomArray) {
						atomNum++;
						molAtom.setAtomMap(atomNum);
					}
					rxnMol.addComponent(molecule, RxnMolecule.PRODUCTS);
				}				
				
				System.out.println(key);
				returnVal.put(key+'\t'+productAtomIndex, rxnMol);								
			}
			return returnVal;
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
		return null;
		
		
		
	}

	private static HashMap<String, Molecule> getCpdSmiles(String cpdfile) {
		HashMap<String, Molecule> returnVal = new HashMap<String,Molecule>();
		BufferedReader br = null;
		String filePath = cpdfile;
		try {

			String ln;
			String key = "";
			br = new BufferedReader(new FileReader(filePath));
			while ((ln = br.readLine()) != null) {
				
				try {
					key = ln.replace("'", "");
					returnVal.put(key, MolImporter.importMol(naomiSmiles.get(key)));
					//System.out.println(key);
				} catch (Exception e) {
					System.out.println(key);
					//e.printStackTrace();
				}
			}
			return returnVal;
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
		return null;
	}

	private static HashMap<String, String> getNaomiSmiles(String folderPath) {
		HashMap<String, String> returnVal = new HashMap<String,String>();
		File dir = new File(folderPath);
		File[] listFiles = dir.listFiles();
		String ideezStr;
		for (File f : listFiles) {
			if (!f.getName().endsWith(".smiles"))
				continue;
			if (!f.getName().contains("naomi"))
				continue;
			try {
				CSVReader reader = new CSVReader(new FileReader(f), '\t');
				String[] nextLine;
				while ((nextLine = reader.readNext()) != null) {
					if (nextLine[3].equalsIgnoreCase("chembl")) continue;
					ideezStr = (nextLine[3]+":").replace("brenda:", "br:").replace("chebi:", "ci:").replace("chembl:", "cl:").replace("ecmdb:", "eb:").replace("ecocyc:", "ey:").replace("hmdb:", "hm:").replace("kegg:", "kg:")
							.replace("metacyc:", "my:")+nextLine[2]+":"+nextLine[1];
					returnVal.put(ideezStr, nextLine[0].replace(" tab", ""));
				}
				

			} catch (FileNotFoundException e) {

				e.printStackTrace();
			} catch (IOException e) {

				e.printStackTrace();
			}
			
		}
		return returnVal;
	}

}
