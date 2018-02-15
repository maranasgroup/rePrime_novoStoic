package primeR;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.derby.tools.sysinfo;
import org.hsqldb.lib.HashSet;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.struc.RxnMolecule;

public class overLayTransitionsIntoCML {

	private static PrintWriter pw;

	public static void main(String[] args) {
		/* read CML, read file, change the INDEXES write it back into a CML */
		/* remove the properties and print to folder. */
		/* generate a stat */
		String rxnId = "BRENDA_29544#2";
		File f = new File("reactions/" + rxnId + ".mrv");
		RxnMolecule rxnMol = gamsFromCML.getRxnMol(f);

		TreeMap<Integer, HashMap<Integer, Integer>> indexMappingAll = getIndexMapping(
				"GAMS/" + rxnId + "/CLCA4.transition");
		
		
		try {
			pw = new PrintWriter("GAMS/" + rxnId + "/CLCA2_mapped.smiles");
			String smiles = rxnMol.getProperty(molProcess.SMILES);
			System.out.println(smiles);
			
			Set<Integer> keySet = indexMappingAll.keySet();
			for (Integer key : keySet) {
				try {
					Molecule mappedMol = getIndicesToAtomMapping(indexMappingAll.get(key), smiles);
					String mappedSmiles = MolExporter.exportToFormat(mappedMol, "smiles");
					int numberOfBondchanges = getBondMismatch(mappedMol);
					int fragmentCount = getFragmentCount(mappedMol);
					
					pw.println(key+"\t"+mappedSmiles+"\t"+numberOfBondchanges+"\t"+fragmentCount);
				} catch (IOException e) {

					e.printStackTrace();
				}
			}
		} catch (FileNotFoundException e) {
			
			e.printStackTrace();
		}
		
		pw.close();

	}

	private static Molecule getIndicesToAtomMapping(HashMap<Integer, Integer> AtomMap, String smiles) throws IOException {

		Molecule rxnMol = MolImporter.importMol(smiles);

		

		MolAtom[] atomArray2 = rxnMol.getAtomArray();
		for (MolAtom molAtom : atomArray2) {
			if (AtomMap.containsKey(molAtom.getAtomMap())) {
				molAtom.setAtomMap(AtomMap.get(molAtom.getAtomMap()));
			}
		}
		
		
		
		return rxnMol;
		
		

	}

	private static int getFragmentCount(Molecule rxnMol) {
		return rxnMol.getFragCount();
		
	}

	private static String getBondIndex(MolBond molBond) {
		if (molBond.getAtom1().getAtomMap() < molBond.getAtom2().getAtomMap()) {
			return molBond.getAtom1().getAtomMap()+"*"+molBond.getAtom2().getAtomMap();
		} else {
			return molBond.getAtom2().getAtomMap()+"*"+molBond.getAtom1().getAtomMap();
		}		
	}

	private static int getBondMismatch(Molecule rxnMol) {
		HashMap<String,Integer> bondCount = new HashMap<String,Integer>();
		MolBond[] bondArray = rxnMol.getBondArray();
		for (MolBond molBond : bondArray) {
			String bondIndex = getBondIndex(molBond);
			if (!bondCount.containsKey(bondIndex)) {
				bondCount.put(bondIndex, 1);
			} else {
				bondCount.put(bondIndex, bondCount.get(bondIndex)+1);
			}
			molBond.putProperty("INDEX", bondIndex);
		}
		LinkedList<MolBond> bondsToRemove = new LinkedList<MolBond>();
		int returnVal = 0;
		for (MolBond molBond : bondArray) {
			if (bondCount.get(molBond.getProperty("INDEX")) == 1 ) {
				bondsToRemove.add(molBond);
				returnVal++;
			}
		}
		for (MolBond molBond : bondsToRemove) {
			rxnMol.removeBond(molBond);
		}
		return returnVal;
	}

	private static void updateMappings(RxnMolecule rxnMol, HashMap<String, String> indexMapping) {
		MolAtom[] atomArray = rxnMol.getAtomArray();
		// System.out.println(rxnMol.getPropertyCount());
		for (MolAtom molAtom : atomArray) {
			if (indexMapping.containsKey(molAtom.getProperty(molProcess.INDEX).toString())) {
				molAtom.putProperty(molProcess.INDEX,
						indexMapping.get(molAtom.getProperty(molProcess.INDEX).toString()));
			}
			Object[] propertyKeySet = molAtom.propertyKeySet().toArray();
			for (Object key : propertyKeySet) {
				if (!key.toString().equalsIgnoreCase(molProcess.INDEX)) {
					molAtom.removeProperty(key.toString());
					// System.out.println("removing prop");
				}
			}
		}

		MolBond[] bondArray = rxnMol.getBondArray();
		for (MolBond molBond : bondArray) {
			Object[] array = molBond.propertyKeySet().toArray();
			for (Object key : array) {
				molBond.removeProperty(key.toString());
			}
		}
	}

	private static TreeMap<Integer, HashMap<Integer, Integer>> getIndexMapping(String filePath) {
		BufferedReader br = null;
		TreeMap<Integer, HashMap<Integer, Integer>> returnVal = new TreeMap<Integer, HashMap<Integer, Integer>>();
		try {

			String ln;

			br = new BufferedReader(new FileReader(filePath));
			
			while ((ln = br.readLine()) != null) {
				String[] split = ln.replace("\"", "").split("\t");
				if (!returnVal.containsKey(Integer.valueOf(split[0])))
					returnVal.put(Integer.valueOf(split[0]), new HashMap<Integer, Integer>());
				try {
					returnVal.get(Integer.valueOf(split[0])).put(Integer.valueOf(split[2]), Integer.valueOf(split[1]));
				} catch (NumberFormatException e) {
				}
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

		return returnVal;
	}

}
