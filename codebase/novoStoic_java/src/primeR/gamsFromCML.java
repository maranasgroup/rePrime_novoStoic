package primeR;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Sets;

import chemaxon.formats.MolFileHandler;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.io.formats.MoleculeImporter;
import chemaxon.marvin.io.formats.mdl.MolfileUtil;
import chemaxon.marvin.util.MolImportUtil;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.struc.RxnMolecule;
import chemaxon.util.MolImporterUtil;

public class gamsFromCML {

	public static void main(String[] args) {

		String cmlfolderPath = "reactions/";
		File[] listOfCMLs = getCMLs(cmlfolderPath);

		for (File file : listOfCMLs) {
			if (!file.getName().endsWith(".mrv"))
				continue;
			try {
				createGAMSEnv(file);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}

		/*
		 * the bond can be mapped only if both the adjoining atoms are mapped.
		 */

	}

	private static File[] getCMLs(String cmlfolderPath) {
		File cmlFolder = new File(cmlfolderPath);
		if (cmlFolder.isDirectory()) {
			return cmlFolder.listFiles();
		} else
			return null;

	}

	private static void createGAMSEnv(File file) throws FileNotFoundException {
		String basePath = createFolderStructure(file);
		RxnMolecule mol = getRxnMol(file);

		createGAMSfiles(mol, basePath);
	}

	private static String createFolderStructure(File file) {
		String basePath = "GAMS/" + file.getName().replaceFirst("[\\.].+", "");
		try {

			new File("GAMS").mkdir();
			new File(basePath).mkdir();
			new File(basePath + "/sets").mkdir();
			new File(basePath + "/parameters").mkdir();
			new File(basePath + "/mapping").mkdir();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return basePath;

	}

	private static void createGAMSfiles(RxnMolecule RxnMol, String basePath) throws FileNotFoundException {
		print_step_1(RxnMol, basePath);
		print_step_2(RxnMol, basePath);
		print_step_3(RxnMol, basePath);

	}

	private static void print_step_3(RxnMolecule rxnMol, String basePath) throws FileNotFoundException {
		Molecule[] molecules = rxnMol.getReactants();
		PrintWriter pw_a = new PrintWriter(basePath+"/mapping/r.nonTerm");
		PrintWriter pw_b = new PrintWriter(basePath+"/mapping/u.nonTerm");
		for (Molecule mol : molecules) {
			MolAtom[] atomArray = mol.getAtomArray();
			
			
			for (MolAtom molAtom : atomArray) {
				if (molAtom.containsPropertyKey(molProcess.TERMINAL)) {
					if (!(boolean)molAtom.getProperty(molProcess.TERMINAL)) {
						pw_a.println("'"+molAtom.getProperty(molProcess.INDEX)+"'");	
					}
				}
			}
			
			
			MolBond[] bondArray = mol.getBondArray();
			for (MolBond molBond : bondArray) {
				if (molBond.containsPropertyKey(molProcess.TERMINAL)) {
					if (!(boolean)molBond.getProperty(molProcess.TERMINAL)) {
						pw_b.println("'"+getBondIndex(molBond)+"'");
					}					
				}
			}
		}
		
		pw_a.close();
		pw_b.close();
		
		molecules = rxnMol.getProducts();
		pw_a = new PrintWriter(basePath+"/mapping/p.nonTerm");
		pw_b = new PrintWriter(basePath+"/mapping/v.nonTerm");
		for (Molecule mol : molecules) {
			MolAtom[] atomArray = mol.getAtomArray();			
			
			for (MolAtom molAtom : atomArray) {
				if (molAtom.containsPropertyKey(molProcess.TERMINAL)) {
					if (!(boolean)molAtom.getProperty(molProcess.TERMINAL)) {
						pw_a.println("'"+molAtom.getProperty(molProcess.INDEX)+"'");	
					}
				}
			}
			
			
			MolBond[] bondArray = mol.getBondArray();
			for (MolBond molBond : bondArray) {
				if (molBond.containsPropertyKey(molProcess.TERMINAL)) {
					if (!(boolean)molBond.getProperty(molProcess.TERMINAL)) {
						pw_b.println("'"+getBondIndex(molBond)+"'");
					}					
				}
			}
		}
		
		pw_a.close();
		pw_b.close();
		
		
	}

	private static void print_step_2(RxnMolecule rxnMol, String basePath) throws FileNotFoundException {
		String mappingPath = basePath + "/mapping";
		String paramPath = basePath + "/parameters";
		Molecule[] molecules = rxnMol.getReactants();
		HashMultimap<String, MolAtom> rAtomsHash = HashMultimap.create();
		HashMultimap<String, MolBond> uBondsHash = HashMultimap.create();

		PrintWriter pw_a_nbhrs = new PrintWriter(mappingPath + "/r.neighbours");
		PrintWriter pw_b_nbhrs = new PrintWriter(mappingPath + "/u.neighbours");
		for (Molecule mol : molecules) {
			MolAtom[] atomArray = mol.getAtomArray();
			for (MolAtom molAtom : atomArray) {
				rAtomsHash.put(molAtom.getAliasstr(), molAtom);
				/*
				 * if (molAtom.isPseudo()) {
				 * rAtomsHash.put(molAtom.getAliasstr(), molAtom); } else {
				 * rAtomsHash.put(molAtom.getSymbol(), molAtom);
				 * 
				 * }
				 */
			}
			MolBond[] bondArray = mol.getBondArray();
			String neighbours;
			for (MolBond molBond : bondArray) {
				uBondsHash.put(molBond.getProperty(molProcess.PRIME + "0").toString(), molBond);
				neighbours = getNeighbours(molBond, "atoms");
				if (null != neighbours)
					pw_a_nbhrs.println(neighbours);

				neighbours = getNeighbours(molBond, "bonds");
				if (null != neighbours)
					pw_b_nbhrs.println(neighbours);

			}
		}

		/*
		 * instead of getProducts, is this a bug in Chemaxon
		 */

		pw_a_nbhrs.close();
		pw_b_nbhrs.close();

		HashMultimap<String, MolAtom> pAtomsHash = HashMultimap.create();
		HashMultimap<String, MolBond> vBondsHash = HashMultimap.create();

		pw_a_nbhrs = new PrintWriter(mappingPath + "/p.neighbours");
		pw_b_nbhrs = new PrintWriter(mappingPath + "/v.neighbours");

		// molecules = rxnMol.getAgents();
		molecules = rxnMol.getProducts();
		for (Molecule mol : molecules) {
			MolAtom[] atomArray = mol.getAtomArray();
			for (MolAtom molAtom : atomArray) {
				pAtomsHash.put(molAtom.getAliasstr(), molAtom);
				/*
				 * if (molAtom.isPseudo()) {
				 * pAtomsHash.put(molAtom.getAliasstr(), molAtom); } else {
				 * pAtomsHash.put(molAtom.getSymbol(), molAtom); }
				 */
			}
			MolBond[] bondArray = mol.getBondArray();
			String neighbours;
			for (MolBond molBond : bondArray) {
				vBondsHash.put(molBond.getProperty(molProcess.PRIME + "0").toString(), molBond);

				neighbours = getNeighbours(molBond, "atoms");
				if (null != neighbours)
					pw_a_nbhrs.println(neighbours);

				neighbours = getNeighbours(molBond, "bonds");
				if (null != neighbours)
					pw_b_nbhrs.println(neighbours);
			}
		}

		pw_a_nbhrs.close();
		pw_b_nbhrs.close();

		PrintWriter pw = new PrintWriter(mappingPath + "/rp_atoms.map");
		PrintWriter pw_indic = new PrintWriter(paramPath + "/rp.indicator");

		Set<String> atomKeys = rAtomsHash.keySet();
		for (String akeys : atomKeys) {
			Set<MolAtom> rAtoms = rAtomsHash.get(akeys);
			Set<MolAtom> pAtoms = pAtomsHash.get(akeys);
			Set<List<MolAtom>> cartesianProduct = Sets.cartesianProduct(rAtoms, pAtoms);
			for (List<MolAtom> list : cartesianProduct) {
				pw.println("'" + list.get(0).getProperty(molProcess.INDEX).toString() + "'.'"
						+ list.get(1).getProperty(molProcess.INDEX).toString() + "'");
				pw_indic.println(getAtomWeightParam(list.get(0), list.get(1)));
			}
		}
		pw.close();
		pw_indic.close();

		pw = new PrintWriter(mappingPath + "/uv_bonds.map");
		pw_indic = new PrintWriter(paramPath + "/uv.indicator");
		Set<String> bondKeys = uBondsHash.keySet();
		for (String bkeys : bondKeys) {
			Set<MolBond> uBonds = uBondsHash.get(bkeys);
			Set<MolBond> vBonds = vBondsHash.get(bkeys);
			Set<List<MolBond>> cartesianProduct = Sets.cartesianProduct(uBonds, vBonds);
			for (List<MolBond> list : cartesianProduct) {
				pw.println("'" + getBondIndex(list.get(0)) + "'.'" + getBondIndex(list.get(1)) + "'");
				pw_indic.println(getBondWeights(list.get(0), list.get(1)));
			}
		}
		pw.close();
		pw_indic.close();
	}

	private static String getBondWeights(MolBond uBond, MolBond vBond) {
		// TODO Auto-generated method stub
		Set<String> propertyKeySet = uBond.propertyKeySet();
		int cnt = 0;
		for (String propKey : propertyKeySet) {
			if (propKey.contains(molProcess.PRIME)) {
				if (uBond.getProperty(propKey).toString().equalsIgnoreCase(vBond.getProperty(propKey).toString())) {
					cnt++;
				}
			}
		}
		return "'" + getBondIndex(uBond) + "'.'" + getBondIndex(vBond) + "' " + cnt;
	}

	private static String getAtomWeightParam(MolAtom rAtom, MolAtom pAtom) {
		Set<String> propertyKeySet = rAtom.propertyKeySet();
		int cnt = 0;
		for (String propKey : propertyKeySet) {
			if (propKey.contains(molProcess.PRIME)) {
				if (rAtom.getProperty(propKey).toString().equalsIgnoreCase(pAtom.getProperty(propKey).toString())) {
					cnt++;
				}
			}
		}
		return "'" + rAtom.getProperty(molProcess.INDEX).toString() + "'.'"
				+ pAtom.getProperty(molProcess.INDEX).toString() + "' " + cnt;
	}

	private static String getNeighbours(MolBond molBond, String ab) {
		if (ab.equalsIgnoreCase("atoms")) {

			/*
			 * if (molBond.containsPropertyKey(molProcess.IGNOREDBYATOM)) { if
			 * (molBond.getProperty(molProcess.IGNOREDBYATOM).toString()
			 * .equalsIgnoreCase(molBond.getAtom2().getProperty(molProcess.INDEX
			 * ).toString())) { return "'" +
			 * molBond.getAtom1().getProperty(molProcess.INDEX).toString() +
			 * "'.'" +
			 * molBond.getAtom2().getProperty(molProcess.INDEX).toString() +
			 * "'"; } else { return "'" +
			 * molBond.getAtom2().getProperty(molProcess.INDEX).toString() +
			 * "'.'" +
			 * molBond.getAtom1().getProperty(molProcess.INDEX).toString() +
			 * "'"; } } else {
			 */
			return "'" + molBond.getAtom1().getProperty(molProcess.INDEX).toString() + "'.'"
					+ molBond.getAtom2().getProperty(molProcess.INDEX).toString() + "'" + "\n" + "'"
					+ molBond.getAtom2().getProperty(molProcess.INDEX).toString() + "'.'"
					+ molBond.getAtom1().getProperty(molProcess.INDEX).toString() + "'";
			// }
		} else {
			return "'" + getBondIndex(molBond) + "'.'" + molBond.getAtom1().getProperty(molProcess.INDEX).toString()
					+ "'" + "\n" + "'" + getBondIndex(molBond) + "'.'"
					+ molBond.getAtom2().getProperty(molProcess.INDEX).toString() + "'";
		}
	}

	private static void print_step_1(RxnMolecule rxnMol, String basePath) {
		try {
			String setPath = basePath + "/sets";
			String paramPath = basePath + "/parameters";
			String mappingPath = basePath + "/mapping";
			Molecule[] reactants = rxnMol.getReactants();
			PrintWriter pw_considerNode = new PrintWriter(paramPath + "/r_atoms.node");

			PrintWriter pw_a = new PrintWriter(setPath + "/r_atoms.index");

			PrintWriter pw_b = new PrintWriter(setPath + "/u_bonds.index");

			PrintWriter pw_gi = new PrintWriter(setPath + "/rgi.index");
			PrintWriter pw_gi_mapping = new PrintWriter(mappingPath + "/r_gi.map");
			TreeSet<String> gi_printed = new TreeSet<String>();
			String GI_INDEX = rxnMol.getProperty("GI_INDEX");
			String gi_key;
			for (int i = 0; i < reactants.length; i++) {
				MolAtom[] atomArray = reactants[i].getAtomArray();
				for (MolAtom molAtom : atomArray) {
					pw_a.println("'" + molAtom.getProperty(molProcess.INDEX) + "'");
					/*
					 * if (molAtom.isPseudo()) { if
					 * (molAtom.getAliasstr().equalsIgnoreCase("e")) {
					 * pw_considerNode.println("'" +
					 * molAtom.getProperty(molProcess.INDEX) + "' 1"); } } else
					 * { if (!molAtom.getSymbol().equalsIgnoreCase("H")) {
					 * pw_considerNode.println("'" +
					 * molAtom.getProperty(molProcess.INDEX) + "' 1"); }
					 * 
					 * }
					 */
					if (!molAtom.isPseudo()) {
						if (!molAtom.getSymbol().equalsIgnoreCase("H")) {
							pw_considerNode.println("'" + molAtom.getProperty(molProcess.INDEX) + "' 1");
							gi_key = "r_" + molAtom.getProperty(GI_INDEX).toString();
							if (!gi_printed.contains(gi_key)) {
								pw_gi.println("'" + gi_key + "'");
							}
							gi_printed.add(gi_key);
							pw_gi_mapping.println("'" + gi_key + "'.'" + molAtom.getProperty(molProcess.INDEX) + "'");
						}
					} else {
						if (molAtom.getAliasstr().equalsIgnoreCase("e")) {
							pw_considerNode.println("'" + molAtom.getProperty(molProcess.INDEX) + "' 1");
						}
					}

				}

				MolBond[] bondArray = reactants[i].getBondArray();

				for (MolBond molBond : bondArray) {
					getBondIndex(molBond);
					pw_b.println("'" + getBondIndex(molBond) + "'");
				}
			}
			pw_considerNode.close();
			pw_a.close();
			pw_b.close();
			pw_gi.close();
			pw_gi_mapping.close();

			Molecule[] products = rxnMol.getProducts();
			// Molecule[] products = rxnMol.getAgents();

			pw_considerNode = new PrintWriter(paramPath + "/p_atoms.node");

			pw_a = new PrintWriter(setPath + "/p_atoms.index");
			pw_b = new PrintWriter(setPath + "/v_bonds.index");

			pw_gi = new PrintWriter(setPath + "/pgi.index");
			pw_gi_mapping = new PrintWriter(mappingPath + "/p_gi.map");
			gi_printed = new TreeSet<String>();

			for (int i = 0; i < products.length; i++) {
				MolAtom[] atomArray = products[i].getAtomArray();
				for (MolAtom molAtom : atomArray) {
					pw_a.println("'" + molAtom.getProperty(molProcess.INDEX) + "'");
					/*
					 * if (molAtom.isPseudo()) { if
					 * (molAtom.getAliasstr().equalsIgnoreCase("e")) {
					 * pw_considerNode.println("'" +
					 * molAtom.getProperty(molProcess.INDEX) + "' 1"); } } else
					 * { pw_considerNode.println("'" +
					 * molAtom.getProperty(molProcess.INDEX) + "' 1"); }
					 */
					if (!molAtom.isPseudo()) {
						if (!molAtom.getSymbol().equalsIgnoreCase("H")) {

							pw_considerNode.println("'" + molAtom.getProperty(molProcess.INDEX) + "' 1");
							gi_key = "p_" + molAtom.getProperty(GI_INDEX).toString();
							if (!gi_printed.contains(gi_key)) {
								pw_gi.println("'" + gi_key + "'");
							}
							gi_printed.add(gi_key);
							pw_gi_mapping.println("'" + gi_key + "'.'" + molAtom.getProperty(molProcess.INDEX) + "'");
						}
					} else {
						if (molAtom.getAliasstr().equalsIgnoreCase("e")) {
							pw_considerNode.println("'" + molAtom.getProperty(molProcess.INDEX) + "' 1");
						}
					}
				}
				MolBond[] bondArray = products[i].getBondArray();
				for (MolBond molBond : bondArray) {
					pw_b.println("'" + getBondIndex(molBond) + "'");
				}
			}
			pw_considerNode.close();
			pw_a.close();
			pw_b.close();
			pw_gi.close();
			pw_gi_mapping.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	private static String getBondIndex(MolBond molBond) {
		TreeSet<String> sorter = new TreeSet<String>();
		sorter.add(molBond.getAtom1().getProperty(molProcess.INDEX).toString());
		sorter.add(molBond.getAtom2().getProperty(molProcess.INDEX).toString());
		return sorter.toString();
	}

	public static RxnMolecule getRxnMol(File file) {
		try {
			RxnMolecule reaction = new RxnMolecule().getReaction(new MolImporter(file).read());
			// reaction.rebuildStructures();
			System.out.println(reaction.getReactantCount() + " " + reaction.getProductCount());
			String indexOfFinalGI = getGIKey(reaction);
			reaction.setProperty("GI_INDEX", indexOfFinalGI);
			return reaction;

		} catch (IOException e) {

			e.printStackTrace();
		}
		return null;

	}

	private static String getGIKey(RxnMolecule reaction) {
		Set<String> propertyKeySet = reaction.getAtom(0).propertyKeySet();
		TreeMap<Integer, String> sorter = new TreeMap<Integer, String>();
		for (String key : propertyKeySet) {
			if (key.contains(molProcess.PRIME)) {
				sorter.put(Integer.valueOf(key.split("[_]")[1]), key);
			}
		}
		return sorter.get(sorter.lastKey());
	}

}
