package nonEssential;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

import org.jgrapht.alg.BellmanFordShortestPath;
import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.alg.FloydWarshallShortestPaths;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleDirectedWeightedGraph;
import org.jgrapht.graph.SimpleWeightedGraph;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.LinkedHashMultimap;
import com.google.common.collect.LinkedListMultimap;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.sketch.swing.actions.NewAction;
import chemaxon.struc.Molecule;

public class readSolutionFile {

	private static final String STEP_0 = "Please enter the file name";
	private static final String SOURCE = "src";
	private static final String DIRECTION = "drxn";
	private static final String TYPE = "type";
	private static final String SUBSTRATES = "s";
	private static final String PRODUCTS = "p";
	private static final String RXNSET = "rxns";
	private static final String MOIETIES = "mtes";
	private static final String MOLECULE = "str";
	private static final String STEP_1 = "Please enter the substrate name";
	// private static Scanner scanIn;
	private static String screenInput;
	private static String targetMoleculeId;
	private static LinkedListMultimap<String, String> names;
	private static String[] split;
	private static HashMap<String, String> metabTypes = new HashMap<String, String>();
	private static Connection conn = Connect.getConnection();
	private static Statement st;
	private static LinkedHashMultimap<String, String> rxnNamesLinkedMap;
	private static LinkedHashMultimap<Object, Object> metabNamesLinkedMap;
	private static SimpleDirectedWeightedGraph<String, DefaultWeightedEdge> wGraph;
	private static HashBasedTable<String, String, Object> rxnTbl;
	private static HashBasedTable<String, String, Object> ruleTbl;
	private static HashBasedTable<String, String, Object> metTbl;
	private static HashMultimap<String, String> edgeNames;
	private static BufferedReader input;;

	public static void main(String[] args) {
		/* Ask which record to start parsing from, or which molecule to read */
		/*
		 * 1) print the name and smiles of the molecule2) print the hypothetical conversion for this molecules3) find the
		 * closest match and print that reaction name and smiles4) Start retro, print the target, using smilarity to
		 * identify the molecule that produced it, 4.1) Show the reaction and name and identify the most similar molecule
		 * to print, ask if that one is correct, else get the correct number to trace 4.2) Remove the visited reaction5)
		 * Terminate or ask to continue to next molecule
		 */

		/* Add data will be stored in a graph */
		/*
		 * get the rxn ids and rules from the file, get the substrate and product and add edges to them. Based on the
		 * moieties, calculate the similarity, also store mols in the node, names in the node and edges
		 */
		input = new BufferedReader(new InputStreamReader(System.in));
		metabTypes.put("ex_wl", "rule_l");
		/* metabTypes.put("ex_ul", "rule_l"); */
		metabTypes.put("intermediate_l", "rule_l");

		metabTypes.put("intermediate_r", "rule_r");
		metabTypes.put("ex_wr", "rule_r");
		/* metabTypes.put("ex_ur", "rule_l"); */

		metabTypes.put("intermediate_iso", "rule_iso");

		File inputFile = (File) prompt(0);
		populateDataStructures(inputFile);
		findPathways();
		prompt(1);

	}

	private static void findPathways() {
		
		String substrate = (String) prompt(1);
		String product = (String) prompt(2);
		List<DefaultWeightedEdge> path = BellmanFordShortestPath.findPathBetween(wGraph, substrate, product); 
		DijkstraShortestPath.findPathBetween(wGraph, substrate, product);
				Set<String> rowKeySet = metTbl.rowKeySet();
		for (String metId : rowKeySet) {
			if (metId.contains(targetMoleculeId)) {
				targetMoleculeId = metId;
			}
		}
		String edgeName;
		double edgeWeight;
		while (null != targetMoleculeId) {
			try {
				System.out.println("target now is "+ targetMoleculeId);
				Set<DefaultWeightedEdge> incomingEdgesOf = wGraph.incomingEdgesOf(targetMoleculeId);
				
				for (DefaultWeightedEdge edge : incomingEdgesOf) {
					edgeName = edge.toString().replaceAll("[\\(\\)]", "");
					edgeWeight = wGraph.getEdgeWeight(edge);
					edgeName = edgeName.split("[ ][:][ ]")[0] + "->" + edgeName.split("[ ][:][ ]")[1];					
					System.out.println(edge + " = " + edgeWeight + ", rxns are "+ edgeNames.get(edgeName));
				}
				/* please enter the most likely substrate */
				/* print the reactions that contains the target and product, name and smiles */
				substrate = (String) prompt(1);
				System.out.println("the SMILES for this metabolite is = " + MolExporter.exportToFormat((Molecule) metTbl.get(substrate, MOLECULE), "smiles"));
				//System.out.println("reactions for this substrate/product is = " + edgeNames.get(substrate + "->" + targetMoleculeId));
				targetMoleculeId = substrate;
			}
			catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}

		/*
		 * String type; for (String start : rowKeySet) { type = (String) metTbl.get(start, TYPE); if
		 * (type.equalsIgnoreCase("ex_ul")||type.equalsIgnoreCase("ex_wl")) {
		 * 
		 * try { List<DefaultWeightedEdge> path = BellmanFordShortestPath.findPathBetween(wGraph,start, targetMoleculeId);
		 * System.out.println(path.toString()); } catch (Exception e) { // TODO: handle exception }
		 * 
		 * } }
		 */

	}

	private static void populateDataStructures(File inputFile) {
		targetMoleculeId = inputFile.getName().replace(".tsv", "").split("[_]")[1];
		rxnTbl = HashBasedTable.create();
		ruleTbl = HashBasedTable.create();
		metTbl = HashBasedTable.create();
		HashSet<String> listOfMetabolites = new HashSet<String>();
		try {
			BufferedReader inputStream = new BufferedReader(new FileReader(inputFile));
			String inLine;
			String id;
			String source;
			String direction;
			while ((inLine = inputStream.readLine()) != null) {
				split = inLine.split("\t");

				if (split[3].startsWith("rxn")) {
					id = split[0];
					source = split[0].split("[_]|[#]")[0];

					if (split[2].equalsIgnoreCase("1")) {
						direction = "fwd";
					}
					else {
						direction = "bck";
					}
					rxnTbl.put(id, SOURCE, source);
					rxnTbl.put(id, DIRECTION, direction);
					rxnTbl.put(id, TYPE, split[3]);
					rxnTbl.put(id, SUBSTRATES, new HashSet<String>());
					rxnTbl.put(id, PRODUCTS, new HashSet<String>());
				}
				else if (split[3].startsWith("rule")) {
					id = split[0];
					source = split[0].split("[_]|[#]")[0];

					if (split[1].equalsIgnoreCase("1")) {
						direction = "fwd";
					}
					else {
						direction = "bck";
					}
					ruleTbl.put(id, SOURCE, source);
					ruleTbl.put(id, DIRECTION, direction);
					ruleTbl.put(id, TYPE, split[3]);
					ruleTbl.put(id, RXNSET, new HashSet<String>());
					ruleTbl.put(id, SUBSTRATES, new HashSet<String>());
					ruleTbl.put(id, PRODUCTS, new HashSet<String>());
				}
				else {
					id = split[0];
					source = split[0].split("[:]")[0];
					source = getSourceId(source);
					if (split[1].equalsIgnoreCase("1")) {
						direction = PRODUCTS;
					}
					else {
						direction = SUBSTRATES;
					}
					metTbl.put(id, SOURCE, source);
					metTbl.put(id, DIRECTION, direction);
					metTbl.put(id, TYPE, split[3]);
					metTbl.put(id, MOIETIES, new HashMap<String, Integer>());
				}
			}
			inputStream.close();

		}
		catch (IOException e) {
			System.err.println("IOException:");
			e.printStackTrace();
		}

		updateRxnTable(rxnTbl, metTbl);
		updateRuleTable(ruleTbl, metTbl);
		updateMetTable(metTbl);
		updateNames(rxnTbl, ruleTbl, metTbl);

		createGraph(rxnTbl, ruleTbl, metTbl);

	}

	private static void createGraph(HashBasedTable<String, String, Object> rxnTbl, HashBasedTable<String, String, Object> ruleTbl, HashBasedTable<String, String, Object> metTbl) {

		wGraph = new SimpleDirectedWeightedGraph<String, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		Set<String> rowKeySet = metTbl.rowKeySet();
		for (String metabId : rowKeySet) {
			wGraph.addVertex(metabId);
		}

		rowKeySet = rxnTbl.rowKeySet();
		edgeNames = HashMultimap.create();
		for (String rxnId : rowKeySet) {
			HashSet<String> substrates = (HashSet<String>) rxnTbl.get(rxnId, SUBSTRATES);
			HashSet<String> products = (HashSet<String>) rxnTbl.get(rxnId, PRODUCTS);
			for (String s : substrates) {
				for (String p : products) {
					wGraph.addEdge(s, p);
					DefaultWeightedEdge edge = wGraph.getEdge(s, p);
					wGraph.setEdgeWeight(edge, getEdgeWeight(s, p, metTbl));
					edgeNames.put(s + "->" + p, rxnId);
				}
			}
		}
		rowKeySet = ruleTbl.rowKeySet();
		for (String ruleId : rowKeySet) {
			HashSet<String> substrates = (HashSet<String>) ruleTbl.get(ruleId, SUBSTRATES);
			HashSet<String> products = (HashSet<String>) ruleTbl.get(ruleId, PRODUCTS);
			for (String s : substrates) {
				for (String p : products) {
					wGraph.addEdge(s, p);
					DefaultWeightedEdge edge = wGraph.getEdge(s, p);
					wGraph.setEdgeWeight(edge, getEdgeWeight(s, p, metTbl));
					edgeNames.put(s + "->" + p, "rule:" + ruleId);
				}
			}
		}
	}

	private static double getEdgeWeight(String s, String p, HashBasedTable<String, String, Object> metTbl) {
		HashMap<String, Integer> sMty = (HashMap<String, Integer>) metTbl.get(s, MOIETIES);
		HashMap<String, Integer> pMty = (HashMap<String, Integer>) metTbl.get(p, MOIETIES);
		double totalS = 0;
		Collection<Integer> values = sMty.values();
		for (Integer cnt : values) {
			totalS = totalS + cnt;
		}
		double totalP = 0;
		values = pMty.values();
		for (Integer cnt : values) {
			totalP = totalP + cnt;
		}

		double similarity = 0;
		Set<String> moietySbst = sMty.keySet();
		Set<String> moietyPrd = pMty.keySet();
		for (String m : moietySbst) {
			if (!moietyPrd.contains(m)) continue;
			similarity = similarity + 2 * Math.min(sMty.get(m), pMty.get(m));
		}

		return similarity / (totalS + totalP - similarity); /* change to dissimilarity */
	}

	private static void updateNames(HashBasedTable<String, String, Object> rxnTbl, HashBasedTable<String, String, Object> ruleTbl, HashBasedTable<String, String, Object> metTbl) {
		rxnNamesLinkedMap = LinkedHashMultimap.create();
		metabNamesLinkedMap = LinkedHashMultimap.create();
		try {
			st = conn.createStatement();
			String query = "";
			Set<String> rowKeySet = rxnTbl.rowKeySet();
			String source;
			for (String rxnId : rowKeySet) {
				source = rxnId.split("[_]|[#]")[0];
				rxnId = rxnId.split("[_]|[#]")[1];
				query = "SELECT name FROM metrxn0200.reactionnames where id = '" + rxnId + "' and source = '" + source + "' order by length(name) asc";
				ResultSet rs = st.executeQuery(query);
				String rxnName;
				while (rs.next()) {
					rxnName = rs.getString("name");
					rxnNamesLinkedMap.put(rxnId, rxnName);
				}
			}
			rowKeySet = ruleTbl.rowKeySet();
			for (String ruleId : rowKeySet) {
				HashSet<String> rowKeySet2 = (HashSet<String>) ruleTbl.get(ruleId, RXNSET);
				for (String rxnId : rowKeySet2) {
					source = rxnId.split("[_]|[#]")[0];
					rxnId = rxnId.split("[_]|[#]")[1];
					query = "SELECT name FROM metrxn0200.reactionnames where id = '" + rxnId + "' and source = '" + source + "' order by length(name) asc";
					ResultSet rs = st.executeQuery(query);
					String rxnName;
					while (rs.next()) {
						rxnName = rs.getString("name");
						rxnNamesLinkedMap.put(rxnId, rxnName);
					}
				}
			}

			rowKeySet = metTbl.rowKeySet();
			for (String metId : rowKeySet) {
				source = getSourceId(metId.split("[:]")[0]);
				metId = metId.split("[:]")[2];
				query = "SELECT name FROM metrxn0200.names where id = '" + metId + "' and source = '" + source + "' order by length(name) asc";
				ResultSet rs = st.executeQuery(query);
				String metName;
				while (rs.next()) {
					metName = rs.getString("name");
					metabNamesLinkedMap.put(metId, metName);
				}
			}
			st.close();
		}
		catch (SQLException ex) {
			System.err.println(ex.getMessage());
		}

	}

	private static void updateMetTable(HashBasedTable<String, String, Object> metTbl) {
		Set<String> rowKeySet = metTbl.rowKeySet();
		/* update moieties */
		try {
			BufferedReader inputStream = new BufferedReader(new FileReader("parameters/atoms.formula"));
			String inLine;
			String moietyId;
			Integer mCount;
			String metabId;
			while ((inLine = inputStream.readLine()) != null) {
				String[] split2 = inLine.split("[']|[.]|[ ]");
				metabId = split2[4];
				if (rowKeySet.contains(metabId)) {
					moietyId = inLine.split("[']|[.]|[ ]")[1];
					mCount = Integer.valueOf(inLine.split("[ ]")[1]);
					HashMap<String, Integer> moietyMap = (HashMap<String, Integer>) metTbl.get(metabId, MOIETIES);
					moietyMap.put(moietyId, mCount);
				}
			}
			inputStream.close();
		}
		catch (IOException e) {
			System.err.println("IOException:");
			e.printStackTrace();
		}

		String source;
		String id;
		String strType;
		Molecule mol;
		for (String metabId : rowKeySet) {
			source = getSourceId(metabId.split("[:]")[0]);
			strType = metabId.split("[:]")[1];
			id = metabId.split("[:]")[2];
			if (strType.equalsIgnoreCase("mol")) {
				mol = getMoleculeFromDirectory(id, source, strType);
			}
			else {
				mol = getMoleculeFromDb(id, source, strType);
			}
			metTbl.put(metabId, MOLECULE, mol);
		}
	}

	private static Molecule getMoleculeFromDb(String id, String source, String strType) {
		try {
			st = conn.createStatement();
			String query = "SELECT structure `str` from `megatable`.`structures` where source = '" + source + "' and id = '" + id + "' and type = '" + strType + "' limit 1";
			ResultSet rs = st.executeQuery(query);
			String str = null;
			while (rs.next()) {
				str = rs.getString("str");
			}
			st.close();
			return MolImporter.importMol(str);
		}
		catch (SQLException ex) {
			System.err.println(ex.getMessage());
		}
		catch (MolFormatException e) {

			e.printStackTrace();
		}
		return null;
	}

	private static Molecule getMoleculeFromDirectory(String id, String source, String strType) {

		try {
			Scanner scanner = new Scanner(new File("structure/" + strType + "/" + source + "/" + id + "." + strType));
			String content = scanner.useDelimiter("\\A").next();
			scanner.close();
			return MolImporter.importMol(content);
		}
		catch (FileNotFoundException e) {

			e.printStackTrace();
		}
		catch (MolFormatException e) {
			e.printStackTrace();
		}
		return null;
	}

	private static void updateRuleTable(HashBasedTable<String, String, Object> ruleTbl, HashBasedTable<String, String, Object> metTbl) {
		try {
			BufferedReader inputStream = new BufferedReader(new FileReader("parameters/ruleToRid.map"));
			String inLine;
			String ruleId;
			String rid;
			Set<String> rowKeySet = ruleTbl.rowKeySet();
			String ruleType;
			while ((inLine = inputStream.readLine()) != null) {
				ruleId = inLine.split("[']|[.]")[1];
				if (!rowKeySet.contains(ruleId)) continue;
				HashSet<String> rxns = (HashSet<String>) ruleTbl.get(ruleId, RXNSET);
				rxns.add(inLine.split("[']|[.]")[4]);
			}

			Set<String> metabKeySet = metTbl.rowKeySet();
			String metType;
			Set<String> metTypekeySet = metabTypes.keySet();
			String direction;
			for (String metabId : metabKeySet) {
				metType = (String) metTbl.get(metabId, TYPE);
				if (metTypekeySet.contains(metType)) {
					for (String rule : rowKeySet) {
						ruleType = (String) ruleTbl.get(rule, TYPE);
						if (!metabTypes.get(metType).equalsIgnoreCase(ruleType)) continue;

						direction = (String) metTbl.get(metabId, DIRECTION);
						if (null == direction) continue;
						if (direction.equalsIgnoreCase(SUBSTRATES)) {
							HashSet<String> substrates = (HashSet<String>) ruleTbl.get(rule, SUBSTRATES);
							substrates.add(metabId);
						}
						if (direction.equalsIgnoreCase(PRODUCTS)) {
							HashSet<String> products = (HashSet<String>) ruleTbl.get(rule, PRODUCTS);
							products.add(metabId);
						}
					}
				}
			}
			inputStream.close();
		}
		catch (IOException e) {
			System.err.println("IOException:");
			e.printStackTrace();
		}

	}

	private static void updateRxnTable(HashBasedTable<String, String, Object> rxnTbl, HashBasedTable<String, String, Object> metTbl) {
		Set<String> rowKeySet = rxnTbl.rowKeySet();
		try {
			BufferedReader inputStream = new BufferedReader(new FileReader("parameters/rule.scj"));
			String inLine;
			String rid;
			String met;
			String source;
			while ((inLine = inputStream.readLine()) != null) {
				rid = inLine.split("[\\']|[\\.]|[ ]")[4];
				// System.out.println(inLine.split("[ ]")[1]);
				int scj = Integer.valueOf(inLine.split("[ ]")[1]);
				if (!rowKeySet.contains(rid)) continue;
				if (rxnTbl.get(rid, DIRECTION).equals("bck")) {
					scj = scj * -1;
				}
				met = inLine.split("[']|[.]|[ ]")[1];
				if (scj < 0) {
					HashSet<String> substrates = (HashSet<String>) rxnTbl.get(rid, SUBSTRATES);
					substrates.add(met);
				}
				else {
					HashSet<String> products = (HashSet<String>) rxnTbl.get(rid, PRODUCTS);
					products.add(met);
				}
				if (!metTbl.containsRow(met)) {
					source = met.split("[:]")[0];
					source = getSourceId(source);
					metTbl.put(met, SOURCE, source);
					metTbl.put(met, TYPE, "met");
					metTbl.put(met, MOIETIES, new HashMap<String, Integer>());
				}

			}
			inputStream.close();
		}
		catch (IOException e) {
			System.err.println("IOException:");
			e.printStackTrace();
		}

	}

	private static String getSourceId(String source) {
		switch (source) {
		case "br":
			return "brenda";
		case "ci":
			return "chebi";
		case "cl":
			return "chembl";
		case "eb":
			return "ecmdb";
		case "ey":
			return "ecocyc";
		case "hm":
			return "hmdb";
		case "kg":
			return "kegg";
		case "my":
			return "metacyc";
		default:
			return source;
		}
	}

	private static Object prompt(int step_id) {
		try {
			if (step_id == 0) {
				System.out.println(STEP_0);
				screenInput = input.readLine();
				return new File(screenInput);
			}
			else if (step_id == 1) {
				System.out.println(STEP_1);
				return input.readLine();
			}
		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

}
