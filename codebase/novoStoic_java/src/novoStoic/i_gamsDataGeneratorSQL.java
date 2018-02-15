package paper6;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.Semaphore;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;






import com.google.common.collect.Sets;
import com.google.common.collect.TreeMultimap;
import com.google.common.hash.Hashing;

public class i_gamsDataGeneratorSQL {
	/* major changes to make in this files and gams.. similarity and dissimilarity in the file.ext name */
	private static final int THREADCOUNT = 7;

	private static int loop = 3;
	private static String type = "stereo";

	private static Connection conn;
	private static String source_qualifier = "NOT IN";
	private static String sources = "('CHEBI')";

	private static SortedMap<String, String> idLookup;
	private static String path;
	private static PrintWriter metabIdMappingStream;
	private static SortedMap<String, SortedSet<String>> moietyhash_id;
	private static SortedMap<String, String> moietyHash_metabHash;
	private static PrintWriter atomRulesStream;
	private static PrintWriter bondRulesStream;
	private static PrintWriter metabStream;
	private static PrintWriter atomStream;
	private static PrintWriter bondStream;
	private static String mapFolder;
	private static String paramFolder;
	private static String setFolder;
	private static String exchangeFolder;
	private static SortedSet<String> atomSet;
	private static SortedSet<String> bondSet;
	private static SortedSet<String> cpdSet;

	private static PrintWriter moietyMappingStream;
	// private static HashMap<String, String> equationMd5;
	private static SortedMap<String, SortedSet<String>> metabhash_id;
	private static PrintWriter SCJStream;
	private static PrintWriter SAJStream;
	private static PrintWriter SBJStream;
	private static HashMap<String, String> rulesMap;
	private static PrintWriter cpdSetStream;
	private static List<String> exchangeMetabolites;
	private static List<String> exchangeFlux;
	private static PrintWriter ruleIdStream;
	private static String level;
	private static Map<String, Integer> ridPrinted;
	private static boolean printAll = true;
	private static boolean printByMoietyHash = false;/* reactions are being printed twice, fix this */

	private static boolean printforBinary = true;/* reactions are being printed twice, fix this */

	private static String sourceName;

	private static HashMap<String, HashMap<String, Integer>> exchangeAtomFrequencyMap;
	private static HashMap<String, HashMap<String, Integer>> exchangeBondFrequencyMap;

	private static Map<String, HashMap<String, Integer>> cpdAtomFrequencyMap;
	private static Map<String, HashMap<String, Integer>> cpdBondFrequencyMap;

	private static SortedMap<String, String> cpdNIdHash;

	private static Set<String> toPrintSimilarity;

	private static SortedSet<String> cpdRid;

	private static Map<String, Double> cpdSimilarity;

	private static HashMap<String, Integer> rulesRidCnt;

	private static PrintWriter SCJStream_rev;

	private static PrintWriter SAJStream_rev;

	private static PrintWriter SBJStream_rev;

	private static PrintWriter cpdSetStream_rev;

	private static PrintWriter ruleIdStream_rev;

	private static PrintWriter rulePairsStream;

	private static PrintWriter ex_scjStream_rev;

	private static PrintWriter ex_setCpdStream_rev;

	private static PrintWriter cpdSimilarityStream_rev;

	private static PrintWriter ruleSimilarityStream_rev;

	private static PrintWriter ex_pairsStream;

	private static SortedSet<String> metabIdSet;

	private static PrintWriter SCJSetStream;

	private static PrintWriter SAJSetStream;

	private static PrintWriter SBJSetStream;

	private static PrintWriter SAJSetStream_rev;

	private static PrintWriter SCJSetStream_rev;

	private static PrintWriter SBJSetStream_rev;

	private static PrintWriter atomRulesSetStream;

	private static PrintWriter bondRulesSetStream;

	private static PrintWriter ex_scjSetStream;

	private static PrintWriter ex_scjSetStream_rev;

	private static PrintWriter metabIdPriceStream;

	private static PrintWriter metabIdCarbonStream;

	private static PrintWriter ruleIdEcStream;
	private static PrintWriter ruleIdEcStream_rev;

	private static TreeMap<String, HashSet<String>> rxnEc;

	// private static SortedSet<String> rulePrinted;

	public static void main(String[] args) {

		{
			/*
			 * Things to do, using, for the conversion that we want, are there other metabolites that are present in the
			 * same way, can we fix those and piggy back on those reactions. i.e .. The conversion we want will be a rule,
			 * now, what is the combination of known metabolites that will form the same rule.
			 */
			// comments

			{
				/*
				 * run these ddls before you start the whole process, set GROUP_CONCAT_MAX_LEN = 999999999; set GLOBAL
				 * max_allowed_packet=1073741824; set global sql_big_selects = 1; (so that the optimizer aborts soon) set
				 * global sql_buffer_result = 1;(so that the data is put into the temporary tables early and results are
				 * buffered sooner) set global max_heap_table_size = 3221225472; set global tmp_table_size = 4294967294; set
				 * global join_buffer_size = 3221225472;(apparently this is the max on stupid windows) set global
				 * sort_buffer_size = 4294967294; set global myisam_sort_buffer_size = 4294967294; set global
				 * myisam_max_sort_file_size = 3221225472; set global myisam_repair_threads = 7;
				 */
			}

			/* start by asking what you want to do, print the gams file or print the exchanges, */

			/* print gams file */
			/*
			 * use the equations which are marked as true, parse them and create the sij.
			 */

			/* get the moieties from the table and make a set out of them */
			/* get the metabolites as different sij file */

			/* 1) for each a, sum(j, R(a,j)*V(j)) =e= 0; */
			/* 2) for each i, sum(j, S(i,j)*V(j)) =e= 0; */

			/* the sets are a, i and j */

			/*
			 * for the set i and set a and metaboliteRules select A.`hash`,B.moiety_hash, group_concat(distinct
			 * concat(source,':',A.`id`,':') order by A.`id`) `ideez`, B.`TF`, B.`primes`, B.`atom_bond` from
			 * metrxn0200.gi_stereo A, `mtable`.`stereo_three` B Where A.source = B.source and A.hash = B.hash and A.id =
			 * B.id and A.source not in ('chembl') and B.source not in ('chembl') group by A.`hash`,B.moiety_hash; we take
			 * the first id to print for now, until the metrxn id's are assigned to them, this will be mutlithreaded
			 */

			/*
			 * next we need to get the R(a,j) and S(i,j)
			 * 
			 * 
			 * 1) for R(a,j) we use the rules table join on the megatable and the reaction rulestable,
			 * 
			 * 2) for S(i,j) we use the balanceInfo table and parse them based on the rules above. 3) two types of S(i,j)
			 * need to be printed, normalized by moiety_hash and just hash
			 */

			/*
			 * exchange files, one for moiety exchange and one for cpd exchange, the cpds only from the reaction will be
			 * printed.
			 */

			/*
			 * then call another proc to print the desired exchanges This will allow me to seperate the files into three
			 * gdx, one for data loading, one for exchange and one set for the solving.
			 */

			/*
			 * we also need to read from mtable.connection_zero, this contains all the formula that will be needed to find
			 * balanced reactions to design starting and ending points in a pathway
			 */

			/* eventually, we would also need directions and predicted kcat/km by the reaction rules */

			/*
			 * use emolecules to get the cost of each molecule, to start with we can randomly assign the cost for each
			 * molecule and test
			 */

		}

		String OS = System.getProperty("os.name");
		if (OS.startsWith("Windows")) {
			path = "D:/paper 6/alchemist/";
			mapFolder = "mapping/";
			paramFolder = "parameters/";
			setFolder = "sets/";
			exchangeFolder = "exchange/";
		}
		else {
			path = "/home/azk172/fromWindows/Akhil/D/paper 6/alchemist/";
			mapFolder = "mapping/";
			paramFolder = "parameters/";
			setFolder = "sets/";
			exchangeFolder = "exchange/";
		}

		if (loop == 0) {
			level = "zero";
		}
		else if (loop == 1) {
			level = "one";
		}
		else if (loop == 2) {
			level = "two";
		}
		else {
			level = "three";
		}

		// exchangeMetabolites = Arrays.asList(new String[] { "18383", "30831", "4877", "16526", "15377" });
		// sourceName = "chebi";
		// exchangeFlux = Arrays.asList(new String[] { "-100", "-200", "100", "100", "400" });

		/*
		 * exchangeMetabolites = Arrays.asList(new String[] { "30953", "16134", "17790", "45285"}); sourceName = "chebi";
		 * exchangeFlux = Arrays.asList(new String[] { "-100", "-100", "-100", "100"});
		 */

//		exchangeMetabolites = Arrays.asList(new String[] { "C00031", "C02752", "C00001" }); /*TAL*/
		exchangeMetabolites = Arrays.asList(new String[] { "HMDB04629", "HMDB03312"});
		sourceName = "HMDB";
		exchangeFlux = Arrays.asList(new String[] { "1", "-1"});

		try {
			conn = Connect.getConnection();

			changeDbSessionParameters();

			if (printAll) {
				metaboliteAndRules();
				reactionsAndRules();
				printExchange(true);
				printSimilarity();
			}
			else {
				printExchange(false);
				printSimilarity();
			}

			conn.close();
		}
		catch (SQLException e) {

			e.printStackTrace();
		}

	}

	private static void printSimilarity() {
		/*
		 * for each metabolite in the exchange, find the similarity to other reactions and also find the similarity to
		 * other metabolites cpd.
		 */
		/* exchange to cpd similarity in HashMap<String,Integer> ex_cpd, */
		try {
			Set<String> keySetNid = cpdNIdHash.keySet();
			toPrintSimilarity = Collections.synchronizedSet(new HashSet<String>());
			cpdSimilarity = Collections.synchronizedMap(new HashMap<String, Double>());

			Set<String> keySetE = exchangeAtomFrequencyMap.keySet();
			final HashMap<String, Integer> network = new HashMap<String, Integer>();
			for (String eMid : keySetE) {
				HashMap<String, Integer> frequencyMapE = exchangeAtomFrequencyMap.get(eMid);
				Set<String> keySet = frequencyMapE.keySet();
				for (String m : keySet) {
					if (!network.containsKey(m)) {
						network.put(m, 0);
					}
					network.put(m, network.get(m) + frequencyMapE.get(m));
				}
			}

			keySetE = exchangeBondFrequencyMap.keySet();

			for (String eMid : keySetE) {
				HashMap<String, Integer> frequencyMapE = exchangeBondFrequencyMap.get(eMid);
				Set<String> keySet = frequencyMapE.keySet();
				for (String m : keySet) {
					if (!network.containsKey(m)) {
						network.put(m, 0);
					}
					network.put(m, network.get(m) + frequencyMapE.get(m));
				}
			}

			final Semaphore s = new Semaphore(THREADCOUNT);
			for (String nid : keySetNid) {
				s.acquireUninterruptibly();
				final String f_nid = nid;
				new Thread() {
					public void run() {
						try {
							Integer min_diff = 100000;
							String m_hash;
							HashMap<String, Integer> frequencyMapA;
							HashMap<String, Integer> frequencyMapB;

							m_hash = cpdNIdHash.get(f_nid);
							frequencyMapA = cpdAtomFrequencyMap.get(m_hash);
							frequencyMapB = cpdBondFrequencyMap.get(m_hash);

							Integer totalCpd = 0;
							Collection<Integer> valuesCpdF = frequencyMapA.values();
							for (Integer sum : valuesCpdF) {
								totalCpd = totalCpd + sum;
							}
							valuesCpdF = frequencyMapB.values();
							for (Integer sum : valuesCpdF) {
								totalCpd = totalCpd + sum;
							}

							Integer difference = 0;
							Integer common = 0;
							Integer totalE = 0;
							Collection<Integer> values = network.values();
							for (Integer sum : values) {
								totalE = totalE + sum;
							}

							Set<String> keySet = network.keySet();
							for (String m : keySet) {
								if (frequencyMapA.containsKey(m)) {
									common = common + 2 * Math.min(frequencyMapA.get(m), network.get(m));
								}
								if (frequencyMapB.containsKey(m)) {
									common = common + 2 * Math.min(frequencyMapB.get(m), network.get(m));
								}
							}

							difference = totalCpd + totalE - common;
							/* choose the least difference */

							if (totalCpd <= 3) {
								difference = 1;
							}

							toPrintSimilarity.add("'" + f_nid + "_x' " + (Math.pow(difference, 2) + 1));
							cpdSimilarity.put("'" + f_nid + "'", (Math.pow(difference, 2) + 1));
						}
						finally {
							s.release();
						}
					}
				}.start();

			}
			/* pause for all threads to close */
			while (s.availablePermits() < THREADCOUNT) {
			}

			PrintWriter cpdSimilarityStream = new PrintWriter(new FileWriter(path + paramFolder + exchangeFolder + "rule_cpd_ex.NetworkSimilarity"));
			PrintWriter ruleSimilarityStream = new PrintWriter(new FileWriter(path + paramFolder + "rules.NetworkSimilarity"));

			if (printforBinary) {
				cpdSimilarityStream_rev = new PrintWriter(new FileWriter(path + paramFolder + exchangeFolder + "rule_cpd_ex_rev.NetworkSimilarity"));
				ruleSimilarityStream_rev = new PrintWriter(new FileWriter(path + paramFolder + "rules_rev.NetworkSimilarity"));
			}

			for (String pr : toPrintSimilarity) {
				cpdSimilarityStream.println(pr);
				if (printforBinary) {
					cpdSimilarityStream_rev.println(pr.replaceFirst("_x'", "_r_x'"));
				}
			}

			HashMap<String, Double> ridSimilarity = new HashMap<String, Double>();

			String rid;

			Double similarity;
			for (String rcpd : cpdRid) {
				rid = "'" + rcpd.split("['][\\.][']")[1];
				if (!ridSimilarity.containsKey(rid)) {
					ridSimilarity.put(rid, 0.00);
				}
				similarity = cpdSimilarity.get(rcpd.split("['][\\.][']")[0] + "'");

				if (ridSimilarity.get(rid) < similarity) {
					ridSimilarity.put(rid, similarity);
				}

			}

			/* each reaction has a its cpd, find the one with least similarity ?? */
			Set<String> keySet = ridSimilarity.keySet();
			for (String r : keySet) {
				ruleSimilarityStream.println(r + " " + ridSimilarity.get(r));
				if (printforBinary) {
					ruleSimilarityStream_rev.println(StringUtils.removeEnd(r, "'") + "_r' " + ridSimilarity.get(r));
				}
			}

			cpdSimilarityStream.close();
			ruleSimilarityStream.close();

			if (printforBinary) {
				cpdSimilarityStream_rev.close();
				ruleSimilarityStream_rev.close();
			}

		}
		catch (IOException e) {
			System.err.println("IOException:");
			e.printStackTrace();
		}

	}

	private static void printExchange(boolean exchangeOfnwIO) {

		HashMap<String, String> exchangefluxLoookup = new HashMap<String, String>();
		for (int i = 0; i < exchangeMetabolites.size(); i++) {
			exchangefluxLoookup.put(exchangeMetabolites.get(i) + "_x", exchangeFlux.get(i));
		}

		exchangeAtomFrequencyMap = new HashMap<String, HashMap<String, Integer>>();
		exchangeBondFrequencyMap = new HashMap<String, HashMap<String, Integer>>();

		String exchange_met_str = "('" + StringUtils.join(exchangeMetabolites, "','") + "')";
		try {
			Statement st = conn.createStatement();
			String query = "SELECT tf,primes, atom_bond, concat(`source`,':',`type`,':',`id`,'_x') `ex_id`,`type` FROM mtable." + type + "_" + level + " WHERE id IN " + exchange_met_str + " " + "	AND source = '" + sourceName + "' "
					+ "	GROUP BY id,source HAVING type = MIN(type);";
			ResultSet rs = st.executeQuery(query);

			String[] tf;
			String[] primes;
			String[] atom_bond;
			String ex_id;
			String key;
			PrintWriter ex_sajStream = new PrintWriter(new FileWriter(path + paramFolder + exchangeFolder + "rule_ex.saj"));
			PrintWriter ex_sbjStream = new PrintWriter(new FileWriter(path + paramFolder + exchangeFolder + "rule_ex.sbj"));
			PrintWriter ex_setStream = new PrintWriter(new FileWriter(path + setFolder + exchangeFolder + "rule_ex.index"));
			PrintWriter ex_SimilarityStream = new PrintWriter(new FileWriter(path + paramFolder + exchangeFolder + "rule_ex.similarity"));
			PrintWriter ex_flux_Stream = new PrintWriter(new FileWriter(path + paramFolder + exchangeFolder + "rule.ex"));

			PrintWriter ex_scjStream;
			PrintWriter ex_setCpdStream;

			PrintWriter ex_sajStream_b = null;
			PrintWriter ex_sbjStream_b = null;
			PrintWriter ex_setStream_b = null;
			PrintWriter ex_SimilarityStream_b = null;
			PrintWriter ex_flux_Stream_b = null;

			if (printforBinary) {
				ex_sajStream_b = new PrintWriter(new FileWriter(path + paramFolder + exchangeFolder + "rule_ex_unity.saj"));
				ex_sbjStream_b = new PrintWriter(new FileWriter(path + paramFolder + exchangeFolder + "rule_ex_unity.sbj"));
				ex_setStream_b = new PrintWriter(new FileWriter(path + setFolder + exchangeFolder + "rule_ex_unity.index"));
				ex_SimilarityStream_b = new PrintWriter(new FileWriter(path + paramFolder + exchangeFolder + "rule_ex_unity.similarity"));
				ex_flux_Stream_b = new PrintWriter(new FileWriter(path + paramFolder + exchangeFolder + "rule_unity.ex"));
			}

			while (rs.next()) {
				HashMap<String, Integer> atomFrequencyMap = new HashMap<String, Integer>();
				HashMap<String, Integer> bondFrequencyMap = new HashMap<String, Integer>();
				tf = rs.getString("tf").replaceAll("[\\[\\]]", "").split("[,][ ]");
				primes = rs.getString("primes").replaceAll("[\\[\\]]", "").split("[,][ ]");
				atom_bond = rs.getString("atom_bond").replaceAll("[\\[\\]]", "").split("[,][ ]");

				ex_id = rs.getString("ex_id");
				ex_flux_Stream.println("'" + ex_id + "' " + exchangefluxLoookup.get(ex_id.split("[:]")[2]));

				if (printforBinary) {
					int metabCount = Integer.parseInt(exchangefluxLoookup.get(ex_id.split("[:]")[2]));
					if (metabCount < 0) {
						for (int i = 0; i < Math.abs(metabCount); i++) {
							ex_flux_Stream_b.println("'" + ex_id + "_s" + i + "' 1");
							ex_SimilarityStream_b.println("'" + ex_id + "_s" + i + "' 1");
							ex_setStream_b.println("'" + ex_id + "_s" + i + "'");
						}
					}
					else {
						for (int i = 0; i < metabCount; i++) {
							ex_flux_Stream_b.println("'" + ex_id + "_p" + i + "' 1");
							ex_SimilarityStream_b.println("'" + ex_id + "_p" + i + "' 1");
							ex_setStream_b.println("'" + ex_id + "_p" + i + "'");
						}
					}
				}

				for (int i = 0; i < tf.length; i++) {
					if (tf[i].equalsIgnoreCase("f")) continue;
					if (atom_bond[i].equalsIgnoreCase("a")) {
						key = atom_bond[i] + ":" + primes[i];
						if (!atomFrequencyMap.containsKey(key)) atomFrequencyMap.put(key, 0);
						atomFrequencyMap.put(key, atomFrequencyMap.get(key) + 1);
					}
					else {
						key = atom_bond[i] + ":" + primes[i];
						if (!bondFrequencyMap.containsKey(key)) bondFrequencyMap.put(key, 0);
						bondFrequencyMap.put(key, bondFrequencyMap.get(key) + 1);
					}
				}
				Set<String> keySet = atomFrequencyMap.keySet();

				for (String a : keySet) {
					ex_sajStream.println("'" + a + "'.'" + ex_id + "' " + atomFrequencyMap.get(a));
				}

				if (printforBinary) {
					int metabCount = Integer.parseInt(exchangefluxLoookup.get(ex_id.split("[:]")[2]));
					if (metabCount < 0) {
						for (int i = 0; i < Math.abs(metabCount); i++) {
							for (String a : keySet) {
								ex_sajStream_b.println("'" + a + "'.'" + ex_id + "_s" + i + "' -" + atomFrequencyMap.get(a));
							}
						}
					}
					else {
						for (int i = 0; i < metabCount; i++) {
							for (String a : keySet) {
								ex_sajStream_b.println("'" + a + "'.'" + ex_id + "_p" + i + "' " + atomFrequencyMap.get(a));
							}
						}
					}
				}

				keySet = bondFrequencyMap.keySet();
				for (String b : keySet) {
					ex_sbjStream.println("'" + b + "'.'" + ex_id + "' " + bondFrequencyMap.get(b));
				}

				if (printforBinary) {
					int metabCount = Integer.parseInt(exchangefluxLoookup.get(ex_id.split("[:]")[2]));
					if (metabCount < 0) {
						for (int i = 0; i < Math.abs(metabCount); i++) {
							for (String b : keySet) {
								ex_sbjStream_b.println("'" + b + "'.'" + ex_id + "_s" + i + "' -" + bondFrequencyMap.get(b));
							}
						}
					}
					else {
						for (int i = 0; i < metabCount; i++) {
							for (String b : keySet) {
								ex_sbjStream_b.println("'" + b + "'.'" + ex_id + "_p" + i + "' " + bondFrequencyMap.get(b));
							}
						}
					}
				}

				ex_setStream.println("'" + ex_id + "'");
				ex_SimilarityStream.println("'" + ex_id + "' 1");

				exchangeAtomFrequencyMap.put(ex_id, atomFrequencyMap);
				exchangeBondFrequencyMap.put(ex_id, bondFrequencyMap);
			}
			st.close();

			ex_sajStream.close();
			ex_sbjStream.close();
			ex_setStream.close();
			ex_SimilarityStream.close();
			ex_flux_Stream.close();

			if (printforBinary) {
				ex_sajStream_b.close();
				ex_sbjStream_b.close();
				ex_setStream_b.close();
				ex_SimilarityStream_b.close();
				ex_flux_Stream_b.close();
			}
			if (exchangeOfnwIO) {
				ex_scjStream = new PrintWriter(new FileWriter(path + paramFolder + exchangeFolder + "rule_ex.scj"));
				ex_scjSetStream = new PrintWriter(new FileWriter(path + setFolder + exchangeFolder + "rule_ex_scj.index"));
				
				ex_setCpdStream = new PrintWriter(new FileWriter(path + setFolder + exchangeFolder + "rule_cpd_ex.index"));

				if (printforBinary) {
					ex_scjStream_rev = new PrintWriter(new FileWriter(path + paramFolder + exchangeFolder + "rule_ex_rev.scj"));
					ex_scjSetStream_rev = new PrintWriter(new FileWriter(path + setFolder + exchangeFolder + "rule_ex_rev_scj.index"));
					
					ex_setCpdStream_rev = new PrintWriter(new FileWriter(path + setFolder + exchangeFolder + "rule_cpd_ex_rev.index"));
					ex_pairsStream = new PrintWriter(new FileWriter(path + paramFolder + exchangeFolder + "rule_cpd_ex.pairs"));
				}

				for (String cpd : cpdSet) {
					ex_scjStream.println("'" + cpd + "'.'" + cpd + "_x' 1");
					ex_setCpdStream.println("'" + cpd + "_x'");
				}

				if (printforBinary) {
					for (String cpd : cpdSet) {
						ex_scjStream_rev.println("'" + cpd + "'.'" + cpd + "_r_x' -1");
						ex_setCpdStream_rev.println("'" + cpd + "_r_x'");
						ex_pairsStream.println("'" + cpd + "_x'.'" + cpd + "_r_x' 1");
					}
				}
				ex_scjStream.close();
				ex_scjSetStream.close();
				ex_setCpdStream.close();

				if (printforBinary) {
					ex_scjStream_rev.close();
					ex_scjSetStream_rev.close();
					ex_setCpdStream_rev.close();
					ex_pairsStream.close();
				}
			}

		}
		catch (SQLException ex) {
			System.err.println(ex.getMessage());
		}
		/* only these reactions need to be printed, all the cpds from the reactions will be printed as full exchange as 1 */
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private static void changeDbSessionParameters() {
		try {
			Statement st = conn.createStatement();

			st.executeUpdate("SET SESSION GROUP_CONCAT_MAX_LEN = 999999999;");
			st.executeUpdate("SET GLOBAL max_allowed_packet=1073741824;");
			st.executeUpdate("SET SESSION sql_big_selects = 1;");
			st.executeUpdate("SET SESSION sql_buffer_result = 1;");
			st.executeUpdate("SET SESSION max_heap_table_size = 3221225472;");
			st.executeUpdate("SET SESSION tmp_table_size = 4294967294;");
			st.executeUpdate("SET SESSION join_buffer_size = 3221225472;");
			st.executeUpdate("SET SESSION sort_buffer_size = 4294967294;");
			st.executeUpdate("SET SESSION myisam_sort_buffer_size = 4294967294;");
			st.executeUpdate("SET GLOBAL myisam_max_sort_file_size = 3221225472;");
			st.executeUpdate("SET SESSION myisam_repair_threads = 7;");

			st.close();
		}
		catch (SQLException ex) {
			System.err.println(ex.getMessage());
		}
	}

	private static void metaboliteAndRules() {

		try {
			final Semaphore s = new Semaphore(THREADCOUNT);
			Statement st = conn.createStatement();
			String query = "SELECT A.`hash`,B.`moiety_hash`, " + "	GROUP_CONCAT(DISTINCT CONCAT(A.`source`,':',A.`type`,':',A.`id`) ORDER BY A.`source`,A.`type`,A.`id` SEPARATOR ',') `ideez`, " + " 	B.`TF`, B.`primes`, B.`atom_bond` "
					+ "	from metrxn0200.gi_" + type + " A, `mtable`.`" + type + "_" + level + "` B " + "	WHERE " + "	A.`source` = B.`source` " + "	AND A.`hash` = B.`hash` " + "	AND A.id =  B.id " + "	AND A.source " + source_qualifier + " " + sources
					+ " " + "	AND B.source " + source_qualifier + " " + sources + " " + "	GROUP BY A.`hash`,B.moiety_hash;";
			ResultSet rs = st.executeQuery(query);
			idLookup = Collections.synchronizedSortedMap(new TreeMap<String, String>());
			moietyhash_id = Collections.synchronizedSortedMap(new TreeMap<String, SortedSet<String>>());
			metabhash_id = Collections.synchronizedSortedMap(new TreeMap<String, SortedSet<String>>());
			moietyHash_metabHash = Collections.synchronizedSortedMap(new TreeMap<String, String>());

			cpdAtomFrequencyMap = Collections.synchronizedMap(new HashMap<String, HashMap<String, Integer>>());
			cpdBondFrequencyMap = Collections.synchronizedMap(new HashMap<String, HashMap<String, Integer>>());

			atomSet = Collections.synchronizedSortedSet(new TreeSet<String>());
			bondSet = Collections.synchronizedSortedSet(new TreeSet<String>());

			metabIdSet = Collections.synchronizedSortedSet(new TreeSet<String>());

			metabIdMappingStream = new PrintWriter(new BufferedWriter(new FileWriter(path + mapFolder + "metabolite.map")));
			moietyMappingStream = new PrintWriter(new BufferedWriter(new FileWriter(path + mapFolder + "moiety_hash.map")));

			atomRulesStream = new PrintWriter(new BufferedWriter(new FileWriter(path + paramFolder + "atoms.formula")));
			bondRulesStream = new PrintWriter(new BufferedWriter(new FileWriter(path + paramFolder + "bonds.formula")));
			
			atomRulesSetStream = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "atoms_formula.index")));
			bondRulesSetStream = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "bonds_formula.index")));

			metabStream = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "metabolite.index")));
			atomStream = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "atoms.index")));
			bondStream = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "bonds.index")));
			
			metabIdPriceStream = new PrintWriter(new BufferedWriter(new FileWriter(path + paramFolder + "metabolite.price")));
			metabIdCarbonStream = new PrintWriter(new BufferedWriter(new FileWriter(path + paramFolder + "metabolite.carbon")));

			while (rs.next()) {
				final String metab_hash = rs.getString("hash");
				final String moiety_hash = rs.getString("moiety_hash");
				String ideezStr = rs.getString("ideez").replace("brenda:", "br:").replace("chebi:", "ci:").replace("chembl:", "cl:").replace("ecmdb:", "eb:").replace("ecocyc:", "ey:").replace("hmdb:", "hm:").replace("kegg:", "kg:")
						.replace("metacyc:", "my:");
				final String[] ideez = ideezStr.split("[,]");
				final String[] primes = rs.getString("primes").replaceAll("[\\[\\]]", "").split("[,][ ]", -1);
				final String[] tf = rs.getString("tf").replaceAll("[\\[\\]]", "").split("[,][ ]", -1);
				final String[] atom_bond = rs.getString("atom_bond").replaceAll("[\\[\\]]", "").split("[,][ ]", -1);

				s.acquireUninterruptibly();

				new Thread() {
					public void run() {
						try {

							if (!moietyhash_id.containsKey(moiety_hash)) {
								SortedSet<String> midezz = Collections.synchronizedSortedSet(new TreeSet<String>());
								moietyhash_id.put(moiety_hash, midezz);
							}
							moietyhash_id.get(moiety_hash).addAll(Arrays.asList(ideez));

							if (!metabhash_id.containsKey(metab_hash)) {
								SortedSet<String> midezz = Collections.synchronizedSortedSet(new TreeSet<String>());
								metabhash_id.put(metab_hash, midezz);
							}
							metabhash_id.get(metab_hash).addAll(Arrays.asList(ideez));
							String midFirst = metabhash_id.get(metab_hash).first();

							moietyHash_metabHash.put(metab_hash, moiety_hash);

							HashMap<String, Integer> atomfrequencyMap = new HashMap<String, Integer>();
							HashMap<String, Integer> bondfrequencyMap = new HashMap<String, Integer>();
							String m_key;
							// String m_m_key;
							for (int i = 0; i < tf.length; i++) {
								if (tf[i].equalsIgnoreCase("f")) continue;
								if (atom_bond[i].equalsIgnoreCase("a")) {
									m_key = atom_bond[i] + ":" + primes[i];
									// m_m_key = "'" + ideez[0] + "'." + m_key;
									if (!atomfrequencyMap.containsKey(m_key)) {
										atomfrequencyMap.put(m_key, 0);
									}
									atomfrequencyMap.put(m_key, atomfrequencyMap.get(m_key) + 1);
									atomSet.add(m_key);
								}
								else {
									m_key = atom_bond[i] + ":" + primes[i];
									// m_m_key = "'" + ideez[0] + "'." + m_key;
									if (!bondfrequencyMap.containsKey(m_key)) {
										bondfrequencyMap.put(m_key, 0);
									}
									bondfrequencyMap.put(m_key, bondfrequencyMap.get(m_key) + 1);
									bondSet.add(m_key);
								}
							}

							metabIdMappingStream.println(midFirst + "\t" + Arrays.deepToString(ideez));
							Set<String> keySet = atomfrequencyMap.keySet();
							for (String key : keySet) {
								atomRulesStream.println("'" + key + "'.'" + midFirst + "' " + atomfrequencyMap.get(key));
								atomRulesSetStream.println("'" + key + "'.'" + midFirst + "'");
							}

							keySet = bondfrequencyMap.keySet();
							for (String key : keySet) {
								bondRulesStream.println("'" + key + "'.'" + midFirst + "' " + bondfrequencyMap.get(key));
								bondRulesSetStream.println("'" + key + "'.'" + midFirst + "'");
							}

							// metabStream.println("'" + ideez[0] + "'");

							metabIdSet.add(midFirst);

							cpdAtomFrequencyMap.put(moiety_hash, atomfrequencyMap);
							cpdBondFrequencyMap.put(moiety_hash, bondfrequencyMap);
						}
						finally {
							s.release();
						}
					}
				}.start();
			}

			/* pause for all threads to close */
			while (s.availablePermits() < THREADCOUNT) {
			}

			for (String mid : metabIdSet) {
				metabStream.println("'" + mid + "'");
			}

			for (String a : atomSet) {
				atomStream.println("'" + a + "'");
			}

			for (String b : bondSet) {
				bondStream.println("'" + b + "'");
			}

			Set<String> keySet = moietyhash_id.keySet();
			for (String str : keySet) {
				moietyMappingStream.println(moietyhash_id.get(str).first() + "\t" + moietyhash_id.get(str).toString());
			}

			st.close();

			metabIdMappingStream.close();
			moietyMappingStream.close();

			atomRulesStream.close();
			bondRulesStream.close();
			
			atomRulesSetStream.close();
			bondRulesSetStream.close();

			metabStream.close();
			atomStream.close();
			bondStream.close();
			
			/*print price*/
			/*based on metabHash*/
			st = conn.createStatement();
			query = "select distinct A.hash,B.price_per_g `price`,A.atomBondNumber from metrxn0200.gi_" + type + " A, "
					+ "emolecules.price_source B where A.id = B.id and A.source = B.source AND A.source " + source_qualifier + " " + sources;
			
			rs = st.executeQuery(query);
			int carbonCount;
			HashSet<String> ifPrinted = new HashSet<String>();
			String metabHash;
			while (rs.next()) {
				metabHash = metabhash_id.get(rs.getString("hash")).first();
				if (ifPrinted.contains(metabHash)) continue;
				ifPrinted.add(metabHash);
				metabIdPriceStream.println("'"+metabHash+"'"+"\t"+rs.getString("price"));
				List<String> asList = Arrays.asList(rs.getString("atomBondNumber").replaceAll("[\\[\\]]", "").split("[,][ ]"));
				carbonCount = Collections.frequency(asList,"6");
				metabIdCarbonStream.println("'"+metabHash+"'"+"\t"+carbonCount);
			}
			metabIdPriceStream.close();
			metabIdCarbonStream.close();
			

		}
		catch (SQLException ex) {
			System.err.println(ex.getMessage());
		}
		catch (IOException e1) {
			e1.printStackTrace();
		}

	}

	private static void reactionsAndRules() {
		// equationMd5 = new HashMap<String, String>();
		rulesRidCnt = new HashMap<String, Integer>();
		try {
			Statement st = conn.createStatement();
			String query = "SELECT DISTINCT rule_" + type + "_" + loop + "_md5 `rule_md5`,rid,source " + "FROM `megatable`.`megatable_myisam` WHERE balanced = 'true';";
			ResultSet rs = st.executeQuery(query);

			String rule_md5;
			while (rs.next()) {
				rule_md5 = rs.getString("rule_md5").split("[!]")[0];
				if (!rulesRidCnt.containsKey(rule_md5)) rulesRidCnt.put(rule_md5, 0);
				rulesRidCnt.put(rule_md5, rulesRidCnt.get(rule_md5) + 1);
			}
			st.close();
		}
		catch (SQLException ex) {
			System.err.println(ex.getMessage());
		}

		rulesMap = new HashMap<String, String>();

		try {
			Statement st = conn.createStatement();
			String query = "SELECT DISTINCT rule_" + type + "_" + loop + "_md5 `rule_md5`,rule_" + type + "_" + loop + " `rule` " + "FROM `metrxn0200`.`reactionrules` WHERE balanced = 'true';";
			ResultSet rs = st.executeQuery(query);

			String rule_md5;
			String rule;
			while (rs.next()) {
				rule_md5 = rs.getString("rule_md5");
				rule = rs.getString("rule");
				rulesMap.put(rule_md5, rule.split("[!]", -1)[0]);
			}
			st.close();
		}
		catch (SQLException ex) {
			System.err.println(ex.getMessage());
		}

		try {
			ridPrinted = Collections.synchronizedMap(new HashMap<String, Integer>());
			cpdRid = Collections.synchronizedSortedSet(new TreeSet<String>());
			
			SCJStream = new PrintWriter(new BufferedWriter(new FileWriter(path + paramFolder + "rule.scj")));
			SCJSetStream = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "rule_scj.index")));
			
			SAJStream = new PrintWriter(new BufferedWriter(new FileWriter(path + paramFolder + "rule.saj")));
			SAJSetStream = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "rule_saj.index")));
			
			SBJStream = new PrintWriter(new BufferedWriter(new FileWriter(path + paramFolder + "rule.sbj")));
			SBJSetStream = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "rule_sbj.index")));
			
			cpdSetStream = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "cpd.index")));
			ruleIdStream = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "rule.index")));
			
			ruleIdEcStream = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "rule.ec")));
			
			
			
			if (printforBinary) {
				SAJStream_rev = new PrintWriter(new BufferedWriter(new FileWriter(path + paramFolder + "rule_rev.saj")));
				SAJSetStream_rev = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "rule_rev_saj.index")));
				
				SCJStream_rev = new PrintWriter(new BufferedWriter(new FileWriter(path + paramFolder + "rule_rev.scj")));
				SCJSetStream_rev = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "rule_rev_scj.index")));
				
				SBJStream_rev = new PrintWriter(new BufferedWriter(new FileWriter(path + paramFolder + "rule_rev.sbj")));
				SBJSetStream_rev = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "rule_rev_sbj.index")));
				
				ruleIdStream_rev = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "rule_rev.index")));
				
				ruleIdEcStream_rev = new PrintWriter(new BufferedWriter(new FileWriter(path + setFolder + "rule_rev.ec")));

				rulePairsStream = new PrintWriter(new BufferedWriter(new FileWriter(path + paramFolder + "rules.pairs")));
			}
			Statement st = conn.createStatement();
			String query = "SELECT A." + type + "_eq_md5 `str_eq_md5`, A.native_Equation_md5 `eq_md5`, A.rule_" + type + "_" + loop + "_md5 `rule_md5`,"
					+ "	GROUP_CONCAT(DISTINCT concat(B.id,':',B.type,':',B.source,':',B.hash)) `mid_type_src_hash`,A.equation,concat(A.source,'_',A.rid) `rid`, coalesce(C.ec,'unknown') `ec`" 
					+ "   FROM megatable.megatable_myisam A left join paper3.ecrxn C on A.rid = C.rxnAbbreviation AND A.source = C.source inner join metrxn0200.gi_" + type + " B " + " on "
					+ "	A.source = B.source " + "	AND A.mid = B.id " + " AND balanced = 'true'	AND A." + type + " = B.hash " 
					+ "   AND (A.stereo <> '35933406de6e4eb07b711f8c5630dd3c' " + " OR A.`order` <> '35933406de6e4eb07b711f8c5630dd3c' "
					+ " 	OR A.`connection` <> '35933406de6e4eb07b711f8c5630dd3c')	" + " AND A.source " + source_qualifier + " " + sources 
					+ "	AND A.native_Equation_md5 IN (SELECT MIN(C.native_Equation_md5) `min_native_equation` "
					+ "	FROM megatable.megatable_myisam C WHERE C.source " + source_qualifier + " " + sources + " AND balanced = 'true' GROUP BY C." + type + "_eq_md5) " + "	GROUP BY A." + type + "_eq_md5, A.native_Equation_md5,A.rule_" + type + "_"
					+ 	loop + "_md5 HAVING COUNT(DISTINCT mid) = COUNT(DISTINCT stereo) " + " OR COUNT(DISTINCT mid) = COUNT(DISTINCT `order`);";
			System.out.println(query);
			ResultSet rs = st.executeQuery(query);

			final Semaphore s = new Semaphore(THREADCOUNT);
			cpdSet = Collections.synchronizedSortedSet(new TreeSet<String>());
			cpdNIdHash = Collections.synchronizedSortedMap(new TreeMap<String, String>());		
			
			// rulePrinted = Collections.synchronizedSortedSet(new TreeSet<String>());
			while (rs.next()) {
				final String str_eq_md5 = rs.getString("str_eq_md5");
				final String eq_md5 = rs.getString("eq_md5");
				final String rule_md5 = rs.getString("rule_md5");
				final String[] mid_type_src_hash = rs.getString("mid_type_src_hash").split("[,]");
				final String equation = rs.getString("equation");
				final String rid_str = rs.getString("rid");
				final String ec = rs.getString("ec");
				// if (rulesRidCnt.get(rule_md5.split("[!]")[0]) == 1) continue;
				

				if (!ridPrinted.containsKey(rid_str.toUpperCase())) {
					ridPrinted.put(rid_str.toUpperCase(), 0);
				} else {
					ridPrinted.put(rid_str.toUpperCase(), ridPrinted.get(rid_str.toUpperCase()) + 1);
				}
				
				
				final String rid = rid_str.toUpperCase()+"#"+ridPrinted.get(rid_str.toUpperCase());
				
				s.acquireUninterruptibly();
				new Thread() {
					public void run() {

						Integer stoich_d;
						HashMap<String, String> mid_N_mid = new HashMap<String, String>();
						String normalized_mid;
						String[] rule_a_b = rulesMap.get(rule_md5.split("[!]")[0]).replaceAll("[\\[\\]]", "").replace("\n", "").split("[,][ ]", -1);

						for (int i = 0; i < mid_type_src_hash.length; i++) {
							if (printByMoietyHash) {
								normalized_mid = moietyhash_id.get(moietyHash_metabHash.get(mid_type_src_hash[i].split("[:]", -1)[3])).first();
							}
							else {
								normalized_mid = metabhash_id.get(mid_type_src_hash[i].split("[:]", -1)[3]).first();
							}

							if (normalized_mid.equalsIgnoreCase("null")) {
								System.err.println(mid_type_src_hash[i]);
							}

							cpdSet.add(normalized_mid);
							mid_N_mid.put(mid_type_src_hash[i].split("[:]", -1)[0].toUpperCase(), normalized_mid);

							cpdNIdHash.put(normalized_mid, moietyHash_metabHash.get(mid_type_src_hash[i].split("[:]", -1)[3]));

						}
						if (rule_md5.split("[!]")[1].equalsIgnoreCase("0"))
							stoich_d = 1;
						else
							stoich_d = -1;

						try {
							if (!equation.split("[ ][<][=][>][ ]")[0].equalsIgnoreCase(equation.split("[ ][<][=][>][ ]")[1])) {

								boolean containsNoUnityStoich = false;

								String cpd_r;
								int s_r;

								TreeMap<String, Integer> frequencyMap = new TreeMap<String, Integer>();
								TreeMap<String, Integer> scjToPrint = new TreeMap<String, Integer>();
								String[] split_r = equation.split("[ ][<][=][>][ ]")[0].split("[ ][+][ ]");
								for (int i = 0; i < split_r.length; i++) {

									cpd_r = split_r[i].trim();
									s_r = 1;

									if (cpd_r.matches("[0-9]+?[ ].+")) {
										s_r = Integer.valueOf(cpd_r.replace(cpd_r.replaceFirst("[0-9]+?[ ]", ""), "").trim());
									}
									cpd_r = cpd_r.replaceFirst("[0-9]+?[ ]", "");

									if (!frequencyMap.containsKey(cpd_r)) frequencyMap.put(cpd_r, 0);

									frequencyMap.put(cpd_r, frequencyMap.get(cpd_r) - (stoich_d * s_r));

								}
								split_r = equation.split("[ ][<][=][>][ ]")[1].split("[ ][+][ ]");
								for (int i = 0; i < split_r.length; i++) {

									cpd_r = split_r[i].trim();
									s_r = 1;

									if (cpd_r.matches("[0-9]+?[ ].+")) {
										s_r = Integer.valueOf(cpd_r.replace(cpd_r.replaceFirst("[0-9]+?[ ]", ""), "").trim());
									}
									cpd_r = cpd_r.replaceFirst("[0-9]+?[ ]", "");
									if (!frequencyMap.containsKey(cpd_r)) frequencyMap.put(cpd_r, 0);

									frequencyMap.put(cpd_r, frequencyMap.get(cpd_r) + (stoich_d * s_r));

								}

								Set<String> keySet = frequencyMap.keySet();
								String key2;
								for (String key : keySet) {
									if (frequencyMap.get(key) != 0) {
										if (null == mid_N_mid.get(key.toUpperCase())) {
											System.out.println(key + " " + rid);
											continue;
										}

										if (mid_N_mid.get(key.toUpperCase()).equalsIgnoreCase("brenda:mol:12")) continue;

										key2 = "'" + mid_N_mid.get(key.toUpperCase()) + "'.'" + rid + "'";

										if (!scjToPrint.containsKey(key2)) {
											scjToPrint.put(key2, 0);
										}
										scjToPrint.put(key2, scjToPrint.get(key2) + frequencyMap.get(key));
									}
								}
								Set<String> keySet2 = scjToPrint.keySet();
								boolean printRid = true;

								for (String key : keySet2) {
									if (Math.abs(scjToPrint.get(key)) > 1) {
										printRid = false; /* so that rules from equations with non unity stoich is not printed */
									}
								}

								/*
								 * if (rulePrinted.contains(rule_md5.split("[!]")[0])) { printRid = false; } else {
								 * rulePrinted.add(rule_md5.split("[!]")[0]); }
								 */
								if (printRid) {

									for (String key : keySet2) {
										cpdRid.add(key);
										if (Math.abs(scjToPrint.get(key)) > 0) {
											SCJStream.println(key + " " + scjToPrint.get(key));
											SCJSetStream.println(key);
											
											SCJStream_rev.println(StringUtils.removeEnd(key, "'") + "_r' " + scjToPrint.get(key) * -1);
											SCJSetStream_rev.println(StringUtils.removeEnd(key, "'") + "_r'");
										}
									}
									ruleIdStream.println("'" + rid + "'");
									
									ruleIdEcStream.println("'" + rid + "'.'"+ec+"'");
									
									if (printforBinary) {
										ruleIdStream_rev.println("'" + StringUtils.removeEnd(rid, "'") + "_r'");
										ruleIdEcStream_rev.println("'" + StringUtils.removeEnd(rid, "'") + "_r'.'"+ec+"'");
										rulePairsStream.println("'" + rid + "'.'" + StringUtils.removeEnd(rid, "'") + "_r' 1");
									}
								}

								/****** sij done ******/

								if (printRid) {
									/******* saj and sbj start ***********/
									TreeMap<String, Integer> atomfrequencyMap = new TreeMap<String, Integer>();
									TreeMap<String, Integer> bondfrequencyMap = new TreeMap<String, Integer>();
									String key;
									for (String m : rule_a_b) {
										if (m.trim().equalsIgnoreCase("")) continue;
										key = m.replace("-", "").trim();

										if (key.startsWith("a:")) {
											key = "'" + key + "'.'" + rid + "'";
											if (!atomfrequencyMap.containsKey(key)) atomfrequencyMap.put(key, 0);

											if (m.startsWith("-")) {
												atomfrequencyMap.put(key, atomfrequencyMap.get(key) - 1);
											}
											else {
												atomfrequencyMap.put(key, atomfrequencyMap.get(key) + 1);
											}
										}
										else {
											key = "'" + key + "'.'" + rid + "'";
											if (!bondfrequencyMap.containsKey(key)) bondfrequencyMap.put(key, 0);

											if (m.startsWith("-")) {
												bondfrequencyMap.put(key, bondfrequencyMap.get(key) - 1);
											}
											else {
												bondfrequencyMap.put(key, bondfrequencyMap.get(key) + 1);
											}
										}
									}
									keySet2 = atomfrequencyMap.keySet();
									for (String key3 : keySet2) {
										SAJStream.println(key3 + " " + atomfrequencyMap.get(key3));										
										SAJStream_rev.println(StringUtils.removeEnd(key3, "'") + "_r' " + atomfrequencyMap.get(key3) * -1);
										SAJSetStream.println(key3);
										SAJSetStream_rev.println(StringUtils.removeEnd(key3, "'") + "_r'");
									}
									keySet2 = bondfrequencyMap.keySet();
									for (String key3 : keySet2) {
										SBJStream.println(key3 + " " + bondfrequencyMap.get(key3));
										SBJStream_rev.println(StringUtils.removeEnd(key3, "'") + "_r' " + bondfrequencyMap.get(key3) * -1);
										SBJSetStream.println(key3);
										SBJSetStream_rev.println(StringUtils.removeEnd(key3, "'") + "_r'");
									}

								}

								/******* saj and sbj end ***********/

							}
						}
						finally {
							s.release();
						}
					}
				}.start();
			}
			/* pause for all threads to close */
			while (s.availablePermits() < THREADCOUNT) {
			}

			for (String cpd : cpdSet) {
				cpdSetStream.println("'" + cpd + "'");
			}

			st.close();
			SCJStream.close();
			SAJStream.close();
			SBJStream.close();
			SCJSetStream.close();
			SAJSetStream.close();
			SBJSetStream.close();
			cpdSetStream.close();
			ruleIdStream.close();
			ruleIdEcStream.close();
			if (printforBinary) {
				SCJStream_rev.close();
				SAJStream_rev.close();
				SBJStream_rev.close();
				SCJSetStream_rev.close();
				SAJSetStream_rev.close();
				SBJSetStream_rev.close();
				ruleIdStream_rev.close();
				rulePairsStream.close();
				ruleIdEcStream_rev.close();
			}

		}
		catch (SQLException ex) {
			System.err.println(ex.getMessage());
		}
		catch (IOException e) {
			e.printStackTrace();
		}

	}

}
