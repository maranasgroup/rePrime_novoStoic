/*TODO : convert most of this to be called from the database on MONDAY 8/10/2015*/

package paper6;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.Semaphore;

import com.google.common.collect.TreeMultimap;
import com.opencsv.CSVReader;

public class i_gamsDataGenerator {

	private static final int THREADCOUNT = 12;
	private static String source = "kegg";
	private static String loop = "3";
	private static String type = "stereo";
	private static TreeMultimap<String, String> ruleHashRid;
	private static TreeMultimap<String, String> ruleHash_Mid;
	private static TreeMultimap<String, String> moietyHashMid;

	private static HashMap<String, String> mid_moietyHash;
	private static HashMap<String, String> rid_ruleHash;
	private static TreeMap<String, TreeMap<String, Integer>> reactionMoiety;
	private static TreeMap<String, TreeMap<String, Integer>> metaboliteMoiety;
	private static HashSet<String> metaboliteIdSet;
	private static boolean ignoreBonds = true;
	private static HashSet<String> moietyIndex;
	private static HashSet<String> moietyIndex_extended;
	private static List<String> exchangeMetabolites;
	private static List<String> exchangeFlux;
	private static Set<String> reactionFilter;
	private static boolean toFilterReactions = false;
	private static HashMap<String, String> exchangeMidMapping;
	private static boolean printReverse = true;

	private static TreeMap<String, TreeMap<String, Float>> ruleHashMetaboliteMoiety;
	private static SortedMap<String, TreeMap<String, Float>> similarity;
	private static int largestk = 0;
	private static int countUsable = 0;
	private static HashMap<String, HashSet<String>> equations;
	private static HashMap<String, String> rid_direction;
	private static TreeSet<String> compoundNotFound;
	private static boolean printSimilarity = false;

	public static void main(String[] args) {

		populateTables();
		printTables();

	}

	private static void populateTables() {

		int offset = 0;

		int level = Integer.valueOf(loop);
		if (type.equalsIgnoreCase("connection")) {
			offset = 5;
		}
		else if (type.equalsIgnoreCase("order")) {
			offset = 9;
		}
		else if (type.equalsIgnoreCase("stereo")) {
			offset = 13;
		}
		offset = offset + level;
		String[] nextLine;
		ruleHashRid = TreeMultimap.create();
		moietyHashMid = TreeMultimap.create();
		reactionMoiety = new TreeMap<String, TreeMap<String, Integer>>();
		metaboliteMoiety = new TreeMap<String, TreeMap<String, Integer>>();
		exchangeMidMapping = new HashMap<String, String>();
		ruleHash_Mid = TreeMultimap.create();

		mid_moietyHash = new HashMap<String, String>();
		rid_ruleHash = new HashMap<String, String>();
		rid_direction = new HashMap<String, String>();

		metaboliteIdSet = new HashSet<String>();
		compoundNotFound = new TreeSet<String>();
		moietyIndex = new HashSet<String>();
		moietyIndex_extended = new HashSet<String>();
		reactionFilter = new TreeSet<String>();
		reactionFilter.addAll(Arrays.asList(new String[] { "R00230", "R00343", "R00354", "R00410", "R00472", "R00479", "R00605", "R00761", "R01056", "R01323", "R01324", "R01361", "R01529", "R01641", "R01644", "R05336", "R05337", "R05338", "R05339",
				"R05576", "R08575", "R09098", "R09280", "R10000", "R09317", "R05595", "R03031", "R01195" }));
		try {
			exchangeMetabolites = Arrays.asList(new String[] { "C16485", "C00014", "C11014", "C00001" });
			exchangeFlux = Arrays.asList(new String[] { "-10", "-10", "10", "10" });

			CSVReader reader = new CSVReader(new FileReader("D:/metrxn0200/input/" + source + ".mbt"), '\t', '"', 1);
			while ((nextLine = reader.readNext()) != null) {
				if (toFilterReactions) {
					if (!reactionFilter.contains(nextLine[4])) continue;
				}
				reactionFilter.remove(nextLine[4]);
				ruleHashRid.put(nextLine[offset].split("[!]")[0], nextLine[4]);
				rid_ruleHash.put(nextLine[4], nextLine[offset].split("[!]")[0]);
				rid_direction.put(nextLine[4], nextLine[offset].split("[!]")[1]);
				ruleHash_Mid.put(nextLine[offset].split("[!]")[0], nextLine[0].toUpperCase());
				metaboliteIdSet.add(nextLine[0].toUpperCase());
			}
			reader.close();
			/*******/
			try {
				HashSet<String> equationMetaboliteAsMoiety;
				reader = new CSVReader(new FileReader("D:/metrxn0200/input/reactions/balanced/" + source + ".tab"), '\t', '"', 1);
				equations = new HashMap<String, HashSet<String>>();
				Integer stoich;
				Integer stoich_d;
				while ((nextLine = reader.readNext()) != null) {
					if (!rid_direction.containsKey(nextLine[0])) continue;

					if (nextLine[2].split("[ ][<][=][>][ ]")[0].equalsIgnoreCase(nextLine[2].split("[ ][<][=][>][ ]")[1])) continue;

					equationMetaboliteAsMoiety = new HashSet<String>();

					if (rid_direction.get(nextLine[0]).equalsIgnoreCase("0")) {
						stoich_d = 1;
					}
					else
						stoich_d = -1;

					String[] split_l = nextLine[2].split("[ ][<][=][>][ ]")[0].split("[ ][+][ ]");
					String[] split_r = nextLine[2].split("[ ][<][=][>][ ]")[1].split("[ ][+][ ]");
					String cpd_l;
					String cpd_r;
					Integer s_l;
					Integer s_r;
					/* also fix when cpds are mentioned multiple times on one side */
					for (int i = 0; i < split_r.length; i++) {
						cpd_r = split_r[i].trim();
						s_r = 1;

						if (cpd_r.matches("[0-9]+?[ ].+")) {
							s_r = Integer.valueOf(cpd_r.replace(cpd_r.replaceFirst("[0-9]+?[ ]", ""), "").trim());
						}
						cpd_r = cpd_r.replaceFirst("[0-9]+?[ ]", "");
						for (int j = 0; j < split_l.length; j++) {
							cpd_l = split_l[j].trim();

							s_l = 1;

							if (cpd_l.matches("[0-9]+?[ ].+")) {
								s_l = Integer.valueOf(cpd_l.replace(cpd_l.replaceFirst("[0-9]+?[ ]", ""), "").trim());
							}
							cpd_l = cpd_l.replaceFirst("[0-9]+?[ ]", "");
							if (cpd_r.equalsIgnoreCase(cpd_l)) {
								if (s_l == s_r) {
									split_r[i] = "*";
									split_l[j] = "*";
								}
								else if (s_l > s_r) {
									split_r[i] = "*";
									split_l[j] = (s_l - s_r) + " " + cpd_l;
								}
								else {
									split_r[i] = (s_r - s_l) + " " + cpd_r;
									split_l[j] = "*";
								}
							}
						}
					}

					for (String cpd : split_l) {
						if (cpd.equalsIgnoreCase("*")) continue;

						if (cpd.matches("[0-9]+?[ ].+")) {
							stoich = -1 * Integer.valueOf(cpd.replaceAll("[ ].+", "").trim()) * stoich_d;
						}
						else
							stoich = -1 * stoich_d;

						cpd = cpd.replaceFirst("[0-9]+?[ ]", "");
						equationMetaboliteAsMoiety.add("'cpd:" + cpd + "'.'" + nextLine[0] + "' " + stoich);
						moietyIndex_extended.add("cpd:" + cpd);
					}

					for (String cpd : split_r) {
						if (cpd.equalsIgnoreCase("*")) continue;

						if (cpd.matches("[0-9]+?[ ].+")) {
							stoich = Integer.valueOf(cpd.replaceAll("[ ].+", "").trim()) * stoich_d;
						}
						else
							stoich = 1 * stoich_d;

						cpd = cpd.replaceFirst("[0-9]+?[ ]", "");

						equationMetaboliteAsMoiety.add("'cpd:" + cpd + "'.'" + nextLine[0] + "' " + stoich);
						moietyIndex_extended.add("cpd:" + cpd);
					}
					equations.put(nextLine[0], equationMetaboliteAsMoiety);
				}
				reader.close();
			}
			catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			/*******/
			if (toFilterReactions) {

				if (reactionFilter.size() > 0) {
					try {
						reader = new CSVReader(new FileReader("D:/metrxn0200/input/reactions/equation/KEGG.tab"), '\t', '"', 1);
						HashMap<String, String> equations = new HashMap<String, String>();
						while ((nextLine = reader.readNext()) != null) {
							if (reactionFilter.contains(nextLine[0])) {
								equations.put(nextLine[0], nextLine[1]);
								ruleHashRid.put(nextLine[0], nextLine[0]);
								String[] split = nextLine[1].split("[ ][<][=][>][ ]")[0].split("[ ][+][ ]");
								for (String cpd : split) {
									metaboliteIdSet.add(cpd.replaceFirst("[0-9]+?[ ]", ""));
								}
							}
						}
						reader.close();
					}
					catch (FileNotFoundException e) {
						e.printStackTrace();
					}
				}
			}
			/**/

			String[] tf;
			String[] primes;
			String[] atom_bond;
			TreeSet<String> sortedPrimes;
			reader = new CSVReader(new FileReader("D:/metrxn0200/input/metabolites/structure/" + loop + "." + type + ".gi.moity"), '\t', '"', 1);
			String key;

			while ((nextLine = reader.readNext()) != null) {
				if (!nextLine[1].equalsIgnoreCase(source)) continue;

				TreeMap<String, Integer> frequencyMap = new TreeMap<String, Integer>();
				tf = nextLine[3].replaceAll("[\\[\\]]", "").split("[,][ ]");
				primes = nextLine[4].replaceAll("[\\[\\]]", "").split("[,][ ]");
				atom_bond = nextLine[8].replaceAll("[\\[\\]]", "").split("[,][ ]");
				sortedPrimes = new TreeSet<String>();
				sortedPrimes.addAll(Arrays.asList(primes));
				if (!metaboliteIdSet.contains(nextLine[0])) {
					compoundNotFound.add(sortedPrimes.size() + "\t" + sortedPrimes.toString() + "\t" + nextLine[0]);
				}

				for (int i = 0; i < atom_bond.length; i++) {
					if (tf[i].equalsIgnoreCase("f")) continue;
					if (ignoreBonds) {
						if (atom_bond[i].equalsIgnoreCase("b")) continue;
					}

					key = atom_bond[i] + ":" + primes[i];

					if (!frequencyMap.containsKey(key)) frequencyMap.put(key, 0);

					frequencyMap.put(key, frequencyMap.get(key) + 1);

					moietyIndex.add(key);
					moietyIndex_extended.add(key);
				}

				metaboliteMoiety.put(nextLine[7], frequencyMap);
				// System.out.println("size of freq map making = "+frequencyMap.size());
				moietyHashMid.put(nextLine[7], nextLine[0]);

				mid_moietyHash.put(nextLine[0], nextLine[7]);

				if (exchangeMetabolites.contains(nextLine[0])) {
					exchangeMidMapping.put(nextLine[0], nextLine[7]);
					ruleHashRid.put(nextLine[7], "RXN_" + nextLine[0]);
					reactionMoiety.put(nextLine[7], frequencyMap);
				}
			}
			reader.close();
			offset = 0;
			if (type.equalsIgnoreCase("connection")) {
				offset = 0;
			}
			else if (type.equalsIgnoreCase("order")) {
				offset = 4;
			}
			else if (type.equalsIgnoreCase("stereo")) {
				offset = 8;
			}
			offset = offset + level;
			reader = new CSVReader(new FileReader("D:/metrxn0200/input/All.umbt"), '\t', '"');
			SortedSet<String> keySet = ruleHashRid.keySet();
			String[] moieties;

			while ((nextLine = reader.readNext()) != null) {
				if (!keySet.contains(nextLine[offset])) continue;
				TreeMap<String, Integer> frequencyMap = new TreeMap<String, Integer>();
				moieties = nextLine[offset + 12].replaceAll("[\\[\\]]", "").split("[,][ ]");
				for (String m : moieties) {
					if (m.equalsIgnoreCase("")) continue;
					key = m.replace("-", "").trim();
					if (ignoreBonds) {
						if (key.startsWith("b:")) continue;
					}
					if (!moietyIndex.contains(key)) {
						System.err.println("INCONSISTENCY for " + key + " hash = " + nextLine[offset]);
						continue;
					}
					if (ignoreBonds) {
						if (m.contains("b:")) continue;
					}

					// moietyIndex.add(key);
					if (!frequencyMap.containsKey(key)) {
						frequencyMap.put(key, 0);
					}
					if (m.startsWith("-")) {
						frequencyMap.put(key, frequencyMap.get(key) - 1);
					}
					else {
						frequencyMap.put(key, frequencyMap.get(key) + 1);
					}
				}
				reactionMoiety.put(nextLine[offset], frequencyMap);
			}
			reader.close();

			if (printSimilarity) {
				/**************/
				System.out.println("calculating similarities");
				final Object[] m_hash_key1 = metaboliteMoiety.keySet().toArray();
				final Object[] m_hash_key2 = metaboliteMoiety.keySet().toArray();
				similarity = Collections.synchronizedSortedMap(new TreeMap<String, TreeMap<String, Float>>());
				final Semaphore s = new Semaphore(THREADCOUNT);
				for (int i = 0; i < m_hash_key1.length; i++) {
					final String m_key1;
					m_key1 = (String) m_hash_key1[i];
					final int i_f = i;
					s.acquireUninterruptibly();
					new Thread() {
						public void run() {
							try {

								String m_key2;

								if (!similarity.containsKey(m_key1)) {
									TreeMap<String, Float> weightMap = new TreeMap<String, Float>();

									weightMap.put(m_key1, 1.00000000f);
									similarity.put(m_key1, weightMap);
								}
								TreeMap<String, Integer> frequencyMap = metaboliteMoiety.get(m_key1);
								Set<String> moiety_key_set1 = new HashSet<String>();
								moiety_key_set1.addAll(frequencyMap.keySet());

								Integer total1 = 0;

								for (Integer sum : (Collection<Integer>) frequencyMap.values()) {
									total1 = total1 + sum;
								}

								for (int j = i_f + 1; j < m_hash_key2.length; j++) {
									m_key2 = (String) m_hash_key2[j];

									if (!similarity.containsKey(m_key2)) {
										TreeMap<String, Float> weightMap = new TreeMap<String, Float>();
										weightMap.put(m_key2, 1.00000000f);
										similarity.put(m_key2, weightMap);
									}

									TreeMap<String, Integer> frequencyMap2 = metaboliteMoiety.get(m_key2);
									Set<String> moiety_key_set2 = new HashSet<String>();
									moiety_key_set2.addAll(frequencyMap2.keySet());

									Integer total2 = 0;
									for (Integer sum : (Collection<Integer>) frequencyMap2.values()) {
										total2 = total2 + sum;
									}

									moiety_key_set2.retainAll(moiety_key_set1);

									int commonCount = 0;

									for (String commonKeys : moiety_key_set2) {
										if (frequencyMap.get(commonKeys) <= frequencyMap2.get(commonKeys)) {
											commonCount = commonCount + 2 * frequencyMap.get(commonKeys);
										}
										else {
											commonCount = commonCount + 2 * frequencyMap2.get(commonKeys);
										}
									}
									float k;
									if ((total1 + total2 - commonCount) > largestk) {
										largestk = total1 + total2 - commonCount;
									}

									if (commonCount == 0) {
										k = 10.0000000f;
									}
									else {
										k = 1.00000000f + (float) (total1 + total2 - commonCount) / (float) (1000);
									}

									similarity.get(m_key1).put(m_key2, k);
									similarity.get(m_key2).put(m_key1, k);

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
				System.out.println("similarity calculation done");
				ruleHashMetaboliteMoiety = new TreeMap<String, TreeMap<String, Float>>();
				SortedSet<String> ruleHash = ruleHash_Mid.keySet();
				SortedSet<String> mid_set;
				String m_hash_from_rule;
				Set<String> m_hash_others_set = metaboliteMoiety.keySet();
				float sim_score = 2.00000000f;
				float prev_sim_score = 2.00000000f;
				for (String r_hash : ruleHash) {
					if (!ruleHashMetaboliteMoiety.containsKey(r_hash)) {
						ruleHashMetaboliteMoiety.put(r_hash, new TreeMap<String, Float>());
					}

					mid_set = ruleHash_Mid.get(r_hash);
					for (String id : mid_set) {
						m_hash_from_rule = mid_moietyHash.get(id);
						ruleHashMetaboliteMoiety.get(r_hash).put(m_hash_from_rule, 1.00000000f);
						for (String m_hash_other : m_hash_others_set) {
							sim_score = 20.00000000f;
							prev_sim_score = 20.00000000f;

							if (similarity.get(m_hash_other).containsKey(m_hash_from_rule)) {
								sim_score = similarity.get(m_hash_other).get(m_hash_from_rule);
							}

							if (ruleHashMetaboliteMoiety.get(r_hash).containsKey(m_hash_other)) {
								prev_sim_score = ruleHashMetaboliteMoiety.get(r_hash).get(m_hash_other);
							}

							if (sim_score < prev_sim_score) {
								ruleHashMetaboliteMoiety.get(r_hash).put(m_hash_other, sim_score);
							}

						}
					}
				}
				/******/
				/* extended */
			}

		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static void printTables() {
		try {
			String setPath = "D:/paper 6/alchemist_" + loop + "/sets/";
			String prmPath = "D:/paper 6/alchemist_" + loop + "/parameters/";
			String mapPath = "D:/paper 6/alchemist_" + loop + "/mapping/";
			// com.google.common.io.Files.copy("D:/paper 6/run_rxnBalance3.gms",
			// "D:/paper 6/alchemist_"+loop+"/run_rxnBalance3.gms");
			PrintWriter outputStream = new PrintWriter(new FileWriter(setPath + "moiety.txt"));
			outputStream.println("/");
			for (String m : moietyIndex) {
				outputStream.println("'" + m + "'");
			}
			outputStream.println("/");
			outputStream.close();

			/********/
			outputStream = new PrintWriter(new FileWriter(setPath + "moiety_extended.txt"));
			outputStream.println("/");
			for (String m : moietyIndex_extended) {
				outputStream.println("'" + m + "'");
			}
			outputStream.println("/");
			outputStream.close();

			/*******/
			SortedSet<String> keySet = moietyHashMid.keySet();
			String mid;
			outputStream = new PrintWriter(new FileWriter(setPath + "metabolite.txt"));
			outputStream.println("/");
			for (String m_hash : keySet) {
				mid = moietyHashMid.get(m_hash).first();
				outputStream.println("'" + mid + "'");
			}
			outputStream.println("/");
			outputStream.close();
			/*******/
			keySet = ruleHashRid.keySet();
			String rid;
			outputStream = new PrintWriter(new FileWriter(setPath + "rules.txt"));
			outputStream.println("/");
			for (String r_hash : keySet) {
				rid = ruleHashRid.get(r_hash).first();
				outputStream.println("'" + rid + "'");
			}
			outputStream.println("/");
			outputStream.close();
			/*******/
			Set<String> keySet3 = rid_ruleHash.keySet();

			outputStream = new PrintWriter(new FileWriter(setPath + "rules_extended.txt"));
			outputStream.println("/");
			for (String r_id : keySet3) {
				outputStream.println("'" + r_id + "'");
			}
			for (String r_id : exchangeMetabolites) {
				outputStream.println("'RXN_" + r_id + "'");
			}
			for (String m_id : moietyIndex_extended) {
				if (m_id.startsWith("cpd:")) {
					outputStream.println("'RXN_" + m_id.replace("cpd:", "") + "_cpd'");
				}

			}
			outputStream.println("/");
			outputStream.close();
			/*******/
			keySet = moietyHashMid.keySet();
			// System.out.println("size of keyset for moietyHashMid = " +
			// moietyHashMid.size());
			outputStream = new PrintWriter(new FileWriter(prmPath + "moietyMetabolite.txt"));
			outputStream.println("/");
			TreeMap<String, Integer> frequencyMap;
			Set<String> keySet2;
			for (String m_hash : keySet) {
				mid = moietyHashMid.get(m_hash).first();
				frequencyMap = metaboliteMoiety.get(m_hash);
				keySet2 = frequencyMap.keySet();
				// System.out.println("size of keyset for frequencyMap = "+keySet2.size());
				for (String m : keySet2) {
					outputStream.println("'" + m + "'.'" + mid + "'" + " " + frequencyMap.get(m));
				}
			}
			outputStream.println("/");
			outputStream.close();
			/*******/
			keySet = ruleHashRid.keySet();
			outputStream = new PrintWriter(new FileWriter(prmPath + "moietyRules.txt"));
			outputStream.println("/");
			for (String r_hash : keySet) {
				rid = ruleHashRid.get(r_hash).first();
				frequencyMap = reactionMoiety.get(r_hash);
				keySet2 = frequencyMap.keySet();
				for (String m : keySet2) {
					outputStream.println("'" + m + "'.'" + rid + "'" + " " + frequencyMap.get(m));
				}
			}
			outputStream.println("/");
			outputStream.close();
			/*******/

			keySet3 = rid_ruleHash.keySet();
			outputStream = new PrintWriter(new FileWriter(prmPath + "moietyRules_extended.txt"));
			outputStream.println("/");
			String r_h;
			for (String r_id : keySet3) {
				r_h = rid_ruleHash.get(r_id);
				frequencyMap = reactionMoiety.get(r_h);
				keySet2 = frequencyMap.keySet();
				for (String m : keySet2) {
					outputStream.println("'" + m + "'.'" + r_id + "'" + " " + frequencyMap.get(m));
				}
				HashSet<String> hashSet = equations.get(r_id);
				for (String out : hashSet) {
					outputStream.println(out);
				}
			}
			Set<String> keySet4 = exchangeMidMapping.keySet();
			for (String r_id : keySet4) {
				r_h = exchangeMidMapping.get(r_id);
				frequencyMap = reactionMoiety.get(r_h);
				keySet2 = frequencyMap.keySet();
				for (String m : keySet2) {
					outputStream.println("'" + m + "'.'RXN_" + r_id + "'" + " " + frequencyMap.get(m));
				}
			}

			for (String m_id : moietyIndex_extended) {
				if (m_id.startsWith("cpd:")) {
					outputStream.println("'" + m_id + "'.'RXN_" + m_id.replace("cpd:", "") + "_cpd' 1");
				}
			}

			outputStream.println("/");
			outputStream.close();
			/*******/
			outputStream = new PrintWriter(new FileWriter(prmPath + "exchangeFlux.txt"));
			outputStream.println("/");
			SortedSet<String> keyset = moietyHashMid.keySet();
			for (int i = 0; i < exchangeMetabolites.size(); i++) {
				mid = moietyHashMid.get(exchangeMidMapping.get(exchangeMetabolites.get(i))).first();
				outputStream.println("'" + mid + "'.'" + "RXN_" + exchangeMetabolites.get(i) + "' " + exchangeFlux.get(i));
				for (String m_hash : keyset) {
					if (!m_hash.equalsIgnoreCase(exchangeMidMapping.get(exchangeMetabolites.get(i)))) {
						mid = moietyHashMid.get(m_hash).first();
						outputStream.println("'" + mid + "'.'" + "RXN_" + exchangeMetabolites.get(i) + "' 0");
					}
				}
			}

			outputStream.println("/");
			outputStream.close();
			/********/

			/***********/
			//
			outputStream = new PrintWriter(new FileWriter(prmPath + "exchangeFluxReactionsOnly.txt"));
			outputStream.println("/");

			for (int i = 0; i < exchangeMetabolites.size(); i++) {
				outputStream.println("'RXN_" + exchangeMetabolites.get(i) + "' " + exchangeFlux.get(i));
			}
			outputStream.println("/");
			outputStream.close();
			/**********/

			outputStream = new PrintWriter(new FileWriter(prmPath + "filter.txt"));
			outputStream.println("/");
			keyset = moietyHashMid.keySet();
			for (String m_hash : keyset) {
				for (String e : exchangeMetabolites) {
					outputStream.println("'" + moietyHashMid.get(m_hash).first() + "'.'" + "RXN_" + e + "' 1");
				}
			}
			outputStream.println("/");
			outputStream.close();

			/*********/
			outputStream = new PrintWriter(new FileWriter(mapPath + "metabolite.txt"));
			keyset = moietyHashMid.keySet();
			for (String m_hash : keyset) {
				outputStream.println(moietyHashMid.get(m_hash).first() + "\t" + moietyHashMid.get(m_hash).toString());
			}
			outputStream.close();
			/********/
			outputStream = new PrintWriter(new FileWriter(mapPath + "reactions.txt"));
			keyset = ruleHashRid.keySet();
			for (String r_hash : keyset) {
				outputStream.println(ruleHashRid.get(r_hash).first() + "\t" + ruleHashRid.get(r_hash).toString());
			}
			outputStream.close();
			/*********/
			if (printReverse) {

				keySet = ruleHashRid.keySet();
				outputStream = new PrintWriter(new FileWriter(setPath + "rules_rev.txt"));
				outputStream.println("/");
				for (String r_hash : keySet) {
					rid = ruleHashRid.get(r_hash).first();
					if (rid.startsWith("RXN_")) {
						outputStream.println("'" + rid + "'");
						continue;
					}

					outputStream.println("'b_" + rid + "'");
				}
				outputStream.println("/");
				outputStream.close();
			}

			/**********/

			outputStream = new PrintWriter(new FileWriter(mapPath + "notFound.txt"));
			for (String m_id : compoundNotFound) {
				outputStream.println(m_id);
			}
			outputStream.close();
			/**********/

			if (printSimilarity) {
				/*********/
				try {
					outputStream = new PrintWriter(new FileWriter(prmPath + "similarity.txt"));

					PrintWriter outputStream2 = new PrintWriter(new FileWriter(prmPath + "prefilter.txt"));
					PrintWriter outputStream3 = new PrintWriter(new FileWriter(mapPath + "similarity_Forsql.txt"));

					outputStream.println("/");
					outputStream2.println("/");
					Set<String> r_hash_set = ruleHashMetaboliteMoiety.keySet();
					Set<String> m_hash_set_all = metaboliteMoiety.keySet();
					Set<String> m_hash_set_rule;
					Float float1;
					for (String r_hash_key : r_hash_set) {
						Set<String> m_hash_set_all_copy = new HashSet<String>();
						m_hash_set_all_copy.addAll(m_hash_set_all);
						rid = ruleHashRid.get(r_hash_key).first();
						m_hash_set_rule = ruleHashMetaboliteMoiety.get(r_hash_key).keySet();
						m_hash_set_all_copy.removeAll(m_hash_set_rule);
						for (String m_hash_key : m_hash_set_rule) {
							mid = moietyHashMid.get(m_hash_key).first();
							float1 = ruleHashMetaboliteMoiety.get(r_hash_key).get(m_hash_key);
							if (float1 < 1.02) {
								countUsable++;
							}
							outputStream.println("'" + mid + "'.'" + rid + "' " + float1);
							outputStream3.println(mid + "\t" + rid + "\t" + float1);
						}

						for (String m_hash_key : m_hash_set_all_copy) {
							mid = moietyHashMid.get(m_hash_key).first();
							outputStream2.println("'" + mid + "'.'" + rid + "' 0");
						}
					}
					outputStream.println("/");
					outputStream.close();
					outputStream2.println("/");
					outputStream2.close();
					outputStream3.close();
				}
				catch (IOException e) {
					System.err.println("IOException:");
					e.printStackTrace();
				}
				try {
					outputStream = new PrintWriter(new FileWriter(mapPath + "metaboliteSimilarity.txt"));
					Set<String> m_hash_set = similarity.keySet();
					TreeMap<String, Float> treeMap;
					for (String m_hash_key : m_hash_set) {
						mid = moietyHashMid.get(m_hash_key).first();
						treeMap = similarity.get(m_hash_key);
						Set<String> m_hash_set2 = treeMap.keySet();
						for (String m_hash_key2 : m_hash_set2) {
							outputStream.println(mid + "\t" + moietyHashMid.get(m_hash_key2).first() + "\t" + treeMap.get(m_hash_key2));
						}
					}
					outputStream.close();
				}
				catch (IOException e) {
					System.err.println("IOException:");
					e.printStackTrace();
				}
			}

		}
		catch (IOException e) {
			System.err.println("IOException:");
			e.printStackTrace();
		}

		System.out.println("largest value = " + largestk);
		System.out.println("usable variables = " + countUsable);
	}

}
