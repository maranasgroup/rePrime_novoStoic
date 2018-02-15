package paper6;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.concurrent.Semaphore;

import org.apache.commons.math3.util.MathUtils;

import com.google.common.collect.TreeMultimap;

public class i_gamsDataSimilarityBetweenMetabolites {

	private static final int D = 50;/*dissimilarity*/
	private static final int THREADCOUNT = 7;
	private static HashSet<String> rulesIndex;
	private static HashMap<String, HashMap<String, Integer>> metabolite_p;
	private static HashMap<String, HashMap<String, Integer>> prime_m;
	// private static SortedMap<String, SortedMap<String, Integer>> m_m;
	private static HashMap<String, Integer> metabolite_t;
	private static HashSet<String> ruleMid;
	private static TreeMultimap<String, String> ruleMidMap;
	private static HashSet<String> rules;
	private static String path;
	private static PrintWriter outputStream;
	/*TODO: Every rule acts on a certain set of moieties, those moieties need to be present on the metabolite for it work, if those are not present, do not */
	public static void main(String[] args) {
		// read saj

		String OS = System.getProperty("os.name");
		if (OS.startsWith("Windows")) {
			path = "D:/paper 6/alchemist/";
		}
		else {
			path = "/home/azk172/fromWindows/Akhil/D/paper 6/alchemist/";
		}

		ruleMid = new HashSet<String>();
		ruleMidMap = TreeMultimap.create();

		String[] split;
		String[] split2;
		String inLine;
		BufferedReader inputStream;
		try {
			inputStream = new BufferedReader(new FileReader(path + "mapping/ruleToMid.map"));

			while ((inLine = inputStream.readLine()) != null) {
				inLine = inLine.replace("[", "").replace("]", "");
				split = inLine.split("\t")[1].split("[,][ ]");
				split2 = inLine.split("\t")[2].split("[,][ ]");

				for (int i = 0; i < split.length; i++) {
					split[i] = split[i].replace("cpd:", "");
					ruleMid.add(split[i]);
					for (int j = 0; j < split2.length; j++) {
						ruleMidMap.put(split2[j].replace("'", ""), split[i]);
					}
				}
			}
			inputStream.close();
		}
		catch (IOException e) {
			System.err.println("IOException:");
			e.printStackTrace();
		}

		rules = new HashSet<String>();

		try {
			inputStream = new BufferedReader(new FileReader(path + "sets/onlyRule.index"));
			while ((inLine = inputStream.readLine()) != null) {
				rules.add(inLine.replace("'", ""));
			}
			inputStream.close();
		}
		catch (IOException e) {
			System.err.println("IOException:");
			e.printStackTrace();
		}

		metabolite_p = new HashMap<String, HashMap<String, Integer>>();
		metabolite_t = new HashMap<String, Integer>();
		prime_m = new HashMap<String, HashMap<String, Integer>>();

		try {
			inputStream = new BufferedReader(new FileReader(path + "parameters/atoms_filtered.formula"));

			while ((inLine = inputStream.readLine()) != null) {
				split = inLine.split("(['][\\.]['])|(['][ ])");

				if (!metabolite_p.containsKey(split[1])) metabolite_p.put(split[1], new HashMap<String, Integer>());

				if (!prime_m.containsKey(split[0])) prime_m.put(split[0], new HashMap<String, Integer>());

				if (!metabolite_t.containsKey(split[1])) metabolite_t.put(split[1], 0);

				metabolite_p.get(split[1]).put(split[0], Integer.valueOf(split[2]));
				prime_m.get(split[0]).put(split[1], Integer.valueOf(split[2]));
				metabolite_t.put(split[1], metabolite_t.get(split[1]) + Integer.valueOf(split[2]));
			}
			inputStream.close();
		}
		catch (IOException e) {
			System.err.println("IOException:");
			e.printStackTrace();
		}

		Set<String> metabIndex1 = metabolite_p.keySet();

		final Semaphore s = new Semaphore(THREADCOUNT);
		int count = 0;
		try {
			outputStream = new PrintWriter(new BufferedWriter(new FileWriter(path + "parameters/metabolite_rule.similarity")));
			for (String metab : metabIndex1) {

				/*count++;
				if (count > 1000) {
					continue;
				}*/

				s.acquireUninterruptibly();

				final String metab_f = metab;
				new Thread() {
					public void run() {
						try {
							HashMap<String, Integer> primeCount1 = metabolite_p.get(metab_f);
							Set<String> primes = primeCount1.keySet();
							HashMap<String, Integer> dissimilarityMap = new HashMap<String, Integer>();
							for (String p : primes) {
								HashMap<String, Integer> metabCount = prime_m.get(p);
								Set<String> metabIndex2 = metabCount.keySet();
								for (String m : metabIndex2) {

									if (ruleMid.contains(m)) {

										if (!dissimilarityMap.containsKey(m)) {
											dissimilarityMap.put(m, 0);
										}
										int common = 2 * Math.min(primeCount1.get(p), metabCount.get(m));
										dissimilarityMap.put(m, dissimilarityMap.get(m) + common);
									}
								}
							}

							SortedSet<String> rxnkeySet = ruleMidMap.keySet();
							for (String rxn : rxnkeySet) {
								if (rules.contains(rxn)) {
									SortedSet<String> sortedSet = ruleMidMap.get(rxn);
									Double leastDissimilar = 1000.00;
									for (String m2 : sortedSet) {
										if (dissimilarityMap.containsKey(m2)) {
											Double dissimilarity = (double) ((double) 100 - ((double) dissimilarityMap.get(m2) / (double) (metabolite_t.get(metab_f) + metabolite_t.get(m2))) * (double) 100);
											if (leastDissimilar > dissimilarity) {
												leastDissimilar = dissimilarity;
											}
										}
									}
									if (leastDissimilar < D) {
										outputStream.println("'"+metab_f + "'.'" + rxn + "' " + Math.floor(leastDissimilar+1));
										outputStream.println("'"+metab_f + "'.'" + rxn + "_r' " + Math.floor(leastDissimilar+1));
									}
								}
							}

						}
						finally {
							s.release();
						}
					}
				}.start();
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		/* pause for all threads to close */
		while (s.availablePermits() < THREADCOUNT) {}
		
		outputStream.close();		

	}

}
