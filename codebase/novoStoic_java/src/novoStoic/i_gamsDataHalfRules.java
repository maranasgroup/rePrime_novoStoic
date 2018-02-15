package paper6;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeMap;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.LinkedHashMultimap;
import com.google.common.collect.TreeMultimap;
import com.google.common.hash.Hashing;


public class i_gamsDataHalfRules {

	public static void main(String[] args) {
		/* read all the rules from existing file, filter and create unique half rules in both directions. */
		/* same direction is twice the positive. */
		/* opposite direction or rev is twice the negative parts */
		/* index i will remain the same */
		/* e will remain the same, */
		/* E will all be positve, no negative */
		/* C will be calculated using old E, the moieties present in the difference will be changed to absolute */
		/*TODO: add cofactor information*/
		TreeMap<String, TreeMap<String, Integer>> store_atoms = new TreeMap<String, TreeMap<String, Integer>>();
		TreeMap<String, TreeMap<String, Integer>> store_bonds = new TreeMap<String, TreeMap<String, Integer>>();/*has not been used to print sbj yet*/
		try {
			BufferedReader inputStream = new BufferedReader(new FileReader("D:/paper 6/alchemist/parameters/rule.saj"));
			String inLine;
			String rid;
			String moieties;
			Integer moietyCount;

			String[] split;
			while ((inLine = inputStream.readLine()) != null) {
				split = inLine.split("['][\\.][']|['][ ]");
				moieties = split[0] + "'";
				rid = "'" + split[1] + "'";
				moietyCount = Integer.valueOf(split[2]);

				if (!store_atoms.containsKey(rid)) store_atoms.put(rid, new TreeMap<String, Integer>());

				store_atoms.get(rid).put(moieties, moietyCount);
			}
			inputStream.close();
			inputStream = new BufferedReader(new FileReader("D:/paper 6/alchemist/parameters/rule.sbj"));
			while ((inLine = inputStream.readLine()) != null) {
				split = inLine.split("['][\\.][']|['][ ]");
				moieties = split[0] + "'";
				rid = "'" + split[1] + "'";
				moietyCount = Integer.valueOf(split[2]);

				if (!store_bonds.containsKey(rid)) store_bonds.put(rid, new TreeMap<String, Integer>());

				store_bonds.get(rid).put(moieties, moietyCount);
			}
			inputStream.close();

			inputStream = new BufferedReader(new FileReader("D:/paper 6/alchemist/parameters/rule.scj"));
			String mid;
			HashMultimap<String, String> rxnCpd = HashMultimap.create();
			HashMap<String,Integer> rxnCpdStoich = new HashMap<String, Integer>();
			while ((inLine = inputStream.readLine()) != null) {
				mid = inLine.split("['][\\.][']|['][ ]")[0].replaceFirst("'", "");
				rid = inLine.split("['][\\.][']|['][ ]")[1];
				moietyCount = Integer.valueOf(inLine.split("['][\\.][']|['][ ]")[2]);
				if (moietyCount != 0) {
					rxnCpd.put(rid, mid);					
					rid = "'"+rid+"'";
					if (!rxnCpdStoich.containsKey(rid))
						rxnCpdStoich.put(rid, moietyCount);
					
					if (rxnCpdStoich.get(rid) < moietyCount) 
						rxnCpdStoich.put(rid,moietyCount);					
				}				
			}
			inputStream.close();
			
			inputStream = new BufferedReader(new FileReader("D:/paper 6/alchemist/sets/rule.ec"));
			String ec;
			HashMultimap<String, String> rxnEc = HashMultimap.create();			
			while ((inLine = inputStream.readLine()) != null) {
				rid = inLine.split("['][\\.][']|['][ ]")[0].replaceFirst("'", "");
				ec = inLine.split("['][\\.][']|['][ ]")[1].replaceFirst("'", "");
				rxnEc.put(rid, ec);								
			}
			inputStream.close();			

			PrintWriter onlyRuleIndexStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/sets/onlyRule.index"));
			PrintWriter onlyRuleSAJStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/onlyRule.saj"));
			PrintWriter onlyRuleSBJStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/onlyRule.sbj"));
			
			PrintWriter onlyRuleEcStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/sets/onlyRule.ec"));

			PrintWriter revOnlyRuleIndexStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/sets/revOnlyRule.index"));
			PrintWriter revOnlyRuleSAJStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/revOnlyRule.saj"));
			PrintWriter revOnlyRuleSBJStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/revOnlyRule.sbj"));

			PrintWriter pairsOnlyRuleIndexStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/onlyRule.pairs"));
			PrintWriter rfrbLinkIndexStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/rfrb.link"));

			PrintWriter halfRuleIndexStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/sets/halfRule.index"));
			PrintWriter halfRuleFwdStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/sets/halfRule_fwd.index"));
			PrintWriter halfRuleRevStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/sets/halfRule_rev.index"));

			PrintWriter halfRulePairStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/halfRule.pairs"));

			PrintWriter halfRuleSAJStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/halfRule.saj"));

			PrintWriter halfRuleSBJStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/halfRule.sbj"));

			
			PrintWriter ruleMidMapStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/mapping/ruleToMid.map"));
			PrintWriter ruleRidMapStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/mapping/ruleToRid.map"));
			PrintWriter ruleMetaboliteCountStream_fwd = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/fixMetabolitePerRule.fix"));
			PrintWriter ruleMetaboliteCountStream_back = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/fixMetabolitePerRule_back.fix"));
			PrintWriter ruleMetaboliteMaxStoichCountStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/fixMaxMetaboliteStoichPerRule.fix"));

			HashSet<String> ifPrinted = new HashSet<String>();
			String hashValue;
			Set<String> keySet = store_atoms.keySet();
			String halFrid;
			LinkedHashMultimap<String, String> ruleMapper = LinkedHashMultimap.create();
			LinkedHashMultimap<String, String> ruleEcMapper = LinkedHashMultimap.create();
			TreeMultimap<String, Integer> ruleMaxMetaboliteCount = TreeMultimap.create();
			HashMultimap<String, String> ruleMidMapper = HashMultimap.create();
			for (String rids : keySet) {
				halFrid = "'h_" + rids.replaceFirst("[']", "");
				TreeMap<String, Integer> treeMap = store_atoms.get(rids);
				hashValue = Hashing.md5().hashString(store_atoms.get(rids).toString(),Charset.defaultCharset() ).toString();
				ruleMapper.put(hashValue, rids);
				Set<String> ecSet = rxnEc.get(rids.replaceAll("[']", ""));
				for (String e_c : ecSet) {
					ruleEcMapper.put(hashValue, e_c);
				}
				
				Set<String> set = rxnCpd.get(rids.replaceAll("'", ""));
				for (String m : set) {
					ruleMidMapper.put(hashValue, m);
					//System.out.println(m);
				}
				ruleMaxMetaboliteCount.put(hashValue, set.size());
				
				if (!ifPrinted.contains(hashValue)) {

					onlyRuleIndexStream.println(rids);
					revOnlyRuleIndexStream.println("'" + rids.replace("'", "") + "_r'");

					pairsOnlyRuleIndexStream.println(rids + "." + "'" + rids.replace("'", "") + "_r'" + " 1");
					rfrbLinkIndexStream.println(rids + "." + "'" + rids.replace("'", "") + "_r'" +"."+ rids + "." + "'" + rids.replace("'", "") + "_r'"+" 1");

					halfRuleIndexStream.println(halFrid);
					halfRuleIndexStream.println("'" + halFrid.replace("'", "") + "_r'");

					halfRuleFwdStream.println(halFrid);
					halfRuleRevStream.println("'" + halFrid.replace("'", "") + "_r'");
					halfRulePairStream.println(halFrid + "." + "'" + halFrid.replace("'", "") + "_r'" + " 1");

					Set<String> keySet2 = treeMap.keySet();
					for (String m : keySet2) {
						if (m.startsWith("'a:")) {
							moietyCount = treeMap.get(m);
							if (moietyCount > 0) {
								halfRuleSAJStream.println(m + "." + halFrid + " " + (treeMap.get(m) * 2));
							}
							else {
								halfRuleSAJStream.println(m + ".'" + halFrid.replace("'", "") + "_r'" + " " + (treeMap.get(m) * -2));
							}

							onlyRuleSAJStream.println(m + "." + rids + " " + treeMap.get(m));
							revOnlyRuleSAJStream.println(m + "." + "'" + rids.replace("'", "") + "_r'" + " " + treeMap.get(m) * -1);

						}
						else {
							moietyCount = treeMap.get(m);
							if (moietyCount > 0) {
								halfRuleSBJStream.println(m + "." + halFrid + " " + (treeMap.get(m) * 2));
							}
							else {
								halfRuleSBJStream.println(m + ".'" + halFrid.replace("'", "") + "_r'" + " " + (treeMap.get(m) * -2));
							}

							onlyRuleSBJStream.println(m + "." + rids + " " + treeMap.get(m));
							revOnlyRuleSBJStream.println(m + "." + "'" + rids.replace("'", "") + "_r'" + " " + treeMap.get(m) * -1);

						}
					}
					ifPrinted.add(hashValue);
				}
			}
			String rids;
			for (String hash : ifPrinted) {
				rids = (String) ruleMapper.get(hash).toArray()[0];
				ruleMetaboliteCountStream_fwd.println(rids+" "+ruleMaxMetaboliteCount.get(hash).last());
				ruleMetaboliteCountStream_back.println("'" + rids.replace("'", "") + "_r'"+" "+ruleMaxMetaboliteCount.get(hash).last());
				
				ruleMetaboliteMaxStoichCountStream.println(rids+" "+rxnCpdStoich.get(rids));
				ruleMetaboliteMaxStoichCountStream.println("'" + rids.replace("'", "") + "_r'"+" "+rxnCpdStoich.get(rids));
				Set<String> set = ruleMapper.get(hash);
				for (String j_index : set) {
					ruleRidMapStream.println(rids+"."+j_index);
				}
				set = ruleEcMapper.get(hash);
				for (String e_c : set) {
					onlyRuleEcStream.println(rids+".'"+e_c+"'");
				}
			}
					
			Set<String> keySet4 = ruleMidMapper.keySet();
			for (String hash : keySet4) {
				ruleMidMapStream.println(hash+"\t"+ruleMidMapper.get(hash).toString()+"\t"+ruleMapper.get(hash).toString());
				
				
			}
			
			ruleMidMapStream.close();
			ruleRidMapStream.close();
			
			halfRuleIndexStream.close();
			halfRuleSAJStream.close();
			halfRuleSBJStream.close();

			halfRuleFwdStream.close();
			halfRuleRevStream.close();
			halfRulePairStream.close();

			onlyRuleIndexStream.close();
			onlyRuleEcStream.close();
			onlyRuleSAJStream.close();
			onlyRuleSBJStream.close();

			revOnlyRuleIndexStream.close();
			revOnlyRuleSAJStream.close();
			revOnlyRuleSBJStream.close();

			pairsOnlyRuleIndexStream.close();
			rfrbLinkIndexStream.close();
			ruleMetaboliteCountStream_fwd.close();
			ruleMetaboliteCountStream_back.close();
			ruleMetaboliteMaxStoichCountStream.close();
			
			try {
				BufferedReader inputStream1 = new BufferedReader(new FileReader("D:/paper 6/alchemist/parameters/atoms.formula"));
				String inLine1;
				HashMap<String, String> filter = new HashMap<String, String>();
				while ((inLine1 = inputStream1.readLine()) != null) {
					filter.put(inLine1.split("['][ ]")[0], inLine1);
				}
				inputStream1.close();

				PrintWriter outputStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/atoms_filtered.formula"));
				Set<String> keySet2 = filter.keySet();
				for (String string : keySet2) {
					outputStream.println(filter.get(string));
				}

				outputStream.close();

			}
			catch (IOException e) {
				System.err.println("IOException:");
				e.printStackTrace();
			}

			try {
				BufferedReader inputStream1 = new BufferedReader(new FileReader("D:/paper 6/alchemist/parameters/bonds.formula"));
				String inLine1;
				HashMap<String, String> filter = new HashMap<String, String>();
				while ((inLine1 = inputStream1.readLine()) != null) {
					filter.put(inLine1.split("['][ ]")[0], inLine1);
				}
				inputStream1.close();

				PrintWriter outputStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/bonds_filtered.formula"));
				Set<String> keySet2 = filter.keySet();
				for (String string : keySet2) {
					outputStream.println(filter.get(string));
				}

				outputStream.close();

			}
			catch (IOException e) {
				System.err.println("IOException:");
				e.printStackTrace();
			}

			HashMap<String, Integer> network = new HashMap<String, Integer>();

			try {
				inputStream = new BufferedReader(new FileReader("D:/paper 6/alchemist/parameters/exchange/rule_ex.saj"));
				while ((inLine = inputStream.readLine()) != null) {
					moieties = inLine.split("['][\\.][']|['][ ]")[0] + "'";
					moietyCount = Integer.valueOf(inLine.split("['][\\.][']|['][ ]")[2]);
					if (!network.containsKey(moieties)) network.put(moieties, 0);
					network.put(moieties, network.get(moieties) + moietyCount);
				}
				inputStream.close();
			}
			catch (IOException e) {
				System.err.println("IOException:");
				e.printStackTrace();
			}

			try {
				inputStream = new BufferedReader(new FileReader("D:/paper 6/alchemist/parameters/exchange/rule_ex.sbj"));
				while ((inLine = inputStream.readLine()) != null) {
					moieties = inLine.split("['][\\.][']|['][ ]")[0] + "'";
					moietyCount = Integer.valueOf(inLine.split("['][\\.][']|['][ ]")[2]);
					if (!network.containsKey(moieties)) network.put(moieties, 0);
					network.put(moieties, network.get(moieties) + moietyCount);
				}
				inputStream.close();
			}
			catch (IOException e) {
				System.err.println("IOException:");
				e.printStackTrace();
			}

			Integer networkTotal = 0;
			Collection<Integer> NetworkValues = network.values();
			for (Integer n : NetworkValues) {
				networkTotal = networkTotal + n;
			}

			HashMap<String, HashMap<String, Integer>> metabolites = new HashMap<String, HashMap<String, Integer>>();

			HashMap<String, HashMap<String, Integer>> halfRules = new HashMap<String, HashMap<String, Integer>>();

			try {
				inputStream = new BufferedReader(new FileReader("D:/paper 6/alchemist/parameters/atoms_filtered.formula"));
				
				while ((inLine = inputStream.readLine()) != null) {
					moieties = inLine.split("['][\\.][']|['][ ]")[0] + "'";
					mid = "'" + inLine.split("['][\\.][']|['][ ]")[1] + "'";
					moietyCount = Integer.valueOf(inLine.split("['][\\.][']|['][ ]")[2]);

					if (!metabolites.containsKey(mid)) metabolites.put(mid, new HashMap<String, Integer>());

					metabolites.get(mid).put(moieties, moietyCount);
				}
				inputStream.close();

				inputStream = new BufferedReader(new FileReader("D:/paper 6/alchemist/parameters/bonds_filtered.formula"));

				while ((inLine = inputStream.readLine()) != null) {
					moieties = inLine.split("['][\\.][']|['][ ]")[0] + "'";
					mid = "'" + inLine.split("['][\\.][']|['][ ]")[1] + "'";
					moietyCount = Integer.valueOf(inLine.split("['][\\.][']|['][ ]")[2]);

					if (!metabolites.containsKey(mid)) metabolites.put(mid, new HashMap<String, Integer>());

					metabolites.get(mid).put(moieties, moietyCount);
				}
				inputStream.close();

				inputStream = new BufferedReader(new FileReader("D:/paper 6/alchemist/parameters/halfRule.saj"));

				while ((inLine = inputStream.readLine()) != null) {
					moieties = inLine.split("['][\\.][']|['][ ]")[0] + "'";
					rid = "'" + inLine.split("['][\\.][']|['][ ]")[1] + "'";
					rid = rid.replace("_r'", "'");
					moietyCount = Integer.valueOf(inLine.split("['][\\.][']|['][ ]")[2]);

					if (!halfRules.containsKey(rid)) halfRules.put(rid, new HashMap<String, Integer>());

					halfRules.get(rid).put(moieties, moietyCount);
				}
				inputStream.close();

				inputStream = new BufferedReader(new FileReader("D:/paper 6/alchemist/parameters/halfRule.sbj"));

				while ((inLine = inputStream.readLine()) != null) {
					moieties = inLine.split("['][\\.][']|['][ ]")[0] + "'";
					rid = "'" + inLine.split("['][\\.][']|['][ ]")[1] + "'";
					rid = rid.replace("_r'", "'");
					moietyCount = Integer.valueOf(inLine.split("['][\\.][']|['][ ]")[2]);

					if (!halfRules.containsKey(rid)) halfRules.put(rid, new HashMap<String, Integer>());

					halfRules.get(rid).put(moieties, moietyCount);
				}
				inputStream.close();
			}
			catch (IOException e) {
				System.err.println("IOException:");
				e.printStackTrace();
			}

			/* calculate similarity with network */

			try {
				inputStream = new BufferedReader(new FileReader("D:/paper 6/alchemist/mapping/metabolite.map"));
				HashSet<String> metabolitesToUse = new HashSet<String>();
				while ((inLine = inputStream.readLine()) != null) {
					if (inLine.contains("kegg")) {
						metabolitesToUse.add(inLine.split("\t")[0]);
					}
				}
				inputStream.close();
				PrintWriter outputStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/sets/keggMetabolites.index"));
				for (String str : metabolitesToUse) {
					outputStream.println("'" + str + "'");
				}
				outputStream.close();

				inputStream = new BufferedReader(new FileReader("D:/paper 6/alchemist/sets/metabolite.index"));
				metabolitesToUse = new HashSet<String>();
				while ((inLine = inputStream.readLine()) != null) {
					metabolitesToUse.add(inLine.trim());
				}
				inputStream.close();

				outputStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/sets/metabolite.index"));
				for (String m : metabolitesToUse) {
					outputStream.println(m);
				}
				outputStream.close();

				inputStream = new BufferedReader(new FileReader("D:/paper 6/alchemist/sets/cpd.index"));
				HashSet<String> cpdToUse = new HashSet<String>();
				while ((inLine = inputStream.readLine()) != null) {
					inLine = inLine.replace("'cpd:", "'");
					if (metabolitesToUse.contains(inLine)) {
						cpdToUse.add(inLine.trim());
					}
				}
				inputStream.close();

				outputStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/sets/metabolitesNoReactions.index"));
				for (String m : metabolitesToUse) {
					if (!cpdToUse.contains(m)) {
						outputStream.println(m);
					}
				}
				outputStream.close();

			}
			catch (IOException e) {
				System.err.println("IOException:");
				e.printStackTrace();
			}

			/* filter for only kegg metabolites, and find the similarity of those metabolites to network metabolites. */
			/* rule_ex.saj and sbj -> dissimilarity exchange */

			try {
				PrintWriter outputStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/metabolites.NetworkSimilarity"));
				Set<String> keySet2 = metabolites.keySet();
				for (String mids : keySet2) {
					try {
						Integer metaboliteTotal = 0;
						HashMap<String, Integer> moietyCountSet = metabolites.get(mids);
						Collection<Integer> moietyCnt = moietyCountSet.values();
						for (Integer cnt : moietyCnt) {
							metaboliteTotal = metaboliteTotal + cnt;
						}

						Set<String> mKeys = network.keySet();
						int common = 0;
						for (String m : mKeys) {

							if (moietyCountSet.containsKey(m)) {
								common = common + 2*Math.min(network.get(m), moietyCountSet.get(m));
							}
						}
						if (metaboliteTotal <= 3) {
							outputStream.println(mids + " " + 1);
						}
						else
							outputStream.println(mids + " " + Math.pow((networkTotal + metaboliteTotal - common), 2));
					}
					catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				outputStream.close();
			}
			catch (IOException e) {
				System.err.println("IOException:");
				e.printStackTrace();
			}

			try {
				PrintWriter outputStream = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/halfRules.NetworkSimilarity"));
				PrintWriter outputStream2 = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/rulesOnly.NetworkSimilarity"));
				PrintWriter outputStream3 = new PrintWriter(new FileWriter("D:/paper 6/alchemist/parameters/rulesOnly_rev.NetworkSimilarity"));
				Set<String> keySet2 = halfRules.keySet();
				for (String r : keySet2) {
					Integer ruleTotal = 0;
					HashMap<String, Integer> moietyCountSet = halfRules.get(r);
					Collection<Integer> moietyCnt = moietyCountSet.values();
					for (Integer cnt : moietyCnt) {
						ruleTotal = ruleTotal + cnt;
					}

					Set<String> mKeys = network.keySet();
					int diff = 0;
					for (String m : mKeys) {
						if (moietyCountSet.containsKey(m)) {
							diff = diff + Math.min(network.get(m), (moietyCountSet.get(m) / 2));
						}
					}
					double disimilarity = Math.pow((networkTotal + ruleTotal - 2 * diff), 2);
					outputStream.println(r + " " + disimilarity);
					outputStream.println((r + " ").replace("' ", "_r'") + " " + disimilarity);

					r = r.replaceFirst("'h_", "'");
					outputStream2.println(r + " " + disimilarity);
					outputStream3.println((r + " ").replace("' ", "_r'") + " " + disimilarity);
				}
				outputStream.close();
				outputStream2.close();
				outputStream3.close();
			}
			catch (IOException e) {
				System.err.println("IOException:");
				e.printStackTrace();
			}

		}

		catch (IOException e) {
			System.err.println("IOException:");
			e.printStackTrace();
		}

	}

}
