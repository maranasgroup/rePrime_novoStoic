package primeR;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigInteger;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.TreeBidiMap;

public class primes {

	public static HashMap<Integer, BigInteger> getPrimesAsHash(String primeFolderPath) throws IOException {
		//TreeMap<Integer, Integer> primesKey = new TreeMap<Integer, Integer>();
		//TreeSet<BigInteger> primesKey = new TreeSet<BigInteger>();
		HashMap<Integer, BigInteger> primeHash = new HashMap<Integer, BigInteger>();
		//HashMap<Integer,BigInteger> primeHash =  new HashMap<Integer, BigInteger>();		
		
		//File f = new File(primeFolderPath);
		
		BufferedReader inputStream;
		int index = 0;
		primeHash.put(0, new BigInteger("1"));
		for (int i = 1; i <= 1; i++) {
			File f = new File(primeFolderPath+"primes"+i+".txt");
			if (f.exists()) {
				inputStream = new BufferedReader(new FileReader(f));
				String inLine;
				int cnt = 0; 
				while ((inLine = inputStream.readLine()) != null) {
				
					cnt++;
					if (cnt < 2) continue;
					String[] split = inLine.trim().split("\\b");
					for (int j = 0; j < split.length; j++) {
						if (!split[j].trim().equalsIgnoreCase("")) {
							index++;							
							primeHash.put(index, new BigInteger(split[j].trim()));							
						}
					}
				}
				inputStream.close();
			}
		}		
		System.out.println("total number of primes in the dictionary " + primeHash.size());
		return primeHash;		
	}
	public static TreeSet<Integer> getSortedPrime(){
		try {
			BufferedReader inputStream = new BufferedReader(new FileReader(
					"primes.txt"));
			String inLine = null;
			TreeSet<Integer> primesKey = new TreeSet<Integer>();
			
			String[] split;
			try {
				while ((inLine = inputStream.readLine().trim()) != null) {					
					split = inLine.split("\\b");
					for (int i = 0; i < split.length; i++) {
						if (!split[i].trim().equalsIgnoreCase("")) {
							primesKey.add(Integer.parseInt(split[i].trim()));
						}						
					}
				}
			} catch (NullPointerException e) {
				
			}
			inputStream.close();
			
			return primesKey;			
		} catch (IOException e) {

			System.out.println("IOException:");
			e.printStackTrace();

		}
		return null;
	}
	
	public static LinkedList<Integer> getSortedPrimeAsList(){
		LinkedList<Integer> returnVal = new LinkedList<>();
		returnVal.addAll(getSortedPrime());		
		return returnVal;
	}
	
	public static HashMap<BigInteger,Integer> getPrimeIndex(HashMap<Integer,BigInteger> primeHash){
		HashMap<BigInteger,Integer> returnVal = new HashMap<BigInteger,Integer>();
		Set<Integer> keySet = primeHash.keySet();
		for (Integer index : keySet) {
			returnVal.put(primeHash.get(index), index);
		}
		return returnVal;
	}
}
