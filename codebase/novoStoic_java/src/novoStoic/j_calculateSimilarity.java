package paper6;

import org.apache.commons.lang3.StringUtils;

public class j_calculateSimilarity {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		/* 1) for all metabolites, 
		 * 	index by moietyids 
		 * 	in a TreeMap<moiety,Treemap <moiety_hash,m_count]>>
		 * 	This will be multithreaded
		 * 
		 * for each moiety, 
		 * 	loop over each moiety_hash 
		 * 		and find the difference, store it as hashmap<hash hashMap<hash:[total:difference]>> 
		 * 		 
		 * 	if already exists,
		 * 		update new sum and difference
		 * 
		 * Since these are treemaps, lexicographic ordering is already guaranteed,
		 * 
		 * for each hashmap, print in a single line, hash$total difference twice.	
		 * 
		 * 		
		 * 		 
		 * 	
		 * 
		 * 	
		 * 	
		 * 
		 *  	
		*/
		
		

	}

}
