package nonEssential;

import java.util.Arrays;

public class testCrossProduct {

	public static void main(String[] args) {
		String[][] arr1 = new String[4][];
		arr1[0] = new String[] { "a1","a2" };
		arr1[1] = new String[] { "b1", "b2" };
		arr1[2] = new String[] { "c1", "c2" };
		arr1[3] = new String[] { "d" };
		String[] arr2 = new String[4];
		for (int j = 0; j < arr1[0].length; j++) {
			recurse(arr1, 0, j, arr2);
		}
	}

	private static void recurse(String[][] arr1, int i, int j, String[] arr2) {
		arr2[i] = arr1[i][j];
		i++;

		if (i == arr1.length) {
			System.out.println(Arrays.deepToString(arr2));
			return;
		}	

		for (int k_j = 0; k_j < arr1[i].length; k_j++) {
			recurse(arr1, i, k_j, arr2);
		}

	}

}
