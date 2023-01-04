package AlgorithmsSolver;

import java.util.ArrayList;
import java.util.Collections;

public class AlgoSolver {

	public static  ArrayList<Integer> increasingArrayList(ArrayList<Integer> numbers) {
		int max = 0;
		ArrayList<Integer> tmp = new ArrayList<>();
        for (int i = 0; i < numbers.size(); i++) {
              tmp.add(0);
        }

        ArrayList<Integer> result = new ArrayList<>();
        for (int i = 1; i < numbers.size(); i++) {
              if (numbers.get(i) > numbers.get(i - 1)) {
                    tmp.set(i, tmp.get(i - 1) + 1);
              } else tmp.set(i, 0);
              {
                    max = Math.max(tmp.get(i), max);
              }
        }
        for (int i = 0; i < numbers.size(); i++) {
              if (tmp.get(i) == max) {
                    for (int j = i - tmp.get(i); j <= i; j++) {
                          result.add(numbers.get(j));
                    }
              }
        }
		return result;
	}
	public static  ArrayList<Integer> dencreasingArrayList(ArrayList<Integer> numbers){
        int max = 0;
        ArrayList<Integer> tmp = new ArrayList<>();
        for (int i = 0; i < numbers.size(); i++) {
              tmp.add(0);
        }
        ArrayList<Integer> result = new ArrayList<>();
        for (int i = 1; i < numbers.size(); i++) {
              if (numbers.get(i) < numbers.get(i - 1)) {
                    tmp.set(i, tmp.get(i - 1) + 1);
              } else tmp.set(i, 0);
              {
                    max = Math.max(tmp.get(i), max);
              }
        }

        if(max == 0) 
        {
        	result = null;
        	}
        else {
        	for (int i = 0; i < numbers.size(); i++) {
        		if (tmp.get(i) == max) {
                    	for (int j = i - tmp.get(i); j <= i; j++) {
                    		result.add(numbers.get(j));
                    	}
        		}
              }
        }
	return result;
	}
	public static float median (ArrayList<Integer> numbers) {
		 float median = 0;
		 Collections.sort(numbers);
         if (numbers.size() % 2 == 0) {
        	 median = ((float) numbers.get(Math.round(((float) numbers.size()) / 2) ) + ((float) numbers.get(Math.round(((float) numbers.size()) / 2) - 1))) / 2;
         } else {
        	 median = numbers.get(numbers.size()/ 2);


         }
		return median;
	}
	public static void main(String[] args)
 {
		ArrayList<Integer> testArrayList = new ArrayList<>();
		for (int i = 0; i<10; i++) {
			testArrayList.add(i);
		}

		ArrayList<Integer> result = increasingArrayList(testArrayList);
		System.out.println(result);
		result = dencreasingArrayList(testArrayList);
		System.out.println(result);
		float median = median(testArrayList);
		System.out.println(median);
}
}