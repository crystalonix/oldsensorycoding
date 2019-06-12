package sensoryCoding.network;

import java.util.Comparator;

public class SpikeIndexComparator2 implements Comparator<Integer> {
	private static final double MIN_COEFFICIENT_VALUE = -1000000;
	private final double[] coefficients;
	private final double[] thresholds;
	public Integer[] indexes = null;

	public SpikeIndexComparator2(double[] coefficients, double[] thresholds) {
		this.coefficients = new double[coefficients.length];
		this.thresholds = new double[thresholds.length];
		System.arraycopy(coefficients, 0, this.coefficients, 0, coefficients.length);
		System.arraycopy(thresholds, 0, this.thresholds, 0, thresholds.length);
	}

	public Integer[] createIndexArray() {
		indexes = new Integer[coefficients.length];
		for (int i = 0; i < coefficients.length; i++) {
			indexes[i] = i; // Autoboxing
		}
		return indexes;
	}

	/**
	 * This method is invoked to update the coefficients of the kernels generating spikes
	 * @param coeffs
	 * @param n
	 */
	public void updateCoefficients(double [] coeffs, int n) {
		for(int i=0; i<n ; i++) {
			coefficients[indexes[i]] = coeffs[i];
		}
		/*for(int i=n; i<coefficients.length ; i++) {
			coefficients[indexes[i]] = MIN_COEFFICIENT_VALUE;
		}*/
	}
	
	@Override
	public int compare(Integer index1, Integer index2) {
		// Autounbox from Integer to int to use as array indexes
		if (coefficients[index1] * thresholds[index1] > coefficients[index2] * thresholds[index2]) {
			return -1;
		} else if (coefficients[index1] * thresholds[index1] == coefficients[index2] * thresholds[index2]) {
			return 0;
		}
		return 1;
	}
}
