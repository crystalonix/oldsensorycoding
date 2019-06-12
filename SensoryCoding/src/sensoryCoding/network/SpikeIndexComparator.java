package sensoryCoding.network;

import java.util.Comparator;

public class SpikeIndexComparator implements Comparator<Integer> {
	private final double[] coefficients;
	private final double[] thresholds;
	public Integer[] indexes = null;

	public SpikeIndexComparator(double[] coefficients, double[] thresholds) {
		this.coefficients = coefficients;
		this.thresholds = thresholds;
	}

	public Integer[] createIndexArray() {
		indexes = new Integer[coefficients.length];
		for (int i = 0; i < coefficients.length; i++) {
			indexes[i] = i; // Autoboxing
		}
		return indexes;
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
