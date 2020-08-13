package sensoryCoding.network;

public class NeoKernelCalculator {
	public static double timeStep = 1;// ConfigurationParameters.TIME_STEP_FOR_CONVOLUTIONS;

	public static void main(String[] args) {
		System.out.println();
	}

	/**
	 * This 2-D array stores the basic b-spline components of the kernel signals
	 */
	public double[][] basicExpandedComponentBSplines;

	/**
	 * Initialize the kernel calculator with the basic components used in the kernel
	 * manager
	 * 
	 * @param basicComps
	 */
	public void initKernelCalculator(double[][] basicComps) {
		this.basicExpandedComponentBSplines = basicComps;
		// TODO: use the time_step factor here
	}

	public double bSplineByBspline(int kernelIndex1, int kernelIndex2, double delta) {
		double conv = 0;
		double[] signal1 = basicExpandedComponentBSplines[kernelIndex1];
		double[] signal2 = basicExpandedComponentBSplines[kernelIndex2];
		if (delta > 0) {
			for (int i = 0; i < signal2.length; i++) {
				int j = i + (int) (delta / timeStep);
				if (j >= signal1.length) {
					break;
				}
				conv += signal1[j] * signal2[i] * timeStep * ConfigurationParameters.TIME_STEP;
			}
		} else {
			for (int i = 0; i < signal1.length; i++) {
				int j = i + (int) (-delta / timeStep);
				if (j >= signal2.length) {
					break;
				}
				conv += signal1[i] * signal2[j] * timeStep * ConfigurationParameters.TIME_STEP;
			}
		}
		return conv;
		// throw new NotImplementedException();
	}

	public double bsplineByDbSpline(int kernelIndex1, int kernelIndex2, double delta) {
		throw new IllegalArgumentException("Not implemented yet");
	}

	public double calculateSignalKernelCompConvolution(Signal thisSignal, int kernelIndex, int startTime) {
		double t = startTime;
		double total = 0;
		double[] thisExpandedBSpline = basicExpandedComponentBSplines[kernelIndex];
		for (int i = 0; i < thisExpandedBSpline.length; i++) {
			total += thisExpandedBSpline[i] * thisSignal.getSignalValue(t) * timeStep;
			t += timeStep;
		}
		return total;
	}
}
