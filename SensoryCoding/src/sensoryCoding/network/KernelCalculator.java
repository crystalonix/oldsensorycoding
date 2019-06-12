package sensoryCoding.network;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class KernelCalculator {
	public static double timeStep = ConfigurationParameters.TIME_STEP_FOR_CONVOLUTIONS;

	public static void main(String[] args) {
		// int i = 51;
		// int z = (int)(i/timeStep);
		// double k = ((double)z)*timeStep;
		// System.out.println(z);
		// System.out.println(k);//*timeStep);
		// TODO Auto-generated method stub
		/*System.out.println(bSplineByBspline(9, 17, -24));
		System.out.println(kernelMultiplication(17, 17, 0, 4, null, null)-4*bSplineByBspline(17, 17, 0)+6*bSplineByBspline(17, 17, 34));
		System.out.println();*/
		System.out.println();//*/(17, 17, 34));
	}

	/*public static double kernelByBSpline(int alpha1, int alpha2, double delta, int numberOfComponents,
			List<Double> coeff){
		double conv = 0;
		if (coeff == null) {
			for (int i = 0; i < numberOfComponents; i++) {
				//for (int j = 0; j < numberOfComponents; j++) {
					double newdelta =delta +(-alpha1*2*i);
					double c1 = Math.pow(-1, i);
					//double c2 =  Math.pow(-1, j);
					//System.out.println("i,j,newdelta"+i+","+j+","+newdelta);
					conv += c1*bSplineByBspline(alpha1, alpha2, newdelta);
				//}
			}
		} else {
			for (int i = 0; i < numberOfComponents; i++) {
				for (int j = 0; j < numberOfComponents; j++) {
					double newdelta =delta +(-alpha1*2*i);
					double c1 = coeff.get(i);
					//double c2 =  coeff2.get(j);
					double val= c1*bSplineByBspline(alpha1, alpha2, newdelta);
					conv+=val;
				}
			}
		}
		return conv;
	}

	public static double kernelMultiplication(int alpha1, int alpha2, double delta, int numberOfComponents,
			List<Double> coeff1, List<Double> coeff2) {
		double conv = 0;
		if (coeff1 == null) {
			for (int i = 0; i < numberOfComponents; i++) {
				for (int j = 0; j < numberOfComponents; j++) {
					double newdelta =delta +(-alpha1*2*i + alpha2*2*j);
					double c1 = Math.pow(-1, i);
					double c2 =  Math.pow(-1, j);
					//System.out.println("i,j,newdelta"+i+","+j+","+newdelta);
					conv += c1*c2*bSplineByBspline(alpha1, alpha2, newdelta);
				}
			}
		} else {
			for (int i = 0; i < numberOfComponents; i++) {
				for (int j = 0; j < numberOfComponents; j++) {
					double newdelta =delta +(-alpha1*2*i + alpha2*2*j);
					double c1 = coeff1.get(i);
					double c2 =  coeff2.get(j);
					double val= c1*c2*bSplineByBspline(alpha1, alpha2, newdelta);
					conv+=val;
				}
			}
		}
		return conv;
	}
	public static double kernelDifferntialMultiplication(int alpha1, int alpha2, double delta, int numberOfComponents,
			List<Double> coeff1, List<Double> coeff2) {
		double conv = 0;
		if (coeff1 == null) {
			for (int i = 0; i < numberOfComponents; i++) {
				for (int j = 0; j < numberOfComponents; j++) {
					double newdelta =delta +(-alpha1*2*i + alpha2*2*j);
					double c1 = Math.pow(-1, i);
					double c2 =  Math.pow(-1, j);
					conv += c1*c2*bsplineByDbSpline(alpha1, alpha2, newdelta);
				}
			}
		} else {
			for (int i = 0; i < numberOfComponents; i++) {
				for (int j = 0; j < numberOfComponents; j++) {
					double newdelta =delta +(-alpha1*2*i + alpha2*2*j);
					double c1 = coeff1.get(i);
					double c2 =  coeff2.get(j);
					conv += c1*c2*bsplineByDbSpline(alpha1, alpha2, newdelta);
				}
			}
		}
		return conv;
	}
*/
	/**
	 * This 2-D array stores the basic b-spline components of the kernel signals
	 */
	public double[][] basicExpandedComponentBSplines;

	public void initKernelCalculator(){
		basicExpandedComponentBSplines = new double[ConfigurationParameters.numberOfKernels][];
		for(int i=0; i< ConfigurationParameters.numberOfKernels; i++){
			double expansionFactor = 1;
			if(ConfigurationParameters.NATURE_OF_KERNEL_SPREAD==1){
				expansionFactor = 1+i*ConfigurationParameters.FREQUENCY_SCALING_FACTOR;
			}else{
				throw new IllegalArgumentException("kernel calculator currently does not support geometric kernels");
			}
			double lengthOfThisBspline = ConfigurationParameters.lengthOfBasicBSpline*expansionFactor/timeStep;
			basicExpandedComponentBSplines[i] = new double[(int)(lengthOfThisBspline) + 1];
			for(int j =0; j<basicExpandedComponentBSplines[i].length; j++){
				basicExpandedComponentBSplines[i][j] = getValueForTheChosenBSpline(j*(3/lengthOfThisBspline));
			}
		}
	}

	/**
	 * This method displays the basic bsplines
	 */
	public void drawBsplines(){
		List<XYSeries> BSplinesDisplayData = new ArrayList<>();
		for(int i =0; i< ConfigurationParameters.numberOfKernels; i++){
			BSplinesDisplayData.add(new Signal(basicExpandedComponentBSplines[i]).getSignalDisplayData("Bspline#:"+i));
		}
		Utilities.displaySetOfSignals(BSplinesDisplayData, "basic bsplines", "time", "value");
	}

	public double bSplineByBspline(int kernelIndex1, int kernelIndex2, double delta) {
		double conv = 0;
		// signal 1
		double[] signal1= basicExpandedComponentBSplines[kernelIndex1];
/*		double[] signal1 = new double[(int) ((3 * alpha1) / timeStep) + 1];
		for (int i = 0; i < signal1.length; i++) {
			signal1[i] = getValueForTheChosenBSpline((i * timeStep) / alpha1);
		}
*/		// signal 2
		double[] signal2 = basicExpandedComponentBSplines[kernelIndex2];//new double[(int) ((3 * alpha2) / timeStep) + 1];
/*		for (int i = 0; i < signal2.length; i++) {
			signal2[i] = getValueForTheChosenBSpline((i * timeStep) / alpha2);
		}
*/
		if (delta > 0) {
			for (int i = 0; i < signal2.length; i++) {
				int j = i + (int) (delta / timeStep);
				if (j >= signal1.length) {
					break;
				}
				conv += signal1[j] * signal2[i] * timeStep*ConfigurationParameters.TIME_STEP;
			}
		} else {
			for (int i = 0; i < signal1.length; i++) {
				int j = i + (int) (-delta / timeStep);
				if (j >= signal2.length) {
					break;
				}
				conv += signal1[i] * signal2[j] * timeStep* ConfigurationParameters.TIME_STEP;
			}
		}
		return conv;
	}

	public double bsplineByDbSpline(int kernelIndex1, int kernelIndex2, double delta) {
		double conv = 0;
		// signal 1
		double[] signal1 = basicExpandedComponentBSplines[kernelIndex1];//new double[(int) ((3 * alpha1) / timeStep) + 1];
/*		for (int i = 0; i < signal1.length; i++) {
			// System.out.println("i here:"+i*timeStep);
			// System.out.println("i here:"+(i*timeStep)/alpha1);
			signal1[i] = getValueForTheChosenBSpline((i * timeStep) / alpha1);
		}
*/		// signal 2
		double[] signalDiff = basicExpandedComponentBSplines[kernelIndex2];//new double[(int) ((3 * alpha2) / timeStep) + 1];
/*		for (int i = 0; i < signal2.length; i++) {
			signal2[i] = getValueForTheChosenBSpline((i * timeStep) / alpha2);
		}
*/		double [] signal2 = new double[signalDiff.length];
		for (int i = 0; i < signal2.length - 1; i++) {
			signal2[i] = (signalDiff[i + 1] - signalDiff[i]) / (timeStep*ConfigurationParameters.TIME_STEP);
		}

		signal2[signal2.length - 1] = signal2[signal2.length - 2];

		if (delta > 0) {
			for (int i = 0; i < signal2.length; i++) {
				int j = i + (int) (delta / timeStep);
				if (j >= signal1.length) {
					break;
				}
				conv += signal1[j] * signal2[i] * timeStep *ConfigurationParameters.TIME_STEP;
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
	}

	/**
	 * This method computes the value of the chosen BSpline basis function at a
	 * given time
	 *
	 * @param i
	 * @return
	 */
	private static double getValueForTheChosenBSpline(double t) {
		if (t < 0 || t > 3) {
			return 0;// throw new IllegalArgumentException("Time provided for
						// computing in the basis function is out of range"+t);
		} else {
			if (t < 1) {
				return t * t / 2;
			} else if (t < 2) {
				return (-2 * t * t + 6 * t - 3) / 2;
			} else {
				return (3 - t) * (3 - t) / 2;
			}
		}
	}

	public double calculateSignalKernelCompConvolution(Signal thisSignal, int kernelIndex, int startTime) {
		double t = startTime;
		double total = 0;
		double []thisExpandedBSpline = basicExpandedComponentBSplines[kernelIndex];
		for(int i=0; i<thisExpandedBSpline.length; i++){
			total += thisExpandedBSpline[i]*thisSignal.getSignalValue(t)*timeStep;
			t+=timeStep;
		}
		return total;
	}
}
