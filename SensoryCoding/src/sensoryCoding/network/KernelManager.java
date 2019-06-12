package sensoryCoding.network;

import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class KernelManager {
	private int numberOfKernels;
	// denotes the number of component basis functions the kernel is composed of
	private int numberOfKernelComponents;
	// this matrix stores the coefficients of the basis functions for each
	// kernel
	public double[][] kernelCoefficients;
	// this matrix stores the threshold values corresponding to each kernel
	private double[][] thresholdValues;
	// this indicates the initial value at which the threshold is set
	private double initialThresHoldValue;
	// this variable denotes the length of each component signal
	private int lengthOfSignal;
	// this stores all the required kernel signals
	private List<Signal> kernelsCache;
	// this stores all the required kernel signals revrted in time
	private List<Signal> invertedKernelsCache;
	// this stores all the required kernel differential signals
	private List<Signal> kernelDifferentialsCache;
	// this stores all the required kernel component differential signals
	private List<Signal> componentDifferentialsCache;
	// this stores the component BSpline functions
	private List<List<Signal>> basisComponentsCache;
	// this stores all the inner product of two time-shifted kernel
	// components(always chooses the 0th component for each kernel)
	public Signal[][] componentByComponentSignalConvolutionCache;
	// this stores all the inner product of one time-shifted kernel component
	// and another time-shifted kernel component differential signal (always
	// chooses the 0th component for each kernel)
	public Signal[][] componentByDiffComponentSignalConvolutionCache;
	// This stores the length of component b-splines
	public double[] componentBsplineLengths;
	// This stores the shifts of component b-splines
	public double[][] componentBsplineShifts;
	// component BSpline function without time shift
	// and without scaling used to generate the kernels
	private Signal basicBSplineSignal;
	// denotes the length of underlying basic BSpline function inclusive of both
	// ends
	private int lengthOfBasicBSpline = ConfigurationParameters.lengthOfBasicBSpline;
	// kernel calculator used to facilitate kernel related computations
	public KernelCalculator kernelCalc = null;
	/**
	 * stores the squared norm of the bspline edges
	 */
	private static final double A_INTEGRAL = 1.0 / 30.0;

	/**
	 * stores the squared norm of the bspline middle
	 */
	private static final double B_INTEGRAL = 3.0 / 20.0;
	/**
	 * stores the squared norm of the bspline cross
	 */
	private static final double C_INTEGRAL = 1.0 / 180.0;

	/**
	 * Stores the length of component b-splines of each kernel
	 */
	public double[] kernelBsplineLengths = new double[ConfigurationParameters.numberOfKernels];

	/**
	 * Stores the lengths of each kernels
	 */
	public double[] kernelLengths = new double[ConfigurationParameters.numberOfKernels];

	/**
	 * Constructor for the kernel manager
	 * 
	 * @param numberOfKernels
	 * @param numberOfKernelComponents
	 * @param lengthOfSignal
	 * @param shouldNormalizeKernels
	 * @throws Exception
	 */
	public KernelManager(int numberOfKernels, int numberOfKernelComponents, int lengthOfSignal,
			boolean shouldNormalizeKernels) throws Exception {

		kernelCalc = new KernelCalculator();
		kernelCalc.initKernelCalculator();
		this.numberOfKernelComponents = numberOfKernelComponents;
		this.numberOfKernels = numberOfKernels;
		this.lengthOfSignal = lengthOfSignal;
		this.initialThresHoldValue = ConfigurationParameters.initialThresHoldValue;
		initializeKernelCoefficients();
		initializeThresholdValues();

		loadBasisComponents();

		if (shouldNormalizeKernels) {
			normalizeAllKernels(false);
		}

		loadKernelCache();
		loadTimeInvertedKernelCache();
		// TODO: check if we really need kernel differentials and kernel comp
		// differentials evaluated
		loadDifferentialKernels();
		loadDifferentialKernelComponents();
		loadComponentByComponentSignalConvolutionCache();
		loadComponentByDifferentialComponentSignalConvolutionCache();
		if (ConfigurationParameters.USER_MODE.equals(MODE.DISPLAY_MODE)) {
			displayKernels();
		}
	}

	/**
	 * Constructor for the kernel manager
	 * 
	 * @param numberOfKernels
	 * @param numberOfKernelComponents
	 * @param lengthOfSignal
	 * @throws Exception
	 */
	public KernelManager(int numberOfKernels, int numberOfKernelComponents, int lengthOfSignal) throws Exception {
		this(numberOfKernels, numberOfKernelComponents, lengthOfSignal, true);
	}

	/**
	 * Constructor for the kernel manager with default parameters
	 * 
	 * @param shouldNormalizeKernels
	 * @throws Exception
	 */
	public KernelManager(boolean shouldNormalizeKernels) throws Exception {
		this(ConfigurationParameters.numberOfKernels, ConfigurationParameters.numberofKernelComponents,
				ConfigurationParameters.lengthOfComponentSignals, shouldNormalizeKernels);
	}

	public void loadTimeInvertedKernelCache() throws Exception {
		invertedKernelsCache = new ArrayList<>();
		for (int i = 0; i < numberOfKernels; i++) {
			double[] data = new double[lengthOfSignal];
			Signal result = new Signal(data, 0, 0);
			for (int j = 0; j < numberOfKernelComponents; j++) {
				Signal comp = getBasisForKernelComp(i, j);
				result = SignalUtils.addTwoSignals(result,
						SignalUtils.scalarMultiply(
								SignalUtils.shiftSignal(comp, -(comp.getEndTime() + comp.getStartTime())),
								kernelCoefficients[i][j]));
			}
			invertedKernelsCache.add(result);
		}
	}

	/**
	 * This method will store all the kernel signals into files
	 * 
	 * @param folderPath
	 * @throws Exception
	 */
	public void saveKernelsToFile(String folderPath) throws Exception {
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			Signal kernelSignal = getKernel(i);
			SignalUtils.storeSignal(folderPath + "kernel-" + i + ".txt", kernelSignal, false);
		}
	}

	public void showBasicBsplines() throws Exception {
		List<XYSeries> allBSplines = new ArrayList<>();
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			allBSplines.add(getBasisForKernelComp(i, 0).getSignalDisplayData("bspline#:" + i));
		}
		Utilities.displaySetOfSignals(allBSplines, "basic bsplines", "time", "value");
	}

	public void initKernelCoefficients(double[][] kernelCoeffs) throws Exception {
		this.kernelCoefficients = kernelCoeffs;
		loadKernelCache();
		loadTimeInvertedKernelCache();
	}

	/**
	 * Is invoked every time kernel coefficients are modified
	 * 
	 * @throws Exception
	 */
	public void updateCache() throws Exception {
		loadKernelCache();
		loadTimeInvertedKernelCache();
	}

	/**
	 * displays the kernels in a chart
	 */
	public void displayKernels(String title) {
		XYSeriesCollection kernelDataSet = new XYSeriesCollection();
		for (int i = 0; i < numberOfKernels; i++) {
			// if(i%10==9||i==0)
			kernelDataSet.addSeries(kernelsCache.get(i).getSignalDisplayData("kernel" + i));

		}
		// JFreeChart chart = ChartFactory.createXYLineChart(, , "value",
		// kernelDataSet);
		JFreeChart linechart = ChartFactory.createXYLineChart(title, "time", "value", kernelDataSet);
		ChartFrame chartFrame = new ChartFrame("plot potentials", linechart);
		chartFrame.setVisible(true);
	}

	/**
	 * displays the kernels in a chart
	 */
	public void displaySelectedKernels(int[] kernelsList) {
		XYSeriesCollection kernelDataSet = new XYSeriesCollection();
		for (int i = 0; i < kernelsList.length; i++) {
			int kernelIndex = kernelsList[i];
			kernelDataSet.addSeries(kernelsCache.get(kernelIndex).getSignalDisplayData("kernel" + kernelIndex));

		}
		// JFreeChart chart = ChartFactory.createXYLineChart(, , "value",
		// kernelDataSet);
		JFreeChart linechart = ChartFactory.createXYLineChart("kernels chart", "time", "value", kernelDataSet);
		ChartFrame chartFrame = new ChartFrame("plot potentials", linechart);
		chartFrame.setVisible(true);
	}

	/**
	 * displays the kernels in a chart
	 */
	public void displayKernels() {
		XYSeriesCollection kernelDataSet = new XYSeriesCollection();
		for (int i = 0; i < numberOfKernels; i++) {
			// if(i%10==9||i==0)
			kernelDataSet.addSeries(kernelsCache.get(i).getSignalDisplayData("kernel" + i));

		}
		// JFreeChart chart = ChartFactory.createXYLineChart(, , "value",
		// kernelDataSet);
		JFreeChart linechart = ChartFactory.createXYLineChart("kernels chart", "time", "value", kernelDataSet);
		ChartFrame chartFrame = new ChartFrame("plot potentials", linechart);
		chartFrame.setVisible(true);
	}

	/**
	 * displays the kernels in a chart
	 */
	public void displayInvertedKernels() {
		XYSeriesCollection kernelDataSet = new XYSeriesCollection();
		for (int i = 0; i < numberOfKernels; i++) {
			kernelDataSet.addSeries(invertedKernelsCache.get(i).getSignalDisplayData("inverted kernel" + i));
		}
		JFreeChart linechart = ChartFactory.createXYLineChart("kernels chart", "time", "value", kernelDataSet);
		ChartFrame chartFrame = new ChartFrame("plot potentials", linechart);
		chartFrame.setVisible(true);
	}

	/**
	 * This method loads all the basis components of the kernels
	 */
	private void loadBasisComponents() {
		kernelLengths = new double[ConfigurationParameters.numberOfKernels];
		componentBsplineLengths = new double[ConfigurationParameters.numberOfKernels];
		componentBsplineShifts = new double[ConfigurationParameters.numberOfKernels][ConfigurationParameters.numberofKernelComponents];
		basisComponentsCache = new ArrayList<>();
		for (int i = 0; i < numberOfKernels; i++) {
			List<Signal> componentSignals = new ArrayList<>();
			double expansionFactor = 1;
			if (ConfigurationParameters.NATURE_OF_KERNEL_SPREAD == 1) {
				expansionFactor = (ConfigurationParameters.FREQUENCY_SCALING_FACTOR * (i) + 1);
			} else if (ConfigurationParameters.NATURE_OF_KERNEL_SPREAD == 2) {
				expansionFactor = Math.pow(ConfigurationParameters.FREQUENCY_SCALING_FACTOR, i);
			}
			double thisBsplineLength = expansionFactor * lengthOfBasicBSpline;
			componentBsplineLengths[i] = thisBsplineLength;
			// put the kernel length into the array
			kernelBsplineLengths[i] = thisBsplineLength;
			int shift = 0;
			for (int j = 0; j < numberOfKernelComponents; j++) {
				shift = (int) ((j * thisBsplineLength) * (1.0 - ConfigurationParameters.OVERLAP_FACTOR));
				componentBsplineShifts[i][j] = shift;
				Signal componentSignal = getScaledBasicBSplineSignal(expansionFactor);
				Signal shiftedSignal = SignalUtils.shiftSignal(componentSignal, shift);
				// Signal scaledSignal = SignalUtils.scaleSignal(shiftedSignal,
				// );
				// Signal debugSignal = SignalUtils.scaleSignal(shiftedSignal,
				// 1/(i+1));
				componentSignals.add(shiftedSignal);
			}
			kernelLengths[i] = shift + thisBsplineLength;
			basisComponentsCache.add(componentSignals);
		}
	}

	/**
	 * This method initializes the kernel cache
	 *
	 * @throws Exception
	 */
	public void loadKernelCache() throws Exception {
		kernelsCache = new ArrayList<>();
		for (int i = 0; i < numberOfKernels; i++) {
			double[] data = new double[lengthOfSignal];
			Signal result = new Signal(data, 0, 0);
			for (int j = 0; j < numberOfKernelComponents; j++) {
				result = SignalUtils.addTwoSignals(result,
						SignalUtils.scalarMultiply(getBasisForKernelComp(i, j), kernelCoefficients[i][j]));
			}
			kernelsCache.add(result);
		}
	}

	/**
	 * TODO: This need not be stored at its entire length This method initializes
	 * the threshold values at which each kernel generates a spike
	 */
	private void initializeThresholdValues() {
		thresholdValues = new double[numberOfKernels][lengthOfSignal];
		for (int i = 0; i < numberOfKernels; i++) {
			for (int j = 0; j < lengthOfSignal; j++) {
				thresholdValues[i][j] = initialThresHoldValue;
			}
		}
	}

	/**
	 * This method initializes the kernel manager with required number of kernels
	 * and kernel components and individual time shifts and so on
	 */
	public void initializeKernelCoefficients() {
		kernelCoefficients = new double[numberOfKernels][numberOfKernelComponents];
		// iterate through the elements of the array and initialize

		for (int j = 0; j < kernelCoefficients.length; j++) {
			double[] kernelArray = kernelCoefficients[j];
			boolean toggle = true;
			for (int i = 0; i < kernelArray.length; i++) {
				double p = 1;
				if (ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS) {
					p = Math.random();
				}
				if (toggle) {
					kernelArray[i] = p * ConfigurationParameters.initialCoefficientValue
							/ (1 + j * ConfigurationParameters.FREQUENCY_SCALING_FACTOR);
				} else {
					kernelArray[i] = -p * ConfigurationParameters.initialCoefficientValue
							/ (1 + j * ConfigurationParameters.FREQUENCY_SCALING_FACTOR);
				}
				toggle = !toggle;
			}
		}
	}

	/**
	 * Given the index of the kernel this method returns the value of the kernel at
	 * a particular instant of time
	 *
	 * @param time        at which we want to calculate the value of the kernel
	 * @param kernelIndex index of the kernel whose value we want to calculate, it
	 *                    starts with 0
	 * @return
	 * @throws Exception
	 */
	public double getKernelValue(int time, int kernelIndex) throws Exception {
		if (kernelsCache != null) {
			return kernelsCache.get(kernelIndex).getSignalValue(time);
		}
		throw new Exception("kernel value has not been loaded yet");
		// double total =0;
		// for(int i=0; i<numberOfKernelComponents; i++){
		// total+=
		// kernelCoefficients[kernelIndex][i]*getBasisValueForKernelComp(time,
		// kernelIndex, i);
		// }
		// return total;
	}

	/**
	 * returns the value of the basis function corresponding to a kernel comp. at a
	 * given time instant
	 *
	 * @param t
	 * @param thisKernelIndex
	 * @param thisKernelComponent
	 * @return
	 * @throws Exception
	 */
	public double getBasisValueForKernelComp(int t, int thisKernelIndex, int thisKernelComponent) throws Exception {
		if (basisComponentsCache != null) {
			return basisComponentsCache.get(thisKernelIndex).get(thisKernelComponent).getSignalValue(t);
		}
		throw new Exception("basis function's value has not been loaded yet");
		// int timeShiftForThisComponent =
		// (thisKernelComponent*3*numberOfKernels)/(thisKernelIndex + 1);
		// if(t<timeShiftForThisComponent || t> timeShiftForThisComponent+
		// 3*numberOfKernels/(thisKernelIndex+1)){
		// return 0;
		// }else{
		// int shiftedTime = t - timeShiftForThisComponent;
		// return
		// getValueForTheChosenBSpline(shiftedTime*(thisKernelIndex+1)/numberOfKernels);
		// }
	}

	/**
	 * This method computes the value of the chosen BSpline basis function at a
	 * given time
	 *
	 * @param i
	 * @return
	 */
	private double getValueForTheChosenBSpline(double t) {
		if (t < 0 || t > 3) {
			throw new IllegalArgumentException("Time provided for computing in the basis function is out of range" + t);
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

	/**
	 * returns the component BSpline signal without time shift and without scaling
	 * used to generate the kernels
	 *
	 * @return
	 */
	public Signal getBasicBSplineSignal() {
		if (basicBSplineSignal != null) {
			return basicBSplineSignal;
		}
		double[] data = new double[lengthOfSignal];
		double scalingFactor = 3.0 / lengthOfBasicBSpline;
		for (int i = 0; i <= lengthOfBasicBSpline; i++) {
			data[i] = getValueForTheChosenBSpline(i * scalingFactor);
		}
		basicBSplineSignal = new Signal(data, 0, lengthOfBasicBSpline, 0);
		return basicBSplineSignal;
	}

	/**
	 * Returns the scaled basic BSpline
	 *
	 * @param alpha factor by which the BSpline is expanded
	 * @return scaled version of the BSpline signal
	 */
	public Signal getScaledBasicBSplineSignal(double alpha) {
		double[] data = new double[lengthOfSignal];
		// double scalingFactor = (3.0)/(((double)lengthOfBasicBSpline)*alpha);
		for (int i = 0; i <= lengthOfBasicBSpline * alpha; i++) {
			data[i] = getValueForTheChosenBSpline((i * 3.0) / (lengthOfBasicBSpline * alpha));
		}
		basicBSplineSignal = new Signal(data, 0, (int) (lengthOfBasicBSpline * alpha), 0);
		return basicBSplineSignal;
	}

	/**
	 * returns the basis function corresponding to a kernel comp. as a signal array
	 *
	 * @param thisKernelIndex
	 * @param thisKernelComponent
	 * @return
	 * @throws Exception
	 */
	public Signal getBasisForKernelComp(int thisKernelIndex, int thisKernelComponent) throws Exception {
		if (basisComponentsCache != null) {
			return basisComponentsCache.get(thisKernelIndex).get(thisKernelComponent);
		}
		// double[] thisBasisSignal = new double[lengthOfSignal];
		// for (int i = 0; i < lengthOfSignal; i++) {
		// thisBasisSignal[i] = getBasisValueForKernelComp(i, thisKernelIndex,
		// thisKernelComponent);
		// }
		// return thisBasisSignal;
		throw new Exception("Kernel component has not been initialized yet");
	}

	/**
	 * returns the basis function corresponding to a kernel comp. as a signal array
	 *
	 * @param thisKernelIndex
	 * @param thisKernelComponent
	 * @return
	 * @throws Exception
	 */
	public Signal getDifferentialKernelComp(int thisKernelIndex) throws Exception {
		if (componentDifferentialsCache != null) {
			return componentDifferentialsCache.get(thisKernelIndex);
		}
		throw new Exception("Kernel component has not been initialized yet");
	}

	/**
	 * This method returns the kernel as a signal array
	 *
	 * @param thisKernelIndex
	 * @return
	 * @throws Exception
	 */
	public Signal getKernel(int thisKernelIndex) throws Exception {
		return kernelsCache.get(thisKernelIndex);
	}

	/**
	 * This method returns the kernel as a signal array
	 *
	 * @param thisKernelIndex
	 * @return
	 * @throws Exception
	 */
	public Signal getInvertedKernel(int thisKernelIndex) throws Exception {
		return invertedKernelsCache.get(thisKernelIndex);
	}

	/**
	 * This method loads all the kernel differential signals
	 *
	 * @throws Exception
	 */
	private void loadDifferentialKernels() throws Exception {
		kernelDifferentialsCache = new ArrayList<>();
		for (int i = 0; i < numberOfKernels; i++) {
			kernelDifferentialsCache.add(SignalUtils.calculateSignalDifferential(getKernel(i)));
		}
	}

	private void loadDifferentialKernelComponents() throws Exception {
		componentDifferentialsCache = new ArrayList<>();
		for (int i = 0; i < numberOfKernels; i++) {
			componentDifferentialsCache.add(SignalUtils.calculateSignalDifferential(getBasisForKernelComp(i, 0)));
		}
	}

	/**
	 * This method returns the differential signal of the kernel signal
	 *
	 * @param kernelIndex
	 * @return
	 * @throws Exception
	 */
	public Signal getDifferentialKernel(int kernelIndex) throws Exception {
		// return
		// SignalUtils.calculateSignalDifferential(getKernel(kernelIndex));
		if (kernelDifferentialsCache != null) {
			return kernelDifferentialsCache.get(kernelIndex);
		}
		throw new Exception("The kernel differentials cache has not been initialized yet");
	}

	/**
	 * This function returns the multiplication of two time shifted kernel component
	 * values
	 *
	 * @param kernelIndex1 index of the first kernel
	 * @param kernelIndex2 index of the second kernel
	 * @param compIndex1   index of the first component
	 * @param compIndex2   index of the second component
	 * @param timeShift1   shift of the first signal
	 * @param timeShift2   shift of the second signal
	 * @return the integration of the multiplied signal
	 * @throws Exception
	 */
	public double getCompByCompMult(int kernelIndex1, int kernelIndex2, int compIndex1, int compIndex2,
			double timeShift1, double timeShift2) throws Exception {
		int i, j, ci, cj;
		double ti, tj;
		if (kernelIndex1 < kernelIndex2) {
			i = kernelIndex1;
			j = kernelIndex2;
			ci = compIndex1;
			cj = compIndex2;
			ti = timeShift1;
			tj = timeShift2;
		} else {
			j = kernelIndex1;
			i = kernelIndex2;
			cj = compIndex1;
			ci = compIndex2;
			tj = timeShift1;
			ti = timeShift2;
		}

		int deli = (int) getBasisForKernelComp(i, ci).getEndTime();
		int delj = (int) getBasisForKernelComp(j, cj).getEndTime();
		// refer to the underlying mathematics to see how below equation is
		// valid
		// if(i==0 && j==0) {
		// System.out.println("faulty condition");
		// }
		/******************/
		/** accuracy check **/
		/*****************/
		// return
		// KernelCalculator.bSplineByBspline((int)getBasisForKernelComp(i,
		// 0).end/3, (int)getBasisForKernelComp(j, 0).end/3, -deli+delj-ti+tj);
		/******************/
		/*** revert back ****/
		/******************/
		return componentByComponentSignalConvolutionCache[i][j].getSignalValue(deli - delj - ti + tj);
	}

	/**
	 * TODO: start from here from tomorrow This method populates the
	 * ComponentByComponentSignalConvolution cache
	 *
	 * @throws Exception
	 */
	private void loadComponentByComponentSignalConvolutionCache() throws Exception {
		componentByComponentSignalConvolutionCache = new Signal[numberOfKernels][numberOfKernels];
		for (int i = 0; i < numberOfKernels; i++) {
			for (int j = i; j < numberOfKernels; j++) {
				Signal sigi = getBasisForKernelComp(i, 0);
				Signal sigj = getBasisForKernelComp(j, 0);
				int e1 = (int) sigi.getEndTime();
				int e2 = (int) sigj.getEndTime();
				double[] innerproducts = new double[e1 + e2];
				for (int delta = -e2; delta < e1; delta++) {
					/**********************/
					/***** new version *****/
					/**********************/
					innerproducts[delta + e2] = kernelCalc.bSplineByBspline(i, j, delta);
					/**********************/
					/***** old version *****/
					/**********************/
					/*
					 * Signal shiftedSignal = SignalUtils.shiftSignal(sigj, delta); Signal
					 * multipliedSignal = SignalUtils.multiplyTwoSignals(sigi, shiftedSignal);
					 * innerproducts[delta + e2] =
					 * SignalUtils.calculateSignalIntegral(multipliedSignal);
					 */
				}
				componentByComponentSignalConvolutionCache[i][j] = new Signal(innerproducts, -e2, e1, -e2);
				// componentByComponentSignalConvolutionCache[i][j].DrawSignal("compbycomp"+i+":"+j);
			}
		}
	}

	/**
	 * TODO: check if all the calculations here are correct This method populates
	 * the ComponentByComponentSignalConvolution cache
	 *
	 * @throws Exception
	 */
	private void loadComponentByDifferentialComponentSignalConvolutionCache() throws Exception {
		componentByDiffComponentSignalConvolutionCache = new Signal[numberOfKernels][numberOfKernels];
		for (int i = 0; i < numberOfKernels; i++) {
			for (int j = 0; j < numberOfKernels; j++) {
				Signal sigi = getBasisForKernelComp(i, 0);
				Signal sigj = getDifferentialKernelComp(j);
				int e1 = (int) sigi.getEndTime();
				int e2 = (int) sigj.getEndTime();
				double[] innerproducts = new double[e1 + e2];
				for (int delta = -e2; delta < e1; delta++) {
					/**********************/
					/***** new version *****/
					/**********************/
					innerproducts[delta + e2] = kernelCalc.bsplineByDbSpline(i, j, delta);
					/**********************/
					/***** old version *****/
					/**********************/
					/*
					 * Signal shiftedSignal = SignalUtils.shiftSignal(sigj, delta); Signal
					 * multipliedSignal = SignalUtils.multiplyTwoSignals(sigi, shiftedSignal);
					 * innerproducts[delta + e2] =
					 * SignalUtils.calculateSignalIntegral(multipliedSignal);
					 */
				}
				componentByDiffComponentSignalConvolutionCache[i][j] = new Signal(innerproducts, -e2, e1, -e2);
			}
		}
	}

	/**
	 * returns the threshold value of a kernel at a particular time instance
	 *
	 * @param kernelIndex
	 * @param time
	 * @return
	 */
	public double getThreshold(int kernelIndex, int time) {
		return thresholdValues[kernelIndex][time];
	}

	/**
	 * returns the inner product of a comp. (1st arg.) with another comp.
	 * differential signal (2nd arg.)
	 *
	 * @param kernelIndex
	 * @param differentialKernelIndex
	 * @param componentindex
	 * @param c
	 * @param ti
	 * @param tj
	 * @return
	 * @throws Exception
	 */
	public double getCompByCompDiff(int kernelIndex, int differentialKernelIndex, int componentindex, int c, double ti,
			double tj) throws Exception {
		// TODO: now we need to fill this out
		int deli = (int) getBasisForKernelComp(kernelIndex, componentindex).getEndTime();
		int delj = (int) getBasisForKernelComp(differentialKernelIndex, c).getEndTime();
		// refer to the underlying mathematics to see how below equation is
		// valid

		return -componentByDiffComponentSignalConvolutionCache[kernelIndex][differentialKernelIndex]
				.getSignalValue(deli - delj - ti + tj);
	}

	/**
	 * This method increments the value of a particular kernel coefficient after a
	 * gradient step has been taken
	 *
	 * @param kernelIndex
	 * @param coefficientIndex
	 * @param increment
	 */
	public void incrementKernelCoefficient(int kernelIndex, int coefficientIndex, double increment) {
		// TODO Auto-generated method stub
		kernelCoefficients[kernelIndex][coefficientIndex] += increment;
	}

	/**
	 * This method updates each kernel coefficient by the amount given in the 2-D
	 * array
	 *
	 * @param updateOnKernelCoeffs
	 */
	public void incrementKernelCoefficient(double[][] updateOnKernelCoeffs) {
		// TODO Auto-generated method stub
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			for (int j = 0; j < ConfigurationParameters.numberofKernelComponents; j++) {
				kernelCoefficients[i][j] += updateOnKernelCoeffs[i][j];
			}
		}
	}

	/**
	 * Method to roll back previous update step
	 *
	 * @param updateOnKernelCoeffs
	 */
	public void rollBackUpdate(double[][] updateOnKernelCoeffs) {
		// TODO Auto-generated method stub
		// TODO Auto-generated method stub
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			for (int j = 0; j < ConfigurationParameters.numberofKernelComponents; j++) {
				kernelCoefficients[i][j] -= updateOnKernelCoeffs[i][j];
			}
		}
	}

	/**
	 * This method normalizes all the kernels to give them a unit l2 norm
	 * 
	 * @param shouldUpdateCache
	 * @throws Exception
	 */
	public void normalizeAllKernels(boolean shouldUpdateCache) throws Exception {
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			normalizeKernel(i, false);
		}
		if(shouldUpdateCache) {
			updateCache();
		}
	}

	/**
	 * This method normalizes a kernel
	 * 
	 * @param index of the kernel
	 * @throws Exception
	 */
	public void normalizeKernel(int index, boolean shouldUpdateCache) throws Exception {
		double bSquares = 0;
		double bCrossSquares = 0;
		for (int i = 0; i < ConfigurationParameters.numberofKernelComponents; i++) {
			bSquares += kernelCoefficients[index][i] * kernelCoefficients[index][i];
			if (i != 0) {
				bCrossSquares += kernelCoefficients[index][i] * kernelCoefficients[index][i - 1];
			}
		}
		double kernelNorm = Math.sqrt(
				((A_INTEGRAL + B_INTEGRAL) * bSquares + (C_INTEGRAL) * bCrossSquares) * kernelBsplineLengths[index]);
		for (int i = 0; i < ConfigurationParameters.numberofKernelComponents; i++) {
			kernelCoefficients[index][i] /= kernelNorm;
		}
		if (shouldUpdateCache) {
			updateCache();
		}
	}

	/**
	 * Returns the normalized gradient for the kernel norm constraint
	 * 
	 * @param kernelIndex
	 * @return
	 */
	public double[] getGradGNormalized(int kernelIndex) {
		double[] gradG = new double[ConfigurationParameters.numberofKernelComponents];
		for (int i = 0; i < ConfigurationParameters.numberofKernelComponents; i++) {
			gradG[i] = (2 * (A_INTEGRAL + B_INTEGRAL) * kernelCoefficients[kernelIndex][i]);
			if (i > 0) {
				gradG[i] += (2 * (C_INTEGRAL) * kernelCoefficients[kernelIndex][i - 1]);
			}
			if (i + 1 < ConfigurationParameters.numberofKernelComponents) {
				gradG[i] += (2 * (C_INTEGRAL) * kernelCoefficients[kernelIndex][i + 1]);
			}
		}
		gradG = Utilities.normalizeVector(gradG);
		return gradG;// fix this has a bug due to wrong kernel norm condition
	}

	/**
	 * Projects the kernel gradient on the surface of the kernel constraint
	 * 
	 * @param kernelIndex
	 * @param deltas
	 * @return returns the delta step projected onto the constraint surface
	 */
	public double[] projectAlongKernelConstraint(int kernelIndex, double[] deltas) {
		double[] kernelGGradNormalized = getGradGNormalized(kernelIndex);
		double innerProductWthGradient = Utilities.vectorInnerProduct(deltas, kernelGGradNormalized);
		kernelGGradNormalized = Utilities.scaleVector(-innerProductWthGradient, kernelGGradNormalized);
		kernelGGradNormalized = Utilities.addTwoVectors(kernelGGradNormalized, deltas);
		return kernelGGradNormalized;
	}

	/**
	 * Given a gradient vector this method gives its projection along the kernel
	 * constraint surface
	 * 
	 * @param gradient
	 * @return
	 */
	public double[][] getGradientAlongKernelConstraints(double[][] gradient) {
		if (gradient.length != numberOfKernels) {
			throw new IllegalArgumentException("Not sufficient arguments provided in gradient update");
		}
		double[][] projectedGradient = new double[numberOfKernels][];
		for (int i = 0; i < numberOfKernels; i++) {
			projectedGradient[i] = projectAlongKernelConstraint(i, gradient[i]);
		}
		return projectedGradient;
	}
}
