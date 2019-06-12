package sensoryCoding.test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.Test;

import sensoryCoding.network.ConfigurationParameters;
import sensoryCoding.network.DebugConfigurationParameters;
import sensoryCoding.network.KernelManager;
import sensoryCoding.network.Network;
import sensoryCoding.network.Signal;
import sensoryCoding.network.SignalUtils;
import sensoryCoding.network.Utilities;

public class NetworkRunnerTest {
	public static double LEARNING_RATE = 0.01;
	public static double INITIAL_LEARNING_RATE = 0.01;
	public static boolean IS_CONSTRAINED_GRADIENT_DESCENT = false;
	public static boolean DYNAMIC_LEARNING_RATE = false;
	public static String OUTPUT_FOLDER_PATH = null;
	private static final String KERNEL_FILE_PATH = "kernelSignal";
	private static final String ERROR_FILE_PATH = "movingErrorAvg.txt";
	private static boolean SHOULD_NORMALIZE_KERNEL_UPDATE = false;
	private static boolean SHOULD_NORMALIZE_ERRORGRAD = false;
	private static boolean ERROR_DEPENDENT_LEARNING = false;

	@Test
	public void testUtilSimple() throws Exception {
		double[] sampleErr = { 0.2, -0.1, 0.5, 0.6 };
		SignalUtils.storeSignal("ExmapleSignal.txt", new Signal(sampleErr));
		sampleErr[0] = 0.3;
		SignalUtils.storeSignal("ExmapleSignal.txt", new Signal(sampleErr));
	}

	/*
	 * @Test public void testUpdateOnKernelNormalization() throws Exception {
	 * ConfigurationParameters.numberofKernelComponents = 6;
	 * ConfigurationParameters.numberOfKernels = 1;
	 * ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = true;
	 * ConfigurationParameters.lengthOfBasicBSpline = 60;
	 * 
	 * KernelManager krMgr = new KernelManager(true); int rows = 1; int cols = 6;
	 * double maxAbsValue = 100; double[][] testSignalRandomCoeffs =
	 * getRandomValuesMatrix(rows, cols, maxAbsValue); krMgr.kernelCoefficients =
	 * testSignalRandomCoeffs; krMgr.normalizeAllKernels(true);
	 * System.out.println(SignalUtils.calculateSquaredNormAccurate(krMgr.
	 * getInvertedKernel(0)));
	 * 
	 * testSignalRandomCoeffs = getRandomValuesMatrix(rows, cols, maxAbsValue);
	 * testSignalRandomCoeffs =
	 * krMgr.getGradientAlongKernelConstraints(testSignalRandomCoeffs);
	 * testSignalRandomCoeffs = Utilities.scale(testSignalRandomCoeffs, 0.001);
	 * krMgr.incrementKernelCoefficient(testSignalRandomCoeffs);
	 * krMgr.updateCache();
	 * System.out.println(SignalUtils.calculateSquaredNormAccurate(krMgr.
	 * getInvertedKernel(0))); }
	 */
	/**
	 * This method runs stochastic gradient descent of signal snippets
	 * 
	 * @throws Exception
	 */
	@Test
	public void testLearningOnKernelSignalsStochastic() throws Exception {
		/**
		 * ======================== Configure here ======================
		 */
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 1;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = true;
		ConfigurationParameters.lengthOfBasicBSpline = 60;
		ConfigurationParameters.AHP_REFRACTORY_PERIOD = 30.0;
		ConfigurationParameters.AHP_HIGH_VALUE = 1000.0;
		ConfigurationParameters.AHP_SLOPE = ConfigurationParameters.AHP_HIGH_VALUE
				/ ConfigurationParameters.AHP_REFRACTORY_PERIOD;
		ConfigurationParameters.initialThresHoldValue = 1;

		INITIAL_LEARNING_RATE = 0.01;
		LEARNING_RATE = 0.01;
		DYNAMIC_LEARNING_RATE = false;
		ERROR_DEPENDENT_LEARNING = true;
		SHOULD_NORMALIZE_ERRORGRAD = true;
		boolean normalizeInBetween = false;
		IS_CONSTRAINED_GRADIENT_DESCENT = false;
		SHOULD_NORMALIZE_KERNEL_UPDATE = true;
		int samplingInterval = 200;
		DebugConfigurationParameters.DEBUG_MODE_ON = false;
		KernelManager krMgr = new KernelManager(true);

		int numberOfIterations = 500000;
		double[] movingAverage = new double[numberOfIterations];
		double[][] testSignalRandomCoeffs = { { 0.1, 0.3, -2.2, 4.1, 1.2, -2.2 } };
		krMgr.kernelCoefficients = testSignalRandomCoeffs;
		krMgr.normalizeAllKernels(true);
		krMgr.getInvertedKernel(0).DrawSignal("Signal Kernel");
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		/**
		 * Set the stage here
		 */
		Signal blank = null;
		Network net = new Network(blank);
		Signal[] initialKernels = new Signal[ConfigurationParameters.numberOfKernels];
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			initialKernels[i] = net.kernelMgr.getInvertedKernel(i);
			initialKernels[i].DrawSignal("Initial Kernels");
		}
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		/**
		 * Do the computation
		 */
		double initialErrorRate = 1.0;
		double updatedErrorRate = 1.0;
		List<Double> errorsRegular = new ArrayList<>();
		double[] errorsFast = new double[numberOfIterations];
		for (int itr = 0; itr < numberOfIterations; itr++) {
			/**
			 * ======================= Set the test signal===============================
			 */
			Signal testSignal = getRandomSignalFromKernels(new Signal[] { krMgr.getInvertedKernel(0) }, 1);
			if (itr % 100 == 0) {
				// testSignal.DrawSignal("Signal for test");
			}

			/**
			 * ======================= Run the network here===============================
			 */
			net.init(testSignal);

			net.calculateSpikeTimesAndReconstructSignal();
			errorsRegular.add(net.calculateError());
			double error = net.calculateErrorFast();
			errorsFast[itr] = error;
			net.calculateErrorGradient();
			if (itr == 1) {
				net.kernelMgr.getInvertedKernel(0).DrawSignal("This kernel");
			}
			double[][] errorGrads = net.errorGradients;
			double errorRate = error / net.signalNormSquare;
			movingAverage[itr] = itr == 0 ? errorRate : (movingAverage[itr - 1] * itr + errorRate) / (itr + 1);

			updatedErrorRate = errorRate;
			if (DYNAMIC_LEARNING_RATE) {
				if (updatedErrorRate > initialErrorRate) {
					// drop the learning rate
					LEARNING_RATE /= 2;
				}
				if (updatedErrorRate <= initialErrorRate && LEARNING_RATE < INITIAL_LEARNING_RATE) {
					LEARNING_RATE *= 1.414;
				}
			}
			initialErrorRate = updatedErrorRate;
			updateNetworkState(net, errorGrads, error, errorRate);
			// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			/**
			 * Display at the end
			 */
			System.out.println("Error computed in two different ways:" + errorsRegular.get(itr) + "," + errorsFast[itr]
					+ " at step#" + itr + "number of spikes:" + net.spikeTimings.size() + " error rate:" + errorRate
					+ ", moving average of errorate" + movingAverage[itr]);
			if (itr % samplingInterval == 0 && itr != 0) {

				if (normalizeInBetween) {
					net.kernelMgr.normalizeAllKernels(true);
				}
				/**
				 * Save the state here
				 */
				Signal[] finalKernels = new Signal[ConfigurationParameters.numberOfKernels];
				for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
					finalKernels[i] = net.kernelMgr.getInvertedKernel(i);
					SignalUtils.storeSignal(KERNEL_FILE_PATH + i + ".txt", finalKernels[i]);
					// finalKernels[i].DrawSignal("final kernel");
				}
				Signal dummyMovingAvgSignal = new Signal(movingAverage);
				SignalUtils.storeSignal(ERROR_FILE_PATH, dummyMovingAvgSignal);
			}

		}

	}
	/*
	 * @Test public void testSignalConv() throws Exception {
	 *//**
		 * ======================== Configure here ======================
		 */
	/*
	 * ConfigurationParameters.numberofKernelComponents = 6;
	 * ConfigurationParameters.numberOfKernels = 10;
	 * ConfigurationParameters.lengthOfBasicBSpline = 60; KernelManager krMgr = new
	 * KernelManager(true); int compIndex = 2;
	 *//**
		 * Set the test components here
		 */
	/*
	 * double[][] testSignalRandomCoeffs =
	 * getRandomValuesMatrix(ConfigurationParameters.numberOfKernels,
	 * ConfigurationParameters.numberofKernelComponents, 1);
	 * krMgr.kernelCoefficients = testSignalRandomCoeffs;
	 * krMgr.normalizeAllKernels(true); // Store the component Snippets Signal[]
	 * compSignals = new Signal[ConfigurationParameters.numberOfKernels]; for (int i
	 * = 0; i < ConfigurationParameters.numberOfKernels; i++) { compSignals[i] =
	 * krMgr.getInvertedKernel(i); } //
	 * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	 * 
	 *//**
		 * Set the stage here before actual run
		 */
	/*
	 * Signal blank = null; Network net = new Network(blank); Signal[]
	 * initialKernels = new Signal[ConfigurationParameters.numberOfKernels]; for
	 * (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
	 * initialKernels[i] = net.kernelMgr.getInvertedKernel(i);
	 * //initialKernels[i].DrawSignal("Initial Kernel-" + i);
	 * //SignalUtils.storeSignal(OUTPUT_FOLDER_PATH + KERNEL_FILE_PATH + i + ".txt",
	 * initialKernels[i]); } //
	 * ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	 *//**
		 * Do the computation
		 */
	/*
	 * Signal testSignal = SignalUtils.shiftSignal(krMgr.getBasisForKernelComp(0,
	 * 0), 100);//getRandomSignalFromKernels(compSignals, 4);
	 * testSignal.DrawSignal("Test Signal here", 0, 400);
	 * 
	 *//**
		 * ======================= Run the network here===============================
		 *//*
			 * net.init(testSignal);
			 * net.signalKernelComponentConvolutions[0].DrawSignal("Comp Conv", -100, 400);
			 * Signal basisComp = net.kernelMgr.getBasisForKernelComp(0, compIndex);
			 * basisComp.DrawSignal("comp :"+compIndex); System.out.println("Checked");
			 * 
			 * }
			 */

	/**
	 * This method runs stochastic gradient descent of signal snippets with multiple
	 * components
	 * 
	 * @throws Exception
	 */
	@Test
	public void testLearningOnKernelSignalsStochasticMultComp() throws Exception {
		/**
		 * ======================== Configure here ======================
		 */
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 10;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = true;
		ConfigurationParameters.lengthOfBasicBSpline = 60;
		ConfigurationParameters.AHP_HIGH_VALUE = 1000.0;
		// ConfigurationParameters.initialThresHoldValue = 1;

		INITIAL_LEARNING_RATE = 0.2;
		LEARNING_RATE = 0.2;
		DYNAMIC_LEARNING_RATE = false;
		// ++++++++++Important+++++++++++++++++++++++//
		ConfigurationParameters.AHP_REFRACTORY_PERIOD = 60.0;
		ConfigurationParameters.AHP_SLOPE = ConfigurationParameters.AHP_HIGH_VALUE
				/ ConfigurationParameters.AHP_REFRACTORY_PERIOD;

		ERROR_DEPENDENT_LEARNING = true;
		SHOULD_NORMALIZE_ERRORGRAD = true;
		SHOULD_NORMALIZE_KERNEL_UPDATE = true;
		OUTPUT_FOLDER_PATH = "testResultsWith10CompAttempt3/";
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 0.2;
		ConfigurationParameters.initialThresHoldValue = 10;
		ConfigurationParameters.lengthOfComponentSignals = 15000;
		int numberOfCompSnippetsUsed = 5;

		// +++++++++++++++++++++++++++++++++++++++++//
		boolean normalizeInBetween = false;
		IS_CONSTRAINED_GRADIENT_DESCENT = false;
		int numberOfIterations = 500000;
		int samplingInterval = 30;
		DebugConfigurationParameters.DEBUG_MODE_ON = false;
		KernelManager krMgr = new KernelManager(true);

		/**
		 * Set the test components here
		 */
		double[] movingAverage = new double[numberOfIterations];
		double[][] testSignalRandomCoeffs = getRandomValuesMatrix(ConfigurationParameters.numberOfKernels,
				ConfigurationParameters.numberofKernelComponents, 1);
		krMgr.kernelCoefficients = testSignalRandomCoeffs;
		krMgr.normalizeAllKernels(true);
		// Store the component Snippets
		Signal[] compSignals = new Signal[ConfigurationParameters.numberOfKernels];
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			compSignals[i] = krMgr.getInvertedKernel(i);
			//compSignals[i].DrawSignal("Signal Kernel");
			SignalUtils.storeSignal(OUTPUT_FOLDER_PATH + "SignalComp" + i + ".txt", krMgr.getInvertedKernel(i), false);
		}
		Utilities.writeMatrixToFileWithoutAnyHeader(OUTPUT_FOLDER_PATH + "signalCompCoeffs.txt",
				krMgr.kernelCoefficients);
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		/**
		 * Set the stage here before actual run
		 */
		Signal blank = null;
		Network net = new Network(blank);
		Signal[] initialKernels = new Signal[ConfigurationParameters.numberOfKernels];
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			initialKernels[i] = net.kernelMgr.getInvertedKernel(i);
			//initialKernels[i].DrawSignal("Initial Kernel-" + i);
			SignalUtils.storeSignal(OUTPUT_FOLDER_PATH + KERNEL_FILE_PATH + "-initial" + i + ".txt", initialKernels[i]);
			Utilities.writeMatrixToFileWithoutAnyHeader(OUTPUT_FOLDER_PATH + "initialKernelCoeffsCoeffs.txt",
					net.kernelMgr.kernelCoefficients);
		}
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		/**
		 * Do the computation
		 */
		double initialErrorRate = 1.0;
		double updatedErrorRate = 1.0;
		List<Double> errorsRegular = new ArrayList<>();
		double[] errorsFast = new double[numberOfIterations];
		for (int itr = 0; itr < numberOfIterations; itr++) {
			/**
			 * ======================= Set the test signal===============================
			 */
			Signal testSignal = getRandomSignalFromKernels(compSignals, numberOfCompSnippetsUsed);
			if (itr % 100 == 0) {
				// testSignal.DrawSignal("Signal for test");
			}

			/**
			 * ======================= Run the network here===============================
			 */
			net.init(testSignal);

			net.calculateSpikeTimesAndReconstructSignal();
			errorsRegular.add(net.calculateError());
			double error = net.calculateErrorFast();
			errorsFast[itr] = error;
			net.calculateErrorGradient();
			if (itr == 1) {
				net.kernelMgr.getInvertedKernel(0).DrawSignal("This kernel");
			}
			double[][] errorGrads = net.errorGradients;
			double errorRate = error / net.signalNormSquare;
			movingAverage[itr] = itr == 0 ? errorRate : (movingAverage[itr - 1] * itr + errorRate) / (itr + 1);

			updatedErrorRate = errorRate;
			if (DYNAMIC_LEARNING_RATE) {
				if (updatedErrorRate > initialErrorRate) {
					// drop the learning rate
					LEARNING_RATE /= 2;
				}
				if (updatedErrorRate <= initialErrorRate && LEARNING_RATE < INITIAL_LEARNING_RATE) {
					LEARNING_RATE *= 1.414;
				}
			}
			initialErrorRate = updatedErrorRate;
			updateNetworkState(net, errorGrads, error, errorRate);
			// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			/**
			 * Display at the end
			 */
			System.out.println("Error computed in two different ways:" + errorsRegular.get(itr) + "," + errorsFast[itr]
					+ " at step#" + itr + "number of spikes:" + net.spikeTimings.size() + " error rate:" + errorRate
					+ ", moving average of errorate" + movingAverage[itr]);
			if (itr % samplingInterval == 0 && itr != 0) {

				if (normalizeInBetween) {
					net.kernelMgr.normalizeAllKernels(true);
				}
				/**
				 * Save the state here
				 */
				Signal[] finalKernels = new Signal[ConfigurationParameters.numberOfKernels];
				for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
					finalKernels[i] = net.kernelMgr.getInvertedKernel(i);
					SignalUtils.storeSignal(OUTPUT_FOLDER_PATH + KERNEL_FILE_PATH + "-final-" + i + ".txt",
							finalKernels[i], false);
					// finalKernels[i].DrawSignal("final kernel");
				}
				Utilities.writeMatrixToFileWithoutAnyHeader(OUTPUT_FOLDER_PATH + "finalKernelCoeffsCoeffs.txt",
						net.kernelMgr.kernelCoefficients);
				Signal dummyMovingAvgSignal = new Signal(movingAverage);
				SignalUtils.storeSignal(OUTPUT_FOLDER_PATH + ERROR_FILE_PATH, dummyMovingAvgSignal, false);
			}

		}

	}
	
	/**
	 * This method runs stochastic gradient descent of signal snippets with multiple
	 * components
	 * 
	 * @throws Exception
	 */
	@Test
	public void trainAndtestOnKernelSignalsStochasticMultComp() throws Exception {
		/**
		 * ======================== Configure here ======================
		 */
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 10;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = true;
		ConfigurationParameters.lengthOfBasicBSpline = 60;
		ConfigurationParameters.AHP_HIGH_VALUE = 1000.0;
		// ConfigurationParameters.initialThresHoldValue = 1;

		INITIAL_LEARNING_RATE = 0.2;
		LEARNING_RATE = 0.2;
		DYNAMIC_LEARNING_RATE = false;
		// ++++++++++Important+++++++++++++++++++++++//
		ConfigurationParameters.AHP_REFRACTORY_PERIOD = 50.0;
		ConfigurationParameters.AHP_SLOPE = ConfigurationParameters.AHP_HIGH_VALUE
				/ ConfigurationParameters.AHP_REFRACTORY_PERIOD;

		ERROR_DEPENDENT_LEARNING = true;
		SHOULD_NORMALIZE_ERRORGRAD = true;
		SHOULD_NORMALIZE_KERNEL_UPDATE = true;
		OUTPUT_FOLDER_PATH = "testResultsWith10CompTestAndTrain/";
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 0.2;
		ConfigurationParameters.initialThresHoldValue = 10;
		ConfigurationParameters.lengthOfComponentSignals = 15000;
		int numberOfCompSnippetsUsed = 5;

		// +++++++++++++++++++++++++++++++++++++++++//
		boolean normalizeInBetween = false;
		IS_CONSTRAINED_GRADIENT_DESCENT = false;
		int numberOfIterations = 500000;
		int samplingInterval = 5000;
		DebugConfigurationParameters.DEBUG_MODE_ON = false;
		KernelManager krMgr = new KernelManager(true);

		/**
		 * Set the test components here
		 */
		double[] movingAverage = new double[numberOfIterations];
		double[][] testSignalRandomCoeffs = getRandomValuesMatrix(ConfigurationParameters.numberOfKernels,
				ConfigurationParameters.numberofKernelComponents, 1);
		krMgr.kernelCoefficients = testSignalRandomCoeffs;
		krMgr.normalizeAllKernels(true);
		// Store the component Snippets
		Signal[] compSignals = new Signal[ConfigurationParameters.numberOfKernels];
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			compSignals[i] = krMgr.getInvertedKernel(i);
			//compSignals[i].DrawSignal("Signal Kernel");
			SignalUtils.storeSignal(OUTPUT_FOLDER_PATH + "SignalComp" + i + ".txt", krMgr.getInvertedKernel(i), false);
		}
		Utilities.writeMatrixToFileWithoutAnyHeader(OUTPUT_FOLDER_PATH + "signalCompCoeffs.txt",
				krMgr.kernelCoefficients);
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		/**
		 * Set the stage here before actual run
		 */
		Signal blank = null;
		Network net = new Network(blank);
		Signal[] initialKernels = new Signal[ConfigurationParameters.numberOfKernels];
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			initialKernels[i] = net.kernelMgr.getInvertedKernel(i);
			//initialKernels[i].DrawSignal("Initial Kernel-" + i);
			SignalUtils.storeSignal(OUTPUT_FOLDER_PATH + KERNEL_FILE_PATH + "-initial" + i + ".txt", initialKernels[i]);
			Utilities.writeMatrixToFileWithoutAnyHeader(OUTPUT_FOLDER_PATH + "initialKernelCoeffsCoeffs.txt",
					net.kernelMgr.kernelCoefficients);
		}
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		/**
		 * Do the computation
		 */
		double initialErrorRate = 1.0;
		double updatedErrorRate = 1.0;
		List<Double> errorsRegular = new ArrayList<>();
		double[] errorsFast = new double[numberOfIterations];
		for (int itr = 0; itr < numberOfIterations; itr++) {
			/**
			 * ======================= Set the test signal===============================
			 */
			Signal testSignal = getRandomSignalFromKernels(compSignals, numberOfCompSnippetsUsed);
			if (itr % 100 == 0) {
				// testSignal.DrawSignal("Signal for test");
			}

			/**
			 * ======================= Run the network here===============================
			 */
			net.init(testSignal);

			net.calculateSpikeTimesAndReconstructSignal();
			errorsRegular.add(net.calculateError());
			double error = net.calculateErrorFast();
			errorsFast[itr] = error;
			net.calculateErrorGradient();
			if (itr == 1) {
				net.kernelMgr.getInvertedKernel(0).DrawSignal("This kernel");
			}
			double[][] errorGrads = net.errorGradients;
			double errorRate = error / net.signalNormSquare;
			movingAverage[itr] = itr == 0 ? errorRate : (movingAverage[itr - 1] * itr + errorRate) / (itr + 1);

			updatedErrorRate = errorRate;
			if (DYNAMIC_LEARNING_RATE) {
				if (updatedErrorRate > initialErrorRate) {
					// drop the learning rate
					LEARNING_RATE /= 2;
				}
				if (updatedErrorRate <= initialErrorRate && LEARNING_RATE < INITIAL_LEARNING_RATE) {
					LEARNING_RATE *= 1.414;
				}
			}
			initialErrorRate = updatedErrorRate;
			updateNetworkState(net, errorGrads, error, errorRate);
			// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			/**
			 * Display at the end
			 */
			System.out.println("Error computed in two different ways:" + errorsRegular.get(itr) + "," + errorsFast[itr]
					+ " at step#" + itr + "number of spikes:" + net.spikeTimings.size() + " error rate:" + errorRate
					+ ", moving average of errorate" + movingAverage[itr]);
			if (itr % samplingInterval == 0 && itr != 0) {

				if (normalizeInBetween) {
					net.kernelMgr.normalizeAllKernels(true);
				}
				/**
				 * Save the state here
				 */
				Signal[] finalKernels = new Signal[ConfigurationParameters.numberOfKernels];
				for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
					finalKernels[i] = net.kernelMgr.getInvertedKernel(i);
					SignalUtils.storeSignal(OUTPUT_FOLDER_PATH + KERNEL_FILE_PATH + "-final-" + i + ".txt",
							finalKernels[i], false);
					// finalKernels[i].DrawSignal("final kernel");
				}
				/*Utilities.writeMatrixToFileWithoutAnyHeader(OUTPUT_FOLDER_PATH + "finalKernelCoeffsCoeffs.txt",
						net.kernelMgr.kernelCoefficients);
*/				Signal dummyMovingAvgSignal = new Signal(movingAverage);
				SignalUtils.storeSignal(OUTPUT_FOLDER_PATH + ERROR_FILE_PATH, dummyMovingAvgSignal, false);
				Utilities.writeMatrixToFileWithoutAnyHeader(OUTPUT_FOLDER_PATH + "kernelCoeffs-"+itr+".txt",
						net.kernelMgr.kernelCoefficients);

			}

		}

	}

	private Signal getRandomSignalFromKernels(Signal[] signals, int numberOfComps) {
		/**
		 * construct a test signal
		 */
		List<Signal> signalComps = new ArrayList<>();
		// initial shift of a component
		double shift = 100;
		// for leaving gap between components
		double gap = 100;
		/**
		 * Change this later to randomize
		 */
		Random r = new Random();
		for (int i = 0; i < numberOfComps; i++) {
			/*
			 * int index = r.nextInt(10); selectedKers[i] = index;
			 */
			int index = r.nextInt(ConfigurationParameters.numberOfKernels);
			shift = shift + signals[index].getSignalSpan() + (r.nextInt(100) - 20);
			Signal comp = SignalUtils.scalarMultiply(SignalUtils.shiftSignal(signals[index], shift),
					r.nextInt(90) + 10);
			signalComps.add(comp);
		}
		return SignalUtils.addNSignals(signalComps);
	}

	/**
	 * This method updates the kernel coeffs and refreshes the cache
	 * 
	 * @param net
	 * @param errorGrads
	 * @param error
	 * @throws Exception
	 */
	private void updateNetworkState(Network net, double[][] errorGrads, double error, double errorRate)
			throws Exception {
		double errorGradsNorm = Math.sqrt(Utilities.calculateMatrixNormSquare(errorGrads));
		double learningRate = ERROR_DEPENDENT_LEARNING ? LEARNING_RATE * errorRate : LEARNING_RATE;
		if (SHOULD_NORMALIZE_ERRORGRAD) {
			learningRate /= errorGradsNorm;
		}
		double[][] updatedErrorGrad = Utilities.scale(errorGrads, -learningRate);
		if (IS_CONSTRAINED_GRADIENT_DESCENT) {
			updatedErrorGrad = net.kernelMgr.getGradientAlongKernelConstraints(updatedErrorGrad);
		}
		net.kernelMgr.incrementKernelCoefficient(updatedErrorGrad);
		if (SHOULD_NORMALIZE_KERNEL_UPDATE) {
			net.kernelMgr.normalizeAllKernels(true);
		} else {
			net.kernelMgr.updateCache();
		}
	}

	/**
	 * This method fills matrix with random values within a certain range and
	 * returns
	 * 
	 * @param rows
	 * @param cols
	 * @param maxAbsValue
	 * @return
	 */
	private double[][] getRandomValuesMatrix(int rows, int cols, double maxAbsValue) {
		double[][] randomValues = new double[rows][cols];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				randomValues[i][j] = (Math.random()-0.5) * maxAbsValue;
			}
		}
		return randomValues;
	}

	/**
	 * This is an overfit case where the same kernel is tested repeatedly
	 */
	@Test
	public void testSingleKernelSignalConvergence() throws Exception {
		/**
		 * ======================== Configure here ======================
		 */
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 1;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = true;
		ConfigurationParameters.lengthOfBasicBSpline = 60;
		ConfigurationParameters.AHP_REFRACTORY_PERIOD = 10.0;
		ConfigurationParameters.AHP_HIGH_VALUE = 1000.0;
		ConfigurationParameters.AHP_SLOPE = ConfigurationParameters.AHP_HIGH_VALUE
				/ ConfigurationParameters.AHP_REFRACTORY_PERIOD;
		ConfigurationParameters.initialThresHoldValue = 1;

		KernelManager krMgr = new KernelManager(true);

		DebugConfigurationParameters.DEBUG_MODE_ON = false;
		double[][] testSignalRandomCoeffs = { { 0.1, 0.3, -2.2, 4.1, 1.2, -2.2 } };
		krMgr.kernelCoefficients = testSignalRandomCoeffs;
		krMgr.normalizeAllKernels(true);
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		/**
		 * ======================= Set the test signal===============================
		 */
		Signal signal = getRandomSignalFromKernels(new Signal[] { krMgr.getInvertedKernel(0) }, 1);
		Signal testSignal = SignalUtils.scalarMultiply(signal, 30);
		testSignal.DrawSignal("Signal for test");

		List<Double> errorsRegular = new ArrayList<>();
		List<Double> errorsFast = new ArrayList<>();

		/**
		 * Set the stage here
		 */
		Signal blank = null;
		Network net = new Network(blank);
		net.init(testSignal);
		Signal[] initialKernels = new Signal[ConfigurationParameters.numberOfKernels];
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			initialKernels[i] = net.kernelMgr.getInvertedKernel(i);
		}
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		/**
		 * Do the computation
		 */
		double initialErrorRate = 1.0;
		double updatedErrorRate = 1.0;
		for (int itr = 0; itr < 100000; itr++) {
			net.calculateSpikeTimesAndReconstructSignal();
			errorsRegular.add(net.calculateError());
			double error = net.calculateErrorFast();
			errorsFast.add(error);
			net.calculateErrorGradient();
			if (itr == 1) {
				net.kernelMgr.getInvertedKernel(0).DrawSignal("This kernel");
			}
			double[][] errorGrads = net.errorGradients;
			double errorRate = error / net.signalNormSquare;
			updatedErrorRate = errorRate;
			if (DYNAMIC_LEARNING_RATE) {
				if (updatedErrorRate > initialErrorRate) {
					// drop the learning rate
					LEARNING_RATE /= 2;
				}
				if (updatedErrorRate <= initialErrorRate && LEARNING_RATE < INITIAL_LEARNING_RATE) {
					LEARNING_RATE *= 1.414;
				}
			}
			initialErrorRate = updatedErrorRate;
			updateNetworkState(net, errorGrads, error, errorRate);
			// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			/**
			 * Display at the end
			 */
			System.out.println("Error computed in two different ways:" + errorsRegular.get(itr) + ","
					+ errorsFast.get(itr) + " at step#" + itr + "number of spikes:" + net.spikeTimings.size()
					+ " error rate:" + errorRate);

		}
		/**
		 * Save the state here
		 */
		Signal[] finalKernels = new Signal[ConfigurationParameters.numberOfKernels];
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			finalKernels[i] = net.kernelMgr.getInvertedKernel(i);
			finalKernels[i].DrawSignal("final kernel");
		}
		// net.kernelMgr.getInvertedKernel(0);
		System.out.println("saving the state here");
	}

}
