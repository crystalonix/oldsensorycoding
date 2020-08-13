package sensoryCoding.test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.List;
import java.util.stream.IntStream;

import org.junit.Test;

// import Jama.Matrix;

//import com.sun.media.sound.WaveFileReader;

import sensoryCoding.network.ConfigurationParameters;
import sensoryCoding.network.KernelManager;
import sensoryCoding.network.Network;
import sensoryCoding.network.Signal;
import sensoryCoding.network.SignalUtils;
import sensoryCoding.network.Utilities;

/**************************************************************************************/
// go to the root /cise/research/compneuro1/anonymous/sensoryOnRealData/SensoryCoding
// find ./src -name "*.java" > sources_list.txt
// javac -cp "./lib/*" @sources_list.txt
/**************************************************************************************/

public class NetworkWithGammatoneTest {

	@Test
	public void testGammaFilterRecons() throws Exception {
		ConfigurationParameters.numberOfKernels = 10;
		ConfigurationParameters.OVERLAP_FACTOR = 2.0 / 3.0;
		ConfigurationParameters.lengthOfBasicBSpline = 30;
		ConfigurationParameters.numberOfKernels = 1;
		ConfigurationParameters.numberofKernelComponents = 10;
		KernelManager kr = new KernelManager(false);
		Signal kernelSignal = kr.getKernel(0);
		kernelSignal.DrawSignal("Gammatone Signal");
		System.out.println("This is the signal we want");
	}

	@Test
	public void testGammaFilterBankRecons() throws Exception {
		ConfigurationParameters.lengthOfComponentSignals = 882000;
		ConfigurationParameters.numberOfKernels = 90;
		ConfigurationParameters.numberofKernelComponents = 100;
		ConfigurationParameters.OVERLAP_FACTOR = 2.0 / 3.0;
		ConfigurationParameters.SAMPLING_FREQUENCY = 44100;
		ConfigurationParameters.initialThresHoldValue = 0.0001;
		ConfigurationParameters.AHP_REFRACTORY_PERIOD = 100.0;
		ConfigurationParameters.numberOfThreads = 40;
		double[] centerFrequencies = { 31.074903570227125, 42.64272211867989, 54.72539400122914, 67.34583399255955,
				80.52797674401583, 94.2968221756443, 108.67848288851539, 123.70023368724425, 139.3905633066296,
				155.77922844050948, 172.8973101753005, 190.77727293524583, 209.45302605116194, 228.9599880694496,
				249.3351539233301, 270.61716509369796, 292.84638289264905, 316.06496500866746, 340.31694545863735,
				365.6483180983111, 392.1071238496081, 419.74354181017424, 448.6099844179877, 478.76119685149365,
				510.2543608537783, 543.1492031776851, 577.5081088575394, 613.3962395223019, 650.8816569745305,
				690.0354522695166, 730.9318805393995, 773.6485018179459, 818.2663281330758, 864.8699771460863,
				913.5478326289608, 964.3922120840992, 1017.4995418243674, 1072.9705398454967, 1130.9104068376619,
				1191.4290256984762, 1254.6411699257976, 1320.6667212855448, 1389.6308971673438, 1461.6644880591762,
				1536.9041055914035, 1615.4924416205856, 1697.578538844434, 1783.3180734611367, 1872.8736504091044,
				1966.4151117470694, 2064.119858759372, 2166.173188397313, 2272.7686446946295, 2384.108385823553,
				2500.4035674875754, 2621.8747433780186, 2748.752283453884, 2881.276810838241, 3019.6996581597227,
				3164.2833442046062, 3315.302071783412, 3473.04224775625, 3637.8030262031357, 3809.8968757693856,
				3989.650172262086, 4177.403817621469, 4373.513886441091, 4578.35230126292, 4792.307537928067,
				5015.785362320801, 5249.209599903136, 5493.022939499375, 5747.687772855005, 6013.687071562171,
				6291.525303014791, 6581.72938713044, 6884.849695653443, 7201.4610959343, 7532.164041165039,
				7877.585709138046, 8238.38119168812, 8615.234737073451, 9008.861047651717, 9420.00663531239,
				9849.451237235739, 10298.009294663583, 10766.53149748628, 11255.906397575207, 11767.062093920447,
				12300.967992769602, 12858.636646105682, 13441.125671950942, 14049.539760138337, 14685.032767354625,
				15348.809905428418, 16042.130027013145, 16766.308012999776, 17522.71726618713, 18312.79231593877,
				19138.031538766543 };
		int[] bsplineLengths = new int[ConfigurationParameters.numberOfKernels];
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			bsplineLengths[i] = (int) Math.ceil(
					(22.0 / centerFrequencies[i]) * (1.0 / (ConfigurationParameters.numberofKernelComponents + 2.0))
							* ConfigurationParameters.SAMPLING_FREQUENCY);
		}
		System.out.println("Check the lengths of the bsplines");
		double sum = 0;
		for (int l : bsplineLengths) {
			System.out.print(l);
			System.out.print(",");
			sum += l * (ConfigurationParameters.numberofKernelComponents + 2.0);
		}
		System.out.println("Total required length:" + sum);
		double[][] kernelCoeffs = new double[ConfigurationParameters.numberOfKernels][];
		// populate the B-spline kernel coefficients
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			System.out.println("kernel number:" + i);
			String file = "/Users/anonymouschattopadhyay/Desktop/research/oldsensorycoding/SensoryCoding/src/sensoryCoding/bsplineCoefficients/solution-"
					+ i + ".values";
			kernelCoeffs[i] = Utilities.getCoefficientsArrayFromFile(file,
					ConfigurationParameters.numberofKernelComponents);
		}

		KernelManager kr = new KernelManager(false, kernelCoeffs, bsplineLengths);
		kr.normalizeAllKernels(true);
		IntStream.range(0, ConfigurationParameters.numberOfKernels).parallel().forEach(i -> {
			Signal sig = null;
			try {
				sig = kr.getKernel(i);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println("Value of the signal norm:" + SignalUtils.calculateSquaredNormAccurate(sig)
					+ " for kernel index:" + i);
		});
		String soundData = "/Users/anonymouschattopadhyay/Desktop/research/audio_text/4341.txt";
		FileReader in = new FileReader(soundData);
		BufferedReader brSound = new BufferedReader(in);
		List<Double> data = Utilities.readListFromFile(brSound);
		double[] dataVals = new double[data.size()];
		for (int i = 0; i < data.size(); i++) {
			dataVals[i] = data.get(i);
		}
		Signal testSignal = new Signal(dataVals);
		testSignal.DrawSignal("Sample Gamma Signal");
		Network net = new Network(testSignal, kr);
		double before = System.currentTimeMillis();
		net.init(testSignal, true, false);
		System.out.println("Time for init:" + (System.currentTimeMillis() - before));
		before = System.currentTimeMillis();
		net.calculateSpikeTimesAndReconstructSignal();
		System.out.println("Time for reconstruction:" + (System.currentTimeMillis() - before));
		Signal reconsSignal = net.getReconstructedSignal();
		reconsSignal.DrawSignal("Reconstructed signal checked");
		System.out.println("Error Rate:" + net.calculateErrorRate());
		int stepSize = 500;
//		for (int i = net.spikeTimings.size() - 1; i > 0; i--) {
//			net.reconstructSignalWithNSpikes(i);
//			double thisError = net.calculateErrorFastWithNSpikes();
//			System.out.println("error with" + i + " spikes is:" + thisError);
//			if (i % stepSize == 0) {
//				Signal lossySignal = net.getReconstructedSignalWithNSpikes();
//				lossySignal.DrawSignal("Lossy Signal after- " + i + " steps");
//			}
//		}

		net.applyFlatCompression(1);
		System.out.println("End of reconstruction");
	}

//	@Test
//	public void testWavRead() throws UnsupportedAudioFileException, IOException {
//		WaveFileReader wfr = new WaveFileReader();
//		java.io.File waveFile = new java.io.File(
//				"/Users/anonymouschattopadhyay/Downloads/freesound-audio-tagging/audio_train/0a0a8d4c.wav");
//		AudioInputStream str = wfr.getAudioInputStream(waveFile);
//		System.out.println(str.getFrameLength());
//		byte[] btr = new byte[2];
//		str.read(btr);
//		System.out.println("The byte read:" + btr[1]);
//	}

	/**
	 * 
	 */
	@Test
	public void testArrayOpsSeq() {
		int length = 1000000;
		double[] sig = new double[length];

		for (int i = 0; i < sig.length; i++) {
			sig[i] = Math.random();
		}
		double total = 0;
		double t1 = System.currentTimeMillis();
		for (int i = 0; i < sig.length; i++) {
			for (int j = i; j < sig.length; j++) {
				total += sig[i] * sig[j];
			}
		}
		double t2 = System.currentTimeMillis();
		System.out.println("Total time taken is:" + (t2 - t1) / 1000 + " and value of the entire array is:" + total);
	}

	/**
	 * 
	 */
	@Test
	public void testArrayOpsPar() {
		int length = 1000000;
		double[] sig = new double[length];

		for (int i = 0; i < sig.length; i++) {
			sig[i] = Math.random();
		}
		double[] total = new double[length];
		double t1 = System.currentTimeMillis();

		for (int i = 0; i < sig.length; i++) {
			for (int j = i; j < sig.length; j++) {
				total[i] += sig[i] * sig[j];
			}
		}
		double t2 = System.currentTimeMillis();
		System.out.println("Total time taken is:" + (t2 - t1) / 1000 + " and value of the entire array is:" + total[0]);
	}

	/**
	 * 
	 */
	@Test
	public void testArrayOpsParWithStreams() {
		int length = 1000000;
		double[] sig = new double[length];

		for (int i = 0; i < sig.length; i++) {
			sig[i] = Math.random();
		}
		double[] total = new double[length];
		double t1 = System.currentTimeMillis();

		IntStream.range(0, length).parallel().forEach(i -> {
			IntStream.range(i, length).parallel().forEach(j -> {
				total[i] += sig[i] * sig[j];
			});
		});
		double t2 = System.currentTimeMillis();
		System.out.println("Total time taken is:" + (t2 - t1) / 1000 + " and value of the entire array is:" + total[0]);
	}

	@Test
	public void testJacobi() {
		double[][] A = { { 2, -2 }, { 2, 2 } };
		double[] b = { -1, 2 };
		double x[] = Utilities.solveByJacobi(A, b);
		System.out.println("Check the values:" + x[0]);
	}

	@Test
	public void testConjugateGradient() {
		ConfigurationParameters.numberOfThreads = 3;
		double[][] A = { { 1, 2, 3 }, { 2, 5, 1 }, { 3, 1, 2 } };
		double[] b = { 6, 8, 6 };
		double x[] = Utilities.applyConjugateGradient(A, b, 1, -1, null);
		System.out.println("Check the values:" + x[0]);
	}

	/**
	 * 
	 */
	@Test
	public void testArrayRead() {
		String file = "/Users/anonymouschattopadhyay/Desktop/research/oldsensorycoding/SensoryCoding/src/sensoryCoding/bsplineCoefficients/solution-"
				+ 0 + ".values";
		System.out.println(ConfigurationParameters.bsplineCoeffsPath + "solution-" + 0 + ".values");
		System.out.println(ConfigurationParameters.bsplineCoeffsPath + "solution-" + 0 + ".values");
		double[] testArr = Utilities.getCoefficientsArrayFromFile(file, 3);
		System.out.println(testArr);
	}

	@Test
	public void streamTest() {
		IntStream.range(0, 10).parallel().forEach(t -> {
			System.out.println("Thread number is:" + t);
		});
	}

//	@Test
//	public void jamaTest() {
//		double[][] m = new double[][] { { 1, 2, }, { 3, 4 } };
//		Matrix M = new Matrix(m);
//		Matrix mInv = M.inverse();
//		for (int i = 0; i < m.length; i++) {
//			for (int j = 0; j < m[0].length; j++) {
//				System.out.println("The entries are:" + mInv.get(i, j));
//			}
//		}
//		System.out.println("determinant value is:" + M.det());
//	}
}
