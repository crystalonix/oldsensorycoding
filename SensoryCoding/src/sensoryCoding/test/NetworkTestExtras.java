package sensoryCoding.test;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import sensoryCoding.network.ConfigurationParameters;
import sensoryCoding.network.DebugConfigurationParameters;
import sensoryCoding.network.KernelManager;
import sensoryCoding.network.Network;
import sensoryCoding.network.Signal;
import sensoryCoding.network.SignalUtils;

public class NetworkTestExtras {
	/**
	 * try the following 2 things: 1. sort based on coeff*threshold; 2. stepwise
	 * remove spikes 
	 * @throws Exception
	 */

	/**
	 * tests the reconstruction of signal at reduced spike rate
	 * @throws Exception
	 */
	@Test
	public void testReconstructionWithSpikeCompression() throws Exception {
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 10;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.lengthOfBasicBSpline = 30;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = false;
		ConfigurationParameters.AHP_REFRACTORY_PERIOD = 5.0;
		ConfigurationParameters.AHP_HIGH_VALUE = 1000.0;
		ConfigurationParameters.AHP_SLOPE = ConfigurationParameters.AHP_HIGH_VALUE
				/ ConfigurationParameters.AHP_REFRACTORY_PERIOD;
		KernelManager krMgr = new KernelManager(false);

		DebugConfigurationParameters.DEBUG_MODE_ON = false;
		/**
		 * construct a test signal
		 */
		List<Signal> signalComps = new ArrayList<>();
		// initial shift of a component
		double shift = 0;
		// for leaving gap between components
		double gap = 100;
		/**
		 * Change this later to randomize
		 */
		int[] selectedKers = { 3, 8, 6 };
		// Construct the test signal
		for (int i = 0; i < 3; i++) {
			/*
			 * int index = r.nextInt(10); selectedKers[i] = index;
			 */
			int index = selectedKers[i];
			shift = shift + gap + krMgr.kernelLengths[index];
			Signal comp = SignalUtils.shiftSignal(krMgr.getInvertedKernel(index), shift);
			signalComps.add(comp);
		}

		// run the network once
		Signal testSignal = SignalUtils.addNSignals(signalComps);
		testSignal = SignalUtils.scalarMultiply(testSignal, 30);
		// testSignal.DrawSignal("Experimental Signal");

		Signal blank = null;
		Network net = new Network(blank);
		net.init(testSignal);
		net.calculateSpikeTimesAndReconstructSignal();
		System.out.println("Actual signal norm square:" + net.signalNormSquare);
		for (int i = net.spikeTimings.size()-1; i > 0; i--) {
			net.reconstructSignalWithNSpikes(i);
			double thisError = net.calculateErrorFastWithNSpikes();
			//Signal reconstructionWithNSpikes = net.getReconstructedSignalWithNSpikes();
			//Signal[] signals = { testSignal, reconstructionWithNSpikes };
			// Utilities.displaySetOfSignals(signals, "actual and reconstruction");
			System.out.println("error with" + i + " spikes is:" + thisError);
		}
		System.out.println("End of reconstruction");
	}
}
