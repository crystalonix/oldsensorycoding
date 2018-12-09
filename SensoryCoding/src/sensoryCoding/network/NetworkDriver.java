package sensoryCoding.network;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class NetworkDriver {
	public static Network net = null;
	public static List<Integer> spikeCounts = new ArrayList<>();
	public static List<Double> errors = new ArrayList<>();
	public static List<Double> movingErrorAvg = null;
	private static int stepNumber = 0;

	public static void main(String[] args) throws Exception {
		int numberOfSelectedKernels =10;
		List<Integer> selectedKernels = new ArrayList<>();
		int numberOfBuckets = ConfigurationParameters.numberOfKernels/numberOfSelectedKernels;
		for(int i=0; i<ConfigurationParameters.numberOfKernels; i++){
			if(i%numberOfBuckets==0){
				selectedKernels.add(i);
			}
		}

		int mode = 1;
		/*******************************/
		/**** test the data partition ****/
		/*******************************/
		if (mode == 1) {
			int testPartitionNumber = 14;
			Signal blankSignal = null;
			net = new Network(selectedKernels);
			/*if (args.length > 0) {
			*/	String kernelCoeffsFileName = "kernelCoefficients-300000.txt";
				FileReader fr = new FileReader(ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + kernelCoeffsFileName);
				BufferedReader br = new BufferedReader(fr);
				double[][] kernelCoeffs = Utilities.readMatrixToBufferedWriter(br);
				net.kernelMgr.initKernelCoefficients(kernelCoeffs);
				br.close();
			testOnSinglePartition(testPartitionNumber);
		}
		/*******************************/
		/**** train the data partition ***/
		/*******************************/
		else {
			int high = 5;
			int low = 1;
			int numberOfEpochs = 2000;
			net = new Network(selectedKernels);
			if (args.length > 0) {
				String kernelCoeffsFileName = args[0];
				FileReader fr = new FileReader(ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + kernelCoeffsFileName);
				BufferedReader br = new BufferedReader(fr);
				double[][] kernelCoeffs = Utilities.readMatrixToBufferedWriter(br);
				net.kernelMgr.initKernelCoefficients(kernelCoeffs);
				br.close();
			}
			stepNumber = 0;
			movingErrorAvg = new ArrayList<>();
			for (int k = 0; k < numberOfEpochs; k++) {
				for (int i = low; i <= high; i++) {
					// skipping 7 because this partition is not ready yet
					if(i==7){
						continue;
					}
					trainOnSinglePartition(i);
					if(ConfigurationParameters.SHOULD_COLLECT_SPIKE_STATISTICS){
						collectSpikeStatistics(net);
					}
				}
			}
		}
	}

	private static void collectSpikeStatistics(Network net) throws IOException {
		// TODO Auto-generated method stub
		FileWriter fw1 = new FileWriter(ConfigurationParameters.avgReconsCoeffsFileName);
		FileWriter fw2 = new FileWriter(ConfigurationParameters.avgSpikeCountsFileName);
		BufferedWriter bw1 = new BufferedWriter(fw1);
		Utilities.writeListToBufferedWriter(bw1, net.averageCoefficientOfReConstruction);
		BufferedWriter bw2 = new BufferedWriter(fw2);
		Utilities.writeIntegerListToBufferedWriter(bw2, net.totalSpikeCount);
		bw1.close();
		bw2.close();
	}

	private static void testInEfficientlyOnData(int partitionNumber) throws Exception {
		List<Double> errorRates = new ArrayList<>();
		List<Integer> spikeCounts = new ArrayList<>();
		String seedfileName = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "partitionSeed.txt";
		Signal parts = SignalUtils.readSignalFromFile(seedfileName);
		double[] indexes = parts.data;
		Signal blankSignal = null;
		int startIndex = (partitionNumber - 1) * ConfigurationParameters.numberOfFilesForPartition;
		int endIndex = (partitionNumber) * ConfigurationParameters.numberOfFilesForPartition;
		Network net = new Network(blankSignal);
		for (int j = startIndex; j < endIndex; j++) {
			int fileIndex = (int) indexes[j];
			String fileName = ConfigurationParameters.SOUND_DATA_FOLDER_PATH + "output" + fileIndex + ".txt";
			// File file = new File(fileName);
			double allSignalsInThisPartition[][] = DataReader.readSignalPieces(fileName);
			for (int k = 0; k < allSignalsInThisPartition.length; k++) {
				Signal input = new Signal(allSignalsInThisPartition[k]);
				net.init(input);
				net.calculateSpikeTimesAndReconstructSignal();
				net.calculateErrorGradient();
				double errorate = net.calculateErrorRate();
				errorRates.add(errorate);
				spikeCounts.add(net.spikeTimings.size());
				// as per requirement add the spike count and error avg
				System.out.println("error rate at step number:" + stepNumber + " is:" + errorate);
				if (stepNumber % 500 == 0) {
					// Utilities.writeMatrixToFile(null,
					// ConfigurationParameters.kernelCoeffFileName
					// +"-"+stepNumber+".txt",
					// net.kernelMgr.kernelCoefficients);
					FileWriter fw1 = new FileWriter(ConfigurationParameters.errorValuesFileName + ".txt", true);
					FileWriter fw2 = new FileWriter(ConfigurationParameters.spikeCountFileName + ".txt", true);
					BufferedWriter bw1 = new BufferedWriter(fw1);
					Utilities.writeListToBufferedWriter(bw1, errorRates);
					BufferedWriter bw2 = new BufferedWriter(fw2);
					Utilities.writeIntegerListToBufferedWriter(bw2, spikeCounts);
					bw1.close();
					bw2.close();
					errorRates = new ArrayList<>();
					spikeCounts = new ArrayList<>();
				}
				stepNumber++;
			}
		}
		System.out.println("done here");
	}

	private static void testOnSinglePartition(int partitionNumber) throws Exception {
		List<Double> errorRates = new ArrayList<>();
		List<Integer> spikeCounts = new ArrayList<>();
		DataReader dtr = new DataReader(net);
		dtr.loadPartitionIntoMemory(partitionNumber, 0);
		for (int k = 0; k < dtr.allSignalPieces.length; k++) {
			net.initFromFileData(dtr.allSignalPieces[k], dtr.allDifferentialSignalPieces[k], dtr.signalByKernelComp[k]);
			net.calculateSpikeTimesAndReconstructSignal();
			double errorate = net.calculateErrorRate();
			errorRates.add(errorate);
			spikeCounts.add(net.spikeTimings.size());
			// as per requirement add the spike count and error avg
			System.out.println("error rate at step number:" + stepNumber + " is:" + errorate);
			if (stepNumber % 500 == 0) {
				// Utilities.writeMatrixToFile(null,
				// ConfigurationParameters.kernelCoeffFileName
				// +"-"+stepNumber+".txt", net.kernelMgr.kernelCoefficients);
				FileWriter fw1 = new FileWriter(ConfigurationParameters.errorValuesFileName + ".txt", true);
				FileWriter fw2 = new FileWriter(ConfigurationParameters.spikeCountFileName + ".txt", true);
				BufferedWriter bw1 = new BufferedWriter(fw1);
				Utilities.writeListToBufferedWriter(bw1, errorRates);
				BufferedWriter bw2 = new BufferedWriter(fw2);
				Utilities.writeIntegerListToBufferedWriter(bw2, spikeCounts);
				bw1.close();
				bw2.close();
				errorRates = new ArrayList<>();
				spikeCounts = new ArrayList<>();
			}
			stepNumber++;
		}
		System.out.println("done here");
	}
	private static void trainInefficientlyOnData() throws Exception {
		// TODO Auto-generated method stub
		//int stepNumber =0;
		List<Double> errorRates = new ArrayList<>();
		List<Integer> spikeCounts = new ArrayList<>();
		String seedfileName = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "partitionSeed.txt";
		Signal parts = SignalUtils.readSignalFromFile(seedfileName);
		double[] indexes = parts.data;
		Signal blankSignal = null;
		Network net = new Network(blankSignal);
		for (int j = 0; j < 1000; j++) {
			int fileIndex = (int) indexes[j];
			String fileName = ConfigurationParameters.SOUND_DATA_FOLDER_PATH + "output" + fileIndex + ".txt";
			//File file = new File(fileName);
			double allSignalsInThisPartition[][] = DataReader.readSignalPieces(fileName);
			for (int k = 0; k < allSignalsInThisPartition.length; k++) {
				Signal input = new Signal(allSignalsInThisPartition[k]);
				net.init(input);
				net.calculateSpikeTimesAndReconstructSignal();
				net.calculateErrorGradient();
				double errorate = net.updateKernelCoefficients();
				errorRates.add(errorate);
				spikeCounts.add(net.spikeTimings.size());
				// as per requirement add the spike count and error avg
				System.out.println("error rate at step number:"+stepNumber+ " is:" + errorate);
				if(stepNumber%500==0){
					Utilities.writeMatrixToFile(null, ConfigurationParameters.kernelCoeffFileName +"-"+stepNumber+".txt", net.kernelMgr.kernelCoefficients);
					FileWriter fw1 = new FileWriter(ConfigurationParameters.errorValuesFileName +".txt", true);
					FileWriter fw2 = new FileWriter(ConfigurationParameters.spikeCountFileName +".txt", true);
					BufferedWriter bw1 = new BufferedWriter(fw1);
					Utilities.writeListToBufferedWriter(bw1, errorRates);
					BufferedWriter bw2 = new BufferedWriter(fw2);
					Utilities.writeIntegerListToBufferedWriter(bw2, spikeCounts);
					bw1.close();
					bw2.close();
					errorRates = new ArrayList<>();
					spikeCounts = new ArrayList<>();
				}
				stepNumber++;
			}
			if(j%10==0){
				Utilities.writeMatrixToFile(null, ConfigurationParameters.kernelCoeffFileName +".txt", net.kernelMgr.kernelCoefficients);
			}
		}
		System.out.println("done here");
	}

	private static void saveMovingErrorAverages(int high, int low, List<Double> movingErrAvg) throws IOException {
		// TODO Auto-generated method stub
		// TODO Auto-generated method stub
		FileWriter fw1 = new FileWriter(
				ConfigurationParameters.movingerroravgFileName + "-" + high + "-" + low + ".txt");
		BufferedWriter bw1 = new BufferedWriter(fw1);
		Utilities.writeListToBufferedWriter(bw1, movingErrAvg);
		bw1.close();
	}

	private static void trainOnSinglePartition(int partitionNumber) throws Exception {
		DataReader dtr = new DataReader(net);
		for (int ck = 0; ck < ConfigurationParameters.numberOfPartitionChunks; ck++) {
			dtr.loadPartitionIntoMemory(partitionNumber, ck);
			for (int i = 0; i < dtr.allSignalPieces.length; i++) {
				net.initFromFileData(dtr.allSignalPieces[i], dtr.allDifferentialSignalPieces[i],
						dtr.signalByKernelComp[i]);
				net.calculateSpikeTimesAndReconstructSignal();
				net.calculateErrorGradient();
				double errorrate = net.updateKernelCoefficients();
				// errorValues.add(net.calculateError());
				errors.add(errorrate);
				spikeCounts.add(net.spikeTimings.size());
				// as per requirement add the spike count and error avg
				System.out.println("error rate at step number:" + stepNumber + " is:" + errorrate);
				if (stepNumber % 1000 == 0) {
					Utilities.writeMatrixToFile(null,
							ConfigurationParameters.kernelCoeffFileName + "-" + stepNumber + ".txt",
							net.kernelMgr.kernelCoefficients);
					FileWriter fw1 = new FileWriter(ConfigurationParameters.errorValuesFileName + ".txt", true);
					FileWriter fw2 = new FileWriter(ConfigurationParameters.spikeCountFileName + ".txt", true);
					BufferedWriter bw1 = new BufferedWriter(fw1);
					Utilities.writeListToBufferedWriter(bw1, errors);
					BufferedWriter bw2 = new BufferedWriter(fw2);
					Utilities.writeIntegerListToBufferedWriter(bw2, spikeCounts);
					bw1.close();
					bw2.close();
					errors = new ArrayList<>();
					spikeCounts = new ArrayList<>();
				}
				stepNumber++;
				System.out.println(stepNumber + ": steps executed");
			}
		}
		Utilities.writeMatrixToFile(null, ConfigurationParameters.kernelCoeffFileName + ".txt",
				net.kernelMgr.kernelCoefficients);
	}

	private static void writeErrorValues(int partitionNumber, List<Double> errorValues, List<Integer> spikeCounts)
			throws IOException {
		// TODO Auto-generated method stub
		FileWriter fw1 = new FileWriter(ConfigurationParameters.errorValuesFileName + "-" + partitionNumber + ".txt");
		FileWriter fw2 = new FileWriter(ConfigurationParameters.spikeCountFileName + "-" + partitionNumber + ".txt");
		BufferedWriter bw1 = new BufferedWriter(fw1);
		Utilities.writeListToBufferedWriter(bw1, errorValues);
		BufferedWriter bw2 = new BufferedWriter(fw2);
		Utilities.writeIntegerListToBufferedWriter(bw2, spikeCounts);
		bw1.close();
		bw2.close();
	}

	private static void saveKernelCoefficients(int partitionNumber, int stepNumber) throws IOException {
		// TODO Auto-generated method stub
		FileWriter fw = new FileWriter(
				ConfigurationParameters.kernelCoeffFileName + /*"-" + partitionNumber + "-" + stepNumber +*/ ".txt");
		BufferedWriter bw = new BufferedWriter(fw);
		Utilities.writeMatrixToBufferedWriter(bw, net.kernelMgr.kernelCoefficients);
		bw.close();
	}

	private static void runOnTrainingSet(int stepNumber) throws Exception {
		// steps denotes the number of steps that has already been executed
		int steps = stepNumber;
		// TODO Auto-generated method stub
		DataReader dtr = new DataReader();
		int numberOfStepsInOnePartition = ConfigurationParameters.numberOfFilesForPartition
				* ConfigurationParameters.numberOfSignalSegmentsInOneFile;
		int startPartitionIndex = ConfigurationParameters.TRAINING_START_PARTITION
				+ stepNumber / numberOfStepsInOnePartition;
		int numberofStepsInOnePartitionChunk = numberOfStepsInOnePartition
				/ ConfigurationParameters.numberOfPartitionChunks;
		int startPartitionChunkIndex = (stepNumber % (numberOfStepsInOnePartition))
				/ (numberofStepsInOnePartitionChunk);
		int startingSignalPieceindex = (stepNumber % (numberOfStepsInOnePartition))
				% (numberofStepsInOnePartitionChunk);
		dtr.loadPartitionIntoMemory(startPartitionIndex, startPartitionChunkIndex);
		for (int i = startingSignalPieceindex; i < dtr.allSignalPieces.length; i++) {
			net.initFromFileData(dtr.allSignalPieces[i], dtr.allDifferentialSignalPieces[i], dtr.signalByKernelComp[i],
					dtr.signalByDiffKernelComp[i]);
			net.calculateSpikeTimesAndReconstructSignal();
			net.calculateErrorGradient();
			net.updateKernelCoefficients();
			errors.add(net.calculateError());
			steps++;
			if (steps % ConfigurationParameters.OUTPUT_STEPS_INTERVAL == 0) {
				net.saveNetworkState(steps);
			}
		}
		for (int chunkPart = startPartitionChunkIndex
				+ 1; chunkPart < ConfigurationParameters.numberOfPartitionChunks; chunkPart++) {
			dtr.loadPartitionIntoMemory(startPartitionIndex, chunkPart);
			for (int i = 0; i < dtr.allSignalPieces.length; i++) {
				net.initFromFileData(dtr.allSignalPieces[i], dtr.allDifferentialSignalPieces[i],
						dtr.signalByKernelComp[i], dtr.signalByDiffKernelComp[i]);
				net.calculateSpikeTimesAndReconstructSignal();
				net.calculateErrorGradient();
				net.updateKernelCoefficients();
				errors.add(net.calculateError());
				steps++;
				if (steps % ConfigurationParameters.OUTPUT_STEPS_INTERVAL == 0) {
					net.saveNetworkState(steps);
				}
			}
		}
		for (int p = startPartitionIndex + 1; p <= ConfigurationParameters.TRAINING_END_PARTITION; p++) {
			for (int pc = 0; pc < ConfigurationParameters.numberOfPartitionChunks; pc++) {
				dtr.loadPartitionIntoMemory(p, pc);
				for (int i = 0; i < dtr.allSignalPieces.length; i++) {
					net.initFromFileData(dtr.allSignalPieces[i], dtr.allDifferentialSignalPieces[i],
							dtr.signalByKernelComp[i], dtr.signalByDiffKernelComp[i]);
					net.calculateSpikeTimesAndReconstructSignal();
					net.calculateErrorGradient();
					net.updateKernelCoefficients();
					errors.add(net.calculateError());
					steps++;
					if (steps % ConfigurationParameters.OUTPUT_STEPS_INTERVAL == 0) {
						net.saveNetworkState(steps);
					}
				}
			}
		}
	}

	private static void runOnTestSet(int stepNumber) throws Exception {
		// TODO Auto-generated method stub
		// steps denotes the number of steps that has already been executed
		int steps = stepNumber;
		// TODO Auto-generated method stub
		DataReader dtr = new DataReader();
		int numberOfStepsInOnePartition = ConfigurationParameters.numberOfFilesForPartition
				* ConfigurationParameters.numberOfSignalSegmentsInOneFile;
		int startPartitionIndex = ConfigurationParameters.TRAINING_START_PARTITION
				+ stepNumber / numberOfStepsInOnePartition;
		int numberofStepsInOnePartitionChunk = numberOfStepsInOnePartition
				/ ConfigurationParameters.numberOfPartitionChunks;
		int startPartitionChunkIndex = (stepNumber % (numberOfStepsInOnePartition))
				/ (numberofStepsInOnePartitionChunk);
		int startingSignalPieceindex = (stepNumber % (numberOfStepsInOnePartition))
				% (numberofStepsInOnePartitionChunk);
		dtr.loadPartitionIntoMemory(startPartitionIndex, startPartitionChunkIndex);
		for (int i = startingSignalPieceindex; i < dtr.allSignalPieces.length; i++) {
			net.initFromFileData(dtr.allSignalPieces[i], dtr.allDifferentialSignalPieces[i], dtr.signalByKernelComp[i],
					dtr.signalByDiffKernelComp[i]);
			net.calculateSpikeTimesAndReconstructSignal();
			// net.calculateErrorGradient();
			// net.updateKernelCoefficients();
			errors.add(net.calculateError());
			steps++;
			if (steps % ConfigurationParameters.OUTPUT_STEPS_INTERVAL == 0) {
				net.saveNetworkState(steps);
			}
		}
		for (int chunkPart = startPartitionChunkIndex
				+ 1; chunkPart < ConfigurationParameters.numberOfPartitionChunks; chunkPart++) {
			dtr.loadPartitionIntoMemory(startPartitionIndex, chunkPart);
			for (int i = 0; i < dtr.allSignalPieces.length; i++) {
				net.initFromFileData(dtr.allSignalPieces[i], dtr.allDifferentialSignalPieces[i],
						dtr.signalByKernelComp[i], dtr.signalByDiffKernelComp[i]);
				net.calculateSpikeTimesAndReconstructSignal();
				// net.calculateErrorGradient();
				// net.updateKernelCoefficients();
				errors.add(net.calculateError());
				steps++;
				if (steps % ConfigurationParameters.OUTPUT_STEPS_INTERVAL == 0) {
					net.saveNetworkState(steps);
				}
			}
		}
		for (int p = startPartitionIndex + 1; p <= ConfigurationParameters.TRAINING_END_PARTITION; p++) {
			for (int pc = 0; pc < ConfigurationParameters.numberOfPartitionChunks; pc++) {
				dtr.loadPartitionIntoMemory(p, pc);
				for (int i = 0; i < dtr.allSignalPieces.length; i++) {
					net.initFromFileData(dtr.allSignalPieces[i], dtr.allDifferentialSignalPieces[i],
							dtr.signalByKernelComp[i], dtr.signalByDiffKernelComp[i]);
					net.calculateSpikeTimesAndReconstructSignal();
					// net.calculateErrorGradient();
					// net.updateKernelCoefficients();
					errors.add(net.calculateError());
					steps++;
					if (steps % ConfigurationParameters.OUTPUT_STEPS_INTERVAL == 0) {
						net.saveNetworkState(steps);
					}
				}
			}
		}
	}
}
