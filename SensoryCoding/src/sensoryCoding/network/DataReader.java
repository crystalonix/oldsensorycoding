package sensoryCoding.network;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeriesCollection;

public class DataReader {

	public DataReader() throws Exception {
		Signal blankSignal = null;
		network = new Network(blankSignal);
	}

	public DataReader(Network net) throws Exception {
		network = net;
	}

	public static void main(String[] args) throws Exception {
		int val = 5;
		if(val ==0){
			Signal blankSignal = null;
			// check kernels
			Network testNet = new Network(blankSignal);
			testNet.kernelMgr.displayInvertedKernels();
		}
		else if (val == 1) {
			// task 1: display kernel after execution
			String filename = "kernelCoefficients-5-2500.txt";
			displayKernelsFromFile(filename, "updated kernel");
			filename = "kernelCoefficients-1-0.txt";
			displayKernelsFromFile(filename, "initial kernel");
		}
		else if(val == 2){
			// read moving error average from file and display
			String fileName = "errorAverages-15-7.txt";
			FileReader fr = new FileReader(ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + fileName);
			BufferedReader br = new BufferedReader(fr);
			List<Double> errorValues = Utilities.readListFromFile(br);
			double [] dataInitialErrors = new double [100];
			double [] dataAfter100Steps = new double [errorValues.size()-100];
			for(int i=1; i<errorValues.size(); i++){
				if(i<=100){
					dataInitialErrors[i-1] = errorValues.get(i);
				}
				else{
					dataAfter100Steps[i-101] = errorValues.get(i);
				}
			}
			(new Signal(dataInitialErrors, 0, 99)).DrawSignal("Error avg. in inital 100 steps experiment 1");
			(new Signal(dataAfter100Steps, 100, errorValues.size(), 100)).DrawSignal("Error avg. after first 100 steps experiment 1");
			/*double [][] allSignalsData = new double[errorValues.size()/1000+1][];

			for(int i=1; i<errorValues.size(); i++){
				if(i %1000 ==0){
					allSignalsData[(i-1)/1000] = data;
					data = new double[1000];
				}
				data[i%1000] = errorValues.get(i);
			}
			allSignalsData[allSignalsData.length-1]= data;
			for(int i = 0;i<allSignalsData.length; i++){
			(new Signal(allSignalsData[i], i*1000, (i+1)*1000, i*1000)).DrawSignal("Moving avg of error Values part"+i);
			}*/
			br.close();
		}
		else if(val==4){
			int k = 900;
			int partitionNumber = 15;
			String f1 = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "originalSignalFile-"+partitionNumber+"-"+k+".txt";
			String f2 = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "reconstructedSignalFile-"+partitionNumber+"-"+k+".txt";
			Signal orig = SignalUtils.readSignalFromFile(f1);
			Signal recon = SignalUtils.readSignalFromFile(f2);
			/*orig.DrawSignal("original signal after "+k+" steps");
			recon.DrawSignal("reconsturcted signal after "+k+" steps");
			*/
			XYSeriesCollection kernelDataSet = new XYSeriesCollection();
			kernelDataSet.addSeries(orig.getSignalDisplayData("Original Signal at sample:"+k+" of partition:"+partitionNumber));
			kernelDataSet.addSeries(recon.getSignalDisplayData("Reconstructed Signal at sample:"+k+" of partition:"+partitionNumber));
			// JFreeChart chart = ChartFactory.createXYLineChart(, , "value",
			// kernelDataSet);
			JFreeChart linechart = ChartFactory.createXYLineChart("Original vs Reconstructed Signal for sample#"+k+" of partition:"+partitionNumber, "time", "value", kernelDataSet);
			String fileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "signalRecons-"+partitionNumber+"-"+k+".png";
			Utilities.saveLineChartToFile(fileName, linechart, null, null);
		}
		else {
			// task 2
			DataReader dtr = new DataReader();
			// generatePartitionSeed();
			dtr.processSignalDataPartition(20, 11);
		}
	}

	private static void displayKernelsFromFile(String filename, String title) throws Exception {
		// TODO Auto-generated method stub
		FileReader fr = new FileReader(ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + filename);
		BufferedReader br = new BufferedReader(fr);
		double[][] kernelCoeffs = Utilities.readMatrixToBufferedWriter(br);
		KernelManager krMg = (new KernelManager(ConfigurationParameters.numberOfKernels,
				ConfigurationParameters.numberofKernelComponents, ConfigurationParameters.lengthOfComponentSignals));
		//krMg.displayKernels();
		krMg.initKernelCoefficients(kernelCoeffs);
		br.close();
		krMg.displayKernels(title);
	}

	/**
	 * This method reads signal data from an input data folder, partitions
	 *
	 * @param partitionNumber
	 * @throws Exception
	 */
	public void processSignalDataPartition(int high, int low) throws Exception {

		List<Integer> randomizedPartitions = new ArrayList<>();
		Signal partitionSeed = SignalUtils.readSignalFromFile(ConfigurationParameters.PARTITION_SEED_FILE);
		for(int i=0; i< partitionSeed.data.length; i++){
			randomizedPartitions.add((int) partitionSeed.data[i]);
		}
		int partitionNumber = 0;
		/*String dataFolder = ConfigurationParameters.SOUND_DATA_FOLDER_PATH;
		File folder = new File(dataFolder);
		File[] listOfFiles = folder.listFiles();*/
		for (partitionNumber = low; partitionNumber <= high; partitionNumber++) {
			int startFileIndex = (partitionNumber-1) * ConfigurationParameters.numberOfFilesForPartition;
			int endFileIndex = startFileIndex + ConfigurationParameters.numberOfFilesForPartition;
			String outputFolderPath = ConfigurationParameters.PARTITION_PATH + partitionNumber + "/";
			for (int i = startFileIndex; i < endFileIndex; i++) {
				String selectedFileName = ConfigurationParameters.SOUND_DATA_FOLDER_PATH+"output"+randomizedPartitions.get(i)+".txt";
				// get the entire chunk of signals from one file
				double[][] signalsInPartition = readSignalPieces(selectedFileName);//listOfFiles[randomizedPartitions.get(i)]);
				// iterate through each signal piece individually
				for (int j = 0; j < signalsInPartition.length; j++) {
					Signal signalPiece = new Signal(signalsInPartition[j]);
					network.init(signalPiece);
					String outputFileName = outputFolderPath + ConfigurationParameters.signalPieceName
							+ (i - startFileIndex) + "-" + j + ".txt";
					String DifferentialSignalFileName = outputFolderPath + ConfigurationParameters.diffSignalPieceName
							+ (i - startFileIndex) + "-" + j + ".txt";
					SignalUtils.storeSignal(outputFileName, signalPiece);
					SignalUtils.storeSignal(DifferentialSignalFileName, network.differentialOfThisSignal);
					for (int k = 0; k < ConfigurationParameters.numberOfKernels; k++) {
						String signalKernelConvFileName = outputFolderPath + ConfigurationParameters.signalByBName
								+ (i - startFileIndex) + "-" + j + "-" + k + ".txt";
						// removing differential kernel convolutions part because it can be calculated from kernelConvolutions
						/*String signalByDiffKernelConvFileName = outputFolderPath
								+ ConfigurationParameters.signalByBDiffName + (i - startFileIndex) + "-" + j + "-" + k
								+ ".txt";*/
						SignalUtils.storeSignal(signalKernelConvFileName, network.signalKernelComponentConvolutions[k]);
						/*SignalUtils.storeSignal(signalByDiffKernelConvFileName,
								network.signalDifferentialKernelComponentConvolutions[k]);*/
					}
				}
			}
		}
	}

	public static void generatePartitionSeed() throws Exception{
		List<Integer> filePartitions = new ArrayList<>();
		String dataFolder = ConfigurationParameters.SOUND_DATA_FOLDER_PATH;
		File folder = new File(dataFolder);
		File[] listOfFiles = folder.listFiles();
		for (int i = 0; i < listOfFiles.length; i++) {
			filePartitions.add(i);
		}
		Collections.shuffle(filePartitions);
		double data[] = new double[filePartitions.size()];
		for(int i=0; i<data.length; i++){
			data[i] = filePartitions.get(i);
		}
		SignalUtils.storeSignal(ConfigurationParameters.PARTITION_SEED_FILE, new Signal(data));
	}

	/**
	 * This method loads the entire data set present in the partition into
	 * memory
	 * @throws Exception
	 */
	public void loadPartitionIntoMemory(int partitionNumber, int partitionChunkNumber) throws Exception {
		int startFileNumber = partitionChunkNumber * ConfigurationParameters.numberOfFIlesInAPartitionChunk;
		int endFileNumber = (partitionChunkNumber + 1) * ConfigurationParameters.numberOfFIlesInAPartitionChunk;
		int numberOfPiecesInOnePartChunk = ConfigurationParameters.numberOfFIlesInAPartitionChunk*ConfigurationParameters.numberOfSignalSegmentsInOneFile;
		allSignalPieces = new Signal[numberOfPiecesInOnePartChunk];
		allDifferentialSignalPieces = new Signal[ConfigurationParameters.numberOfFIlesInAPartitionChunk*ConfigurationParameters.numberOfSignalSegmentsInOneFile];
		signalByKernelComp = new Signal[numberOfPiecesInOnePartChunk][ConfigurationParameters.numberOfKernels];
		//signalByDiffKernelComp = new Signal[numberOfPiecesInOnePartChunk][ConfigurationParameters.numberOfKernels];
		int pieceIndex = 0;
		for (int fileNumber = startFileNumber; fileNumber < endFileNumber; fileNumber++) {
			for (int segmentNumber = 0; segmentNumber < ConfigurationParameters.numberOfSignalSegmentsInOneFile; segmentNumber++) {
				String inputFile = ConfigurationParameters.PARTITION_PATH + partitionNumber + "/"
						+ ConfigurationParameters.signalPieceName + fileNumber + "-" + segmentNumber + ".txt";
				String diffInputFile = ConfigurationParameters.PARTITION_PATH + partitionNumber + "/"
						+ ConfigurationParameters.diffSignalPieceName + fileNumber + "-" + segmentNumber + ".txt";
				allSignalPieces[pieceIndex] = SignalUtils.readSignalFromFile(inputFile);
				if(allSignalPieces[pieceIndex]==null){
					System.out.println("Null signal was encountered for file:"+ inputFile);
					return;
				}
				allDifferentialSignalPieces[pieceIndex] = SignalUtils.readSignalFromFile(diffInputFile);
				Signal[] signalByB = new Signal[ConfigurationParameters.numberOfKernels];
				//Signal[] signalByDiffB = new Signal[ConfigurationParameters.numberOfKernels];
				for (int k = 0; k < ConfigurationParameters.numberOfKernels; k++) {
					String sigByB = ConfigurationParameters.PARTITION_PATH + partitionNumber + "/"
							+ ConfigurationParameters.signalByBName + fileNumber + "-" + segmentNumber + "-"+k + ".txt";
					//String sigByDiffB = ConfigurationParameters.PARTITION_PATH + partitionNumber + "/"
					//		+ ConfigurationParameters.signalByBDiffName + fileNumber + "-" + segmentNumber +"-"+ k + ".txt";
					signalByB[k] = SignalUtils.readSignalFromFile(sigByB);
					//signalByDiffB[k] = SignalUtils.readSignalFromFile(sigByDiffB);
				}
				signalByKernelComp[pieceIndex] = signalByB;
				//signalByDiffKernelComp[pieceIndex] = signalByDiffB;
				pieceIndex++;
			}
		}
	}
	// network object used to perform computations with the signals accessing
	// the bspline components
	private Network network = null;
	/**
	 * Number of signal pieces present in a partition
	 */
	private int numberOfSignalPiecesInAPartition = ConfigurationParameters.numberOfSignalSegmentsInOneFile
			* ConfigurationParameters.numberOfFilesForPartition;
	/**
	 * List of all signal pieces to be loaded in memory
	 */
	public Signal[] allSignalPieces = null;
	/**
	 * List of all differential signal pieces to be loaded in memory
	 */
	public Signal[] allDifferentialSignalPieces = null;
	/**
	 * List of all signal by bspline convolutions to be loaded in memory
	 */
	public Signal[][] signalByKernelComp = null;//new Signal[numberOfSignalPiecesInAPartition][ConfigurationParameters.numberOfDataPartitions];
	/**
	 * List of all signal by differential bspline convolutions to be loaded in
	 * memory
	 */
	public Signal[][] signalByDiffKernelComp = null;//new Signal[numberOfSignalPiecesInAPartition][ConfigurationParameters.numberOfDataPartitions];

	/**
	 * This method reads the data corresponding to signal pieces from a .txt
	 * file
	 *
	 * @param selectedFile
	 * @return
	 */
	public static double[][] readSignalPieces(String selectedFileName) {
		double[][] data = new double[ConfigurationParameters.numberOfSignalSegmentsInOneFile][ConfigurationParameters.lengthOfComponentSignals];
		System.out.println("check the file name here:" + selectedFileName);
		BufferedReader bstream = null;
		try {
			bstream = new BufferedReader(new FileReader(selectedFileName));
			int counter = 0;
			String value = null;
			while ((value = bstream.readLine()) != null) {
				data[counter / (ConfigurationParameters.lengthOfComponentSignals)][counter
						% (ConfigurationParameters.lengthOfComponentSignals)] = Double.parseDouble(value);
				counter++;
			}
		} catch (Exception e) {
			throw new IllegalArgumentException("error in reading signal from file:"+ selectedFileName);
		} finally {
			try {
				bstream.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return data;
	}
}
