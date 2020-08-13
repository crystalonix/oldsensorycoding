package sensoryCoding.network;

public class ConfigurationParameters {
	/*******************************/
	/***** these parameters need ****/
	/***** to be checked before ******/
	/********** every run ************/
	// indicates the time step for the signals in ms time
	public static final double TIME_STEP = 1;
	// nature how the kernel frequencies are separated: 1 for AP and 2 for GP
	public static int NATURE_OF_KERNEL_SPREAD = 1;
	// This denotes the number of kernels used for the experiment
	public static int numberOfKernels = 50;
	// This denotes the number of kernel components
	public static int numberofKernelComponents = 6;
	// This denotes the length of the basic spline function used to generate the
	// kernels
	public static int lengthOfBasicBSpline = (int) (3 / TIME_STEP);
	// denotes the factor to which the frequency of the next kernel should be
	// incremented
	public static double FREQUENCY_SCALING_FACTOR = 1;
	// This indicates the frequency at which all the signals are sampled (in kHz)
	public static int SAMPLING_FREQUENCY = 1;
	// This denotes the length of each component signal
	public static int lengthOfComponentSignals = 8820;
	// This denotes the number of threads spawned in different instances of the code
	public static int numberOfThreads = 100;
	// number of chunks into which a single partition has to be broken before
	// loading into memory
	public static int numberOfPartitionChunks = 1;
	// This is the absolute path to the root folder specific to the machine
	// Users⁩ ▸ ⁨anonymouschattopadhyay⁩ ▸ ⁨Desktop⁩ ▸ ⁨research⁩ ▸ ⁨oldsensorycoding⁩ ▸
	// ⁨SensoryCoding⁩ ▸ ⁨src⁩ ▸ ⁨sensoryCoding⁩ ▸
	public static String ABSOLUTE_MACHINE_PATH = "/Users/anonymouschattopadhyay⁩/Desktop/research/oldsensorycoding/SensoryCoding⁩/⁨src⁩/⁨sensoryCoding⁩/";// "/cise/research/compneuro3/anonymous/version2/setup1/proc3/SensoryCoding/src/sensoryCoding/";
	// "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding
	// (2)/SensoryCoding/src/sensoryCoding/";
	// Absolute path to the folder where at each step output is stored
	public static String STEP_OUTPUT_FOLDER_PATH = ABSOLUTE_MACHINE_PATH + "stepOutPut/testAfterTrainOnDataPart14/";
	/*******************************/
	/*******************************/
	/*******************************/

	/****** play with the *********/
	/******* threshold values *******/
	/*******************************/
	// coeffcients to initialize the kernel with
	public static double initialCoefficientValue = 1;
	// This denotes the value at which the threshold for each spike generation is
	// initially set to
	public static double initialThresHoldValue = 0.1;
	// rate at which the coefficients need to be updated
	public static double INITIAL_LEARNING_RATE = 0.001;
	// this specifies the threshold value of the gradient vector length beyond which
	// we tune the gradient values
	public static double INITIAL_THRESHOLD_GRADIENT_STEP_LENGTH = 0.01;
	// this is high value at which the threshold is kicked up when a spike occurs
	public static double AHP_HIGH_VALUE = 1000.0;
	// refractory period within which the threshold comes down to its original value
	public static double AHP_REFRACTORY_PERIOD = 50.0;
	// stores the positive slope of the ahp graph
	public static double AHP_SLOPE = AHP_HIGH_VALUE / AHP_REFRACTORY_PERIOD;
	// constant factor with which the ahp function is multiplied with
	public static int AHP_CONSTANT = 10;
	// time constant for ahp function
	public static double TIME_CONSTANT = 1.5 * 44.1;
	/*******************************/
	/*******************************/
	/*******************************/

	/*******************************/
	/******** Mode parameters ********/
	/*******************************/
	// defines the mode of running the application
	public static final MODE USER_MODE = MODE.APPLICATION_MODE;
	// directs if the fail-safe option need to be put in the gradient update
	public static boolean FAIL_SAFE_ON = true;
	// denotes whether the values for this signal can be picked from the cache
	public static boolean IS_SIGNAL_CACHED = false;
	// indicates if the kernel coefficients has to be randomized
	public static boolean SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = false;
	// for testing on a large corpora we set this to 0 because since we are applying
	// stochastic gradient descent this removes all the local minima problems
	public static int NUMBER_OF_KICKS = 0;
	// directs if for each iteration of network we should collect spike statistics
	public static boolean SHOULD_COLLECT_SPIKE_STATISTICS = true;
	// when applying conjugate gradient it checks the residuals
	public static boolean SHOULD_CHECK_RESIDUALS = false;
	// this flag denotes whether a pre-conditioner should be applied
	public static boolean USE_PRECONDITIONER = false;
	/*******************************/
	/*******************************/
	/*******************************/

	// denotes the factor to which two component BSplines should overlap
	public static double OVERLAP_FACTOR = 1.0 / 3.0;
	// number of signal segments present in file
	public static int numberOfSignalSegmentsInOneFile = 220500 / lengthOfComponentSignals;
	// number of partitions into which the entire data set has to be broken at
	// processing stage
	public static int numberOfDataPartitions = 20;
	// number of data files read for each partition
	public static int numberOfFilesForPartition = 100;
	// name attached to the text file for a signal piece
	public static String signalPieceName = "signalpiece-";
	public static String diffSignalPieceName = "diffSignalpiece-";
	public static String signalByBName = "signalByB-";
	public static String signalByBDiffName = "signalByBDiff-";
	public static int numberOfFIlesInAPartitionChunk = numberOfFilesForPartition / numberOfPartitionChunks;
	public static String originalFileName = "originalSignalFile";
	public static String reconstructedFileName = "reconstructedSignalFile";

	// This title is for display on any graph being plotted
	public static final String APPLICATION_TITLE = "Sensory Coding Signals Plot";
	// labels the time axis
	public static final String TIME_LABEL = "Time in ms";
	// labels the signal values axis
	public static final String VALUE_LABEL = "Signal value";
	// file where kernel coeffs are stored
	public static String kernelCoeffFileName = STEP_OUTPUT_FOLDER_PATH + "kernelCoefficients";
	// path to the folder where the b-spline coefficients are stored
	public static String bsplineCoeffsPath = ABSOLUTE_MACHINE_PATH + "bsplineCoefficients⁩/";
	// file where kernel coeffs are stored
	public static final String movingerroravgFileName = STEP_OUTPUT_FOLDER_PATH + "errorAverages";
	// file where kernel coeffs are stored
	public static String errorValuesFileName = STEP_OUTPUT_FOLDER_PATH + "errorValues";
	// file where kernel coeffs are stored
	public static String spikeCountFileName = STEP_OUTPUT_FOLDER_PATH + "spikeCounts";
	// This stores the absolute path where the application is running
	public static String ABSOLUTE_FILE_PATH = ABSOLUTE_MACHINE_PATH + "output/";
	// This is the path for test data folder
	public static String SOUND_DATA_FOLDER_PATH = ABSOLUTE_MACHINE_PATH + "soundData/";
	// This is the path to the folder to processed signal data partitions
	public static String PARTITION_PATH = ABSOLUTE_MACHINE_PATH + "processedData/dataPartition-";
	// "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/dataBackUpForAllPartitions/6k6cDataParts/dataPartition-";
	// This is the path to the folder to processed signal data partitions
	public static String PARTITION_SEED_FILE = ABSOLUTE_MACHINE_PATH + "partitionSeed.txt";
	// This is the path to the file where the state of the network is stored
	public static final String STATE_OF_NETWORK = ABSOLUTE_MACHINE_PATH + "state.txt";
	// file path of where the log file will be stored
	public static final String ERRORGRAD_LOG_FILE_NAME = ABSOLUTE_FILE_PATH + "sensoryLogFileErrorGard";
	// file path of where the log file will be stored
	public static final String KERNELCOEFF_LOG_FILE_NAME = ABSOLUTE_FILE_PATH + "sensoryLogFileKernelCoeffs";
	// file path for where the signal image files will be stored
	public static final String SIGNAL_IMAGE_FILE_NAME = ABSOLUTE_FILE_PATH + "ReconstructedSignalStep-";
	// Number of steps after which statistics will be calculated
	public static final int SAMPLING_STEPS = 1;
	// Total number of iterations
	public static final int NUMBER_OF_ITERATIONS = 10000000;
	// Denotes the height of the signal image being stored in the png file
	public static final int IMAGE_WIDTH = 1000;
	// Denotes the width of the signal image being stored in the png file
	public static final int IMAGE_HEIGHT = 500;
	// Denotes the mode of operation
	public static final OP_MODE OPERATION_MDOE = OP_MODE.TRAIN_MODE;
	// Index of the first partition for training set
	public static final int TRAINING_START_PARTITION = 1;
	// Index of the last partition for training set
	public static final int TRAINING_END_PARTITION = 12;
	// Index of the first partition for training set
	public static final int TEST_START_PARTITION = 13;
	// Index of the last partition for training set
	public static final int TEST_END_PARTITION = 20;
	// Interval of number of steps after which we need to store the state of the
	// network
	public static final int OUTPUT_STEPS_INTERVAL = 10000;
	public static final boolean LOCAL_MINIMA_MODE = false;
	// This is the time step for accurately computing the convolutions
	public static final double TIME_STEP_FOR_CONVOLUTIONS = 0.1;
	// This is the field that tells whether the differential comps should be loaded
	public static final boolean loadDifferentialComps = false;
	// file where average value of the coefficients of reconstruction are stored
	public static String avgReconsCoeffsFileName = STEP_OUTPUT_FOLDER_PATH + "avgCoeffsOfRecons.txt";
	// file where avgerage spike count values are stored
	public static String avgSpikeCountsFileName = STEP_OUTPUT_FOLDER_PATH + "avgSpikeCounts.txt";
	// This field decides whether least square method has to be applied to find the
	// coefficients
	public static boolean USE_LEAST_SQUARE = false;
	// tells whether the signal has to be zero padded for calculating signal kernel
	// convolution
	public static boolean ZERO_PADDING_NEEDED = true;
}

enum MODE {
	DEBUG_MODE, DISPLAY_MODE, APPLICATION_MODE, TEST_MODE
}

enum OP_MODE {
	TEST_MODE, TRAIN_MODE
}
