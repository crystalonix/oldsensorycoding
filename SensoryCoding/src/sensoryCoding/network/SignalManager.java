package sensoryCoding.network;

public class SignalManager {
	// this is the original sound signal that we break and piecewise feed that into our algorithm
	private double[] timelessSignal;
	private boolean isZeroPadded = false;
	private int lengthOfSignalPiece;
	private static int pieceCount = 0;

	/**
	 *
	 * @param timelessSignal
	 */
	public SignalManager(double [] timelessSignal){
		this.timelessSignal = timelessSignal;
		this.lengthOfSignalPiece = ConfigurationParameters.lengthOfComponentSignals;
	}

	/**
	 * Constructor of Signal Manager
	 * @param timelessSignal
	 * @param lengthOfSignalPiece
	 */
	public SignalManager(double [] timelessSignal, int lengthOfSignalPiece){
		this.timelessSignal = timelessSignal;
		this.lengthOfSignalPiece = ConfigurationParameters.lengthOfComponentSignals;
		this.lengthOfSignalPiece = lengthOfSignalPiece;
	}
	/**
	 * This method initializes the signal manager
	 */
	public void initializeSignalManager(){
		// TODO: this function should initialize the signal manager with
	}

	/**
	 * This method returns a piece of signal from the original signal
	 * without zero padding at both the ends
	 * @return
	 */
	public double[] getPieceOfSignalWithOutZeroPadding(){
		if((pieceCount+1)*lengthOfSignalPiece > timelessSignal.length){
			pieceCount = 0;
		}
		double signalPiece[] = new double[lengthOfSignalPiece];
		// copy the values from the timelessSignal into the signalPiece
		for(int i= pieceCount*lengthOfSignalPiece, j=0; i< (pieceCount+1)*lengthOfSignalPiece; i++, j++){
			signalPiece[j] = timelessSignal[i];
		}
		pieceCount++;
		return signalPiece;
	}
}
