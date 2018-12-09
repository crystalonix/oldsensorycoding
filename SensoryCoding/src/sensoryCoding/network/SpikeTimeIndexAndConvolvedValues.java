package sensoryCoding.network;

public class SpikeTimeIndexAndConvolvedValues {
	public double spikeTime;
	public int kernelIndex;
	public double convolvedValue;
	public SpikeTimeIndexAndConvolvedValues(double spikeTime, int kernelIndex, double convolvedValue){
		this.kernelIndex = kernelIndex;
		this.spikeTime = spikeTime;
		this.convolvedValue = convolvedValue;
	}
}
