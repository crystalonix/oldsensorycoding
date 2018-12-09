package sensoryCoding.network;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYDataset;

public class SignalUtils {

	public static double[] multiply() {
		// TODO Auto-generated method stub
		return null;
	}

	/**
 	 * This method calculates the differential signal for the given input signal
 	 * @param inputSignal // TODO: move to signal utils
 	 * @return differential signal
 	 */
	public static Signal calculateSignalDifferential(Signal inputSignal) {
		double [] differentialSignal = new double[inputSignal.getLength()];
		double shift = inputSignal.getTimeShift();
		for(double i=inputSignal.start+1; i<=inputSignal.end; i++){
			// here I am taking the left difference
			differentialSignal[(int)(i-shift)]=(inputSignal.getSignalValue(i)-inputSignal.getSignalValue(i-1))/ConfigurationParameters.TIME_STEP;
		}
		// the first point should get the right difference
		differentialSignal[(int)(inputSignal.start-shift)]= (inputSignal.getSignalValue(inputSignal.start+1)-inputSignal.getSignalValue(inputSignal.start))/ ConfigurationParameters.TIME_STEP;
		return new Signal(differentialSignal, inputSignal.start, inputSignal.end, inputSignal.getTimeShift());
	}

	/**
	 * This method calculates the integral of the given signal
	 *
	 * @return
	 */
	public static double calculateSignalIntegral(Signal input) {
		double summation = 0;
		for (double i = input.start; i <= input.end; i++) {
			summation += input.getSignalValue(i);
		}
		return summation*ConfigurationParameters.TIME_STEP;
	}

	/**
	 * This method displays the signal as comma separated values
	 * @param signal
	 */
	public static void outputSignal(Signal signal){
		System.out.println("***see the signal***");
		if(signal == null) {
			System.out.println("Signal is null and cannot be displayed");
			return;
		}
		for (int i = 0; i < signal.getLength(); i++) {
			System.out.print(signal.data[i]);
			if(i!=(signal.getLength()-1)){
				System.out.print(",");
			}
		}
		System.out.println("***see the signal***");
	}

	/**
	 * Checks whether two signals overlap based on the start and end time
	 * @param signal1
	 * @param signal2
	 * @return
	 */
	public static boolean doesOverlap(Signal signal1, Signal signal2) {
		if(signal1.start > signal2.end || signal1.end < signal2.start){
			return false;
		}
		return true;
	}


	/**
	 * This method returns a signal shifted by delta in time
	 * @param originalSignal
	 * @param timeShift
	 * @return time shifted signal
	 */
	public static Signal shiftSignal(Signal originalSignal, double delta){
		return new Signal(originalSignal.data, originalSignal.start + delta, originalSignal.end + delta, originalSignal.getTimeShift() + delta);
	}

	/**
	 * takes a signal B(t) and returns the signal B(t*alpha)
	 *
	 * @param originalSignal
	 * @param alpha
	 * @return scaled signal
	 */
	public static Signal scaleSignal(Signal originalSignal, double alpha) {
		int newStart = (int) (originalSignal.start / alpha);
		int newEnd = (int) (originalSignal.end / alpha);
		double[] data = new double[originalSignal.getLength()];
		for (int i = newStart; i <= newEnd; i++) {
			data[i] = originalSignal.getSignalValue((int) (i * alpha));
		}
		return new Signal(data, newStart, newEnd);
	}

	/**
	 * Returns the summation of two signals
	 * @param signal1
	 * @param signal2
	 * @return
	 */
	public static Signal addTwoSignals(Signal signal1, Signal signal2){
		int dataLength = signal1.getLength();
		double start = Math.min(signal1.start, signal2.start);
		double end = Math.max(signal1.end, signal2.end);
		double data[] = new double[signal1.getLength()];
		for(int i = (int)(signal1.start-start); i <= (signal1.end-start) && i<dataLength; i++){
			data[i]+= signal1.getSignalValue(i+start);
		}
		for(int i = (int)(signal2.start-start); i <= (signal2.end-start) && i<dataLength; i++){
			data[i]+= signal2.getSignalValue(i+start);
		}
		return new Signal(data, start, Math.min(end, dataLength-1), start);
	}

	/**
	 * Returns the multiplication of the two signals
	 * p.s. there would be no time shift
	 * @param signal1
	 * @param signal2
	 * @return
	 */
	public static Signal multiplyTwoSignals(Signal signal1, Signal signal2){
		double start = Math.max(signal1.start, signal2.start);
		double end = Math.min(signal1.end, signal2.end);
		double timeshift = start;
		double data[] = new double[signal1.getLength()];
		for(double i = start; i <= end; i++){
			data[(int)(i-timeshift)] = (signal1.getSignalValue(i))*(signal2.getSignalValue(i));
		}
		return new Signal(data, start, end, timeshift);
	}

	/**
	 * returns a signal which is the result of scalar multiplication of the given signal with b
	 * @param originalSignal
	 * @param b scalar multiplier
	 * @return resulting signal which is obtained by multiplying the original signal with b
	 */
	public static Signal scalarMultiply(Signal originalSignal, double b){
		double [] newData = new double[originalSignal.getLength()];
		for(double i=originalSignal.start; i<=originalSignal.end; i++){
			newData[(int)(i-originalSignal.getTimeShift())] = originalSignal.getSignalValue(i)*b;
		}
		return new Signal(newData, originalSignal.start, originalSignal.end, originalSignal.getTimeShift());
	}

	/**
	 * This function displays the line chart for a given set of signals
	 * @param dataset
	 * @param chartTitle
	 * @param panelTitle
	 */
	public static void displaySignal(XYDataset dataset, String chartTitle) {
		JFreeChart lineChart = ChartFactory.createXYLineChart(chartTitle, "time", "value", dataset);
		ChartPanel cpanel = new ChartPanel(lineChart);
		cpanel.setVisible(true);
	}

	public static Signal readSignalFromFile(String fileName) throws Exception{
		double data[] = null;
		BufferedReader bstream = null;
		double start =0;
		double end =0;
		double timeShift =0;
		try {
			bstream = new BufferedReader(new FileReader(fileName));
			start = Double.parseDouble(bstream.readLine());
			end = Double.parseDouble(bstream.readLine());
			timeShift = Double.parseDouble(bstream.readLine());
			int dataLength = Integer.parseInt(bstream.readLine());
			data = new double[dataLength];
			for(int i=0; i<dataLength; i++){
				data[i] = Double.parseDouble(bstream.readLine());
			}
		} catch (Exception e) {
			throw new IllegalArgumentException("error reading a file:"+fileName);
		} finally {
			try {
				bstream.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return  new Signal(data, start, end, timeShift);
	}

	public static void storeSignal(String outputFileName, Signal signalPiece) throws Exception {
		// TODO Auto-generated method stub
		FileWriter fw = null;
		BufferedWriter bw = null;
		try {
			fw = new FileWriter(outputFileName);
			bw = new BufferedWriter(fw);
			bw.write(String.valueOf(signalPiece.start));
			bw.write('\n');
			bw.write(String.valueOf(signalPiece.end));
			bw.write('\n');
			bw.write(String.valueOf(signalPiece.getTimeShift()));
			bw.write('\n');
			bw.write(String.valueOf(signalPiece.data.length));
			bw.write('\n');
			for (int time = (int) (signalPiece.start - signalPiece.getTimeShift()); time <= (int) (signalPiece.end
					- signalPiece.getTimeShift()); time++) {
				bw.write(String.valueOf(signalPiece.data[time]));
				bw.write('\n');
			}
		} catch (IOException e) {
			throw new IllegalArgumentException("error in writing the signal to data file:" + outputFileName);
		} finally {
			try {
				bw.close();

				// bw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	public static void storeSignal(String outputFileName, Signal signalPiece, boolean includeHeader) throws Exception {
		// TODO Auto-generated method stub
		FileWriter fw = null;
		BufferedWriter bw = null;
		try {
			fw = new FileWriter(outputFileName);
			bw = new BufferedWriter(fw);
			if (includeHeader) {
				bw.write(String.valueOf(signalPiece.start));
				bw.write('\n');
				bw.write(String.valueOf(signalPiece.end));
				bw.write('\n');
				bw.write(String.valueOf(signalPiece.getTimeShift()));
				bw.write('\n');
				bw.write(String.valueOf(signalPiece.data.length));
				bw.write('\n');
			}
			for (int time = (int) (signalPiece.start - signalPiece.getTimeShift()); time <= (int) (signalPiece.end
					- signalPiece.getTimeShift()); time++) {
				bw.write(String.valueOf(signalPiece.data[time]));
				bw.write('\n');
			}
		} catch (IOException e) {
			throw new IllegalArgumentException("error in writing the signal to data file:" + outputFileName);
		} finally {
			try {
				bw.close();

				// bw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	/**
	 * Function that returns the squared norm of the signal
	 * @param input
	 * @return
	 */
	public static double calculateSquaredNorm(Signal input) {
		double summation = 0;
		for (double i = input.start; i <= input.end; i++) {
			double value = input.getSignalValue(i);
			summation += value*value;
		}
		return summation*ConfigurationParameters.TIME_STEP;
	}

	/**
	 * This method computes the convolution at single point considering the time
	 * step into account
	 *
	 * @param kernelSignal
	 * @param inputSignal
	 * @param timeStep
	 * @return
	 */
	public static double calculateAccurateConvolutionWithTimeStep(Signal kernelSignal, Signal inputSignal,
			double timeStep) {
		double total = 0;
		for (double t = kernelSignal.start; t <= kernelSignal.end; t += timeStep) {
			total += kernelSignal.getSignalValue(t) * inputSignal.getSignalValue(t) * timeStep;
		}
		return total;
	}
}
