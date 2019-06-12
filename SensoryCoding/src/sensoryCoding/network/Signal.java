package sensoryCoding.network;

import java.util.ArrayList;
import java.util.List;

import org.jfree.data.xy.XYSeries;

public class Signal {
	public double[] data;
	public double start;
	public double end;
	private double timeShift;

	public Signal(double[] data) {
		this.data = data;
		start = 0;
		end = data.length - 1;
		timeShift = 0;
	}

	public Signal(double[] data, double start, double end) {
		this.data = data;
		this.start = start;
		this.end = end;
		timeShift = 0;
	}

	public Signal(double[] data, double start, double end, double timeShift) {
		this.data = data;
		this.start = start;
		this.end = Math.min(end, getLength() - 1 + timeShift);
		this.timeShift = timeShift;
	}

	/**
	 * Constructor that defensively copies the elements in signal
	 * 
	 * @param sig
	 */
	public Signal(Signal sig) {

		this.data = new double[sig.data.length];
		for (int i = 0; i < data.length; i++) {
			data[i] = sig.data[i];
		}
		this.start = sig.start;
		this.end = sig.end;
		this.timeShift = sig.timeShift;
	}
	// deprecated now
	// /**
	// * Returns the signal value at time
	// * @param t
	// * @return
	// */
	// public double getSignalValue(int t) {
	// if(t< start || t > end || isOutOfRange(t-timeShift)){
	// return 0;
	// }
	// else{
	// return data[t - timeShift];
	//
	// }
	// }

	public double getSignalValue(double t) {
		if (t < start || t > end || isOutOfRange(t - timeShift)) {
			return 0;
		}
		double arrayLocation = t - timeShift;
		int left = (int) arrayLocation;
		if (arrayLocation == left) {
			return data[left];
		}
		int right = left + 1;
		double leftValue = data[left];
		double rightValue = data[right];
		return leftValue + (rightValue - leftValue) * (arrayLocation - left);
		// else{
		// return data[t - timeShift];
		//
		// }
	}

	/**
	 * Check if the time is out of range for the signal
	 * 
	 * @param t
	 * @return
	 */
	private boolean isOutOfRange(double t) {
		if (t < 0 || t > (data.length - 1)) {
			return true;
		}
		return false;
	}

	/**
	 * This method returns the length of the signal
	 * 
	 * @return
	 */
	public int getLength() {
		return data.length;
	}

	/**
	 * This would get invoked when we are updating the timeshift
	 */
	public void setTimeShift(int timeShift) {
		this.timeShift = timeShift;
		this.end = Math.min(end, getLength() - 1 - timeShift);
	}

	/**
	 * Returns the time shift value
	 * 
	 * @return
	 */
	public double getTimeShift() {
		return timeShift;
	}

	/**
	 * This method returns the span of the signal during which the signal assumes
	 * non-zero values
	 * 
	 * @return
	 */
	public double getSignalSpan() {
		return end - start;
	}

	/**
	 * returns the start time of the signal
	 * 
	 * @return
	 */
	public double getStartTime() {
		return start;
	}

	/**
	 * returns the end time of the signal
	 * 
	 * @return
	 */
	public double getEndTime() {
		return end;
	}

	/**
	 * populates the signal data to be displayed on a line chart
	 */
	public XYSeries getSignalDisplayData(String signalTitle) {
		// DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		XYSeries series = new XYSeries(signalTitle);
		for (int i = (int) start; i <= end; i++) {
			series.add(i, getSignalValue(i));// addValue(getSignalValue(i), signalTitle, Integer.toString(i));
		}
		return series;
	}
	
	/**
	 * populates the signal data to be displayed on a line chart
	 */
	public XYSeries getSignalDisplayData(String signalTitle, int start, int end) {
		// DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		XYSeries series = new XYSeries(signalTitle);
		for (int i = start; i <= end; i++) {
			series.add(i, getSignalValue(i));// addValue(getSignalValue(i), signalTitle, Integer.toString(i));
		}
		return series;
	}

	/**
	 * This method will draw the given signal
	 */
	public void DrawSignal(String chartTitle) {
		List<Double> xValues = new ArrayList<>();
		List<Double> yValues = new ArrayList<>();
		for (int i = (int) start; i <= end; i++) {
			xValues.add((double) i);
			yValues.add(getSignalValue(i));
		}
		Utilities.drawLineChart(xValues, yValues, chartTitle, ConfigurationParameters.APPLICATION_TITLE,
				ConfigurationParameters.TIME_LABEL, ConfigurationParameters.VALUE_LABEL);
	}

	/**
	 * This method will draw the given signal
	 */
	public void DrawSignal(String chartTitle, int start, int end) {
		List<Double> xValues = new ArrayList<>();
		List<Double> yValues = new ArrayList<>();
		for (int i = start; i <= end; i++) {
			xValues.add((double) i);
			yValues.add(getSignalValue(i));
		}
		Utilities.drawLineChart(xValues, yValues, chartTitle, ConfigurationParameters.APPLICATION_TITLE,
				ConfigurationParameters.TIME_LABEL, ConfigurationParameters.VALUE_LABEL);
	}

	/**
	 * This method saves the plot of the signal into a file
	 *
	 * @param chartTitle
	 * @param fileName
	 */
	public void saveSignalIntoFile(String chartTitle, String fileName) {
		List<Double> xValues = new ArrayList<>();
		List<Double> yValues = new ArrayList<>();
		for (int i = (int) start; i <= end; i++) {
			xValues.add((double) i);
			yValues.add(getSignalValue(i));
		}
		Utilities.saveLineChartToFile(fileName, xValues, yValues, chartTitle, ConfigurationParameters.APPLICATION_TITLE,
				ConfigurationParameters.TIME_LABEL, ConfigurationParameters.VALUE_LABEL);
	}
}
