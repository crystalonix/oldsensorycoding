package sensoryCoding.network;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

import org.ejml.simple.SimpleBase;
import org.ejml.simple.SimpleMatrix;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class Utilities {
	public static void displaySetOfSignals(List<XYSeries> allSignalData, String title, String xLabel, String yLabel) {
		XYSeriesCollection signalData = new XYSeriesCollection();
		for (XYSeries xySeries : allSignalData) {
			signalData.addSeries(xySeries);
		}
		JFreeChart linechart = ChartFactory.createXYLineChart(title, xLabel, yLabel, signalData);
		ChartFrame chartFrame = new ChartFrame(title, linechart);
		chartFrame.setVisible(true);
	}

	public static void saveSetOfSignalsToFile(String fileName, List<XYSeries> allSignalData, String title, String xLabel, String yLabel) {
		XYSeriesCollection signalData = new XYSeriesCollection();
		for (XYSeries xySeries : allSignalData) {
			signalData.addSeries(xySeries);
		}
		JFreeChart linechart = ChartFactory.createXYLineChart(title, xLabel, yLabel, signalData);
		OutputStream out = null;
		try {
			out = new FileOutputStream(new File(fileName));
			ChartUtilities.writeChartAsPNG(out, linechart, ConfigurationParameters.IMAGE_WIDTH,
					ConfigurationParameters.IMAGE_HEIGHT);
		} catch (Exception e) {
			System.out.println("error in saving image to file");
		} finally {
			try {
				if (out != null)
					out.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
	}

	/**
	 * Returns a 2d array from a simplematrix
	 *
	 * @param matrix
	 * @return
	 */
	public static double[][] getArrayFromMatrix(SimpleMatrix matrix) {
		double[][] result = new double[matrix.numRows()][matrix.numCols()];
		for (int i = 0; i < matrix.numRows(); i++) {
			for (int j = 0; j < matrix.numCols(); j++) {
				result[i][j] = matrix.get(i, j);
			}
		}
		return result;
	}

	/**
	 * This method multiplies each element of the matrix by the scaling factor
	 *
	 * @param arrayFromMatrix
	 * @param scalingFactor
	 * @return
	 */
	public static double[][] scale(double[][] arrayFromMatrix, int scalingFactor) {
		double[][] result = new double[arrayFromMatrix.length][arrayFromMatrix[0].length];
		for (int i = 0; i < arrayFromMatrix.length; i++) {
			for (int j = 0; j < arrayFromMatrix[0].length; j++) {
				result[i][j] = (arrayFromMatrix[i][j]) * scalingFactor;
			}
		}
		return result;
	}

	/**
	 *
	 * @param data
	 * @return
	 */
	public static double[][] calculateInverse(double[][] data) {
		SimpleMatrix mat = new SimpleMatrix(data);
		SimpleBase base = mat.invert();
		double[][] invMat = new double[base.numRows()][base.numCols()];
		for (int i = 0; i < base.numRows(); i++)
			for (int j = 0; j < base.numRows(); j++) {
				invMat[i][j] = base.get(i, j);
			}
		return invMat;
	}

	/**
	 * Given a x-y data set this method draws the dataset as a line chart
	 *
	 * @param xValues
	 * @param yValues
	 * @param chartTitle
	 * @param applicationTitle
	 * @param xLabel
	 * @param yLabel
	 */
	public static void drawLineChart(List<Double> xValues, List<Double> yValues, String chartTitle,
			String applicationTitle, String xLabel, String yLabel) {
		JFreeChart linechart = ChartFactory.createXYLineChart(chartTitle, xLabel, yLabel,
				createDataSet(xValues, yValues, chartTitle));
		ChartFrame chartFrame = new ChartFrame(chartTitle, linechart);
		chartFrame.setVisible(true);
	}

	/**
	 * Given a x-y data set this method saves the graph into a png file
	 *
	 * @param fileName
	 * @param xValues
	 * @param yValues
	 * @param chartTitle
	 * @param applicationTitle
	 * @param xLabel
	 * @param yLabel
	 */
	public static void saveLineChartToFile(String fileName, List<Double> xValues, List<Double> yValues,
			String chartTitle, String applicationTitle, String xLabel, String yLabel) {
		OutputStream out = null;
		JFreeChart linechart = ChartFactory.createXYLineChart(chartTitle, xLabel, yLabel,
				createDataSet(xValues, yValues, chartTitle));
		try {
			out = new FileOutputStream(new File(fileName));
			ChartUtilities.writeChartAsPNG(out, linechart, ConfigurationParameters.IMAGE_WIDTH,
					ConfigurationParameters.IMAGE_HEIGHT);
		} catch (Exception e) {
			System.out.println("error in saving image to file");
		} finally {
			try {
				if (out != null)
					out.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
	}

	public static void saveLineChartToFile(String fileName, JFreeChart lineChart, Integer width, Integer height) {
		OutputStream out = null;
		try {
			out = new FileOutputStream(new File(fileName));
			ChartUtilities.writeChartAsPNG(out, lineChart, ConfigurationParameters.IMAGE_WIDTH,
					ConfigurationParameters.IMAGE_HEIGHT);
		} catch (Exception e) {
			System.out.println("error in saving image to file");
		} finally {
			try {
				if (out != null)
					out.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
	}

	private static XYDataset createDataSet(List<Double> xValues, List<Double> yValues, String chartTitle) {
		// TODO Auto-generated method stub
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries xySeries = new XYSeries(chartTitle);
		for (int i = 0; i < xValues.size(); i++)
			xySeries.add(xValues.get(i), yValues.get(i));
		dataset.addSeries(xySeries);
		return dataset;
	}

	private static XYDataset createDataSet(List<Double> xValues, List<Double> yValues) {
		// TODO Auto-generated method stub
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries xySeries = new XYSeries("potential at time");
		for (int i = 0; i < xValues.size(); i++)
			xySeries.add(xValues.get(i), yValues.get(i));
		dataset.addSeries(xySeries);
		return dataset;
	}

	/**
	 * This method writes matrix into a log file
	 *
	 * @param title
	 *            title for the log
	 * @param logFileName
	 *            name of the log file
	 * @param matrix
	 *            input matrix
	 */
	public static void writeMatrixToFile(String title, String logFileName, double[][] matrix) {
		FileWriter fw = null;
		BufferedWriter bw = null;
		try {
			StringBuilder sb = new StringBuilder();
			//sb.append(title).append("\n");
			sb.append(Integer.toString(matrix.length)).append("\n");
			sb.append(Integer.toString(matrix[0].length)).append("\n");
			for (int i = 0; i < matrix.length; i++) {
				for (int j = 0; j < matrix[0].length; j++) {
					sb.append(matrix[i][j]).append(", ");
				}
				sb.append("\n");
			}
			fw = new FileWriter(logFileName);
			bw = new BufferedWriter(fw);
			bw.write(sb.toString());
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (bw != null)
					bw.close();
				if (fw != null)
					fw.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
	}

	public static void outputMatrix(double[][] matrix) {
		if (matrix == null) {
			System.out.println("matrix is null and cannot be displayed");
			return;
		}
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				System.out.print(matrix[i][j] + ",");
			}
			System.out.println();
		}
	}

	public static void writeMatrixToBufferedWriter(BufferedWriter bw, double[][] kernelCoefficients)
			throws IOException {
		// TODO Auto-generated method stub
		int rows = kernelCoefficients.length;
		bw.write(Integer.toString(rows));
		bw.write("\n");
		int cols = kernelCoefficients[0].length;
		bw.write(Integer.toString(cols));
		bw.write("\n");
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				bw.write(Double.toString(kernelCoefficients[i][j]));
				bw.write("\n");
			}
		}
	}

	public static double[][] readMatrixToBufferedWriter(BufferedReader br) throws IOException {
		// TODO Auto-generated method stub
		int rows = Integer.parseInt(br.readLine());
		int cols = Integer.parseInt(br.readLine());
		double[][] kernelCoefficients = new double[rows][cols];
		for (int i = 0; i < rows; i++) {
			String rowNumbers = br.readLine();
			String []numbers = rowNumbers.split(",");
			for (int j = 0; j < cols; j++) {
				double value = Double.parseDouble(numbers[j]);
				kernelCoefficients[i][j] = value;
			}
		}
		return kernelCoefficients;
	}

	public static void writeListToBufferedWriter(BufferedWriter bw, List<Double> errorValues) throws IOException {
		// TODO Auto-generated method stub
		for (int i = 0; i < errorValues.size(); i++) {
			bw.write(Double.toString(errorValues.get(i)));
			bw.write("\n");

		}
	}
	public static void writeListToBufferedWriter(BufferedWriter bw, double[] errorValues) throws IOException {
		// TODO Auto-generated method stub
		for (int i = 0; i < errorValues.length; i++) {
			bw.write(Double.toString(errorValues[i]));
			bw.write("\n");

		}
	}
	public static void writeIntegerListToBufferedWriter(BufferedWriter bw, List<Integer> spikeCounts)
			throws IOException {
		// TODO Auto-generated method stub
		for (int i = 0; i < spikeCounts.size(); i++) {
			bw.write(Integer.toString(spikeCounts.get(i)));
			bw.write("\n");

		}
	}

	public static void writeIntegerListToBufferedWriter(BufferedWriter bw, int[] spikeCounts)
			throws IOException {
		// TODO Auto-generated method stub
		for (int i = 0; i < spikeCounts.length; i++) {
			bw.write(Integer.toString(spikeCounts[i]));
			bw.write("\n");
		}
	}

	public static List<Double> readListFromFile(BufferedReader br) throws IOException {
		String s = "";
		List<Double> values = new ArrayList<>();
		while ((s = br.readLine()) != null) {
			values.add(Double.parseDouble(s));
		}
		return values;
	}

	/**
	 * This method returns the moving average of a given list of values
	 * @param values
	 * @return
	 */
	public static List<Double> convertToMovingErrorAverage(List<Double> values, boolean skipNoReconstructions){
		List<Double> movingAvgs = new ArrayList<>();
		for(int i=0; i< values.size(); i++){
			int len = movingAvgs.size();
			if(values.get(i)<=0){
				continue;
			}
			if(values.get(i)>=1 && skipNoReconstructions){
				continue;
			}
			if(len==0){
				movingAvgs.add(values.get(i));
			}
			else{
				movingAvgs.add((movingAvgs.get(len-1)*len+values.get(i))/(len+1));
			}
		}
		return movingAvgs;
	}

	/**
	 * This method returns the moving average of a given list of values
	 * @param values
	 * @return
	 */
	public static List<Double> convertToMovingSpikeCountAverage(List<Double> values, boolean shouldSkipZeros){
		List<Double> movingAvgs = new ArrayList<>();
		for(int i=0; i< values.size(); i++){
			int len = movingAvgs.size();
			if(values.get(i)<=0/* && values.get(i)>=1*/ && shouldSkipZeros){
				continue;
			}
			if(len==0){
				movingAvgs.add(values.get(i));
			}
			else{
				movingAvgs.add((movingAvgs.get(len-1)*len+values.get(i))/(len+1));
			}
		}
		return movingAvgs;
	}
}
