package sensoryCoding.test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import sensoryCoding.network.Signal;
import sensoryCoding.network.SignalUtils;
import sensoryCoding.network.Utilities;

public class CollateResults {

	@Test
	public void testArray() throws IOException {
		String filePath = "/Users/anonymouschattopadhyay/Desktop/research/oldsensorycoding/SensoryCoding/Reports/output/";
		String tempWritePath = "/Users/anonymouschattopadhyay/Desktop/tempWrite/";
		double spRateLimit = -1;
		double spRateQuery = 25;
		int startIndex = 1;
		int endIndex = 400;
		List<Double> avgSnrs = new ArrayList<Double>();
		double signal_len = 500;
		int steps_interval = 2000;
		List<Integer> counts = new ArrayList<Integer>();
		for (int i = startIndex; i <= endIndex; i++) {
			String fileName = filePath + "report-" + i + ".txt";
			List<Double> errorRates = Utilities.readListByFileName(fileName);
			System.out.println("len of these files:" + errorRates.size());
			for (int j = errorRates.size() - 1; j >= 0; j--) {

				int index = errorRates.size() - j - 1;
				double thiserror = 0;
				if (errorRates.get(j) > 0) {
					thiserror = Utilities.getDBValue(errorRates.get(j));
				} else {
					thiserror = 50;	
				}
				if (avgSnrs.size() - 1 < index) {
					avgSnrs.add(thiserror);
					counts.add(1);
				} else {
					int n = counts.get(index);
					avgSnrs.set(index, (thiserror + n * avgSnrs.get(index)) / (n + 1));
					counts.set(index, n + 1);
				}
			}
		}
		double[] snrValues = new double[avgSnrs.size()];
		for (int i = 0; i < avgSnrs.size(); i++) {
			snrValues[i] = avgSnrs.get(i);
		}
		List<Double> spikeRates = new ArrayList<Double>();
		int counter = 0;

		for (int i = 10; i <= 1000 && counter < avgSnrs.size(); i += 10) {
			spikeRates.add(i * 44.1 / signal_len);
			counter++;
		}
		for (int i = 1050; i <= 5000 && counter < avgSnrs.size(); i += 50) {
			spikeRates.add(i * 44.1 / signal_len);
			counter++;
		}

		for (int i = 6000; counter < avgSnrs.size(); i += steps_interval) {
			spikeRates.add(i * 44.1 / signal_len);
			counter++;
		}
		
		int limitIndex =0;
		int queryIndex =0; 
		if(spRateLimit >0) {
			while(spikeRates.get(limitIndex)<=spRateLimit) {
				limitIndex++;
			}
		}
		
		if(spRateQuery >0) {
			while(spikeRates.get(queryIndex)<=spRateQuery) {
				queryIndex++;
			}
		}
		
		double [] rates = new double[endIndex- startIndex+1];
		for (int i = startIndex; i <= endIndex; i++) {
			String fileName = filePath + "report-" + i + ".txt";
			List<Double> errorRates = Utilities.readListByFileName(fileName);			
					double thiserror = Utilities.getDBValue(errorRates.get(errorRates.size()-queryIndex));
					rates[i-startIndex] = thiserror;
			}
		double avg =0;
		for (int i = startIndex; i <= endIndex; i++) {
			avg+= rates[i-startIndex];
		}
		System.out.println("avg:"+avg/(endIndex-startIndex+1));
		Utilities.writeListToFile(tempWritePath + "x-values.txt", spikeRates);
		Utilities.writeListToFile(tempWritePath + "y-values.txt", avgSnrs);
		if(limitIndex >0)
		{Utilities.drawLineChart(spikeRates.subList(0, limitIndex), avgSnrs.subList(0, limitIndex), "average SNR plot", "average SNR plot", "spike rate in kHz",
				"SNR in DB");}
		else {
			Utilities.drawLineChart(spikeRates, avgSnrs, "average SNR plot", "average SNR plot", "spike rate in kHz",
					"SNR in DB");
		}
		Utilities.writeListToFile(tempWritePath+ "snrValues-"+spRateQuery+".txt", rates);
		Utilities.writeListToFile(tempWritePath+ "avgSnrValues.txt", avgSnrs);
		System.out.println("Done");
	}
}
