package sensoryCoding.test;

import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.ejml.simple.SimpleMatrix;
import org.jfree.data.xy.XYSeries;
import org.junit.Assert;
import org.junit.Test;

import com.sun.corba.se.impl.javax.rmi.CORBA.Util;

import jdk.nashorn.internal.runtime.regexp.joni.Config;
import sensoryCoding.network.ConfigurationParameters;
import sensoryCoding.network.DataReader;
import sensoryCoding.network.DebugConfigurationParameters;
import sensoryCoding.network.KernelCalculator;
import sensoryCoding.network.KernelManager;
import sensoryCoding.network.Network;
import sensoryCoding.network.Signal;
import sensoryCoding.network.SignalUtils;
import sensoryCoding.network.Utilities;
import sun.nio.ch.Net;

public class NetworkTest {
	@Test
	public void checkTimeDifferential() throws Exception {
		for (int index = 0; index < 14; index++) {
			KernelManager kernelMgr = new KernelManager(ConfigurationParameters.numberOfKernels,
					ConfigurationParameters.numberofKernelComponents, ConfigurationParameters.lengthOfComponentSignals);
			int kernelindex = 0;
			int componentIndex = 0;
			double eps = 0.005;
			Signal kernelSignal = getKernelSignal();
			Signal blankSignal = null;
			Network nt = new Network(blankSignal);
			nt.init(kernelSignal);
			nt.calculateSpikeTimesAndReconstructSignal();
			nt.calculateErrorGradient();
			List<Integer> indexes = nt.kernelIndexesForSpikes;
			/*
			 * int index = 0; for (int i = 0; i < nt.kernelIndexesForSpikes.size(); i++) {
			 * if (nt.kernelIndexesForSpikes.get(i) == kernelindex) { index = i; break; } }
			 */
			double t0 = nt.spikeTimings.get(index);
			double timeDifferential = nt.timeDifferentials[index][componentIndex];

			double error1 = nt.calculateError();

			double actualError = nt.errorGradients[kernelindex][componentIndex] * eps;
			nt.kernelMgr.kernelCoefficients[kernelindex][componentIndex] += eps;
			nt.init(kernelSignal);
			nt.calculateSpikeTimesAndReconstructSignal();

			/*
			 * for (int i = 0; i < nt.kernelIndexesForSpikes.size(); i++) { if
			 * (nt.kernelIndexesForSpikes.get(i) == kernelindex) { index = i; break; } }
			 */

			double t0New = nt.spikeTimings.get(index);

			double error2 = nt.calculateError();
			System.out.println("real:" + (t0New - t0));
			System.out.println("ideal:" + timeDifferential * eps);
		}
		// nt.kernelMgr.dis
	}

	/*
	 * @Test public void checkCoeffDifferential() throws Exception { // TODO: next
	 * check for 16, 1 int rowIndex = 12; int colIndex = 24; int kernelindex = 6;
	 * int componentIndex = 0; double eps = 0.005;
	 *
	 * KernelManager kernelMgr = new
	 * KernelManager(ConfigurationParameters.numberOfKernels,
	 * ConfigurationParameters.numberofKernelComponents,
	 * ConfigurationParameters.lengthOfComponentSignals); Signal kernelSignal =
	 * getKernelSignal(); Network nt = new Network(kernelSignal);
	 *//************************/
	/*
	*//******* first run *******/
	/*
	*//***********************/
	/*
	 * nt.calculateSpikeTimesAndReconstructSignal(); nt.calculateErrorGradient(); //
	 * System.out.println("spike time12:"+nt.spikeTimings.get(12)+" spike // time //
	 * 28:"+nt.spikeTimings.get(28)+"timediff:"+(nt.spikeTimings.get(12)-nt.
	 * spikeTimings.get(28))); int rowSpikeIndex =
	 * nt.kernelIndexesForSpikes.get(rowIndex); double rowSpikeTime =
	 * nt.spikeTimings.get(rowIndex); int colSpikeIndex =
	 * nt.kernelIndexesForSpikes.get(colIndex); double colSpikeTime =
	 * nt.spikeTimings.get(colIndex);
	 *
	 * System.out.println("check out spike time: index, kernel index, time::" +
	 * rowIndex + "," + rowSpikeIndex + "," + rowSpikeTime + "," + colIndex + "," +
	 * colSpikeIndex + "," + colSpikeTime); // System.out.println(); double[][]
	 * pdiffs = nt.pdiffs[kernelindex][componentIndex]; double[][] pmatOld =
	 * nt.pMatrix; double t0Old = nt.spikeTimings.get(rowIndex); double timeDiff =
	 * nt.timeDifferentials[rowIndex][componentIndex];
	 *
	 * System.out.println("see the p value:" + pmatOld[rowIndex][colIndex]);
	 *
	 * int alpha1 = (int) kernelMgr.getBasisForKernelComp(rowSpikeIndex,
	 * 0).getEndTime() / 3; int alpha2 = (int)
	 * kernelMgr.getBasisForKernelComp(colSpikeIndex, 0).getEndTime() / 3; double
	 * delta = nt.spikeTimings.get(colIndex) - nt.spikeTimings.get(rowIndex);
	 * System.out.println("ideal p value:" +
	 * KernelCalculator.kernelMultiplication(alpha1, alpha2, delta,
	 * ConfigurationParameters.numberofKernelComponents, null, null));
	 *//************************/
	/*
	*//** ideal p differential **/
	/*
	*//************************/
	/*
	 * double idealPdiff = KernelCalculator.kernelByBSpline(alpha2, alpha1, -delta,
	 * ConfigurationParameters.numberofKernelComponents, null) - timeDiff *
	 * KernelCalculator.kernelDifferntialMultiplication(alpha2, alpha1, -delta,
	 * ConfigurationParameters.numberofKernelComponents, null, null);
	 * System.out.println("ideal p diff*eps:" + idealPdiff * eps); //
	 * System.out.println("ideal p //
	 * diff:"+(KernelCalculator.kernelByBSpline(alpha1, alpha2, delta, //
	 * ConfigurationParameters.numberofKernelComponents, null)- //
	 * timeDiff*KernelCalculator.kernelDifferntialMultiplication(alpha1, // alpha2,
	 * delta,ConfigurationParameters.numberofKernelComponents, null, // null)));
	 *
	 * // nt.displayPByPInv();
	 *
	 * // int index=0; // for(int i=0; i<nt.kernelIndexesForSpikes.size(); i++) { //
	 * if(nt.kernelIndexesForSpikes.get(i)==kernelindex) { // rowIndex = i; //
	 * break; // } // } // colIndex = rowIndex;
	 *
	 * // Signal sg = nt.kernelMgr.getBasisForKernelComp(6, 1); //
	 * nt.kernelMgr.componentByComponentSignalConvolutionCache[6][6].DrawSignal(
	 * ""); Signal compBYDiffComp =
	 * nt.kernelMgr.componentByDiffComponentSignalConvolutionCache[2][6];
	 * compBYDiffComp.DrawSignal("2-6 comp by diff comp");
	 * nt.kernelMgr.getKernel(2).DrawSignal("kernelSignal"); double[] cOld =
	 * nt.coefficientsOfReconstructedSignal; double[] cDiffs =
	 * nt.coefficientDifferentials[kernelindex][componentIndex]; double cDiff =
	 * nt.coefficientDifferentials[kernelindex][componentIndex][rowIndex]; double c0
	 * = nt.coefficientsOfReconstructedSignal[rowIndex]; double p0 =
	 * pmatOld[rowIndex][colIndex]; double pdiff = pdiffs[rowIndex][colIndex];
	 * System.out.println("p diff real:" + pdiff); System.out.println("pvalue:" +
	 * pmatOld[rowIndex][colIndex]);
	 *
	 *//************************/
	/*
	*//******* second run ******/

	/*
	*//***********************//*
								 * nt.init(kernelSignal); nt.kernelMgr.kernelCoefficients[kernelindex][
								 * componentIndex] += eps; nt.kernelMgr.loadKernelCache();
								 * nt.calculateSpikeTimesAndReconstructSignal();
								 * System.out.println("spike time12 new:" + nt.spikeTimings.get(12) +
								 * " spike time 28:" + nt.spikeTimings.get(28) + " diff:" +
								 * (nt.spikeTimings.get(12) - nt.spikeTimings.get(28))); double t0New =
								 * nt.spikeTimings.get(rowIndex); // nt.displayPByPInv(); rowSpikeIndex =
								 * nt.kernelIndexesForSpikes.get(rowIndex); rowSpikeTime =
								 * nt.spikeTimings.get(rowIndex); colSpikeIndex =
								 * nt.kernelIndexesForSpikes.get(colIndex); colSpikeTime =
								 * nt.spikeTimings.get(colIndex); System.out.
								 * println("check out spike time: index, kernel index, time::" + rowIndex + ","
								 * + rowSpikeIndex + "," + rowSpikeTime + "," + colIndex + "," + colSpikeIndex +
								 * "," + colSpikeTime);
								 *
								 * double[][] pmatNew = nt.pMatrix; double c0New =
								 * nt.coefficientsOfReconstructedSignal[rowIndex ]; double p0New =
								 * pmatNew[rowIndex][colIndex];
								 *
								 * alpha1 = (int) kernelMgr.getBasisForKernelComp( rowSpikeIndex,
								 * 0).getEndTime() / 3; alpha2 = (int) kernelMgr.getBasisForKernelComp(
								 * colSpikeIndex, 0).getEndTime() / 3; delta = nt.spikeTimings.get(colIndex) -
								 * nt.spikeTimings.get(rowIndex); System.out.println("ideal p value:" +
								 * KernelCalculator.kernelMultiplication(alpha1, alpha2, delta,
								 * ConfigurationParameters. numberofKernelComponents, null, null));
								 *
								 * // check the pmatrix errors for (int i = 0; i < nt.spikeTimings.size(); i++)
								 * { for (int j = 0; j < nt.spikeTimings.size(); j++) { if (pdiffs[i][j] != 0) {
								 * double errorRate = (pmatNew[i][j] - pmatOld[i][j]) / ((pdiffs[i][j]) * eps);
								 * System.out.println( "real p value:" + pmatNew[i][j] + "Accuracy:" + errorRate
								 * + " for " + i + ", " + j); } } } // check coeff errors double[] cNew =
								 * nt.coefficientsOfReconstructedSignal; for (int i = 0; i < cOld.length; i++) {
								 * System.out.println("value:" + cNew[i] + "accuracy:" + (cNew[i] - cOld[i]) /
								 * (cDiffs[i] * eps)); } System.out.println("ideal spike time diff:" + timeDiff
								 * * eps); System.out.println("real spike time diff:" + (t0New - t0Old));
								 * System.out.println("see the p value now:" + pmatNew[rowIndex][colIndex]);
								 *
								 * alpha1 = (int) kernelMgr.getBasisForKernelComp( rowSpikeIndex,
								 * 0).getEndTime() / 3; alpha2 = (int) kernelMgr.getBasisForKernelComp(
								 * colSpikeIndex, 0).getEndTime() / 3; delta = nt.spikeTimings.get(colIndex) -
								 * nt.spikeTimings.get(rowIndex); System.out.println("ideal p value now:" +
								 * KernelCalculator.kernelMultiplication(alpha1, alpha2, delta,
								 * ConfigurationParameters. numberofKernelComponents, null, null));
								 *
								 * System.out.println("real:" + (c0New - c0)); System.out.println("ideal:" +
								 * cDiff * eps); System.out.println("real pdiff:" + (p0New - p0));
								 * System.out.println("ideal pdiff:" + pdiff * eps); // nt.kernelMgr.dis }
								 */

	@Test
	public void checkErrorDifferential() throws Exception {
		// TODO: next check for 16, 1
		/*
		 * int rowIndex = 12; int colIndex = 24;
		 */
		Signal blankSignal = null;
		Network nt = new Network(blankSignal);
		Signal kernelSignal = getKernelSignal();
		kernelSignal.DrawSignal("original Signal");
		double eps = 0.00005;
		for (int kernelindex = 0; kernelindex < ConfigurationParameters.numberOfKernels; kernelindex++) {
			for (int componentIndex = 0; componentIndex < ConfigurationParameters.numberofKernelComponents; componentIndex++) {

				System.out.println("you are seeing kernel index:" + kernelindex + ",component index:" + componentIndex);

				nt.init(kernelSignal);
				/************************/
				/******* first run *******/
				/***********************/
				nt.calculateSpikeTimesAndReconstructSignal();
				nt.calculateErrorGradient();
				nt.getReconstructedSignal().DrawSignal("recon sig old");
				// double oldErrorGradients[][] = nt.errorGradients;
				double initialError = nt.calculateError();
				double idealErrorRate = nt.errorGradients[kernelindex][componentIndex];

				/************************/
				/******* second run ******/
				/***********************/
				nt.init(kernelSignal);
				nt.kernelMgr.kernelCoefficients[kernelindex][componentIndex] += eps;
				nt.kernelMgr.loadKernelCache();
				nt.calculateSpikeTimesAndReconstructSignal();
				double newError = nt.calculateError();
				double actualError = (newError - initialError);
				nt.getReconstructedSignal().DrawSignal("recon sig new");
				double idealError = idealErrorRate * eps;
				double errorRate = 100 * (actualError - idealError) / idealError;
				System.out.println("(" + actualError + "," + idealError + "," + initialError + "," + errorRate + "%)");
				/*
				 * System.out.println("real error diff:" + actualError);
				 * System.out.println("ideal error:" + idealError);
				 */ }
		}
	}

	/*
	 * @Test public void checkPDifferentials() throws Exception { double eps =
	 * 0.005; for (int kernelIndex = 0; kernelIndex <
	 * ConfigurationParameters.numberOfKernels; kernelIndex++) { for (int compIndex
	 * = 0; compIndex < ConfigurationParameters.numberofKernelComponents;
	 * compIndex++) { System.out.println("you are looking for kernelIndex:" +
	 * kernelIndex + " component index:" + compIndex); KernelManager kernelMgr = new
	 * KernelManager(ConfigurationParameters.numberOfKernels,
	 * ConfigurationParameters.numberofKernelComponents,
	 * ConfigurationParameters.lengthOfComponentSignals); Signal kernelSignal =
	 * getKernelSignal(); Network nt = new Network(kernelSignal);
	 *//************************/
	/*
	*//******* first run ******/
	/*
	*//***********************/
	/*
	 * nt.calculateSpikeTimesAndReconstructSignal(); nt.calculateErrorGradient();
	 * double[][] pOld = nt.pMatrix; double pdiffs[][] =
	 * nt.pdiffs[kernelIndex][compIndex];
	 *//************************/
	/*
	*//******* second run ******/

	/*
	*//***********************//*
								 * nt.init(kernelSignal); nt.kernelMgr.kernelCoefficients[kernelIndex][
								 * compIndex] += eps; nt.kernelMgr.loadKernelCache();
								 * nt.calculateSpikeTimesAndReconstructSignal(); double[][] pNew = nt.pMatrix;
								 * for (int i = 0; i < nt.spikeTimings.size(); i++) { for (int j = 0; j <
								 * nt.spikeTimings.size(); j++) { if ((pNew[i][j] == pOld[i][j]) ||
								 * Math.abs((pNew[i][j] - pOld[i][j]) / (pdiffs[i][j] * eps) - 1) < 0.1) { } //
								 * System.out.println("0"); else { double errorFrac = Math.abs((pNew[i][j] -
								 * pOld[i][j]) / (pdiffs[i][j] * eps) - 1) * 100; System.out.println("entry" + i
								 * + "," + j + ":" + "(" + pdiffs[i][j] * eps + "," + (pNew[i][j] - pOld[i][j])
								 * + ") error rate:" + errorFrac + "%"); }
								 *
								 * System.out.println("entry number:"+i+", "+j);
								 * System.out.println("ideal pdiff:"+pdiffs[i][j ]*eps);
								 * System.out.println("real pdiff:"+(pNew[i][j]- pOld[i][j]));
								 * System.out.println( "+++++++++++++++++++++++++++++++++"); }
								 * System.out.println(); } System.out.println(
								 * "======================================"); } } }
								 */

	@Test
	public void sensitiveDependenceCheckOnPMatrix() throws Exception {
		double eps = 0.1;
		KernelManager kernelMgr = new KernelManager(ConfigurationParameters.numberOfKernels,
				ConfigurationParameters.numberofKernelComponents, ConfigurationParameters.lengthOfComponentSignals);
		Signal kernelSignal = getKernelSignal();
		Network nt = new Network(kernelSignal);
		/************************/
		/******* first run *****/
		/***********************/
		nt.calculateSpikeTimesAndReconstructSignal();
		nt.calculateErrorGradient();
		double[][] pmat = nt.pMatrix;
		double[] oldCoeffs = nt.coefficientsOfReconstructedSignal;

		for (int i = 0; i < pmat.length; i++) {
			for (int j = 0; j < pmat.length; j++) {
				double del = (pmat[i][j] * eps) * (Math.random() - 0.5);
				System.out.println(del);
				pmat[i][j] += del;
			}
		}
		SimpleMatrix pInvMatrix = new SimpleMatrix(nt.pMatrix).invert();
		double[] newCoeffs = Utilities.getArrayFromMatrix((pInvMatrix.mult(nt.convolvedValuesMatrix)).transpose())[0];
		for (int i = 0; i < newCoeffs.length; i++) {
			System.out.println("error rate for " + i + "th coeff is:"
					+ ((newCoeffs[i] - oldCoeffs[i]) / oldCoeffs[i]) * 100 + "% ," + oldCoeffs[i]);
			// System.out.print("("+kernelindexesOld.get(i)+","+spiketimesOld.get(i)+",
			// "+nt.coefficientsOfReconstructedSignal[i]+") ");
		}
		System.out.println();
		// Utilities.outputMatrix(matrix);
		// Utilities.outputMatrix(Utilities.getArrayFromMatrix(mat));
		System.out.println("++++++++++++++++++++++++++++++++++++++++++");
	}
	/*
	 * @Test public void checkCoeffcientErrors() throws Exception { double eps =
	 * 0.005;
	 *
	 * int kernelIndex = 0; int compIndex = 0;
	 *
	 * for (int kernelIndex = 0; kernelIndex <
	 * ConfigurationParameters.numberOfKernels; kernelIndex++) { for (int compIndex
	 * = 0; compIndex < ConfigurationParameters.numberofKernelComponents;
	 * compIndex++) { System.out.println("seeing kernelindex:" + kernelIndex +
	 * " comp index:" + compIndex); KernelManager kernelMgr = new
	 * KernelManager(ConfigurationParameters.numberOfKernels,
	 * ConfigurationParameters.numberofKernelComponents,
	 * ConfigurationParameters.lengthOfComponentSignals); Signal kernelSignal =
	 * getKernelSignal(); Network nt = new Network(kernelSignal);
	 *
	 *//************************/
	/*
	*//******* first run *****/
	/*
	*//***********************/
	/*
	 * nt.calculateSpikeTimesAndReconstructSignal(); nt.calculateErrorGradient();
	 * double[] oldCoeffs = nt.coefficientsOfReconstructedSignal; double pdiffs[][]
	 * = nt.pdiffs[kernelIndex][compIndex]; double[][] alphaArray = new
	 * double[1][nt.coefficientsOfReconstructedSignal.length]; alphaArray[0] =
	 * nt.coefficientsOfReconstructedSignal; SimpleMatrix alphaColumn = (new
	 * SimpleMatrix(alphaArray)).transpose(); SimpleMatrix pdiffByAlpha = new
	 * SimpleMatrix(pdiffs).mult(alphaColumn); double coeffDiffMat[][] =
	 * Utilities.getArrayFromMatrix(nt.pInvMatrix.mult(pdiffByAlpha).transpose() );
	 * double idealCoeffDiffs[] = (Utilities.scale(coeffDiffMat, -1))[0];
	 *
	 *//************************/
	/*
	*//******* second run ******/

	/*
	*//***********************//*
								 * nt.init(kernelSignal); nt.kernelMgr.kernelCoefficients[kernelIndex][
								 * compIndex] += eps; nt.kernelMgr.loadKernelCache();
								 * nt.calculateSpikeTimesAndReconstructSignal(); double newCoeffs[] =
								 * nt.coefficientsOfReconstructedSignal; for (int i = 0; i < newCoeffs.length;
								 * i++) { if ((newCoeffs[i] == oldCoeffs[i]) || ((((newCoeffs[i] - oldCoeffs[i])
								 * / (idealCoeffDiffs[i] * eps) - 1) < 0.1) && (((newCoeffs[i] - oldCoeffs[i]) /
								 * (idealCoeffDiffs[i] * eps) - 1) > 0))) { } // System.out.println("0"); else {
								 * double errorFrac = ((newCoeffs[i] - oldCoeffs[i]) / (idealCoeffDiffs[i] *
								 * eps) - 1) * 100; System.out.println("entry" + i + ":" + oldCoeffs[i] + " (" +
								 * idealCoeffDiffs[i] * eps + "," + (newCoeffs[i] - oldCoeffs[i]) +
								 * ") error rate:" + errorFrac + "%"); } } } } }
								 */
	@Test
	public void checkSpikeTimeKernelIndex() throws Exception {
		int kernelIndex = 0;
		int compIndex = 0;
		double eps = 0.005;
		KernelManager kernelMgr = new KernelManager(ConfigurationParameters.numberOfKernels,
				ConfigurationParameters.numberofKernelComponents, ConfigurationParameters.lengthOfComponentSignals);
		Signal kernelSignal = getKernelSignal();
		Network nt = new Network(kernelSignal);
		/************************/
		/******* first run *****/
		/***********************/
		nt.calculateSpikeTimesAndReconstructSignal();
		nt.calculateErrorGradient();
		double[][] pmat = nt.pMatrix;
		double thresholds[] = new double[nt.pMatrix.length];

		for (int i = 0; i < pmat.length; i++) {
			thresholds[i] = ConfigurationParameters.initialThresHoldValue;
		}
		// double [] calculatedCoeffs =
		Utilities.outputMatrix(Utilities.getArrayFromMatrix((new SimpleMatrix(pmat)).invert()));
		System.out.println("++++++++++++++++++++++");
		Utilities.outputMatrix(pmat);
		List<Integer> kernelindexesOld = nt.kernelIndexesForSpikes;
		List<Double> spiketimesOld = nt.spikeTimings;
		// System.out.println("old
		// sizes:("+kernelindexesOld.size()+","+spiketimesOld.size()+") ");

		for (int i = 0; i < kernelindexesOld.size(); i++) {
			System.out.print("(" + kernelindexesOld.get(i) + "," + spiketimesOld.get(i) + ", "
					+ nt.coefficientsOfReconstructedSignal[i] + ")  ");
		}
		System.out.println();
		// Utilities.outputMatrix(matrix);
		// Utilities.outputMatrix(Utilities.getArrayFromMatrix(mat));
		System.out.println("++++++++++++++++++++++++++++++++++++++++++");

		/************************/
		/******* second run ******/
		/***********************/
		nt.init(kernelSignal);
		nt.kernelMgr.kernelCoefficients[kernelIndex][compIndex] += eps;
		nt.kernelMgr.loadKernelCache();
		nt.calculateSpikeTimesAndReconstructSignal();
		List<Integer> kernelindexesNew = nt.kernelIndexesForSpikes;
		List<Double> spiketimesNew = nt.spikeTimings;
		// System.out.println("new
		// sizes:("+kernelindexesNew.size()+","+spiketimesNew.size()+") ");
		for (int i = 0; i < kernelindexesNew.size(); i++) {
			System.out.print("(" + kernelindexesNew.get(i) + "," + spiketimesNew.get(i) + ","
					+ nt.coefficientsOfReconstructedSignal[i] + ")  ");
		}
		System.out.println();
		// Utilities.outputMatrix(pmat);
	}

	@Test
	public void pinvTest() throws Exception {
		int kernelIndex = 0;
		int compIndex = 0;
		double eps = 0.005;
		KernelManager kernelMgr = new KernelManager(ConfigurationParameters.numberOfKernels,
				ConfigurationParameters.numberofKernelComponents, ConfigurationParameters.lengthOfComponentSignals);
		Signal kernelSignal = getKernelSignal();
		Network nt = new Network(kernelSignal);
		/************************/
		/******* first run *****/
		/***********************/
		nt.calculateSpikeTimesAndReconstructSignal();
		nt.calculateErrorGradient();
		SimpleMatrix mat = (nt.pInvMatrix).mult(new SimpleMatrix(nt.pMatrix));
		double oldInv[][] = Utilities.getArrayFromMatrix(nt.pInvMatrix);
		Utilities.outputMatrix(Utilities.getArrayFromMatrix(nt.pInvMatrix));
		// Utilities.outputMatrix(matrix);
		// Utilities.outputMatrix(Utilities.getArrayFromMatrix(mat));
		System.out.println("++++++++++++++++++++++++++++++++++++++++++");

		/************************/
		/******* second run ******/
		/***********************/
		nt.init(kernelSignal);
		nt.kernelMgr.kernelCoefficients[kernelIndex][compIndex] += eps;
		nt.kernelMgr.loadKernelCache();
		nt.calculateSpikeTimesAndReconstructSignal();
		mat = (nt.pInvMatrix).mult(new SimpleMatrix(nt.pMatrix));
		Utilities.outputMatrix(Utilities.getArrayFromMatrix(nt.pInvMatrix));
		double newInv[][] = Utilities.getArrayFromMatrix(nt.pInvMatrix);
		System.out.println("++++++++++++++++++++++++++++++++++++++++++");
		for (int i = 0; i < newInv.length; i++) {
			for (int j = 0; j < newInv[0].length; j++) {
				if ((newInv[i][j] == 0 && oldInv[i][j] == 0)
						|| Math.abs((newInv[i][j] - oldInv[i][j]) / newInv[i][j]) < 0.1) {
				} // System.out.println("0");
				else {
					System.out.println(
							"entry" + i + "," + j + ":" + Math.abs((newInv[i][j] - oldInv[i][j]) / newInv[i][j]));
				}
			}
		}
	}

	/*
	 * @Test public void checkAlphas() throws Exception { // TODO: next check for
	 * 16, 1
	 *
	 * int rowIndex = 12; int colIndex = 24;
	 *
	 * KernelManager kernelMgr = new
	 * KernelManager(ConfigurationParameters.numberOfKernels,
	 * ConfigurationParameters.numberofKernelComponents,
	 * ConfigurationParameters.lengthOfComponentSignals); Signal kernelSignal =
	 * getKernelSignal(); Network nt = new Network(kernelSignal);
	 *//************************/
	/*
	*//******* first run *******/

	/*
	*//***********************//*
								 * nt.calculateSpikeTimesAndReconstructSignal(); nt.calculateErrorGradient();
								 * for (int i = 0; i < nt.spikeTimings.size(); i++) { int ji =
								 * nt.kernelIndexesForSpikes.get(i); double ti = nt.spikeTimings.get(i); double
								 * value = 0; for (int k = 0; k < nt.spikeTimings.size(); k++) { int jk =
								 * nt.kernelIndexesForSpikes.get(k); int alpha1 = (int)
								 * kernelMgr.getBasisForKernelComp(jk, 0).getEndTime() / 3; int alpha2 = (int)
								 * kernelMgr.getBasisForKernelComp(ji, 0).getEndTime() / 3; double tk =
								 * nt.spikeTimings.get(k); value += nt.coefficientsOfReconstructedSignal[k]
								 * KernelCalculator.kernelMultiplication(alpha1, alpha2, ti - tk, 4, null,
								 * null); } System.out.println("multiplied value:" + value); }
								 *
								 * }
								 */

	@Test
	public void ahpTest() throws Exception {
		List<Integer> selectedKernels = new ArrayList<>();
		selectedKernels.add(0);
		Network net = new Network(selectedKernels);
		Signal input = SignalUtils.shiftSignal(net.kernelMgr.getKernel(0), 100);
		input = SignalUtils.addTwoSignals(input, SignalUtils.shiftSignal(input, 800));
		net.init(input);
		net.calculateSpikeTimesAndReconstructSignal();
		Signal reconsSignal = net.getReconstructedSignal();
		input.DrawSignal("input signal");
		reconsSignal.DrawSignal("reconstructed signal");
		Signal thresholdSig = new Signal(net.thresholdForTesting, 50, 1000, 0);
		thresholdSig.DrawSignal("thresholds");
		System.out.println("Ahp Function tested");
	}

	@Test
	public void checkKernels() throws Exception {
		Signal blankSig = null;
		Network net = new Network(blankSig);
		net.kernelMgr.displayKernels();
		System.out.println("kernels are checked");
	}

	@Test
	public void checkErrorRateOnASingleOne() throws Exception {
		// List<Double> errors = new ArrayList<>();
		List<Integer> selectedKernels = new ArrayList<>();
		selectedKernels.add(0);
		Network net = new Network(selectedKernels);
		Signal input = SignalUtils.shiftSignal(net.kernelMgr.getKernel(0), 100);
		input = SignalUtils.addTwoSignals(input, SignalUtils.shiftSignal(input, 800));
		for (int i = 0; i < 100; i++) {
			net.init(input);
			net.calculateSpikeTimesAndReconstructSignal();
			net.calculateErrorGradient();
			net.updateKernelCoefficients();
		}
	}

	@Test
	public void testPartitionedSignals() throws Exception {
		DataReader dtr = new DataReader(null);
		dtr.loadPartitionIntoMemory(10, 1);
		ConfigurationParameters.numberOfKernels = 20;
		ConfigurationParameters.numberofKernelComponents = 4;
		List<Integer> selectedKernels = new ArrayList<>();
		selectedKernels.add(0);
		selectedKernels.add(5);
		Signal blankSignal = null;
		Network net = new Network(blankSignal);
		net.kernelMgr.displayKernels();
		int k = 0;
		Signal inputSignal = dtr.allSignalPieces[k];
		inputSignal.saveSignalIntoFile("input signal",
				ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "inputSignal.png");
		net.init(dtr.allSignalPieces[0]);
		for (int i = 0; i < 10000; i++) {
			long t0 = System.currentTimeMillis();
			net.init(dtr.allSignalPieces[k]);
			System.out.println("time for update:" + (System.currentTimeMillis() - t0));
			long t1 = System.currentTimeMillis();
			net.calculateSpikeTimesAndReconstructSignal();
			long t2 = System.currentTimeMillis();
			net.calculateErrorGradient();
			long t3 = System.currentTimeMillis();
			System.out.println("Recons time:" + (t2 - t1));
			System.out.println("error gard time:" + (t3 - t2));
			double errorrate = net.updateKernelCoefficients();
			long t4 = System.currentTimeMillis();
			System.out.println("update time:" + (t4 - t3));
			System.out.println(i + " steps executed");
			if (i == 0) {
				net.getReconstructedSignal().saveSignalIntoFile("reconsturcted signal after" + i + " steps",
						ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "recons-" + i + ".png");
			}
			if (i % 50 == 0) {
				net.getReconstructedSignal().saveSignalIntoFile("reconsturcted signal after" + i + " steps",
						ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "recons-" + i + "-" + errorrate + ".png");
			}
		}
		net.getReconstructedSignal().DrawSignal("final reconstruction");
		System.out.println("this is done");
	}

	@Test
	public void trainOnDataSet() throws Exception {
		for (int i = 0; i < ConfigurationParameters.numberOfPartitionChunks; i++) {
			DataReader dtr = new DataReader(null);
			dtr.loadPartitionIntoMemory(10, i);
			ConfigurationParameters.numberOfKernels = 20;
			ConfigurationParameters.numberofKernelComponents = 4;
			Signal blankSignal = null;
			Network net = new Network(blankSignal);
			for (int j = 0; j < dtr.allSignalPieces.length; j++) {
				net.init(dtr.allSignalPieces[j]);
				net.calculateSpikeTimesAndReconstructSignal();
				net.calculateErrorGradient();
				double errorrate = net.updateKernelCoefficients();
				if (j % 50 == 0) {
					dtr.allSignalPieces[j].saveSignalIntoFile("input signal",
							ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "input-" + j + "-" + errorrate + ".png");
					net.getReconstructedSignal().saveSignalIntoFile("reconsturcted signal after" + j + " steps",
							ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "recons-" + j + "-" + errorrate + ".png");
				}

			}
		}
		System.out.println("this is done");
	}

	private static Signal getKernelSignal() throws Exception {
		KernelManager krnmgr = new KernelManager(ConfigurationParameters.numberOfKernels,
				ConfigurationParameters.numberofKernelComponents, ConfigurationParameters.lengthOfComponentSignals);
		// return SignalUtils.shiftSignal(krnmgr.getInvertedKernel(0),
		// 1000/ConfigurationParameters.TIME_STEP);
		return SignalUtils.addTwoSignals(
				SignalUtils.shiftSignal(krnmgr.getKernel(0), 1000 / ConfigurationParameters.TIME_STEP),
				SignalUtils.shiftSignal(krnmgr.getKernel(0), 10 / ConfigurationParameters.TIME_STEP));
	}

	@Test
	public void singleReconstructionImprovement() throws Exception {
		/************************/
		/**** set up for this ****/
		/******** experiment ******/
		ConfigurationParameters.NATURE_OF_KERNEL_SPREAD = 1;
		ConfigurationParameters.numberOfKernels = 20;
		ConfigurationParameters.numberofKernelComponents = 4;

		ConfigurationParameters.initialCoefficientValue = 10;
		ConfigurationParameters.initialThresHoldValue = 6;
		ConfigurationParameters.INITIAL_LEARNING_RATE = 0.01;
		ConfigurationParameters.INITIAL_THRESHOLD_GRADIENT_STEP_LENGTH = 0.01;
		ConfigurationParameters.NUMBER_OF_KICKS = 100;
		/************************/
		/************************/
		/************************/
		Signal nullSignal = null;
		Network net = new Network(nullSignal);
		Signal input = SignalUtils.addTwoSignals(
				SignalUtils.addTwoSignals(SignalUtils.shiftSignal(net.kernelMgr.getInvertedKernel(5), 500),
						SignalUtils.shiftSignal(net.kernelMgr.getInvertedKernel(0), 700)),
				SignalUtils.shiftSignal(net.kernelMgr.getInvertedKernel(3), 800));
		input.DrawSignal("input signal");
		net.init(input);
		String logFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "initialKernel.txt";
		Utilities.writeMatrixToFile("inital kernelCoeffs", logFileName, net.kernelMgr.kernelCoefficients);
		XYSeries inputDisplay = input.getSignalDisplayData("input signal");

		double eps = 0.05;
		/************** check for the gradients ***********/
		/************************************************/
		int kernelIndex = 2;
		int compIndex = 0;
		net.calculateSpikeTimesAndReconstructSignal();
		net.calculateErrorGradient();
		double error1 = net.calculateError();
		double errorideal = net.errorGradients[kernelIndex][compIndex] * eps;
		net.kernelMgr.kernelCoefficients[kernelIndex][compIndex] += eps;

		net.init(input);
		net.calculateSpikeTimesAndReconstructSignal();
		net.calculateErrorGradient();
		double error2 = net.calculateError();
		double errorActual = error2 - error1;
		/************************************************/
		/************************************************/

		for (int i = 0; i < 1000000; i++) {
			net.init(input, false);
			net.calculateSpikeTimesAndReconstructSignal();
			net.calculateErrorGradient();
			net.updateKernelCoefficients();
			if (i % 50 == 0) {
				String logKernelFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "initialKernel-" + i
						+ ".txt";
				Utilities.writeMatrixToFile("inital kernelCoeffs", logKernelFileName, net.kernelMgr.kernelCoefficients);
				net.init(input, false);
				net.calculateSpikeTimesAndReconstructSignal();
				List<XYSeries> signals = new ArrayList<>();
				String fileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "reconstructedSignal-" + i
						+ "-spikeCount-" + net.spikeTimings.size() + ".png";
				XYSeries reconsSignalData = net.getReconstructedSignal()
						.getSignalDisplayData("reconstructed signal after" + i + " steps");
				signals.add(inputDisplay);
				signals.add(reconsSignalData);
				Utilities.saveSetOfSignalsToFile(fileName, signals, "reconstruction after" + i + " steps", "time in ms",
						"signal value");
			}
		}
		System.out.println("done here");
	}

	@Test
	public void testSignalConvs() throws Exception {
		Signal input = getKernelSignal();
		input.DrawSignal("input");
		Signal nullSig = null;
		Network net = new Network(nullSig);
		net.init(input);
		List<XYSeries> bBybs = new ArrayList<>();
		List<XYSeries> bByDbs = new ArrayList<>();
		List<XYSeries> signalByKernels = new ArrayList<>();
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			bBybs.add(net.kernelMgr.componentByComponentSignalConvolutionCache[0][i]
					.getSignalDisplayData("kernel by kernel comp:" + i));
			bByDbs.add(net.kernelMgr.componentByDiffComponentSignalConvolutionCache[0][i]
					.getSignalDisplayData("kernel by diff kernel comp:" + i));
			signalByKernels
					.add(net.signalKernelComponentConvolutions[i].getSignalDisplayData("signal by kernel comp:" + i));
		}
		Utilities.displaySetOfSignals(signalByKernels, "signaByKernels", "time", "value");
		Utilities.displaySetOfSignals(bBybs, "bBybs", "time", "value");
		Utilities.displaySetOfSignals(bByDbs, "bByDbs", "time", "value");
		System.out.println("seen here");
	}

	@Test
	public void improvementOnASignalAudioSignal() throws Exception {
		/**********************************/
		/******* define the mode params *****/
		/**********************************/
		ConfigurationParameters.ABSOLUTE_MACHINE_PATH = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/";
		ConfigurationParameters.SOUND_DATA_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "soundData/";
		ConfigurationParameters.numberOfKernels = 10;
		ConfigurationParameters.numberofKernelComponents = 100;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (50.0 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.initialCoefficientValue = 1;
		ConfigurationParameters.initialThresHoldValue = 0.01;
		ConfigurationParameters.INITIAL_LEARNING_RATE = 0.0001;
		ConfigurationParameters.INITIAL_THRESHOLD_GRADIENT_STEP_LENGTH = 0.00001;
		// directs if the fail-safe option need to be put in the gradient update
		ConfigurationParameters.FAIL_SAFE_ON = true;
		ConfigurationParameters.NUMBER_OF_KICKS = 100;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = false;
		ConfigurationParameters.lengthOfComponentSignals = 8820;
		ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH
				+ "stepOutPut/testOnASingleSignalSnippet/";
		// ConfigurationParameters.TIME_CONSTANT = 20;
		// ConfigurationParameters.AHP_CONSTANT = 100;
		int fileIndex = 3;
		String fileName = ConfigurationParameters.SOUND_DATA_FOLDER_PATH + "output" + fileIndex + ".txt";
		// File file = new File(fileName);
		double signaldata[][] = DataReader.readSignalPieces(fileName);
		Signal input = new Signal(signaldata[0]);
		Signal blankSignal = null;
		Network net = new Network(blankSignal);
		XYSeries inputDisplay = input.getSignalDisplayData("input signal");
		net.init(input);
		SignalUtils.storeSignal(ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "input.txt", input, false);
		List<Double> timers = new ArrayList<>();
		// Initialize a total of 6 timers
		for (int i = 0; i < 6; i++) {
			timers.add(0.0);
		}
		double currentTime = System.currentTimeMillis();
		for (int i = 0; i < 10000; i++) {
			net.init(input, false);
			timers.add(0, timers.get(0) + (System.currentTimeMillis() - currentTime));
			currentTime = System.currentTimeMillis();

			net.calculateSpikeTimesAndReconstructSignal();
			timers.add(1, timers.get(1) + (System.currentTimeMillis() - currentTime));
			currentTime = System.currentTimeMillis();

			Signal reconstructedSignal = net.getReconstructedSignal();
			List<List<Double>> allsptimes = net.kernelSpikeTimes;
			List<Double> spiktimes = net.spikeTimings;
			currentTime = System.currentTimeMillis();
			net.calculateErrorGradient();
			timers.add(2, timers.get(2) + (System.currentTimeMillis() - currentTime));
			currentTime = System.currentTimeMillis();

			net.updateKernelCoefficients();
			System.out.println("done with step#" + i);
			if (i > 100) {
				writeKernelSpikeTimes(allsptimes);
				break;
			}
			if (i % 10 == 0) {
				SignalUtils.storeSignal(
						ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "reconstructedSignal-" + i + ".txt",
						reconstructedSignal);
				FileWriter fw1 = new FileWriter(
						ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "spikeTimes-" + i + ".txt");
				BufferedWriter bw1 = new BufferedWriter(fw1);
				Utilities.writeListToBufferedWriter(bw1, spiktimes);
				bw1.close();
				/*
				 * String logKernelFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH +
				 * "initialKernel-" + i + ".txt";
				 * Utilities.writeMatrixToFile("inital kernelCoeffs", logKernelFileName,
				 * net.kernelMgr.kernelCoefficients); net.init(input, false);
				 * net.calculateSpikeTimesAndReconstructSignal(); List<XYSeries> signals = new
				 * ArrayList<>(); String reconsSignalFileName =
				 * ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "reconstructedSignal-" + i
				 * + "-spikeCount-" + net.spikeTimings.size() + ".png"; XYSeries
				 * reconsSignalData = net.getReconstructedSignal()
				 * .getSignalDisplayData("reconstructed signal after" + i + " steps");
				 * signals.add(inputDisplay); signals.add(reconsSignalData);
				 * Utilities.saveSetOfSignalsToFile(reconsSignalFileName, signals,
				 * "reconstruction after" + i + " steps", "time in ms", "signal value");
				 */

			}
		}
		System.out.println("done here");
	}

	/**
	 * This is testing sound snippets broken frequency wise
	 * 
	 * @throws Exception
	 */
	@Test
	public void improvementOnPartialAudioSignal() throws Exception {
		/**********************************/
		/******* define the mode params *****/
		/**********************************/// C:\Users\crystalonix\Downloads\compNeuroScience\researchProj\sensoryCodingWithPreciseSpikeTime\SensoryCoding
											// (2)\frequencyTestBridge
		ConfigurationParameters.ABSOLUTE_MACHINE_PATH = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/frequencyTestBridge/";
		// ConfigurationParameters.SOUND_DATA_FOLDER_PATH =
		// ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "soundData/";

		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 10;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.lengthOfBasicBSpline = 6;

		ConfigurationParameters.initialCoefficientValue = 1;
		ConfigurationParameters.initialThresHoldValue = 0.01;
		ConfigurationParameters.INITIAL_LEARNING_RATE = 0.01;
		ConfigurationParameters.INITIAL_THRESHOLD_GRADIENT_STEP_LENGTH = 0.001;
		// directs if the fail-safe option need to be put in the gradient update
		ConfigurationParameters.FAIL_SAFE_ON = true;
		ConfigurationParameters.NUMBER_OF_KICKS = 100;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = false;
		ConfigurationParameters.lengthOfComponentSignals = 8820;
		ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH
				+ "stepOutPut/testOnASingleSignalSnippet/";
		// ConfigurationParameters.TIME_CONSTANT = 20;
		// ConfigurationParameters.AHP_CONSTANT = 100;
		// int fileIndex = 3;
		String fileName = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "partialSignalSnippet" + ".txt";
		// File file = new File(fileName);
		ConfigurationParameters.numberOfSignalSegmentsInOneFile = 1;
		double signaldata[][] = DataReader.readSignalPieces(fileName);
		Signal input = new Signal(signaldata[0]);
		Signal blankSignal = null;
		Network net = new Network(blankSignal);
		// XYSeries inputDisplay = input.getSignalDisplayData("input signal");
		net.init(input);
		// SignalUtils.storeSignal(ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"input.txt",
		// input, false);
		List<Double> timers = new ArrayList<>();
		// Initialize a total of 6 timers
		for (int i = 0; i < 6; i++) {
			timers.add(0.0);
		}
		double currentTime = System.currentTimeMillis();
		for (int i = 0; i < 10000; i++) {
			net.init(input, false);
			timers.add(0, timers.get(0) + (System.currentTimeMillis() - currentTime));
			currentTime = System.currentTimeMillis();

			net.calculateSpikeTimesAndReconstructSignal();
			timers.add(1, timers.get(1) + (System.currentTimeMillis() - currentTime));
			currentTime = System.currentTimeMillis();

			Signal reconstructedSignal = net.getReconstructedSignal();
			List<List<Double>> allsptimes = net.kernelSpikeTimes;
			List<Double> spiktimes = net.spikeTimings;
			currentTime = System.currentTimeMillis();
			net.calculateErrorGradient();
			timers.add(2, timers.get(2) + (System.currentTimeMillis() - currentTime));
			currentTime = System.currentTimeMillis();

			double errorRate = net.updateKernelCoefficients();
			System.out.println("done with step#" + i);
			if (i % 10 == 0) {
				SignalUtils.storeSignal(
						"C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/frequencyTestBridge/reconsSignal-"
								+ i + ".txt",
						reconstructedSignal);
			}
			if (i > 100) {
				writeKernelSpikeTimes(allsptimes);
				break;
			}
//			if (i % 10 == 0) {
//				SignalUtils.storeSignal(ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"reconstructedSignal-"+i+".txt", reconstructedSignal);
//				FileWriter fw1 = new FileWriter(
//						ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"spikeTimes-"+i+".txt");
//				BufferedWriter bw1 = new BufferedWriter(fw1);
//				Utilities.writeListToBufferedWriter(bw1, spiktimes);
//				bw1.close();
//
//			}
		}
		System.out.println("done here");
	}

	private static void writeKernelSpikeTimes(List<List<Double>> kernelSpikeTimes) throws IOException {
		String octavePath = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/octaveWorkSpace/rasterplotFinal";
		for (int i = 0; i < kernelSpikeTimes.size(); i++) {
			List<Double> spikeTimes = kernelSpikeTimes.get(i);
			FileWriter fw1 = new FileWriter(octavePath + "spikeTimes-" + i + ".txt");
			BufferedWriter bw1 = new BufferedWriter(fw1);
			Utilities.writeListToBufferedWriter(bw1, spikeTimes);
			bw1.close();

		}
	}

	@Test
	public void testOnCorporaOfSignals() throws Exception {
		/**********************************/
		/******* define the mode params *****/
		/**********************************/
		ConfigurationParameters.numberOfKernels = 50;
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (50.0 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.initialCoefficientValue = 1;
		ConfigurationParameters.initialThresHoldValue = 0.05;
		ConfigurationParameters.INITIAL_LEARNING_RATE = 0.0001;
		ConfigurationParameters.INITIAL_THRESHOLD_GRADIENT_STEP_LENGTH = 0.0001;
		// directs if the fail-safe option need to be put in the gradient update
		ConfigurationParameters.FAIL_SAFE_ON = true;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = false;
		ConfigurationParameters.lengthOfComponentSignals = 8820;
		ConfigurationParameters.ABSOLUTE_MACHINE_PATH = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/";
		ConfigurationParameters.SOUND_DATA_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "soundData/";
		ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH
				+ "stepOutPut/testOnCorpora/";
		List<Integer> selectedFiles = new ArrayList<>();
		int stepnumber = 0;
		/*******************
		 * make a selection of files here
		 ********************/
		for (int i = 100; i >= 1; i--) {
			selectedFiles.add(i);
		}
		/*********************************************************************/
		Signal blankSignal = null;
		Network net = new Network(blankSignal);
		for (int j = 0; j < selectedFiles.size(); j++) {
			int fileIndex = selectedFiles.get(j);
			String fileName = ConfigurationParameters.SOUND_DATA_FOLDER_PATH + "output" + fileIndex + ".txt";
			// File file = new File(fileName);
			double signaldata[][] = DataReader.readSignalPieces(fileName);
			for (int k = 0; k < signaldata.length; k++) {
				Signal input = new Signal(signaldata[k]);
				XYSeries inputDisplay = input.getSignalDisplayData("input signal");
				net.init(input);
				net.calculateSpikeTimesAndReconstructSignal();
				net.calculateErrorGradient();
				double errorate = net.calculateErrorRate();// net.updateKernelCoefficients();
				System.out.println("error rate at step number:" + stepnumber + " is:" + errorate);
				if (stepnumber % 10 == 0) {
					String logKernelFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "initialKernel-"
							+ stepnumber + ".txt";
					Utilities.writeMatrixToFile("inital kernelCoeffs", logKernelFileName,
							net.kernelMgr.kernelCoefficients);
					net.init(input, false);
					net.calculateSpikeTimesAndReconstructSignal();
					List<XYSeries> signals = new ArrayList<>();
					String reconsSignalFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH
							+ "reconstructedSignal-" + stepnumber + "-spikeCount-" + net.spikeTimings.size() + ".png";
					XYSeries reconsSignalData = net.getReconstructedSignal()
							.getSignalDisplayData("reconstructed signal after" + stepnumber + " steps");
					signals.add(inputDisplay);
					signals.add(reconsSignalData);
					Utilities.saveSetOfSignalsToFile(reconsSignalFileName, signals,
							"reconstruction after" + stepnumber + " steps", "time in ms", "signal value");

				}
				stepnumber++;
			}
		}
		System.out.println("done here");
	}

	@Test
	public void checkKernelShape() throws Exception {

		int[] selectedKernels = new int[50 - 5]; // {7,16,26,36,46,49};
		for (int i = 0; i < selectedKernels.length; i++) {
			selectedKernels[i] = i + 5;
		}
		int[][] sets = new int[10][];
		sets[0] = new int[] { 0, 9, 19, 29, 39, 49 };
		sets[1] = new int[] { 1, 10, 20, 30, 40, 49 };
		sets[2] = new int[] { 2, 11, 21, 31, 41, 49 };
		sets[3] = new int[] { 3, 12, 22, 32, 42, 49 };
		sets[4] = new int[] { 4, 13, 23, 33, 43, 49 };
		sets[5] = new int[] { 5, 14, 24, 34, 44, 49 };
		sets[6] = new int[] { 6, 15, 25, 35, 45, 49 };
		sets[7] = new int[] { 7, 16, 26, 36, 46, 49 };
		sets[8] = new int[] { 8, 17, 27, 37, 47, 49 };
		sets[9] = new int[] { 9, 18, 28, 38, 48, 49 };
		Signal blankSignal = null;
		Network net = new Network(blankSignal);
		net.kernelMgr.displaySelectedKernels(selectedKernels);
		// int i = 2490;

		String folderName = "sample50kernels-11-8pm/";
		String kernelFileName = "kernelCoefficients-395000.txt";
		String kernelFile = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/stepOutPut/finalResultsFromServerExecution/"
				+ folderName + kernelFileName;// ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"testOnCorpora/"+"initialKernel-"+i+".txt";
		FileReader fr = new FileReader(kernelFile);
		BufferedReader br = new BufferedReader(fr);
		double[][] kernelCoeffs = Utilities.readMatrixToBufferedWriter(br);
		net.kernelMgr.initKernelCoefficients(kernelCoeffs);
		br.close();
		net.kernelMgr.displayKernels();
		net.kernelMgr.displaySelectedKernels(selectedKernels);
		/*
		 * for(int i=0; i<10; i++) net.kernelMgr.displaySelectedKernels(sets[i]);
		 * System.out.println("displayed the kernels");
		 */
	}

	@Test
	public void testKernelData() throws Exception {
		ConfigurationParameters.numberOfKernels = 10;
		ConfigurationParameters.numberofKernelComponents = 6;
		Signal blankSignal = null;
		Network net = new Network(blankSignal);
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			Signal tempSignal = net.kernelMgr.getKernel(i);
		}
	}

	@Test
	public void testFromServer() throws Exception {

		ConfigurationParameters.ABSOLUTE_MACHINE_PATH = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/";
		ConfigurationParameters.SOUND_DATA_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "soundData/";
		ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH
				+ "stepOutPut/testForServerRun";
		ConfigurationParameters.kernelCoeffFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH
				+ "kernelCoefficients";
		ConfigurationParameters.errorValuesFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "errorValues";
		ConfigurationParameters.spikeCountFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "spikeCounts";

		int stepNumber = 0;
		// TODO Auto-generated method stub
		// int stepNumber =0;
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
			// File file = new File(fileName);
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
				System.out.println("error rate at step number:" + stepNumber + " is:" + errorate + "first kernel coeff:"
						+ net.kernelMgr.kernelCoefficients[0][0]);
				if (stepNumber % 500 == 0) {
					Utilities.writeMatrixToFile(null,
							ConfigurationParameters.kernelCoeffFileName + "-" + stepNumber + ".txt",
							net.kernelMgr.kernelCoefficients);
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
			/*
			 * if(j%10==0){ Utilities.writeMatrixToFile(null,
			 * ConfigurationParameters.kernelCoeffFileName +".txt",
			 * net.kernelMgr.kernelCoefficients); }
			 */
		}
		System.out.println("done here");
	}

	@Test
	public void testDataPreProcessing() throws Exception {
		ConfigurationParameters.ABSOLUTE_MACHINE_PATH = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/";
		ConfigurationParameters.SOUND_DATA_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "soundData/";
		ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH
				+ "stepOutPut/testForServerRun";
		ConfigurationParameters.kernelCoeffFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH
				+ "kernelCoefficients";
		ConfigurationParameters.errorValuesFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "errorValues";
		ConfigurationParameters.spikeCountFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "spikeCounts";
		int partitionNumber = 5;
		int signalIndex = 0;
		int pieceIndex = 6;
		int kernelIndex = 40;
		String seedfileName = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "partitionSeed.txt";
		Signal parts = SignalUtils.readSignalFromFile(seedfileName);
		double[] indexes = parts.data;

		String path = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/testDataPart/dataPartition-"
				+ partitionNumber + "/";

		String processedfileName = path + "signalpiece-" + signalIndex + "-" + pieceIndex
				+ /* "-"+kernelIndex+ */".txt";

		int fileIndex = (int) indexes[(partitionNumber - 1) * ConfigurationParameters.numberOfFilesForPartition
				+ signalIndex];
		String fileName = ConfigurationParameters.SOUND_DATA_FOLDER_PATH + "output" + fileIndex + ".txt";
		// File file = new File(fileName);
		double allSignalsInThisPartition[][] = DataReader.readSignalPieces(fileName);

		Signal emptySignal = null;
		Network net = new Network(emptySignal);
		net.init(new Signal(allSignalsInThisPartition[pieceIndex]));

		Signal orgSignal = net.thisSignal;// signalKernelComponentConvolutions[kernelIndex];
		orgSignal.DrawSignal("orginial signal");

		SignalUtils.readSignalFromFile(processedfileName).DrawSignal("prepocessed Signal");
		System.out.println("I am done");
	}

	@Test
	public void writeDataPart() throws Exception {
		int partitionNumber = 5;
		ConfigurationParameters.ABSOLUTE_MACHINE_PATH = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/";
		ConfigurationParameters.SOUND_DATA_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "soundData/";
		ConfigurationParameters.PARTITION_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH
				+ "testDataPart/dataPartition-";
		ConfigurationParameters.PARTITION_SEED_FILE = ConfigurationParameters.ABSOLUTE_MACHINE_PATH
				+ "partitionSeed.txt";
		Signal emptySignal = null;
		Network net = new Network(emptySignal);
		DataReader dtr = new DataReader(net);
		dtr.processSignalDataPartition(5, 5);
	}

	@Test
	public void testWriting() throws IOException {
		List<Double> errorRates = Arrays.asList(new Double[] { 1.0, 2.0, 3.0, 4.0 });
		String fileName = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/stepOutPut/anikTest.txt";
		FileWriter fw1 = new FileWriter(fileName, true);
		// FileWriter fw2 = new FileWriter(ConfigurationParameters.spikeCountFileName +
		// "anikTest.txt", true);
		BufferedWriter bw1 = new BufferedWriter(fw1);
		Utilities.writeListToBufferedWriter(bw1, errorRates);

		// BufferedWriter bw2 = new BufferedWriter(fw2);
		// Utilities.writeIntegerListToBufferedWriter(bw2, spikeCounts);
		bw1.close();
		// bw2.close();

		errorRates = Arrays.asList(new Double[] { 10.0, 20.0, 30.0, 4.0 });
		fw1 = new FileWriter(fileName, true);
		// FileWriter fw2 = new FileWriter(ConfigurationParameters.spikeCountFileName +
		// "anikTest.txt", true);
		bw1 = new BufferedWriter(fw1);
		Utilities.writeListToBufferedWriter(bw1, errorRates);

		// BufferedWriter bw2 = new BufferedWriter(fw2);
		// Utilities.writeIntegerListToBufferedWriter(bw2, spikeCounts);
		bw1.close();
	}

	@Test
	public void testErrorAndSpikeCountAverages() throws IOException {
		String foldreName = "sample10kernels-11-11pm/";
		String filePath = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/stepOutPut/finalResultsFromServerExecution/"
				+ foldreName;
		String errorValuesFile = filePath + "errorValues.txt";
		String spikeCountsFile = filePath + "spikeCounts.txt";
		FileReader fr1 = new FileReader(errorValuesFile);
		BufferedReader br1 = new BufferedReader(fr1);
		List<Double> errorValues = Utilities.readListFromFile(br1);
		br1.close();
		List<Double> movingErrorAvgs = Utilities.convertToMovingErrorAverage(errorValues, true);

		FileReader fr2 = new FileReader(spikeCountsFile);
		BufferedReader br2 = new BufferedReader(fr2);
		List<Double> spikeCounts = Utilities.readListFromFile(br2);
		br2.close();
		List<Double> movingSpikeAvgs = Utilities.convertToMovingSpikeCountAverage(spikeCounts, true);

		double finalErrorValue = movingErrorAvgs.get(movingErrorAvgs.size() - 1);
		double finalSpikeValue = movingSpikeAvgs.get(movingSpikeAvgs.size() - 1);
		System.out.println("done");
	}

	@Test
	public void testSpikeCounts() throws IOException {
		String fileName = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/stepOutPut/finalResultsFromServerExecution/testAfterTrainOnDataPart12/"
				+ "spikeCounts.txt";
		FileReader fr = new FileReader(fileName);
		BufferedReader br = new BufferedReader(fr);
		List<Double> spikeCounts = Utilities.readListFromFile(br);
		br.close();
		List<Double> movingAvgs = Utilities.convertToMovingSpikeCountAverage(spikeCounts, true);
		double finalValue = movingAvgs.get(movingAvgs.size() - 1);
		System.out.println("done");
	}

	@Test
	public void initializeAndSaveKernelsToFile() throws Exception {
		int[] selectedKernels = { 0, 1, 2, 3 };
		String folderPath = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/octaveWorkSpace/kernelsFolderFinal2/";
		Signal nullSignal = null;
		Network net = new Network(nullSignal);

		String folderName = "sample50kernels-11-8-30pm/";
		String kernelFileName = "kernelCoefficients-402000.txt";
		String kernelFile = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/stepOutPut/finalResultsFromServerExecution/"
				+ folderName + kernelFileName;// ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"testOnCorpora/"+"initialKernel-"+i+".txt";
		FileReader fr = new FileReader(kernelFile);
		BufferedReader br = new BufferedReader(fr);
		double[][] kernelCoeffs = Utilities.readMatrixToBufferedWriter(br);
		net.kernelMgr.initKernelCoefficients(kernelCoeffs);

		net.kernelMgr.saveKernelsToFile(folderPath);
		net.kernelMgr.displaySelectedKernels(selectedKernels);
		System.out.println("done");
	}

	@Test
	public void collatePictureForReconstructionRepresentation() throws Exception {
		ConfigurationParameters.initialThresHoldValue = 0.5;
		Signal blankSignal = null;
		Network net = new Network(blankSignal);
		Signal inputSignal = SignalUtils.addTwoSignals(
				SignalUtils.addTwoSignals(net.kernelMgr.getKernel(49),
						SignalUtils.shiftSignal(net.kernelMgr.getKernel(1), 40)),
				SignalUtils.shiftSignal(net.kernelMgr.getKernel(3), 2));
		inputSignal = SignalUtils.addTwoSignals(inputSignal, SignalUtils.shiftSignal(net.kernelMgr.getKernel(20), 30));
		inputSignal = SignalUtils.addTwoSignals(inputSignal, SignalUtils.shiftSignal(net.kernelMgr.getKernel(5), 30));
		inputSignal = SignalUtils.addTwoSignals(inputSignal, SignalUtils.shiftSignal(net.kernelMgr.getKernel(15), 60));
		inputSignal.DrawSignal("input signal");
		List<Integer> selectedKernels = new ArrayList<>();
		selectedKernels.add(0);
		selectedKernels.add(2);
		selectedKernels.add(4);
		net = new Network(selectedKernels);
		net.init(inputSignal);
		net.calculateSpikeTimesAndReconstructSignal();
		Signal reconsSignal = net.getReconstructedSignal();
		reconsSignal.DrawSignal("reconstructedSignal");
		String octavePath = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/octaveWorkSpace/demoReconstruction/";
		SignalUtils.storeSignal(octavePath + "sampleInput.txt", inputSignal);
		SignalUtils.storeSignal(octavePath + "sampleReconstruction.txt", reconsSignal);
		for (int i = 0; i < net.spikeTimings.size(); i++) {
			double time = net.spikeTimings.get(i);
			int kernelIndex = net.kernelIndexesForSpikes.get(i);
			SignalUtils.storeSignal(octavePath + "invertedKernel-" + i + ".txt",
					SignalUtils.shiftSignal(net.kernelMgr.getInvertedKernel(kernelIndex), time));
		}
		FileWriter fw1 = new FileWriter(octavePath + "spikeTimes.txt");
		BufferedWriter bw1 = new BufferedWriter(fw1);
		Utilities.writeListToBufferedWriter(bw1, net.spikeTimings);
		bw1.close();
		FileWriter fw2 = new FileWriter(octavePath + "reconstructionCoeffs.txt");
		BufferedWriter bw2 = new BufferedWriter(fw2);
		Utilities.writeListToBufferedWriter(bw2, net.coefficientsOfReconstructedSignal);
		bw2.close();
		System.out.println("done");
	}

	@Test
	public void GenerateSignal() throws Exception {
		ConfigurationParameters.numberofKernelComponents = 60;
		ConfigurationParameters.numberOfKernels = 1;
		ConfigurationParameters.lengthOfBasicBSpline = 50;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = true;
		ConfigurationParameters.TIME_CONSTANT = 10;
		Signal blank = null;
		Network net = new Network(blank);
		net.kernelMgr.getKernel(0).DrawSignal("input");
		String fileName = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/octaveWorkSpace/demoReconstruction/inputsignice.txt";
		SignalUtils.storeSignal(fileName, net.kernelMgr.getKernel(0));
		System.out.println("input here");
	}

	@Test
	public void checkErrorRate() {
		int[] errorsini = { 69, 70, 72, 79, 87 };
		int[] numspikesIni = { 524, 434, 324, 216, 108 };
		int[] errorsfinal = { 32, 39, 55, 65, 74 };
		int[] numspikesfinal = { 372, 316, 198, 83, 46 };
		double[] ratesini = new double[5];
		double[] ratesfinal = new double[5];
		for (int i = 4; i >= 0; i--) {
			ratesfinal[i] = errorsfinal[i] * errorsfinal[i] * (numspikesfinal[i] / (i + 1)) / 200000.0;
			ratesini[i] = errorsini[i] * errorsini[i] * (numspikesIni[i] / (i + 1)) / 200000.0;
		}
		System.out.println("done");
	}

	@Test
	public void writeKernelForFrequencyAnalysis() throws Exception {
		/********************************/
		/** modify this part for expt. ****/
		/********************************/
		ConfigurationParameters.numberofKernelComponents = 20;
		ConfigurationParameters.numberOfKernels = 400;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.lengthOfBasicBSpline = 6;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = true;
		Signal blank = null;
		Network net = new Network(blank);
		int kernelIndex = 300;
		Signal kernelSignal = net.kernelMgr.getKernel(kernelIndex);
		kernelSignal.DrawSignal("input");
		String fileName = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/frequencyTestBridge/kernelSignal.txt";
		SignalUtils.storeSignal(fileName, kernelSignal);
		System.out.println("input here");
	}

	@Test
	public void testReconstructionWithMultipleKernelInstances() {
		/**
		 * Generate a test signal
		 */
		// Signal testSignal =
	}

	@Test
	public void testReconstructionWithNoiseOnSignalComp() throws Exception {
		/**
		 * parameter selection
		 */
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 10;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.lengthOfBasicBSpline = 6;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = true;

		/**
		 * network initialization
		 */
		Signal blank = null;
		Network net = new Network(blank);
		Signal mainSignal = null;

		int totalShift = 100;
		int numberOfSignalComps = 3;
		int kernelIndex = 0;
		List<Integer> shifts = new ArrayList<>();
		List<Integer> kernelIndexes = new ArrayList<>();
		/**
		 * iteratively add the components to the main signal
		 */
		for (int i = 0; i < numberOfSignalComps; i++) {
			shifts.add(totalShift);
			kernelIndex = new Random().nextInt(ConfigurationParameters.numberOfKernels);
			kernelIndexes.add(kernelIndex);
			Signal signalComp = net.kernelMgr.getInvertedKernel(kernelIndex);
			Signal noisyComp = SignalUtils.addNoiseToSignal(signalComp, 0.05);
			double compKernelInnerProd = SignalUtils
					.calculateSignalIntegral(SignalUtils.multiplyTwoSignals(signalComp, noisyComp));
			// calculate the scale factor now
			double deltaThres = (ConfigurationParameters.initialThresHoldValue) * 1.1;
			double scaleFactor = 1;
			if (mainSignal != null) {
				double signalKernelInnerProd = SignalUtils.calculateSignalIntegral(
						SignalUtils.multiplyTwoSignals(mainSignal, SignalUtils.shiftSignal(signalComp, totalShift)));
				scaleFactor = (deltaThres - signalKernelInnerProd) / compKernelInnerProd;
			} else {
				scaleFactor = deltaThres / compKernelInnerProd;
			}

			/**
			 * Scale the component and add it to the signal
			 */
			Signal scaledNoisyComp = SignalUtils.scalarMultiply(noisyComp, scaleFactor);
			scaledNoisyComp = SignalUtils.shiftSignal(scaledNoisyComp, totalShift);
			mainSignal = SignalUtils.addTwoSignals(scaledNoisyComp, mainSignal);
			totalShift += new Random().nextInt(100) + 50;

		}

		/**
		 * reconstruct the signal
		 */
		net.init(mainSignal);
		// mainSignal.DrawSignal("mainSignal");
		// net.kernelMgr.getInvertedKernel(kernelIndex).DrawSignal("kernel");;
		net.calculateSpikeTimesAndReconstructSignal();
		List<Double> timings = net.spikeTimings;
		List<Integer> realKernelIndexes = net.kernelIndexesForSpikes;
		System.out.println(timings.size());

		/**
		 * check for the spurious spike coefficients
		 */
		double[] coefficientsOfRecons = net.coefficientsOfReconstructedSignal;

		/**
		 * iterate through all ideal spikes and check for spurious spikes
		 */
		List<Integer> spuriousSpikes = new ArrayList<>();
		for (int i = 0, j = 0; j < coefficientsOfRecons.length; j++) {
			double idealSpikeTime = 0;
			if (i < shifts.size()) {
				idealSpikeTime = shifts.get(i);
			}
			double realSpikeTime = timings.get(j);

			if (i >= shifts.size() || Math.abs(realSpikeTime - idealSpikeTime) > 2
					|| kernelIndexes.get(i) != (int) realKernelIndexes.get(j)) {
				System.out.println("got a spurious spike here with coeff:" + coefficientsOfRecons[j]);
				// Assert.assertTrue(coefficientsOfRecons[j]<0.002);
				spuriousSpikes.add(j);
				continue;
			}
			i++;
		}
		double actualErrorRate = net.calculateErrorRate();
		System.out.println("actual error rate is:" + actualErrorRate);
		net.removeSpikes(spuriousSpikes);
		net.calculateSpikeTimesAndReconstructSignal();
		double idealErrorRate = net.calculateErrorRate();
		System.out.println("ideal error rate is:" + idealErrorRate);
		assertTrue(idealErrorRate >= actualErrorRate);
	}

	@Test
	public void normalizeKernelTest() throws Exception {
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 10;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.lengthOfBasicBSpline = 6;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = true;

		KernelManager krMgr = new KernelManager(ConfigurationParameters.numberOfKernels,
				ConfigurationParameters.numberofKernelComponents, ConfigurationParameters.lengthOfComponentSignals);

		List<Signal> signalComps = new ArrayList<>();

		double shift = 0;
		double gap = 10;
		// Construct the test signal
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			shift = shift + gap + krMgr.kernelLengths[i];
			Signal comp = SignalUtils.shiftSignal(krMgr.getInvertedKernel(i), shift);
			signalComps.add(comp);
		}

//		Signal comp = SignalUtils.shiftSignal(krMgr.getKernel(3),10);
//        signalComps.add(comp);   
//       
//        comp = SignalUtils.shiftSignal(krMgr.getKernel(9),60);
//        signalComps.add(comp);   

		Signal testSignal = SignalUtils.addNSignals(signalComps);
		testSignal.DrawSignal("Experimental Signal");
		krMgr.kernelCoefficients[0][0] = 10;
		krMgr.normalizeKernel(0, true);
		(krMgr.getKernel(0)).DrawSignal("normalized signal");
		krMgr.normalizeKernel(0, true);
		Signal normalSignal = krMgr.getKernel(0);
		(normalSignal).DrawSignal("renormalized signal");
		Signal sqSignal = SignalUtils.multiplyTwoSignals(normalSignal, normalSignal);
		double norm = SignalUtils.calculateSignalIntegral(sqSignal);
		assertTrue(norm < 1.1 && norm > 0.9);
		System.out.println();
	}

	@Test
	public void testForSpikeGen() throws Exception {
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 10;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.lengthOfBasicBSpline = 30;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = true;

		KernelManager krMgr = new KernelManager(ConfigurationParameters.numberOfKernels,
				ConfigurationParameters.numberofKernelComponents, ConfigurationParameters.lengthOfComponentSignals);

		List<Signal> signalComps = new ArrayList<>();

		double shift = 0;
		// for leaving gap between components
		double gap = 100;
		Random r = new Random();
		// Construct the test signal
		for (int i = 0; i < 3; i++) {
			int index = r.nextInt(10);
			shift = shift + gap + krMgr.kernelLengths[index];
			Signal comp = SignalUtils.shiftSignal(krMgr.getInvertedKernel(index), shift);
			signalComps.add(comp);
		}

		Signal testSignal = SignalUtils.addNSignals(signalComps);
		testSignal.DrawSignal("Experimental Signal");

		Signal blank = null;
		Network net = new Network(blank);
		net.init(testSignal);
		net.calculateSpikeTimesAndReconstructSignal();

		net.signalKernelConvolutionCache[0].DrawSignal("signal kernel conv", 0, 500);
		net.kernelMgr.getInvertedKernel(1).DrawSignal("kernel");
		// List<Double> spikeTimes = net.kernelSpikeTimes.get(1);
		net.kernelThresholds[1].DrawSignal("kernel thresholds", 0, 1000);
		System.out.println("The End");
	}

	/**
	 * Unit test for spike time differential with linear ahp
	 * 
	 * @throws Exception
	 */
	@Test
	public void testSpikeTimeDiffLinearAhp() throws Exception {
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 10;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.lengthOfBasicBSpline = 30;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = false;
		ConfigurationParameters.AHP_REFRACTORY_PERIOD = 15.0;
		ConfigurationParameters.AHP_HIGH_VALUE = 1000.0;
		ConfigurationParameters.AHP_SLOPE = ConfigurationParameters.AHP_HIGH_VALUE
				/ ConfigurationParameters.AHP_REFRACTORY_PERIOD;
		KernelManager krMgr = new KernelManager(false);

		DebugConfigurationParameters.DEBUG_MODE_ON = true;
		DebugConfigurationParameters.TEST_KERNEL_INDEX = 3;
		DebugConfigurationParameters.TEST_COMPONENT_INDEX = 0;
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
		Random r = new Random();
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
		System.out.println("Comp kernels" + selectedKers[0] + selectedKers[1] + selectedKers[2]);

		// run the network once
		Signal testSignal = SignalUtils.addNSignals(signalComps);
		testSignal = SignalUtils.scalarMultiply(testSignal, 30);
		testSignal.DrawSignal("Experimental Signal");

		Signal blank = null;
		Network net = new Network(blank);
		net.init(testSignal);
		net.calculateSpikeTimesAndReconstructSignal();

		/**
		 * Do a little bit of display here
		 */
		net.kernelMgr.getInvertedKernel(DebugConfigurationParameters.TEST_KERNEL_INDEX).DrawSignal("kernel");
		// make a slight perturbation with the kernel coeffs.

		double delB = 0.01;
		net.calculateTimeDifferentials();
		int spikeIndex = -1;
		double lastspiketime = 0;
		int toggle = 0;

		for (int i = 0; i < net.kernelIndexesForSpikes.size(); i++) {
			/*
			 * if(net.kernelIndexesForSpikes.get(i)==kernelIndex) { spikeIndex = i; break; }
			 */
			// comment: revert to this
			if (net.kernelIndexesForSpikes.get(i) == DebugConfigurationParameters.TEST_KERNEL_INDEX && toggle == 1
					&& (net.spikeTimings.get(i) - lastspiketime) < (ConfigurationParameters.AHP_REFRACTORY_PERIOD
							- 2)) {
				spikeIndex = i;
				break;
			} else if (net.kernelIndexesForSpikes.get(i) == DebugConfigurationParameters.TEST_KERNEL_INDEX) {
				toggle = 1;
				lastspiketime = net.spikeTimings.get(i);
			}
		}
		// reveal the spike index
		System.out.println("test spike index:" + spikeIndex);

		double initialTime = net.spikeTimings.get(spikeIndex);
		double timeDiff = net.timeDifferentials[spikeIndex][DebugConfigurationParameters.TEST_COMPONENT_INDEX];
		double idealTimeDiff = timeDiff * delB;

		net.kernelMgr.kernelCoefficients[DebugConfigurationParameters.TEST_KERNEL_INDEX][DebugConfigurationParameters.TEST_COMPONENT_INDEX] += delB;
		net.kernelMgr.updateCache();
		net.init(testSignal);

		// rerun the entire setup
		net.calculateSpikeTimesAndReconstructSignal();

		double finalTime = net.spikeTimings.get(spikeIndex);
		double actualTimeDiff = finalTime - initialTime;
		System.out.println("real/ideal:" + actualTimeDiff / idealTimeDiff + "actualdiff:" + actualTimeDiff);
		System.out.println("The End" + selectedKers[0] + selectedKers[1] + selectedKers[2]);
	}

	/**
	 * This is for checking if the reconstruction is good enough
	 * 
	 * @throws Exception
	 */
	@Test
	public void testReconstruction() throws Exception {
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 10;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.lengthOfBasicBSpline = 30;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = false;
		ConfigurationParameters.AHP_REFRACTORY_PERIOD = 15.0;
		ConfigurationParameters.AHP_HIGH_VALUE = 1000.0;
		ConfigurationParameters.AHP_SLOPE = ConfigurationParameters.AHP_HIGH_VALUE
				/ ConfigurationParameters.AHP_REFRACTORY_PERIOD;
		KernelManager krMgr = new KernelManager(false);

		DebugConfigurationParameters.DEBUG_MODE_ON = false;
		DebugConfigurationParameters.TEST_KERNEL_INDEX = 3;
		DebugConfigurationParameters.TEST_COMPONENT_INDEX = 0;
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
		Random r = new Random();
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
		System.out.println("Comp kernels" + selectedKers[0] + selectedKers[1] + selectedKers[2]);

		// run the network once
		Signal testSignal = SignalUtils.addNSignals(signalComps);
		testSignal = SignalUtils.scalarMultiply(testSignal, 30);
		// testSignal.DrawSignal("Experimental Signal");

		Signal blank = null;
		Network net = new Network(blank);
		net.init(testSignal);
		net.calculateSpikeTimesAndReconstructSignal();
		Signal reconstruction = net.getReconstructedSignal();
		Signal reconstructionInaccurate = net.getReconstructedSignalInaccurate();
		Signal[] signals = {  testSignal,  reconstructionInaccurate, reconstruction };
		Utilities.displaySetOfSignals(signals, "actual and reconstruction");
		System.out.println("End of reconstruction");
	}

	@Test
	public void testVectorNormalizationAndInnerProduct() {
		double[] vec = { 1, 1, 1 };
		double[] vec2 = { 0, -1, 1 };
		vec = Utilities.normalizeVector(vec);
		double val = Utilities.calculateVectorNorm(vec);
		System.out.println(vec[0]);
		Assert.assertTrue(val == 1);
		double ip = Utilities.vectorInnerProduct(vec, vec2);
		System.out.println(ip);
		assertTrue(ip == 0);
	}

	@Test
	public void testNormalizeKernel() throws Exception {
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = true;
		KernelManager krMgr = new KernelManager(true);
		int kernelIndex = 8;
		double norm = Math.sqrt(SignalUtils.calculateSquaredNorm(krMgr.getKernel(kernelIndex)));
		System.out.println(norm);
		Assert.assertTrue(norm > 0.9 && norm < 1.1);
		krMgr.normalizeKernel(kernelIndex, true);
		norm = Math.sqrt(SignalUtils.calculateSquaredNorm(krMgr.getKernel(kernelIndex)));
		System.out.println(norm);
		Assert.assertTrue(norm > 0.9 && norm < 1.1);
	}

	@Test
	public void testKernelGradG() throws Exception {
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = true;
		KernelManager krMgr = new KernelManager(true);

		int kernelIndex = 1;
		double eps = 0.1;
		double[] thisKerCoeffs = krMgr.kernelCoefficients[kernelIndex];
		double[] deltas = new double[ConfigurationParameters.numberofKernelComponents];
		for (int i = 0; i < ConfigurationParameters.numberofKernelComponents; i++) {
			deltas[i] = Math.random() * eps;
		}
		deltas = krMgr.projectAlongKernelConstraint(kernelIndex, deltas);
		thisKerCoeffs = Utilities.addTwoVectors(thisKerCoeffs, deltas);
		krMgr.kernelCoefficients[kernelIndex] = thisKerCoeffs;
		// krMgr.kernelCoefficients[kernelIndex] = deltas;
		krMgr.updateCache();
		double updatedKrNorm = SignalUtils.calculateSquaredNorm(krMgr.getKernel(kernelIndex));
		System.out.println("Updated Norm:" + updatedKrNorm);
		Assert.assertTrue(updatedKrNorm < 1.1 && updatedKrNorm > 0.9);
	}

	@Test
	public void testVectorScaleAndAdd() {
		double[] vec = { 1, 1, 1 };
		double[] vec2 = { 0, -1, 1 };
		double[] result = Utilities.addTwoVectors(vec, vec2);
		assertTrue(result[1] == 0);
		result = Utilities.scaleVector(10, vec);
		assertTrue(result[1] == 10);
		result = Utilities.scaleVector(0.1, vec2);
		assertTrue(result[1] == -0.1);
	}

	@Test
	public void testErrorGradient() throws Exception {
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 1;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.lengthOfBasicBSpline = 30;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = false;
		ConfigurationParameters.AHP_REFRACTORY_PERIOD = 150.0;
		ConfigurationParameters.AHP_HIGH_VALUE = 1000.0;
		ConfigurationParameters.AHP_SLOPE = ConfigurationParameters.AHP_HIGH_VALUE
				/ ConfigurationParameters.AHP_REFRACTORY_PERIOD;
		ConfigurationParameters.initialThresHoldValue = 120;
		KernelManager krMgr = new KernelManager(false);

		DebugConfigurationParameters.DEBUG_MODE_ON = true;
		DebugConfigurationParameters.TEST_KERNEL_INDEX = 0;
		DebugConfigurationParameters.TEST_COMPONENT_INDEX = 0;
		double delB = 0.00001;

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
		 *//*
			 * Random r = new Random(); int[] selectedKers = { 3, 8, 6 }; // Construct the
			 * test signal for (int i = 0; i < 3; i++) {
			 * 
			 * int index = r.nextInt(10); selectedKers[i] = index;
			 * 
			 * int index = selectedKers[i]; shift = shift + gap +
			 * krMgr.kernelLengths[index]; Signal comp =
			 * SignalUtils.shiftSignal(krMgr.getInvertedKernel(index), shift);
			 * signalComps.add(comp); } System.out.println("Comp kernels" + selectedKers[0]
			 * + selectedKers[1] + selectedKers[2]); Signal testSignal =
			 * SignalUtils.addNSignals(signalComps);
			 */
		Signal testSignal = SignalUtils.shiftSignal(krMgr.getInvertedKernel(0), 500);
		testSignal = SignalUtils.scalarMultiply(testSignal, 30);
		testSignal.DrawSignal("Experimental Signal");

		/**
		 * Now run the network once
		 */
		Signal blank = null;
		Network net = new Network(blank);
		net.init(testSignal);
		net.calculateSpikeTimesAndReconstructSignal();
		Signal reconstruction = net.getReconstructedSignal();
		Signal differenceSignal = SignalUtils.addTwoSignals(testSignal, SignalUtils.scalarMultiply(reconstruction, -1));
		double initialError = SignalUtils.calculateSquaredNorm(differenceSignal);

		// make a slight perturbation with the kernel coeffs.

		net.calculateErrorGradient();
		double expected = net.errorGradients[DebugConfigurationParameters.TEST_KERNEL_INDEX][DebugConfigurationParameters.TEST_COMPONENT_INDEX]
				* delB;

		// update the network state now
		net.kernelMgr.kernelCoefficients[DebugConfigurationParameters.TEST_KERNEL_INDEX][DebugConfigurationParameters.TEST_COMPONENT_INDEX] += delB;
		net.kernelMgr.updateCache();
		net.init(testSignal);

		// rerun the entire setup
		net.calculateSpikeTimesAndReconstructSignal();

		// recalculate the error
		reconstruction = net.getReconstructedSignal();
		differenceSignal = SignalUtils.addTwoSignals(testSignal, SignalUtils.scalarMultiply(reconstruction, -1));
		double finalError = SignalUtils.calculateSquaredNorm(differenceSignal);
		double real = finalError - initialError;
		double accuracy = expected / real;
		System.out.println("real:" + real + " ideal:" + expected + " accuracy:" + accuracy);
		Assert.assertTrue(accuracy > 0.9 && accuracy < 1.1);
	}

	/*	*//**
			 * This is for testing signal kernel convolution caches
			 * 
			 * @throws Exception
			 */
	/*
	 * @Test public void testSignalKernelConvCaches() throws Exception {
	 * ConfigurationParameters.numberofKernelComponents = 6;
	 * ConfigurationParameters.numberOfKernels = 10;
	 * ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 /
	 * ConfigurationParameters.numberOfKernels);
	 * ConfigurationParameters.lengthOfBasicBSpline = 30;
	 * ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = false;
	 * ConfigurationParameters.AHP_REFRACTORY_PERIOD = 15.0;
	 * ConfigurationParameters.AHP_HIGH_VALUE = 1000.0;
	 * ConfigurationParameters.AHP_SLOPE = ConfigurationParameters.AHP_HIGH_VALUE /
	 * ConfigurationParameters.AHP_REFRACTORY_PERIOD; KernelManager krMgr = new
	 * KernelManager(false);
	 * 
	 * DebugConfigurationParameters.DEBUG_MODE_ON = true;
	 * DebugConfigurationParameters.TEST_KERNEL_INDEX = 1;
	 * DebugConfigurationParameters.TEST_COMPONENT_INDEX = 1;
	 *//**
		 * construct a test signal
		 */
	/*
	 * List<Signal> signalComps = new ArrayList<>(); // initial shift of a component
	 * double shift = 0; // for leaving gap between components double gap = 100;
	 *//**
		 * Change this later to randomize
		 */

	/*
	 * Random r = new Random(); int[] selectedKers = { 3, 8, 6 }; // Construct the
	 * test signal for (int i = 0; i < 3; i++) {
	 * 
	 * int index = r.nextInt(10); selectedKers[i] = index;
	 * 
	 * int index = selectedKers[i]; shift = shift + gap +
	 * krMgr.kernelLengths[index]; Signal comp =
	 * SignalUtils.shiftSignal(krMgr.getInvertedKernel(index), shift);
	 * signalComps.add(comp); }
	 * 
	 * System.out.println("Comp kernels" + selectedKers[0] + selectedKers[1] +
	 * selectedKers[2]); Signal testSignal = SignalUtils.addNSignals(signalComps);
	 * testSignal = SignalUtils.scalarMultiply(testSignal, 30);
	 * testSignal.DrawSignal("Experimental Signal");
	 * 
	 *//**
		 * Now run the network once
		 *//*
			 * Signal blank = null; Network net = new Network(blank); net.init(testSignal);
			 * net.calculateSpikeTimesAndReconstructSignal();
			 * net.getReconstructedSignal().DrawSignal("reconstruction", 10 ,1000); Signal
			 * basicComp = net.kernelMgr.getBasisForKernelComp(DebugConfigurationParameters.
			 * TEST_KERNEL_INDEX, 0); basicComp.DrawSignal("basic comp"); Signal
			 * kernelCalcComp = new Signal(
			 * net.kernelMgr.kernelCalc.basicExpandedComponentBSplines[
			 * DebugConfigurationParameters.TEST_KERNEL_INDEX]);
			 * kernelCalcComp.DrawSignal("kernelCalcComp"); Signal conv =
			 * net.signalKernelComponentConvolutions[DebugConfigurationParameters.
			 * TEST_KERNEL_INDEX]; conv.DrawSignal("kernel Convolution", 0, 1500); double
			 * time = 150; double val = SignalUtils.calculateSignalIntegral(
			 * SignalUtils.multiplyTwoSignals(testSignal, SignalUtils.shiftSignal(basicComp,
			 * time))); System.out.println(val); }
			 */

	@Test
	public void testForSecondAndThirdTerm() throws Exception {
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 10;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.lengthOfBasicBSpline = 30;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = false;
		ConfigurationParameters.AHP_REFRACTORY_PERIOD = 200.0;
		ConfigurationParameters.AHP_HIGH_VALUE = 1000.0;
		ConfigurationParameters.AHP_SLOPE = ConfigurationParameters.AHP_HIGH_VALUE
				/ ConfigurationParameters.AHP_REFRACTORY_PERIOD;
		ConfigurationParameters.initialThresHoldValue = 10;
		KernelManager krMgr = new KernelManager(false);

		DebugConfigurationParameters.DEBUG_MODE_ON = false;
		/**
		 * construct a test signal========================================
		 */
		List<Signal> signalComps = new ArrayList<>();
		// initial shift of a component
		double shift = 0;
		// for leaving gap between components
		double gap = 20;

		Random r = new Random();
		int[] selectedKers = { 2, 1, 0, /* 2, 0 */ };
		for (int i = 0; i < selectedKers.length; i++) {
			int index = selectedKers[i];
			shift = shift + gap + krMgr.kernelLengths[index];
			Signal comp = SignalUtils.shiftSignal(krMgr.getInvertedKernel(index), shift);
			signalComps.add(comp);
		}
		Signal testSignal = SignalUtils.shiftSignal(krMgr.getKernel(0), 200);
		/* SignalUtils.addNSignals(signalComps) */;
		testSignal = SignalUtils.scalarMultiply(testSignal, 30);
		testSignal.DrawSignal("Experimental Signal");

		/**
		 * run the actual tests here========================================
		 */
		double[][] accuracies = new double[ConfigurationParameters.numberOfKernels][ConfigurationParameters.numberofKernelComponents];
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			for (int j = 0; j < ConfigurationParameters.numberofKernelComponents; j++) {
				DebugConfigurationParameters.TEST_KERNEL_INDEX = i;
				DebugConfigurationParameters.TEST_COMPONENT_INDEX = j;
				double delB = 0.00001;

				/**
				 * Now run the network once
				 */
				Signal blank = null;
				Network net = new Network(blank);
				net.init(testSignal);
				net.calculateSpikeTimesAndReconstructSignal();
				net.calculateTimeDifferentials();

				// Signal kInvertedInitial = net.kernelMgr.getInvertedKernel(0);

				Signal reconstruction = net.getReconstructedSignal();
				Signal differenceSignal = SignalUtils.addTwoSignals(testSignal,
						SignalUtils.scalarMultiply(reconstruction, -1));
				double initialError = SignalUtils.calculateSquaredNorm(differenceSignal);
				net.calculateErrorGradient();
				double expected = net.errorGradients[DebugConfigurationParameters.TEST_KERNEL_INDEX][DebugConfigurationParameters.TEST_COMPONENT_INDEX]
						* delB;
				// double dtdb = net.timeDifferentials[0][0];

				/**
				 * =========================about to rerun=================================
				 */
				// update the network state now and rerun the entire setup
				net.kernelMgr.kernelCoefficients[DebugConfigurationParameters.TEST_KERNEL_INDEX][DebugConfigurationParameters.TEST_COMPONENT_INDEX] += delB;
				net.kernelMgr.updateCache();
				net.init(testSignal);
				// Signal kInvertedWithDelB = net.kernelMgr.getInvertedKernel(0);

				net.calculateSpikeTimesAndReconstructSignal();
				// double spikeTime2 = net.spikeTimings.get(0);

				// recalculate the error
				Signal reconstruction2 = net.getReconstructedSignal();
				Signal differenceSignal2 = SignalUtils.addTwoSignals(testSignal,
						SignalUtils.scalarMultiply(reconstruction2, -1));
				// double newCoeff = net.coefficientsOfReconstructedSignal[0];
				double finalError = SignalUtils.calculateSquaredNorm(differenceSignal2);
				double real = finalError - initialError;
				double accuracy = expected / real;
				accuracies[i][j] = accuracy;
				/**
				 * ================================= show the real and ideal
				 * error==================
				 */
			}
		}
		/**
		 * ==========================================
		 */
		System.out.println("Check the accuracies");
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			System.out.println("For kernel:" + i + ":---------------->");
			for (int j = 0; j < ConfigurationParameters.numberofKernelComponents; j++) {
				System.out.println("Accuracy:" + accuracies[i][j]);
			}
		}
	}

	@Test
	public void testErrorForSingleInstance() throws Exception {
		ConfigurationParameters.numberofKernelComponents = 1;
		ConfigurationParameters.numberOfKernels = 5;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.lengthOfBasicBSpline = 30;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = false;
		ConfigurationParameters.AHP_REFRACTORY_PERIOD = 200.0;
		ConfigurationParameters.AHP_HIGH_VALUE = 1000.0;
		ConfigurationParameters.AHP_SLOPE = ConfigurationParameters.AHP_HIGH_VALUE
				/ ConfigurationParameters.AHP_REFRACTORY_PERIOD;
		ConfigurationParameters.initialThresHoldValue = 10;
		KernelManager krMgr = new KernelManager(false);

		DebugConfigurationParameters.DEBUG_MODE_ON = false;
		/**
		 * construct a test signal========================================
		 */
		List<Signal> signalComps = new ArrayList<>();
		// initial shift of a component
		double shift = 200;
		// for leaving gap between components
		double gap = 20;

		Random r = new Random();
		int[] selectedKers = { 2/* , 1, 0, 2, 0 */ };
		for (int i = 0; i < selectedKers.length; i++) {
			int index = selectedKers[i];
			shift = shift + gap + krMgr.kernelLengths[index];
			Signal comp = SignalUtils.shiftSignal(krMgr.getInvertedKernel(index), shift);
			signalComps.add(comp);
		}
		Signal testSignal = SignalUtils.addNSignals(signalComps);
		testSignal = SignalUtils.scalarMultiply(testSignal, 30);
		// testSignal.DrawSignal("Experimental Signal");

		/**
		 * run the actual tests here========================================
		 */
		DebugConfigurationParameters.TEST_KERNEL_INDEX = 3;
		DebugConfigurationParameters.TEST_COMPONENT_INDEX = 0;
		double delB = 0.00001;

		List<Integer> selectedKernels = new ArrayList<>();
		selectedKernels.add(DebugConfigurationParameters.TEST_KERNEL_INDEX);
		/**
		 * Now run the network once
		 */
		Network net = new Network(selectedKernels);
		net.init(testSignal);
		net.calculateSpikeTimesAndReconstructSignal();
		Signal kInvertedInitial = net.kernelMgr.getInvertedKernel(DebugConfigurationParameters.TEST_KERNEL_INDEX);

		Signal reconstruction = net.getReconstructedSignal();
		Signal unshiftedReconsSignal = reconstruction;
		Signal differenceSignal = SignalUtils.subtractSignal(testSignal, reconstruction);
		differenceSignal.DrawSignal("Reconstruction gap with shift=0");
		Signal unshiftedReconsGap = differenceSignal;
		double initialError = SignalUtils.calculateSquaredNorm(differenceSignal);

		/**
		 * =============store the values after first run here=======================
		 */
		double coeff = net.coefficientsOfReconstructedSignal[0];
		double spikeTime1 = net.spikeTimings.get(0);
		System.out.println("Initial error in reconstruction:" + initialError + ", spike time:" + spikeTime1
				+ ", coefficient:" + coeff);
		net.calculateErrorGradient();
		Signal ktiOld = SignalUtils.shiftSignal(
				net.kernelMgr.getInvertedKernel(DebugConfigurationParameters.TEST_KERNEL_INDEX), spikeTime1);
		double idealT2 = net.calculateSecondTermInefficiently(DebugConfigurationParameters.TEST_KERNEL_INDEX,
				DebugConfigurationParameters.TEST_COMPONENT_INDEX);
		double idealT3 = net.calculateThirdTermInefficiently(DebugConfigurationParameters.TEST_KERNEL_INDEX,
				DebugConfigurationParameters.TEST_COMPONENT_INDEX);
		System.out.println("check the ideal T2 and T3:" + idealT2 + "," + idealT3);
		double dtdb = net.timeDifferentials[0][DebugConfigurationParameters.TEST_COMPONENT_INDEX];
		double expected = net.errorGradients[DebugConfigurationParameters.TEST_KERNEL_INDEX][DebugConfigurationParameters.TEST_COMPONENT_INDEX]
				* delB;
		System.out.println("Expected error:" + expected);

		/**
		 * =========================about to rerun=================================
		 */
		net.kernelMgr.kernelCoefficients[DebugConfigurationParameters.TEST_KERNEL_INDEX][DebugConfigurationParameters.TEST_COMPONENT_INDEX] += delB;
		net.kernelMgr.updateCache();
		net.init(testSignal);
		Signal kInvertedWithDelB = net.kernelMgr.getInvertedKernel(DebugConfigurationParameters.TEST_KERNEL_INDEX);

		net.calculateSpikeTimesAndReconstructSignal();
		double spikeTime2 = net.spikeTimings.get(0);
		Signal ktiNew = SignalUtils.shiftSignal(
				net.kernelMgr.getInvertedKernel(DebugConfigurationParameters.TEST_KERNEL_INDEX), spikeTime2);
		double newCoeff = net.coefficientsOfReconstructedSignal[0];
		// recalculate the error
		Signal reconstruction2 = net.getReconstructedSignal();
		Signal differenceSignal2 = SignalUtils.subtractSignal(testSignal, reconstruction2);
		double finalError = SignalUtils.calculateSquaredNorm(differenceSignal2);
		System.out.println("final error in reconstruction:" + finalError + ", spike time:" + spikeTime2
				+ ", coefficient:" + newCoeff);
		double real = finalError - initialError;
		double accuracy = expected / real;
		System.out.println("========accuracy info here================");
		System.out.println("accuracy for kernel#" + DebugConfigurationParameters.TEST_KERNEL_INDEX + ", componentindex#"
				+ DebugConfigurationParameters.TEST_COMPONENT_INDEX + " is:" + accuracy + ", real error:" + real
				+ " ,ideal error:" + expected);
		System.out.println("========accuracy info here================");
		System.out.println("Difference in spike time:" + (spikeTime2 - spikeTime1));

		/**
		 * =========================do a little bit of postmortem here
		 * ==============================
		 */
		Signal dkForT3 = SignalUtils.addTwoSignals(SignalUtils.shiftSignal(kInvertedWithDelB, spikeTime2),
				SignalUtils.scalarMultiply(SignalUtils.shiftSignal(kInvertedWithDelB, spikeTime1), -1));
		Signal dkForT2 = SignalUtils.addTwoSignals(SignalUtils.shiftSignal(kInvertedWithDelB, spikeTime1),
				SignalUtils.scalarMultiply(SignalUtils.shiftSignal(kInvertedInitial, spikeTime1), -1));
		// dkForT3.DrawSignal("T3 diff signal", 300, 500);
		/*
		 * Signal dkCalculatedForT3 = SignalUtils.scalarMultiply(dKSignal, delB * dtdb);
		 * dkCalculatedForT3.DrawSignal("T3 diff in code", 300, 500);
		 */
		System.out.println("time diff:" + dtdb);
		double t2Numerical = -2 * coeff
				* SignalUtils.calculateSignalIntegral(SignalUtils.multiplyTwoSignals(dkForT2, differenceSignal));
		double t3Numerical = -2 * coeff
				* SignalUtils.calculateSignalIntegral(SignalUtils.multiplyTwoSignals(dkForT3, differenceSignal));
		System.out.println("Numerically computed t2, t3:" + t2Numerical + "," + t3Numerical);
		ktiNew.DrawSignal("k new shifted");
		Signal diffK = SignalUtils.subtractSignal(ktiNew, ktiOld);
		Signal[] signals = { ktiNew, ktiOld };
		// Utilities.displaySetOfSignals(signals, "shifted ks");
		// diffK.DrawSignal("diff in ks");
		Signal diffKMultByCoeffs = SignalUtils.subtractSignal(SignalUtils.scalarMultiply(ktiNew, newCoeff),
				SignalUtils.scalarMultiply(ktiOld, coeff));
		// diffKMultByCoeffs.DrawSignal("diff in coeff multiplied ks");
		SignalUtils.subtractSignal(reconstruction2, reconstruction).DrawSignal("diff in recons");
		double numericaldE = -2 * coeff
				* SignalUtils.calculateInnerProduct(differenceSignal, SignalUtils.subtractSignal(ktiNew, ktiOld));
		double actualdE = -2 * SignalUtils.calculateInnerProduct(differenceSignal,
				SignalUtils.subtractSignal(reconstruction2, reconstruction));
		double deltaSq = SignalUtils.calculateSquaredNorm(SignalUtils.subtractSignal(reconstruction2, reconstruction));
		System.out.println("Difference due to error term:" + deltaSq);
		System.out.println("Numerically computed dE:" + numericaldE + ", actual dE:" + actualdE);
		double deviation = SignalUtils
				.calculateSquaredNorm(SignalUtils.subtractSignal(reconstruction2, reconstruction));
		System.out.println("Deviation from actual reconstruction is:" + deviation);
		Signal kernelShiftedRecons = SignalUtils.scalarMultiply(SignalUtils.shiftSignal(kInvertedInitial, spikeTime1),
				coeff);
		Signal diffToNote = SignalUtils.subtractSignal(kernelShiftedRecons, reconstruction);
		double idealIp = coeff * coeff * SignalUtils.calculateSquaredNorm(kInvertedInitial);
		double reconsSqNorm = SignalUtils.calculateSquaredNorm(reconstruction);
		diffToNote.DrawSignal("Look at the difference");
		double diffInValues = SignalUtils.calculateInnerProduct(differenceSignal, diffToNote);
		System.out.println("The difference in values is:" + diffInValues + ", ideal iP:" + idealIp + ", real recons ip:"
				+ reconsSqNorm);
	}

	@Test
	public void criticalUtilsTest() throws Exception {
		KernelManager krMgr = new KernelManager(false);
		Signal sg = krMgr.getInvertedKernel(5);
		Signal test = krMgr.getKernel(3);
		System.out.println("length of inverted kernel:" + sg.getSignalSpan() + "start point:" + sg.start);
		Signal shiftedkernel1 = SignalUtils.shiftSignal(sg, 8.2);
		Signal shiftedkernel2 = SignalUtils.shiftSignal(sg, 208.2);
		Signal composedSignal = SignalUtils.addTwoSignals(shiftedkernel1, SignalUtils.shiftSignal(test, 10));
		Signal composedSignal2 = SignalUtils.addTwoSignals(shiftedkernel2, SignalUtils.shiftSignal(test, 210));
		double sgNorm1 = SignalUtils.calculateSquaredNorm(shiftedkernel1);
		double sgNorm2 = SignalUtils.calculateSquaredNorm(shiftedkernel2);
		double norm1 = SignalUtils.calculateSquaredNorm(composedSignal);
		double norm2 = SignalUtils.calculateSquaredNorm(composedSignal2);
		assertTrue(sgNorm1 == sgNorm2);
		System.out.println("first one passed");
		System.out.println((norm2 - norm1) + ", accuracy:" + (norm2 - norm1) / norm1 + " ,norm:" + norm1);
		assertTrue(norm1 == norm2);
	}

	@Test
	public void testForSpikeTimesThorough() throws Exception {
		ConfigurationParameters.numberofKernelComponents = 6;
		ConfigurationParameters.numberOfKernels = 20;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.lengthOfBasicBSpline = 30;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = false;
		ConfigurationParameters.AHP_REFRACTORY_PERIOD = 10.0;
		ConfigurationParameters.AHP_HIGH_VALUE = 1000.0;
		ConfigurationParameters.AHP_SLOPE = ConfigurationParameters.AHP_HIGH_VALUE
				/ ConfigurationParameters.AHP_REFRACTORY_PERIOD;
		ConfigurationParameters.initialThresHoldValue = 10;
		KernelManager krMgr = new KernelManager(false);

		DebugConfigurationParameters.DEBUG_MODE_ON = true;
		double[] realTimeDeltas = null;
		double[] idealTimeDeltas = null;
		Integer[] kernelIndexes = null;
		/**
		 * construct a test signal========================================
		 */
		List<Signal> signalComps = new ArrayList<>();
		// initial shift of a component
		double shift = 0;
		// for leaving gap between components
		double gap = 20;

		Random r = new Random();
		int[] selectedKers = { 2, 1, 0, 2, 0 };
		for (int i = 0; i < 3; i++) {
			int index = selectedKers[i];
			shift = shift + gap + krMgr.kernelLengths[index];
			Signal comp = SignalUtils.shiftSignal(krMgr.getInvertedKernel(index), shift);
			signalComps.add(comp);
		}
		Signal testSignal = SignalUtils.addNSignals(signalComps);
		testSignal = SignalUtils.scalarMultiply(testSignal, 30);
		// testSignal.DrawSignal("Experimental Signal");

		/**
		 * ====================== now start with the actual
		 * experiment====================
		 */
		for (int i = 0; i < ConfigurationParameters.numberOfKernels; i++) {
			for (int j = 0; j < ConfigurationParameters.numberofKernelComponents; j++) {
				DebugConfigurationParameters.TEST_KERNEL_INDEX = i;
				DebugConfigurationParameters.TEST_COMPONENT_INDEX = j;
				double delB = 0.00001;

				/**
				 * Now run the network once
				 */
				Signal blank = null;
				Network net = new Network(blank);
				net.init(testSignal);
				net.calculateSpikeTimesAndReconstructSignal();
				net.calculateTimeDifferentials();
				List<Double> initialSpikeTimes = net.spikeTimings;
				double[][] timeDifferentials = net.timeDifferentials;
				kernelIndexes = net.kernelIndexesForSpikes.toArray(new Integer[0]);
				idealTimeDeltas = new double[initialSpikeTimes.size()];
				for (int k = 0; k < initialSpikeTimes.size(); k++) {
					if (net.kernelIndexesForSpikes.get(k) == DebugConfigurationParameters.TEST_KERNEL_INDEX) {
						idealTimeDeltas[k] = timeDifferentials[k][DebugConfigurationParameters.TEST_COMPONENT_INDEX]
								* delB;
					} else {
						idealTimeDeltas[k] = 0;
					}
				}
				/**
				 * ======piece of debugging
				 */
				if (i == 0 && j == 0) {
					double ti = net.spikeTimings.get(2);
					Signal shiftedBSpline = net.kernelMgr.getBasisForKernelComp(0, 0);
					shiftedBSpline = SignalUtils.shiftSignal(shiftedBSpline,
							ti - net.kernelMgr.kernelBsplineLengths[0]);
					double num = -SignalUtils
							.calculateSignalIntegral(SignalUtils.multiplyTwoSignals(testSignal, shiftedBSpline));
					double den = SignalUtils.calculateSignalIntegral(
							SignalUtils.multiplyTwoSignals(SignalUtils.calculateSignalDifferential(testSignal),
									SignalUtils.shiftSignal(net.kernelMgr.getInvertedKernel(0), ti)));
					System.out.println("numerator and denominator:" + num + "," + den);
					System.out.println("kernel spike times:" + net.kernelSpikeTimes.get(0).get(0) + ","
							+ net.kernelSpikeTimes.get(0).get(1));
				}
				Signal initialConv = net.signalKernelConvolutionCache[0];
				// ++++++++++++++++++++++++++++++++++++++++

				/**
				 * =========================about to rerun=================================
				 */
				// update the network state now and rerun the entire setup
				net.kernelMgr.kernelCoefficients[DebugConfigurationParameters.TEST_KERNEL_INDEX][DebugConfigurationParameters.TEST_COMPONENT_INDEX] += delB;
				net.kernelMgr.updateCache();
				net.init(testSignal);
				// Signal kInvertedWithDelB = net.kernelMgr.getInvertedKernel(0);

				net.calculateSpikeTimesAndReconstructSignal();
				// double spikeTime2 = net.spikeTimings.get(0);
				Signal finalConv = net.signalKernelConvolutionCache[0];
				// Signal[] convs = { initialConv, finalConv };
				// Signal convDifference = SignalUtils.addTwoSignals(finalConv,
				// SignalUtils.scalarMultiply(initialConv, -1));
				// convDifference.DrawSignal("convolution difference", 10, 200);
				// Utilities.displaySetOfSignals(convs, "Signal Convolutions", 10, 200);

				realTimeDeltas = new double[net.spikeTimings.size()];
				for (int k = 0; k < Math.min(net.spikeTimings.size(), initialSpikeTimes.size()); k++) {
					realTimeDeltas[k] = net.spikeTimings.get(k) - initialSpikeTimes.get(k);
				}
				/**
				 * ==================Print out the accuracy values ========================
				 */
				System.out.println("Check the accuracies for kernel#" + i + ", comp index#" + j);
				for (int k = 0; k < realTimeDeltas.length; k++) {
					double real = realTimeDeltas[k];
					double ideal = idealTimeDeltas[k];
					if (real != 0 || ideal != 0)
						System.out.println("real time deltas:" + real + " ,ideal time deltas:" + ideal + ", accuracy:"
								+ ((real == 0 && ideal == 0) ? 1 : real / ideal) + ", generating kernel:"
								+ kernelIndexes[k] + ", initial time:" + initialSpikeTimes.get(k));
				}
			}
		}
	}

	/**
	 * Tests some basic signal util methods like integral, shift, scaling, norm,
	 * etc.
	 * 
	 * @throws Exception
	 */
	@Test
	public void testSignalUtilShiftAndIntegral() throws Exception {
		ConfigurationParameters.numberOfKernels = 10;
		KernelManager kr = new KernelManager(false);
		Signal testSignal1 = kr.getKernel(5);
		Signal testSignal2 = kr.getKernel(9);
		Signal shiftedSignal1 = SignalUtils.shiftSignal(testSignal1, 400.2);
		Signal shiftedSignal2 = SignalUtils.shiftSignal(testSignal2, 400.2);
		Signal scaledSignal = SignalUtils.scalarMultiply(testSignal1, 3.2);

		double norm1 = SignalUtils.calculateSquaredNorm(testSignal1);
		double normshifted = SignalUtils.calculateSquaredNorm(shiftedSignal1);
		double scaledNorm = SignalUtils.calculateSquaredNorm(scaledSignal);

		double innerProduct = SignalUtils
				.calculateSignalIntegral(SignalUtils.multiplyTwoSignals(testSignal1, testSignal2));
		double innerProductShifted = SignalUtils
				.calculateSignalIntegral(SignalUtils.multiplyTwoSignals(shiftedSignal1, shiftedSignal2));

		Assert.assertTrue(norm1 == normshifted);
		Assert.assertTrue(innerProduct == innerProductShifted);
		Assert.assertTrue(scaledNorm == norm1 * 3.2 * 3.2);// This will be approximately equal
	}

	/**
	 * Simple test for checking the error gradient
	 * 
	 * @throws Exception
	 */
	//********** This works only for 1 spike case ***********//
	@Test
	public void simplestErrorGradTest() throws Exception {
		/**
		 * ======================== Configure here ======================
		 */
		ConfigurationParameters.numberofKernelComponents = 3;
		ConfigurationParameters.numberOfKernels = 1;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = false;
		ConfigurationParameters.lengthOfBasicBSpline = 60;
		ConfigurationParameters.AHP_REFRACTORY_PERIOD = 100.0;
		ConfigurationParameters.AHP_HIGH_VALUE = 1000.0;
		ConfigurationParameters.AHP_SLOPE = ConfigurationParameters.AHP_HIGH_VALUE
				/ ConfigurationParameters.AHP_REFRACTORY_PERIOD;
		ConfigurationParameters.initialThresHoldValue = 1;

		KernelManager krMgr = new KernelManager(false);

		DebugConfigurationParameters.DEBUG_MODE_ON = true;
		DebugConfigurationParameters.TEST_KERNEL_INDEX = 0;
		DebugConfigurationParameters.TEST_COMPONENT_INDEX = 2;
		double delB = 0.00001;
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		/**
		 * ======================= Set the test signal===============================
		 */
		Signal signal = SignalUtils.shiftSignal(krMgr.getInvertedKernel(0), 500);
		Signal testSignal = SignalUtils.scalarMultiply(signal, 30);
		testSignal = SignalUtils.addTwoSignals(testSignal, SignalUtils.scalarMultiply(SignalUtils.shiftSignal(signal, 100),10));
		testSignal.DrawSignal("Signal for test");
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		/**
		 * Now run the network once
		 */
		Signal blank = null;
		Network net = new Network(blank);
		net.init(testSignal);
		net.calculateSpikeTimesAndReconstructSignal();
		net.calculateErrorGradient();
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		/**
		 * ===============Get the statistics of first run
		 * here===========================
		 */
		double spikeTimeOld = net.spikeTimings.get(0);
		List<Double> spikeTimesOld = net.spikeTimings;
		double coeffOld = net.coefficientsOfReconstructedSignal[0];
		double [] coeffsOld = net.coefficientsOfReconstructedSignal;
		Signal reconsOld = net.getReconstructedSignal();
		Signal shiftedKernelOld = SignalUtils.shiftSignal(
				net.kernelMgr.getInvertedKernel(DebugConfigurationParameters.TEST_KERNEL_INDEX), spikeTimeOld);
		Signal [] oldInvertedKernels = new Signal[spikeTimesOld.size()];
		for(int i=0;i<spikeTimesOld.size();i++) {
			oldInvertedKernels[i] = net.kernelMgr.getInvertedKernel(net.kernelIndexesForSpikes.get(i));
		}
		Signal errorSignalOld = SignalUtils.subtractSignal(testSignal, reconsOld);

		double dtdb = net.timeDifferentials[0][DebugConfigurationParameters.TEST_COMPONENT_INDEX];
		double dEdb = net.errorGradients[DebugConfigurationParameters.TEST_KERNEL_INDEX][DebugConfigurationParameters.TEST_COMPONENT_INDEX];
		double dET2 = net.calculateSecondTermInefficiently(DebugConfigurationParameters.TEST_KERNEL_INDEX,
				DebugConfigurationParameters.TEST_COMPONENT_INDEX);
		double dET2Accurate = net.calculateSecondTerm(DebugConfigurationParameters.TEST_KERNEL_INDEX, DebugConfigurationParameters.TEST_COMPONENT_INDEX);
		double dET3 = net.calculateThirdTermInefficiently(DebugConfigurationParameters.TEST_KERNEL_INDEX,
				DebugConfigurationParameters.TEST_COMPONENT_INDEX);
		double dET3Accurate = net.calculateThirdTerm(DebugConfigurationParameters.TEST_KERNEL_INDEX, DebugConfigurationParameters.TEST_COMPONENT_INDEX);

		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		/**
		 * ================= update the network state and rerun
		 * =======================================
		 */
		net.kernelMgr.kernelCoefficients[DebugConfigurationParameters.TEST_KERNEL_INDEX][DebugConfigurationParameters.TEST_COMPONENT_INDEX] += delB;
		net.kernelMgr.updateCache();
		net.init(testSignal);
		net.calculateSpikeTimesAndReconstructSignal();
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		/**
		 * ======================= Get the updated statistics
		 * here===================================
		 */
		double spikeTimeNew = net.spikeTimings.get(0);
		List<Double> spikeTimesNew = net.spikeTimings;
		double[] coeffsNew = net.coefficientsOfReconstructedSignal;
		Signal reconsNew = net.getReconstructedSignal();
		Signal shiftedKernelNew = SignalUtils.shiftSignal(
				net.kernelMgr.getInvertedKernel(DebugConfigurationParameters.TEST_KERNEL_INDEX), spikeTimeNew);
		Signal[] newInvertedKernels = new Signal[spikeTimesNew.size()];
		for (int i = 0; i < spikeTimesNew.size(); i++) {
			newInvertedKernels[i] = net.kernelMgr.getInvertedKernel(net.kernelIndexesForSpikes.get(i));
		}
		Signal errorSignalNew = SignalUtils.subtractSignal(testSignal, reconsNew);		
		

		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		/**
		 * ======================= Do the analytics here
		 * ===================================
		 */
		double reconsErrorOld = SignalUtils.calculateSquaredNorm(errorSignalOld);
		double reconsErrorNew = SignalUtils.calculateSquaredNorm(errorSignalNew);
		double realError = reconsErrorNew - reconsErrorOld;
		double idealError = dEdb * delB;
		double accuracy = idealError / realError;

		double timeDiffReal = spikeTimeNew - spikeTimeOld;
		double timeDiffIdeal = dtdb * delB;
		System.out.println("============accuracy info ================");
		System.out.println(
				"real error diff:" + realError + " ideal error diff:" + idealError + " accuracy in error:" + accuracy);
		System.out.println("real time diff:" + timeDiffReal + " ideal time diff:" + timeDiffIdeal + " accuracy in time:"
				+ timeDiffIdeal/timeDiffReal);
		System.out.println("============accuracy info ================");

		Signal diffRecons = SignalUtils.subtractSignal(reconsNew, reconsOld);
		double manualError = -2 * SignalUtils.calculateInnerProduct(diffRecons, errorSignalOld);
		System.out.println("manual error diff:" + manualError);

		System.out.println("new spike time:" + net.spikeTimings.get(0));
		double numericaldE = -2 * coeffOld * SignalUtils.calculateInnerProduct(errorSignalOld,
				SignalUtils.subtractSignal(shiftedKernelNew, shiftedKernelOld));
		System.out.println("================ sensitive================");
		System.out.println("Numerical dE:" + numericaldE);
		System.out.println("================ sensitive================");

		// compute T2 numerically
		List<Signal> shiftedKs = new ArrayList<>();
		for(int i =0; i<spikeTimesOld.size();i++) {
			Signal comp = SignalUtils.scalarMultiply(SignalUtils.subtractSignal(SignalUtils.shiftSignal(newInvertedKernels[i], spikeTimesNew.get(i)), SignalUtils.shiftSignal(oldInvertedKernels[i], spikeTimesNew.get(i))), coeffsOld[i]);
			shiftedKs.add(comp);
		}
		Signal dXb = SignalUtils.addNSignals(shiftedKs);
		double numericalT2 = -2*SignalUtils.calculateInnerProduct(errorSignalOld, dXb);
		System.out.println("======accuracy T2,t3==========");
		System.out.println("Numerical T2:"+ numericalT2+", accurate T2:"+dET2Accurate+ ", T2 used:"+ dET2);
		// compute T3 numerically
		List<Signal> shiftedKts = new ArrayList<>();
		for(int i =0; i<spikeTimesOld.size();i++) {
			Signal comp = SignalUtils.scalarMultiply(SignalUtils.subtractSignal(SignalUtils.shiftSignal(oldInvertedKernels[i], spikeTimesNew.get(i)), SignalUtils.shiftSignal(oldInvertedKernels[i], spikeTimesOld.get(i))), coeffsOld[i]);
			shiftedKts.add(comp);
		}
		Signal dXt = SignalUtils.addNSignals(shiftedKts);
		double numericalT3 = -2*SignalUtils.calculateInnerProduct(errorSignalOld, dXt);
		System.out.println("Numerical T3:"+ numericalT3+", accurate T3:"+dET3Accurate+ ", T3 used:"+ dET3);
		System.out.println("+++++++++++++++++++++++++++++++++++");
		/**
		 * ======================= put some display here
		 * ===================================
		 */
		Signal errorDiffSignal = SignalUtils.subtractSignal(reconsNew, reconsOld);
		Signal[] reconstructionDiffs = { reconsNew, reconsOld, errorDiffSignal, testSignal };
		Utilities.displaySetOfSignals(reconstructionDiffs, "check the reconstruction diffs");

		Assert.assertTrue(accuracy > 0.9 && accuracy < 1.1);
	}

	@Test
	public void simplestErrorGradTestCase() throws Exception {
		ConfigurationParameters.numberofKernelComponents = 1;
		ConfigurationParameters.numberOfKernels = 1;
		ConfigurationParameters.lengthOfBasicBSpline = 30;
		double delB = 0.1;
		int kernelIndex = 0;
		int compIndex = 0;
		/**
		 * ========================== create simplest signal
		 */
		Signal blank = null;
		Network net = new Network(blank);
		Signal testSignal = SignalUtils.shiftSignal(net.kernelMgr.getInvertedKernel(0), 300);

		/**
		 * ========================== run the network once
		 */
		net.init(testSignal);
		net.calculateSpikeTimesAndReconstructSignal();
		net.calculateErrorGradient();
		Signal recons1 = net.getReconstructedSignal();
		double spTime1 = net.spikeTimings.get(0);
		double coeff1 = net.coefficientsOfReconstructedSignal[0];
		double dtdb = net.timeDifferentials[0][compIndex];
		/**
		 * ========================== update the kernel coeffs.
		 */
		net.kernelMgr.kernelCoefficients[kernelIndex][compIndex] += delB;
		net.kernelMgr.updateCache();
		net.init(testSignal);
		net.calculateSpikeTimesAndReconstructSignal();
		/**
		 * ========================== rerun the net.
		 */
		Signal recons2 = net.getReconstructedSignal();
		double spTime2 = net.spikeTimings.get(0);
		double coeff2 = net.coefficientsOfReconstructedSignal[0];
		/**
		 * ========================== display here
		 */
		Signal[] signalsToDisplay = { testSignal, recons1, recons2 };
		Utilities.displaySetOfSignals(signalsToDisplay, "displayable signals");
		System.out.println("spike times" + spTime1 + "," + spTime2 + "real time diff:" + (spTime2 - spTime1)
				+ "ideal time diff:" + delB * dtdb + ", coefficients: " + coeff1 + "," + coeff2);
		/**
		 * ========================== compute the error terms here
		 */
		// TODO: fill this part
		double error = 0;
		System.out.println("spike times" + spTime1 + "," + spTime2 + ", coefficients: " + coeff1 + "," + coeff2);

	}	
}
