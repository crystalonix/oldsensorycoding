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
			 * int index = 0; for (int i = 0; i <
			 * nt.kernelIndexesForSpikes.size(); i++) { if
			 * (nt.kernelIndexesForSpikes.get(i) == kernelindex) { index = i;
			 * break; } }
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
			 * (nt.kernelIndexesForSpikes.get(i) == kernelindex) { index = i;
			 * break; } }
			 */

			double t0New = nt.spikeTimings.get(index);

			double error2 = nt.calculateError();
			System.out.println("real:" + (t0New - t0));
			System.out.println("ideal:" + timeDifferential * eps);
		}
		// nt.kernelMgr.dis
	}

	/*
	 * @Test public void checkCoeffDifferential() throws Exception { // TODO:
	 * next check for 16, 1 int rowIndex = 12; int colIndex = 24; int
	 * kernelindex = 6; int componentIndex = 0; double eps = 0.005;
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
	 * nt.calculateSpikeTimesAndReconstructSignal();
	 * nt.calculateErrorGradient(); //
	 * System.out.println("spike time12:"+nt.spikeTimings.get(12)+" spike //
	 * time //
	 * 28:"+nt.spikeTimings.get(28)+"timediff:"+(nt.spikeTimings.get(12)-nt.
	 * spikeTimings.get(28))); int rowSpikeIndex =
	 * nt.kernelIndexesForSpikes.get(rowIndex); double rowSpikeTime =
	 * nt.spikeTimings.get(rowIndex); int colSpikeIndex =
	 * nt.kernelIndexesForSpikes.get(colIndex); double colSpikeTime =
	 * nt.spikeTimings.get(colIndex);
	 *
	 * System.out.println("check out spike time: index, kernel index, time::" +
	 * rowIndex + "," + rowSpikeIndex + "," + rowSpikeTime + "," + colIndex +
	 * "," + colSpikeIndex + "," + colSpikeTime); // System.out.println();
	 * double[][] pdiffs = nt.pdiffs[kernelindex][componentIndex]; double[][]
	 * pmatOld = nt.pMatrix; double t0Old = nt.spikeTimings.get(rowIndex);
	 * double timeDiff = nt.timeDifferentials[rowIndex][componentIndex];
	 *
	 * System.out.println("see the p value:" + pmatOld[rowIndex][colIndex]);
	 *
	 * int alpha1 = (int) kernelMgr.getBasisForKernelComp(rowSpikeIndex,
	 * 0).getEndTime() / 3; int alpha2 = (int)
	 * kernelMgr.getBasisForKernelComp(colSpikeIndex, 0).getEndTime() / 3;
	 * double delta = nt.spikeTimings.get(colIndex) -
	 * nt.spikeTimings.get(rowIndex); System.out.println("ideal p value:" +
	 * KernelCalculator.kernelMultiplication(alpha1, alpha2, delta,
	 * ConfigurationParameters.numberofKernelComponents, null, null));
	 *//************************/
	/*
	*//** ideal p differential **/
	/*
	*//************************/
	/*
	 * double idealPdiff = KernelCalculator.kernelByBSpline(alpha2, alpha1,
	 * -delta, ConfigurationParameters.numberofKernelComponents, null) -
	 * timeDiff * KernelCalculator.kernelDifferntialMultiplication(alpha2,
	 * alpha1, -delta, ConfigurationParameters.numberofKernelComponents, null,
	 * null); System.out.println("ideal p diff*eps:" + idealPdiff * eps); //
	 * System.out.println("ideal p //
	 * diff:"+(KernelCalculator.kernelByBSpline(alpha1, alpha2, delta, //
	 * ConfigurationParameters.numberofKernelComponents, null)- //
	 * timeDiff*KernelCalculator.kernelDifferntialMultiplication(alpha1, //
	 * alpha2, delta,ConfigurationParameters.numberofKernelComponents, null, //
	 * null)));
	 *
	 * // nt.displayPByPInv();
	 *
	 * // int index=0; // for(int i=0; i<nt.kernelIndexesForSpikes.size(); i++)
	 * { // if(nt.kernelIndexesForSpikes.get(i)==kernelindex) { // rowIndex = i;
	 * // break; // } // } // colIndex = rowIndex;
	 *
	 * // Signal sg = nt.kernelMgr.getBasisForKernelComp(6, 1); //
	 * nt.kernelMgr.componentByComponentSignalConvolutionCache[6][6].DrawSignal(
	 * ""); Signal compBYDiffComp =
	 * nt.kernelMgr.componentByDiffComponentSignalConvolutionCache[2][6];
	 * compBYDiffComp.DrawSignal("2-6 comp by diff comp");
	 * nt.kernelMgr.getKernel(2).DrawSignal("kernelSignal"); double[] cOld =
	 * nt.coefficientsOfReconstructedSignal; double[] cDiffs =
	 * nt.coefficientDifferentials[kernelindex][componentIndex]; double cDiff =
	 * nt.coefficientDifferentials[kernelindex][componentIndex][rowIndex];
	 * double c0 = nt.coefficientsOfReconstructedSignal[rowIndex]; double p0 =
	 * pmatOld[rowIndex][colIndex]; double pdiff = pdiffs[rowIndex][colIndex];
	 * System.out.println("p diff real:" + pdiff); System.out.println("pvalue:"
	 * + pmatOld[rowIndex][colIndex]);
	 *
	 *//************************/
	/*
	*//******* second run ******/

	/*
	*//***********************//*
								 * nt.init(kernelSignal);
								 * nt.kernelMgr.kernelCoefficients[kernelindex][
								 * componentIndex] += eps;
								 * nt.kernelMgr.loadKernelCache();
								 * nt.calculateSpikeTimesAndReconstructSignal();
								 * System.out.println("spike time12 new:" +
								 * nt.spikeTimings.get(12) + " spike time 28:" +
								 * nt.spikeTimings.get(28) + " diff:" +
								 * (nt.spikeTimings.get(12) -
								 * nt.spikeTimings.get(28))); double t0New =
								 * nt.spikeTimings.get(rowIndex); //
								 * nt.displayPByPInv(); rowSpikeIndex =
								 * nt.kernelIndexesForSpikes.get(rowIndex);
								 * rowSpikeTime = nt.spikeTimings.get(rowIndex);
								 * colSpikeIndex =
								 * nt.kernelIndexesForSpikes.get(colIndex);
								 * colSpikeTime = nt.spikeTimings.get(colIndex);
								 * System.out.
								 * println("check out spike time: index, kernel index, time::"
								 * + rowIndex + "," + rowSpikeIndex + "," +
								 * rowSpikeTime + "," + colIndex + "," +
								 * colSpikeIndex + "," + colSpikeTime);
								 *
								 * double[][] pmatNew = nt.pMatrix; double c0New
								 * =
								 * nt.coefficientsOfReconstructedSignal[rowIndex
								 * ]; double p0New =
								 * pmatNew[rowIndex][colIndex];
								 *
								 * alpha1 = (int)
								 * kernelMgr.getBasisForKernelComp(
								 * rowSpikeIndex, 0).getEndTime() / 3; alpha2 =
								 * (int) kernelMgr.getBasisForKernelComp(
								 * colSpikeIndex, 0).getEndTime() / 3; delta =
								 * nt.spikeTimings.get(colIndex) -
								 * nt.spikeTimings.get(rowIndex);
								 * System.out.println("ideal p value:" +
								 * KernelCalculator.kernelMultiplication(alpha1,
								 * alpha2, delta, ConfigurationParameters.
								 * numberofKernelComponents, null, null));
								 *
								 * // check the pmatrix errors for (int i = 0; i
								 * < nt.spikeTimings.size(); i++) { for (int j =
								 * 0; j < nt.spikeTimings.size(); j++) { if
								 * (pdiffs[i][j] != 0) { double errorRate =
								 * (pmatNew[i][j] - pmatOld[i][j]) /
								 * ((pdiffs[i][j]) * eps); System.out.println(
								 * "real p value:" + pmatNew[i][j] + "Accuracy:"
								 * + errorRate + " for " + i + ", " + j); } } }
								 * // check coeff errors double[] cNew =
								 * nt.coefficientsOfReconstructedSignal; for
								 * (int i = 0; i < cOld.length; i++) {
								 * System.out.println("value:" + cNew[i] +
								 * "accuracy:" + (cNew[i] - cOld[i]) /
								 * (cDiffs[i] * eps)); }
								 * System.out.println("ideal spike time diff:" +
								 * timeDiff * eps);
								 * System.out.println("real spike time diff:" +
								 * (t0New - t0Old));
								 * System.out.println("see the p value now:" +
								 * pmatNew[rowIndex][colIndex]);
								 *
								 * alpha1 = (int)
								 * kernelMgr.getBasisForKernelComp(
								 * rowSpikeIndex, 0).getEndTime() / 3; alpha2 =
								 * (int) kernelMgr.getBasisForKernelComp(
								 * colSpikeIndex, 0).getEndTime() / 3; delta =
								 * nt.spikeTimings.get(colIndex) -
								 * nt.spikeTimings.get(rowIndex);
								 * System.out.println("ideal p value now:" +
								 * KernelCalculator.kernelMultiplication(alpha1,
								 * alpha2, delta, ConfigurationParameters.
								 * numberofKernelComponents, null, null));
								 *
								 * System.out.println("real:" + (c0New - c0));
								 * System.out.println("ideal:" + cDiff * eps);
								 * System.out.println("real pdiff:" + (p0New -
								 * p0)); System.out.println("ideal pdiff:" +
								 * pdiff * eps); // nt.kernelMgr.dis }
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
	 * ConfigurationParameters.numberOfKernels; kernelIndex++) { for (int
	 * compIndex = 0; compIndex <
	 * ConfigurationParameters.numberofKernelComponents; compIndex++) {
	 * System.out.println("you are looking for kernelIndex:" + kernelIndex +
	 * " component index:" + compIndex); KernelManager kernelMgr = new
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
	 * nt.calculateSpikeTimesAndReconstructSignal();
	 * nt.calculateErrorGradient(); double[][] pOld = nt.pMatrix; double
	 * pdiffs[][] = nt.pdiffs[kernelIndex][compIndex];
	 *//************************/
	/*
	*//******* second run ******/

	/*
	*//***********************//*
								 * nt.init(kernelSignal);
								 * nt.kernelMgr.kernelCoefficients[kernelIndex][
								 * compIndex] += eps;
								 * nt.kernelMgr.loadKernelCache();
								 * nt.calculateSpikeTimesAndReconstructSignal();
								 * double[][] pNew = nt.pMatrix; for (int i = 0;
								 * i < nt.spikeTimings.size(); i++) { for (int j
								 * = 0; j < nt.spikeTimings.size(); j++) { if
								 * ((pNew[i][j] == pOld[i][j]) ||
								 * Math.abs((pNew[i][j] - pOld[i][j]) /
								 * (pdiffs[i][j] * eps) - 1) < 0.1) { } //
								 * System.out.println("0"); else { double
								 * errorFrac = Math.abs((pNew[i][j] -
								 * pOld[i][j]) / (pdiffs[i][j] * eps) - 1) *
								 * 100; System.out.println("entry" + i + "," + j
								 * + ":" + "(" + pdiffs[i][j] * eps + "," +
								 * (pNew[i][j] - pOld[i][j]) + ") error rate:" +
								 * errorFrac + "%"); }
								 *
								 * System.out.println("entry number:"+i+", "+j);
								 * System.out.println("ideal pdiff:"+pdiffs[i][j
								 * ]*eps);
								 * System.out.println("real pdiff:"+(pNew[i][j]-
								 * pOld[i][j])); System.out.println(
								 * "+++++++++++++++++++++++++++++++++"); }
								 * System.out.println(); } System.out.println(
								 * "======================================"); }
								 * } }
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
	 * ConfigurationParameters.numberOfKernels; kernelIndex++) { for (int
	 * compIndex = 0; compIndex <
	 * ConfigurationParameters.numberofKernelComponents; compIndex++) {
	 * System.out.println("seeing kernelindex:" + kernelIndex + " comp index:" +
	 * compIndex); KernelManager kernelMgr = new
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
	 * nt.calculateSpikeTimesAndReconstructSignal();
	 * nt.calculateErrorGradient(); double[] oldCoeffs =
	 * nt.coefficientsOfReconstructedSignal; double pdiffs[][] =
	 * nt.pdiffs[kernelIndex][compIndex]; double[][] alphaArray = new
	 * double[1][nt.coefficientsOfReconstructedSignal.length]; alphaArray[0] =
	 * nt.coefficientsOfReconstructedSignal; SimpleMatrix alphaColumn = (new
	 * SimpleMatrix(alphaArray)).transpose(); SimpleMatrix pdiffByAlpha = new
	 * SimpleMatrix(pdiffs).mult(alphaColumn); double coeffDiffMat[][] =
	 * Utilities.getArrayFromMatrix(nt.pInvMatrix.mult(pdiffByAlpha).transpose()
	 * ); double idealCoeffDiffs[] = (Utilities.scale(coeffDiffMat, -1))[0];
	 *
	 *//************************/
	/*
	*//******* second run ******/

	/*
	*//***********************//*
								 * nt.init(kernelSignal);
								 * nt.kernelMgr.kernelCoefficients[kernelIndex][
								 * compIndex] += eps;
								 * nt.kernelMgr.loadKernelCache();
								 * nt.calculateSpikeTimesAndReconstructSignal();
								 * double newCoeffs[] =
								 * nt.coefficientsOfReconstructedSignal; for
								 * (int i = 0; i < newCoeffs.length; i++) { if
								 * ((newCoeffs[i] == oldCoeffs[i]) ||
								 * ((((newCoeffs[i] - oldCoeffs[i]) /
								 * (idealCoeffDiffs[i] * eps) - 1) < 0.1) &&
								 * (((newCoeffs[i] - oldCoeffs[i]) /
								 * (idealCoeffDiffs[i] * eps) - 1) > 0))) { } //
								 * System.out.println("0"); else { double
								 * errorFrac = ((newCoeffs[i] - oldCoeffs[i]) /
								 * (idealCoeffDiffs[i] * eps) - 1) * 100;
								 * System.out.println("entry" + i + ":" +
								 * oldCoeffs[i] + " (" + idealCoeffDiffs[i] *
								 * eps + "," + (newCoeffs[i] - oldCoeffs[i]) +
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
	 * @Test public void checkAlphas() throws Exception { // TODO: next check
	 * for 16, 1
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
								 * nt.calculateSpikeTimesAndReconstructSignal();
								 * nt.calculateErrorGradient(); for (int i = 0;
								 * i < nt.spikeTimings.size(); i++) { int ji =
								 * nt.kernelIndexesForSpikes.get(i); double ti =
								 * nt.spikeTimings.get(i); double value = 0; for
								 * (int k = 0; k < nt.spikeTimings.size(); k++)
								 * { int jk = nt.kernelIndexesForSpikes.get(k);
								 * int alpha1 = (int)
								 * kernelMgr.getBasisForKernelComp(jk,
								 * 0).getEndTime() / 3; int alpha2 = (int)
								 * kernelMgr.getBasisForKernelComp(ji,
								 * 0).getEndTime() / 3; double tk =
								 * nt.spikeTimings.get(k); value +=
								 * nt.coefficientsOfReconstructedSignal[k]
								 * KernelCalculator.kernelMultiplication(alpha1,
								 * alpha2, ti - tk, 4, null, null); }
								 * System.out.println("multiplied value:" +
								 * value); }
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
		ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "stepOutPut/testOnASingleSignalSnippet/";
		//ConfigurationParameters.TIME_CONSTANT = 20;
		//ConfigurationParameters.AHP_CONSTANT = 100;
		int fileIndex = 3;
		String fileName = ConfigurationParameters.SOUND_DATA_FOLDER_PATH + "output" + fileIndex + ".txt";
		//File file = new File(fileName);
		double signaldata[][] = DataReader.readSignalPieces(fileName);
		Signal input = new Signal(signaldata[0]);
		Signal blankSignal = null;
		Network net = new Network(blankSignal);
		XYSeries inputDisplay = input.getSignalDisplayData("input signal");
		net.init(input);
		SignalUtils.storeSignal(ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"input.txt", input, false);
		List<Double> timers = new ArrayList<>();
		// Initialize a total of 6 timers
		for(int i=0; i<6; i++){
			timers.add(0.0);
		}
		double currentTime = System.currentTimeMillis();
		for (int i = 0; i < 10000; i++) {
			net.init(input, false);
			timers.add(0,timers.get(0)+(System.currentTimeMillis()-currentTime));
			currentTime = System.currentTimeMillis();

			net.calculateSpikeTimesAndReconstructSignal();
			timers.add(1,timers.get(1)+(System.currentTimeMillis()-currentTime));
			currentTime = System.currentTimeMillis();

			Signal reconstructedSignal = net.getReconstructedSignal();
			List<List<Double>> allsptimes = net.kernelSpikeTimes;
			List<Double> spiktimes = net.spikeTimings;
			currentTime = System.currentTimeMillis();
			net.calculateErrorGradient();
			timers.add(2,timers.get(2)+(System.currentTimeMillis()-currentTime));
			currentTime = System.currentTimeMillis();

			net.updateKernelCoefficients();
			System.out.println("done with step#" + i);
			if(i>100){
				writeKernelSpikeTimes(allsptimes);
				break;
			}
			if (i % 10 == 0) {
				SignalUtils.storeSignal(ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"reconstructedSignal-"+i+".txt", reconstructedSignal);
				FileWriter fw1 = new FileWriter(
						ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"spikeTimes-"+i+".txt");
				BufferedWriter bw1 = new BufferedWriter(fw1);
				Utilities.writeListToBufferedWriter(bw1, spiktimes);
				bw1.close();
				/*String logKernelFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "initialKernel-" + i
						+ ".txt";
				Utilities.writeMatrixToFile("inital kernelCoeffs", logKernelFileName, net.kernelMgr.kernelCoefficients);
				net.init(input, false);
				net.calculateSpikeTimesAndReconstructSignal();
				List<XYSeries> signals = new ArrayList<>();
				String reconsSignalFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "reconstructedSignal-"
						+ i + "-spikeCount-" + net.spikeTimings.size() + ".png";
				XYSeries reconsSignalData = net.getReconstructedSignal()
						.getSignalDisplayData("reconstructed signal after" + i + " steps");
				signals.add(inputDisplay);
				signals.add(reconsSignalData);
				Utilities.saveSetOfSignalsToFile(reconsSignalFileName, signals, "reconstruction after" + i + " steps",
						"time in ms", "signal value");*/

			}
		}
		System.out.println("done here");
	}

	/**
	 * This is testing sound snippets broken frequency wise
	 * @throws Exception
	 */
	@Test
	public void improvementOnPartialAudioSignal() throws Exception {
		/**********************************/
		/******* define the mode params *****/
		/**********************************///C:\Users\crystalonix\Downloads\compNeuroScience\researchProj\sensoryCodingWithPreciseSpikeTime\SensoryCoding (2)\frequencyTestBridge
		ConfigurationParameters.ABSOLUTE_MACHINE_PATH = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/frequencyTestBridge/";
		//ConfigurationParameters.SOUND_DATA_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "soundData/";

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
		ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "stepOutPut/testOnASingleSignalSnippet/";
		//ConfigurationParameters.TIME_CONSTANT = 20;
		//ConfigurationParameters.AHP_CONSTANT = 100;
		//int fileIndex = 3;
		String fileName = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "partialSignalSnippet" + ".txt";
		//File file = new File(fileName);
		ConfigurationParameters.numberOfSignalSegmentsInOneFile = 1;
		double signaldata[][] = DataReader.readSignalPieces(fileName);
		Signal input = new Signal(signaldata[0]);
		Signal blankSignal = null;
		Network net = new Network(blankSignal);
		//XYSeries inputDisplay = input.getSignalDisplayData("input signal");
		net.init(input);
		//SignalUtils.storeSignal(ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"input.txt", input, false);
		List<Double> timers = new ArrayList<>();
		// Initialize a total of 6 timers
		for(int i=0; i<6; i++){
			timers.add(0.0);
		}
		double currentTime = System.currentTimeMillis();
		for (int i = 0; i < 10000; i++) {
			net.init(input, false);
			timers.add(0,timers.get(0)+(System.currentTimeMillis()-currentTime));
			currentTime = System.currentTimeMillis();

			net.calculateSpikeTimesAndReconstructSignal();
			timers.add(1,timers.get(1)+(System.currentTimeMillis()-currentTime));
			currentTime = System.currentTimeMillis();

			Signal reconstructedSignal = net.getReconstructedSignal();
			List<List<Double>> allsptimes = net.kernelSpikeTimes;
			List<Double> spiktimes = net.spikeTimings;
			currentTime = System.currentTimeMillis();
			net.calculateErrorGradient();
			timers.add(2,timers.get(2)+(System.currentTimeMillis()-currentTime));
			currentTime = System.currentTimeMillis();

			double errorRate = net.updateKernelCoefficients();
			System.out.println("done with step#" + i);
			if(i%10==0){
				SignalUtils.storeSignal("C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/frequencyTestBridge/reconsSignal-"+i+".txt", reconstructedSignal);
			}
			if(i>100){
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
	private static void writeKernelSpikeTimes(List<List<Double>> kernelSpikeTimes) throws IOException{
		String octavePath = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/octaveWorkSpace/rasterplotFinal";
		for(int i=0; i<kernelSpikeTimes.size(); i++){
			List<Double> spikeTimes = kernelSpikeTimes.get(i);
			FileWriter fw1 = new FileWriter(octavePath+"spikeTimes-"+i+".txt");
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
		ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH +"stepOutPut/testOnCorpora/";
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

		int []selectedKernels = new int[50-5]; //{7,16,26,36,46,49};
		for(int i=0; i<selectedKernels.length; i++){
			selectedKernels[i] = i+5;
		}
		int [][]sets = new int[10][];
		sets[0]= new int[]{0,9,19,29,39,49};
		sets[1]= new int[]{1,10,20,30,40,49};
		sets[2]= new int[]{2,11,21,31,41,49};
		sets[3]= new int[]{3,12,22,32,42,49};
		sets[4]= new int[]{4,13,23,33,43,49};
		sets[5]= new int[]{5,14,24,34,44,49};
		sets[6]= new int[]{6,15,25,35,45,49};
		sets[7]= new int[]{7,16,26,36,46,49};
		sets[8]= new int[]{8,17,27,37,47,49};
		sets[9]= new int[]{9,18,28,38,48,49};
		Signal blankSignal = null;
		Network net = new Network(blankSignal);
		net.kernelMgr.displaySelectedKernels(selectedKernels);
		// int i = 2490;

		String folderName = "sample50kernels-11-8pm/";
		String kernelFileName = "kernelCoefficients-395000.txt";
		String kernelFile = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/stepOutPut/finalResultsFromServerExecution/"+folderName+ kernelFileName;// ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"testOnCorpora/"+"initialKernel-"+i+".txt";
		FileReader fr = new FileReader(kernelFile);
		BufferedReader br = new BufferedReader(fr);
		double[][] kernelCoeffs = Utilities.readMatrixToBufferedWriter(br);
		net.kernelMgr.initKernelCoefficients(kernelCoeffs);
		br.close();
		net.kernelMgr.displayKernels();
		net.kernelMgr.displaySelectedKernels(selectedKernels);
		/*for(int i=0; i<10; i++)
		net.kernelMgr.displaySelectedKernels(sets[i]);
		System.out.println("displayed the kernels");*/
	}

	@Test
    public void testKernelData() throws Exception{
    	ConfigurationParameters.numberOfKernels = 10;
    	ConfigurationParameters.numberofKernelComponents = 6;
    	Signal blankSignal = null;
    	Network net = new Network(blankSignal);
    	for(int i=0 ; i<ConfigurationParameters.numberOfKernels; i++)
    		{Signal tempSignal = net.kernelMgr.getKernel(i);}
    }

	@Test
	public void testFromServer() throws Exception {

		ConfigurationParameters.ABSOLUTE_MACHINE_PATH = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/";
		ConfigurationParameters.SOUND_DATA_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "soundData/";
		ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH+ "stepOutPut/testForServerRun";
		ConfigurationParameters.kernelCoeffFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "kernelCoefficients";
		ConfigurationParameters.errorValuesFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"errorValues";
		ConfigurationParameters.spikeCountFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"spikeCounts";


		int stepNumber =0;
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
						System.out.println("error rate at step number:"+stepNumber+ " is:" + errorate+"first kernel coeff:"+net.kernelMgr.kernelCoefficients[0][0]);
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
					/*if(j%10==0){
						Utilities.writeMatrixToFile(null, ConfigurationParameters.kernelCoeffFileName +".txt", net.kernelMgr.kernelCoefficients);
					}*/
				}
				System.out.println("done here");
	}

	@Test
	public void testDataPreProcessing() throws Exception{
		ConfigurationParameters.ABSOLUTE_MACHINE_PATH = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/";
		ConfigurationParameters.SOUND_DATA_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "soundData/";
		ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH+ "stepOutPut/testForServerRun";
		ConfigurationParameters.kernelCoeffFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH + "kernelCoefficients";
		ConfigurationParameters.errorValuesFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"errorValues";
		ConfigurationParameters.spikeCountFileName = ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"spikeCounts";
		int partitionNumber = 5;
		int signalIndex = 0;
		int pieceIndex = 6;
		int kernelIndex = 40;
		String seedfileName = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "partitionSeed.txt";
		Signal parts = SignalUtils.readSignalFromFile(seedfileName);
		double[] indexes = parts.data;

		String path = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/testDataPart/dataPartition-"+partitionNumber+"/";

		String processedfileName = path + "signalpiece-"+signalIndex+"-"+pieceIndex+/*"-"+kernelIndex+*/".txt";

		int fileIndex = (int) indexes[(partitionNumber-1)*ConfigurationParameters.numberOfFilesForPartition+signalIndex];
		String fileName = ConfigurationParameters.SOUND_DATA_FOLDER_PATH + "output" + fileIndex + ".txt";
		//File file = new File(fileName);
		double allSignalsInThisPartition[][] = DataReader.readSignalPieces(fileName);

		Signal emptySignal = null;
		Network net = new Network(emptySignal);
		net.init(new Signal(allSignalsInThisPartition[pieceIndex]));

		Signal orgSignal = net.thisSignal;//signalKernelComponentConvolutions[kernelIndex];
		orgSignal.DrawSignal("orginial signal");

		SignalUtils.readSignalFromFile(processedfileName).DrawSignal("prepocessed Signal");
		System.out.println("I am done");
	}

	@Test
	public void writeDataPart() throws Exception {
		int partitionNumber = 5;
		ConfigurationParameters.ABSOLUTE_MACHINE_PATH = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/";
		ConfigurationParameters.SOUND_DATA_FOLDER_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "soundData/";
		ConfigurationParameters.PARTITION_PATH = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "testDataPart/dataPartition-";
		ConfigurationParameters.PARTITION_SEED_FILE = ConfigurationParameters.ABSOLUTE_MACHINE_PATH + "partitionSeed.txt";
		Signal emptySignal = null;
		Network net = new Network(emptySignal);
		DataReader dtr = new DataReader(net);
		dtr.processSignalDataPartition(5, 5);
	}
	@Test
	public void testWriting() throws IOException{
		List<Double> errorRates = Arrays.asList(new Double[]{1.0,2.0,3.0,4.0});
		String fileName = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/stepOutPut/anikTest.txt";
		FileWriter fw1 = new FileWriter(fileName, true);
		//FileWriter fw2 = new FileWriter(ConfigurationParameters.spikeCountFileName + "anikTest.txt", true);
		BufferedWriter bw1 = new BufferedWriter(fw1);
		Utilities.writeListToBufferedWriter(bw1, errorRates);

		//BufferedWriter bw2 = new BufferedWriter(fw2);
		//Utilities.writeIntegerListToBufferedWriter(bw2, spikeCounts);
		bw1.close();
		//bw2.close();

		 errorRates = Arrays.asList(new Double[]{10.0,20.0,30.0,4.0});
		fw1 = new FileWriter(fileName, true);
		//FileWriter fw2 = new FileWriter(ConfigurationParameters.spikeCountFileName + "anikTest.txt", true);
		bw1 = new BufferedWriter(fw1);
		Utilities.writeListToBufferedWriter(bw1, errorRates);

		//BufferedWriter bw2 = new BufferedWriter(fw2);
		//Utilities.writeIntegerListToBufferedWriter(bw2, spikeCounts);
		bw1.close();
	}

	@Test
	public void testErrorAndSpikeCountAverages() throws IOException{
		String foldreName = "sample10kernels-11-11pm/";
		String filePath = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/stepOutPut/finalResultsFromServerExecution/"+foldreName;
		String errorValuesFile = filePath+"errorValues.txt";
		String spikeCountsFile = filePath+"spikeCounts.txt";
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

		double finalErrorValue = movingErrorAvgs.get(movingErrorAvgs.size()-1);
		double finalSpikeValue = movingSpikeAvgs.get(movingSpikeAvgs.size()-1);
		System.out.println("done");
	}

	@Test
	public void testSpikeCounts() throws IOException{
		String fileName = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/stepOutPut/finalResultsFromServerExecution/testAfterTrainOnDataPart12/"+
				"spikeCounts.txt";
		FileReader fr = new FileReader(fileName);
		BufferedReader br = new BufferedReader(fr);
		List<Double> spikeCounts = Utilities.readListFromFile(br);
		br.close();
		List<Double> movingAvgs = Utilities.convertToMovingSpikeCountAverage(spikeCounts, true);
		double finalValue = movingAvgs.get(movingAvgs.size()-1);
		System.out.println("done");
	}

	@Test
	public void initializeAndSaveKernelsToFile() throws Exception{
		int [] selectedKernels = {0,1,2,3};
		String folderPath = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/octaveWorkSpace/kernelsFolderFinal2/";
		Signal nullSignal = null;
		Network net = new Network(nullSignal);


		String folderName = "sample50kernels-11-8-30pm/";
		String kernelFileName = "kernelCoefficients-402000.txt";
		String kernelFile = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/SensoryCoding/src/sensoryCoding/stepOutPut/finalResultsFromServerExecution/"+folderName+ kernelFileName;// ConfigurationParameters.STEP_OUTPUT_FOLDER_PATH+"testOnCorpora/"+"initialKernel-"+i+".txt";
		FileReader fr = new FileReader(kernelFile);
		BufferedReader br = new BufferedReader(fr);
		double[][] kernelCoeffs = Utilities.readMatrixToBufferedWriter(br);
		net.kernelMgr.initKernelCoefficients(kernelCoeffs);




		net.kernelMgr.saveKernelsToFile(folderPath);
		net.kernelMgr.displaySelectedKernels(selectedKernels);
		System.out.println("done");
	}


	@Test
	public void collatePictureForReconstructionRepresentation() throws Exception{
		ConfigurationParameters.initialThresHoldValue= 0.5;
		Signal blankSignal = null;
		Network net = new Network(blankSignal);
		Signal inputSignal = SignalUtils.addTwoSignals(SignalUtils.addTwoSignals(net.kernelMgr.getKernel(49), SignalUtils.shiftSignal(net.kernelMgr.getKernel(1), 40)),  SignalUtils.shiftSignal(net.kernelMgr.getKernel(3), 2));
		inputSignal = SignalUtils.addTwoSignals(inputSignal, SignalUtils.shiftSignal(net.kernelMgr.getKernel(20),30));
		inputSignal = SignalUtils.addTwoSignals(inputSignal, SignalUtils.shiftSignal(net.kernelMgr.getKernel(5),30));
		inputSignal = SignalUtils.addTwoSignals(inputSignal, SignalUtils.shiftSignal(net.kernelMgr.getKernel(15), 60));
		inputSignal.DrawSignal("input signal");
		List<Integer> selectedKernels = new ArrayList<>();
		selectedKernels.add(0);
		selectedKernels.add(2);
		selectedKernels.add(4);
		net  = new Network(selectedKernels);
		net.init(inputSignal);
		net.calculateSpikeTimesAndReconstructSignal();
		Signal reconsSignal = net.getReconstructedSignal();
		reconsSignal.DrawSignal("reconstructedSignal");
		String octavePath = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/octaveWorkSpace/demoReconstruction/";
		SignalUtils.storeSignal(octavePath+ "sampleInput.txt",inputSignal);
		SignalUtils.storeSignal(octavePath+ "sampleReconstruction.txt",reconsSignal);
		for(int i=0; i<net.spikeTimings.size(); i++){
			double time = net.spikeTimings.get(i);
			int kernelIndex = net.kernelIndexesForSpikes.get(i);
			SignalUtils.storeSignal(octavePath+"invertedKernel-"+i+".txt", SignalUtils.shiftSignal(net.kernelMgr.getInvertedKernel(kernelIndex),time));
		}
		FileWriter fw1 = new FileWriter(octavePath+"spikeTimes.txt");
		BufferedWriter bw1 = new BufferedWriter(fw1);
		Utilities.writeListToBufferedWriter(bw1, net.spikeTimings);
		bw1.close();
		FileWriter fw2 = new FileWriter(octavePath+"reconstructionCoeffs.txt");
		BufferedWriter bw2 = new BufferedWriter(fw2);
		Utilities.writeListToBufferedWriter(bw2, net.coefficientsOfReconstructedSignal);
		bw2.close();
		System.out.println("done");
	}

	@Test
	public void GenerateSignal() throws Exception{
		ConfigurationParameters.numberofKernelComponents = 60;
		ConfigurationParameters.numberOfKernels =1 ;
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
	public void checkErrorRate(){
		int [] errorsini = {69, 70, 72, 79, 87};
		int [] numspikesIni = {524, 434, 324, 216, 108};
		int [] errorsfinal = {32, 39, 55, 65, 74};
		int [] numspikesfinal = {372, 316, 198, 83, 46};
		double [] ratesini = new double[5];
		double [] ratesfinal = new double[5];
		for(int i=4; i>=0; i--){
			ratesfinal[i] = errorsfinal[i]*errorsfinal[i]*(numspikesfinal[i]/(i+1))/200000.0;
			ratesini[i]	= errorsini[i]*errorsini[i]*(numspikesIni[i]/(i+1))/200000.0;
		}
		System.out.println("done");
	}

	@Test
	public void writeKernelForFrequencyAnalysis() throws Exception{
		/********************************/
		/**modify this part for expt.****/
		/********************************/
		ConfigurationParameters.numberofKernelComponents = 20;
		ConfigurationParameters.numberOfKernels = 400;
		ConfigurationParameters.FREQUENCY_SCALING_FACTOR = 1.0 * (1.5 / ConfigurationParameters.numberOfKernels);
		ConfigurationParameters.lengthOfBasicBSpline = 6;
		ConfigurationParameters.SHOULD_RANDOMIZE_KERNEL_COEFFICIENTS = true;
		Signal blank = null;
		Network net = new Network(blank);
		int kernelIndex =300;
		Signal kernelSignal = net.kernelMgr.getKernel(kernelIndex);
		kernelSignal.DrawSignal("input");
		String fileName = "C:/Users/crystalonix/Downloads/compNeuroScience/researchProj/sensoryCodingWithPreciseSpikeTime/SensoryCoding (2)/frequencyTestBridge/kernelSignal.txt";
		SignalUtils.storeSignal(fileName, kernelSignal);
		System.out.println("input here");
	}

	@Test
	public void testReconstructionWithMultipleKernelInstances(){
		/**
		 * Generate a test signal
		 */
		//Signal testSignal =
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
		for(int i=0; i<numberOfSignalComps; i++) {
			shifts.add(totalShift);
			kernelIndex = new Random().nextInt(ConfigurationParameters.numberOfKernels);
			kernelIndexes.add(kernelIndex);
			Signal signalComp = net.kernelMgr.getInvertedKernel(kernelIndex);
			Signal noisyComp = SignalUtils.addNoiseToSignal(signalComp, 0.05);
			double compKernelInnerProd = SignalUtils.calculateSignalIntegral(SignalUtils.multiplyTwoSignals(signalComp, noisyComp));
			// calculate the scale factor now
			double deltaThres = (ConfigurationParameters.initialThresHoldValue)*1.1;
			double scaleFactor = 1;
			if(mainSignal!= null) {
				double signalKernelInnerProd = SignalUtils.calculateSignalIntegral(SignalUtils.multiplyTwoSignals(mainSignal, SignalUtils.shiftSignal(signalComp, totalShift)));
				scaleFactor = (deltaThres- signalKernelInnerProd)/compKernelInnerProd;
			}
			else {
				scaleFactor = deltaThres/compKernelInnerProd;
			}
			
			/**
			 * Scale the component and add it to the signal
			 */
			Signal scaledNoisyComp = SignalUtils.scalarMultiply(noisyComp, scaleFactor);		
			scaledNoisyComp = SignalUtils.shiftSignal(scaledNoisyComp, totalShift);			
			mainSignal = SignalUtils.addTwoSignals(scaledNoisyComp, mainSignal);			
			totalShift+= new Random().nextInt(100) +50;
			
		}
		
		/**
		 * reconstruct the signal
		 */
		net.init(mainSignal);
		// mainSignal.DrawSignal("mainSignal");
		// net.kernelMgr.getInvertedKernel(kernelIndex).DrawSignal("kernel");;
		net.calculateSpikeTimesAndReconstructSignal();
		List<Double> timings= net.spikeTimings;
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
		for(int i=0,j=0; j<coefficientsOfRecons.length; j++) {
			double idealSpikeTime=0;
			if (i < shifts.size()) {
				idealSpikeTime = shifts.get(i);
			}
			double realSpikeTime = timings.get(j);
			
			if(i>=shifts.size()||Math.abs(realSpikeTime-idealSpikeTime)>2||kernelIndexes.get(i)!=(int)realKernelIndexes.get(j)) {
				System.out.println("got a spurious spike here with coeff:"+coefficientsOfRecons[j]);
				// Assert.assertTrue(coefficientsOfRecons[j]<0.002);
				spuriousSpikes.add(j);
				continue;
			}
			i++;
		}
		double actualErrorRate = net.calculateErrorRate();
		System.out.println("actual error rate is:"+ actualErrorRate);
		net.removeSpikes(spuriousSpikes);
		net.reconstructSignal();
		double idealErrorRate = net.calculateErrorRate();		
		System.out.println("ideal error rate is:"+ idealErrorRate);
		assertTrue(idealErrorRate>=actualErrorRate);
	}
}
