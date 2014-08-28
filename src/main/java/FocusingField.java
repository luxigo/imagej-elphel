/**
**
** FocusingField.jave - save/restore/process sagittal/tangential PSF width
** over FOV, together with related data
**
** Copyright (C) 2014 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**
** CalibrationHardwareInterface.java is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program. If not, see <http://www.gnu.org/licenses/>.
** -----------------------------------------------------------------------------**
**
*/

import ij.IJ;
import ij.gui.GenericDialog;
import ij.text.TextWindow;

import java.awt.Point;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Properties;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;








//import Distortions.LMAArrays; // may still reuse?
import Jama.LUDecomposition;
import Jama.Matrix;


public class FocusingField {
//    public String path;
// restored from properties	
	FieldFitting fieldFitting=null;
	public double pX0_distortions;
	public double pY0_distortions;
	public double currentPX0;
	public double currentPY0;
	boolean sagittalMaster; // center data is the same, when true sagittal fitting only may change r=0 coefficients,
	boolean parallelOnly; // only process measurements for parallel moves
	boolean filterInput;
	double  filterInputMotorDiff;
	double  filterInputDiff; // um
	boolean filterInputFirstLast; 
	boolean filterInputTooFar; // filter samples that are too far from the "center of mass" of other samples 
	double  filterInputFarRatio; // remove samples that are farther than this ration of average distance
	boolean filterInputConcave;
	double  filterInputConcaveSigma;
	boolean filterInputConcaveRemoveFew;
	int filterInputConcaveMinSeries;
	double  filterInputConcaveScale;
	boolean    filterZ;    // (adjustment mode)filter samples by Z
	double  filterByValueScale;
	int        minLeftSamples; // minimal number of samples (channel/dir/location) for adjustment
    double zMin;
    double zMax;
    double zStep;
    double tMin;
    double tMax;
    double tStep;
	
	double targetRelFocalShift; // target focal shift relative to best composite "sharpness" 

	// when false - tangential is master
	double [] minMeas; // pixels
	double [] maxMeas; // pixels
	double [] thresholdMax; // pixels
	boolean useMinMeas;
	boolean useMaxMeas;
	boolean useThresholdMax;
	int weightMode; // 0; // 0 - same weight, 1 - linear threshold difference, 2 - quadratic thershold difference
	double weightRadius; //2.0; // Gaussian sigma in mm
	private double k_red;
	private double k_blue;
	private double qb_scan_below; // um
	private double qb_scan_above; // um
	private double qb_scan_step; // um
	private boolean qb_use_corrected;
	private boolean qb_invert;
	private boolean z_relative;  // focal distance relative to center greeen
	private boolean rslt_show_z_axial;
	private boolean rslt_show_z_smooth;
	private boolean rslt_show_z_individual;
	private boolean rslt_show_f_axial;
	private boolean rslt_show_f_smooth;
	private boolean rslt_show_f_individual;
	private double  rslt_show_smooth_sigma;
	private double rslt_scan_below;
	private double rslt_scan_above;
	private double rslt_scan_step;
	private boolean rslt_mtf50_mode;
	private boolean rslt_solve; // find z for minimum of f, if false - use parameter z0
	private boolean [] rslt_show_chn;

// not saved/restored
	private double [] z0_estimates=null; // initial estimates for z0 from the lowest value on data
	private double lambdaStepUp; // multiply lambda by this if result is worse
	private double lambdaStepDown; // multiply lambda by this if result is better
	private double thresholdFinish; // (copied from series) stop iterations if 2 last steps had less improvement (but not worsening )
	private int numIterations; // maximal number of iterations
	private double maxLambda; // max lambda to fail
	private double lambda; // copied from series
	private String strategyComment;
	private boolean lastInSeries;
	private int currentStrategyStep; // -1 do not read from strategies
	private boolean stopEachStep; // open dialog after each fitting step
	private boolean stopEachSeries; // stop after each series
	private boolean stopOnFailure; // open dialog when fitting series failed
	private boolean showParams; // show modified parameters
	private boolean showDisabledParams;
	private boolean showCorrectionParams;
	private boolean keepCorrectionParameters;
	private boolean resetCenter; // use distortion center
	private boolean saveSeries; // just for the dialog
	private boolean showMotors;
	private boolean [] showMeasCalc;
	private boolean [] showColors;
	private boolean [] showDirs;
	private boolean [] showSamples;
	private boolean showAllSamples;
	private boolean showIgnoredData;
	private boolean showRad;
	private boolean correct_measurement_ST;
	private boolean updateWeightWhileFitting;
    private int debugPoint;
    private int debugParameter;
    // not reset to defaults
	private boolean [][][][][] sampleMask=null;
	public int debugLevel;
	public boolean debugDerivatives;
	public boolean debugDerivativesFxDxDy;
	private Properties savedProperties=null; // to-be applied
	private String propertiesPrefix=null;
	public double fwhm_to_mtf50=2*Math.log(2.0)/Math.PI*1000; //pi/0.004
    public boolean updateStatus=true;
    public String [] debugParameterNames=null;
	private double [] lastImprovements= {-1.0,-1.0}; // {last improvement, previous improvement}. If both >0 and < thresholdFinish - done
	private int iterationStepNumber=0;
	private long startTime=0;
	private AtomicInteger stopRequested=null; // 1 - stop now, 2 - when convenient
	public static final double PIXEL_SIZE=0.0022; // mm
	public static final String sep = " ";
	public static final String regSep = "\\s";
	public String serialNumber;
	public String lensSerial; // if null - do not add average
	public String comment;
	
	public double [][][] sampleCoord;
	public ArrayList<FocusingFieldMeasurement> measurements;
	double [] weightReference=null; // calculated per-channel (6) array of maximal PSF FWHM after applying min/max correction
	MeasuredSample [] dataVector;
	double [][][] zRanges; // min/max "reliable" z for each channel/sample - will be used during adjustment
	boolean [] prevEnable; // used in adjustment mode to save previous result of filterByZRanges()
//	boolean changedEnable; // used in adjustment mode to signal if new result of filterByZRanges() differes from the previous one
	double [] dataValues;
	double [] dataWeights;
//	int [][][] dataIndex=null; // [measurement][channel][sample] - index in dataValues (and  dataWeights) or -1
	//    double sumWeights=0.0;
	double [][] jacobian=null; // rows - parameters, columns - samples
	double [] currentVector=null;
	double [] nextVector=null;
	double [] savedVector=null;
	
	private LMAArrays lMAArrays=null;
	private LMAArrays savedLMAArrays=null;
	// temporarily changing visibility of currentfX
	double [] currentfX=null; // array of "f(x)" - simulated data for all images, combining pixel-X and pixel-Y (odd/even)
	private double [] nextfX=null; // array of "f(x)" - simulated data for all images, combining pixel-X and pixel-Y (odd/even)
	private double currentRMS=-1.0; // calculated RMS for the currentVector->currentfX
	private double currentRMSPure=-1.0; // calculated RMS for the currentVector->currentfX
	private double nextRMS=-1.0; // calculated RMS for the nextVector->nextfX
	private double nextRMSPure=-1.0; // calculated RMS for the nextVector->nextfX

	private double firstRMS=-1.0; // RMS before current series of LMA started
	private double firstRMSPure=-1.0; // RMS before current series of LMA started


    public void setDefaults(){
    	zRanges=null;
    	prevEnable=null;
//    	changedEnable=true;
    	z0_estimates=null;
    	sagittalMaster=false; // center data is the same, when true sagittal fitting only may change r=0 coefficients,
    	parallelOnly = true; // only process measurements for parallel moves
    	filterInput = true;
    	filterInputMotorDiff = 500.0;
    	filterInputDiff = 2.0; // um
    	filterInputFirstLast = true; 
    	filterInputTooFar = true; // filter samples that are too far from the "center of mass" of other samples 
    	filterInputFarRatio = 3.0; // remove samples that are farther than this ration of average distance
    	filterInputConcave = true; //um
    	filterInputConcaveSigma = 8.0; //um
    	filterInputConcaveRemoveFew=true;
    	filterInputConcaveMinSeries=5;
    	filterInputConcaveScale=0.9;
    	filterZ=true;      // (adjustment mode)filter samples by Z
    	filterByValueScale=1.5; // (adjustment mode)filter samples by value - remove higher than scaled best FWHM
    	minLeftSamples=10;  // minimal number of samples (channel/dir/location) for adjustment

    	zMin=-40.0;
        zMax= 40.0;
        zStep=2.0;
        tMin=-15.0;
        tMax= 15.0;
        tStep=1.0;
    	
    	targetRelFocalShift=8.0; // target focal shift relative to best composite "sharpness" 
    	// when false - tangential is master
    	double [] minMeasDflt= {0.5,0.5,0.5,0.5,0.5,0.5}; // pixels
    	minMeas= minMeasDflt; // pixels
    	double [] maxMeasDflt= {4.5,4.5,4.5,4.5,4.5,4.5}; // pixels
    	maxMeas= maxMeasDflt; // pixels
//    	double [] thresholdMaxDflt= {2.4,3.0,2.6,3.0,3.1,3.0}; // pixels
    	double [] thresholdMaxDflt= {3.5,3.5,3.5,3.5,3.5,3.5}; // pixels
    	thresholdMax= thresholdMaxDflt; // pixels
    	useMinMeas= true;
    	useMaxMeas= true;
    	useThresholdMax=true;
    	weightMode=2; // 1; // 0 - same weight, 1 - linear threshold difference, >1 - power of PSF radius
    	weightRadius=0.0; //2.0; // Gaussian sigma in mm
    	k_red=0.7;
    	k_blue=0.4;
    	qb_scan_below=-40.0; // um
    	qb_scan_above= 80.0; // um
    	qb_scan_step= 0.5; // um
    	qb_use_corrected=true;
    	qb_invert=true;
    	z_relative=true;  // focal distance relative to best overall (false - center green )
    	rslt_show_z_axial=true;
    	rslt_show_z_smooth=false;
    	rslt_show_z_individual=true;
    	rslt_show_f_axial=true;
    	rslt_show_f_smooth=false;
    	rslt_show_f_individual=true;
    	rslt_show_smooth_sigma=0.3; // mm
    	rslt_scan_below=-10.0;
    	rslt_scan_above= 10.0;
    	rslt_scan_step= 5.0;
    	rslt_mtf50_mode= true;
    	rslt_solve = false;
    	boolean [] rslt_show_chnDflt={true,true,true,true,true,true};
    	rslt_show_chn=rslt_show_chnDflt.clone();
    	// not saved/restored
    	lambdaStepUp= 8.0; // multiply lambda by this if result is worse
    	lambdaStepDown= 0.5; // multiply lambda by this if result is better
    	thresholdFinish=0.001; // (copied from series) stop iterations if 2 last steps had less improvement (but not worsening )
    	numIterations= 100; // maximal number of iterations
    	maxLambda= 100.0; // max lambda to fail
    	lambda=0.001; // copied from series
    	stopEachStep= true; // open dialog after each fitting step
    	stopEachSeries= false;
    	stopOnFailure= true; // open dialog when fitting series failed
    	strategyComment="";
    	lastInSeries=true;
    	currentStrategyStep=-1; // -1 do not read from strategies
    	
    	
    	showParams= false; // show modified parameters
    	showDisabledParams = false;
    	showCorrectionParams = false;
    	keepCorrectionParameters = true;
    	resetCenter=false;
    	saveSeries=false; // just for the dialog

    	showMotors = true;
    	boolean [] showMeasCalcDflt={true,true,true};
    	showMeasCalc=showMeasCalcDflt.clone();
    	boolean [] showColorsDflt = {true,true,true};
    	showColors = showColorsDflt.clone();
    	boolean [] showDirsDflt = {true,true};
    	showDirs = showDirsDflt.clone();
    	showSamples = null;
    	showAllSamples = true;
    	showIgnoredData= false;
    	showRad = true;
    	sampleMask=null;

    	correct_measurement_ST=true;
    	updateWeightWhileFitting=false;
        debugPoint=-1;
        debugParameter=-1;
        
        currentPX0=pX0_distortions;
        currentPY0=pY0_distortions;
    	
    }
    
    
	public void setProperties(String prefix,Properties properties){
		if (debugLevel>1) System.out.println("FocusingField: setProperties()");
		if (fieldFitting == null) {
			System.out.println("fieldFitting is not initialized, nothing to save");
			return;
		}

		fieldFitting.setProperties(prefix+"fieldFitting.",properties);
		properties.setProperty(prefix+"pX0_distortions",pX0_distortions+"");
		properties.setProperty(prefix+"pY0_distortions",pY0_distortions+"");
		properties.setProperty(prefix+"currentPX0",currentPX0+"");
		properties.setProperty(prefix+"currentPY0",currentPY0+"");
		properties.setProperty(prefix+"sagittalMaster",sagittalMaster+"");
		properties.setProperty(prefix+"parallelOnly",parallelOnly+"");
		properties.setProperty(prefix+"filterInput",filterInput+"");
		properties.setProperty(prefix+"filterInputMotorDiff",filterInputMotorDiff+"");
		properties.setProperty(prefix+"filterInputDiff",filterInputDiff+"");
		properties.setProperty(prefix+"filterInputFirstLast",filterInputFirstLast+"");
		properties.setProperty(prefix+"filterInputTooFar",filterInputTooFar+"");
		properties.setProperty(prefix+"filterInputFarRatio",filterInputFarRatio+"");
		properties.setProperty(prefix+"filterInputConcave",filterInputConcave+"");
		properties.setProperty(prefix+"filterInputConcaveSigma",filterInputConcaveSigma+"");
		properties.setProperty(prefix+"filterInputConcaveRemoveFew",filterInputConcaveRemoveFew+"");
		properties.setProperty(prefix+"filterInputConcaveMinSeries",filterInputConcaveMinSeries+"");
		properties.setProperty(prefix+"filterInputConcaveScale",filterInputConcaveScale+"");
		properties.setProperty(prefix+"filterZ",filterZ+"");
		properties.setProperty(prefix+"filterByValueScale",filterByValueScale+"");
		
		properties.setProperty(prefix+"minLeftSamples",minLeftSamples+"");
		
		properties.setProperty(prefix+"zMin",zMin+"");
		properties.setProperty(prefix+"zMax",zMax+"");
		properties.setProperty(prefix+"zStep",zStep+"");
		properties.setProperty(prefix+"tMin",tMin+"");
		properties.setProperty(prefix+"tMax",tMax+"");
		properties.setProperty(prefix+"tStep",tStep+"");
		
		properties.setProperty(prefix+"targetRelFocalShift",targetRelFocalShift+"");
		for (int chn=0; chn<minMeas.length; chn++) properties.setProperty(prefix+"minMeas_"+chn,minMeas[chn]+"");
		for (int chn=0; chn<maxMeas.length; chn++) properties.setProperty(prefix+"maxMeas_"+chn,maxMeas[chn]+"");
		for (int chn=0; chn<thresholdMax.length; chn++) properties.setProperty(prefix+"thresholdMax_"+chn,thresholdMax[chn]+"");
		properties.setProperty(prefix+"useMinMeas",useMinMeas+"");
		properties.setProperty(prefix+"useMaxMeas",useMaxMeas+"");
		properties.setProperty(prefix+"useThresholdMax",useThresholdMax+"");
		properties.setProperty(prefix+"weightMode",weightMode+"");
		properties.setProperty(prefix+"weightRadius",weightRadius+"");
		properties.setProperty(prefix+"k_red",k_red+"");
		properties.setProperty(prefix+"k_blue",k_blue+"");
		properties.setProperty(prefix+"qb_scan_below",qb_scan_below+"");
		properties.setProperty(prefix+"qb_scan_above",qb_scan_above+"");
		properties.setProperty(prefix+"qb_scan_step",qb_scan_step+"");
		properties.setProperty(prefix+"qb_use_corrected",qb_use_corrected+"");
		properties.setProperty(prefix+"qb_invert",qb_invert+"");
		properties.setProperty(prefix+"z_relative",z_relative+"");
		properties.setProperty(prefix+"rslt_show_z_axial",rslt_show_z_axial+"");
		properties.setProperty(prefix+"rslt_show_z_smooth",rslt_show_z_smooth+"");
		properties.setProperty(prefix+"rslt_show_z_individual",rslt_show_z_individual+"");
		properties.setProperty(prefix+"rslt_show_f_axial",rslt_show_f_axial+"");
		properties.setProperty(prefix+"rslt_show_f_smooth",rslt_show_f_smooth+"");
		properties.setProperty(prefix+"rslt_show_f_individual",rslt_show_f_individual+"");
		properties.setProperty(prefix+"rslt_show_smooth_sigma",rslt_show_smooth_sigma+"");
		properties.setProperty(prefix+"rslt_scan_below",rslt_scan_below+"");
		properties.setProperty(prefix+"rslt_scan_above",rslt_scan_above+"");
		properties.setProperty(prefix+"rslt_scan_step",rslt_scan_step+"");
		properties.setProperty(prefix+"rslt_mtf50_mode",rslt_mtf50_mode+"");
		properties.setProperty(prefix+"rslt_solve",rslt_solve+"");
		for (int chn=0; chn<rslt_show_chn.length; chn++) properties.setProperty(prefix+"rslt_show_chn_"+chn,rslt_show_chn[chn]+"");
// always re-calculate here?		
		zRanges=calcZRanges(dataWeightsToBoolean());
		if (zRanges!=null){
			properties.setProperty(prefix+"zRanges_length",zRanges.length+"");
			for (int chn=0;chn<zRanges.length;chn++) if (zRanges[chn]!=null) {
				properties.setProperty(prefix+"zRanges_"+chn+"_length",zRanges[chn].length+"");
				for (int sample=0;sample<zRanges[chn].length;sample++) if (zRanges[chn][sample]!=null) {
					properties.setProperty(prefix+"zRanges_"+chn+"_"+sample,zRanges[chn][sample][0]+","+zRanges[chn][sample][1]);
				}
			}
		}
	}

	public void getProperties(String prefix,Properties properties){
		savedProperties=properties;
		propertiesPrefix=prefix;
		if (debugLevel>1) System.out.println("FocusingField: getProperties()");
		if (fieldFitting == null) {
			System.out.println("fieldFitting is not initialized, will apply properties later");
			return; //fieldFitting=new FieldFitting();
		}
		fieldFitting.getProperties(prefix+"fieldFitting.",properties);

		if (properties.getProperty(prefix+"pX0_distortions")!=null)
			pX0_distortions=Double.parseDouble(properties.getProperty(prefix+"pX0_distortions"));
		if (properties.getProperty(prefix+"pY0_distortions")!=null)
			pY0_distortions=Double.parseDouble(properties.getProperty(prefix+"pY0_distortions"));
		if (properties.getProperty(prefix+"currentPX0")!=null)
			currentPX0=Double.parseDouble(properties.getProperty(prefix+"currentPX0"));
		if (properties.getProperty(prefix+"currentPY0")!=null)
			currentPY0=Double.parseDouble(properties.getProperty(prefix+"currentPY0"));
		if (properties.getProperty(prefix+"sagittalMaster")!=null)
			sagittalMaster=Boolean.parseBoolean(properties.getProperty(prefix+"sagittalMaster"));
		if (properties.getProperty(prefix+"parallelOnly")!=null)
			parallelOnly=Boolean.parseBoolean(properties.getProperty(prefix+"parallelOnly"));
		if (properties.getProperty(prefix+"filterInput")!=null)
			filterInput=Boolean.parseBoolean(properties.getProperty(prefix+"filterInput"));
		if (properties.getProperty(prefix+"filterInputMotorDiff")!=null)
			filterInputMotorDiff=Double.parseDouble(properties.getProperty(prefix+"filterInputMotorDiff"));
		if (properties.getProperty(prefix+"filterInputDiff")!=null)
			filterInputDiff=Double.parseDouble(properties.getProperty(prefix+"filterInputDiff"));
		if (properties.getProperty(prefix+"filterInputFirstLast")!=null)
			filterInputFirstLast=Boolean.parseBoolean(properties.getProperty(prefix+"filterInputFirstLast"));
		if (properties.getProperty(prefix+"filterInputTooFar")!=null)
			filterInputTooFar=Boolean.parseBoolean(properties.getProperty(prefix+"filterInputTooFar"));
		if (properties.getProperty(prefix+"filterInputFarRatio")!=null)
			filterInputFarRatio=Double.parseDouble(properties.getProperty(prefix+"filterInputFarRatio"));
		
		if (properties.getProperty(prefix+"filterInputConcave")!=null)
			filterInputConcave=Boolean.parseBoolean(properties.getProperty(prefix+"filterInputConcave"));
		if (properties.getProperty(prefix+"filterInputConcaveSigma")!=null)
			filterInputConcaveSigma=Double.parseDouble(properties.getProperty(prefix+"filterInputConcaveSigma"));
		
		if (properties.getProperty(prefix+"filterInputConcaveRemoveFew")!=null)
			filterInputConcaveRemoveFew=Boolean.parseBoolean(properties.getProperty(prefix+"filterInputConcaveRemoveFew"));
		if (properties.getProperty(prefix+"filterInputConcaveMinSeries")!=null)
			filterInputConcaveMinSeries=Integer.parseInt(properties.getProperty(prefix+"filterInputConcaveMinSeries"));
		if (properties.getProperty(prefix+"filterInputConcaveScale")!=null)
			filterInputConcaveScale=Double.parseDouble(properties.getProperty(prefix+"filterInputConcaveScale"));
		if (properties.getProperty(prefix+"filterZ")!=null)
			filterZ=Boolean.parseBoolean(properties.getProperty(prefix+"filterZ"));
		if (properties.getProperty(prefix+"filterByValueScale")!=null)
			filterByValueScale=Double.parseDouble(properties.getProperty(prefix+"filterByValueScale"));
		if (properties.getProperty(prefix+"minLeftSamples")!=null)
			minLeftSamples=Integer.parseInt(properties.getProperty(prefix+"minLeftSamples"));
		if (properties.getProperty(prefix+"zMin")!=null)
			zMin=Double.parseDouble(properties.getProperty(prefix+"zMin"));
		if (properties.getProperty(prefix+"zMax")!=null)
			zMax=Double.parseDouble(properties.getProperty(prefix+"zMax"));
		if (properties.getProperty(prefix+"zStep")!=null)
			zStep=Double.parseDouble(properties.getProperty(prefix+"zStep"));
		if (properties.getProperty(prefix+"tMin")!=null)
			tMin=Double.parseDouble(properties.getProperty(prefix+"tMin"));
		if (properties.getProperty(prefix+"tMax")!=null)
			tMax=Double.parseDouble(properties.getProperty(prefix+"tMax"));
		if (properties.getProperty(prefix+"tStep")!=null)
			tStep=Double.parseDouble(properties.getProperty(prefix+"tStep"));
		if (properties.getProperty(prefix+"targetRelFocalShift")!=null)
			targetRelFocalShift=Double.parseDouble(properties.getProperty(prefix+"targetRelFocalShift"));
		for (int chn=0; chn<minMeas.length; chn++) if (properties.getProperty(prefix+"minMeas_"+chn)!=null)
			minMeas[chn]=Double.parseDouble(properties.getProperty(prefix+"minMeas_"+chn));
		for (int chn=0; chn<maxMeas.length; chn++) if (properties.getProperty(prefix+"maxMeas_"+chn)!=null)
			maxMeas[chn]=Double.parseDouble(properties.getProperty(prefix+"maxMeas_"+chn));
		for (int chn=0; chn<thresholdMax.length; chn++) if (properties.getProperty(prefix+"thresholdMax_"+chn)!=null)
			thresholdMax[chn]=Double.parseDouble(properties.getProperty(prefix+"thresholdMax_"+chn));
		if (properties.getProperty(prefix+"useMinMeas")!=null)
			useMinMeas=Boolean.parseBoolean(properties.getProperty(prefix+"useMinMeas"));
		if (properties.getProperty(prefix+"useMaxMeas")!=null)
			useMaxMeas=Boolean.parseBoolean(properties.getProperty(prefix+"useMaxMeas"));
		if (properties.getProperty(prefix+"useThresholdMax")!=null)
			useThresholdMax=Boolean.parseBoolean(properties.getProperty(prefix+"useThresholdMax"));
		if (properties.getProperty(prefix+"weightMode")!=null)
			weightMode=Integer.parseInt(properties.getProperty(prefix+"weightMode"));
		if (properties.getProperty(prefix+"weightRadius")!=null)
			weightRadius=Double.parseDouble(properties.getProperty(prefix+"weightRadius"));
		if (properties.getProperty(prefix+"k_red")!=null)
			k_red=Double.parseDouble(properties.getProperty(prefix+"k_red"));
		if (properties.getProperty(prefix+"k_blue")!=null)
			k_blue=Double.parseDouble(properties.getProperty(prefix+"k_blue"));
		if (properties.getProperty(prefix+"qb_scan_below")!=null)
			qb_scan_below=Double.parseDouble(properties.getProperty(prefix+"qb_scan_below"));
		if (properties.getProperty(prefix+"qb_scan_above")!=null)
			qb_scan_above=Double.parseDouble(properties.getProperty(prefix+"qb_scan_above"));
		if (properties.getProperty(prefix+"qb_scan_step")!=null)
			qb_scan_step=Double.parseDouble(properties.getProperty(prefix+"qb_scan_step"));
		if (properties.getProperty(prefix+"qb_use_corrected")!=null)
			qb_use_corrected=Boolean.parseBoolean(properties.getProperty(prefix+"qb_use_corrected"));
		if (properties.getProperty(prefix+"qb_invert")!=null)
			qb_invert=Boolean.parseBoolean(properties.getProperty(prefix+"qb_invert"));
		if (properties.getProperty(prefix+"z_relative")!=null)
			z_relative=Boolean.parseBoolean(properties.getProperty(prefix+"z_relative"));
		if (properties.getProperty(prefix+"rslt_show_z_axial")!=null)
			rslt_show_z_axial=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_show_z_axial"));
		if (properties.getProperty(prefix+"rslt_show_z_smooth")!=null)
			rslt_show_z_smooth=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_show_z_smooth"));
		if (properties.getProperty(prefix+"rslt_show_z_individual")!=null)
			rslt_show_z_individual=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_show_z_individual"));
		if (properties.getProperty(prefix+"rslt_show_f_axial")!=null)
			rslt_show_f_axial=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_show_f_axial"));
		if (properties.getProperty(prefix+"rslt_show_f_smooth")!=null)
			rslt_show_f_smooth=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_show_f_smooth"));
		if (properties.getProperty(prefix+"rslt_show_f_individual")!=null)
			rslt_show_f_individual=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_show_f_individual"));
		if (properties.getProperty(prefix+"rslt_show_smooth_sigma")!=null)
			rslt_show_smooth_sigma=Double.parseDouble(properties.getProperty(prefix+"rslt_show_smooth_sigma"));
		if (properties.getProperty(prefix+"rslt_scan_below")!=null)
			rslt_scan_below=Double.parseDouble(properties.getProperty(prefix+"rslt_scan_below"));
		if (properties.getProperty(prefix+"rslt_scan_above")!=null)
			rslt_scan_above=Double.parseDouble(properties.getProperty(prefix+"rslt_scan_above"));
		if (properties.getProperty(prefix+"rslt_scan_step")!=null)
			rslt_scan_step=Double.parseDouble(properties.getProperty(prefix+"rslt_scan_step"));
		if (properties.getProperty(prefix+"rslt_mtf50_mode")!=null)
			rslt_mtf50_mode=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_mtf50_mode"));
		if (properties.getProperty(prefix+"rslt_solve")!=null)
			rslt_solve=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_solve"));
		for (int chn=0; chn<rslt_show_chn.length; chn++) if (properties.getProperty(prefix+"rslt_show_chn_"+chn)!=null)
			rslt_show_chn[chn]=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_show_chn_"+chn));
		zRanges=null;
		if (properties.getProperty(prefix+"zRanges_length")!=null){
			zRanges=new double [Integer.parseInt(properties.getProperty(prefix+"zRanges_length"))][][];
			for (int chn=0;chn<zRanges.length;chn++) {
				zRanges[chn]=null;
				if (properties.getProperty(prefix+"zRanges_"+chn+"_length")!=null){
					zRanges[chn]=new double [Integer.parseInt(properties.getProperty(prefix+"zRanges_"+chn+"_length"))][];
					for (int sample=0;sample<zRanges[chn].length;sample++) {
						zRanges[chn][sample]=null;
						String s=properties.getProperty(prefix+"zRanges_"+chn+"_"+sample);
						if (s!=null){
							zRanges[chn][sample]=new double[2];
							String [] ss=s.split(",");
							zRanges[chn][sample][0]=Double.parseDouble(ss[0]);
							zRanges[chn][sample][1]=Double.parseDouble(ss[1]);
						}
					}
				}
			}
		}
	}
	public void setDebugLevel(int debugLevel){
		this.debugLevel=debugLevel;
	}
	public void setAdjustMode(boolean mode){
		if ((fieldFitting!=null) && (fieldFitting.mechanicalFocusingModel!=null)) {
			fieldFitting.mechanicalFocusingModel.setAdjustMode(mode);
		}
	}

	public class LMAArrays { // reuse from Distortions?
		public double [][] jTByJ= null; // jacobian multiplied by Jacobian transposed
		public double [] jTByDiff=null; // jacobian multiplied difference vector
		public LMAArrays clone() {
			LMAArrays lma=new LMAArrays();
			lma.jTByJ = this.jTByJ.clone();
			for (int i=0;i<this.jTByJ.length;i++) lma.jTByJ[i]=this.jTByJ[i].clone();
			lma.jTByDiff=this.jTByDiff.clone();
			return lma;
		}
	}

public enum MECH_PAR {
        K0, // Average motor center travel","um/step","0.0124"},
        KD1, // M1 and M2 travel disbalance","um/step","0.0"},
        KD3, // M3 to average of M1 and M2 travel disbalance","um/step","0.0"},
        sM1, // M1: sin component amplitude, relative to tread pitch","","0.0"},
        cM1, // M1: cos component amplitude, relative to tread pitch","","0.0"},
        sM2, // M2: sin component amplitude, relative to tread pitch","","0.0"},
        cM2, // M2: cos component amplitude, relative to tread pitch","","0.0"},
        sM3, // M3: sin component amplitude, relative to tread pitch","","0.0"},
        cM3, // M3: cos component amplitude, relative to tread pitch","","0.0"},
        Lx, // Half horizontal distance between M3 and and M2 supports", "mm","21.0"},
        Ly, // Half vertical distance between M1 and M2 supports", "mm","10.0"},
        mpX0, // pixel X coordinate of mechanical center","px","1296.0"},
        mpY0, // pixel Y coordinate of mechanical center","px","968.0"},
        z0, // center shift, positive away form the lens","um","0.0"},
        tx, // horizontal tilt", "um/pix","0.0"},
        ty // vertical tilt", "um/pix","0.0"}};
}
public static double getPixelMM(){return PIXEL_SIZE;}
public static double getPixelUM(){return PIXEL_SIZE*1000;}
public int flattenIndex(int i, int j){return j+i*sampleCoord[0].length;}
public int getNumSamples(){return sampleCoord.length*sampleCoord[0].length;}
public int getNumChannels(){return 6;}

public int getSampleWidth(){return sampleCoord[0].length;}
public double [][] flattenSampleCoord(){
///sampleCoord
    double [][] flatSampleCoord=new double [sampleCoord.length*sampleCoord[0].length][2];
    int index=0;
    for (int i=0;i<sampleCoord.length;i++) for (int j=0;j<sampleCoord[0].length;j++) flatSampleCoord[index++]= sampleCoord[i][j];
    return flatSampleCoord; // last dimension is not cloned
}
public class MeasuredSample{
    public int [] motors = new int[3];
    public String timestamp;
    public double px;
    public double py;
    public int sampleIndex=0;
    public int channel;
    public double value;
    public double [] dPxyc=new double[2]; // derivative of the value by optical (aberration) center pixel X,Y

    //     public double weight;
//     public MeasuredSample(){}
    public MeasuredSample(
            int [] motors,
            String timestamp,
            double px,
            double py,
            int sampleIndex,
            int channel,
            double value,
            double dPxc,
            double dPyc
            ){
        this.motors = motors;
        this.timestamp=timestamp;
        this.px = px;
        this.py = py;
        this.channel = channel;
        this.value = value;
        this.dPxyc[0]=dPxc;
        this.dPxyc[1]=dPyc;
        this.sampleIndex=sampleIndex;
    }
}
public boolean configureDataVector(String title, boolean forcenew, boolean enableReset){
    if ((fieldFitting == null) && !forcenew){
        forcenew=true;
    }
    boolean setupMasks=true,setupParameters=true,showDisabled=true;
    FieldFitting tmpFieldFitting=fieldFitting;
    if (tmpFieldFitting==null) tmpFieldFitting=    new FieldFitting(); // just to get field description
    int [] numCurvPars=tmpFieldFitting.getNumCurvars();
    GenericDialog gd = new GenericDialog(title+(forcenew?" RESETTING DATA":""));
    gd.addCheckbox("Only use measurements acquired during parallel moves (false - use all)",parallelOnly);
    
    gd.addCheckbox("Remove \"crazy\" input data (samll motor move causing large variations of FWHM)",filterInput);
    gd.addNumericField("Maximal motor move to be considered small",filterInputMotorDiff,0,5,"steps (~90um/step)");
    gd.addNumericField("Maximal allowed PSF FWHM variations fro the move above",filterInputDiff,3,5,"um");
    gd.addCheckbox("Remove first/last in a series of measuremnts separated by small (see above) steps",filterInputFirstLast);
    gd.addCheckbox("Remove measurements taken too far from the rest for the same channel/sample",filterInputTooFar);
    gd.addNumericField("\"Too far\" ratio to the average distance to the center of measurements",filterInputFarRatio,3,5,"um");

    gd.addCheckbox("Filter non-concave areas from best focus for each sample",filterInputConcave);
    gd.addNumericField("Concave filter sigma",filterInputConcaveSigma,3,5,"um");
    gd.addCheckbox("Remove small series ",filterInputConcaveRemoveFew);
    gd.addNumericField("Minimal number of samples (to remove / apply concave filter) ",filterInputConcaveMinSeries,3,5,"samples");
    gd.addNumericField("Concave filter scale",filterInputConcaveScale,3,5,"<=1.0");
    
    gd.addCheckbox("Sagittal channels are master channels (false - tangential are masters)",sagittalMaster);
    gd.addMessage("=== Setting minimal measured PSF radius for different colors/directions ===");
    
        
    for (int i=0;i<minMeas.length;i++){
        gd.addNumericField(tmpFieldFitting.getDescription(i),this.minMeas[i],3,5,"pix");
    }
    gd.addCheckbox("Use minimal radius",useMinMeas);
    gd.addMessage("=== Setting maximal measured PSF radius for different colors/directions ===");
    for (int i=0;i<maxMeas.length;i++){
        gd.addNumericField(tmpFieldFitting.getDescription(i),this.maxMeas[i],3,5,"pix");
    }
    gd.addCheckbox("Use maximal radius",useMaxMeas);
    gd.addMessage("=== Setting maximal usable PSF radius for different colors/directions ===");
    for (int i=0;i<thresholdMax.length;i++){
        gd.addNumericField(tmpFieldFitting.getDescription(i),this.thresholdMax[i],3,5,"pix");
    }
    gd.addCheckbox("Discard measurements with PSF radius above specified above threshold",useThresholdMax);
    if (forcenew) {
        gd.addMessage("=== Setting number of parameters for approximation of the PSF dimensions ===");
        gd.addNumericField("Number of parameters for psf(z) approximation (>=3)",numCurvPars[0],0);
        gd.addNumericField("Number of parameters for radial dependence of PSF curves (>=1)",numCurvPars[1],0);
    }
    gd.addMessage("");
    gd.addNumericField("Data weight mode (0 - equal mode, 1 -linear treshold diff, 2 - quadratic threshold diff)",weightMode,0);
    gd.addNumericField("Data weight radius (multiply weight by Gaussian), 0 - no dependence on radius",weightRadius,3,5,"mm");
    gd.addCheckbox("Setup parameter masks?",setupMasks);
    gd.addCheckbox("Setup parameter values?",setupParameters);
    gd.addCheckbox("Show/modify disabled for auto-adjustment parameters?",showDisabled);
    gd.addCheckbox("Debug feature: update measurements and /dxc, /dyc if center is being fitted",correct_measurement_ST);
    gd.addCheckbox("Debug feature: update sample weights during fitting",updateWeightWhileFitting);
    if (enableReset) gd.enableYesNoCancel("OK","Reset to defaults, re-open"); // default OK (on enter) - "Apply"
    WindowTools.addScrollBars(gd);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    if (!gd.wasOKed()) {
    	savedProperties=null;
    	setDefaults();
    	if (!configureDataVector(title, true,false)) return false;
    	return true;
    }
    
 // boolean configureDataVector(String title, boolean forcenew, boolean moreset)   
    
    parallelOnly=            gd.getNextBoolean();
    filterInput=             gd.getNextBoolean();
    filterInputMotorDiff=    gd.getNextNumber();
    filterInputDiff=         gd.getNextNumber();
    
    filterInputFirstLast=    gd.getNextBoolean();
    filterInputTooFar=       gd.getNextBoolean();
    filterInputFarRatio=     gd.getNextNumber();
	filterInputConcave=      gd.getNextBoolean();
	filterInputConcaveSigma= gd.getNextNumber();
	filterInputConcaveRemoveFew=       gd.getNextBoolean();
	filterInputConcaveMinSeries= (int) gd.getNextNumber();
	filterInputConcaveScale=           gd.getNextNumber();
    
    sagittalMaster= gd.getNextBoolean();
    for (int i=0;i<minMeas.length;i++)this.minMeas[i]= gd.getNextNumber();
    useMinMeas= gd.getNextBoolean();
    for (int i=0;i<maxMeas.length;i++) this.maxMeas[i]= gd.getNextNumber();
    useMaxMeas= gd.getNextBoolean();
    for (int i=0;i<thresholdMax.length;i++) this.thresholdMax[i] = gd.getNextNumber();
    useThresholdMax= gd.getNextBoolean();
    if (forcenew) {
        numCurvPars[0] = (int) gd.getNextNumber();
        numCurvPars[1] = (int) gd.getNextNumber();
    }
    weightMode = (int) gd.getNextNumber();
        weightRadius = gd.getNextNumber();

    setupMasks= gd.getNextBoolean();
    setupParameters= gd.getNextBoolean();
    showDisabled= gd.getNextBoolean();
    correct_measurement_ST= gd.getNextBoolean();
    updateWeightWhileFitting= gd.getNextBoolean();
    if (forcenew) {
        this.fieldFitting= new FieldFitting(
                currentPX0,
                currentPY0,
                numCurvPars[0],
                numCurvPars[1]);
        if (savedProperties!=null){
        	if (debugLevel>0) System.out.println("configureDataVector(): Applying properties");
        	getProperties(propertiesPrefix,savedProperties); // overwrites parallelOnly!
        }
    }
	fieldFitting.setCenterXY(currentPX0,currentPY0);
    if (setupMasks) {
        if (!fieldFitting.maskSetDialog("Setup parameter masks")) return false;
    }
    if (setupParameters) {
        if (!fieldFitting.showModifyParameterValues("Setup parameter values",showDisabled)) return false;
    }
    double [] centerXY=fieldFitting.getCenterXY();
    currentPX0=centerXY[0];
    currentPY0=centerXY[1];
    this.savedVector=fieldFitting.createParameterVector(sagittalMaster);
//     initialVector     
    return true;
}

private boolean [] dataWeightsToBoolean(){
	boolean [] enable = new boolean [dataVector.length];
	for (int i=0;i<enable.length;i++){
		enable[i]=dataWeights[i]>0.0;
	}
	return enable;
}

public int numEnabled(boolean [] en){
	int num=0;
	if (en!=null) for (boolean e:en) if (e) num++;
	return num;
}
private void maskDataWeights(boolean [] enable){
	for (int i=0;i<enable.length;i++){
		if (!enable[i]) dataWeights[i]=0.0;
	}
}
public double [][] getSeriesWeights(){
	double [][] seriesWeights=new double [getNumChannels()][getNumSamples()];
	for (int chn=0;chn<seriesWeights.length;chn++)  for (int sample=0;sample<seriesWeights[chn].length;sample++) seriesWeights[chn][sample]=0.0;
	for (int index=0;index<dataVector.length;index++) if (dataWeights[index]>0.0){
		seriesWeights[dataVector[index].channel][dataVector[index].sampleIndex]+=dataWeights[index];
	}	
	if (debugLevel>1){
		System.out.println("==== getSeriesWeights():");
		for (int chn=0;chn<seriesWeights.length;chn++)  for (int sample=0;sample<seriesWeights[chn].length;sample++){
			System.out.println("chn="+chn+" sample="+sample+" weight="+IJ.d2s(seriesWeights[chn][sample],3));
		}
	}
	return seriesWeights;
}

private double [][][] calcZRanges(
		boolean [] enable){
	double [][][] zRanges=new double[getNumChannels()][getNumSamples()][];
	for (int chn=0;chn<zRanges.length;chn++) for (int sample=0;sample<zRanges[chn].length;sample++) zRanges[chn][sample]=null;
	double [][] sCoord=	flattenSampleCoord();

	for (int index=0;index<dataVector.length;index++) if ((index>=enable.length) ||enable[index]){
		int chn=dataVector[index].channel;
		int sample=dataVector[index].sampleIndex;
		double z=     fieldFitting.getMotorsZ(
				dataVector[index].motors, // 3 motor coordinates
				sCoord[sample][0], // pixel x
				sCoord[sample][1]); // pixel y
		if (zRanges[chn][sample]==null){
			zRanges[chn][sample]=new double[2];
			zRanges[chn][sample][0]=z;
			zRanges[chn][sample][1]=z;
		} else {
			if (z<zRanges[chn][sample][0]) zRanges[chn][sample][0]=z;
			if (z>zRanges[chn][sample][1]) zRanges[chn][sample][1]=z;
		}
	}
	if (debugLevel>0) System.out.println("calcZRanges()");
	return zRanges;
}

private boolean [] filterByZRanges (
		double [][][] zRanges,
		boolean [] enable_in){
	boolean [] enable_out=enable_in.clone();
	double [][] sCoord=	flattenSampleCoord();
	int numFiltered=0;
	int numLeft=0;

	if (zRanges!=null) {
		for (int index=0;index<dataVector.length;index++) if ((index>=enable_in.length) || enable_in[index]){
			int chn=dataVector[index].channel;
			int sample=dataVector[index].sampleIndex;
			double z=     fieldFitting.getMotorsZ(
					dataVector[index].motors, // 3 motor coordinates
					sCoord[sample][0], // pixel x
					sCoord[sample][1]); // pixel y
			if ((zRanges[chn]!=null) && (zRanges[chn][sample]!=null)){
				if ((z<zRanges[chn][sample][0]) || (z>zRanges[chn][sample][1])) {
					enable_out[index]=false;
					numFiltered++;
				} else {
					numLeft++;
				}
			}
		}
	}
	if (debugLevel>1) System.out.println("filterByZRanges(): Filtered "+numFiltered+" samples, left "+numLeft+" samples");
	return enable_out;
}

private boolean [] filterByValue (
		double scale, // scale to best FWHM - larger are ignored
		boolean [] enable_in){
	boolean [] enable_out=enable_in.clone();
	double [][] fwhm = fieldFitting.getFWHM(
    		true, // boolean corrected,
    		true //boolean allChannels 
    		);
	int numFiltered=0;
	int numLeft=0;
	if (scale>0.0) {
		for (int index=0;index<dataVector.length;index++) if ((index>=enable_in.length) || enable_in[index]){
			int chn=dataVector[index].channel;
			int sample=dataVector[index].sampleIndex;
			if (dataVector[index].value > scale*fwhm[chn][sample]){
				enable_out[index]=false;
				numFiltered++;
			} else {
				numLeft++;
			}
		}
	}
	if (debugLevel>0) System.out.println("filterByValue(): Filtered "+numFiltered+" samples, left "+numLeft+" samples");
	return enable_out;
}

private int getNumEnabledSamples(
		boolean [] enable){
	int num_en=0;
	for (int index=0;index<dataVector.length;index++) if ((index>=enable.length) || enable[index]) num_en++;
	return num_en;
}

private boolean [] filterConcave(
		boolean [] scanMask, // do not filter if false
		double sigma,
		boolean removeInsufficient,
		int minSeries,
		double  concaveScale,
		boolean [] enable_in){
	if (enable_in==null) {
		enable_in=new boolean [dataVector.length];
		for (int i=0;i<enable_in.length;i++)enable_in[i]=true;
	}
	if (scanMask==null) {
		scanMask=new boolean [dataVector.length];
		for (int i=0;i<scanMask.length;i++) scanMask[i]=true;
	}
	boolean [] enable_masked=enable_in.clone();
	for (int i=0;i<enable_masked.length;i++) if ((i<scanMask.length) && !scanMask[i]) enable_masked[i]=false;
	boolean [] enable_out=enable_masked.clone();

	
	int debugThreshold=1;
	double maxGap=sigma; // this point has this gap towards minimal
	double kexp=-0.5/(sigma*sigma);
//	boolean [] enable_out=enable_in.clone();
	double keepNearMin=sigma; // when removing non-concave points around min, skip very close ones
	
	double [][] flatSampleCoordinates=fieldFitting.getSampleCoordinates();
	int numFilteredInsufficient = 0;
	int numFiltered = 0;
	int [][] numPoints=new int [getNumChannels()][getNumSamples()];
	double [][][] z0EstData=new double[getNumChannels()][getNumSamples()][2];
	for (int chn=0;chn<numPoints.length;chn++)  for (int sample=0;sample<numPoints[chn].length;sample++) {
		numPoints[chn][sample]=0;
		z0EstData[chn][sample][0]=0.0;
		z0EstData[chn][sample][1]=0.0;
	}
	for (int index=0;index<dataVector.length;index++) if ((index>=enable_masked.length) ||enable_masked[index]){
		numPoints[dataVector[index].channel][dataVector[index].sampleIndex]++;
	}
	int [][][] indices=new int[numPoints.length][numPoints[0].length][];
	for (int chn=0;chn<numPoints.length;chn++)  for (int sample=0;sample<numPoints[chn].length;sample++){
		indices[chn][sample]=new int [numPoints[chn][sample]]; // may be 0 length
		numPoints[chn][sample]=0; // will be used as a counter
	}
	for (int index=0;index<dataVector.length;index++) if ((index>=enable_masked.length) ||enable_masked[index]){
		int chn=dataVector[index].channel;
		int sample=dataVector[index].sampleIndex;
//		numPoints[dataVector[index].channel][dataVector[index].sampleIndex]++;
		indices[chn][sample][numPoints[chn][sample]++]=index;
	}
	for (int chn=0;chn<numPoints.length;chn++) for (int sample=0;sample<numPoints[chn].length;sample++){
		if (indices[chn][sample].length<minSeries){
			if (indices[chn][sample].length>0) {
				if (debugLevel>0) System.out.println("filterConcave(): Channel "+chn+" sample "+sample+" has too few points - "+indices[chn][sample].length+" < "+minSeries);
				if (removeInsufficient){
					for (int i=0;i<indices[chn][sample].length;i++){
						enable_out[indices[chn][sample][i]]=false;
						numFilteredInsufficient++;
					}
				}
			}
		} else {
			int [] thisIndices=indices[chn][sample];
			double [] point_z = new double[thisIndices.length];
			double [] point_v = new double[thisIndices.length];
			double [] point_filt = new double[thisIndices.length];
			double [] point_vdz = new double[thisIndices.length];
			double [] point_slope = new double[thisIndices.length];
			boolean [] nonConcave=new boolean[thisIndices.length];
			for (int i=0;i<thisIndices.length;i++){
//if ((chn==0) && (sample==7) && (i>=16)){
//   System.out.println("DEBUG00");				
//}
				point_z[i]=fieldFitting.getMotorsZ(
						dataVector[thisIndices[i]].motors, // 3 motor coordinates
						flatSampleCoordinates[sample][0], // pixel x
						flatSampleCoordinates[sample][1]); // pixel y
				point_v[i]=dataVector[thisIndices[i]].value;
				nonConcave[i]=false;
			}			
			for (int i=0;i<thisIndices.length;i++){
				point_filt[i]=0.0;
				double weight=0.0;
				for (int j=0;j<thisIndices.length;j++){
					double r=point_z[i]-point_z[j];
					double w=Math.exp(kexp*r*r);
					weight+=w;
					point_filt[i]+=w*point_v[j];
				}
				point_filt[i]/=weight;
			}
			for (int i=0;i<thisIndices.length;i++){
				point_vdz[i]=0.0;
				double S0=0.0,SX=0.0,SY=0.0,SX2=0.0,SXY=0.0;
				for (int j=0;j<thisIndices.length;j++){
					double x=point_z[j]-point_z[i];
					double v=point_filt[j];
					double w=Math.exp(kexp*x*x);
					S0+=w;
					SX+=w*x;
					SY+=w*v;
					SX2+=w*x*x;
					SXY+=w*x*v;
				}
				point_vdz[i]=(SXY*S0-SX*SY)/(SX2*S0-SX*SX);
			}
			// find min on filtered
			int minIndex=0;
			for (int i=1;i<thisIndices.length;i++) if (point_filt[i]<point_filt[minIndex]) minIndex=i;
			for (int i=0;i<thisIndices.length;i++) {
				if (i == minIndex) point_slope[i]=0.0;
				else point_slope[i]= (point_filt[i]-point_filt[minIndex])/(point_z[i]-point_z[minIndex]);
			}
			//concaveScale
			for (int i=0;i<thisIndices.length;i++) {
				if (
						(((point_z[i]-point_z[minIndex])>keepNearMin) && (concaveScale*point_slope[i]>point_vdz[i])) ||
						(((point_z[minIndex]-point_z[i])>keepNearMin) && (concaveScale*point_slope[i]<point_vdz[i]))){
					nonConcave[i]=true;
				}
			}
			// find gaps 	double maxGap=sigma; // this point has this gap towards minimal
			for (int i=0;i<thisIndices.length;i++) if (!nonConcave[i]){
				if ((point_z[i]-point_z[minIndex])>keepNearMin){
					boolean goodPoint=false;
					for (int j=0;j<thisIndices.length;j++) if (
							(point_z[j]>point_z[minIndex]) &&
							(point_z[j]<point_z[i]) &&
							((point_z[i]-point_z[j]) < maxGap )	){
						goodPoint=true;
						break;
					}
					if (!goodPoint) nonConcave[i]=true;
				} else if ((point_z[minIndex]-point_z[i])>keepNearMin) {
					boolean goodPoint=false;
					for (int j=0;j<thisIndices.length;j++) if (
							(point_z[j]<point_z[minIndex]) &&
							(point_z[j]>point_z[i]) &&
							((point_z[j]-point_z[i]) < maxGap )	){
						goodPoint=true;
						break;
					}
					if (!goodPoint) nonConcave[i]=true;
				}
			}
			
			// propagate
			for (int i=0;i<thisIndices.length;i++) if (!nonConcave[i]){
				if ((point_z[i]-point_z[minIndex])>keepNearMin){
					for (int j=0;j<thisIndices.length;j++) if (
							nonConcave[j] &&
							(point_z[j]>point_z[minIndex]) &&
							(point_z[j]<point_z[i])){
						nonConcave[i]=true;
						break;
					}
				} else if ((point_z[minIndex]-point_z[i])>keepNearMin) {
					for (int j=0;j<thisIndices.length;j++) if (
							nonConcave[j] &&
							(point_z[j]<point_z[minIndex]) &&
							(point_z[j]>point_z[i])){
						nonConcave[i]=true;
						break;
					}
				}
			}
			for (int i=0;i<thisIndices.length;i++) if (nonConcave[i]){
				enable_out[thisIndices[i]]=false;
				numFiltered++;
			}
			
			// See if too few are left - remove them
			int numPointsLeft=0;
			for (int i=0;i<thisIndices.length;i++) if (enable_out[thisIndices[i]]){
				numPointsLeft++;
			}
			if (numPointsLeft<minSeries){
				if (debugLevel>0) System.out.println("filterConcave(): Channel "+chn+" sample "+sample+" has too few points left after filter - "+numPointsLeft+" < "+minSeries);
				if (removeInsufficient){
					for (int i=0;i<indices[chn][sample].length;i++){
						if (enable_out[thisIndices[i]]) {
							enable_out[indices[chn][sample][i]]=false;
							numFilteredInsufficient++;
						}
					}
				}
			}
			

			if (debugLevel>debugThreshold) {
				System.out.println("filterConcave(), chn="+chn+", sample="+sample);

				for (int i=0;i<thisIndices.length;i++){
					System.out.println(i+": z="+ IJ.d2s(point_z[i],3)+", v="+ IJ.d2s(point_v[i],3)+
							", filt="+ IJ.d2s(point_filt[i],3)+", vdz="+ IJ.d2s(100*point_vdz[i],3)+
							", slope="+ IJ.d2s(100*point_slope[i],3)+ ", concave="+(nonConcave[i]?0.0:1.0));
				}
			}
			// contribute to z0 calculation
			for (int i=0;i<thisIndices.length;i++) if (!nonConcave[i]){
				z0EstData[chn][sample][1]+=dataWeights[thisIndices[i]];
			}
			z0EstData[chn][sample][0]=point_z[minIndex];
		}
	}
	if (debugLevel>0) System.out.println("filterConcave(): removed for too few points "+numFilteredInsufficient+" samples");
	if (debugLevel>0) System.out.println("filterConcave(): removed for non-concave "+numFiltered+" samples");
	
	//	for (int chn=0;chn<numPoints.length;chn++) for (int sample=0;sample<numPoints[chn].length;sample++){
	z0_estimates=new double[getNumChannels()];
	for (int chn=0;chn<z0_estimates.length;chn++){
		double z=0;
		double w=0;
		for (int sample=0;sample<numPoints[chn].length;sample++){
			z+=z0EstData[chn][sample][0]*z0EstData[chn][sample][1];
			w+=z0EstData[chn][sample][1];
		}
		z0_estimates[chn]= (w>0.0)?z/w:Double.NaN;
	}
	// restore masked out data
	for (int i=0;i<enable_out.length;i++) if (
			(i<scanMask.length) &&
			!scanMask[i]  &&
			enable_in[i]) enable_out[i]=true;
	return enable_out;
}

private boolean [] filterTooFar(
		boolean [] scanMask, // do not filter if false
		double ratio,
		boolean [] enable_in){
	if (enable_in==null) {
		enable_in=new boolean [dataVector.length];
		for (int i=0;i<enable_in.length;i++)enable_in[i]=true;
	}
	if (scanMask==null) {
		scanMask=new boolean [dataVector.length];
		for (int i=0;i<scanMask.length;i++) scanMask[i]=true;
	}
	boolean [] enable_masked=enable_in.clone();
	for (int i=0;i<enable_masked.length;i++) if ((i<scanMask.length) && !scanMask[i]) enable_masked[i]=false;
	boolean [] enable_out=enable_masked.clone();
	int numFiltered = 0;
	double [][][] data=new double [getNumChannels()][getNumSamples()][3];
	double [] z_sample=new double [dataVector.length];
	for (int chn=0;chn<data.length;chn++)
		for (int sample=0;sample<data[chn].length;sample++)
		for (int i=0;i<data[chn][sample].length;i++)data[chn][sample][i]=0;
	for (int index=0;index<dataVector.length;index++) if (((index>=enable_masked.length) ||enable_masked[index])){
		int chn=dataVector[index].channel;
		int sample=dataVector[index].sampleIndex;
		double z= fieldFitting.getMotorsZ(
				dataVector[index].motors, //int [] motors, // 3 motor coordinates
				dataVector[index].px, //        double px, // pixel x
				dataVector[index].px); //        double py)
		double w=weightReference[chn]-dataVector[index].value; // maybe use square of w?
		data[chn][sample][0]+=w;
		data[chn][sample][1]+=w*z;
		data[chn][sample][2]+=w*z*z;
		z_sample[index]=z;
	}
	for (int chn=0;chn<data.length;chn++) for (int sample=0;sample<data[chn].length;sample++) {
		if (data[chn][sample][0]>0.0) {
			double z_av=data[chn][sample][1]/data[chn][sample][0];
			double z2 =(data[chn][sample][2]*data[chn][sample][0]-data[chn][sample][1]*data[chn][sample][1])/
					(data[chn][sample][0]*data[chn][sample][0]);
			data[chn][sample][2]=z2*ratio*ratio; // squared radius
			data[chn][sample][1]=z_av; // center z
			if (debugLevel>1) System.out.println("filterTooFar(): chn="+chn+", sample="+
			sample+", r_av="+IJ.d2s(Math.sqrt(z2),3)+" z_av="+IJ.d2s(z_av,3));
		}
	}
	for (int index=0;index<dataVector.length;index++) if ((index>=enable_masked.length) || enable_masked[index]){
		int chn=dataVector[index].channel;
		int sample=dataVector[index].sampleIndex;
		if (data[chn][sample][0]>0.0) {
			double diff_z=z_sample[index]-data[chn][sample][1];
			if ((diff_z*diff_z > data[chn][sample][2]) &&  enable_out[index]){
				enable_out[index]=false;
				numFiltered++;
			}
		}
	}
	
	if (debugLevel>0) System.out.println("filterTooFar(): removed "+numFiltered+" samples");
	// restore masked out data
	for (int i=0;i<enable_out.length;i++) if (
			(i<scanMask.length) &&
			!scanMask[i]  &&
			enable_in[i]) enable_out[i]=true;
	return enable_out;
}


private boolean [] filterCrazyInput(
		boolean [] scanMask, // do not filter if false
		boolean [] enable_in, // [meas][cjn][sample] (or null) // can be shorter or longer than dataVector
		double maxMotDiff,
		double diff,
		boolean removeFirstLast // very first, very last in all samples (or after big move) - OK
		){
	if (enable_in==null) {
		enable_in=new boolean [dataVector.length];
		for (int i=0;i<enable_in.length;i++)enable_in[i]=true;
	}
	if (scanMask==null) {
		scanMask=new boolean [dataVector.length];
		for (int i=0;i<scanMask.length;i++)scanMask[i]=true;
	}
	boolean [] enable_out=enable_in.clone();
//	int lastIndex=-1;
	int numFiltered=0;
	int [][] lastIndex=new int [getNumChannels()][getNumSamples()];
	int [][] lastTimestampIndex=new int [getNumChannels()][getNumSamples()];
	for (int i=0;i<lastIndex.length;i++) for (int j=0;j<lastIndex[i].length;j++) {
		lastIndex[i][j]=-1;
		lastTimestampIndex[i][j]=-1;
	}
	String lastTimestamp=null;
	int thisTimestampIndex=0;
	int lastIndexAny=-1;
	boolean smallMove=false;
//	for (int index=0;index<dataVector.length;index++) if ((index>=enable_in.length) ||enable_in[index]){
	for (int index=0;index<dataVector.length;index++) if (scanMask[index]){ // crazy neighbor still kills even if is ignored itself - needed
		int chn=dataVector[index].channel;
		int sample=dataVector[index].sampleIndex;
		if (lastTimestamp==null) lastTimestamp=dataVector[index].timestamp;
		if (!dataVector[index].timestamp.equals(lastTimestamp)){
			if (smallMove){ // see if any of the samples/channel did not move last time
				for (int i=0;i<lastIndex.length;i++) for (int j=0;j<lastIndex[i].length;j++) {
					if ((lastIndex[i][j]>=0) && (lastTimestampIndex[i][j]<thisTimestampIndex)){
						if(removeFirstLast){
							if (enable_out[lastIndex[i][j]]) numFiltered++;
							enable_out[lastIndex[i][j]]=false;
						}
						lastIndex[i][j]=-1;
					}
				}
			}
			smallMove=((Math.abs(dataVector[index].motors[0]-dataVector[lastIndexAny].motors[0]))<=maxMotDiff) &&
					((Math.abs(dataVector[index].motors[1]-dataVector[lastIndexAny].motors[1]))<=maxMotDiff) &&
					((Math.abs(dataVector[index].motors[2]-dataVector[lastIndexAny].motors[2]))<=maxMotDiff);
			thisTimestampIndex++;
		}
		// is it a first enabled sample after small move?
		if (smallMove){
			if ((lastIndex[chn][sample]<0) || (lastTimestampIndex[chn][sample]<(thisTimestampIndex-1))){
				if (removeFirstLast) {
					if (enable_out[index]) numFiltered++;
					enable_out[index]=false;
				}
			} else { // large difference?
				if ((Math.abs(dataVector[index].value-dataVector[lastIndex[chn][sample]].value))>=diff){ 
					// yes, remove both this and previous
					if (enable_out[lastIndex[chn][sample]]) numFiltered++;
					enable_out[lastIndex[chn][sample]]=false;
					if (enable_out[index]) numFiltered++;
					enable_out[index]=false;
				}
			}
		}
		lastIndex[chn][sample]=index;
		lastTimestampIndex[chn][sample]=thisTimestampIndex;
		lastTimestamp=dataVector[index].timestamp;
		lastIndexAny=index;
	}
	if (debugLevel>0) System.out.println("filterCrazyInput(): removed "+numFiltered+" samples");
	return enable_out;
}


private int [] getParallelDiff(MeasuredSample [] vector){
	HashMap<Point,AtomicInteger> map=new HashMap<Point,AtomicInteger>();
	for (MeasuredSample ms:vector) {
		Point diff=new Point (ms.motors[1]-ms.motors[0],ms.motors[2]-ms.motors[0]);
		if (map.containsKey(diff)) map.get(diff).incrementAndGet();
		else map.put(diff, new AtomicInteger(1));
	}
	Point parallelDiff=new Point(0,0);
	int maxRun=0;
	for (Point diff:map.keySet()){
		if (map.get(diff).get()>maxRun){
			maxRun=map.get(diff).get();
			parallelDiff=diff;
		}
	}
	if (debugLevel>1) System.out.println("getParallelDiff(): maximal number of parallel measurements is "+
	maxRun+", for M2-M1="+parallelDiff.x+", M3-M2="+parallelDiff.y);
	int [] result= {parallelDiff.x,parallelDiff.y};
	return result;
}
private boolean [] createScanMask(MeasuredSample [] vector){
	int [] diffs=getParallelDiff(vector);
	int longestStart=0;
	int longestRun=0;
	int thisStart=0;
	int thisRun=0;
//	boolean [] chanSel=fieldFitting.getSelectedChannels();
//	int numSamples=0;
	for (int i=0;i<vector.length+1;i++) { // if (chanSel[vector[i].channel]) {
		boolean diffMatch= (i>=vector.length)? false:
			(
					((vector[i].motors[1]-vector[i].motors[0]) == diffs[0]) &&
					((vector[i].motors[2]-vector[i].motors[0]) == diffs[1]));
		if (diffMatch){
			if (thisRun>0){
				thisRun++;
			} else {
				thisStart=i;
				thisRun=1;
			}
		} else {
			if (thisRun>0){
				if (thisRun>longestRun){
					longestRun=thisRun;
					longestStart=thisStart;
				}
				thisRun=0;
			}
		}
//		numSamples++;
	}
//	numSamples--;
	boolean []scanMask=new boolean[vector.length]; // numSamples];
	int index=0;
	for (int i=0;i<vector.length;i++) { //if (chanSel[vector[i].channel]) {
		scanMask[index]=(index>=longestStart) && ((index<(longestStart+longestRun)));
		index++;
	}
	return scanMask;
}

// includes deselected channels
public void setDataVector(
		boolean calibrateMode,
		MeasuredSample [] vector){ // remove unused channels if any. vector is already corrected from input data, FWHM psf
	if (debugLevel>1) System.out.println("+++++ (Re)calculating sample weights +++++");
//	int [] diffs=null;
//	if (calibrateMode && parallelOnly) diffs=getParallelDiff(vector);
	boolean [] chanSel=fieldFitting.getSelectedChannels();
	boolean [] fullScanMask=createScanMask(vector);
	int numSamples=0;
//	int index=0;
	for (int i=0;i<vector.length;i++) if (chanSel[vector[i].channel]){
		
		if (calibrateMode && parallelOnly && !fullScanMask[i]) continue; // skip non-scan
//		if ((diffs!=null) && ( 
//				((vector[i].motors[1]-vector[i].motors[0]) != diffs[0]) ||
//				((vector[i].motors[2]-vector[i].motors[0]) != diffs[1]))) continue;
		numSamples++;
	}
	dataVector=new MeasuredSample [numSamples];
	boolean [] scanMask=new boolean[numSamples];
	int n=0;
	for (int i=0;i<vector.length;i++) if (chanSel[vector[i].channel]) {
		if (calibrateMode && parallelOnly && !fullScanMask[i]) continue;
//		if ((diffs!=null) && ( 
//				((vector[i].motors[1]-vector[i].motors[0]) != diffs[0]) ||
//				((vector[i].motors[2]-vector[i].motors[0]) != diffs[1]))) continue;
		scanMask[n]=fullScanMask[i];
		dataVector[n++]=vector[i];
	}
//	if (calibrateMode) {
//		if (parallelOnly){
//			scanMask=new boolean [numSamples];
//			for (int i=0;i<scanMask.length;i++) scanMask[i]=true;
//		}
//	} else {
//		scanMask=new boolean [numSamples];
//		for (int i=0;i<scanMask.length;i++) scanMask[i]=false;
//	}
	int corrLength=fieldFitting.getNumberOfCorrParameters();
	dataValues = new double [dataVector.length+corrLength];
	dataWeights = new double [dataVector.length+corrLength];
	double kw= (weightRadius>0.0)?(-0.5*getPixelMM()*getPixelMM()/(weightRadius*weightRadius)):0;
	for (int i=0;i<dataVector.length;i++){
		MeasuredSample ms=dataVector[i];
		dataValues[i]=ms.value;
		dataWeights[i]=1.0/Math.pow(ms.value,weightMode);
		if (weightRadius>0.0){
			double r2=(ms.px-currentPX0)*(ms.px-currentPX0)+(ms.py-currentPY0)*(ms.py-currentPY0);
			dataWeights[i]*=Math.exp(kw*r2);
		}
	}
	for (int i=0;i<corrLength;i++){
		dataValues[i+dataVector.length]=0.0; // correction target is always 0
		dataWeights[i+dataVector.length]=1.0; // improve?
	}
	if (calibrateMode && filterInput){
		boolean [] en=dataWeightsToBoolean();
		en= filterCrazyInput(
				scanMask,
				en, // [meas][cjn][sample] (or null) // can be shorter or longer than dataVector
				filterInputMotorDiff,
				filterInputDiff,
				filterInputFirstLast
				);
		maskDataWeights(en);
	}
	if (calibrateMode && filterInputTooFar){
		boolean [] en=dataWeightsToBoolean();
		en= filterTooFar(
				scanMask,
				filterInputFarRatio,
				en);
		maskDataWeights(en);
	}

	if (calibrateMode && filterInputConcave){
		boolean [] en=dataWeightsToBoolean();
		en= filterConcave(
				scanMask,
				filterInputConcaveSigma,
				filterInputConcaveRemoveFew,
				filterInputConcaveMinSeries,
				filterInputConcaveScale,
				en);
		maskDataWeights(en);

	}
// TODO: add filtering for tilt motor calibration
	
	fieldFitting.initSampleCorrVector(
			flattenSampleCoord(), //double [][] sampleCoordinates,
			getSeriesWeights()); //double [][] sampleSeriesWeights);
}

// for compatibility with Distortions class\

public void commitParameterVector(double [] vector){
	fieldFitting.commitParameterVector(vector,sagittalMaster);
	// recalculate measured S,T (depend on center) if center is among fitted parameters
	boolean [] centerSelect=fieldFitting.getCenterSelect();
	if (centerSelect[0] ||centerSelect[1]){ // do not do that if XC, YC are not modified
		// recalculate data vector
		double [] pXY=fieldFitting.getCenterXY();
		if (debugLevel>0) System.out.println("Updated currentPX0="+pXY[0]+"("+currentPX0+")"+", currentPY0="+pXY[1]+"("+currentPY0+")");
		currentPX0=pXY[0];
		currentPY0=pXY[1];
		if (correct_measurement_ST && updateWeightWhileFitting) {
			setDataVector(
					true,
					createDataVector(
					false, // boolean updateSelection,
					pXY[0], //double centerPX,
					pXY[1])); //double centerPY
		}
	}
}

public double [] createFXandJacobian( double [] vector, boolean createJacobian){
    commitParameterVector(vector);
    return createFXandJacobian(createJacobian);
}


public double [] createFXandJacobian(boolean createJacobian){
	int numCorrPar=fieldFitting.getNumberOfCorrParameters();
	double [] fx=new double[dataVector.length + numCorrPar ];
	double [][] derivs=null;
	double [] subData=null;
	boolean [] selChannels=fieldFitting.getSelectedChannels();
	int [] selChanIndices= new int[selChannels.length];
	selChanIndices[0]=0;
	for (int i=1;i<selChanIndices.length;i++){
		selChanIndices[i]= selChanIndices[i-1]+(selChannels[i-1]?1:0);
	}
	int numPars=fieldFitting.getNumberOfParameters(sagittalMaster);
	int numRegPars=fieldFitting.getNumberOfRegularParameters(sagittalMaster);
	if (createJacobian) {
		jacobian=new double [numPars][dataVector.length+numCorrPar];
		for (double [] row : jacobian)
			Arrays.fill(row, 0.0);
		derivs=new double [fieldFitting.getNumberOfChannels()][];
	}
	String prevTimeStamp="";
	double prevPx=-1,prevPy=-1;
	for (int n=0;n<dataVector.length;n++){
		MeasuredSample ms=dataVector[n];
//		int saveDebugLevel=debugLevel;
//		if (n==9) {
//			System.out.println("createFXandJacobian(): n="+n);
//			debugLevel=10;
//		}
		if (!ms.timestamp.equals(prevTimeStamp) || (ms.px!=prevPx) || (ms.py!=prevPy)){
			subData=fieldFitting.getValsDerivatives(
					ms.sampleIndex,
					sagittalMaster,
					ms.motors, // 3 motor coordinates
					ms.px, // pixel x
					ms.py, // pixel y
					derivs);
			prevTimeStamp=ms.timestamp;
			prevPx=ms.px;
			prevPy=ms.py;
		}
//		debugLevel=saveDebugLevel; // restore debugLevel
		fx[n]=subData[selChanIndices[ms.channel]];
		if (createJacobian) {
			double [] thisDerivs=derivs[selChanIndices[ms.channel]];
			//         for (int i=0;i<numRegPars;i++){
			if ((debugLevel>1) && (debugParameter>=0)){
				if ((debugParameter<thisDerivs.length) && (thisDerivs[debugParameter]!=0.0)){
					System.out.println("createFXandJacobian(): n="+n+" channel="+ms.channel+" chn. index="+selChanIndices[ms.channel]+
							", sample="+ms.sampleIndex+" timestamp="+ms.timestamp+" derivative="+thisDerivs[debugParameter]);
				}
				
			}

			// contains derivatives for normal and correction parameters
			for (int i=0;i<numPars;i++){
				//             jacobian[i][n]=derivs[selChanIndices[ms.channel]][i];
				jacobian[i][n]=thisDerivs[i];
			}
			//TODO: correct /dpx, /dpy to compensate for measured S,T calcualtion
			boolean [] centerSelect=fieldFitting.getCenterSelect();
			if (correct_measurement_ST && (centerSelect[0] || centerSelect[1])){ // do not do that if both X and Y are disabled
				int np=0;
				for (int i=0;i<2;i++) if (centerSelect[i]){
					jacobian[np++][n]-=ms.dPxyc[i]; // subtract, as effect is opposite to fX
				}
			}
		}
		// add mutual dependence of correction parameters. first - values (fx)
//		System.out.println("Using sampleCorrVector 1");
		int index=dataVector.length; // add to the end of vector
		if (fieldFitting.sampleCorrChnParIndex!=null){
			int numSamples=getNumSamples();
			for (int chn=0;chn<fieldFitting.sampleCorrChnParIndex.length;chn++) if (fieldFitting.sampleCorrChnParIndex[chn]!=null) {
				for (int np=0;np<fieldFitting.sampleCorrChnParIndex[chn].length;np++){
					int pindex=fieldFitting.sampleCorrChnParIndex[chn][np];
					if (pindex>=0) {
						for (int i=0;i<numSamples;i++){
							double f=0.0;
							for (int j=0;j<numSamples;j++){
								f+=fieldFitting.sampleCorrVector[pindex+j]*fieldFitting.sampleCorrCrossWeights[chn][np][i][j];
							}
							fx[index]=f;
							//                                 f+=fieldFitting.sampleCorrVector[pindex+i]
							if (createJacobian) {
								for (int j=0;j<numSamples;j++){
									jacobian[numRegPars+pindex+j][index]=fieldFitting.sampleCorrCrossWeights[chn][np][i][j];
								}
							}                         
							index++;
						}
					}
				}
			}
		}
	}
	if (createJacobian && (debugLevel>1)){
		if (debugPoint>=0) debugJacobianPoint(debugPoint);
		if (debugParameter>=0) debugJacobianParameter(debugParameter);
	}
	return fx;
}

public void debugJacobianPoint(int nPoint){
	System.out.println("==== Non-zero parameters on which point #"+nPoint+" depends:");
	if (nPoint>=jacobian[0].length){
		System.out.println("Jacobian is defined for "+jacobian[0].length+" points only");
		return;
	}
	for (int i=0;i<jacobian.length;i++){
		if (jacobian[i][nPoint]!=0.0){
			String name=(debugParameterNames==null)?"":debugParameterNames[i];
			System.out.println(i+": "+name+" = "+ jacobian[i][nPoint]);
		}
	}
}
public void debugJacobianParameter(int nPar){
	String name=(debugParameterNames==null)?"":debugParameterNames[nPar];
	System.out.println("==== points that depend on parameter #"+nPar+": "+name+" :");
	if (nPar>=jacobian.length){
		System.out.println("Jacobian is defined for "+jacobian.length+" parameters only");
		return;
	}
	String [] corrNames=fieldFitting.getCorrNames();
	for (int i=0;i<jacobian[nPar].length;i++){
		if (jacobian[nPar][i]!=0.0){
			int nMeasPoints=dataVector.length;
			String pointName="";
			if (i<nMeasPoints){
				pointName="chn"+dataVector[i].channel+":"+dataVector[i].sampleIndex+"__"+dataVector[i].timestamp;
			} else {
				int overData=i-nMeasPoints;
				if ((corrNames!=null) && (overData< corrNames.length)){
					pointName="correction_parameter-"+corrNames[overData];
				} else {
					pointName="unknown-point_+"+overData;
				}
			}
			System.out.println(i+": "+ jacobian[nPar][i]+"    "+pointName);
		}
	}
}

public double getRMS(double [] fx, boolean pure){
    int len=pure?dataVector.length:fx.length;
    double sum=0.0;
    double sum_w=0.0;
    if (dataWeights!=null){
        for (int i=0;i<len;i++){
            double d=fx[i]-dataValues[i];
            sum+=dataWeights[i]*d*d;
            sum_w+=dataWeights[i];
        }
    } else {
        for (int i=0;i<len;i++){
            double d=fx[i]-dataValues[i];
            sum+=d*d;
            sum_w+=1.0;
        }
    }
    if (sum_w>0) {
        sum/=sum_w;
    }
    return Math.sqrt(sum);
}

public MeasuredSample [] createDataVector(FocusingFieldMeasurement measurement){
	 ArrayList<FocusingFieldMeasurement> singleMeasurement=new ArrayList<FocusingFieldMeasurement>();
	 singleMeasurement.add(measurement);
    return createDataVector(
    		singleMeasurement,
    		false, // calibrate
            true,  // update selection
            currentPX0, // ignored
            currentPY0, // ignored
            (this.useMinMeas?this.minMeas:null), // pixels
            (this.useMaxMeas?this.maxMeas:null), // pixels
            (this.useThresholdMax?this.thresholdMax:null)); // pixels
}



public MeasuredSample [] createDataVector(){
    return createDataVector(
            true, // boolean updateSelection,
            currentPX0, // double centerPX,
            currentPY0); //double centerPY
}
public MeasuredSample [] createDataVector(
        boolean updateSelection,
        double centerPX,
        double centerPY
        ){ // use this data
return createDataVector(
		measurements,
		true, // calibrate
        updateSelection,
        centerPX,
        centerPY,
        (this.useMinMeas?this.minMeas:null), // pixels
        (this.useMaxMeas?this.maxMeas:null), // pixels
        (this.useThresholdMax?this.thresholdMax:null)); // pixels
}

/**
* Generate of usable measurement samples
* @param minMeas minimal measurement PSF radius (in pixels) - correction that increases result
* resolution in "sharp" areas (to compensate for measurement error). Individual for each of
* the 6 color/direction components.
* @param updateSelection when true - updates selection of "good" samples, when false - reuses existent one
* @param maxMeas maximal measurement PSF radius (in pixels) - correction that decreases result
* resolution out-of-focus areas to compensate for the limited size of the PSF window.
* Individual for each of the 6 color/direction components.
* @param thresholdMax maximal PSF radius to consider data usable
* @return array of the MeasuredSample instances, including motors, PSF radius, channel and value
*/
/*
Need to find partial derivatives of each of the 3 coefficients: c2, s2 and cs for both px0 and py0
r2 = (x-x0)^2+ (y-y0)^2

c2=(x-x0)^2/r2
s2=(y-y0)^2/r2
cs=(x-x0)*(y-y0)/r2
(p/q)'=(p'*q-q'*p)/q^2

d_r2/d_x0=-2*(x-x0)
d_r2/d_y0=-2*(y-y0)

d_c2/d_x0=(2*(x-x0)*(-1)*r2 - ((d_r2/d_x0)*(x-x0)^2))/r2^2=
(-2*(x-x0)*r2 + (2*(x-x0)*(x-x0)^2))/r2^2=
2*delta_x*(delta_x^2 - r2)/r2^2

d_c2/d_y0= (0-((d_r2/d_y0)*(x-x0)^2))/r2^2=
((2*(y-y0))*(x-x0)^2)/r2^2=
2*delta_y*delta_x^2/r2^2

d_cs/dx0=(-(y-y0)*r2- ((d_r2/d_x0)*(x-x0)*(y-y0)))/r2^2=
(-(y-y0)*r2 + 2*(x-x0)*(x-x0)*(y-y0))/r2^2=
(-delta_y*r2 + 2*delta_x^2*delta_y)/r2^2=
delta_y*(2*delta_x^2-r2)/r2^2


d_c2/d_x0= 2*delta_x*(delta_x^2 - r2)/r2^2
d_c2/d_y0= 2*delta_y*delta_x^2/r2^2

d_s2/d_y0= 2*delta_y*(delta_y^2 - r2)/r2^2
d_s2/d_x0= 2*delta_x*delta_y^2/r2^2

d_cs/dx0= delta_y*(2*delta_x^2-r2)/r2^2
d_cs/dy0= delta_x*(2*delta_y^2-r2)/r2^2

*/
public MeasuredSample [] createDataVector(
		ArrayList<FocusingFieldMeasurement> measurements,
		boolean calibrate, // false - adjust, should have updateSelection==true and a single-element measurements list
		boolean updateSelection,
		double centerPX,
		double centerPY,
		double [] minMeas, // pixels
		double [] maxMeas, // pixels
		double [] thresholdMax){ // pixels
	debugDerivatives=debugLevel==3;
	if (calibrate) {
		currentPX0=centerPX;
		currentPY0=centerPY;
	}
	final int numColors=3;
	final int numDirs=2;
	if (sampleMask== null) updateSelection=true;
	if (updateSelection){
		if (debugLevel>1) System.out.println("createDataVector(true,"+centerPX+","+centerPY+"...)");
		sampleMask= new boolean[measurements.size()] [sampleCoord.length][sampleCoord[0].length][numColors][numDirs];
		for (int n=0;n<sampleMask.length;n++)
			for (int i=0;i<sampleMask[n].length;i++)
				for (int j=0;j<sampleMask[n][i].length;j++)
					for (int c=0;c<numColors;c++)
						for (int d=0;d<numDirs;d++) sampleMask[n][i][j][c][d]=false;

	}
	/*
d_c2/d_x0= 2*delta_x*(delta_x^2 - r2)/r2^2
d_c2/d_y0= 2*delta_y*delta_x^2/r2^2

d_s2/d_y0= 2*delta_y*(delta_y^2 - r2)/r2^2
d_s2/d_x0= 2*delta_x*delta_y^2/r2^2

2*d_cs/dx0= 2*delta_y*(2*delta_x^2-r2)/r2^2
2*d_cs/dy0= 2*delta_x*(2*delta_y^2-r2)/r2^2
	 */
	double [][][] cosSin2Tab=new double[sampleCoord.length][sampleCoord[0].length][3];
	double [][][] debugCosSin2Tab_dx=null;
	double [][][] debugCosSin2Tab_dy=null;
	double debugDelta_x_dx=0,debugDelta_y_dy=0,debugR2_dx=0,debugR2_dy=0;
	if (debugDerivatives){
		debugCosSin2Tab_dx=new double[sampleCoord.length][sampleCoord[0].length][3];
		debugCosSin2Tab_dy=new double[sampleCoord.length][sampleCoord[0].length][3];

	}
	double [][][] cosSin2Tab_dx0=new double[sampleCoord.length][sampleCoord[0].length][3];
	double [][][] cosSin2Tab_dy0=new double[sampleCoord.length][sampleCoord[0].length][3];
	for (int i=0;i<sampleCoord.length;i++) for (int j=0;j<sampleCoord[i].length;j++){
		double delta_x=sampleCoord[i][j][0]-currentPX0;
		double delta_y=sampleCoord[i][j][1]-currentPY0;
		double r2=delta_x*delta_x+delta_y*delta_y;
		double r4=r2*r2;
		if (debugDerivatives){
			debugDelta_x_dx=sampleCoord[i][j][0]-currentPX0-1;
			debugDelta_y_dy=sampleCoord[i][j][1]-currentPY0-1;
			debugR2_dx=debugDelta_x_dx*debugDelta_x_dx+delta_y*delta_y;
			debugR2_dy=delta_x*delta_x+debugDelta_y_dy*debugDelta_y_dy;

		}
		if (r2>0.0) {
			cosSin2Tab[i][j][0]= delta_x*delta_x/r2; // cos^2
			cosSin2Tab[i][j][1]= delta_y*delta_y/r2; // sin^2
			cosSin2Tab[i][j][2]=2*delta_x*delta_y/r2; // 2*cos*sin

			cosSin2Tab_dx0[i][j][0]= 2*delta_x*(delta_x*delta_x - r2)/r4; // d(cos^2)/d(x0)
			cosSin2Tab_dx0[i][j][1]= 2*delta_x*delta_y*delta_y/r4; // d(sin^2)/d(x0)
			cosSin2Tab_dx0[i][j][2]= 2*delta_y*(2*delta_x*delta_x-r2)/r4; // d(2*cos*sin)/d(x0)

			cosSin2Tab_dy0[i][j][0]= 2*delta_y*delta_x*delta_x/r4; // d(cos^2)/d(y0)
			cosSin2Tab_dy0[i][j][1]= 2*delta_y*(delta_y*delta_y - r2)/r4; // d(sin^2)/d(y0)
			cosSin2Tab_dy0[i][j][2]= 2*delta_x*(2*delta_y*delta_y-r2)/r4; // d(2*cos*sin)/d(y0)
			if (debugDerivatives){
				debugCosSin2Tab_dx[i][j][0]= debugDelta_x_dx*debugDelta_x_dx/debugR2_dx; // cos^2
				debugCosSin2Tab_dx[i][j][1]= delta_y*delta_y/debugR2_dx; // sin^2
				debugCosSin2Tab_dx[i][j][2]=2*debugDelta_x_dx*delta_y/debugR2_dx; // 2*cos*sin

				debugCosSin2Tab_dy[i][j][0]= delta_x*delta_x/debugR2_dy; // cos^2
				debugCosSin2Tab_dy[i][j][1]= debugDelta_y_dy*debugDelta_y_dy/debugR2_dy; // sin^2
				debugCosSin2Tab_dy[i][j][2]=2*delta_x*debugDelta_y_dy/debugR2_dy; // 2*cos*sin
			}

		} else {
			cosSin2Tab[i][j][0]=1.0;
			cosSin2Tab[i][j][1]=0.0;
			cosSin2Tab[i][j][2]=0.0;

			cosSin2Tab_dx0[i][j][0]=0.0;
			cosSin2Tab_dx0[i][j][1]=0.0;
			cosSin2Tab_dx0[i][j][2]=0.0;

			cosSin2Tab_dy0[i][j][0]=0.0;
			cosSin2Tab_dy0[i][j][1]=0.0;
			cosSin2Tab_dy0[i][j][2]=0.0;
			if (debugDerivatives){
				debugCosSin2Tab_dx[i][j][0]= 0.0;
				debugCosSin2Tab_dx[i][j][1]= 0.0;
				debugCosSin2Tab_dx[i][j][2]= 0.0;

				debugCosSin2Tab_dy[i][j][0]= 0.0;
				debugCosSin2Tab_dy[i][j][1]= 0.0;
				debugCosSin2Tab_dy[i][j][2]= 0.0;
			}
		}
	}

	ArrayList<MeasuredSample> sampleList=new ArrayList<MeasuredSample>();
	if (thresholdMax != null){
		weightReference=thresholdMax.clone();
		for (int c=0;c<weightReference.length;c++){
			// correct for for minimal measurement;
			if (minMeas != null){
				if (weightReference[c]<minMeas[c]){
					weightReference[c]=0; // do not use
					System.out.println ("Weight reference calculation failed (below minimal), all samples may only have the same weight");
					weightReference=null;
					break;
				}
				weightReference[c]=Math.sqrt(weightReference[c]*weightReference[c]-minMeas[c]*minMeas[c]);
			}
			// correct for for maximal measurement;
			if (maxMeas != null) {
				if (weightReference[c] >= maxMeas[c]){
					weightReference[c]=0; // do not use
					System.out.println ("Weight reference calculation failed (above maximal), all samples may only have the same weight");
					weightReference=null;
					break;
				}
				weightReference[c] = 1.0/Math.sqrt(1.0/(weightReference[c]*weightReference[c])-1.0/(maxMeas[c]*maxMeas[c]));
			}
			// convert to microns from pixels
			weightReference[c]*=getPixelUM();
			if (debugLevel>1) System.out.println("==== weightReference["+c+"]="+weightReference[c]);
		}
	} else {
		thresholdMax=null;
	}
	int nMeas=0;
	for (FocusingFieldMeasurement ffm:measurements){
		double [][][][] samples=ffm.samples;
		if (samples!=null) for (int i=0;i<sampleCoord.length;i++){
			if ((i<samples.length) && (samples[i]!=null)) for (int j=0;j<sampleCoord[i].length;j++){
				if ((j<samples[i].length) && (samples[i][j]!=null)) for (int c=0;c<numColors;c++){

					if ((c<samples[i][j].length) && (samples[i][j][c]!=null)){
						double [] sagTan = {
								Math.sqrt(
										cosSin2Tab[i][j][0]*samples[i][j][c][0]+
										cosSin2Tab[i][j][1]*samples[i][j][c][1]+
										cosSin2Tab[i][j][2]*samples[i][j][c][2]),
										Math.sqrt(
												cosSin2Tab[i][j][1]*samples[i][j][c][0]+
												cosSin2Tab[i][j][0]*samples[i][j][c][1]-
												cosSin2Tab[i][j][2]*samples[i][j][c][2])};
						double [] debugSagTan_dx={0.0,0.0},debugSagTan_dy={0.0,0.0};
						if (debugDerivatives){
							debugSagTan_dx[0]= Math.sqrt(
									debugCosSin2Tab_dx[i][j][0]*samples[i][j][c][0]+
									debugCosSin2Tab_dx[i][j][1]*samples[i][j][c][1]+
									debugCosSin2Tab_dx[i][j][2]*samples[i][j][c][2]);
							debugSagTan_dx[1]= Math.sqrt(
									debugCosSin2Tab_dx[i][j][1]*samples[i][j][c][0]+
									debugCosSin2Tab_dx[i][j][0]*samples[i][j][c][1]-
									debugCosSin2Tab_dx[i][j][2]*samples[i][j][c][2]);
							debugSagTan_dy[0]= Math.sqrt(
									debugCosSin2Tab_dy[i][j][0]*samples[i][j][c][0]+
									debugCosSin2Tab_dy[i][j][1]*samples[i][j][c][1]+
									debugCosSin2Tab_dy[i][j][2]*samples[i][j][c][2]);
							debugSagTan_dy[1]= Math.sqrt(
									debugCosSin2Tab_dy[i][j][1]*samples[i][j][c][0]+
									debugCosSin2Tab_dy[i][j][0]*samples[i][j][c][1]-
									debugCosSin2Tab_dy[i][j][2]*samples[i][j][c][2]);
						}
						double [] sagTan_dx0 = {
								(0.5*(
										cosSin2Tab_dx0[i][j][0]*samples[i][j][c][0]+
										cosSin2Tab_dx0[i][j][1]*samples[i][j][c][1]+
										cosSin2Tab_dx0[i][j][2]*samples[i][j][c][2])/sagTan[0]),
										(0.5*(
												cosSin2Tab_dx0[i][j][1]*samples[i][j][c][0]+
												cosSin2Tab_dx0[i][j][0]*samples[i][j][c][1]-
												cosSin2Tab_dx0[i][j][2]*samples[i][j][c][2])/sagTan[1])
						};

						double [] sagTan_dy0 = {
								(0.5*(
										cosSin2Tab_dy0[i][j][0]*samples[i][j][c][0]+
										cosSin2Tab_dy0[i][j][1]*samples[i][j][c][1]+
										cosSin2Tab_dy0[i][j][2]*samples[i][j][c][2])/sagTan[0]),
										(0.5*(
												cosSin2Tab_dy0[i][j][1]*samples[i][j][c][0]+
												cosSin2Tab_dy0[i][j][0]*samples[i][j][c][1]-
												cosSin2Tab_dy0[i][j][2]*samples[i][j][c][2])/sagTan[1])
						};

						if (debugLevel>3) System.out.print("\n"+ffm.motors[2]+" i= "+i+" j= "+j+" c= "+c+" sagTan= "+sagTan[0]+" "+sagTan[1]+" ");
						for (int d=0;d<numDirs;d++){
							if (debugLevel>3) System.out.print(" d= "+d+" ");

							if (!updateSelection && !sampleMask[nMeas][i][j][c][d]) continue;
							int chn=d+numDirs*c;
							//                     System.out.println("i="+i+", j="+j+", c="+c+", d="+d);
							// saved values are PSF radius, convert to FWHM by multiplying by 2.0
							double value=sagTan[d]*2.0;
							double value_dx0=sagTan_dx0[d]*2.0;
							double value_dy0=sagTan_dy0[d]*2.0;
							double debugValue_dx0=0.0,debugValue_dy0=0.0;
							if (debugDerivatives){
								debugValue_dx0=debugSagTan_dx[d]*2.0;
								debugValue_dy0=debugSagTan_dy[d]*2.0;
							}

							boolean dbg=(debugLevel==3) && (i==1) && (j==3);
							double dbg_delta_x=sampleCoord[i][j][0]-currentPX0;
							double dbg_delta_y=sampleCoord[i][j][1]-currentPY0;

							if (dbg) System.out.print("mot="+ffm.motors[2]+" dx="+dbg_delta_x+" dy="+dbg_delta_y);
							if (dbg) System.out.println(" value="+value+" value_dx0="+value_dx0+" value_dy0="+value_dy0);
							if (dbg) System.out.println(" value(dx)="+debugValue_dx0+" value(dy)="+debugValue_dy0+
									" debugDiff_dx0="+(debugValue_dx0-value)+" debugDiff_dy0="+(debugValue_dy0-value));
							// discard above threshold (in pixels, raw FWHM data)
							if (Double.isNaN(value)) {
								if (debugLevel>3) System.out.println("samples["+i+"]["+j+"]["+c+"]["+d+"] = Double.NaN, motors[0]="+ffm.motors[0]);
								if (updateSelection) continue; // bad measurement
							}
							if (debugLevel>3) System.out.print(" A "+value);
							if (thresholdMax != null){
								if (value >= thresholdMax[chn]) {
									if (debugLevel>3) System.out.print(" > "+thresholdMax[chn]);
									if (updateSelection) continue; // bad measurement (above threshold)
								}
							}
							if (debugLevel>3) System.out.print(" B "+value);
							// correct for for minimal measurement;
							if (minMeas != null){
								if (value<minMeas[chn]) {
									if (debugLevel>3) System.out.print(" < "+minMeas[chn]);
									if (updateSelection) continue; // bad measurement (smaller than correction)
								}
								double f=value;
								value=Math.sqrt(value*value-minMeas[chn]*minMeas[chn]);
								value_dx0*=f/value;
								value_dy0*=f/value;
								if (dbg) {
									System.out.println("2. value="+value+" value_dx0="+value_dx0+" value_dy0="+value_dy0);
									if (debugDerivatives){
										debugValue_dx0=Math.sqrt(debugValue_dx0*debugValue_dx0-minMeas[chn]*minMeas[chn]);
										debugValue_dy0=Math.sqrt(debugValue_dy0*debugValue_dy0-minMeas[chn]*minMeas[chn]);
										if (dbg) System.out.println("2. value(dx)="+debugValue_dx0+" value(dy)="+debugValue_dy0+
												" debugDiff_dx0="+(debugValue_dx0-value)+" debugDiff_dy0="+(debugValue_dy0-value));
									}
								}
							}
							if (debugLevel>3) System.out.print(" C "+value);
							// correct for for maximal measurement;
							if (maxMeas != null) {
								if (value >= maxMeas[chn]){
									if (debugLevel>3) System.out.print(" > "+maxMeas[chn]);
									if (updateSelection) continue; // bad measurement (larger than correction)
								}
								double f=value;
								value = 1.0/Math.sqrt(1.0/(value*value)-1.0/(maxMeas[chn]*maxMeas[chn]));
								//                             value_dx0*=1.0/(value*f*f*f);
								//                             value_dy0*=1.0/(value*f*f*f);
								f=value/f;
								f*=f*f;
								value_dx0*=f;
								value_dy0*=f;

								if (dbg) {
									System.out.println("3. value="+value+" value_dx0="+value_dx0+" value_dy0="+value_dy0);
									if (debugDerivatives){
										debugValue_dx0 = 1.0/Math.sqrt(1.0/(debugValue_dx0*debugValue_dx0)-1.0/(maxMeas[chn]*maxMeas[chn]));
										debugValue_dy0 = 1.0/Math.sqrt(1.0/(debugValue_dy0*debugValue_dy0)-1.0/(maxMeas[chn]*maxMeas[chn]));
										if (dbg) System.out.println("3. value(dx)="+debugValue_dx0+" value(dy)="+debugValue_dy0+
												" debugDiff_dx0="+(debugValue_dx0-value)+" debugDiff_dy0="+(debugValue_dy0-value));
									}
								}
							}
							if (debugLevel>3) System.out.print(" D "+value);
							// convert to microns from pixels
							value*=getPixelUM();
							value_dx0*=getPixelUM();
							value_dy0*=getPixelUM();
							if (dbg) {
								System.out.println("4. value="+value+" value_dx0="+value_dx0+" value_dy0="+value_dy0);
								if (debugDerivatives){
									debugValue_dx0*=getPixelUM();
									debugValue_dy0*=getPixelUM();
									if (dbg) System.out.println("4. value(dx)="+debugValue_dx0+" value(dy)="+debugValue_dy0+
											" debugDiff_dx0="+(debugValue_dx0-value)+" debugDiff_dy0="+(debugValue_dy0-value));
								}
							}

							sampleList.add(new MeasuredSample(
									ffm.motors,
									ffm.timestamp,
									sampleCoord[i][j][0], // px,
									sampleCoord[i][j][1], // py,
									flattenIndex(i,j),
									chn,
									value,
									value_dx0, //double dPxc; // derivative of the value by optical (aberration) center pixel X
									value_dy0 //double dPyc; // derivative of the value by optical (aberration) center pixel Y
									));
							if (debugLevel>3) System.out.print(" E "+value);
							if (updateSelection) sampleMask[nMeas][i][j][c][d]=true;
						}
					}
				}
			}
		}
		nMeas++;
	}
	if (debugLevel>3) System.out.println();
	if (debugLevel>1) System.out.println("createDataVector -> "+sampleList.size()+" elements");
	return sampleList.toArray(new MeasuredSample[0]);
}
    /**
     * Calculate differences vector
     * @param fX vector of calculated pixelX,pixelY on the sensors
     * @return same dimension vector of differences from this.Y (measured grid pixelxX, pixelY)
     */
    public double [] calcYminusFx(double [] fX){
        double [] result=this.dataValues.clone();
        for (int i=0;i<result.length;i++) result[i]-=fX[i];
     return result;    
    }
    public double calcErrorDiffY(double [] fX, boolean pure){
        int len=pure?dataVector.length:fX.length;
        double result=0;
        double sumWeights=0;
        if (this.dataWeights!=null) {
            for (int i=0;i<len;i++){ //
                double diff=this.dataValues[i]-fX[i];
                result+=diff*diff*this.dataWeights[i];
                if (debugLevel>2) System.out.println(""+i+" fx="+fX[i]+" data="+this.dataValues[i]+" diff="+diff+" w="+this.dataWeights[i]);
                if (debugLevel>1){
                	if (pure) {
                		int chn=dataVector[i].channel;
                		int sample=dataVector[i].sampleIndex;
                		if ((chn==2) && (sample==16)) {
                			System.out.println(""+i+" fx="+fX[i]+" data="+this.dataValues[i]+" diff="+diff+" w="+this.dataWeights[i]);
                		}
                	}
                }
                sumWeights+=this.dataWeights[i];
            }
            if (sumWeights>0) result/=sumWeights;
        } else {
            for (int i=0;i<len;i++){
                double diff=this.dataValues[i]-fX[i];
                result+=diff*diff;
            }
            result/=fX.length;
        }
        return Math.sqrt(result);
    }
    

    public LMAArrays calculateJacobianArrays(double [] fX){
        // calculate JtJ
        double [] diff=calcYminusFx(fX);
        int numPars=this.jacobian.length; // number of parameters to be adjusted
        int length=diff.length; // should be the same as this.jacobian[0].length
     double [][] JtByJmod=new double [numPars][numPars]; //Transposed Jacobian multiplied by Jacobian
     double [] JtByDiff=new double [numPars];
     for (int i=0;i<numPars;i++) for (int j=i;j<numPars;j++){
         JtByJmod[i][j]=0.0;
         if (this.dataWeights!=null)
             for (int k=0;k<length;k++) JtByJmod[i][j]+=this.jacobian[i][k]*this.jacobian[j][k]*this.dataWeights[k];
         else
             for (int k=0;k<length;k++) JtByJmod[i][j]+=this.jacobian[i][k]*this.jacobian[j][k];
     }
     for (int i=0;i<numPars;i++) { // subtract lambda*diagonal , fill the symmetrical half below the diagonal
//         JtByJmod[i][i]+=lambda*JtByJmod[i][i]; //Marquardt mod
         for (int j=0;j<i;j++) JtByJmod[i][j]= JtByJmod[j][i]; // it is symmetrical matrix, just copy
     }
     for (int i=0;i<numPars;i++) {
         JtByDiff[i]=0.0;
         if (this.dataWeights!=null)
             for (int k=0;k<length;k++) JtByDiff[i]+=this.jacobian[i][k]*diff[k]*this.dataWeights[k];
         else
             for (int k=0;k<length;k++) JtByDiff[i]+=this.jacobian[i][k]*diff[k];

     }
    
        LMAArrays lMAArrays = new LMAArrays();
        lMAArrays.jTByJ=JtByJmod;
        lMAArrays.jTByDiff=JtByDiff;
        return lMAArrays;
    }
    public double [] solveLMA(
            LMAArrays lMAArrays,
            double lambda,
            int debugLevel){
        this.debugLevel=debugLevel;
        double [][] JtByJmod= lMAArrays.jTByJ.clone();
        int numPars=JtByJmod.length;
        for (int i=0;i<numPars;i++){
            JtByJmod[i]=lMAArrays.jTByJ[i].clone();
            JtByJmod[i][i]+=lambda*JtByJmod[i][i]; //Marquardt mod
        }
//     M*Ma=Mb
     Matrix M=new Matrix(JtByJmod);
        if (debugLevel>2) {
            System.out.println("Jt*J -lambda* diag(Jt*J), lambda="+lambda+":");
            M.print(10, 5);
        }

     Matrix Mb=new Matrix(lMAArrays.jTByDiff,numPars); // single column
     if (!(new LUDecomposition(M)).isNonsingular()){
         double [][] arr=M.getArray();
            System.out.println("Singular Matrix "+arr.length+"x"+arr[0].length);
            // any rowsx off all 0.0?
            for (int n=0;n<arr.length;n++){
                boolean zeroRow=true;
                for (int i=0;i<arr[n].length;i++) if (arr[n][i]!=0.0){
                    zeroRow=false;
                    break;
                }
                if (zeroRow){
                    System.out.println("Row of all zeros: "+n);
                }
            }
//            M.print(10, 5);
         return null;
     }
     Matrix Ma=M.solve(Mb); // singular
     return Ma.getColumnPackedCopy();
    }

    public void compareDrDerivatives(double [] vector){
    	double delta=0.00010; // make configurable
    	boolean [] centerSelect=fieldFitting.getCenterSelect();
    	if (centerSelect[0] && centerSelect[1]) delta=1.0;
    	if (debugParameter>=0){
    		String parName="";
    		if ((debugParameterNames!=null) && (debugParameterNames.length>debugParameter)) parName=debugParameterNames[debugParameter];
    		System.out.println("Debugging derivatives for parameter #"+debugParameter+" ("+parName+")");
    		//debugParameterNames
            double [] vector_dp=vector.clone();
            vector_dp[debugParameter]+=delta;
            double [] fx_dp=createFXandJacobian(vector_dp,false);
            double [] fx= createFXandJacobian(vector,true);
            for (int i=0;i<fx.length;i++){
            	if ((debugPoint>=0) && (debugPoint!=i)) continue; // debug only single point
    			int nMeasPoints=dataVector.length;
    			String pointName="";
    			if (i<nMeasPoints){
    				pointName="chn"+dataVector[i].channel+":"+dataVector[i].sampleIndex+"__"+dataVector[i].timestamp;
//    			} else {
//    				int overData=i-nMeasPoints;
//    				if ((corrNames!=null) && (overData< corrNames.length)){
//    					pointName="correction_parameter-"+corrNames[overData];
//    				} else {
//    					pointName="unknown-point_+"+overData;
//    				}
    			}
                System.out.println(i+": "+pointName+" fx= "+fx[i]+" delta_fx= "+((fx_dp[i]-fx[i])/delta)+" df/dp= "+jacobian[debugParameter][i]);
            }
    	}
/*   	
        boolean [] centerSelect=fieldFitting.getCenterSelect();
        if (!centerSelect[0] || !centerSelect[1]){
            System.out.println("compareDrDerivatives(): Both px0 and px1 parameters should be enabled, aborting");
            return;
        }
        double [] vector_dx=vector.clone();
        double [] vector_dy=vector.clone();
        vector_dx[0]+=1.0;
        vector_dy[1]+=1.0;
     double [] fx_dx=createFXandJacobian(vector_dx,false);
     double [] fx_dy=createFXandJacobian(vector_dy,false);
     double [] fx= createFXandJacobian(vector,true);
     for (int i=0;i<fx.length;i++){
         System.out.println(i+" fx= "+fx[i]+" delta_fx= "+(fx_dx[i]-fx[i])+" delta_fy= "+(fx_dy[i]-fx[i]));
         System.out.println(" df/dx= "+jacobian[0][i]+" df/dy= "+jacobian[1][i]);
         
     }
*/     
    }
    
    
    /**
     * Calculates next parameters vector, holds some arrays
     * @param numSeries
     * @return array of two booleans: { improved, finished}
     */
    public boolean [] stepLevenbergMarquardtFirst(int debugLevel){
        double [] deltas=null;
        if (this.currentVector==null) {
            this.currentVector=this.savedVector.clone();
            this.currentRMS=-1;
            this.currentfX=null; // invalidate
            this.jacobian=null; // invalidate
            this.lMAArrays=null;
            lastImprovements[0]=-1.0;
            lastImprovements[1]=-1.0;
        }
        this.debugLevel=debugLevel;
        // calculate this.currentfX, this.jacobian if needed
        if (debugLevel>2) {
            System.out.println("this.currentVector");
            for (int i=0;i<this.currentVector.length;i++){
                System.out.println(i+": "+ this.currentVector[i]);
            }
        }
        //     if ((this.currentfX==null)|| ((this.jacobian==null) && !this.threadedLMA )) {
        if ((this.currentfX==null)|| (this.lMAArrays==null)) {
			String msg="initial Jacobian matrix calculation. Points:"+this.dataValues.length+" Parameters:"+this.currentVector.length;
			if (debugLevel>1) System.out.println(msg);
			if (this.updateStatus) IJ.showStatus(msg);
            this.currentfX=createFXandJacobian(this.currentVector, true); // is it always true here (this.jacobian==null)
            this.lMAArrays=calculateJacobianArrays(this.currentfX);
            this.currentRMS= calcErrorDiffY(this.currentfX,false);
            this.currentRMSPure=calcErrorDiffY(this.currentfX, true);
			msg=this.currentStrategyStep+": initial RMS="+IJ.d2s(this.currentRMS,8)+" (pure RMS="+IJ.d2s(this.currentRMSPure,8)+")"+
					". Calculating next Jacobian. Points:"+this.dataValues.length+" Parameters:"+this.currentVector.length;
			if (debugLevel>1) System.out.println(msg);
			if (this.updateStatus) IJ.showStatus(msg);
        } else {
            // continuing, but need to re-build Y vector to match this.currentVector
            commitParameterVector(this.currentVector); // may be extra job?
            this.currentRMS= calcErrorDiffY(this.currentfX, false); // Verify - currentfX is correct, but Y vector is modified! Yes, it is so
            this.currentRMSPure=calcErrorDiffY(this.currentfX, true);
    }
    if (this.firstRMS<0) {
        this.firstRMS=this.currentRMS;
        this.firstRMSPure=this.currentRMSPure;
    }
    
    
        deltas=solveLMA(this.lMAArrays,    this.lambda, debugLevel);
        
    boolean matrixNonSingular=true;
    if (deltas==null) {
        deltas=new double[this.currentVector.length];
        for (int i=0;i<deltas.length;i++) deltas[i]=0.0;
        matrixNonSingular=false;
    }
        if (debugLevel>1) {
            System.out.println("deltas");
            for (int i=0;i<deltas.length;i++){
                System.out.println(i+": "+ deltas[i]);
            }
        }
// apply deltas     
    this.nextVector=this.currentVector.clone();
    for (int i=0;i<this.nextVector.length;i++) this.nextVector[i]+=deltas[i];
// another option - do not calculate J now, just fX. and late - calculate both if it was improvement     
//     save current Jacobian
    
        if (debugLevel>2) {
            System.out.println("this.nextVector");
            for (int i=0;i<this.nextVector.length;i++){
                System.out.println(i+": "+ this.nextVector[i]);
            }
        }

// this.savedJacobian=this.jacobian;
this.savedLMAArrays=lMAArrays.clone();
this.jacobian=null; // not needed, just to catch bugs
this.nextfX=createFXandJacobian(this.nextVector,true);
this.lMAArrays=calculateJacobianArrays(this.nextfX);
this.nextRMS= calcErrorDiffY(this.nextfX,false);
this.nextRMSPure= calcErrorDiffY(this.nextfX,true);
        this.lastImprovements[1]=this.lastImprovements[0];
        this.lastImprovements[0]=this.currentRMS-this.nextRMS;
		String msg="currentRMS="+this.currentRMS+
				", currentRMSPure="+this.currentRMSPure+
				", nextRMS="+this.nextRMS+
				", nextRMSPure="+this.nextRMSPure+
				", delta="+(this.currentRMS-this.nextRMS);
		if (debugLevel>2) System.out.println("stepLMA "+msg);
		if (this.updateStatus) IJ.showStatus(msg);
        boolean [] status={matrixNonSingular && (this.nextRMS<=this.currentRMS),!matrixNonSingular};
        // additional test if "worse" but the difference is too small, it was be caused by computation error, like here:
        //stepLevenbergMarquardtAction() step=27, this.currentRMS=0.17068403807026408, this.nextRMS=0.1706840380702647
        
        if (!status[0] && matrixNonSingular) {
            if (this.nextRMS<(this.currentRMS+this.currentRMS*this.thresholdFinish*0.01)) {
                this.nextRMS=this.currentRMS;
                status[0]=true;
                status[1]=true;
                this.lastImprovements[0]=0.0;
                if (debugLevel>1) {
                    System.out.println("New RMS error is larger than the old one, but the difference is too small to be trusted ");
                    System.out.println(
                            "stepLMA this.currentRMS="+this.currentRMS+
                            ", this.currentRMSPure="+this.currentRMSPure+
                            ", this.nextRMS="+this.nextRMS+
                            ", this.nextRMSPure="+this.nextRMSPure+
                            ", delta="+(this.currentRMS-this.nextRMS));
                }
                
            }
        }
    if (status[0] && matrixNonSingular) { //improved
        status[1]=(this.iterationStepNumber>this.numIterations) || ( // done
                (this.lastImprovements[0]>=0.0) &&
                (this.lastImprovements[0]<this.thresholdFinish*this.currentRMS) &&
                (this.lastImprovements[1]>=0.0) &&
                (this.lastImprovements[1]<this.thresholdFinish*this.currentRMS));
    } else if (matrixNonSingular){
//         this.jacobian=this.savedJacobian;// restore saved Jacobian
        this.lMAArrays=this.savedLMAArrays; // restore Jt*J and Jt*diff
        
        status[1]=(this.iterationStepNumber>this.numIterations) || // failed
        ((this.lambda*this.lambdaStepUp)>this.maxLambda);
    }
///this.currentRMS     
//TODO: add other failures leading to result failure?     
        if (debugLevel>2) {
            System.out.println("stepLevenbergMarquardtFirst("+debugLevel+")=>"+status[0]+","+status[1]);
        }
        return status;
}
    
public void stepLevenbergMarquardtAction(int debugLevel){//
	this.iterationStepNumber++;
	// apply/revert,modify lambda 
	String msg="currentRMS="+this.currentRMS+
			", currentRMSPure="+this.currentRMSPure+
			", nextRMS="+this.nextRMS+
			", nextRMSPure="+this.nextRMSPure+
			", delta="+(this.currentRMS-this.nextRMS)+
			", lambda="+this.lambda;
	if (debugLevel>1) System.out.println("stepLevenbergMarquardtAction() "+msg);
//	if (this.updateStatus) IJ.showStatus(msg);
	if (this.nextRMS<this.currentRMS) { //improved
		this.lambda*=this.lambdaStepDown;
		this.currentRMS=this.nextRMS;
		this.currentfX=this.nextfX;
		this.currentVector=this.nextVector;
	} else {
		this.lambda*=this.lambdaStepUp;
		this.lMAArrays=this.savedLMAArrays; // restore Jt*J and Jt*diff
	}
}

/**
* Dialog to select Levenberg-Marquardt algorithm and related parameters
* @return true if OK, false if canceled
*/
public boolean selectLMAParameters(){
//     int numSeries=fittingStrategy.getNumSeries();
//    boolean resetCorrections=false;
     GenericDialog gd = new GenericDialog("Levenberg-Marquardt algorithm parameters for cameras distortions/locations");
     
     //TODO: change to selection using series comments
//     	gd.addNumericField("Fitting series number", this.currentStrategyStep, 0, 3," (-1 - current)");
    	FieldStrategies fs=fieldFitting.fieldStrategies;
        String [] indices=new String[fs.getNumStrategies()+1];
        indices[0]="current strategy";
        for (int i=0;i<fs.getNumStrategies();i++) {
        	indices[i+1]=i+": "+fs.getComment(i)+" ("+(fs.isStopAfterThis(i)?"STOP":"CONTINUE")+")";
        }
        if (this.currentStrategyStep>=(indices.length-1)) this.currentStrategyStep=indices.length-2;
        gd.addChoice("Fitting series", indices,indices[this.currentStrategyStep+1]);
     	
        gd.addCheckbox("Debug df/dX0, df/dY0", false);
        gd.addNumericField("Debug Jacobian for point number", this.debugPoint, 0, 5,"(-1 - none)");
        gd.addNumericField("Debug Jacobian for parameter number", this.debugParameter, 0, 5,"(-1 - none)");
        
//        gd.addCheckbox("Keep current correction parameters (do not reset)", this.keepCorrectionParameters);
        gd.addNumericField("Initial LMA Lambda ", 0.0, 5, 8, "0 - keep, last was "+this.lambda);
        gd.addNumericField("Multiply lambda on success", this.lambdaStepDown, 5);
        gd.addNumericField("Threshold RMS to exit LMA", this.thresholdFinish, 7,9,"pix");
        gd.addNumericField("Multiply lambda on failure", this.lambdaStepUp, 5);
        gd.addNumericField("Threshold lambda to fail", this.maxLambda, 5);
        gd.addNumericField("Maximal number of iterations", this.numIterations, 0);

        gd.addCheckbox("Dialog after each iteration step", this.stopEachStep);
        gd.addCheckbox("Dialog after each iteration series", this.stopEachSeries);
        gd.addCheckbox("Dialog after each failure", this.stopOnFailure);
        gd.addCheckbox("Show modified parameters", this.showParams);
        gd.addCheckbox("Show disabled parameters", this.showDisabledParams);
        gd.addCheckbox("Show per-sample correction parameters", this.showCorrectionParams);

        //        gd.addCheckbox("Reset all per-sample corrections to zero", resetCorrections);

        //        gd.addCheckbox("Show debug images before correction",this.showThisImages);
        //        gd.addCheckbox("Show debug images after correction", this.showNextImages);
        //        gd.addNumericField("Maximal number of threads", this.threadsMax, 0);
        //        gd.addCheckbox("Use memory-saving/multithreaded version", this.threadedLMA);
        gd.showDialog();
        if (gd.wasCanceled()) return false;
        this.currentStrategyStep= gd.getNextChoiceIndex()-1; //(int) gd.getNextNumber();
        if (this.currentStrategyStep>=0){
        	getStrategy(this.currentStrategyStep);
        }
        this.debugDerivativesFxDxDy=gd.getNextBoolean();

        debugPoint=     (int) gd.getNextNumber();
        debugParameter= (int) gd.getNextNumber();

//        this.keepCorrectionParameters = gd.getNextBoolean();
        double preLambda=gd.getNextNumber();
        if (preLambda>0.0) this.lambda= preLambda;
        this.lambdaStepDown= gd.getNextNumber();
        this.thresholdFinish= gd.getNextNumber();
        this.lambdaStepUp= gd.getNextNumber();
        this.maxLambda= gd.getNextNumber();
        this.numIterations= (int) gd.getNextNumber();
        this.stopEachStep= gd.getNextBoolean();
        this.stopEachSeries= gd.getNextBoolean();
        this.stopOnFailure= gd.getNextBoolean();
        this.showParams= gd.getNextBoolean();
        this.showDisabledParams= gd.getNextBoolean();
        this.showCorrectionParams= gd.getNextBoolean();
//    if (!keepCorrectionParameters) fieldFitting.resetSampleCorr();
     return true;
}

public void listParameters(String title, String path){
    String header="State\tDescription\tValue\tUnits";
        StringBuffer sb = new StringBuffer();
        for (String s:this.fieldFitting.getParameterValueStrings(true,true)){ //this.showDisabledParams)){ //
        sb.append(s+"\n");
        }
        if (path!=null) {
            CalibrationFileManagement.saveStringToFile (
                    path,
                    header+"\n"+sb.toString());
        } else {
            new TextWindow(title, header, sb.toString(), 800,1000);
        }
}

public void listData(String title, String path){
    if ((showSamples==null) || (showSamples.length!=sampleCoord.length*sampleCoord[0].length)){
        showSamples=new boolean [sampleCoord.length*sampleCoord[0].length];
        for (int i=0;i<showSamples.length;i++) showSamples[i]=false;
    }
     GenericDialog gd = new GenericDialog(title);
        gd.addCheckbox("Show motor positions", this.showMotors);
        gd.addCheckbox("Show measured PSF FWHM", this.showMeasCalc[0]);
        gd.addCheckbox("Show calculated PSF FWHM", this.showMeasCalc[1]);
        gd.addCheckbox("Show mechanical distance", this.showMeasCalc[2]);
        gd.addCheckbox("Show RED data", this.showColors[0]);
        gd.addCheckbox("Show GREEN data", this.showColors[1]);
        gd.addCheckbox("Show BLUE data", this.showColors[2]);
        gd.addCheckbox("Show SAGITAL data", this.showDirs[0]);
        gd.addCheckbox("Show TANGENTIAL data", this.showDirs[1]);
    gd.addMessage("Select field samples");
    for (int i=0; i<sampleCoord.length;i++) for (int j=0; j<sampleCoord[0].length;j++){
        gd.addCheckbox("Sample X"+j+"Y"+i+" (x="+sampleCoord[i][j][0]+
        " y="+sampleCoord[i][j][1]+")",showSamples[j+i*sampleCoord[0].length]);
    }
    
        gd.addCheckbox("Show all samples", this.showAllSamples);
        gd.addCheckbox("Show ignored data", this.showIgnoredData);
        
        gd.addCheckbox("Show sample distance fom the optical center", this.showRad);
        WindowTools.addScrollBars(gd);
     gd.showDialog();
     if (gd.wasCanceled()) return;
        this.showMotors=gd.getNextBoolean();
        this.showMeasCalc[0]=gd.getNextBoolean();
        this.showMeasCalc[1]=gd.getNextBoolean();
        this.showMeasCalc[2]=gd.getNextBoolean();
        this.showColors[0]=gd.getNextBoolean();
        this.showColors[1]=gd.getNextBoolean();
        this.showColors[2]=gd.getNextBoolean();
        this.showDirs[0]=gd.getNextBoolean();
        this.showDirs[1]=gd.getNextBoolean();
    for (int i=0; i<sampleCoord.length;i++) for (int j=0; j<sampleCoord[0].length;j++){
        showSamples[j+i*sampleCoord[0].length]=gd.getNextBoolean();
    }
        this.showAllSamples=gd.getNextBoolean();
        this.showIgnoredData=gd.getNextBoolean();
        this.showRad=gd.getNextBoolean();
     boolean [] localShowSamples=showSamples.clone();
     if (showAllSamples){
         for (int i=0;i<localShowSamples.length;i++) localShowSamples[i]=true;
     }
     listData(title,
    		 path,
    		 showMotors,
    		 showMeasCalc,
    		 showColors,
    		 showDirs,
    		 localShowSamples,
    		 showIgnoredData,
    		 showRad);
}

public void listData(String title,
         String path,
         boolean showMotors,
         boolean []showMeasCalc,
         boolean [] showColors,
         boolean [] showDirs,
         boolean [] showSamples,
         boolean showIgnoredData,
         boolean showRad
         ){
    String header="";
    String [] sMeasCalc={"M","C","Z"};
    String [] sDirs={"S","T"};
    String [] sColors={"R","G","B"};
    int numColors=3;
    int numDirs=2;
    boolean [] selChannels=fieldFitting.getSelectedChannels();
    int [] selChanIndices= new int[selChannels.length];
    selChanIndices[0]=0;
    for (int i=1;i<selChanIndices.length;i++){
        selChanIndices[i]= selChanIndices[i-1]+(selChannels[i-1]?1:0);
    }
    double [][][] cosSin2Tab=new double[sampleCoord.length][sampleCoord[0].length][3];
        for (int i=0;i<sampleCoord.length;i++) for (int j=0;j<sampleCoord[i].length;j++){
            double delta_x=sampleCoord[i][j][0]-currentPX0;
            double delta_y=sampleCoord[i][j][1]-currentPY0;
            double r2=delta_x*delta_x+delta_y*delta_y;
            if (r2>0.0) {
                cosSin2Tab[i][j][0]= delta_x*delta_x/r2; // cos^2
                cosSin2Tab[i][j][1]= delta_y*delta_y/r2; // sin^2
                cosSin2Tab[i][j][2]=2*delta_x*delta_y/r2; // 2*cos*sin
            } else {
                cosSin2Tab[i][j][0]=1.0;
                cosSin2Tab[i][j][1]=0.0;
                cosSin2Tab[i][j][2]=0.0;
            }
        }

        StringBuffer sb = new StringBuffer();
    if (showMotors) header +="M1\tM2\tM3\t";
    boolean first=true;
    for (int i=0;i<sampleCoord.length;i++) for (int j=0;j<sampleCoord[i].length;j++) if (showSamples[j+i*sampleCoord[i].length]){
    	if (showMeasCalc[2]) {
    		if (!first) header+="\t";
    		first=false;
    		header+="Y"+i+"X"+j+sMeasCalc[2];
    	}
    	for (int c=0;c<numColors;c++) if (showColors[c]){
    		for (int d=0;d<numDirs;d++) if (showDirs[d]){
    			for (int m=0;m<2;m++) if (showMeasCalc[m]){
    				if (!first) header+="\t";
    				first=false;
    				header+="Y"+i+"X"+j+sColors[c]+sDirs[d]+sMeasCalc[m];
    			}                 
    		}
    	}
    }
    for (FocusingFieldMeasurement ffm:measurements){
//        double [][][][] samples=ffm.samples;
        double [][][][] samplesFull=new double [sampleCoord.length][sampleCoord[0].length][numColors][numDirs];
        double [][][][] calcSamples=new double [sampleCoord.length][sampleCoord[0].length][numColors][numDirs];
        double [][] calcZ= new double [sampleCoord.length][sampleCoord[0].length];

        for (int i=0;i<sampleCoord.length;i++)
        	for (int j=0;j<sampleCoord[0].length;j++) {
        		calcZ[i][j]=Double.NaN;
        		for (int c=0;c<numColors;c++)
        			for (int d=0;d<numDirs;d++) {
        				samplesFull[i][j][c][d]=Double.NaN;
        				calcSamples[i][j][c][d]=Double.NaN;
        			}
        	}
//showIgnoredData MeasuredSample [] dataVector
        boolean [] enable=dataWeightsToBoolean();
        for (int index=0;index<dataVector.length;index++) if (ffm.timestamp.equals(dataVector[index].timestamp) && (showIgnoredData || (index>=enable.length) || enable[index])){
    		int chn=dataVector[index].channel;
    		int sample=dataVector[index].sampleIndex;
    		int sampleRow=sample / sampleCoord[0].length;
    		int sampleCol=sample % sampleCoord[0].length;
    		int color = chn / 2;
    		int dir = chn % 2;
    		samplesFull[sampleRow][sampleCol][color][dir]=dataVector[index].value;
        }
// still needed if     showIgnoredData!    
/*        
        if (showMeasCalc[0] && (samples!=null)) for (int i=0;i<sampleCoord.length;i++){
        	if ((i<samples.length) && (samples[i]!=null)) for (int j=0;j<sampleCoord[i].length;j++){
        		if ((j<samples[i].length) && (samples[i][j]!=null)) for (int c=0;c<numColors;c++){
        			if ((c<samples[i][j].length) && (samples[i][j][c]!=null)) {
        				double sagTan[] = {
        						Math.sqrt(
        								cosSin2Tab[i][j][0]*samples[i][j][c][0]+
        								cosSin2Tab[i][j][1]*samples[i][j][c][1]+
        								cosSin2Tab[i][j][2]*samples[i][j][c][2]),
        								Math.sqrt(
        										cosSin2Tab[i][j][1]*samples[i][j][c][0]+
        										cosSin2Tab[i][j][0]*samples[i][j][c][1]-
        										cosSin2Tab[i][j][2]*samples[i][j][c][2])};
        				if (debugLevel>3) System.out.print("\n"+ffm.motors[2]+" i= "+i+" j= "+j+" c= "+c+" sagTan= "+sagTan[0]+" "+sagTan[1]+" ");
        				for (int d=0;d<numDirs;d++){
        					if (debugLevel>3) System.out.print(" d= "+d+" ");
        					int chn=d+numDirs*c;
        					//                     System.out.println("i="+i+", j="+j+", c="+c+", d="+d);
        					// saved values are PSF radius, convert to FWHM by multiplying by 2.0
        					double value=sagTan[d]*2.0;
        					// discard above threshold (in pixels, raw FWHM data)
        					if (Double.isNaN(value)) continue; // bad measurement
        					if (debugLevel>3) System.out.print(" A "+value);
        					if (!showIgnoredData && (thresholdMax != null)){
        						if (value >= thresholdMax[chn]){
        							if (debugLevel>3) System.out.print(" > "+thresholdMax[chn]);
        							continue; // bad measurement (above threshold)
        						}
        					}
        					if (debugLevel>3) System.out.print(" B "+value);
        					// correct for for minimal measurement;
        					if (useMinMeas && (minMeas != null)){
        						if (value<minMeas[chn]) {
        							if (debugLevel>3) System.out.print(" < "+minMeas[chn]);
        							continue; // bad measurement (smaller than correction)
        						}
        						value=Math.sqrt(value*value-minMeas[chn]*minMeas[chn]);
        					}
        					if (debugLevel>3) System.out.print(" C "+value);
        					// correct for for maximal measurement;
        					if (useMaxMeas && (maxMeas != null)) {
        						if (value >= maxMeas[chn]) {
        							if (debugLevel>3) System.out.print(" > "+maxMeas[chn]);
        							continue; // bad measurement (larger than correction)
        						}
        						value = 1.0/Math.sqrt(1.0/(value*value)-1.0/(maxMeas[chn]*maxMeas[chn]));
        					}
        					if (debugLevel>3) System.out.print(" D "+value);
        					// convert to microns from pixels
        					value*=getPixelUM();
        					samplesFull[i][j][c][d]=value;
        					if (debugLevel>3) System.out.print(" E "+value);
        				}
        			}
        		}
        	}
        }
*/

        // Now calculate values for the same samples
        if (showMeasCalc[1]){
        	for (int i=0;i<sampleCoord.length;i++)
        		for (int j=0;j<sampleCoord[0].length;j++) {
        			
        			if ((i==0) &&  (j==3) && (ffm.motors[0]==2209)){
        				System.out.println("listData(), i="+i+", j="+j);
        			}
        			
        			double [] subData=fieldFitting.getValsDerivatives(
        					flattenIndex(i,j),
        					sagittalMaster, // dependent channel does not have center parameters, but that is only used for derivs.
        					ffm.motors, // 3 motor coordinates
        					sampleCoord[i][j][0], // pixel x
        					sampleCoord[i][j][1], // pixel y
        					null);
        			for (int c=0;c<numColors;c++)
        				for (int d=0;d<numDirs;d++) {
        					int index=d+c*numDirs;
        					if (selChannels[index])
        						calcSamples[i][j][c][d]=subData[selChanIndices[index]];
        				}
        		}
        }
        // calculate Z from the mechanical
        if (showMeasCalc[2]){
        	for (int i=0;i<sampleCoord.length;i++)
        		for (int j=0;j<sampleCoord[0].length;j++) {
        			calcZ[i][j]=     fieldFitting.getMotorsZ(
        					ffm.motors, // 3 motor coordinates
        					sampleCoord[i][j][0], // pixel x
        					sampleCoord[i][j][1]); // pixel y
        		}
        }
        // combine line
        if (showMotors) sb.append(""+ffm.motors[0]+"\t"+ffm.motors[1]+"\t"+ffm.motors[2]+"\t");
        first=true;
        for (int i=0;i<sampleCoord.length;i++) for (int j=0;j<sampleCoord[i].length;j++) if (showSamples[j+i*sampleCoord[i].length]){
        	if (showMeasCalc[2]) {
        		if (!first) sb.append("\t");
        		first=false;
        		sb.append(calcZ[i][j]);
        	}

        	for (int c=0;c<showColors.length;c++) if (showColors[c]){
        		for (int d=0;d<showDirs.length;d++) if (showDirs[d]){
        			for (int m=0;m<2;m++) if (showMeasCalc[m]){
        				if (!first) sb.append("\t");
        				first=false;
        				if (m==0) sb.append(samplesFull[i][j][c][d]);
        				else sb.append(calcSamples[i][j][c][d]);
        			}                 
        		}
        	}
        }
        sb.append("\n");
    }

    if (debugLevel>3) System.out.println();

    if (showRad){
    	sb.append(header+"\n");
    	if (showMotors) sb.append("Sample\tradius\t(mm)\t");
    	first=true;
    	for (int i=0;i<sampleCoord.length;i++) for (int j=0;j<sampleCoord[i].length;j++) if (showSamples[j+i*sampleCoord[i].length]){
    		double rad=fieldFitting.getRadiusMM(
    				sampleCoord[i][j][0], // pixel x
    				sampleCoord[i][j][1]); // pixel y
    		if (showMeasCalc[2]) {
    			if (!first) sb.append("\t");
    			first=false;
    			sb.append(rad);
    		}
    		for (int c=0;c<numColors;c++) if (showColors[c]){
    			for (int d=0;d<numDirs;d++) if (showDirs[d]){
    				for (int m=0;m<2;m++) if (showMeasCalc[m]){
    					if (!first) sb.append("\t");
    					first=false;
    					sb.append(rad);
    				}                 
    			}
    		}
    	}
    	sb.append("\n");
    }
    if (path!=null) {
    	CalibrationFileManagement.saveStringToFile (
    			path,
    			header+"\n"+sb.toString());
    } else {
    	new TextWindow(title, header, sb.toString(), 800,1000);
    }

}
public void showCurvCorr(){
    fieldFitting.showCurvCorr("curv_corr");
}

public void listCombinedResults(){
// private boolean rslt_show_z_axial=true;
// private boolean rslt_show_z_individual=true;
// private boolean rslt_show_f_axial=true;
// private boolean rslt_show_f_individual=true;
// private double rslt_scan_below=-10.0;
// private double rslt_scan_above= 10.0;
// private double rslt_scan_step= 5.0;
// private boolean rslt_mtf50_mode= true;
// public double fwhm_to_mtf50=500.0; // put actual number
    
    double [] center_z=fieldFitting.getZCenters(false); // do not solve, use z0 coefficient
    double [] centerFWHM={
    		fieldFitting.getCalcValuesForZ(center_z[0],0.0,null)[1],
    		fieldFitting.getCalcValuesForZ(center_z[1],0.0,null)[3],
    		fieldFitting.getCalcValuesForZ(center_z[2],0.0,null)[5]
    };
    double [] best_qb_axial= fieldFitting.getBestQualB(
            k_red,
            k_blue,
            false);
    double [] best_qb_corr= fieldFitting.getBestQualB(
            k_red,
            k_blue,
            true);
    GenericDialog gd = new GenericDialog("Setup results table FWHM="+IJ.d2s(best_qb_corr[1],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/best_qb_corr[1],2)+" lp/mm");
    gd.addMessage("Best center focus for Red "+ IJ.d2s(center_z[0],3)+" um"+
    		", FWHM="+IJ.d2s(centerFWHM[0],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/centerFWHM[0],2)+" lp/mm");
    gd.addMessage("Best center focus for Green "+ IJ.d2s(center_z[1],3)+" um"+
    		", FWHM="+IJ.d2s(centerFWHM[1],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/centerFWHM[1],2)+" lp/mm");
    gd.addMessage("Best center focus for Blue "+ IJ.d2s(center_z[2],3)+" um"+
    		", FWHM="+IJ.d2s(centerFWHM[2],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/centerFWHM[2],2)+" lp/mm");
    gd.addMessage("Best composite distance for FWHM^4, axial model "+ IJ.d2s(best_qb_axial[0],3)+" um"+
    		", FWHM="+IJ.d2s(best_qb_axial[1],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/best_qb_axial[1],2)+" lp/mm");
    gd.addMessage("Best composite distance for FWHM^4, individual "+ IJ.d2s(best_qb_corr[0],3)+"  um"+
    		", FWHM="+IJ.d2s(best_qb_corr[1],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/best_qb_corr[1],2)+" lp/mm");
    for (int i=0;i<rslt_show_chn.length;i++){
    	gd.addCheckbox("Show results for "+fieldFitting.getDescription(i), this.rslt_show_chn[i]);
    }
    gd.addCheckbox("Show best focus distance (axial model)", this.rslt_show_z_axial);
    gd.addCheckbox("Show ring-averaged focus distance", this.rslt_show_z_smooth);
    gd.addCheckbox("Show best focus distance (per-sample adjusted)", this.rslt_show_z_individual);
    gd.addCheckbox("Show constant-z sections (axial model)", this.rslt_show_f_axial);
    gd.addCheckbox("Show ring-averaged per-sample adjusted section data", this.rslt_show_f_smooth);
    gd.addCheckbox("Show constant-z sections (per-sample adjusted)", this.rslt_show_f_individual);
    gd.addNumericField("Ring averaging radial sigma", this.rslt_show_smooth_sigma, 3,5,"mm");
    gd.addCheckbox("Show mtf50 (false - PSF FWHM)", this.rslt_mtf50_mode);
    gd.addCheckbox("Find z for minimum, unchecked - use parameter", this.rslt_solve);
        
    gd.addCheckbox("Show focal distance relative to best composite focus (false - to center green )", this.z_relative);
    
    gd.addMessage("Multiple section setup:");
        gd.addNumericField("Scan from (relative to green center)", this.rslt_scan_below, 3,7,"um");
        gd.addNumericField("Scan to (relative to green center)", this.rslt_scan_above, 3,7,"um");
        gd.addNumericField("Scan step", this.rslt_scan_step, 3,7,"um");
    gd.showDialog();
    if (gd.wasCanceled()) {
         return;
     }
    
    for (int i=0;i<rslt_show_chn.length;i++){
    	this.rslt_show_chn[i]=gd.getNextBoolean();
    }
    this.rslt_show_z_axial= gd.getNextBoolean();
    this.rslt_show_z_smooth= gd.getNextBoolean();
    this.rslt_show_z_individual= gd.getNextBoolean();
    this.rslt_show_f_axial= gd.getNextBoolean();
    this.rslt_show_f_smooth= gd.getNextBoolean();
    this.rslt_show_f_individual= gd.getNextBoolean();
    this.rslt_show_smooth_sigma= gd.getNextNumber();
    this.rslt_mtf50_mode= gd.getNextBoolean();
    this.rslt_solve=gd.getNextBoolean();
    this.z_relative= gd.getNextBoolean();
    this.rslt_scan_below= gd.getNextNumber();
    this.rslt_scan_above= gd.getNextNumber();
    this.rslt_scan_step= gd.getNextNumber();

    listCombinedResults(
    		"Field curvature measuremnts results", //String title,
    		null, //String path,
            z_relative?best_qb_corr[0]:center_z[1],
    		rslt_show_chn, //boolean [] show_chn,
    		rslt_show_z_axial, //boolean show_z_axial,
    		rslt_show_z_smooth, // boolean rslt_show_z_smooth;
    		rslt_show_z_individual, //boolean show_z_individual,
    		rslt_show_f_axial, //boolean show_f_axial,
    		rslt_show_f_smooth, //boolean rslt_show_f_smooth;
    		rslt_show_f_individual, //boolean show_f_individual,
    		rslt_show_smooth_sigma, // double rslt_show_smooth_sigma;
    		(z_relative?best_qb_corr[0]:center_z[1])+rslt_scan_below, // double scan_below,
    		(z_relative?best_qb_corr[0]:center_z[1])+rslt_scan_above, //double scan_above,
    		rslt_scan_step, //double scan_step,
    		rslt_mtf50_mode, //boolean freq_mode)
    		rslt_solve); // currently if not using z^(>=2) no numeric solution is required - z0 is the minimum 
}

public double [][] filterListSamples(
		double sigma,
		double [] radiuses, // numsaples+1 long!
		double [] saggital, // numsamples!
		double []tangential,
		double [] saggitalWeights, // numsamples long
		double [] tangentialWeights){// numsamples long
	double [][] data={saggital,tangential};
	double [][] weights={saggitalWeights,tangentialWeights};
	
	int len=saggital.length+1;
	double [][] result = new double[2][len];
	double kexp=-0.5/(sigma*sigma);
	for (int dir=0;dir<2;dir++) {
		if (data[dir]!=null) {
			result[dir] = new double[len];
			for (int i=0;i<len;i++){
				double sum_w=0.0;
				double sum_v=0.0;
				for (int otherDir=0;otherDir<2;otherDir++) {
					if (data[otherDir]!=null)
						for (int j=1;j<len;j++){
							double r= (otherDir==dir)?(radiuses[i]-radiuses[j]):(radiuses[i]+radiuses[j]);
							double w=Math.exp(kexp*r*r);
							if (weights[otherDir]!=null){
								w*=weights[otherDir][j-1];
							}
							if (!Double.isNaN(data[otherDir][j-1])) {
								sum_v+=data[otherDir][j-1]*w;
								sum_w+=w;
							}
						}
				}
				result[dir][i]=(sum_w>0.0)?(sum_v/sum_w):0.0;
			}
		} else {
			result[dir]=null;
		}
	}
	return result;
}


public void listCombinedResults(
		String title,
		String path,
		double z0, // subtract from z
		boolean [] show_chn,
		boolean show_z_axial,
		boolean show_z_smooth,
		boolean show_z_individual,
		boolean show_f_axial,
		boolean show_f_smooth,
		boolean show_f_individual,
		double smooth_sigma,
		double scan_below,
		double scan_above,
		double scan_step,
		boolean freq_mode,
		boolean solveZ){ // currently if not using z^(>=2) no numeric solution is required - z0 is the minimum 

	String [] chnNames={"RS","RT","GS","GT","BS","BT"};
	// Calculate weights of each channel/sample
	double [][] sampleWeights=new double[getNumChannels()][getNumSamples()];
	for (int chn=0;chn<sampleWeights.length;chn++)  for (int sample=0;sample<sampleWeights[chn].length;sample++) {
		sampleWeights[chn][sample]=0.0;
	}
	for (int index=0;index<dataVector.length;index++) {
		sampleWeights[dataVector[index].channel][dataVector[index].sampleIndex] += dataWeights[index];
	}

	double k=this.fwhm_to_mtf50; //TODO: correct psf fwhm to mtf50 conversion
	//     double k=2*Math.log(2.0)/Math.PI*1000;
	//     String header="Z(um)\tComposite\tRed\tGreen\tBlue";
	String header="radius(mm)";
	if (show_z_axial){
		for (int i=0;i<show_chn.length;i++){
			if (show_chn[i])header+="\tZ"+chnNames[i];
		}
	}
	if (show_z_smooth){
		for (int i=0;i<show_chn.length;i++){
			if (show_chn[i])header+="\tZ"+chnNames[i]+"~";
		}
	}
	
	if (show_z_individual){
		for (int i=0;i<show_chn.length;i++){
			if (show_chn[i])header+="\tZ"+chnNames[i]+"*";
		}
	}
	int numSect=0;
	if (show_f_axial || show_f_individual) for (double z=scan_below;z<=scan_above;z+=scan_step){
		if (show_f_axial) {
			for (int i=0;i<show_chn.length;i++){
				if (show_chn[i])header+="\t"+chnNames[i]+IJ.d2s(z-z0,1);
			}
		}
		if (show_f_smooth) {
			for (int i=0;i<show_chn.length;i++){
				if (show_chn[i])header+="\t"+chnNames[i]+IJ.d2s(z-z0,1)+"~";
			}
		}
		if (show_f_individual) {
			for (int i=0;i<show_chn.length;i++){
				if (show_chn[i])header+="\t"+chnNames[i]+IJ.d2s(z-z0,1)+"*";
			}
		}
		numSect++;
	}
	double [] radiuses0=fieldFitting.getSampleRadiuses();
	int numSamples=radiuses0.length;
	double [] radiuses=new double [numSamples+1];
	radiuses[0]=0.0;
	for (int i=0;i<numSamples;i++) radiuses[i+1]=radiuses0[i];
	double [][][][] f_values=new double [numSect][3][][];
	int sect=0;
	for (double z=scan_below;z<=scan_above;z+=scan_step){
		if (show_f_axial){
			double [][] f=fieldFitting.getCalcValuesForZ(z, false,true);
			double [] f0=fieldFitting.getCalcValuesForZ(z, 0.0,null);
			f_values[sect][0]=new double [f.length][];
			for (int chn=0;chn<f.length;chn++){
				if (f[chn]!=null){
					f_values[sect][0][chn]=new double[f[chn].length+1];
					f_values[sect][0][chn][0]=f0[chn];
					for (int i=0;i<f[chn].length;i++){
						f_values[sect][0][chn][i+1]=f[chn][i];
					}
				} else f_values[sect][0][chn]=null;
			}
		} else {
			f_values[sect][0]=null;
		}
		if (show_f_individual){
			double [][] f=fieldFitting.getCalcValuesForZ(z, true,true);
			f_values[sect][1]=new double [f.length][];
			for (int chn=0;chn<f.length;chn++){
				if (f[chn]!=null){
					f_values[sect][1][chn]=new double[f[chn].length+1];
					f_values[sect][1][chn][0]=Double.NaN;
					for (int i=0;i<f[chn].length;i++){
						f_values[sect][1][chn][i+1]=f[chn][i];
					}
				} else f_values[sect][1][chn]=null;
			}
		} else {
			f_values[sect][1]=null;
		}
		if (show_f_smooth){
			double [][] f=fieldFitting.getCalcValuesForZ(z, true,true); // corrected (individual), 
			f_values[sect][2]=new double [f.length][];
			for (int chn=0;chn<f.length;chn+=2){
				double [][] smooth=filterListSamples(
						smooth_sigma,
						radiuses,  // all 3 arrays should be numsaples+1 long
						f[chn] ,   //double [] saggital, // not ordered, but [0]==0.0
						f[chn+1] , //double []tangential,
						sampleWeights[chn], //double [] saggitalWeights, // numsamples long
						sampleWeights[chn+1]); //double [] tangentialWeights){// numsamples long
				f_values[sect][2][chn]=smooth[0];
				f_values[sect][2][chn+1]=smooth[1];
			}
		} else {
			f_values[sect][2]=null;
		}
		sect++;
	}
	double [][][] z_values=new double [3][][];
	if (show_z_axial){
		double [][] zai=fieldFitting.getCalcZ(false,true,solveZ);
		double [] zai0=fieldFitting.getCalcZ(0.0,solveZ);
		z_values[0]=new double [zai.length][];
		for (int chn=0;chn<zai.length;chn++){
			if (zai[chn]!=null){
				z_values[0][chn]=new double[zai[chn].length+1];
				z_values[0][chn][0]=zai0[chn];
				for (int i=0;i<zai[chn].length;i++){
					z_values[0][chn][i+1]=zai[chn][i];
				}
				
			} else z_values[0][chn]=null;
		}
	} else z_values[0] = null;
	if (show_z_individual){
		double [][] zai=fieldFitting.getCalcZ(true,true,solveZ);
		z_values[1]=new double [zai.length][];
		for (int chn=0;chn<zai.length;chn++){
			if (zai[chn]!=null){
				z_values[1][chn]=new double[zai[chn].length+1];
				z_values[1][chn][0]=Double.NaN;
				for (int i=0;i<zai[chn].length;i++){
					z_values[1][chn][i+1]=zai[chn][i];
				}
				
			} else z_values[1][chn]=null;
		}
	} else z_values[1] = null;
	if (show_f_smooth){
		double [][] zai=fieldFitting.getCalcZ(true,true,solveZ);
		z_values[2]=new double [zai.length][];
		for (int chn=0;chn<zai.length;chn+=2){
			double [][] smooth=filterListSamples(
					smooth_sigma,
					radiuses,  // all 3 arrays should be numsaples+1 long
					zai[chn] ,   //double [] saggital, // not ordered, but [0]==0.0
					zai[chn+1] , //double []tangential,
					sampleWeights[chn], //double [] saggitalWeights, // numsamples long
					sampleWeights[chn+1]); //double [] tangentialWeights){// numsamples long
			z_values[2][chn]=smooth[0];
			z_values[2][chn+1]=smooth[1];
		}
	} else {
		z_values[2]=null;
	}
	
	
	
	
	
	StringBuffer sb = new StringBuffer();
	for (int n=0;n<radiuses.length;n++){
		sb.append(IJ.d2s(radiuses[n],3));
		if (show_z_axial){
			for (int i=0;i<show_chn.length;i++){
				if (show_chn[i]) sb.append("\t"+IJ.d2s(z_values[0][i][n]-z0,3));
			}
		}
		if (show_z_smooth){
			for (int i=0;i<show_chn.length;i++){
				if (show_chn[i]) sb.append("\t"+IJ.d2s(z_values[2][i][n]-z0,3));
			}
		}
		if (show_z_individual){
			for (int i=0;i<show_chn.length;i++){
				if (show_chn[i]) sb.append("\t"+IJ.d2s(z_values[1][i][n]-z0,3));
			}
		}
		if (show_f_axial || show_f_individual) for (int s=0;s<numSect;s++){
			if (show_f_axial) for (int i=0;i<show_chn.length;i++) if (show_chn[i]) {
				if (freq_mode){
					sb.append("\t"+IJ.d2s(k/f_values[s][0][i][n],3));
				} else {
					sb.append("\t"+IJ.d2s(f_values[s][0][i][n],3));
				}
			}
			if (show_f_smooth) for (int i=0;i<show_chn.length;i++) if (show_chn[i]) {
				if (freq_mode){
					sb.append("\t"+IJ.d2s(k/f_values[s][2][i][n],3));
				} else {
					sb.append("\t"+IJ.d2s(f_values[s][2][i][n],3));
				}
			}
			if (show_f_individual) for (int i=0;i<show_chn.length;i++) if (show_chn[i]) {
				if (freq_mode){
					sb.append("\t"+IJ.d2s(k/f_values[s][1][i][n],3));
				} else {
					sb.append("\t"+IJ.d2s(f_values[s][1][i][n],3));
				}
			}
		}
		sb.append("\n");
	}
	if (path!=null) {
		CalibrationFileManagement.saveStringToFile (
				path,
				header+"\n"+sb.toString());
	} else {
		new TextWindow(title, header, sb.toString(), 800,1000);
	}
}


public void listScanQB(){
    
    double [] center_z=fieldFitting.getZCenters(false); // do not solve, use z0 coefficient
    double [] centerFWHM={
    		fieldFitting.getCalcValuesForZ(center_z[0],0.0,null)[1],
    		fieldFitting.getCalcValuesForZ(center_z[1],0.0,null)[3],
    		fieldFitting.getCalcValuesForZ(center_z[2],0.0,null)[5]
    };
    double [] best_qb_axial= fieldFitting.getBestQualB(
                k_red,
                k_blue,
                false);
    double [] best_qb_corr= fieldFitting.getBestQualB( //best_qb_corr[0] - distance (motorZ)
                k_red,
                k_blue,
                true);
     GenericDialog gd = new GenericDialog("Setup quality-B table FWHM="+IJ.d2s(best_qb_corr[1],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/best_qb_corr[1],2)+" lp/mm");
     gd.addMessage("Best center focus for Red "+ IJ.d2s(center_z[0],3)+" um"+
     		", FWHM="+IJ.d2s(centerFWHM[0],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/centerFWHM[0],2)+" lp/mm");
     gd.addMessage("Best center focus for Green "+ IJ.d2s(center_z[1],3)+" um"+
     		", FWHM="+IJ.d2s(centerFWHM[1],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/centerFWHM[1],2)+" lp/mm");
     gd.addMessage("Best center focus for Blue "+ IJ.d2s(center_z[2],3)+" um"+
     		", FWHM="+IJ.d2s(centerFWHM[2],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/centerFWHM[2],2)+" lp/mm");
     gd.addMessage("Best composite distance for FWHM^4, axial model "+ IJ.d2s(best_qb_axial[0],3)+" um"+
     		", FWHM="+IJ.d2s(best_qb_axial[1],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/best_qb_axial[1],2)+" lp/mm");
     gd.addMessage("Best composite distance for FWHM^4, individual "+ IJ.d2s(best_qb_corr[0],3)+"  um"+
     		", FWHM="+IJ.d2s(best_qb_corr[1],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/best_qb_corr[1],2)+" lp/mm");
        gd.addNumericField("Scan from (relative to green center)", this.qb_scan_below, 3,7,"um");
        gd.addNumericField("Scan to (relative to green center)", this.qb_scan_above, 3,7,"um");
        gd.addNumericField("Scan step", this.qb_scan_step, 3,7,"um");
        gd.addNumericField("Relative (to green) weight of red channel",100* this.k_red, 3,7,"%");
        gd.addNumericField("Relative (to green) weight of blue channel",100* this.k_blue, 3,7,"%");
        gd.addCheckbox("Use per-sample location corrected data", this.qb_use_corrected);
        gd.addCheckbox("Show mtf50 (false - PSF FWHM)", this.qb_invert);
        gd.addCheckbox("Show focal distance relative to best composite focus (false - to center green )", this.z_relative);
        
     gd.showDialog();
     if (gd.wasCanceled()) {
         return;
     }
        this.qb_scan_below= gd.getNextNumber();
        this.qb_scan_above= gd.getNextNumber();
        this.qb_scan_step= gd.getNextNumber();
        this.k_red= 0.01*gd.getNextNumber();
        this.k_blue= 0.01*gd.getNextNumber();
        this.qb_use_corrected= gd.getNextBoolean();
        this.qb_invert= gd.getNextBoolean();
        this.z_relative= gd.getNextBoolean();
        
        
        listScanQB(
             "Focus quality vs. distance", // String title,
             null, //String path,
             z_relative?best_qb_corr[0]:center_z[1],
             center_z[1]+this.qb_scan_below, //double min_z, // absolute
             center_z[1]+this.qb_scan_above, //double max_z,
             this.qb_scan_step , //double step_z,
             this.k_red, //double kr,
             this.k_blue, //double kb,
             this.qb_use_corrected, //boolean corrected,
             this.qb_invert); //boolean freq_mode)
}

public void listScanQB(
        String title,
        String path,
        double z0,
        double min_z, // absolute
        double max_z,
        double step_z,
            double kr,
            double kb,
            boolean corrected,
            boolean freq_mode){
    double k=this.fwhm_to_mtf50; //TODO: correct psf fwhm to mtf50 conversion
    String header="Z(um)\tComposite\tRed\tGreen\tBlue";
        StringBuffer sb = new StringBuffer();
        for (double z=min_z; z<=max_z;z+=step_z){
            double qb_w= fieldFitting.getQualB(
                    z,
                    kr,
                    kb,
                    corrected);
            double [] qb_rgb=fieldFitting.getQualB(
                    z,
                    corrected);
            if (freq_mode){
                sb.append(IJ.d2s(z-z0,3)+"\t"+IJ.d2s(k/qb_w,3)+"\t"+IJ.d2s(k/qb_rgb[0],3)+"\t"+IJ.d2s(k/qb_rgb[1],3)+"\t"+IJ.d2s(k/qb_rgb[2],3) +"\n");
            } else {
                sb.append(IJ.d2s(z-z0,3)+"\t"+IJ.d2s(qb_w,3)+"\t"+IJ.d2s(qb_rgb[0],3)+"\t"+IJ.d2s(qb_rgb[1],3)+"\t"+IJ.d2s(qb_rgb[2],3) +"\n");
            }
        }
        
//        for (String s:this.fieldFitting.getParameterValueStrings(true,true)){ //this.showDisabledParams)){ //
//         sb.append(s+"\n");
//        }
        
// Add result to the bottom of the file
        double [] center_z=fieldFitting.getZCenters(false); // z0 coefficient, do not find minimum
        double [] centerFWHM={
        		fieldFitting.getCalcValuesForZ(center_z[0],0.0,null)[1],
        		fieldFitting.getCalcValuesForZ(center_z[1],0.0,null)[3],
        		fieldFitting.getCalcValuesForZ(center_z[2],0.0,null)[5]
        };
//        double [] best_qb_axial= fieldFitting.getBestQualB(
//                k_red,
//                k_blue,
//                false);
        double [] best_qb_corr= fieldFitting.getBestQualB(
                k_red,
                k_blue,
                true);
        sb.append("Results:\t---\t---\t---\t---\n");
        sb.append("Red center"+"\t"+  IJ.d2s(center_z[0],3)+"\t"+IJ.d2s(centerFWHM[0],3)+"\t"+IJ.d2s(fwhm_to_mtf50/centerFWHM[0],2)+"\t"+"lp/mm" +"\n");
        sb.append("Green center"+"\t"+IJ.d2s(center_z[1],3)+"\t"+IJ.d2s(centerFWHM[1],3)+"\t"+IJ.d2s(fwhm_to_mtf50/centerFWHM[1],2)+"\t"+"lp/mm" +"\n");
        sb.append("Blue center"+"\t"+ IJ.d2s(center_z[2],3)+"\t"+IJ.d2s(centerFWHM[2],3)+"\t"+IJ.d2s(fwhm_to_mtf50/centerFWHM[2],2)+"\t"+"lp/mm" +"\n");
        sb.append("Composite ^4"+"\t"+ IJ.d2s(best_qb_corr[0],3)+"\t"+IJ.d2s(best_qb_corr[1],3)+"\t"+IJ.d2s(fwhm_to_mtf50/best_qb_corr[1],2)+"\t"+"lp/mm" +"\n");
        
        sb.append("Lens center "+"\t"+  "used" +             "\t"+IJ.d2s(currentPX0,1)+"\t"+ IJ.d2s(currentPY0,1)+"\t"+"pix" +"\n");
        sb.append("Lens center "+"\t"+  "distortions" +      "\t"+IJ.d2s(pX0_distortions,1)+"\t"+IJ.d2s(pY0_distortions,1)+"\t"+"pix" +"\n");
        
        if (path!=null) {
            CalibrationFileManagement.saveStringToFile (
                    path,
                    header+"\n"+sb.toString());
        } else {
            new TextWindow(title, header, sb.toString(), 800,1000);
        }
}


public boolean dialogLMAStep(boolean [] state){
    String [] states={
            "Worse, increase lambda",
            "Better, decrease lambda",
            "Failed to fit",
            "Fitting Successful"};
    int iState=(state[0]?1:0)+(state[1]?2:0);

     GenericDialog gd = new GenericDialog("Levenberg-Marquardt algorithm step");
//     String [][] parameterDescriptions=fittingStrategy.distortionCalibrationData.parameterDescriptions;
    gd.addMessage("Current state="+states[iState]);
    gd.addMessage("Current series="+this.currentStrategyStep);
    gd.addMessage("Iteration step="+this.iterationStepNumber);
    
    gd.addMessage("Initial RMS="+IJ.d2s(this.firstRMS,6)+", Current RMS="+IJ.d2s(this.currentRMS,6)+", new RMS="+IJ.d2s(this.nextRMS,6));
    gd.addMessage("Pure initial RMS="+IJ.d2s(this.firstRMSPure,6)+", Current RMS="+IJ.d2s(this.currentRMSPure,6)+", new RMS="+IJ.d2s(this.nextRMSPure,6));
    
    if (this.showParams) {
        gd.addMessage("==== Current parameter values ===");
        for (String s:this.fieldFitting.getParameterValueStrings(this.showDisabledParams, this.showCorrectionParams)){ //
        gd.addMessage(s);
        }
        gd.addMessage("");

    }
    
        gd.addNumericField("Lambda ", this.lambda, 5);
        gd.addNumericField("Multiply lambda on success", this.lambdaStepDown, 5);
        gd.addNumericField("Threshold RMS to exit LMA", this.thresholdFinish, 7,9,"pix");
        gd.addNumericField("Multiply lambda on failure", this.lambdaStepUp, 5);
        gd.addNumericField("Threshold lambda to fail", this.maxLambda, 5);
        gd.addNumericField("Maximal number of iterations", this.numIterations, 0);

        gd.addCheckbox("Dialog after each iteration step", this.stopEachStep);
        gd.addCheckbox("Dialog after each iteration series", this.stopEachSeries);
        gd.addCheckbox("Dialog after each failure", this.stopOnFailure);
        gd.addCheckbox("Show modified parameters", this.showParams);
        gd.addCheckbox("Show disabled parameters", this.showDisabledParams);
        gd.addCheckbox("Show per-sample correction parameters", this.showCorrectionParams);

//        gd.addCheckbox("Show debug images before correction",this.showThisImages);
//        gd.addCheckbox("Show debug images after correction", this.showNextImages);
        gd.addMessage("Done will save the current (not new!) state and exit, Continue will proceed according to LMA");
        gd.enableYesNoCancel("Continue", "Done");
        WindowTools.addScrollBars(gd);

     gd.showDialog();
     if (gd.wasCanceled()) {
         this.saveSeries=false;
         return false;
     }
        this.lambda= gd.getNextNumber();
        this.lambdaStepDown= gd.getNextNumber();
        this.thresholdFinish= gd.getNextNumber();
        this.lambdaStepUp= gd.getNextNumber();
        this.maxLambda= gd.getNextNumber();
        this.numIterations= (int) gd.getNextNumber();
        this.stopEachStep= gd.getNextBoolean();
        this.stopEachSeries= gd.getNextBoolean();
        this.stopOnFailure= gd.getNextBoolean();
        this.showParams= gd.getNextBoolean();
        this.showDisabledParams= gd.getNextBoolean();
        this.showCorrectionParams= gd.getNextBoolean();

//        this.showThisImages= gd.getNextBoolean();
//        this.showNextImages= gd.getNextBoolean();
     this.saveSeries=true;
     return gd.wasOKed();
}
public double getAdjustRMS(
		FocusingFieldMeasurement measurement,
		boolean filterZ,
		double filterByValueScale,
		double z,
		double tx,
		double ty){
	fieldFitting.selectZTilt(false);
	fieldFitting.mechanicalFocusingModel.setZTxTy(z,tx,ty);
	double [] sv= fieldFitting.createParameterVector(sagittalMaster);
	setDataVector(
			false, // calibrate mode
			createDataVector(measurement));
	if (filterZ) {
		boolean [] en=dataWeightsToBoolean();
		en= filterByZRanges(
				zRanges,
				en);
		maskDataWeights(en);
		prevEnable=en;
    	int numEn=getNumEnabledSamples(en);
    	if (numEn<minLeftSamples) return Double.NaN;
	}
	if (filterByValueScale>0.0){
		boolean [] en=dataWeightsToBoolean();
		en= filterByValue(
				filterByValueScale,
				en);
		maskDataWeights(en);
		prevEnable=en;
    	int numEn=getNumEnabledSamples(en);
    	if (numEn<minLeftSamples) return Double.NaN;
	}
	double [] focusing_fx= createFXandJacobian(sv, false);
	double rms_pure=       calcErrorDiffY(focusing_fx, true);
	//	System.out.println("rms_pure="+rms_pure);
	return rms_pure;
}
public double [] findAdjustZ(
		FocusingFieldMeasurement measurement,
		boolean filterZ,
		double filterByValueScale,
		double zMin,
		double zMax,
		double zStep,
        double tMin,
        double tMax,
        double tStep){
	double zBest=Double.NaN;
	double tXBest=Double.NaN;
	double tYBest=Double.NaN;
	double bestRMS=Double.NaN;
	int bestEn=0;
	for (double z=zMin;z<=zMax;z+=zStep) for (double tx=tMin;tx<=tMax;tx+=tStep) for (double ty=tMin;ty<=tMax;ty+=tStep){
		double rms=getAdjustRMS(
				measurement,
				filterZ,
				filterByValueScale,
				z,
				tx,
				ty);
		if (((Double.isNaN(bestRMS) || (bestRMS>=rms)) && !Double.isNaN(rms) && (rms>0.0))){
			zBest=z;
			tXBest=tx;
			tYBest=ty;
			bestRMS=rms;
			bestEn=numEnabled(dataWeightsToBoolean());
		}

		if (debugLevel>1) {
			int numEn=numEnabled(dataWeightsToBoolean());
			System.out.println("findAdjustZ(): z="+z+" tx="+tx+" ty="+ty+" rms="+rms+" used "+numEn+" samples");
		}
	}
	if (debugLevel>0) System.out.println("findAdjustZ()-> z="+zBest+" tx="+tXBest+" ty="+tYBest+" (best RMS = "+bestRMS+" used "+bestEn+" samples)");
	double [] result={zBest,tXBest,tYBest};
	return result;
}



public boolean LevenbergMarquardt(
		FocusingFieldMeasurement measurement, // null in calibrate mode
		boolean openDialog,
//		boolean filterZ, // for adjust mode
		int debugLevel){
	boolean calibrate=measurement==null;
	double savedLambda=this.lambda;
	this.debugLevel=debugLevel;
	if (openDialog && !selectLMAParameters()) return false;
	this.startTime=System.nanoTime();
	// create savedVector (it depends on parameter masks), restore from it if aborted
//	fieldFitting.initSampleCorrVector(
//			flattenSampleCoord(), //double [][] sampleCoordinates,
//			getSeriesWeights()); //double [][] sampleSeriesWeights);
//	fieldFitting.setEstimatedZ0( z0_estimates, false); // boolean force)

//	this.savedVector=this.fieldFitting.createParameterVector(sagittalMaster);
//	if (debugDerivativesFxDxDy){
//		compareDrDerivatives(this.savedVector);
//	}
	if (!calibrate) {
		this.currentStrategyStep=-1;
		fieldFitting.selectZTilt(false);
		keepCorrectionParameters=true;
		resetCenter=false;
		if (!openDialog) stopEachStep=false;
	}
	this.iterationStepNumber=0;
	this.firstRMS=-1; //undefined
	while (true) { // loop for all series
		
        if (this.currentStrategyStep>=0){
        	if (!getStrategy(this.currentStrategyStep)) break; //invalid strategy
        }
        if (!keepCorrectionParameters) fieldFitting.resetSampleCorr();
        if (resetCenter){
			if (debugLevel>0) System.out.println("Resetting center: X "+IJ.d2s(currentPX0,2)+" -> "+IJ.d2s(pX0_distortions,2));
			if (debugLevel>0) System.out.println("Resetting center: Y "+IJ.d2s(currentPY0,2)+" -> "+IJ.d2s(pY0_distortions,2));
            currentPX0=pX0_distortions;
            currentPY0=pY0_distortions;
    		fieldFitting.setCenterXY(currentPX0,currentPY0);
        }
//        setDataVector(createDataVector()); //new
       	fieldFitting.initSampleCorrChnParIndex(flattenSampleCoord());
       	if (calibrate) {
       		setDataVector(
       				true, // calibrate mode
       				createDataVector()); // Make it different for adjustment mode
       		fieldFitting.initSampleCorrVector(
       				flattenSampleCoord(), //double [][] sampleCoordinates,
       				getSeriesWeights()); //double [][] sampleSeriesWeights);
       		fieldFitting.setEstimatedZ0( z0_estimates, false); // boolean force)
       	} else {
       		setDataVector(
       				false, // calibrate mode
       				createDataVector(measurement)); // Make it different for adjustment mode
       		if (filterZ) {
       	    	boolean [] en=dataWeightsToBoolean();
       	    	en= filterByZRanges(
       	    			zRanges,
       	    			en);
       	    	maskDataWeights(en);
       	    	prevEnable=en;
       	    	int numEn=getNumEnabledSamples(en);
       	    	if (numEn<minLeftSamples) return false;
       		}
       		if (filterByValueScale>0.0){
       			boolean [] en=dataWeightsToBoolean();
       			en= filterByValue(
       					filterByValueScale,
       					en);
       			maskDataWeights(en);
       			prevEnable=en;
       	    	int numEn=getNumEnabledSamples(en);
       	    	if (numEn<minLeftSamples) return false;
       		}       		
       		
       		fieldFitting.initSampleCorrVector(
       				flattenSampleCoord(), //double [][] sampleCoordinates,
       				null); //getSeriesWeights()); //double [][] sampleSeriesWeights);
//       		fieldFitting.setEstimatedZ0( z0_estimates, false); // boolean force)

       	}

    	this.savedVector=this.fieldFitting.createParameterVector(sagittalMaster);
    	if (debugDerivativesFxDxDy){
    		compareDrDerivatives(this.savedVector);
    	}
       	
    	
    	
    	
		//     while (this.fittingStrategy.isSeriesValid(this.seriesNumber)){ // TODO: Add "stop" tag to series
		this.currentVector=null; // invalidate for the new series
		//         boolean wasLastSeries=false;
		while (true) { // loop for the same series
	        
			boolean [] state=stepLevenbergMarquardtFirst(debugLevel);
			if (state==null) {
				String msg="Calculation aborted by user request, restoring saved parameter vector";
				IJ.showMessage(msg);
				System.out.println(msg);
				commitParameterVector(this.savedVector);
				this.lambda=savedLambda;
				return false;
			}

			if (debugLevel>1) System.out.println(this.currentStrategyStep+":"+this.iterationStepNumber+": stepLevenbergMarquardtFirst("+debugLevel+")==>"+state[1]+":"+state[0]);
			boolean cont=true;
			// Make it success if this.currentRMS<this.firstRMS even if LMA failed to converge
			if (state[1] && !state[0] && (this.firstRMS>this.currentRMS)){
				if (debugLevel>1) System.out.println("LMA failed to converge, but RMS improved from the initial value ("+
						this.currentRMS+" < "+this.firstRMS+"), currentRMSPure="+currentRMSPure+", firstRMSPure="+firstRMSPure);
				state[0]=true;
			}
			if (
					(this.stopRequested.get()>0) || // graceful stop requested
					(this.stopEachStep) ||
					(this.stopEachSeries && state[1]) ||
					(this.stopOnFailure && state[1] && !state[0])){
				if (state[1] && !state[0] && !calibrate){
					return false;
				}

				if (debugLevel>0){
					if (this.stopRequested.get()>0) System.out.println("User requested stop");
					System.out.println("LevenbergMarquardt(): step ="+this.currentStrategyStep+":"+this.iterationStepNumber+
							", RMS="+IJ.d2s(this.currentRMS,8)+
							" ("+IJ.d2s(this.firstRMS,8)+") "+
							") at "+ IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3));
				}
				long startDialogTime=System.nanoTime();
				cont=dialogLMAStep(state);
				this.stopRequested.set(0); // Will not stop each run
				this.startTime+=(System.nanoTime()-startDialogTime); // do not count time used by the User.
				//                 if (this.showThisImages) showDiff (this.currentfX, "fit-"+this.iterationStepNumber);
				//                 if (this.showNextImages) showDiff (this.nextfX, "fit-"+(this.iterationStepNumber+1));
			}
			stepLevenbergMarquardtAction(debugLevel); // apply step - in any case?
			if (this.updateStatus){
				IJ.showStatus("Step #"+this.currentStrategyStep+":"+this.iterationStepNumber+
						" RMS="+IJ.d2s(this.currentRMS,8)+
						" ("+IJ.d2s(this.firstRMS,8)+")"+
						" RMSPure="+IJ.d2s(this.currentRMSPure,8)+
						" ("+IJ.d2s(this.firstRMSPure,8)+")"+
						" ");
			}
			if (!cont){
				if (this.saveSeries) {
					savedLambda=this.lambda;

					this.savedVector=this.currentVector.clone();
					//                        saveFittingSeries(); // will save series even if it ended in failure, vector will be only updated
					//                        updateCameraParametersFromCalculated(true); // update camera parameters from all (even disabled) images
					//                        updateCameraParametersFromCalculated(false); // update camera parameters from enabled only images (may overwrite some of the above)
				}
				// if RMS was decreased. this.saveSeries==false after dialogLMAStep(state) only if "cancel" was pressed
				commitParameterVector(this.savedVector); // either new or original
				this.lambda=savedLambda;
				return this.saveSeries; // TODO: Maybe change result?
			}
			//stepLevenbergMarquardtAction();             
			if (state[1]) {
				if (!state[0]) {
					commitParameterVector(this.savedVector);
					this.lambda=savedLambda;
					return false; // sequence failed
				}
				this.savedVector=this.currentVector.clone();
				//                 saveFittingSeries();
				//                    updateCameraParametersFromCalculated(true); // update camera parameters from all (even disabled) images
				//                    updateCameraParametersFromCalculated(false); // update camera parameters from enabled only images (may overwrite some of the above)
				//                    wasLastSeries=this.fittingStrategy.isLastSeries(this.seriesNumber);
				//                 this.seriesNumber++;
				break; // while (true), proceed to the next series
			}
		} // while true - same series
		//         if (wasLastSeries) break;
		//     } // while (this.fittingStrategy.isSeriesValid(this.seriesNumber)){ // TODO: Add "stop" tag to series
		if (fieldFitting.fieldStrategies.isLast(this.currentStrategyStep)) break;
		String msg="LMA series="+this.currentStrategyStep+ " RMS="+this.currentRMS+" ("+this.firstRMS+") "+
				", pure RMS="+this.currentRMSPure+" ("+this.firstRMSPure+") "+
				" at "+ IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3);
		if (debugLevel>0) System.out.println("stepLevenbergMarquardtAction() "+msg);
		this.currentStrategyStep++;
		this.iterationStepNumber=0;
	} // for all series
	String msg="RMS="+this.currentRMS+" ("+this.firstRMS+") "+
			", pure RMS="+this.currentRMSPure+" ("+this.firstRMSPure+") "+
			" at "+ IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3);
	if (debugLevel>0) System.out.println("stepLevenbergMarquardtAction() "+msg);
	//    	if (this.updateStatus) IJ.showStatus(msg);
	if (this.updateStatus){
		IJ.showStatus("Done: Step #"+this.currentStrategyStep+":"+this.iterationStepNumber+
				" RMS="+IJ.d2s(this.currentRMS,8)+
				" ("+IJ.d2s(this.firstRMS,8)+")"+
				" RMSPure="+IJ.d2s(this.currentRMSPure,8)+
				" ("+IJ.d2s(this.firstRMSPure,8)+")"+
				" ");
	}
	this.savedVector=this.currentVector.clone();
	commitParameterVector(this.savedVector);
	if (calibrate) zRanges=calcZRanges(dataWeightsToBoolean());
	return true; // all series done
}

    public class FocusingFieldMeasurement{
        public String timestamp;
        public double temperature;
        public int [] motors;
        double [][][][] samples = null; // last - psf radius in pixels: {x2,y2, xy}
        
        public FocusingFieldMeasurement(
                String timestamp,
                double temperature,
                int [] motors,
                double [][][][] samples
                ){
            this.timestamp=timestamp;
            this.temperature=temperature;
            this.motors = new int[motors.length];
            for (int i=0;i<motors.length;i++) this.motors[i]= motors[i];
            if (samples !=null) {
//                this.samples=new double[samples.length][samples[0].length][samples[0][0].length][2];
                try {
                    this.samples=new double[samples.length][samples[0].length][samples[0][0].length][3];
                } catch (Exception e){
                    return;
                }
                for (int i=0;i<this.samples.length;i++) for (int j=0;j<this.samples[i].length;j++) for (int c=0;c<this.samples[i][j].length;c++){
                    try {
//                        double rt=samples[i][j][c][0]/Math.sqrt((1+samples[i][j][c][1]*samples[i][j][c][1])/2.0); // tangential;
//                        double rs=rt*samples[i][j][c][1]; // sagittal
//                        this.samples[i][j][c][0]=rs;
//                        this.samples[i][j][c][1]=rt;
                        this.samples[i][j][c]=samples[i][j][c].clone();
                    } catch (Exception e){
                        for (int ii=0;ii<this.samples[i][j][c].length;ii++) this.samples[i][j][c][ii]=Double.NaN;
                    }
                }
            }
        }
        public FocusingFieldMeasurement(
                String timestamp,
                double temperature,
                int [] motors,
//                ArrayList<String> sampleStrings
                String [] sampleStrings
                ){
            this.timestamp=timestamp;
            this.temperature=temperature;
            this.motors = new int[motors.length];
            for (int i=0;i<motors.length;i++) this.motors[i]= motors[i];
            
//            if ((sampleStrings!=null) && (sampleStrings.size()>0)) {
            if ((sampleStrings!=null) && (sampleStrings.length>0)) {
                int maxi=0,maxj=0;
                for (String s:sampleStrings){
                    String [] ps=s.split(regSep);
                    int i=Integer.parseInt(ps[0]);
                    int j=Integer.parseInt(ps[1]);
                    if (i>maxi) maxi=i;
                    if (j>maxj) maxj=j;
                }
                int rows=maxi+1;
                int cols=maxj+1;
                this.samples=new double [rows][cols][][];
                for (String s:sampleStrings){
                    String [] ps=s.split(regSep);
                    int i=Integer.parseInt(ps[0]);
                    int j=Integer.parseInt(ps[1]);
                    int colors=(ps.length-2)/3; // 2;
                    this.samples[i][j]=new double [colors][3]; //[2];
                    for (int c=0;c<colors;c++){
//                        this.samples[i][j][c][0]=Double.parseDouble(ps[2*c+2]);
//                        this.samples[i][j][c][1]=Double.parseDouble(ps[2*c+3]);
                        this.samples[i][j][c][0]=Double.parseDouble(ps[3*c+2]);
                        this.samples[i][j][c][1]=Double.parseDouble(ps[3*c+3]);
                        this.samples[i][j][c][2]=Double.parseDouble(ps[3*c+4]);
                    }
                }                
            }
        }

        public ArrayList<String> asListString(){
        ArrayList<String> nodeList=new ArrayList<String>();
        if (this.samples!=null){
            for (int i=0;i<this.samples.length;i++) for (int j=0;j<this.samples[i].length;j++){
                     String sdata=i+sep+j;
                     for (int c=0;c<samples[i][j].length;c++) {
//                         sdata += sep+this.samples[i][j][c][0]+ // sagittal
//                                 sep+this.samples[i][j][c][1]; // tangential
                         sdata += sep+this.samples[i][j][c][0]+ // x2
                                 sep+this.samples[i][j][c][1]+ // y2
                                 sep+this.samples[i][j][c][2]; // xy
                     }
                     nodeList.add(sdata);
            }
        }
            return nodeList;
        }
    }
    public FocusingField(
        String serialNumber,
            String lensSerial, // if null - do not add average
        String comment,
            double pX0,
            double pY0,
            double [][][] sampleCoord, //){ // x,y,r
            AtomicInteger stopRequested){
    	setDefaults();
        this.serialNumber=serialNumber;
        this.lensSerial=lensSerial;
        this.comment=comment;
        this.pX0_distortions=pX0;
        this.pY0_distortions=pY0;
        // copy distortions to current PX0/PY0
//        this.currentPX0=pX0_distortions;
//        this.currentPY0=pY0_distortions;
        this.sampleCoord=sampleCoord;
        this.measurements=new ArrayList<FocusingFieldMeasurement>();
        this.stopRequested=stopRequested;
    }
    
    public FocusingField(
    		boolean smart, // do not open dialog if default matches
    		String defaultPath, //){
    		AtomicInteger stopRequested){
    	setDefaults();
    	this.stopRequested=stopRequested;
    	loadXML(smart,defaultPath);
    }
    public void addSample(
    		String timestamp,
    		double temperature,
    		int [] motors,
    		double [][][][] samples
    		)
    {
    	measurements.add(new FocusingFieldMeasurement(
    			timestamp,
    			temperature,
    			motors,
    			samples
    			));
    }
    public boolean loadXML(
            boolean smart, // do not open dialog if default matches
            String defaultPath){ // x,y,r
        String [] extensions={".history-xml"};
        CalibrationFileManagement.MultipleExtensionsFileFilter parFilter = new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,"*.history-xml files");
        String pathname=CalibrationFileManagement.selectFile(
                smart,
                false,
                "Restore focusing field measurement data",
                "Restore",
                parFilter,
                defaultPath); //String defaultPath
        if ((pathname==null) || (pathname=="")) return false;
    XMLConfiguration hConfig=null;
    try {
            hConfig=new XMLConfiguration(pathname);
        } catch (ConfigurationException e) {
            return false;
        }
        hConfig.setThrowExceptionOnMissing(false); // default value, will return null on missing
        comment= hConfig.getString("comment","no comments");
//        if ((comment.length()>10) && comment.substring(0,9).equals("<![CDATA[")) comment=comment.substring(9,comment.length()-3);
        
        serialNumber= hConfig.getString("serialNumber","???");
        lensSerial= hConfig.getString("lensSerial","???");
        pX0_distortions=Double.parseDouble(hConfig.getString("lens_center_x","0.0"));
        pY0_distortions=Double.parseDouble(hConfig.getString("lens_center_y","0.0"));
        // copy distortions to current PX0/PY0
        this.currentPX0=pX0_distortions;
        this.currentPY0=pY0_distortions;
        int rows=Integer.parseInt(hConfig.getString("samples_y","0"));
        int cols=Integer.parseInt(hConfig.getString("samples_x","0"));
        sampleCoord=new double [rows][cols][2];
        for (int i=0;i<rows;i++) for (int j=0;j<cols;j++){
            String [] coords= hConfig.getString("sample_"+i+"_"+j,"0 0").split(regSep);
            sampleCoord[i][j][0]=Double.parseDouble(coords[0]);
            sampleCoord[i][j][1]=Double.parseDouble(coords[1]);
        }
        int numMeasurements=Integer.parseInt(hConfig.getString("measurements","0"));
        measurements=new ArrayList<FocusingFieldMeasurement>();
        for (int m=0;m<numMeasurements;m++){
         String prefix="measurement_"+m+".";

            String timestamp= hConfig.getString(prefix+"timestamp","0");
            double temperature=Double.parseDouble(hConfig.getString(prefix+"temperature","0.0"));
            String [] sMotors=hConfig.getString(prefix+"motors","0 0 0").split(regSep);
            int [] motors=new int [sMotors.length];
            for (int i=0;i<sMotors.length;i++) motors[i]=Integer.parseInt(sMotors[i]);
            String [] sampleStrings = hConfig. getStringArray(prefix+"sample");
            measurements.add(new FocusingFieldMeasurement(
                    timestamp,
                    temperature,
                    motors,
                    sampleStrings));
        }
        return true;
    }
    public void saveXML(
            String path){ // x,y,r
    XMLConfiguration hConfig=new XMLConfiguration();
    hConfig.setRootElementName("focusingHistory");
    if (comment!=null){
        String comment_esc=comment.replace(",","\\,");
//          hConfig.addProperty("comment","<![CDATA["+comment_esc+ "]]>");
         hConfig.addProperty("comment",comment_esc);
    }
    if (serialNumber!=null) hConfig.addProperty("serialNumber",serialNumber);
    if (lensSerial!=null) hConfig.addProperty("lensSerial",lensSerial);
     hConfig.addProperty("lens_center_x",pX0_distortions); // distortions center, not aberrations!
     hConfig.addProperty("lens_center_y",pY0_distortions);
     if ((sampleCoord!=null) && (sampleCoord.length>0) && (sampleCoord[0] != null) && (sampleCoord[0].length>0)){
         hConfig.addProperty("samples_x",sampleCoord[0].length);
         hConfig.addProperty("samples_y",sampleCoord.length);
         for (int i=0;i<sampleCoord.length;i++)
             for (int j=0;j<sampleCoord[i].length;j++){
//          double coord[] = {sampleCoord[i][j][0],sampleCoord[i][j][1]};
         hConfig.addProperty("sample_"+i+"_"+j,sampleCoord[i][j][0]+sep+sampleCoord[i][j][1]);
             }
     }
     hConfig.addProperty("measurements",this.measurements.size());
     for (int i=0;i<this.measurements.size();i++){
         FocusingFieldMeasurement meas=this.measurements.get(i);
         String prefix="measurement_"+i+".";
         if (meas.timestamp!=null) hConfig.addProperty(prefix+"timestamp",meas.timestamp);
         hConfig.addProperty(prefix+"temperature",meas.temperature);
         hConfig.addProperty(prefix+"motors",meas.motors[0]+sep+meas.motors[1]+sep+meas.motors[2]);
             hConfig.addProperty(prefix+"sample",meas.asListString());
     }     
     File file=new File (path);
     BufferedWriter writer;
     try {
         writer = new BufferedWriter(new FileWriter(file));
         hConfig.save(writer);
     } catch (IOException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
     } catch (ConfigurationException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
     }
    }
    public void testMeasurement(){
        GenericDialog gd = new GenericDialog("Select measurement");
        int nMeas=measurements.size()/2;
//        double zMin=-40.0;
//        double zMax= 40.0;
//        double zStep=2.0;
        double targetTiltX=0.0; // for testing, normally should be 0 um/mm
        double targetTiltY=0.0; // for testing, normally should be 0 um/mm
        
//    	fieldFitting.mechanicalFocusingModel.setZTxTy(0.0,0.0,0.0); // to correctly find Z centers,
    	fieldFitting.mechanicalFocusingModel.setAdjustMode(false); // to correctly find Z centers,
        
        double [] center_z=fieldFitting.getZCenters(false);
        double [] centerFWHM={
        		fieldFitting.getCalcValuesForZ(center_z[0],0.0,null)[1],
        		fieldFitting.getCalcValuesForZ(center_z[1],0.0,null)[3],
        		fieldFitting.getCalcValuesForZ(center_z[2],0.0,null)[5]
        };
//        String path=null;
        String title="Test adjustment results";
        double [] best_qb_corr= fieldFitting.getBestQualB(
                k_red,
                k_blue,
                true);
        gd.addMessage("Best composite distance for FWHM^4 "+ IJ.d2s(best_qb_corr[0],3)+"  um"+
        		", FWHM="+IJ.d2s(best_qb_corr[1],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/best_qb_corr[1],2)+" lp/mm");
        gd.addMessage("Best center focus for Red (relative to best composite) = "+ IJ.d2s(center_z[0]-best_qb_corr[0],3)+" um"+
        		", FWHM="+IJ.d2s(centerFWHM[0],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/centerFWHM[0],2)+" lp/mm");
        gd.addMessage("Best center focus for Green (relative to best composite) = "+ IJ.d2s(center_z[1]-best_qb_corr[0],3)+" um"+
        		", FWHM="+IJ.d2s(centerFWHM[1],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/centerFWHM[1],2)+" lp/mm");
        gd.addMessage("Best center focus for Blue (relative to best composite) = "+ IJ.d2s(center_z[2]-best_qb_corr[0],3)+" um"+
        		", FWHM="+IJ.d2s(centerFWHM[2],3)+"um, MTF50="+IJ.d2s(fwhm_to_mtf50/centerFWHM[2],2)+" lp/mm");
        
        gd.addNumericField("Measurement number",nMeas,0,5,"0.."+(measurements.size()-1));
    	if (fieldFitting.channelSelect==null) fieldFitting.channelSelect=fieldFitting.getDefaultMask();
    	for (int i=0;i<fieldFitting.channelSelect.length;i++) {
    		gd.addCheckbox(fieldFitting.getDescription(i), fieldFitting.channelSelect[i]);                    
    	}
        
//    	filterZ=true;      // (adjustment mode)filter samples by Z
//    	minLeftSamples=10;  // minimal number of samples (channel/dir/location) for adjustment
        
        
        gd.addCheckbox("Filter samples/channels by Z",filterZ);
        gd.addNumericField("Filter by value (remove samples above scaled best FWHM for channel/location)",filterByValueScale,2,5,"x");
        gd.addNumericField("Minimal required number of channels/samples",minLeftSamples,0,3,"samples");

        gd.addNumericField("Z min",zMin,2,5,"um");
        gd.addNumericField("Z max",zMax,2,5,"um");
        gd.addNumericField("Z step",zStep,2,5,"um");

        gd.addNumericField("Tilt min",tMin,2,5,"um/mm");
        gd.addNumericField("Tilt max",tMax,2,5,"um/mm");
        gd.addNumericField("Tilt step",tStep,2,5,"um/mm");
        
        gd.addNumericField("Target focus (relative to best composirte)",targetRelFocalShift,2,5,"um");
        gd.addNumericField("Target horizontal tilt (normally 0)",targetTiltX,2,5,"um/mm");
        gd.addNumericField("Target vertical tilt (normally 0)",targetTiltY,2,5,"um/mm");
    	WindowTools.addScrollBars(gd);
    	gd.showDialog();
    	if (gd.wasCanceled()) return;
    	
    	nMeas=(int) gd.getNextNumber();

    	for (int i=0;i<fieldFitting.channelSelect.length;i++) {
    		fieldFitting.channelSelect[i]=gd.getNextBoolean();
    	}

    	filterZ=gd.getNextBoolean();
    	filterByValueScale=gd.getNextNumber();
        minLeftSamples=(int) gd.getNextNumber();
    	
        zMin= gd.getNextNumber();
        zMax= gd.getNextNumber();
        zStep=gd.getNextNumber();

        tMin= gd.getNextNumber();
        tMax= gd.getNextNumber();
        tStep=gd.getNextNumber();

        targetRelFocalShift=gd.getNextNumber();
        targetTiltX=gd.getNextNumber(); // for testing, normally should be 0 um/mm
        targetTiltY=gd.getNextNumber(); // for testing, normally should be 0 um/mm
       
        boolean OK;
    	fieldFitting.mechanicalFocusingModel.setAdjustMode(true); // to correctly find Z centers,
    	fieldFitting.mechanicalSelect=fieldFitting.mechanicalFocusingModel.maskSetZTxTy();

        String header="# measurement";
		for (int i=0;i<fieldFitting.mechanicalFocusingModel.paramValues.length;i++){
			if ((fieldFitting.mechanicalSelect==null) || fieldFitting.mechanicalSelect[i] ) {
				header+="\t"+fieldFitting.mechanicalFocusingModel.getName(i)+" ("+fieldFitting.mechanicalFocusingModel.getUnits(i)+")";
			}
		}
        for (int i=0;i<3;i++) header+="\tmz"+i;
        for (int i=0;i<3;i++) header+="\tm"+i;
        StringBuffer sb = new StringBuffer();


        boolean single= (nMeas>=0);
    	if (single){
    		OK=testMeasurement(
    				measurements.get(nMeas),    				
//    				nMeas,
    				zMin, //+best_qb_corr[0],
			        zMax, //+best_qb_corr[0],
			        zStep, 
    				tMin,
			        tMax,
			        tStep 
			        );
    		if (!OK){
    			if (debugLevel>0) System.out.println("testMeasurement("+nMeas+") failed");
    		} else {
    			
        		if (debugLevel>0) System.out.print("======== testMeasurement("+nMeas+") ========");
    			for (int i=0;i<fieldFitting.mechanicalFocusingModel.paramValues.length;i++){
    				if ((fieldFitting.mechanicalSelect==null) || fieldFitting.mechanicalSelect[i] ) {
    					System.out.println(
    							fieldFitting.mechanicalFocusingModel.getDescription(i)+": "+
    									IJ.d2s(fieldFitting.mechanicalFocusingModel.paramValues[i],3)+" "+
    									fieldFitting.mechanicalFocusingModel.getUnits(i));
    				}
    			}
    			double [] dmz= getAdjustedMotors(
    		            targetRelFocalShift+best_qb_corr[0],
    		            targetTiltX, // for testing, normally should be 0 um/mm
    		            targetTiltY,
    		            false);
    			if ((dmz!=null) && (debugLevel>0)){
    				System.out.println("Suggested motor linearized positions: "+IJ.d2s(dmz[0],2)+":"+IJ.d2s(dmz[1],2)+":"+IJ.d2s(dmz[2],2));
    			}
    			double [] dm= getAdjustedMotors(
    		            targetRelFocalShift+best_qb_corr[0],
    		            targetTiltX, // for testing, normally should be 0 um/mm
    		            targetTiltY,
    		            true);
    			if ((dm!=null) && (debugLevel>0)){
    				System.out.println("Suggested motor positions: "+IJ.d2s(dm[0],0)+":"+IJ.d2s(dm[1],0)+":"+IJ.d2s(dm[2],0));
    			}
    			sb.append(nMeas);
    			for (int i=0;i<fieldFitting.mechanicalFocusingModel.paramValues.length;i++){
    				if ((fieldFitting.mechanicalSelect==null) || fieldFitting.mechanicalSelect[i] ) {
    					sb.append("\t"+IJ.d2s(fieldFitting.mechanicalFocusingModel.paramValues[i],3));
    				}
    			}
    			if (dmz!=null){
    				for (int i=0;i<dmz.length;i++) sb.append("\t"+IJ.d2s(dmz[i],1)); 
    			} else {
    				sb.append("\t---\t---\t---");
    			}
    			if (dm!=null){
    				for (int i=0;i<dm.length;i++) sb.append("\t"+IJ.d2s(dm[i],1)); 
    			} else {
    				sb.append("\t---\t---\t---");
    			}
    			sb.append("\n");
    		}
    	} else {
    		for (nMeas=0;nMeas<measurements.size();nMeas++){
        		if (debugLevel>0) System.out.print("======== testMeasurement("+nMeas+") ========      ");
    			OK=testMeasurement(
    					measurements.get(nMeas),
    			        zMin, //+best_qb_corr[0],
    			        zMax, // +best_qb_corr[0],
    			        zStep,
        				tMin,
    			        tMax,
    			        tStep);
        		if (!OK){
        			if (debugLevel>0) System.out.println("testMeasurement("+nMeas+") failed");
        		} else {
        			for (int i=0;i<fieldFitting.mechanicalFocusingModel.paramValues.length;i++){
        				if ((fieldFitting.mechanicalSelect==null) || fieldFitting.mechanicalSelect[i] ) {
        					System.out.println(
        							fieldFitting.mechanicalFocusingModel.getDescription(i)+": "+
        									IJ.d2s(fieldFitting.mechanicalFocusingModel.paramValues[i],3)+" "+
        									fieldFitting.mechanicalFocusingModel.getUnits(i));
        				}
        			}
        			double [] dmz= getAdjustedMotors(
        		            targetRelFocalShift+best_qb_corr[0],
        		            targetTiltX, // for testing, normally should be 0 um/mm
        		            targetTiltY,
        		            false);
        			if ((dmz!=null) && (debugLevel>0)){
        				System.out.println("Suggested motor linearized positions: "+IJ.d2s(dmz[0],2)+":"+IJ.d2s(dmz[1],2)+":"+IJ.d2s(dmz[2],2));
        			}
        			double [] dm= getAdjustedMotors(
        		            targetRelFocalShift+best_qb_corr[0],
        		            targetTiltX, // for testing, normally should be 0 um/mm
        		            targetTiltY,
        		            true);
        			if ((dm!=null) && (debugLevel>0)){
        				System.out.println("Suggested motor positions: "+IJ.d2s(dm[0],0)+":"+IJ.d2s(dm[1],0)+":"+IJ.d2s(dm[2],0));
        			}
        			sb.append(nMeas);
        			for (int i=0;i<fieldFitting.mechanicalFocusingModel.paramValues.length;i++){
        				if ((fieldFitting.mechanicalSelect==null) || fieldFitting.mechanicalSelect[i] ) {
        					sb.append("\t"+IJ.d2s(fieldFitting.mechanicalFocusingModel.paramValues[i],3));
        				}
        			}
        			if (dmz!=null){
        				for (int i=0;i<dmz.length;i++) sb.append("\t"+IJ.d2s(dmz[i],1)); 
        			} else {
        				sb.append("\t---\t---\t---");
        			}
        			if (dm!=null){
        				for (int i=0;i<dm.length;i++) sb.append("\t"+IJ.d2s(dm[i],1)); 
        			} else {
        				sb.append("\t---\t---\t---");
        			}
        			sb.append("\n");
        		}
    		}
    	}
    	if (!single) {
//    		if (path!=null) {
//    			CalibrationFileManagement.saveStringToFile (
//    					path,
//    					header+"\n"+sb.toString());
//    		} else {
    			new TextWindow(title, header, sb.toString(), 800,1000);
 //   		}
    	}
//    	fieldFitting.mechanicalFocusingModel.setZTxTy(0.0,0.0,0.0); // restore zeros to correctly find Z centers,
    	fieldFitting.mechanicalFocusingModel.setAdjustMode(false); // to correctly find Z centers,
    }
    public double [] getAdjustedMotors(
            double targetRelFocalShift,
            double targetTiltX, // for testing, normally should be 0 um/mm
            double targetTiltY,
            boolean motorSteps){ // for testing, normally should be 0 um/mm
    	double [] zM=fieldFitting.mechanicalFocusingModel.getZM(
    			currentPX0, //fieldFitting.getCenterXY()[0],
    			currentPY0, //fieldFitting.getCenterXY()[1],
    			targetRelFocalShift,
    			targetTiltX,
    			targetTiltY);
    	if (zM==null) return null; // not yet used
    	if (!motorSteps) return zM;
//		if (debugLevel>0) System.out.println("Suggested motor linearized positions: "+IJ.d2s(zM[0],2)+":"+IJ.d2s(zM[1],2)+":"+IJ.d2s(zM[2],2));
    	double [] dm=new double[zM.length];
    	for (int index=0;index<dm.length;index++){
    		dm[index]=fieldFitting.mechanicalFocusingModel.zmToM(
        			zM[index],
        			index);
    	}
    	return dm;
    }
    
    public boolean testMeasurement(
    		FocusingFieldMeasurement measurement, // null in calibrate mode
//    		int nMeas,
            double zMin,
            double zMax,
            double zStep,
            double tMin,
            double tMax,
            double tStep
            ){
    	int retryLimit=20;
    	fieldFitting.mechanicalFocusingModel.setAdjustMode(true);
    	setDataVector(
    			false,
    			createDataVector(measurement)); //measurements.get(nMeas)));
// TODO: Adjust by center samples  only?
    	double [] zTxTy=findAdjustZ(
//    			measurements.get(nMeas),
    			measurement,
    			filterZ, //boolean filterZ,
    			filterByValueScale,
    			zMin,
    			zMax,
    			zStep,
                tMin,
                tMax,
                tStep);
    	if (Double.isNaN(zTxTy[0])) {
    		if (debugLevel>1) System.out.println("testMeasurement() failed: insufficient data to get z initial estimation");
    		return false;
    	}
    	fieldFitting.mechanicalFocusingModel.setZTxTy(zTxTy[0],zTxTy[1],zTxTy[2]); //z,0.0,0.0);// z,tx,ty
    	boolean [] wasPrevEnable=null;
    	for (int n=0;n<retryLimit;n++) { // TODO: Watch for the mask remain stable
    		zTxTy[0]=fieldFitting.mechanicalFocusingModel.getValue(MECH_PAR.z0);
    		zTxTy[1]=fieldFitting.mechanicalFocusingModel.getValue(MECH_PAR.tx);
    		zTxTy[2]=fieldFitting.mechanicalFocusingModel.getValue(MECH_PAR.ty);
    		if (debugLevel>0) System.out.println("testMeasurement(), run "+n+" (z="+zTxTy[0]+" tx="+zTxTy[1]+" ty="+zTxTy[2]+")");
    		boolean [] was2PrevEnable=(wasPrevEnable==null)?null:wasPrevEnable.clone();
    		wasPrevEnable=(prevEnable==null)?null:prevEnable.clone();
    		boolean OK=LevenbergMarquardt(
    				measurement, 
    				false, // true, // open dialog
//    				filterZ, // filterZ
    				debugLevel);
    		if (!OK){
//        		if (debugLevel>0) System.out.println("testMeasurement("+nMeas+") failed");
        		if (debugLevel>1) System.out.println("testMeasurement() failed: LMA failed");
        		return false;
    		}
/*    		
    		for (int i=0;i<fieldFitting.mechanicalFocusingModel.paramValues.length;i++){
    			if ((fieldFitting.mechanicalSelect==null) || fieldFitting.mechanicalSelect[i] ) {
    				System.out.println(
    						fieldFitting.mechanicalFocusingModel.getDescription(i)+": "+
    								IJ.d2s(fieldFitting.mechanicalFocusingModel.paramValues[i],3)+" "+
    								fieldFitting.mechanicalFocusingModel.getUnits(i));
    			}
    		}
*/    		
    		if ((wasPrevEnable!=null) && (prevEnable!=null) && (wasPrevEnable.length==prevEnable.length)){
    			boolean changedEnable=false;
    			for (int i=0;i<prevEnable.length;i++) if (prevEnable[i]!=wasPrevEnable[i]){
    				changedEnable=true;
    				break;
    			}
    			if (!changedEnable) {
    				if (debugLevel>0) System.out.println("No filter chnage, finished in "+(n+1)+" steps");
    				return true;
    			} else {
    				if ((was2PrevEnable!=null) && (prevEnable!=null) && (was2PrevEnable.length==prevEnable.length)){
    					changedEnable=false;
    					for (int i=0;i<prevEnable.length;i++) if (prevEnable[i]!=was2PrevEnable[i]){
    						changedEnable=true;
    						break;
    					}
    					if (!changedEnable) {
    						if (debugLevel>0) System.out.println("Filter repeats one before previous, finished in "+(n+1)+" steps");
    						return true;
    					}
    				}
    			}
    		}
    	}
		if (debugLevel>0) System.out.println("Maximal retries exceeded in "+retryLimit+" steps");
    	return true; //?
    }
    public class FieldFitting{
 //   	private Properties savedProperties=null;
    	private double [] pXY=null;
    	private boolean [] centerSelect=null;
    	private boolean [] centerSelectDefault={true,true};
    	private MechanicalFocusingModel mechanicalFocusingModel;
    	private CurvatureModel [] curvatureModel=new CurvatureModel[6]; // 3 colors, sagittal+tangential each
    	private boolean [] channelSelect=null;
    	private boolean [] mechanicalSelect=null;
    	private boolean [][] curvatureSelect=new boolean[6][];

    	private boolean [][] sampleCorrSelect= new boolean[6][]; // enable individual (per sample coordinates) correction of parameters
    	private double [][] sampleCorrCost= new double[6][]; // equivalent cost of one unit of parameter value (in result units, um)
    	private double [][] sampleCorrSigma= new double[6][]; // sigma (in mm) for neighbors influence
    	private double [][] sampleCorrPullZero=new double[6][]; // 1.0 - only difference from neighbors matters, 0.0 - only difference from 0
//    	private String strategyComment="";
//    	private boolean lastInSeries=true;
//    	private double lambda=0.001; // synchronize with top?
    	public FieldStrategies fieldStrategies;
    	
//    	private double [] sampleCorrRadius=null;
    	private double [][] sampleCoordinates=null;
    	private double [][][][] sampleCorrCrossWeights= new double[6][][][];
    	private double [] sampleCorrVector=null; // currently adjusted per-sample parameters
    	private double [][][] correctionParameters=new double[6][][]; // all
    	public int numberOfLocations=0;

    	private int [][] sampleCorrChnParIndex=null;
    	private boolean [] dflt_sampleCorrSelect= {false,false,false,false};
//    	private double [] dflt_sampleCorrCost= {1.0,50.0,1.0,1.0};
//    	private double [] dflt_sampleCorrCost= {0.3,20.0,0.3,1.0};
//    	private double [] dflt_sampleCorrCost= {0.05,10.0,2.0,1.0};
//    	private double [] dflt_sampleCorrCost= {0.05,1.0,1.0,1.0};
//    	private double [] dflt_sampleCorrCost= {0.5,0.5,2.0,1.0};
//    	private double [] dflt_sampleCorrCost= {0.1,0.5,2.0,1.0};
//    	private double [] dflt_sampleCorrCost= {0.1,1.0,1.0,0.5};
//    	private double [] dflt_sampleCorrCost= {0.2,2.0,2.0,1.0,1.0};
//    	private double [] dflt_sampleCorrCost= {0.1,1.0,1.0,1.0,1.0};
    	private double [] dflt_sampleCorrCost= {0.01,1.0,1.0,1.0,1.0};
    	private double dflt_sampleCorrSigma= 2.0; // mm
    	private double dflt_sampleCorrPullZero= 0.75; // fraction
    	public final String [] channelDescriptions={
    			"Red, sagittal","Red, tangential",
    			"Green, sagittal","Green, tangential",
    			"Blue, sagittal","Blue, tangential"};

        public void setProperties(String prefix,Properties properties){
        	if (mechanicalFocusingModel==null){
        		if (debugLevel>1) System.out.println ("Mechanical properties not yet initialized, will save properties later");
        		return;
        	}
        	properties.setProperty(prefix+"numberOfLocations",numberOfLocations+"");
        	properties.setProperty(prefix+"centerSelect_X",centerSelect[0]+"");
        	properties.setProperty(prefix+"centerSelect_Y",centerSelect[1]+"");
        	mechanicalFocusingModel.setProperties(prefix+"mechanicalFocusingModel.",properties);
        	for (int i=0;i<curvatureModel.length;i++){
        		if (curvatureModel[i]!=null) curvatureModel[i].setProperties(prefix+"curvatureModel_"+i+".",properties);
        	}
        	if (channelSelect!=null) for (int i=0;i<channelSelect.length;i++){
        		properties.setProperty(prefix+"channelSelect_"+i,channelSelect[i]+"");
        	}
        	if (mechanicalSelect!=null) for (int i=0;i<mechanicalSelect.length;i++){
        		properties.setProperty(prefix+"mechanicalSelect_"+i,mechanicalSelect[i]+"");
        	}
        	for (int chn=0;chn<curvatureSelect.length;chn++) if (curvatureSelect[chn]!=null) for (int i=0;i<curvatureSelect[chn].length;i++){
        		properties.setProperty(prefix+"curvatureSelect_"+chn+"_"+i,curvatureSelect[chn][i]+"");
        	}
        	for (int chn=0;chn<sampleCorrSelect.length;chn++) if (sampleCorrSelect[chn]!=null) for (int i=0;i<sampleCorrSelect[chn].length;i++){
        		properties.setProperty(prefix+"sampleCorrSelect_"+chn+"_"+i,sampleCorrSelect[chn][i]+"");
        	}
        	for (int chn=0;chn<sampleCorrCost.length;chn++) if (sampleCorrCost[chn]!=null) for (int i=0;i<sampleCorrCost[chn].length;i++){
        		properties.setProperty(prefix+"sampleCorrCost_"+chn+"_"+i,sampleCorrCost[chn][i]+"");
        	}
        	for (int chn=0;chn<sampleCorrSigma.length;chn++) if (sampleCorrSigma[chn]!=null) for (int i=0;i<sampleCorrSigma[chn].length;i++){
        		properties.setProperty(prefix+"sampleCorrSigma_"+chn+"_"+i,sampleCorrSigma[chn][i]+"");
        	}
        	for (int chn=0;chn<sampleCorrPullZero.length;chn++) if (sampleCorrPullZero[chn]!=null) for (int i=0;i<sampleCorrPullZero[chn].length;i++){
        		properties.setProperty(prefix+"sampleCorrPullZero_"+chn+"_"+i,sampleCorrPullZero[chn][i]+"");
        	}
        	// save correction parameters values
//        	private double [][][] correctionParameters=new double[6][][]; // all
        	if (correctionParameters!=null){
        		for (int chn=0;chn<correctionParameters.length; chn++) if (correctionParameters[chn]!=null){
        			for (int np=0;np<correctionParameters[chn].length;np++) if (correctionParameters[chn][np]!=null){
        				for (int i=0;i<correctionParameters[chn][np].length;i++){
        	        		properties.setProperty(prefix+"correctionParameters_"+chn+"_"+np+"_"+i,correctionParameters[chn][np][i]+"");
        				}
        			}
        		}
        	}
            fieldStrategies.setProperties(prefix+"fieldStrategies.",properties);
        }
        public void getProperties(String prefix,Properties properties){
        	if (properties.getProperty(prefix+"numberOfLocations")!=null)
        		numberOfLocations=Integer.parseInt(properties.getProperty(prefix+"numberOfLocations"));

        	if (properties.getProperty(prefix+"centerSelect_X")!=null)
        		centerSelect[0]=Boolean.parseBoolean(properties.getProperty(prefix+"centerSelect_X"));
        	if (properties.getProperty(prefix+"centerSelect_Y")!=null)
        		centerSelect[1]=Boolean.parseBoolean(properties.getProperty(prefix+"centerSelect_Y"));
        	if (mechanicalFocusingModel==null){
        		if (debugLevel>1) System.out.println ("Mechanical properties not yet initialized, will apply properties later");
        		return;
        	}
        	mechanicalFocusingModel.getProperties(prefix+"mechanicalFocusingModel.",properties);
        	for (int i=0;i<curvatureModel.length;i++){
        		if (curvatureModel[i]!=null) curvatureModel[i].getProperties(prefix+"curvatureModel_"+i+".",properties);
        	}
        	if (channelSelect==null) {
        		channelSelect=new boolean [6];            	 
        		for (int i=0;i<channelSelect.length;i++)channelSelect[i]=true;
        	}
        	for (int i=0;i<channelSelect.length;i++) if (properties.getProperty(prefix+"channelSelect_"+i)!=null) {
        		channelSelect[i]=Boolean.parseBoolean(properties.getProperty(prefix+"channelSelect_"+i));
        	}
        	if ((mechanicalSelect==null) || (mechanicalSelect.length!=mechanicalFocusingModel.getDefaultMask().length)){
        		mechanicalSelect=mechanicalFocusingModel.getDefaultMask();
        	}
        	for (int i=0;i<mechanicalSelect.length;i++) if (properties.getProperty(prefix+"mechanicalSelect_"+i)!=null) {
        		mechanicalSelect[i]=Boolean.parseBoolean(properties.getProperty(prefix+"mechanicalSelect_"+i));
        	}
        	// curvature parameter selection: first index : channel, inner index - parameter number (radial fast, z - outer)
        	if (curvatureSelect==null) {
        		curvatureSelect=new boolean [curvatureModel.length][];
        		for (int i=0;i<curvatureSelect.length;i++) curvatureSelect[i]=null;
        	}
        	for (int chn=0;chn<correctionParameters.length; chn++) if (curvatureModel[chn]!=null){
        		if ((curvatureSelect[chn]==null) || (curvatureSelect[chn].length!=curvatureModel[chn].getDefaultMask().length)){
        			curvatureSelect[chn]=curvatureModel[chn].getDefaultMask();
        		}
        		for (int i=0;i<curvatureSelect[chn].length;i++) if (properties.getProperty(prefix+"curvatureSelect_"+chn+"_"+i)!=null){
        			curvatureSelect[chn][i]=Boolean.parseBoolean(properties.getProperty(prefix+"curvatureSelect_"+chn+"_"+i));
        		}
        	}
        	
// get correction setup parameters:
        	if (sampleCorrSelect==null){
        		sampleCorrSelect=new boolean [curvatureModel.length][];
        		for (int i=0;i<sampleCorrSelect.length;i++) sampleCorrSelect[i]=null;
        	}
        	for (int chn=0;chn<curvatureModel.length;chn++) if (curvatureModel[chn]!=null){
        		if ((sampleCorrSelect[chn]==null) || (sampleCorrSelect[chn].length!=getDefaultSampleCorrSelect().length)){
        			sampleCorrSelect[chn]=getDefaultSampleCorrSelect();
        		}
        		for (int i=0;i<sampleCorrSelect[chn].length;i++) if (properties.getProperty(prefix+"sampleCorrSelect_"+chn+"_"+i)!=null){
        			sampleCorrSelect[chn][i]=Boolean.parseBoolean(properties.getProperty(prefix+"sampleCorrSelect_"+chn+"_"+i));
        		}
        	}
        		
        	if (sampleCorrCost==null){
        		sampleCorrCost=new double[curvatureModel.length][];
        		for (int i=0;i<sampleCorrCost.length;i++) sampleCorrCost[i]=null;
        	}
        	for (int chn=0;chn<curvatureModel.length;chn++) if (curvatureModel[chn]!=null){
        		if ((sampleCorrCost[chn]==null) || (sampleCorrCost[chn].length!=getDefaultSampleCorrCost().length)){
        			sampleCorrCost[chn]=getDefaultSampleCorrCost();
        		}
        		for (int i=0;i<sampleCorrCost[chn].length;i++) if (properties.getProperty(prefix+"sampleCorrCost_"+chn+"_"+i)!=null){
        			sampleCorrCost[chn][i]=Double.parseDouble(properties.getProperty(prefix+"sampleCorrCost_"+chn+"_"+i));
        		}
        	}

        	if (sampleCorrSigma==null){
        		sampleCorrSigma=new double[curvatureModel.length][];
        		for (int i=0;i<sampleCorrSigma.length;i++) sampleCorrSigma[i]=null;
        	}
        	for (int chn=0;chn<curvatureModel.length;chn++) if (curvatureModel[chn]!=null){
        		if ((sampleCorrSigma[chn]==null) || (sampleCorrSigma[chn].length!=getDefaultSampleCorrSigma().length)){
        			sampleCorrSigma[chn]=getDefaultSampleCorrSigma();
        		}
        		for (int i=0;i<sampleCorrSigma[chn].length;i++) if (properties.getProperty(prefix+"sampleCorrSigma_"+chn+"_"+i)!=null){
        			sampleCorrSigma[chn][i]=Double.parseDouble(properties.getProperty(prefix+"sampleCorrSigma_"+chn+"_"+i));
        		}
        	}

        	if (sampleCorrPullZero==null){
        		sampleCorrPullZero=new double[curvatureModel.length][];
        		for (int i=0;i<sampleCorrPullZero.length;i++) sampleCorrPullZero[i]=null;
        	}
        	for (int chn=0;chn<curvatureModel.length;chn++) if (curvatureModel[chn]!=null){
        		if ((sampleCorrPullZero[chn]==null) || (sampleCorrPullZero[chn].length!=getDefaultCorrPullZero().length)){
        			sampleCorrPullZero[chn]=getDefaultCorrPullZero();
        		}
        		for (int i=0;i<sampleCorrPullZero[chn].length;i++) if (properties.getProperty(prefix+"sampleCorrPullZero_"+chn+"_"+i)!=null){
        			sampleCorrPullZero[chn][i]=Double.parseDouble(properties.getProperty(prefix+"sampleCorrPullZero_"+chn+"_"+i));
        		}
        	}
        	
        	//  read/restore correction parameters values
        	if (correctionParameters==null){
        		correctionParameters=new double[6][][];
        		for (int i=0;i<correctionParameters.length;i++) correctionParameters[i]=null;
        	}
        	if (numberOfLocations>0) {
        		for (int chn=0;chn<correctionParameters.length; chn++) if (curvatureModel[chn]!=null){
        			int numPars=curvatureModel[chn].getNumPars()[0]; // number of Z-parameters
        			if ((correctionParameters[chn]==null) || (correctionParameters[chn].length!=numPars)){
        				correctionParameters[chn]=new double [numPars][numberOfLocations]; // numberOfLocations==0 ?
        				for (int np=0;np<numPars;np++) for (int i=0;i<numberOfLocations;i++)
        					correctionParameters[chn][np][i]=0.0;
        			}
        			for (int np=0;np<numPars;np++) {
        				if ((correctionParameters[chn][np]==null) || (correctionParameters[chn][np].length!=numberOfLocations)){
        					correctionParameters[chn][np]=new double [numberOfLocations];
        					for (int i=0;i<numberOfLocations;i++) correctionParameters[chn][np][i]=0.0;
        				}
        				for (int i=0;i<numberOfLocations;i++)
        					if (properties.getProperty(prefix+"correctionParameters_"+chn+"_"+np+"_"+i)!=null) {
        						correctionParameters[chn][np][i]=Double.parseDouble(properties.getProperty(prefix+"correctionParameters_"+chn+"_"+np+"_"+i));
        					}
        			}
        		}
        	} else {
        		if (debugLevel>1) System.out.println("numberOfLocations==0, can not restore");
        	}

            fieldStrategies= new FieldStrategies(); // reset old
            fieldStrategies.getProperties(prefix+"fieldStrategies.",properties);
        }        
        
//        public double [] getSampleRadiuses(){ // distance from the current center to each each sample
//            return sampleCorrRadius;
//        }
        
        public double [] getSampleRadiuses(){
            double [] sampleCorrRadius=new double [numberOfLocations];
            //pXY
            for (int i=0;i<numberOfLocations;i++){
                double dx=sampleCoordinates[i][0]-pXY[0];
                double dy=sampleCoordinates[i][1]-pXY[1];
                sampleCorrRadius[i]=getPixelMM()*Math.sqrt(dx*dx+dy*dy);
            }
            return sampleCorrRadius;
        }

        public double getSampleRadius(int sample){
        	double dx=sampleCoordinates[sample][0]-pXY[0];
        	double dy=sampleCoordinates[sample][1]-pXY[1];
        	return getPixelMM()*Math.sqrt(dx*dx+dy*dy);
        }

        
        public void showCurvCorr(String title){
            int width=getSampleWidth();
            int numSamples=getNumSamples();
            String [] chnShortNames={"RS","RT","GS","GT","BS","BT"};
            int numCorrPar=0;
            int maxNumPars=0;
            for (int chn=0;chn<correctionParameters.length;chn++)
                if (correctionParameters[chn]!=null)
                    for (int np=0;np<correctionParameters[chn].length;np++)
                        if (correctionParameters[chn][np]!=null) {
                            numCorrPar++;
                            if (np>maxNumPars) maxNumPars=np;
                        }
            maxNumPars++;
            if (numCorrPar==0){
                String msg="No correction parameters are enabled";
                IJ.showMessage(msg);
                if (debugLevel>1) System.out.println(msg);
                return;
            }
            double [][] pixels=new double [numCorrPar][numSamples];
            String [] titles=new String[numCorrPar];
            int index=0;
            for (int np=0;np<maxNumPars;np++)
                for (int chn=0;chn<correctionParameters.length;chn++)
                    if ((correctionParameters[chn]!=null) && (correctionParameters[chn].length>np) && (correctionParameters[chn][np]!=null)) {
                        titles[index]=chnShortNames[chn]+"-"+curvatureModel[chn].getZName(np);
                        for (int nSample=0;nSample<numSamples;nSample++) {
                            pixels[index][nSample]=correctionParameters[chn][np][nSample];
                        }
                        index++;
                    }
             (new showDoubleFloatArrays()). showArrays(
                     pixels,
                     width,
                     numSamples/width,
                     true, //boolean asStack,
                     title,
                     titles);
            
        }
        public double [] getCalcValuesForZ(double z, double r, double [] corrPars){
        	double [] result=new double [6];
        	for (int chn=0;chn<result.length;chn++) {
        		if (curvatureModel[chn]!=null){
        			result[chn]=curvatureModel[chn].getFdF(
        					corrPars,
        					r, // in mm,
        					Double.NaN, // py,
        					z, //mot_z,
        					null); //deriv_curv[c]);
        		} else {
        			result[chn]=Double.NaN;
        		}
        	}
        	return result;
        }
        
        /**
         * Calculate values (sagital and tangential PSF FWHM in um for each of color channels) for specified z
         * @param z distance from lens (from some zero point), image plane is perpendicular to the axis
         * @param corrected when false - provide averaged (axial model) for radius, if true - with individual correction applied
         * @param allChannels calculate for all (even disabled) channels, false - only for currently selected
         * @return outer dimension - number of channel, inner - number of sample (use getSampleRadiuses for radius of each)
         */
        public double [][] getCalcValuesForZ(double z, boolean corrected, boolean allChannels){
        	double [] sampleCorrRadius=getSampleRadiuses();
            int numSamples=sampleCorrRadius.length;
            double [][] result=new double [6][];
            for (int chn=0;chn<result.length;chn++) {
                if ((curvatureModel[chn]!=null) && (allChannels || channelSelect[chn])){
                    result[chn]=new double [numSamples];
                    for (int sampleIndex=0;sampleIndex<numSamples;sampleIndex++) {
//if ((chn==3) &&  (sampleIndex==23)){
//	System.out.println("getCalcValuesForZ(), chn="+chn+", sampleIndex="+sampleIndex);
//}
                    	
                        result[chn][sampleIndex]=curvatureModel[chn].getFdF(
                                corrected?getCorrPar(chn,sampleIndex):null,
                                        sampleCorrRadius[sampleIndex], //px,
                                        Double.NaN, // py,
                                        z, //mot_z,
                                        null); //deriv_curv[c]);
                    }
                } else {
                    result[chn]=null;
                }
            }
            return result;
        }
        
        public double getChannelBestFWHM(
        		int channel,
        		int sampleIndex, // -1 for center
        		boolean corrected // 
        		){
        	int r0Index=2; // index of "r0" parameter (fwhm is twice that)
        	double [] corrPars=corrected?getCorrPar(channel,sampleIndex):null;
        	double r=(sampleIndex>=0)?getSampleRadius(sampleIndex):0.0;
        	double fwhm=2.0 * curvatureModel[channel].getAr( r, corrPars)[r0Index];
        	return fwhm;
        }

        public double [] getCalcZ(double r,
        		boolean solve // currently if not using z^(>=2) no numeric solution is required - z0 is the minimum 
        		){
            double [] result=new double [6];
            for (int chn=0;chn<result.length;chn++) {
                if (curvatureModel[chn]!=null){
//                    result[chn]=curvatureModel[chn].getAr(r, null)[0];
                    result[chn]=findBestZ(chn, -1, false,solve);
                } else {
                    result[chn]=Double.NaN;
                }
            }
            return result;
        }
        public double findBestZ(
        		int channel,
        		int sampleIndex, // -1 for center
        		boolean corrected,
        		boolean solve // currently if not using z^(>=2) no numeric solution is required - z0 is the minimum 
        		){
        	return findBestZ(
            		channel,
            		sampleIndex, // -1 for center
            		corrected, //
            		solve, // currently if not using z^(>=2) no numeric solution is required - z0 is the minimum 
                    1.0, // double iniStep,
                    0.0001); //double precision)
        }

        public double findBestZ(
        		int channel,
        		int sampleIndex, // -1 for center
        		boolean corrected, // 
        		boolean solve, // currently if not using z^(>=2) no numeric solution is required - z0 is the minimum 
                double iniStep,
                double precision){
            int maxSteps=100;
        	double [] corrPars=corrected?getCorrPar(channel,sampleIndex):null;
        	double r=(sampleIndex>=0)?getSampleRadius(sampleIndex):0.0;
        	double z0=curvatureModel[channel].getAr( r, corrPars)[0];
        	if (!solve) return z0;
        	if (Double.isNaN(z0)) return z0;
        	double f0=getCalcValuesForZ(z0, r,corrPars)[channel];
            double z1=z0+iniStep;
            double f1=getCalcValuesForZ(z1,r,corrPars)[channel];
            double dir = (f1<f0)?1.0:-1.0;
            double z_prev,f_prev;
            if (dir>0) {
                z_prev=z0;
                f_prev=f0;
                z0=z1;
                f0=f1;
            } else {
                z_prev=z1;
                f_prev=f1;
            }
            int step;
            for (step=0;step<maxSteps;step++) {
                z1=z0+dir*iniStep;
                f1=getCalcValuesForZ(z1,r,corrPars)[channel];
                if (f1>f0) break;
                z_prev=z0;
                f_prev=f0;
                z0=z1;
                f0=f1;
            }
            if (step>=maxSteps){
                System.out.println("Failed to find minimum in "+maxSteps+" steps");
                return Double.NaN;
            }
            // now dividing z_prev - z0 - z1 range
            for (step=0;step<maxSteps;step++) {
                if (f_prev>f1){
                    z_prev=z0;
                    f_prev=f0;
                } else {
                    z1=z0;
                    f1=f0;
                }
                z0=(z_prev+z1)/2;
//                f0=getCalcValuesForZ(z1,r)[channel]; // ????
                f0=getCalcValuesForZ(z0,r,corrPars)[channel];
                if (Math.abs(z0-z_prev)<precision) break;
            }
            return z0;
        }
        
        /**
         * calculate distance to "best focus" for each channel (color and S/T) for each sample
         * @param corrected when false - provide averaged (axial model) for radius, if true - with individual correction applied
         * @param allChannels calculate for all (even disabled) channels, false - only for currently selected
         * @return outer dimension - number of channel, inner - number of sample (use getSampleRadiuses for radius of each)
         */
        public double [][] getCalcZ(
        		boolean corrected,
        		boolean allChannels,
        		boolean solve // currently if not using z^(>=2) no numeric solution is required - z0 is the minimum 
        		){
        	double [] sampleCorrRadius=getSampleRadiuses();
            int numSamples=sampleCorrRadius.length;
            boolean [][] goodSamples=new boolean[getNumChannels()][getNumSamples()];
            for (int i=0;i<goodSamples.length;i++) for (int j=0;j<goodSamples[0].length;j++) goodSamples[i][j]=false;
            for (int n=0;n<dataVector.length;n++) if (dataWeights[n]>0.0){
            	goodSamples[dataVector[n].channel][dataVector[n].sampleIndex]=true;
            }
            double [][] result=new double [6][];
            for (int chn=0;chn<result.length;chn++) {
                if ((curvatureModel[chn]!=null) && (allChannels || channelSelect[chn])){
                    result[chn]=new double [numSamples];
                    for (int sampleIndex=0;sampleIndex<numSamples;sampleIndex++) {
                    	if (goodSamples[chn][sampleIndex]) {
/*                    		
                        result[chn][sampleIndex]=curvatureModel[chn].getAr(
                                sampleCorrRadius[sampleIndex],
                                corrected?getCorrPar(chn,sampleIndex):null
                                )[0];
//That was just the initial estimation, true only if z^1... are all 0                        
*/                                

                        result[chn][sampleIndex]=findBestZ(
                        		chn,         // int channel,
                        		sampleIndex, // int sampleIndex,
                        		corrected,   //boolean corrected,
                        		solve);
                    	} else {
                    		result[chn][sampleIndex]=Double.NaN;
                    	}
                    }
                } else {
                    result[chn]=null;
                }
            }
            return result;
        }        

        /**
         * calculate FWHM of the  "best focus" for each channel (color and S/T) for each sample
         * @param corrected when false - provide averaged (axial model) for radius, if true - with individual correction applied
         * @param allChannels calculate for all (even disabled) channels, false - only for currently selected
         * @return outer dimension - number of channel, inner - number of sample (use getSampleRadiuses for radius of each)
         */
        public double [][] getFWHM(
        		boolean corrected,
        		boolean allChannels 
        		){
        	double [] sampleCorrRadius=getSampleRadiuses();
            int numSamples=sampleCorrRadius.length;
            boolean [][] goodSamples=new boolean[getNumChannels()][getNumSamples()];
            for (int i=0;i<goodSamples.length;i++) for (int j=0;j<goodSamples[0].length;j++) goodSamples[i][j]=false;
            for (int n=0;n<dataVector.length;n++) if (dataWeights[n]>0.0){
            	goodSamples[dataVector[n].channel][dataVector[n].sampleIndex]=true;
            }
            double [][] result=new double [6][];
            for (int chn=0;chn<result.length;chn++) {
                if ((curvatureModel[chn]!=null) && (allChannels || channelSelect[chn])){
                    result[chn]=new double [numSamples];
                    for (int sampleIndex=0;sampleIndex<numSamples;sampleIndex++) {
                    	if (goodSamples[chn][sampleIndex]) {
                        result[chn][sampleIndex]=getChannelBestFWHM(
                        		chn,         // int channel,
                        		sampleIndex, // int sampleIndex,
                        		corrected);   //boolean corrected,
                    	} else {
                    		result[chn][sampleIndex]=Double.NaN;
                    	}
                    }
                } else {
                    result[chn]=null;
                }
            }
            return result;
        }        
        
        public double getGreenZCenter(){
            int chn=3; // Green, Tangential
            return curvatureModel[chn].getCenterVector()[0];
        }
        public double [] getZCenters(boolean solve){
//            int chn=3; // Green, Tangential
            double [] result = {
                    curvatureModel[1].getCenterVector()[0], // Red, Tangential
                    curvatureModel[3].getCenterVector()[0], // Green, Tangential
                    curvatureModel[5].getCenterVector()[0]}; // Blue, Tangential
            if (solve) {
            	double [] result_1 = {
            			findBestZ(1, -1, false,true),
            			findBestZ(3, -1, false,true),
            			findBestZ(5, -1, false,true),
            	};
            	return solve?result_1:result;
            }
            return result;
        }
        public double [] getQualB(double z, boolean corrected){
            double [][] data=getCalcValuesForZ(z,corrected, true );
            double [] qualB = {0.0,0.0,0.0};
        	double [] sampleCorrRadius=getSampleRadiuses();
            int numSamples=sampleCorrRadius.length;
            boolean [][] goodSamples=new boolean[getNumChannels()][getNumSamples()];
            for (int i=0;i<goodSamples.length;i++) for (int j=0;j<goodSamples[0].length;j++) goodSamples[i][j]=false;
            for (int n=0;n<dataVector.length;n++) if (dataWeights[n]>0.0){
            	goodSamples[dataVector[n].channel][dataVector[n].sampleIndex]=true;
            }
            for (int c=0;c<3;c++) {
                if ((data[2*c]!=null) && (data[2*c+1]!=null)){
                	int nSamp=0;
                    qualB[c]=0.0;
                    for (int i=0;i<numSamples;i++){
                    	for (int dir=0;dir<2;dir++) {
                        	if (goodSamples[2*c+dir][i]){
                                qualB[c]+=data[2*c+dir][i]*data[2*c+dir][i]*data[2*c+dir][i]*data[2*c+dir][i];
                                nSamp++;
                        	}
                    	}
                    }
                    if (nSamp>0){
                        qualB[c]/=nSamp;
                    }
                    qualB[c]=Math.sqrt(Math.sqrt(qualB[c]));
                } else {
                    qualB[c]=Double.NaN;
                }
            }
            //TODO: Move to a separate function
            int [] numBad={0,0,0,0,0,0};
            boolean hasBad=false;
            for (int i=0;i<goodSamples.length;i++) for (int j=0;j<goodSamples[0].length;j++) if (!goodSamples[i][j]){
            	numBad[i]++;
            	hasBad=true;
            }
            if ((debugLevel>1) && hasBad){ // was 0
            	for (int i=0;i<numBad.length;i++) if (numBad[i]>0){
            		System.out.println(numBad[i]+" sample locations are missing data for "+fieldFitting.getDescription(i));
            	}
            }
            return qualB;
        }

        public double getQualB(
                double z,
                double kr,
                double kb,
                boolean corrected){
            double [] k={kr,1.0,kb};
            double [] qualB=getQualB(z,corrected);
            for (int i=0;i<qualB.length;i++) if (Double.isNaN(qualB[i])){
                qualB[i]=0.0;
                k[i]=0.0;
            }
            return Math.sqrt(Math.sqrt((
                    k[0]*qualB[0]*qualB[0]*qualB[0]*qualB[0]+
                    k[1]*qualB[1]*qualB[1]*qualB[1]*qualB[1]+
                    k[2]*qualB[2]*qualB[2]*qualB[2]*qualB[2])/
                    (k[0]+k[1]+k[2])));
        }

        public double [] getBestQualB(
                double kr,
                double kb,
                boolean corrected){
            return getBestQualB(kr,kb,corrected,1.0,0.001);
            
        }
        
        public double [] getBestQualB( // find minimum
                double kr,
                double kb,
                boolean corrected,
                double iniStep,
                double precision){
            int maxSteps=100;
            double z0=getGreenZCenter();
            double qb0=getQualB(z0,kr,kb,corrected);
            double z1=z0+iniStep;
            double qb1=getQualB(z1,kr,kb,corrected);
            double dir = (qb1<qb0)?1.0:-1.0;
            double z_prev,qb_prev;
            double [] result={Double.NaN,Double.NaN};
            if (dir>0) {
            	z_prev=z0;
            	qb_prev=qb0;
            	z0=z1;
            	qb0=qb1;
            } else {
            	z_prev=z1;
            	qb_prev=qb1;
            }
            int step;
            for (step=0;step<maxSteps;step++) {
            	z1=z0+dir*iniStep;
            	qb1=getQualB(z1,kr,kb,corrected);
            	if (qb1>qb0) break;
            	z_prev=z0;
            	qb_prev=qb0;
            	z0=z1;
            	qb0=qb1;
            }
            if (step>=maxSteps){
            	System.out.println("Failed to find minimum in "+maxSteps+" steps");
            	return result;
            }
            // now dividing z_prev - z0 - z1 range
            for (step=0;step<maxSteps;step++) {
            	if (qb_prev>qb1){
            		z_prev=z0;
            		qb_prev=qb0;
            	} else {
            		z1=z0;
            		qb1=qb0;
            	}
            	z0=(z_prev+z1)/2;
//            	qb0=getQualB(z1,kr,kb,corrected); // ???????????????????
            	qb0=getQualB(z0,kr,kb,corrected);
            	if (Math.abs(z0-z_prev)<precision) break;
            }
            result[0]=z0;
            result[1]=qb0;
            		return result; // z0;
        }
        
        
        public String getDescription(int i){
            return channelDescriptions[i];
        }

        public boolean [] getCenterSelect(){
            return centerSelect;
        }
        public double[] getCenterXY(){
            return pXY;
        }
        public void setCenterXY(double px, double py){
            pXY[0]=px;
            pXY[1]=py;
        }
        public void setDefaultSampleCorr(){
//            int numPars= getNumCurvars()[0]; // number of Z parameters ( [1] - numbnr of radial parameters).
            for (int n=0;n<channelDescriptions.length;n++){
                sampleCorrSelect[n]=getDefaultSampleCorrSelect(); //new boolean[numPars];
                sampleCorrCost[n]=getDefaultSampleCorrCost(); // new double [numPars];
                sampleCorrSigma[n]=getDefaultSampleCorrSigma(); //new double [numPars];
                sampleCorrPullZero[n]=getDefaultCorrPullZero(); //new double [numPars];
//                for (int i=0;i<numPars;i++){
//                    sampleCorrSelect[n][i]=(i<dflt_sampleCorrSelect.length)?dflt_sampleCorrSelect[i]:false;
//                    sampleCorrCost[n][i]=(i<dflt_sampleCorrCost.length)?dflt_sampleCorrCost[i]:1.0;
//                    sampleCorrSigma[n][i]=dflt_sampleCorrSigma;
//                    sampleCorrPullZero[n][i]=dflt_sampleCorrPullZero;
//                }
            }
        }
        public boolean [] getDefaultSampleCorrSelect(){
            boolean [] dflt=new boolean[getNumCurvars()[0]]; // number of Z parameters ( [1] - numbnr of radial parameters).
            for (int i=0;i<dflt.length;i++){
            	dflt[i]=(i<dflt_sampleCorrSelect.length)?dflt_sampleCorrSelect[i]:false;
            }
        	return dflt;
        }

        public double [] getDefaultSampleCorrCost(){
            double [] dflt=new double[getNumCurvars()[0]]; // number of Z parameters ( [1] - numbnr of radial parameters).
            for (int i=0;i<dflt.length;i++){
            	dflt[i]=(i<dflt_sampleCorrCost.length)?dflt_sampleCorrCost[i]:1.0;
            }
        	return dflt;
        }

        public double [] getDefaultSampleCorrSigma(){
            double [] dflt=new double[getNumCurvars()[0]]; // number of Z parameters ( [1] - numbnr of radial parameters).
            for (int i=0;i<dflt.length;i++){
            	dflt[i]=dflt_sampleCorrSigma;
            }
        	return dflt;
        }

        public double [] getDefaultCorrPullZero(){
            double [] dflt=new double[getNumCurvars()[0]]; // number of Z parameters ( [1] - numbnr of radial parameters).
            for (int i=0;i<dflt.length;i++){
            	dflt[i]=dflt_sampleCorrPullZero;
            }
        	return dflt;
        }
/**
         * Dialog to setup per-sample (coordinate) corrections
         * @param title Dialog title
         * @param individualChannels configure each of the six color/dir channels separately (false - apply to all)
         * @param disabledPars enable use of the parameters that are currently disabled from fitting
         * @return true if OK was pressed, false - if cancel
         */
        public boolean setupSampleCorr(String title, boolean individualChannels, boolean disabledPars){
        	int firstChn=0;
        	for (int i=0;i<channelSelect.length;i++) if (channelSelect[i]){
        		firstChn=i;
        		break;
        	}
        	if (!channelSelect[firstChn]){
        		String msg="No channels selected, please select at least one";
        		IJ.showMessage(msg);
        		System.out.println(msg);
        	}
        	int fromChn=individualChannels?0:firstChn;
        	int toChn=individualChannels?channelSelect.length:(firstChn+1);
        	boolean resetCorrections=false;
        	GenericDialog gd = new GenericDialog(title);
        	gd.addCheckbox("Reset all per-sample corrections to zero", resetCorrections);

        	int numParsZR[] = getNumCurvars(); // [0] - Z, [1] - r,
        	for (int nChn=fromChn;nChn<toChn;nChn++) if (channelSelect[nChn]){

        		String chnName=individualChannels?channelDescriptions[nChn]:"All selected Channels";
        		for (int i=0;i<sampleCorrSelect[nChn].length;i++) if (disabledPars || curvatureSelect[nChn][i*numParsZR[1]]){
        			gd.addMessage("===== "+chnName+", "+curvatureModel[nChn].getZDescription(i)+" =====");
        			gd.addCheckbox("Enable per-sample correction for \""+curvatureModel[nChn].getZDescription(i)+"\"", sampleCorrSelect[nChn][i]);
        			gd.addNumericField("Correction cost",sampleCorrCost[nChn][i],5,8,"um");
        			gd.addNumericField("Correction sigma",sampleCorrSigma[nChn][i],5,8,"mm");
        			gd.addNumericField("Pull to zero fraction",100*sampleCorrPullZero[nChn][i],4,8,"%");
        		}
        	}
        	gd.enableYesNoCancel("Apply","Keep"); // default OK (on enter) - "Apply"
        	WindowTools.addScrollBars(gd);
        	gd.showDialog();
        	if (gd.wasCanceled()) return false;
        	if (gd.wasOKed()) { // selected non-default "Apply"
        		resetCorrections= gd.getNextBoolean();
        		for (int nChn=fromChn;nChn<toChn;nChn++) if (channelSelect[nChn]){
        			for (int i=0;i<sampleCorrSelect[nChn].length;i++) if (disabledPars || curvatureSelect[nChn][i*numParsZR[1]]){
        				sampleCorrSelect[nChn][i]= gd.getNextBoolean();
        				sampleCorrCost[nChn][i]= gd.getNextNumber();
        				sampleCorrSigma[nChn][i]= gd.getNextNumber();
        				sampleCorrPullZero[nChn][i]= 0.01*gd.getNextNumber();
        			} else {
        				sampleCorrSelect[nChn][i]= false;
        			}
        		}
        		if (!individualChannels){ // copy settings to other channels
        			for (int nChn=0;nChn<channelSelect.length;nChn++) if (channelSelect[nChn] && (nChn!=fromChn) ){
        				for (int i=0;i<sampleCorrSelect[nChn].length;i++){
        					sampleCorrSelect[nChn][i]= sampleCorrSelect[fromChn][i];
        					sampleCorrCost[nChn][i]= sampleCorrCost[fromChn][i];
        					sampleCorrSigma[nChn][i]= sampleCorrSigma[fromChn][i];
        					sampleCorrPullZero[nChn][i]= sampleCorrPullZero[fromChn][i];
        				}
        			}
        		}
        	}
        	if (resetCorrections) resetSampleCorr();
        	return true;
        }
        // once per data set
        public void resetSampleCorr(){
        	if (debugLevel>0) System.out.println("---resetSampleCorr()---");
            for (int i=0; i<correctionParameters.length;i++) correctionParameters[i]=null;
        }
        
        public double [] getCorrVector(){
            //numberOfLocations
            int numPars=0;
            for (int nChn=0; nChn< sampleCorrChnParIndex.length;nChn++) {
                if (sampleCorrChnParIndex[nChn]!=null){
                    for (int nPar=0;nPar< sampleCorrChnParIndex[nChn].length;nPar++) {
                        if (sampleCorrChnParIndex[nChn][nPar]>=0){
                            numPars+=numberOfLocations;
                        }
                    }
                }
            }
            if (debugLevel>1) System.out.println("getCorrVector 1");
            sampleCorrVector=new double [numPars];
            for (int nChn=0; nChn< sampleCorrChnParIndex.length;nChn++) {
                if (sampleCorrChnParIndex[nChn]!=null){
                    for (int nPar=0;nPar< sampleCorrChnParIndex[nChn].length;nPar++) {
                        if (sampleCorrChnParIndex[nChn][nPar]>=0){
                            for (int i=0;i<numberOfLocations;i++){
                                if ((correctionParameters[nChn]!=null) && (correctionParameters[nChn][nPar]!=null) && (correctionParameters[nChn][nPar].length>i))
                                    sampleCorrVector[sampleCorrChnParIndex[nChn][nPar]+i]=correctionParameters[nChn][nPar][i];
                                else {
                                    sampleCorrVector[sampleCorrChnParIndex[nChn][nPar]+i]=0.0;
                                    if ((correctionParameters[nChn]!=null) && (correctionParameters[nChn][nPar]!=null) && (correctionParameters[nChn][nPar].length<i)){
                                    	if (debugLevel>1) System.out.println("correctionParameters["+nChn+"]["+nPar+"].length < "+i);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return sampleCorrVector;
        }
        
        public String [] getCorrNames(){
            //numberOfLocations
            int numPars=0;
            for (int nChn=0; nChn< sampleCorrChnParIndex.length;nChn++) {
                if (sampleCorrChnParIndex[nChn]!=null){
                    for (int nPar=0;nPar< sampleCorrChnParIndex[nChn].length;nPar++) {
                        if (sampleCorrChnParIndex[nChn][nPar]>=0){
                            numPars+=numberOfLocations;
                        }
                    }
                }
            }
            if (debugLevel>1) System.out.println("getCorrVector 1");
            String [] corrNames=new String [numPars];
            for (int nChn=0; nChn< sampleCorrChnParIndex.length;nChn++) {
                if (sampleCorrChnParIndex[nChn]!=null){
                    for (int nPar=0;nPar< sampleCorrChnParIndex[nChn].length;nPar++) {
                        if (sampleCorrChnParIndex[nChn][nPar]>=0){
                            for (int i=0;i<numberOfLocations;i++){
                                corrNames[sampleCorrChnParIndex[nChn][nPar]+i]="chn"+nChn+":"+nPar+"-"+i;
                            }
                        }
                    }
                }
            }
            return corrNames;
        }
        
        
        public void commitCorrVector(){
        	if (debugLevel>1) System.out.println("commitCorrVector()");
            commitCorrVector(sampleCorrVector);
        }
        
        public void commitCorrVector(double [] vector){
            for (int nChn=0; nChn< sampleCorrChnParIndex.length;nChn++) {
                if (sampleCorrChnParIndex[nChn]!=null){
                    for (int nPar=0;nPar< sampleCorrChnParIndex[nChn].length;nPar++) {
                        if (sampleCorrChnParIndex[nChn][nPar]>=0){
                            if(correctionParameters[nChn]==null)
                                correctionParameters[nChn]=new double [sampleCorrChnParIndex[nChn].length][];
                            if ((correctionParameters[nChn][nPar]==null) || (correctionParameters[nChn][nPar].length!=numberOfLocations)){
                            	if (debugLevel>1) System.out.println("commitCorrVector(): correctionParameters["+nChn+"]["+nPar+"].length < "+numberOfLocations);

                                correctionParameters[nChn][nPar]=new double [numberOfLocations];
                            }
                            for (int i=0;i<numberOfLocations;i++){
                                correctionParameters[nChn][nPar][i]=vector[sampleCorrChnParIndex[nChn][nPar]+i];
                            }
                        }
                    }
                }
            }
        }
        
        public void setEstimatedZ0(
        		double [] z0,
        		boolean force){
        	if ((z0==null)|| (curvatureModel==null)) return; // no estimation available
        	for (int chn=0;chn<curvatureModel.length;chn++){
        		if (!Double.isNaN(z0[chn]) && (!curvatureModel[chn].z0IsValid() || force)){
        			curvatureModel[chn].set_z0(z0[chn]);
        			if (debugLevel>1) System.out.println("Setting initial (estimated) best focal position for channel "+chn+" = "+z0[chn]);
        		}
        	}
        }
        
        /**
         * create matrix of weights of the other parameters influence
         * @param sampleCoordinates [sample number]{x,y} - flattened array of sample coordinates
         * Run in the beginning of fitting series (zeroes the values)
         */
        // once per fitting series (or parameter change
        public void initSampleCorrVector(
        		double [][] sampleCoordinates,
        		double [][] sampleSeriesWeights){
        	if (debugLevel>1) System.out.println("initSampleCorrVector()");
            numberOfLocations=sampleCoordinates.length;
            this.sampleCoordinates=new double[sampleCoordinates.length][];
            for (int i=0;i<sampleCoordinates.length;i++) this.sampleCoordinates[i]=sampleCoordinates[i].clone();
//            int numSamples=sampleCoordinates.length;
            for (int nChn=0; nChn< sampleCorrSelect.length;nChn++) {
                if (channelSelect[nChn]){
                    int numPars=sampleCorrSelect[nChn].length;
                    sampleCorrCrossWeights[nChn]=new double[numPars][][];
                    for (int nPar=0;nPar<numPars;nPar++){
                        if (sampleCorrSelect[nChn][nPar]){
                            sampleCorrCrossWeights[nChn][nPar]=new double[numberOfLocations][numberOfLocations];
                            double k=-getPixelMM()*getPixelMM()/(sampleCorrSigma[nChn][nPar]*sampleCorrSigma[nChn][nPar]);
                            for (int i=0;i<numberOfLocations;i++){
                                double sw=0.0;
                                for (int j=0;j<numberOfLocations;j++) {
                                    if (i!=j){
                                        double dx=sampleCoordinates[i][0]-sampleCoordinates[i][0];
                                        double dy=sampleCoordinates[i][1]-sampleCoordinates[i][1];
                                        double a= Math.exp(k*(dx*dx+dy*dy));
                                        if (a<0.01) a=0.0; // reduce numbner of non-zero matrix elements
                                        sampleCorrCrossWeights[nChn][nPar][i][j]=a;
                                        sw+=a;
                                    }
                                }
                                double normalizedCost=sampleCorrCost[nChn][nPar];
                                if ((sampleSeriesWeights!=null) && ((sampleSeriesWeights[nChn][i]!=0.0))) normalizedCost*=sampleSeriesWeights[nChn][i];
                                if ((sampleCorrPullZero[nChn][nPar]==0) ||(sampleCorrCost[nChn][nPar]==0)) sw=0.0;
                                else if (sw!=0.0) sw=-normalizedCost*sampleCorrPullZero[nChn][nPar]/sw;
                                for (int j=0;j<numberOfLocations;j++) {
                                    if (i!=j){
                                        sampleCorrCrossWeights[nChn][nPar][i][j]*=sw;
                                    } else {
                                        sampleCorrCrossWeights[nChn][nPar][i][j]=normalizedCost;
                                    }
                                }
                            }
                        } else {
                            sampleCorrCrossWeights[nChn][nPar]=null;
                        }
                    }
                } else {
                    sampleCorrCrossWeights[nChn]=null;
                }
            }
        	getCorrVector();
        }
        
 
        public void initSampleCorrChnParIndex(
        		double [][] sampleCoordinates){
        	numberOfLocations=sampleCoordinates.length;
        	this.sampleCoordinates=new double[sampleCoordinates.length][];
        	for (int i=0;i<sampleCoordinates.length;i++) this.sampleCoordinates[i]=sampleCoordinates[i].clone();
        	sampleCorrChnParIndex=new int [sampleCorrSelect.length][];
        	int numPars=0;
        	for (int nChn=0; nChn< sampleCorrSelect.length;nChn++) {
        		if (channelSelect[nChn]) {
        			sampleCorrChnParIndex[nChn]=new int [sampleCorrSelect[nChn].length];
        			for (int nPar=0;nPar< sampleCorrChnParIndex[nChn].length;nPar++) {
        				if (sampleCorrSelect[nChn][nPar]) {
        					sampleCorrChnParIndex[nChn][nPar]=numPars; // pointer to the first sample
        					numPars+=numberOfLocations;
        				} else {
        					sampleCorrChnParIndex[nChn][nPar]=-1;
        				}
        			}                
        		} else {
        			sampleCorrChnParIndex[nChn]=null;
        		}
        	}
        	if (debugLevel>1) System.out.println("initSampleCorrChnParIndex()");
        	// currently all correction parameters are initialized as zeros.
        	getCorrVector();
        }
        
        
        
        
        public double [][] getSampleCoordinates(){
        	return sampleCoordinates;
        }
        public double [] getCorrPar(int chn, int sampleIndex){
            if (correctionParameters[chn]==null) return null;
            double [] corr =new double [correctionParameters[chn].length];
            for (int i=0;i<corr.length;i++){
            	if ((correctionParameters[chn][i] !=null) && (correctionParameters[chn][i].length<=i)){
            		if (debugLevel>1) System.out.println("getCorrPar(): correctionParameters["+chn+"]["+i+"].length="+correctionParameters[chn][i].length);
            	}
                if ((correctionParameters[chn][i] !=null) && (correctionParameters[chn][i].length>i)) corr[i]=correctionParameters[chn][i][sampleIndex];
                else corr[i]=0.0;
            }

            return corr;
        }
        
        public double [][] getCorrPar(int sampleIndex){
            if (sampleCorrChnParIndex==null) return null;
            double [][] result=new double [sampleCorrChnParIndex.length][];
            boolean non_null=false;
            for (int i=0;i<sampleCorrChnParIndex.length;i++){
                result[i]=getCorrPar(i, sampleIndex);
                non_null |= (result[i]!=null);
            }
            return non_null?result:null;
        }
        public boolean[] getDefaultMask(){
            boolean [] mask = {true,true,true,true,true,true};
            return mask;
        }
        public FieldFitting(){} // just to get descriptions

        public FieldFitting(
                double pX0,
                double pY0,
                int distanceParametersNumber,
                int radialParametersNumber)
        {
            fieldStrategies= new FieldStrategies();
            pXY=new double [2];
            pXY[0]=pX0;
            pXY[1]=pY0;
            channelSelect=getDefaultMask();
            mechanicalFocusingModel=new MechanicalFocusingModel();
            mechanicalSelect=mechanicalFocusingModel.getDefaultMask();
            centerSelect=centerSelectDefault.clone();
            for (int i=0;i<curvatureModel.length;i++){
                curvatureModel[i]= new CurvatureModel(
                        pX0,
                        pY0,
                        distanceParametersNumber,
                        radialParametersNumber);
                curvatureSelect[i]=curvatureModel[i].getDefaultMask();
            }
            setDefaultSampleCorr(); // should be after curvatureModel
        }
        
        public boolean [] getSelectedChannels(){
            return this.channelSelect;
        }
        public int [] getNumCurvars(){
            try {
                return curvatureModel[0].getNumPars();
            }catch (Exception e) {
                int [] dflt_numPars={CurvatureModel.dflt_distanceParametersNumber,CurvatureModel.dflt_radialParametersNumber};
                return dflt_numPars;

            }
        }
        public boolean maskSetDialog(
        		String title//,
//         		String strategyComment,
//        		double lambda,
//        		boolean lastInSeries
        		){
        	GenericDialog gd = new GenericDialog(title);
        	boolean editMechMask=false;
        	boolean editCurvMask=false;
        	boolean commonCurvMask=true;
        	boolean detailedCurvMask=false;
        	boolean setupCorrectionPars=false;
        	boolean commonCorrectionPars=true;
        	boolean disabledCorrectionPars=false;
        	if (centerSelect==null) centerSelect=centerSelectDefault.clone();
        	gd.addCheckbox("Adjust aberration center (pX0)", centerSelect[0]);
        	gd.addCheckbox("Adjust aberration center (pY0)", centerSelect[1]);

        	if (channelSelect==null) channelSelect=getDefaultMask();
        	for (int i=0;i<channelSelect.length;i++) {
        		gd.addCheckbox(getDescription(i), channelSelect[i]);                    
        	}
        	gd.addCheckbox("Edit mechanical parameters masks", editMechMask);
        	gd.addCheckbox("Edit curvature model parameters mask(s)", editCurvMask);
        	gd.addCheckbox("Apply same curvature model parameters mask to all channels", commonCurvMask);
        	gd.addCheckbox("Edit full matrix of the curvature model parameters masks", detailedCurvMask);
        	gd.addMessage("");
        	gd.addCheckbox("Setup per-sample correction", setupCorrectionPars);
        	gd.addCheckbox("Apply same per-sample corrections to all channels", commonCorrectionPars);
        	gd.addCheckbox("Setup correction parameters when the parameter itself is disabled", disabledCorrectionPars);
        	gd.addMessage("---");
        	
     		gd.addStringField("Strategy comment",strategyComment,60);
     		gd.addNumericField("Initial LMA lambda",lambda,3,5,"");
        	gd.addCheckbox("Reset optical center to distortions center", resetCenter);
        	gd.addCheckbox("Reset correction parameters before this LMA step", !keepCorrectionParameters);
        	gd.addCheckbox("Stop after this LMA step", lastInSeries);

        	//         gd.enableYesNoCancel("Keep","Apply"); // default OK (on enter) - "Keep"
        	gd.showDialog();
        	if (gd.wasCanceled()) return false;
        	centerSelect[0]=gd.getNextBoolean();
        	centerSelect[1]=gd.getNextBoolean();
        	for (int i=0;i<channelSelect.length;i++) {
        		channelSelect[i]=gd.getNextBoolean();
        	}
        	editMechMask=gd.getNextBoolean();
        	editCurvMask=gd.getNextBoolean();
        	commonCurvMask=gd.getNextBoolean();
        	detailedCurvMask=gd.getNextBoolean();
        	setupCorrectionPars=gd.getNextBoolean();
        	commonCorrectionPars=gd.getNextBoolean();
        	disabledCorrectionPars=gd.getNextBoolean();

        	strategyComment=gd.getNextString();
     		lambda=gd.getNextNumber();
     		resetCenter=gd.getNextBoolean();
        	keepCorrectionParameters=!gd.getNextBoolean();
        	lastInSeries=gd.getNextBoolean();
        	//         boolean OK;
        	if (editMechMask){
        		boolean [] mask=mechanicalFocusingModel.maskSetDialog("Focusing mechanical parameters mask", mechanicalSelect);
        		if (mask!=null) mechanicalSelect=mask;
        		else return false; // canceled
        	}
        	if (editCurvMask) {
        		if (commonCurvMask){
        			boolean [] mask=new boolean [curvatureSelect[0].length];
        			for (int i=0;i<mask.length;i++){
        				mask[i]=curvatureSelect[0][i];
        				for (int j=1;j<curvatureSelect.length;j++){
        					mask[i] |= curvatureSelect[j][i];
        				}
        			}
        			mask=curvatureModel[0].maskSetDialog(
        					"All channels mask for all curvature models (colors,S/T)",
        					detailedCurvMask,
        					mask);
        			if (mask==null) return false; // canceled
        			for (int i=0;i<curvatureSelect.length;i++){
        				for (int j=0;j<mask.length;j++){
        					curvatureSelect[i][j] = mask[j];
        				}
        			}
        		} else {
        			for (int i=0;i<channelSelect.length;i++) if (channelSelect[i]){
        				boolean [] mask=curvatureSelect[i];
        				mask=curvatureModel[0].maskSetDialog(
//        						"Parameter mask for curvature model, channel \""+getDescription(i)+"\"",
        						"Channel \""+getDescription(i)+"\" parameter mask for curvature model",
        						detailedCurvMask,
        						mask);
        				if (mask==null) return false; // canceled
        				curvatureSelect[i]=mask;
        			}
        		}
        	}
        	if (setupCorrectionPars) {
        		if (!setupSampleCorr("Setup per-sample correction parameters",
        				!commonCorrectionPars,
        				disabledCorrectionPars)) return false;
        		//             initSampleCorr(flattenSampleCoord());
        	}
        	// will modify
        	initSampleCorrChnParIndex(flattenSampleCoord()); // run always regardless of configured or not (to create zero-length array of corr)
        	return true;
        }
        
        public void selectZTilt(boolean allChannels){
        	mechanicalSelect=mechanicalFocusingModel.maskSetZTxTy(); // enable z0, tx, ty
        	// enable all color/dir channels (add separate selection dialog?)
        	for (int i=0;i<channelSelect.length;i++) {
        		if (allChannels) channelSelect[i]=true;
				curvatureSelect[i]=curvatureModel[0].maskAllDisabled();
        	}
        	if (sampleCorrSelect!=null){
            	for (int i=0;i<sampleCorrSelect.length;i++) if (sampleCorrSelect[i]!=null) {
                	for (int j=0;j<sampleCorrSelect[i].length;j++) sampleCorrSelect[i][j]=false;
            	}
        	}
        	initSampleCorrChnParIndex(flattenSampleCoord());
        }
        
        ArrayList<String> getParameterValueStrings(boolean showDisabled, boolean showCorrection){
            ArrayList<String> parList=new ArrayList<String>();
            parList.add("\t ===== Aberrations center =====\t\t");
            String [] centerDescriptions={"Aberrations center X","Aberrations center Y"};
            for (int i=0;i<2;i++){
                if (centerSelect[i] ) {
                    parList.add("\t"+centerDescriptions[i]+":\t"+pXY[i]+"\tpix");
                } else if (showDisabled){
                    parList.add("(disabled)\t"+centerDescriptions[i]+":\t"+pXY[i]+"\tpix");
                }
            }
            parList.add("\t ===== Mechanical model parameters =====\t\t");
            for (int i=0;i<mechanicalFocusingModel.paramValues.length;i++){
                if ((mechanicalSelect==null) || mechanicalSelect[i] ) {
                    parList.add("\t"+mechanicalFocusingModel.getDescription(i)+":\t"+mechanicalFocusingModel.paramValues[i]+
                            "\t"+mechanicalFocusingModel.getUnits(i));
                } else if (showDisabled){
                    parList.add("(disabled)\t"+mechanicalFocusingModel.getDescription(i)+":\t"+mechanicalFocusingModel.paramValues[i]+
                            "\t"+mechanicalFocusingModel.getUnits(i));
                }
            }
            for (int n=0;n<channelSelect.length;n++) if (channelSelect[n]){
                boolean isMasterDir=(n%2) == (sagittalMaster?0:1);
                parList.add("\t ===== Curvature model parameters for \""+ getDescription(n)+"\"=====\t\t");
                int n1= n ^ 1;
             int index=0;
                for (int i=0;i<curvatureModel[n].modelParams.length;i++) for (int j=0;j<curvatureModel[n].modelParams[0].length;j++){
                    String name=curvatureModel[n].getZDescription(i)+", "+curvatureModel[n].getRadialDecription(j);
                    if ((j==0) && !isMasterDir){ // dependent, copy center coefficients from the master
                        if ((curvatureSelect[n1]==null) || curvatureSelect[n1][index] ) {
                            parList.add("(master)\t"+name+":\t"+curvatureModel[n1].modelParams[i][j]+"\t");
                            // debugging
                            parList.add("(this)\t"+name+":\t"+curvatureModel[n].modelParams[i][j]+"\t");
                        } else if (showDisabled) {
                            parList.add("(master_disabled)\t"+name+":\t"+curvatureModel[n1].modelParams[i][j]+"\t");
                            // debugging
                            parList.add("(this_disabled)\t"+name+":\t"+curvatureModel[n1].modelParams[i][j]+"\t");
                        }
                    } else {
                        if ((curvatureSelect[n]==null) || curvatureSelect[n][index] ) {
                            parList.add("\t"+name+":\t"+curvatureModel[n].modelParams[i][j]+"\t");
                        } else if (showDisabled) {
                            parList.add("(disabled)\t"+name+":\t"+curvatureModel[n].modelParams[i][j]+"\t");
                        }
                    }
                    index++;
                }
            }
//            if (showCorrection && (getNumberOfCorrParameters()>0)){
            if (showCorrection){
                parList.add("\t ===== Per-sample correction parameters =====\t\t");
                for (int n=0;n<correctionParameters.length;n++) if (correctionParameters[n]!=null){
                	for (int np=0;np<correctionParameters[n].length;np++)
                		if ((correctionParameters[n][np]!=null) && (correctionParameters[n][np].length==numberOfLocations)){
                			//                        int numSamples=sampleCorrCrossWeights[n][np].length;
                			parList.add("\t ----- correction parameters for \""+ getDescription(n)+" "+curvatureModel[n].getZDescription(np)+"\" -----\t\t");
                			for (int i=0;i<numberOfLocations;i++){
                				parList.add(i+"\t"+curvatureModel[n].getZDescription(np)+":\t"+correctionParameters[n][np][i]+"\t");
                			}
                		} else {
                			if ((correctionParameters[n][np]!=null) && (correctionParameters[n][np].length!=numberOfLocations)){
                				if (debugLevel>1) System.out.println("getParameterValueStrings(): correctionParameters["+n+"]["+np+"].length="+correctionParameters[n][np].length);
                			}
                		}
                }
            }
            return parList;
        }
        // Show/modify all parameters in a single window.
        public boolean showModifyParameterValues(String title, boolean showDisabled){
         GenericDialog gd = new GenericDialog(title);
         gd.addMessage("===== Aberrations center =====");
            String [] centerDescriptions={"Aberrations center X","Aberrations center Y"};
            for (int i=0;i<2;i++){
            	double distXY=(i==0)?pX0_distortions:pY0_distortions;
                if (centerSelect[i] ) {
                    gd.addNumericField(centerDescriptions[i],pXY[i],5,10,"pix ("+IJ.d2s(distXY,1)+")");
                } else if (showDisabled){
                    gd.addNumericField("(disabled) "+centerDescriptions[i],pXY[i],5,10,"pix ("+IJ.d2s(distXY,1)+")");
                }
            }
         
         gd.addMessage("===== Mechanical model parameters =====");
            for (int i=0;i<mechanicalFocusingModel.paramValues.length;i++){
                if ((mechanicalSelect==null) || mechanicalSelect[i] ) {
                    gd.addNumericField(mechanicalFocusingModel.getDescription(i),mechanicalFocusingModel.paramValues[i],5,8,
                            mechanicalFocusingModel.getUnits(i));
                } else if (showDisabled){
                    gd.addNumericField("(disabled) "+mechanicalFocusingModel.getDescription(i),mechanicalFocusingModel.paramValues[i],5,8,
                            mechanicalFocusingModel.getUnits(i));
                }
            }
            for (int n=0;n<channelSelect.length;n++) if (channelSelect[n]){
                boolean isMasterDir=(n%2) == (sagittalMaster?0:1);
                int n1= n ^ 1;

             gd.addMessage("===== Curvature model parameters for \""+ getDescription(n)+"\"=====");
             int index=0;
                for (int i=0;i<curvatureModel[n].modelParams.length;i++) for (int j=0;j<curvatureModel[n].modelParams[0].length;j++){
                    String name=curvatureModel[n].getZDescription(i)+", "+curvatureModel[n].getRadialDecription(j);
                    if ((j==0) && !isMasterDir){ // dependent, copy center coefficients from the master
                        if ((curvatureSelect[n1]==null) || curvatureSelect[n1][index] ) {
                            gd.addNumericField("(copied from master) "+name,curvatureModel[n1].modelParams[i][j],5,8,"");
                        } else if (showDisabled) {
                            gd.addNumericField("(copied from disabled master) "+name,curvatureModel[n1].modelParams[i][j],5,8,"");
                        }
                    } else {
                        if ((curvatureSelect[n]==null) || curvatureSelect[n][index] ) {
                            gd.addNumericField(name,curvatureModel[n].modelParams[i][j],5,8,"");
                        } else if (showDisabled) {
                            gd.addNumericField("(disabled) "+name,curvatureModel[n].modelParams[i][j],5,8,"");
                        }
                    }
                    index++;
                }
            }
         
         gd.enableYesNoCancel("Apply","Keep"); // default OK (on enter) - "Apply"
     WindowTools.addScrollBars(gd);
    gd.showDialog();
         if (gd.wasCanceled()) return false;
         if (gd.wasOKed()) { // selected default "Apply"
                for (int i=0;i<2;i++){
                    if (centerSelect[i] || showDisabled) {
                        pXY[i]=gd.getNextNumber();
                    }
                }
                for (int i=0;i<mechanicalFocusingModel.paramValues.length;i++){
                    if ((mechanicalSelect==null) || mechanicalSelect[i] || showDisabled) {
                        mechanicalFocusingModel.paramValues[i]=gd.getNextNumber();
                    }
                }
                for (int n=0;n<channelSelect.length;n++) if (channelSelect[n]){
                    int index=0;
                    for (int i=0;i<curvatureModel[n].modelParams.length;i++) for (int j=0;j<curvatureModel[n].modelParams[0].length;j++){
                        if ((curvatureSelect[n]==null) || curvatureSelect[n][index] || showDisabled) {
                            curvatureModel[n].modelParams[i][j]=gd.getNextNumber();
                        }
                    }
                    index++;
                }
         }         
         return true;
        }

        public int getNumberOfCorrParameters(){ // selected for fitting
            return (sampleCorrVector!=null)?sampleCorrVector.length:0;
        }
        
        public int getNumberOfParameters(boolean sagittalMaster){
            return getNumberOfRegularParameters(sagittalMaster)+getNumberOfCorrParameters();
        }        
        /**
         * @return number of selected parameters (including center, mechanical and each selected - up to 6 - curvature)
         */
        public int getNumberOfRegularParameters(boolean sagittalMaster){
            int np=0;
            for (int i=0;i<2;i++){
                if ( centerSelect[i]) np++;
            }
            for (int i=0;i<mechanicalFocusingModel.paramValues.length;i++){
                if ((mechanicalSelect==null) || mechanicalSelect[i] ) np++;
            }
            for (int n=0;n<channelSelect.length;n++) if (channelSelect[n]){
                boolean isMasterDir=(n%2) == (sagittalMaster?0:1);
             int index=0;
                for (int i=0;i<curvatureModel[n].modelParams.length;i++) for (int j=0;j<curvatureModel[n].modelParams[0].length;j++){
                    if ((isMasterDir || (j!=0)) && ((curvatureSelect[n]==null) || curvatureSelect[n][index] )) np++;
                    index++;
                }
            }
            return np;
        }
        
        /**
         * @return number of selected channels (up to 6 - colors and S/T)
         */
        public int getNumberOfChannels(){
            int nc=0;    
            for (int n=0;n<channelSelect.length;n++) if (channelSelect[n]) nc++;
            return nc;
        }
        
        /**
         * @return vector of the current selected parameters values
         */
        public double [] createParameterVector(boolean sagittalMaster){
        	int debugThreshold=0;
            double [] pars = new double [getNumberOfParameters(sagittalMaster)];
            int np=0;
            if (debugLevel>debugThreshold) debugParameterNames=new String [pars.length];
            for (int i=0;i<2;i++){
                if ( centerSelect[i]) {
                    if (debugLevel>debugThreshold) debugParameterNames[np]="pXY"+i;
                	pars[np++]=pXY[i];
                }
            }

            for (int i=0;i<mechanicalFocusingModel.paramValues.length;i++){
                if ((mechanicalSelect==null) || mechanicalSelect[i] ) {
                    if (debugLevel>debugThreshold) debugParameterNames[np]=mechanicalFocusingModel.getName(i);
                	pars[np++]=mechanicalFocusingModel.paramValues[i];
                }
            }
            for (int n=0;n<channelSelect.length;n++) if (channelSelect[n]){
                boolean isMasterDir=(n%2) == (sagittalMaster?0:1);
             int index=0;
                for (int i=0;i<curvatureModel[n].modelParams.length;i++) for (int j=0;j<curvatureModel[n].modelParams[0].length;j++){
                    if ((isMasterDir || (j!=0)) && ((curvatureSelect[n]==null) || curvatureSelect[n][index] )) {
                        if (debugLevel>debugThreshold) debugParameterNames[np]="chn"+n+"-"+curvatureModel[n].getZName(i)+":"+curvatureModel[n].getRadialName(j);
                        pars[np++]=curvatureModel[n].modelParams[i][j];
                    }
                    index++;
                }
            }
            if (debugLevel>1) System.out.println("createParameterVector(): using sampleCorrVector - do we need to create it first?");
            getCorrVector(); // do we need that?            
            int nCorrPars=getNumberOfCorrParameters();
            String [] corrNames=(debugLevel>debugThreshold)?getCorrNames():null;
            for (int i=0;i<nCorrPars;i++) {
                if (debugLevel>debugThreshold) debugParameterNames[np]="corr_par-"+corrNames[i];
            	pars[np++]=sampleCorrVector[i];
            }
            if (debugLevel>1){
            	for (int i=0;i<pars.length;i++){
            		System.out.println(i+" "+debugParameterNames[i]+" = "+pars[i]);
            	}
            }
            return pars;
        }
        
        /**
         * Apply (modified) parameter values to selected ones
         * @param pars vector corresponding to selected parameters
         */
        public void commitParameterVector(double [] pars, boolean sagittalMaster){
            int np=0;
            for (int i=0;i<2;i++){
                if ( centerSelect[i]) pXY[i]=pars[np++];
            }
            for (int i=0;i<mechanicalFocusingModel.paramValues.length;i++){
                if ((mechanicalSelect==null) || mechanicalSelect[i] ) mechanicalFocusingModel.paramValues[i] = pars[np++];;
            }
            for (int n=0;n<channelSelect.length;n++) if (channelSelect[n]){
                boolean isMasterDir=(n%2) == (sagittalMaster?0:1);
             int index=0;
                for (int i=0;i<curvatureModel[n].modelParams.length;i++) for (int j=0;j<curvatureModel[n].modelParams[0].length;j++){
                    if ((isMasterDir || (j!=0)) && ((curvatureSelect[n]==null) || curvatureSelect[n][index] )) {
                        curvatureModel[n].modelParams[i][j]=pars[np++];
                    }
                    index++;
                }
            }
            // copy correction parameters
            if (debugLevel>1) System.out.println("commitParameterVector():  Creating and committing sampleCorrVector");
            int nCorrPars=getNumberOfCorrParameters();
            for (int i=0;i<nCorrPars;i++) sampleCorrVector[i]=pars[np++];
            commitCorrVector();

// copy center parameters to dependent
// copy if master is selected, regardless of is dependent selected or not            
            for (int n=0;n<channelSelect.length;n++) if (channelSelect[n]){
                boolean isMasterDir = (n%2) == (sagittalMaster?0:1);
                if (isMasterDir){
                    curvatureModel[n ^ 1].setCenterVector(curvatureModel[n].getCenterVector());
                }
            }
// propagate pXY to each channel (even disabled)            
            for (int n=0;n<curvatureModel.length;n++){
                curvatureModel[n].setCenter(pXY[0],pXY[1]);
            }
        }
        public double getMotorsZ(
                int [] motors, // 3 motor coordinates
                double px, // pixel x
                double py) {// pixel y
            return mechanicalFocusingModel.calc_ZdZ(
                    motors,
                    px,
                    py,
                    null);

        }
        public double getRadiusMM(
                double px, // pixel x
                double py) {// pixel y
            double pd=Math.sqrt((px-currentPX0)*(px-currentPX0)+(py-currentPY0)*(py-currentPY0));
            return pd*getPixelMM();
            
        }

        public double getRadiusMM_distortions(
                double px, // pixel x
                double py) {// pixel y
            double pd=Math.sqrt((px-pX0_distortions)*(px-pX0_distortions)+(py-pY0_distortions)*(py-pY0_distortions));
            return pd*getPixelMM();
        }
        /**
         * Generate vector of function values (up to 6 - 1 per selected channel), corresponding to motor positions and pixel
         * coordinates, optionally generate partial derivatives for each channel and parameter
         * @param motors array of 3 motors positions
         * @param px pixel X-coordinate of the sample center
         * @param py pixel Y-coordinate of the sample center
         * @param deriv 2d array with outer dimension correspond to number of selected channels (getNumberOfChannels())
         * or null if derivatives are not needed, just values
         * @return array of [getNumberOfChannels()] calculated function values
         */
        public double [] getValsDerivatives(
                int sampleIndex, // double [][] corrPars, // [6][nParZ]
                boolean sagittalMaster, // false - tangential master, true - sagittal master (for center coefficients)
                int [] motors, // 3 motor coordinates
                double px, // pixel x
                double py, // pixel y
                double [][] deriv // array of (1..6[][], matching getNumberOfChannels) or null if derivatives are not required
                ){
            
            double [][] corrPars=getCorrPar(sampleIndex);

            double [] motorDerivs=(deriv==null)? null:(new double [mechanicalFocusingModel.getNumPars()]);
            double [] chnValues=new double [getNumberOfChannels()];
            double mot_z=mechanicalFocusingModel.calc_ZdZ(
                    motors,
                    px,
                    py,
                    motorDerivs);
            int nChn=0;
            double [][] deriv_curv = new double [channelSelect.length][];
            for (int c=0;c<channelSelect.length;c++) deriv_curv[c]=null;
            for (int c=0;c<channelSelect.length;c++) if (channelSelect[c]){
                deriv_curv[c]=(deriv==null)?null:(new double [curvatureModel[c].getSize()]); // nr*nz+1
                chnValues[nChn++]=curvatureModel[c].getFdF(
                        (corrPars==null)?null:corrPars[c], // param_corr
                        px,
                        py,
                        mot_z,
                        deriv_curv[c]);
            }
            if (deriv!=null){
                double dX=px-pXY[0];
                double dY=py-pXY[1];
                double r=Math.sqrt(dX*dX+dY*dY);
                double [] dr_dxy={1.0,1.0};
                if (r>0.0){
                    dr_dxy[0]=-dX/r;
                    dr_dxy[1]=-dY/r;
                }
                nChn=0;
                for (int i=0;i<2;i++){
                    dr_dxy[i]*=getPixelMM(); // radius in mm, dx, dy - in pixels
                }
                for (int c=0;c<channelSelect.length;c++) if (channelSelect[c]){ // c -full, nChn - disabled skipped - dependent COMPONET 
                    boolean isMasterDir=(c%2) == (sagittalMaster?0:1);
                    int otherChannel= c ^ 1;
//                    deriv[nChn]=new double [getNumberOfRegularParameters(sagittalMaster)]; //???????????
                    deriv[nChn]=new double [getNumberOfParameters(sagittalMaster)]; //???????????
                    int np=0;
                    for (int i=0;i<2;i++){
                        if (centerSelect[i] ) {
                            deriv[nChn][np++]=dr_dxy[i]*deriv_curv[c][deriv_curv[c].length-1]; // needs correction from measurement conversion
                        }
                    }
                    // For dependent/master channels no difference for mechanical parameters
                    // TODO: verify they are the same?
                    // calculate derivatives for the center variations (/dXc, /dYc). Will need to be modified for Y vector too as
                    // Y-vector two values (S,T) are calculated from 3 (X2,Y2,XY) and depend on the center position.

                    for (int i=0;i<mechanicalFocusingModel.paramValues.length;i++){
                        if ((mechanicalSelect==null) || mechanicalSelect[i] ) {
                            deriv[nChn][np++]=-motorDerivs[i]*deriv_curv[c][0]; // minus d/dz0 const part
                        }
                    }
                    // other parameters - for dependent - skip center (j==0), for master add dependent
                    for (int n=0;n<channelSelect.length;n++) if (channelSelect[n]){ // n - group of channel parameters
                        int [] ncp=curvatureModel[n].getNumPars(); // {(z),(r)}
                        boolean isDependMasterDir=(n%2) == (sagittalMaster?0:1); // n is master
                        for (int i=0;i<curvatureSelect[n].length; i++) if (curvatureSelect[n][i] ){ // i - parameter number
                            if (((i%ncp[1])!=0) || isDependMasterDir) { // non center or master
                                int dependOnChannel=(((i%ncp[1])==0) && !isMasterDir)?otherChannel:c;
                                deriv[nChn][np++]=(n==dependOnChannel)?(deriv_curv[c][i]):0.0; // deriv[nChn][np++]=(n==dependOnChannel)?(deriv_curv[n][i]):0.0;

                            }
                        }
                    }
                    if (debugLevel==10){
                    	System.out.println("getValsDerivatives(), #="+np+" (of "+deriv[nChn].length+ ")");
                    	for (int ii=0;ii<np;ii++) if (deriv[nChn][ii]!=0.0){
                    		System.out.println("getValsDerivatives(), #="+ii+" (of "+deriv[nChn].length+ ") c="+c+" nChn="+nChn+" deriv="+deriv[nChn][ii]);
                    	}
                    }
                    // add correction parameters?
                    // now np points to the first correction parameter
                    // correction parameters do not depend on sagittalMaster - each mayy have own shift
                    if ((corrPars!=null) && (sampleCorrChnParIndex!=null)){
                        for (int n=0;n<channelSelect.length;n++) if (sampleCorrChnParIndex[n]!=null){
                            int [] ncp=curvatureModel[n].getNumPars(); // {(z),(r)}
                            for (int i=0;i<ncp[0];i++) if (sampleCorrChnParIndex[n][i]>=0){
                                // np now points to the first corr parameter
                                deriv[nChn][np+sampleCorrChnParIndex[n][i]+sampleIndex]=(nChn==n)?deriv_curv[n][i*ncp[1]]:0.0;
                            }
                        }
                    }
                    nChn++;
                }
            }
            return chnValues;
        }
        
        
        
        
        
    }
    
    public class MechanicalFocusingModel{

    	public final String [][] descriptions={
    			{"K0", "Average motor center travel","um/step","0.0124"},
    			{"KD1","M1 and M2 travel disbalance","um/step","0.0"},
    			{"KD3","M3 to average of M1 and M2 travel disbalance","um/step","0.0"},
    			{"sM1","M1: sin component amplitude, relative to tread pitch","","0.0"},
    			{"cM1","M1: cos component amplitude, relative to tread pitch","","0.0"},
    			{"sM2","M2: sin component amplitude, relative to tread pitch","","0.0"},
    			{"cM2","M2: cos component amplitude, relative to tread pitch","","0.0"},
    			{"sM3","M3: sin component amplitude, relative to tread pitch","","0.0"},
    			{"cM3","M3: cos component amplitude, relative to tread pitch","","0.0"},
    			{"Lx", "Half horizontal distance between M3 and and M2 supports", "mm","21.0"},
    			{"Ly", "Half vertical distance between M1 and M2 supports", "mm","10.0"},
    			{"mpX0","pixel X coordinate of mechanical center","px","1296.0"},
    			{"mpY0","pixel Y coordinate of mechanical center","px","968.0"},
    			{"z0", "center shift, positive away from the lens","um","0.0"},
    			{"tx", "horizontal tilt", "um/mm","0.0"},
    			{"ty", "vertical tilt", "um/mm","0.0"}};
    	public double PERIOD=3584.0; // steps/revolution
    	//        public double PIXEL_SIZE=0.0022; // mm
    	public double [] paramValues=new double [descriptions.length];
    	public double [] paramValuesCalibrate=null; // save calibration mode parameters while adjusting
    	public boolean adjustMode=false;

    	public MechanicalFocusingModel(){ // add arguments?
    		initDefaults();
    	}
    	public void setAdjustMode(boolean mode){
    		if (mode==adjustMode) return;
    		if (adjustMode) restoreCalibration();
    		else saveCalibration();
    		adjustMode=mode;
    	}
    	private void saveCalibration(boolean force){ // do not save if in adjust mode
    		if (force) adjustMode=false;
    		saveCalibration();
    	}
    	private void saveCalibration(){ // do not save if in adjust mode
    		if (paramValues!=null){
    			if (!adjustMode) paramValuesCalibrate=paramValues.clone();
    			else System.out.println("Possible BUG - tried to save mechanical parameters while in adjust mode");
    		}
    	}
    	private void restoreCalibration(){
    		if (paramValuesCalibrate!=null)	paramValues=paramValuesCalibrate.clone();
    		adjustMode=false;
    	}

    	public void setProperties(String prefix,Properties properties){
    		setAdjustMode(false); // restore if needed
    		if (paramValues!=null) {
    			for (int i=0;i<paramValues.length;i++){
    				properties.setProperty(prefix+descriptions[i][0],paramValues[i]+"");
    			}
    		}
    	}
    	public void getProperties(String prefix,Properties properties){
    		if ((paramValues==null) || (paramValues.length!=descriptions.length))
    			initDefaults();
    		for (int i=0;i<paramValues.length;i++){
    			if (properties.getProperty(prefix+descriptions[i][0])!=null)
    				paramValues[i]=Double.parseDouble(properties.getProperty(prefix+descriptions[i][0]));
    		}
    		saveCalibration(true);
    	}
    	public void initDefaults(){
    		paramValues=new double [descriptions.length];
    		for (int i=0;i<descriptions.length;i++) paramValues[i]=Double.parseDouble(descriptions[i][3]);
    		saveCalibration(true);
    	}
    	public void setVector(double[] vector){
    		paramValues=new double [descriptions.length];
    		for (int i=0;i<vector.length;i++) paramValues[i]=vector[i];
    	}
    	public void setZTxTy(double z, double tx, double ty){
    		paramValues[getIndex(MECH_PAR.z0)]=z;
    		paramValues[getIndex(MECH_PAR.tx)]=tx;
    		paramValues[getIndex(MECH_PAR.ty)]=ty;
    	}

    	public void setVector(double[] vector, boolean [] mask){
    		for (int i=0;i<vector.length;i++) if (mask[i]) paramValues[i]=vector[i];
    	}

    	public double [] getVector() {return paramValues;}

    	public int getNumPars(){return descriptions.length;}
    	public int getIndex(String parName){
    		for (int i=0;i<descriptions.length;i++) if (descriptions[i].equals(parName)) return i;
    		return -1;
    	}
    	public int getIndex(MECH_PAR mech_par){
    		return mech_par.ordinal();
    	}
    	public String getDescription(int i){
    		return descriptions[i][1];
    	}
    	public String getName(int i){
    		return descriptions[i][0];
    	}
    	public String getDescription(MECH_PAR mech_par){
    		return descriptions[mech_par.ordinal()][1];
    	}
    	public String getUnits(int i){
    		return descriptions[i][2];
    	}
    	public String getUnits(MECH_PAR mech_par){
    		return descriptions[mech_par.ordinal()][2];
    	}
    	public Double getValue(int i){
    		return paramValues[i];
    	}

    	public Double getValue(MECH_PAR mech_par){
    		return paramValues[mech_par.ordinal()];
    	}

    	public double [] debugDeriv_ZdZ(
    			double scale,
    			int [] motors,
    			double px,
    			double py){
    		double [] derivSteps ={
    				0.001, // [0]	K0 843.6001052876973	
    				0.001, // [1]	KD1 1600.0723700390713	
    				0.001, // [2]	KD3 -2428.8998947123027	
    				0.1,   // [3]	sM1 -1.1720600074585734	
    				0.1,   // [4]	cM1 0.6081060292866338	
    				0.1,   // [5]	sM2 -2.3717384278673035	
    				0.1,   // [6]	cM2 1.2305414643438266	
    				0.1,   // [7]	sM3 2.7345656468251707	
    				0.1,   // [8]	cM3 1.4187890097199358	
    				1.0,   // [9]	Lx -0.535718420298317	
    				1.0,   // [10]  Ly 	0.0	
    				10.0,  // [11]  mpX0	-2.816702366108815E-5	
    				10.0,  // [12]  mpY0	-8.351458550253758E-5	
    				1.0,   // [13]  z0	 1.0	
    				0.1,   // [14]  tx -2.5168000000000004	
    				0.1    // [15]  ty -1.76	
    		};
    		double [] initialVector=paramValues.clone();
    		double [] derivs=new double [paramValues.length];
    		for (int i=0;i<derivs.length;i++){
    			paramValues[i]=initialVector[i]-scale*derivSteps[i];
    			double zm=calc_ZdZ(motors,px,py,null);
    			paramValues[i]=initialVector[i]+scale*derivSteps[i];
    			double zp=calc_ZdZ(motors,px,py,null);
    			paramValues[i]=initialVector[i];
    			derivs[i]=(zp-zm)/(2.0*scale*derivSteps[i]);
    		}
    		return derivs;

    	}


    	/**
    	 * Calculate distance from the selected sensor pixel to the common "focal" plane
    	 * and optionally partial derivatives of the pixel distance from the focal plane by selected parameters
    	 * @param motors array of 3 motor positions
    	 * @param px horizontal sensor position
    	 * @param py vertical sensor pixel position
    	 * @param deriv returns partial derivatives if array is provided. If null - does not calculate derivatives
    	 * @return array of partial derivatives
    	 */
    	public double calc_ZdZ(
    			int [] motors,
    			double px,
    			double py,
    			double[] deriv){
//    			double[][] adjData){
    		double debugMot=6545;
    		int debugThreshold=2;
    		boolean dbg = (debugLevel>debugThreshold);
    		/*
             kM3=K0+KD3
             kM1=K0+KD1-KD3
             kM2=K0-KD1-KD3
             K0 may be fixed (overall scale), kM1, kM2 - variable

            dZc(m1)= 0.25* kM1 * (m1 + sM1*P/(2*pi)*sin(2pi*m1/P) + cM1*P/(2*pi)*cos(2pi*m1/P))
            dZc(m2)= 0.25* kM2 * (m2 + sM2*P/(2*pi)*sin(2pi*m2/P) + cM2*P/(2*pi)*cos(2pi*m2/P))
            dZc(m3)= 0.5 * kM3 *( m3 + sM3*P/(2*pi)*sin(2pi*m3/P) + cM3*P/(2*pi)*cos(2pi*m3/P))

            Assuming X - towards M3, Y - towards M1

            d2Z/dX/dm1=-ps/(4*Lx)* kM1 * (m1 + sM1*P/(2*pi)*sin(2pi*m1/P) + cM1*P/(2*pi)*cos(2pi*m1/P))
            d2Z/dX/dm2=-ps/(4*Lx)* kM2 *( m2 + sM2*P/(2*pi)*sin(2pi*m2/P) + cM2*P/(2*pi)*cos(2pi*m2/P))
            d2Z/dX/dm3= ps/(2*Lx)* kM3 * (m3 + sM3*P/(2*pi)*sin(2pi*m3/P) + cM3*P/(2*pi)*cos(2pi*m3/P))

            d2Z/dY/dm1=+ps/(2*Ly)* kM1 * (m1 + sM1*P/(2*pi)*sin(2pi*m1/P) + cM1*P/(2*pi)*cos(2pi*m1/P))
            d2Z/dY/dm2= -ps/(2*Ly)* kM2 * (m2 + sM2*P/(2*pi)*sin(2pi*m2/P) + cM2*P/(2*pi)*cos(2pi*m2/P))
            d2Z/dY/dm3= 0

    		 */
    		//            double [] deriv=new double [paramValues.length];
    		double kM1=    getValue(MECH_PAR.K0)+getValue(MECH_PAR.KD1)-getValue(MECH_PAR.KD3);
    		double kM2=    getValue(MECH_PAR.K0)-getValue(MECH_PAR.KD1)-getValue(MECH_PAR.KD3);
    		double kM3=    getValue(MECH_PAR.K0)+getValue(MECH_PAR.KD3);
    		double p2pi= PERIOD/2/Math.PI;
    		double m1=motors[0],m2=motors[1],m3=motors[2];
    		double aM1=(m1 + getValue(MECH_PAR.sM1)*p2pi*Math.sin(m1/p2pi) + getValue(MECH_PAR.cM1)*p2pi*Math.cos(m1/p2pi));
    		double aM2=(m2 + getValue(MECH_PAR.sM2)*p2pi*Math.sin(m2/p2pi) + getValue(MECH_PAR.cM2)*p2pi*Math.cos(m2/p2pi));
    		double aM3=(m3 + getValue(MECH_PAR.sM3)*p2pi*Math.sin(m3/p2pi) + getValue(MECH_PAR.cM3)*p2pi*Math.cos(m3/p2pi));
    		double zM1=kM1 * aM1;
    		double zM2=kM2 * aM2;
    		double zM3=kM3 * aM3;

    		double zc= 0.25* zM1+ 0.25* zM2+ 0.5 * zM3+getValue(MECH_PAR.z0);
    		double dx=PIXEL_SIZE*(px-getValue(MECH_PAR.mpX0));
    		double dy=PIXEL_SIZE*(py-getValue(MECH_PAR.mpY0));
    		double zx=dx*(getValue(MECH_PAR.tx)+(2*zM3-zM1-zM2)/(4*getValue(MECH_PAR.Lx))) ;
    		//            double zy=dy*(getValue(MECH_PAR.ty)+(zM1-zM2)/(2*getValue(MECH_PAR.Ly)));
    		double zy=dy*(getValue(MECH_PAR.ty)+(zM2-zM1)/(2*getValue(MECH_PAR.Ly)));
    		double z=zc+zx+zy;
    		if (dbg) if ((Math.abs(m1)==debugMot)&& (Math.abs(m2)==debugMot)){
    			System.out.print ("M: "+((int)m1)+":"+((int)m2)+":"+((int)m3)+
    					" dxy="+IJ.d2s(dx,3)+":"+IJ.d2s(dy,3)+" zcxy="+IJ.d2s(zc,3)+":"+IJ.d2s(zx,3)+":"+IJ.d2s(zy,3)+
    					" zxy(t)="+IJ.d2s(dx*getValue(MECH_PAR.tx),3)+":"+IJ.d2s(dy*getValue(MECH_PAR.ty),3));
    		}
    		if (deriv==null) {
    			if (dbg) if ((Math.abs(m1)==debugMot)&& (Math.abs(m2)==debugMot)){
    				System.out.println();
    			}
    			return z;
    		}
    		for (int i=0;i<deriv.length;i++) deriv[i]=0.0;
    		// Above same as calc_Z
    		double dx_mpX0=-PIXEL_SIZE;
    		double dy_mpY0=-PIXEL_SIZE;
    		double zM1_K0= aM1;
    		double zM1_KD1= aM1;
    		double zM1_KD3=-aM1;
    		double zM2_K0= aM2;
    		double zM2_KD1=-aM2;
    		double zM2_KD3=-aM2;
    		double zM3_K0= aM3;
    		//            double zM3_KD1=-aM3;
    		//            double zM3_KD3= 0.0;
    		double zM3_KD1= 0.0;
    		double zM3_KD3= aM3;
    		double zM1_sM1= kM1*p2pi*Math.sin(m1/p2pi);
    		double zM1_cM1= kM1*p2pi*Math.cos(m1/p2pi);
    		double zM2_sM2= kM2*p2pi*Math.sin(m2/p2pi);
    		double zM2_cM2= kM2*p2pi*Math.cos(m2/p2pi);
    		double zM3_sM3= kM3*p2pi*Math.sin(m3/p2pi);
    		double zM3_cM3= kM3*p2pi*Math.cos(m3/p2pi);

    		double zc_K0= 0.25* zM1_K0+ 0.25* zM2_K0+ 0.5 * zM3_K0;
    		double zc_KD1= 0.25* zM1_KD1+ 0.25* zM2_KD1+ 0.5 * zM3_KD1;
    		double zc_KD3= 0.25* zM1_KD3+ 0.25* zM2_KD3+ 0.5 * zM3_KD3;
    		double zc_sM1= 0.25* zM1_sM1;
    		double zc_cM1= 0.25* zM1_cM1;
    		double zc_sM2= 0.25* zM2_sM2;
    		double zc_cM2= 0.25* zM2_cM2;
    		double zc_sM3= 0.5* zM3_sM3;
    		double zc_cM3= 0.5* zM3_cM3;

    		//            double zx_K0=(2*zM3-zM1-zM2)* dx/(4*getValue(MECH_PAR.Lx));
    		double zx_a=dx/(4*getValue(MECH_PAR.Lx));
    		double zx_K0= (2*zM3_K0-zM1_K0-zM2_K0)*zx_a;
    		double zx_KD1=(2*zM3_KD1-zM1_KD1-zM2_KD1)*zx_a;
    		double zx_KD3=(2*zM3_KD3-zM1_KD3-zM2_KD3)*zx_a;
    		double zx_sM1= (-zM1_sM1)*zx_a;
    		double zx_cM1= (-zM1_cM1)*zx_a;
    		double zx_sM2= (-zM2_sM2)*zx_a;
    		double zx_cM2= (-zM2_cM2)*zx_a;
    		double zx_sM3= (2*zM3_sM3)*zx_a;
    		double zx_cM3= (2*zM3_cM3)*zx_a;
    		double zx_mpX0=dx_mpX0*(getValue(MECH_PAR.tx)+(2*zM3-zM1-zM2)/(4*getValue(MECH_PAR.Lx))); //  double zx_mpX0=dx_mpX0/(4*getValue(MECH_PAR.Lx));
    		double zx_tx= dx;
    		double zx_Lx= -dx*(2*zM3-zM1-zM2)/(4*getValue(MECH_PAR.Lx)*getValue(MECH_PAR.Lx));
    		//          double zy=dy*(getValue(MECH_PAR.ty)+(zM2-zM1)/(2*getValue(MECH_PAR.Ly)));
    		double zy_a=dy/(2*getValue(MECH_PAR.Ly));
    		double zy_K0=  (zM2_K0- zM1_K0) *zy_a;
    		double zy_KD1= (zM2_KD1-zM1_KD1)*zy_a;
    		double zy_KD3= (zM2_KD3-zM1_KD3)*zy_a;
    		double zy_sM1= (-zM1_sM1)*zy_a;
    		double zy_cM1= (-zM1_cM1)*zy_a;
    		double zy_sM2= (zM2_sM2)*zy_a;
    		double zy_cM2= (zM2_cM2)*zy_a;
    		double zy_sM3= 0.0;
    		double zy_cM3= 0.0;
    		double zy_mpY0=dy_mpY0*(getValue(MECH_PAR.ty)+(zM2-zM1)/(2*getValue(MECH_PAR.Ly))); // double zy_mpY0=dy_mpY0/(2*getValue(MECH_PAR.Ly));
    		double zy_ty= dy;
    		//            double zy_Ly= -dy*(zM1-zM2)/(2*getValue(MECH_PAR.Ly)*getValue(MECH_PAR.Ly));
    		double zy_Ly= -dy*(zM2-zM1)/(2*getValue(MECH_PAR.Ly)*getValue(MECH_PAR.Ly));

    		deriv[getIndex(MECH_PAR.K0)]= zc_K0+zx_K0+zy_K0;
    		deriv[getIndex(MECH_PAR.KD1)]=zc_KD1+zx_KD1+zy_KD1;
    		deriv[getIndex(MECH_PAR.KD3)]=zc_KD3+zx_KD3+zy_KD3;
    		deriv[getIndex(MECH_PAR.sM1)]=zc_sM1+zx_sM1+zy_sM1;
    		deriv[getIndex(MECH_PAR.cM1)]=zc_cM1+zx_cM1+zy_cM1;
    		deriv[getIndex(MECH_PAR.sM2)]=zc_sM2+zx_sM2+zy_sM2;
    		deriv[getIndex(MECH_PAR.cM2)]=zc_cM2+zx_cM2+zy_cM2;
    		deriv[getIndex(MECH_PAR.sM3)]=zc_sM3+zx_sM3+zy_sM3;
    		deriv[getIndex(MECH_PAR.cM3)]=zc_cM3+zx_cM3+zy_cM3;
    		deriv[getIndex(MECH_PAR.Lx)] = zx_Lx;
    		deriv[getIndex(MECH_PAR.Ly)] = zy_Ly;
    		deriv[getIndex(MECH_PAR.mpX0)] = zx_mpX0;
    		deriv[getIndex(MECH_PAR.mpY0)] = zy_mpY0;
    		deriv[getIndex(MECH_PAR.z0)] = 1.0;
    		deriv[getIndex(MECH_PAR.tx)] = zx_tx;
    		deriv[getIndex(MECH_PAR.ty)] = zy_ty;
    		if (dbg) if ((Math.abs(m1)==debugMot)&& (Math.abs(m2)==debugMot)){
    			if (m1*m2>0){
    				System.out.println("same sign");
    			} else {
    				System.out.println("opposite sign");
    			}
    			System.out.println (" zxy_txy="+IJ.d2s(zx_tx,3)+":"+IJ.d2s(zy_ty,3)+" zxy_Lxy="+IJ.d2s(zx_Lx,5)+":"+IJ.d2s(zy_Ly,5));
    			double [] dbg_derivs= debugDeriv_ZdZ(
    					0.1, //scale,
    					motors,
    					px,
    					py);
    			for (int i=0;i<deriv.length;i++){
    				System.out.println(i+": "+descriptions[i][0]+" deriv="+deriv[i]+", dbg_derivs="+dbg_derivs[i]);
    			}
    		}
    		return z;
    	}
    	
    	/**
    	 * Calculate linearized mount (motor) displacement from motor position in steps
    	 * @param m motor position in steps
    	 * @param index motor number (0..2)
    	 * @return mount displacement (in microns)
    	 */
    	public double mToZm(
    			double m,
    			int index) {
    		double p2pi= PERIOD/2/Math.PI;
    		double kM=Double.NaN,aM=Double.NaN;
    		switch (index){
    		case 0:
        		kM=    getValue(MECH_PAR.K0)+getValue(MECH_PAR.KD1)-getValue(MECH_PAR.KD3);
        		aM=(m + getValue(MECH_PAR.sM1)*p2pi*Math.sin(m/p2pi) + getValue(MECH_PAR.cM1)*p2pi*Math.cos(m/p2pi));
        		break;
    		case 1:
        		kM=    getValue(MECH_PAR.K0)-getValue(MECH_PAR.KD1)-getValue(MECH_PAR.KD3);
        		aM=(m + getValue(MECH_PAR.sM2)*p2pi*Math.sin(m/p2pi) + getValue(MECH_PAR.cM2)*p2pi*Math.cos(m/p2pi));
        		break;
    		case 2:
        		kM=    getValue(MECH_PAR.K0)+getValue(MECH_PAR.KD3);
        		aM=(m + getValue(MECH_PAR.sM3)*p2pi*Math.sin(m/p2pi) + getValue(MECH_PAR.cM3)*p2pi*Math.cos(m/p2pi));
        		break;
    		}
    		return kM*aM;
    	}
    	
    	private double getDzmDm(
    			double m,
    			double kM,
    			double s,
    			double c) {
    		double p2pi= PERIOD/2/Math.PI;
    		return kM*(1.0+ s*Math.cos(m/p2pi)-c*Math.sin(m/p2pi));
    		
    	}

    	/**
    	 * Convert linearized motor position to motor steps using current mechanical parameter values
    	 * @param zm linear motor position
    	 * @param index motor index (0..2)
    	 * @return motor position in steps
    	 */
    	public double zmToM(
    			double zm,
    			int index) {
    		double p2pi= PERIOD/2/Math.PI;
    		double kM=Double.NaN, m=0,s=0.0,c=1.0;
    		switch (index){
    		case 0:
        		kM=    getValue(MECH_PAR.K0)+getValue(MECH_PAR.KD1)-getValue(MECH_PAR.KD3);
        		s=getValue(MECH_PAR.sM1);
        		c=getValue(MECH_PAR.cM1);
        		break;
    		case 1:
        		kM=    getValue(MECH_PAR.K0)-getValue(MECH_PAR.KD1)-getValue(MECH_PAR.KD3);
        		s=getValue(MECH_PAR.sM2);
        		c=getValue(MECH_PAR.cM2);
        		break;
    		case 2:
        		kM=    getValue(MECH_PAR.K0)+getValue(MECH_PAR.KD3);
        		s=getValue(MECH_PAR.sM3);
        		c=getValue(MECH_PAR.cM3);
        		break;
    		}
    		m=zm/kM; // without thread sin/cos
    		double eps=0.1;
    		int maxRetries=100;
//     		aM=(m + getValue(MECH_PAR.sM3)*p2pi*Math.sin(m/p2pi) + getValue(MECH_PAR.cM3)*p2pi*Math.cos(m/p2pi));
//	double dzM1_dm1=kM1*(1.0+getValue(MECH_PAR.sM1)*Math.cos(m1/p2pi)-getValue(MECH_PAR.cM1)*Math.sin(m1/p2pi));
    		double a=Math.sqrt(s*s+c*c);
    		double phase=0.0;
    		double extrenPhase=0.0;
    		if (a>1.0){
    			phase=Math.atan2(s,c);
//    			double rd=1.0+ s*Math.cos(m/p2pi)-c*Math.sin(m/p2pi);
//    			double rd=1.0+ a*Math.cos((m+p2pi*phase)/p2pi);
// rd==0 => Math.cos((m+p2pi*phase)/p2pi)=-1.0/a
    			extrenPhase=phase-Math.asin(1.0/a);
//    			minPhase-=(2*Math.PI)*Math.floor(minPhase/(2*Math.PI));
    			extrenPhase-=Math.PI*Math.floor(extrenPhase/Math.PI);  // min/max 0<=hase<PI of the min/max zm(m)
    		}
    		for (int retry=0;retry<maxRetries;retry++){
    			double zm1=mToZm(m, index);
    			if (Math.abs(zm1-zm)<eps) break;
        		double dzm_dm=getDzmDm(m,kM,s,c);
        		if (dzm_dm<=0.0){
        			if (debugLevel>0) System.out.println("Negative derivative dzm_dm");
        			double mPhase=m/p2pi;
        			double halfPeriods=Math.floor(mPhase/Math.PI);
        			double fractPhase=mPhase-Math.PI*halfPeriods;
        			double mirrExtrenFractPhase=extrenPhase;
        			if (mirrExtrenFractPhase<fractPhase) mirrExtrenFractPhase+=Math.PI;
        			if (zm1>zm) mirrExtrenFractPhase-=Math.PI;
        			double mirrPhase=Math.PI*halfPeriods+mirrExtrenFractPhase;
        			double mirrM=p2pi*mirrPhase;
        			m=2*mirrM-m; // symmetrical aqround extrenum, up if needed more, down - if needed less
        		}
        		double l=1.0;
        		double m2=m;
    			for (int retry2=0;retry2<maxRetries;retry2++){
        			m2 = m+l*((zm-zm1)/dzm_dm);
        			double zm2= mToZm(m2, index);
    				if (Math.abs(zm2-zm)<Math.abs(zm1-zm)) break;
    				l*=0.5;
    			}
    			m=m2;
    		}
    		return m;
    	}

    	/**
    	 * Calculate three linearized values of motor positions for current parameters, target center focal shift and tilt 
    	 * @param px lens center X (pixels)
    	 * @param py lens center Y (pixels)
    	 * @param targetZ target focal shift uin the center, microns, positive - away
    	 * @param targetTx target horizontal tilt (normally 0)
    	 * @param targetTy target vertical tilt (normallty 0) 
    	 * @return array of 3 linearized motor positions (microns) 
    	 */
    	public double [] getZM(
    			double px,
    			double py,
    			double targetZ,
    			double targetTx,
    			double targetTy){
    		double dx=PIXEL_SIZE*(px-getValue(MECH_PAR.mpX0));
    		double dy=PIXEL_SIZE*(py-getValue(MECH_PAR.mpY0));
    		//    		double zc= 0.25* zM1+ 0.25* zM2+ 0.5 * zM3+getValue(MECH_PAR.z0);
    		//    		double zx=dx*(getValue(MECH_PAR.tx)+(2*zM3-zM1-zM2)/(4*getValue(MECH_PAR.Lx))) ;
    		//    		double zy=dy*(getValue(MECH_PAR.ty)+(zM2-zM1)/(2*getValue(MECH_PAR.Ly)));
    		//          double z=zc+zx+zy    		
    		//			A*{zM1,zM2,zM3}={targetZ,targetTx,targetTy}
    		//    		A*{zM1,zM2,zM3}={targetZ-getValue(MECH_PAR.z0),targetTx-dx*getValue(MECH_PAR.tx),targetTy-dy*getValue(MECH_PAR.ty)}
    		double [][] A={
    				{
    					0.25 - dx/(4*getValue(MECH_PAR.Lx)) -dy/(2*getValue(MECH_PAR.Ly)),
    					0.25 - dx/(4*getValue(MECH_PAR.Lx)) +dy/(2*getValue(MECH_PAR.Ly)),
    					0.5  + dx/(2*getValue(MECH_PAR.Lx))
    				} , {
    					-1.0/(4*getValue(MECH_PAR.Lx)),
    					-1.0/(4*getValue(MECH_PAR.Lx)),
    					1.0/ (2*getValue(MECH_PAR.Lx))
    				} , {
    					-1.0/(2*getValue(MECH_PAR.Ly)),
    					1.0/ (2*getValue(MECH_PAR.Ly)),
    					0.0
    				}
    		};
    		double [][] B={
    			{targetZ-getValue(MECH_PAR.z0)},
    			{targetTx-dx*getValue(MECH_PAR.tx)},
    			{targetTy-dy*getValue(MECH_PAR.ty)}
    		};
    		Matrix MA=new Matrix(A);
    		Matrix MB=new Matrix(B);
    		Matrix S=MA.solve(MB);
    		return S.getColumnPackedCopy();
    	}
    	
    	public boolean[] getDefaultMask(){
    		boolean [] mask = new boolean[this.paramValues.length];
    		for (int i=0;i<mask.length;i++) mask[i]=false;
    		mask[getIndex(MECH_PAR.K0)]= false;
    		mask[getIndex(MECH_PAR.KD1)]=false; //true;
    		mask[getIndex(MECH_PAR.KD3)]=false; //true;
    		mask[getIndex(MECH_PAR.sM1)]=false;
    		mask[getIndex(MECH_PAR.cM1)]=false;
    		mask[getIndex(MECH_PAR.sM2)]=false;
    		mask[getIndex(MECH_PAR.cM2)]=false;
    		mask[getIndex(MECH_PAR.sM3)]=false;
    		mask[getIndex(MECH_PAR.cM3)]=false;
    		mask[getIndex(MECH_PAR.Lx)] =false;
    		mask[getIndex(MECH_PAR.Ly)] =false; //true;
    		mask[getIndex(MECH_PAR.mpX0)] = false; //true;
    		mask[getIndex(MECH_PAR.mpY0)] = false; //true;
    		mask[getIndex(MECH_PAR.tx)] = true;
    		mask[getIndex(MECH_PAR.ty)] = true;
    		return mask;
    	}

    	public boolean[] maskSetDialog(String title, boolean [] currentMask){
    		GenericDialog gd = new GenericDialog(title);
    		boolean [] mask = new boolean[this.paramValues.length];
    		if (currentMask==null) currentMask=getDefaultMask();
    		for (int i=0;i<mask.length;i++) {
    			mask[i]=currentMask[i];
    			gd.addCheckbox(getDescription(i), mask[i]);                    
    		}
    		gd.enableYesNoCancel("Apply","Keep"); // default OK (on enter) - "Apply"
    		gd.showDialog();
    		if (gd.wasCanceled()) return null;
    		if (gd.wasOKed()) { // selected default "Apply"
    			for (int i=0;i<mask.length;i++) {
    				mask[i]=gd.getNextBoolean();
    			}
    		}
    		return mask;
    	}
    	public boolean [] maskSetZTxTy(){
    		boolean [] mask = new boolean[this.paramValues.length];
    		for (int i=0;i<mask.length;i++)	mask[i]=false;
    		mask[getIndex(MECH_PAR.z0)]=true;
    		mask[getIndex(MECH_PAR.tx)]=true;
    		mask[getIndex(MECH_PAR.ty)]=true;
    		return mask;
    	}
    	public boolean showModifyParameterValues(String title, boolean showDisabled, boolean [] mask){
    		GenericDialog gd = new GenericDialog(title);
    		for (int i=0;i<this.paramValues.length;i++){
    			if ((mask==null) || mask[i] ) {
    				gd.addNumericField(getDescription(i),this.paramValues[i],5,8,getUnits(i));
    			} else if (showDisabled){
    				//                    gd.addMessage(getDescription(i) +": "+this.paramValues[i]+" ("+getUnits(i)+")");
    				gd.addNumericField("(disabled) "+getDescription(i),this.paramValues[i],5,8,getUnits(i));
    			}
    		}
    		gd.enableYesNoCancel("Apply","Keep"); // default OK (on enter) - "Apply"
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		if (gd.wasOKed()) { // selected default "Apply"
    			for (int i=0;i<this.paramValues.length;i++){
    				if ((mask==null) || mask[i] || showDisabled) {
    					this.paramValues[i]=gd.getNextNumber();
    				}
    			}
    		}         
    		return true;    
    	}        
    }
    
    public class CurvatureModel{
        private double dflt_na=Math.log(0.15); // um/um
        private double dflt_r0=Math.log(4.0); //3.3; // um (1.5 pix)
        private double [][] modelParams=null;
        public static final int dflt_distanceParametersNumber=5;
        public static final int dflt_radialParametersNumber=4;
        private static final int min_distanceParametersNumber=4;
        private static final int min_radialParametersNumber=1;
        private double pX0,pY0;
        public CurvatureModel(
                double pX0,
                double pY0,
                int distanceParametersNumber,
                int radialParametersNumber){
            this.pX0=pX0;
            this.pY0=pY0;
            if (distanceParametersNumber<min_distanceParametersNumber) distanceParametersNumber=min_distanceParametersNumber;
            if (radialParametersNumber<min_radialParametersNumber) radialParametersNumber=min_radialParametersNumber;
            this.modelParams=new double [distanceParametersNumber][radialParametersNumber];
            setDefaults();
        }
    public void setProperties(String prefix,Properties properties){
            properties.setProperty(prefix+"distanceParametersNumber",modelParams.length+"");
            properties.setProperty(prefix+"radialParametersNumber",modelParams[0].length+"");
            properties.setProperty(prefix+"pX0",pX0+"");
            properties.setProperty(prefix+"pY0",pY0+"");
            for (int i=0;i<modelParams.length;i++) for (int j=0;j<modelParams[i].length;j++){
                properties.setProperty(prefix+"modelParams"+i+"_"+j,modelParams[i][j]+"");
            }
    }
    public void getProperties(String prefix,Properties properties){
        int numZPars,numRPars;
        if (properties.getProperty(prefix+"pX0")!=null)
            pX0=Double.parseDouble(properties.getProperty(prefix+"pX0"));
        if (properties.getProperty(prefix+"pY0")!=null)
            pX0=Double.parseDouble(properties.getProperty(prefix+"pY0"));
        if (modelParams!=null){
            numZPars=modelParams.length;
            numRPars=modelParams[0].length;
        } else {
            numZPars=dflt_distanceParametersNumber;
            numRPars=dflt_radialParametersNumber;
        }
        if (properties.getProperty(prefix+"distanceParametersNumber")!=null)
            numZPars=Integer.parseInt(properties.getProperty(prefix+"distanceParametersNumber"));
        if (properties.getProperty(prefix+"radialParametersNumber")!=null)
            numRPars=Integer.parseInt(properties.getProperty(prefix+"radialParametersNumber"));
        if ((modelParams!=null) || (modelParams.length!=numZPars)|| (modelParams[0].length!=numRPars)){
            modelParams=new double [numZPars][numRPars];
            setDefaults();
        }

        for (int i=0;i<modelParams.length;i++) for (int j=0;j<modelParams[i].length;j++){
            if (properties.getProperty(prefix+"modelParams"+i+"_"+j)!=null) {
                modelParams[i][j]=Double.parseDouble(properties.getProperty(prefix+"modelParams"+i+"_"+j));
            }
        }
    }

    
        
        public void setCenter(double pX, double pY){
            pX0=pX;
            pY0=pY;
        }
        public int [] getNumPars(){
            if (modelParams==null) return null;
            int [] numPars={modelParams.length,modelParams[0].length};
            return numPars;
        }
        public void setDefaults(){
            for (int i=0;i<this.modelParams.length;i++) for (int j=0;j<this.modelParams[0].length;j++){
                this.modelParams[i][j]=0.0;
            }
            this.modelParams[1][0]=dflt_na;
            this.modelParams[2][0]=dflt_r0;
            this.modelParams[0][0]=Double.NaN;
        }
        public boolean z0IsValid(){
        	return !Double.isNaN(this.modelParams[0][0]);
        }
        
        public void set_z0(double z0){
            this.modelParams[0][0]=z0;
        }
        
        public double [] getCenterVector(){
            double [] vector=new double [this.modelParams.length];
            for (int i=0;i<this.modelParams.length;i++){
                vector[i]=this.modelParams[i][0];
            }
            return vector;
        }
        public void setCenterVector(double [] vector){
            for (int i=0;i<this.modelParams.length;i++){
                this.modelParams[i][0]=vector[i];
            }
        }

        public double [] getVector(){
            double [] vector=new double [this.modelParams.length*this.modelParams[0].length];
            int index=0;
            for (int i=0;i<this.modelParams.length;i++) for (int j=0;j<this.modelParams[0].length;j++){
                vector[index++]=this.modelParams[i][j];
            }
            return vector;
        }
        public void setVector(double [] vector, boolean[] mask){ // mask may be null
            int index=0;
            for (int i=0;i<this.modelParams.length;i++) for (int j=0;j<this.modelParams[0].length;j++) if ((mask==null) || mask[index]){
                this.modelParams[i][j]=vector[index++];
            }
        }
        public int getSize(){
            try {
                return this.modelParams.length*this.modelParams[0].length+1; // last element df/dr
            } catch (Exception e){
                return 0;
            }
        }
        double [] getAr(
                double r,
                double [] param_corr
                ){
            double [] ar= new double [this.modelParams.length];
            for (int i=0;i<this.modelParams.length;i++){
                ar[i]=this.modelParams[i][0];
                if (param_corr!=null) ar[i]+=param_corr[i];
//                double rp=1.0;
                double rp=r;
                for (int j=1;j<this.modelParams[i].length;j++){
//                    rp*=r2;
                    rp*=r;
                    ar[i]+=this.modelParams[i][j]*rp;
                }
            }
            return ar;
        }
        /**
         * Calculate function value and (optionally) partial derivatives over all parameters as 1d array
         * inner scanning by radial coefficient (for even powers of radius), outer - for f(z) parameters
         * first parameter - z0 (shift), then scale (related to numeric aperture), then r0 (related to minimal PSF radius),
         * other parameters - just polynomial coefficients for (z_in-z0)^n
         * @param param_corr - per-sample corrections, added to this.modelParams[i][0] - or null
         * @param pX - sample pixel x coordinate (currently just for radius calculation) or radius in mm (if pY is NaN)
         * @param pY - sample pixel y coordinate (currently just for radius calculation) or NaN, then pX - radius in mm
         * @param z_in - z (from mechanical) before subtracting of z0
         * @param deriv array to accommodate all partial derivatives or null (should have modelParams.length*modelParams[0].length+1
         * @return function value
         */
        public double getFdF(
                double [] param_corr,
                double pX, // if isNaN(pY), then pX is radius in mm
                double pY,
                double z_in,
                double [] deriv){
/*
f=sqrt(( a*(zin-z0))^2 + r0^2)+a0+ a1*(zin-z0)+...aN*(zin-z0)^N
each of the z0,z0,a,a[i] is polynomial of even powers of r (r^0, r^2, r^4...)
z=z_in-z0

modified parameters, r0 - PSF FWHM at z=0, k (instead of a0), so that old r0 now reff=
r0*exp(-k), old a0= r0*(1-exp(-k)).

f=sqrt((a*(zin-z0))^2 + (r0*(exp(-k))^2)+r0*(1-exp(-k))+ a1*(zin-z0)+...aN*(zin-z0)^N

z0 - ar[0]
a - ar[1]
r0 - ar[2]
k - ar[3]
ar1 - ar[4]

modified to avoid zero-crossing for a and r0:
a=exp(ln_a)
r0=exp(ln_r0)

z0    - ar[0]
ln_a  - ar[1]
ln_r0 - ar[2]
k     - ar[3]
ar1   - ar[4]
=================
f(x)=sqrt((ax)^2+r^2)+kx-f_corr
f(z)=sqrt((a*(z-z_corr)^2+r^2)+k*(z-z_corr)-f_corr
f'(x)=0 -> x=sqrt(1/r^2-a^2)
z_corr=(kx*rc/a)/*sqrt(a^2-kx^2)
f_corr=sqrt((a*z_corr)^2+r_eff^2)-kx*(z_corr) – r_eff
k=a*func(tilt]), func(-inf)=-1, func(0)=0, func(+inf)=1 
k=a*2/pi*atan(tilt);
d_k/d_tilt = 2/pi*(1/tilt^2)

x=z_in-(z0+z_corr)

Ideally I should match second derivative near minimum to isolate tilt from ar[3] - adjust r_eff and then shift
Maybe not, just keep iot as it is (matching f" at minimum while tilting will cause a sharp turn nearby)
a=exp(ln_a)
r0=exp(ln_r0)
z0    - ar[0]
ln_a  - ar[1]
ln_r0 - ar[2]
k     - ar[3]
a1    - ar[4]

kx=a*2/pi*atan(a1);
z_corr=(kx*rc/a)/sqrt(a^2-kx^2)
f_corr=sqrt((a*z_corr)^2+r_eff^2)-kx*(z_corr) – r_eff
f=sqrt((a*(zin-z0-z_corr))^2 + (r0*(exp(-k))^2)+r0*(1-exp(-k))+ kx*(zin-z0-z_corr) {+...aN*(zin-z0-z_corr)^N} - {} are not likely to be ever used 

dependence:
kx: a,a1 (ar[1], ar[4]
z_corr: kx,r_eff,a - ar[1], ar[4], ar[2], ar[3]
f_corr: d_fcorr/d_zcorr=0, other: a, reff, kx ->  ar[1], ar[2], ar[3],  ar[4]
*/
            double r=Double.isNaN(pY)?pX:Math.sqrt((pX-pX0)*(pX-pX0)+(pY-pY0)*(pY-pY0))*PIXEL_SIZE; // in mm
            double [] ar= new double [this.modelParams.length];
            for (int i=0;i<this.modelParams.length;i++){
                ar[i]=this.modelParams[i][0];
                if (param_corr!=null) ar[i]+=param_corr[i];
                double rp=r;
                for (int j=1;j<this.modelParams[i].length;j++){
                    rp*=r;
                    ar[i]+=this.modelParams[i][j]*rp;
                }
            }
            double exp_a=Math.exp(ar[1]);
            double exp_r=Math.exp(ar[2]);
            double exp_mk=Math.exp(-ar[3]);
            double reff=exp_r*exp_mk;
            double reff2=reff*reff;
            double exp_a2=exp_a*exp_a;
            double a1=ar[4];
//            kx=a*2/pi*atan(a1);
            double kx=exp_a*2.0/Math.PI*Math.atan(a1);
//            z_corr=(kx*rc/a)/*sqrt(a^2-kx^2)
            double z_corr=(kx*reff/exp_a)/Math.sqrt(exp_a2 - kx*kx);
            double z_corr2=z_corr*z_corr;
//          f_corr=sqrt((a*z_corr)^2+r_eff^2)-kx*(z_corr) – r_eff
            double f_corr=Math.sqrt(exp_a2*z_corr2 + reff2)-kx*z_corr - reff;
            double z=z_in-ar[0]-z_corr;
            double sqrt=Math.sqrt(exp_a2*z*z + reff2);
//            f=sqrt((a*(zin-z0-z_corr))^2 + (r0*(exp(-k))^2)+r0*(1-exp(-k))+ kx*(zin-z0-z_corr)-f_corr {+...aN*(zin-z0-z_corr)^N} - {} are not likely to be ever used 
            double f=sqrt+exp_r*(1-exp_mk) +kx*z -f_corr;
            double zp=z;
            for (int i=5;i<ar.length;i++){
                zp*=z;
                f+=ar[i]*zp;
            }

            if (deriv==null) return f; // only value, no derivatives
            
            
// derivatives calculation independent of z - move to a separate function that can be called once for channel/sample, stored asnd then applied to multiple z measurements        
//  double kx=exp_a*2.0/Math.PI*Math.atan(a1);
            double dkx_dar1=kx; // exp_a*2.0/Math.PI*Math.atan(a1);
            double dkx_dar4=exp_a*2.0/Math.PI/(1.0+ a1*a1);
            // z_corr=(kx*reff/exp_a)/Math.sqrt(exp_a*exp_a - kx*kx)  //kx,r_eff,a - ar[1], ar[4], ar[2], ar[3]
            //dzcorr_dkx=rc*a/(a^2-kx^2)/sqrt(a^2-kx^2)
            double kx2=kx*kx;
            double exp_a2_kx2=exp_a2-kx2;
            double sqrt_exp_a2_kx2=Math.sqrt(exp_a2_kx2);
            
            double dzcorr_dkx=reff*exp_a/(exp_a2_kx2*sqrt_exp_a2_kx2); //(exp_a2-kx2)/Math.sqrt(exp_a2-kx2);
            double dzcorr_dar4=dzcorr_dkx*dkx_dar4; // d_tilt
//d_zcorr/d_tilt=d_kx/d_tilt*rc*a/(a^2-kx^2)/sqrt(a^2-kx^2)
            double dzcorr_dar1=-z_corr; // Nice!
// d_reff/dar2==reff; d_reff/dar3=-reff   
// d_zcorr/dar2=z_corr d_zcorr/dar3=-z_corr verified
//z_corr: kx,r_eff,a - ar[1], ar[4], ar[2], ar[3]
            double sqr_azr=Math.sqrt(exp_a2*z_corr2+reff2);

//            double dfcorr_da=exp_a*z_corr2/sqr_azr;
            double dfcorr_dreff=reff/sqr_azr-1;
            double dfcorr_dkx=-z_corr;
         // dfcorr_dar1==0
            double dfcorr_dar2=dfcorr_dreff*exp_r*exp_mk;
//                 dfcorr_dar3=dfcorr_dar2
            double dfcorr_dar4=dfcorr_dkx*dkx_dar4;
            

            
//            double z=z_in-ar[0]-z_corr;
//            double sqrt=Math.sqrt(exp_a2*z*z + reff2);
//            f=sqrt((a*(zin-z0-z_corr))^2 + (r0*(exp(-k))^2)+r0*(1-exp(-k))+ kx*(zin-z0-z_corr) -f_corr {+...aN*(zin-z0-z_corr)^N} - {} are not likely to be ever used 
            
            
            // Dependent on z:
            double z2=z*z;
            double [] df_da=new double[this.modelParams.length]; // last element - derivative for dz
            // derivative for z0 (shift) - ar[0]
            df_da[0]=-1.0/sqrt*exp_a2*z;
            df_da[0]-=kx;
            zp=z;
            for (int i=5;i<this.modelParams.length;i++){
                df_da[0]-=ar[i]*(i-3)*zp; // ar[i] calculated coefficients for current radius
                zp*=z;
            }
            double df_dz_corr=df_da[0];
//            double df_df_corr=-1;//, so subtract each of dfcorr_dar* from df_dar*
            
//          double z=z_in-ar[0]-z_corr;
//          double sqrt=Math.sqrt(exp_a2*z*z + reff2);
//          f=sqrt((a*(zin-z0-z_corr))^2 + (r0*(exp(-k))^2) +r0*(1-exp(-k)) + kx*(zin-z0-z_corr) -f_corr{+...aN*(zin-z0-z_corr)^N} - {} are not likely to be ever used 

            
        // derivative for a (related to numeric aperture) - ar[1]
//            df_da[1]=1.0/sqrt*exp_a*z*z  *exp_a; // d(f)/d(exp_a) *exp_a
// first - calculate derivatives w/O f_corr, z_corr - then apply them            
            df_da[1]=1.0/sqrt*exp_a2*z2; // d(f)/d(exp_a) *exp_a
       // derivative for a (related to lowest PSF radius) - ar[2]
            df_da[2]=(1.0/sqrt*reff*exp_mk + (1-exp_mk)) * exp_r; // d(f)/d(exp_r) *exp_r
       // derivative for k (ar[3]
            df_da[3]=1.0/sqrt*reff*exp_r*exp_mk*(-1) + exp_r*exp_mk;
       // derivative for tilt (ar[4])
            df_da[4]=z*dkx_dar4;
            // derivatives for rest (polynomial) coefficients, probably not ever needed
            zp=z;
            for (int i=5;i<this.modelParams.length;i++){
                zp*=z;
                df_da[i]=zp;
            }
         // new extra term dependent on ar[1] f=...+ kx*(zin-z0-z_corr) +...            
            df_da[1]+=dkx_dar1*z;
         // now apply corrections for z_corr and f_corr
            df_da[1]+=df_dz_corr*dzcorr_dar1; // bad 
            df_da[2]+=df_dz_corr*z_corr;      // good
            df_da[3]+=df_dz_corr*(-z_corr);   // good 
            df_da[4]+=df_dz_corr*dzcorr_dar4; // good
            
            df_da[2]-=dfcorr_dar2; // good
            df_da[3]+=dfcorr_dar2; // good
            df_da[4]-=dfcorr_dar4; // good
            
            // derivative for z (to be combined with mechanical is just negative of derivative for z0, no need to calcualate separately
            // calculate even powers of radius
            double [] dar= new double [this.modelParams[0].length];
            dar[0]=1;
            for (int j=1;j<dar.length;j++){
                dar[j]=dar[j-1]*r;
                if (j==1) dar[j]*=r; // 0,2,3,4,5...
            }
            int index=0;
            for (int i=0;i<this.modelParams.length;i++) for (int j=0;j<this.modelParams[0].length;j++){
                deriv[index++]=df_da[i]*dar[j];
            }
//            calculate d(f)/d(r)
            double df_dr=0.0;
            for (int i=0;i<this.modelParams.length;i++){
                double da_dr=0.0;
                double rp=1.0;
                for (int j=1; j<this.modelParams[0].length;j++) {
                    rp*=r;
                    da_dr+=this.modelParams[i][j]*(j+1)*rp;
                }
                df_dr+=df_da[i]*da_dr;
            }
            deriv[this.modelParams.length*this.modelParams[0].length]=df_dr; // last element
            return f; // function value
        }
        public String getRadialName(int i){
            return "ar_"+(i+1); //TODO: chnage ar_1-> ar_0 (or ar_c),but that will break configuration files
        }
        public String getRadialDecription(int i){
            if (i==0) return "Radial constant coefficient";
            return "Radial polynomial coefficient for r^"+(i+1);
        }
        public String getZName(int i){
            if (i==0) return "z0";
            if (i==1) return "ln(na)";
            if (i==2) return "ln(r0)";
            if (i==3) return "ln(k)";
            if (i==4) return "tilt";
            else return "az_"+(i-3);
        }
        public String getZDescription(int i){
            if (i==0) return "Focal shift";
            if (i==1) return "Defocus/focus shift (~NA), ln()";
            if (i==2) return "Best PSF radius, ln()";
            if (i==3) return "Cross shift, ln()";
            if (i==4) return "Tilt (asymmetry)";
            else return "Polynomial coefficient for z^"+(i-3);
        }
        
        public boolean[] getDefaultMask(){
         boolean [] mask = new boolean[this.modelParams.length*this.modelParams[0].length];
        int index=0;
            for (int i=0;i<this.modelParams.length;i++) for (int j=0;j<this.modelParams[0].length;j++){
                mask[index++]= (i<4) && (j<2);
            }
            return mask;
        }
        public boolean[] maskAllDisabled(){
            boolean [] mask = new boolean[this.modelParams.length*this.modelParams[0].length];
            for (int i=0;i<mask.length;i++) mask[i]=false;
            return mask;
        }
        public boolean[] maskSetDialog(String title, boolean detailed, boolean [] currentMask){
         GenericDialog gd = new GenericDialog(title);
         boolean [] mask = new boolean[this.modelParams.length*this.modelParams[0].length];
         if (currentMask==null) currentMask=getDefaultMask();
         for (int i=0;i<mask.length;i++) {
             mask[i]=currentMask[i];
         }
        boolean [] zMask=new boolean [this.modelParams.length];
        boolean [] rMask=new boolean [this.modelParams[0].length];
            for (int i=0;i<zMask.length;i++) zMask[i]=false;
            for (int i=0;i<rMask.length;i++) rMask[i]=false;

         if (detailed){
             int index=0;
                for (int i=0;i<this.modelParams.length;i++) for (int j=0;j<this.modelParams[0].length;j++){
                    gd.addCheckbox(getZDescription(i)+", "+getRadialDecription(j), mask[index++]);                    
                }
         } else {
             int index=0;
                for (int i=0;i<this.modelParams.length;i++) for (int j=0;j<this.modelParams[0].length;j++){
                    zMask[i] |=mask[index];
                    rMask[j] |=mask[index];
                    index++;
             }            
             
             gd.addMessage("===== Focal parameters =====");
                for (int i=0;i<this.modelParams.length;i++){
                    gd.addCheckbox(getZDescription(i), zMask[i]);                    
                }
             gd.addMessage("===== Radial dependence parameters =====");
                for (int j=0;j<this.modelParams[0].length;j++){
                    gd.addCheckbox(getRadialDecription(j), rMask[j]);                    
                }

         }
         gd.enableYesNoCancel("Apply","Keep"); // default OK (on enter) - "Apply"
     WindowTools.addScrollBars(gd);
    gd.showDialog();
         if (gd.wasCanceled()) return null;
         if (gd.wasOKed()) { // selected non-default "Apply"
             if (detailed){
                 int index=0;
                    for (int i=0;i<this.modelParams.length;i++) for (int j=0;j<this.modelParams[0].length;j++){
                        mask[index++]=gd.getNextBoolean();
                    }
             } else {
                    for (int i=0;i<zMask.length;i++){
                        zMask[i]=gd.getNextBoolean();
                    }
                    for (int i=0;i<rMask.length;i++){
                        rMask[i]=gd.getNextBoolean();
                    }
                 int index=0;
                    for (int i=0;i<this.modelParams.length;i++) for (int j=0;j<this.modelParams[0].length;j++){
                        mask[index++]= zMask[i] && rMask[j];
                 }            
             }             
         }
            return mask;
        }
        
        public boolean showModifyParameterValues(String title, boolean showDisabled, boolean [] mask, boolean isMaster){
         GenericDialog gd = new GenericDialog(title);
        int index=0;
            for (int i=0;i<this.modelParams.length;i++) for (int j=0;j<this.modelParams[0].length;j++){
                String name=getZDescription(i)+", "+getRadialDecription(j);
                boolean dependent=!isMaster && (j==0);
                if (((mask==null) || mask[index] ) && !dependent) {
                    gd.addNumericField(name,this.modelParams[i][j],5,8,"");
                } else if (showDisabled) {
                    if (dependent) {
                        gd.addNumericField("(from master) "+name,this.modelParams[i][j],5,8,"");
                    } else {
                        gd.addNumericField("(disabled) "+name,this.modelParams[i][j],5,8,"");
                    }
//                    gd.addMessage(name +": "+this.modelParams[i][j]);
                }
                index++;
            }
         gd.enableYesNoCancel("Apply","Keep"); // default OK (on enter) - "Apply"
     WindowTools.addScrollBars(gd);
    gd.showDialog();
         if (gd.wasCanceled()) return false;
         if (gd.wasOKed()) { // selected default "Apply"
             index=0;
                for (int i=0;i<this.modelParams.length;i++) for (int j=0;j<this.modelParams[0].length;j++){
                    if ((mask==null) || mask[index] || showDisabled) {
                        this.modelParams[i][j]=gd.getNextNumber();
                    }
                }
                index++;
         }         
         return true;    
        }
    }
    
    public boolean getStrategy(int strategyIndex){
    	FieldStrategies fs=fieldFitting.fieldStrategies;
    	if ((strategyIndex>=0) && (strategyIndex<fs.getNumStrategies())) {
    		fs.getFromStrategy( strategyIndex, fieldFitting);
    		this.strategyComment=fs.getComment(strategyIndex);
    		this.lambda=fs.getInitialLambda(strategyIndex);
    		this.lastInSeries=fs.isStopAfterThis(strategyIndex);
    		this.keepCorrectionParameters=!fs.isResetCorrection(strategyIndex);
    		this.resetCenter=fs.isResetCenter(strategyIndex);
    		return true;
    	} else return false;
    }

    public int organizeStrategies(String title){
    	String [] actions={
    			"<select action>",                      // 0
    			"Restore strategy",                     // 1
    			"Save (replace) strategy",              // 2
    			"Save (insert before/append) strategy", // 3
    			"Remove strategy",                      // 4
    			"Edit strategy (restore-edit-save)"};   // 5
    	FieldStrategies fs=fieldFitting.fieldStrategies;
    	int selectedActionIndex=0;
    	int selectedStrategyIndex=fs.getNumStrategies();
    	boolean editStrategy=false;
        GenericDialog gd = new GenericDialog(title);
        gd.addMessage("Current strategies:");
        String [] indices=new String[fs.getNumStrategies()+1];
        for (int i=0;i<fs.getNumStrategies();i++) {
        	indices[i]=i+": "+fs.getComment(i)+
        			" ("+(fs.isStopAfterThis(i)?"STOP":"CONTINUE")+
        			(fs.isResetCenter(i)?", RESET CENTER":"")+
        			(fs.isResetCorrection(i)?", RESET CORRECTIONS":"")+
        			")";
        }
        indices[fs.getNumStrategies()]="very end";
        for (int i=0;i<fs.getNumStrategies();i++){
            gd.addMessage(i+": "+fs.getComment(i)+
            		" ("+(fs.isStopAfterThis(i)?"STOP":"CONTINUE")+
        			(fs.isResetCenter(i)?", RESET CENTER":"")+
            		(fs.isResetCorrection(i)?", RESET CORRECTIONS":"")+
            		")");
        }
        gd.addMessage("=======================");
        gd.addChoice("Action",actions,actions[selectedActionIndex]);
        gd.addChoice("Index", indices,indices[selectedStrategyIndex]);
        gd.addMessage("=======================");
        gd.addCheckbox("Edit strategy",editStrategy);
        
        gd.enableYesNoCancel("OK","Done");
        WindowTools.addScrollBars(gd);
        gd.showDialog();
        if (gd.wasCanceled()) return -1;
        selectedActionIndex=gd.getNextChoiceIndex();
        selectedStrategyIndex=gd.getNextChoiceIndex();
        editStrategy=gd.getNextBoolean();
        if ((selectedActionIndex!=3) && (selectedStrategyIndex>=fs.getNumStrategies())){
        	selectedStrategyIndex=fs.getNumStrategies()-1; // last
        }
        if (selectedStrategyIndex>=0){
        	switch (selectedActionIndex){
        	case 1:
        		getStrategy(selectedStrategyIndex);
        		if (editStrategy){
        	        if (!fieldFitting.maskSetDialog("Setup restored strategy "+selectedStrategyIndex)) break;
        		}
        		break;
        	case 3:
        		if (editStrategy){
        	        if (!fieldFitting.maskSetDialog("Setup strategy "+selectedStrategyIndex)) break;
        		}
        		if (selectedStrategyIndex>=fs.getNumStrategies()) fs.addStrategy();
        		else fs.insertStrategy(selectedStrategyIndex);
        		// fall through to the next case
        	case 2:
        		if ((selectedActionIndex!=3) && editStrategy){
        	        if (!fieldFitting.maskSetDialog("Setup strategy "+selectedStrategyIndex)) break;
        		}
        		fs.setStrategy(selectedStrategyIndex,fieldFitting);
        		fs.setComment(selectedStrategyIndex,this.strategyComment);
        		fs.setInitialLambda(selectedStrategyIndex,this.lambda);
        		fs.setStopAfterThis(selectedStrategyIndex,this.lastInSeries);
        		fs.setResetCorrection(selectedStrategyIndex,!this.keepCorrectionParameters);
        		fs.setResetCenter(selectedStrategyIndex,this.resetCenter);
        		break;
        	case 4:
        		fs.removeStrategy(selectedStrategyIndex);
        		break;
        	case 0:
        		if (editStrategy){
        	        if (!fieldFitting.maskSetDialog("Setup current strategy")) break;
        		}
        		break;
        	case 5:
        		getStrategy(selectedStrategyIndex);
        		if (!fieldFitting.maskSetDialog("Edit strategy "+selectedStrategyIndex)) break;
        		fs.setStrategy(selectedStrategyIndex,fieldFitting);
        		fs.setComment(selectedStrategyIndex,this.strategyComment);
        		fs.setInitialLambda(selectedStrategyIndex,this.lambda);
        		fs.setStopAfterThis(selectedStrategyIndex,this.lastInSeries);
        		fs.setResetCorrection(selectedStrategyIndex,!this.keepCorrectionParameters);
        		fs.setResetCenter(selectedStrategyIndex,this.resetCenter);
        		break;
        	}
        	
        }
        if (gd.wasOKed()) return 0;
        return 1; // "Done" selected
    }
    
    
    
    
    public class FieldStrategies{
    	ArrayList<FieldSrategy> strategies=new ArrayList<FieldSrategy>();
        public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"strategies_length",strategies.size()+"");
			for (int i=0;i<strategies.size();i++){
				FieldSrategy strategy=strategies.get(i);
				if (strategy!=null) strategy.setProperties(prefix+"strategy_"+i+"_",properties);
			}
        }
    	
        public void getProperties(String prefix,Properties properties){
        	strategies=new ArrayList<FieldSrategy>();
        	String s=properties.getProperty(prefix+"strategies_length");
        	if (s!=null) {
        		int len=Integer.parseInt(s);
        		for (int i=0;i<len;i++){
    				FieldSrategy strategy=new FieldSrategy();
    				// compatibility with old version
    				if (properties.getProperty(prefix+"strategy_"+i+"_"+"centerSelect")!=null){
    					if (debugLevel>0) System.out.println("Restoring new format strategy #"+i);
        				strategy.getProperties(prefix+"strategy_"+i+"_", properties);
    				} else if (properties.getProperty(prefix+"_"+i+"_"+"centerSelect")!=null){
    					if (debugLevel>0) System.out.println("Restoring old format strategy #"+i);
        				strategy.getProperties(prefix+"_"+i+"_", properties);
    				} else {
    					if (debugLevel>0) System.out.println("No info for the field LMA strategy #"+i);
    				}
    				strategies.add(strategy);
        		}
        	}
        }
        public int getNumStrategies(){
        	return strategies.size();
        }

        public void getFromStrategy( // any of the arguments can be null - do not set this array
    			int strategyIndex,
    			FieldFitting fieldFitting){
    		strategies.get(strategyIndex).getFromStrategy( // any of the arguments can be null - do not set this array
    				fieldFitting.centerSelect,
    				fieldFitting.channelSelect,
    				fieldFitting.mechanicalSelect,
    				fieldFitting.curvatureSelect,
    				fieldFitting.sampleCorrSelect,
    				fieldFitting.sampleCorrCost,
    				fieldFitting.sampleCorrSigma,
    				fieldFitting.sampleCorrPullZero
        			);
    	}

        public void setStrategy( // any of the arguments can be null - do not set this array
    			int strategyIndex,
    			FieldFitting fieldFitting){
    		strategies.get(strategyIndex).setStrategy( // any of the arguments can be null - do not set this array
    				fieldFitting.centerSelect,
    				fieldFitting.channelSelect,
    				fieldFitting.mechanicalSelect,
    				fieldFitting.curvatureSelect,
    				fieldFitting.sampleCorrSelect,
    				fieldFitting.sampleCorrCost,
    				fieldFitting.sampleCorrSigma,
    				fieldFitting.sampleCorrPullZero
        			);
    	}
		public double getInitialLambda(
				int strategyIndex) {
			return strategies.get(strategyIndex).getInitialLambda();
		}
		public void setInitialLambda(
				int strategyIndex,
				double initialLambda) {
			strategies.get(strategyIndex).setInitialLambda(initialLambda);
		}
		public boolean isStopAfterThis(
				int strategyIndex) {
			return strategies.get(strategyIndex).isStopAfterThis();
		}
		public boolean isResetCorrection(
				int strategyIndex) {
			return strategies.get(strategyIndex).isResetCorrection();
		}

		public boolean isResetCenter(
				int strategyIndex) {
			return strategies.get(strategyIndex).isResetCenter();
		}
				
		public boolean isLast(
				int strategyIndex) {
			if (strategyIndex < 0) return true;
			if (strategyIndex>=(getNumStrategies()-1)) return true;// last
			return isStopAfterThis(strategyIndex);
		}
		
		public void setStopAfterThis(
				int strategyIndex,
				boolean stopAfterThis) {
			strategies.get(strategyIndex).setStopAfterThis(stopAfterThis);
		}
		public void setResetCorrection(
				int strategyIndex,
				boolean resetCorrection) {
			strategies.get(strategyIndex).setResetCorrection(resetCorrection);
		}
		public void setResetCenter(
				int strategyIndex,
				boolean resetCenter) {
			strategies.get(strategyIndex).setResetCenter(resetCenter);
		}
		public String getComment(
				int strategyIndex){
			return strategies.get(strategyIndex).getComment();
		}

		public void setComment(
				int strategyIndex,
				String comment){
			strategies.get(strategyIndex).setComment(comment);
		}
		
		
        public void insertStrategy(
        		int strategyIndex){
        	strategies.add(strategyIndex,new FieldSrategy());
        }
        public void removeStrategy(
        		int strategyIndex){
        	strategies.remove(strategyIndex);
        }
        public void addStrategy(){
        	strategies.add(new FieldSrategy());
        }

        public void saveStrategies(
        		String path,      // full path w/o extension or null
        		String directory ){
        	String [] patterns= {".fstg-xml",".xml"};
        	if (path==null) {
        		path= CalibrationFileManagement.selectFile(true, // save  
        				"Save Field LMA Strategy selection", // title
        				"Select Field LMA Strategy file", // button
        				new CalibrationFileManagement.MultipleExtensionsFileFilter(patterns,
        						"Strategy files (*.fstg-xml)"), // filter
        						directory); // may be ""
        	} else path+=patterns[0];
        	if (path==null) return;
        	Properties properties=new Properties();
        	setProperties("",properties); // no prefix
        	OutputStream os;
        	try {
        		os = new FileOutputStream(path);
        	} catch (FileNotFoundException e1) {
        		IJ.showMessage("Error","Failed to open field LMA strategy file for writing: "+path);
        		return;
        	}
        	try {
        		properties.storeToXML(os,
        				"last updated " + new java.util.Date(), "UTF8");

        	} catch (IOException e) {
        		IJ.showMessage("Error","Failed to write XML configuration file: "+path);
        		return;
        	}
        	try {
        		os.close();
        	} catch (IOException e) {
        		// TODO Auto-generated catch block
        		e.printStackTrace();
        	}
        	if (debugLevel>0) System.out.println("Field LMA strategy parameters are saved to "+path);
        }

        public void loadStrategies(
        		String path,      // full path w/o extension or null
        		String directory ){
        	String [] patterns= {".fstg-xml",".xml"};
        	if (path==null) {
        		path= CalibrationFileManagement.selectFile(false, // save  
        				"Field LMA Strategy selection", // title
        				"Select Field LMA Strategy file", // button
        				new CalibrationFileManagement.MultipleExtensionsFileFilter(patterns,
        						"Strategy files (*.fstg-xml)"), // filter
        						directory); // may be ""
        	} else {
  	    	  // do not add extension if it already exists
  	    	  if ((path.length()<patterns[0].length()) || (!path.substring(path.length()-patterns[0].length()).equals(patterns[0]))){
  	    		  path+=patterns[0];
  	    	  }
        	}
        	if (path==null) return;
        	InputStream is;
        	try {
        		is = new FileInputStream(path);
        	} catch (FileNotFoundException e) {
        		IJ.showMessage("Error","Failed to open field LMA strategy file: "+path);
        		return;
        	}
        	Properties properties=new Properties();
        	try {
        		properties.loadFromXML(is);

        	} catch (IOException e) {
        		IJ.showMessage("Error","Failed to read field LMA strategy file: "+path);
        		return;
        	}
        	try {
        		is.close();
        	} catch (IOException e) {
        		// TODO Auto-generated catch block
        		e.printStackTrace();
        	}      
        	getProperties("",properties); // no prefix
        	if (debugLevel>0) System.out.println("Field LMA strategy parameters are restored from "+path);
        }
        
    	class FieldSrategy{
    		private boolean [] centerSelect=null;
        	private boolean [] channelSelect=null;
        	private boolean [] mechanicalSelect=null;
        	private boolean [][] curvatureSelect=new boolean[6][];
        	private boolean [][] sampleCorrSelect= new boolean[6][]; // enable individual (per sample coordinates) correction of parameters
        	private double [][] sampleCorrCost= new double[6][]; // equivalent cost of one unit of parameter value (in result units, um)
        	private double [][] sampleCorrSigma= new double[6][]; // sigma (in mm) for neighbors influence
        	private double [][] sampleCorrPullZero=new double[6][]; // 1.0 - only difference from neighbors matters, 0.0 - only difference from 0
        	// TODO: add LMA-specific (initial lambda, stop after this
        	private double initialLambda=0.001;
			private boolean stopAfterThis=true;
			private boolean resetCorrection=false;
			private boolean resetCenter=false;
			private String strategyComment="";
			private Properties properties=null;
			private String prefix=null;
			public FieldSrategy(){
				setDefaults();
			}
			public FieldSrategy(
					String strategyComment,
					double lambda,
					boolean lastInSeries,
					boolean resetCorrection,
					boolean resetCenter
					){
				this.strategyComment=strategyComment;
				initialLambda=lambda;
				stopAfterThis=lastInSeries;
				this.resetCorrection=resetCorrection;
				this.resetCenter=resetCenter;
				setDefaults();
			}
			private void setDefaults(){
				boolean [][] booleanNull6={null,null,null,null,null,null};
				double [][]  doubleNull6= {null,null,null,null,null,null};
				curvatureSelect=    booleanNull6.clone();
				sampleCorrSelect=   booleanNull6.clone();
				sampleCorrCost=     doubleNull6.clone();
				sampleCorrSigma=    doubleNull6.clone();
				sampleCorrPullZero= doubleNull6.clone();
			}
			private String boolToStr(boolean [] ba){
				if (ba==null) return "";
				String result="";
				for (boolean b : ba) result+=b?"+":"-";
				return result;
			}
			private void setPropBool(boolean [] arr, String name){
            	if (arr!=null) properties.setProperty(prefix+name,boolToStr(arr));
        	}
			private void setPropBool(boolean [][] arr, String name){
				if (arr!=null) {
					properties.setProperty(prefix+name+"_length",arr.length+"");
					for (int i=0;i<arr.length;i++) {
						setPropBool(arr[i], name+"_"+i);
					}
				}
			}
			private void setPropDouble(double [] arr, String name){
            	if (arr!=null) properties.setProperty(prefix+name,doubleToStr(arr));
        	}
			private void setPropDouble(double [][] arr, String name){
				if (arr!=null) {
					properties.setProperty(prefix+name+"_length",arr.length+"");
					for (int i=0;i<arr.length;i++) {
						setPropDouble(arr[i], name+"_"+i);
					}
				}
			}
			private boolean [] strToBool(String s){
				if (s==null) return new boolean [0];
				boolean [] result= new boolean [s.length()];
				for (int i=0;i<result.length;i++) result[i]=(s.charAt(i)=='+');
				return result;
			}
			private String doubleToStr(double [] da){
				if (da==null) return "";
				String result="";
				for (double d : da) result+=","+d;
				if (result.length()>0) result=result.substring(1);
				return result;
			}
			private double [] strToDouble(String s){
				if (s==null) return new double [0];
				String [] sa=s.split(",");
				double [] result= new double [sa.length];
				for (int i=0;i<result.length;i++) result[i]=Double.parseDouble(sa[i]);
				return result;
			}
			
			private boolean [] getPropBool(boolean[] arr, String name){
				String s=properties.getProperty(prefix+name);
				if (s!=null) return strToBool(s);
				else return arr;
        	}

			private boolean[][] getPropBool(boolean[][] arr, String name){
				if (arr==null){
					String s=properties.getProperty(prefix+name+"_length");
					if (s==null) return null;
					arr=new boolean[Integer.parseInt(s)][];
					for (int i=0;i<arr.length;i++) arr[i]=null;
				}
				for (int i=0;i<arr.length;i++){
					boolean [] a=getPropBool((boolean[]) null,name+"_"+i);
					if (a!=null) arr[i]=a;
				}
				return arr;
        	}

			private double [] getPropDouble(double [] arr, String name){
				String s=properties.getProperty(prefix+name);
				if (s!=null) return strToDouble(s);
				else return arr;
        	}
			
			private double[][] getPropDouble(double[][] arr, String name){
				if (arr==null){
					String s=properties.getProperty(prefix+name+"_length");
					if (s==null) return null;
					arr=new double[Integer.parseInt(s)][];
					for (int i=0;i<arr.length;i++) arr[i]=null;
				}
				for (int i=0;i<arr.length;i++){
					double [] a=getPropDouble((double []) null,name+"_"+i);
					if (a!=null) arr[i]=a;
				}
				return arr;
        	}

			public String getComment(){
				return strategyComment;
			}
			public void setComment(String comment){
				strategyComment=comment;
			}

			public double getInitialLambda() {
				return initialLambda;
			}
			public void setInitialLambda(double initialLambda) {
				this.initialLambda = initialLambda;
			}
			public boolean isStopAfterThis() {
				return stopAfterThis;
			}
			public void setStopAfterThis(boolean stopAfterThis) {
				this.stopAfterThis = stopAfterThis;
			}
			public boolean isResetCorrection() {
				return resetCorrection;
			}
			public void setResetCorrection(boolean resetCorrection) {
				this.resetCorrection = resetCorrection;
			}
			public boolean isResetCenter() {
				return resetCenter;
			}
			public void setResetCenter(boolean resetCenter) {
				this.resetCenter = resetCenter;
			}

			public void setStrategy( // any of the arguments can be null - do not set this array
            		boolean [] centerSelect,
                	boolean [] channelSelect,
                	boolean [] mechanicalSelect,
                	boolean [][] curvatureSelect,
                	boolean [][] sampleCorrSelect,
                	double [][] sampleCorrCost,
                	double [][] sampleCorrSigma,
                	double [][] sampleCorrPullZero
        			){
        		if (centerSelect!=null)    this.centerSelect=      centerSelect.clone();
        		if (channelSelect!=null)    this.channelSelect=      channelSelect.clone();
        		if (mechanicalSelect!=null) this.mechanicalSelect=   mechanicalSelect.clone();
        		if (curvatureSelect!=null) {
        			this.curvatureSelect=    new boolean[curvatureSelect.length][];
        			for (int i=0;i<curvatureSelect.length;i++)    this.curvatureSelect[i]=curvatureSelect[i].clone();
        		}
        		if (sampleCorrSelect!=null) {
        			this.sampleCorrSelect=   new boolean[sampleCorrSelect.length][];
        			for (int i=0;i<sampleCorrSelect.length;i++)   this.sampleCorrSelect[i]=sampleCorrSelect[i].clone();
        		}
        		if (sampleCorrCost!=null) {
        			this.sampleCorrCost=     new double[sampleCorrCost.length][];
        			for (int i=0;i<sampleCorrCost.length;i++)     this.sampleCorrCost[i]=sampleCorrCost[i].clone();
        		}
        		if (sampleCorrSigma!=null) {
        			this.sampleCorrSigma=    new double[sampleCorrSigma.length][];
        			for (int i=0;i<sampleCorrSigma.length;i++)    this.sampleCorrSigma[i]=sampleCorrSigma[i].clone();
        		}
        		if (sampleCorrPullZero!=null) {
        			this.sampleCorrPullZero= new double[sampleCorrPullZero.length][];
        			for (int i=0;i<sampleCorrPullZero.length;i++) this.sampleCorrPullZero[i]=sampleCorrPullZero[i].clone();
        		}
        	}
        	public void getFromStrategy( // any of the arguments can be null - do not set this array
            		boolean [] centerSelect,
                	boolean [] channelSelect,
                	boolean [] mechanicalSelect,
                	boolean [][] curvatureSelect,
                	boolean [][] sampleCorrSelect,
                	double [][] sampleCorrCost,
                	double [][] sampleCorrSigma,
                	double [][] sampleCorrPullZero
        			){
        		if (centerSelect!=null) for (int i=0;i<centerSelect.length;i++) centerSelect[i]=this.centerSelect[i];
        		if (channelSelect!=null) for (int i=0;i<channelSelect.length;i++) channelSelect[i]=this.channelSelect[i];
        		if (mechanicalSelect!=null) for (int i=0;i<mechanicalSelect.length;i++) mechanicalSelect[i]=this.mechanicalSelect[i];
        		if (curvatureSelect!=null) {
        			for (int i=0;i<curvatureSelect.length;i++)    curvatureSelect[i]=this.curvatureSelect[i].clone();
        		}
        		if (sampleCorrSelect!=null) {
        			for (int i=0;i<sampleCorrSelect.length;i++)   sampleCorrSelect[i]=this.sampleCorrSelect[i].clone();
        		}
        		if (sampleCorrCost!=null) {
        			for (int i=0;i<sampleCorrCost.length;i++)     sampleCorrCost[i]=this.sampleCorrCost[i].clone();
        		}
        		if (sampleCorrSigma!=null) {
        			for (int i=0;i<sampleCorrSigma.length;i++)    sampleCorrSigma[i]=this.sampleCorrSigma[i].clone();
        		}
        		if (sampleCorrPullZero!=null) {
        			for (int i=0;i<sampleCorrPullZero.length;i++) sampleCorrPullZero[i]=this.sampleCorrPullZero[i].clone();
        		}
        	}
            public void setProperties(String prefix,Properties properties){
            	this.prefix=prefix;
            	this.properties=properties;
            	setPropBool(centerSelect, "centerSelect");
            	setPropBool(channelSelect, "channelSelect");
            	setPropBool(mechanicalSelect, "mechanicalSelect");
            	setPropBool(curvatureSelect, "curvatureSelect");
            	setPropBool(sampleCorrSelect, "sampleCorrSelect");
            	setPropDouble(sampleCorrCost, "sampleCorrCost");
            	setPropDouble(sampleCorrSigma, "sampleCorrSigma");
            	setPropDouble(sampleCorrPullZero, "sampleCorrPullZero");
				properties.setProperty(prefix+"initialLambda",getInitialLambda()+"");
				properties.setProperty(prefix+"stopAfterThis",isStopAfterThis()+"");
				properties.setProperty(prefix+"resetCorrection",isResetCorrection()+"");
				properties.setProperty(prefix+"resetCenter",isResetCenter()+"");
				properties.setProperty(prefix+"strategyComment","<![CDATA["+strategyComment+ "]]>");
            }
            public void getProperties(String prefix,Properties properties){
            	this.prefix=prefix;
            	this.properties=properties;
            	centerSelect=      getPropBool(centerSelect, "centerSelect");
            	channelSelect=     getPropBool(channelSelect, "channelSelect");
            	mechanicalSelect=  getPropBool(mechanicalSelect, "mechanicalSelect");
            	curvatureSelect=   getPropBool(curvatureSelect, "curvatureSelect");
            	sampleCorrSelect=  getPropBool(sampleCorrSelect, "sampleCorrSelect");
            	sampleCorrCost=    getPropDouble(sampleCorrCost, "sampleCorrCost");
            	sampleCorrSigma=   getPropDouble(sampleCorrSigma, "sampleCorrSigma");
            	sampleCorrPullZero=getPropDouble(sampleCorrPullZero, "sampleCorrPullZero");
            	String s=properties.getProperty(prefix+"initialLambda");
            	if (s!=null) initialLambda=Double.parseDouble(s);
            	s=properties.getProperty(prefix+"stopAfterThis");
            	if (s!=null) stopAfterThis=Boolean.parseBoolean(s);
            	s=properties.getProperty(prefix+"resetCorrection");
            	if (s!=null) resetCorrection=Boolean.parseBoolean(s);
            	s=properties.getProperty(prefix+"resetCenter");
            	if (s!=null) resetCenter=Boolean.parseBoolean(s);
            	s=properties.getProperty(prefix+"strategyComment");
            	if (s!=null){
            		strategyComment=s;
        			if ((strategyComment.length()>10) && strategyComment.substring(0,9).equals("<![CDATA[")) {
        				strategyComment=strategyComment.substring(9,strategyComment.length()-3); 
        			}
            	}
            }
    	}
    }
}



