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
import java.io.FileWriter;
import java.io.IOException;
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
	double  filterInpuFarRatio; // remove samples that are farther than this ration of average distance
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
	private boolean rslt_show_z_individual;
	private boolean rslt_show_f_axial;
	private boolean rslt_show_f_individual;
	private double rslt_scan_below;
	private double rslt_scan_above;
	private double rslt_scan_step;
	private boolean rslt_mtf50_mode;
	private boolean [] rslt_show_chn;

// not saved/restored
	
	private double lambdaStepUp; // multiply lambda by this if result is worse
	private double lambdaStepDown; // multiply lambda by this if result is better
	private double thresholdFinish; // (copied from series) stop iterations if 2 last steps had less improvement (but not worsening )
	private int numIterations; // maximal number of iterations
	private double maxLambda; // max lambda to fail
	private double lambda; // copied from series
	private boolean stopEachStep; // open dialog after each fitting step
	private boolean stopOnFailure; // open dialog when fitting series failed
	private boolean showParams; // show modified parameters
	private boolean showDisabledParams;
	private boolean showCorrectionParams;
	private boolean keepCorrectionParameters;
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
	private double [] currentfX=null; // array of "f(x)" - simulated data for all images, combining pixel-X and pixel-Y (odd/even)
	private double [] nextfX=null; // array of "f(x)" - simulated data for all images, combining pixel-X and pixel-Y (odd/even)
	private double currentRMS=-1.0; // calculated RMS for the currentVector->currentfX
	private double currentRMSPure=-1.0; // calculated RMS for the currentVector->currentfX
	private double nextRMS=-1.0; // calculated RMS for the nextVector->nextfX
	private double nextRMSPure=-1.0; // calculated RMS for the nextVector->nextfX

	private double firstRMS=-1.0; // RMS before current series of LMA started
	private double firstRMSPure=-1.0; // RMS before current series of LMA started


    public void setDefaults(){
    	sagittalMaster=false; // center data is the same, when true sagittal fitting only may change r=0 coefficients,
    	parallelOnly = true; // only process measurements for parallel moves
    	filterInput = true;
    	filterInputMotorDiff = 500.0;
    	filterInputDiff = 2.0; // um
    	filterInputFirstLast = true; 
    	filterInputTooFar = true; // filter samples that are too far from the "center of mass" of other samples 
    	filterInpuFarRatio = 3.0; // remove samples that are farther than this ration of average distance
    	// when false - tangential is master
    	double [] minMeasDflt= {0.5,0.5,0.5,0.5,0.5,0.5}; // pixels
    	minMeas= minMeasDflt; // pixels
    	double [] maxMeasDflt= {4.5,4.5,4.5,4.5,4.5,4.5}; // pixels
    	maxMeas= maxMeasDflt; // pixels
    	double [] thresholdMaxDflt= {2.4,3.0,2.6,3.0,3.1,3.0}; // pixels
    	thresholdMax= thresholdMaxDflt; // pixels
    	useMinMeas= true;
    	useMaxMeas= true;
    	useThresholdMax=true;
    	weightMode=1; // 0; // 0 - same weight, 1 - linear threshold difference, 2 - quadratic thershold difference
    	weightRadius=0.0; //2.0; // Gaussian sigma in mm
    	k_red=0.7;
    	k_blue=0.4;
    	qb_scan_below=-40.0; // um
    	qb_scan_above= 80.0; // um
    	qb_scan_step= 0.5; // um
    	qb_use_corrected=true;
    	qb_invert=true;
    	z_relative=true;  // focal distance relative to center greeen
    	rslt_show_z_axial=true;
    	rslt_show_z_individual=true;
    	rslt_show_f_axial=true;
    	rslt_show_f_individual=true;
    	rslt_scan_below=-10.0;
    	rslt_scan_above= 10.0;
    	rslt_scan_step= 5.0;
    	rslt_mtf50_mode= true;
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
    	stopOnFailure= true; // open dialog when fitting series failed
    	showParams= false; // show modified parameters
    	showDisabledParams = false;
    	showCorrectionParams = false;
    	keepCorrectionParameters = true;
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

    	
    	
    	
    }
    
    
	public void setProperties(String prefix,Properties properties){
		if (debugLevel>1) System.out.println("FocusingField: setProperties()");
		if (fieldFitting == null) {
			System.out.println("fieldFitting is not initioalized, nothing to save");
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
		properties.setProperty(prefix+"filterInpuFarRatio",filterInpuFarRatio+"");
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
		properties.setProperty(prefix+"rslt_show_z_individual",rslt_show_z_individual+"");
		properties.setProperty(prefix+"rslt_show_f_axial",rslt_show_f_axial+"");
		properties.setProperty(prefix+"rslt_show_f_individual",rslt_show_f_individual+"");
		properties.setProperty(prefix+"rslt_scan_below",rslt_scan_below+"");
		properties.setProperty(prefix+"rslt_scan_above",rslt_scan_above+"");
		properties.setProperty(prefix+"rslt_scan_step",rslt_scan_step+"");
		properties.setProperty(prefix+"rslt_mtf50_mode",rslt_mtf50_mode+"");
		for (int chn=0; chn<rslt_show_chn.length; chn++) properties.setProperty(prefix+"rslt_show_chn_"+chn,rslt_show_chn[chn]+"");
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
		if (properties.getProperty(prefix+"filterInpuFarRatio")!=null)
			filterInpuFarRatio=Double.parseDouble(properties.getProperty(prefix+"filterInpuFarRatio"));
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
		if (properties.getProperty(prefix+"rslt_show_z_individual")!=null)
			rslt_show_z_individual=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_show_z_individual"));
		if (properties.getProperty(prefix+"rslt_show_f_axial")!=null)
			rslt_show_f_axial=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_show_f_axial"));
		if (properties.getProperty(prefix+"rslt_show_f_individual")!=null)
			rslt_show_f_individual=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_show_f_individual"));
		if (properties.getProperty(prefix+"rslt_scan_below")!=null)
			rslt_scan_below=Double.parseDouble(properties.getProperty(prefix+"rslt_scan_below"));
		if (properties.getProperty(prefix+"rslt_scan_above")!=null)
			rslt_scan_above=Double.parseDouble(properties.getProperty(prefix+"rslt_scan_above"));
		if (properties.getProperty(prefix+"rslt_scan_step")!=null)
			rslt_scan_step=Double.parseDouble(properties.getProperty(prefix+"rslt_scan_step"));
		if (properties.getProperty(prefix+"rslt_mtf50_mode")!=null)
			rslt_mtf50_mode=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_mtf50_mode"));
		for (int chn=0; chn<rslt_show_chn.length; chn++) if (properties.getProperty(prefix+"rslt_show_chn_"+chn)!=null)
			rslt_show_chn[chn]=Boolean.parseBoolean(properties.getProperty(prefix+"rslt_show_chn_"+chn));
	}
	public void setDebugLevel(int debugLevel){
		this.debugLevel=debugLevel;
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
    gd.addNumericField("\"Too far\" ratio to the average distance to the center of measurements",filterInpuFarRatio,3,5,"um");
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
    
    parallelOnly=         gd.getNextBoolean();
    filterInput=          gd.getNextBoolean();
    filterInputMotorDiff= gd.getNextNumber();
    filterInputDiff=      gd.getNextNumber();
    
    filterInputFirstLast= gd.getNextBoolean();
    filterInputTooFar=    gd.getNextBoolean();
    filterInpuFarRatio=   gd.getNextNumber();
    
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
        	if (debugLevel>1) System.out.println("configureDataVector(): Applying properties");
        	getProperties(propertiesPrefix,savedProperties);
        }
    }
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

private void maskDataWeights(boolean [] enable){
	for (int i=0;i<enable.length;i++){
		if (!enable[i]) dataWeights[i]=0.0;
	}
}

private boolean [] filterTooFar(double ratio,boolean [] enable_in){
	if (enable_in==null) {
		enable_in=new boolean [dataVector.length];
		for (int i=0;i<enable_in.length;i++)enable_in[i]=true;
	}
	boolean [] enable_out=enable_in.clone();
	int numFiltered = 0;
	double [][][] data=new double [getNumChannels()][getNumSamples()][3];
	double [] z_sample=new double [dataVector.length];
	for (int chn=0;chn<data.length;chn++)
		for (int sample=0;sample<data[chn].length;sample++)
		for (int i=0;i<data[chn][sample].length;i++)data[chn][sample][i]=0;
	for (int index=0;index<dataVector.length;index++) if ((index>=enable_in.length) ||enable_in[index]){
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
	for (int index=0;index<dataVector.length;index++) if ((index>=enable_in.length) ||enable_in[index]){
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
	return enable_out;
}


private boolean [] filterCrazyInput(
		boolean [] enable_in, // [meas][cjn][sample] (or null) // can be shorter or longer than dataVector
		double maxMotDiff,
		double diff,
		boolean removeFirstLast // very first, very last in all samples (or after big move) - OK
		){
	if (enable_in==null) {
		enable_in=new boolean [dataVector.length];
		for (int i=0;i<enable_in.length;i++)enable_in[i]=true;
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
	for (int index=0;index<dataVector.length;index++) { // crazy neighbor still kills even if is ignored itself - needed
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


// includes deselected channels
public void setDataVector(MeasuredSample [] vector){ // remove unused channels if any. vector is already corrected from input data, FWHM psf
    if (debugLevel>1) System.out.println("+++++ (Re)calculating sample weights +++++");
    int [] diffs=null;
    if (parallelOnly) diffs=getParallelDiff(vector);
    boolean [] chanSel=fieldFitting.getSelectedChannels();
    int numSamples=0;
    for (int i=0;i<vector.length;i++) if (chanSel[vector[i].channel]){
    	if ((diffs!=null) && ( 
    					((vector[i].motors[1]-vector[i].motors[0]) != diffs[0]) ||
    					((vector[i].motors[2]-vector[i].motors[0]) != diffs[1]))) continue;
    	numSamples++;
    }
    dataVector=new MeasuredSample [numSamples];
    int n=0;
    for (int i=0;i<vector.length;i++) if (chanSel[vector[i].channel]) {
    	if ((diffs!=null) && ( 
				((vector[i].motors[1]-vector[i].motors[0]) != diffs[0]) ||
				((vector[i].motors[2]-vector[i].motors[0]) != diffs[1]))) continue;
    	dataVector[n++]=vector[i];
    }
    int corrLength=fieldFitting.getNumberOfCorrParameters();
    dataValues = new double [dataVector.length+corrLength];
    dataWeights = new double [dataVector.length+corrLength];
//     sumWeights=0.0;
    int mode=weightMode;
    double kw= (weightRadius>0.0)?(-0.5*getPixelMM()*getPixelMM()/(weightRadius*weightRadius)):0;
//weightRadius     
    if (weightReference==null)mode=0;
    for (int i=0;i<dataVector.length;i++){
        MeasuredSample ms=dataVector[i];
        dataValues[i]=ms.value;
        double diff=weightReference[ms.channel]-ms.value;
        if (diff<0.0) diff=0;
        switch (mode){
        case 0:    dataWeights[i]=1.0; break;
        case 1:    dataWeights[i]=diff; break;
        case 2:    dataWeights[i]=diff*diff; break;
        default: dataWeights[i]=1.0;
        }
        if (weightRadius>0.0){
            double r2=(ms.px-currentPX0)*(ms.px-currentPX0)+(ms.py-currentPY0)*(ms.py-currentPY0);
            dataWeights[i]*=Math.exp(kw*r2);
        }
//         sumWeights+=dataWeights[i];
    }
    for (int i=0;i<corrLength;i++){
        dataValues[i+dataVector.length]=0.0; // correction target is always 0
        dataWeights[i+dataVector.length]=1.0; // improve?
    }
    if (filterInput){
    	boolean [] en=dataWeightsToBoolean();
    	en= filterCrazyInput(
    			en, // [meas][cjn][sample] (or null) // can be shorter or longer than dataVector
    			filterInputMotorDiff,
    			filterInputDiff,
    			filterInputFirstLast
    			);
    	maskDataWeights(en);
    }
    if (filterInputTooFar){
    	boolean [] en=dataWeightsToBoolean();
    	en= filterTooFar(
    			filterInpuFarRatio,
    			en);
    	maskDataWeights(en);
    }
    
}

// for compatibility with Distortions class\

public void commitParameterVector(double [] vector){
    fieldFitting.commitParameterVector(vector,sagittalMaster);
    // recalculate measured S,T (depend on center) if center is among fitted parameters
    boolean [] centerSelect=fieldFitting.getCenterSelect();
    if (centerSelect[0] ||centerSelect[1]){ // do not do that if XC, YC are not modified
        // recalculate data vector
        double [] pXY=fieldFitting.getCenterXY();
    currentPX0=pXY[0];
    currentPY0=pXY[1];
        if (debugLevel>1) System.out.println("Updated currentPX0="+currentPX0+", currentPY0="+currentPY0);
    if (correct_measurement_ST && updateWeightWhileFitting) {
        setDataVector(createDataVector(
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
        boolean updateSelection,
        double centerPX,
        double centerPY,
        double [] minMeas, // pixels
        double [] maxMeas, // pixels
        double [] thresholdMax){ // pixels
    debugDerivatives=debugLevel==3;
    currentPX0=centerPX;
    currentPY0=centerPY;
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
			msg=": initial RMS="+IJ.d2s(this.currentRMS,8)+" (pure RMS"+IJ.d2s(this.currentRMSPure,8)+")"+
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
        gd.addCheckbox("Debug df/dX0, df/dY0", false);
        gd.addNumericField("Debug Jacobian for point number", this.debugPoint, 0, 5,"(-1 - none)");
        gd.addNumericField("Debug Jacobian for parameter number", this.debugParameter, 0, 5,"(-1 - none)");
        
        gd.addCheckbox("Keep current correction parameters (do not reset)", this.keepCorrectionParameters);
//        gd.addNumericField("Iteration number to start (0.."+(numSeries-1)+")", this.seriesNumber, 0);
        gd.addNumericField("Initial LMA Lambda ", this.lambda, 5);
        gd.addNumericField("Multiply lambda on success", this.lambdaStepDown, 5);
        gd.addNumericField("Threshold RMS to exit LMA", this.thresholdFinish, 7,9,"pix");
        gd.addNumericField("Multiply lambda on failure", this.lambdaStepUp, 5);
        gd.addNumericField("Threshold lambda to fail", this.maxLambda, 5);
        gd.addNumericField("Maximal number of iterations", this.numIterations, 0);

        gd.addCheckbox("Dialog after each iteration step", this.stopEachStep);
//        gd.addCheckbox("Dialog after each iteration series", this.stopEachSeries);
        gd.addCheckbox("Dialog after each failure", this.stopOnFailure);
//        gd.addCheckbox("Ask for weight function filter", this.askFilter);
        
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
     this.debugDerivativesFxDxDy=gd.getNextBoolean();
     
		debugPoint=     (int) gd.getNextNumber();
		debugParameter= (int) gd.getNextNumber();

        this.keepCorrectionParameters = gd.getNextBoolean();
//     this.seriesNumber= (int) gd.getNextNumber();
        this.lambda= gd.getNextNumber();
        this.lambdaStepDown= gd.getNextNumber();
        this.thresholdFinish= gd.getNextNumber();
        this.lambdaStepUp= gd.getNextNumber();
        this.maxLambda= gd.getNextNumber();
        this.numIterations= (int) gd.getNextNumber();
        this.stopEachStep= gd.getNextBoolean();
//        this.stopEachSeries= gd.getNextBoolean();
        this.stopOnFailure= gd.getNextBoolean();
//        this.askFilter= gd.getNextBoolean();
        this.showParams= gd.getNextBoolean();
        this.showDisabledParams= gd.getNextBoolean();
        this.showCorrectionParams= gd.getNextBoolean();
//        this.showThisImages= gd.getNextBoolean();
//        this.showNextImages= gd.getNextBoolean();
//        this.threadsMax= (int) gd.getNextNumber();
//        this.threadedLMA= gd.getNextBoolean();
//        resetCorrections= gd.getNextBoolean();
    if (!keepCorrectionParameters) fieldFitting.resetSampleCorr();
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
    
    double [] center_z=fieldFitting.getZCenters();
     GenericDialog gd = new GenericDialog("Setup results table");
    double best_qb_axial= fieldFitting.getBestQualB(
                k_red,
                k_blue,
                false);
    double best_qb_corr= fieldFitting.getBestQualB(
            k_red,
            k_blue,
            true);
    gd.addMessage("Best center focus for Red "+ IJ.d2s(center_z[0],3)+" um");
    gd.addMessage("Best center focus for Green "+ IJ.d2s(center_z[1],3)+" um");
    gd.addMessage("Best center focus for Blue "+ IJ.d2s(center_z[2],3)+" um");
    gd.addMessage("Best composite distance for FWHM^4, axial model "+ IJ.d2s(best_qb_axial,3)+" um");
    gd.addMessage("Best composite distance for FWHM^4, individual "+ IJ.d2s(best_qb_corr,3)+" um");
    for (int i=0;i<rslt_show_chn.length;i++){
    gd.addCheckbox("Show results for "+fieldFitting.getDescription(i), this.rslt_show_chn[i]);
    }
    gd.addCheckbox("Show best focus distance (axial model)", this.rslt_show_z_axial);
    gd.addCheckbox("Show best focus distance (per-sample adjusted)", this.rslt_show_z_individual);
    gd.addCheckbox("Show constant-z sections (axial model)", this.rslt_show_f_axial);
    gd.addCheckbox("Show constant-z sections (per-sample adjusted)", this.rslt_show_f_individual);
    gd.addCheckbox("Show mtf50 (false - PSF FWHM)", this.rslt_mtf50_mode);
    gd.addCheckbox("Show focal distance relative to center green (false - motor absolute)", this.z_relative);
    
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
    this.rslt_show_z_individual= gd.getNextBoolean();
    this.rslt_show_f_axial= gd.getNextBoolean();
    this.rslt_show_f_individual= gd.getNextBoolean();
    this.rslt_mtf50_mode= gd.getNextBoolean();
    this.z_relative= gd.getNextBoolean();
    this.rslt_scan_below= gd.getNextNumber();
    this.rslt_scan_above= gd.getNextNumber();
    this.rslt_scan_step= gd.getNextNumber();

    listCombinedResults(
    		"Field curvature measuremnts results", //String title,
    		null, //String path,
            z_relative?center_z[1]:0.0,
    		rslt_show_chn, //boolean [] show_chn,
    		rslt_show_z_axial, //boolean show_z_axial,
    		rslt_show_z_individual, //boolean show_z_individual,
    		rslt_show_f_axial, //boolean show_f_axial,
    		rslt_show_f_individual, //boolean show_f_individual,
    		center_z[1]+rslt_scan_below, // double scan_below,
    		center_z[1]+rslt_scan_above, //double scan_above,
    		rslt_scan_step, //double scan_step,
    		rslt_mtf50_mode); //boolean freq_mode)
}

public void listCombinedResults(
		String title,
		String path,
		double z0, // subtract from z
		boolean [] show_chn,
		boolean show_z_axial,
		boolean show_z_individual,
		boolean show_f_axial,
		boolean show_f_individual,
		double scan_below,
		double scan_above,
		double scan_step,
		boolean freq_mode){
	String [] chnNames={"RS","RT","GS","GT","BS","BT"};

	double k=this.fwhm_to_mtf50; //TODO: correct psf fwhm to mtf50 conversion
	//     double k=2*Math.log(2.0)/Math.PI*1000;
	//     String header="Z(um)\tComposite\tRed\tGreen\tBlue";
	String header="radius(mm)";
	if (show_z_axial){
		for (int i=0;i<show_chn.length;i++){
			if (show_chn[i])header+="\tZ"+chnNames[i];
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
	double [][][][] f_values=new double [numSect][2][][];
	int sect=0;
	for (double z=scan_below;z<=scan_above;z+=scan_step){
		if (show_f_axial){
			double [][] f=fieldFitting.getCalcValuesForZ(z, false,true);
			double [] f0=fieldFitting.getCalcValuesForZ(z, 0.0);
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
		sect++;
	}
	double [][][] z_values=new double [2][][];
	if (show_z_axial){
		double [][] zai=fieldFitting.getCalcZ(false,true);
		double [] zai0=fieldFitting.getCalcZ(0.0);
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
		double [][] zai=fieldFitting.getCalcZ(true,true);
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
	StringBuffer sb = new StringBuffer();
	for (int n=0;n<radiuses.length;n++){
		sb.append(IJ.d2s(radiuses[n],3));
		if (show_z_axial){
			for (int i=0;i<show_chn.length;i++){
				if (show_chn[i]) sb.append("\t"+IJ.d2s(z_values[0][i][n]-z0,3));
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
    
    double [] center_z=fieldFitting.getZCenters();
    double best_qb_axial= fieldFitting.getBestQualB(
                k_red,
                k_blue,
                false);
    double best_qb_corr= fieldFitting.getBestQualB(
                k_red,
                k_blue,
                true);
     GenericDialog gd = new GenericDialog("Setup quality-B table");
    gd.addMessage("Best center focus for Red "+ IJ.d2s(center_z[0],3)+" um");
    gd.addMessage("Best center focus for Green "+ IJ.d2s(center_z[1],3)+" um");
    gd.addMessage("Best center focus for Blue "+ IJ.d2s(center_z[2],3)+" um");
    gd.addMessage("Best composite distance for FWHM^4, axial model "+ IJ.d2s(best_qb_axial,3)+" um");
    gd.addMessage("Best composite distance for FWHM^4, individual "+ IJ.d2s(best_qb_corr,3)+" um");
        gd.addNumericField("Scan from (relative to green center)", this.qb_scan_below, 3,7,"um");
        gd.addNumericField("Scan to (relative to green center)", this.qb_scan_above, 3,7,"um");
        gd.addNumericField("Scan step", this.qb_scan_step, 3,7,"um");
        gd.addNumericField("Relative (to green) weight of red channel",100* this.k_red, 3,7,"%");
        gd.addNumericField("Relative (to green) weight of blue channel",100* this.k_blue, 3,7,"%");
        gd.addCheckbox("Use per-sample location corrected data", this.qb_use_corrected);
        gd.addCheckbox("Show mtf50 (false - PSF FWHM)", this.qb_invert);
        gd.addCheckbox("Show focal distance relative to center green (false - motor absolute)", this.z_relative);
        
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
             z_relative?center_z[1]:0.0,
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
//        gd.addCheckbox("Dialog after each iteration series", this.stopEachSeries);
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
//        this.stopEachSeries= gd.getNextBoolean();
        this.stopOnFailure= gd.getNextBoolean();
        this.showParams= gd.getNextBoolean();
        this.showDisabledParams= gd.getNextBoolean();
        this.showCorrectionParams= gd.getNextBoolean();

//        this.showThisImages= gd.getNextBoolean();
//        this.showNextImages= gd.getNextBoolean();
     this.saveSeries=true;
     return gd.wasOKed();
}



public boolean LevenbergMarquardt(boolean openDialog, int debugLevel){
    double savedLambda=this.lambda;
    this.debugLevel=debugLevel;
    if (openDialog && !selectLMAParameters()) return false;
    this.startTime=System.nanoTime();
// create savedVector (it depends on parameter masks), restore from it if aborted     
    this.savedVector=this.fieldFitting.createParameterVector(sagittalMaster);
    if (debugDerivativesFxDxDy){
        compareDrDerivatives(this.savedVector);
    }
    this.iterationStepNumber=0;
    this.firstRMS=-1; //undefined
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
            
            if (debugLevel>1) System.out.println(":"+this.iterationStepNumber+": stepLevenbergMarquardtFirst("+debugLevel+")==>"+state[1]+":"+state[0]);
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
//                     (this.stopEachSeries && state[1]) ||
                    (this.stopOnFailure && state[1] && !state[0])){
                
                if (debugLevel>0){
                    if (this.stopRequested.get()>0) System.out.println("User requested stop");
                    System.out.println("LevenbergMarquardt(): step ="+this.iterationStepNumber+
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
                     IJ.showStatus("Step #"+this.iterationStepNumber+
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
        }
//         if (wasLastSeries) break;
//     } // while (this.fittingStrategy.isSeriesValid(this.seriesNumber)){ // TODO: Add "stop" tag to series
    	String msg="RMS="+this.currentRMS+" ("+this.firstRMS+") "+
    	           ", pure RMS="+this.currentRMSPure+" ("+this.firstRMSPure+") "+
    	           " at "+ IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3);
    	if (debugLevel>0) System.out.println("stepLevenbergMarquardtAction() "+msg);
    	//    	if (this.updateStatus) IJ.showStatus(msg);
    	if (this.updateStatus){
    		IJ.showStatus("Done: Step #"+this.iterationStepNumber+
    				" RMS="+IJ.d2s(this.currentRMS,8)+
    				" ("+IJ.d2s(this.firstRMS,8)+")"+
    				" RMSPure="+IJ.d2s(this.currentRMSPure,8)+
    				" ("+IJ.d2s(this.firstRMSPure,8)+")"+
    				" ");
    	}
    	this.savedVector=this.currentVector.clone();
    	commitParameterVector(this.savedVector);
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
        this.currentPX0=pX0_distortions;
        this.currentPY0=pY0_distortions;
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
//    	private double [] sampleCorrRadius=null;
    	private double [][] sampleCoordinates=null;
    	private double [][][][] sampleCorrCrossWeights= new double[6][][][];
    	private double [] sampleCorrVector=null; // currently adjusted per-sample parameters
    	private double [][][] correctionParameters=new double[6][][]; // all
    	public int numberOfLocations=0;

    	private int [][] sampleCorrChnParIndex=null;
    	private boolean [] dflt_sampleCorrSelect= {false,false,false,false};
    	private double [] dflt_sampleCorrCost= {1.0,50.0,1.0,1.0};
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
        public double [] getCalcValuesForZ(double z, double r){
        	double [] result=new double [6];
        	for (int chn=0;chn<result.length;chn++) {
        		if (curvatureModel[chn]!=null){
        			result[chn]=curvatureModel[chn].getFdF(
        					null,
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
if ((chn==4) &&  (sampleIndex==3)){
	System.out.println("getCalcValuesForZ(), chn="+chn+", sampleIndex="+sampleIndex);
}
                    	
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

        public double [] getCalcZ(double r){
            double [] result=new double [6];
            for (int chn=0;chn<result.length;chn++) {
                if (curvatureModel[chn]!=null){
                    result[chn]=curvatureModel[chn].getAr(r, null)[0];
                } else {
                    result[chn]=Double.NaN;
                }
            }
            return result;
        }        
        
        /**
         * calculate distance to "best focus" for each channel (color and S/T) for each sample
         * @param corrected when false - provide averaged (axial model) for radius, if true - with individual correction applied
         * @param allChannels calculate for all (even disabled) channels, false - only for currently selected
         * @return outer dimension - number of channel, inner - number of sample (use getSampleRadiuses for radius of each)
         */
        public double [][] getCalcZ(boolean corrected, boolean allChannels){
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
                        result[chn][sampleIndex]=curvatureModel[chn].getAr(
                                sampleCorrRadius[sampleIndex],
                                corrected?getCorrPar(chn,sampleIndex):null
                                )[0];
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
        public double [] getZCenters(){
//            int chn=3; // Green, Tangential
            double [] result = {
                    curvatureModel[1].getCenterVector()[0], // Red, Tangential
                    curvatureModel[3].getCenterVector()[0], // Green, Tangential
                    curvatureModel[5].getCenterVector()[0]}; // Blue, Tangential
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
            if ((debugLevel>0) && hasBad){
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

        public double getBestQualB(
                double kr,
                double kb,
                boolean corrected){
            return getBestQualB(kr,kb,corrected,1.0,0.001);
            
        }
        
        public double getBestQualB( // find minimum
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
                return Double.NaN;
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
                qb0=getQualB(z1,kr,kb,corrected);
                if (Math.abs(z0-z_prev)<precision) break;
            }
            return z0;
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
        
        /**
         * create matrix of weights of the other parameters influence
         * @param sampleCoordinates [sample number]{x,y} - flattened array of sample coordinates
         * Run in the beginning of fitting series (zeroes the values)
         */
        // once per fitting series (or parameter change
        public void initSampleCorrVector(double [][] sampleCoordinates){
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
                                        sampleCorrCrossWeights[nChn][nPar][i][j]=a;
                                        sw+=a;
                                    }
                                }
                                if ((sampleCorrPullZero[nChn][nPar]==0) ||(sampleCorrCost[nChn][nPar]==0)) sw=0.0;
                                else if (sw!=0.0) sw=-sampleCorrCost[nChn][nPar]*sampleCorrPullZero[nChn][nPar]/sw;
                                for (int j=0;j<numberOfLocations;j++) {
                                    if (i!=j){
                                        sampleCorrCrossWeights[nChn][nPar][i][j]*=sw;
                                    } else {
                                        sampleCorrCrossWeights[nChn][nPar][i][j]=sampleCorrCost[nChn][nPar];
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
            sampleCorrChnParIndex=new int [sampleCorrSelect.length][];
            int numPars=0;
            for (int nChn=0; nChn< sampleCorrCrossWeights.length;nChn++) {
                if (sampleCorrCrossWeights[nChn]!=null){
                    sampleCorrChnParIndex[nChn]=new int [sampleCorrCrossWeights[nChn].length];
                    for (int nPar=0;nPar< sampleCorrCrossWeights[nChn].length;nPar++) {
                        if (sampleCorrCrossWeights[nChn][nPar]!=null){
                            sampleCorrChnParIndex[nChn][nPar]=numPars; // pointer to the first sample
                            numPars+=sampleCorrCrossWeights[nChn][nPar].length;
                        } else {
                            sampleCorrChnParIndex[nChn][nPar]=-1;
                        }
                    }                
                } else {
                    sampleCorrChnParIndex[nChn]=null;
                }
            }
            // currently all correction parameters are initialized as zeros.
            getCorrVector();
            if (debugLevel>1) System.out.println("was resetting sampleCorrVector here");
            
//            sampleCorrVector=new double [numPars];
//            for (int i=0;i<numPars;i++)sampleCorrVector[i]=0.0;
/*            
            sampleCorrRadius=new double [numberOfLocations];
            //pXY
            for (int i=0;i<numberOfLocations;i++){
                double dx=sampleCoordinates[i][0]-pXY[0];
                double dy=sampleCoordinates[i][1]-pXY[0];
                sampleCorrRadius[i]=getPixelMM()*Math.sqrt(dx*dx+dy*dy);
            }
*/            
        }
        
        public double [] getCorrPar(int chn, int sampleIndex){
/*            if ((sampleCorrChnParIndex==null) || (sampleCorrChnParIndex[chn]==null)) return null;
            double [] corr =new double [sampleCorrChnParIndex[chn].length];
            for (int i=0;i<corr.length;i++){
                if (sampleCorrChnParIndex[chn][i]<0) corr[i]=0.0;
                else corr[i]=sampleCorrVector[sampleCorrChnParIndex[chn][i]+sampleIndex];
            }
*/
//    		System.out.println("used sampleCorrVector here, now  correctionParameters");

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
        public boolean maskSetDialog(String title){
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
        	initSampleCorrVector(flattenSampleCoord()); // run always regardless of configured or not (to create zero-length array of corr)
        	return true;
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
                if (centerSelect[i] ) {
                    gd.addNumericField(centerDescriptions[i],pXY[i],5,10,"pix");
                } else if (showDisabled){
                    gd.addNumericField("(disabled) "+centerDescriptions[i],pXY[i],5,10,"pix");
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
                for (int c=0;c<channelSelect.length;c++) if (channelSelect[c]){
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
                    for (int n=0;n<channelSelect.length;n++) if (channelSelect[n]){
                        int [] ncp=curvatureModel[n].getNumPars(); // {(z),(r)}
                        boolean isDependMasterDir=(n%2) == (sagittalMaster?0:1);
                        for (int i=0;i<curvatureSelect[n].length; i++) if (curvatureSelect[n][i] ){
                            if (((i%ncp[1])!=0) || isDependMasterDir) {
                                int dependOnChannel=(((i%ncp[1])==0) && !isMasterDir)?otherChannel:c;
                                deriv[nChn][np++]=(n==dependOnChannel)?(deriv_curv[n][i]):0.0;
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
/*
        public double [] getValsDerivativesOld(
                boolean sagittalMaster,
                int [] motors, // 3 motor coordinates
                double px, // pixel x
                double py, // pixel y
                double [][] deriv // array of (1..6[][], matching getNumberOfChannels) or null if derivatives are not required
                ){
            double [] motorDerivs=(deriv==null)? null:(new double [mechanicalFocusingModel.getNumPars()]);
            double [] chnValues=new double [getNumberOfChannels()];
            double mot_z=mechanicalFocusingModel.calc_ZdZ(
                    motors,
                    px,
                    py,
                    motorDerivs);
            int nChn=0;
            for (int c=0;c<channelSelect.length;c++) if (channelSelect[c]){
                double [] deriv_curv=(deriv==null)?null:(new double [curvatureModel[c].getSize()]);
                chnValues[nChn]=curvatureModel[c].getFdF(
                        null, // param_corr
                        px,
                        py,
                        mot_z,
                        deriv_curv);
                if (deriv!=null){
                    deriv[nChn]=new double [getNumberOfParameters(sagittalMaster)];
                    int np=0;

                    for (int i=0;i<mechanicalFocusingModel.paramValues.length;i++){
                        if ((mechanicalSelect==null) || mechanicalSelect[i] ) deriv[nChn][np++]=-motorDerivs[i]*deriv_curv[0]; // minus d/dz0 const part
                    }
                    
                    for (int n=0;n<channelSelect.length;n++) if (channelSelect[n]){
                        for (int i=0;i<curvatureSelect[n].length; i++) if (curvatureSelect[n][i] ){
                            deriv[nChn][np++]=(n==c)?(deriv_curv[i]):0.0;
                        }
                    }
                }
                nChn++;
            }
            return chnValues;
        }
*/        
        
        
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
                {"z0", "center shift, positive away form the lens","um","0.0"},
                {"tx", "horizontal tilt", "um/mm","0.0"},
                {"ty", "vertical tilt", "um/mm","0.0"}};
        public double PERIOD=3584.0; // steps/revolution
//        public double PIXEL_SIZE=0.0022; // mm
        public double [] paramValues=new double [descriptions.length];
        
        public MechanicalFocusingModel(){ // add arguments?
            initDefaults();
        }
        
    public void setProperties(String prefix,Properties properties){
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
    }
        public void initDefaults(){
            paramValues=new double [descriptions.length];
            for (int i=0;i<descriptions.length;i++) paramValues[i]=Double.parseDouble(descriptions[i][3]);
        }
        public void setVector(double[] vector){
            paramValues=new double [descriptions.length];
            for (int i=0;i<vector.length;i++) paramValues[i]=vector[i];
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
        private double dflt_na=0.15; // um/um
        private double dflt_r0=4.0; //3.3; // um (1.5 pix)
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
f=sqrt((a*(zin-z0))^2 + r0^2)+a0+ a1*(zin-z0)+...aN*(zin-z0)^N
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

*/
            double r=Double.isNaN(pY)?pX:Math.sqrt((pX-pX0)*(pX-pX0)+(pY-pY0)*(pY-pY0))*PIXEL_SIZE; // in mm
//            double r2=r*r;
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
            double z=z_in-ar[0];
            double exp=Math.exp(-ar[3]);
            double reff=ar[2]*exp;
//            double sqrt=Math.sqrt((ar[1]*z)*(ar[1]*z) + ar[2]*ar[2]);
            double sqrt=Math.sqrt((ar[1]*z)*(ar[1]*z) + reff*reff);
            
            double f=sqrt+ar[2]*(1-exp);
            double zp=1.0;
/*            
            for (int i=3;i<ar.length;i++){
                f+=ar[i]*zp;
                zp*=z;
            }
*/            
            for (int i=4;i<ar.length;i++){
                zp*=z;
                f+=ar[i]*zp;
            }

            if (deriv==null) return f; // only value, no derivatives
            double [] df_da=new double[this.modelParams.length]; // last element - derivative for dz
            // derivative for z0 (shift) - ar[0]
            df_da[0]=-1.0/sqrt*ar[1]*ar[1]*z;
            zp=1.0;
            for (int i=4;i<this.modelParams.length;i++){
                df_da[0]-=ar[i]*(i-3)*zp; // ar[i] calculated coefficients for current radius
                zp*=z;
            }
            // derivative for a (related to numeric aperture) - ar[1]
            df_da[1]=1.0/sqrt*ar[1]*z*z;
            // derivative for a (related to lowest PSF radius) - ar[2]
//            df_da[2]=1.0/sqrt*ar[2];
            df_da[2]=1.0/sqrt*reff*exp + (1-exp); // * exp(-k)
            // derivative for k (ar[3]
            
            df_da[3]=1.0/sqrt*reff*ar[2]*exp*(-1) + ar[2]*exp;
            
            // derivatives for rest (polynomial) coefficients
            zp=1.0;
/*            
            for (int i=3;i<this.modelParams.length;i++){
                df_da[i]=zp;
                zp*=z;
            }
*/            
            for (int i=4;i<this.modelParams.length;i++){
                zp*=z;
                df_da[i]=zp;
            }
            // derivative for z (to be combined with mechanical is just negative of derivative for z0, no need to calcualate separately
            // calculate even powers of radius
            double [] dar= new double [this.modelParams[0].length];
            dar[0]=1;
            for (int j=1;j<dar.length;j++){
//                dar[j]=dar[j-1]*r2;
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
            return "ar_"+(i+1);
        }
        public String getRadialDecription(int i){
            if (i==0) return "Radial constant coefficient";
            return "Radial polynomial coefficient for r^"+(i+1);
        }
        public String getZName(int i){
            if (i==0) return "z0";
            if (i==1) return "na";
            if (i==2) return "r0";
            else return "az_"+(i-3);
        }
        public String getZDescription(int i){
            if (i==0) return "Focal shift";
            if (i==1) return "Defocus/focus shift (~NA)";
            if (i==2) return "Best PSF radius";
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
}



