/**
 **
 ** FocusingField.jave - save/restore/process sagittal/tangential PSF width
 ** over FOV, together with related data
 **
 ** Copyright (C) 2014 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  CalibrationHardwareInterface.java is free software: you can redistribute it and/or modify
 **  it under the terms of the GNU General Public License as published by
 **  the Free Software Foundation, either version 3 of the License, or
 **  (at your option) any later version.
 **
 **  This program is distributed in the hope that it will be useful,
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **  GNU General Public License for more details.
 **
 **  You should have received a copy of the GNU General Public License
 **  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ** -----------------------------------------------------------------------------**
 **
 */

import ij.IJ;
import ij.gui.GenericDialog;
import ij.text.TextWindow;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;






//import Distortions.LMAArrays; // may still reuse?
import Jama.LUDecomposition;
import Jama.Matrix;


public class FocusingField {
//	public String path;
	public static final double PIXEL_SIZE=0.0022; // mm
	public static final String sep = " ";
	public static final String regSep = "\\s";
	public String serialNumber;
	public String lensSerial; // if null - do not add average
	public String comment;
	public double pX0;
	public double pY0;
	public double [][][] sampleCoord;
	public ArrayList<FocusingFieldMeasurement> measurements;
	boolean sagittalMaster=false; // center data is the same, when true sagittal fitting only may change r=0 coefficients,
	                              // when false - tangential is master
	double [] minMeas=      {1.0,1.0,1.0,1.0,1.0,1.0}; // pixels
	double [] maxMeas=      {4.5,4.5,4.5,4.5,4.5,4.5}; // pixels
	double [] thresholdMax= {3.0,3.0,3.0,3.0,4.0,3.0}; // pixels
	double [] weightReference=null;
	boolean useMinMeas=     true;
	boolean useMaxMeas=     true;
	boolean useThresholdMax=true;
	FieldFitting fieldFitting=null;
	int weighMode=0; // 0 - same weight, 1 - linear threshold difference, 2 - quadratic thershold difference
	double weightRadius=2.0; // Gaussian sigma in mm
	MeasuredSample [] dataVector;
	double [] dataValues;
	double [] dataWeights;
//	double sumWeights=0.0;
	double [][] jacobian=null; // rows - parameters, columns - samples
	double [] currentVector=null;
	double [] nextVector=null;
	double [] savedVector=null;
    private LMAArrays lMAArrays=null;
    private LMAArrays  savedLMAArrays=null;
    private double []   currentfX=null; // array of "f(x)" - simulated data for all images, combining pixel-X and pixel-Y (odd/even)
    private double []   nextfX=null; // array of "f(x)" - simulated data for all images, combining pixel-X and pixel-Y (odd/even)
    private double      currentRMS=-1.0; // calculated RMS for the currentVector->currentfX
    private double      nextRMS=-1.0; // calculated RMS for the nextVector->nextfX
    private double      firstRMS=-1.0; // RMS before current series of LMA started
    private double lambdaStepUp=   8.0; // multiply lambda by this if result is worse
    private double lambdaStepDown= 0.5; // multiply lambda by this if result is better
    private double thresholdFinish=0.001; // (copied from series) stop iterations if 2 last steps had less improvement (but not worsening ) 
    private int    numIterations=  100; // maximal number of iterations
    private double maxLambda=      100.0;  // max lambda to fail
    
    private double lambda=0.001;        // copied from series
    private double [] lastImprovements= {-1.0,-1.0}; // {last improvement, previous improvement}. If both >0 and < thresholdFinish - done
    private int    iterationStepNumber=0;
    private boolean stopEachStep=  true;  // open dialog after each fitting step
    private boolean stopOnFailure= true;  // open dialog when fitting series failed
    private boolean showParams=   false;   // show modified parameters
    private boolean showDisabledParams = false;
    private long startTime=0;
    private AtomicInteger stopRequested=null; // 1 - stop now, 2 - when convenient
    private boolean saveSeries=false;   // just for the dialog
    
    private boolean showMotors = true;
    private boolean [] showMeasCalc = {true,true,true};
    private boolean [] showColors =   {true,true,true};
    private boolean [] showDirs =     {true,true};
    private boolean [] showSamples =   null;
    private boolean    showAllSamples = true;
    private boolean    showIgnoredData= false;
    private boolean showRad = true;

	
    public class LMAArrays { // reuse from Distortions?
        public double [][] jTByJ=  null; // jacobian multiplied by Jacobian transposed
        public double []   jTByDiff=null; // jacobian multiplied difference vector
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
    
    public class MeasuredSample{
    	public int [] motors = new int[3];
    	public String timestamp;
    	public double px;
    	public double py;
    	public int channel;
    	public double value;
    	//    	public double weight;
//    	public MeasuredSample(){}
    	public MeasuredSample(
    			int [] motors,
    			String timestamp,
    			double px,
    			double py,
    			int channel,
    			double value
    			){
    		this.motors = motors;
    		this.timestamp=timestamp;
    		this.px = px;
    		this.py = py;
    		this.channel = channel;
    		this.value = value;
    	}
    }
    public boolean configureDataVector(String title, boolean forcenew){
    	if ((fieldFitting == null) && !forcenew){
    		forcenew=true;
    	}
    	boolean setupMasks=true,setupParameters=true,showDisabled=true;
    	FieldFitting tmpFieldFitting=fieldFitting;
    	if (tmpFieldFitting==null) tmpFieldFitting=	new FieldFitting(); // just to get field description
    	int [] numCurvPars=tmpFieldFitting.getNumCurvars();
    	GenericDialog gd = new GenericDialog(title+(forcenew?" RESETTING DATA":""));
    	gd.addMessage("=== Setting minimal measured PSF radius for different colors/directions ===");
    	gd.addCheckbox("Sagittal channels are master channels (false - tangential are masters)",sagittalMaster);

    	
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
    		gd.addNumericField("Number of parameters for pzf(z) approximation (>=3)",numCurvPars[0],0);
    		gd.addNumericField("Number of parameters for radial dependence of PSF curves (>=1)",numCurvPars[1],0);
    	}
		gd.addMessage("");
		gd.addNumericField("Data weight mode (0 - equal mode, 1 -linear treshold diff, 2 - quadratic threshold diff)",weighMode,0);
		gd.addNumericField("Data weight radius (multiply weight by Gaussian), 0 - no dependence on radius",weightRadius,3,5,"mm");
    	gd.addCheckbox("Setup parameter masks?",setupMasks);
    	gd.addCheckbox("Setup parameter values?",setupParameters);
    	gd.addCheckbox("Show/modify disabled for auto-adjustment parameters?",showDisabled);
	    WindowTools.addScrollBars(gd);
    	gd.showDialog();
    	if (gd.wasCanceled()) return false;
    	sagittalMaster=                                                gd.getNextBoolean();
    	for (int i=0;i<minMeas.length;i++)this.minMeas[i]=             gd.getNextNumber();
    	useMinMeas=                                                    gd.getNextBoolean();
    	for (int i=0;i<maxMeas.length;i++) this.maxMeas[i]=            gd.getNextNumber();
    	useMaxMeas=                                                    gd.getNextBoolean();
    	for (int i=0;i<thresholdMax.length;i++) this.thresholdMax[i] = gd.getNextNumber();
    	useThresholdMax=                                               gd.getNextBoolean();
    	if (forcenew) {
    		numCurvPars[0] =                                         (int) gd.getNextNumber();
    		numCurvPars[1] =                                         (int) gd.getNextNumber();
    	}
    	weighMode =                                              (int) gd.getNextNumber();
		weightRadius =                                                 gd.getNextNumber();

    	setupMasks=                                                    gd.getNextBoolean();
    	setupParameters=                                               gd.getNextBoolean();
    	showDisabled=                                                  gd.getNextBoolean();
    	if (forcenew) {
    		this.fieldFitting= new FieldFitting(
    				pX0,
    				pY0,
    				numCurvPars[0],
    				numCurvPars[1]);
    	}
    	if (setupMasks) this.fieldFitting.maskSetDialog("Setup parameter masks");
    	if (setupParameters) this.fieldFitting.showModifyParameterValues("Setup parameter values",showDisabled);
    	this.savedVector=this.fieldFitting.createParameterVector(sagittalMaster);
//    	initialVector    	
    	return true;
    }
    

    
    
    // includes deselected channels
    public void setDataVector(MeasuredSample [] vector){ // remove unused channels if any
    	boolean [] chanSel=fieldFitting.getSelectedChannels();
    	int numSamples=0;
    	for (int i=0;i<vector.length;i++) if (chanSel[vector[i].channel]) numSamples++;
    	dataVector=new MeasuredSample [numSamples];
    	int n=0;
    	for (int i=0;i<vector.length;i++) if (chanSel[vector[i].channel]) dataVector[n++]=vector[i];
    	dataValues = new double [dataVector.length];
    	dataWeights = new double [dataVector.length];
 //   	sumWeights=0.0;
    	int mode=weighMode;
    	double kw= (weightRadius>0.0)?(-0.5*getPixelMM()*getPixelMM()/(weightRadius*weightRadius)):0;
//weightRadius    	
    	if (weightReference==null)mode=0;
    	for (int i=0;i<dataVector.length;i++){
    		MeasuredSample ms=dataVector[i];
    		dataValues[i]=ms.value;
    		double diff=weightReference[ms.channel]-ms.value;
    		if (diff<0.0) diff=0;
    		switch (mode){
    		case 0:	dataWeights[i]=1.0; break;
    		case 1:	dataWeights[i]=diff; break;
    		case 2:	dataWeights[i]=diff*diff; break;
    		default: dataWeights[i]=1.0;
    		}
    		if (weightRadius>0.0){
    			double r2=(ms.px-pX0)*(ms.px-pX0)+(ms.py-pY0)*(ms.py-pY0);
    			dataWeights[i]*=Math.exp(kw*r2);
    		}
//    		sumWeights+=dataWeights[i];
    	}
    }
    
// for compatibility with Distortions class    
    public double [] createFXandJacobian( double [] vector, boolean createJacobian){
    	fieldFitting.commitParameterVector(vector,sagittalMaster);
    	return createFXandJacobian(createJacobian);
    }

    
    public double [] createFXandJacobian(boolean createJacobian){
    	double [] fx=new double[dataVector.length];
    	double [][] derivs=null;
    	double [] subData=null;
    	boolean [] selChannels=fieldFitting.getSelectedChannels();
    	int [] selChanIndices= new int[selChannels.length];
    	selChanIndices[0]=0;
    	for (int i=1;i<selChanIndices.length;i++){
    		selChanIndices[i]= selChanIndices[i-1]+(selChannels[i-1]?1:0);
    	}
    	int numPars=fieldFitting.getNumberOfParameters(sagittalMaster);
    	if (createJacobian) {
    		jacobian=new double [numPars][dataVector.length];
    		derivs=new double [fieldFitting.getNumberOfChannels()][];
    	}
    	String prevTimeStamp="";
    	double prevPx=-1,prevPy=-1;
    	for (int n=0;n<dataVector.length;n++){
    		MeasuredSample ms=dataVector[n];
    		if (!ms.timestamp.equals(prevTimeStamp) || (ms.px!=prevPx) || (ms.py!=prevPy)){
    			subData=fieldFitting.getValsDerivatives(
    					sagittalMaster,
    					ms.motors, // 3 motor coordinates
    					ms.px, // pixel x 
    					ms.py, // pixel y
    					derivs);
    			prevTimeStamp=ms.timestamp;
    			prevPx=ms.px;
    			prevPy=ms.py;
    		}
    		fx[n]=subData[selChanIndices[ms.channel]];
        	if (createJacobian) {
        		for (int i=0;i<numPars;i++){
        			jacobian[i][n]=derivs[selChanIndices[ms.channel]][i];
        		}
        	}    		
    		
    	}
    	return fx;
    }
    
    public double getRMS(double [] fx){
    	double sum=0.0;
    	double sum_w=0.0;
    	for (int i=0;i<fx.length;i++){
    		double d=fx[i]-dataValues[i];
    		sum+=dataWeights[i]*d*d;
    		sum_w+=dataWeights[i];
    	}
    	if (sum_w>0) {
    		sum/=sum_w;
    	}
    	return Math.sqrt(sum);
    }
    
    public MeasuredSample [] createDataVector(){ // use this data
        return createDataVector(
        		(this.useMinMeas?this.minMeas:null), // pixels
        		(this.useMaxMeas?this.maxMeas:null), // pixels
        		(this.useThresholdMax?this.thresholdMax:null)); // pixels
    }
    
    /**
     * Generate of usable measurement samples
     * @param minMeas minimal measurement PSF radius (in pixels) - correction that increases result
     * resolution in "sharp" areas (to compensate for measurement error). Individual for each of
     *  the 6 color/direction components.
     * @param maxMeas maximal measurement PSF radius (in pixels) - correction that decreases result
     * resolution out-of-focus areas to compensate for the limited size of the PSF window.
     * Individual for each of the 6 color/direction components.
     * @param thresholdMax maximal PSF radius to consider data usable
     * @return array of the MeasuredSample instances, including motors, PSF radius, channel and value
     */
    public MeasuredSample [] createDataVector(
    		double [] minMeas, // pixels
    		double [] maxMeas, // pixels
    		double [] thresholdMax){ // pixels
    	final int numColors=3;
    	final int numDirs=2;
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
	    	}
		} else {
			thresholdMax=null;
		}
    	for (FocusingFieldMeasurement ffm:measurements){
    		double [][][][] samples=ffm.samples;
    		if (samples!=null) for (int i=0;i<sampleCoord.length;i++){
    			if ((i<samples.length) && (samples[i]!=null)) for (int j=0;j<sampleCoord[i].length;j++){
    				if ((j<samples[i].length) && (samples[i][j]!=null)) for (int c=0;c<numColors;c++){
        				if ((c<samples[i][j].length) && (samples[i][j][c]!=null)) for (int d=0;d<numDirs;d++){
        					int chn=d+numDirs*c;
//        					System.out.println("i="+i+", j="+j+", c="+c+", d="+d);
        					// saved values are PSF radius, convert to FWHM by multiplying by 2.0
        					double value=samples[i][j][c][d]*2.0;
        					// discard above threshold (in pixels, raw FWHM data)
        					if (Double.isNaN(value)) {
        						System.out.println("samples["+i+"]["+j+"]["+c+"]["+d+"] = Double.NaN, motors[0]="+ffm.motors[0]);
        						continue; // bad measurement
        					}
        					if (thresholdMax != null){
        						if (value >= thresholdMax[chn]) continue; // bad measurement (above threshold)
        					}
        					// correct for for minimal measurement; 
        					if (minMeas != null){
        						if (value<minMeas[chn]) continue; // bad measurement (smaller than correction)
        						value=Math.sqrt(value*value-minMeas[chn]*minMeas[chn]);
        					}
        					// correct for for maximal measurement; 
        					if (maxMeas != null) {
        						if (value >= maxMeas[chn]) continue; // bad measurement (larger than correction)
        						value = 1.0/Math.sqrt(1.0/(value*value)-1.0/(maxMeas[chn]*maxMeas[chn]));
        					}
        					// convert to microns from pixels
        					value*=getPixelUM();
        					sampleList.add(new MeasuredSample(
        							ffm.motors,
        							ffm.timestamp,
        							sampleCoord[i][j][0], // px,
        							sampleCoord[i][j][1], // py,
        							chn,
        							value
        							));
        				}    					
    				}
    			}
    		}
    	}
    	System.out.println("createDataVector -> "+sampleList.size()+" elements");
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
	public double calcErrorDiffY(double [] fX){
		double result=0;
		double sumWeights=0;
		if (this.dataWeights!=null) {
			for (int i=0;i<fX.length;i++){
				double diff=this.dataValues[i]-fX[i];
				result+=diff*diff*this.dataWeights[i]; System.out.println(""+i+"fx="+fX[i]+" data="+this.dataValues[i]+" diff="+diff+" w="+this.dataWeights[i]);
				sumWeights+=this.dataWeights[i];
			}
			if (sumWeights>0) result/=sumWeights;
		} else {
			for (int i=0;i<fX.length;i++){
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
//	    	JtByJmod[i][i]+=lambda*JtByJmod[i][i]; //Marquardt mod
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
		double [][] JtByJmod= lMAArrays.jTByJ.clone();
		int numPars=JtByJmod.length;
		for (int i=0;i<numPars;i++){
			JtByJmod[i]=lMAArrays.jTByJ[i].clone();
   			JtByJmod[i][i]+=lambda*JtByJmod[i][i]; //Marquardt mod
		}
//	    M*Ma=Mb
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
//			M.print(10, 5);
	    	return null;
	    }
	    Matrix Ma=M.solve(Mb); // singular
	    return Ma.getColumnPackedCopy();
	}

	/**
	 * Calculates  next parameters vector, holds some arrays
	 * @param numSeries
	 * @return array of two booleans: { improved, finished}
	 */
	public boolean [] stepLevenbergMarquardtFirst(int debugLevel){
		double [] deltas=null;
		if (this.currentVector==null) {
			this.currentVector=this.savedVector.clone();
			this.currentRMS=-1;
			this.currentfX=null; // invalidate
			this.jacobian=null;  // invalidate
			this.lMAArrays=null;
			lastImprovements[0]=-1.0;
			lastImprovements[1]=-1.0;
		}
		// calculate  this.currentfX, this.jacobian if needed
		if (debugLevel>2) {
			System.out.println("this.currentVector");
			for (int i=0;i<this.currentVector.length;i++){
				System.out.println(i+": "+ this.currentVector[i]);
			}
		}
		//    	if ((this.currentfX==null)|| ((this.jacobian==null) && !this.threadedLMA )) {
		if ((this.currentfX==null)|| (this.lMAArrays==null)) {
			if (debugLevel>1) {
				System.out.println(": initial Jacobian matrix calculation. Points:"+this.dataValues.length+" Parameters:"+this.currentVector.length);
			}
			this.currentfX=createFXandJacobian(this.currentVector, true); // is it always true here (this.jacobian==null)
			this.lMAArrays=calculateJacobianArrays(this.currentfX);
			this.currentRMS=calcErrorDiffY(this.currentfX);
			if (debugLevel>1) {
				System.out.println(": initial RMS="+IJ.d2s(this.currentRMS,8)+
						". Calculating next Jacobian. Points:"+this.dataValues.length+" Parameters:"+this.currentVector.length);
			}
		} else {
			this.currentRMS= calcErrorDiffY(this.currentfX);
    	}
    	if (this.firstRMS<0) {
    		this.firstRMS=this.currentRMS;
    	}
	
    	
		deltas=solveLMA(this.lMAArrays,	this.lambda, debugLevel);
		
    	boolean matrixNonSingular=true;
    	if (deltas==null) {
    		deltas=new double[this.currentVector.length];
    		for (int i=0;i<deltas.length;i++) deltas[i]=0.0;
    		matrixNonSingular=false;
    	}
		if (debugLevel>2) {
			System.out.println("deltas");
			for (int i=0;i<deltas.length;i++){
				System.out.println(i+": "+ deltas[i]);
			}
		}
// apply deltas    	
    	this.nextVector=this.currentVector.clone();
    	for (int i=0;i<this.nextVector.length;i++) this.nextVector[i]+=deltas[i];
// another option - do not calculate J now, just fX. and late - calculate both if it was improvement    	
//    	save current Jacobian
		if (debugLevel>2) {
			System.out.println("this.nextVector");
			for (int i=0;i<this.nextVector.length;i++){
				System.out.println(i+": "+ this.nextVector[i]);
			}
		}

//        this.savedJacobian=this.jacobian;
        this.savedLMAArrays=lMAArrays.clone();
        this.jacobian=null; // not needed, just to catch bugs
	    	this.nextfX=createFXandJacobian(this.nextVector,true);
			this.lMAArrays=calculateJacobianArrays(this.nextfX);
  		this.nextRMS= calcErrorDiffY(this.nextfX);
		this.lastImprovements[1]=this.lastImprovements[0];
		this.lastImprovements[0]=this.currentRMS-this.nextRMS;
		if (debugLevel>2) {
			System.out.println("stepLMA this.currentRMS="+this.currentRMS+
					", this.nextRMS="+this.nextRMS+
					", delta="+(this.currentRMS-this.nextRMS));
		}
		boolean [] status={matrixNonSingular && (this.nextRMS<=this.currentRMS),!matrixNonSingular};
		// additional test if "worse" but the difference is too small, it was be caused by computation error, like here:
		//stepLevenbergMarquardtAction() step=27, this.currentRMS=0.17068403807026408,   this.nextRMS=0.1706840380702647
		
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
							", this.nextRMS="+this.nextRMS+
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
//    		this.jacobian=this.savedJacobian;// restore saved Jacobian
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
		if (debugLevel>1) {
			System.out.println(
					"stepLevenbergMarquardtAction() step="+this.iterationStepNumber+
					", this.currentRMS="+this.currentRMS+
					", this.nextRMS="+this.nextRMS+
					" lambda="+this.lambda); //+" at "+IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec");
		}
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
//    	int numSeries=fittingStrategy.getNumSeries();
	    GenericDialog gd = new GenericDialog("Levenberg-Marquardt algorithm parameters for cameras distortions/locations");
//		gd.addNumericField("Iteration number to start (0.."+(numSeries-1)+")", this.seriesNumber, 0);
		gd.addNumericField("Initial LMA Lambda ",            this.lambda, 5);
		gd.addNumericField("Multiply lambda on success",     this.lambdaStepDown, 5);
		gd.addNumericField("Threshold RMS to exit LMA",      this.thresholdFinish, 7,9,"pix");
		gd.addNumericField("Multiply lambda on failure",     this.lambdaStepUp, 5);
		gd.addNumericField("Threshold lambda to fail",       this.maxLambda, 5);
		gd.addNumericField("Maximal number of iterations",   this.numIterations, 0);

		gd.addCheckbox("Dialog after each iteration step",   this.stopEachStep);
//		gd.addCheckbox("Dialog after each iteration series", this.stopEachSeries);
		gd.addCheckbox("Dialog after each failure",          this.stopOnFailure);
//		gd.addCheckbox("Ask for weight function filter",     this.askFilter);
		
		gd.addCheckbox("Show modified parameters",           this.showParams);
		gd.addCheckbox("Show disabled parameters",           this.showDisabledParams);
		
//		gd.addCheckbox("Show debug images before correction",this.showThisImages);
//		gd.addCheckbox("Show debug images after correction", this.showNextImages);
//		gd.addNumericField("Maximal number of threads",   this.threadsMax, 0);
//		gd.addCheckbox("Use memory-saving/multithreaded version", this.threadedLMA);
	    gd.showDialog();
	    if (gd.wasCanceled()) return false;
//	    this.seriesNumber=     (int) gd.getNextNumber();
		this.lambda=                 gd.getNextNumber();
		this.lambdaStepDown=         gd.getNextNumber();
		this.thresholdFinish=        gd.getNextNumber();
		this.lambdaStepUp=           gd.getNextNumber();
		this.maxLambda=              gd.getNextNumber();
		this.numIterations=    (int) gd.getNextNumber();
		this.stopEachStep=           gd.getNextBoolean();
//		this.stopEachSeries=         gd.getNextBoolean();
		this.stopOnFailure=          gd.getNextBoolean();
//		this.askFilter=              gd.getNextBoolean();
		this.showParams=             gd.getNextBoolean();
		this.showDisabledParams=     gd.getNextBoolean();
//		this.showThisImages=         gd.getNextBoolean();
//		this.showNextImages=         gd.getNextBoolean();
//		this.threadsMax=       (int) gd.getNextNumber();
//		this.threadedLMA=            gd.getNextBoolean();
	    return true;
    }

    public void listParameters(String title, String path){
    	String header="State\tDescription\tValue\tUnits";
 		StringBuffer sb = new StringBuffer();
		for (String s:this.fieldFitting.getParameterValueStrings(true)){ //this.showDisabledParams)){ //
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
 		gd.addCheckbox("Show measured PSF FWHM",    this.showMeasCalc[0]);
 		gd.addCheckbox("Show calculated  PSF FWHM", this.showMeasCalc[1]);
 		gd.addCheckbox("Show mechanical distance",  this.showMeasCalc[2]);
 		gd.addCheckbox("Show RED data",             this.showColors[0]);
 		gd.addCheckbox("Show GREEN data",           this.showColors[1]);
 		gd.addCheckbox("Show BLUE data",            this.showColors[2]);
 		gd.addCheckbox("Show SAGGITAL data",        this.showDirs[0]);
 		gd.addCheckbox("Show TANGENTIAL data",      this.showDirs[1]);
    	gd.addMessage("Select field samples");
    	for (int i=0; i<sampleCoord.length;i++) for (int j=0; j<sampleCoord[0].length;j++){
     		gd.addCheckbox("Sample X"+j+"Y"+i+" (x="+sampleCoord[i][j][0]+
     		" y="+sampleCoord[i][j][1]+")",showSamples[j+i*sampleCoord[0].length]);
    	}
    	
 		gd.addCheckbox("Show all samples",          this.showAllSamples);
 		gd.addCheckbox("Show ignored data",         this.showIgnoredData);
 		
 		gd.addCheckbox("Show sample distance fom the optical center",  this.showRad);
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
    		double [][][][] samples=ffm.samples;
    		double [][][][] samplesFull=new double [sampleCoord.length][sampleCoord[0].length][numColors][numDirs];
    		double [][][][] calcSamples=new double [sampleCoord.length][sampleCoord[0].length][numColors][numDirs];
    		double [][]     calcZ=      new double [sampleCoord.length][sampleCoord[0].length];

    		for (int i=0;i<sampleCoord.length;i++)
    			for (int j=0;j<sampleCoord[0].length;j++) {
    				calcZ[i][j]=Double.NaN;
        			for (int c=0;c<numColors;c++)
            			for (int d=0;d<numDirs;d++) {
            				samplesFull[i][j][c][d]=Double.NaN;
            				calcSamples[i][j][c][d]=Double.NaN;
            			}
    			}
    		
    		if (showMeasCalc[0] && (samples!=null)) for (int i=0;i<sampleCoord.length;i++){
    			if ((i<samples.length) && (samples[i]!=null)) for (int j=0;j<sampleCoord[i].length;j++){
    				if ((j<samples[i].length) && (samples[i][j]!=null)) for (int c=0;c<numColors;c++){
        				if ((c<samples[i][j].length) && (samples[i][j][c]!=null)) for (int d=0;d<numDirs;d++){
        					int chn=d+numDirs*c;
//        					System.out.println("i="+i+", j="+j+", c="+c+", d="+d);
        					// saved values are PSF radius, convert to FWHM by multiplying by 2.0
        					double value=samples[i][j][c][d]*2.0;
        					// discard above threshold (in pixels, raw FWHM data)
        					if (Double.isNaN(value)) continue; // bad measurement
        					if (!showIgnoredData && (thresholdMax != null)){
        						if (value >= thresholdMax[chn]) continue; // bad measurement (above threshold)
        					}
        					// correct for for minimal measurement; 
        					if (useMinMeas && (minMeas != null)){
        						if (value<minMeas[chn]) continue; // bad measurement (smaller than correction)
        						value=Math.sqrt(value*value-minMeas[chn]*minMeas[chn]);
        					}
        					// correct for for maximal measurement; 
        					if (useMaxMeas && (maxMeas != null)) {
        						if (value >= maxMeas[chn]) continue; // bad measurement (larger than correction)
        						value = 1.0/Math.sqrt(1.0/(value*value)-1.0/(maxMeas[chn]*maxMeas[chn]));
        					}
        					// convert to microns from pixels
        					value*=getPixelUM();
        					samplesFull[i][j][c][d]=value;
        				}
    				}
    			}
    		}
    		// Now calculate values for the same samples
    		if (showMeasCalc[1]){
    			for (int i=0;i<sampleCoord.length;i++)
    				for (int j=0;j<sampleCoord[0].length;j++) {
    					double [] subData=fieldFitting.getValsDerivatives(
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
    					calcZ[i][j]= 	fieldFitting.getMotorsZ(
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
    
    
    public boolean dialogLMAStep(boolean [] state){
    	String [] states={
    			"Worse, increase lambda",
    			"Better, decrease lambda",
    			"Failed to fit",
    			"Fitting Successful"};
    	int iState=(state[0]?1:0)+(state[1]?2:0);
    
	    GenericDialog gd = new GenericDialog("Levenberg-Marquardt algorithm step");
//    	String [][] parameterDescriptions=fittingStrategy.distortionCalibrationData.parameterDescriptions;
    	gd.addMessage("Current state="+states[iState]);
    	gd.addMessage("Iteration step="+this.iterationStepNumber);
    	
    	gd.addMessage("Initial RMS="+IJ.d2s(this.firstRMS,6)+", Current RMS="+IJ.d2s(this.currentRMS,6)+", new RMS="+IJ.d2s(this.nextRMS,6));
//   	gd.addMessage("Pure initial RMS="+IJ.d2s(this.firstRMSPure,6)+", Current RMS="+IJ.d2s(this.currentRMSPure,6)+", new RMS="+IJ.d2s(this.nextRMSPure,6));
    	
    	if (this.showParams) {
    		gd.addMessage("==== Current parameter values ===");
    		for (String s:this.fieldFitting.getParameterValueStrings(this.showDisabledParams)){ //
        		gd.addMessage(s);
    		}
    		gd.addMessage("");

    	}
    	
		gd.addNumericField("Lambda ",                        this.lambda, 5);
		gd.addNumericField("Multiply lambda on success",     this.lambdaStepDown, 5);
		gd.addNumericField("Threshold RMS to exit LMA",      this.thresholdFinish, 7,9,"pix");
		gd.addNumericField("Multiply lambda on failure",     this.lambdaStepUp, 5);
		gd.addNumericField("Threshold lambda to fail",       this.maxLambda, 5);
		gd.addNumericField("Maximal number of iterations",   this.numIterations, 0);

		gd.addCheckbox("Dialog after each iteration step",   this.stopEachStep);
//		gd.addCheckbox("Dialog after each iteration series", this.stopEachSeries);
		gd.addCheckbox("Dialog after each failure",          this.stopOnFailure);
		gd.addCheckbox("Show modified parameters",           this.showParams);
		gd.addCheckbox("Show disabled parameters",           this.showDisabledParams);
//		gd.addCheckbox("Show debug images before correction",this.showThisImages);
//		gd.addCheckbox("Show debug images after correction", this.showNextImages);
		gd.addMessage("Done will save the current (not new!) state and exit, Continue will proceed according to LMA");
		gd.enableYesNoCancel("Continue", "Done");
		WindowTools.addScrollBars(gd);

	    gd.showDialog();
	    if (gd.wasCanceled()) {
	    	this.saveSeries=false;
	    	return false;
	    }
		this.lambda=                 gd.getNextNumber();
		this.lambdaStepDown=         gd.getNextNumber();
		this.thresholdFinish=        gd.getNextNumber();
		this.lambdaStepUp=           gd.getNextNumber();
		this.maxLambda=              gd.getNextNumber();
		this.numIterations=    (int) gd.getNextNumber();
		this.stopEachStep=           gd.getNextBoolean();
//		this.stopEachSeries=         gd.getNextBoolean();
		this.stopOnFailure=          gd.getNextBoolean();
		this.showParams=             gd.getNextBoolean();
		this.showDisabledParams=     gd.getNextBoolean();

//		this.showThisImages=         gd.getNextBoolean();
//		this.showNextImages=         gd.getNextBoolean();
	    this.saveSeries=true;
	    return gd.wasOKed();
    }
    
    
    
    public boolean LevenbergMarquardt(boolean openDialog, int debugLevel){
    	if (openDialog && !selectLMAParameters()) return false;
    	this.startTime=System.nanoTime();
// create savedVector (it depends on parameter masks), restore from it if aborted     	
    	this.savedVector=this.fieldFitting.createParameterVector(sagittalMaster);
    	this.iterationStepNumber=0;
    	this.firstRMS=-1; //undefined
//    	while (this.fittingStrategy.isSeriesValid(this.seriesNumber)){ // TODO: Add "stop" tag to series
    		this.currentVector=null; // invalidate for the new series
//    		boolean wasLastSeries=false;
    		while (true) { // loop for the same series
    			
    			boolean [] state=stepLevenbergMarquardtFirst(debugLevel);
    			if (state==null) {
    				String msg="Calculation aborted by user request, restoring saved parameter vector";
    				IJ.showMessage(msg);
    				System.out.println(msg);
    				this.fieldFitting.commitParameterVector(this.savedVector,sagittalMaster);
    				return false; 
    			}
    			
    			if (debugLevel>1) System.out.println(":"+this.iterationStepNumber+": stepLevenbergMarquardtFirst("+debugLevel+")==>"+state[1]+":"+state[0]);
				boolean cont=true;
				// Make it success if this.currentRMS<this.firstRMS even if LMA failed to converge
				if (state[1] && !state[0] && (this.firstRMS>this.currentRMS)){
					if (debugLevel>1) System.out.println("LMA failed to converge, but RMS improved from the initial value ("+this.currentRMS+" < "+this.firstRMS+")");
					state[0]=true;
				}
    			if (
    					(this.stopRequested.get()>0) || // graceful stop requested
    					(this.stopEachStep) ||
//    					(this.stopEachSeries && state[1]) ||
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
//    				if (this.showThisImages) showDiff (this.currentfX, "fit-"+this.iterationStepNumber);
//    				if (this.showNextImages) showDiff (this.nextfX,    "fit-"+(this.iterationStepNumber+1));
    			}
				stepLevenbergMarquardtAction(debugLevel); // apply step - in any case?
/*				
				if (this.updateStatus){
   	   				IJ.showStatus(this.seriesNumber+": "+"Step #"+this.iterationStepNumber+
   	   						" RMS="+IJ.d2s(this.currentRMS,8)+
   	   						" ("+IJ.d2s(this.firstRMS,8)+")"+
   	   						" RMSPure="+IJ.d2s(this.currentRMSPure,8)+
   	   						" ("+IJ.d2s(this.firstRMSPure,8)+")"+
   	   						" ");
//   	   				showStatus(this.seriesNumber+": "+"Step #"+this.iterationStepNumber+" RMS="+IJ.d2s(this.currentRMS,8)+ " ("+IJ.d2s(this.firstRMS,8)+")",0);
				}
*/				
				if (!cont){
					if (this.saveSeries) {
						this.savedVector=this.currentVector.clone();
//						saveFittingSeries(); // will save series even if it ended in failure, vector will be only updated
//						updateCameraParametersFromCalculated(true); // update camera parameters from all (even disabled) images
//						updateCameraParametersFromCalculated(false); // update camera parameters from enabled only images (may overwrite some of the above)
					}
					// if RMS was decreased. this.saveSeries==false after dialogLMAStep(state) only if "cancel" was pressed 
					this.fieldFitting.commitParameterVector(this.savedVector,sagittalMaster); // either new or original
					return this.saveSeries; // TODO: Maybe change result?
				}
//stepLevenbergMarquardtAction();    			
    			if (state[1]) {
    				if (!state[0]) {
    					this.fieldFitting.commitParameterVector(this.savedVector,sagittalMaster);
    					return false; // sequence failed
    				}
					this.savedVector=this.currentVector.clone();
//    				saveFittingSeries();
//					updateCameraParametersFromCalculated(true); // update camera parameters from all (even disabled) images
//					updateCameraParametersFromCalculated(false); // update camera parameters from enabled only images (may overwrite some of the above)
//					wasLastSeries=this.fittingStrategy.isLastSeries(this.seriesNumber);
//    				this.seriesNumber++;
    				break; // while (true), proceed to the next series
    			}
    		}
//    		if (wasLastSeries) break;
//    	} // while (this.fittingStrategy.isSeriesValid(this.seriesNumber)){ // TODO: Add "stop" tag to series
		if (debugLevel>0) System.out.println("LevenbergMarquardt(): all done"+
                ", step="+    	this.iterationStepNumber+
				", RMS="+this.currentRMS+
				" ("+this.firstRMS+") "+
				" at "+ IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3));
		this.savedVector=this.currentVector.clone();
		this.fieldFitting.commitParameterVector(this.savedVector,sagittalMaster);
    	return true; // all series done
    	
    }
    
    
    
    
	
	
	
    
	public class FocusingFieldMeasurement{
		public String timestamp;
		public double temperature;
		public int [] motors;
		double [][][][] samples = null;
		
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
				this.samples=new double[samples.length][samples[0].length][samples[0][0].length][2];
				for (int i=0;i<this.samples.length;i++) for (int j=0;j<this.samples[i].length;j++) for (int c=0;c<this.samples[i][j].length;c++){
					try {
						double rt=samples[i][j][c][0]/Math.sqrt((1+samples[i][j][c][1]*samples[i][j][c][1])/2.0); // tangential;
						double rs=rt*samples[i][j][c][1]; // sagittal
						this.samples[i][j][c][0]=rs;
						this.samples[i][j][c][1]=rt;
					} catch (Exception e){
						this.samples[i][j][c][0]=Double.NaN;
						this.samples[i][j][c][1]=Double.NaN;
					}
				}
			}
		}
		public FocusingFieldMeasurement(
				String timestamp,
				double temperature,
				int [] motors,
//				ArrayList<String> sampleStrings
				String [] sampleStrings
				){
			this.timestamp=timestamp;
			this.temperature=temperature;
			this.motors = new int[motors.length];
			for (int i=0;i<motors.length;i++) this.motors[i]= motors[i];
			
//			if ((sampleStrings!=null) && (sampleStrings.size()>0)) {
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
					int colors=(ps.length-2)/2;
					this.samples[i][j]=new double [colors][2];
					for (int c=0;c<colors;c++){
						this.samples[i][j][c][0]=Double.parseDouble(ps[2*c+2]);
						this.samples[i][j][c][1]=Double.parseDouble(ps[2*c+3]);
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
						 sdata += sep+this.samples[i][j][c][0]+ // sagittal
								  sep+this.samples[i][j][c][1]; // tangential
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
		this.serialNumber=serialNumber;
		this.lensSerial=lensSerial;
		this.comment=comment;
		this.pX0=pX0;
		this.pY0=pY0;
		this.sampleCoord=sampleCoord;
		this.measurements=new ArrayList<FocusingFieldMeasurement>();
		this.stopRequested=stopRequested;
	}
	
	public FocusingField(
			boolean smart,       // do not open dialog if default matches 
			String defaultPath, //){
			AtomicInteger stopRequested){
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
			boolean smart,       // do not open dialog if default matches 
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
		comment=      hConfig.getString("comment","no comments");
		serialNumber= hConfig.getString("serialNumber","???");
		lensSerial=   hConfig.getString("lensSerial","???");
		pX0=Double.parseDouble(hConfig.getString("lens_center_x","0.0"));
		pY0=Double.parseDouble(hConfig.getString("lens_center_y","0.0"));
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
    	if (comment!=null)      hConfig.addProperty("comment",comment);
    	if (serialNumber!=null) hConfig.addProperty("serialNumber",serialNumber);
    	if (lensSerial!=null)   hConfig.addProperty("lensSerial",lensSerial);
    	 hConfig.addProperty("lens_center_x",pX0);
    	 hConfig.addProperty("lens_center_y",pY0);
    	 if ((sampleCoord!=null) && (sampleCoord.length>0) && (sampleCoord[0] != null) && (sampleCoord[0].length>0)){
    		 hConfig.addProperty("samples_x",sampleCoord[0].length);
    		 hConfig.addProperty("samples_y",sampleCoord.length);
    		 for (int i=0;i<sampleCoord.length;i++)
    			 for (int j=0;j<sampleCoord[i].length;j++){
//            		 double coord[] = {sampleCoord[i][j][0],sampleCoord[i][j][1]};
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
		MechanicalFocusingModel mechanicalFocusingModel;
		CurvatureModel [] curvatureModel=new CurvatureModel[6]; // 3 colors, sagittal+tangential each
		boolean [] channelSelect=null;
		boolean [] mechanicalSelect=null;
		boolean [][] curvatureSelect=new boolean[6][];
		public final String [] channelDescriptions={
			"Red, sagittal","Red, tangential",
			"Green, sagittal","Green, tangential",
			"Blue, sagittal","Blue, tangential"};
		public String getDescription(int i){
			return channelDescriptions[i];
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
			channelSelect=getDefaultMask();
			mechanicalFocusingModel=new MechanicalFocusingModel();
			mechanicalSelect=mechanicalFocusingModel.getDefaultMask();
			for (int i=0;i<curvatureModel.length;i++){
				curvatureModel[i]= new CurvatureModel(
						pX0,
						pY0,
						distanceParametersNumber,
						radialParametersNumber);
				curvatureSelect[i]=curvatureModel[i].getDefaultMask();
			}
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
	    	
	    	if (channelSelect==null) channelSelect=getDefaultMask();
    		for (int i=0;i<channelSelect.length;i++) {
				gd.addCheckbox(getDescription(i), channelSelect[i]);					
	    	}
			gd.addCheckbox("Edit mechanical parameters masks", editMechMask);
			gd.addCheckbox("Edit curvature model parameters mask(s)", editCurvMask);
			gd.addCheckbox("Apply same curvature model parameters mask to all channels", commonCurvMask);
			gd.addCheckbox("Edit full matrix of the curvature model parameters masks", detailedCurvMask);
    		
//	    	gd.enableYesNoCancel("Keep","Apply"); // default OK (on enter) - "Keep"
        	gd.showDialog();
	    	if (gd.wasCanceled()) return false;
    		for (int i=0;i<channelSelect.length;i++) {
    			channelSelect[i]=gd.getNextBoolean();
    		}
    		editMechMask=gd.getNextBoolean();
    		editCurvMask=gd.getNextBoolean();
    		commonCurvMask=gd.getNextBoolean();
    		detailedCurvMask=gd.getNextBoolean();
//    		boolean OK;
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
    						"Parameter mask for all curvature model channels (colors,S/T)",
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
        						"Parameter mask for curvature model, channel \""+getDescription(i)+"\"",
        						detailedCurvMask,
        						mask);
        				if (mask==null) return false; // canceled
        				curvatureSelect[i]=mask;
    	    		}
    			}
    		}
	    	return true;
		}
		ArrayList<String> getParameterValueStrings(boolean showDisabled){
			ArrayList<String> parList=new ArrayList<String>();
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
			return parList;
		}
		// Show/modify all parameters in a single window.
		public boolean showModifyParameterValues(String title, boolean showDisabled){
	    	GenericDialog gd = new GenericDialog(title);
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
	    	
	    	gd.enableYesNoCancel("Keep","Apply"); // default OK (on enter) - "Keep"
    	    WindowTools.addScrollBars(gd);
        	gd.showDialog();
	    	if (gd.wasCanceled()) return false;
	    	if (!gd.wasOKed()) { // selected non-default "Apply"
				for (int i=0;i<mechanicalFocusingModel.paramValues.length;i++){
					if ((mechanicalSelect==null) || mechanicalSelect[i]  || showDisabled) {
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
		
		/**
		 * @return number of selected parameters (including mechanical and each selected - up to 6 - curvature)
		 */
		public int getNumberOfParameters(boolean sagittalMaster){
			int np=0;
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
			double [] pars = new double [getNumberOfParameters(sagittalMaster)];
			int np=0;
			for (int i=0;i<mechanicalFocusingModel.paramValues.length;i++){
				if ((mechanicalSelect==null) || mechanicalSelect[i] ) pars[np++]=mechanicalFocusingModel.paramValues[i];
			}
			for (int n=0;n<channelSelect.length;n++) if (channelSelect[n]){
				boolean isMasterDir=(n%2) == (sagittalMaster?0:1);
	    		int index=0;
				for (int i=0;i<curvatureModel[n].modelParams.length;i++) for (int j=0;j<curvatureModel[n].modelParams[0].length;j++){
					if ((isMasterDir || (j!=0)) && ((curvatureSelect[n]==null) || curvatureSelect[n][index] )) {
						pars[np++]=curvatureModel[n].modelParams[i][j];
					}
					index++;
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
// copy center parameters to dependent
// copy if master is selected, regardless of is dependent selected or not			
			for (int n=0;n<channelSelect.length;n++) if (channelSelect[n]){
				boolean isMasterDir = (n%2) == (sagittalMaster?0:1);
				if (isMasterDir){
					curvatureModel[n ^ 1].setCenterVector(curvatureModel[n].getCenterVector());
				}
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
			double pd=Math.sqrt((px-pX0)*(px-pX0)+(py-pY0)*(py-pY0));
			return pd*getPixelMM();
			
		}
		/**
		 * Generate vector of function values (up to 6 - 1 per selected channel), corresponding to motor positions and pixel
		 * coordinates, optionally generate partial derivatives for each channel and parameter
		 * @param motors array of 3 motors positions
		 * @param px pixel X-coordinate of the sample center
		 * @param py pixel Y-coordinate of the sample center
		 * @param deriv 2d array with outer dimension correspond to number of selected channels (getNumberOfChannels())
		 *  or null if derivatives are not needed, just values
		 * @return array of [getNumberOfChannels()] calculated function values
		 */
		public double [] getValsDerivatives(
				boolean sagittalMaster, // false - tangential master, true - sagittal master (for center coefficients)
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
			double [][] deriv_curv = new double [channelSelect.length][];
			for (int c=0;c<channelSelect.length;c++) deriv_curv[c]=null; 
			for (int c=0;c<channelSelect.length;c++) if (channelSelect[c]){
				deriv_curv[c]=(deriv==null)?null:(new double [curvatureModel[c].getSize()]);
				chnValues[nChn++]=curvatureModel[c].getFdF(
						px,
						py,
						mot_z,
						deriv_curv[c]);
			}
			if (deriv!=null){
				nChn=0;
				for (int c=0;c<channelSelect.length;c++) if (channelSelect[c]){
					boolean isMasterDir=(c%2) == (sagittalMaster?0:1);
					int otherChannel= c ^ 1;
					deriv[nChn]=new double [getNumberOfParameters(sagittalMaster)];
					int np=0;
					// For dependent/master channels no difference for mechanical parameters
					// TOSO: verify they are the same?

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
					nChn++;
				}
			}
			return chnValues;
		}

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
				{"z0",  "center shift, positive away form the lens","um","0.0"},
				{"tx",  "horizontal tilt", "um/mm","0.0"},
				{"ty",  "vertical tilt", "um/mm","0.0"}};
		public double PERIOD=3584.0; // steps/revolution
//		public double PIXEL_SIZE=0.0022; // mm
		public double [] paramValues=new double [descriptions.length];
		
		public MechanicalFocusingModel(){ // add arguments?
			initDefaults();
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
			/*
			 kM3=K0+KD3
			 kM1=K0+KD1-KD3
			 kM2=K0-KD1-KD3
			 K0 may be fixed (overall scale), kM1, kM2 - variable

			dZc(m1)=  0.25* kM1 * (m1 + sM1*P/(2*pi)*sin(2pi*m1/P) + cM1*P/(2*pi)*cos(2pi*m1/P))
			dZc(m2)= 0.25* kM2 * (m2 + sM2*P/(2*pi)*sin(2pi*m2/P) + cM2*P/(2*pi)*cos(2pi*m2/P))
			dZc(m3)= 0.5 *  kM3 *( m3 + sM3*P/(2*pi)*sin(2pi*m3/P) + cM3*P/(2*pi)*cos(2pi*m3/P))

			Assuming X - towards M3, Y - towards M1

			d2Z/dX/dm1=-ps/(4*Lx)* kM1 * (m1 + sM1*P/(2*pi)*sin(2pi*m1/P) + cM1*P/(2*pi)*cos(2pi*m1/P))
			d2Z/dX/dm2=-ps/(4*Lx)* kM2 *( m2 + sM2*P/(2*pi)*sin(2pi*m2/P) + cM2*P/(2*pi)*cos(2pi*m2/P))
			d2Z/dX/dm3= ps/(2*Lx)* kM3 * (m3 + sM3*P/(2*pi)*sin(2pi*m3/P) + cM3*P/(2*pi)*cos(2pi*m3/P))

			d2Z/dY/dm1=+ps/(2*Ly)* kM1 * (m1 + sM1*P/(2*pi)*sin(2pi*m1/P) + cM1*P/(2*pi)*cos(2pi*m1/P))
			d2Z/dY/dm2= -ps/(2*Ly)* kM2 * (m2 + sM2*P/(2*pi)*sin(2pi*m2/P) + cM2*P/(2*pi)*cos(2pi*m2/P))
			d2Z/dY/dm3= 0

			*/
//			double [] deriv=new double [paramValues.length];
			double kM1=	getValue(MECH_PAR.K0)+getValue(MECH_PAR.KD1)-getValue(MECH_PAR.KD3);
			double kM2=	getValue(MECH_PAR.K0)-getValue(MECH_PAR.KD1)-getValue(MECH_PAR.KD3);
			double kM3=	getValue(MECH_PAR.K0)+getValue(MECH_PAR.KD3);
			double p2pi= PERIOD/2/Math.PI;
			double m1=motors[0],m2=motors[1],m3=motors[2];
			double aM1=(m1 + getValue(MECH_PAR.sM1)*p2pi*Math.sin(m1/p2pi) + getValue(MECH_PAR.cM1)*p2pi*Math.cos(m1/p2pi));
			double aM2=(m2 + getValue(MECH_PAR.sM2)*p2pi*Math.sin(m2/p2pi) + getValue(MECH_PAR.cM2)*p2pi*Math.cos(m2/p2pi));
			double aM3=(m3 + getValue(MECH_PAR.sM3)*p2pi*Math.sin(m3/p2pi) + getValue(MECH_PAR.cM3)*p2pi*Math.cos(m3/p2pi));
			double zM1=kM1 * aM1;
			double zM2=kM2 * aM2;
			double zM3=kM3 * aM3;

			double zc= 0.25* zM1+ 0.25* zM2+ 0.5 * zM3;
			double dx=PIXEL_SIZE*(px-getValue(MECH_PAR.mpX0));
			double dy=PIXEL_SIZE*(py-getValue(MECH_PAR.mpY0));
			double zx=dx*(getValue(MECH_PAR.tx)+(2*zM3-zM1-zM2)/(4*getValue(MECH_PAR.Lx))) ;
			double zy=dy*(getValue(MECH_PAR.ty)+(zM1-zM2)/(2*getValue(MECH_PAR.Ly)));
			double z=zc+zx+zy;
			if (deriv==null) return z;
			for (int i=0;i<deriv.length;i++) deriv[i]=0.0;
			// Above same as calc_Z
			double dx_mpX0=-PIXEL_SIZE;
			double dy_mpY0=-PIXEL_SIZE;
			double zM1_K0=  aM1;
			double zM1_KD1= aM1;
			double zM1_KD3=-aM1;
			double zM2_K0=  aM2;
			double zM2_KD1=-aM2;
			double zM2_KD3=-aM2;
			double zM3_K0=  aM3;
			double zM3_KD1=-aM3;
			double zM3_KD3= 0.0;
			double zM1_sM1= kM1*p2pi*Math.sin(m1/p2pi);
			double zM1_cM1= kM1*p2pi*Math.cos(m1/p2pi);
			double zM2_sM2= kM2*p2pi*Math.sin(m2/p2pi);
			double zM2_cM2= kM2*p2pi*Math.cos(m2/p2pi);
			double zM3_sM3= kM3*p2pi*Math.sin(m3/p2pi);
			double zM3_cM3= kM3*p2pi*Math.cos(m3/p2pi);
			
			double zc_K0=  0.25* zM1_K0+ 0.25* zM2_K0+ 0.5 * zM3_K0;
			double zc_KD1= 0.25* zM1_KD1+ 0.25* zM2_KD1+ 0.5 * zM3_KD1;
			double zc_KD3= 0.25* zM1_KD3+ 0.25* zM2_KD3+ 0.5 * zM3_KD3;
			double zc_sM1= 0.25* zM1_sM1;
			double zc_cM1= 0.25* zM1_cM1;
			double zc_sM2= 0.25* zM2_sM2;
			double zc_cM2= 0.25* zM2_cM2;
			double zc_sM3= 0.5*  zM3_sM3;
			double zc_cM3= 0.5*  zM3_cM3;
			
//			double zx_K0=(2*zM3-zM1-zM2)* dx/(4*getValue(MECH_PAR.Lx));
			double zx_a=dx/(4*getValue(MECH_PAR.Lx));
			double zx_K0= (2*zM3_K0-zM1_K0-zM2_K0)*zx_a;
			double zx_KD1=(2*zM3_KD1-zM1_KD1-zM2_KD1)*zx_a;
			double zx_KD3=(2*zM3_K0-zM1_K0-zM2_K0)*zx_a;
			double zx_sM1= (-zM1_sM1)*zx_a;
			double zx_cM1= (-zM1_cM1)*zx_a;
			double zx_sM2= (-zM2_sM2)*zx_a;
			double zx_cM2= (-zM2_cM2)*zx_a;
			double zx_sM3= (2*zM3_sM3)*zx_a;
			double zx_cM3= (2*zM3_cM3)*zx_a;
			double zx_mpX0=dx_mpX0/(4*getValue(MECH_PAR.Lx));
			double zx_tx=  dx;
			double zx_Lx= -dx*(2*zM3-zM1-zM2)/(4*getValue(MECH_PAR.Lx)*getValue(MECH_PAR.Lx));

			double zy_a=dx/(2*getValue(MECH_PAR.Ly));
			double zy_K0=  (zM1_K0- zM2_K0) *zy_a;
			double zy_KD1= (zM1_KD1-zM2_KD1)*zy_a;
			double zy_KD3= (zM1_KD3-zM2_KD3)*zy_a;
			double zy_sM1= (zM1_sM1)*zy_a;
			double zy_cM1= (zM1_cM1)*zy_a;
			double zy_sM2= (-zM2_sM2)*zy_a;
			double zy_cM2= (-zM2_cM2)*zy_a;
			double zy_sM3= 0.0;
			double zy_cM3= 0.0;
			double zy_mpY0=dy_mpY0/(2*getValue(MECH_PAR.Ly));
			double zy_ty=  dy;
			double zy_Ly= -dy*(zM1-zM2)/(2*getValue(MECH_PAR.Ly)*getValue(MECH_PAR.Ly));

			deriv[getIndex(MECH_PAR.K0)]= zc_K0+zx_K0+zy_K0;
			deriv[getIndex(MECH_PAR.KD1)]=zc_KD1+zx_KD1+zy_KD1;
			deriv[getIndex(MECH_PAR.KD3)]=zc_KD3+zx_KD3+zy_KD3;
			deriv[getIndex(MECH_PAR.sM1)]=zc_sM1+zx_sM1+zy_sM1;
			deriv[getIndex(MECH_PAR.cM1)]=zc_cM1+zx_cM1+zy_cM1;
			deriv[getIndex(MECH_PAR.sM2)]=zc_sM2+zx_sM2+zy_sM2;
			deriv[getIndex(MECH_PAR.cM2)]=zc_cM2+zx_cM2+zy_cM2;
			deriv[getIndex(MECH_PAR.sM3)]=zc_sM3+zx_sM3+zy_sM3;
			deriv[getIndex(MECH_PAR.cM3)]=zc_cM3+zx_cM3+zy_cM3;
			deriv[getIndex(MECH_PAR.Lx)] =       zx_Lx;
			deriv[getIndex(MECH_PAR.Ly)] =              zy_Ly;
			deriv[getIndex(MECH_PAR.mpX0)] =     zx_mpX0;
			deriv[getIndex(MECH_PAR.mpY0)] =            zy_mpY0;
			deriv[getIndex(MECH_PAR.tx)] =       zx_tx;
			deriv[getIndex(MECH_PAR.ty)] =              zy_ty;
			return  z;
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
	    	gd.enableYesNoCancel("Keep","Apply"); // default OK (on enter) - "Keep"
        	gd.showDialog();
	    	if (gd.wasCanceled()) return null;
	    	if (!gd.wasOKed()) { // selected non-default "Apply"
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
//					gd.addMessage(getDescription(i) +":   "+this.paramValues[i]+" ("+getUnits(i)+")");
					gd.addNumericField("(disabled) "+getDescription(i),this.paramValues[i],5,8,getUnits(i));
				}
			}
	    	gd.enableYesNoCancel("Keep","Apply"); // default OK (on enter) - "Keep"
        	gd.showDialog();
	    	if (gd.wasCanceled()) return false;
	    	if (!gd.wasOKed()) { // selected non-default "Apply"
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
		public static final int dflt_distanceParametersNumber=4;
		public static final int dflt_radialParametersNumber=3;
		private static final int min_distanceParametersNumber=3;
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
				return this.modelParams.length*this.modelParams[0].length;
			} catch (Exception e){
				return 0;
			}
		}
		/**
		 * Calculate function value and (optionally) partial derivatives over all parameters as 1d array
		 * inner scanning by radial coefficient (for even powers of radius), outer - for f(z) parameters
		 * first parameter - z0 (shift), then scale (related to numeric aperture), then r0 (related to minimal PSF radius),
		 * other parameters - just polynomial coefficients for (z_in-z0)^n
		 * @param pX - sample pixel x coordinate (currently just for radius calculation)
		 * @param pY - sample pixel y coordinate  (currently just for radius calculation)
		 * @param z_in - z (from mechanical) before subtracting of  z0
		 * @param deriv array to accommodate all partial derivatives or null
		 * @return function value
		 */
		public double getFdF(
				double pX,
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
z0 -  ar[0]
a -   ar[1]
r0 -  ar[2]
k  -  ar[3]
ar1 - ar[4]

*/
			double r=Math.sqrt((pX-pX0)*(pX-pX0)+(pY-pY0)*(pY-pY0))*PIXEL_SIZE; // in mm
			double r2=r*r;
			double [] ar= new double [this.modelParams.length];
			for (int i=0;i<this.modelParams.length;i++){
				ar[i]=this.modelParams[i][0];
				double rp=1.0;
				for (int j=1;j<this.modelParams[i].length;j++){
					rp*=r2;
					ar[i]+=this.modelParams[i][j]*rp;
				}
			}
			double z=z_in-ar[0];
			double exp=Math.exp(-ar[3]);
			double reff=ar[2]*exp;
//			double sqrt=Math.sqrt((ar[1]*z)*(ar[1]*z) + ar[2]*ar[2]);
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
//			df_da[2]=1.0/sqrt*ar[2];
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
				dar[j]=dar[j-1]*r2;
			}
			int index=0;
			for (int i=0;i<this.modelParams.length;i++) for (int j=0;j<this.modelParams[0].length;j++){
				deriv[index++]=df_da[i]*dar[j];
			}
			return f; // function value
		}
		public String getRadialName(int i){
			return "ar_"+(2*i);
		}
		public String getRadialDecription(int i){
			if (i==0) return "Radial constant coefficient";
			return "Radial polynomial coefficient for r^"+(2*i);
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
	    	gd.enableYesNoCancel("Keep","Apply"); // default OK (on enter) - "Keep"
    	    WindowTools.addScrollBars(gd);
        	gd.showDialog();
	    	if (gd.wasCanceled()) return null;
	    	if (!gd.wasOKed()) { // selected non-default "Apply"
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
//					gd.addMessage(name +":   "+this.modelParams[i][j]);
				}
				index++;
			}
	    	gd.enableYesNoCancel("Keep","Apply"); // default OK (on enter) - "Keep"
    	    WindowTools.addScrollBars(gd);
        	gd.showDialog();
	    	if (gd.wasCanceled()) return false;
	    	if (!gd.wasOKed()) { // selected non-default "Apply"
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
