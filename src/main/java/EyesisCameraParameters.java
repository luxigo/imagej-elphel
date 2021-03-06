/*
 **
 ** EyesisCameraParameters.java
 **
 ** Copyright (C) 2011-2014 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  EyesisCameraParameters.java is free software: you can redistribute it and/or modify
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
import ij.gui.GenericDialog;

import java.util.Properties;

import org.apache.commons.configuration.XMLConfiguration;


    public  class EyesisCameraParameters{
    	final public String [] distortionModelDescriptions= {
    		"Radial model",
    		"Non radial with shift/elongation, non cummulative",
    		"Non radial with shift/elongation, cummulative",
    		"With non-radial polynomial terms"
    	};
    	final int [] distortionModels={0,100,101,200};
    	public int defaultLensDistortionModel=200;
    	public double [] goniometerHorizontal; // goniometer rotation around "horizontal" axis (tilting from the target - positive)
    	public double [] goniometerAxial; // goniometer rotation around Eyesis axis (clockwise in plan - positive 
    	public EyesisSubCameraParameters [][] eyesisSubCameras=null;
    	public double [] interAxisDistance; // distance in mm between two goniometer axes
    	public double [] interAxisAngle;    // angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
    	                            // clockwise when eyesis is in 'normal' position, looking to the target
    	public double [] horAxisErrPhi;   // angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)
    	public double [] horAxisErrPsi;   // angle in degrees "horizontal" goniometer axis is rotated around moving X axis (up)
    	public double [] entrancePupilForward; // common to all lenses - distance from the sensor to the lens entrance pupil
    	public double [] centerAboveHorizontal; // camera center distance along camera axis above the closest point to horizontal rotation axis (adds to height of each SFE)
    	public double [][] GXYZ=null; // [numStations]{x,y,z}  coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system

    	// non-adjustable parameters, not parts of vector
    	public int     numStations;
    	public double [] stationWeight; // reprojection error weights (close station - relax errors)
    	public boolean isTripod=false; // when true - make goniometerHorizontal rotation around "vertical" axis and "goniometerAxial" - around 
        // rotated horizontal.
    	public int sensorWidth=      2592;
    	public int sensorHeight=     1936;
    	public int    shrinkGridForMask=4; //2; //shrink detected grids by one point for/vert this number of times before calculating masks
    	public double maskBlurSigma=    -3; //2.0;   // blur sensor masks (>0 - pixels, <0 - in grid units)
    	public int    decimateMasks=     1;
    	public double badNodeThreshold=0.1; // filter out grid nodes with difference from quadratically predicted from 8 neighbors in pixels
    	public int    maxBadNeighb=      1; // maximal number of bad nodes around the corrected one to fix
    	public int    minimalValidNodes=50; // do not use images with less than this number of non-zero nodes (after all applicable weight masks) 
    	public int    weightMultiImageMode=1; // increase weight for multi-image sets (0 - do not increase, 1 - multiply by number of images in a set to weightMultiExponent power)
    	public double weightMultiExponent= 1.0;
    	public double weightDiameterExponent=1.0; // if( >0) use grid diameter to scale weights of this image
    	public double weightYtoX=1.0; // relative Y-to-X errors weight (to somewhat compensate for rectabular shape of the sensor)
    	public double minimalGridContrast=0.4; // (normally max ~0.8)
    	public double shrinkBlurSigma = 4.0;
    	public double shrinkBlurLevel = 0.5;
    	public double balanceChannelWeightsMode=-1.0; // <0 - use defaults, 0 - keep, >0 balance by number of points to this power 
    	public double removeOverRMS=2.0;            // error is multiplied by weight function before comparison (more permissive on the borders
    	public double removeOverRMSNonweighted=4.0; // error is not multiplied (no more permissions on tyhe borders
        public int [] extrinsicIndices={6,7,8,9,10,11,12,13,14,15,16}; //support variations
        public double [][] variationsDefaults={
        		null,               // 0
        		null,               // 1
        		null,               // 2
        		null,               // 3
        		null,               // 4
        		null,               // 5
        		{10.0,0.1,0.0,1.0}, // 6 goniometerHorizontal
        		{10.0,0.1,0.0,1.0}, // 7 goniometerAxial
        		{10.0,2.0,0.0,1.0}, // 8 interAxisDistance
        		{10.0,0.2,0.0,1.0}, // 9 interAxisAngle
        		{10.0,0.2,0.0,1.0}, // 10 horAxisErrPhi
        		{10.0,0.2,0.0,1.0}, // 11 horAxisErrPsi
        		{10.0,2.0,0.0,1.0}, // 12 entrancePupilForward
        		{10.0,2.0,0.0,1.0}, // 13 centerAboveHorizontal
        		{10.0,2.0,0.0,1.0}, // 14 GXYZ0
        		{10.0,2.0,0.0,1.0}, // 15 GXYZ1
        		{10.0,2.0,0.0,1.0}  // 16 GXYZ2
        };
        public int    tiltIndex=6;
        private ParameterVariationCosts []  parameterVariationCosts=null;
        
        public boolean isExtrinsic (int index){
        	for (int i=0;i<this.extrinsicIndices.length; i++) if (this.extrinsicIndices[i]==index) return true;
        	return false;
        }
        public boolean isTilt (int index){
        	return (index==this.tiltIndex);
        }
    	public boolean isVarianceCostSet(int index){
    		return (parameterVariationCosts!=null) && (index< parameterVariationCosts.length) && (parameterVariationCosts[index]!=null);
    	}
    	public double varianceCostScale(int index){
    		if (!isVarianceCostSet(index)) return -1.0;
    		return parameterVariationCosts[index].scale;
    	}
    	public double varianceCostVariationAbs(int index){
    		if (!isVarianceCostSet(index)) return -1.0;
    		return parameterVariationCosts[index].variationAbs;
    	}
    	public double varianceCostVariationDiff(int index){
    		if (!isVarianceCostSet(index)) return -1.0;
    		return parameterVariationCosts[index].variationDiff;
    	}
    	public double varianceCostVariationExponent(int index){
    		if (!isVarianceCostSet(index)) return -1.0;
    		return parameterVariationCosts[index].exponent;
    	}
    	public String varianceCostVariationName(int index){
    		if (!isVarianceCostSet(index)) return null;
    		return parameterVariationCosts[index].parName;
    	}
    	public int getVarParsLength(){
    		return parameterVariationCosts.length;
    	}
    	public boolean editCostProperties(
    			int index,
       			String parameterName,
       			String parameterDescription,
       			String parameterUnits) {
    		if (!isExtrinsic(index)) return false; // does not have varaince parameters
    		getVariationCosts();
    		if (parameterVariationCosts[index]==null) parameterVariationCosts[index]=new ParameterVariationCosts(
					this.variationsDefaults[index][0],	
					this.variationsDefaults[index][1],	
					this.variationsDefaults[index][2],	
					this.variationsDefaults[index][3]);
    		return parameterVariationCosts[index].showVarianceDialog(
           			parameterName,
           			parameterDescription,
           			parameterUnits);
    	}
        private ParameterVariationCosts [] getVariationCosts(){
        	if (this.parameterVariationCosts==null){
            	int max=0;
            	for (int i=0;i<this.extrinsicIndices.length; i++) if (this.extrinsicIndices[i]>max) max=this.extrinsicIndices[i];
        		this.parameterVariationCosts=new ParameterVariationCosts [max+1];
        		for (int i=0;i<this.parameterVariationCosts.length;i++) this.parameterVariationCosts[i]=null;
        	}
        	return this.parameterVariationCosts;
        }
        public void getCostsPropertiesXML(String prefix,XMLConfiguration hConfig){
    		getVariationCosts();
    		for (int i=0;i<this.extrinsicIndices.length;i++){
    			int index=this.extrinsicIndices[i];
//    			System.out.println("getCostsPropertiesXML("+prefix+",hconfig) index="+index);
            		if (hConfig.configurationsAt(prefix+"varianceCosts_"+index).size()!=0){
//    				System.out.println("hConfig.configurationAt(prefix+\"varianceCosts_\"+index).isEmpty()=false");
    				this.parameterVariationCosts[index]=new ParameterVariationCosts(
    						this.variationsDefaults[index][0],	
    						this.variationsDefaults[index][1],	
    						this.variationsDefaults[index][2],	
    						this.variationsDefaults[index][3]);
    				boolean isSet= this.parameterVariationCosts[index].getPropertiesXML(prefix+"varianceCosts_"+index+".", hConfig);
    				if (!isSet) this.parameterVariationCosts[index]=null;
    			}
    		}
    	}      
        public void setCostsPropertiesXML(String prefix,XMLConfiguration hConfig){
    		if (this.parameterVariationCosts==null) return;
    		for (int i=0;i<this.extrinsicIndices.length;i++){
    			int index=this.extrinsicIndices[i];
    			if (this.parameterVariationCosts[index]!=null) {
    				hConfig.addProperty(prefix+"varianceCosts_"+index,"");
    				this.parameterVariationCosts[index].setPropertiesXML(prefix+"varianceCosts_"+index+".", hConfig);
    			}
    		}
    	}      
        private class ParameterVariationCosts{
    		public double scale = 1.0; // 1 pixel
    		public double variationAbs  =0.0; // variation of the parameter to cost 1 pixel
    		public double variationDiff =0.0; // variation of the parameter to cost 1 pixel
    		public double exponent=         1.0; // 1.0 - square diff
    		public String parName=null;
 //   		public ParameterVariationCosts(){}
    		public ParameterVariationCosts(
    				double scale,
    				double variationAbs,
    				double variationDiff,
    				double exponent,
    				String parName){
				this.scale=scale;
				this.variationAbs=variationAbs;
				this.variationDiff=variationDiff;
				this.exponent=exponent;
				this.parName=parName;
    		}
    		public ParameterVariationCosts(
    				double scale,
    				double variationAbs,
    				double variationDiff,
    				double exponent){
				this.scale=scale;
				this.variationAbs=variationAbs;
				this.variationDiff=variationDiff;
				this.exponent=exponent;
    		}
    		public ParameterVariationCosts clone(){
    			return new ParameterVariationCosts(
    					this.scale,
    					this.variationAbs,
    					this.variationDiff,
    					this.exponent,
    					this.parName);
    		}

        	public void setPropertiesXML(String prefix,XMLConfiguration hConfig){
        		hConfig.addProperty(prefix+"scale",this.scale+"");
        		hConfig.addProperty(prefix+"variationAbs",this.variationAbs+"");
        		hConfig.addProperty(prefix+"variationDiff",this.variationDiff+"");
        		hConfig.addProperty(prefix+"exponent",this.exponent+"");
        		if (this.parName!=null) hConfig.addProperty(prefix+"parName",this.parName+"");
        	}
        	public boolean getPropertiesXML(String prefix,XMLConfiguration hConfig){
//    			System.out.println("getPropertiesXML("+prefix+",hconfig)");
        		boolean isSet=false;
        		
        		if (hConfig.getString(prefix+"scale")!=null){
//        			System.out.println("getPropertiesXML("+prefix+",hconfig), hConfig.getString(prefix+\"scale\")!=null");
        			this.scale=Double.parseDouble(hConfig.getString(prefix+"scale"));
        			isSet=true;
        		}
        		if (hConfig.getString(prefix+"variationAbs")!=null){
//        			System.out.println("getPropertiesXML("+prefix+",hconfig), hConfig.getString(prefix+\"variationAbs\")!=null");
        			this.variationAbs=Double.parseDouble(hConfig.getString(prefix+"variationAbs"));
        			isSet=true;
        		}
        		if (hConfig.getString(prefix+"variationDiff")!=null){
//       			System.out.println("getPropertiesXML("+prefix+",hconfig), hConfig.getString(prefix+\"variationDiff\")!=null");
        			this.variationDiff=Double.parseDouble(hConfig.getString(prefix+"variationDiff"));
        			isSet=true;
        		}
        		if (hConfig.getString(prefix+"exponent")!=null){
//        			System.out.println("getPropertiesXML("+prefix+",hconfig), hConfig.getString(prefix+\"exponent\")!=null");
        			this.exponent=Double.parseDouble(hConfig.getString(prefix+"exponent"));
        			isSet=true;
        		}
        		if (hConfig.getString(prefix+"parName")!=null){
//        			System.out.println("getPropertiesXML("+prefix+",hconfig), hConfig.getString(prefix+\"parName\")!=null");
        			this.parName=hConfig.getString(prefix+"parName");
        			isSet=true;
        		}
        		return isSet;
        	}

           	public boolean showVarianceDialog(
           			String parameterName,
           			String parameterDescription,
           			String parameterUnits) {
           		String title="Setup costs for image set variance of "+parameterName;
        		GenericDialog gd = new GenericDialog(title);
        		gd.addMessage("Parameter: "+parameterName+" - "+parameterDescription);
				gd.addNumericField("Effective cost of the image set when the "+parameterName+" variance reaches value below",this.scale, 2,4,"pix");
				gd.addNumericField("Variance amount from the average of "+parameterName+" to result in the specified cost", this.variationAbs, 3,8,parameterUnits);
				gd.addNumericField("Variance amount from the p to result in the specified cost", this.variationDiff, 3,8,parameterUnits);
				gd.addNumericField("Exponent of the cost vs. variance (1.0 - squared err "+parameterName+" to result in the specified cost", this.exponent, 3,8,"");
        	    WindowTools.addScrollBars(gd);
        		gd.showDialog();
        		if (gd.wasCanceled()) return false;
				this.scale=         gd.getNextNumber();
				this.variationAbs=  gd.getNextNumber();
				this.variationDiff= gd.getNextNumber();
				this.exponent= gd.getNextNumber();
				this.parName=parameterName;
				return true;
           	}
    	}

    	public int getNumStations() {return this.numStations;}
    	public int getNumChannels() {
    		return getNumChannels(0);
    	}
    	public int getNumChannels(int numStation) {
    		return this.eyesisSubCameras[numStation].length;
    	}
//    	this.eyesisSubCameras[numStation].length
    	public EyesisCameraParameters () {} // just create new instance, all parameters data will be provided additionally
    	public EyesisCameraParameters (
    			int numStations,
    			boolean isTripod,
    	    	double goniometerHorizontal, // goniometer rotation around "horizontal" axis (tilting from the target - positive)
    	    	double goniometerAxial, // goniometer rotation around Eyesis axis (clockwise in plan - positive 
    			int numSubCameras,
    	    	double interAxisDistance, // distance in mm between two goniometer axes, positive if the vertical axis (when Eyesis is head up) is closer to the target
    	    	double interAxisAngle,    // angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
    	    	                            // clockwise when eyesis is in 'normal' position, looking to the target
    	    	double horAxisErrPhi,   // angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)
    	    	double horAxisErrPsi,   // angle in degrees "horizontal" goniometer axis is rotated around moving Z axis (CW looking at target)
    	    	double entrancePupilForward, // common to all lenses - distance from the sensor to the lens entrance pupil
    	    	double centerAboveHorizontal, // camera center distance along camera axis above the closest point to horizontal rotation axis (adds to height of each 
    	    	double GXYZ_0, // coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system
    	    	double GXYZ_1, // coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system
    	    	double GXYZ_2, // coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system
    	    	int sensorWidth,
    	    	int sensorHeight,
    	    	int    shrinkGridForMask, //shrink detected grids by one point for/vert this number of times before calculating masks
    	    	double maskBlurSigma,      // blur sensor masks (in grid units)
    	    	int    decimateMasks,       // reduce masks resolution
    	    	double badNodeThreshold, // filter out grid nodes with difference from quadratically predicted from 8 neighbors in pixels
    	    	int    maxBadNeighb, // maximal number of bad nodes around the corrected one to fix
    	    	int minimalValidNodes,
    	    	int    weightMultiImageMode, // increase weight for multi-image sets (0 - do not increase, 1 - multiply by number of images in a set)
    	    	double  weightMultiExponent, 
    	    	double  weightDiameterExponent, // if( >0) use grid diameter to scale weights of this image
    	    	double weightYtoX,
    	    	double minimalGridContrast,
    	    	double shrinkBlurSigma,
    	    	double shrinkBlurLevel,
    	    	double balanceChannelWeightsMode,
    	    	double removeOverRMS,
    	    	double removeOverRMSNonweighted
    			){
    		double [] GXYZ={GXYZ_0,GXYZ_1,GXYZ_2};
    		setSameEyesisCameraParameters (
    				numStations,
    				isTripod,
    				goniometerHorizontal, // goniometer rotation around "horizontal" axis (tilting from the target - positive)
    				goniometerAxial, // goniometer rotation around Eyesis axis (clockwise in plan - positive 
    				numSubCameras,
    				interAxisDistance, // distance in mm between two goniometer axes
    				interAxisAngle,    // angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
    				// clockwise when eyesis is in 'normal' position, looking to the target
    				horAxisErrPhi,   // angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)
    				horAxisErrPsi,   // angle in degrees "horizontal" goniometer axis is rotated moving Z axis (CW looking at target)
    				entrancePupilForward, // common to all lenses - distance from the sensor to the lens entrance pupil
    				centerAboveHorizontal, // camera center distance along camera axis above the closest point to horizontal rotation axis (adds to height of each 
    				GXYZ, // coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system
    				sensorWidth,
    				sensorHeight,
    				shrinkGridForMask, //shrink detected grids by one point for/vert this number of times before calculating masks
    				maskBlurSigma,      // blur sensor masks (in grid units)
    				decimateMasks,       // reduce masks resolution
    				badNodeThreshold, // filter out grid nodes with difference from quadratically predicted from 8 neighbors in pixels
    				maxBadNeighb, // maximal number of bad nodes around the corrected one to fix
    				minimalValidNodes,
    				weightMultiImageMode, // increase weight for multi-image sets (0 - do not increase, 1 - multiply by number of images in a set)
        	    	weightMultiExponent, 
        	    	weightDiameterExponent, // if( >0) use grid diameter to scale weights of this image
        	    	weightYtoX,
        	    	minimalGridContrast,
        	    	shrinkBlurSigma,
    		    	shrinkBlurLevel,
    		    	balanceChannelWeightsMode,
        	    	removeOverRMS,
        	    	removeOverRMSNonweighted
    		);
    	}
    	public EyesisCameraParameters (
    			int numStations,
    			boolean isTripod,
    	    	double goniometerHorizontal, // goniometer rotation around "horizontal" axis (tilting from the target - positive)
    	    	double goniometerAxial, // goniometer rotation around Eyesis axis (clockwise in plan - positive 
    			int numSubCameras,
    	    	double interAxisDistance, // distance in mm between two goniometer axes
    	    	double interAxisAngle,    // angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
    	    	                            // clockwise when eyesis is in 'normal' position, looking to the target
    	    	double horAxisErrPhi,   // angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)
    	    	double horAxisErrPsi,   // angle in degrees "horizontal" goniometer axis is rotated moving Z axis (CW looking at target)
    	    	double entrancePupilForward, // common to all lenses - distance from the sensor to the lens entrance pupil
    	    	double centerAboveHorizontal, // camera center distance along camera axis above the closest point to horizontal rotation axis (adds to height of each 
    	    	double [] GXYZ, // coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system
    	    	int sensorWidth,
    	    	int sensorHeight,
    	    	int    shrinkGridForMask, //shrink detected grids by one point for/vert this number of times before calculating masks
    	    	double maskBlurSigma,      // blur sensor masks (in grid units)
    	    	int    decimateMasks,       // reduce masks resolution
    	    	double badNodeThreshold, // filter out grid nodes with difference from quadratically predicted from 8 neighbors in pixels
    	    	int    maxBadNeighb, // maximal number of bad nodes around the corrected one to fix
    	    	int    minimalValidNodes,
    	    	int    weightMultiImageMode, // increase weight for multi-image sets (0 - do not increase, 1 - multiply by number of images in a set)
    	    	double  weightMultiExponent, 
    	    	double  weightDiameterExponent, // if( >0) use grid diameter to scale weights of this image
    	    	double weightYtoX,
    	    	double minimalGridContrast,
    	    	double shrinkBlurSigma,
    	    	double shrinkBlurLevel,
    	    	double balanceChannelWeightsMode,
    	    	double removeOverRMS,
    	    	double removeOverRMSNonweighted
    			){
    		setSameEyesisCameraParameters (
    				numStations,
    				isTripod,
    				goniometerHorizontal, // goniometer rotation around "horizontal" axis (tilting from the target - positive)
    				goniometerAxial, // goniometer rotation around Eyesis axis (clockwise in plan - positive 
    				numSubCameras,
    				interAxisDistance, // distance in mm between two goniometer axes
    				interAxisAngle,    // angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
    				// clockwise when eyesis is in 'normal' position, looking to the target
    				horAxisErrPhi,   // angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)
    				horAxisErrPsi,   // angle in degrees "horizontal" goniometer axis is rotated moving Z axis (CW looking at target)
    				entrancePupilForward, // common to all lenses - distance from the sensor to the lens entrance pupil
    				centerAboveHorizontal, // camera center distance along camera axis above the closest point to horizontal rotation axis (adds to height of each 
    				GXYZ, // coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system
    				sensorWidth,
    				sensorHeight,
    				shrinkGridForMask, //shrink detected grids by one point for/vert this number of times before calculating masks
    				maskBlurSigma,      // blur sensor masks (in grid units)
    				decimateMasks,       // reduce masks resolution
    				badNodeThreshold, // filter out grid nodes with difference from quadratically predicted from 8 neighbors in pixels
    				maxBadNeighb, // maximal number of bad nodes around the corrected one to fix
    				minimalValidNodes,
    				weightMultiImageMode, // increase weight for multi-image sets (0 - do not increase, 1 - multiply by number of images in a set)
        	    	weightMultiExponent, 
        	    	weightDiameterExponent, // if( >0) use grid diameter to scale weights of this image
        	    	weightYtoX,
        	    	minimalGridContrast,
        	    	shrinkBlurSigma,
    		    	shrinkBlurLevel,
    		    	balanceChannelWeightsMode,
        	    	removeOverRMS,
        	    	removeOverRMSNonweighted
    		);
    	}
    	void setSameEyesisCameraParameters (
    			int numStations,
    			boolean isTripod,
    	    	double goniometerHorizontal, // goniometer rotation around "horizontal" axis (tilting from the target - positive)
    	    	double goniometerAxial, // goniometer rotation around Eyesis axis (clockwise in plan - positive 
    			int numSubCameras,
    	    	double interAxisDistance, // distance in mm between two goniometer axes
    	    	double interAxisAngle,    // angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
    	    	                            // clockwise when eyesis is in 'normal' position, looking to the target
    	    	double horAxisErrPhi,   // angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)
    	    	double horAxisErrPsi,   // angle in degrees "horizontal" goniometer axis is rotated moving Z axis (CW looking at target)
    	    	double entrancePupilForward, // common to all lenses - distance from the sensor to the lens entrance pupil
    	    	double centerAboveHorizontal, // camera center distance along camera axis above the closest point to horizontal rotation axis (adds to height of each 
    	    	double [] GXYZ, // coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system
    	    	int sensorWidth,
    	    	int sensorHeight,
    	    	int    shrinkGridForMask, //shrink detected grids by one point for/vert this number of times before calculating masks
    	    	double maskBlurSigma,      // blur sensor masks (in grid units)
    	    	int    decimateMasks,       // reduce masks resolution
    	    	double badNodeThreshold, // filter out grid nodes with difference from quadratically predicted from 8 neighbors in pixels
    	    	int    maxBadNeighb, // maximal number of bad nodes around the corrected one to fix
    	    	int    minimalValidNodes,
    	    	int    weightMultiImageMode, // increase weight for multi-image sets (0 - do not increase, 1 - multiply by number of images in a set)
    	    	double  weightMultiExponent, 
    	    	double  weightDiameterExponent, // if( >0) use grid diameter to scale weights of this image
    	    	double weightYtoX,
    	    	double minimalGridContrast,
    	    	double shrinkBlurSigma,
    	    	double shrinkBlurLevel,
    	    	double balanceChannelWeightsMode,
    	    	double removeOverRMS,
    	    	double removeOverRMSNonweighted
    			){
    		this.numStations=numStations;
    		this.isTripod=isTripod;
	    	this.sensorWidth=sensorWidth;
	    	this.sensorHeight=sensorHeight;
	    	this.shrinkGridForMask=shrinkGridForMask; //shrink detected grids by one point for/vert this number of times before calculating masks
	    	this.maskBlurSigma=maskBlurSigma;      // blur sensor masks (in grid units)
	    	this.decimateMasks=decimateMasks;
	    	this.badNodeThreshold=badNodeThreshold; // filter out grid nodes with difference from quadratically predicted from 8 neighbors in pixels
	    	this.maxBadNeighb=maxBadNeighb; // maximal number of bad nodes around the corrected one to fix
	    	this.minimalValidNodes=minimalValidNodes;
	    	this.weightMultiImageMode=weightMultiImageMode; // increase weight for multi-image sets (0 - do not increase, 1 - multiply by number of images in a set)
	    	this.weightMultiExponent=weightMultiExponent; 
	    	this.weightDiameterExponent=weightDiameterExponent; // if( >0) use grid diameter to scale weights of this image
	    	this.weightYtoX=weightYtoX;
	    	this.minimalGridContrast=minimalGridContrast;
	    	this.goniometerHorizontal=new double[numStations];
    		this.goniometerAxial=new double[numStations]; 
	    	this.interAxisDistance=new double[numStations];
	    	this.interAxisAngle=new double[numStations];
	    	this.horAxisErrPhi=new double[numStations];
	    	this.horAxisErrPsi=new double[numStations];
	    	this.entrancePupilForward=new double[numStations]; // common to all lenses - distance from the sensor to the lens entrance pupil
	    	this.centerAboveHorizontal=new double[numStations]; // camera center distance along camera axis above the closest point to horizontal rotation axis (adds to height of each 
	    	this.GXYZ=new double[numStations][];
	    	if (numSubCameras>0) this.eyesisSubCameras=new EyesisSubCameraParameters[numStations][];
	    	this.stationWeight=new double[numStations];
	    	for (int numStation=0;numStation<numStations;numStation++) {
	    		this.stationWeight[numStation]=1.0;
	    		this.goniometerHorizontal[numStation]=goniometerHorizontal;
	    		this.goniometerAxial[numStation]=goniometerAxial; 
		    	this.interAxisDistance[numStation]=interAxisDistance;
		    	this.interAxisAngle[numStation]=interAxisAngle;
		    	this.horAxisErrPhi[numStation]=horAxisErrPhi;
		    	this.horAxisErrPsi[numStation]=horAxisErrPsi;
		    	this.entrancePupilForward[numStation]=entrancePupilForward; // common to all lenses - distance from the sensor to the lens entrance pupil
		    	this.centerAboveHorizontal[numStation]=centerAboveHorizontal; // camera center distance along camera axis above the closest point to horizontal rotation axis (adds to height of each
		    	this.GXYZ[numStation]=new double [3];
		    	for (int i=0;i<3;i++) this.GXYZ[numStation][i]=GXYZ[i];
		    	if (numSubCameras>0) initSubCameras(numStation,numSubCameras);
	    	}
    	}
    	
    	public EyesisCameraParameters clone() {
    		EyesisCameraParameters result= new EyesisCameraParameters ();
    		copyData(
        			this.numStations,
        			this,
        			result);
    		return result;
    	}

    	
    	/**
    	 * Capoy parameters from source EysesisCamerParameters to destination, trimming/expanding nu,ber of stations
    	 * @param newNumStations new number of stations
    	 * @param source source EysesisCamerParameters
    	 * @param destination destination EysesisCamerParameters
    	 */
    	public void copyData(
    			int newNumStations,
    			EyesisCameraParameters source,
    			EyesisCameraParameters destination) {
    		destination.numStations=newNumStations;
    		destination.isTripod=source.isTripod;
    		destination.sensorWidth=source.sensorWidth;
    		destination.sensorHeight=source.sensorHeight;
    		destination.shrinkGridForMask=source.shrinkGridForMask; //shrink detected grids by one point for/vert this number of times before calculating masks
    		destination.maskBlurSigma=source.maskBlurSigma;      // blur sensor masks (in grid units)
    		destination.decimateMasks=source.decimateMasks;
    		destination.badNodeThreshold=source.badNodeThreshold; // filter out grid nodes with difference from quadratically predicted from 8 neighbors in pixels
    		destination.maxBadNeighb=source.maxBadNeighb; // maximal number of bad nodes around the corrected one to fix
    		destination.minimalValidNodes=source.minimalValidNodes;
    		destination.weightMultiImageMode=source.weightMultiImageMode; // increase weight for multi-image sets (0 - do not increase, 1 - multiply by number of images in a set)
    		destination.weightMultiExponent= source.weightMultiExponent; // if( >0) use grid diameter to scale weights of this image
    		destination.weightDiameterExponent=source.weightDiameterExponent;
    		destination.weightYtoX=source.weightYtoX;
	    	destination.minimalGridContrast=source.minimalGridContrast;
	    	destination.shrinkBlurSigma=source.shrinkBlurSigma;
	    	destination.shrinkBlurLevel=source.shrinkBlurLevel;
	    	destination.balanceChannelWeightsMode=source.balanceChannelWeightsMode;
	    	destination.removeOverRMSNonweighted=source.removeOverRMSNonweighted;
    		destination.goniometerHorizontal=new double[destination.numStations];
    		destination.goniometerAxial=new double[destination.numStations]; 
    		destination.interAxisDistance=new double[destination.numStations];
    		destination.interAxisAngle=new double[destination.numStations];
    		destination.horAxisErrPhi=new double[destination.numStations];
    		destination.horAxisErrPsi=new double[destination.numStations];
    		destination.entrancePupilForward=new double[destination.numStations]; // common to all lenses - distance from the sensor to the lens entrance pupil
    		destination.centerAboveHorizontal=new double[destination.numStations]; // camera center distance along camera axis above the closest point to horizontal rotation axis (adds to height of each 
    		destination.GXYZ=new double[destination.numStations][];
    		if (source.eyesisSubCameras!=null) destination.eyesisSubCameras=new EyesisSubCameraParameters[destination.numStations][];
    		else destination.eyesisSubCameras=null;
    		destination.stationWeight=new double[destination.numStations];
    		for (int numStation=0;numStation<destination.numStations;numStation++) {
    			int srcNumStation=(numStation<source.numStations)?numStation:(source.numStations-1);
    			destination.stationWeight[numStation]=source.stationWeight[srcNumStation];
    			destination.goniometerHorizontal[numStation]=source.goniometerHorizontal[srcNumStation];
    			destination.goniometerAxial[numStation]=source.goniometerAxial[srcNumStation]; 
    			destination.interAxisDistance[numStation]=source.interAxisDistance[srcNumStation];
    			destination.interAxisAngle[numStation]=source.interAxisAngle[srcNumStation];
    			destination.horAxisErrPhi[numStation]=source.horAxisErrPhi[srcNumStation];
    			destination.horAxisErrPsi[numStation]=source.horAxisErrPsi[srcNumStation];
    			destination.entrancePupilForward[numStation]=source.entrancePupilForward[srcNumStation]; // common to all lenses - distance from the sensor to the lens entrance pupil
    			destination.centerAboveHorizontal[numStation]=source.centerAboveHorizontal[srcNumStation]; // camera center distance along camera axis above the closest point to horizontal rotation axis (adds to height of each
    			destination.GXYZ[numStation]=new double [3];
    			for (int i=0;i<3;i++) destination.GXYZ[numStation][i]=source.GXYZ[srcNumStation][i];
    			if (destination.eyesisSubCameras!=null){
    				destination.eyesisSubCameras[numStation]=new EyesisSubCameraParameters[source.eyesisSubCameras[srcNumStation].length];
    				for (int i=0; i<destination.eyesisSubCameras[numStation].length;i++) if (source.eyesisSubCameras[srcNumStation][i]!=null){
    					destination.eyesisSubCameras[numStation][i]=source.eyesisSubCameras[srcNumStation][i].clone();
    				}
    			}
    		}
    	}
    	
    	public void setProperties(String prefix,Properties properties){
    		properties.setProperty(prefix+"isTripod",this.isTripod+"");
    		properties.setProperty(prefix+"sensorWidth",this.sensorWidth+"");
    		properties.setProperty(prefix+"sensorHeight",this.sensorHeight+"");
    		properties.setProperty(prefix+"shrinkGridForMask",this.shrinkGridForMask+"");
    		properties.setProperty(prefix+"maskBlurSigma",this.maskBlurSigma+"");
    		properties.setProperty(prefix+"decimateMasks",this.decimateMasks+"");
    		properties.setProperty(prefix+"badNodeThreshold",this.badNodeThreshold+"");
    		properties.setProperty(prefix+"maxBadNeighb",this.maxBadNeighb+"");
    		properties.setProperty(prefix+"minimalValidNodes",this.minimalValidNodes+"");
    		properties.setProperty(prefix+"weightMultiImageMode",this.weightMultiImageMode+"");
    		properties.setProperty(prefix+"weightMultiExponent",this.weightMultiExponent+"");
    		properties.setProperty(prefix+"weightDiameterExponent",this.weightDiameterExponent+"");
    		properties.setProperty(prefix+"weightYtoX",this.weightYtoX+"");
    		properties.setProperty(prefix+"minimalGridContrast",this.minimalGridContrast+"");
    		properties.setProperty(prefix+"shrinkBlurSigma",this.shrinkBlurSigma+"");
    		properties.setProperty(prefix+"shrinkBlurLevel",this.shrinkBlurLevel+"");
    		properties.setProperty(prefix+"balanceChannelWeightsMode",this.balanceChannelWeightsMode+"");
    		properties.setProperty(prefix+"removeOverRMS",this.removeOverRMS+"");
    		properties.setProperty(prefix+"removeOverRMSNonweighted",this.removeOverRMSNonweighted+"");
			properties.setProperty(prefix+"numSubCameras",this.eyesisSubCameras[0].length+"");
    		properties.setProperty(prefix+"numStations",this.numStations+"");
    		for (int numStation=0;numStation<this.numStations;numStation++){
    			properties.setProperty(prefix+"stationWeight_"+numStation,this.stationWeight[numStation]+"");
    			properties.setProperty(prefix+"goniometerHorizontal_"+numStation,this.goniometerHorizontal[numStation]+"");
    			properties.setProperty(prefix+"goniometerAxial_"+numStation,this.goniometerAxial[numStation]+"");
    			properties.setProperty(prefix+"interAxisDistance_"+numStation,this.interAxisDistance[numStation]+"");
    			properties.setProperty(prefix+"interAxisAngle_"+numStation,this.interAxisAngle[numStation]+"");
    			properties.setProperty(prefix+"horAxisErrPhi_"+numStation,this.horAxisErrPhi[numStation]+"");
    			properties.setProperty(prefix+"horAxisErrPsi_"+numStation,this.horAxisErrPsi[numStation]+"");
    			properties.setProperty(prefix+"entrancePupilForward_"+numStation,this.entrancePupilForward[numStation]+"");
    			properties.setProperty(prefix+"centerAboveHorizontal_"+numStation,this.centerAboveHorizontal[numStation]+"");
    			properties.setProperty(prefix+"GXYZ_0_"+numStation,this.GXYZ[numStation][0]+"");
    			properties.setProperty(prefix+"GXYZ_1_"+numStation,this.GXYZ[numStation][1]+"");
    			properties.setProperty(prefix+"GXYZ_2_"+numStation,this.GXYZ[numStation][2]+"");
    			for (int i=0;i<this.eyesisSubCameras[numStation].length;i++) {
    				this.eyesisSubCameras[numStation][i].setProperties(prefix+numStation+"_subCamera_"+i+'.',properties);
    			}
    		}
 //   		setCostsProperties(prefix,properties);
    	}
    	public void getProperties(String prefix,Properties properties){
    		if (properties.getProperty(prefix+"isTripod")!=null)
    			this.isTripod=Boolean.parseBoolean(properties.getProperty(prefix+"isTripod"));
    		if (properties.getProperty(prefix+"sensorWidth")!=null)
    			this.sensorWidth=Integer.parseInt(properties.getProperty(prefix+"sensorWidth"));
    		if (properties.getProperty(prefix+"sensorHeight")!=null)
    			this.sensorHeight=Integer.parseInt(properties.getProperty(prefix+"sensorHeight"));
    		if (properties.getProperty(prefix+"shrinkGridForMask")!=null)
    			this.shrinkGridForMask=Integer.parseInt(properties.getProperty(prefix+"shrinkGridForMask"));
    		if (properties.getProperty(prefix+"maskBlurSigma")!=null)
    			this.maskBlurSigma=Double.parseDouble(properties.getProperty(prefix+"maskBlurSigma"));
    		if (properties.getProperty(prefix+"decimateMasks")!=null)
    			this.decimateMasks=Integer.parseInt(properties.getProperty(prefix+"decimateMasks"));

    		if (properties.getProperty(prefix+"badNodeThreshold")!=null)
    			this.badNodeThreshold=Double.parseDouble(properties.getProperty(prefix+"badNodeThreshold"));
    		if (properties.getProperty(prefix+"maxBadNeighb")!=null)
    			this.maxBadNeighb=Integer.parseInt(properties.getProperty(prefix+"maxBadNeighb"));
    		if (properties.getProperty(prefix+"minimalValidNodes")!=null)
    			this.minimalValidNodes=Integer.parseInt(properties.getProperty(prefix+"minimalValidNodes"));
    		if (properties.getProperty(prefix+"weightMultiImageMode")!=null)
    			this.weightMultiImageMode=Integer.parseInt(properties.getProperty(prefix+"weightMultiImageMode"));
    		if (properties.getProperty(prefix+"weightMultiExponent")!=null)
    			this.weightMultiExponent=Double.parseDouble(properties.getProperty(prefix+"weightMultiExponent"));
    		if (properties.getProperty(prefix+"weightDiameterExponent")!=null)
    			this.weightDiameterExponent=Double.parseDouble(properties.getProperty(prefix+"weightDiameterExponent"));
    		if (properties.getProperty(prefix+"weightYtoX")!=null)
    			this.weightYtoX=Double.parseDouble(properties.getProperty(prefix+"weightYtoX"));
    		if (properties.getProperty(prefix+"minimalGridContrast")!=null)
    			this.minimalGridContrast=Double.parseDouble(properties.getProperty(prefix+"minimalGridContrast"));
    		if (properties.getProperty(prefix+"shrinkBlurSigma")!=null)
    			this.shrinkBlurSigma=Double.parseDouble(properties.getProperty(prefix+"shrinkBlurSigma"));
    		if (properties.getProperty(prefix+"shrinkBlurLevel")!=null)
    			this.shrinkBlurLevel=Double.parseDouble(properties.getProperty(prefix+"shrinkBlurLevel"));
    		if (properties.getProperty(prefix+"balanceChannelWeightsMode")!=null)
    			this.balanceChannelWeightsMode=Double.parseDouble(properties.getProperty(prefix+"balanceChannelWeightsMode"));
    		if (properties.getProperty(prefix+"removeOverRMS")!=null)
    			this.removeOverRMS=Double.parseDouble(properties.getProperty(prefix+"removeOverRMS"));
    		if (properties.getProperty(prefix+"removeOverRMSNonweighted")!=null)
    			this.removeOverRMSNonweighted=Double.parseDouble(properties.getProperty(prefix+"removeOverRMSNonweighted"));
    		boolean multiStation=true; // new default
    		int newNumStations=1;
    		int numSubCameras=0;
    		if (properties.getProperty(prefix+"numSubCameras")!=null) numSubCameras=Integer.parseInt(properties.getProperty(prefix+"numSubCameras"));
    		if (properties.getProperty(prefix+"numStations")!=null){
    			newNumStations=Integer.parseInt(properties.getProperty(prefix+"numStations"));
    		}else {
    			multiStation=false; // old config format
    		}
// TODO: trim/expand stations    		
    		updateNumstations (newNumStations);
//    		this.numStations
    		
// read old/new format data
//    		this.numStations=newNumStations;
    		if (multiStation){
    			for (int numStation=0;numStation<this.numStations;numStation++){
        			if (properties.getProperty(prefix+"stationWeight_"+numStation)!=null)
        				this.stationWeight[numStation]=Double.parseDouble(properties.getProperty(prefix+"stationWeight_"+numStation));
        			if (properties.getProperty(prefix+"goniometerHorizontal_"+numStation)!=null)
        				this.goniometerHorizontal[numStation]=Double.parseDouble(properties.getProperty(prefix+"goniometerHorizontal_"+numStation));
        			if (properties.getProperty(prefix+"goniometerAxial_"+numStation)!=null)
        				this.goniometerAxial[numStation]=Double.parseDouble(properties.getProperty(prefix+"goniometerAxial_"+numStation));
        			if (properties.getProperty(prefix+"interAxisDistance_"+numStation)!=null)
        				this.interAxisDistance[numStation]=Double.parseDouble(properties.getProperty(prefix+"interAxisDistance_"+numStation));
        			if (properties.getProperty(prefix+"interAxisAngle_"+numStation)!=null)
        				this.interAxisAngle[numStation]=Double.parseDouble(properties.getProperty(prefix+"interAxisAngle_"+numStation));
        			if (properties.getProperty(prefix+"horAxisErrPhi_"+numStation)!=null)
        				this.horAxisErrPhi[numStation]=Double.parseDouble(properties.getProperty(prefix+"horAxisErrPhi_"+numStation));
        			if (properties.getProperty(prefix+"horAxisErrPsi_"+numStation)!=null)
        				this.horAxisErrPsi[numStation]=Double.parseDouble(properties.getProperty(prefix+"horAxisErrPsi_"+numStation));
        			if (properties.getProperty(prefix+"entrancePupilForward_"+numStation)!=null)
        				this.entrancePupilForward[numStation]=Double.parseDouble(properties.getProperty(prefix+"entrancePupilForward_"+numStation));
        			if (properties.getProperty(prefix+"centerAboveHorizontal_"+numStation)!=null)
        				this.centerAboveHorizontal[numStation]=Double.parseDouble(properties.getProperty(prefix+"centerAboveHorizontal_"+numStation));
        			if (properties.getProperty(prefix+"GXYZ_0_"+numStation)!=null)
        				this.GXYZ[numStation][0]=Double.parseDouble(properties.getProperty(prefix+"GXYZ_0_"+numStation));
        			if (properties.getProperty(prefix+"GXYZ_1_"+numStation)!=null)
        				this.GXYZ[numStation][1]=Double.parseDouble(properties.getProperty(prefix+"GXYZ_1_"+numStation));
        			if (properties.getProperty(prefix+"GXYZ_2_"+numStation)!=null)
        				this.GXYZ[numStation][2]=Double.parseDouble(properties.getProperty(prefix+"GXYZ_2_"+numStation));
        			if (numSubCameras>0) {
        				initSubCameras(numStation, numSubCameras); // set array with default parameters
        				for (int i=0;i<numSubCameras;i++){
        					this.eyesisSubCameras[numStation][i].getProperties(prefix+numStation+"_subCamera_"+i+'.',properties,i);
        				}
        			}
    			}
    		} else { // read old format
    			if (properties.getProperty(prefix+"goniometerHorizontal")!=null)
    				this.goniometerHorizontal[0]=Double.parseDouble(properties.getProperty(prefix+"goniometerHorizontal"));
    			if (properties.getProperty(prefix+"goniometerAxial")!=null)
    				this.goniometerAxial[0]=Double.parseDouble(properties.getProperty(prefix+"goniometerAxial"));
    			if (properties.getProperty(prefix+"interAxisDistance")!=null)
    				this.interAxisDistance[0]=Double.parseDouble(properties.getProperty(prefix+"interAxisDistance"));
    			if (properties.getProperty(prefix+"interAxisAngle")!=null)
    				this.interAxisAngle[0]=Double.parseDouble(properties.getProperty(prefix+"interAxisAngle"));
    			if (properties.getProperty(prefix+"horAxisErrPhi")!=null)
    				this.horAxisErrPhi[0]=Double.parseDouble(properties.getProperty(prefix+"horAxisErrPhi"));
    			if (properties.getProperty(prefix+"horAxisErrPsi")!=null)
    				this.horAxisErrPsi[0]=Double.parseDouble(properties.getProperty(prefix+"horAxisErrPsi"));
    			if (properties.getProperty(prefix+"entrancePupilForward")!=null)
    				this.entrancePupilForward[0]=Double.parseDouble(properties.getProperty(prefix+"entrancePupilForward"));
    			if (properties.getProperty(prefix+"centerAboveHorizontal")!=null)
    				this.centerAboveHorizontal[0]=Double.parseDouble(properties.getProperty(prefix+"centerAboveHorizontal"));
    			if (properties.getProperty(prefix+"GXYZ_0")!=null)
    				this.GXYZ[0][0]=Double.parseDouble(properties.getProperty(prefix+"GXYZ_0"));
    			if (properties.getProperty(prefix+"GXYZ_1")!=null)
    				this.GXYZ[0][1]=Double.parseDouble(properties.getProperty(prefix+"GXYZ_1"));
    			if (properties.getProperty(prefix+"GXYZ_2")!=null)
    				this.GXYZ[0][2]=Double.parseDouble(properties.getProperty(prefix+"GXYZ_2"));
    			System.out.println(" === numSubCameras="+numSubCameras);
    			if (numSubCameras>0) {
    				initSubCameras(0, numSubCameras); // set array with default parameters
    				for (int i=0;i<numSubCameras;i++){
    					System.out.println("this.eyesisSubCameras[0]["+i+"].getProperties("+prefix+"subCamera_"+i+".,properties);");
    					this.eyesisSubCameras[0][i].getProperties(prefix+"subCamera_"+i+'.',properties);
    				}
    			}
    		}
//    		getCostsProperties(prefix,properties);
    	}
    	void  updateNumstations (int newNumStations){
//    		System.out.println("updateNumstations("+newNumStations+"), was "+this.numStations);
    		if (newNumStations==this.numStations) return;
    		System.out.println("updateNumstations("+newNumStations+"), was "+this.numStations);
    		EyesisCameraParameters newData=this.clone();
    		copyData(
    				newNumStations,
    				newData,
        			this);
    		System.out.println("updateNumstations() after copyData() this.numStations="+this.numStations);
    	}
    	
    	// returns -1 - canceled, 0 - done, >0 - number of sub-camera to edit

    	private int subShowDialog(String title, int nextSubCamera) {
    		GenericDialog gd = new GenericDialog(title);
//        	public double goniometerHorizontal; // goniometer rotation around "horizontal" axis (tilting from the target - positive)
//        	public double goniometerAxial; // goniometer rotation around Eyesis axis (clockwise in plan - positive 
//this.isTripod
    		String [] modelChoice=new String [distortionModelDescriptions.length+1];
    		modelChoice[0]="--- keep current ---";
    		for (int i=0;i<distortionModelDescriptions.length;i++) modelChoice[i+1]=distortionModelDescriptions[i];
    		gd.addChoice("Change camera distortion model for all channels", modelChoice, modelChoice[0]);
    		gd.addCheckbox("Tripod mode (first vertical axis, then horizontal), changes meaning of the next 2 fields",this.isTripod);
    		gd.addMessage("=== Camera parameters to be fitted ===");
    		for (int numStation=0;numStation<this.numStations;numStation++) {
    			if (this.numStations>1){
    				gd.addMessage("--- Station number "+numStation+" ---");
    				gd.addNumericField("Station reprojection errors weight",      100*this.stationWeight[numStation], 1,5,"%");
    			}
    			if (this.isTripod) {
    				gd.addNumericField("Tripod rotation around vertical axis (clockwise from top - positive)",this.goniometerHorizontal[numStation], 3,7,"degrees");
    				gd.addNumericField("Tripod rotation around horizontal axis (camera up - positive)",       this.goniometerAxial[numStation], 3,7,"degrees");
    			} else {
    				gd.addNumericField("Goniometer rotation around 'horizontal' axis (tilting from the target - positive)",this.goniometerHorizontal[numStation], 3,7,"degrees");
    				gd.addNumericField("Rotation around Eyesis main axis (clockwise in plan - positive)",      this.goniometerAxial[numStation], 3,7,"degrees");
    			}

    			gd.addNumericField("Distance between goniometer axes",                                     this.interAxisDistance[numStation], 3,7,"mm");
    			gd.addNumericField("Angle error between goniometer axes (<0 if vertical axis rotated CW )",this.interAxisAngle[numStation],  3,7,"degrees");
    			if (this.isTripod) {
    				gd.addNumericField("Vertical tripod axis tilt  from true vertical",                    this.horAxisErrPhi[numStation],  3,7,"degrees");
    				gd.addNumericField("Vertical tripod axis roll error from true vertical",               this.horAxisErrPsi[numStation],  3,7,"degrees");
    			} else {
    				gd.addNumericField("Horizontal axis azimuth error (CW in plan)",                       this.horAxisErrPhi[numStation],  3,7,"degrees");
    				gd.addNumericField("Horizontal axis roll error (CW looking to target)",                this.horAxisErrPsi[numStation],  3,7,"degrees");
    			}

    			gd.addNumericField("Distance between the lens entrace pupil and the sensor",               this.entrancePupilForward[numStation],  3,7,"mm");
    			gd.addNumericField("Camera center above goniometer horizontal axis",                       this.centerAboveHorizontal[numStation],  3,7,"mm");

    			gd.addNumericField("Goniometer reference point position X (target coordinates, left)",     this.GXYZ[numStation][0], 3,7,"mm");
    			gd.addNumericField("Goniometer reference point position Y (target coordinates, up)",       this.GXYZ[numStation][1], 3,7,"mm");
    			gd.addNumericField("Goniometer reference point position Z (target coordinates, away)",     this.GXYZ[numStation][2], 3,7,"mm");
    		}
    		gd.addMessage("=== Other parameters ===");

    		gd.addNumericField("Image sensor width (maximal if different)",                            this.sensorWidth, 0,4,"pix");
    		gd.addNumericField("Image sensor height (maximal if different)",                           this.sensorHeight, 0,4,"pix");
    		gd.addNumericField("Shrink detected grid by this number of nodes (half/periods) for masks",this.shrinkGridForMask, 0,4,"grid nodes");
    		gd.addNumericField("Gaussian blur masks for the sensors (positive - pixels, negative - grid half-periods)", this.maskBlurSigma, 2,6,"pix");
    		gd.addNumericField("Reduce sensor resolution when calculating masks",                      this.decimateMasks, 0);

    		gd.addNumericField("Filter out grid nodes with difference from quadratically predicted from 8 neighbors", this.badNodeThreshold, 2,6,"pix");
    		gd.addNumericField("Maximal number of bad nodes around the corrected one to fix",          this.maxBadNeighb, 0);
    		gd.addNumericField("Minimal number of valid (with all filters applied) nodes in each image",this.minimalValidNodes, 0);
    		gd.addNumericField("Increase weight of the multi-image sets (0 - do not increase, 1 - multiply by number of images in a set (to power ), 2 - same but remove single-image ",this.weightMultiImageMode, 0);
    		gd.addNumericField("Increase weight of the multi-image sets power (used with mode above)", this.weightMultiExponent, 2,6,"");
    		gd.addNumericField("Increase weight of the images by the power of their diameters", this.weightDiameterExponent, 2,6,"");
    		gd.addNumericField("Increase weight Y-error with respect to X-error",                100.0*this.weightYtoX, 2,6,"%");
    		gd.addNumericField("Filter out grid nodes with the contrast less than this value (maximal is ~0.8) ", this.minimalGridContrast, 2,6,"");
    		gd.addNumericField("Shrink-blur detected grids alpha",                                     this.shrinkBlurSigma, 2,6,"grid nodes");
    		gd.addNumericField("Shrink-blur detected grids  level (-1..+1)",                           this.shrinkBlurLevel, 2,6,"");
    		gd.addNumericField("Channel balace mode: <0 - use specified defaults, 0 - keep curent, >0 - exponent for correction (1.0 - precise equalization)", this.balanceChannelWeightsMode, 3,6,"");
    		gd.addNumericField("Remove nodes with error greater than scaled RMS in that image, weighted",  this.removeOverRMS, 2,6,"xRMS");
    		gd.addNumericField("Same, not weghted (not more permissive near the borders with low weight)", this.removeOverRMSNonweighted, 2,6,"xRMS");
    		gd.addNumericField("Number of sub-camera modules",                                         this.eyesisSubCameras[0].length,   0,2,"");
    		gd.addNumericField("Number of sub-camera module to edit (<=0  - none)",                    nextSubCamera,   0,2,"");
   	        gd.enableYesNoCancel("OK", "Done");
    	    WindowTools.addScrollBars(gd);
    		gd.showDialog();
    		if (gd.wasCanceled()) return -1;
    		int modelIndex=gd.getNextChoiceIndex()-1;
    		if (modelIndex>=0){
    			for (EyesisSubCameraParameters [] esps:this.eyesisSubCameras){
    				for (EyesisSubCameraParameters esp:esps){
    					esp.lensDistortionModel=distortionModels[modelIndex];
    				}
    			}
    		}

    		this.isTripod=                 gd.getNextBoolean();
    		for (int numStation=0;numStation<this.numStations;numStation++) {
    			if (this.numStations>1){
    				this.stationWeight[numStation]=0.01*gd.getNextNumber();
    			}
    			this.goniometerHorizontal[numStation]=     gd.getNextNumber();
    			this.goniometerAxial[numStation]=          gd.getNextNumber();
    			this.interAxisDistance[numStation]=        gd.getNextNumber();
    			this.interAxisAngle[numStation]=           gd.getNextNumber();
    			this.horAxisErrPhi[numStation]=            gd.getNextNumber();
    			this.horAxisErrPsi[numStation]=            gd.getNextNumber();
    			this.entrancePupilForward[numStation]=     gd.getNextNumber();
    			this.centerAboveHorizontal[numStation]=    gd.getNextNumber();

    			this.GXYZ[numStation][0]=                  gd.getNextNumber();
    			this.GXYZ[numStation][1]=                  gd.getNextNumber();
    			this.GXYZ[numStation][2]=                  gd.getNextNumber();
    		}
    		this.sensorWidth=         (int) gd.getNextNumber();
    		this.sensorHeight=        (int) gd.getNextNumber();
    		this.shrinkGridForMask=   (int) gd.getNextNumber();
    		this.maskBlurSigma=             gd.getNextNumber();
    		this.decimateMasks=       (int) gd.getNextNumber();
	    	this.badNodeThreshold=          gd.getNextNumber();
	    	this.maxBadNeighb=        (int) gd.getNextNumber();
	    	this.minimalValidNodes=   (int) gd.getNextNumber();
	    	this.weightMultiImageMode=(int) gd.getNextNumber();
    		this.weightMultiExponent=       gd.getNextNumber();
    		this.weightDiameterExponent=    gd.getNextNumber();
    		this.weightYtoX=           0.01*gd.getNextNumber();
	    	this.minimalGridContrast=       gd.getNextNumber();
	    	this.shrinkBlurSigma =          gd.getNextNumber();
	    	this.shrinkBlurLevel =          gd.getNextNumber();
	    	this.balanceChannelWeightsMode= gd.getNextNumber();
	    	this.removeOverRMS =            gd.getNextNumber();
	    	this.removeOverRMSNonweighted=  gd.getNextNumber();
    		int numSubCams=           (int) gd.getNextNumber();
    		int numSubCam=            (int) gd.getNextNumber();
    		if (numSubCams!=this.eyesisSubCameras[0].length){
    			this.eyesisSubCameras=new EyesisSubCameraParameters[this.numStations][];
    	    	for (int numStation=0;numStation<numStations;numStation++) {
    	    		initSubCameras(numStation,numSubCams); // re-initialize from defaults, discard current!
    	    	}
    		}
    		if (numSubCam<0)numSubCam=0;
    		return gd.wasOKed()?numSubCam:-1;
    	}

    	public boolean showSubcameraDialog(
    			int numSubCam,
    			String title) {
    		GenericDialog gd = new GenericDialog(title);
    		for (int numStation=0;numStation<this.numStations;numStation++) {
    			EyesisSubCameraParameters subCam=this.eyesisSubCameras[numStation][numSubCam];
    			if (subCam!=null) {
		    		gd.addMessage("Channel weight "+subCam.channelWeightCurrent);
    				if (numStation==0){
        				gd.addNumericField("Channel "+numSubCam+" default weight",      subCam.channelWeightDefault, 5,8,"");
    				}
    				if (this.numStations>1) gd.addMessage("--- Station number "+numStation+" ---");
    				gd.addNumericField("Subcamera lens distortion model",               subCam.lensDistortionModel, 5,0,"");
    				gd.addCheckbox    ("Enable matching w/o laser pointers",            subCam.enableNoLaser);
    				gd.addNumericField("Subcamera azimuth",                             subCam.azimuth, 5,9,"degrees");
    				gd.addNumericField("Subcamera distance from the axis",              subCam.radius,  5,9,"mm");
    				gd.addNumericField("Subcamera height from the 'equator'",           subCam.height,  5,9,"mm");
    				gd.addNumericField("Optical axis heading (relative to azimuth)",    subCam.phi,     5,9,"degrees");
    				gd.addNumericField("Optical axis elevation (up from equator)",      subCam.theta,   5,9,"degrees");
    				gd.addNumericField("Camera roll, positive CW looking to the target",subCam.psi,     5,9,"degrees");
    				gd.addNumericField("Lens focal length",               subCam.focalLength, 5,8,"mm");
    				gd.addNumericField("Sensor pixel period",             subCam.pixelSize, 5,8,"um");
    				gd.addNumericField("Distortion radius (half width)",  subCam.distortionRadius, 6,8,"mm");
    				gd.addNumericField("Distortion A8 (r^8)",             subCam.distortionA8, 8,10,"");
    				gd.addNumericField("Distortion A7 (r^7)",             subCam.distortionA7, 8,10,"");
    				gd.addNumericField("Distortion A6 (r^6)",             subCam.distortionA6, 8,10,"");
    				gd.addNumericField("Distortion A5 (r^5)",             subCam.distortionA5, 8,10,"");
    				gd.addNumericField("Distortion A (r^4)",              subCam.distortionA, 8,10,"");
    				gd.addNumericField("Distortion B (r^3)",              subCam.distortionB, 8,10,"");
    				gd.addNumericField("Distortion C (r^2)",              subCam.distortionC, 8,10,"");
    				gd.addNumericField("Lens axis on the sensor (horizontal, from left edge)", subCam.px0, 4,9,"pixels");
    				gd.addNumericField("Lens axis on the sensor (vertical, from top  edge)",   subCam.py0, 4,9,"pixels");
    				
    				gd.addMessage("=== non-radial model parameters ===");
    				gd.addMessage("For r^2 (Distortion C):");
    				gd.addNumericField("Orthogonal elongation for r^2",   100*subCam.r_od[0][0], 8,10,"%");
    				gd.addNumericField("Diagonal elongation for r^2",     100*subCam.r_od[0][1], 8,10,"%");
    				gd.addMessage("For r^3 (Distortion B):");
    				gd.addNumericField("Distortion center shift X for r^3", 100*subCam.r_xy[0][0], 8,10,"%");
    				gd.addNumericField("Distortion center shift Y for r^3", 100*subCam.r_xy[0][1], 8,10,"%");
    				gd.addNumericField("Orthogonal elongation for r^3",  100*subCam.r_od[1][0], 8,10,"%");
    				gd.addNumericField("Diagonal elongation for r^3",    100*subCam.r_od[1][1], 8,10,"%");
    				gd.addMessage("For r^4 (Distortion A):");
    				gd.addNumericField("Distortion center shift X for r^4", 100*subCam.r_xy[1][0], 8,10,"%");
    				gd.addNumericField("Distortion center shift Y for r^4", 100*subCam.r_xy[1][1], 8,10,"%");
    				gd.addNumericField("Orthogonal elongation for r^4",  100*subCam.r_od[2][0], 8,10,"%");
    				gd.addNumericField("Diagonal elongation for r^4",    100*subCam.r_od[2][1], 8,10,"%");
    				gd.addMessage("For r^5 (Distortion A5):");
    				gd.addNumericField("Distortion center shift X for r^5", 100*subCam.r_xy[2][0], 8,10,"%");
    				gd.addNumericField("Distortion center shift Y for r^5", 100*subCam.r_xy[2][1], 8,10,"%");
    				gd.addNumericField("Orthogonal elongation for r^5",  100*subCam.r_od[3][0], 8,10,"%");
    				gd.addNumericField("Diagonal elongation for r^5",    100*subCam.r_od[3][1], 8,10,"%");
    				gd.addMessage("For r^6 (Distortion A6:");
    				gd.addNumericField("Distortion center shift X for r^6", 100*subCam.r_xy[3][0], 8,10,"%");
    				gd.addNumericField("Distortion center shift Y for r^6", 100*subCam.r_xy[3][1], 8,10,"%");
    				gd.addNumericField("Orthogonal elongation for r^6",  100*subCam.r_od[4][0], 8,10,"%");
    				gd.addNumericField("Diagonal elongation for r^6",    100*subCam.r_od[4][1], 8,10,"%");
    				gd.addMessage("For r^7 (Distortion A7):");
    				gd.addNumericField("Distortion center shift X for r^7", 100*subCam.r_xy[4][0], 8,10,"%");
    				gd.addNumericField("Distortion center shift Y for r^7", 100*subCam.r_xy[4][1], 8,10,"%");
    				gd.addNumericField("Orthogonal elongation for r^7",  100*subCam.r_od[5][0], 8,10,"%");
    				gd.addNumericField("Diagonal elongation for r^7",    100*subCam.r_od[5][1], 8,10,"%");
    				gd.addMessage("For r^8 (Distortion A8):");
    				gd.addNumericField("Distortion center shift X for r^8", 100*subCam.r_xy[5][0], 8,10,"%");
    				gd.addNumericField("Distortion center shift Y for r^8", 100*subCam.r_xy[5][1], 8,10,"%");
    				gd.addNumericField("Orthogonal elongation for r^8",  100*subCam.r_od[6][0], 8,10,"%");
    				gd.addNumericField("Diagonal elongation for r^8",    100*subCam.r_od[6][1], 8,10,"%");
    				
    			}
    		}
   	        gd.enableYesNoCancel("OK", "Done");
    	    WindowTools.addScrollBars(gd);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		double channelWeightDefault=1.0;
    		for (int numStation=0;numStation<this.numStations;numStation++) {
    			EyesisSubCameraParameters subCam=this.eyesisSubCameras[numStation][numSubCam];
    			if (subCam!=null) {
    				if (numStation==0){
    		    		gd.addMessage("Channel weight "+subCam.channelWeightCurrent);
    		    		channelWeightDefault=gd.getNextNumber();
    				}
    				subCam.channelWeightDefault= channelWeightDefault; // assign to all stations
    				subCam.lensDistortionModel= (int) gd.getNextNumber();
    				subCam.enableNoLaser =            gd.getNextBoolean();
    				subCam.azimuth=         gd.getNextNumber();
    				subCam.radius=          gd.getNextNumber();
    				subCam.height=          gd.getNextNumber();
    				subCam.phi=             gd.getNextNumber();
    				subCam.theta=           gd.getNextNumber();
    				subCam.psi=             gd.getNextNumber();
    				subCam.focalLength=     gd.getNextNumber();
    				subCam.pixelSize=       gd.getNextNumber();
    				subCam.distortionRadius=gd.getNextNumber();
    				subCam.distortionA8=    gd.getNextNumber();
    				subCam.distortionA7=    gd.getNextNumber();
    				subCam.distortionA6=    gd.getNextNumber();
    				subCam.distortionA5=    gd.getNextNumber();
    				subCam.distortionA=     gd.getNextNumber();
    				subCam.distortionB=     gd.getNextNumber();
    				subCam.distortionC=     gd.getNextNumber();
    				subCam.px0=             gd.getNextNumber();
    				subCam.py0=             gd.getNextNumber();
    				
    				subCam.r_od[0][0]= 0.01*gd.getNextNumber();
    				subCam.r_od[0][1]= 0.01*gd.getNextNumber();
    				subCam.r_xy[0][0]= 0.01*gd.getNextNumber();
    				subCam.r_xy[0][1]= 0.01*gd.getNextNumber();
    				subCam.r_od[1][0]= 0.01*gd.getNextNumber();
    				subCam.r_od[1][1]= 0.01*gd.getNextNumber();
    				subCam.r_xy[1][0]= 0.01*gd.getNextNumber();
    				subCam.r_xy[1][1]= 0.01*gd.getNextNumber();
    				subCam.r_od[2][0]= 0.01*gd.getNextNumber();
    				subCam.r_od[2][1]= 0.01*gd.getNextNumber();
    				subCam.r_xy[2][0]= 0.01*gd.getNextNumber();
    				subCam.r_xy[2][1]= 0.01*gd.getNextNumber();
    				subCam.r_od[3][0]= 0.01*gd.getNextNumber();
    				subCam.r_od[3][1]= 0.01*gd.getNextNumber();
    				subCam.r_xy[3][0]= 0.01*gd.getNextNumber();
    				subCam.r_xy[3][1]= 0.01*gd.getNextNumber();
    				subCam.r_od[4][0]= 0.01*gd.getNextNumber();
    				subCam.r_od[4][1]= 0.01*gd.getNextNumber();
    				subCam.r_xy[4][0]= 0.01*gd.getNextNumber();
    				subCam.r_xy[4][1]= 0.01*gd.getNextNumber();
    				subCam.r_od[5][0]= 0.01*gd.getNextNumber();
    				subCam.r_od[5][1]= 0.01*gd.getNextNumber();
    				subCam.r_xy[5][0]= 0.01*gd.getNextNumber();
    				subCam.r_xy[5][1]= 0.01*gd.getNextNumber();
    				subCam.r_od[6][0]= 0.01*gd.getNextNumber();
    				subCam.r_od[6][1]= 0.01*gd.getNextNumber();
    			}
    		}
    		return gd.wasOKed();
    	}

    	// TODO: make subcameras show numbers from 0
    	public boolean showDialog(String title) {
    		System.out.println("this.eyesisSubCameras.length="+this.eyesisSubCameras.length);
    		System.out.println("this.eyesisSubCameras[0].length="+this.eyesisSubCameras[0].length);
    		int subCam=1;
    		while (subCam<=this.eyesisSubCameras[0].length) {
        		System.out.println("subCam="+subCam);
    			subCam=subShowDialog(title, subCam);
    			if (subCam<0) return false;
    			if (subCam==0) return true;
    			while (subCam<=this.eyesisSubCameras[0].length) {
//    				if (!this.eyesisSubCameras[subCam-1].showDialog(title+": subcamera "+subCam)) break;
    				if (!showSubcameraDialog(subCam-1,title+": subcamera "+subCam+" ("+(subCam-1)+")")) break;
    				subCam++;
    			}
    		}
    		return true;
    	}
    	public int getLensDistortionModel(int stationNumber,int subCamNumber){
        	if (
        			(this.eyesisSubCameras==null) ||
        			(this.numStations<=stationNumber) ||
        			(this.eyesisSubCameras.length<=stationNumber) ||
        			(this.eyesisSubCameras[stationNumber].length<=subCamNumber)) throw new IllegalArgumentException
            ("Nonexistent subcamera "+subCamNumber+ " and/or station number="+stationNumber+" this.numStations="+this.numStations+" this.eyesisSubCameras.length="+this.eyesisSubCameras.length);
        	EyesisSubCameraParameters subCam=this.eyesisSubCameras[stationNumber][subCamNumber];
        	return subCam.lensDistortionModel;
    	}
    	boolean getEnableNoLaser(int stationNumber,int subCamNumber){
        	if (
        			(this.eyesisSubCameras==null) ||
        			(this.numStations<=stationNumber) ||
        			(this.eyesisSubCameras.length<=stationNumber) ||
        			(this.eyesisSubCameras[stationNumber].length<=subCamNumber)) throw new IllegalArgumentException
            ("Nonexistent subcamera "+subCamNumber+ " and/or station number="+stationNumber+" this.numStations="+this.numStations+" this.eyesisSubCameras.length="+this.eyesisSubCameras.length);
        	EyesisSubCameraParameters subCam=this.eyesisSubCameras[stationNumber][subCamNumber];
        	return subCam.enableNoLaser;
    	}
        /**
         * 
         * @param eyesisCameraParameters current parameters of the Eyesis camera, subcameras and goniometer
         * @param subCamNumber number of sub-camera (from 0)
         * @return array of the parameters (both individual sub-camera and common to all sub-cameras)
         *
         */
        public double [] getParametersVector(
        		int stationNumber,
        		int subCamNumber // 
                ){
        	if (
        			(this.eyesisSubCameras==null) ||
        			(this.numStations<=stationNumber) ||
        			(this.eyesisSubCameras.length<=stationNumber) ||
        			(this.eyesisSubCameras[stationNumber].length<=subCamNumber)) throw new IllegalArgumentException
            ("Nonexistent subcamera "+subCamNumber+ " and/or station number="+stationNumber+" this.numStations="+this.numStations+" this.eyesisSubCameras.length="+this.eyesisSubCameras.length);
        	
        	EyesisSubCameraParameters subCam=this.eyesisSubCameras[stationNumber][subCamNumber];
//        	System.out.println("getParametersVector("+stationNumber+","+subCamNumber+"), subCam is "+((subCam==null)?"null":"NOT null"));
        	double [] parVect={
        			subCam.azimuth,            // 0 azimuth of the lens entrance pupil center, degrees, clockwise looking from top
        			subCam.radius,             // 1 mm, distance from the rotation axis
        			subCam.height,             // 2 mm, up (was downwards?) - from the origin point
        			subCam.phi,                // 3 degrees, optical axis from azimuth/r vector, clockwise
        			subCam.theta,              // 4 degrees, optical axis from the eyesis horizon, positive - up
        			subCam.psi,                // 5 degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
        			this.goniometerHorizontal[stationNumber], // 6 goniometer rotation around "horizontal" axis (tilting from the target - positive)
        			this.goniometerAxial[stationNumber],      // 7 goniometer rotation around Eyesis axis (clockwise in plan - positive 
        			this.interAxisDistance[stationNumber],    // 8 distance in mm between two goniometer axes
        			this.interAxisAngle[stationNumber],       // 9 angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
        			this.horAxisErrPhi[stationNumber],        //10 angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)
        			this.horAxisErrPsi[stationNumber],        //11 angle in degrees "horizontal" goniometer axis is rotated around moving X axis (up)
            		this.entrancePupilForward[stationNumber], //12 common to all lenses - distance from the sensor to the lens entrance pupil
            		this.centerAboveHorizontal[stationNumber],//13 camera center distance along camera axis above the closest point to horizontal rotation axis (adds to height of each
        			
        			this.GXYZ[stationNumber][0],              //14 (12) coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system 
        			this.GXYZ[stationNumber][1],              //15 (13)  y
        			this.GXYZ[stationNumber][2],              //16 (14)  z
        			subCam.focalLength,        //17 (15)  lens focal length
        			subCam.px0,                //18 (16) center of the lens on the sensor, pixels
        			subCam.py0,                //19 (17) center of the lens on the sensor, pixels
        			subCam.distortionA8,       //20 (18)  r^8 (normalized to focal length or to sensor half width?)
        			subCam.distortionA7,       //21 (19)  r^7 (normalized to focal length or to sensor half width?)
        			subCam.distortionA6,       //22 (20)  r^6 (normalized to focal length or to sensor half width?)
        			subCam.distortionA5,       //23 (21)  r^5 (normalized to focal length or to sensor half width?)
        			subCam.distortionA,        //24 (22)  r^4 (normalized to focal length or to sensor half width?)
        			subCam.distortionB,        //25 (23)  r^3
        			subCam.distortionC,        //26 (24)  r^2
        			subCam.r_od[0][0],
        			subCam.r_od[0][1],
        			subCam.r_xy[0][0],
        			subCam.r_xy[0][1],
        			subCam.r_od[1][0],
        			subCam.r_od[1][1],
        			subCam.r_xy[1][0],
        			subCam.r_xy[1][1],
        			subCam.r_od[2][0],
        			subCam.r_od[2][1],
        			subCam.r_xy[2][0],
        			subCam.r_xy[2][1],
        			subCam.r_od[3][0],
        			subCam.r_od[3][1],
        			subCam.r_xy[3][0],
        			subCam.r_xy[3][1],
        			subCam.r_od[4][0],
        			subCam.r_od[4][1],
        			subCam.r_xy[4][0],
        			subCam.r_xy[4][1],
        			subCam.r_od[5][0],
        			subCam.r_od[5][1],
        			subCam.r_xy[5][0],
        			subCam.r_xy[5][1],
        			subCam.r_od[6][0],
        			subCam.r_od[6][1]
        	};
        	// Global parameters, not adjusted - just copied once when camera is selected
    // or should they stay fixed and not copied at all?    	
//    		this.lensDistortionParameters.pixelSize=subCam.pixelSize; // has to be set separately
//    		this.lensDistortionParameters.distortionRadius=subCam.distortionRadius;
        	return parVect;
        }
        
 //       public int getNumSubCameras (){return (this.eyesisSubCameras==null)?0:this.eyesisSubCameras.length;}
        public int getGoniometerHorizontalIndex(){return 6;}
        public int getGoniometerAxialIndex(){return 7;}
        public int getInterAxisAngleIndex(){return 9;}
        public int getSensorWidth() { return this.sensorWidth;}
        public int getSensorHeight() { return this.sensorHeight;}
        public int getSensorWidth(int subCam) { return this.sensorWidth;} // for the future? different sensors
        public int getSensorHeight(int subCam) { return this.sensorHeight;}// for the future? different sensors
        public double getPixelSize(int subCamNumber){return  this.eyesisSubCameras[0][subCamNumber].pixelSize;} // use station 0's pixel size
        public double getDistortionRadius(int subCamNumber){return  this.eyesisSubCameras[0][subCamNumber].distortionRadius;}

        public void setParametersVectorAllStations(
        		double [] parVect,
        		boolean [] update,
        		int subCamNumber // 
                ){
        	for (int stationNumber=0;stationNumber<this.numStations;stationNumber++){
        		setParametersVector(
                		parVect,
                		update,
                		stationNumber,
                		subCamNumber );
        	}
        	
        }

        /**
         * Set camera/subcamera parameters from the parameters vector
         * @param parVect array of the parameters (both individual sub-camera and common to all sub-cameras)
         * @param update - which parameter of the vector to update
         * @param eyesisCameraParameters current parameters of the Eyesis camera, subcameras and goniometer
         * @param subCamNumber number of sub-camera (from 0)
         */
        public void setParametersVector(
        		double [] parVect,
        		boolean [] update,
        		int stationNumber,
        		int subCamNumber // 
                ){
//        	if (parVect.length!=27) throw new IllegalArgumentException ("Wrong length of the parameters vector: "+parVect.length+"(should be 27)");
        	if (parVect.length!=53) throw new IllegalArgumentException ("Wrong length of the parameters vector: "+parVect.length+"(should be 53)");
        	if (
        			(this.eyesisSubCameras==null) ||
        			(this.numStations<=stationNumber) ||
        			(this.eyesisSubCameras.length<=stationNumber) ||
        			(this.eyesisSubCameras[stationNumber].length<=subCamNumber)) throw new IllegalArgumentException
            ("Nonexistent subcamera "+subCamNumber+ " and/or station number="+stationNumber);
        	EyesisSubCameraParameters subCam=this.eyesisSubCameras[stationNumber][subCamNumber];
        	if (update[0]) subCam.azimuth=parVect[0];            // 0 azimuth of the lens entrance pupil center, degrees, clockwise looking from top
        	if (update[1]) subCam.radius=parVect[1];             // 1 mm, distance from the rotation axis
        	if (update[2]) subCam.height=parVect[2];             // 2 mm, up (was downwards?) - from the origin point
        	if (update[2]) subCam.phi=parVect[3];                // 3 degrees, optical axis from azimuth/r vector, clockwise
        	if (update[4]) subCam.theta=parVect[4];              // 4 degrees, optical axis from the eyesis horizon, positive - up
        	if (update[5]) subCam.psi=parVect[5];                // 5 degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
        	if (update[6]) this.goniometerHorizontal[stationNumber]=parVect[6]; // 6 goniometer rotation around "horizontal" axis (tilting from the target - positive)
        	if (update[7]) this.goniometerAxial[stationNumber]=parVect[7];      // 7 goniometer rotation around Eyesis axis (clockwise in plan - positive 
        	if (update[8]) this.interAxisDistance[stationNumber]=parVect[8];    // 8 distance in mm between two goniometer axes
        	if (update[9]) this.interAxisAngle[stationNumber]=parVect[9];       // 9 angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
        	if (update[10]) this.horAxisErrPhi[stationNumber]=parVect[10];       //10 angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)
        	if (update[11]) this.horAxisErrPsi[stationNumber]=parVect[11];       //11 angle in degrees "horizontal" goniometer axis is rotated around moving X axis (up)
        	if (update[12]) this.entrancePupilForward[stationNumber]=parVect[12];             //12 coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system 
        	if (update[13]) this.centerAboveHorizontal[stationNumber]=parVect[13];             //13  y
        	if (update[14]) this.GXYZ[stationNumber][0]=parVect[14];             //14 coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system 
        	if (update[15]) this.GXYZ[stationNumber][1]=parVect[15];             //15  y
        	if (update[16]) this.GXYZ[stationNumber][2]=parVect[16];             //16  z
        	if (update[17]) subCam.focalLength=parVect[17];       //17  lens focal length
        	if (update[18]) subCam.px0=parVect[18];               //18 center of the lens on the sensor, pixels
        	if (update[19]) subCam.py0=parVect[19];               //19 center of the lens on the sensor, pixels
        	if (update[20]) subCam.distortionA8=parVect[20];      //20  r^8 (normalized to focal length or to sensor half width?)
        	if (update[21]) subCam.distortionA7=parVect[21];      //21  r^7 (normalized to focal length or to sensor half width?)
        	if (update[22]) subCam.distortionA6=parVect[22];      //22  r^6 (normalized to focal length or to sensor half width?)
        	if (update[23]) subCam.distortionA5=parVect[23];      //23  r^5 (normalized to focal length or to sensor half width?)
        	if (update[24]) subCam.distortionA= parVect[24];      //24  r^4 (normalized to focal length or to sensor half width?)
        	if (update[25]) subCam.distortionB= parVect[25];      //25  r^3
        	if (update[26]) subCam.distortionC= parVect[26];      //26  r^2
        	if (update[27]) subCam.r_od[0][0]= parVect[27];
        	if (update[28]) subCam.r_od[0][1]= parVect[28];
        	if (update[29]) subCam.r_xy[0][0]= parVect[29];
        	if (update[30]) subCam.r_xy[0][1]= parVect[30];
        	if (update[31]) subCam.r_od[1][0]= parVect[31];
        	if (update[32]) subCam.r_od[1][1]= parVect[32];
        	if (update[33]) subCam.r_xy[1][0]= parVect[33];
        	if (update[34]) subCam.r_xy[1][1]= parVect[34];
        	if (update[35]) subCam.r_od[2][0]= parVect[35];
        	if (update[36]) subCam.r_od[2][1]= parVect[36];
        	if (update[37]) subCam.r_xy[2][0]= parVect[37];
        	if (update[38]) subCam.r_xy[2][1]= parVect[38];
        	if (update[39]) subCam.r_od[3][0]= parVect[39];
        	if (update[40]) subCam.r_od[3][1]= parVect[40];
        	if (update[41]) subCam.r_xy[3][0]= parVect[41];
        	if (update[42]) subCam.r_xy[3][1]= parVect[42];
        	if (update[43]) subCam.r_od[4][0]= parVect[43];
        	if (update[44]) subCam.r_od[4][1]= parVect[44];
        	if (update[45]) subCam.r_xy[4][0]= parVect[45];
        	if (update[46]) subCam.r_xy[4][1]= parVect[46];
        	if (update[47]) subCam.r_od[5][0]= parVect[47];
        	if (update[48]) subCam.r_od[5][1]= parVect[48];
        	if (update[49]) subCam.r_xy[5][0]= parVect[49];
        	if (update[50]) subCam.r_xy[5][1]= parVect[50];
        	if (update[51]) subCam.r_od[6][0]= parVect[51];
        	if (update[52]) subCam.r_od[6][1]= parVect[52];
        }
        
    	public void initSubCameras(
    			int numStation,
    			int numSubCameras){
    		System.out.println("initSubCameras("+numStation+","+numSubCameras+")");
    		this.eyesisSubCameras[numStation]=new EyesisSubCameraParameters[numSubCameras];
    		for (int i=0;i<numSubCameras;i++)  this.eyesisSubCameras[numStation][i]=null;
    		if (numSubCameras==3) {
    			this.eyesisSubCameras[numStation][0]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					0.0, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					52.53, // double radius,  // mm, distance from the rotation axis
    					34.64, // double height,  // mm, up (was downwards) - from the origin point
    					0.0, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					180.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		4.5, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][1]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					30.0, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					60.0, // double radius,  // mm, distance from the rotation axis
    					-17.32, // double height,  // mm, up (was downwards) - from the origin point
    					-30.0, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					180.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		4.5, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][2]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					-30.0, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					60.0, // double radius,  // mm, distance from the rotation axis
    					-17.32, // double height,  // mm, up (was downwards) - from the origin point
    					30.0, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					180.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		4.5, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    		} else if (numSubCameras==1) {
    			this.eyesisSubCameras[numStation][0]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					0.0, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					0.0, // double radius,  // mm, distance from the rotation axis
    					0.0, // double height,  // mm, up (was downwards) - from the origin point
    					0.0, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					0.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		4.5, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    		} else if (numSubCameras == 21) {
    			// ================
    			// PHG21 parameters
    			//
    			this.eyesisSubCameras[numStation][0]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					0.0, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					46.57, // double radius,  // mm, distance from the rotation axis
    					0.0, // double height,  // mm, up (was downwards) - from the origin point
    					0.0, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					-90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][1]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					21.0, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					50.36, // double radius,  // mm, distance from the rotation axis
    					-15.0, // double height,  // mm, up (was downwards) - from the origin point
    					-53.0, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][2]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					-21.0, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					50.36, // double radius,  // mm, distance from the rotation axis
    					-15.0, // double height,  // mm, up (was downwards) - from the origin point
    					53.0, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][3]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					0.0, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					46.57, // double radius,  // mm, distance from the rotation axis
    					70.0, // double height,  // mm, up (was downwards) - from the origin point
    					0.0, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					-90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][4]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					21.0, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					50.36, // double radius,  // mm, distance from the rotation axis
    					55.0, // double height,  // mm, up (was downwards) - from the origin point
    					-53.0, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][5]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					-21.0, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					50.36, // double radius,  // mm, distance from the rotation axis
    					55.0, // double height,  // mm, up (was downwards) - from the origin point
    					53.0, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][6]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					52.47, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					76.45, // double radius,  // mm, distance from the rotation axis
    					35.0, // double height,  // mm, up (was downwards) - from the origin point
    					-52.47, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					-90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][7]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					59.13, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					91.65, // double radius,  // mm, distance from the rotation axis
    					20.0, // double height,  // mm, up (was downwards) - from the origin point
    					-91.13, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][8]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					42.16, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					63.43, // double radius,  // mm, distance from the rotation axis
    					20.0, // double height,  // mm, up (was downwards) - from the origin point
    					-10.16, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][9]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					52.47, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					76.45, // double radius,  // mm, distance from the rotation axis
    					-35.0, // double height,  // mm, up (was downwards) - from the origin point
    					-52.47, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					-90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][10]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					59.13, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					91.65, // double radius,  // mm, distance from the rotation axis
    					-50.0, // double height,  // mm, up (was downwards) - from the origin point
    					-91.13, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][11]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					42.16, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					63.43, // double radius,  // mm, distance from the rotation axis
    					-50.0, // double height,  // mm, up (was downwards) - from the origin point
    					-10.16, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][12]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					0.0, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					46.57, // double radius,  // mm, distance from the rotation axis
    					-70.0, // double height,  // mm, up (was downwards) - from the origin point
    					0.0, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					-90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][13]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					21.0, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					50.36, // double radius,  // mm, distance from the rotation axis
    					-85.0, // double height,  // mm, up (was downwards) - from the origin point
    					-53.0, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][14]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					-21.0, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					50.36, // double radius,  // mm, distance from the rotation axis
    					-85.0, // double height,  // mm, up (was downwards) - from the origin point
    					53.0, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][15]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					-52.47, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					76.45, // double radius,  // mm, distance from the rotation axis
    					-35.0, // double height,  // mm, up (was downwards) - from the origin point
    					52.47, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					-90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][16]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					-42.16, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					63.43, // double radius,  // mm, distance from the rotation axis
    					-50.0, // double height,  // mm, up (was downwards) - from the origin point
    					10.16, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][17]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					-59.13, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					91.65, // double radius,  // mm, distance from the rotation axis
    					-50.0, // double height,  // mm, up (was downwards) - from the origin point
    					91.13, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][18]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					-52.47, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					76.45, // double radius,  // mm, distance from the rotation axis
    					35.0, // double height,  // mm, up (was downwards) - from the origin point
    					52.47, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					-90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][19]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					-42.16, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					63.43, // double radius,  // mm, distance from the rotation axis
    					20.0, // double height,  // mm, up (was downwards) - from the origin point
    					10.16, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			this.eyesisSubCameras[numStation][20]=new EyesisSubCameraParameters( //TODO:  modify for lens adjustment defaults?
    					defaultLensDistortionModel,
    					true,
    					-59.13, // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					91.65, // double radius,  // mm, distance from the rotation axis
    					20.0, // double height,  // mm, up (was downwards) - from the origin point
    					91.13, // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0, //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		5.4, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			//
    			// end of PHG21 parameters
    			// =======================
    			} else {
    			// default setup for the 26 sub-cameras	    		
    			for (int i=0;i<8;i++) if (i<numSubCameras) 	this.eyesisSubCameras[numStation][i]=new EyesisSubCameraParameters( // top 8 cameras
    					defaultLensDistortionModel,
    					true,
    					45.0*i,      // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					41.540,      // double radius,  // mm, distance from the rotation axis
    					42.883,      // double height,  // mm, up (was downwards?) - from the origin point
    					0.0,         // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					60.0,        //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		4.5, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			for (int i=8;i<16;i++) if (i<numSubCameras) 	this.eyesisSubCameras[numStation][i]=new EyesisSubCameraParameters( // middle 8 cameras
    					defaultLensDistortionModel,
    					true,
    					45.0*(i-8),  // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					54.525,      // double radius,  // mm, distance from the rotation axis
    					0.0, // double height,  // mm, up (was downwards) - from the origin point
    					0.0,         // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0,         //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		4.5, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			for (int i=16;i<24;i++) if (i<numSubCameras) 	this.eyesisSubCameras[numStation][i]=new EyesisSubCameraParameters( // bottom eight cameras
    					defaultLensDistortionModel,
    					true,
    					45.0*(i-16), // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					41.540,      // double radius,  // mm, distance from the rotation axis
    					-42.883,     // double height,  // mm, up (was downwards?) - from the origin point
    					0.0,         // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					-60.0,       //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					-90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		4.5, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					1.0); //channelWeightDefault
    			if (24<numSubCameras) 	this.eyesisSubCameras[numStation][24]=new EyesisSubCameraParameters(
    					defaultLensDistortionModel,
    					false,
    					90,      // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					12.025,  // double radius,  // mm, distance from the rotation axis
    					-807.0,  // double height,  // mm, up - from the origin point
    					0.0,     // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0,     //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					90.0,  //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		4.5, // double focalLength
    		    		2.2, // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					8.0); //channelWeightDefault (was 4)
    			if (25<numSubCameras) 	this.eyesisSubCameras[numStation][25]=new EyesisSubCameraParameters(
    					defaultLensDistortionModel,
    					false,
    					270,     // double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    					12.025,  // double radius,  // mm, distance from the rotation axis
    					-841.0,  // double height,  // mm, up - from the origin point
    					0.0,     // double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    					0.0,     //double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    					-90.0,    //double psi;      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    		    		4.5,    // double focalLength
    		    		2.2,    // double pixelSize (um)
    		    		2.8512, //double distortionRadius mm - half width of the sensor
    		    		0.0, // double distortionA8 // r^8 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA7 // r^7 (normalized to focal length or to sensor half width?)
    		    		0.0, // double distortionA6 // r^6 (normalized to focal length or to sensor half width?)
    		    		0.0,    // double distortionA5 // r^5 (normalized to focal length or to sensor half width?)
    		    		0.0,    // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		    		0.0,    // double distortionB // r^3
    		    		0.0,    // double distortionC // r^2
    		    		1296.0, // double px0=1296.0;          // center of the lens on the sensor, pixels
    					968.0, // double py0=968.0;           // center of the lens on the sensor, pixels
    					null,  // eccentricity for b,a,a5,a6,a7,a8
    					null,  // elongation for c,b,a,a5,a6,a7,a8
    					8.0); //channelWeightDefault (was 4)
    		}
    	}
    	public void recenterVertically(boolean [] subcams, boolean [] stations){
    		boolean sameCenterAboveHorizontal=true;
    		double cah=Double.NaN;
    		for (int i=0;i<stations.length;i++){
    			if (stations[i]) {
    				if (Double.isNaN(cah)) cah=this.centerAboveHorizontal[i];
    				else if (cah!=this.centerAboveHorizontal[i]) sameCenterAboveHorizontal=false;
    			}
    		}
    		if (Double.isNaN(cah)){
    			System.out.println("No stations enabled, nothing to do for vertical centering");
    			return;
    		}
    		System.out.println("centerAboveHorizontal "+
    		(sameCenterAboveHorizontal?"is common for all stations":"differs between stations"));
    		
    		if (sameCenterAboveHorizontal){
    	   		System.out.println("Re-centering vertically for centerAboveHorizontal common for all stations ");
	
    			double sumWeightedHeights=0.0;
    			double sumWeights=0.0;
    					
    			for (int i=0;i<stations.length;i++) if (stations[i]) {
    				for (int subIndex=0; subIndex<subcams.length;subIndex++) if (subcams[subIndex]){
    					System.out.println("Averaging station "+i+", subcamera "+subIndex);
    					sumWeights+=stationWeight[i];
    					sumWeightedHeights+=stationWeight[i]*eyesisSubCameras[i][subIndex].height;
    				}
    			}
        		if (sumWeights==0.0){
        			System.out.println("No subcams are enabled, nothing to do for vertical centering");
        			return;
        		}
    			sumWeightedHeights/=sumWeights;
    			for (int i=0;i<stations.length;i++) if (stations[i]) {
    				for (int subIndex=0; subIndex<subcams.length;subIndex++) { // need to update all channels, not only averaged
    					eyesisSubCameras[i][subIndex].height-=sumWeightedHeights;
    				}
					this.centerAboveHorizontal[i]+=sumWeightedHeights;
    			}
    		} else {
    	   		System.out.println("Re-centering vertically for centerAboveHorizontal individual for each station");
    			for (int i=0;i<stations.length;i++) if (stations[i]) {
        			double sumWeightedHeights=0.0;
        			double sumWeights=0.0;
    				for (int subIndex=0; subIndex<subcams.length;subIndex++) if (subcams[subIndex]){
    					System.out.println("Averaging station "+i+", subcamera "+subIndex);
    					sumWeights+=stationWeight[i];
    					sumWeightedHeights+=stationWeight[i]*eyesisSubCameras[i][subIndex].height;
    				}
            		if (sumWeights==0.0){
            			System.out.println("No subcams are enabled, nothing to do for vertical centering");
            			return;
            		}
        			sumWeightedHeights/=sumWeights;
    				for (int subIndex=0; subIndex<subcams.length;subIndex++) {  // need to update all channels, not only averaged
    					eyesisSubCameras[i][subIndex].height-=sumWeightedHeights;
    				}
					this.centerAboveHorizontal[i]+=sumWeightedHeights;
    			}

    		}
    	}
    }
 
    