import java.util.Properties;

import ij.IJ;
import ij.gui.GenericDialog;
import Jama.Matrix;

//import EyesisCameraParameters;
//import Distortions.EyesisSubCameraParameters;
//import LensDistortionParameters;
//import PatternParameters;
//import DistortionCalibrationData.GridImageParameters;
//import DistortionCalibrationData.GridImageSet;

	public class LensDistortionParameters {
		/*
		 * Hugin Rsrc=a*Rdest^4+b*Rdest^3+c*Rdest^2+d*Rdest; d=1-(a+b+c)
		 */
		// lens parameters (add more later)

		final int numInputs=27; // with A8...// 24;   // parameters in subcamera+...
	    final int numOutputs=16; // with A8...//13;  // parameters in a single camera
		
		public double focalLength=4.5;
		public double pixelSize=  2.2; //um
		public double distortionRadius=  2.8512; // mm - half width of the sensor
		public double distortionA8=0.0; // r^8 (normalized to focal length or to sensor half width?)
		public double distortionA7=0.0; // r^7 (normalized to focal length or to sensor half width?)
		public double distortionA6=0.0; // r^6 (normalized to focal length or to sensor half width?)
		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
		public double distortionB=0.0; // r^3
		public double distortionC=0.0; // r^2
		// orientation/position parameters
		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
		public double x0=0;      // lens axis from pattern center, mm (to the right)
		public double y0=0;      // lens axis from pattern center, mm (down)
		public double z0=0;      // lens axis from pattern center, mm (away from the camera)
		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
		public double px0=1296.0;           // center of the lens on the sensor, pixels
		public double py0=968.0;           // center of the lens on the sensor, pixels
		public boolean flipVertical; // acquired image is mirrored vertically (mirror used)
		public int debugLevel=1; // was 2

// intermediate values
		public double phi, theta,psi,cPH,sPH,cTH,sTH,cPS,sPS;
		public double [][] rotMatrix=new double[3][3]; // includes mirroring for Y (target coordinates y- down, camera - y  up)
		
	    public LensDistortionParameters(
//	    		LensDistortionParameters lensDistortionParameters,
	    		boolean isTripod,
	            double [][] interParameterDerivatives, //partial derivative matrix from subcamera-camera-goniometer to single camera (12x21) if null - just values, no derivatives
	    		double [] parVect,
	    		boolean [] mask, // calculate only selected derivatives (all parVect values are still
	    		int debugLevel
//	    		boolean calculateDerivatives // calculate this.interParameterDerivatives -derivatives array (false - just this.values)
	    		){
	    	this.debugLevel=debugLevel;
	    	lensCalcInterParamers( // changed name to move calcInterParamers method from enclosing class
	    			this,
		    		isTripod,
		            interParameterDerivatives, //partial derivative matrix from subcamera-camera-goniometer to single camera (12x21) if null - just values, no derivatives
		    		parVect,
		    		mask // calculate only selected derivatives (all parVect values are still
//		    		boolean calculateDerivatives // calculate this.interParameterDerivatives -derivatives array (false - just this.values)
		    		);
	    }

		public LensDistortionParameters(
				double focalLength,
				double pixelSize,  //um
				double distortionRadius, // mm
				double distortionA8, // r^8
				double distortionA7, // r^7
				double distortionA6, // r^6
				double distortionA5, // r^5
				double distortionA, // r^4
				double distortionB, // r^3
				double distortionC, // r^2
				// orientation/position parameters
				double yaw,   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
				double pitch, // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
				double roll,  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
				double x0,      // lens axis from pattern center, mm (to the right)
				double y0,      // lens axis from pattern center, mm (down)
				double z0,      // lens axis from pattern center, mm (away from the camera)
				double distance, // distance from the lens input pupil to the pattern plane along the camera axis, mm
				double px0,           // center of the lens on the senosr, pixels
				double py0,           // center of the lens on the senosr, pixels
				boolean flipVertical // acquired image is mirrored vertically (mirror used)
		){
			setLensDistortionParameters(
			focalLength,
			pixelSize,
			distortionRadius,
			distortionA8,
			distortionA7,
			distortionA6,
			distortionA5,
			distortionA,
			distortionB,
			distortionC,
			yaw,
			pitch,
			roll,
			x0,
			y0,
			z0,
			distance,
			px0,
			py0,
			flipVertical);
		}
		public LensDistortionParameters(){
			this.focalLength=4.5;
			this.pixelSize=2.2;
			this.distortionRadius=2.8512;
			this.distortionA8=0.0;
			this.distortionA7=0.0;
			this.distortionA6=0.0;
			this.distortionA5=0.0;
			this.distortionA=0.0;
			this.distortionB=0.0;
			this.distortionC=0.0;
			this.yaw=0.0;
			this.pitch=0.0;
			this.roll=0.0;
			this.x0=0.0;
			this.y0=0.0;
			this.z0=0.0;
			this.distance=2360.0;
			this.px0=1296.0;
			this.py0=968.0;
			this.flipVertical=true;
		}
		public LensDistortionParameters clone() {
			return new LensDistortionParameters(
					this.focalLength,
					this.pixelSize,
					this.distortionRadius,
					this.distortionA8,
					this.distortionA7,
					this.distortionA6,
					this.distortionA5,
					this.distortionA,
					this.distortionB,
					this.distortionC,
					this.yaw,
					this.pitch,
					this.roll,
					this.x0,
					this.y0,
					this.z0,
					this.distance,
					this.px0,
					this.py0,
					this.flipVertical
			);
		}
		public void setLensDistortionParameters(
				double focalLength,
				double pixelSize,  //um
				double distortionRadius, // mm
				double distortionA8, // r^7
				double distortionA7, // r^6
				double distortionA6, // r^5
				double distortionA5, // r^4
				double distortionA, // r^4
				double distortionB, // r^3
				double distortionC, // r^2
				// orientation/position parameters
				double yaw,   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
				double pitch, // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
				double roll,  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
				double x0,      // lens axis from pattern center, mm (to the right)
				double y0,      // lens axis from pattern center, mm (down)
				double z0,      // lens axis from pattern center, mm (away)
				double distance, // distance from the lens input pupil to the pattern plane along the camera axis, mm
				double px0,           // center of the lens on the sensor, pixels
				double py0,           // center of the lens on the sensor, pixels
				boolean flipVertical // acquired image is mirrored vertically (mirror used)
		){
			this.focalLength=focalLength;
			this.pixelSize=pixelSize;
			this.distortionRadius=distortionRadius;
			this.distortionA8=distortionA8;
			this.distortionA7=distortionA7;
			this.distortionA6=distortionA6;
			this.distortionA5=distortionA5;
			this.distortionA=distortionA;
			this.distortionB=distortionB;
			this.distortionC=distortionC;
			this.yaw=yaw;
			this.pitch=pitch;
			this.roll=roll;
			this.x0=x0;
			this.y0=y0;
			this.z0=z0;
			this.distance=distance;
			this.px0=px0;
			this.py0=py0;
			this.flipVertical=flipVertical;
			recalcCommons();
		}
		public void setIntrincicFromSubcamera(EyesisSubCameraParameters pars){
			setLensDistortionParameters(
					pars.focalLength,
					pars.pixelSize,  //um
					pars.distortionRadius, // mm
					pars.distortionA8, // r^7
					pars.distortionA7, // r^6
					pars.distortionA6, // r^5
					pars.distortionA5, // r^4
					pars.distortionA,  // r^4
					pars.distortionB,  // r^3
					pars.distortionC,  // r^2
					// orientation/position parameters
					this.yaw,          // (keep) angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
					this.pitch,        // (keep) angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
					this.roll,         // (keep) angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
					this.x0,           // (keep) lens axis from pattern center, mm (to the right)
					this.y0,           // (keep) lens axis from pattern center, mm (down)
					this.z0,           // (keep) lens axis from pattern center, mm (away)
					this.distance,     // (keep) distance from the lens input pupil to the pattern plane along the camera axis, mm
					pars.px0,          //        center of the lens on the sensor, pixels
					pars.py0,          // center of the lens on the sensor, pixels
					this.flipVertical  // (keep)  acquired image is mirrored vertically (mirror used)
			);
		}
		
		
		public void setLensDistortionParameters(LensDistortionParameters ldp
		){
			this.focalLength=ldp.focalLength;
			this.pixelSize=ldp.pixelSize;
			this.distortionRadius=ldp.distortionRadius;
			this.distortionA8=ldp.distortionA8;
			this.distortionA7=ldp.distortionA7;
			this.distortionA6=ldp.distortionA6;
			this.distortionA5=ldp.distortionA5;
			this.distortionA=ldp.distortionA;
			this.distortionB=ldp.distortionB;
			this.distortionC=ldp.distortionC;
			this.yaw=ldp.yaw;
			this.pitch=ldp.pitch;
			this.roll=ldp.roll;
			this.x0=ldp.x0;
			this.y0=ldp.y0;
			this.z0=ldp.z0;
			this.distance=ldp.distance;
			this.px0=ldp.px0;
			this.py0=ldp.py0;
			this.flipVertical=ldp.flipVertical;
			recalcCommons();
		}
		// just for debugging				
		public void setLensDistortionParameters(LensDistortionParameters ldp,
				int index, // parameter to add delta, 1..13->14->17
				double delta
		){
			this.focalLength=ldp.focalLength+((index==7)?delta:0);
			this.pixelSize=ldp.pixelSize;
			this.distortionRadius=ldp.distortionRadius;
			this.distortionA8=ldp.distortionA8+((index==9)?delta:0);
			this.distortionA7=ldp.distortionA7+((index==10)?delta:0);
			this.distortionA6=ldp.distortionA6+((index==11)?delta:0);
			this.distortionA5=ldp.distortionA5+((index==12)?delta:0);
			this.distortionA=ldp.distortionA+((index==13)?delta:0);
			this.distortionB=ldp.distortionB+((index==14)?delta:0);
			this.distortionC=ldp.distortionC+((index==15)?delta:0);
			this.yaw=ldp.yaw+((index==1)?delta:0);
			this.pitch=ldp.pitch+((index==2)?delta:0);
			this.roll=ldp.roll+((index==3)?delta:0);
			this.x0=ldp.x0+((index==4)?delta:0);
			this.y0=ldp.y0+((index==5)?delta:0);
			this.z0=ldp.z0+((index==6)?delta:0);
			this.distance=ldp.distance+((index==8)?delta:0);
			this.px0=ldp.px0+((index==16)?delta:0);
			this.py0=ldp.py0+((index==17)?delta:0);
			this.flipVertical=ldp.flipVertical;
			recalcCommons();
		}
		
// recalculate common (point-invariant) intermediate values (cos, sin, rotation matrix) 		
        public void recalcCommons(){
//		public double phi, theta,psi,cPH,sPH,cTH,sTH,cPS,sPS;
//    		public double [][] rotMatrix=new double[3][3];
        	this.phi=  this.yaw*Math.PI/180;
        	this.theta=this.pitch*Math.PI/180;
        	this.psi=  this.roll*Math.PI/180;
        	this.sPH=Math.sin(this.phi);
        	this.cPH=Math.cos(this.phi);
        	this.sTH=Math.sin(this.theta);
        	this.cTH=Math.cos(this.theta);
        	this.sPS=Math.sin(this.psi);
        	this.cPS=Math.cos(this.psi);
/*
| Xe |   |   0  |   | cPS*cPH+sPS*sTH*sPH   -sPS*cTH  -cPS*sPH+sPS*sTH*cPH |   | Xp |
| Ye | = |   0  | + | sPS*cPH-cPS*sTH*sPH    cPS*cTH  -sPS*sPH-cPS*sTH*cPH | * |-Yp |
| Ze |   | dist |   | cTH*sPH                    sTH       cTH*cPH         |   | Zp |

| PX | =(1000*f)/(Ze*Psz) * |  Xe | + | PX0 |
| PY | =                    | -Ye |   | PY0 |


 Xe =   (cPS*cPH+sPS*sTH*sPH)*Xp   +sPS*cTH*Yp  +(-cPS*sPH+sPS*sTH*cPH)*Zp
 Ye =   (sPS*cPH-cPS*sTH*sPH)*Xp   -cPS*cTH*Yp  +(-sPS*sPH-cPS*sTH*cPH)*Zp
 Ze =   (cTH*sPH)*Xp               -sTH*Yp      +( cTH*cPH)* Zp + dist

        	
theta==0, psi==0:        	
 Xe =   (cPH)*Xp 
 Ye =   Yp
 Ze =   cPH* Zp + dist

(4) PXmmc =f/(cPH* Zp + dist)*  (cPH)*Xp  // mm, left from the lens axis intersection with the sensor
        	
dPXmmc/dphi=
        	
 */
        	this.rotMatrix[0][0]= cPS*cPH+sPS*sTH*sPH;
        	this.rotMatrix[0][1]= sPS*cTH;
           	this.rotMatrix[0][2]=-cPS*sPH+sPS*sTH*cPH;
           	this.rotMatrix[1][0]= sPS*cPH-cPS*sTH*sPH;
           	this.rotMatrix[1][1]=-cPS*cTH;
           	this.rotMatrix[1][2]=-sPS*sPH-cPS*sTH*cPH;
           	this.rotMatrix[2][0]= cTH*sPH;
           	this.rotMatrix[2][1]=-sTH;
           	this.rotMatrix[2][2]= cTH*cPH;
           	if (this.debugLevel>2){
           		System.out.println("recalcCommons():this.rotMatrix:");
        		(new Matrix(this.rotMatrix)).print(10, 5);
           	}
        	
        }
        
        private String [][] descriptions=
        {
        		{"distance",    "Distance from the intersection of the lens axis with z=0 target plane to the camera lens entrance pupil", "mm", "e"},
        		{"x0",          "Lens axis from pattern center, (to the right)", "mm","e"},
        		{"y0",          "Lens axis from pattern center, (down)", "mm","e"},
        		{"yaw",         "Angle from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top", "degrees","e"},
        		{"pitch",       "Angle from perpendicular to the pattern, 0 - towards wall, positive - up", "degrees","e"},
        		{"roll",        "Angle around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise", "degrees","e"},
        		{"focalLength", "Lens focal length", "mm","i"},
        		{"px0",         "horizontal (left to right) pixel number of the lens axis intersection with the sensor", "pix","i"},
        		{"py0",         "vertical (up to down) pixel number of the lens axis intersection with the sensor", "pix","i"},
        		{"distortionA8","Lens distortion coefficient for r^8 (r0=half sensor width)", "r0","i"},
        		{"distortionA7","Lens distortion coefficient for r^7 (r0=half sensor width)", "r0","i"},
        		{"distortionA6","Lens distortion coefficient for r^6 (r0=half sensor width)", "r0","i"},
        		{"distortionA5","Lens distortion coefficient for r^5 (r0=half sensor width)", "r0","i"},
        		{"distortionA", "Lens distortion coefficient for r^4 (r0=half sensor width)", "r0","i"},
        		{"distortionB", "Lens distortion coefficient for r^3 (r0=half sensor width)", "r0","i"},
        		{"distortionC", "Lens distortion coefficient for r^2 (r0=half sensor width)", "r0","i"}
        };
        private int numberIntExtrisic(String type){
        	int num=0;
        	for (int i=0;i<this.descriptions.length;i++) if (type.indexOf(descriptions[i][3].charAt(0))>=0) num++;
        	return num;
        }
        /**
         * Verifies that the camera is looking towards the target
         * @return true if looking tio the target, false - if away
         */
        public boolean isTargetVisible(boolean verbose){
        	if (verbose) System.out.println("isTargetVisible(): this.distance="+this.distance+", this.yaw="+this.yaw+", this.pitch="+this.pitch);
        	if (this.distance <=0.0) return false;
        	if (Math.cos(this.yaw*Math.PI/180)<0.0) return false;
        	if (Math.cos(this.pitch*Math.PI/180)<0.0) return false;
        	return true;
        }

        public double [] getExtrinsicVector(){
        	double [] extVector = {
        			this.distance,
        			this.x0,
        			this.y0,
        			this.yaw,
        			this.pitch,
        			this.roll
        			};
        	return extVector;
        }
        public double [] getIntrinsicVector(){
        	double [] extVector = {
        			this.focalLength,
        			this.px0,
        			this.py0,
        			this.distortionA8,
        			this.distortionA7,
        			this.distortionA6,
        			this.distortionA5,
        			this.distortionA,
        			this.distortionB,
        			this.distortionC
        	};
        	return extVector;
        }
        public double [] getAllVector(){
        	double [] allVector = new double[getExtrinsicVector().length+getIntrinsicVector().length];
        	int index=0;
        	double [] extVector=getExtrinsicVector();
        	double [] intVector=getIntrinsicVector();
        	for (int i=0;i<extVector.length;i++) allVector[index++]= extVector[i];
        	for (int i=0;i<intVector.length;i++) allVector[index++]= intVector[i];
        	return allVector;
        }
        public void setAllVector(double [] vector){
        	if (vector.length!=(getExtrinsicVector().length+getIntrinsicVector().length)){
        		String msg="Parameter vector should have exactly"+(getExtrinsicVector().length+getIntrinsicVector().length)+" elements";
        		IJ.showMessage("Error",msg); 
        		throw new IllegalArgumentException (msg);
        	}
        	this.distance=    vector[ 0];
        	this.x0=          vector[ 1];
        	this.y0=          vector[ 2];
        	this.yaw=         vector[ 3];
        	this.pitch=       vector[ 4];
        	this.roll=        vector[ 5];
        	this.focalLength= vector[ 6];
        	this.px0=         vector[ 7];
        	this.py0=         vector[ 8];
        	this.distortionA8=vector[ 9];
        	this.distortionA7=vector[10];
        	this.distortionA6=vector[11];
        	this.distortionA5=vector[12];
        	this.distortionA= vector[13];
        	this.distortionB= vector[14];
        	this.distortionC= vector[15];
        	recalcCommons();
        }
        
        public String [] getExtrinsicNames()       {return getDescriptionStrings("e", 0);}
        public String [] getExtrinsicDescriptions(){return getDescriptionStrings("e", 1);}
        public String [] getExtrinsicUnits()       {return getDescriptionStrings("e", 2);}
        public String [] getIntrinsicNames()       {return getDescriptionStrings("i", 0);}
        public String [] getIntrinsicDescriptions(){return getDescriptionStrings("i", 1);}
        public String [] getIntrinsicUnits()       {return getDescriptionStrings("i", 2);}
        public String [] getAllNames()             {return getDescriptionStrings("ei", 0);}
        public String [] getAllDescriptions()      {return getDescriptionStrings("ei", 1);}
        public String [] getAllUnits()             {return getDescriptionStrings("ei", 2);}
        public String [] getAllFlags()             {return getDescriptionStrings("ei", 3);}

        private String [] getDescriptionStrings(String type, int var){
        	String [] s=new String [numberIntExtrisic(type)];
        	int num=0;
        	for (int i=0;i<this.descriptions.length;i++) if (type.indexOf(descriptions[i][3].charAt(0))>=0) s[num++]=this.descriptions[i][var];
        	return s;
        }
        
        
/*
 * Calculate pixel value of projection from the pattern point [xp,yp,zp] using current distortion/position parameters
(1) Xe =   (cPS*cPH+sPS*sTH*sPH)*(Xp-X0)   +sPS*cTH*(Yp-Y0)  +(-cPS*sPH+sPS*sTH*cPH)*(Zp-Z0)
(2) Ye =   (sPS*cPH-cPS*sTH*sPH)*(Xp-X0)   -cPS*cTH*(Yp-Y0)  +(-sPS*sPH-cPS*sTH*cPH)*(Zp-Z0)
(3) Ze =   (cTH*sPH)*(Xp-X0)               -sTH*(Yp-Y0)      +( cTH*cPH)* (Zp-Z0) + dist

(4) PXmmc =f/Ze*  Xe  // mm, left from the lens axis intersection with the sensor
(5) PYmmc =f/Ze*  Ye  // mm, up from the lens axis intersection with the sensor

(6) r=sqrt(PXmmc^2+PYmmc^2) // distance from the image point to the lens axis intersection with the sensor (pinhole model)
(7) kD=(Da*(r/r0)^3+Db*(r/r0)^2+Dc*(r/r0)^1+(1-Da-Db-Dc)) correction to the actual distance from the image point to the lens axis due to distortion
(8) xDist = kD *  PXmmc // horisontal distance (mm) from the lens axis on the sensor to the image point, mm (positive - right)
(9) yDist = kD *  PYmmc // vertical distance (mm) from the lens axis on the sensor to the image point, mm (positive - up)


(10) PX = 1000/Psz*( xDist) + PX0 // horizontal pixel of the image (positive - right)
(11) PY = 1000/Psz*(-yDist) + PY0 // vertical pixel of the image (positive - down)
        
 */
        public double [] patternToPixels(
        		double xp, // target point horizontal, positive - right,  mm
        		double yp, // target point vertical,   positive - down,  mm
        		double zp // target point horizontal, positive - away from camera,  mm
        		){
        	return calcPartialDerivatives(xp,yp,zp,false)[0];
        }
/*
0		public double x0=0;      // lens axis from pattern center, mm (to the right)
1		public double y0=0;      // lens axis from pattern center, mm (down)
2		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
3		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
4		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
5		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
6		public double focalLength=4.5;
7		public double px0=1296.0;           // center of the lens on the sensor, pixels
8		public double py0=968.0;           // center of the lens on the sensor, pixels
 		public double distortionRadius=  2.8512; // mm - half width of the sensor
9		public double distortionA8=0.0; // r^8 (normalized to focal length or to sensor half width?)
10		public double distortionA7=0.0; // r^7 (normalized to focal length or to sensor half width?)
11		public double distortionA6=0.0; // r^6 (normalized to focal length or to sensor half width?)
12		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
13		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
14		public double distortionB=0.0; // r^3
15		public double distortionC=0.0; // r^2

 */
        /**
         *  extract needed ones, and reorder partial derivatives to match calcInterParameters
         *  @param srcDerivatives - values and 15 derivatives for px, py
         */
        public double [][] reorderPartialDerivatives (double [][] srcDerivatives){
        	int [] order={
        			 4, // 0		public double x0=0;      // lens axis from pattern center, mm (to the right)
        			 5, // 1		public double y0=0;      // lens axis from pattern center, mm (down)
        			 8, // 2		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
        			 1, // 3		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
        			 2, // 4		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
        			 3, // 5		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
        			 7, // 6		public double focalLength=4.5;
        			16, // 7		public double px0=1296.0;           // center of the lens on the sensor, pixels
        			17, // 8		public double py0=968.0;           // center of the lens on the sensor, pixels
        			    // 	public double distortionRadius=  2.8512; // mm - half width of the sensor
       			     9, // 9		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
    			    10, // 10		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
    			    11, // 11		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
        			12, // 12		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
        			13, // 13		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
        			14, // 14		public double distortionB=0.0; // r^3
        			15  // 15		public double distortionC=0.0; // r^2
        	};
        	double [][] result = new double [order.length][2];
        	for (int i=0; i<order.length;i++){
        		result[i][0]=srcDerivatives[order[i]][0];
        		result[i][1]=srcDerivatives[order[i]][1];
        	}
        	return result;
        	
        }
/**
 * Reorder to match the sequence of names - seems to be different :-(        
 * @param srcDerivatives  values and 15 derivatives for px, py
 * @return
 */
        public double [][] reorderPartialDerivativesAsNames (double [][] srcDerivatives){
        	int [] order={
        			 8, // 2		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
        			 4, // 0		public double x0=0;      // lens axis from pattern center, mm (to the right)
        			 5, // 1		public double y0=0;      // lens axis from pattern center, mm (down)
        			 1, // 3		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
        			 2, // 4		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
        			 3, // 5		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
        			 7, // 6		public double focalLength=4.5;
        			16, // 7		public double px0=1296.0;           // center of the lens on the sensor, pixels
        			17, // 8		public double py0=968.0;           // center of the lens on the sensor, pixels
        			 9, // 9		public double distortionA8=0.0; // r^5 (normalized to focal length or to sensor half width?)
        			10, // 10		public double distortionA7=0.0; // r^5 (normalized to focal length or to sensor half width?)
        			11, // 11		public double distortionA6=0.0; // r^5 (normalized to focal length or to sensor half width?)
        			12, // 12		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
        			13, // 13		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
        			14, // 14		public double distortionB=0.0; // r^3
        			15  // 15		public double distortionC=0.0; // r^2
        	};
        	double [][] result = new double [order.length][2];
        	for (int i=0; i<order.length;i++){
        		result[i][0]=srcDerivatives[order[i]][0];
        		result[i][1]=srcDerivatives[order[i]][1];
        	}
        	return result;
        	
        }
/**
 * Calculate lens center location in target coordinates
 * @return lens center coordinates
 */
        public double [] getLensCenterCoordinates(){
        	double [] p={
        			 this.x0-this.distance*this.rotMatrix[2][0],
        			-this.y0-this.distance*this.rotMatrix[2][1], // this.y0 - up?
        			 this.z0-this.distance*this.rotMatrix[2][2]};
/*        	
        	Matrix MR=new Matrix(this.rotMatrix);
        	Matrix MRT=MR.transpose();
        	Matrix E=MR.times(MRT);
        	System.out.println("===MR==");  	MR.print(8, 5);System.out.println("");
        	System.out.println("===MRT==");  	MRT.print(8, 5);System.out.println("");
        	System.out.println("===E===");  	E.print(8, 5);System.out.println("");
        	
        	System.out.println("x0="+IJ.d2s(this.x0,1)+" y0="+IJ.d2s(this.y0,1)+" z0="+IJ.d2s(this.z0,1)+" this.distance="+IJ.d2s(this.distance,1));
        	System.out.println("phi="+IJ.d2s(this.phi,4)+" theta="+IJ.d2s(this.theta,4)+" psi="+IJ.d2s(this.psi,4));
        	
        	System.out.println("|"+IJ.d2s(this.rotMatrix[0][0],4)+" "+IJ.d2s(this.rotMatrix[0][1],4)+" "+IJ.d2s(this.rotMatrix[0][2],4)+"|");
        	System.out.println("|"+IJ.d2s(this.rotMatrix[1][0],4)+" "+IJ.d2s(this.rotMatrix[1][1],4)+" "+IJ.d2s(this.rotMatrix[1][2],4)+"|");
        	System.out.println("|"+IJ.d2s(this.rotMatrix[2][0],4)+" "+IJ.d2s(this.rotMatrix[2][1],4)+" "+IJ.d2s(this.rotMatrix[2][2],4)+"|");
*/        	
        	return p;
        }

/*
 * Calculate pixel value and partial derivatives for different parameters
 * of projection from the pattern point [xp,yp,zp] using current distortion/position parameters
 * [0][0] - pixel x - NaN if looking away
 * [1][ 0] - pixel y- NaN if looking away
 * [*][ 1] - pixel x[or y] partial derivative for phi    (yaw)
 * [*][ 2] - pixel x[or y] partial derivative for theta  (pitch)
 * [*][ 3] - pixel x[or y] partial derivative for psi    (roll)
 * [*][ 4] - pixel x[or y] partial derivative for X0 (intersection of the lens axis with zp==Z0 plane of the target)
 * [*][ 5] - pixel x[or y] partial derivative for Y0 (intersection of the lens axis with zp==Z0 plane of the target)
 * [*][ 6] - pixel x[or y] partial derivative for Z0 (intersection of the lens axis with zp==Z0 plane of the target) - not used
 * [*][ 7] - pixel x[or y] partial derivative for f  (focal length)
 * [*][ 8] - pixel x[or y] partial derivative for dist (distance from [X0,Y0,Z0] to the lens entrance pupil
 * [*][ 9] - pixel x[or y] partial derivative for Da8 (distortion coefficient for r^8)
 * [*][10] - pixel x[or y] partial derivative for Da7 (distortion coefficient for r^7)
 * [*][11] - pixel x[or y] partial derivative for Da6 (distortion coefficient for r^6)
 * [*][12] - pixel x[or y] partial derivative for Da5 (distortion coefficient for r^5)
 * [*][13] - pixel x[or y] partial derivative for Da (distortion coefficient for r^4)
 * [*][14] - pixel x[or y] partial derivative for Db (distortion coefficient for r^3)
 * [*][15] - pixel x[or y] partial derivative for Dc (distortion coefficient for r^2)
 * [*][16] - pixel x[or y] partial derivative for PX0 (lens axis on the sensor, horizontal, right,in pixels)
 * [*][17] - pixel x[or y] partial derivative for PY0 (lens axis on the sensor, vertical, down, in pixels)
 */
        /*
         * TODO: minimaize calculations for individual {xp,yp,zp}
         */
        public double [][] calcPartialDerivatives(
        		double xp, // target point horizontal, positive - right,  mm
        		double yp, // target point vertical,   positive - down,  mm
        		double zp, // target point horizontal, positive - away from camera,  mm
        		boolean calculateAll){ // calculate derivatives, false - values only
        	double partDeriv[][] = new double [calculateAll?18:1][2];
//        	double [] XYZ= {xp-this.x0, yp-this.y0, zp-this.z0};
        	double [] XYZ= {xp-this.x0, yp+this.y0, zp-this.z0};
        	double [] XeYeZe={
        			this.rotMatrix[0][0]*XYZ[0] + this.rotMatrix[0][1]*XYZ[1] + this.rotMatrix[0][2]*XYZ[2],
        			this.rotMatrix[1][0]*XYZ[0] + this.rotMatrix[1][1]*XYZ[1] + this.rotMatrix[1][2]*XYZ[2],
        			this.rotMatrix[2][0]*XYZ[0] + this.rotMatrix[2][1]*XYZ[1] + this.rotMatrix[2][2]*XYZ[2]+this.distance
        	};
        	double [] PXYmmc={this.focalLength/XeYeZe[2]*XeYeZe[0],this.focalLength/XeYeZe[2]*XeYeZe[1]};
        	double r= Math.sqrt(PXYmmc[0]*PXYmmc[0]+PXYmmc[1]*PXYmmc[1]);
        	double rr=r/this.distortionRadius;
//        	double kD=((this.distortionA*rr+this.distortionB)*rr+this.distortionC)*rr + 1.0-this.distortionA-this.distortionB-this.distortionC;
        	double kD=((((((this.distortionA8*rr+this.distortionA7)*rr+ this.distortionA6)*rr+ this.distortionA5)*rr+this.distortionA)*rr+this.distortionB)*rr+this.distortionC)*rr +
        	1.0-this.distortionA8-this.distortionA7-this.distortionA6-this.distortionA5-this.distortionA-this.distortionB-this.distortionC;
        	double [] xyDist={kD*PXYmmc[0],kD*PXYmmc[1]};
        	partDeriv[0][0]=  1000.0/this.pixelSize*xyDist[0] + this.px0;
        	partDeriv[0][1]= -1000.0/this.pixelSize*xyDist[1] + this.py0;
        	if (!calculateAll) {
        		// TODO: Looking away from the target, trying only with no dervatives. Do the same for derivatives too?
            	if (XeYeZe[2]<0.0) {
            		partDeriv[0][0]=Double.NaN;
            		partDeriv[0][1]=Double.NaN;
            	}
        		return partDeriv;
        	}
        	
            double [][] dXeYeZe=new double[9][3]; //[14];
            dXeYeZe[0][0]=XeYeZe[0];
            dXeYeZe[0][1]=XeYeZe[1];
            dXeYeZe[0][2]=XeYeZe[2];
// /dphi            
            dXeYeZe[1][0]=(-cPS*sPH + sPS*sTH*cPH)*XYZ[0]                    +(-cPS*cPH-sPS*sTH*sPH)*XYZ[2];
            dXeYeZe[1][1]=(-sPS*sPH-cPS*sTH*cPH)  *XYZ[0]                    +(-sPS*cPH+cPS*sTH*sPH)*XYZ[2];
            dXeYeZe[1][2]=(cTH*cPH)*XYZ[0]                                     -cTH*sPH* XYZ[2];
// /dtheta
            dXeYeZe[2][0]=(sPS*cTH*sPH)*XYZ[0]              -sPS*sTH*XYZ[1]  +(sPS*cTH*cPH)*XYZ[2];
            dXeYeZe[2][1]=(-cPS*cTH*sPH)*XYZ[0]             +cPS*sTH*XYZ[1]  +(-cPS*cTH*cPH)*XYZ[2];
            dXeYeZe[2][2]=(-sTH*sPH)*XYZ[0]                 -cTH*XYZ[1]      -sTH*cPH* XYZ[2];
// /dpsi            
            dXeYeZe[3][0]=(-sPS*cPH+cPS*sTH*sPH)*XYZ[0]     +cPS*cTH*XYZ[1]  +(sPS*sPH+cPS*sTH*cPH)*XYZ[2];
            dXeYeZe[3][1]=(cPS*cPH+ sPS*sTH*sPH)*XYZ[0]     +sPS*cTH*XYZ[1]  +(-cPS*sPH+sPS*sTH*cPH)*XYZ[2];
            dXeYeZe[3][2]=0.0;
// /dX0
//          dXeYeZe[4][0]=-cPS*cPH+sPS*sTH*sPH; // bad?
            dXeYeZe[4][0]=-cPS*cPH-sPS*sTH*sPH; // bad?
            dXeYeZe[4][1]=-sPS*cPH+cPS*sTH*sPH;
            dXeYeZe[4][2]=-cTH*sPH;
// /dY0
//            dXeYeZe[5][0]=-sPS*cTH;
//            dXeYeZe[5][1]= cPS*cTH;
//            dXeYeZe[5][2]= sTH;
            dXeYeZe[5][0]=+sPS*cTH;
            dXeYeZe[5][1]=-cPS*cTH;
            dXeYeZe[5][2]=-sTH;
// /dZ0
//          dXeYeZe[6][0]= cPS*sPH+sPS*sTH*cPH; //bad?
            dXeYeZe[6][0]= cPS*sPH-sPS*sTH*cPH; //bad?
            dXeYeZe[6][1]= sPS*sPH+cPS*sTH*cPH;
            dXeYeZe[6][2]=-cTH*cPH;
// /df
            dXeYeZe[7][0]=0.0;
            dXeYeZe[7][1]=0.0;
            dXeYeZe[7][2]=0.0;
// /ddist
            dXeYeZe[8][0]=0.0;
            dXeYeZe[8][1]=0.0;
            dXeYeZe[8][2]=1.0;

            double [][] dPXYmmc=new double[9][2]; //[14];
            dPXYmmc[0][0]=PXYmmc[0];
            dPXYmmc[0][1]=PXYmmc[1];
//(4) PXmmc =f/Ze*  Xe  // mm, left from the lens axis intersection with the sensor
//dPXmmc/dphi   = f/Ze * dXe/dphi - f*Xe/Ze^2 * dZe/dphi
            dPXYmmc[1][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[1][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[1][2]);
//dPXmmc/dtheta = f/Ze * dXe/dtheta - f*Xe/Ze^2 * dZe/dtheta
            dPXYmmc[2][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[2][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[2][2]);
//dPXmmc/dpsi   = f/Ze * dXe/dpsi - f*Xe/Ze^2 * dZe/dpsi
            dPXYmmc[3][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[3][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[3][2]);
//dPXmmc/dX0    = f/Ze * dXe/dX0 - f*Xe/Ze^2 * dZe/dX0
            dPXYmmc[4][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[4][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[4][2]);
//dPXmmc/dY0    = f/Ze * dXe/dY0 - f*Xe/Ze^2 * dZe/dY0
            dPXYmmc[5][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[5][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[5][2]);
//dPXmmc/dZ0    = f/Ze * dXe/dZ0 - f*Xe/Ze^2 * dZe/dZ0
            dPXYmmc[6][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[6][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[6][2]); //bad?
//dPXmmc/df     = Xe/Ze
            dPXYmmc[7][0]=dXeYeZe[0][0]/dXeYeZe[0][2];
//dPXmmc/ddist =  - f*Xe/Ze^2
            dPXYmmc[8][0]=-this.focalLength*dXeYeZe[0][0]/(dXeYeZe[0][2]*dXeYeZe[0][2]);

//(5) PYmmc =f/Ze*  Ye  // mm, up from the lens axis intersection with the sensor
//dPYmmc/dphi   = f/Ze * dYe/dphi - f*Ye/Ze^2 * dZe/dphi
            dPXYmmc[1][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[1][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[1][2]);
//dPYmmc/dtheta = f/Ze * dYe/dtheta - f*Ye/Ze^2 * dZe/dtheta
            dPXYmmc[2][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[2][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[2][2]);
//dPYmmc/dpsi   = f/Ze * dYe/dpsi - f*Ye/Ze^2 * dZe/dpsi
            dPXYmmc[3][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[3][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[3][2]);
//dPYmmc/dX0    = f/Ze * dYe/dX0 - f*Ye/Ze^2 * dZe/dX0
            dPXYmmc[4][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[4][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[4][2]);
//dPYmmc/dY0    = f/Ze * dYe/dY0 - f*Ye/Ze^2 * dZe/dY0
            dPXYmmc[5][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[5][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[5][2]);
//dPYmmc/dZ0    = f/Ze * dYe/dZ0 - f*Ye/Ze^2 * dZe/dZ0
            dPXYmmc[6][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[6][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[6][2]); // good?
//dPYmmc/df     = Ye/Ze
            dPXYmmc[7][1]=dXeYeZe[0][1]/dXeYeZe[0][2];
//dPYmmc/ddist =  - f*Ye/Ze^2
            dPXYmmc[8][1]=-this.focalLength*dXeYeZe[0][1]/(dXeYeZe[0][2]*dXeYeZe[0][2]);
            double[] dr=new double [9];
//(6) r=sqrt(PXmmc^2+PYmmc^2) // distance from the image point to the lens axis intersection with the sensor (pinhole model)
            dr[0]=r;
            // Added 3 zeros
            if (dr[0]<0.00000001*this.distortionRadius) dr[0]=0.00000001*this.distortionRadius; // avoid 1/0.0
// dr/dphi   = (PXmmc*dPXmmc/dphi+PYmmc*dPYmmc/dphi)/r
            dr[1]=(PXYmmc[0]*dPXYmmc[1][0]+PXYmmc[1]*dPXYmmc[1][1])/dr[0];
// dr/dtheta = (PXmmc*dPXmmc/dtheta+PYmmc*dPYmmc/dtheta)/r
            dr[2]=(PXYmmc[0]*dPXYmmc[2][0]+PXYmmc[1]*dPXYmmc[2][1])/dr[0];
// dr/dpsi   = (PXmmc*dPXmmc/dpsi+PYmmc*dPYmmc/dpsi)/r
            dr[3]=(PXYmmc[0]*dPXYmmc[3][0]+PXYmmc[1]*dPXYmmc[3][1])/dr[0];
// dr/dX0    = (PXmmc*dPXmmc/dX0+PYmmc*dPYmmc/dX0)/r
            dr[4]=(PXYmmc[0]*dPXYmmc[4][0]+PXYmmc[1]*dPXYmmc[4][1])/dr[0];
// dr/dY0    = (PXmmc*dPXmmc/dY0+PYmmc*dPYmmc/dY0)/r
            dr[5]=(PXYmmc[0]*dPXYmmc[5][0]+PXYmmc[1]*dPXYmmc[5][1])/dr[0];
// dr/dZ0    = (PXmmc*dPXmmc/dZ0+PYmmc*dPYmmc/dZ0)/r
            dr[6]=(PXYmmc[0]*dPXYmmc[6][0]+PXYmmc[1]*dPXYmmc[6][1])/dr[0];
// dr/df     = (PXmmc*dPXmmc/df+PYmmc*dPYmmc/df)/r
            dr[7]=(PXYmmc[0]*dPXYmmc[7][0]+PXYmmc[1]*dPXYmmc[7][1])/dr[0];
// dr/ddist =  (PXmmc*dPXmmc/ddist+PYmmc*dPYmmc/ddist)/r
            dr[8]=(PXYmmc[0]*dPXYmmc[8][0]+PXYmmc[1]*dPXYmmc[8][1])/dr[0];
//            double[] dkD=new double [12];
            double[] dkD=new double [16];
// (7) kD=(Da*(r/r0)^3+Db*(r/r0)^2+Dc*(r/r0)^1+(1-Da-Db-Dc)) correction to the actual distance from the image point to the lens axis due to distortion
// (7) kD=(Da5*(r/r0)^4+(Da*(r/r0)^3+Db*(r/r0)^2+Dc*(r/r0)^1+(1-Da-Db-Dc-Da5)) correction to the actual distance from the image point to the lens axis due to distortion
            dkD[0]=kD;
// dkD/dr=  Da/r0^3 * 3*r^2 + Db/r0^2*2*r + Dc/r0
// dkD/dr=  Da5/r0^4 * 4*r^3 + Da/r0^3 * 3*r^2 + Db/r0^2*2*r + Dc/r0
            
// dkD/dr=  1/r0*(3*Da*(r/r0)^2 + 2*Db*(r/r0) + Dc)
// dkD/dr=  1/r0*(4*Da5*(r/r0)^3 + 3*Da*(r/r0)^2 + 2*Db*(r/r0) + Dc)            
            
            double dkDdr=1.0/this.distortionRadius*((((((8*this.distortionA8*rr+6*this.distortionA7)*rr+5*this.distortionA6)*rr+4*this.distortionA5)*rr+3*this.distortionA)*rr+2*this.distortionB)*rr + this.distortionC);
// dkD/dphi   = dkD/dr * dr/dphi
            dkD[1]=dkDdr*dr[1];
// dkD/dtheta = dkD/dr * dr/dtheta
            dkD[2]=dkDdr*dr[2];
// dkD/dpsi   = dkD/dr * dr/dpsi
            dkD[3]=dkDdr*dr[3];
// dkD/dX0    = dkD/dr * dr/dX0
            dkD[4]=dkDdr*dr[4];
// dkD/dY0    = dkD/dr * dr/dY0
            dkD[5]=dkDdr*dr[5];
// dkD/dZ0    = dkD/dr * dr/dZ0
            dkD[6]=dkDdr*dr[6];
// dkD/df     = dkD/dr * dr/df
            dkD[7]=dkDdr*dr[7];
// dkD/ddist  = dkD/dr * dr/ddist
            dkD[8]=dkDdr*dr[8];

// dkD/dDa8   = (r/r0)^7 - 1
            dkD[ 9]=rr*rr*rr*rr*rr*rr*rr-1.0;
// dkD/dDa7   = (r/r0)^6 - 1
            dkD[10]=rr*rr*rr*rr*rr*rr-1.0;
// dkD/dDa6   = (r/r0)^5 - 1
            dkD[11]=rr*rr*rr*rr*rr-1.0;
            
            
// dkD/dDa5   = (r/r0)^4 - 1
            dkD[12]=rr*rr*rr*rr-1.0;
// dkD/dDa    = (r/r0)^3 - 1
            dkD[13]=rr*rr*rr-1.0;
// dkD/dDb    = (r/r0)^2 - 1
            dkD[14]=rr*rr-1.0;
// dkD/dDc    = (r/r0)   - 1
            dkD[15]=rr-1.0;
            
//            double[][] dxyDist=new double [12][2];
            double[][] dxyDist=new double [16][2];
// (8) xDist = kD *  PXmmc // horisontal distance (mm) from the lens axis on the sensor to the image point, mm (positive - right)
// (9) yDist = kD *  PYmmc // vertical distance (mm) from the lens axis on the sensor to the image point, mm (positive - up)
            dxyDist[0][0]=xyDist[0];
            dxyDist[0][1]=xyDist[1];
// dxDist/dphi   = dkD/dphi*PXmmc   + kD*dPXmmc/dphi
// dyDist/dphi   = dkD/dphi*PYmmc   + kD*dPYmmc/dphi
            dxyDist[ 1][0]=dkD[ 1]*dPXYmmc[0][0]+dkD[0]*dPXYmmc[ 1][0];
            dxyDist[ 1][1]=dkD[ 1]*dPXYmmc[0][1]+dkD[0]*dPXYmmc[ 1][1];
// dxDist/dtheta = dkD/dtheta*PXmmc + kD*dPXmmc/dtheta
            dxyDist[ 2][0]=dkD[ 2]*dPXYmmc[0][0]+dkD[0]*dPXYmmc[ 2][0];
            dxyDist[ 2][1]=dkD[ 2]*dPXYmmc[0][1]+dkD[0]*dPXYmmc[ 2][1];
// dxDist/dpsi   = dkD/dpsi*PXmmc   + kD*dPXmmc/dpsi
            dxyDist[ 3][0]=dkD[ 3]*dPXYmmc[0][0]+dkD[0]*dPXYmmc[ 3][0];
            dxyDist[ 3][1]=dkD[ 3]*dPXYmmc[0][1]+dkD[0]*dPXYmmc[ 3][1];
// dxDist/dX0    = dkD/dX0*PXmmc    + kD*dPXmmc/dX0
            dxyDist[ 4][0]=dkD[ 4]*dPXYmmc[0][0]+dkD[0]*dPXYmmc[ 4][0]; // large error
            dxyDist[ 4][1]=dkD[ 4]*dPXYmmc[0][1]+dkD[0]*dPXYmmc[ 4][1];
// dxDist/dY0    = dkD/dY0*PXmmc    + kD*dPXmmc/dY0
            dxyDist[ 5][0]=dkD[ 5]*dPXYmmc[0][0]+dkD[0]*dPXYmmc[ 5][0];
            dxyDist[ 5][1]=dkD[ 5]*dPXYmmc[0][1]+dkD[0]*dPXYmmc[ 5][1];
// dxDist/dZ0    = dkD/dZ0*PXmmc    + kD*dPXmmc/dZ0
            dxyDist[ 6][0]=dkD[ 6]*dPXYmmc[0][0]+dkD[0]*dPXYmmc[ 6][0]; // large error
            dxyDist[ 6][1]=dkD[ 6]*dPXYmmc[0][1]+dkD[0]*dPXYmmc[ 6][1];
// dxDist/df     = dkD/df*PXmmc     + kD*dPXmmc/df
            dxyDist[ 7][0]=dkD[ 7]*dPXYmmc[0][0]+dkD[0]*dPXYmmc[ 7][0];
            dxyDist[ 7][1]=dkD[ 7]*dPXYmmc[0][1]+dkD[0]*dPXYmmc[ 7][1];
// dxDist/ddist  = dkD/ddist*PXmmc  + kD*dPXmmc/ddist
            dxyDist[ 8][0]=dkD[ 8]*dPXYmmc[0][0]+dkD[0]*dPXYmmc[ 8][0];
            dxyDist[ 8][1]=dkD[ 8]*dPXYmmc[0][1]+dkD[0]*dPXYmmc[ 8][1];

// dxDist/dDa5    = dkD/dDa8*PXmmc
            dxyDist[ 9][0]=dkD[ 9]*dPXYmmc[0][0];
            dxyDist[ 9][1]=dkD[ 9]*dPXYmmc[0][1];
// dxDist/dDa5    = dkD/dDa7*PXmmc
            dxyDist[10][0]=dkD[10]*dPXYmmc[0][0];
            dxyDist[10][1]=dkD[10]*dPXYmmc[0][1];
// dxDist/dDa5    = dkD/dD6a*PXmmc
            dxyDist[11][0]=dkD[11]*dPXYmmc[0][0];
            dxyDist[11][1]=dkD[11]*dPXYmmc[0][1];
            
            
// dxDist/dDa5    = dkD/dDa5*PXmmc`
            dxyDist[12][0]=dkD[12]*dPXYmmc[0][0];
            dxyDist[12][1]=dkD[12]*dPXYmmc[0][1];
// dxDist/dDa    = dkD/dDa*PXmmc
            dxyDist[13][0]=dkD[13]*dPXYmmc[0][0];
            dxyDist[13][1]=dkD[13]*dPXYmmc[0][1];
// dxDist/dDb    = dkD/dDb*PXmmc
            dxyDist[14][0]=dkD[14]*dPXYmmc[0][0];
            dxyDist[14][1]=dkD[14]*dPXYmmc[0][1];
// dxDist/dDc    = dkD/dDc*PXmmc
            dxyDist[15][0]=dkD[15]*dPXYmmc[0][0];
            dxyDist[15][1]=dkD[15]*dPXYmmc[0][1];

            double K=Math.PI/180; // multiply all derivatives my angles 
// (10) PX = 1000/Psz*( xDist) + PX0 // horizontal pixel of the image (positive - right)
// dPX/dphi   =  1000/Psz* dxDist/dphi
        	partDeriv[ 1][0]=  K*1000.0/this.pixelSize*dxyDist[ 1][0];
        	partDeriv[ 1][1]= -K*1000.0/this.pixelSize*dxyDist[ 1][1];
// dPX/dtheta =  1000/Psz* dxDist/dtheta
        	partDeriv[ 2][0]=  K*1000.0/this.pixelSize*dxyDist[ 2][0];
        	partDeriv[ 2][1]= -K*1000.0/this.pixelSize*dxyDist[ 2][1];
// dPX/dpsi   =  1000/Psz* dxDist/dpsi
        	partDeriv[ 3][0]=  K*1000.0/this.pixelSize*dxyDist[ 3][0];
        	partDeriv[ 3][1]= -K*1000.0/this.pixelSize*dxyDist[ 3][1];
// dPX/dX0    =  1000/Psz* dxDist/dX0
        	partDeriv[ 4][0]=  1000.0/this.pixelSize*dxyDist[ 4][0]; // large error
        	partDeriv[ 4][1]= -1000.0/this.pixelSize*dxyDist[ 4][1];
// dPX/dY0    =  1000/Psz* dxDist/dY0
        	partDeriv[ 5][0]=  1000.0/this.pixelSize*dxyDist[ 5][0];
        	partDeriv[ 5][1]= -1000.0/this.pixelSize*dxyDist[ 5][1];
// dPX/dZ0    =  1000/Psz* dxDist/dZ0
        	partDeriv[ 6][0]=  1000.0/this.pixelSize*dxyDist[ 6][0]; // large error, including DC (sig -0.8..+0.5, err -0.2..-0.1
        	partDeriv[ 6][1]= -1000.0/this.pixelSize*dxyDist[ 6][1];
// dPX/df     =  1000/Psz* dxDist/df
        	partDeriv[ 7][0]=  1000.0/this.pixelSize*dxyDist[ 7][0];
        	partDeriv[ 7][1]= -1000.0/this.pixelSize*dxyDist[ 7][1];
// dPX/ddist  =  1000/Psz* dxDist/ddist
        	partDeriv[ 8][0]=  1000.0/this.pixelSize*dxyDist[ 8][0];
        	partDeriv[ 8][1]= -1000.0/this.pixelSize*dxyDist[ 8][1];

// dPX/dDa8    =  1000/Psz* dxDist/dDa5
        	partDeriv[ 9][0]=  1000.0/this.pixelSize*dxyDist[ 9][0];
        	partDeriv[ 9][1]= -1000.0/this.pixelSize*dxyDist[ 9][1];
// dPX/dDa7    =  1000/Psz* dxDist/dDa5
        	partDeriv[10][0]=  1000.0/this.pixelSize*dxyDist[10][0];
        	partDeriv[10][1]= -1000.0/this.pixelSize*dxyDist[10][1];
// dPX/dDa6    =  1000/Psz* dxDist/dDa5
        	partDeriv[11][0]=  1000.0/this.pixelSize*dxyDist[11][0];
        	partDeriv[11][1]= -1000.0/this.pixelSize*dxyDist[11][1];
        	
// dPX/dDa5    =  1000/Psz* dxDist/dDa5
        	partDeriv[12][0]=  1000.0/this.pixelSize*dxyDist[12][0];
        	partDeriv[12][1]= -1000.0/this.pixelSize*dxyDist[12][1];
// dPX/dDa    =  1000/Psz* dxDist/dDa
        	partDeriv[13][0]=  1000.0/this.pixelSize*dxyDist[13][0];
        	partDeriv[13][1]= -1000.0/this.pixelSize*dxyDist[13][1];
// dPX/dDb    =  1000/Psz* dxDist/dDb
        	partDeriv[14][0]=  1000.0/this.pixelSize*dxyDist[14][0];
        	partDeriv[14][1]= -1000.0/this.pixelSize*dxyDist[14][1];
// dPX/dDc    =  1000/Psz* dxDist/dDc
        	partDeriv[15][0]=  1000.0/this.pixelSize*dxyDist[15][0];
        	partDeriv[15][1]= -1000.0/this.pixelSize*dxyDist[15][1];
// dPX/dPX0   =  1
// dPY/dPX0   =  0
        	partDeriv[16][0]=  1.0;
        	partDeriv[16][1]=  0.0;
// dPX/dPY0   =  0
// dPY/dPY0   =  1
        	partDeriv[17][0]=  0.0;
        	partDeriv[17][1]=  1.0;
            return partDeriv;
        }
        
        
        public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"focalLength",this.focalLength+"");
			properties.setProperty(prefix+"pixelSize",this.pixelSize+"");
			properties.setProperty(prefix+"distortionRadius",this.distortionRadius+"");
			properties.setProperty(prefix+"distortionA8",this.distortionA8+"");
			properties.setProperty(prefix+"distortionA7",this.distortionA7+"");
			properties.setProperty(prefix+"distortionA6",this.distortionA6+"");
			properties.setProperty(prefix+"distortionA5",this.distortionA5+"");
			properties.setProperty(prefix+"distortionA",this.distortionA+"");
			properties.setProperty(prefix+"distortionB",this.distortionB+"");
			properties.setProperty(prefix+"distortionC",this.distortionC+"");
			properties.setProperty(prefix+"yaw",this.yaw+"");
			properties.setProperty(prefix+"pitch",this.pitch+"");
			properties.setProperty(prefix+"roll",this.roll+"");
			properties.setProperty(prefix+"x0",this.x0+"");
			properties.setProperty(prefix+"y0",this.y0+"");
			properties.setProperty(prefix+"z0",this.z0+"");
			properties.setProperty(prefix+"distance",this.distance+"");
			properties.setProperty(prefix+"px0",this.px0+"");
			properties.setProperty(prefix+"py0",this.py0+"");
			properties.setProperty(prefix+"flipVertical",this.flipVertical+"");
		}
		public void getProperties(String prefix,Properties properties){
			if (properties.getProperty(prefix+"focalLength")!=null)
				this.focalLength=Double.parseDouble(properties.getProperty(prefix+"focalLength"));
			if (properties.getProperty(prefix+"pixelSize")!=null)
				this.pixelSize=Double.parseDouble(properties.getProperty(prefix+"pixelSize"));
			if (properties.getProperty(prefix+"distortionRadius")!=null)
				this.distortionRadius=Double.parseDouble(properties.getProperty(prefix+"distortionRadius"));
			if (properties.getProperty(prefix+"distortionA8")!=null)
				this.distortionA8=Double.parseDouble(properties.getProperty(prefix+"distortionA8"));
			if (properties.getProperty(prefix+"distortionA7")!=null)
				this.distortionA7=Double.parseDouble(properties.getProperty(prefix+"distortionA7"));
			if (properties.getProperty(prefix+"distortionA6")!=null)
				this.distortionA6=Double.parseDouble(properties.getProperty(prefix+"distortionA6"));
			if (properties.getProperty(prefix+"distortionA5")!=null)
				this.distortionA5=Double.parseDouble(properties.getProperty(prefix+"distortionA5"));
			if (properties.getProperty(prefix+"distortionA")!=null)
				this.distortionA=Double.parseDouble(properties.getProperty(prefix+"distortionA"));
			if (properties.getProperty(prefix+"distortionB")!=null)
				this.distortionB=Double.parseDouble(properties.getProperty(prefix+"distortionB"));
			if (properties.getProperty(prefix+"distortionC")!=null)
				this.distortionC=Double.parseDouble(properties.getProperty(prefix+"distortionC"));
			if (properties.getProperty(prefix+"yaw")!=null)
				this.yaw=Double.parseDouble(properties.getProperty(prefix+"yaw"));
			if (properties.getProperty(prefix+"pitch")!=null)
				this.pitch=Double.parseDouble(properties.getProperty(prefix+"pitch"));
			if (properties.getProperty(prefix+"roll")!=null)
				this.roll=Double.parseDouble(properties.getProperty(prefix+"roll"));
			if (properties.getProperty(prefix+"x0")!=null)
				this.x0=Double.parseDouble(properties.getProperty(prefix+"x0"));
			if (properties.getProperty(prefix+"y0")!=null)
				this.y0=Double.parseDouble(properties.getProperty(prefix+"y0"));
			if (properties.getProperty(prefix+"z0")!=null)
				this.z0=Double.parseDouble(properties.getProperty(prefix+"z0"));
			if (properties.getProperty(prefix+"distance")!=null)
				this.distance=Double.parseDouble(properties.getProperty(prefix+"distance"));
			if (properties.getProperty(prefix+"px0")!=null)
				this.px0=Double.parseDouble(properties.getProperty(prefix+"px0"));
			if (properties.getProperty(prefix+"py0")!=null)
				this.py0=Double.parseDouble(properties.getProperty(prefix+"py0"));
			if (properties.getProperty(prefix+"flipVertical")!=null)
				this.flipVertical=Boolean.parseBoolean(properties.getProperty(prefix+"flipVertical"));
			recalcCommons();
		}
		public boolean showDialog() {
			GenericDialog gd = new GenericDialog("Lens distortion, location and orientation");
			gd.addNumericField("Lens focal length",               this.focalLength, 3,6,"mm");
			gd.addNumericField("Sensor pixel period",             this.pixelSize, 3,6,"um");
			gd.addNumericField("Distortion radius (halw width)",  this.distortionRadius, 5,8,"mm");
			gd.addNumericField("Distortion A8(r^5)",              this.distortionA8, 6,8,"");
			gd.addNumericField("Distortion A7(r^5)",              this.distortionA7, 6,8,"");
			gd.addNumericField("Distortion A6(r^5)",              this.distortionA6, 6,8,"");
			gd.addNumericField("Distortion A5(r^5)",              this.distortionA5, 6,8,"");
			gd.addNumericField("Distortion A (r^4)",              this.distortionA, 6,8,"");
			gd.addNumericField("Distortion B (r^3)",              this.distortionB, 6,8,"");
			gd.addNumericField("Distortion C (r^2)",              this.distortionC, 6,8,"");
			gd.addNumericField("Lens axis from perpendicular to the pattern, positive - clockwise (from top)", this.yaw, 2,6,"degrees");
			gd.addNumericField("Lens axis from perpendicular to the pattern, positive - up",                   this.pitch, 2,6,"degrees");
			gd.addNumericField("Rotation around lens axis, positive - clockwise (looking to pattern)",         this.roll, 2,6,"degrees");
			gd.addNumericField("Lens axis from the pattern center, (to the right)",                            this.x0, 1,6,"mm");
			gd.addNumericField("Lens axis from the pattern center, (down)",                                    this.y0, 1,6,"mm");
			gd.addNumericField("Lens axis from the pattern center, (away from camera, normally 0.0)",          this.z0, 1,6,"mm");
			gd.addNumericField("Distance from the lens input pupil to the pattern plane along the camera axis",this.distance, 1,6,"mm");
			gd.addNumericField("Lens axis on the sensor (horizontal, from left edge)",                         this.px0, 1,6,"pixels");
			gd.addNumericField("Lens axis on the sensor (vertical, from top edge)",                            this.py0, 1,6,"pixels");
			gd.addCheckbox    ("Camera looks through the mirror",                                              this.flipVertical);
		
			gd.showDialog();
			if (gd.wasCanceled()) return false;
			this.focalLength=     gd.getNextNumber();
			this.pixelSize=       gd.getNextNumber();
			this.distortionRadius=gd.getNextNumber();
			this.distortionA8=    gd.getNextNumber();
			this.distortionA7=    gd.getNextNumber();
			this.distortionA6=    gd.getNextNumber();
			this.distortionA5=    gd.getNextNumber();
			this.distortionA=     gd.getNextNumber();
			this.distortionB=     gd.getNextNumber();
			this.distortionC=     gd.getNextNumber();
			this.yaw=             gd.getNextNumber();
			this.pitch=           gd.getNextNumber();
			this.roll=            gd.getNextNumber();
			this.x0=              gd.getNextNumber();
			this.y0=              gd.getNextNumber();
			this.z0=              gd.getNextNumber();
			this.distance=        gd.getNextNumber(); 
			this.px0=             gd.getNextNumber();
			this.py0=             gd.getNextNumber();
			this.flipVertical=   gd.getNextBoolean();
			return true;
		}
	    /**
	     * Calculate/set  this.lensDistortionParameters and this.interParameterDerivatives
     * UPDATE - Modifies lensDistortionParameters, not "this" formulti-threaded 
	     * @param parVect 21-element vector for eyesis sub-camera, including common and individual parameters
	     * @param mask -mask - which partial derivatives are needed to be calculated (others will be null)
	     * @param calculateDerivatives calculate array of partial derivatives, if false - just the values
	     */
	    public void lensCalcInterParamers(
	    		LensDistortionParameters lensDistortionParameters,
	    		boolean isTripod,
	            double [][] interParameterDerivatives, //partial derivative matrix from subcamera-camera-goniometer to single camera (12x21) if null - just values, no derivatives
	    		double [] parVect,
	    		boolean [] mask // calculate only selected derivatives (all parVect values are still
//	    		boolean calculateDerivatives // calculate this.interParameterDerivatives -derivatives array (false - just this.values)
	    		){
//	    	LensDistortionParameters lensDistortionParameters=this;
	    	boolean calculateDerivatives=(interParameterDerivatives!=null);  // calculate this.interParameterDerivatives -derivatives array (false - just this.values) 
	    	 // change meaning of goniometerHorizontal (tripod vertical) and goniometerAxial (tripod horizontal) 
//	    	boolean isTripod=this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.isTripod;
	    	double azimuth=parVect[0];
	    	double radius= parVect[1];
	    	double height= parVect[2];
	    	double phi=    parVect[3];
	    	double theta=  parVect[4];
	    	double psi=    parVect[5];
	    	double goniometerHorizontal=parVect[6];
	    	double goniometerAxial=parVect[7];
	    	double interAxisDistance=parVect[8];
	    	double interAxisAngle=parVect[9];
	    	double horAxisErrPhi=parVect[10];
	    	double horAxisErrPsi=parVect[11];
	    	double entrancePupilForward=parVect[12];
	    	double centerAboveHorizontal=parVect[13];
	    	double GXYZ0=parVect[14];
	    	double GXYZ1=parVect[15];
	    	double GXYZ2=parVect[16];
	    	
	    	double cPS=   Math.cos(psi*Math.PI/180); // subCam.psi
	    	double sPS=   Math.sin(psi*Math.PI/180); // subCam.psi
	    	double cTH=   Math.cos(theta*Math.PI/180); // subCam.theta
	    	double sTH=   Math.sin(theta*Math.PI/180); // subCam.theta
	    	double cAZP=  Math.cos((azimuth+phi)*Math.PI/180); //subCam.azimuth+subCam.phi
	    	double sAZP=  Math.sin((azimuth+phi)*Math.PI/180); //subCam.azimuth+subCam.phi
	    	double cAZ=   Math.cos(azimuth*Math.PI/180); //subCam.azimuth
	    	double sAZ=   Math.sin(azimuth*Math.PI/180); //subCam.azimuth
	    	double cGA=   Math.cos(goniometerAxial*Math.PI/180); //eyesisCameraParameters.goniometerAxial
	    	double sGA=   Math.sin(goniometerAxial*Math.PI/180); //eyesisCameraParameters.goniometerAxial
	    	double cGH=   Math.cos(goniometerHorizontal*Math.PI/180); //eyesisCameraParameters.goniometerHorizontal
	    	double sGH=   Math.sin(goniometerHorizontal*Math.PI/180); //eyesisCameraParameters.goniometerHorizontal
	    	double cGIAA= Math.cos(interAxisAngle*Math.PI/180); // eyesisCameraParameters.interAxisAngle
	    	double sGIAA= Math.sin(interAxisAngle*Math.PI/180); //eyesisCameraParameters.interAxisAngle
	    	double cHAEPH=Math.cos(horAxisErrPhi*Math.PI/180); //eyesisCameraParameters.horAxisErrPhi
	    	double sHAEPH=Math.sin(horAxisErrPhi*Math.PI/180); //eyesisCameraParameters.horAxisErrPhi
	    	double cHAEPS=Math.cos(horAxisErrPsi*Math.PI/180); //eyesisCameraParameters.horAxisErrPsi
	    	double sHAEPS=Math.sin(horAxisErrPsi*Math.PI/180); //eyesisCameraParameters.horAxisErrPsi\
	    	
	/*
	 0) Translate by distance to entrance pupil (lens center)
	    	| Xc0 |   | 0                     |   |Xc|
	    	| Yc0 | = | 0                     | + |Yc|
	    	| Zc0 |   | entrancePupilForward  |   |Zc|
	    	*/
	    	double [][] aT0={{0.0},{0.0},{entrancePupilForward}};
	    	Matrix T0=new Matrix(aT0);
	    	
	/*
	Converting from the sub-camera coordinates to the target coordinates
	1) rotate by -psi around CZ: Vc1= R1*Vc
	| Xc1 |   | cos(psi)  sin(psi)    0         |   |Xc0|
	| Yc1 | = |-sin(psi)  cos(psi)    0         | * |Yc0|
	| Zc1 |   |    0         0        1         |   |Zc0|
	*/
	    	double [][] aR1={{cPS,sPS,0.0},{-sPS,cPS,0.0},{0.0,0.0,1.0}};
	    	Matrix R1=new Matrix(aR1);
	/*    	
	2) rotate by - theta around C1X:Vc2= R2*Vc1
	| Xc2 |   |    1         0         0        |   |Xc1|
	| Yc2 | = |    0    cos(theta)   sin(theta) | * |Yc1|
	| Zc2 |   |    0   -sin(theta)   cos(theta) |   |Zc1|
	*/
	    	double [][] aR2={{1.0,0.0,0.0},{0.0,cTH,sTH},{0.0,-sTH,cTH}};
	    	Matrix R2=new Matrix(aR2);
	/*    	
	3) rotate by -(azimuth+phi) around C2Y:Vc3= R3*Vc2
	| Xc3 |   | cos(azimuth+phi)    0   sin(azimuth+phi)   |   |Xc2|
	| Yc3 | = |     0               1         0            | * |Yc2|
	| Zc3 |   | -sin(azimuth+phi)   0   cos(azimuth+phi)   |   |Zc2|
	*/
	    	double [][] aR3={{cAZP,0.0,sAZP},{0.0,1.0,0.0},{-sAZP,0.0,cAZP}};
	    	Matrix R3=new Matrix(aR3);
	/*    	
	4) Now axes are aligned, just translate to get to eyesis coordinates: Vey= T1+Vc3
	| Xey |   |      r * sin (azimuth)       |   |Xc3|
	| Yey | = | height+centerAboveHorizontal | + |Yc3|
	| Zey |   |      r * cos (azimuth)       |   |Zc3|
	*/
	    	double [][] aT1={{radius*sAZ},{(height+centerAboveHorizontal)},{radius*cAZ}}; // {{subCam.radius*sAZ},{subCam.height},{subCam.radius*cAZ}};
	    	Matrix T1=new Matrix(aT1);
	    	
	/**
	5) rotate around moving goniometer axis, by (-goniometerAxial) - same as around EY: Vgm1=R4*Vey
	| Xgm1 |   |    1           0                     0              |   |Xey|
	| Ygm1 | = |    0   cos(goniometerAxial)  sin(goniometerAxial)   | * |Yey|
	| Zgm1 |   |    0   -sin(goniometerAxial) cos(goniometerAxial)   |   |Zey|
	*/
	    	double [][] aR4_tripod=    {{1.0,0.0,0.0},{0.0,cGA,sGA},{0.0,-sGA,cGA}};
	/*    	
	5) rotate around moving goniometer axis, by (-goniometerAxial) - same as around EY: Vgm1=R4*Vey
	| Xgm1 |   | cos(goniometerAxial)    0   sin(goniometerAxial)   |   |Xey|
	| Ygm1 | = |        0                1             0            | * |Yey|
	| Zgm1 |   |-sin(goniometerAxial)    0   cos(goniometerAxial)   |   |Zey|
	*/
			double [][] aR4_goniometer={{cGA,0.0,sGA},{0.0,1.0,0.0},{-sGA,0.0,cGA}};
			Matrix R4=new Matrix(isTripod?aR4_tripod:aR4_goniometer);
	/*    	
	6) move to the goniometer horizontal axis:Vgm2=T2+Vgm1
	| Xgm2 |   |                 0  |   |Xgm1|
	| Ygm2 | = |                 0  | + |Ygm1|
	| Zgm2 |   | interAxisDistance  |   |Zgm1|
	*/
	    	double [][] aT2={{0.0},{0.0},{interAxisDistance}}; //eyesisCameraParameters.interAxisDistance
	    	Matrix T2=new Matrix(aT2);
	/*    	
	7) rotate around Zgm2 by -interAxisAngle, so Xgm3 axis is the same as horizontal goniometer axis: Vgm3=R5*Vgm2
	| Xgm3 |   | cos(interAxisAngle)  sin(interAxisAngle)    0         |   |Xgm2|
	| Ygm3 | = |-sin(interAxisAngle)  cos(interAxisAngle)    0         | * |Ygm2|
	| Zgm3 |   |    0                         0              1         |   |Zgm2|
	*/
	    	double [][] aR5={{cGIAA,sGIAA,0.0},{-sGIAA,cGIAA,0.0},{0.0,0.0,1.0}};
	    	Matrix R5=new Matrix(aR5);
	/**
	8)  rotate around goniometer horizontal axis (Xgm3) by -goniometerHorizontal: Vgm4=R6*Vgm3
	| Xgm4 |   |    cos(goniometerHorizontal)  0   sin(goniometerHorizontal) |   |Xgm3|
	| Ygm4 | = |               0               1               0             | * |Ygm3|
	| Zgm4 |   |   -sin(goniometerHorizontal)  0   cos(goniometerHorizontal) |   |Zgm3|
	*/
	    	double [][] aR6_tripod=    {{cGH,0.0,sGH},{0.0,1.0,0.0},{-sGH,0.0,cGH}};
	/*    	
	8)  rotate around goniometer horizontal axis (Xgm3) by -goniometerHorizontal: Vgm4=R6*Vgm3
	| Xgm4 |   |    1                 0                           0            |   |Xgm3|
	| Ygm4 | = |    0    cos(goniometerHorizontal)   sin(goniometerHorizontal) | * |Ygm3|
	| Zgm4 |   |    0   -sin(goniometerHorizontal)   cos(goniometerHorizontal) |   |Zgm3|
	*/
	    	double [][] aR6_goniometer={{1.0,0.0,0.0},{0.0,cGH,sGH},{0.0,-sGH,cGH}};
			Matrix R6=new Matrix(isTripod?aR6_tripod:aR6_goniometer);
	/*    	
	9) Correct roll error of the goniometer horizontal axis - rotate by -horAxisErrPsi around Zgm4: Vgm5=R7*Vgm4
	| Xgm5 |   | cos(horAxisErrPsi)  sin(horAxisErrPsi)    0         |   |Xgm4|
	| Ygm5 | = |-sin(horAxisErrPsi)  cos(horAxisErrPsi)    0         | * |Ygm4|
	| Zgm5 |   |         0                   0             1         |   |Zgm4|
	*/
	    	double [][] aR7={{cHAEPS,sHAEPS,0.0},{-sHAEPS,cHAEPS,0.0},{0.0,0.0,1.0}};
	    	Matrix R7=new Matrix(aR7);
	/*    	
	10) Correct azimuth error of the goniometer hoirizontal axis - rotate by -horAxisErrPhi around Ygm5: Vgm6=R8*Vgm5
	 
	| Xgm6 |   | cos(horAxisErrPhi)      0   sin(horAxisErrPhi)   |   |Xgm5|
	| Ygm6 | = |        0                1             0          | * |Ygm5|
	| Zgm6 |   |-sin(horAxisErrPhi)      0   cos(horAxisErrPhi)   |   |Zgm5|
	For Tripod - rotate around X-axis (like theta)
	| Xgm6 |   |   1           0                       0      |   |Xgm5|
	| Ygm6 | = |   0  cos(horAxisErrPhi)  sin(horAxisErrPhi)  | * |Ygm5|
	| Zgm6 |   |   0 -sin(horAxisErrPhi)  cos(horAxisErrPhi)  |   |Zgm5|

	*/
	    	double [][] aR8_tripod=    {{1.0,   0.0,   0.0},{0.0,cHAEPH,sHAEPH},{0.0,  -sHAEPH,cHAEPH}};
	    	double [][] aR8_goniometer={{cHAEPH,0.0,sHAEPH},{0.0,   1.0,   0.0},{-sHAEPH,  0.0,cHAEPH}};
	    	Matrix R8=new Matrix(isTripod?aR8_tripod:aR8_goniometer);
	/*    	
	11) translate to the target zero point: Vt=  T3+Vgm6
	| Xt |   | GXYZ[0]  |   |Xgm6|
	| Yt | = |-GXYZ[1]  | + |Ygm6| // Y - up positive
	| Zt |   |-GXYZ[2]  |   |Zgm6| // Z - away positive
	    	 */
//	    	double [][] aT3={{parVect[12]},{-parVect[13]},{-parVect[14]}};//{{eyesisCameraParameters.GXYZ[0]},{eyesisCameraParameters.GXYZ[1]},{eyesisCameraParameters.GXYZ[2]}};
	   	
	    	double [][] aT3={{GXYZ0},{-GXYZ1},{-GXYZ2}};//{{eyesisCameraParameters.GXYZ[0]},{eyesisCameraParameters.GXYZ[1]},{eyesisCameraParameters.GXYZ[2]}};
//	    	double [][] aT3={{parVect[12]},{ parVect[13]},{-parVect[14]}};//{{eyesisCameraParameters.GXYZ[0]},{eyesisCameraParameters.GXYZ[1]},{eyesisCameraParameters.GXYZ[2]}};
	    	Matrix T3=new Matrix(aT3);

	// MA=R8*R7*R6*R5*R4*R3*R2*R1;
	// MB=T3+(R8*R7*R6*R5*(T2+R4*T1)); - old
	// MB=T3+(R8*R7*R6*R5*(T2+R4*(T1+R3*R2*R1*T0)));
	        Matrix MA=R8.times(R7.times(R6.times(R5.times(R4.times(R3.times(R2.times(R1)))))));
//	        Matrix MB=T3.plus(R8.times(R7.times(R6.times(R5.times(T2.plus(R4.times(T1)))))));
	        Matrix MB=T3.plus(R8.times(R7.times(R6.times(R5.times(T2.plus(R4.times(T1.plus(R3.times(R2.times(R1.times(T0)))) )))))));
			if (this.debugLevel>2) {
				System.out.println("MA:");
				MA.print(10, 5);
				System.out.println("MB:");
				MB.print(10, 5);
				System.out.println("T3:");
				T3.print(10, 5);
//				System.out.println("interParameterDerivatives[0]="+sprintfArray(interParameterDerivatives[0]));
			}
	        double [] extrinsicParams=parametersFromMAMB(MA,MB); // all after 6 are 0;
	        lensDistortionParameters.distance=extrinsicParams[2];
	        lensDistortionParameters.x0=extrinsicParams[0];
	        lensDistortionParameters.y0=extrinsicParams[1];
	        lensDistortionParameters.z0=0.0; // used here
	        lensDistortionParameters.pitch=extrinsicParams[4];
	        lensDistortionParameters.yaw=  extrinsicParams[3];
	        lensDistortionParameters.roll= extrinsicParams[5];
//			lensDistortionParameters.focalLength=parVect[15]; //subCam.focalLength;
//			lensDistortionParameters.px0=parVect[16]; //subCam.px0;
//			lensDistortionParameters.py0=parVect[17]; //subCam.py0;
//			lensDistortionParameters.distortionA5=parVect[18]; //subCam.distortion5;
//			lensDistortionParameters.distortionA=parVect[19]; //subCam.distortionA;
//			lensDistortionParameters.distortionB=parVect[20]; //subCam.distortionB;
//			lensDistortionParameters.distortionC=parVect[21]; //subCam.distortionC;

			lensDistortionParameters.focalLength=parVect[17]; //subCam.focalLength;
			lensDistortionParameters.px0=parVect[18]; //subCam.px0;
			lensDistortionParameters.py0=parVect[19]; //subCam.py0;
			lensDistortionParameters.distortionA8=parVect[20]; //subCam.distortion5;
			lensDistortionParameters.distortionA7=parVect[21]; //subCam.distortion5;
			lensDistortionParameters.distortionA6=parVect[22]; //subCam.distortion5;
			lensDistortionParameters.distortionA5=parVect[23]; //subCam.distortion5;
			lensDistortionParameters.distortionA=parVect[24]; //subCam.distortionA;
			lensDistortionParameters.distortionB=parVect[25]; //subCam.distortionB;
			lensDistortionParameters.distortionC=parVect[26]; //subCam.distortionC;

			
			lensDistortionParameters.recalcCommons();
	        
	        
	        if (this.debugLevel>2){
	        	System.out.println("lensDistortionParameters.recalcCommons()");
	        	System.out.println("lensDistortionParameters.distance="+IJ.d2s(lensDistortionParameters.distance, 3));
	        	System.out.println("lensDistortionParameters.x0="+      IJ.d2s(lensDistortionParameters.x0, 3));
	        	System.out.println("lensDistortionParameters.y0="+      IJ.d2s(lensDistortionParameters.y0, 3));
	        	System.out.println("lensDistortionParameters.z0="+      IJ.d2s(lensDistortionParameters.z0, 3));
	        	System.out.println("lensDistortionParameters.pitch="+   IJ.d2s(lensDistortionParameters.pitch, 3));
	        	System.out.println("lensDistortionParameters.yaw="+IJ.d2s(lensDistortionParameters.yaw, 3));
	        	System.out.println("lensDistortionParameters.roll="+IJ.d2s(lensDistortionParameters.roll, 3));
	        	System.out.println("lensDistortionParameters.focalLength="+IJ.d2s(lensDistortionParameters.focalLength, 3));
	        	System.out.println("lensDistortionParameters.px0="+IJ.d2s(lensDistortionParameters.px0, 3));
	        	System.out.println("lensDistortionParameters.py0="+IJ.d2s(lensDistortionParameters.py0, 3));
	        	System.out.println("lensDistortionParameters.distortionA8="+IJ.d2s(lensDistortionParameters.distortionA8, 5));
	        	System.out.println("lensDistortionParameters.distortionA7="+IJ.d2s(lensDistortionParameters.distortionA7, 5));
	        	System.out.println("lensDistortionParameters.distortionA6="+IJ.d2s(lensDistortionParameters.distortionA6, 5));
	        	System.out.println("lensDistortionParameters.distortionA5="+IJ.d2s(lensDistortionParameters.distortionA5, 5));
	        	System.out.println("lensDistortionParameters.distortionA="+IJ.d2s(lensDistortionParameters.distortionA, 5));
	        	System.out.println("lensDistortionParameters.distortionB="+IJ.d2s(lensDistortionParameters.distortionB, 5));
	        	System.out.println("lensDistortionParameters.distortionC="+IJ.d2s(lensDistortionParameters.distortionC, 5));
	        }
	        if (!calculateDerivatives) return; 
	/* Calculate all derivativs as a matrix.
	 *  Input parameters (columns):
	0    	public double azimuth; // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
	1    	public double radius;  // mm, distance from the rotation axis
	2   	public double height;       // mm, up (was downwards?) - from the origin point
	3   	public double phi;     // degrees, optical axis from azimuth/r vector, clockwise
	4   	public double theta;   // degrees, optical axis from the eyesis horizon, positive - up
	5   	public double psi;     // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
	6   	public double goniometerHorizontal; // goniometer rotation around "horizontal" axis (tilting from the target - positive)
	7   	public double goniometerAxial; // goniometer rotation around Eyesis axis (clockwise in plan - positive 
	8   	public double interAxisDistance; // distance in mm between two goniometer axes
	9    	public double interAxisAngle;    // angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
	10   	public double horAxisErrPhi;   // angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)
	11   	public double horAxisErrPsi;   // angle in degrees "horizontal" goniometer axis is rotated around moving X axis (up)
	Two new parameters
	12    	public double entrancePupilForward; // common to all lenses - distance from the sensor to the lens entrance pupil
	13    	public double centerAboveHorizontal; // camera center distance along camera axis above the closest point to horizontal rotation axis (adds to 
	14(12)  x	public double [] GXYZ=null; // coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system 
	15(13)  y
	16(14)  z
	17(15)	public double focalLength=4.5;
	18(16)	public double px0=1296.0;          // center of the lens on the sensor, pixels
	19(17)		public double py0=968.0;           // center of the lens on the sensor, pixels
	20(18)		public double distortionA8=0.0; // r^8 (normalized to focal length or to sensor half width?)
	21(19)		public double distortionA7=0.0; // r^7 (normalized to focal length or to sensor half width?)
	22(20)		public double distortionA6=0.0; // r^6 (normalized to focal length or to sensor half width?)
	23(21)		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
	24(22)		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
	25(23)		public double distortionB=0.0; // r^3
	26(24)		public double distortionC=0.0; // r^2

	 * Output parameters (rows):
	0		public double x0=0;      // lens axis from pattern center, mm (to the right)
	1		public double y0=0;      // lens axis from pattern center, mm (down)
	2		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
	3		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
	4		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
	5		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
	6		public double focalLength=4.5;
	7		public double px0=1296.0;           // center of the lens on the sensor, pixels
	8		public double py0=968.0;           // center of the lens on the sensor, pixels
	 		public double distortionRadius=  2.8512; // mm - half width of the sensor
	9		public double distortionA8=0.0; // r^8 (normalized to focal length or to sensor half width?)
	10		public double distortionA7=0.0; // r^7 (normalized to focal length or to sensor half width?)
	11		public double distortionA6=0.0; // r^6 (normalized to focal length or to sensor half width?)
	12		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
	13		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
	14		public double distortionB=0.0; // r^3
	15		public double distortionC=0.0; // r^2
	 */
//	        interParameterDerivatives=new double[getNumInputs()][]; //partial derivative matrix from subcamera-camera-goniometer to single camera (12x21)
	        
	//calculate dMA_azimuth        
	//calculate dMB_azimuth        
	/*
	// MA=R8*R7*R6*R5*R4*R3*R2*R1;
	// MB=T3+(R8*R7*R6*R5*(T2+R4*T1)); - old
	// MB=T3+(R8*R7*R6*R5*(T2+R4*(T1+R3*R2*R1*T0)));
	        Matrix MA=R8.times(R7.times(R6.times(R5.times(R4.times(R3.times(R2.times(R1)))))));
	        Matrix MB=T3.plus(R8.times(R7.times(R6.times(R5.times(T2.plus(R4.times(T1)))))));  old
	        Matrix MB=T3.plus(R8.times(R7.times(R6.times(R5.times(T2.plus(R4.times(T1.plus(R3.times(R2.times(R1.times(T0)))) )))))));
	3) rotate by -(azimuth+phi) around C2Y:Vc3= R3*Vc2
	| Xc3 |   | cos(azimuth+phi)    0   sin(azimuth+phi)   |   |Xc2|
	| Yc3 | = |     0               1         0            | * |Yc2|
	| Zc3 |   | -sin(azimuth+phi)   0   cos(azimuth+phi)   |   |Zc2|
	    	double [][] aR3={{cAZP,0.0,sAZP},{0.0,1.0,0.0},{-sAZP,0.0,cAZP}};
	    	Matrix R3=new Matrix(aR3);
	4) Now axes are aligned, just translate to get to eyesis coordinates: Vey= T1+Vc3
	| Xey |   |  r * sin (azimuth)  |   |Xc3|
	| Yey | = |             height  | + |Yc3|
	| Zey |   |  r * cos (azimuth)  |   |Zc3|
	    	double [][] aT1={{subCam.radius*sAZ},{subCam.height},{subCam.radius*cAZ}};
	    	Matrix T1=new Matrix(aT1);

	        Matrix 
	// Make a function MA, MB - >parameters (column) - reuse it above and for each    interParameterDerivatives row
	Matrix dMA_azimuth=R8*R7*R6*R5*R4*dR3_azimuth*R2*R1;
	Matrix dMB_azimuth=T3+(R8*R7*R6*R5*(R4*dT1_azimuth));

	Use extrinsicParams=parametersFromMAMB(dMA_azimuth,dMB_azimuth);
		       		if (this.debugLevel>3) {
		    			System.out.println("calculateFxAndJacobian->calcPartialDerivatives("+IJ.d2s(targetXYZ[fullIndex][0],2)+","+
		    					IJ.d2s(targetXYZ[fullIndex][1],2)+","+
		    					IJ.d2s(targetXYZ[fullIndex][2],2)+" ("+calcJacobian+") -> "+
		    					IJ.d2s(derivatives14[0][0],2)+"/"+IJ.d2s(derivatives14[0][1],2));
		    		}

	Which parameters affect which matrices
	                                               R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 || T0 | T1 | T2 | T3 |
	0    	public double azimuth;              //    |    | +  |    |    |    |    |    ||    | +  |    |    |
	1    	public double radius;               //    |    |    |    |    |    |    |    ||    | +  |    |    |
	2   	public double height;               //    |    |    |    |    |    |    |    ||    | +  |    |    |
	3   	public double phi;                  //    |    | +  |    |    |    |    |    ||    |    |    |    |
	4   	public double theta;                //    | +  |    |    |    |    |    |    ||    |    |    |    |
	5   	public double psi;                  // +  |    |    |    |    |    |    |    ||    |    |    |    |
	6   	public double goniometerHorizontal; //    |    |    |    |    | +  |    |    ||    |    |    |    |
	7   	public double goniometerAxial;      //    |    |    | +  |    |    |    |    ||    |    |    |    |
	8   	public double interAxisDistance;    //    |    |    |    |    |    |    |    ||    |    | +  |    |
	9    	public double interAxisAngle;       //    |    |    |    | +  |    |    |    ||    |    |    |    |
	10   	public double horAxisErrPhi;        //    |    |    |    |    |    |    | +  ||    |    |    |    |
	11   	public double horAxisErrPsi;        //    |    |    |    |    |    | +  |    ||    |    |    |    |
	12      public double entrancePupilForward  //    |    |    |    |    |    |    |    || +  |    |    |    |
	13      public double centerAboveHorizontal //    |    |    |    |    |    |    |    ||    | +  |    |    |
	14  x	public double [] GXYZ=null;         //    |    |    |    |    |    |    |    ||    |    |    | +  |
	15  y                                       //    |    |    |    |    |    |    |    ||    |    |    | +  |
	16  z                                       //    |    |    |    |    |    |    |    ||    |    |    | +  |

	 */
	// Below can be optimized with common intermediate results
	   		if (this.debugLevel>2) {
	   			for (int i=0;i<parVect.length;i++){
	   				System.out.println("calcInterParamers(): parVect["+i+"]="+parVect[i]);
	   			}
	   		}
	//0    	public double azimuth; // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
	        if (mask[0]) {
	        	double [][] adR3_azimuth={{-sAZP,0.0,cAZP},{0.0,0.0,0.0},{-cAZP,0.0,-sAZP}};
	        	Matrix dR3_azimuth=new Matrix(adR3_azimuth);
//	        	double [][] adT1_azimuth={{radius*cAZ},{height},{-radius*sAZ}}; //{{subCam.radius*cAZ},{subCam.height},{-subCam.radius*sAZ}}
	        	double [][] adT1_azimuth={{radius*cAZ},{0.0},{-radius*sAZ}}; //{{subCam.radius*cAZ},{subCam.height},{-subCam.radius*sAZ}}
	        	Matrix dT1_azimuth=new Matrix(adT1_azimuth);
	        	
	        	Matrix dMA_azimuth=R8.times(R7.times(R6.times(R5.times(R4.times(dR3_azimuth.times(R2.times(R1)))))));
	        	Matrix dMB0_azimuth=R8.times(R7.times(R6.times(R5.times(R4.times(dT1_azimuth)))));
	        	Matrix dMB_azimuth=dMB0_azimuth.plus(dMA_azimuth.times(T0)); // new term
	        	interParameterDerivatives[0]=d_parametersFromMAMB(dMA_azimuth,dMB_azimuth,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_azimuth:");
	    			dMA_azimuth.print(10, 5);
	    			System.out.println("dMB_azimuth:");
	    			dMB_azimuth.print(10, 5);
	    			System.out.println("interParameterDerivatives[0]="+sprintfArray(interParameterDerivatives[0]));
	    		}
	        	
	        } else interParameterDerivatives[0]=null;
	//1    	public double radius;  // mm, distance from the rotation axis
	        if (mask[1]) {
	        	double [][] adT1_radius={{sAZ},{0.0},{cAZ}}; //{{subCam.radius*sAZ},{0.0},{subCam.radius*cAZ}}
	        	Matrix dT1_radius=new Matrix(adT1_radius);
	        	Matrix dMA_radius=new Matrix(3,3,0.0);
	        	Matrix dMB_radius=R8.times(R7.times(R6.times(R5.times(R4.times(dT1_radius)))));
	        	interParameterDerivatives[1]=d_parametersFromMAMB(dMA_radius,dMB_radius,MA,MB,false); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_radius:");
	    			dMA_radius.print(10, 5);
	    			System.out.println("dMB_radius:");
	    			dMB_radius.print(10, 5);
	    			System.out.println("interParameterDerivatives[1]="+sprintfArray(interParameterDerivatives[1]));
	    		}
	        } else interParameterDerivatives[1]=null; 
	//2   	public double height;       // mm, downwards - from the origin point
	        if (mask[2]) {
	        	double [][] adT1_height={{0.0},{1.0},{0.0}};
	        	Matrix dT1_height=new Matrix(adT1_height);
	        	Matrix dMA_height=new Matrix(3,3,0.0);
	        	Matrix dMB_height=R8.times(R7.times(R6.times(R5.times(R4.times(dT1_height)))));
	        	interParameterDerivatives[2]=d_parametersFromMAMB(dMA_height,dMB_height,MA,MB,false); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_height:");
	    			dMA_height.print(10, 5);
	    			System.out.println("dMB_height:");
	    			dMB_height.print(10, 5);
	    			System.out.println("interParameterDerivatives[2]="+sprintfArray(interParameterDerivatives[2]));
	    		}
	        } else interParameterDerivatives[2]=null; 
	//3   	public double phi;     // degrees, optical axis from azimuth/r vector, clockwise
	        if (mask[3]) {
	        	double [][] adR3_phi={{-sAZP,0.0,cAZP},{0.0,0.0,0.0},{-cAZP,0.0,-sAZP}}; // same as adR3_azimuth
	        	Matrix dR3_phi=new Matrix(adR3_phi); // same as dR3_azimuth
	        	Matrix dMA_phi=R8.times(R7.times(R6.times(R5.times(R4.times(dR3_phi.times(R2.times(R1))))))); //same as dMA_azimuth
//	        	Matrix dMB_phi=new Matrix(3,1,0.0);
	        	Matrix dMB_phi=dMA_phi.times(T0); // new term
	        	interParameterDerivatives[3]=d_parametersFromMAMB(dMA_phi,dMB_phi,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_phi:");
	    			dMA_phi.print(10, 5);
	    			System.out.println("dMB_phi:");
	    			dMB_phi.print(10, 5);
	    			System.out.println("interParameterDerivatives[3]="+sprintfArray(interParameterDerivatives[3]));
	    		}
	        } else interParameterDerivatives[3]=null; 
	//4   	public double theta;   // degrees, optical axis from the eyesis horizon, positive - up
	        if (mask[4]) {
	        	double [][] adR2_theta={{0.0,0.0,0.0},{0.0,-sTH,cTH},{0.0,-cTH,-sTH}};
	        	Matrix dR2_theta=new Matrix(adR2_theta);
	        	Matrix dMA_theta=R8.times(R7.times(R6.times(R5.times(R4.times(R3.times(dR2_theta.times(R1)))))));
//	        	Matrix dMB_theta=new Matrix(3,1,0.0);
	        	Matrix dMB_theta=dMA_theta.times(T0); // new term
	        	interParameterDerivatives[4]=d_parametersFromMAMB(dMA_theta,dMB_theta,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_theta:");
	    			dMA_theta.print(10, 5);
	    			System.out.println("dMB_theta:");
	    			dMB_theta.print(10, 5);
	    			System.out.println("interParameterDerivatives[4]="+sprintfArray(interParameterDerivatives[4]));
	    		}
	        } else interParameterDerivatives[4]=null; 
	//5   	public double psi;     // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
	        if (mask[5]) {
	        	double [][] adR1_psi={{-sPS,cPS,0.0},{-cPS,-sPS,0.0},{0.0,0.0,0.0}};
	        	Matrix dR1_psi=new Matrix(adR1_psi);
	        	Matrix dMA_psi=R8.times(R7.times(R6.times(R5.times(R4.times(R3.times(R2.times(dR1_psi)))))));
//	        	Matrix dMB_psi=new Matrix(3,1,0.0);
	        	Matrix dMB_psi=dMA_psi.times(T0); // new term
	        	interParameterDerivatives[5]=d_parametersFromMAMB(dMA_psi,dMB_psi,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	/*    			
	    			System.out.print("R1:");
	    			R1.print(10, 5);
	    			System.out.print("dR1_psi:");
	    			dR1_psi.print(10, 5);
	    			Matrix R82_psi=R8.times(R7.times(R6.times(R5.times(R4.times(R3.times(R2))))));
	    			System.out.print("R82_psi:");
	    			R82_psi.print(10, 5);
	    			*/
	    			System.out.print("dMA_psi:");
	    			dMA_psi.print(10, 5);
	    			System.out.print("dMB_psi:");
	    			dMB_psi.print(10, 5);
	    			System.out.println("interParameterDerivatives[5]="+sprintfArray(interParameterDerivatives[5]));
	    		}
	        } else interParameterDerivatives[5]=null; 
	//6   	public double goniometerHorizontal; // goniometer rotation around "horizontal" axis (tilting from the target - positive)
	        if (mask[6]) {
	/* define for isTripod */
	        	double [][] adR6_goniometerHorizontal_tripod=    {{-sGH,0.0,cGH},{0.0, 0.0,0.0},{-cGH, 0.0,-sGH}};
	        	double [][] adR6_goniometerHorizontal_goniometer={{ 0.0,0.0,0.0},{0.0,-sGH,cGH},{ 0.0,-cGH,-sGH}};
	        	Matrix dR6_goniometerHorizontal=new Matrix(isTripod?adR6_goniometerHorizontal_tripod:adR6_goniometerHorizontal_goniometer);
	        	Matrix dMA_goniometerHorizontal=R8.times(R7.times(dR6_goniometerHorizontal.times(R5.times(R4.times(R3.times(R2.times(R1)))))));
	        	Matrix dMB_goniometerHorizontal=R8.times(R7.times(dR6_goniometerHorizontal.times(R5.times(T2.plus(R4.times(T1))))));
	        	interParameterDerivatives[6]=d_parametersFromMAMB(dMA_goniometerHorizontal,dMB_goniometerHorizontal,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_goniometerHorizontal:");
	    			dMA_goniometerHorizontal.print(10, 5);
	    			System.out.println("dMB_goniometerHorizontal:");
	    			dMB_goniometerHorizontal.print(10, 5);
	    			System.out.println("interParameterDerivatives[6]="+sprintfArray(interParameterDerivatives[6]));
	    		}
	        } else interParameterDerivatives[6]=null; 
	//7   	public double goniometerAxial; // goniometer rotation around Eyesis axis (clockwise in plan - positive 
	        if (mask[7]) {
	// define for isTripod
	        	double [][] adR4_goniometerAxial_tripod=    {{ 0.0,0.0,0.0},{0.0,-sGA,cGA},{ 0.0,-cGA,-sGA}};
	        	double [][] adR4_goniometerAxial_goniometer={{-sGA,0.0,cGA},{0.0, 0.0,0.0},{-cGA, 0.0,-sGA}};
	        	Matrix dR4_goniometerAxial=new Matrix(isTripod?adR4_goniometerAxial_tripod:adR4_goniometerAxial_goniometer);
	        	Matrix dMA_goniometerAxial=R8.times(R7.times(R6.times(R5.times(dR4_goniometerAxial.times(R3.times(R2.times(R1)))))));
	        	Matrix dMB_goniometerAxial=R8.times(R7.times(R6.times(R5.times(dR4_goniometerAxial.times(T1)))));
	        	interParameterDerivatives[7]=d_parametersFromMAMB(dMA_goniometerAxial,dMB_goniometerAxial,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_goniometerAxial:");
	    			dMA_goniometerAxial.print(10, 5);
	    			System.out.println("dMB_goniometerAxial:");
	    			dMB_goniometerAxial.print(10, 5);
	    			System.out.println("interParameterDerivatives[7]="+sprintfArray(interParameterDerivatives[7]));
	    		}
	        } else interParameterDerivatives[7]=null; 
	//8   	public double interAxisDistance; // distance in mm between two goniometer axes
	        if (mask[8]) {
	        	double [][] adT2_interAxisDistance={{0.0},{0.0},{1.0}};
	        	Matrix dT2_interAxisDistance=new Matrix(adT2_interAxisDistance);
	        	Matrix dMA_interAxisDistance=new Matrix(3,3,0.0);
	        	Matrix dMB_interAxisDistance=R8.times(R7.times(R6.times(R5.times(dT2_interAxisDistance))));
	        	interParameterDerivatives[8]=d_parametersFromMAMB(dMA_interAxisDistance,dMB_interAxisDistance,MA,MB,false); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_interAxisDistance:");
	    			dMA_interAxisDistance.print(10, 5);
	    			System.out.println("dMB_interAxisDistance:");
	    			dMB_interAxisDistance.print(10, 5);
	    			System.out.println("interParameterDerivatives[8]="+sprintfArray(interParameterDerivatives[8]));
	    		}
	        } else interParameterDerivatives[8]=null; 
	//9    	public double interAxisAngle;    // angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
	        if (mask[9]) {
	        	double [][] adR5_interAxisAngle={{-sGIAA,cGIAA,0.0},{-cGIAA,-sGIAA,0.0},{0.0,0.0,0.0}};
	        	Matrix dR5_interAxisAngle=new Matrix(adR5_interAxisAngle);
	        	Matrix dMA_interAxisAngle=R8.times(R7.times(R6.times(dR5_interAxisAngle.times(R4.times(R3.times(R2.times(R1)))))));
	        	Matrix dMB_interAxisAngle=R8.times(R7.times(R6.times(dR5_interAxisAngle.times(T2.plus(R4.times(T1))))));    	
	        	interParameterDerivatives[9]=d_parametersFromMAMB(dMA_interAxisAngle,dMB_interAxisAngle,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_interAxisAngle:");
	    			dMA_interAxisAngle.print(10, 5);
	    			System.out.println("dMB_interAxisAngle:");
	    			dMB_interAxisAngle.print(10, 5);
	    			System.out.println("interParameterDerivatives[9]="+sprintfArray(interParameterDerivatives[9]));
	    		}
	        } else interParameterDerivatives[9]=null; 
	//10   	public double horAxisErrPhi;   // angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)
	        
	// change sHAEPH to rotate like theta /1 0 0 /0 c s/ 0 -s c/ (derivative -/0 0 0 /0 -s c/ 0 -c -s/         
	        if (mask[10]) {
	        	double [][] adR8_horAxisErrPhi_tripod=    {{0.0,   0.0,    0.0},{0.0,-sHAEPH,cHAEPH},{    0.0,-cHAEPH,-sHAEPH}};
	        	double [][] adR8_horAxisErrPhi_goniometer={{-sHAEPH,0.0,cHAEPH},{0.0,    0.0,   0.0},{-cHAEPH,    0.0,-sHAEPH}};
	        	Matrix dR8_horAxisErrPhi=new Matrix(isTripod?adR8_horAxisErrPhi_tripod:adR8_horAxisErrPhi_goniometer);
	        	Matrix dMA_horAxisErrPhi=dR8_horAxisErrPhi.times(R7.times(R6.times(R5.times(R4.times(R3.times(R2.times(R1)))))));
	        	Matrix dMB_horAxisErrPhi=dR8_horAxisErrPhi.times(R7.times(R6.times(R5.times(T2.plus(R4.times(T1))))));
	        	interParameterDerivatives[10]=d_parametersFromMAMB(dMA_horAxisErrPhi,dMB_horAxisErrPhi,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_horAxisErrPhi:");
	    			dMA_horAxisErrPhi.print(10, 5);
	    			System.out.println("dMB_horAxisErrPhi:");
	    			dMB_horAxisErrPhi.print(10, 5);
	    			System.out.println("interParameterDerivatives[10]="+sprintfArray(interParameterDerivatives[10]));
	    		}
	        } else interParameterDerivatives[10]=null; 
	//11   	public double horAxisErrPsi;   // angle in degrees "horizontal" goniometer axis is rotated around moving X axis (up)
	        if (mask[11]) {
	        	double [][] adR7_horAxisErrPsi={{-sHAEPS,cHAEPS,0.0},{-cHAEPS,-sHAEPS,0.0},{0.0,0.0,0.0}};
	        	Matrix dR7_horAxisErrPsi=new Matrix(adR7_horAxisErrPsi);
	        	Matrix dMA_horAxisErrPsi=R8.times(dR7_horAxisErrPsi.times(R6.times(R5.times(R4.times(R3.times(R2.times(R1)))))));
	        	Matrix dMB_horAxisErrPsi=R8.times(dR7_horAxisErrPsi.times(R6.times(R5.times(T2.plus(R4.times(T1))))));
	        	interParameterDerivatives[11]=d_parametersFromMAMB(dMA_horAxisErrPsi,dMB_horAxisErrPsi,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_horAxisErrPsi:");
	    			dMA_horAxisErrPsi.print(10, 5);
	    			System.out.println("dMB_horAxisErrPsi:");
	    			dMB_horAxisErrPsi.print(10, 5);
	    			System.out.println("interParameterDerivatives[11]="+sprintfArray(interParameterDerivatives[11]));
	    		}
	        } else interParameterDerivatives[11]=null; 
//	                                          // R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 || T0 | T1 | T2 | T3 |
	//12   	public double entrancePupilForward  //    |    |    |    |    |    |    |    || +  |    |    |    |
	//13   	public double centerAboveHorizontal //    |    |    |    |    |    |    |    ||    | +  |    |    |

	        if (mask[12]) {
	        	double [][] adT0_entrancePupilForward={{0.0},{0.0},{1.0}};
	        	Matrix dT0_entrancePupilForward=new Matrix(adT0_entrancePupilForward);
	        	Matrix dMA_entrancePupilForward=new Matrix(3,3,0.0);
	        	Matrix dMB_entrancePupilForward=MA.times(dT0_entrancePupilForward);
	        	interParameterDerivatives[12]=d_parametersFromMAMB(dMA_entrancePupilForward,dMB_entrancePupilForward,MA,MB,false); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_centerAboveHorizontal:");
	    			dMA_entrancePupilForward.print(10, 5);
	    			System.out.println("dMB_centerAboveHorizontal:");
	    			dMB_entrancePupilForward.print(10, 5);
	    			System.out.println("interParameterDerivatives[12]="+sprintfArray(interParameterDerivatives[12]));
	    		}
	        } else interParameterDerivatives[12]=null; 
	        
	        if (mask[13]) {
	        	double [][] adT1_centerAboveHorizontal={{0.0},{1.0},{0.0}};
	        	Matrix dT1_centerAboveHorizontal=new Matrix(adT1_centerAboveHorizontal);
	        	Matrix dMA_centerAboveHorizontal=new Matrix(3,3,0.0);
	        	Matrix dMB_centerAboveHorizontal=R8.times(R7.times(R6.times(R5.times(R4.times(dT1_centerAboveHorizontal)))));
	        	interParameterDerivatives[13]=d_parametersFromMAMB(dMA_centerAboveHorizontal,dMB_centerAboveHorizontal,MA,MB,false); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_centerAboveHorizontal:");
	    			dMA_centerAboveHorizontal.print(10, 5);
	    			System.out.println("dMB_centerAboveHorizontal:");
	    			dMB_centerAboveHorizontal.print(10, 5);
	    			System.out.println("interParameterDerivatives[13]="+sprintfArray(interParameterDerivatives[13]));
	    		}
	        } else interParameterDerivatives[13]=null; 
	        
	        
	//14  x	public double [] GXYZ=null; // coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system
	        if (mask[14]) {
	        	double [][] adT3_GXYZ0={{1.0},{0.0},{0.0}};
	        	Matrix dMA_GXYZ0=new Matrix(3,3,0.0);
	        	Matrix dMB_GXYZ0=new Matrix(adT3_GXYZ0);
	        	interParameterDerivatives[14]=d_parametersFromMAMB(dMA_GXYZ0,dMB_GXYZ0,MA,MB,false); // all after 6 are 0;
	        	if (this.debugLevel>2) {
	        		System.out.println("dMA_GXYZ0:");
	        		dMA_GXYZ0.print(10, 5);
	        		System.out.println("dMA_GXYZ0:");
	        		dMB_GXYZ0.print(10, 5);
	        		System.out.println("interParameterDerivatives[14]="+sprintfArray(interParameterDerivatives[14]));
	        	}
	        } else interParameterDerivatives[14]=null; 
	//13  y
	        if (mask[15]) {
//	    	double [][] adT3_GXYZ1={{0.0},{1.0},{0.0}};
	        	double [][] adT3_GXYZ1={{0.0},{-1.0},{0.0}}; // up - positive
//	        	double [][] adT3_GXYZ1={{0.0},{ 1.0},{0.0}}; // up - positive
	        	Matrix dMA_GXYZ1=new Matrix(3,3,0.0);
	        	Matrix dMB_GXYZ1=new Matrix(adT3_GXYZ1);
	        	interParameterDerivatives[15]=d_parametersFromMAMB(dMA_GXYZ1,dMB_GXYZ1,MA,MB,false); // all after 6 are 0;
	        	if (this.debugLevel>2) {
	        		System.out.println("dMA_GXYZ1:");
	        		dMA_GXYZ1.print(10, 5);
	        		System.out.println("dMB_GXYZ1:");
	        		dMB_GXYZ1.print(10, 5);
	        		System.out.println("interParameterDerivatives[15]="+sprintfArray(interParameterDerivatives[15]));
	        	}
	        } else interParameterDerivatives[15]=null; 
	//14  z
	        if (mask[16]) {
//	    	double [][] adT3_GXYZ2={{0.0},{0.0},{1.0}};
	        	double [][] adT3_GXYZ2={{0.0},{0.0},{-1.0}}; // away - positive
	        	Matrix dMA_GXYZ2=new Matrix(3,3,0.0);
	        	Matrix dMB_GXYZ2=new Matrix(adT3_GXYZ2);
	        	interParameterDerivatives[16]=d_parametersFromMAMB(dMA_GXYZ2,dMB_GXYZ2,MA,MB,false); // all after 6 are 0;
	        	if (this.debugLevel>2) {
	        		System.out.println("dMA_GXYZ2:");
	        		dMA_GXYZ2.print(10, 5);
	        		System.out.println("dMB_GXYZ2:");
	        		dMB_GXYZ2.print(10, 5);
	        		System.out.println("interParameterDerivatives[16]="+sprintfArray(interParameterDerivatives[16]));
	        	}
	        } else interParameterDerivatives[16]=null; 
	// now fill the rest, unchanged parameters
	/*
	    int numInputs=22;   //was 21  parameters in subcamera+...
	    int numOutputs=13;  //was 12  parameters in a single camera

	 */
	//17 (15)		public double focalLength=4.5;
	//18 (16)		public double px0=1296.0;          // center of the lens on the sensor, pixels
	//19 (17)		public double py0=968.0;           // center of the lens on the sensor, pixels
	//20 (18)		public double distortionA8=0.0; // r^8 (normalized to focal length or to sensor half width?)
	//21 (19)		public double distortionA7=0.0; // r^7 (normalized to focal length or to sensor half width?)
	//22 (20)		public double distortionA6=0.0; // r^6 (normalized to focal length or to sensor half width?)
	//23 (21)		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
	//24 (22)		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
	//25 (23)		public double distortionB=0.0; // r^3
	//26 (24)		public double distortionC=0.0; // r^2

//	        for (int inpPar=15; inpPar<getNumInputs(); inpPar++){
	        for (int inpPar=17; inpPar<getNumInputs(); inpPar++){ // OK with A8
	        	if (mask[inpPar]){
	        		interParameterDerivatives[inpPar]=new double[getNumOutputs()];
	        		for (int outPar=0; outPar<getNumOutputs(); outPar++){
	        			interParameterDerivatives[inpPar][outPar]=((getNumOutputs()-outPar)==(getNumInputs()-inpPar))?1.0:0.0;
	            		if (this.debugLevel>2) {
	            			System.out.println("interParameterDerivatives["+inpPar+"]["+outPar+"]="+interParameterDerivatives[inpPar][outPar]);
	            		}
	        		}
	        	} else interParameterDerivatives[inpPar]=null;
	        }
	    }

	    
	     public double [] parametersFromMAMB(Matrix MA, Matrix MB){
	    	double [] result= new double [getNumOutputs()]; // only first 6 are used
	    	for (int i=6; i<getNumOutputs();i++) result[i]=0.0;
	/*
	0		public double x0=0;      // lens axis from pattern center, mm (to the right)
	1		public double y0=0;      // lens axis from pattern center, mm (down)
	2		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
	3		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
	4		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
	5		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise

	        Vt=MB+MA*Vc
	        calculate X0,Y0 (Z0==0), dist, phi, theta, psi from MB and MA
	        | MA00 MA01 MA02 |
	        | MA10 MA11 MA12 |
	        | MA20 MA21 MA22 |
	        ===================
	        | MB0 |
	        | MB1 |
	        | MB2 |
	        Take point Pc on a sub-camera axis and on the target Zt=0 plane: [0,0,dist]
	        MA*Pc+MB= dist*[MA02,MA12,MA22]+MB =[Tx0,Ty0,0]
	        MA02*dist+MB0=Tx0
	        MA12*dist+MB1=Ty0
	        MA22*dist+MB2=0
	        dist=-MB2/MA22;
	        Tx0=MA02*(-MB2/MA22)+MB0
	        Ty0=MA12*(-MB2/MA22)+MB1
	        Tx0=MB0-MB2*MA02/MA22
	        Ty0=MB1-MB2*MA12/MA22
	        */
	    	result[2]=-MB.get(2,0)/MA.get(2,2); // distance
	    	result[0]= MB.get(0,0)-MB.get(2,0)*MA.get(0,2)/MA.get(2,2); // x0
	    	result[1]= MB.get(1,0)-MB.get(2,0)*MA.get(1,2)/MA.get(2,2); // y0
	        /*
	        // now find phi, theta, psi
	        MA - rotation from camera to target, transp(MA) - rotation from target to camera - same as rot(phi, theta,psi)
	        MT=transp(MA),
	        */
	        /*
	        MA[1,2]= sin(theta)
	        MA[0,2]= cos(theta)*sin(phi)
	        MA[2,2]= cos(theta)*cos(phi)
	        MA[1,0]=-cos(theta)*sin(psi)
	        MA[1,1]= cos(theta)*cos(psi)
	        theta=arcsin(MA[1,2])
	        phi=  atan2(MA[0,2],MA[2,2])
	        psi= -atan2(MA[1,0],MA[1,1])
	         */
	        result[4]=180.0/Math.PI* Math.asin(MA.get(1, 2));  //pitch
	        result[3]=  180.0/Math.PI* Math.atan2(MA.get(0,2),MA.get(2, 2));//yaw
	        result[5]=-180.0/Math.PI* Math.atan2(MA.get(1, 0),MA.get(1, 1));//roll
	    	return result;
	    }
	    /**
	     * 
	     * @param d_MA differential of the rotational matrix MA by some parameter
	     * @param d_MB differential of the translational matrix MB by some parameter
	     * @param MA - rotation matrix
	     * @param MB - translation matrix
	     * @param isAngle - when true, the partial derivative is for angles, d_MA, d_MB should be divided by 180/pi
	     * @return differentials of the {x0,y0,dist,phi,theta,psi} by that parameter
	     */
	    public double [] d_parametersFromMAMB(Matrix d_MA, Matrix d_MB, Matrix MA, Matrix MB, boolean isAngle){
	    	double [] result= new double [getNumOutputs()]; // only first 6 are used
	    	for (int i=6; i<getNumOutputs();i++) result[i]=0.0;
	/*
	0		public double x0=0;      // lens axis from pattern center, mm (to the right)
	1		public double y0=0;      // lens axis from pattern center, mm (down)
	2		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
	3		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
	4		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
	5		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise

	        Vt=MB+MA*Vc
	        calculate X0,Y0 (Z0==0), dist, phi, theta, psi from MB and MA
	        | MA00 MA01 MA02 |
	        | MA10 MA11 MA12 |
	        | MA20 MA21 MA22 |
	        ===================
	        | MB0 |
	        | MB1 |
	        | MB2 |
	        Take point Pc on a sub-camera axis and on the target Zt=0 plane: [0,0,dist]
	        MA*Pc+MB= dist*[MA02,MA12,MA22]+MB =[Tx0,Ty0,0]
	        MA02*dist+MB0=Tx0
	        MA12*dist+MB1=Ty0
	        MA22*dist+MB2=0
	        dist=-MB2/MA22;
	        Tx0=MB0-MB2*MA02/MA22
	        Ty0=MB1-MB2*MA12/MA22
	        */
//	    	result[2]=-MB.get(2,0)/MA.get(2,2); // distance
	/*
	        d_dist=-(d_MB2/MA22 - d_MA22*MB2/(MA22^2) ;
	 */
	    	double K=isAngle?(Math.PI/180):1.0;
	    	result[2]=K*(-d_MB.get(2,0)/MA.get(2,2)+
	    	d_MA.get(2,2)*MB.get(2,0)/(MA.get(2,2)*MA.get(2,2))); // d_distance
	    	/*
	        Tx0=MB0-MB2*MA02/MA22
	        d_Tx0=d_MB0 -  d_MB2*MA02/MA22 -d_MA02*MB2/MA22 +d_MA22*MB2*MA02/(MA22^2)
	    	 */
	    	result[0]= K*(d_MB.get(0,0) -
	    	d_MB.get(2,0)*MA.get(0,2)/MA.get(2,2) -
	    	d_MA.get(0,2)*MB.get(2,0)/MA.get(2,2) +
	    	d_MA.get(2,2)*MB.get(2,0)*MA.get(0,2)/(MA.get(2,2)*MA.get(2,2))); // x0

	    	/*
	        Ty0=MB0-MB2*MA02/MA22
	        d_Ty0=d_MB1 -  d_MB2*MA12/MA22 -d_MA12*MB2/MA22 +d_MA22*MB2*MA12/(MA22^2)
	    	 */
	    	result[1]= K*(d_MB.get(1,0) -
	    	d_MB.get(2,0)*MA.get(1,2)/MA.get(2,2) -
	    	d_MA.get(1,2)*MB.get(2,0)/MA.get(2,2) +
	    	d_MA.get(2,2)*MB.get(2,0)*MA.get(1,2)/(MA.get(2,2)*MA.get(2,2))); // y0
	        /*
	        // now find phi, theta, psi
	        MA - rotation from camera to target, transp(MA) - rotation from target to camera - same as rot(phi, theta,psi)
	        MT=transp(MA),
	        */
	        /*
	        MA[1,2]= sin(theta)
	        MA[0,2]= cos(theta)*sin(phi)
	        MA[2,2]= cos(theta)*cos(phi)
	        MA[1,0]=-cos(theta)*sin(psi)
	        MA[1,1]= cos(theta)*cos(psi)
	        theta=arcsin(MA[1,2])
	        phi=  atan2(MA[0,2],MA[2,2])
	        psi= -atan2(MA[1,0],MA[1,1])
	        arcsin(x)'= 1/sqrt(1-x^2)
	        arccos(x)'=-1/sqrt(1-x^2)
	        arctan(x)'= 1/sqrt(1+x^2)

	         */
	/*
	        theta=arcsin(MA12)
	        d_theta=d_MA12/sqrt(1-MA12^2)
	 */
	        result[4]=K*d_MA.get(1, 2)*180.0/Math.PI/Math.sqrt(1.0-MA.get(1, 2)*MA.get(1, 2));  //pitch
	/*
	        phi=  atan2(MA02,MA22)
	        d_phi=(d_MA02*MA22-d_MA22*MA02) / (MA22^2+MA02^2)
	 */

	        result[3]=  K*180.0/Math.PI*(d_MA.get(0,2)*MA.get(2,2) - d_MA.get(2,2)*MA.get(0,2)) /
	          ((MA.get(2,2)*MA.get(2,2)) + (MA.get(0,2)*MA.get(0,2)));//yaw
	/*
	        psi=  -atan2(MA10,MA11)
	        d_psi=-(d_MA10*MA11-d_MA11*MA10) / (MA10^2+MA11^2)
	*/
	        result[5]=  -K* 180.0/Math.PI*(d_MA.get(1,0)*MA.get(1,1) - d_MA.get(1,1)*MA.get(1,0)) /
	          ((MA.get(1,1)*MA.get(1,1)) + (MA.get(1,0)*MA.get(1,0)));//roll
	    	return result;
	    }

	    public String sprintfArray(double []arr){
	    	String result="";
	    	for (int i=0;i<arr.length;i++) result+= ((i>0)?", ":"")+arr[i];
	    	return "["+result+"]";
	    }
		
		public int getNumInputs(){return numInputs;}
		public int getNumOutputs(){return numOutputs;}
		
		
		
	}
 