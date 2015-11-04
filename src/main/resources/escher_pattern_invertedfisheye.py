#!/usr/bin/env python3
'''
# Copyright (C) 2015, Elphel.inc.
# Usage: known
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http:#www.gnu.org/licenses/>.
@author:     Oleg K Dzhimiev
@copyright:  2015 Elphel, Inc.
@license:    GPLv3.0+
@contact:    oleg@elphel.com
@deffield    updated: unknown
'''

__author__ = "Oleg K Dzhimiev"
__copyright__ = "Copyright 2015, Elphel, Inc."
__license__ = "GPL"
__version__ = "3.0+"
__maintainer__ = "Oleg K Dzhimiev"
__email__ = "oleg@elphel.com"
__status__ = "Silent Development"

import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(
    prog="escher_pattern_invertedfisheye.py",
    description="* To change the parameters of the lens & pattern edit the source file.\n\
* Python dependencies: numpy,scipy,pyx: \n\
 NumPy,SciPy - Ubuntu: sudo apt-get install python3-numpy python3-scipy\n\
 PyX - Download and install from: http://pyx.sourceforge.net/",
    usage="./escher_pattern_invertedfisheye.py or python3 escher_pattern_invertedfisheye.py",
    formatter_class=RawTextHelpFormatter
    )
args = parser.parse_args()


import numpy as np
from pyx import *
import math

import scipy as sc
from scipy.optimize import leastsq

#TODO:
# DONE 1. draw red arcs on cylinder pattern - don't break the filling - off/on
# DONE 2. approximate curves with more arcs thans just 1
# NOT CRITICAL 3. split floor pattern into 2 halves
# NOT CRITICAL 4. try creating horizontal/vertical lines
# DONE 5. Compare how extra points affect the calculations

#CONSTANTS

DEBUG = 0

#PRIMARY
LENS_FOV = 190
LENS_FOCAL_LENGTH = 2.0
PATTERN_CYLINDER_HEIGHT_M = 0.8

PATTERN_FLOOR_OVERLAP = 0.01 #1cm

PATTERN_SCALE = unit.t_m
#PATTERN_SCALE = unit.t_cm

#lens
LENS_POLYNOME_COEFFICIENTS = [LENS_FOCAL_LENGTH,0,0,0,0]
#pattern (object plane side)
#PATTERN_CYLINDER_RADIUS_M = 1.5
#Black spot size on the sensor in pixels
PATTERN_PERIOD_IMG_PLANE_NPX = 50
PATTERN_E = 2 #some constant - affects the curves
PATTERN_ANGLE = 0 #departure angle from vertical

#production values
PATTERN_NUMBER_OF_POINTS_FOR_CURVE_APPROXIMATION = 50
PATTERN_NUMBER_OF_EXTRA_END_POINTS_APPROXIMATION = 100
PATTERN_EPSILON_FOR_LEAST_SQUARE_APPROXIMATION = 0.00001 #in scale parts

#debug values
#PATTERN_NUMBER_OF_POINTS_FOR_CURVE_APPROXIMATION = 20
#PATTERN_NUMBER_OF_EXTRA_END_POINTS_APPROXIMATION = 10
#PATTERN_EPSILON_FOR_LEAST_SQUARE_APPROXIMATION = 0.0001 #in scale parts

#LEAST IMPORTANT OR NOT USED
#sensor (image plane)
SENSOR_WIDTH_MM = 5.70
SENSOR_HEIGHT_MM = 4.28
SENSOR_PXSIZE_UM = 2.2

#CLASSES
class Lens:
    def __init__(self,k=[2,0,0,0,0],fov=180):
        try: 
            #let the first coefficient be focal length
            self.f = k[0]
            self.k = k
            #field of view
            self.fov = fov
        except IndexError:
            print("ERROR: Lens init")
            
    #forward translation
    #theta - angle between principal axis and incoming ray, theta_max = fov/2
    def r_mm (self,theta_rad):
        result = 0
        for i in range(len(self.k)):
            result += self.k[i]*(theta_rad**(2*i+1))
        return result
    
    #backward translation
    def theta_rad (self,r_mm,tolerance_rad=0.0001):
        theta1 = 0
        theta2 = 2*math.pi
        while (abs(theta2-theta1)>tolerance_rad):
            r1  = self.r_mm(theta1)-r_mm
            r2  = self.r_mm(theta2)-r_mm
            theta  = theta1 - (theta2-theta1)*r1/(r2-r1)
            theta1 = theta2
            theta2 = theta
        return theta
        
    def theta_deg (self,r_mm,tolerance_rad=0.0001):
        return (self.theta_rad(r_mm,tolerance_rad)*180/math.pi)
        
class Sensor:
    def __init__(self,w_mm=5.70,h_mm=4.28,px_um=2.2):
        self.w_mm=w_mm
        self.h_mm=h_mm
        self.d_mm = math.sqrt(w_mm**2+h_mm*2)
        self.px_um=px_um
        
class Pattern:
    def __init__(self,H_m,lens_fov,period_n=40,e=2,angle=0,N=50,N_extra=10,epsilon=0.0001,scale=10,floor_overlap=0.01):
        self.H_m = H_m
        self.lens_fov = lens_fov
        tg = math.tan(lens_fov/2*math.pi/180)
        self.R_m = H_m*tg/(tg-1)
        self.period_n = period_n
        #pattern parameters
        self.angle = angle
        self.e = e
        self.size = 1
        self.hsize = self.size/2
        self.qsize = self.size/4
        a = e*(math.sqrt(2)-1)
        self.a = a
        self.r = (a*a+1)/(2*a)*self.qsize
        self.h = math.sqrt(self.r*self.r-self.qsize*self.qsize)
        self.dc = 2*self.qsize-self.h
        self.halfangle = 180/math.pi*math.atan(self.qsize/self.h)
        self.N = N
        self.N_extra = N_extra
        self.epsilon = epsilon
        self.scale = scale
        self.floor_overlap = floor_overlap

# init objects
L = Lens(
            LENS_POLYNOME_COEFFICIENTS,
            LENS_FOV
        )

S = Sensor(
            SENSOR_WIDTH_MM,
            SENSOR_HEIGHT_MM,
            SENSOR_PXSIZE_UM
            )

P = Pattern(
            PATTERN_CYLINDER_HEIGHT_M,
            L.fov,
            PATTERN_PERIOD_IMG_PLANE_NPX,PATTERN_E,
            PATTERN_ANGLE,
            PATTERN_NUMBER_OF_POINTS_FOR_CURVE_APPROXIMATION,
            PATTERN_NUMBER_OF_EXTRA_END_POINTS_APPROXIMATION,
            PATTERN_EPSILON_FOR_LEAST_SQUARE_APPROXIMATION,
            PATTERN_SCALE,
            PATTERN_FLOOR_OVERLAP
            )

print("CHECK 0: The pattern scale units are:")
print(P.scale)

#FUNCTIONS

def CircleCenterBy2PointsAndRadius(x1,y1,x2,y2,R):
    k = (y2-y1)/(x2-x1)
    k = -1/k
    xm = (x1+x2)/2
    ym = (y1+y2)/2
    
    b = ym - k*xm
    
    A = 1+(k**2)
    B = -2*x1-2*k*(y1-b)
    C = (x1**2)+(y1-b)**2-(R**2)
    
    D = B*B-4*A*C
    xc1 = (-B+math.sqrt(D))/(2*A)
    xc2 = (-B-math.sqrt(D))/(2*A)
    
    yc1 = k*xc1 + b
    yc2 = k*xc2 + b
    
    return [[xc1,yc1],[xc2,yc2]]

def MoveCircle(x1,y1,x2,y2,xc,yc,R):
    #move circle so the end points would lie to it
    points = CircleCenterBy2PointsAndRadius(x1,y1,x2,y2,R)
    d1 = math.sqrt((points[0][0]-xc)*(points[0][0]-xc)+(points[0][1]-yc)*(points[0][1]-yc))
    d2 = math.sqrt((points[1][0]-xc)*(points[1][0]-xc)+(points[1][1]-yc)*(points[1][1]-yc))
    if (d1<d2):
        return points[0]
    else:
        return points[1]

def ArcsFit(points,spl):
    FunCircle = lambda p,x: p[1]-np.sqrt((p[2])*(p[2])-(x-p[0])*(x-p[0]))
    Fun = lambda p,x: -p[0]/p[1]*x-p[2]/p[1]
    ErrorFun = lambda p,x,y: (y*y+x*x)-(p[0]*x+p[1]*y+p[2])
    
    result = []
    
    xdata = np.array([])
    ydata = np.array([])
    for point in points:
        xdata = np.append(xdata,point[0])
        ydata = np.append(ydata,point[1])
    
    p0 = (1,1,1)
    x0w = np.array([])
    y0w = np.array([])
    xnw = np.array([])
    ynw = np.array([])
    #add weights like 1000x replicate end points
    for i in range(P.N_extra):
        x0w = np.append(x0w,xdata[0])
        y0w = np.append(y0w,ydata[0])
        xnw = np.append(xnw,xdata[len(xdata)-1])
        ynw = np.append(ynw,ydata[len(xdata)-1])
    
    xdata_leastsq = np.concatenate((x0w,xdata,xnw),axis=0)
    ydata_leastsq = np.concatenate((y0w,ydata,ynw),axis=0)
    if (DEBUG==1): print("LEASTSQ RESULTS:")
    abc,cov_x = leastsq(ErrorFun,p0,args=(xdata,ydata))
    x0 = abc[0]/2
    y0 = abc[1]/2
    if (DEBUG==1): print(x0,y0,math.sqrt(abc[2]+x0*x0+y0*y0))
    abc,cov_x = leastsq(ErrorFun,p0,args=(xdata_leastsq,ydata_leastsq))
    x0 = abc[0]/2
    y0 = abc[1]/2
    if (DEBUG==1): print(x0,y0,math.sqrt(abc[2]+x0*x0+y0*y0))
    x0 = abc[0]/2
    y0 = abc[1]/2
    r = math.sqrt(abc[2]+x0*x0+y0*y0)

    #place circle through border points
    # sometimes the distance between points is bigger than 2*r and this is ok
    if (r>=(math.sqrt((xdata[0]-xdata[len(xdata)-1])**2+(ydata[0]-ydata[len(xdata)-1])**2))/2):
        x0,y0 = MoveCircle(xdata[0],ydata[0],xdata[len(xdata)-1],ydata[len(xdata)-1],x0,y0,r)

    if (DEBUG==1): print(x0,y0,r)
    s = np.sqrt(np.sum((np.sqrt((xdata-x0)**2+(ydata-y0)**2)-r)**2))
    if (DEBUG==1): print(s)
    
    a0 = 180/math.pi*math.atan2(xdata[0]-x0,ydata[0]-y0)-90
    an = 180/math.pi*math.atan2(xdata[len(xdata)-1]-x0,ydata[len(xdata)-1]-y0)-90
    
    if (s>P.epsilon):
        #print("Split=%d RMS=%f"%(spl,s))
        #cut in half
        if (len(points)>=6):
            newpoints = np.array_split(points, 2)
            #arrays should overlap
            newpoints[0] = np.append(newpoints[0],[newpoints[1][0]],axis=0)
            tmp1 = ArcsFit(newpoints[0],spl+1)
            tmp2 = ArcsFit(newpoints[1],spl+1)
            result.extend(tmp1)
            result.extend(tmp2)
        else:
            result.extend([[x0,y0,r,-a0,-an]])
    else:
        result.extend([[x0,y0,r,-a0,-an]])
    
    return result

def cylinder_curve_points(x,y,xc,yc,a0,an,n,delta):
    global P,L,S
    adata = np.linspace(a0,an,n)
    xdata = (xc)+P.r*np.cos(adata*math.pi/180)
    ydata = (yc)+P.r*np.sin(adata*math.pi/180)
    
    phidata = np.arctan2(ydata,xdata)
    #change pi to -pi for atan2
    boundary_case1 = (phidata[0]>0)and(phidata[len(adata)-1]<0)
    boundary_case2 = (phidata[0]<0)and(phidata[len(adata)-1]>0)
    if (boundary_case1): phidata[0]=-phidata[0]
    if (boundary_case2): phidata[len(adata)-1]=-phidata[len(adata)-1]
       
    rdata = np.sqrt(xdata**2+ydata**2)*img_pp_mm
    tdata = np.fromiter((L.theta_rad(d) for d in rdata),dtype=float)
        
    udata = []
    for i in range(len(adata)):
        if (x<0 and y==0):
            if (phidata[i]<0): 
                udata.append(P.R_m*(-(phidata[i]+math.pi)+delta))
            else:           
                udata.append(P.R_m*(2*math.pi-(phidata[i]+math.pi)+delta))
        else:
            udata.append(P.R_m*(2*math.pi-(phidata[i]+math.pi)))
        
    vdata = P.H_m+P.R_m*np.tan(tdata-math.pi/2)
    
    pth = []
    for j in range(len(udata)):
        pth.append([udata[j],vdata[j]])
    return pth
            

def floor_curve_points(xc,yc,a0,an,n):
    global P,L,S
    adata = np.linspace(a0,an,n)
    xdata = (xc)+P.r*np.cos(adata*math.pi/180)
    ydata = (yc)+P.r*np.sin(adata*math.pi/180)
    rdata = np.sqrt(xdata**2+ydata**2)*img_pp_mm
    cosdata = xdata/rdata*img_pp_mm
    sindata = ydata/rdata*img_pp_mm
    tdata = np.fromiter((L.theta_rad(d) for d in rdata),dtype=float)
    R_m = P.H_m*np.tan(tdata)
    udata = R_m*cosdata
    vdata = R_m*sindata
    pth = []
    for j in range(len(udata)):
        pth.append([udata[j],vdata[j]])
    return pth

#MAIN

# Array for pattern check arcs - shifts and angles
D = [
    [-P.size+P.dc,-P.qsize    ,  0-P.halfangle,   0+P.halfangle],
    [-P.dc       , P.qsize    ,180+P.halfangle, 180-P.halfangle],
    [-P.qsize    , P.size-P.dc,270-P.halfangle, 270+P.halfangle],
    [ P.qsize    , P.dc       , 90+P.halfangle,  90-P.halfangle],
    [ P.size-P.dc, P.qsize    ,180-P.halfangle, 180+P.halfangle],
    [ P.dc       ,-P.qsize    ,  0+P.halfangle,   0-P.halfangle],
    [ P.qsize    ,-P.size+P.dc, 90-P.halfangle,  90+P.halfangle],
    [-P.qsize    ,-P.dc       ,270+P.halfangle, 270-P.halfangle]
]

# find max theta for the sensor
img_rmin_mm = S.h_mm
img_rmax_mm = S.d_mm
#print("CHECK 1: Min/Max radii on image plane: r_min = %0.2fmm, r_max = %0.2fmm"%(img_rmin_mm,img_rmax_mm))

img_pattern_period_mm = S.px_um*P.period_n/1000
img_pp_mm = img_pattern_period_mm

print("CHECK 2: Pattern cell size on image plane: %fmm"%img_pattern_period_mm)

# next... draw an ideal image plane filled with straight pattern - it's not needed - so, jus tt ocheck

#L.r_mm(L.fov/2*math.pi/180)

img_w_n = 2*L.r_mm(L.fov/2*math.pi/180)/img_pattern_period_mm
img_h_n = 2*L.r_mm(L.fov/2*math.pi/180)/img_pattern_period_mm
#thinking circles
img_r = L.r_mm(L.fov/2*math.pi/180)

print("CHECK 3: Image plane size in periods is: %d x %d"%(img_w_n,img_h_n))

print("CHECK 4: Floor Pattern Radius is: %f m"%(P.R_m))

print("CHECK 5: Cylinder Pattern Radius= %f m, Height= %f m"%(P.R_m,P.H_m))

print("NOTE: Do not pay attention to 'content exceeds the papersize' - it's ok")

# ip - Image Plane - the image of the pattern on the "image plane"
ip = canvas.canvas()

# opf - Object Plane Floor - pattern wallpaper that will go to the floor
# Draw within a circle extended by "floor_overlap" for overlapping.
opf_clippath = path.circle(0,0,(P.R_m+P.floor_overlap)*P.scale)
opf = canvas.canvas([canvas.clip(opf_clippath)])

# opc - Object Plane Cylinder - pattern wallpaper that will go to the cylinder wall
# Draw within a rectangle 2piR x H
opc_clippath = path.rect(0*P.scale,0*P.scale,2*math.pi*P.R_m*P.scale,P.H_m*P.scale)
opc = canvas.canvas([canvas.clip(opc_clippath)])

#add more periods to ends to compensate for tilt. do we need tilt at all? set tilt to 0 degrees for now
wn = img_w_n+img_h_n*math.sin(P.angle*math.pi/180)
hn = img_h_n+img_w_n*math.sin(P.angle*math.pi/180)

# the (optical) center is 0,0
for y in range(int(-hn/2),int(hn/2),2):
    for x in range(int(-wn/2),int(wn/2),2):
        # draw 2 black spots at once - top-left + bottom-right
        for i in range(2):
            tmp_x = x + i
            tmp_y = y + i
            
            # IMAGE PLANE BEGIN
            # nothing to debug
            # tmp_x & tmp_y - center of the spot, the spot is drawn within center +/- (r+dc) radius
            # 8 paths, scaling to mm, draw a single black shape
            # a tile of image plane
            p = path.path()
            for j in range(4):
                # D[][] is the array that contains shifts and angles for the base arcs
                p += path.path(path.arc ((tmp_x+D[2*j+0][0])*img_pp_mm,(tmp_y+D[2*j+0][1])*img_pp_mm, P.r*img_pp_mm,D[2*j+0][2],D[2*j+0][3]))
                p += path.path(path.arcn((tmp_x+D[2*j+1][0])*img_pp_mm,(tmp_y+D[2*j+1][1])*img_pp_mm, P.r*img_pp_mm,D[2*j+1][2],D[2*j+1][3]))
             
            #### drawing ###
            ip.stroke(p,[style.linewidth(0.000*img_pp_mm),deco.filled([color.rgb.black]),trafo.rotate(-P.angle)])
            # IMAGE PLANE END
            
            # OBJECT PLANE: FLOOR & CYLINDER
            # spot center distance from the pattern center @(0,0)
            tmp_r = math.sqrt(tmp_x**2+tmp_y**2)*img_pp_mm
            # img_r is radius for the cylinder - anything closer is on the floor
            # img_r is calculated long ago.
            if (tmp_r<img_r):
                # tmp_r can be a zero
                if (DEBUG==1):
                    if (tmp_r==0): print("tmp_r is 0")
                
                # if tmp_r==0 - cos(phi) and sin(phi) do not matter - so, some previous ones can be safely taken
                if (tmp_r>0): tmp_cos_phi = tmp_x/tmp_r*img_pp_mm
                if (tmp_r>0): tmp_sin_phi = tmp_y/tmp_r*img_pp_mm
                
                # max values are not used
                # (P.r+P.dc) is pattern specific
                tmp_r_min = tmp_r - (P.r+P.dc)*img_pp_mm
                if (tmp_r_min<0): tmp_r_min = 0
                tmp_r_max = tmp_r + (P.r+P.dc)*img_pp_mm
                
                # theta - angle between principal axis and incoming ray
                tmp_theta_min = L.theta_rad(tmp_r_min)
                tmp_theta_max = L.theta_rad(tmp_r_max)
                tmp_theta_c   = L.theta_rad(tmp_r)
                    
                #tangent can go negative
                R_min_m = P.H_m*math.tan(tmp_theta_min)
                R_max_m = P.H_m*math.tan(tmp_theta_max)
                
                # for floor
                tmp_u_min = R_min_m*tmp_cos_phi
                tmp_u_max = R_max_m*tmp_cos_phi
                # for floor
                tmp_v_min = R_min_m*tmp_sin_phi
                tmp_v_max = R_max_m*tmp_sin_phi
                
                # not used anymore?
                R_c_m = P.H_m*math.tan(tmp_theta_c)
                tmp_u_c = R_c_m*tmp_cos_phi
                tmp_v_c = R_c_m*tmp_sin_phi
                            
                #### FLOOR PATTERN CHECK ###
                # ATTENTION: 3 chunks of code are basically the same - only some of the called functions are different
                if ((R_min_m>=0)and(abs(tmp_u_min)<=P.R_m)and(abs(tmp_v_min)<=P.R_m)):
                    # reinitialize
                    pth_arcs,pth_lines = [],[]
                    # using 'arc' and 'arcn' - 8/2 = 4
                    for j in range(4):
                        # step 1 - using arc
                        pth_points_array = floor_curve_points(tmp_x+D[2*j+0][0],tmp_y+D[2*j+0][1],D[2*j+0][2],D[2*j+0][3],P.N)
                        pth_lines.extend(pth_points_array)
                        
                        pth_arcs_array = ArcsFit(pth_points_array,0)
                        for pth_arc_single in pth_arcs_array:
                            pth_arc = path.path(path.arc(pth_arc_single[0]*P.scale,pth_arc_single[1]*P.scale, pth_arc_single[2]*P.scale, pth_arc_single[3], pth_arc_single[4]))
                            pth_arcs.append(pth_arc)
                            
                        # step 2 - using arcn (the only difference)
                        pth_points_array = floor_curve_points(tmp_x+D[2*j+1][0],tmp_y+D[2*j+1][1],D[2*j+1][2],D[2*j+1][3],P.N)
                        pth_lines.extend(pth_points_array)
                        
                        pth_arcs_array = ArcsFit(pth_points_array,0)
                        for pth_arc_single in pth_arcs_array:
                            pth_arc = path.path(path.arcn(pth_arc_single[0]*P.scale,pth_arc_single[1]*P.scale, pth_arc_single[2]*P.scale, pth_arc_single[3], pth_arc_single[4]))
                            pth_arcs.append(pth_arc)
                    
                    p_sum_arcs = pth_arcs[0]
                    for k in range(1,len(pth_arcs)): p_sum_arcs = p_sum_arcs + pth_arcs[k]
                    p_sum_arcs.append(path.closepath())
                    
                    #### drawing ###
                    opf.stroke(p_sum_arcs,[style.linewidth(0*P.scale),deco.filled([color.rgb.black])])
                    
                    if (DEBUG==1):
                        p_sum_lines = path.path(path.moveto(pth_lines[0][0]*P.scale,pth_lines[0][1]*P.scale))
                        for k in range(1,len(pth_lines)): p_sum_lines.append(path.lineto(pth_lines[k][0]*P.scale,pth_lines[k][1]*P.scale))
                        p_sum_lines.append(path.closepath())
                        ####
                        opf.stroke(p_sum_lines,[style.linewidth(0.000*P.scale),color.rgb.red])
                        ####
                
                #### CYLINDER WALL PATTERN CHECK ###
                if (tmp_theta_max>=math.atan(P.R_m/P.H_m)):
                    #if the center of the spot is on the cylinder. 
                    #take those 8 curves
                    delta = 0
                    # reinitialize
                    pth_arcs, pth_lines = [],[]
                    for j in range(4):
                        # step 1 - using arc
                        pth_points_array = cylinder_curve_points(x,y,tmp_x+D[2*j+0][0],tmp_y+D[2*j+0][1],D[2*j+0][2],D[2*j+0][3],P.N,delta)
                        pth_lines.extend(pth_points_array)
                        
                        pth_arcs_array = ArcsFit(pth_points_array,0)
                        for pth_arc_single in pth_arcs_array:
                            pth_arc = path.path(path.arc(pth_arc_single[0]*P.scale,pth_arc_single[1]*P.scale, pth_arc_single[2]*P.scale, pth_arc_single[3], pth_arc_single[4]))
                            pth_arcs.append(pth_arc)
                            
                        # step 2 - using arcn (the only difference)
                        pth_points_array = cylinder_curve_points(x,y,tmp_x+D[2*j+1][0],tmp_y+D[2*j+1][1],D[2*j+1][2],D[2*j+1][3],P.N,delta)
                        pth_lines.extend(pth_points_array)
                        
                        pth_arcs_array = ArcsFit(pth_points_array,0)
                        for pth_arc_single in pth_arcs_array:
                            pth_arc = path.path(path.arcn(pth_arc_single[0]*P.scale,pth_arc_single[1]*P.scale, pth_arc_single[2]*P.scale, pth_arc_single[3], pth_arc_single[4]))
                            pth_arcs.append(pth_arc)
                        
                    p_sum_arcs = pth_arcs[0]
                    for k in range(1,len(pth_arcs)): p_sum_arcs = p_sum_arcs + pth_arcs[k]
                    p_sum_arcs.append(path.closepath())
                    
                    #### drawing ###
                    opc.stroke(p_sum_arcs,[style.linewidth(0*P.scale),deco.filled([color.rgb.black])])
                        
                    if (DEBUG==1):
                        p_sum = path.path(path.moveto(pth_lines[0][0]*P.scale,pth_lines[0][1]*P.scale))
                        for k in range(1,len(pth_lines)): p_sum.append(path.lineto(pth_lines[k][0]*P.scale,pth_lines[k][1]*P.scale))    
                        ####
                        opc.stroke(p_sum,[style.linewidth(0.000*P.scale),color.rgb.red])
                        ####
                        
                    #add extra point for 2piR
                    if (x<0 and y==0):
                        delta = 2*math.pi
                        # reinitialize
                        pth_arcs, pth_lines = [],[]
                        for j in range(4):
                            # step 1 - using arc
                            pth_points_array = cylinder_curve_points(x,y,tmp_x+D[2*j+0][0],tmp_y+D[2*j+0][1],D[2*j+0][2],D[2*j+0][3],P.N,delta)
                            pth_lines.extend(pth_points_array)
                            
                            pth_arcs_array = ArcsFit(pth_points_array,0)
                            for pth_arc_single in pth_arcs_array:
                                pth_arc = path.path(path.arc(pth_arc_single[0]*P.scale,pth_arc_single[1]*P.scale, pth_arc_single[2]*P.scale, pth_arc_single[3], pth_arc_single[4]))
                                pth_arcs.append(pth_arc)
                                
                            # step 2 - using arcn (the only difference)
                            pth_points_array = cylinder_curve_points(x,y,tmp_x+D[2*j+1][0],tmp_y+D[2*j+1][1],D[2*j+1][2],D[2*j+1][3],P.N,delta)
                            pth_lines.extend(pth_points_array)
                            
                            pth_arcs_array = ArcsFit(pth_points_array,0)
                            for pth_arc_single in pth_arcs_array:
                                pth_arc = path.path(path.arcn(pth_arc_single[0]*P.scale,pth_arc_single[1]*P.scale, pth_arc_single[2]*P.scale, pth_arc_single[3], pth_arc_single[4]))
                                pth_arcs.append(pth_arc)
                           
                        p_sum_arcs = pth_arcs[0]
                        for k in range(1,len(pth_arcs)): p_sum_arcs = p_sum_arcs + pth_arcs[k]
                        p_sum_arcs.append(path.closepath())
                        
                        #### drawing ###
                        opc.stroke(p_sum_arcs,[style.linewidth(0*P.scale),deco.filled([color.rgb.black])])   
                        
                        if (DEBUG==1):
                            p_sum = path.path(path.moveto(pth_lines[0][0]*P.scale,pth_lines[0][1]*P.scale))
                            for k in range(1,len(pth_lines)): p_sum.append(path.lineto(pth_lines[k][0]*P.scale,pth_lines[k][1]*P.scale))    
                            ####
                            opc.stroke(p_sum,[style.linewidth(0.000*P.scale),color.rgb.red])
                            ####

# option #1 for creating pdf
#ip.writePDFfile("image_plane.pdf")
# option #2
d_ip = document.document()
# centimeter scale for image/sensor plane
d_ip.append(document.page(ip,paperformat=document.paperformat(S.w_mm*unit.t_cm,S.h_mm*unit.t_cm)))
d_ip.writePDFfile("image_plane.pdf",
                 title="\"Escher\" Pattern, Image Plane",
                 subject="Lens calibration pattern",
                 keywords="camera, lens, pattern, checker board, Escher",
                 author="Fabian Gottlieb von Bellingshausen")

d_opf = document.document()
d_opf.append(document.page(opf,paperformat=document.paperformat(2*(P.R_m+P.floor_overlap)*P.scale,2*(P.R_m+P.floor_overlap)*P.scale)))
d_opf.writePDFfile("object_plane_floor.pdf",
                 title="\"Escher\" Pattern, Object Plane (Floor)",
                 subject="Lens calibration pattern",
                 keywords="camera, lens, pattern, checker board, Escher",
                 author="Fabian Gottlieb von Bellingshausen")

d_opc = document.document()
d_opc.append(document.page(opc,paperformat=document.paperformat(2*math.pi*P.R_m*P.scale,P.H_m*P.scale)))
d_opc.writePDFfile("object_plane_cylinder.pdf",
                 title="\"Escher\" Pattern, Object Plane (Cylinder)",
                 subject="Lens calibration pattern",
                 keywords="camera, lens, pattern, checker board, Escher",
                 author="Fabian Gottlieb von Bellingshausen")



# used to draw just a circle - not anymore 
# ip.stroke(path.circle(tmp_x*img_pp_mm,tmp_y*img_pp_mm,0.5*img_pp_mm),[style.linewidth.Thin,deco.filled([color.rgb.black])])
