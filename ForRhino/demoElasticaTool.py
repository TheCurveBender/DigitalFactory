import Rhino
import math
import rhinoscriptsyntax as rs
from Rhino.Geometry import *
import Rhino.Geometry.Point3d as Point3D
import scriptcontext
import sys
from System import Array
from System import Double
import clr
clr.AddReferenceToFileAndPath("DyLibs/ElasticaForRhino_GITHUBV2")
from ElasticaForRhino_GITHUBV2 import Elastica
import cPickle as pickle


### Demostration of elastic tools features
###

#Constants
N_segments = 100
N_samples = 1000
e = Elastica()

#Function
def add_curve_to_rhino(P, N):
    pts_out = []
    for i in range(0,N*2,2):
        pt = Point3d(P[i],P[i+1],0)
        pts_out.append(pt)
    pl  = Rhino.Geometry.Polyline(pts_out)
    plc = Rhino.Geometry.PolylineCurve(pl)
    scriptcontext.doc.Objects.AddCurve(plc)

####Feature 1: discrete elastica solver (discrete_elastica)
###

curve = rs.GetObject("Choose curve to get boundary conditions",4,True)
c_len = rs.CurveLength(curve)
p_s   = rs.CurveStartPoint(curve)
p_e   = rs.CurveEndPoint(curve)
dom   = rs.CurveDomain(curve)
t_s   = rs.CurveTangent(curve,dom[0])
t_e   = rs.CurveTangent(curve,dom[1])

# Compute the discrete elastic curve
P_discrete_elastica = Array.CreateInstance(Double,(N_segments+1)*2)
iter = e.discrete_elastica(p_s.X,p_s.Y, p_e.X,p_e.Y, t_s.X,t_s.Y, t_e.X,t_e.Y, c_len, N_segments, P_discrete_elastica)
add_curve_to_rhino(P_discrete_elastica, N_segments+1)
print("iterations: ", iter) # Print number of iterations needed to compute the discrete elastic curve


### Feature 2: evaluate analytic elastic curve (sample_elastica)
###

#Sample the analytic elastic curve
analytic_elastica_params = Array[float]([0.65, 0, 20,2, 1.4, 1, 1,1]) #k, s0, l, scale, angle, x0,y0, infl/ninfl
P_samples = Array.CreateInstance(Double, 2*N_samples)
e.sample_elastica(analytic_elastica_params, -0.15, 1.15, N_samples, P_samples)
add_curve_to_rhino(P_samples, N_samples)

#### Feature 3: design of numerical elastic splines
####


## Data
with open('objs.pickle') as f:  # Python 3: open(..., 'rb')
    AngError, AngErrorG12, G11x1m,G11y1m,G11x2m,G11y2m, G12x1m,G12y1m,G12x2m,G12y2m = pickle.load(f)
AngErrorV = Array[int](AngError);
G11x1mV = Array[float](G11x1m);
G11y1mV = Array[float](G11y1m);
G11x2mV = Array[float](G11x2m);
G11y2mV = Array[float](G11y2m);
AngErrorG12V = Array[int](AngErrorG12);
G12x1mV = Array[float](G12x1m);
G12y1mV = Array[float](G12y1m);
G12x2mV = Array[float](G12x2m);
G12y2mV = Array[float](G12y2m);


## Initial curve with four elastic curve segments

kn = [0,0,0,1/4,1/4,1/2,1/2,3/4,3/4,1,1,1]
points = [[-6.0000,0.0001,0],[-5.0000,0,0],[   -4.0000,0.0001 ,0],[-2.0000,0,0],[-1.000,0,0],[1.0000,0.0001,0],[2.000,0,0],[4.0000,0.0001,0],[5.00,0,0],[6,0.0001,0]]
crv = rs.AddNurbsCurve(points,kn,3);
N = rs.CurvePointCount(crv);

rs.EnableObjectGrips(crv)
actvpt = rs.GetObjectGrip("Choose active point")
actvpt0 = actvpt[0]
actvpt = actvpt[1]+1;

# Interactive part
while True:
    Cpxvec = [];
    Cpyvec = []
    for i in range(0,N):
        Cpxvec.append(points[i][0])
        Cpyvec.append(points[i][1])
    Cpx = Array[float](Cpxvec);
    Cpy = Array[float](Cpyvec);

    #Determine good locations for the active point
    GoodLocation = Array.CreateInstance(Double,1000000)
    val = e.GoodPointLocations(Cpx,Cpy,N,actvpt,AngErrorV,G11x1mV,G11y1mV,G11x2mV,G11y2mV,AngErrorG12V,G12x1mV,G12y1mV,G12x2mV,G12y2mV,GoodLocation)

    ptcld = [];
    ptcld2 = [];

    if GoodLocation[0] == 0: #Only one curve is affected
        xvals = GoodLocation[2:int(GoodLocation[1])+2]
        yvals = GoodLocation[int(GoodLocation[1])+2:int(GoodLocation[1])*2+2]
        for j in range(0,int(GoodLocation[1])):
            ptcld.append([xvals[j],yvals[j],0])
        ptcld = rs.AddPointCloud(ptcld)
        rs.ObjectColor(ptcld,color=(255,0,255))

    if GoodLocation[0] == 1: #Two curves are affected
        valpc2 = int((val-2-2*GoodLocation[1])/2.)
        xvals1 = GoodLocation[2:int(GoodLocation[1])+2]
        xvals2 = GoodLocation[int(GoodLocation[1])+2:int(GoodLocation[1])+2+valpc2]
        yvals1 = GoodLocation[int(GoodLocation[1])+2+valpc2:int(GoodLocation[1])+2+valpc2+int(GoodLocation[1])]
        yvals2 = GoodLocation[int(GoodLocation[1])+2+valpc2+int(GoodLocation[1]):valpc2+int(GoodLocation[1])+2+valpc2+int(GoodLocation[1])]
        for j in range(0,int(GoodLocation[1])):
            ptcld.append([xvals1[j],yvals1[j],0])
        ptcld = rs.AddPointCloud(ptcld)
        rs.ObjectColor(ptcld,color=(255,0,255))
        for j in range(0,valpc2):
            ptcld2.append([xvals2[j],yvals2[j],0])
        ptcld2 = rs.AddPointCloud(ptcld2)
        rs.ObjectColor(ptcld2,color=(0,255,255))

    #drag active point
    T = rs.GetPoint("Move point")
    rs.DeleteObject(crv);

    points[actvpt-1] = T;
    crv = rs.AddNurbsCurve(points,kn,3);
    rs.EnableObjectGrips(crv)
    actvpt = rs.GetObjectGrip("Choose active point")
    if actvpt is None:
        break
    actvpt0 = actvpt[0]
    actvpt = actvpt[1]+1;
    rs.DeleteObject(ptcld)
    if GoodLocation[0] == 1:
        rs.DeleteObject(ptcld2)
    if not actvpt:
        break

