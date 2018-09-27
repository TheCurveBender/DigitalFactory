# A library for designing with elastic curves (Euler's Elastica)

This GITHUB repo contrains algorithms for designing with an elastic curve, elastic splines, and surfaces foliated by elastic splines (https://en.wikipedia.org/wiki/Elastica_theory).

Because elastic curves describe the shape of a flexible blade, the algorithms can be employed by architects/designers to create form work made with robotic hot-blade cutting (https://odico.dk/en/technologies#robotic-hot-blade-cutting) or wooden strip designs.  

## Getting started

In the folder "For Rhino" we include "demoElasticaTools.py" that contains demonstration of the features from the dll files. In the folder "MATLAB prototypes" we include prototypes of features developed in MATLAB.

### Prerequisites

The code in the folder "For Rhino" was developed for Rhinoceros 5 64 bit. The MATLAB code was developed using MATLAB R2016b.

### Installing

To test the features for Rhino: download the folder "For Rhino", open Rhinoceros and the Python editor (command: "EditPythonScript"). Once the Python editor is open, open the file "demoElasticaTools.py" and run the script (click on the green arrow).

To test the MATLAB prototypes download the MATLAB files in the folder for the prototypes and run the prototype directly from MATLAB.

## Features and demostrations

### Rhinoceros

In the folder "For Rhino" we include the demoElasticaTools.py that contains demonstration of the features from the dll files.
At the top of the file we write:

```
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
clr.AddReferenceToFileAndPath("DyLibs/ElasticaForRhino")
from ElasticaForRhino import Elastica
import cPickle as pickle
```

In the demonstration we test:

- a discrete Elastica boundary value solver. Given endpoints, end tangents, and curve length we compute a dicrete elastica satisfying these conditions.

```
def add_curve_to_rhino(P, N):
    pts_out = []
    for i in range(0,N*2,2):
        pt = Point3d(P[i],P[i+1],0)
        pts_out.append(pt)
    pl  = Rhino.Geometry.Polyline(pts_out)
    plc = Rhino.Geometry.PolylineCurve(pl)
    scriptcontext.doc.Objects.AddCurve(plc)

e = Elastica()
P_discrete_elastica = Array.CreateInstance(Double,(101)*2)
iter = e.discrete_elastica(0,0, 1,0, 0,1,1,1, 2, 100, P_discrete_elastica)
add_curve_to_rhino(P_discrete_elastica, 100+1)
```

![Boundary solver](boundarysolver.png)

- Estimation of an analytical elastic curve given a set of points.

```
e = Elastica()
curve = rs.GetObject("Choose planar curve (z=0) to estimate elastica",4,True)
domain = rs.CurveDomain(curve)
P_samples = Array.CreateInstance(Double,(101)*2)
for i in range(0,101*2,2):
    t = domain[0]*(1-(i/202))+(i/202)*domain[1]
    P = rs.EvaluateCurve(curve,t)
    P_samples[i] = P[0]
    P_samples[i+1] = P[1]

options = Array.CreateInstance(Double,1)
analytic_elastica_params = Array.CreateInstance(Double,9)
options = Array.CreateInstance(Double,1)
e.compute_elastica_parameters(P_samples, 101, options, analytic_elastica_params)
```

- Evaluation of elastica given the seven parameters analytic_elastica_params.

```
samples = Array.CreateInstance(Double, 2*100)
e.sample_elastica(analytic_elastica_params, -0.15, 1.15, 100, samples)
add_curve_to_rhino(samples, 100)
```

![Estimation and evaluation](estimation_evaluation.png)

- Data-driven design tool of elastic spline
```
e = Elastica()

with open('objs.pickle') as f:  
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

# Initial curve with four elastic curve segments

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
```
![Data-driven elastic spline tool with four segments](elastic_spline_data.png)


### MATLAB

In the folder "MATLAB prototypes" we include MATLAB code for the following features:

- Projection tool in the fo

## Authors

David Brander, J. Andreas Bærentzen, Ann-Sofie Fisker, Jens Gravesen.
Technical University of Denmark.

### References

Approximation by planar elastic curves
David Brander, Jens Gravesen, Toke Bjerge Nørbjerg.
Advances in Computational Mathematics, 2017.

Bézier curves that are close to elastica
David Brander, J. Andreas Bærentzen, Ann-Sofie Fisker, Jens Gravesen.
Computer-Aided Design, 2018.

Designing interactively with elastic splines
David Brander, J. Andreas Bærentzen, Ann-Sofie Fisker, Jens Gravesen.
Computer-Aided Geometric Design, 2018. 

## License
This project is licensed under the MIT License.


