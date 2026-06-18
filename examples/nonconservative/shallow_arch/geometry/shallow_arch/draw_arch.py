import math
unit_factor = 1e6
num_points = 201

rise = 3.84e-6*unit_factor
length = 640e-6*unit_factor
thickness = 6.4e-6*unit_factor
width = 32e-6*unit_factor

#---
amplitude = rise/2
frequency = 2*math.pi/length

phase = -math.pi/2

y = lambda x: amplitude*math.sin(frequency*x + phase)
#----

x = [i*length/(num_points-1) for i in range(num_points)]

xy = tuple((x[i], y(x[i])) for i in range(num_points))




# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2024 replay file
# Internal Version: 2023_09_21-13.55.25 RELr426 190762
# Run by gg19546 on Tue Jun 16 16:26:41 2026
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=309.612487792969, 
    height=201.170379638672)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
openMdb('shallow_arch.cae')
#: The model database "C:\Users\gg19546\Documents\PhD\Year 2\Software\Applied Force Reduction\examples\nonconservative\shallow_arch\geometry\shallow_arch\shallow_arch.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=length)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
#splines
s.Spline(points=xy)

#copy splines
s.copyMove(vector=(0.0, thickness), objectList=(g[2], ))

#add connecting lines
s.Line(point1=(0.0, y(x[0])), point2=(0.0, y(x[0])+thickness))
s.VerticalConstraint(entity=g[4], addUndoState=False)

s.Line(point1=(length, y(x[-1])), point2=(length, y(x[-1])+thickness))
s.VerticalConstraint(entity=g[5], addUndoState=False)

#extrude
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseSolidExtrude(sketch=s, depth=width)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
mdb.save()
