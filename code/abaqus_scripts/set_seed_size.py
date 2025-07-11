import sys
import os

model_path = sys.argv[-2]
seed_size = float(sys.argv[-1])

pwd = os.getcwd().split('\\')
pwd = '\\'.join(pwd[:-1])

model_path = repr(pwd).strip("'") + model_path
# This should not have to be this complicated

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import *


executeOnCaeStartup()

# Model
openMdb(pathName= model_path)
model_name = list(mdb.models.keys())[0]
part_name = list(mdb.models[model_name].parts.keys())[0]
part = mdb.models[model_name].parts[part_name]

# Mesh
#elem = mesh.ElemType(elemCode=C3D15)
#part.setElementType(regions=(part.cells, ), elemTypes=(elem,))
part.seedPart(size=seed_size)
part.generateMesh()

# Job
job = mdb.Job(name='remesh', model=model_name)
job.writeInput()


