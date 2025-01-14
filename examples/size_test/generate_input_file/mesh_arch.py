#call from console as 'python mesh_arch.py SEED_SIZE'

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import *

import sys
from os import getcwd

seed_size = float(sys.argv[-1])

current_directory = getcwd()
model_path = current_directory + "\mems arch.cae"
model_path = model_path.replace("\\", "/")

executeOnCaeStartup()

# Model
openMdb(pathName= model_path)
part = mdb.models['Model-1'].parts['Arch']

# Mesh
elem = mesh.ElemType(elemCode=C3D15)
part.setElementType(regions=(part.cells, ), elemTypes=(elem,))
part.seedPart(size=seed_size)
part.generateMesh()

# Job
job = mdb.Job(name='mems_arch', model='Model-1')
job.writeInput()
 