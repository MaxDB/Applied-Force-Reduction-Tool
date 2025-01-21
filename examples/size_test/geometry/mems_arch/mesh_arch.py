#call from console as 'python mesh_arch.py SEED_SIZE'

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import *

import sys
from os import getcwd

seed_size = float(sys.argv[-1])

current_directory = getcwd()
current_directory = current_directory.replace("\\", "/")
file_path = current_directory.split("/")
file_path.pop(-1)
file_path.append("geometry/mems_arch/mems arch.cae")
model_path = "/".join(file_path)


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
 