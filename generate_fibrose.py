import gmsh
import math
import os
import sys

args = sys.argv[1:]

if len(args) < 4:
    print(f"Número de argumentos inválidos. Esperado: 4 Recebido: {len(args)}")
    print(f"Execução: python generate_fibrose.py x_c y_c z_c r")
    sys.exit()

filename = "Patient_7"

x_coord = float(args[0])
y_coord = float(args[1])
z_coord = float(args[2])
raio = float(args[3]) 

gmsh.initialize()
gmsh.clear()

gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 1/3)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 1.0)

gmsh.model.add("FibroseModel")


fibrose_tag = gmsh.model.occ.addSphere(x_coord, y_coord, z_coord, raio, 1)

gmsh.model.occ.synchronize()


loop_fibrose = gmsh.model.geo.addSurfaceLoop([fibrose_tag], 2)


gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)

gmsh.write("fibrose.stl")

gmsh.finalize()
