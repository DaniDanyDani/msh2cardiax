import gmsh
import math
import os
import sys


args = sys.argv[1:]

if len(args[1:4]) >= 4:
    print(f"Número de argumentos inválidos. Esperado: 4 Recebido: {len(args[1:4])}")
    print(f"Execução: python main.py x_c y_c z_c r -nopopup")
    sys.exit()


if float(args[0]) >= 0:
    x_c = "-"+str(float(args[0])) 
else:
    x_c = "+"+str(-1 * float(args[0])) 
if float(args[1]) >= 0:
    y_c = "-"+str(args[1]) 
else:
    y_c = "+"+str(-1 * float(args[1])) 
if float(args[2]) >= 0:
    z_c = "-"+str(args[2]) 
else:
    z_c = "+"+str(-1 * float(args[2]))
if float(args[3]) >= 0:
    r = str(args[3]) 
else:
    r = str(-1 * float(args[3]))


filename = "Patient_7"

gmsh.initialize()
gmsh.clear()

gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 1)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 3)

gmsh.model.add("TestModel")

path = os.path.dirname(os.path.abspath(__file__))
gmsh.merge(os.path.join(path, 'epi0.stl'))
gmsh.merge(os.path.join(path, 'lv0.stl'))
gmsh.merge(os.path.join(path, 'rv0.stl'))
gmsh.merge(os.path.join(path, 'base0.stl'))


gmsh.model.mesh.removeDuplicateNodes()
gmsh.model.mesh.refine()

includeBoundary = True
forceParametrizablePatches=False
angle = 50

curveAngle = 180

gmsh.model.mesh.classifySurfaces(angle * math.pi / 180., includeBoundary,
                                forceParametrizablePatches,
                                curveAngle * math.pi / 180.)

gmsh.model.mesh.createGeometry()


surfaces = gmsh.model.getEntities(dim=2)

surface_tags = [s[1] for s in surfaces]
loop = gmsh.model.geo.addSurfaceLoop(surface_tags, 1)
volume = gmsh.model.geo.addVolume([loop], 1)

gmsh.model.geo.synchronize()


sphere_equation = f"2*Sin((x+y)/5) + 3"
# sphere_equation = f"(x-195.7)*(x-195.7) + (y-180.4)*(y-180.4) + (z+15)*(z+15)"


print(f"{sphere_equation=}")

fibrose_tag_field = gmsh.model.mesh.field.add("MathEval", 3)
gmsh.model.mesh.field.setString(fibrose_tag_field, "F",
                                sphere_equation)

gmsh.model.mesh.field.setAsBackgroundMesh(fibrose_tag_field)

# gmsh.model.addPhysicalGroup(2, [surface_tags[0]], 40, name="epi")
# gmsh.model.addPhysicalGroup(2, [surface_tags[1]], 20, name="ve")
# gmsh.model.addPhysicalGroup(2, [surface_tags[2]], 30, name="vd")
# gmsh.model.addPhysicalGroup(2, [surface_tags[3]], 10, name="base")
# gmsh.model.addPhysicalGroup(3, [volume], 1, name="healthy")

# gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
# gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
# gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

gmsh.model.mesh.generate(3)
gmsh.write(filename + ".msh")



if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
