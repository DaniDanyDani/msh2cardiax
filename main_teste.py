import gmsh
import math
import os
import sys


args = sys.argv[1:]

x_c, y_c, z_c, r = 195.7, 180.4, -15.0, 5.0

if len(args[1:4]) >= 4:
    print(f"Número de argumentos inválidos. Esperado: 4 Recebido: {len(args[1:4])}")
    print(f"Execução: python main.py x_c y_c z_c r -nopopup")
    sys.exit()


x_c = float(args[0])
y_c = float(args[1])
z_c = float(args[2])
r = float(args[3]) 

os.system(f"python generate_fibrose.py {x_c} {y_c} {z_c} {r}")

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


gmsh.merge(os.path.join(path, 'fibrose.stl'))

surfaces = gmsh.model.getEntities(dim=2)

surface_tags = [s[1] for s in surfaces]
loop = gmsh.model.geo.addSurfaceLoop(surface_tags, 1)
loop_fibrose = gmsh.model.geo.addSurfaceLoop([surface_tags[4]], 2)
volume = gmsh.model.geo.addVolume([loop], 1)
volume_fibrose = gmsh.model.geo.addVolume([loop_fibrose], 2)

gmsh.model.geo.synchronize()

if len(surface_tags) == 5:
    print(f"{surface_tags=}")
    gmsh.model.addPhysicalGroup(2, [surface_tags[0]], 40, name="epi")
    gmsh.model.addPhysicalGroup(2, [surface_tags[1]], 20, name="ve")
    gmsh.model.addPhysicalGroup(2, [surface_tags[2]], 30, name="vd")
    gmsh.model.addPhysicalGroup(2, [surface_tags[3]], 10, name="base")
    gmsh.model.addPhysicalGroup(3, [volume], 1, name="healthy")
    gmsh.model.addPhysicalGroup(3, [volume_fibrose], 2, name="fibrose")
    gmsh.model.mesh.generate(3)
    gmsh.write(filename + ".msh")
else:
    print("Atenção, número de superfícies está errado, gerando volume sem physical groups")
    print(f"{surface_tags=}")
    gmsh.model.mesh.generate(3)
    gmsh.write(filename + ".msh")


if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
