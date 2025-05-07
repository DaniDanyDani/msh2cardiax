import gmsh
import math
import os
import sys

filename = "Patient_7"
x_coord = 195.7
y_coord = 180.4
z_coord = -15.0
raio = 1.0

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


fibrose_tag = gmsh.model.occ.addSphere(x_coord, y_coord, z_coord, raio, 1)

gmsh.model.occ.synchronize()

surfaces = gmsh.model.getEntities(dim=2)
surface_tags = [s[1] for s in surfaces]
loop = gmsh.model.geo.addSurfaceLoop(surface_tags, 1)
loop_fibrose = gmsh.model.geo.addSurfaceLoop([fibrose_tag], 2)
volume = gmsh.model.geo.addVolume([loop], 1)
volume_fibrose = gmsh.model.geo.addVolume([loop_fibrose], 2)

gmsh.model.geo.synchronize()



if len(surface_tags) == 5:
    print(f"{surface_tags=}")
    gmsh.model.addPhysicalGroup(2, [surface_tags[1]], 40, name="epi")
    gmsh.model.addPhysicalGroup(2, [surface_tags[2]], 20, name="ve")
    gmsh.model.addPhysicalGroup(2, [surface_tags[3]], 30, name="vd")
    gmsh.model.addPhysicalGroup(2, [surface_tags[4]], 10, name="base")
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
