import gmsh
import math
import os
import sys

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
volume = gmsh.model.geo.addVolume([loop], 2)

volumes = gmsh.model.getEntities(dim=3)
volume_tags = [v[1] for v in volumes]

gmsh.model.geo.synchronize()


if len(surface_tags) == 4:
    print(f"{surface_tags=}")
    print(f"{volume_tags=}")
    gmsh.model.addPhysicalGroup(2, [surface_tags[0]], 40, name="epi")
    gmsh.model.addPhysicalGroup(2, [surface_tags[1]], 20, name="ve")
    gmsh.model.addPhysicalGroup(2, [surface_tags[2]], 30, name="vd")
    gmsh.model.addPhysicalGroup(2, [surface_tags[3]], 10, name="base")
    gmsh.model.addPhysicalGroup(3, [volume], 1, name="healthy")
    gmsh.model.mesh.generate(3)
    gmsh.write(filename + ".msh")
else:
    print("Atenção, número de superfícies está errado, gerando volume sem physical groups")
    print(f"{surface_tags=}")
    print(f"{volume_tags=}")
    gmsh.model.mesh.generate(3)
    gmsh.write(filename + ".msh")


if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
