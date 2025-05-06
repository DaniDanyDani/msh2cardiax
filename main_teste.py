import gmsh
import math
import os
import sys

filename = "Patient_10"

gmsh.initialize()
gmsh.clear()

# gmsh.option.setNumber("General.Terminal", 1)
# gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.75)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 2.25)
# gmsh.option.setNumber("Mesh.Optimize",1)
# gmsh.option.setNumber("Mesh.QualityType",2)

gmsh.model.add("TestModel")

path = os.path.dirname(os.path.abspath(__file__))
gmsh.merge(os.path.join(path, 'epi.stl'))
gmsh.merge(os.path.join(path, 'lv.stl'))
gmsh.merge(os.path.join(path, 'rv.stl'))
gmsh.merge(os.path.join(path, 'base.stl'))

gmsh.model.mesh.removeDuplicateNodes()
gmsh.model.mesh.refine()

includeBoundary = True
forceParametrizablePatches=False
angle = 40

curveAngle = 180

gmsh.model.mesh.classifySurfaces(angle * math.pi / 180., includeBoundary,
                                forceParametrizablePatches,
                                curveAngle * math.pi / 180.)
# curveAngle * math.pi / 180.

gmsh.model.mesh.createGeometry()

surfaces = gmsh.model.getEntities(dim=2)
print(f"Superfícies detectadas: {len(surfaces)}")

surface_tags = [s[1] for s in surfaces]
loop = gmsh.model.geo.addSurfaceLoop(surface_tags)
volume = gmsh.model.geo.addVolume([loop])


gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(3)
gmsh.write(filename + ".msh")



# # gmsh.model.add("t1")

# # lc = 1e-2
# # gmsh.model.geo.addPoint(0, 0, 0, lc, 1)


# # gmsh.model.geo.addPoint(.1, 0, 0, lc, 2)
# # gmsh.model.geo.addPoint(.1, .3, 0, lc, 3)

# # p4 = gmsh.model.geo.addPoint(0, .3, 0, lc)

# # gmsh.model.geo.addLine(1, 2, 1)
# # gmsh.model.geo.addLine(3, 2, 2)
# # gmsh.model.geo.addLine(3, p4, 3)
# # gmsh.model.geo.addLine(4, 1, p4)


# # gmsh.model.geo.addCurveLoop([4, 1, -2, 3], 1)


# # gmsh.model.geo.addPlaneSurface([1], 1)

# # gmsh.model.geo.synchronize()

# # gmsh.model.addPhysicalGroup(1, [1, 2, 4], 5)
# # gmsh.model.addPhysicalGroup(2, [1], name="My surface")


# # gmsh.model.mesh.generate(2)


# # gmsh.write("t1.msh")


if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
