import sys
import os
import time
from math import sin, cos, pi, tan, asin, acos

def criaFibra(baseName):

	from dolfin import *
	from fiberrules import *
	from dolfin_utils.meshconvert import meshconvert

	parameters["form_compiler"]["quadrature_degree"] = 2
	parameters.allow_extrapolation = True

	#convert mesh to dolfin format in order to use in Fiberrules
	ifilename = baseName + '.msh'
	ofilename = baseName + '_fenics.xml'
	iformat = 'gmsh'
	meshconvert.convert2xml(ifilename, ofilename, iformat=iformat)

	#Load mesh for fiberrules
	mesh = Mesh(ofilename)
	domains = MeshFunction("size_t", mesh, baseName + '_fenics_facet_region.xml')

	for facet in facets(mesh):
		mesh.domains().set_marker((facet.index(), domains[facet]), 2)

	#Fibers varying transmurally -60 to 60
	#fiber_angle_epi = 150 #90 - 150 = -60
	#fiber_angle_endo = 30 #90 - 30 = 60
    

    #default values
        fiber_angle_epi = 150. #50.
        fiber_angle_endo = 30. #40.
        sheet_angle_endo = 65.
        sheet_angle_epi = 25.

	fiber_space = FunctionSpace(mesh, 'DG',0)
	#fiber_space = FunctionSpace(mesh, "Lagrange", 1)

	fibers, sheet_normals, cross_sheet = dolfin_fiberrules(mesh, fiber_space, fiber_angle_epi, fiber_angle_endo, sheet_angle_epi, sheet_angle_endo)

	dolfin_to_vtk(fibers, baseName)
	dolfin_to_carpfile(mesh, baseName)
	fibers_to_carpfile(mesh, fibers, baseName)


def criaMalhaTetraedro(baseName):

    cmdstr = '/home/joventino/gmsh/bin/gmsh ' + baseName + '.geo -3 -o ' + baseName + '.msh'

    os.system(cmdstr)

# ------------------------------------------------------------------------------

def convert2cardiax(material_params, pvloop_params, baseName):
	import cardiax_gmsh2xml as c

	c.gmsh2xml(baseName + '.msh', baseName + '_cardiax.xml', material_params, pvloop_params, 'bullseye/pvloop_data.txt')



# ------------------------------------------------------------------------------



def run(a1, c1, bull, material_params, pvloop_params, baseName):


    print baseName, bull

    fileName = baseName + '.geo'


    criaMalhaTetraedro(baseName)

    criaFibra(baseName)

    convert2cardiax(material_params, pvloop_params, baseName)
