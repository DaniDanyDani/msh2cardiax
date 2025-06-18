import cardiax_gmsh2xml as cg
import sys, time, os
import numpy as np


material_params = {'problemtyp': 'THREE_DIM', 
'material_type': 'Guccione', 
'material_coef': [1.1, 6.6, 4.0, 2.6, 0, 300.], 
'fiber_type': 'fiber_transversely_isotropic', 
'num_increments': 1, 
'pressure_marker': 30, 
'pressure_value': 2.0}

material_params = {'problemtyp': 'THREE_DIM', 
        'material_type': 'HolzapfelOgden', 
        'material_coef': [2000, 3000, 2000, 1000, 9.242, 15.972, 10.446, 11.602, 0, 500000.0], 
        'fiber_type': 'fiber_orthotropic', 
        'num_increments': 100, 
        'pressure_marker': 30, 'pressure_value': 2000.0}

pvloop_params = {'C_art': 0.001, 
'R_per': 20000.0, 
'P_o': 1000.0, 
'p_art': 10800.0, 
'stroke_volume': 0.0}


if (len(sys.argv) < 3):
    print("\n Usage: gmsh2xml <gmsh_mesh> <output_xml> <pvloop_data>\n")
    sys.exit(-1)

gmsh_mesh = sys.argv[1]
output_xml = sys.argv[2]
pvloop_data = sys.argv[3]

if (not os.path.isfile(gmsh_mesh)):
   print("\n Error: the input gmsh %s does not exist.\n" % (gmsh_mesh))
   sys.exit(-1)

if (not os.path.isfile(pvloop_data)):
   print("\n Error: the input pvloop_data %s does not exist.\n" % (pvloop_data))
   sys.exit(-1)


#convert from mm to m
factor = 1e-3
biv = False
cg.gmsh2xml(gmsh_mesh, output_xml + '_cardiax.xml', factor, material_params,
        pvloop_params, pvloop_data, biv)



