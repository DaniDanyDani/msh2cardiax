import cardiax_gmsh2xml as cg
import sys, time, os
import numpy as np


baseName = "malha-fibrose-v2"

material_params = {'problemtyp': 'THREE_DIM', 
'material_type': 'Guccione', 
'material_coef': [1.1, 6.6, 4.0, 2.6, 0, 300.], 
'fiber_type': 'transversely_isotropic', 
'num_increments': 1, 
'pressure_marker': 30, 
'pressure_value': 2.0}

pvloop_params = {'C_art': 0.001, 
'R_per': 20000.0, 
'P_o': 1000.0, 
'p_art': 10800.0, 
'stroke_volume': 0.0}



cg.gmsh2xml(baseName + '.msh', baseName + '_cardiax.xml', material_params, markers_boundbox, pvloop_params, './pvloop_data.txt')