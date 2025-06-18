
import bullseye.general_lvspline_gmsh_17 as lvg 


def model(material_params, bull_params, pvloop_params, baseName):

    lvg.run(bull_params[0], bull_params[1], bull_params[2:], material_params, pvloop_params, baseName)




pvloop_params = {'C_art': 0.001, 'R_per': 20000.0, 'P_o': 500.0, 'p_art': 7000.0, 'stroke_volume': 0.0}
#pvloop_params = {'C_art': 0.001, 'R_per': 20000.0, 'P_o': 1000.0, 'p_art': 10800.0, 'stroke_volume': 0.0}

bull_henrik = [0.065, 0.0826, 0.0058, 0.0064, 0.0102, 0.0108, 0.0084, 0.0061, 0.0092, 0.0081, 0.0102, 0.0108, 0.0082, 0.0081, 0.0157, 0.0157, 0.0102, 0.0127, 0.0157]
#bull_params_bai = [0.04, 0.06, 0.00621, 0.00638, 0.00623, 0.00582, 0.00538, 0.00666, 0.00640, 0.00706, 0.00841, 0.00664, 0.00594, 0.00692, 0.00547, 0.00619, 0.00619, 0.00540, 0.00587, 0.00587, 0.00437]
bull_params_bai = [0.04, 0.06, 0.00621, 0.00638, 0.00623, 0.00582, 0.00538, 0.00666, 0.00640, 0.00706, 0.00841, 0.00664, 0.00594, 0.00692, 0.00547, 0.00619, 0.00540, 0.00587, 0.00437]

bull_params_p2 = [0.055, 0.063, 0.007, 0.006, 0.006, 0.007, 0.008, 0.008, 0.007, 0.006, 0.006, 0.007, 0.008, 0.008, 0.007, 0.006, 0.007, 0.008, 0.005]

#bull_params_bai = bull_henrik

#bull_params_bai = [0.04, 0.06]
#for i in range(17):
#    bull_params_bai.append(0.00640)

#bull_params_bai[1+2] *= 0.8
#bull_params_bai[2+2] *= 0.8
#bull_params_bai[7+2] *= 1.3
#bull_params_bai[8+2] *= 1.3
#bull_params_bai[13+2] *= 1.3
#bull_params_bai[4+2] *= 1.3
#bull_params_bai[5+2] *= 1.3
#bull_params_bai[11+2] *= 0.8
#bull_params_bai[10+2] *= 0.8
#bull_params_bai[15+2] *= 1.3




#bull_params_bai = [0.04, 0.06, 0.00621, 0.00666, 0.00538, 0.00582, 0.00623, 0.00638, 0.00640, 0.00692, 0.00594, 0.00664, 0.00841, 0.00706, 0.00547, 0.00587, 0.00540, 0.00619, 0.00437]
material_params = {'problemtyp': 'THREE_DIM', 'material_type': 'Guccione', 'material_coef': [1.1, 6.6, 4.0, 2.6, 0, 300.], 'fiber_type': 'transversely_isotropic', 'num_increments': 1, 'pressure_marker': 30, 'pressure_value': 2.0}
#material_params = {'problemtyp': 'THREE_DIM', 'material_type': 'HolzapfelOgden', 'material_coef': [0.330, 18.535, 2.564, 0.417, 9.242, 15.972, 10.446, 11.602, 0.0, 33333.0], 'fiber_type': 'orthotropic', 'num_increments': 10, 'pressure_marker': 30, 'pressure_value': 2.0}



model(material_params, bull_params_p2, pvloop_params, 'LV_mesh')

