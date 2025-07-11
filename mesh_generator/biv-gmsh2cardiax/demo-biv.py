import cardiax_gmsh2xml as cg
import sys, time, os
import numpy as np


def read_parcellation_from_vtu(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    parcellation_values = []
    inside_parcellation = False

    for line in lines:
        if 'Name="parcellation"' in line:
            inside_parcellation = True
            continue
        if inside_parcellation:
            if "</DataArray>" in line:
                break  # Fim do array
            # Extrair valores numéricos
            numbers = [float(n) for n in line.strip().split()]
            parcellation_values.extend(numbers)

    return parcellation_values

# Edite essa região com os parâmetros do xml =========================
# ====================================================================
# healthy
material_healthy_params = {
   'problemtyp': 'THREE_DIM', 
   'material_type': 'Guccione', 
   'material_coef': [650, 6.62, 3.65, 2.65, 0, 200000.0], 
   'fiber_type': 'fiber_orthotropic', #'fiber_transversely_isotropic', 
   'num_increments': 1
   }

# Fibrotic
material_fibrotic_params = {
   'problemtyp': 'THREE_DIM', 
   'material_type': 'Guccione', 
   'material_coef': [650, 6.62, 3.65, 2.65, 0, 200000.0], 
   'fiber_type': 'fiber_orthotropic', #'fiber_transversely_isotropic', 
   'num_increments': 1
   }

pressure_bc = {
   '10': .0,  # base
   '20': .0,  # lv
   '30': .0,  # rv
   # '40': .0   # epi
}

spring_bc = {
   # '10': .0,  # base
   # '20': .0,  # lv
   # '30': .0,  # rv
   '40': .0   # epi
}

bc_conditions = [pressure_bc, spring_bc]
material_params = [material_healthy_params, material_fibrotic_params]

pvloop_params = {
    'size': 901,
    'total_time': 0.900000,
    'C_art': 0.014,
    'C_ven': 0.3,
    'R_ao': 3850.0,
    'V_ven_zero': 3300.0,
    'P_o': 500.0,
    'B_LA': 0.049,
    'V_art_zero': 580.0,
    'tau': 25.0,
    'R_ven': 1400.0,
    'A_LA': 58.67,
    'R_mv': 1750.0,
    'E_es_LA': 60.0,
    'p_art': 10800.0,
    'p_ven': 1600.0,
    'T_max': 200.0,
    'R_per': 140000.0,
    'stroke_volume': 0.0,
    'T_ref': 1.59
    }

# ====================================================================
# ====================================================================

if (len(sys.argv) < 3):
    print("\n Usage: gmsh2xml <gmsh_mesh> <output_xml> <pvloop_data>\n")
    sys.exit(-1)

gmsh_mesh = sys.argv[1]
output_xml = sys.argv[2]
pvloop_data = sys.argv[3]
aha_mesh = sys.argv[4]

if (not os.path.isfile(gmsh_mesh)):
   print("\n Error: the input gmsh %s does not exist.\n" % (gmsh_mesh))
   sys.exit(-1)

if (not os.path.isfile(pvloop_data)):
   print("\n Error: the input pvloop_data %s does not exist.\n" % (pvloop_data))
   sys.exit(-1)


# Exemplo de uso
aha_list = read_parcellation_from_vtu(aha_mesh)
print(f"{aha_list[:20]=}")  # Mostra os primeiros 20 valo



# convert from mm to m
factor = 1e-3
biv = True
cg.gmsh2xml(gmsh_mesh, output_xml + '_cardiax.xml', aha_list, factor, material_params, bc_conditions,
        pvloop_params, pvloop_data, biv)



