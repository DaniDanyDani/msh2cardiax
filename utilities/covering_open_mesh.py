import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import vtk
from vtk.util import numpy_support as VN
import pandas

def vector_angle(x, start_angle, flat_x):
    if x <= flat_x:
        return start_angle/(flat_x)**2 * (x-flat_x)**2 # This is a quadratic function that begins at start angle and goes to zero at end_x
    else:
        return 0

def plot_3D_vector(ax, coord, vector, color):
    ax.quiver(coord[0], coord[1], coord[2], vector[0], vector[1], vector[2], color=color)

def save_data(save_name, data):
    print("Salvando nuvem de pontos em txt")
    np.savetxt(save_name, data)
    
    
    print("Salvando nuvem de pontos em vtk")
    vtk_save_name = save_name.replace('.txt', '.vtk')
    vtk_points = vtk.vtkPoints()
    vtk_points.SetData(VN.numpy_to_vtk(data, deep=True))

    output_polydata = vtk.vtkPolyData()
    output_polydata.SetPoints(vtk_points)

    verts = vtk.vtkCellArray()
    for i in range(output_polydata.GetNumberOfPoints()):
        verts.InsertNextCell(1)      # Célula com 1 ponto
        verts.InsertCellPoint(i) # Adiciona o ponto de índice 'i'
    output_polydata.SetVerts(verts)

    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(vtk_save_name)
    writer.SetInputData(output_polydata)
    writer.Write()


def generate_cap_surface_data(edge_nodes, start_vectors, longitudinal_vector, lid_rise, flat_percentage, resolution_along_radius):
    # Create curves from each edge node to the middle
    middle_coord = np.mean(edge_nodes, axis=0)
    middle_lid_coord = middle_coord + longitudinal_vector*lid_rise

    data = np.array([[0, 0, 0]])
    for edge_node, start_vector in zip(edge_nodes, start_vectors):
        x_vector = middle_coord - edge_node
        y_axis = np.cross(x_vector/np.linalg.norm(x_vector), longitudinal_vector)
        y_axis = y_axis/np.linalg.norm(y_axis)
        x_axis = np.cross(longitudinal_vector, y_axis)
        x_axis = x_axis/np.linalg.norm(x_axis)
        dx = resolution_along_radius * (x_vector)/np.linalg.norm(x_vector)

        n_points = np.floor(np.linalg.norm(x_vector)/ resolution_along_radius).astype(int)
        curve = np.zeros((n_points, 3))
        curve[0, :] = edge_node
        flat_x = flat_percentage * np.linalg.norm(x_vector)
        end_vector = x_vector / np.linalg.norm(x_vector)
        start_angle = np.arccos(np.dot(start_vector, end_vector))

        current_coord = edge_node
        for i in range(1, n_points):
            angle = vector_angle(i * resolution_along_radius, start_angle, flat_x)
            basis_transformation = np.vstack((x_axis, y_axis, longitudinal_vector)).transpose()
            R = np.array([[np.cos(angle), 0, -np.sin(angle)], [0, 1, 0], [np.sin(angle), 0, np.cos(angle)]])
            new_vector = np.matmul(R, np.matmul(np.linalg.inv(basis_transformation),end_vector))
            current_coord = current_coord + dx
            new_vector_global = np.matmul(basis_transformation, new_vector)
            new_vector_global = new_vector_global/np.linalg.norm(new_vector_global)

            curve[i, :] = curve[i-1,:] + new_vector_global * resolution_along_radius

        lid_scale = lid_rise / (np.dot(longitudinal_vector, curve[-1,:] - edge_node))
        x_scale =  np.linalg.norm(x_vector)/np.dot(x_axis, curve[-1,:] - edge_node)
        for i in range(0, n_points):
            curve_x = np.dot((curve[i,:]-edge_node), x_axis)
            curve_y = np.dot((curve[i,:]-edge_node), y_axis)
            curve_z = np.dot((curve[i,:]-edge_node), longitudinal_vector)
            curve[i, :] = (curve_x* x_axis*x_scale + curve_y*y_axis + curve_z * longitudinal_vector * lid_scale) + edge_node
        curve = curve[1:, :]


        data = np.concatenate((data,curve))
    data = data[1:,:]
    return data
# def generate_cap_surface_data(edge_nodes, start_vectors, longitudinal_vector,
#                               lid_rise, flat_percentage, n_points_per_curve = 100):
#     """
#     Gera uma nuvem de pontos para a tampa, controlada pelo
#     NÚMERO de pontos por raio, em vez da resolução.
#     """
    
#     middle_coord = np.mean(edge_nodes, axis=0)
#     data = np.array([[0, 0, 0]]) # Array placeholder para os pontos da tampa

#     # Garante que o número de pontos é um inteiro
#     n_points = int(n_points_per_curve)
    
#     # Precisamos de pelo menos 2 pontos (início e fim) para fazer uma curva
#     if n_points <= 1:
#         print("AVISO: n_points_per_curve deve ser pelo menos 2.")
#         return np.array([[]]) # Retorna dados vazios

#     for edge_node, start_vector in zip(edge_nodes, start_vectors):
#         x_vector = middle_coord - edge_node
#         radius = np.linalg.norm(x_vector) # Comprimento total da curva (raio)
        
#         # Pula se o ponto da borda estiver (quase) no centro
#         if radius < 1e-6:
#             continue

#         # --- MUDANÇA PRINCIPAL ---
#         # Agora, calculamos o 'resolution' com base no 'n_points' desejado
#         # Há (n_points - 1) segmentos na curva
#         resolution_along_radius = radius / (n_points - 1)
#         # --- FIM DA MUDANÇA ---
        
#         # Cria eixos locais
#         end_vector = x_vector / radius
#         y_axis = np.cross(end_vector, longitudinal_vector)
#         y_axis = y_axis / np.linalg.norm(y_axis)
#         x_axis = np.cross(longitudinal_vector, y_axis)
#         x_axis = x_axis / np.linalg.norm(x_axis)

#         # n_points agora vem do argumento da função
#         curve = np.zeros((n_points, 3))
#         curve[0, :] = edge_node
        
#         # Distância do centro onde a curva começa a ficar plana
#         flat_x = flat_percentage * radius
        
#         # Adiciona np.clip para segurança (evita erros de arccos com float)
#         dot_product = np.clip(np.dot(start_vector, end_vector), -1.0, 1.0)
#         start_angle = np.arccos(dot_product)

#         for i in range(1, n_points):
#             # A função vector_angle espera a distância do *centro*.
#             # 'i' vai de 1 (perto da borda) a n_points-1 (perto do centro)
#             dist_from_edge = i * resolution_along_radius
#             dist_from_center = radius - dist_from_edge
            
#             angle = vector_angle(dist_from_center, start_angle, flat_x)
            
#             basis_transformation = np.vstack((x_axis, y_axis, longitudinal_vector)).transpose()
#             R = np.array([[np.cos(angle), 0, -np.sin(angle)], [0, 1, 0], [np.sin(angle), 0, np.cos(angle)]])
            
#             # Aplica a rotação ao vetor que aponta para o centro
#             new_vector_local = np.matmul(np.linalg.inv(basis_transformation), end_vector)
#             rotated_vector_local = np.matmul(R, new_vector_local)
#             new_vector_global = np.matmul(basis_transformation, rotated_vector_local)
#             new_vector_global = new_vector_global / np.linalg.norm(new_vector_global)

#             # Adiciona o novo ponto ao longo da direção do vetor
#             curve[i, :] = curve[i-1,:] + new_vector_global * resolution_along_radius

#         # --- Lógica de Escala (para garantir altura e raio exatos) ---
#         # Adiciona verificação de segurança para divisão por zero
        
#         delta_z = np.dot(longitudinal_vector, curve[-1,:] - edge_node)
#         delta_x = np.dot(x_axis, curve[-1,:] - edge_node)
        
#         if abs(delta_z) < 1e-6: delta_z = 1e-6
#         if abs(delta_x) < 1e-6: delta_x = 1e-6
            
#         lid_scale = lid_rise / delta_z
#         x_scale =  radius / delta_x
        
#         for i in range(0, n_points):
#             curve_x = np.dot((curve[i,:]-edge_node), x_axis)
#             curve_y = np.dot((curve[i,:]-edge_node), y_axis)
#             curve_z = np.dot((curve[i,:]-edge_node), longitudinal_vector)
#             curve[i, :] = (curve_x* x_axis*x_scale + curve_y*y_axis + curve_z * longitudinal_vector * lid_scale) + edge_node
        
#         # Exclui o primeiro ponto (que já está na borda)
#         curve = curve[1:, :]

#         data = np.concatenate((data,curve))
        
#     data = data[1:,:] # Remove o placeholder inicial
#     return data

if __name__ == '__main__':
    vtk_names = ['lv.vtk', 'rv.vtk', 'epi.vtk']
    edge_names = ['lv_edge_node_start_vectors.csv', 'rv_edge_node_start_vectors.csv', 'epi_edge_node_start_vectors.csv']
    save_names = ['lv_cap.txt', 'rv_cap.txt', 'epi_cap.txt']
    lid_rise_percentages = [0.015, 0.015, 0.02]
    normal_vectors = [np.array([0, 0, 1]), np.array([0, 0, 1]), np.array([0, 0, 1])]


    for vtk_name, edge_name, save_name, lid_rise_percentage, normal_vector in zip(vtk_names, edge_names, save_names, lid_rise_percentages, normal_vectors):
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(vtk_name)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()

        nodes_xyz = VN.vtk_to_numpy(data.GetPoints().GetData())  # Convert from mm to cm
        # fig = plt.figure(figsize=(8, 8))
        # ax = fig.add_subplot(projection='3d')
        # ax.scatter(nodes_xyz[:, 0], nodes_xyz[:, 1], nodes_xyz[:, 2], marker='*', c='b')

        num_points = data.GetNumberOfPoints()
        node_ids = np.arange(num_points)
        print(node_ids)
        existing_lv_length = np.amax(np.abs(np.sum(normal_vector * nodes_xyz, axis=1, keepdims=True)))
        lid_rise = existing_lv_length * lid_rise_percentage
        data = pandas.read_csv(edge_name)
        edge_node_idx = data['PointIds'].values
        start_vectors = np.zeros((edge_node_idx.shape[0], 3))
        start_vectors[:, 0] = data['ScalarGradient:0'].values
        start_vectors[:, 1] = data['ScalarGradient:1'].values
        start_vectors[:, 2] = data['ScalarGradient:2'].values

        edge_nodes = np.zeros((edge_node_idx.shape[0], 3))
        for i, idx in enumerate(edge_node_idx):
            edge_nodes[i, :] = nodes_xyz[np.where(node_ids==idx)[0],:]

        cap_data = generate_cap_surface_data(edge_nodes=edge_nodes, start_vectors=start_vectors * 0.1,
                                             longitudinal_vector=normal_vector, lid_rise=lid_rise,
                                             flat_percentage=0.9,resolution_along_radius=0.5)

        # ax.scatter(cap_data[:, 0], cap_data[:, 1], cap_data[:, 2], marker='*', c='m')
        # plt.show()

        save_data(save_name, cap_data)

        reader = None

# Use MeshLab to generate the surface mesh from this point on
