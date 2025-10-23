# -*- coding: utf-8 -*-
"""
Este script lê um arquivo de superfície VTK (ex: endocárdio do VE),
extrai os nós da borda e calcula um gradiente para cada nó.
O gradiente é baseado no eixo longo do modelo, que é calculado
usando Análise de Componentes Principais (PCA).

O resultado é salvo em um arquivo CSV, pronto para ser usado pelo
script de geração de "tampa" (cap).
"""
import vtk
from vtk.util import numpy_support as VN
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

def extract_border_and_gradient(vtk_input_filename, csv_output_filename):
    """
    Função principal que executa todo o processo.

    Args:
        vtk_input_filename (str): Caminho para o arquivo VTK da superfície.
        csv_output_filename (str): Caminho onde o arquivo CSV de saída será salvo.
    """

    print(f"Lendo o arquivo: {vtk_input_filename}")
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(vtk_input_filename)
    reader.Update()
    polydata = reader.GetOutput()

    # Extrai todos os pontos da malha para um array NumPy
    all_points_vtk = polydata.GetPoints()
    all_points_np = VN.vtk_to_numpy(all_points_vtk.GetData())

    print("Extraindo as bordas da malha...")
    feature_edges = vtk.vtkFeatureEdges()
    feature_edges.SetInputData(polydata)
    feature_edges.BoundaryEdgesOn()  # Queremos APENAS as bordas de fronteira
    feature_edges.FeatureEdgesOff()
    feature_edges.ManifoldEdgesOff()
    feature_edges.NonManifoldEdgesOff()
    feature_edges.Update()

    border_points_np = VN.vtk_to_numpy(feature_edges.GetOutput().GetPoints().GetData())
    print(f"Encontrados {len(border_points_np)} pontos na borda.")

    print("Calculando o eixo longo via PCA...")
    center = np.mean(all_points_np, axis=0)
    centered_points = all_points_np - center
    
    covariance_matrix = np.cov(centered_points, rowvar=False)
    
    eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)
    
    long_axis_vector = eigenvectors[:, -1]
    print(f"Vetor do eixo longo: {long_axis_vector}")

    scalar_values = np.dot(all_points_np, long_axis_vector)

    scalar_array = VN.numpy_to_vtk(scalar_values, deep=True)
    scalar_array.SetName("LongAxisProjection")
    polydata.GetPointData().SetScalars(scalar_array)

    print("Calculando o gradiente do campo escalar...")
    gradient_filter = vtk.vtkGradientFilter()
    gradient_filter.SetInputData(polydata)
    gradient_filter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "LongAxisProjection")
    gradient_filter.SetResultArrayName("GradientVectors")
    gradient_filter.Update()

    output_with_gradients = gradient_filter.GetOutput()
    all_gradients_np = VN.vtk_to_numpy(output_with_gradients.GetPointData().GetArray("GradientVectors"))

    print("Mapeando pontos da borda aos seus IDs e gradientes...")
    kdtree = cKDTree(all_points_np)
    distances, original_indices = kdtree.query(border_points_np)

    border_gradients_np = all_gradients_np[original_indices]

    print(f"Salvando dados da borda em: {csv_output_filename}")
    df = pd.DataFrame({
        'PointIds': original_indices,
        'X': border_points_np[:, 0],
        'Y': border_points_np[:, 1],
        'Z': border_points_np[:, 2],
        'ScalarGradient:0': border_gradients_np[:, 0],
        'ScalarGradient:1': border_gradients_np[:, 1],
        'ScalarGradient:2': border_gradients_np[:, 2]
    })
    df.to_csv(csv_output_filename, index=False)
    print("Processo concluído com sucesso!")

    # fig = plt.figure(figsize=(8, 8))
    # ax = fig.add_subplot(projection='3d')
    
    # ax.scatter(all_points_np[::20, 0], all_points_np[::20, 1], all_points_np[::20, 2], c='gray', alpha=0.2, label='Superfície VE')
    
    # ax.scatter(border_points_np[:, 0], border_points_np[:, 1], border_points_np[:, 2], c='r', marker='o', label='Pontos da Borda')
    
    # ax.quiver(
    #     border_points_np[:, 0], border_points_np[:, 1], border_points_np[:, 2],
    #     border_gradients_np[:, 0], border_gradients_np[:, 1], border_gradients_np[:, 2],
    #     length=2.0,  # Ajuste o comprimento para melhor visualização
    #     normalize=True,
    #     color='b',
    #     label='Gradientes (Eixo Longo)'
    # )
    
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    # ax.legend()
    # plt.title("Visualização da Borda e Gradientes")
    # plt.show()
    


if __name__ == '__main__':

    print("\nGerando csv da borda para o lv")
    INPUT_VTK_FILE = 'lv.vtk' 
    OUTPUT_CSV_FILE = 'lv_edge_node_start_vectors.csv'
    extract_border_and_gradient(INPUT_VTK_FILE, OUTPUT_CSV_FILE)

    print("\nGerando csv da borda para o rv")
    INPUT_VTK_FILE = 'rv.vtk' 
    OUTPUT_CSV_FILE = 'rv_edge_node_start_vectors.csv'
    extract_border_and_gradient(INPUT_VTK_FILE, OUTPUT_CSV_FILE)

    print("\nGerando csv da borda para o epi")
    INPUT_VTK_FILE = 'epi.vtk' 
    OUTPUT_CSV_FILE = 'epi_edge_node_start_vectors.csv'
    extract_border_and_gradient(INPUT_VTK_FILE, OUTPUT_CSV_FILE)
