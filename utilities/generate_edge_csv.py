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
    # --- 1. Ler o arquivo VTK ---
    # Usamos vtkPolyDataReader, que é comum para arquivos de superfície.
    # Se o seu arquivo for UnstructuredGrid, use vtkUnstructuredGridReader.
    print(f"Lendo o arquivo: {vtk_input_filename}")
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(vtk_input_filename)
    reader.Update()
    polydata = reader.GetOutput()

    # Extrai todos os pontos da malha para um array NumPy
    all_points_vtk = polydata.GetPoints()
    all_points_np = VN.vtk_to_numpy(all_points_vtk.GetData())

    # --- 2. Extrair os Pontos da Borda ---
    # vtkFeatureEdges é o filtro perfeito para encontrar bordas de uma malha.
    print("Extraindo as bordas da malha...")
    feature_edges = vtk.vtkFeatureEdges()
    feature_edges.SetInputData(polydata)
    feature_edges.BoundaryEdgesOn()  # Queremos APENAS as bordas de fronteira
    feature_edges.FeatureEdgesOff()
    feature_edges.ManifoldEdgesOff()
    feature_edges.NonManifoldEdgesOff()
    feature_edges.Update()

    # O resultado contém as linhas que formam a borda.
    # Vamos extrair as coordenadas dos pontos dessas linhas.
    border_points_np = VN.vtk_to_numpy(feature_edges.GetOutput().GetPoints().GetData())
    print(f"Encontrados {len(border_points_np)} pontos na borda.")

    # --- 3. Calcular o Eixo Longo (PCA) e Criar Campo Escalar ---
    print("Calculando o eixo longo via PCA...")
    # Centraliza os pontos (subtrai a média)
    center = np.mean(all_points_np, axis=0)
    centered_points = all_points_np - center
    
    # Calcula a matriz de covariância
    covariance_matrix = np.cov(centered_points, rowvar=False)
    
    # Calcula os autovetores e autovalores
    eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)
    
    # O eixo longo é o autovetor correspondente ao maior autovalor
    long_axis_vector = eigenvectors[:, -1]
    print(f"Vetor do eixo longo: {long_axis_vector}")

    # Cria o campo escalar: projeção de cada ponto no eixo longo
    scalar_values = np.dot(all_points_np, long_axis_vector)

    # Adiciona este novo campo escalar de volta ao objeto polydata do VTK
    scalar_array = VN.numpy_to_vtk(scalar_values, deep=True)
    scalar_array.SetName("LongAxisProjection")
    polydata.GetPointData().SetScalars(scalar_array)

    # --- 4. Calcular o Gradiente do Campo Escalar ---
    print("Calculando o gradiente do campo escalar...")
    gradient_filter = vtk.vtkGradientFilter()
    gradient_filter.SetInputData(polydata)
    gradient_filter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "LongAxisProjection")
    gradient_filter.SetResultArrayName("GradientVectors")
    gradient_filter.Update()

    # Extrai os vetores de gradiente para TODOS os pontos
    output_with_gradients = gradient_filter.GetOutput()
    all_gradients_np = VN.vtk_to_numpy(output_with_gradients.GetPointData().GetArray("GradientVectors"))

    # --- 5. Mapear os Pontos da Borda aos seus Gradientes e IDs Originais ---
    print("Mapeando pontos da borda aos seus IDs e gradientes...")
    # Usamos uma k-d tree para encontrar rapidamente o índice original de cada ponto da borda
    kdtree = cKDTree(all_points_np)
    distances, original_indices = kdtree.query(border_points_np)

    # Seleciona os gradientes que correspondem apenas aos pontos da borda
    border_gradients_np = all_gradients_np[original_indices]

    # --- 6. Salvar os resultados no formato CSV ---
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

    # --- 7. (Opcional) Visualização para Verificação ---
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection='3d')
    
    # Plota um subconjunto de todos os pontos para dar contexto
    ax.scatter(all_points_np[::20, 0], all_points_np[::20, 1], all_points_np[::20, 2], c='gray', alpha=0.2, label='Superfície VE')
    
    # Plota os pontos da borda
    ax.scatter(border_points_np[:, 0], border_points_np[:, 1], border_points_np[:, 2], c='r', marker='o', label='Pontos da Borda')
    
    # Plota os vetores de gradiente (setas) em cada ponto da borda
    ax.quiver(
        border_points_np[:, 0], border_points_np[:, 1], border_points_np[:, 2],
        border_gradients_np[:, 0], border_gradients_np[:, 1], border_gradients_np[:, 2],
        length=2.0,  # Ajuste o comprimento para melhor visualização
        normalize=True,
        color='b',
        label='Gradientes (Eixo Longo)'
    )
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    plt.title("Visualização da Borda e Gradientes")
    # Para ajustar a proporção dos eixos e ter uma visualização correta
    # ax.set_aspect('equal', 'box')
    plt.show()
    


if __name__ == '__main__':
    # --- CONFIGURAÇÃO ---
    # Coloque o nome do seu arquivo VTK aqui
    INPUT_VTK_FILE = 'lv.vtk' 
    # Nome do arquivo CSV que será gerado
    OUTPUT_CSV_FILE = 'lv_edge_node_start_vectors.csv'

    # Para testar, vamos criar um arquivo VTK falso em formato de cilindro aberto
    # print("Criando um arquivo VTK de teste ('lv.vtk')...")
    # cylinder = vtk.vtkCylinderSource()
    # cylinder.SetResolution(50)
    # cylinder.SetHeight(10.0)
    # cylinder.SetRadius(5.0)
    # cylinder.CappingOff() # Importante: sem as tampas para ter bordas
    # cylinder.Update()
    
    # writer = vtk.vtkPolyDataWriter()
    # writer.SetFileName(INPUT_VTK_FILE)
    # writer.SetInputData(cylinder.GetOutput())
    # writer.Write()
    # print("Arquivo de teste criado.")

    # Executa a função principal
    extract_border_and_gradient(INPUT_VTK_FILE, OUTPUT_CSV_FILE)
