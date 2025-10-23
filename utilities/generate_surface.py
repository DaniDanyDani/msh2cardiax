# -*- coding: utf-8 -*-
"""
Este script lê uma nuvem de pontos de um arquivo VTK
(como o gerado no script anterior) e usa o algoritmo
vtkDelaunay3D + vtkGeometryFilter para reconstruir uma
malha de superfície 3D (triangulada) a partir desses pontos.
"""
import vtk
import numpy as np
from vtk.util import numpy_support as VN

def create_test_point_cloud_vtk(filename, num_points=1000, radius=10):
    """
    Cria um arquivo VTK de nuvem de pontos de teste (uma esfera).
    Isso nos dá um arquivo de entrada para que o script possa ser
    executado e testado.
    """
    print(f"Criando arquivo de teste: {filename}")
    
    # Gera pontos aleatórios na superfície de uma esfera
    # Sorteia 'phi' (ângulo azimutal) e 'costheta' (cosseno do ângulo polar)
    phi = np.random.uniform(0, 2 * np.pi, num_points)
    costheta = np.random.uniform(-1, 1, num_points)
    
    theta = np.arccos(costheta)
    
    # Converte de coordenadas esféricas para cartesianas
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)
    
    # Combina em um array NumPy (N_pontos x 3)
    points_np = np.vstack((x, y, z)).T
    
    # --- Converte o array NumPy para um vtkPolyData ---
    
    # 1. Cria o vtkPoints
    vtk_points = vtk.vtkPoints()
    vtk_points.SetData(VN.numpy_to_vtk(points_np, deep=True))

    # 2. Cria o vtkPolyData
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)

    # 3. Adiciona um vértice (célula) para cada ponto
    verts = vtk.vtkCellArray()
    for i in range(polydata.GetNumberOfPoints()):
        verts.InsertNextCell(1)
        verts.InsertCellPoint(i)
    polydata.SetVerts(verts)

    # 4. Salva o arquivo VTK
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(polydata)
    writer.Write()
    print("Arquivo de teste criado com sucesso.")


def reconstruct_surface(input_point_cloud_vtk, output_surface_mesh_vtk):
    """
    Lê uma nuvem de pontos VTK e gera uma malha de superfície.

    Args:
        input_point_cloud_vtk (str): Caminho para o .vtk da nuvem de pontos.
        output_surface_mesh_vtk (str): Caminho para salvar o .vtk da malha.
    """
    
    # --- 1. Ler a Nuvem de Pontos ---
    print(f"Lendo a nuvem de pontos de: {input_point_cloud_vtk}")
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(input_point_cloud_vtk)
    reader.Update()
    
    point_cloud_polydata = reader.GetOutput()
    
    if point_cloud_polydata is None or point_cloud_polydata.GetNumberOfPoints() == 0:
        print(f"Erro: Não foi possível ler a nuvem de pontos do arquivo {input_point_cloud_vtk}")
        return

    print(f"Nuvem de pontos lida com {point_cloud_polydata.GetNumberOfPoints()} pontos.")

    # --- 2. Reconstruir a Superfície com vtkDelaunay3D ---
    # vtkPowerCrust não está na instalação padrão. Usamos Delaunay3D + GeometryFilter.
    print("Iniciando a reconstrução da superfície (vtkDelaunay2D)...")
    
    # vtkDelaunay3D cria uma malha de volume (tetraedros) a partir dos pontos
    delaunay = vtk.vtkDelaunay2D()
    delaunay.SetInputData(point_cloud_polydata)
    delaunay.Update()
    
    # vtkGeometryFilter extrai a superfície externa ("pele") da malha de volume
    # print("Extraindo superfície externa da malha de volume...")
    # geometry_filter = vtk.vtkGeometryFilter()
    # geometry_filter.SetInputConnection(delaunay.GetOutputPort())
    # geometry_filter.Update()
    
    print("Reconstrução concluída.")

    # --- 3. Salvar a Malha de Superfície ---
    # O output do geometry_filter é um vtkPolyData com triângulos
    # surface_mesh = geometry_filter.GetOutput()
    surface_mesh = delaunay.GetOutput()
    print(f"Number of triangles: {surface_mesh.GetNumberOfCells()}")
    

    print(f"Salvando malha de superfície em: {output_surface_mesh_vtk}")
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(output_surface_mesh_vtk)
    writer.SetInputData(surface_mesh)
    writer.Write()
    
    print("--- Processo Finalizado ---")
    print(f"Malha final salva. Você pode abri-la no ParaView: {output_surface_mesh_vtk}")


if __name__ == '__main__':
    
    # --- CONFIGURAÇÃO ---
    # Arquivo de entrada: A nuvem de pontos que você gerou
    print("\nGerando csv da borda para o lv")
    INPUT_VTK_POINTS = 'lv_cap.vtk'
    OUTPUT_VTK_MESH = 'lv_cap_surface.vtk'
    reconstruct_surface(INPUT_VTK_POINTS, OUTPUT_VTK_MESH)
    
    print("\nGerando csv da borda para o rv")
    INPUT_VTK_POINTS = 'rv_cap.vtk'
    OUTPUT_VTK_MESH = 'rv_cap_surface.vtk'
    reconstruct_surface(INPUT_VTK_POINTS, OUTPUT_VTK_MESH)
    
    print("\nGerando csv da borda para o epi")
    INPUT_VTK_POINTS = 'epi_cap.vtk'
    OUTPUT_VTK_MESH = 'epi_cap_surface.vtk'
    reconstruct_surface(INPUT_VTK_POINTS, OUTPUT_VTK_MESH)
