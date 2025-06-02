from pathlib import Path
from vtkmodules.vtkIOLegacy import vtkPolyDataReader
from vtkmodules.vtkIOPLY import vtkPLYReader
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader

from vtkmodules.vtkCommonCore import (
    vtkIdList
)

from tqdm import tqdm

import vtk
import os
import sys
import meshio

def ReadUnstructuredGrid(file_name):
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()
    return reader.GetOutput()

def SaveUnstructuredGrid(mesh, file_name):
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(file_name)
    writer.SetInputData(mesh)
    writer.Write()

markers = {
    "epi":40,
    "ve":20,
    "vd":30,
    "base":10,
    "healthy":1,
    "fibrose":2
}


input_file = "Patient_7_mesh.vtk"
output_file = "Patient_7_mesh_editado.vtk"
n_camadas = -3
aumentar = False
diminuir = False

if n_camadas>0:
    aumentar = True
    print(f"{n_camadas=}: Aumentando fibrose.")
elif n_camadas<0:
    print(f"{n_camadas=}: Diminuindo fibrose.")
    diminuir = True
else:
    sys.exit(1)

# Malha de entrada
input_mesh = ReadUnstructuredGrid(input_file)
if not input_mesh.GetCellData().HasArray("CellEntityIds"):
    print("ERRO: 'CellEntityIds' não encontrado em input_mesh.")
    sys.exit(1)

cell_array_input = input_mesh.GetCellData().GetArray("CellEntityIds")
if cell_array_input.GetNumberOfTuples() != input_mesh.GetNumberOfCells():
    print("ERRO: Tamanho do array 'CellEntityIds' não corresponde ao número de células.")
    sys.exit(1)

valid_values = set(markers.values())
for i in range(cell_array_input.GetNumberOfTuples()):
    val = cell_array_input.GetValue(i)
    if val not in valid_values:
        print(f"Aviso: Célula {i} tem valor inesperado: {val}")



# Cópia inicial da malha de entrada
output_mesh = vtk.vtkUnstructuredGrid()
output_mesh.DeepCopy(input_mesh)

output_cell_array = vtk.vtkIntArray()
output_cell_array.DeepCopy(input_mesh.GetCellData().GetArray("CellEntityIds"))
output_cell_array.SetName("CellEntityIds")
output_mesh.GetCellData().SetScalars(output_cell_array)


cell_array = output_mesh.GetCellData().GetArray("CellEntityIds")
if cell_array is None:
    print("ERRO: 'CellEntityIds' não encontrado em output_mesh.")
    sys.exit(1)


# Número total de células
qntCells = input_mesh.GetNumberOfCells()

for camada in range(abs(n_camadas)):
    print(f" Iteração {camada+1}/{abs(n_camadas)}:")

    temp_mesh = vtk.vtkUnstructuredGrid()
    temp_mesh.DeepCopy(output_mesh)

    for idCell in tqdm(range(qntCells), desc="  Processando células", unit="célula"):

        cell_data = output_mesh.GetCellData().GetArray("CellEntityIds").GetValue(idCell)
        
        if cell_data == markers["fibrose"]:

            cell_point_ids = vtk.vtkIdList()
            input_mesh.GetCellPoints(idCell, cell_point_ids)
            
            neighborCellIds = vtkIdList()
            visited_neighbors = set()

            for i in range(cell_point_ids.GetNumberOfIds()):
                point_id_list = vtk.vtkIdList()
                point_id_list.InsertNextId(cell_point_ids.GetId(i))

                local_neighbors = vtk.vtkIdList()
                input_mesh.GetCellNeighbors(idCell, point_id_list, local_neighbors)

                for j in range(local_neighbors.GetNumberOfIds()):
                    neighbor_id = local_neighbors.GetId(j)
                    if neighbor_id != idCell:
                        visited_neighbors.add(neighbor_id)


            for neighbor_id in visited_neighbors:
 
                if 0 <= neighbor_id < qntCells:
                    neighbor_data = output_mesh.GetCellData().GetArray("CellEntityIds").GetValue(neighbor_id)
                else:
                    print("ID INVÁLIDO APRA NEIGHBOR")
                    sys.exit(1)
                
                if aumentar and neighbor_data == markers["healthy"]:
                    temp_mesh.GetCellData().GetArray("CellEntityIds").SetValue(neighbor_id, markers["fibrose"])

                elif diminuir:
                    vizinho_nao_fibrose = any(
                        output_mesh.GetCellData().GetArray("CellEntityIds").GetValue(neighbor_id) != markers["fibrose"]
                        for neighbor_id in visited_neighbors
                    )
                    if vizinho_nao_fibrose:
                        temp_mesh.GetCellData().GetArray("CellEntityIds").SetValue(idCell, markers["healthy"])


    output_mesh.DeepCopy(temp_mesh)

SaveUnstructuredGrid(output_mesh, output_file)
meshio.write("Patient_7_mesh_editado.msh", output_mesh)
print(f"Malha salva em: {output_file}")