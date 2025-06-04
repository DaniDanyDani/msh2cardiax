import argparse
from tqdm import tqdm
import meshio

from vtkmodules.vtkCommonCore import vtkIdList
from vtkmodules.vtkIOLegacy import vtkPolyDataReader
from vtkmodules.vtkIOPLY import vtkPLYReader
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridWriter

import vtk
import sys

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

def write_vtk(mesh, output_file):
    SaveUnstructuredGrid(mesh, output_file+'.vtk')
    print(f"Mesh .vtk write in: {output_file}.vtk")


def write_msh_from_vtk(mesh, filename_msh, field_name="CellEntityIds"):
    with open(filename_msh+'.msh', "w") as f:
        f.write("$MeshFormat\n2.0 0 8\n$EndMeshFormat\n")

        f.write("$PhysicalNames\n6\n")
        f.write('2 10 "base"\n2 20 "ve"\n2 30 "vd"\n2 40 "epi"\n3 1 "healthy"\n3 2 "fibrose"\n')
        f.write("$EndPhysicalNames\n")

        # Nodes
        points = mesh.GetPoints()
        num_points = points.GetNumberOfPoints()
        f.write("$Nodes\n")
        f.write(f"{num_points}\n")
        for i in range(num_points):
            x, y, z = points.GetPoint(i)
            f.write(f"{i+1} {x:.16f} {y:.16f} {z:.16f}\n")
        f.write("$EndNodes\n")

        # Elements
        entity_array = mesh.GetCellData().GetArray(field_name)
        if entity_array is None:
            raise RuntimeError(f"Field '{field_name}' not found in cells.")

        f.write("$Elements\n")
        written_elements = 0
        element_lines = []

        geom_by_phys = {
            1: 1, 2: 2, 40: 1, 30: 2, 20: 3, 10: 4
        }

        for i in range(mesh.GetNumberOfCells()):
            cell = mesh.GetCell(i)
            cell_type = cell.GetCellType()
            node_ids = [cell.GetPointId(j)+1 for j in range(cell.GetNumberOfPoints())]
            phys_id = entity_array.GetValue(i)
            geom_id = geom_by_phys.get(phys_id, 1)

            if cell_type == vtk.VTK_TETRA:
                gmsh_type = 4
            elif cell_type == vtk.VTK_HEXAHEDRON:
                gmsh_type = 5
            elif cell_type == vtk.VTK_TRIANGLE:
                gmsh_type = 2
            elif cell_type == vtk.VTK_QUAD:
                gmsh_type = 3
            else:
                continue

            tag_str = f"{gmsh_type} 2 {phys_id} {geom_id}"
            node_str = " ".join(str(nid) for nid in node_ids)
            element_lines.append(f"{written_elements+1} {tag_str} {node_str}")
            written_elements += 1

        f.write(f"{written_elements}\n")
        for line in element_lines:
            f.write(f"{line}\n")
        f.write("$EndElements\n")

        print(f"Mesh .msh write in: {filename_msh}.msh")

parser = argparse.ArgumentParser(description="Expand or shrink fibrosis in a VTK heart mesh.")
parser.add_argument("-i", "--input", required=True, help="Input .vtk file")
parser.add_argument("-o", "--output", required=False, help="Output file name")
parser.add_argument("-n", "--n_layers", type=int, required=True, help="Number of layers: positive to expand, negative to shrink")

args = parser.parse_args()

input_file = args.input
output_file = args.output
n_layers = args.n_layers

expand = False
shrink = False

if n_layers > 0:
    expand = True
    print(f"{n_layers=}: Expanding fibrosis.")
elif n_layers < 0:
    shrink = True
    print(f"{n_layers=}: Shrinking fibrosis.")
else:
    sys.exit(1)


# Input mesh
input_mesh = ReadUnstructuredGrid(input_file)
if not input_mesh.GetCellData().HasArray("CellEntityIds"):
    print("ERROR: 'CellEntityIds' not found in input_mesh.")
    sys.exit(1)

cell_array_input = input_mesh.GetCellData().GetArray("CellEntityIds")
if cell_array_input.GetNumberOfTuples() != input_mesh.GetNumberOfCells():
    print("ERROR: Size of 'CellEntityIds' array does not match number of cells.")
    sys.exit(1)

valid_values = set(markers.values())
for i in range(cell_array_input.GetNumberOfTuples()):
    val = cell_array_input.GetValue(i)
    if val not in valid_values:
        print(f"Warning: Cell {i} has unexpected value: {val}")



# Initial copy of the input mesh
output_mesh = vtk.vtkUnstructuredGrid()
output_mesh.DeepCopy(input_mesh)

output_cell_array = vtk.vtkIntArray()
output_cell_array.DeepCopy(input_mesh.GetCellData().GetArray("CellEntityIds"))
output_cell_array.SetName("CellEntityIds")
output_mesh.GetCellData().SetScalars(output_cell_array)


cell_array = output_mesh.GetCellData().GetArray("CellEntityIds")
if cell_array is None:
    print("ERROR: 'CellEntityIds' not found in output_mesh.")
    sys.exit(1)


# Total number of cells
total_cells = input_mesh.GetNumberOfCells()

for layer in range(abs(n_layers)):
    print(f" Iteration {layer+1}/{abs(n_layers)}:")

    temp_mesh = vtk.vtkUnstructuredGrid()
    temp_mesh.DeepCopy(output_mesh)
    has_fibrosis = False
    has_healthy = False

    for cell_id in tqdm(range(total_cells), desc="  Processing cells", ascii=True):
    # for cell_id in range(total_cells):

        cell_data = output_mesh.GetCellData().GetArray("CellEntityIds").GetValue(cell_id)
        
        
        if cell_data == markers["healthy"]:
            has_healthy = True
        
        if cell_data == markers["fibrose"]:
            has_fibrosis = True

            cell_point_ids = vtk.vtkIdList()
            input_mesh.GetCellPoints(cell_id, cell_point_ids)
            
            visited_neighbors = set()

            for i in range(cell_point_ids.GetNumberOfIds()):
                point_id_list = vtk.vtkIdList()
                point_id_list.InsertNextId(cell_point_ids.GetId(i))

                local_neighbors = vtk.vtkIdList()
                input_mesh.GetCellNeighbors(cell_id, point_id_list, local_neighbors)

                for j in range(local_neighbors.GetNumberOfIds()):
                    neighbor_id = local_neighbors.GetId(j)
                    if neighbor_id != cell_id:
                        visited_neighbors.add(neighbor_id)

        
            for neighbor_id in visited_neighbors:
 
                if 0 <= neighbor_id < total_cells:
                    neighbor_data = output_mesh.GetCellData().GetArray("CellEntityIds").GetValue(neighbor_id)
                else:
                    print("INVALID NEIGHBOR ID")
                    sys.exit(1)
                
                if expand and neighbor_data == markers["healthy"]:
                    temp_mesh.GetCellData().GetArray("CellEntityIds").SetValue(neighbor_id, markers["fibrose"])

                elif shrink:
                    neighbor_not_fibrosis = any(
                        output_mesh.GetCellData().GetArray("CellEntityIds").GetValue(neighbor_id) != markers["fibrose"]
                        for neighbor_id in visited_neighbors
                    )
                    if neighbor_not_fibrosis:
                        temp_mesh.GetCellData().GetArray("CellEntityIds").SetValue(cell_id, markers["healthy"])

    output_mesh.DeepCopy(temp_mesh)
    if expand and not has_healthy:
        print(f"In iteration {layer + 1}, there is no more healthy tissue. Saving mesh and exiting.")
        write_vtk(output_mesh, output_file)
        write_msh_from_vtk(output_mesh, output_file)
        sys.exit(0)

    if shrink and not has_fibrosis:
        print(f"In iteration {layer + 1}, there is no more fibrosis. Saving mesh and exiting.")
        write_vtk(output_mesh, output_file)
        write_msh_from_vtk(output_mesh, output_file)
        sys.exit(0)


write_vtk(output_mesh, output_file)
write_msh_from_vtk(output_mesh, output_file)