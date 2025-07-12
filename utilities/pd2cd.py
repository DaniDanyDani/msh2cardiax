import meshio
from collections import Counter
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Adiciona parcellation dominante por célula em uma malha VTU.")
parser.add_argument("input_vtu", help="Arquivo .vtu de entrada com point_data['parcellation']")
parser.add_argument("output_vtu", help="Arquivo .vtu de saída com cell_data['parcellation_cells']")

args = parser.parse_args()

mesh = meshio.read(args.input_vtu)

points = mesh.points
parcellation = mesh.point_data["parcellation"]

cells = mesh.cells[0]
cells_connectivity = cells.data

name_parcellation_cells = "parcellation_cells"
parcellation_cells = []

for cell_id, cell_points in enumerate(cells_connectivity):
    parcellation_list = [int(parcellation[pid]) for pid in cell_points]

    contagem = Counter(parcellation_list)
    mais_frequente, freq = contagem.most_common(1)[0]
    parcellation_cells.append(mais_frequente)

    if cell_id == 0:
        print("Formato das conectividades das células")
        print(f"{cell_id=}")
        print(f"    cells_connectivity[{cell_id}] = {cell_points}")
        for pid in cell_points:
            print(f"    parcellation[{pid}] = {int(parcellation[pid])}")
        print("    Frequência de cada número:", contagem)
        print("    Números mais frequentes (ordem):", contagem.most_common())

mesh.cell_data = {
    name_parcellation_cells: [np.array(parcellation_cells)]
}

mesh.write(args.output_vtu, file_format="vtu", binary=False)

