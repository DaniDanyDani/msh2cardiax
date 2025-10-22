# # Creating fibers on a simple BiV ellipsoid
#
# In this demo we create a simple BiV ellipsoid using `mshr`.
# You can install `mshr` using `conda`. It also possible to create ellipsoidal geometries using gmsh, see e.g https://github.com/finsberg/pulse/blob/0d7b5995f62f41df4eec9f5df761fa03da725f69/pulse/geometries.py#L160
#
#

import dolfin as df
import numpy as np

import ldrb

import argparse
import cardiac_geometries
import meshio


def get_physical_region(fname):
    #generate physical groups to xdmf
    import pathlib as pl
    tetra_mesh_name = pl.Path(f"mesh_{fname}.xdmf")
    auxmesh = df.Mesh()
    with df.XDMFFile(tetra_mesh_name.as_posix()) as infile:
        infile.read(auxmesh)

    tecido = df.MeshFunction("size_t", auxmesh, 3)
    with df.XDMFFile(pl.Path(tetra_mesh_name).as_posix()) as f:
        f.read(tecido, "name_to_read")

    tecido.array()[:] = tecido.array()==2 
    tecido.rename("tissue", "tissue")

    import os 

    os.system(f"rm mesh_{fname}.*")
    os.system(f"rm triangle_mesh_{fname}.*")

    return tecido

def solve_laplace(mesh, boundary_markers, boundary_values, ldrb_markers):
    V = df.FunctionSpace(mesh, 'P', 1)

    u_rv, u_lv, u_epi = boundary_values

    bc1 = df.DirichletBC(V, u_rv, boundary_markers, ldrb_markers["rv"]) 
    bc2 = df.DirichletBC(V, u_lv, boundary_markers, ldrb_markers["lv"])
    bc3 = df.DirichletBC(V, u_epi, boundary_markers, ldrb_markers["epi"])

    bcs=[bc1, bc2 ,bc3]

    ds = df.Measure('ds', domain=mesh, subdomain_data=boundary_markers)
    dx = df.Measure('dx', domain=mesh)

    # Define variational problem
    u = df.mesh_generator/biv-gmsh2cardiax/fibers_ldrb.pyTrialFunction(V)
    v = df.TestFunction(V)
    f = df.Constant(0.0)   
    a = df.dot(df.grad(u), df.grad(v))*dx  
    L = f*v*dx

    # Compute solution
    u = df.Function(V)
    df.solve(a == L, u, bcs, solver_parameters=dict(linear_solver='gmres', preconditioner='hypre_amg')) 

    return u


#mesh, ffun, markers = ldrb.gmsh2dolfin("malha3.msh")

#for m in markers:
#    print(m)

# Update the markers which are stored within the mesh


def create_fiber(meshname, biv):

    mesh, markers, markers_functions = cardiac_geometries.gmsh2dolfin(meshname)
    # msh_mesh = meshio.read(meshname)

    # tetra_cells = msh_mesh.get_cells_type("tetra")
    # triangle_cells = msh_mesh.get_cells_type("triangle")

    # tetra_data = msh_mesh.get_cell_data("gmsh:physical", "tetra")
    # triangle_data = msh_mesh.get_cell_data("gmsh:physical", "triangle")

    # tetra_mesh = meshio.Mesh(points=msh_mesh.points, cells=[("tetra", tetra_cells)],
    #                          cell_data={"name_to_read": [tetra_data]})

    # triangle_mesh = meshio.Mesh(points=msh_mesh.points, cells=[("triangle", triangle_cells)],
    #                             cell_data={"name_to_read": [triangle_data]})

    # xml_base = "teste_meshio_diretocodigo"
    # meshio.write(f"{xml_base}.xml", tetra_mesh)
    # meshio.write(f"{xml_base}_physical_region.xml", tetra_mesh)
    # meshio.write(f"{xml_base}_facet_region.xml", triangle_mesh)

    # mesh = df.Mesh()

    # with XDMFFile(f"{xml_base}.xml") as infile: # Usar XDMFFile é mais moderno e robusto que o leitor XML direto
    #     infile.read(mesh)
    # with XDMFFile(f"{xml_base}_physical_region.xml") as infile: # Usar XDMFFile é mais moderno e robusto que o leitor XML direto
    #     infile.read(mesh)
    # with XDMFFile(f"{xml_base}_facet_region.xml") as infile: # Usar XDMFFile é mais moderno e robusto que o leitor XML direto
    #     infile.read(mesh)

    # MIOCARDIO_TAG = 1
    # mesh = SubMesh(mesh, volume_markers, MIOCARDIO_TAG)
    # mesh.init(2,3)
    # parent_cell_indices = mesh.data().array('parent_cell_indices', 3)
    # parent_facet_indices = np.array([], dtype='uintp')

    # for c in df.cells(mesh):
    #     parent_cell = df.Cell(Mesh(f"{xml_base}.xml"), parent_cell_indices[c.index()])
    #     for f in df.facets(parent_cell):
    #         if f.exterior():
    #             parent_facet_indices = np.append(parent_facet_indices,f.index())

    # ffun = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
    # ffun.array()[:] = facet_markers.array()[parent_facet_indices][mesh.data().array('parent_facet_indices', 2)]

    # class MarkersFunctions:
    #     pass

    # markers_functions = MarkersFunctions()
    # markers_functions.ffun = ffun 

    print(markers)
    print(np.unique(markers_functions.ffun.array()[:]))


    ldrb_markers = {
        "base": 10,
        "rv": 20,
        "lv": 30,
        "epi": 40
    }

    mesh_markers = {
        "baseLV": 60,
        "baseRV": 50
    }

    # with df.XDMFFile(mesh.mpi_comm(), "teste.xdmf") as xdmf:
    #     xdmf.parameters.update(
    #     {
    #         "functions_share_mesh": True,
    #         "rewrite_function_mesh": False
    #     })
    #     xdmf.write(markers_functions.ffun)


    # markers_functions.ffun.array()[markers_functions.ffun.array() == mesh_markers["baseLV"]] = ldrb_markers["base"]
    # markers_functions.ffun.array()[markers_functions.ffun.array() == mesh_markers["baseRV"]] = ldrb_markers["base"]

    # Choose space for the fiber fields
    # This is a string on the form {family}_{degree}
    fiber_space = "DG_0"

    fiber, sheet, sheet_normal = None, None, None

    # Compute the microstructure
    print('Calculando fiber, sheet, sheet_normal...')
    if biv:
        fiber, sheet, sheet_normal = ldrb.dolfin_ldrb(
            mesh=mesh,
            fiber_space=fiber_space,
            ffun=markers_functions.ffun,
            markers=ldrb_markers,
            alpha_endo_lv=60,#30,  # Fiber angle on the LV endocardium
            alpha_epi_lv=-60,#-30,  # Fiber angle on the LV epicardium
            beta_endo_lv=0,  # Sheet angle on the LV endocardium
            beta_epi_lv=0,  # Sheet angle on the LV epicardium
            alpha_endo_sept=60,#60,  # Fiber angle on the Septum endocardium
            alpha_epi_sept=-60,#-60,  # Fiber angle on the Septum epicardium
            beta_endo_sept=0,  # Sheet angle on the Septum endocardium
            beta_epi_sept=0,  # Sheet angle on the Septum epicardium
            alpha_endo_rv=60,#80,  # Fiber angle on the RV endocardium
            alpha_epi_rv=-60,#-80,  # Fiber angle on the RV epicardium
            beta_endo_rv=0,  # Sheet angle on the RV endocardium
            beta_epi_rv=0,
        )

    else:
        fiber, sheet, sheet_normal = ldrb.dolfin_ldrb(
            mesh=mesh,
            fiber_space=fiber_space,
            ffun=markers_functions.ffun,
            markers=ldrb_markers,
            alpha_endo_lv=60,  # Fiber angle on the LV endocardium
            alpha_epi_lv=-60,  # Fiber angle on the LV epicardium
            beta_endo_lv=0,  # Sheet angle on the LV endocardium
            beta_epi_lv=-0,  # Sheet angle on the LV epicardium
        )

    print('Calculando long, circ e rad...')
    #long, circ and rad directions
    lv_ffun = df.MeshFunction("size_t", mesh, 2)
    lv_ffun.array()[:] = markers_functions.ffun.array().copy()
    lv_ffun.array()[markers_functions.ffun.array() == ldrb_markers["rv"]] = ldrb_markers["lv"]
    lv_markers = ldrb_markers.copy()
    lv_markers.pop("rv")

    longi, _, _ = ldrb.dolfin_ldrb(
    mesh=mesh,
    fiber_space=fiber_space,
    ffun=lv_ffun,
    markers=lv_markers,
    alpha_endo_lv=-90,
    alpha_epi_lv=-90,
    beta_endo_lv=0,
    beta_epi_lv=0,)

    circ, _, _ = ldrb.dolfin_ldrb(
    mesh=mesh,
    fiber_space=fiber_space,
    ffun=markers_functions.ffun,
    markers=ldrb_markers,
    alpha_endo_lv=0,
    alpha_epi_lv=0,
    beta_endo_lv=0,
    beta_epi_lv=0,)

    rad = df.project(df.cross(sheet, fiber), sheet.function_space())


    fiber.rename("f_0","f_0")
    sheet.rename("s_0","s_0")
    sheet_normal.rename("n_0","n_0")
    longi.rename("long","long")
    circ.rename("circ","circ")
    rad.rename("rad","rad")


    print("Salvando...")

    with df.XDMFFile(mesh.mpi_comm(), meshname + ".xdmf") as xdmf:
        xdmf.parameters.update(
        {
            "functions_share_mesh": True,
            "rewrite_function_mesh": False
        })
        xdmf.write(mesh)
        xdmf.write(fiber, 0)
        xdmf.write(sheet, 0)
        xdmf.write(sheet_normal,0)
        xdmf.write(longi, 0)
        xdmf.write(circ, 0)
        xdmf.write(rad,0)

    #print("Done.")

    f0 = fiber.vector().get_local()
    f0 = f0.reshape(int(len(f0)/3), 3)
    s0 = sheet.vector().get_local()
    s0 = s0.reshape(int(len(s0)/3), 3)
    n0 = sheet_normal.vector().get_local()
    n0 = n0.reshape(int(len(n0)/3), 3)
    l0 = longi.vector().get_local()
    l0 = l0.reshape(int(len(l0)/3), 3)
    c0 = circ.vector().get_local()
    c0 = c0.reshape(int(len(c0)/3), 3)
    r0 = rad.vector().get_local()
    r0 = r0.reshape(int(len(r0)/3), 3)

    return f0, s0, n0, l0, c0, r0


#fib = create_fiber('malha-fibrose-grossa.msh')

#fib = fib.reshape(int(len(fib)/3), 3)

#for i, f in enumerate(fib):
#    print(i, f)
#    print('')


# Store the results
#with df.HDF5File(mesh.mpi_comm(), "biv.h5", "w") as h5file:
#    h5file.write(fiber, "/fiber")
#    h5file.write(sheet, "/sheet")
#    h5file.write(sheet_normal, "/sheet_normal")


# You can also store files in XDMF which will also compute the fiber angle as scalars on the glyph to be visualised in Paraview. Note that these functions don't work (yet) using mpirun

# (These function are not tested in parallel)
#ldrb.fiber_to_xdmf(fiber, "biv_fiber_new1")
#ldrb.fiber_to_xdmf(sheet, "biv_sheet")
#ldrb.fiber_to_xdmf(sheet_normal, "biv_sheet_normal")

# ![_](_static/figures/biv_fiber.png)
# [Link to source code](https://github.com/finsberg/ldrb/blob/main/demos/demo_biv.py)
