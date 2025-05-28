import gmsh
import sys

def read_model(model):

    if len(model) > 1 and model[1][0] != '-':
        gmsh.open(model[1])
    else:
        gmsh.model.occ.addCone(1, 0, 0, 1, 0, 0, 0.5, 0.1)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate()

    print('\nModel ' + gmsh.model.getCurrent() + '.msh (' +
      str(gmsh.model.getDimension()) + 'D)')


    entities = gmsh.model.getEntities()
    # print(len(entities))

    for e in entities:

        dim = e[0]
        tag = e[1]

        # Pega os nós da entidade (dim,tag):
        nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(dim, tag)

        # Pega os elementos da entidade (dim,tag):
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)

        # Tipo e nome da entidade:
        type = gmsh.model.getType(dim, tag)
        name = gmsh.model.getEntityName(dim, tag)
        if len(name): name += ' '
        # print("\nEntity " + name + str(e) + " of type " + type)

        # # Número de nós e elementos
        # numElem = sum(len(i) for i in elemTags)
        # print(" - Mesh has " + str(len(nodeTags)) + " nodes and " + str(numElem) +
        #     " elements")

        # Physical group que a entidade faz parte.
        physicalTags = gmsh.model.getPhysicalGroupsForEntity(dim, tag)
        if len(physicalTags):
            s = ''
            for p in physicalTags:
                n = gmsh.model.getPhysicalName(dim, p)
                if n: n += ''
                s += n + '(' + str(dim) + ', ' + str(p) + ') '
            # print(" - Physical groups: " + s)

        # # Imprime o tipo do elemento
        # for t in elemTypes:
        #     name, dim, order, numv, parv, _ = gmsh.model.mesh.getElementProperties(
        #         t)
        #     print(" - Element type: " + name + ", order " + str(order) + "\n")
        

        # print(f"{n=}")
        if n.strip() == "epi":
            epi = e

        elif n.strip() == "ve":
            # print(f"{e=}")
            ve = e

        elif n.strip() == "vd":
            # print(f"{e=}")
            vd = e
        
        elif n.strip() == "base":
            # print(f"{e=}")
            base = e

        elif n.strip() == "healthy":
            # print(f"{e=}")
            healthy = e
        
        elif n.strip() == "fibrose":
            # print(f"{e=}")
            fibrose = e
        
    modelo_original = {
        "epi": epi,
        "ve": ve,
        "vd": vd,
        "base": base,
        "healthy": healthy,
        "fibrose": fibrose
    }

    print(f"Entidades salvas em: \n    {modelo_original=}\n")
    return modelo_original




gmsh.initialize()
model = sys.argv

modelo_original = read_model(model)

elemenTypes_fibrose, elemenTags_fibrose, elemenNodeTags_fibrose = gmsh.model.mesh.getElements(modelo_original["fibrose"][0], modelo_original["fibrose"][1])
elemenTypes_healthy, elemenTags_healthy, elemenNodeTags_healthy = gmsh.model.mesh.getElements(modelo_original["healthy"][0], modelo_original["healthy"][1])

# print(type(elemenNodeTags_fibrose[0]),elemenNodeTags_fibrose[0],"\n", elemenTags_fibrose[0])


connectado = []
if elemenTypes_fibrose == 4:
    
    i = 0
    for tag_fibrose in elemenTags_fibrose:
        nodes_fibroses = elemenNodeTags_fibrose[i:i+4]    
        # print(nodes_fibroses)

        j = 0
        for tag_healthy in elemenTags_healthy:
            nodes_healthy = elemenNodeTags_healthy[j:j+4]   

            if j%2==0:
                print(f"Tentando {tag_fibrose=} e {tag_healthy}") 

            for nodes_fibroses in nodes_fibroses:
                # print(f'{nodes_healthy=}')
                for node_fibroses in nodes_fibroses:
                    # print(f'{node_fibroses=}')
                    for node_healthy in nodes_healthy:
                        # print(f'{nodes_healthy=}')
                        if node_fibroses == node_healthy.any():
                            connectado.append((tag_fibrose, tag_healthy))
                            i += 4
                            j += 4
                            print(f"conectividade entre {tag_fibrose=} e {tag_healthy=}")

            



gmsh.clear()









# ============================ Final ============================ 
# Visualização
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

# Limpar data
gmsh.clear()

# Finalizar
gmsh.finalize()