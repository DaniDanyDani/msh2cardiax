# Fibrosis Layer Modifier

Este projeto fornece um script Python para expandir ou reduzir regiões de fibrose em um modelo cardíaco tridimensional no formato VTK. A fibrose é identificada por marcações (`CellEntityIds`) presentes na malha.

## Funcionalidades

- Expansão da fibrose por um número definido de camadas celulares adjacentes.
- Redução da fibrose, revertendo células de fibrose para tecido saudável se estiverem adjacentes a regiões não fibróticas.
- Exportação dos resultados nos formatos `.vtk` e `.msh` (compatível com Gmsh).

## Estrutura de Marcações

| Região    | ID  |
| --------- | --- |
| `base`    | 10  |
| `ve`      | 20  |
| `vd`      | 30  |
| `epi`     | 40  |
| `healthy` | 1   |
| `fibrose` | 2   |

## Uso

### Requisitos

- Python 3.8+
- Bibliotecas:
  - `vtk`
  - `tqdm`
  - `argparse`
  - `meshio`

Instale as dependências com:

```bash
pip install vtk tqdm meshio
```

### Execução

```bash
python main.py -i input.vtk -o output -n N
```

- `-i`, `--input`: Caminho para o arquivo .vtk de entrada contendo a malha original.  
- `-o`, `--output`: Nome base do arquivo de saída (sem extensão).  
- `-n`, `--n_layers`: Número de camadas:  
  - Positivo → expande a fibrose.  
  - Negativo → reduz a fibrose.  
  - Zero não é permitido.  

 
## Saídas
- output.vtk: malha modificada com o campo CellEntityIds atualizado.
- output.msh: versão .msh compatível com Gmsh com os nomes físicos das regiões.

## Exemplos
Expandir a fibrose em 3 camadas:

```bash
python main.py -i Example.vtk -o heart_expanded -n 3
```

Reduzir a fibrose em 2 camadas:

```bash
python main.py -i Example.vtk -o heart_reduced -n -2
```
