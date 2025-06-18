# Fibrosis Layer Modifier

Este projeto fornece um script Python para expandir ou reduzir regiões de fibrose em um modelo cardíaco tridimensional no formato MSH. A fibrose é identificada por marcações presentes na malha volumétrica.

## Funcionalidades

- Expansão da fibrose por um número definido de camadas celulares adjacentes.
- Redução da fibrose, revertendo células de fibrose para tecido saudável se estiverem adjacentes a regiões não fibróticas.
- Exportação dos resultados nos formatos `.vtk`, `.msh` (compatível com Gmsh) e `.xml` (compatível com o Cardiax).
- É possível adaptar o formato do xml editando o arquivo `demo-biv.py` em `./mesh_generator/biv-gmsh2cardiax/demo-biv.py`

## Estrutura de Marcações

| Região    | ID  |
| --------- | --- |
| `base`    | 10  |
| `ve`      | 30  |
| `vd`      | 20  |
| `epi`     | 40  |
| `healthy` | 1   |
| `fibrose` | 2   |

## Uso

### Requisitos
- python 3+
- Conda/Miniconda

### Instale as dependências com:

```bash
conda env create -f Dependencies.yml
```

### Execução

Ative o ambiente conda:

```bash
conda activate fibrosis_layer_modifier
```

Execução:
```bash
python main.py -i input_file_name -o output_file_name -n N
```

- `-i`, `--input`: Caminho para o arquivo .msh de entrada contendo a malha original (sem extensão).  
- `-o`, `--output`: Nome base do arquivo de saída (sem extensão).  
- `-n`, `--n_layers`: Número de camadas:  
  - Positivo → expande a fibrose.  
  - Negativo → reduz a fibrose.  
  - Zero não é permitido.  
 
## Saídas
- output.vtk: malha modificada com o campo CellEntityIds atualizado.
- output.msh: versão .msh compatível com Gmsh com os nomes físicos das regiões.
- output.xml: versão xml do modelo

## Exemplos
Expandir a fibrose em 3 camadas:

```bash
python main.py -i Example -o Example_expanded -n 3
```

Reduzir a fibrose em 2 camadas:

```bash
python main.py -i Example -o Example_reduced -n -2
```
