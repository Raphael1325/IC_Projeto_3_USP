🧩 Gerador e Manipulador de Geometrias 3D para Modelagem Estrutural

Este repositório contém um conjunto de funções e classes em Python para:

geração de coordenadas de diferentes geometrias em 3D,

cálculo de centro de massa,

visualização gráfica,

manipulação e exportação de modelos no formato .dat.

O objetivo principal é auxiliar na modelagem computacional de partículas/nanopartículas em simulações físicas e biológicas, com suporte para diferentes formas geométricas e análise de dados experimentais de espalhamento.

✨ Funcionalidades
Geometrias 3D suportadas

Esfera sólida

Cilindro oco (tubo)

Elipsoide sólido

Prisma retangular

Casca elipsoidal (ou esférica)

Meia casca elipsoidal (hemisfério)

Círculo (anel 2D)

Processamento

Geração de coordenadas em grade (d_grid) com pesos associados

Cálculo do centro de massa (CM)

Exportação para arquivos .dat

Visualização 3D colorida por peso

Análise de espalhamento

Plotagem de intensidade I(q) vs q

Função de distribuição de distâncias P(r)

Normalizações com raio de giro Rg

Suporte para arquivos .NIQ e .POR

📦 Dependências

Instale as bibliotecas necessárias com:

pip install numpy numba scipy matplotlib

🚀 Como Usar
1. Geração de Geometrias

Exemplo: gerar uma esfera de raio 10 Å centrada na origem:

from geometria import GerenciadorGeometrias

ger = GerenciadorGeometrias()
ger.adicionar_geometria(
    tipo=1,                      # 1 = esfera
    centro=[0, 0, 0],
    parametros={"R": 10, "d_grid": 1.0}
)

coordenadas = ger.obter_coordenadas_totais()
print("Centro de Massa:", ger.calcular_CM(coordenadas))

2. Salvar em arquivo .dat
ger.salvar_modelo_dat("modelo.dat")


Formato do arquivo:

# Formato: X Y Z Peso
   0.0000    1.0000    2.0000    1.0000
   1.0000    0.0000    2.0000    1.0000
   ...

3. Visualização Gráfica

O código possui funções prontas para análise de dados experimentais:

plota_intensidade("dados.NIQ")
plot_POR("dados.POR")
plota_intensidade_vs_qRg("dados.NIQ", "rg.txt")
plota_POR_vs_rRg("dados.POR", "rg.txt")

📊 Exemplos de Geometrias

Esfera:

parametros = {"R": 20, "d_grid": 2.0}
ger.adicionar_geometria(tipo=1, centro=[0, 0, 0], parametros=parametros)


Cilindro:

parametros = {"raio_externo": 15, "raio_interno": 5, "Altura": 40, "d_grid": 1.0}
ger.adicionar_geometria(tipo=2, centro=[10, 0, 0], parametros=parametros)


Elipsoide:

parametros = {"a": 15, "b": 10, "c": 5, "d_grid": 1.0}
ger.adicionar_geometria(tipo=3, centro=[0, 0, 0], parametros=parametros)

📂 Estrutura do Código

Funções Numba (@njit) → geração eficiente de pontos para diferentes geometrias

Classe Geometria → representa uma geometria individual (tipo, centro, parâmetros, pesos)

Classe GerenciadorGeometrias → conjunto de geometrias, manipulação de coordenadas e exportação

Funções de Plotagem → análise gráfica de resultados experimentais

📝 Observações

O código usa Numba para aceleração JIT; a primeira execução pode ser mais lenta.

Os pesos podem ser definidos em parametros['peso'] (default = 1.0).

Arquivos .NIQ e .POR devem estar no formato padrão de softwares de SAXS.

📜 Licença

Este projeto é distribuído sob a licença MIT.
