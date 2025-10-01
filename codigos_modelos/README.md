# 🧩 Gerador e Manipulador de Geometrias 3D para Modelagem Estrutural

Este repositório contém um conjunto de funções e classes em **Python** para:  
- geração de coordenadas de diferentes geometrias em 3D,  
- cálculo de centro de massa,  
- visualização gráfica,  
- manipulação e exportação de modelos no formato `.dat`.  

O objetivo principal é **auxiliar na modelagem computacional de partículas/nanopartículas** em simulações físicas e biológicas, com suporte para diferentes formas geométricas e análise de dados experimentais de espalhamento.

---

## ✨ Funcionalidades

### Geometrias 3D suportadas
- Esfera sólida  
- Cilindro oco (tubo)  
- Elipsoide sólido  
- Prisma retangular  
- Casca elipsoidal (ou esférica)  
- Meia casca elipsoidal (hemisfério)  
- Círculo (anel 2D)  

### Processamento
- Geração de coordenadas em grade (`d_grid`) com pesos associados  
- Cálculo do centro de massa (CM)  
- Exportação para arquivos `.dat`  
- Visualização 3D colorida por peso  

### Análise de espalhamento
- Plotagem de intensidade `I(q)` vs `q`  
- Função de distribuição de distâncias `P(r)`  
- Normalizações com raio de giro `Rg`  
- Suporte para arquivos `.NIQ` e `.POR`  

---

## 📦 Dependências

Instale as bibliotecas necessárias com:  

bash
pip install numpy numba scipy matplotlib


1. Geração de Geometrias

# Exemplo: gerar uma esfera de raio 10 Å centrada na origem:

```python
from geometria import GerenciadorGeometrias

ger = GerenciadorGeometrias()
ger.adicionar_geometria(
    tipo=1,                      # 1 = esfera
    centro=[0, 0, 0],
    parametros={"R": 10, "d_grid": 1.0}
)

coordenadas = ger.obter_coordenadas_totais()
print("Centro de Massa:", ger.calcular_CM(coordenadas))
