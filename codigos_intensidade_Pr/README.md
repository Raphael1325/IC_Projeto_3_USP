# 📐 Simulação de Intensidade de Espalhamento e Função P(r)  

Este código implementa um algoritmo em **Python (com suporte a Numba para aceleração)** que permite calcular propriedades estruturais de modelos atômicos, a partir de arquivos no formato **PDB**.  

O programa realiza:  
- Leitura de coordenadas atômicas.  
- Construção de histogramas de distâncias entre pares.  
- Cálculo do **Raio de Giro (Rg)** e da **Dimensão Máxima (Dmax)**.  
- Cálculo da **intensidade de espalhamento I(q)**.  
- Cálculo e plotagem da **função de distribuição de distâncias P(r)**.  
- Geração de gráficos e arquivos de saída com os resultados.  

---

## ⚙️ Funcionalidades  

1. **Leitura de arquivo PDB**  
   - Extrai coordenadas (x, y, z) e pesos atômicos.  

2. **Cálculo de propriedades estruturais**  
   - Raio de Giro (**Rg**).  
   - Dimensão máxima (**Dmax**).  

3. **Histograma de distâncias**  
   - Geração do histograma ponderado de todas as distâncias atômicas.  

4. **Cálculo da intensidade I(q)**  
   - Implementação eficiente com **histograma de distâncias**.  
   - Normalização automática da curva de intensidade.  

5. **Função P(r)**  
   - Implementação baseada nas fórmulas de **Glatter**.  
   - Cálculo analítico e plotagem da curva P(r).  

6. **Visualização**  
   - Gráficos de:  
     - Intensidade normalizada `I(q)`  
     - Função de distribuição `P(r)`  
     - Histograma de distâncias  

7. **Exportação dos resultados**  
   - Intensidade → `intensidade.dat`  
   - Função P(r) → `pr_function_analytical.dat`  

---

## 📊 Exemplo de Saída  

- **Gráficos**  
  - `I(q) vs q` (em escala log-log).  
  - `P(r) vs r`.  
  - Histograma de distâncias atômicas.  

- **Arquivos**  
  - `intensidade.dat` → Intensidade calculada.  
  - `pr_function_analytical.dat` → Função P(r) normalizada.  

---

## 🛠️ Tecnologias Utilizadas  

- **Python 3**  
- **NumPy** – operações numéricas.  
- **Matplotlib** – visualização.  
- **Numba** – otimização com JIT compilation.  

---

## ▶️ Como Executar  

1. Instale as dependências:  

  - pip install numpy matplotlib numba
  - Edite no código o caminho para seu arquivo .pdb:
  - path_pdb = r"C:\Interface-python\teste_criacao_modelo_py\Model.pdb"
  - Execute o script
