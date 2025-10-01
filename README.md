# Simulação e Modelagem de Dados: Desenvolvimento de um Algoritmo para Estruturas de Nanopartículas  

Este projeto demonstra minha capacidade de desenvolver ferramentas de análise e simulação, traduzindo problemas complexos em soluções computacionais eficientes e escaláveis. Focado no desafio de analisar a estrutura de nanopartículas, este trabalho se concentra na modelagem de dados, na otimização de algoritmos e na liberdade criativa para a criação de modelos.  

---

## 📌 Introdução  

A intensidade de espalhamento de raios X é uma poderosa ferramenta para a análise de nanopartículas. O problema, entretanto, reside em obter a estrutura 3D de uma partícula a partir dos dados de espalhamento. Para contornar essa dificuldade, o projeto foca no problema inverso:  

**Construir a nanopartícula e simular o perfil de espalhamento de geometrias complexas para auxiliar na interpretação de dados experimentais.**  

O trabalho resultou no desenvolvimento de um algoritmo robusto em Python que utiliza o **Método dos Elementos Finitos (FEM)** para simular o comportamento de espalhamento de partículas com formas complexas, superando as limitações dos modelos analíticos tradicionais.  

---

## ⚙️ Metodologia  

O projeto consistiu em três etapas principais:  

### 1. Modelagem de Dados 3D  
- Desenvolvimento de um script em Python para construir modelos tridimensionais de partículas a partir de combinações de geometrias simples, como esferas e cilindros.  
- Uso da metodologia de *“Grade e Verificação”* para a construção de sólidos.  
- Aplicação da *“Discretização por Comprimento de Arco”* para a modelagem de cascas ocas, garantindo a distribuição uniforme de subunidades.  
- Baseado em um programa desenvolvido em **Turbo Pascal 7** pelo Prof. Dr. Cristiano Oliveira em sua dissertação de mestrado.  

### 2. Otimização de Algoritmos  
- A lógica inicial foi herdada do programa em **Turbo Pascal 7**, que só roda em arquiteturas antigas (ex.: Windows XP).  
- A proposta foi modernizar a abordagem para construção de modelos e cálculos de intensidade, utilizando bibliotecas otimizadas em Python.  
- Uso do **Numba** para acelerar cálculos intensivos de histogramas e intensidade de espalhamento.  
- A otimização resultou em grande redução do tempo de processamento, especialmente nos cálculos de intensidade e na função **P(r)** (*função de distribuição de pares de distância*), crucial para fornecer insights sobre a forma da partícula.  
- Comparado à versão em Fortran, o novo algoritmo demonstrou ganhos expressivos de performance.  

### 3. Validação e Comunicação de Resultados  
- Validação por meio de modelos teóricos.  
- Garantia de cálculos extensivos da intensidade de espalhamento e da função **P(r)**.  
- Mais detalhes disponíveis na seção de relatórios do projeto.  

---

## 📊 Resultados  

- **Flexibilidade e Escalabilidade:**  
  A abordagem por elementos finitos permitiu a modelagem de geometrias complexas e realistas, como halteres, dímeros e outras formas comuns em contextos biológicos. Essa flexibilidade é essencial para lidar com sistemas de alta complexidade.  

- **Inovação e Acessibilidade:**  
  Criação de um protótipo de interface web, representando um passo importante para transformar uma ferramenta técnica em um recurso acessível, capaz de alcançar tanto especialistas quanto não especialistas.  

---

## 🛠️ Tecnologias e Ferramentas  

- **Programação e Otimização:** Python, Numba, NumPy, Matplotlib, SciPy, Fortran 77  
- **Modelagem e Análise de Dados:** Método dos Elementos Finitos (FEM)  
- **Documentação e Prototipagem:** LaTeX, Google Colaboratory  

---

## 👨‍🔬 Informações do Projeto  

- **Autor:** Raphael Lima Alves  
- **Orientador(a):** Prof. Dr. Cristiano Luis Pinto Oliveira  
- **Universidade:** Universidade de São Paulo (USP)  
- **Período:** 2024 - 2025  
- **Projeto FAPESP:** 2024/13288-7  
# Simulação e Modelagem de Dados: Desenvolvimento de um Algoritmo para Estruturas de Nanopartículas  

Este projeto demonstra minha capacidade de desenvolver ferramentas de análise e simulação, traduzindo problemas complexos em soluções computacionais eficientes e escaláveis. Focado no desafio de analisar a estrutura de nanopartículas, este trabalho se concentra na modelagem de dados, na otimização de algoritmos e na liberdade criativa para a criação de modelos.  

---

## 📌 Introdução  

A intensidade de espalhamento de raios X é uma poderosa ferramenta para a análise de nanopartículas. O problema, entretanto, reside em obter a estrutura 3D de uma partícula a partir dos dados de espalhamento. Para contornar essa dificuldade, o projeto foca no problema inverso:  

**Construir a nanopartícula e simular o perfil de espalhamento de geometrias complexas para auxiliar na interpretação de dados experimentais.**  

O trabalho resultou no desenvolvimento de um algoritmo robusto em Python que utiliza o **Método dos Elementos Finitos (FEM)** para simular o comportamento de espalhamento de partículas com formas complexas, superando as limitações dos modelos analíticos tradicionais.  

---

## ⚙️ Metodologia  

O projeto consistiu em três etapas principais:  

### 1. Modelagem de Dados 3D  
- Desenvolvimento de um script em Python para construir modelos tridimensionais de partículas a partir de combinações de geometrias simples, como esferas e cilindros.  
- Uso da metodologia de *“Grade e Verificação”* para a construção de sólidos.  
- Aplicação da *“Discretização por Comprimento de Arco”* para a modelagem de cascas ocas, garantindo a distribuição uniforme de subunidades.  
- Baseado em um programa desenvolvido em **Turbo Pascal 7** pelo Prof. Dr. Cristiano Oliveira em sua dissertação de mestrado.  

### 2. Otimização de Algoritmos  
- A lógica inicial foi herdada do programa em **Turbo Pascal 7**, que só roda em arquiteturas antigas (ex.: Windows XP).  
- A proposta foi modernizar a abordagem para construção de modelos e cálculos de intensidade, utilizando bibliotecas otimizadas em Python.  
- Uso do **Numba** para acelerar cálculos intensivos de histogramas e intensidade de espalhamento.  
- A otimização resultou em grande redução do tempo de processamento, especialmente nos cálculos de intensidade e na função **P(r)** (*função de distribuição de pares de distância*), crucial para fornecer insights sobre a forma da partícula.  
- Comparado à versão em Fortran, o novo algoritmo demonstrou ganhos expressivos de performance.  

### 3. Validação e Comunicação de Resultados  
- Validação por meio de modelos teóricos.  
- Garantia de cálculos extensivos da intensidade de espalhamento e da função **P(r)**.  
- Mais detalhes disponíveis na seção de relatórios do projeto.  

---

## 📊 Resultados  

- **Flexibilidade e Escalabilidade:**  
  A abordagem por elementos finitos permitiu a modelagem de geometrias complexas e realistas, como halteres, dímeros e outras formas comuns em contextos biológicos. Essa flexibilidade é essencial para lidar com sistemas de alta complexidade.  

- **Inovação e Acessibilidade:**  
  Criação de um protótipo de interface web, representando um passo importante para transformar uma ferramenta técnica em um recurso acessível, capaz de alcançar tanto especialistas quanto não especialistas.  

---

## 🛠️ Tecnologias e Ferramentas  

- **Programação e Otimização:** Python, Numba, NumPy, Matplotlib, SciPy, Fortran 77  
- **Modelagem e Análise de Dados:** Método dos Elementos Finitos (FEM)  
- **Documentação e Prototipagem:** LaTeX, Google Colaboratory  

---

## 👨‍🔬 Informações do Projeto  

- **Autor:** Raphael Lima Alves  
- **Orientador(a):** Prof. Dr. Cristiano Luis Pinto Oliveira  
- **Universidade:** Universidade de São Paulo (USP)  
- **Período:** 2024 - 2025  
- **Projeto FAPESP:** 2024/13288-7  
