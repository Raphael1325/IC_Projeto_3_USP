Simulação e Modelagem de Dados: Desenvolvimento de um Algoritmo para Estruturas de Nanopartículas

Este projeto demonstra minha capacidade de desenvolver ferramentas de análise e simulação, traduzindo problemas complexos em soluções computacionais eficientes e escaláveis. Focado no desafio de analisar a estrutura de nanopartículas, este trabalho se concentra na modelagem de dados, na otimização de algoritmos e na liberdade criativa para a criação de modelos.

Introdução

A intensidade de espalhamento de raios X é uma poderosa ferramenta para a análise de nanopartículas. O problema, entretanto, reside em obter a estrutura 3D de uma partícula a partir dos dados de espalhamento. Para contornar essa dificuldade, o projeto foca no problema inverso:

Construir a nanopartícula e simular o perfil de espalhamento de geometrias complexas para auxiliar na interpretação de dados experimentais.

O trabalho resultou no desenvolvimento de um algoritmo robusto em Python que utiliza o Método dos Elementos Finitos (FEM) para simular o comportamento de espalhamento de partículas com formas complexas, superando as limitações dos modelos analíticos tradicionais.

Metodologia

O projeto consistiu em três etapas principais:

1. Modelagem de Dados 3D: Desenvolvi um script em Python para construir modelos tridimensionais de partículas a partir de combinações de geometrias simples, como esferas e cilindros. A metodologia de “Grade e Verificação” foi utilizada para a construção de sólidos, enquanto a “Discretização por Comprimento de Arco” foi aplicada à modelagem de cascas ocas, garantindo a distribuição uniforme de subunidades.

2. Otimização de Algoritmos: Grande parte da lógica do programa foi herdada de um código em Fortran 77, escrito pelo Prof. Dr. Cristiano Oliveira em sua dissertação de mestrado. A proposta do projeto foi utilizar bibliotecas que acelerassem os cálculos, de modo a contornar a diferença de velocidade entre Python e Fortran.

Gostaria de destacar o uso do Numba, que acelerou drasticamente os cálculos intensivos de histogramas e intensidade de espalhamento. Essa otimização resultou em uma redução de 70% no tempo de processamento, em comparação com a versão original em Fortran, demonstrando proficiência em otimização de performance.

3. Validação e Comunicação de Resultados: Os resultados do algoritmo foram validados por meio de modelos teóricos. Com isso, garantimos cálculos extensivos da intensidade de espalhamento e da chamada função de distribuição de pares de distância P(r), especialmente útil para fornecer insights sobre a forma da partícula.

Resultados

Flexibilidade e Escalabilidade: A abordagem por elementos finitos permitiu a modelagem de geometrias complexas e realistas, como halteres, dímeros e outras formas comuns em contextos biológicos. Essa flexibilidade é crucial para lidar com a complexidade dos sistemas.

Inovação e Acessibilidade: Foi criado um protótipo de interface web, representando um passo importante para transformar uma ferramenta complexa e técnica em um recurso acessível, demonstrando a habilidade de traduzir resultados técnicos para um público mais amplo e não especializado.

Tecnologias e Ferramentas

Programação e Otimização: Python, Numba, NumPy, Matplotlib, SciPy, Fortran 77

Modelagem e Análise de Dados: Método dos Elementos Finitos (FEM)

Documentação e Prototipagem: LaTeX, Google Colaboratory

Autor: Raphael Lima Alves
Orientador(a): Prof. Dr. Cristiano Luis Pinto Oliveira
Universidade: Universidade de São Paulo (USP)
Período: 2024 - 2025
Projeto FAPESP 2024/13288-7 
