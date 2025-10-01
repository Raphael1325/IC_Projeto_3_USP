# Simulação e Modelagem de Dados: Desenvolvimento de um Algoritmo de Alta Performance em Python

Este projeto demonstra minha capacidade de desenvolver ferramentas de análise e simulação, traduzindo problemas complexos em soluções computacionais eficientes e escaláveis. Focado no desafio de analisar a estrutura de nanopartículas, este trabalho se concentra na modelagem de dados, otimização de algoritmos e democratização do acesso à informação.

## Introdução

A intensidade de espalhamento de raios X é uma poderosa ferramenta para a análise de nanopartículas. O problema, entretanto, reside em obter a estrutura 3D de uma partícula a partir dos dados de espalhamento. Para contornar essa dificuldade, o projeto foca no problema inverso: 


simular o perfil de espalhamento de geometrias complexas para auxiliar na interpretação de dados experimentais.

O trabalho resultou no desenvolvimento de um algoritmo robusto em Python que utiliza o 

Método de Elementos Finitos (FEM) para simular o comportamento de espalhamento de partículas com formas complexas, superando as limitações dos modelos analíticos tradicionais.



## Metodologia

O projeto consistiu em três etapas principais:


Modelagem de Dados 3D: Implementei um script em Python para construir modelos 3D de partículas a partir de combinações de geometrias simples, como esferas e cilindros. A metodologia de "Grade e Verificação" foi utilizada para a construção de sólidos, enquanto a "Discretização por Comprimento de Arco" foi aplicada para a modelagem de cascas ocas, garantindo a distribuição uniforme de subunidades.






Otimização de Algoritmos: A lógica do programa foi reescrita de uma linguagem legada (Fortran) para Python , utilizando a biblioteca 


Numba para acelerar drasticamente os cálculos intensivos de histogramas e intensidade de espalhamento. 


Essa otimização resultou em uma redução de 70% no tempo de processamento em comparação com a versão original em Fortran, demonstrando proficiência em otimização de performance.


Validação e Comunicação de Resultados: Os resultados do algoritmo foram validados através da comparação com um código de referência em Fortran. Além disso, iniciei o desenvolvimento de um 

protótipo de interface web para democratizar o uso da ferramenta, permitindo que usuários sem conhecimento de programação possam criar modelos e visualizar simulações diretamente no navegador.

## Resultados


Flexibilidade e Escalabilidade: A abordagem por elementos finitos permitiu a modelagem de geometrias complexas e realistas, como o modelo de um haltere e até mesmo o dróide R2-D2. Essa flexibilidade é crucial para lidar com a complexidade do mundo real.


Inovação e Acessibilidade: A criação de um protótipo de interface web representa um passo importante para transformar uma ferramenta complexa e técnica em um recurso acessível, demonstrando a habilidade de traduzir resultados técnicos para um público mais amplo e não-especialista.


## Tecnologias e Ferramentas

Programação e Otimização: Python, Numba, Numpy, Matplotlib, Scipy, Fortran 77.

Modelagem e Análise de Dados: Método de Elementos Finitos (FEM).

Documentação e Prototipagem: LaTeX, GoogleColaboratory.

Autor: Raphael Lima Alves
Orientador(a): Prof. Dr. Cristiano Luis Pinto Oliveira
Universidade: Universidade de São Paulo (USP)
Período: 2023 - 2024
