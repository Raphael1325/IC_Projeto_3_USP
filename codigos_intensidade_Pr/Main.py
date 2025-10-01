import numpy as np
import numba
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import subprocess
from scipy.integrate import quad
from scipy.interpolate import interp1d
# Adicione estas importações no topo do seu arquivo .py
from numba.typed import List
from numba.core import types

@numba.njit
def gerar_esfera(R_esfera, d_grid, peso=1.0):
    """Gera coordenadas e pesos para uma esfera sólida."""
    # --- CORREÇÃO AQUI ---
    x_vals = List.empty_list(types.float64)
    y_vals = List.empty_list(types.float64)
    z_vals = List.empty_list(types.float64)
    pesos = List.empty_list(types.float64)
    
    n_passos = int(R_esfera / d_grid) + 1
    
    for i in range(-n_passos, n_passos + 1):
        x = i * d_grid
        for j in range(-n_passos, n_passos + 1):
            y = j * d_grid
            for k in range(-n_passos, n_passos + 1):
                z = k * d_grid
                if (x**2 + y**2 + z**2) < R_esfera**2:
                    x_vals.append(x)
                    y_vals.append(y)
                    z_vals.append(z)
                    pesos.append(peso)
                    
    return x_vals, y_vals, z_vals, pesos

@numba.njit
def gerar_cilindro(raio_externo, raio_interno, altura, d_grid, peso=1.0):
    """Gera coordenadas e pesos para um tubo (cilindro oco)."""
    # --- CORREÇÃO AQUI ---
    x_vals = List.empty_list(types.float64)
    y_vals = List.empty_list(types.float64)
    z_vals = List.empty_list(types.float64)
    pesos = List.empty_list(types.float64)
    
    limite_xy = int(raio_externo / d_grid) + 1
    limite_z = int((altura / 2) / d_grid) + 1
    
    r_ext_sq = raio_externo**2
    r_int_sq = raio_interno**2
    
    for i in range(-limite_xy, limite_xy + 1):
        x = i * d_grid
        for j in range(-limite_xy, limite_xy + 1):
            y = j * d_grid
            dist_sq = x**2 + y**2
            
            if (dist_sq >= r_int_sq) and (dist_sq < r_ext_sq):
                for k in range(-limite_z, limite_z + 1):
                    z = k * d_grid
                    if abs(z) < altura / 2:
                        x_vals.append(x)
                        y_vals.append(y)
                        z_vals.append(z)
                        pesos.append(peso)
                        
    return x_vals, y_vals, z_vals, pesos

@numba.njit
def gerar_elipsoide(a, b, c, d_grid, peso=1.0):
    """Gera coordenadas e pesos para um elipsoide sólido."""
    # --- CORREÇÃO AQUI ---
    x_vals = List.empty_list(types.float64)
    y_vals = List.empty_list(types.float64)
    z_vals = List.empty_list(types.float64)
    pesos = List.empty_list(types.float64)
    
    limite_x = int(a / d_grid) + 1
    limite_y = int(b / d_grid) + 1
    limite_z = int(c / d_grid) + 1
    
    for i in range(-limite_x, limite_x + 1):
        x = i * d_grid
        for j in range(-limite_y, limite_y + 1):
            y = j * d_grid
            for k in range(-limite_z, limite_z + 1):
                z = k * d_grid
                if ((x/a)**2 + (y/b)**2 + (z/c)**2) < 1:
                    x_vals.append(x)
                    y_vals.append(y)
                    z_vals.append(z)
                    pesos.append(peso)
                        
    return x_vals, y_vals, z_vals, pesos

@numba.njit
def gerar_prisma(a, b, c, d_grid, peso=1.0):
    """Gera coordenadas e pesos para um prisma sólido."""
    # --- CORREÇÃO AQUI ---
    x_vals = List.empty_list(types.float64)
    y_vals = List.empty_list(types.float64)
    z_vals = List.empty_list(types.float64)
    pesos = List.empty_list(types.float64)
    
    limite_x = int((a / 2) / d_grid) + 1
    limite_y = int((b / 2) / d_grid) + 1
    limite_z = int((c / 2) / d_grid) + 1
    
    for i in range(-limite_x, limite_x + 1):
        x = i * d_grid
        for j in range(-limite_y, limite_y + 1):
            y = j * d_grid
            for k in range(-limite_z, limite_z + 1):
                z = k * d_grid
                if (abs(x) < a / 2) and (abs(y) < b / 2) and (abs(z) < c / 2):
                    x_vals.append(x)
                    y_vals.append(y)
                    z_vals.append(z)
                    pesos.append(peso)
                        
    return x_vals, y_vals, z_vals, pesos

@numba.njit
def _gerar_pontos_no_anel(raio, z, num_pontos):
    """Função auxiliar Numba para gerar pontos em um anel 3D."""
    x_vals = np.zeros(num_pontos, dtype=np.float64)
    y_vals = np.zeros(num_pontos, dtype=np.float64)
    z_vals = np.full(num_pontos, z, dtype=np.float64)
    
    for i in range(num_pontos):
        angulo = i * (2 * np.pi / num_pontos)
        x_vals[i] = raio * np.cos(angulo)
        y_vals[i] = raio * np.sin(angulo)
        
    return x_vals, y_vals, z_vals

def gerar_casca_elipsoidal(a, b, d_grid, peso=1.0):
    """
    Gera coordenadas e pesos para uma casca elipsoidal (ou esférica se a=b).
    A lógica é baseada na distribuição uniforme de anéis ao longo do comprimento do arco da superfície.
    """
    if a == 0 or b == 0:
        return [], [], [], []

    # Passo 1: Calcular o comprimento do arco da elipse (de 0 a pi/2) de forma precisa
    # A função a ser integrada para o comprimento do arco de uma elipse
    integrando = lambda teta: np.sqrt((a * np.sin(teta))**2 + (b * np.cos(teta))**2)
    comprimento_arco_total, _ = quad(integrando, 0, np.pi / 2)

    # Passo 2: Criar um mapa (interpolação) entre ângulo e comprimento do arco
    angulos_mapa = np.linspace(0, np.pi / 2, 200) # 200 pontos para um mapa preciso
    arcos_mapa = np.array([quad(integrando, 0, t)[0] for t in angulos_mapa])
    # Invertemos o mapa: dado um comprimento de arco, qual é o ângulo?
    arco_para_angulo = interp1d(arcos_mapa, angulos_mapa, fill_value="extrapolate")

    # Passo 3: Determinar as posições dos anéis ao longo do arco
    num_aneis_quadrante = int(comprimento_arco_total / d_grid)
    if num_aneis_quadrante == 0: num_aneis_quadrante = 1
    
    # Gera os comprimentos de arco onde cada anel deve estar
    arcos_desejados = np.linspace(0, comprimento_arco_total, num_aneis_quadrante)
    
    # Usa nosso mapa para encontrar os ângulos correspondentes
    angulos_para_aneis = arco_para_angulo(arcos_desejados)

    # Listas para armazenar todas as coordenadas e pesos
    all_x, all_y, all_z, all_pesos = [], [], [], []

    # Passo 4: Gerar os pontos para cada anel
    for angulo in angulos_para_aneis:
        # Raio e altura (z) do anel para este ângulo
        raio_anel = b * np.cos(angulo)
        z_pos = a * np.sin(angulo)

        if raio_anel < d_grid / 2: # Evita anéis muito pequenos ou sobrepostos no polo
            # Adiciona apenas um ponto no polo
            if not any(np.isclose(z_pos, np.array(all_z))): # Evita duplicatas
                all_x.append(0.0)
                all_y.append(0.0)
                all_z.append(z_pos)
                all_pesos.append(peso)
                if z_pos > 1e-6: # Adiciona o polo sul se não estiver no equador
                    all_x.append(0.0)
                    all_y.append(0.0)
                    all_z.append(-z_pos)
                    all_pesos.append(peso)
            continue
            
        # Calcula o número de esferas para preencher a circunferência do anel
        num_pontos_anel = int(2 * np.pi * raio_anel / d_grid)
        if num_pontos_anel < 1: continue

        # Gera os pontos para o hemisfério norte
        x_n, y_n, z_n = _gerar_pontos_no_anel(raio_anel, z_pos, num_pontos_anel)
        all_x.extend(x_n)
        all_y.extend(y_n)
        all_z.extend(z_n)
        all_pesos.extend([peso] * num_pontos_anel)

        # Gera os pontos para o hemisfério sul (se não for o equador)
        if z_pos > 1e-6: # Evita duplicar o anel do equador
            x_s, y_s, z_s = _gerar_pontos_no_anel(raio_anel, -z_pos, num_pontos_anel)
            all_x.extend(x_s)
            all_y.extend(y_s)
            all_z.extend(z_s)
            all_pesos.extend([peso] * num_pontos_anel)

    return all_x, all_y, all_z, all_pesos
@numba.njit
def gerar_circulo(raio, d_grid, peso=1.0):
    """Gera coordenadas e pesos para um círculo (anel 2D) no plano XY."""
    x_vals = List.empty_list(types.float64)
    y_vals = List.empty_list(types.float64)
    z_vals = List.empty_list(types.float64)
    pesos = List.empty_list(types.float64)

    # Se o raio for muito pequeno, retorna vazio
    if raio <= 0:
        return x_vals, y_vals, z_vals, pesos

    # Calcula quantos pontos são necessários para preencher a circunferência com o espaçamento d_grid
    circunferencia = 2 * np.pi * raio
    num_pontos = int(round(circunferencia / d_grid))
    
    # Garante que tenhamos pelo menos um ponto se o raio for maior que zero
    if num_pontos == 0:
        num_pontos = 1

    for i in range(num_pontos):
        angulo = i * (2 * np.pi / num_pontos)
        x = raio * np.cos(angulo)
        y = raio * np.sin(angulo)
        
        x_vals.append(x)
        y_vals.append(y)
        z_vals.append(0.0) # Z é sempre 0 para um círculo 2D
        pesos.append(peso)
            
    return x_vals, y_vals, z_vals, pesos
def gerar_meia_casca_elipsoidal(a, b, d_grid, peso=1.0):
    """
    Gera coordenadas e pesos para uma MEIA casca elipsoidal (hemisfério superior).
    Para uma meia casca esférica, basta fazer a = b.
    """
    if a == 0 or b == 0:
        return [], [], [], []

    # Passos 1, 2 e 3 são idênticos à função da casca completa
    integrando = lambda teta: np.sqrt((a * np.sin(teta))**2 + (b * np.cos(teta))**2)
    comprimento_arco_total, _ = quad(integrando, 0, np.pi / 2)
    angulos_mapa = np.linspace(0, np.pi / 2, 200)
    arcos_mapa = np.array([quad(integrando, 0, t)[0] for t in angulos_mapa])
    arco_para_angulo = interp1d(arcos_mapa, angulos_mapa, fill_value="extrapolate")
    num_aneis_quadrante = int(comprimento_arco_total / d_grid)
    if num_aneis_quadrante == 0: num_aneis_quadrante = 1
    arcos_desejados = np.linspace(0, comprimento_arco_total, num_aneis_quadrante)
    angulos_para_aneis = arco_para_angulo(arcos_desejados)

    all_x, all_y, all_z, all_pesos = [], [], [], []

    # Passo 4: Gerar os pontos para cada anel (APENAS HEMISFÉRIO NORTE)
    for angulo in angulos_para_aneis:
        raio_anel = b * np.cos(angulo)
        z_pos = a * np.sin(angulo)

        if raio_anel < d_grid / 2:
            if not any(np.isclose(z_pos, np.array(all_z))):
                all_x.append(0.0)
                all_y.append(0.0)
                all_z.append(z_pos)
                all_pesos.append(peso)
            continue
            
        num_pontos_anel = int(2 * np.pi * raio_anel / d_grid)
        if num_pontos_anel < 1: continue

        # Gera os pontos APENAS para o hemisfério norte (Z positivo)
        x_n, y_n, z_n = _gerar_pontos_no_anel(raio_anel, z_pos, num_pontos_anel)
        all_x.extend(x_n)
        all_y.extend(y_n)
        all_z.extend(z_n)
        all_pesos.extend([peso] * num_pontos_anel)
        
        # A parte que gerava o hemisfério sul (com Z negativo) foi REMOVIDA.

    return all_x, all_y, all_z, all_pesos


def plota_intensidade (arquivo):        
    q_vals, I_vals, erro_q,erro_I = [], [], [],[]

    with open(arquivo, "r") as file:
        linhas = file.readlines()
    for linha in linhas [3:]:
        valores = linha.split()
        if len(valores)==4:
            q,Iq,dq,dI = map(float,valores)
            q_vals.append(q)
            I_vals.append(Iq)
            erro_q.append(dq)
            erro_I.append(dI)
    q_vals = np.array(q_vals)
    I_vals = np.array(I_vals)
    erro_q = np.array(erro_q)
    erro_I = np.array(erro_I)
    plt.figure(figsize=(8,6))
    plt.errorbar(q_vals,I_vals,yerr=None,xerr=None,fmt='o', capsize=3)
    plt.xlabel(r'$q$ $(\AA^{-1})$',fontsize = 14)
    plt.ylabel(r"$I(q)$ (A.U.)",fontsize =14)
    plt.xscale("log")
    plt.yscale("log")
    plt.title("Scattered Intensity",fontsize = 14)
    plt.legend()
    plt.grid(True,which='both', linestyle = '--',alpha = 0.5)
    plt.show ()


def plot_POR(arquivo):
    p_vals,r_vals,erro_p,erro_r = [],[],[],[]
    with open (arquivo,'r') as file:
        linhas = file.readlines()
    for linha  in linhas [3:]:
        valores = linha.split()
        if len(valores) ==4 :
            p,r,dp,dr = map(float,valores)
            p_vals.append(p)
            r_vals.append(r)
            erro_p.append(dp)
            erro_r.append(dr)
    p_vals = np.array(p_vals)
    r_vals = np.array(r_vals)
    erro_p = np.array(erro_p)
    erro_r = np.array(erro_r)
    plt.figure(figsize=(8,6))
    # --- CORREÇÃO AQUI: eixos invertidos para r no x e p no y ---
    plt.errorbar(r_vals, p_vals, yerr=erro_p, xerr=erro_r, fmt='o', capsize=3)
    plt.xlabel(r'r (Å)', fontsize = 14) # Unidade mais comum é Angstrom
    plt.ylabel(r"P(r)", fontsize =14)
    plt.title("Função de Distribuição de Distâncias", fontsize = 14)
    plt.grid(True,which='both', linestyle = '--',alpha = 0.5)
    plt.show ()


def plota_intensidade_vs_qRg(arquivo_niq, arquivo_rg):
    # Lê a intensidade do arquivo .NIQ
    q_vals, I_vals = [], []
    with open(arquivo_niq, "r") as file:
        linhas = file.readlines()[3:]  # Pula cabeçalho

    for linha in linhas:
        valores = linha.split()
        if len(valores) == 4:
            q, Iq, _, _ = map(float, valores)
            q_vals.append(q)
            I_vals.append(Iq)

    q_vals = np.array(q_vals)
    I_vals = np.array(I_vals)

    # Lê o Rg da primeira coluna da primeira linha útil
    with open(arquivo_rg, "r") as file:
        for linha in file:
            if linha.strip():
                rg_str = linha.strip().split()[0]
                rg = float(rg_str)
                break

    # Plota I(q) vs q * Rg
    plt.figure(figsize=(8, 6))
    plt.plot(q_vals * rg, I_vals, '-', color='purple' ,label=r"$I(q)$ vs $q \cdot Rg$")
    plt.xlabel(r'$q \cdot Rg$', fontsize=14)
    plt.ylabel(r'$I(q)$', fontsize=14)
    plt.xscale('log')
    plt.yscale('log')
    plt.title("Scattered Intensity", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.show()
            
def plota_POR_vs_rRg(arquivo_por, arquivo_rg):
    # Lê o Rg da primeira coluna da primeira linha útil
    with open(arquivo_rg, "r") as file:
        for linha in file:
            if linha.strip():
                rg_str = linha.strip().split()[0]
                rg = float(rg_str)
                break

    # Lê os dados do arquivo .POR
    r_vals = []
    p_vals = []
    with open(arquivo_por, 'r') as file:
        linhas = file.readlines()[3:]  # Pula cabeçalho
        for linha in linhas:
            valores = linha.split()
            if len(valores) >= 2:
                r = float(valores[0])
                p = float(valores[1])
                r_vals.append(r)
                p_vals.append(p)

    r_vals = np.array(r_vals)
    p_vals = np.array(p_vals)
    r_div_rg = r_vals / rg  # Normaliza por Rg

    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(r_div_rg, p_vals, 'o-', color='darkgreen', label=r"$p(r)$")
    plt.xlabel(r"$r / R_g$", fontsize=14)
    plt.ylabel(r"$p(r)$", fontsize=14)
    plt.title(r"p(r)", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.show()

class Geometria:
    def __init__(self,tipo,centro,parametros):
        self.tipo = tipo
        self.centro = np.array(centro)
        self.parametros = parametros
        self.coordenadas = []
        self.pesos = []
        self.CM = None  # Inicializa o atributo CM
    def _ajustar_coordenadas(self,x_vals,y_vals,z_vals):
        return [(x+self.centro[0],y+ self.centro[1], z+self.centro[2]) for x,y,z in zip(x_vals,y_vals,z_vals)]
    def calcular_CM(self, coordenadas):
        coord_array = np.array(coordenadas)
        centro_de_massa = np.mean(coord_array, axis=0)
        return centro_de_massa    
    def gerar_coordenadas(self):
        if self.tipo == 1: 
            self.coordenadas,self.pesos = self._gera_esfera()
        elif self.tipo ==2:
            self.coordenadas,self.pesos= self._gera_cilindro()
        elif self.tipo ==3:
            self.coordenadas,self.pesos = self._gera_elipsoide()
        elif self.tipo ==4: 
            self.coordenadas,self.pesos = self._gera_prisma()
        elif self.tipo ==5:
            self.coordenadas,self.pesos = self._gerar_casca_elipsoidal()
        elif self.tipo ==6: 
            self.coordenadas, self.pesos = self._gerar_circulo()
        elif self.tipo ==7: 
            self.coordenadas,self.pesos = self._gerar_meia_casca_elipsoidal()
        else:
            raise ValueError('Tipo de geometria desconhecido')
        if self.coordenadas:   
            self.CM = self.calcular_CM(self.coordenadas)
    
    def _get_peso(self):
        # Método auxiliar para pegar o peso do dicionário ou usar o padrão 1.0
        return self.parametros.get('peso', 1.0)

    def _gera_esfera(self):
        R, d_grid = self.parametros['R'], self.parametros['d_grid']
        x, y, z, p = gerar_esfera(R, d_grid, self._get_peso())
        return self._ajustar_coordenadas(x, y, z), p
    def _gera_cilindro(self):
        re, ri, alt, d_grid = self.parametros['raio_externo'], self.parametros['raio_interno'], self.parametros['Altura'], self.parametros['d_grid']
        x, y, z, p = gerar_cilindro(re, ri, alt, d_grid, self._get_peso())
        return self._ajustar_coordenadas(x, y, z), p
    def _gera_elipsoide(self):
        a, b, c, d_grid = self.parametros['a'], self.parametros['b'], self.parametros['c'], self.parametros['d_grid']
        x, y, z, p = gerar_elipsoide(a, b, c, d_grid, self._get_peso())
        return self._ajustar_coordenadas(x, y, z), p
    def _gera_prisma(self):
        a, b, c, d_grid = self.parametros['a'], self.parametros['b'], self.parametros['c'], self.parametros['d_grid']
        x, y, z, p = gerar_prisma(a, b, c, d_grid, self._get_peso())
        return self._ajustar_coordenadas(x, y, z), p
    def _gerar_casca_elipsoidal(self):
        a, b, d_grid = self.parametros['a'], self.parametros['b'], self.parametros['d_grid']
        x, y, z, p = gerar_casca_elipsoidal(a, b, d_grid, self._get_peso())
        return self._ajustar_coordenadas(x,y,z),p


    def _gerar_circulo(self):
        raio, d_grid= self.parametros['Raio'],self.parametros['d_grid']

        x,y,z,p = gerar_circulo(raio,d_grid,self._get_peso())
    def _gerar_meia_casca_elipsoidal(self):
        a,b, d_grid = self.parametros['a'],self.parametros['b'], self.parametros['d_grid']
        x,y,z,p = gerar_meia_casca_elipsoidal(a,b,d_grid)
        return self._ajustar_coordenadas(x,y,z),p
class GerenciadorGeometrias:
    def __init__(self):
        self.geometrias =[]
    
    def adicionar_geometria(self,tipo,centro,parametros):
        '''Adiciona uma nova geometria ao sistema'''
        nova_geometria = Geometria(tipo,centro,parametros)
        nova_geometria.gerar_coordenadas()
        self.geometrias.append(nova_geometria)
        
    def obter_coordenadas_totais(self):
        coordenadas=[]
        for geo in self.geometrias:
            coordenadas.extend(geo.coordenadas)
            
        return coordenadas
    def calcular_CM(self, coordenadas):
        coord_array = np.array(coordenadas)
        centro_de_massa = np.mean(coord_array, axis=0)
        return centro_de_massa
    def salvar_modelo_dat(self, arquivo):
        print(f"Salvando modelo no formato .dat em '{arquivo}'...")
        with open(arquivo, "w") as file:
            # Escreve o cabeçalho (opcional, mas bom para referência)
            file.write("# Formato: X Y Z Peso\n")
            
            total_coords = 0
            for geo in self.geometrias:
                for i, coord in enumerate(geo.coordenadas):
                    x, y, z = coord
                    # Pega o peso correspondente da lista de pesos
                    peso = geo.pesos[i]
                    file.write(f"{x:10.4f} {y:10.4f} {z:10.4f} {peso:10.4f}\n")
                total_coords += len(geo.coordenadas)
        print(f"Modelo com {total_coords} subunidades salvo com sucesso.")
    def visualizar_com_pesos(self, dat_filename):
        """
        Lê um arquivo .dat no formato 'X Y Z Peso' e plota um gráfico 3D,
        usando os pesos para colorir as subunidades.
        """
        x_vals, y_vals, z_vals, pesos = [], [], [], []

        print(f"Lendo dados de '{dat_filename}' para visualização colorida...")
        with open(dat_filename, "r") as file:
            for line in file:
                # Ignora linhas de comentário ou vazias
                if line.strip().startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split()
                if len(parts) == 4:
                    x, y, z, peso = map(float, parts)
                    x_vals.append(x)
                    y_vals.append(y)
                    z_vals.append(z)
                    pesos.append(peso)

        if not x_vals:
            print("Nenhum dado encontrado para visualização.")
            return

        # Criando a figura 3D
        fig = plt.figure(figsize=(10, 8)) # Um pouco maior para a barra de cores
        ax = fig.add_subplot(111, projection='3d')
        
        scatter_plot = ax.scatter(x_vals, y_vals, z_vals, marker='o', s=20, c=pesos, cmap='viridis')

        # Adicionando a barra de cores como legenda
        cbar = fig.colorbar(scatter_plot, ax=ax, shrink=0.6)
        cbar.set_label('Peso da Subunidade')

        # Configurando rótulos dos eixos
        ax.set_xlabel("Eixo X")
        ax.set_ylabel("Eixo Y")
        ax.set_zlabel("Eixo Z")
        ax.set_title("Visualização 3D com Pesos")
        
        # Melhora a visualização para garantir que todos os pontos apareçam
        max_range = np.array([max(x_vals)-min(x_vals), max(y_vals)-min(y_vals), max(z_vals)-min(z_vals)]).max() / 2.0
        mid_x = (max(x_vals)+min(x_vals)) * 0.5
        mid_y = (max(y_vals)+min(y_vals)) * 0.5
        mid_z = (max(z_vals)+min(z_vals)) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)


        plt.show()

    def exibir_resumo(self):
        for i, geo in enumerate(self.geometrias):
            print(f"Geometrias{i+1}: {geo.tipo}, Origem: {geo.centro}, quantidade de elementos: {len(geo.coordenadas)}, Centro_de_Massa = {geo.CM}")

           
    def conversao_dat_pdb(self, dat_filename, output_dir):
        """
        Converte um arquivo .dat no formato 'X Y Z Peso' para o formato .pdb,
        lidando corretamente com um número de átomos > 99999.
        """
        os.makedirs(output_dir, exist_ok=True)
        pdb_filename = os.path.join(output_dir, "Model.pdb")

        print(f"Convertendo '{dat_filename}' para '{pdb_filename}'...")
        
        with open(dat_filename, 'r') as dat_file, open(pdb_filename, 'w') as pdb_file:
            lines = dat_file.readlines()
            
            atom_index = 1
            for line in lines:
                if line.strip().startswith('#') or not line.strip():
                    continue

                parts = line.strip().split()
                
                if len(parts) == 4:
                    x, y, z, peso = map(float, parts)

                    # --- LÓGICA PARA LIDAR COM NÚMEROS GRANDES ---
                    # Garante que o número do átomo não passe de 5 dígitos
                    atom_serial_pdb = atom_index % 100000
                    # Garante que o número do resíduo não passe de 4 dígitos
                    residue_seq_pdb = atom_index % 10000
                    
                    # --- LINHA COM FORMATAÇÃO ROBUSTA ---
                    # Usamos as novas variáveis para garantir o alinhamento
                    pdb_line = f"ATOM  {atom_serial_pdb:5d}  CA  ASP  {residue_seq_pdb:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  {peso:6.2f}  1.00 0 2 201 C  \n"
                    
                    pdb_file.write(pdb_line)
                    atom_index += 1

        print(f"Arquivo PDB salvo com sucesso em: {pdb_filename}")
    def rotacionar(self, geometria, angle_degrees, axis, ponto_rotacao=None):
        """
        Rotaciona uma geometria em torno de um ponto específico ou do centro de massa.
        
        Args:
            geometria: Objeto Geometria a ser rotacionado
            angle_degrees: Ângulo de rotação em graus
            axis: Eixo de rotação ('x', 'y' ou 'z')
            ponto_rotacao: Ponto em torno do qual rotacionar (None para usar o centro de massa)
        """
        angle_radians = np.radians(angle_degrees)
        
        # Verificar se existem coordenadas
        if not geometria.coordenadas:
            raise ValueError("A geometria não possui coordenadas definidas")
        
        # Determinar o ponto de rotação (CM se não for fornecido)
        if ponto_rotacao is None:
            if geometria.CM is None:
                geometria.CM = self.calcular_CM(geometria.coordenadas)
            ponto_ref = np.array(geometria.CM)
        else:
            ponto_ref = np.array(ponto_rotacao)
        
        # Verificar formato do ponto de referência
        if ponto_ref.shape != (3,):
            raise ValueError(f"Ponto de rotação deve ter formato (3,), recebido {ponto_ref.shape}")
        
        # Matrizes de rotação
        if axis == "x":
            rotation_matrix = np.array([
                [1, 0, 0],
                [0, np.cos(angle_radians), -np.sin(angle_radians)],
                [0, np.sin(angle_radians), np.cos(angle_radians)]
            ])
        elif axis == "y":
            rotation_matrix = np.array([
                [np.cos(angle_radians), 0, np.sin(angle_radians)],
                [0, 1, 0],
                [-np.sin(angle_radians), 0, np.cos(angle_radians)]
            ])
        elif axis == "z":
            rotation_matrix = np.array([
                [np.cos(angle_radians), -np.sin(angle_radians), 0],
                [np.sin(angle_radians), np.cos(angle_radians), 0],
                [0, 0, 1]
            ])
        else:
            raise ValueError("Eixo inválido. Escolha 'x', 'y' ou 'z'.")

        # Converter coordenadas para array numpy
        coords = np.array(geometria.coordenadas)
        
        # Verificar formato das coordenadas
        if coords.size == 0:
            raise ValueError("Nenhuma coordenada encontrada para rotacionar")
        
        # 1. Transladar pontos para o sistema de referência do ponto de rotação
        pontos_relativos = coords - ponto_ref
        
        # 2. Aplicar rotação
        pontos_rotacionados = np.dot(pontos_relativos, rotation_matrix.T)
        
        # 3. Transladar de volta para o sistema original
        novas_coords = pontos_rotacionados + ponto_ref
        
        # Atualizar geometria
        geometria.coordenadas = [tuple(coord) for coord in novas_coords]
        
        # Recalcular CM se estivermos rotacionando em torno dele
        if ponto_rotacao is None:
            geometria.CM = self.calcular_CM(geometria.coordenadas)
        
        return geometria.coordenadas



gerenciador = GerenciadorGeometrias()

# Parâmetros Globais
d_grid_resolucao = 3
R_subunidade = 0.62035 * d_grid_resolucao
print(f"Usando d_grid = {d_grid_resolucao}, o raio da subunidade para o cálculo será: {R_subunidade:.4f}")


parametros_esfera = {
    'R': 50, 
    'd_grid': 4, 
    'peso': 1.0  # Peso definido como 3.0
}
gerenciador.adicionar_geometria(tipo=1, centro=(0,0,0), parametros=parametros_esfera)


'''
parametros_meia_casca = {'a': 30,'b':30 ,'d_grid': d_grid_resolucao}
gerenciador.adicionar_geometria(tipo=7, centro=(0, 0, 50), parametros=parametros_meia_casca)
parametros_cilindro = {
    'raio_externo': 4, 
    'raio_interno': 0,      # Cilindro sólido
    'Altura': 20,          # Comprimento da barra
    'd_grid': d_grid_resolucao, 
    'peso': 1.0             # Peso definido como 1.0
}
gerenciador.adicionar_geometria(tipo=2, centro=(-6, 0, 80), parametros=parametros_cilindro)
parametros_cilindro_2= {
    'raio_externo': 4, 
    'raio_interno': 0,      # Cilindro sólido
    'Altura': 40,          # Comprimento da barra
    'd_grid': d_grid_resolucao, 
    'peso': 1.0             # Peso definido como 1.0
}
gerenciador.adicionar_geometria(tipo=2, centro=(7, 0, 90), parametros=parametros_cilindro_2)
'''




'''  # 1. Parâmetros para as Esferas (com peso 3.0)
parametros_esfera_pesada = {
    'R': 30, 
    'd_grid': d_grid_resolucao, 
    'peso': 3.0  # Peso definido como 3.0
}

# 2. Parâmetros para o Cilindro (com peso 1.0 e chaves corretas)
parametros_cilindro_leve = {
    'raio_externo': 15, 
    'raio_interno': 0,      # Cilindro sólido
    'Altura': 100,          # Comprimento da barra
    'd_grid': d_grid_resolucao, 
    'peso': 1.0             # Peso definido como 1.0
}

# 3. Adicionar a primeira esfera (à esquerda)
print("Criando esfera 1 (peso 3.0)...")
gerenciador.adicionar_geometria(
    tipo=1, 
    centro=(-50, 0, 0), # Posiciona no final da barra do cilindro
    parametros=parametros_esfera_pesada
)

# 4. Adicionar o cilindro (no centro)
print("Criando cilindro central (peso 1.0)...")
gerenciador.adicionar_geometria(
    tipo=2, 
    centro=(0, 0, 0), # Cria na origem para rotacionar
    parametros=parametros_cilindro_leve
)
# Para rotacionar, primeiro pegamos uma referência ao último objeto adicionado
cilindro_central = gerenciador.geometrias[-1] 
# Rotaciona 90 graus no eixo 'y' para deitá-lo ao longo do eixo 'x'
gerenciador.rotacionar(cilindro_central, 90, 'y')

# 5. Adicionar a segunda esfera (à direita)
print("Criando esfera 2 (peso 3.0)...")
gerenciador.adicionar_geometria(
    tipo=1, 
    centro=(50, 0, 0), # Posiciona na outra ponta da barra
    parametros=parametros_esfera_pesada
)

'''
'''
parametros_meia_casca = {'a': 14,'b':14 ,'d_grid': d_grid_resolucao}
gerenciador.adicionar_geometria(tipo=7, centro=(0, 0, 25), parametros=parametros_meia_casca)

parametros_cilindro = {'raio_externo':15, 'raio_interno':0,'Altura':50,'d_grid':d_grid_resolucao}
gerenciador.adicionar_geometria(tipo=2, centro=(0, 0, 0), parametros=parametros_cilindro)

parametros_cilindro = {'raio_externo':5, 'raio_interno':0,'Altura':35,'d_grid':d_grid_resolucao}
gerenciador.adicionar_geometria(tipo=2, centro=(18, 0, -10), parametros=parametros_cilindro)
parametros_cilindro = {'raio_externo':5, 'raio_interno':0,'Altura':35,'d_grid':d_grid_resolucao}
gerenciador.adicionar_geometria(tipo=2, centro=(-18, 0, -10), parametros=parametros_cilindro)


parametros_prisma = {'a':15, 'b':25,'c':10,'d_grid':d_grid_resolucao}
gerenciador.adicionar_geometria(tipo=4, centro=(18, 0, -30), parametros=parametros_prisma)
parametros_prisma = {'a':15, 'b':25,'c':10,'d_grid':d_grid_resolucao}
gerenciador.adicionar_geometria(tipo=4, centro=(-18, 0, -30), parametros=parametros_prisma)
r2d2_coords = gerenciador.obter_coordenadas_totais()

pivo_rotacao = gerenciador.calcular_CM(r2d2_coords)
print(f"Ponto de rotação (CM do modelo): {np.round(pivo_rotacao, 2)}")

# 3. Faz um loop por cada peça do R2D2 e rotaciona individualmente
for peca_do_r2d2 in gerenciador.geometrias:
    gerenciador.rotacionar(peca_do_r2d2, -9, 'z', ponto_rotacao=pivo_rotacao)
'''



# --- Salvando e Visualizando ---

gerenciador.exibir_resumo()

# Salva o modelo final no formato DAT (X Y Z Peso)
arquivo_final = "model.dat"
gerenciador.salvar_modelo_dat(arquivo_final)

# Chama a NOVA função de visualização!
gerenciador.visualizar_com_pesos(arquivo_final)

# O resto do seu fluxo de trabalho (conversão para PDB, cálculo, etc.)

gerenciador.conversao_dat_pdb("Model.dat", r'C:\Interface-python\teste_criacao_modelo_py')

# Execução dos cálculos adicionais
caminho_executavel = r"C:\Interface-python\teste_criacao_modelo_py\calcIQPR_DM.exe"
subprocess.run(caminho_executavel)  
intensidade = r"C:\Interface-python\teste_criacao_modelo_py\MODEL.NIQ"
raio_de_giro = r"C:\Interface-python\teste_criacao_modelo_py\RG_DMAX.DAT"
POR = r"C:\Interface-python\teste_criacao_modelo_py\MODEL.POR"

plota_intensidade_vs_qRg(intensidade, raio_de_giro)
plota_POR_vs_rRg(POR, raio_de_giro)
plota_intensidade(intensidade)
plot_POR(POR)