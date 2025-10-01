import sys
import os
import subprocess
import numpy as np
import numba
from numba.typed import List
from numba.core import types
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import quad
from scipy.interpolate import interp1d

from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                             QPushButton, QLabel, QLineEdit, QComboBox, QFileDialog, 
                             QStackedWidget, QGroupBox, QDoubleSpinBox, QTextEdit)
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt

# ==============================================================================
# ===== SUAS FUNÇÕES ORIGINAIS DO Main.py (sem nenhuma alteração) ================
# ==============================================================================

@numba.njit
def gerar_esfera(R_esfera, d_grid, peso=1.0):
    """Gera coordenadas e pesos para uma esfera sólida."""
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
    if a == 0 or b == 0:
        return [], [], [], []
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
    for angulo in angulos_para_aneis:
        raio_anel = b * np.cos(angulo)
        z_pos = a * np.sin(angulo)
        if raio_anel < d_grid / 2:
            if not any(np.isclose(z_pos, np.array(all_z))):
                all_x.append(0.0); all_y.append(0.0); all_z.append(z_pos); all_pesos.append(peso)
                if z_pos > 1e-6:
                    all_x.append(0.0); all_y.append(0.0); all_z.append(-z_pos); all_pesos.append(peso)
            continue
        num_pontos_anel = int(2 * np.pi * raio_anel / d_grid)
        if num_pontos_anel < 1: continue
        x_n, y_n, z_n = _gerar_pontos_no_anel(raio_anel, z_pos, num_pontos_anel)
        all_x.extend(x_n); all_y.extend(y_n); all_z.extend(z_n); all_pesos.extend([peso] * num_pontos_anel)
        if z_pos > 1e-6:
            x_s, y_s, z_s = _gerar_pontos_no_anel(raio_anel, -z_pos, num_pontos_anel)
            all_x.extend(x_s); all_y.extend(y_s); all_z.extend(z_s); all_pesos.extend([peso] * num_pontos_anel)
    return all_x, all_y, all_z, all_pesos

@numba.njit
def gerar_circulo(raio, d_grid, peso=1.0):
    x_vals = List.empty_list(types.float64)
    y_vals = List.empty_list(types.float64)
    z_vals = List.empty_list(types.float64)
    pesos = List.empty_list(types.float64)
    if raio <= 0:
        return x_vals, y_vals, z_vals, pesos
    circunferencia = 2 * np.pi * raio
    num_pontos = int(round(circunferencia / d_grid))
    if num_pontos == 0: num_pontos = 1
    for i in range(num_pontos):
        angulo = i * (2 * np.pi / num_pontos)
        x = raio * np.cos(angulo)
        y = raio * np.sin(angulo)
        x_vals.append(x); y_vals.append(y); z_vals.append(0.0); pesos.append(peso)
    return x_vals, y_vals, z_vals, pesos

# Funções de plotagem
def plota_intensidade(arquivo):        
    q_vals, I_vals, erro_q,erro_I = [], [], [],[]
    with open(arquivo, "r") as file:
        linhas = file.readlines()
    for linha in linhas [3:]:
        valores = linha.split()
        if len(valores)==4:
            q,Iq,dq,dI = map(float,valores)
            q_vals.append(q); I_vals.append(Iq); erro_q.append(dq); erro_I.append(dI)
    q_vals, I_vals = np.array(q_vals), np.array(I_vals)
    plt.figure(figsize=(8,6))
    plt.errorbar(q_vals,I_vals,yerr=None,xerr=None,fmt='o', capsize=3)
    plt.xlabel(r'$q$ $(\AA^{-1})$',fontsize = 14)
    plt.ylabel(r"$I(q)$ (A.U.)",fontsize =14)
    plt.xscale("log"); plt.yscale("log")
    plt.title("Scattered Intensity",fontsize = 14)
    plt.grid(True,which='both', linestyle = '--',alpha = 0.5); plt.show()

def plot_POR(arquivo):
    p_vals,r_vals,erro_p,erro_r = [],[],[],[]
    with open (arquivo,'r') as file:
        linhas = file.readlines()
    for linha  in linhas [3:]:
        valores = linha.split()
        if len(valores) ==4:
            p,r,dp,dr = map(float,valores)
            p_vals.append(p); r_vals.append(r); erro_p.append(dp); erro_r.append(dr)
    p_vals, r_vals = np.array(p_vals), np.array(r_vals)
    plt.figure(figsize=(8,6))
    plt.errorbar(r_vals, p_vals, yerr=erro_p, xerr=erro_r, fmt='o', capsize=3)
    plt.xlabel(r'r (Å)', fontsize = 14)
    plt.ylabel(r"P(r)", fontsize =14)
    plt.title("Função de Distribuição de Distâncias", fontsize = 14)
    plt.grid(True,which='both', linestyle = '--',alpha = 0.5); plt.show()

def plota_intensidade_vs_qRg(arquivo_niq, arquivo_rg):
    q_vals, I_vals = [], []
    with open(arquivo_niq, "r") as file:
        for linha in file.readlines()[3:]:
            valores = linha.split()
            if len(valores) == 4:
                q, Iq, _, _ = map(float, valores)
                q_vals.append(q); I_vals.append(Iq)
    q_vals, I_vals = np.array(q_vals), np.array(I_vals)
    with open(arquivo_rg, "r") as file:
        rg = float(file.readline().strip().split()[0])
    plt.figure(figsize=(8, 6)); plt.plot(q_vals * rg, I_vals, '-', color='purple' ,label=r"$I(q)$ vs $q \cdot Rg$")
    plt.xlabel(r'$q \cdot Rg$', fontsize=14); plt.ylabel(r'$I(q)$', fontsize=14)
    plt.xscale('log'); plt.yscale('log')
    plt.title("Scattered Intensity", fontsize=14); plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(); plt.tight_layout(); plt.show()
            
def plota_POR_vs_rRg(arquivo_por, arquivo_rg):
    with open(arquivo_rg, "r") as file:
        rg = float(file.readline().strip().split()[0])
    r_vals, p_vals = [], []
    with open(arquivo_por, 'r') as file:
        for linha in file.readlines()[3:]:
            valores = linha.split()
            if len(valores) >= 2:
                r, p = float(valores[0]), float(valores[1])
                r_vals.append(r); p_vals.append(p)
    r_vals, p_vals = np.array(r_vals), np.array(p_vals)
    plt.figure(figsize=(8, 6)); plt.plot(r_vals / rg, p_vals, 'o-', color='darkgreen', label=r"$p(r)$")
    plt.xlabel(r"$r / R_g$", fontsize=14); plt.ylabel(r"$p(r)$", fontsize=14)
    plt.title(r"p(r)", fontsize=14); plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(); plt.tight_layout(); plt.show()

# Classes de Geometria
class Geometria:
    def __init__(self,tipo,centro,parametros):
        self.tipo = tipo
        self.centro = np.array(centro)
        self.parametros = parametros
        self.coordenadas, self.pesos, self.CM = [], [], None
    def _ajustar_coordenadas(self,x_vals,y_vals,z_vals):
        return [(x+self.centro[0],y+ self.centro[1], z+self.centro[2]) for x,y,z in zip(x_vals,y_vals,z_vals)]
    def calcular_CM(self, coordenadas):
        return np.mean(np.array(coordenadas), axis=0)
    def gerar_coordenadas(self):
        func_map = {1: self._gera_esfera, 2: self._gera_cilindro, 3: self._gera_elipsoide, 
                    4: self._gera_prisma, 5: self._gerar_casca_elipsoidal, 6: self._gerar_circulo}
        if self.tipo in func_map:
            self.coordenadas, self.pesos = func_map[self.tipo]()
        else: raise ValueError('Tipo de geometria desconhecido')
        if self.coordenadas: self.CM = self.calcular_CM(self.coordenadas)
    def _get_peso(self): return self.parametros.get('peso', 1.0)
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
        raio, d_grid = self.parametros['Raio'], self.parametros['d_grid']
        x,y,z,p = gerar_circulo(raio,d_grid,self._get_peso())
        return self._ajustar_coordenadas(x,y,z),p

class GerenciadorGeometrias:
    def __init__(self):
        self.geometrias =[]
    def adicionar_geometria(self,tipo,centro,parametros):
        nova_geometria = Geometria(tipo,centro,parametros)
        nova_geometria.gerar_coordenadas()
        self.geometrias.append(nova_geometria)
    def salvar_modelo_dat(self, arquivo):
        with open(arquivo, "w") as file:
            file.write("# Formato: X Y Z Peso\n")
            total_coords = 0
            for geo in self.geometrias:
                for i, coord in enumerate(geo.coordenadas):
                    peso = geo.pesos[i]
                    file.write(f"{coord[0]:10.4f} {coord[1]:10.4f} {coord[2]:10.4f} {peso:10.4f}\n")
                total_coords += len(geo.coordenadas)
        print(f"Modelo com {total_coords} subunidades salvo em '{arquivo}'.")
    def visualizar_com_pesos(self, dat_filename):
        x_vals, y_vals, z_vals, pesos = [], [], [], []
        with open(dat_filename, "r") as file:
            for line in file:
                if line.strip().startswith('#') or not line.strip(): continue
                parts = line.strip().split()
                if len(parts) == 4:
                    x, y, z, peso = map(float, parts)
                    x_vals.append(x); y_vals.append(y); z_vals.append(z); pesos.append(peso)
        if not x_vals: return
        fig = plt.figure(figsize=(10, 8)); ax = fig.add_subplot(111, projection='3d')
        scatter_plot = ax.scatter(x_vals, y_vals, z_vals, marker='o', s=20, c=pesos, cmap='viridis')
        cbar = fig.colorbar(scatter_plot, ax=ax, shrink=0.6); cbar.set_label('Peso da Subunidade')
        ax.set_xlabel("Eixo X"); ax.set_ylabel("Eixo Y"); ax.set_zlabel("Eixo Z")
        ax.set_title("Visualização 3D com Pesos"); plt.show()
    def conversao_dat_pdb(self, dat_filename, output_dir):
        os.makedirs(output_dir, exist_ok=True)
        pdb_filename = os.path.join(output_dir, "Model.pdb")
        with open(dat_filename, 'r') as dat_file, open(pdb_filename, 'w') as pdb_file:
            atom_index = 1
            for line in dat_file:
                if line.strip().startswith('#') or not line.strip(): continue
                parts = line.strip().split()
                if len(parts) == 4:
                    x, y, z, peso = map(float, parts)
                    atom_serial_pdb = atom_index % 100000
                    residue_seq_pdb = atom_index % 10000
                    pdb_file.write(f"ATOM  {atom_serial_pdb:5d}  CA  ASP  {residue_seq_pdb:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  {peso:6.2f}  1.00 0 2 201 C  \n")
                    atom_index += 1
        print(f"Arquivo PDB salvo em: {pdb_filename}")


# ==============================================================================
# ===== NOVA INTERFACE GRÁFICA (PyQt5) =========================================
# ==============================================================================

class AppJanela(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Gerador de Modelos e Análise SAXS")
        self.setGeometry(100, 100, 800, 600)

        # Layout principal
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QHBoxLayout(self.central_widget)

        # Painel de controle (esquerda)
        self.control_panel = QGroupBox("Controles")
        self.control_layout = QVBoxLayout()
        self.control_panel.setLayout(self.control_layout)
        
        # Painel de Log (direita)
        self.log_panel = QGroupBox("Log do Processo")
        self.log_layout = QVBoxLayout()
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_layout.addWidget(self.log_text)
        self.log_panel.setLayout(self.log_layout)

        self.main_layout.addWidget(self.control_panel, 1) # 1/3 do espaço
        self.main_layout.addWidget(self.log_panel, 2)     # 2/3 do espaço

        self.init_ui()

    def log(self, message):
        self.log_text.append(message)
        QApplication.processEvents() # Atualiza a UI

    def init_ui(self):
        # 1. Seleção de Geometria
        self.control_layout.addWidget(QLabel("1. Escolha a Geometria:"))
        self.geo_combo = QComboBox()
        self.geo_map = {
            "Esfera": 1, "Cilindro": 2, "Elipsoide": 3,
            "Prisma": 4, "Casca Elipsoidal": 5, "Círculo (2D)": 6
        }
        self.geo_combo.addItems(self.geo_map.keys())
        self.control_layout.addWidget(self.geo_combo)

        # 2. Parâmetros Específicos (com StackedWidget)
        self.param_stack = QStackedWidget()
        self.create_param_widgets()
        self.control_layout.addWidget(self.param_stack)
        self.geo_combo.currentIndexChanged.connect(self.param_stack.setCurrentIndex)

        # 3. Parâmetros Gerais
        geral_group = QGroupBox("Parâmetros Gerais")
        geral_layout = QVBoxLayout()
        
        # d_grid
        d_grid_layout = QHBoxLayout()
        d_grid_layout.addWidget(QLabel("d_grid (resolução):"))
        self.d_grid_input = QDoubleSpinBox(); self.d_grid_input.setValue(3.0); self.d_grid_input.setSingleStep(0.1)
        d_grid_layout.addWidget(self.d_grid_input)
        geral_layout.addLayout(d_grid_layout)

        # Centro
        centro_layout = QHBoxLayout()
        centro_layout.addWidget(QLabel("Centro (X, Y, Z):"))
        self.cx = QDoubleSpinBox(); self.cx.setRange(-1000, 1000)
        self.cy = QDoubleSpinBox(); self.cy.setRange(-1000, 1000)
        self.cz = QDoubleSpinBox(); self.cz.setRange(-1000, 1000)
        centro_layout.addWidget(self.cx); centro_layout.addWidget(self.cy); centro_layout.addWidget(self.cz)
        geral_layout.addLayout(centro_layout)
        
        # Peso
        peso_layout = QHBoxLayout()
        peso_layout.addWidget(QLabel("Peso (contraste):"))
        self.peso_input = QDoubleSpinBox(); self.peso_input.setValue(1.0); self.peso_input.setSingleStep(0.1)
        peso_layout.addWidget(self.peso_input)
        geral_layout.addLayout(peso_layout)

        geral_group.setLayout(geral_layout)
        self.control_layout.addWidget(geral_group)

        # 4. Caminhos de Arquivos
        paths_group = QGroupBox("Caminhos")
        paths_layout = QVBoxLayout()
        
        # Diretório de Saída
        self.output_dir_input = QLineEdit(r'C:\Interface-python\teste_criacao_modelo_py')
        btn_output = QPushButton("Selecionar Pasta de Saída")
        btn_output.clicked.connect(self.select_output_dir)
        paths_layout.addWidget(btn_output); paths_layout.addWidget(self.output_dir_input)
        
        # Executável
        self.exe_path_input = QLineEdit(r'C:\Interface-python\teste_criacao_modelo_py\calcIQPR_DM.exe')
        btn_exe = QPushButton("Selecionar Executável")
        btn_exe.clicked.connect(self.select_exe_path)
        paths_layout.addWidget(btn_exe); paths_layout.addWidget(self.exe_path_input)
        
        paths_group.setLayout(paths_layout)
        self.control_layout.addWidget(paths_group)
        
        self.control_layout.addStretch()

        # 5. Botão de Execução
        self.run_button = QPushButton("Gerar Modelo e Analisar")
        self.run_button.setFont(QFont("Arial", 12, QFont.Bold))
        self.run_button.clicked.connect(self.run_analysis)
        self.control_layout.addWidget(self.run_button)

    def create_param_widgets(self):
        """Cria os widgets de parâmetros para cada geometria no QStackedWidget."""
        # Esfera
        widget = QWidget(); layout = QVBoxLayout(widget)
        self.esfera_R = QLineEdit("50"); layout.addWidget(QLabel("Raio (R):")); layout.addWidget(self.esfera_R)
        self.param_stack.addWidget(widget)
        
        # Cilindro
        widget = QWidget(); layout = QVBoxLayout(widget)
        self.cil_re = QLineEdit("50"); layout.addWidget(QLabel("Raio Externo:")); layout.addWidget(self.cil_re)
        self.cil_ri = QLineEdit("25"); layout.addWidget(QLabel("Raio Interno:")); layout.addWidget(self.cil_ri)
        self.cil_alt = QLineEdit("100"); layout.addWidget(QLabel("Altura:")); layout.addWidget(self.cil_alt)
        self.param_stack.addWidget(widget)

        # Elipsoide
        widget = QWidget(); layout = QVBoxLayout(widget)
        self.elip_a = QLineEdit("60"); layout.addWidget(QLabel("Semieixo a:")); layout.addWidget(self.elip_a)
        self.elip_b = QLineEdit("30"); layout.addWidget(QLabel("Semieixo b:")); layout.addWidget(self.elip_b)
        self.elip_c = QLineEdit("30"); layout.addWidget(QLabel("Semieixo c:")); layout.addWidget(self.elip_c)
        self.param_stack.addWidget(widget)

        # Prisma
        widget = QWidget(); layout = QVBoxLayout(widget)
        self.prisma_a = QLineEdit("50"); layout.addWidget(QLabel("Lado a:")); layout.addWidget(self.prisma_a)
        self.prisma_b = QLineEdit("50"); layout.addWidget(QLabel("Lado b:")); layout.addWidget(self.prisma_b)
        self.prisma_c = QLineEdit("50"); layout.addWidget(QLabel("Lado c:")); layout.addWidget(self.prisma_c)
        self.param_stack.addWidget(widget)
        
        # Casca Elipsoidal
        widget = QWidget(); layout = QVBoxLayout(widget)
        self.casca_a = QLineEdit("60"); layout.addWidget(QLabel("Semieixo a (vertical):")); layout.addWidget(self.casca_a)
        self.casca_b = QLineEdit("30"); layout.addWidget(QLabel("Semieixo b (raio equatorial):")); layout.addWidget(self.casca_b)
        self.param_stack.addWidget(widget)
        
        # Círculo
        widget = QWidget(); layout = QVBoxLayout(widget)
        self.circ_r = QLineEdit("50"); layout.addWidget(QLabel("Raio:")); layout.addWidget(self.circ_r)
        self.param_stack.addWidget(widget)

    def select_output_dir(self):
        dir = QFileDialog.getExistingDirectory(self, "Selecionar Pasta de Saída")
        if dir: self.output_dir_input.setText(dir)

    def select_exe_path(self):
        path, _ = QFileDialog.getOpenFileName(self, "Selecionar Executável", filter="Executable (*.exe)")
        if path: self.exe_path_input.setText(path)

    def get_params(self):
        """Coleta os parâmetros da UI."""
        params = {}
        tipo_str = self.geo_combo.currentText()
        tipo_int = self.geo_map[tipo_str]

        try:
            params['d_grid'] = self.d_grid_input.value()
            params['peso'] = self.peso_input.value()
            
            if tipo_int == 1: # Esfera
                params['R'] = float(self.esfera_R.text())
            elif tipo_int == 2: # Cilindro
                params['raio_externo'] = float(self.cil_re.text())
                params['raio_interno'] = float(self.cil_ri.text())
                params['Altura'] = float(self.cil_alt.text())
            elif tipo_int == 3: # Elipsoide
                params['a'] = float(self.elip_a.text())
                params['b'] = float(self.elip_b.text())
                params['c'] = float(self.elip_c.text())
            elif tipo_int == 4: # Prisma
                params['a'] = float(self.prisma_a.text())
                params['b'] = float(self.prisma_b.text())
                params['c'] = float(self.prisma_c.text())
            elif tipo_int == 5: # Casca Elipsoidal
                params['a'] = float(self.casca_a.text())
                params['b'] = float(self.casca_b.text())
            elif tipo_int == 6: # Círculo
                params['Raio'] = float(self.circ_r.text())

            return params
        except ValueError as e:
            self.log(f"ERRO: Parâmetro inválido! Verifique os números digitados. Detalhe: {e}")
            return None

    def run_analysis(self):
        self.log_text.clear()
        self.log("="*40)
        self.log(f"Iniciando processo para: {self.geo_combo.currentText()}")
        self.log("="*40)

        # 1. Coletar parâmetros
        self.log("[1/7] Coletando parâmetros da interface...")
        parametros = self.get_params()
        if not parametros: return
        
        tipo_selecionado = self.geo_map[self.geo_combo.currentText()]
        centro_selecionado = (self.cx.value(), self.cy.value(), self.cz.value())
        output_dir = self.output_dir_input.text()
        exe_path = self.exe_path_input.text()

        # Verifica se os caminhos existem
        if not os.path.isdir(output_dir):
            self.log(f"ERRO: O diretório de saída '{output_dir}' não existe!")
            return
        if not os.path.isfile(exe_path):
            self.log(f"ERRO: O executável '{exe_path}' não foi encontrado!")
            return

        # 2. Gerar geometria
        self.log("\n[2/7] Gerando a geometria...")
        gerenciador = GerenciadorGeometrias()
        gerenciador.adicionar_geometria(
            tipo=tipo_selecionado,
            centro=centro_selecionado,
            parametros=parametros
        )
        geo = gerenciador.geometrias[0]
        self.log(f"Geometria gerada com {len(geo.coordenadas)} subunidades.")

        # 3. Salvar .dat e visualizar
        arquivo_dat = os.path.join(output_dir, "model.dat")
        self.log(f"\n[3/7] Salvando modelo em '{arquivo_dat}'...")
        gerenciador.salvar_modelo_dat(arquivo_dat)
        self.log("... Salvo com sucesso.")
        
        self.log("\n[4/7] Mostrando visualização 3D (feche a janela para continuar)...")
        gerenciador.visualizar_com_pesos(arquivo_dat)

        # 4. Converter para .pdb
        self.log(f"\n[5/7] Convertendo para PDB em '{output_dir}'...")
        gerenciador.conversao_dat_pdb(arquivo_dat, output_dir)
        self.log("... Conversão concluída.")

        # 5. Executar cálculo externo
        self.log(f"\n[6/7] Executando cálculo externo: '{exe_path}'...")
        try:
            # Roda o executável no diretório de saída para que ele encontre os arquivos
            resultado = subprocess.run(exe_path, cwd=output_dir, capture_output=True, text=True, check=True)
            self.log("... Execução finalizada com sucesso.")
            self.log("Saída do processo:\n" + resultado.stdout)
        except subprocess.CalledProcessError as e:
            self.log(f"ERRO ao executar o processo externo!\nRetorno: {e.returncode}\nOutput: {e.stdout}\nErro: {e.stderr}")
            return
        except FileNotFoundError:
            self.log(f"ERRO: Executável não encontrado em '{exe_path}'")
            return

        # 6. Plotar resultados
        self.log("\n[7/7] Gerando gráficos (feche cada gráfico para ver o próximo)...")
        intensidade = os.path.join(output_dir, "MODEL.NIQ")
        raio_de_giro = os.path.join(output_dir, "RG_DMAX.DAT")
        por_file = os.path.join(output_dir, "MODEL.POR")
        
        try:
            plota_intensidade_vs_qRg(intensidade, raio_de_giro)
            plota_POR_vs_rRg(por_file, raio_de_giro)
            plota_intensidade(intensidade)
            plot_POR(por_file)
        except FileNotFoundError as e:
            self.log(f"ERRO: Não foi possível encontrar um dos arquivos de resultado para plotar: {e.filename}")
            return

        self.log("\nPROCESSO CONCLUÍDO COM SUCESSO!")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    janela = AppJanela()
    janela.show()
    sys.exit(app.exec_())