import sys
import os
import time
import numpy as np
import matplotlib.pyplot as plt

# Tenta importar o Numba. Se não estiver instalado, o código ainda funciona, mas mais lentamente.
try:
    from numba import jit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    # Cria um "decorador falso" se o Numba não estiver disponível
    def jit(func=None, **kwargs):
        if func:
            return func
        else:
            return lambda f: f

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QFileDialog, QLabel, QTextEdit, QGroupBox, QFormLayout, QDoubleSpinBox,
    QRadioButton, QMessageBox
)
from PyQt5.QtCore import QObject, QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# ==============================================================================
# SUAS FUNÇÕES DE CÁLCULO ORIGINAIS (sem modificações)
# ==============================================================================

def define_r(q_max):
    """Calcula um raio de subunidade razoável a partir de um q_max."""
    if q_max > 0:
        return 4.0 / q_max
    return 1.0

def ler_pdb(filepath):
    """Lê um ficheiro PDB usando fatiamento de colunas fixas e padrão."""
    coords, weights = [], []
    try:
        with open(filepath, 'r') as pdb_file:
            for linha in pdb_file:
                if linha.startswith(('ATOM', 'HETATM')):
                    x = float(linha[30:38].strip())
                    y = float(linha[38:46].strip())
                    z = float(linha[46:54].strip())
                    try:
                        weight = float(linha[54:63].strip())
                    except (ValueError, IndexError):
                        weight = 1.0
                    weights.append(weight)
                    coords.append([x, y, z])
    except FileNotFoundError:
        return None, None, f'ERRO: Ficheiro não encontrado em {filepath}'
    except Exception as e:
        return None, None, f'ERRO ao ler o PDB: {e}'
    if not coords:
        return None, None, 'AVISO: Nenhuma coordenada de átomo encontrada no ficheiro.'
    return np.array(coords), np.array(weights), None # Success

@jit(nopython=True, fastmath=True)
def calc_rg_dmax(coords, weight, radius):
    """Calcula Rg e Dmax, incluindo a correção do Teorema dos Eixos Paralelos."""
    n_atoms = len(coords)
    if n_atoms < 2: return 0.0, 0.0
    xs, ys, zs, wsum = 0.0, 0.0, 0.0, 0.0
    for j in range(n_atoms):
        xs += coords[j, 0] * weight[j]
        ys += coords[j, 1] * weight[j]
        zs += coords[j, 2] * weight[j]
        wsum += weight[j]
    if wsum == 0: return 0.0, 0.0
    x_cm, y_cm, z_cm = xs / wsum, ys / wsum, zs / wsum
    acumular_dist = 0.0
    for k in range(n_atoms):
        dx, dy, dz = coords[k, 0] - x_cm, coords[k, 1] - y_cm, coords[k, 2] - z_cm
        acumular_dist += (dx*dx + dy*dy + dz*dz) * weight[k]
    rg_sq = (0.6) * radius**2 + acumular_dist / wsum
    rg = np.sqrt(rg_sq)
    dmax2 = 0.0
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            dist2 = ((coords[i, 0] - coords[j, 0])**2 +
                     (coords[i, 1] - coords[j, 1])**2 +
                     (coords[i, 2] - coords[j, 2])**2)
            if dist2 > dmax2:
                dmax2 = dist2
    Dmax = np.sqrt(dmax2)
    return rg, Dmax

# ... (todas as outras funções de cálculo: create_distance_histogram, calc_intensidade, p2, p3, etc. permanecem aqui) ...
# O restante das suas funções de cálculo vai aqui, sem alterações.
# Para manter a resposta concisa, elas foram omitidas, mas você deve colá-las aqui.
@jit(nopython=True, fastmath=True)
def create_distance_histogram(coords, weights, d_max, num_bins=20000):
    n_atom = len(coords)
    histogram = np.zeros(num_bins, dtype=np.float64)
    f_quadrado = np.sum(weights**2)
    if d_max == 0: return histogram, f_quadrado
    fator_de_conversao = num_bins / d_max
    for i in range(n_atom):
        for j in range(i + 1, n_atom):
            dist = np.sqrt((coords[i, 0] - coords[j, 0])**2 + (coords[i, 1] - coords[j, 1])**2 + (coords[i, 2] - coords[j, 2])**2)
            bin_index = int(fator_de_conversao * dist)
            if bin_index < num_bins:
                histogram[bin_index] += weights[i] * weights[j]
    return histogram, f_quadrado

@jit(nopython=True, fastmath=True)
def calc_intensidade(q_values, histogram, f_quadrado, d_max, raio_do_atomo, num_bins=20000):
    intensities = np.zeros_like(q_values, dtype=np.float64)
    bin_width = d_max / num_bins
    r_j = (np.arange(num_bins) + 0.5) * bin_width
    for i in range(len(q_values)):
        q = q_values[i]
        x_form = q * raio_do_atomo
        if x_form < 0.05:
            fator_forma = 1.0 - 0.1 * x_form**2
        else:
            fator_forma = 3 * (np.sin(x_form) - x_form * np.cos(x_form)) / x_form**3
        qr_j = q * r_j
        sinc_values = np.ones_like(qr_j)
        for k in range(len(qr_j)):
            if qr_j[k] != 0:
                sinc_values[k] = np.sin(qr_j[k]) / qr_j[k]
        hist_sum = np.dot(histogram, sinc_values)
        fator_estrutura = f_quadrado + 2 * hist_sum
        intensities[i] = (fator_forma**2) * fator_estrutura
    return intensities

@jit(nopython=True, fastmath=True)
def p2_simplified(r, R):
    pi = np.pi; R2 = R*R; R3 = R2*R
    return (4.0*pi/3.0)*R3*r**2 - pi*R2*r**3 + (pi/12.0)*r**5

@jit(nopython=True, fastmath=True)
def p3_simplified(r, a, R):
    pi = np.pi; R2, R3 = R*R, R*R*R
    if a == 0: return 0.0
    term1 = (2.0/3.0) * R3 * ((2*R)**2 - (a - r)**2)
    term2 = -(2.0/6.0) * R2 * ((2*R)**3 - abs(a - r)**3)
    term3 = (1.0/60.0) * ((2*R)**5 - abs(a - r)**5)
    return (pi * r / (2.0 * a)) * (term1 + term2 + term3)

@jit(nopython=True, fastmath=True)
def p4_simplified(r, a, R):
    pi = np.pi; R3, R5 = R*R*R, R*R*R*R*R
    if a == 0: return 0.0
    term1 = (8.0/15.0) * R5
    term2 = -(2.0/3.0) * R3 * (a - r)**2
    return (pi * r / (2.0 * a)) * (term1 + term2)

@jit(nopython=True, fastmath=True)
def p5_simplified(r, a, R):
    pi = np.pi; R2, R3 = R*R, R*R*R
    if a == 0: return 0.0
    term1 = (8.0/3.0) * R3 * a * r
    term2 = -(2.0/6.0) * R2 * ((a + r)**3 - abs(a - r)**3)
    term3 = (1.0/60.0) * ((a + r)**5 - abs(a - r)**5)
    return (pi * r / (2.0 * a)) * (term1 + term2 + term3)

@jit(nopython=True, fastmath=True)
def p6_simplified(r, a, R):
    pi = np.pi; R2, R3 = R*R, R*R*R
    if a == 0: return 0.0
    term1 = -(2.0/3.0) * R3 * (a - r)**2
    term2 = (2.0/3.0) * R3 * (r + a)**2
    term3 = -(2.0/6.0) * R2 * (r + a)**3
    term4 = (1.0/60.0) * (r + a)**5
    return (pi * r / (2.0 * a)) * (term1 + term2 + term3 + term4)

@jit(nopython=True, fastmath=True)
def calculate_pr_analytical(histogram, f_quadrado, d_max, R, d_max_for_pr, num_bins=20000, nr=501):
    r_axis = np.linspace(0, d_max_for_pr, nr)
    pr_values = np.zeros(nr, dtype=np.float64)
    r_max_p0 = 2 * R
    for k in range(nr):
        r = r_axis[k]
        if r < r_max_p0:
            pr_values[k] += p2_simplified(r, R) * f_quadrado
    bin_width = d_max / num_bins
    r_hist = (np.arange(num_bins) + 0.5) * bin_width
    for i in range(num_bins):
        hist_val = histogram[i]
        if hist_val == 0: continue
        a = r_hist[i]
        for k in range(nr):
            r = r_axis[k]
            p_cross = 0.0
            if a > 0 and a <= R:
                if r <= 2*R - a: p_cross = p6_simplified(r, a, R)
                elif r <= a: p_cross = p4_simplified(r, a, R)
                elif r <= 2*R + a: p_cross = p3_simplified(r, a, R)
            elif a > R and a <= 2*R:
                if r <= 2*R - a: p_cross = p5_simplified(r, a, R)
                elif r <= a: p_cross = p3_simplified(r, a, R)
                elif r <= 2*R + a: p_cross = p3_simplified(r, a, R)
            elif a > 2*R:
                if r >= a - 2*R and r <= 2*R + a: p_cross = p3_simplified(r, a, R)
            pr_values[k] += 2 * p_cross * hist_val
    pr_values[r_axis > d_max] = 0
    return r_axis, pr_values


# ==============================================================================
# CLASSE WORKER: Executa os cálculos em uma thread separada para não travar a GUI
# ==============================================================================
class CalculationWorker(QObject):
    finished = pyqtSignal(object)  # Sinal emitido quando tudo termina, carrega os resultados
    progress = pyqtSignal(str)     # Sinal para enviar mensagens de texto para a GUI
    error = pyqtSignal(str)        # Sinal para reportar erros

    def __init__(self, pdb_path, raio_subunidade, dmax_pr):
        super().__init__()
        self.pdb_path = pdb_path
        self.raio_subunidade = raio_subunidade
        self.dmax_pr = dmax_pr

    def run(self):
        try:
            # --- PARTE 1: LEITURA E PREPARAÇÃO ---
            self.progress.emit(f"Lendo arquivo: {self.pdb_path}...")
            coords, weights, error_msg = ler_pdb(self.pdb_path)
            if error_msg:
                self.error.emit(error_msg)
                return

            model_name = os.path.basename(self.pdb_path)
            self.progress.emit(f"Arquivo lido com sucesso! {len(coords)} átomos encontrados.")
            self.progress.emit(f"Raio da subunidade definido como: {self.raio_subunidade:.3f} Å")

            # --- PARTE 2: CÁLCULOS PRINCIPAIS ---
            self.progress.emit("\n-------- A Calcular Propriedades do Modelo ------")
            start_time = time.time()
            rg, dmax = calc_rg_dmax(coords, weights, self.raio_subunidade)
            self.progress.emit(f"  Tempo de execução: {time.time() - start_time:.4f} segundos")
            self.progress.emit(f"  Raio de giro (Rg) calculado: {rg:.3f} Å")
            self.progress.emit(f"  Dimensão máxima (Dmax) calculada: {dmax:.3f} Å")

            self.progress.emit("\n-------- A Gerar o Histograma de Distâncias ---------")
            start_time = time.time()
            hist, f_quadrado = create_distance_histogram(coords, weights, dmax)
            self.progress.emit(f"  Histograma construído com sucesso!")
            self.progress.emit(f"  Tempo de execução: {time.time() - start_time:.4f} segundos")

            self.progress.emit('\n------ Calculo da intensidade I(q) ---------')
            q_vals = np.linspace(0.4 / 1000, 0.4, 1000)
            start_time = time.time()
            intensidades_brutas = calc_intensidade(q_vals, hist, f_quadrado, dmax, self.raio_subunidade)
            intensities = intensidades_brutas / intensidades_brutas[0] if intensidades_brutas.size > 0 and intensidades_brutas[0] != 0 else intensidades_brutas
            self.progress.emit(f"  Cálculo de intensidade finalizado!")
            self.progress.emit(f"  Tempo de execução: {time.time() - start_time:.4f} segundos")
            
            # --- PARTE 3: CÁLCULO P(r) ---
            self.progress.emit('\n-------- A Calcular P(r) (Método Analítico) --------')
            dmax_for_pr_calc = dmax if self.dmax_pr <= 0 else self.dmax_pr
            self.progress.emit(f"  Usando Dmax para P(r) de: {dmax_for_pr_calc:.3f} Å")
            start_time = time.time()
            r_axis, pr_vals = calculate_pr_analytical(hist, f_quadrado, dmax, self.raio_subunidade, dmax_for_pr_calc)
            self.progress.emit(f"  Cálculo da P(r) finalizado!")
            self.progress.emit(f"  Tempo de execução: {time.time() - start_time:.4f} segundos")

            self.progress.emit("\nCálculos concluídos. Gerando gráficos...")
            
            results = {
                'q_vals': q_vals, 'intensities': intensities,
                'r_axis': r_axis, 'pr_vals': pr_vals,
                'model_name': model_name
            }
            self.finished.emit(results)

        except Exception as e:
            self.error.emit(f"Ocorreu um erro inesperado durante o cálculo: {e}")

# ==============================================================================
# CLASSE DA JANELA PRINCIPAL (GUI)
# ==============================================================================
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Calculadora de Espalhamento a partir de PDB")
        self.setGeometry(100, 100, 800, 600)
        self.selected_pdb_path = None

        # --- Layout Principal ---
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)

        # --- Painel Esquerdo (Controles) ---
        controls_layout = QVBoxLayout()
        main_layout.addLayout(controls_layout, 1) # 1/3 do espaço

        # --- Grupo de Seleção de Arquivo ---
        file_group = QGroupBox("1. Seleção de Arquivo")
        file_layout = QVBoxLayout()
        self.select_file_btn = QPushButton("Selecionar Arquivo PDB...")
        self.select_file_btn.clicked.connect(self.select_pdb_file)
        self.file_path_label = QLabel("Nenhum arquivo selecionado.")
        self.file_path_label.setWordWrap(True)
        file_layout.addWidget(self.select_file_btn)
        file_layout.addWidget(self.file_path_label)
        file_group.setLayout(file_layout)
        controls_layout.addWidget(file_group)
        
        # --- Grupo de Parâmetros ---
        params_group = QGroupBox("2. Parâmetros de Cálculo")
        params_layout = QFormLayout()
        
        self.radio_raio_direto = QRadioButton("Informar raio diretamente (Å):")
        self.radio_raio_qmax = QRadioButton("Calcular raio a partir de q_max (Å⁻¹):")
        self.spin_raio_direto = QDoubleSpinBox(value=3.0, minimum=0.1, maximum=100.0, singleStep=0.1)
        self.spin_raio_qmax = QDoubleSpinBox(value=0.5, minimum=0.01, maximum=10.0, singleStep=0.05)
        
        self.radio_raio_direto.toggled.connect(self.update_raio_input_state)
        self.radio_raio_qmax.setChecked(True) # Padrão
        
        params_layout.addRow(self.radio_raio_direto, self.spin_raio_direto)
        params_layout.addRow(self.radio_raio_qmax, self.spin_raio_qmax)
        
        self.spin_dmax_pr = QDoubleSpinBox(value=0, minimum=0, maximum=1000.0, singleStep=10)
        params_layout.addRow("Dmax para P(r) (0 para auto):", self.spin_dmax_pr)
        params_group.setLayout(params_layout)
        controls_layout.addWidget(params_group)

        # --- Botão de Execução ---
        self.run_button = QPushButton("Iniciar Cálculo")
        self.run_button.setStyleSheet("font-size: 16px; padding: 10px; background-color: #4CAF50; color: white;")
        self.run_button.clicked.connect(self.start_calculation)
        self.run_button.setEnabled(False) # Desabilitado até um arquivo ser selecionado
        controls_layout.addWidget(self.run_button)

        controls_layout.addStretch() # Empurra tudo para cima

        # --- Painel Direito (Log) ---
        log_group = QGroupBox("Log de Progresso")
        log_layout = QVBoxLayout()
        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setFontFamily("Courier")
        log_layout.addWidget(self.log_output)
        log_group.setLayout(log_layout)
        main_layout.addWidget(log_group, 2) # 2/3 do espaço
        
        if NUMBA_AVAILABLE:
            self.log_output.append("Numba detectado. Funções de cálculo serão aceleradas.")
        else:
            self.log_output.append("AVISO: Numba não encontrado. O código funcionará em modo Python puro (mais lento).")

    def update_raio_input_state(self):
        self.spin_raio_direto.setEnabled(self.radio_raio_direto.isChecked())
        self.spin_raio_qmax.setEnabled(self.radio_raio_qmax.isChecked())

    def select_pdb_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Selecione o arquivo PDB", "", "PDB Files (*.pdb);;All Files (*)")
        if path:
            self.selected_pdb_path = path
            self.file_path_label.setText(os.path.basename(path))
            self.run_button.setEnabled(True)
            self.log_output.append(f"\nArquivo selecionado: {path}")

    def start_calculation(self):
        if not self.selected_pdb_path:
            self.show_error("Nenhum arquivo PDB foi selecionado.")
            return

        # Coletar parâmetros da GUI
        if self.radio_raio_direto.isChecked():
            raio = self.spin_raio_direto.value()
        else:
            raio = define_r(self.spin_raio_qmax.value())
        
        dmax_pr = self.spin_dmax_pr.value()

        # Desabilitar botões para evitar cliques múltiplos
        self.run_button.setEnabled(False)
        self.run_button.setText("Calculando...")
        self.log_output.clear() # Limpa o log para a nova execução
        
        # Configurar e iniciar a thread de cálculo
        self.thread = QThread()
        self.worker = CalculationWorker(self.selected_pdb_path, raio, dmax_pr)
        self.worker.moveToThread(self.thread)

        # Conectar sinais
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.on_calculation_finish)
        self.worker.error.connect(self.on_calculation_error)
        self.worker.progress.connect(self.append_log)
        
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)

        self.thread.start()

    def append_log(self, message):
        self.log_output.append(message)
        self.log_output.verticalScrollBar().setValue(self.log_output.verticalScrollBar().maximum())

    def on_calculation_finish(self, results):
        self.append_log("\nProcesso finalizado com sucesso!")
        self.run_button.setEnabled(True)
        self.run_button.setText("Iniciar Cálculo")
        self.show_plots(results)
    
    def on_calculation_error(self, message):
        self.show_error(message)
        self.run_button.setEnabled(True)
        self.run_button.setText("Iniciar Cálculo")

    def show_error(self, message):
        msg_box = QMessageBox(self)
        msg_box.setIcon(QMessageBox.Critical)
        msg_box.setText("Erro")
        msg_box.setInformativeText(message)
        msg_box.setWindowTitle("Erro na Execução")
        msg_box.exec_()
        self.append_log(f"ERRO: {message}")
    
    def show_plots(self, results):
        self.plot_window = PlotWindow(results, self)
        self.plot_window.show()


class PlotWindow(QMainWindow):
    """Uma nova janela para exibir os gráficos Matplotlib."""
    def __init__(self, results, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Resultados para: {results['model_name']}")
        self.setGeometry(200, 200, 1000, 500)
        
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        layout = QHBoxLayout(main_widget)
        
        # --- Gráfico I(q) ---
        fig_iq = Figure()
        canvas_iq = FigureCanvas(fig_iq)
        ax_iq = fig_iq.add_subplot(111)
        ax_iq.plot(results['q_vals'], results['intensities'], label=f'Modelo: {results["model_name"]}', color='blue')
        ax_iq.set_title('Intensidade de Espalhamento Normalizada')
        ax_iq.set_xlabel(r'log q $ \quad (\AA^{-1})$')
        ax_iq.set_ylabel(r'Log  $I(q) \quad [Arb. U]$')
        ax_iq.set_xscale("log")
        ax_iq.set_yscale("log")
        ax_iq.grid(True, which='both', linestyle='--', alpha=0.6)
        ax_iq.legend()
        fig_iq.tight_layout()
        layout.addWidget(canvas_iq)

        # --- Gráfico P(r) ---
        fig_pr = Figure()
        canvas_pr = FigureCanvas(fig_pr)
        ax_pr = fig_pr.add_subplot(111)
        pr_vals = results['pr_vals']
        pr_norm = pr_vals / np.max(pr_vals) if np.max(pr_vals) > 0 else pr_vals
        ax_pr.plot(results['r_axis'], pr_norm, label=f'Modelo: {results["model_name"]}', color='red')
        ax_pr.set_title("Função de Distribuição de Distâncias P(r)")
        ax_pr.set_xlabel(r'$r \quad (\AA)$')
        ax_pr.set_ylabel(r'$P(r) $')
        ax_pr.axhline(0, color='black', linewidth=0.5)
        ax_pr.grid(True, linestyle='--', alpha=0.6)
        ax_pr.legend()
        fig_pr.tight_layout()
        layout.addWidget(canvas_pr)

# ==============================================================================
# PONTO DE ENTRADA DA APLICAÇÃO
# ==============================================================================
if __name__ == '__main__':
    # O seu 'main()' original não é mais necessário, a GUI controla o fluxo.
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())