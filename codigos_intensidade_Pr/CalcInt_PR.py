import numpy as np 
import os
import time
import matplotlib.pyplot as plt
import numba
# Tenta importar o Numba. Se não estiver instalado, o código ainda funciona, mas mais lentamente.
try:
    from numba import jit
    print("Numba detectado. Funções de cálculo serão aceleradas.")
except ImportError:
    print("AVISO: Numba não encontrado. O código funcionará em modo Python puro (mais lento).")
    # Cria um "decorador falso" se o Numba não estiver disponível
    def jit(func=None, **kwargs):
        if func:
            return func
        else:
            return lambda f: f

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
                    # CORREÇÃO: Usando as colunas padrão do formato PDB
                    # X: 31-38, Y: 39-46, Z: 47-54, Ocupância: 55-60
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
        print(f'ERRO: Ficheiro não encontrado em {filepath}'); return None, None
    except Exception as e:
        print(f'ERRO ao ler o PDB: {e}'); return None, None
    if not coords:
        print('AVISO: Nenhuma coordenada de átomo encontrada no ficheiro.'); return None, None
    return np.array(coords), np.array(weights)


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

@jit(nopython=True, fastmath=True)
def create_distance_histogram(coords, weights, d_max, num_bins=20000):
    """Constrói o histograma de distâncias par-a-par de forma eficiente com um laço duplo explícito, ideal para Numba."""
    n_atom = len(coords)
    histogram = np.zeros(num_bins, dtype=np.float64)
    f_quadrado = np.sum(weights**2)
    if d_max == 0: return histogram, f_quadrado
    fator_de_conversao = num_bins / d_max
    for i in range(n_atom):
        for j in range(i + 1, n_atom):
            dist = np.sqrt((coords[i, 0] - coords[j, 0])**2 + 
                             (coords[i, 1] - coords[j, 1])**2 + 
                             (coords[i, 2] - coords[j, 2])**2)
            bin_index = int(fator_de_conversao * dist)
            if bin_index < num_bins:
                histogram[bin_index] += weights[i] * weights[j]
    return histogram, f_quadrado

@jit(nopython=True, fastmath=True)
def calc_intensidade(q_values, histogram, f_quadrado, d_max, raio_do_atomo, num_bins=20000):
    """Calcula a intensidade "bruta" I(q) usando o histograma."""
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
    """Implementa p2(r) de Glatter para R1=R2=R. Corresponde a p0(r)."""
    pi = np.pi
    R2 = R*R
    R3 = R2*R
    return (4.0*pi/3.0)*R3*r**2 - pi*R2*r**3 + (pi/12.0)*r**5  #check

@jit(nopython=True, fastmath=True)
def p3_simplified(r, a, R):
    """Implementa a fórmula p3(r) de Glatter para R1=R2=R."""
    pi = np.pi
    R2, R3, R5 = R*R, R*R*R, R*R*R*R*R
    if a == 0: return 0.0
    
    term1 = (2.0/3.0) * R3 * ((2*R)**2 - (a - r)**2) #check
    term2 = -(2.0/6.0) * R2 * ((2*R)**3 - abs(a - r)**3) #check
    term3 = (1.0/60.0) * ((2*R)**5 - abs(a - r)**5) #check
    return (pi * r / (2.0 * a)) * (term1 + term2 + term3) #check

@jit(nopython=True, fastmath=True)
def p4_simplified(r, a, R):
    """Implementa a fórmula p4(r) de Glatter para R1=R2=R."""
    pi = np.pi
    R3, R5 = R*R*R, R*R*R*R*R
    if a == 0: return 0.0

    term1 = (8.0/15.0) * R5 #check
    term2 = -(2.0/3.0) * R3 * (a - r)**2 # check
    return (pi * r / (2.0 * a)) * (term1 + term2) #check

@jit(nopython=True, fastmath=True)
def p5_simplified(r, a, R):
    """Implementa a fórmula p5(r) de Glatter para R1=R2=R."""
    pi = np.pi
    R2, R3 = R*R, R*R*R
    if a == 0: return 0.0

    term1 = (8.0/3.0) * R3 * a * r #check
    term2 = -(2.0/6.0) * R2 * ((a + r)**3 - abs(a - r)**3) #check
    term3 = (1.0/60.0) * ((a + r)**5 - abs(a - r)**5) #check
    return (pi * r / (2.0 * a)) * (term1 + term2 + term3) #check

@jit(nopython=True, fastmath=True)
def p6_simplified(r, a, R):
    """Implementa a fórmula p6(r) de Glatter para R1=R2=R."""
    pi = np.pi
    R2, R3, R5 = R*R, R*R*R, R*R*R*R*R
    if a == 0: return 0.0

    term1 = -(2.0/3.0) * R3 * (a - r)**2 #check
    term2 = (2.0/3.0) * R3 * (r + a)**2  #check
    term3 = -(2.0/6.0) * R2 * (r + a)**3 # check
    term4 = (1.0/60.0) * (r + a)**5  #check
    return (pi * r / (2.0 * a)) * (term1 + term2 + term3 + term4) #check

@jit(nopython=True, fastmath=True)
def calculate_pr_analytical(histogram, f_quadrado, d_max, R, d_max_for_pr, num_bins=20000, nr=501):
    r_axis = np.linspace(0, d_max_for_pr, nr)
    pr_values = np.zeros(nr, dtype=np.float64)
    
    # 1. Contribuição da auto-correlação (p0(r), que é igual a p2(r) simplificada)
    r_max_p0 = 2 * R
    for k in range(nr):
        r = r_axis[k]
        if r < r_max_p0:
            pr_values[k] += p2_simplified(r, R) * f_quadrado

    # 2. Contribuição dos pares (correlação-cruzada)
    bin_width = d_max / num_bins
    r_hist = (np.arange(num_bins) + 0.5) * bin_width
    
    for i in range(num_bins):
        hist_val = histogram[i]
        if hist_val == 0: continue
        
        a = r_hist[i]
        
        for k in range(nr):
            r = r_axis[k]
            p_cross = 0.0
            
            if a > 0 and a <= R: # Casos II
                if r <= 2*R - a:
                    p_cross = p6_simplified(r, a, R)
                elif r <= a:
                    p_cross = p4_simplified(r, a, R)
                elif r <= 2*R + a:
                    p_cross = p3_simplified(r, a, R)
            
            elif a > R and a <= 2*R: # Caso VI
                if r <= 2*R - a:
                    p_cross = p5_simplified(r, a, R)
                elif r <= a:
                    p_cross = p3_simplified(r, a, R)
                elif r <= 2*R + a:
                    p_cross = p3_simplified(r, a, R)

            elif a > 2*R: # Caso VII
                if r >= a - 2*R and r <= 2*R + a:
                    p_cross = p3_simplified(r, a, R)

            pr_values[k] += 2 * p_cross * hist_val
            
    pr_values[r_axis > d_max] = 0
    return r_axis, pr_values

def plota_intensidade(q_values, intensities, model_name=""):
    """Plota a curva de I(q) vs q."""
    plt.figure(figsize=(10, 7))
    plt.plot(q_values, intensities, label=f'Modelo: {model_name}', color='blue')
    plt.title('Intensidade de Espalhamento Normalizada', fontsize=16)
    plt.xlabel(r'$Log(q) \quad (\AA^{-1})$', fontsize=14)
    plt.ylabel(r'$Log$ I(q) $\quad [Arb.U]$', fontsize=14)
    plt.xscale("log"); plt.yscale("log")
    plt.grid(True, which='both', linestyle='--', alpha=0.6)
    plt.legend(); plt.tight_layout()

def plota_pr(r_vals, pr_vals, model_name=""):
    """Plota a curva P(r) vs. r."""
    plt.figure(figsize=(10, 7))
    if np.max(pr_vals) > 0:
        pr_vals_normalized = pr_vals / np.max(pr_vals)
    else:
        pr_vals_normalized = pr_vals
    plt.plot(r_vals, pr_vals_normalized, label=f'Modelo: {model_name}', color='red')
    plt.title("Função de Distribuição de Distâncias P(r)", fontsize=16)
    plt.xlabel(r'$r \quad (\AA)$', fontsize=14)
    plt.ylabel(r'$P(r)$', fontsize=14)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(); plt.tight_layout()
def plota_histograma(histogram,dmax, model_name, numbins = 20000):
    bin_width = dmax / numbins
    r_hist= (np.arange(numbins) + 0.5) * bin_width #pegando o centro do bin,ou o ponto médio
    plt.figure(figsize=(10,7))
    plt.plot(r_hist,histogram, label=f'Modelo: {model_name}', color='green')
    plt.title('Histograma de distância',fontsize=14)
    plt.xlabel(r'Distância entre pares $r\quad(\AA)$', fontsize=14)
    plt.ylabel('Frequência ponderada pelos pesos',fontsize=14)
    plt.legend()
    plt.tight_layout()

class ModelParameters:
    def __init__(self):
        self.atom_coords, self.weights, self.beads = None, None, 0
        self.raio_atomo = None
        self.hist_dis, self.f_quadrado = None, None
        self.rg, self.d_max = None, None
        self.intensities = None
        self.r_values_for_pr, self.pr_values = None, None

def main():
    Model = ModelParameters()
    
    q_max_default = 1
    n_points = 1000
    q_vals = np.linspace(q_max_default / n_points, q_max_default, n_points)
    
    raio_subunidade = None
    try:
        entrada_subunidade = int(input('Como deseja definir o raio da subunidade?\n  1 - Informar o raio diretamente\n  2 - Calcular a partir de um q_max\n  ---> '))
        if entrada_subunidade == 1:
            raio_subunidade = float(input("Informe o raio em Angstroms -> "))
        elif entrada_subunidade == 2:
            q_max = float(input('Informe o máximo valor de q (q_max) em A-1: '))
            raio_subunidade = define_r(q_max)
        else:
            print("Opção inválida. A usar raio padrão de 1.0 A.")
            raio_subunidade = 1.0
    except ValueError:
        print("Entrada inválida. A usar raio padrão de 1.0 A.")
        raio_subunidade = 1.0

    path_pdb = r"C:\Interface-python\teste_criacao_modelo_py\Model.pdb"
    
    
    coords, weights = ler_pdb(path_pdb)
    if coords is None:
        print("Encerrando o programa devido a erro na leitura do PDB.")
        return

    Model.atom_coords, Model.weights, Model.raio_atomo = coords, weights, raio_subunidade
    print(f'\nFicheiro lido com sucesso! {len(Model.atom_coords)} átomos encontrados.')
    print(coords, weights)
    
    print('\n-------- A Calcular Propriedades do Modelo ------')
    start_time = time.time()
    rg, dmax = calc_rg_dmax(Model.atom_coords, Model.weights, Model.raio_atomo)
    end_time = time.time()
    Model.rg, Model.d_max = rg, dmax
    print(f"   Tempo de execução: {end_time - start_time:.4f} segundos")
    print(f'   Raio de giro (Rg) calculado: {Model.rg:.3f} A')
    print(f'   Dimensão máxima (Dmax) calculada: {Model.d_max:.3f} A')

    print("\n-------- A Gerar o Histograma de Distâncias ---------")
    start_time = time.time()
    hist, f_quadrado = create_distance_histogram(Model.atom_coords, Model.weights, Model.d_max)
    end_time = time.time()
    Model.hist_dis, Model.f_quadrado = hist, f_quadrado
    print(f"   Histograma construído com sucesso!")
    print(f"   Tempo de execução: {end_time - start_time:.4f} segundos")
    
    print('\n------Calculo da intensidade---------')
    start_time = time.time()
    intensidades_brutas = calc_intensidade(q_vals, Model.hist_dis, Model.f_quadrado, Model.d_max, Model.raio_atomo)
    end_time = time.time()
    
    if intensidades_brutas.size > 0 and intensidades_brutas[0] != 0:
        Model.intensities = intensidades_brutas / intensidades_brutas[0]
    else:
        Model.intensities = intensidades_brutas
    
    print(f'   Cálculo de intensidade finalizado!')
    print(f'   Tempo de execução: {end_time - start_time:.4f} segundos')
    # Define o nome da pasta e do arquivo
    intensity_folder = 'C:\\Interface-python\\teste_criacao_modelo_py\\INTPOR.PY'
    intensity_filename = 'intensidade.dat'

    # Cria a pasta (se ela não existir)
    os.makedirs(intensity_folder, exist_ok=True)

    # Cria o caminho completo para o arquivo
    intensity_filepath = os.path.join(intensity_folder, intensity_filename)

    # Salva o arquivo no caminho especificado
    # Lembre-se de usar sua variável de dados correta no lugar de 'dados_de_intensidade'
    np.savetxt(intensity_filepath,np.c_[q_vals, Model.intensities] , fmt='%.6e')

    print(f"Arquivo salvo em: {intensity_filepath}")
    
    # --- CÁLCULO E PLOTAGEM DA P(r) ANALÍTICA ---
    print('\n-------- A Calcular P(r) (Método Analítico) --------')
    try:
        dmax_pr_input = float(input("Digite a Dmax para o gráfico P(r) (ou <0 para auto): "))
    except ValueError:
        print("   Entrada inválida. A usar Dmax automática.")
        dmax_pr_input = -1.0
    dmax_for_pr_calc = dmax if dmax_pr_input < 0 else dmax_pr_input

    start_time = time.time()
    r_axis, pr_vals = calculate_pr_analytical(
        Model.hist_dis, Model.f_quadrado, Model.d_max, Model.raio_atomo, dmax_for_pr_calc
    )
    end_time = time.time()
    Model.r_values_for_pr, Model.pr_values = r_axis, pr_vals
    print(f"   Cálculo da P(r) finalizado!")
    print(f"   Tempo de execução: {end_time - start_time:.4f} segundos")
    # Define o nome da pasta e do arquivo
    pr_folder = 'C:\\Interface-python\\teste_criacao_modelo_py\\INTPOR.PY'
    pr_filename = 'pr_function_analytical.dat'

    # Cria a pasta (se ela não existir)
    os.makedirs(pr_folder, exist_ok=True)

    # Cria o caminho completo para o arquivo (ex: 'p(r)/pr_function_analytical.dat')
    pr_filepath = os.path.join(pr_folder, pr_filename)

    # Salva o arquivo no caminho especificado
    np.savetxt(pr_filepath, np.c_[r_axis, pr_vals], fmt='%.6e')

    print(f"Arquivo salvo em: {pr_filepath}")


    # --- VISUALIZAÇÃO FINAL ---
    print("\n-------- A Gerar Gráficos --------")
    plota_intensidade(q_vals, Model.intensities, os.path.basename(path_pdb))
    plota_pr(Model.r_values_for_pr, Model.pr_values, os.path.basename(path_pdb))
    plota_histograma(Model.hist_dis, Model.d_max, model_name='')
    plt.show()

    print('\nFim do programa.')

if __name__ == '__main__':
    main()
