import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.lines as mlines

# Ler os dados capturados na ferramenta online WebPlotDigitizer
def read_data(file_name):
    x, y = [], []
    with open(file_name, 'r') as file:
        for line in file:
            line = line.replace(',', '.')
            values = line.split()
            if len(values) == 2:
                x.append(float(values[0]))
                y.append(float(values[1]))
    return np.array(x), np.array(y)

# Carregar dados do arquivo
x_data, y_data = read_data('data_.txt')

# Normalizar o eixo x dos dados experimentais
x_min = np.min(x_data)
x_max = np.max(x_data)
x_norm = (x_data - x_min) / (x_max - x_min)  # Normalização de x

# Função para ajustar duas distribuições gaussianas
def gaussian(x, amplitude1, mean1, stddev1, amplitude2, mean2, stddev2):
    return ((amplitude1 * (1/stddev1 * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean1) / stddev1)**2)) + \
           ((amplitude2 * (1/stddev2 * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean2) / stddev2)**2))

# Estimativas iniciais para o ajuste das gaussianas
initial_guess = [1, 0.3, 0.1, 1, 0.7, 0.1]

# Ajustar as gaussianas aos dados
params, _ = curve_fit(gaussian, x_norm, y_data, p0=initial_guess)

# Obter os parâmetros das gaussianas ajustadas
amplitude1, mean1, stddev1, amplitude2, mean2, stddev2 = params

# Gerar os valores ajustados da função gaussiana
y_gaussian = gaussian(x_norm, *params)

# Criar o gráfico principal com eixo x secundário
fig, ax1 = plt.subplots(figsize=(8, 5))
ax2 = ax1.twiny()

# Plotar os dados experimentais
ax1.plot(x_norm, y_data, 'o-', color="black", markersize=4, label="Dados Experimentais")
ax1.plot(x_norm, y_gaussian, '--', color="gray", label="Ajuste Gaussiano")

# Preencher áreas sob as gaussianas
ax1.fill_between(x_norm, 0, ((amplitude1 * (1/stddev1 * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x_norm - mean1) / stddev1)**2)),
                 color="red", alpha=0.8)
ax1.fill_between(x_norm, 0, ((amplitude2 * (1/stddev2 * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x_norm - mean2) / stddev2)**2)),
                 color="blue", alpha=0.5)

# Destacar a região azul
highlight_region = (x_norm >= 0.786) & (x_norm <= 1)
ax1.fill_between(x_norm[highlight_region], 0, (amplitude2 * (1/stddev2 * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x_norm[highlight_region] - mean2) / stddev2)**2),
                 color="blue")

# Adicionar linha vertical no meio das gaussianas
middle_x = (mean1 + mean2) / 2
ax1.axvline(x=middle_x, color='green', linestyle='--', label='$x_{med}$')

# Configuração do eixo secundário (valores reais de x)
ax2.set_xlim(x_min, x_max)

# Sincronizar ticks entre os dois eixos
real_ticks = np.linspace(x_min, x_max, 6)
norm_ticks = (real_ticks - x_min) / (x_max - x_min)
ax1.set_xticks(norm_ticks)
ax2.set_xticks(real_ticks)

# Adicionar setas e anotações
ax1.annotate('', xy=(0.786, 0.4), xytext=(middle_x, 0.4),
             arrowprops=dict(facecolor='blue', edgecolor='blue', alpha=0.5, arrowstyle='<-', lw=3))
ax1.text(0.61, 0.42, "$T_P$$^{baixo}$", color='blue', alpha=0.5, fontsize=18)

ax1.annotate('', xy=(1, 0.4), xytext=(0.786, 0.4),
             arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', lw=3))
ax1.text(0.83, 0.42, "$T_P$$^{alto}$", color='blue', fontsize=18)

ax1.annotate('', xy=(0, 0.4), xytext=(middle_x, 0.4),
             arrowprops=dict(facecolor='red', edgecolor='red', alpha=0.8, arrowstyle='<->', lw=3))
ax1.text(0.22, 0.42, "$T_N$", color='red', alpha=0.8, fontsize=18)

# Configurações dos eixos
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 0.5)

# Adicionar rótulos aos eixos
ax1.set_ylabel("PDF", fontsize=18)
ax1.set_xlabel("Expressão de Antígeno $x_{norm}$", fontsize=18)
ax2.set_xlabel("Expressão de Antígeno $x_{data}$", fontsize=18)

# Ajustes finais e exibição do gráfico
plt.tight_layout()
plt.grid(True, linestyle=':')
ax1.tick_params(axis='both', which='major', labelsize=15)
ax2.tick_params(axis='both', which='major', labelsize=15)

plt.show()
