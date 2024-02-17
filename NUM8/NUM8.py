import numpy as np
import matplotlib.pyplot as plt

# Pobranie danych z pliku
def getData(filename="data.txt"):
    data = np.loadtxt(filename, delimiter=',')
    x_values = data[:, 0]
    y_values = data[:, 1]
    return x_values, y_values

# Funkcja modelująca punkty
def F(x, parameters):
    return parameters[0]* x**2 + parameters[1] * np.sin(x) + parameters[2] * np.cos(5*x) + parameters[3] * np.exp(-x)
def G(x, parameters):
    return parameters[0]*(x-np.sin(x))**2 + parameters[1] * np.exp(np.sin(x)) + parameters[2] * np.sin(x) + parameters[3]*x

def FStack(x):
    return np.column_stack([x ** 2, np.sin(x), np.cos(5 * x), np.exp(-x)])

def GStack(x):
    return np.column_stack([(x-np.sin(x)**3), x**2, np.sin(x), x])

# Prawdziwy przebieg funkcji G
def realG(x, g_realParameters):
    plt.figure(figsize=(10, 6))
    plt.title("Prawdziwy przebieg funkcji G")
    plt.plot(x, G(x, g_realParameters), color="orange")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.savefig("Prawdziwy przebieg funkcji G.svg")
    plt.show()

# Metoda najmniejszych kwadratów za pomocą SVD
def findParameters(x,y, vsatckfunction):
    X = vsatckfunction(x)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    parameters = Vt.T @ np.linalg.inv(np.diag(S)) @ U.T @ y
    return parameters

def drawPlot(f, x, y, p, title):
    plt.figure(figsize=(10,6))
    plt.title(title)
    plt.scatter(x, y, color="blue")
    plt.plot(x, f(x,p), color="orange")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.savefig(f"{title}.svg")
    plt.show()

def main():
    x_values, y_values = getData()
    parameters = findParameters(x_values, y_values, FStack)
    print(f"Parametry funkcji F: {parameters}")
    drawPlot(F, x_values, y_values, parameters, "Dopasowana funkcja F")
    g_realParameters = [-0.025, 0.75, 12, -0.75]
    realG(x_values, g_realParameters)
    y_errorG = G(x_values, g_realParameters) + np.random.normal(0, 1e-4, size=len(x_values))
    g_parameters = findParameters(x_values, y_errorG, GStack)
    print(f"Parametry funkcji G: {g_parameters}")
    drawPlot(G, x_values, y_errorG, g_parameters, "Dopasowana funkcja G")


main()