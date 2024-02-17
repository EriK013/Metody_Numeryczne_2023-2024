import numpy as np
import matplotlib.pyplot as plt
N = 124
max_iterations = list(range(1, 200))
b = np.arange(1, N + 1)

def Dokladne_rozwiazania(N, b):
    AP = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i == j:
                AP[i, j] = 3
            elif j == i + 1:
                AP[i, j] = 1
            elif i == j + 1:
                AP[i, j] = 1
            elif j == i + 2:
                AP[i, j] = 0.15
            elif i == j + 2:
                AP[i, j] = 0.15
    return np.linalg.solve(AP,b)


bP = Dokladne_rozwiazania(N,b)


def A(i,j):  # Funkcja przechowująca wartości macierzy A
    if (i == j):
        return 3
    elif np.abs(i - j) == 1:
        return 1
    elif np.abs(i - j) == 2:
        return 0.15
    else:
        return 0

# Funkcja rozwiązująca układ równań za pomocą metody jacobiego
def Jacobi(N, x0, estop=1e-15, max_iter=400):
    x = x0.copy() # Wektor rozwiązań
    x_prev = x.copy()  # Wektor wartości poprzednich iteracji
    l = len(x)
    errors = []

    for k in range(max_iter):
        for i in range(N):
            # Obliczamy sumę aijxj(k) z wyjeciem j == i
            sum = 0
            for j in range(max(0, i-2),min(l,i+3)):
                if j != i:
                    sum += A(i,j) * x_prev[j]
            x[i] = (b[i] - sum)/A(i,i)

        error = np.linalg.norm(bP - x)
        errors.append(error)
        if error < estop:
            break

        x_prev = x.copy()

    return x, errors


def GaussSeidel(N, x0, estop=1e-15, max_iter=400):
    x = x0.copy()
    x_prev = x.copy()
    errors = []

    for k in range(max_iter):
        for i in range(N):
            sum = 0
            for j in range(max(0,i-2), min(len(x), i + 3)):
                if j < i:
                    sum += A(i,j) * x[j]
                elif j > i:
                    sum += A(i,j) * x_prev[j]
            x[i] = (b[i] - sum)/A(i,i)

        error = np.linalg.norm(bP - x)
        errors.append(error)
        if error < estop:
            break

        x_prev = x.copy()

    return x, errors


# Rozwiązanie
def solution(x0, title=""):
    x_jacobi, errors_jacobi = Jacobi(N, x0)
    x_gs, errors_gauss = GaussSeidel(N, x0)
    print("#######Jacobi#########")
    print(x_jacobi)
    print("#######GaussSeidel#########")
    print(x_gs)
    print("################")
    print("Poprawne rozwiązania")
    print(bP)

    plt.figure(figsize=(10, 6))
    plt.plot(errors_jacobi, label="Jacobi")
    plt.plot(errors_gauss, label="Gauss")
    plt.yscale('log')
    plt.xlabel('Iteracje')
    plt.ylabel('Błąd')
    plt.title(f"Błędy w kolejnych iteracjach {title}")
    plt.legend()
    #plt.savefig("Error_t4(N).svg")
    plt.show()

solution(np.zeros(N))

# Testowanie różnych punktów startowych
solution(np.add(np.zeros(N), 1e15), "[0+1e15,..., 0+1e15]")
solution(np.add(np.zeros(N), -1e7), "[0-1e7,..., 0-1e7]")
solution(np.ones(N), "punkt startowego [1,...,1]")
solution(np.add(np.zeros(N), 100000), "[100000,..., 100000]")