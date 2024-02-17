import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt

# Zadana macierz M
M = np.array([[8, 1, 0, 0],
              [1, 7, 2, 0],
              [0, 2, 6, 3],
              [0, 0, 3, 5]])

max_iterations = list(range(1,200))
tolerancja = 1e-15

# Obliczanie dokładnych wartości własnych
d_lambda, d_vectorlambda = la.eig(M)

# A
def ANajwiekszyModul(A, p=0):
    lambdy = []
    zbierznosci = []
    y_prev = np.array([1,0,0,0]) # Przykładowy wektor początkowy o normie 1
    for k in max_iterations:
        lambda_k = np.dot(y_prev.T, np.dot(A, y_prev)) / np.dot(y_prev.T, y_prev) + p
        y_k = np.dot(A, y_prev) / np.linalg.norm(np.dot(A, y_prev))

        zbierznosc = np.abs(np.max(d_lambda) - lambda_k)
        lambdy.append(lambda_k)

        y_prev = y_k/lambda_k

        zbierznosci.append(zbierznosc)
        if (zbierznosc < tolerancja):
            break

    print(f'"Największa co do modułu wartość własna: {lambdy[-1]}"')
    print(f'"Odpowiadający jej wektor własny {y_prev}"')
    # Wykres zbieżności w skali logarytmicznej
    plt.figure(figsize=(10, 6))
    plt.plot(max_iterations, zbierznosci, label='||Lambda_dokładna - Lambda_k||')
    plt.yscale('log')
    plt.xlabel('Liczba iteracji')
    plt.ylabel('Błąd wartości własnej')
    plt.title('Zbieżność metody potęgowej')
    plt.legend()
    plt.grid(True)
    plt.savefig('zbierznosc1.svg', format='svg')
    plt.show()
    return lambdy

# B
def BQR(A, p=0):
    diagonale = []
    diagonale.append(np.diag(A))
    suma_modulow = []
    for i in max_iterations:
        Q, R = np.linalg.qr(A)
        A = np.dot(R,Q)

        diagonale.append(np.diag(A))
        suma_modulu = np.sum(np.abs(np.tril(A, k=-1))) # Suma modułów elementów pod diagonalą
        suma_modulow.append(suma_modulu)
        if (suma_modulu < tolerancja):
            break


    print(f"Wartości własne macierzy M: {np.diag(A) + p}")
    # Wykres pokazujący ewolucję elementów diagonalnych w zależności od iteracji
    plt.figure(figsize=(10,6))
    plt.plot([a11[0] for a11 in diagonale], label=f'Pierwszy element na diagonali')
    plt.plot([a22[1] for a22 in diagonale], label=f'Drugi element na diagonali')
    plt.plot([a33[2] for a33 in diagonale], label=f'Trzeci element na diagonali')
    plt.plot([a44[3] for a44 in diagonale], label=f'Czwarty element na diagonali')
    plt.xlabel('Iteracja')
    plt.ylabel('Wartości elementów na diagonali')
    plt.title('Ewolucja elementów na diagonali')
    plt.legend()
    plt.grid(True)
    plt.savefig('diagonale.svg', format='svg')
    plt.show()
    # Wykres pokazujący zmianę elementów pod diagonalą w każdej iteracji
    plt.figure(figsize=(10, 6))
    plt.plot(suma_modulow, label=f'Suma modułów elementów pod diagonalą')
    plt.yscale('log')
    plt.xlabel('Iteracja')
    plt.ylabel('Wartości elementów pod diagonalą')
    plt.title('Zanikanie elementów pod diagonalą')
    plt.legend()
    plt.grid(True)
    plt.savefig('ewolucja.svg', format='svg')
    plt.show()
    return diagonale[-1]

def ApplyWilkison(A, lambd):
    p = 1/2*(lambd[1]+lambd[-1])
    A1 = A - p*np.eye(np.size(A[0]))
    return A1,p

print(f'"Wartości własne obliczone za pomocą funkcji bibliotecznej: {d_lambda}"')
ANajwiekszyModul(M)
l = BQR(M)
M1, p = ApplyWilkison(M,l)
ANajwiekszyModul(M1, p)
