import numpy as np
import scipy.linalg as sc
import timeit
import matplotlib.pyplot as plt


N = 80
times = [] # Przechowywanie czasu
sizes = range(0, 1000, 50)


'''
Korzystamy ze wzoru Shermana Morrisona.
x = y - ((v^t*y)*z/(1+v^t*z))
Macierz A możemy rozłożyć na A = uv^T + macierz rzadka(A'), 
która ma elementy niezerowe(11) na diagonali i (7) na superdiagonali
A'y = b
A'z = u
u i v są wektorami wypełnionymi na każdej pozycji 1
'''

################################
#         Rozwiązanie          #
################################
def build(N):
    global A, b, u, v, x
    b = []  # Wektor rozwiązań b
    for i in range(0, N):
        b.append(5)
    u = []  # Wektor u
    for i in range(0, N):
        u.append(1)
    v = []  # Transponowany wektor v
    for i in range(0, N):
        v.append(1)
    # Macierz A' zapisana jako 2N
    ku = 1  # Liczba superdiagonali
    A = [[0 for x in range(N)] for y in range(ku + 1)]
    x = []
    for i in range(0, N):
        x.append(0)

def fillA(N): # O(n)
    for i in range (N):
        # Diagonala
        A[1][i] = 11
        # Superdiagonala
        if (i < (N - 1)):
            A[0][i] = 7

def filly(N): # O(n)
    y = []
    for i in range(N):
        y.append(0)
    for i in range(N-1, -1, -1):
        if (i == N-1):
            y[i] = b[i] / A[1][i]
        else:
            y[i] = (b[i] - A[0][i] * y[i+1]) / A[1][i]
    return y

def fillz(N):
    z = []
    for i in range(N):
        z.append(0)
    for i in range(N - 1, -1, -1):
        if (i == N - 1):
            z[i] = u[i] / A[1][i]
        else:
            z[i] = (u[i] - A[0][i] * z[i + 1]) / A[1][i]
    return z

def printvector(a):
    for i in range(len(a)):
        print(a[i])

def solve(N):
    build(N)
    fillA(N)
    y = filly(N)
    z = fillz(N)
    licznik = 0
    mianownik = 0
    for i in range(len(v)):
        licznik += v[i] * y[i]
        mianownik += v[i] * z[i]
    mianownik += 1
    for j in range(len(y)):
        x[j] = (y[j] - (licznik/mianownik)*z[j])

solve(N)
print("Wektor x: ")
printvector(x)

################################
#TESTOWANIE ZA POMOCA BIBLIOTEK#
################################
# ROZWIAZANIE OGOLNE ZA POMOCA SOLVE
A = np.ones((N, N), dtype=int) # Stworzenie macierzy A
np.fill_diagonal(A, 12)
for i in range(N - 1):
    A[i, i + 1] = 8
b = np.full(N, 5, dtype=int) # Stworzenie wektora z samą wartością 5
b = np.array([b]).T
np.transpose(b)
print("Wektor x obliczony za pomocą bibliotek: ")
print(sc.solve(A, b))

################################
#       Obliczanie czasu        #
################################
for N in sizes:
    if (N < 20):
        czas = timeit.timeit(lambda: solve(N), number=400)
        sredni_czas = czas/500
    else:
        czas = timeit.timeit(lambda: solve(N), number=1)
        sredni_czas = czas
    times.append((N, round(sredni_czas, 10)))


wartosci_N, srednie_czasy = zip(*times)
# Rysowanie wykresu
plt.figure(figsize=(10, 6))
plt.plot(wartosci_N, srednie_czasy)
plt.xlabel('Rozmiar macierzy N')
plt.ylabel('Czas trwania operacji [s]')
plt.grid(True)
plt.savefig('wykres.svg', format='svg')
plt.show()