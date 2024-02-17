import numpy as np
import scipy.linalg as sc
import timeit
import matplotlib.pyplot as plt

times = [] # Przechowywanie czasu
sizes = range(0, 1000, 50)

################################################################
#              Funkcje pomocnicze w rozwiazaniu
################################################################
# Uzupełnianie reprezentacji 4xN macierzy A
def fillA(N): # O(N)
    for i in range (N):
        # Diagonala
        A[1][i] = 1.2

        # Subdiagonala
        if (i < (N - 1)):
            A[0][i] = 0.2

        # Pierwsza superdiagonala
        if (i < (N - 1)):
            A[2][i] = (0.1 / (i + 1))

        # Druga superdiagonala
        if (i < (N - 2)):
            A[3][i] = (0.15 / ((i + 1) ** 2))

# Faktoryzacja LU zapisana w postaci 4xN
def fill_LU(N):  # O(N)
    for i in range(N):
        LU[2][i] = A[1][i] - LU[3][i-1] * LU[1][i-1]  # Diagonala U uii = aii - li(i-1) * u(i-1)i
        LU[1][i] = A[2][i] - LU[3][i-1] * LU[0][i-1 ]  # Pierwsza superdiagonala U ui(i+1) = ai(i+1) - li(i-1) * u(i-1)(i+1)
        LU[0][i] = A[3][i]  # Druga superdiagonala U ui(i+1) = ai(i+2)
        LU[3][i] = A[0][i] / LU[2][i]  # Pasmo pod diagonalą L l(i+1)i = a(i+1)i / uii

# Obliczanie wyznacznika poprzez mnożenie elementów diagonalnych U
def determinantA(N):  # O(N)
    det = 1
    for i in range (N):
        det *= LU[2][i]
    return det

# Rozwiązanie Lz = x
def solveLZ(N):
    z = []  # Wektor rozwiązań z
    for i in range(N):
        z.append(0)
    for i in range(N):
        z[i] = x[i] - LU[3][i-1] * z[i-1]
    return z

# Rozwiązanie y = A^-1 * x poprzez rozwiązanie Uy = z
def solveA(N):
    y = []
    z = solveLZ(N)
    for i in range(N):
        y.append(0)
    for i in range(N-1, -1, -1):
        if (i+2 < N):
            y[i] = (z[i] - (LU[1][i] * y[i+1] + LU[0][i] * y[i+2]))/LU[2][i]
        elif(i+1 < N):
            y[i] = (z[i] - (LU[1][i] * y[i+1])) / LU[2][i]
        else:
            y[i] = (z[i]/LU[2][i])
    return y

################################################################
#              Obliczanie czasu rozwiazania
################################################################
def zbuduj_macierz(N):
    global LU, A, x
    # INICJALIZACJA MACIERZY A i wektora rozwiązań
    kl = 1  # Liczba subdiagonali
    ku = 2  # Liczba superdiagonali
    A = [[0 for x in range(N)] for y in
         range(ku + kl + 1)]  # Zamiast calej macierzy przechowujemy tylko elementy na diagonalach rozmiar 4xN
    '''
    Macierz L będzie miała elementy niezerowe tylko na diagonali i na jednym paśmie pod diagonala. Macierz U tylko na 
    diagonali oraz na dwóch pasmach nad diagonalą. Możemy je, więc rozpisać w postaci 4xN
    U13 U24 ... U(N-2)(N) 0 0
    U12 U23 ... U(N-1)(N) 0
    U11 U12 ... UNN
    L21 L32 ... L(N)(N-1) 0
    Pamiętając, że Lii = 1
    '''
    LU = [[0 for x in range(N)] for y in range(ku + kl + 1)]
    x = []  # Wektor rozwiązań x
    for i in range(1, N + 1):
        x.append(i)

def printvector(a):
    for i in range(len(a)):
        print(a[i])

def SOLUTION(N):
    fillA(N)
    fill_LU(N)
    solveA(N)
    determinantA(N)

################################################################
#                      Testowanie
################################################################
#   Funkcja wypełniająca niezerowe elementy macierzy A

def TEST(N):
    global l, u, TA, p
    TA = np.array(np.zeros((N, N)))  # T: Cała macierz NxN w celach testowych
    l = np.zeros((N, N))
    p = np.zeros((N, N,))
    u = np.zeros((N,N))
    x = []  # Wektor rozwiązań x
    for i in range(1, N + 1):
        x.append(i)
    z = []  # Wektor rozwiązań z
    for i in range(1, N + 1):
        z.append(0)
    fill_TA(N)
    p, l, u = sc.lu(TA)
    print("------------MACIERZ A------------")
    print(TA)
    print("------------MACIERZ L------------")
    print(l)
    print("------------MACIERZ U------------")
    print(u)
    print("------------Rozwiazanie Lz = x --------")
    z = sc.solve(l, x)
    print(z)
    print("------------Rozwiazanie Uy = z --------")
    print(sc.solve(u, z))
    print("------------Rozwiazanie Ay = x --------")
    print(sc.solve(TA,x))
    print(f'"Wyznacznik macierzy A: {determinantTA(N)} "')


def fill_TA(N):
    for i in range(N):
        TA[i, i] = 1.2  # Elementy na diagonali
        if i < N - 1:
            TA[i, i + 1] = 0.1/(i+1)  # Elementy na pierwszej superdiagonali
            TA[i + 1, i] = 0.2  # Elementy na subdiagonali
        if i < N - 2:
            TA[i, i + 2] = 0.15/((i+1)**2)  # Elementy na drugiej superdiagonali

def determinantTA(N):
    d = 1
    for i in range(N):
        d *= u[i][i]
    return d

################################
#    Rozwiazanie dla N=124     #
################################
zbuduj_macierz(124)
fillA(124)
fill_LU(124)
print("Rozwiązanie za pomocą algorytmu")
printvector(solveA(124))
print(f'"Wyznacznik macierzy A: {determinantA(124)} "')
print("Rozwiązanie za pomocą bibliotek")
print(TEST(124))

################################
#        Obliczanie czasu      #
################################
for N in sizes:
    zbuduj_macierz(N)
    if (N < 50):
        czas = timeit.timeit(lambda: SOLUTION(N), number=500)
        sredni_czas = czas/500
    else:
        czas = timeit.timeit(lambda: SOLUTION(N), number=1)
        sredni_czas = czas
    times.append((N, round(sredni_czas, 5)))

wartosci_N, srednie_czasy = zip(*times)
# Rysowanie wykresu
plt.figure(figsize=(10,6))
plt.plot(wartosci_N, srednie_czasy)
plt.xlabel('Rozmiar macierzy N')
plt.ylabel('Czas trwania operacji [s]')
plt.grid(True)
plt.savefig('wykres.svg', format='svg')
plt.show()
