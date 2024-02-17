import numpy as np
import matplotlib.pyplot as plt
def f(x):
    return np.sin(x) - 0.4

def g(x):
    return (f(x))**2

def df(x): # Analitycznie obliczona pochodna f(x)
    return np.cos(x)

def dg(x): # Analitycznie obliczona pochodna g(X)
    return 2 * np.cos(x) * (np.sin(x) - 0.4)

def u(x):
    return g(x) / dg(x)

def du(x): # Analitycznie obliczona pochodna u(x)
    return (((np.sin(x)-0.4)*np.sin(x))/(2*np.cos(x)**2)) + 0.5
def bisekcja(funkcja, przedzial, tol=1e-15):
    if (funkcja(przedzial[0]) * funkcja(przedzial[1]) < 0):  # Sprawdzenie czy możemy zastosować metodę dla funkcji w przedziale
        a = przedzial[0]
        b = przedzial[1]
        roznicaRozwiazania = []
        kroki = 0
        c = a
        while(np.abs(c - np.arcsin(0.4))) > tol:  #
            c = (a+b)/2  # Obliczamy punkt środkowy
            if (funkcja(a) * funkcja(c) < 0): # Sprawdzamy, którą część wykresu mamy pobrać
                b = c
            elif (funkcja(c) * funkcja(b) < 0):
                a = c
            kroki += 1
            roznicaRozwiazania.append(np.abs(c - np.arcsin(0.4)))
        return c, roznicaRozwiazania, kroki
    print("Nie można zastosować metody dla tej funkcji")

def falsi(funkcja, przedzial, tol=1e-15):
    if (funkcja(przedzial[0]) * funkcja(
            przedzial[1]) < 0):
        a = przedzial[0]
        b = przedzial[1]
        roznicaRozwiazania = []
        kroki = 0
        c = a
        while (np.abs(c - np.arcsin(0.4))) > tol:
            c = (-b*funkcja(a) + a*funkcja(b))/(funkcja(b) - funkcja(a))
            if (funkcja(a) * funkcja(c) < 0):
                b = c
            elif (funkcja(c) * funkcja(b) < 0):
                a = c
            kroki += 1
            roznicaRozwiazania.append(np.abs(c - np.arcsin(0.4)))
        return c, roznicaRozwiazania, kroki
    print("Nie można zastosować metody dla tej funkcji")

def siecznych(funkcja, przedzial, tol=1e-15, maxIter=100):
    a = przedzial[0]
    b = przedzial[1]
    roznicaRozwiazan = []
    c = 0
    kroki = 0
    if (funkcja(a) != funkcja(b)):
        while((np.abs(c - np.arcsin(0.4))) > tol and kroki < maxIter):
            c = (a * funkcja(b) - b*funkcja(a))/(funkcja(b) - funkcja(a))
            b = a
            a = c
            roznicaRozwiazan.append(np.abs(c - np.arcsin(0.4)))
            kroki += 1
        return c, roznicaRozwiazan, kroki
    print("Niespełniony warunek f(x1) != f(x2)")

def newton(f, df, przedzial, tol=1e-15, maxIter=100):
    kroki = 0
    x0 = przedzial[0]
    roznicaRozwiazan = []
    while(np.abs(x0 - np.arcsin(0.4))) > tol and kroki < maxIter:
        x0 = x0 - f(x0) / df(x0)
        roznicaRozwiazan.append(np.abs(x0 - np.arcsin(0.4)))
        kroki += 1
    return x0, roznicaRozwiazan, kroki

def main():
    przedzial = [0, np.pi / 2]

    #****PODPUNKT A*****
    print("#### sin(x) - 0.4 ####")
    bx0, bRRozwiazan, bKroki = bisekcja(f, przedzial)
    print(f"Miejsce zerowe uzyskane metodą bisekcji={bx0} prawdziwe miejsce zerowe={np.arcsin(0.4)}. Uzyskano po {bKroki} krokach")

    fx0, fRRozwiazan, fKroki = falsi(f, przedzial)
    print(f"Miejsce zerowe uzyskane metodą falsi={fx0} prawdziwe miejsce zerowe={np.arcsin(0.4)}. Uzyskano po {fKroki} krokach")

    sx0, sRRozwiazan, sKroki = siecznych(f,przedzial)
    print(f"Miejsce zerowe uzyskane metodą siecznych={sx0} prawdziwe miejsce zerowe={np.arcsin(0.4)}. Uzyskano po {sKroki} krokach")

    nx0, nRRozwiazan, nKroki = newton(f, df, przedzial)
    print(f"Miejsce zerowe uzyskane metodą newtona={nx0} prawdziwe miejsce zerowe={np.arcsin(0.4)}. Uzyskano po {nKroki} krokach")

    #****PODPUNKT B*****
    print("#### (sin(x) - 0.4)^2 ####")
    bisekcja(g, przedzial)
    falsi(g, przedzial)

    sgx0, sgRRozwiazan, sgKroki = siecznych(g, przedzial)
    print(f"Miejsce zerowe uzyskane metodą siecznych={sgx0} prawdziwe miejsce zerowe={np.arcsin(0.4)}. Uzyskano po {sgKroki} krokach")

    ngx0, ngRRozwiazan, ngKroki = newton(g, dg, przedzial)
    print(f"Miejsce zerowe uzyskane metodą newtona={ngx0} prawdziwe miejsce zerowe={np.arcsin(0.4)}. Uzyskano po {ngKroki} krokach")

    #****USPRAWNIENIE*****
    print("#### Funkcja usprawniona ####")
    sux0, suRRozwiazan, suKroki = siecznych(u, przedzial)
    print(
        f"Miejsce zerowe uzyskane metodą siecznych={sux0} prawdziwe miejsce zerowe={np.arcsin(0.4)}. Uzyskano po {suKroki} krokach")

    nux0, nuRRozwiazan, nuKroki = newton(u, du, przedzial)
    print(
        f"Miejsce zerowe uzyskane metodą newtona={nux0} prawdziwe miejsce zerowe={np.arcsin(0.4)}. Uzyskano po {nuKroki} krokach")

    # Rysowanie wykresu f
    plt.figure(figsize=(10, 6))
    plt.title("Zachowane xi - x* dla funkcji f")
    plt.plot(bRRozwiazan,  label="Bisekcja")
    plt.plot(fRRozwiazan, label="Falsi")
    plt.plot(sRRozwiazan, label="Siecznych")
    plt.plot(nRRozwiazan, label="Newtona")
    plt.yscale('log')
    plt.xlabel("Liczba kroków")
    plt.ylabel("Wielkość przedziału")
    plt.grid()
    plt.legend()
    #plt.savefig("Zbierznosci_F.svg")
    plt.show()

    # Rysowanie wykresu g
    plt.figure(figsize=(10, 6))
    plt.title("Zachowane xi - x* dla funkcji g")
    plt.plot(sgRRozwiazan, label="Siecznych")
    plt.plot(ngRRozwiazan, label="Newtona")
    plt.yscale('log')
    plt.xlabel("Liczba kroków")
    plt.ylabel("Wielkość przedziału")
    plt.grid()
    plt.legend()
    #plt.savefig("Zbierznosci_G.svg")
    plt.show()

    # Rysowanie wykresu u
    plt.figure(figsize=(10, 6))
    plt.title("Zachowane xi - x* dla funkcji u")
    plt.plot(suRRozwiazan, label="Siecznych")
    plt.plot(nuRRozwiazan, label="Newtona")
    plt.yscale('log')
    plt.xlabel("Liczba kroków")
    plt.ylabel("Wielkość przedziału")
    plt.grid()
    plt.legend()
    #plt.savefig("Zbierznosci_U.svg")
    plt.show()

main()