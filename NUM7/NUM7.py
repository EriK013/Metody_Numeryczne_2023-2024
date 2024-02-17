import numpy as np
import matplotlib.pyplot as plt

def funkcja_y(x):
    return 1 / (1 + 50 * x**2)

def funkcja_y1(x):
    return np.abs(-5*x**7)

def jednorodne_wezly(n):
    wezly = []
    for i in range (n + 1):
        wezly.append(-1 + 2*i/(n+1))
    return wezly


def czebyszew_wezly(n):
    wezly = []
    for i in range (n + 1):
        wezly.append(np.cos((2 * i + 1) * np.pi / (2 * (n+1))))
    return wezly


def interporacja_lagrange(funkcja, wezly, x):
    n = len(wezly) - 1
    def l(x,i):
        result1 = 1
        for j in range(n + 1):
            if j != i:
                result1 *= (x - wezly[j]) / (wezly[i] - wezly[j])
        return result1

    result = 0
    for i in range(n + 1):
        result += funkcja(wezly[i]) * l(x, i)
    return result


def rysuj_wykres(funkcja, wezly, siatka):
    x_values = np.linspace(-1,1,1000)
    y_values = [interporacja_lagrange(funkcja, wezly, xi) for xi in x_values]
    plt.figure(figsize=(10,6))
    plt.title(f"Interpolacja dla {len(wezly)} węzłów {siatka}")
    plt.plot(x_values, y_values, label=f"Wielomian interpolacyjny", color="orange")
    plt.plot(x_values, [funkcja(x) for x in x_values], label=f"Przebieg funkcji", color="green")
    plt.scatter(wezly, [funkcja(x) for x in wezly], label="węzły")
    plt.grid()
    plt.legend()
    #plt.savefig("iczebyszew_71.svg")
    plt.show()

nwartosci = [3,5,12,30,70]
for n in nwartosci:
    rysuj_wykres(funkcja_y, jednorodne_wezly(n), "Siatka jednorodna")
    rysuj_wykres(funkcja_y, czebyszew_wezly(n), "Siatka Czebyszewa")