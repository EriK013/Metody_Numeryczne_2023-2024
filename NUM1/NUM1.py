import numpy as np
import matplotlib.pyplot as plt

x = 0.2
h_skala = np.logspace(-16,0,100) #skala logarytmiczna h od 10^-16 do 10^0
a_blad_float = []
b_blad_float = []
a_blad_double = []
b_blad_double = []


def f(x): #Funkcja f(x)
    return np.sin(x**2)

def f2(x): #Druga funkcja
    return np.cos(x**2)

def df(x): #Pochodna funkcji f(x)
    return 2*x*np.cos(x**2)

def df2(x): #Pochodna drugiej funkcji f(x)
    return -2*x*np.sin(x**2)

#Wzory na przyblizenie funkcji
def a_Dhf(x,h,f):
    return (f(x+h) - f(x)) / h

def b_Dhf(x,h,f):
    return (f(x+h) - f(x-h)) / (2*h)

#Funkcja na obliczenie błędu dla pojedyńczego h
def blad(x, h, f, Df, df):
    przyblizone_f = Df(x, h, f)
    dokladne_f = df(x)
    blad = np.abs(przyblizone_f - dokladne_f) #|Dhf(x) - f'(x)|
    return blad

#Funkcja obliczająca błędy dla h w skali logarytmicznej
def oblicz_bledy(x, skala, f, df):
    a_blad_float.clear()
    a_blad_double.clear()
    b_blad_float.clear()
    b_blad_double.clear()
    for h in skala:
        a_blad_float.append(blad(np.float32(x), np.float32(h), f, a_Dhf, df))
        a_blad_double.append(blad(np.float64(x), np.float64(h), f, a_Dhf, df))
        b_blad_float.append(blad(np.float32(x), np.float32(h), f, b_Dhf, df))
        b_blad_double.append(blad(np.float64(x), np.float64(h), f, b_Dhf, df))


#Ustawienia rysowania wykresu
def rysuj(title):
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('h')
    plt.ylabel('|Dhf(x) - f\'(x)|')
    plt.title(title)
    plt.legend()
    plt.show()

def wykresy_bledu(funkcja, x):
    #Dla float i pierwszego wzoru
    plt.figure(figsize=(10, 6))
    plt.plot(h_skala, a_blad_float, label="(f(x+h) - f(x)) / h")
    rysuj(f'Bład w zależności od h dla float {funkcja} {x}')

    #Dla float i drugiego wzoru
    plt.figure(figsize=(10, 6))
    plt.plot(h_skala, b_blad_float, label="(f(x+h) - f(x-h)) / (2*h)")
    rysuj(f'Bład w zależności od h dla float {funkcja} {x}')

    #Dla float
    plt.figure(figsize=(10, 6))
    plt.plot(h_skala, a_blad_float, label="(f(x+h) - f(x)) / h")
    plt.plot(h_skala, b_blad_float, label="(f(x+h) - f(x-h)) / (2*h)")
    rysuj(f'Bład w zależności od h dla float {funkcja} {x}')


    #Dla double i pierwszego wzoru
    plt.figure(figsize=(10, 6))
    plt.plot(h_skala, a_blad_double, label="(f(x+h) - f(x)) / h")
    rysuj(f'Bład w zależności od h dla double {funkcja} {x}')

    #Dla double i drugiego wzoru
    plt.figure(figsize=(10, 6))
    plt.plot(h_skala, b_blad_double, label="(f(x+h) - f(x-h)) / (2*h)")
    rysuj(f'Bład w zależności od h dla double {funkcja} {x}')

    #Dla double
    plt.figure(figsize=(10, 6))
    plt.plot(h_skala, a_blad_double, label="(f(x+h) - f(x)) / h")
    plt.plot(h_skala, b_blad_double, label="(f(x+h) - f(x-h)) / (2*h)")
    rysuj(f'Bład w zależności od h dla double {funkcja} {x}')

    #Dla wszystkich przypadków
    plt.figure(figsize=(10, 6))
    plt.plot(h_skala, a_blad_float, label="Dhf1 float")
    plt.plot(h_skala, b_blad_float, label="Dhf2 float")
    plt.plot(h_skala, a_blad_double, label="Dhf1 double")
    plt.plot(h_skala, b_blad_double, label="Dhf2 double")
    rysuj(f'Bład w zależności od h dla float i double {funkcja} {x}')

oblicz_bledy(x, h_skala, f, df)
wykresy_bledu("sin(x^2)", x)
x = 100.0
oblicz_bledy(x, h_skala, f2, df2)
wykresy_bledu("cos(x^2)", x)
