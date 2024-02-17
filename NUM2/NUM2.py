import numpy as np
import scipy.linalg as sc

#Macierze
A1 = np.array([
    [2.554219275, 0.871733993, 0.052575899, 0.240740262, 0.316022841],
    [0.871733993, 0.553460938, -0.070921727, 0.255463951, 0.707334556],
    [0.052575899, -0.070921727, 3.409888776, 0.293510439, 0.847758171],
    [0.240740262, 0.255463951, 0.293510439, 1.108336850, -0.206925123],
    [0.316022841, 0.707334556, 0.847758171, -0.206925123, 2.374094162]
])

A2 = np.array([
    [2.645152285, 0.544589368, 0.009976745, 0.327869824, 0.424193304],
    [0.544589368, 1.730410927, 0.082334875, -0.057997220, 0.318175706],
    [0.009976745, 0.082334875, 3.429845092, 0.252693077, 0.797083832],
    [0.327869824, -0.057997220, 0.252693077, 1.191822050, -0.103279098],
    [0.424193304, 0.318175706, 0.797083832, -0.103279098, 2.502769647]
])

#Wektor wyrazów wolnych
b = np.array([[-0.642912346, -1.408195475, 4.595622394, -5.073473196, 2.178020609]]).T
y1 = sc.solve(A1, b) #A1y = b
y2 = sc.solve(A2, b) #A2y = b

#Zaburzenie delta b
delta_b = np.random.normal(0, 1e-6, len(b)).reshape(-1, 1)
#Aiy = b + delta_b
y1_delta = sc.solve(A1, (b + delta_b))
y2_delta = sc.solve(A2, (b + delta_b))

print("*******")
print("ZABURZENIE WEKTORA b")
print(delta_b)
print("*******")
print("A1y = b")
print("y1:")
print(y1)
print("\n")
print("A2y = b")
print("y2:")
print(y2)
print("\n")
print("A1y = b + delta_b")
print("y1_delta:")
print(y1_delta)
print("\n")
print("A2y = b + delta_b")
print("y2_delta:")
print(y2_delta)

k1 = np.linalg.cond(A1)
k2 = np.linalg.cond(A2)
print("K dla A1: ")
print(k1)
print("K dla A2: ")
print(k2)