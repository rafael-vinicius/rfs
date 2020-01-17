import numpy as np
import sympy as sp
from sympy import Eq, solve_linear_system, Matrix, simplify
import matplotlib.pyplot as plt


np.set_printoptions(formatter={'float': '{: 0.6f}'.format})
print("number of nodes: ")
nodes= 30 ## quant nodes
real_nodes = int(nodes)+1
y = np.array(sp.symbols("y0:"+ str(real_nodes)))
VA = sp.symbols("VA")
print(y)


y_1 = np.array(y)
y_1 = np.insert(y_1,-1,VA)
print(y_1)

print("Max value of x: ")
x_max = 4
print("Min value of x: ")
x_min = 0

t = float((float(x_max) - float(x_min))/(float(nodes) - 1))
print("Step size: ")
print(t)

x = np.arange(0,float(x_max), t)
x = np.append(x, float(x_max))
print("X values created: ")
print(x)


print("\nCondition 1")
x1_ab = 1
y_ab1 = 0

print("\nCondition 2")
x2_ab = nodes ## quant nodes
y_ab2 = 0

print("\nCondition 3")
x3_ab = 0
dy_ab = 0
print("\n")

M1EI = 0.00008329

print(y[4])
print(y[3])
print(y[1])
print(x[2])


equations = np.array([])
print("Equations created: ")

for i in range(1, int(real_nodes)):
  if y[i] == sp.symbols("y1") or y[i] == sp.symbols("y30"): ## y6 é sempre o valor da quant nodes
    print("pass")
  else:
    eq = Eq(((y[i+1] - 2*y[i] + y[i-1])/t**2) - VA*(x[i-1]) + 4*VA,((- 15*(x[i-1]**2) + 260)*M1EI))
    eq = eq.subs({y[int(x1_ab)]:int(y_ab1)})
    eq = eq.subs({y[int(x2_ab)]:int(y_ab2)})
    equations = np.append(equations, eq)
    print(equations[i-2])


eq_Deriv_CC2 = Eq((((y[3] - y[1])/(2*t)) - 0.5*VA*(x[2]**2) + 4*VA*(x[2])),((-5*(x[2]**3) +240*(x[2]))*M1EI))
eq_Deriv_CC2 = eq_Deriv_CC2.subs({y[int(x3_ab)]:int(dy_ab)})
eq_Deriv_CC2 = eq_Deriv_CC2.subs({y[int(x1_ab)]:int(y_ab1)})
eq_Deriv_CC2 = eq_Deriv_CC2.subs({y[int(x2_ab)]:int(y_ab2)})
equations = np.append(equations, eq_Deriv_CC2)
print(equations[i-2])

system, equality = sp.linear_eq_to_matrix(equations, y_1[2:31]) ## ALTERAR aqui também, sempre nodes+1 ## Sempre a quantidade desejada de variaveis + 1 !!!
system_array = np.array(system).astype(np.float64)
##system_array2 = np.delete(system_array, 2, 1)
system_equality = np.array(equality).astype(np.float64)

print("\nCoefficient Matrix")
print(system_array)
print("\nEquality Matrix")
print(system_equality)


ans = np.linalg.solve(system_array, system_equality)
ans2 = np.insert(ans, 0, y_ab1)
ans3 = np.insert(ans2,(int(nodes)-1), y_ab2) ## O valor de VA está sendo subst. pela cond contorno.
print("\nVA valor Teste")
print(ans)##teste
ans4 = np.resize(ans3, (int(nodes), 1))

print("\nAnswer column Matrix")
print(ans4)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.axhline(linewidth=2, color='black', alpha=0.4)
plt.axvline(linewidth=2, color='black', alpha=0.4)
plt.grid(color="gray", linewidth=0.5)
plt.plot(x, ans4)
plt.show()



