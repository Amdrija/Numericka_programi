import matplotlib.pyplot as plt
from matplotlib import style
import sympy

style.use("classic")

def f1(x):
	return sympy.exp(x) - 1

def f(x):
	return 2 * x * (1 - x**2 + x) * sympy.log(x) + 1 - x**2

def df(f):
	return sympy.diff(f)

#racuna niz vrednosti za funkciju sa gustinom alpha
#da bi je nacrtali sa plt.plot(niz_x,niz_y)
def calculate_function(function, xaxis_min, xaxis_max, alpha):
	i = 0

	x = []
	y = []

	d = (xaxis_max - xaxis_min) / alpha

	while i <= d:
		x.append(xaxis_min + i * alpha)
		y.append(function(x[i]))
		i += 1

	return x, y

def tangent(x,x0,f):
	function = sympy.lambdify(x, f)
	fprime = sympy.lambdify(x, sympy.diff(f))

	return fprime(x0) * (x - x0) + function(x0)

def secica(x,x0,x1,f):
	function = sympy.lambdify(x, f)

	return abs((function(x0)-function(x1))/(x0-x1)) * (x - x0) + function(x0)


#precision: koliko precizno da nadje resenje; alpha: koliko gusto da crta grafik
#x,yaxis: crta grafik od min do max, plot_graphs: da crta ili samo u terminalu da ispise
def newton(function, precision, alpha, m1, M2, x_starting, xaxis_min, xaxis_max, yaxis_min, yaxis_max, plot_graphs):

	xs = sympy.Symbol('x')

	x0 = []
	x0.append(x_starting)

	#za function ne treba da se uradi lambdify jer je tako zadata
	x, y = calculate_function(function, xaxis_min, xaxis_max, alpha)

	#funkcije tangenti
	xt = []
	yt = []

	i = 0

	curr = 1
	while curr > (precision*abs(m1)/abs(M2))**(1/2):

		#tangent je sa simbolom 'x', pa mora da se lambdify pre poziva calculate_function
		tangent_function = sympy.lambdify(xs, tangent(xs,x0[i],function(xs)))
		xt1, yt1 = calculate_function(tangent_function, xaxis_min, xaxis_max, 0.1)

		xt.append(xt1)
		yt.append(yt1)

		#zasto mora f_lambd kada function radi sa promenljivama????
		f_lambd = sympy.lambdify(xs, function(xs))
		df_lambd = sympy.lambdify(xs, sympy.diff(function(xs)))

		x0.append(x0[i] - f_lambd(x0[i]) / df_lambd(x0[i]))

		curr = abs(x0[i+1] - x0[i])

		i+=1

	#za poslenju tacku
	tangent_function = sympy.lambdify(xs, tangent(xs,x0[i],function(xs)))
	xt1, yt1 = calculate_function(tangent_function, xaxis_min, xaxis_max, 0.1)

	xt.append(xt1)
	yt.append(yt1)

	if plot_graphs:
		round_acc = 5
		lab = "f(x)=0 za x=" + str(round(x0[-1],round_acc))

		#grafik funckije
		plt.title("Newton's method",fontsize=20,pad=25)
		plt.plot(x,y,label=lab)
		plt.axis([xaxis_min,xaxis_max,yaxis_min,yaxis_max])

		#x i y axis
		plt.plot([xaxis_min,xaxis_max],[0,0], color='black')
		plt.plot([0,0],[yaxis_min, yaxis_max], color='black')

		for i in range(0,len(xt)):
		
			#da brojevi x0 ne idu preko grafika
			#nego da budu sa "druge strane" grafika
			pos = 0
			if function(x0[i]) >=0:
				pos = -0.1
			else:
				pos = 0.1

			#tangente, uspravne linije, tacke
			plt.plot(xt[i],yt[i],color='g')
			plt.plot([x0[i],x0[i]],[pos,function(x0[i])],'r--')
			plt.plot(x0[i],function(x0[i]),'bo')

			#text z x0[i]
			plt.text(x0[i], pos, str(round(x0[i],round_acc)),bbox=dict(facecolor='gray', alpha=0.7), horizontalalignment='center',verticalalignment='center')

		plt.legend(loc=0)

		mng = plt.get_current_fig_manager()
		#mng.window.state('zoomed')

		plt.show()
	else:
		for i in range(0,len(x0)):
			print(x0[i])

#TODO: secice, iteration, bisection methods

def main():

	newton(f, 10**(-4), 10**(-3), -7.39198, 186.0319, 0.01, 0, 1, -0.5, 1.1, True)

main()