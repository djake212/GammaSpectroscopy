import tkinter as tk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
from tkinter import filedialog
import sympy as sp
from scipy.optimize import curve_fit

def Gauss(x, N,b,sigma,M,C):
    return N *(1/(sigma*np.sqrt(2*np.pi))) * np.exp(-0.5*((x - b)/sigma)**2) + M*x + C

def linear(x, M, B):
    y = M*x + B
    return y

def quadratic(x, a, b, c):
    return a * x**2 + b * x + c

def r_squared(y_actual, y_predicted):
    y_mean = np.mean(y_actual)
    ss_tot = np.sum((y_actual - y_mean) ** 2)
    ss_res = np.sum((y_actual - y_predicted) ** 2)
    return 1 - (ss_res / ss_tot)

def error_propagation(f, var_values, var_errors):
    """
    f: function with sympy symbols as input variables
    var_values: dictionary with keys as the sympy symbols and values as the variable values
    var_errors: dictionary with keys as the sympy symbols and values as the variable errors
    """

    # Calculate the partial derivatives of the function for each variable
    partial_derivatives = [sp.diff(f, var) for var in var_values.keys()]

    # Calculate the propagated error using partial derivatives
    sigma_f_squared = sum((partial_derivative * var_errors[var])**2 for partial_derivative, var in zip(partial_derivatives, var_values.keys()))

    # Substitute the variable values
    sigma_f_val = sp.sqrt(sigma_f_squared.subs(var_values))

    return sigma_f_val

def on_select(eclick, erelease):
    global channel, counts, parameters2
    
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    if x1 > x2:
        x1, x2 = x2, x1

    start_idx = np.searchsorted(channel, x1)
    end_idx = np.searchsorted(channel, x2)
    
    x_fit = channel[start_idx:end_idx]
    y_fit = counts[start_idx:end_idx]

    N_init = y_fit.max()
    b_init = x_fit[np.argmax(y_fit)]
    sigma_init = (x2 - x1) / 4
    yerror = np.sqrt(y_fit)

    try:
        popt, pcov = curve_fit(Gauss, x_fit, y_fit, p0=[N_init, b_init, sigma_init,0,0], sigma = yerror, absolute_sigma=True)
        gaussian_y_fit = Gauss(x_fit, *popt)
        errormatrix = np.sqrt(np.diag(pcov))

        FWHM = 2.35482 * popt[2]
        Rsquare = r_squared(gaussian_y_fit,y_fit)
        x = sp.symbols('x')
        KeV = parameters2[0] * x**2 + parameters2[1] * x + parameters2[2]
        var_values = {x:popt[1]}
        var_errors = {x:errormatrix[1]}

      # Calculate propagated error for findcurrentactivity function
        error_KeV = error_propagation(KeV, var_values, var_errors)

        graphlabel = f"peakchannel = {popt[1]:.0f}(error:+-{errormatrix[1]:.3f}), KeV = {KeV.subs(x, popt[1]):.2f}(error:+-{error_KeV:.3f}) \n| integral value = {abs(popt[0]):.2f}(error:+-{errormatrix[0]:.3f}) | standard dev = {abs(popt[2]):.2f}(error:+-{errormatrix[2]:.3f}) R^2 = {Rsquare:.2f} | FWHM = {abs(FWHM):.3f}(error: +- {(errormatrix[2]*2.35482):.3f})"
        print(graphlabel)
        plt.plot(x_fit, gaussian_y_fit, label=graphlabel)
        plt.legend()
        plt.draw()
    except:
        pass

#grab graph title from the filename
file_path = filedialog.askopenfilename()
graphtitle = str(file_path.split('/')[-1]).rstrip('.csv')

#create grab the csv data and ignore and comments and the heading info
data = np.loadtxt(file_path,comments='#',skiprows=25,delimiter=',')
data = data.T
channel = data[0]
energy = data[1]
counts = data[2]

#plotting the channel vs counts data
fig,ax = plt.subplots()
ax.scatter(channel,counts,s=10)
ax.set_xlabel('channel number')
ax.set_ylabel('Counts')
ax.set_yscale('log')
ax.set_title(graphtitle)
fig.set_size_inches(16, 8)

# fig2,ax2 = plt.subplots()
# ax2.scatter(channel,energy,s=5)
# ax2.set_ylabel('Gamma Energy KeV')
# ax2.set_xlabel('Channel Number')
# ax2.set_title(f"{graphtitle} energy vs channel")
# ax2.legend()

parameters1, pcov1 = curve_fit(linear, channel, energy)
parameters2, pcov2 = curve_fit(quadratic, channel, energy)
linear_fit = linear(channel, *parameters1)
quadratic_fit = quadratic(channel, *parameters2)
r_squared_linear = r_squared(energy, linear_fit)
r_squared_quadratic = r_squared(energy, quadratic_fit)
linear_label = f"Linear fit: y = {parameters1[0]:.2f}x + {parameters1[1]:.2f}, R^2 = {r_squared_linear:.6f}"
quadratic_label = f"Quadratic fit: y = {parameters2[0]:.2e}x^2 + {parameters2[1]:.2f}x + {parameters2[2]:.2f}, R^2 = {r_squared_quadratic:.6f}"

# ax2.plot(channel, linear(channel, *parameters1).T, label=linear_label, color='red')  
# ax2.plot(channel, quadratic(channel, *parameters2).T, label=quadratic_label, color='green')
# ax2.legend()

# Register mouse event for selecting range
rect_select = RectangleSelector(ax, on_select, useblit=True)

plt.show()
