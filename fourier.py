import PySimpleGUI as sg
import numpy as np
import numexpr as ne
import matplotlib.pyplot as plt
import math, cmath
import time as tm

#TODO Implement real fourier series

def complex_fourier(c, N, t, T=math.pi*2):
    """Calculates value of f(t) with the help of its complex fourier series approximation.

    Arguments:
        c {function/lambda} -- function expecting n as parameter and returning n-th coefficient (c_n * exp(i*2π/T*nt))
        N {int} -- upper limit of sum, the higher N the better the approximation 
        t {np.array} -- the points where the series should be evaluated
        T {double} -- period of function (default: {2π})

    Returns:
        double -- returns fourier series approximation of f(t)
    """ 
    w = 2 * math.pi / T # calculate angular frequency

    sum = np.zeros(t.shape,dtype=np.complex128)
    c_0 = c(0) # get first coeffiecient
    sum = sum + c_0

    for n in range(1, N+1):
        # Get coefficients for e^inwt and e^-inwt
        c_n = c(n)
        c_nn = c(-n)

        # Vectorize calculation with numpy to improve performance
        exponent = 1j * n * w * t
        sum = sum + c_n * np.exp(exponent)
        sum = sum + c_nn * np.exp(-exponent)

    return np.real(sum)

def real_fourier(a, b, N, t):
    """Calculates value of f(t) with the help of its real fourier series approximation.

    Arguments:
        a {function/lambda} -- function expecting n as parameter and returning n-th coefficient (a_n * cos(n*2π/T*x))
        b {function/lambda} -- function expecting n as parameter and returning n-th coefficient (b_n * sin(n*2π/T*x))
        N {int} -- upper limit of sum, the higher N the better the approximation 
        t {double} -- the point where the series should be evaluated
        T {double} -- period of function (default: {2π})


    Returns:
        double -- returns fourier series approximation of f(t)
    """
    w = 2 * math.pi / T # calculate angular frequency

    a_0 = a(0) # get first coefficient
    sum = a_0/2

    for n in range(1, N+1):
        a_n = a(n)
        b_n = b(n)

        sum += a_n * math.cos(n*w*t) 
        sum += b_n * math.sin(n*w*t)
    
    return sum


def check_input(values):
    """Validates input and prepares for evaluation

    Arguments:
        values {list} -- list of input returned from window.read() call on PySimpleGui layout, directly modified
    """
    for i in values:
        for pi in ["pi", "PI", "Pi", "π"]:
            if pi in values[i]:
                values[i] = values[i].replace(pi, str(math.pi))
        for im in ["*i", "*j"]:
            if im in values[i]:
                values[i] = values[i].replace(im, "*1j")
        for im in ["j", "i"]:
            if im in values[i]:
                values[i] = values[i].replace(im, "1j")


def main():
    # Create layout
    layout = [[sg.Text('Please input data for Fourier Series visualization:')],
            [sg.Text('Coefficient c_0'), sg.InputText()],
            [sg.Text('Coefficients c_n'), sg.InputText()],
            [sg.Text('Period T'), sg.InputText()],
            [sg.Text('Interval'), sg.InputText(), sg.InputText()],
            [sg.Button('Ok'), sg.Button('Cancel')]]

    window = sg.Window('Fourier Series', layout)    
    event, values = window.read()    
    check_input(values)
    window.close()

    # Assign data to variables
    c_0 = ne.evaluate(values[0])
    c_n = values[1]
    T = ne.evaluate(values[2])
    interval = (ne.evaluate(values[3]), ne.evaluate(values[4]))
    N = 100

    # Create lambda for coefficient calculation
    c = lambda n: ne.evaluate(c_n) if n != 0 else c_0 

    # Create input list, calculate output and plot data
    time = np.linspace(interval[0], interval[1], 5000)
    fig, axs = plt.subplots(3, 2, num="Fourier Series")

    for i, N in enumerate([1, 2, 5, 10, 100, 1000]):
        y = complex_fourier(c, N, time)

        col = i % 2
        row = i//2 
        axs[row, col].plot(time, y, linewidth=1.5)
        axs[row, col].set_title(f'N = {N}')

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()