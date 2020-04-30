import PySimpleGUI as sg
import numpy as np
import numexpr as ne
import matplotlib.pyplot as plt
import math, cmath
import time as tm

def complex_fourier(c, N, t, T=math.pi*2):
    """Calculates value of f(t) with the help of its complex fourier series approximation.

    Arguments:
        c {function/lambda} -- function expecting n as parameter and returning n-th coefficient (c_n * exp(i*2π/T*nt))
        N {int} -- upper limit of sum, the higher N the better the approximation 
        t {np.array} -- the points where the series should be evaluated
        T {double} -- period of function (default: {2π})

    Returns:
        np.array -- returns array with fourier series approximation of f(t)
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

def real_fourier(a, b, N, t, T=math.pi*2):
    """Calculates value of f(t) with the help of its real fourier series approximation.

    Arguments:
        a {function/lambda} -- function expecting n as parameter and returning n-th coefficient (a_n * cos(n*2π/T*x))
        b {function/lambda} -- function expecting n as parameter and returning n-th coefficient (b_n * sin(n*2π/T*x))
        N {int} -- upper limit of sum, the higher N the better the approximation 
        t {np.array} -- the points where the series should be evaluated
        T {double} -- period of function (default: {2π})


    Returns:
        np.array -- returns array with fourier series approximation of f(t)
    """
    w = 2 * math.pi / T # calculate angular frequency

    sum = np.zeros(t.shape)
    a_0 = a(0) # get first coefficient
    sum = sum + a_0/2

    for n in range(1, N+1):
        # Get coefficients for cos(nwt) and sin(nwt)
        a_n = a(n)
        b_n = b(n)

        # Vectorize calculation with numpy to improve performance
        argument = n*w*t
        sum = sum + a_n * np.cos(argument) 
        sum = sum + b_n * np.sin(argument)
    
    return sum


def check_input(values):
    """Validates input and prepares for evaluation

    Arguments:
        values {list} -- list of input returned from window.read() call on PySimpleGui layout, directly modified
    Returns:
        string -- returns type of input {real, complex}
    """
    # Replace symbols with mathematical values
    for i in values:
        for pi in ['pi', 'PI', 'Pi', 'π']:
            if pi in values[i]:
                values[i] = values[i].replace(pi, str(math.pi))
        for im in ['*i', '*j']:
            if im in values[i]:
                values[i] = values[i].replace(im, '*1j')
        for im in ['j', 'i']:
            if im in values[i]:
                values[i] = values[i].replace(im, '1j')

    # Check mode
    mode = values['tab_group']

    # Check if all inputs provided and valid
    if mode == 'Real':
        for key in ['a_0', 'a_n', 'b_n', 'T_r', 'I0_r', 'I1_r']:
            if key != '':
                try:
                    if key == 'a_n' or key == 'b_n':
                        n = 1 # define n for evaluation
                    ne.evaluate(values[key])
                except:
                    return None

    if mode == 'Complex':
        for key in ['c_0', 'c_n', 'T_c', 'I0_c', 'I1_c']:
            if key != '':
                try:
                    if key == 'c_n':
                        n = 1 # define n for evaluation
                    ne.evaluate(values[key])
                except:
                    print(key)
                    return None
    return mode

def create_layout():
    """Creates layout of GUI.

    Returns:
        sg.Layout -- returns layout object
    """
    complex_layout = [
            [sg.Text('Please input data for Fourier Series visualization:')],
            [sg.Text('Coefficient c_0'), sg.InputText(key='c_0')],
            [sg.Text('Coefficients c_n'), sg.InputText(key='c_n')],
            [sg.Text('Period T'), sg.InputText(key='T_c')],
            [sg.Text('Interval'), sg.InputText(key='I0_c'), sg.InputText(key='I1_c')]
            ]
    
    real_layout = [
            [sg.Text('Please input data for Fourier Series visualization:')],
            [sg.Text('Coefficient a_0'), sg.InputText(key='a_0')],
            [sg.Text('Coefficients a_n'), sg.InputText(key='a_n')],
            [sg.Text('Coefficients b_n'), sg.InputText(key='b_n')],
            [sg.Text('Period T'), sg.InputText(key='T_r')],
            [sg.Text('Interval'), sg.InputText(key='I0_r'), sg.InputText(key='I1_r')]
            ]

    layout = [
            [sg.TabGroup([[sg.Tab('Real', real_layout, tooltip='Real Fourier Series'), sg.Tab('Complex', complex_layout, tooltip='Complex Fourier Series')]], key="tab_group")],
            [sg.Button('Ok'), sg.Button('Cancel')]
            ]    
        
    return layout


def read_input(window):
    """Reads input from window.

    Arguments:
        window {sg.Window} -- main window of GUI
    Returns:
        (string, dict) -- returns tuple with mode {'Complex', 'Real'} and inputted values
    """
    # Read input values
    while True:
        event, values = window.read() 
        if event is None or event == 'Cancel': #exit    
            return  
        
        # Check and process input   
        mode = check_input(values) 
        if mode is None:
            invalid = sg.PopupOKCancel('Please check your input and look for invalid symbols.', title='Invalid input')
            if invalid == 'Cancel':
                return
            else:
                continue
        else:
            break
    return (mode, values)

def main():
    # Create layout
    layout = create_layout()
    window = sg.Window('Fourier Series', layout) 
    mode, values = read_input(window)
    window.close()

    # Assign data to variables
    if mode == 'Real':
        a_0 = ne.evaluate(values['a_0'])
        a_n = values['a_n']
        b_n = values['b_n']
        T = ne.evaluate(values['T_r'])
        interval = (ne.evaluate(values['I0_r']), ne.evaluate(values['I1_r']))
        N = 100

        # Create lambdas for coefficient calculation
        a = lambda n: ne.evaluate(a_n) if n != 0 else a_0 
        b = lambda n: ne.evaluate(b_n) 

    
    elif mode == 'Complex':
        
        c_0 = ne.evaluate(values['c_0'])
        c_n = values['c_n']
        T = ne.evaluate(values['T_c'])
        interval = (ne.evaluate(values['I0_c']), ne.evaluate(values['I1_c']))
        N = 100

        # Create lambda for coefficient calculation
        c = lambda n: ne.evaluate(c_n) if n != 0 else c_0 

    # Create input list, calculate output and plot data
    time = np.linspace(interval[0], interval[1], 5000)
    fig, axs = plt.subplots(3, 2, num='Fourier Series')

    for i, N in enumerate([1, 2, 5, 10, 100, 1000]):
        if mode == 'Real':
            y = real_fourier(a, b, N, time, T)
        elif mode == 'Complex':
            y = complex_fourier(c, N, time, T)

        col = i % 2
        row = i//2 
        axs[row, col].plot(time, y, linewidth=1.5)
        axs[row, col].set_title(f'N = {N}')

    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()