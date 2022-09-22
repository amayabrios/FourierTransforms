from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Poisson2D:
    #solve Poisson's equation in 2D: 0 = c^2(u_xx + u_yy) + f(x, y)
    def __init__(self, N1):
        self.N1 = N1
        self.N2 = N1
        self.h = 1/N1
        self.c = 1

        self.utotal = zeros((self.N1 + 1, self.N2 + 1), dtype=complex)    # row:x, col:y

        self.x = linspace(0, 1 - self.h, self.N1)
        self.y = linspace(0, 1 - self.h, self.N2)

    def particularFFT(self):
        '''fft for f(x, y) = 1/pi**2 * sin(2*pi*x) * sin(2*pi*y)'''

        f = zeros((self.N1, self.N2))  # initial size of f, not including far right and top boundaries
        for i in range(self.N1):
            for j in range(self.N2):
                f[i, j] = pi ** 2 * sin(2 * pi * self.x[i]) * sin(2 * pi * self.y[j])

        return fft.fft2(f) #fhat

    def calc_u(self):
        fhat = self.particularFFT()
        uhat = zeros((self.N1, self.N2), dtype=complex)

        for n1 in range(self.N1):
            for n2 in range(self.N2):
                W1 = 2 * pi * 1j * n1 / self.N1
                W2 = 2 * pi * 1j * n2 / self.N2
                denom = (exp(-W1) - 4 + exp(W1) + exp(-W2) + exp(W2))
                if denom != 0:
                    uhat[n1, n2] = -self.h ** 2 * fhat[n1, n2] / denom

        return fft.ifft2(uhat) #u

    def uexact(self, x, y):
        # this was the most exact answer I have yet to figure out
        return pi ** 2 / 75 * sin(2 * pi * x) * sin(2 * pi * y)

    def calc_utotal(self):
        u = self.calc_u()
        # partially fill utotal with all of u
        self.utotal[0:self.N1, 0:self.N2] = u
        # fill in last row of utotal
        self.utotal[-1, 0:self.N1] = u[1, 0:self.N1]
        # fill in last column of utotal
        self.utotal[:, -1] = self.utotal[:, 1]

    def displayFFT(self):
        self.calc_utotal()
        xtotal = linspace(0, 1, self.N1 + 1)
        ytotal = linspace(0, 1, self.N2 + 1)

        X, Y = meshgrid(xtotal, ytotal)
        ax = plt.axes(projection='3d')
        ax.plot_surface(X, Y, real(self.utotal.transpose()), rstride=1, cstride=1, cmap='cool')

        ax.set_title('u(x,y)')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.view_init(azim=None);

    def displayExact(self):
        xtotal = linspace(0, 1, self.N1 + 1)
        ytotal = linspace(0, 1, self.N2 + 1)

        X, Y = meshgrid(xtotal, ytotal)
        ax = plt.axes(projection = '3d')
        ax.plot_surface(X, Y, self.uexact(X,Y), rstride = 1, cstride = 1, cmap = 'cool')

        ax.set_title('uexact(x,y)')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.view_init(azim = None);