from numpy import *
from matplotlib import pyplot as plt


class Alias:
    def __init__(self, L, N):
        self.L = L
        self.N = N
        self.dt = L/N

    def f(self):
        t = linspace(0, self.L-self.dt, self.N)
        return 6 * cos(3 * t) + 7 * cos(t)

    def shift(self):
        fhat = fft.fft(self.f())
        freqs = self.N * fft.fftfreq(self.N)

        freqs_shift = fft.fftshift(freqs)

        fhat_shift = fft.fftshift(fhat)
        omega_shift = 2 * pi * freqs_shift / (self.N * self.dt)
        return omega_shift, fhat_shift

    def display(self):
        omega_shift, fhat_shift = self.shift()[0], self.shift()[1]

        plt.plot(omega_shift, self.dt * fhat_shift, '.')

        omega_fine = linspace(omega_shift[0], omega_shift[self.N - 1], 3000)
        ghat = (self.L / 2) * 6 * (sinc(self.L * (omega_fine - 3)) + sinc(self.L * (omega_fine + 3))) + (self.L / 2) * 7 * (
                    sinc(self.L * (omega_fine - 1)) + sinc(self.L * (omega_fine + 1)))

        plt.plot(omega_fine, ghat)
        plt.show()

#alias = Alais(4*pi, 2**6)
#alias.display()