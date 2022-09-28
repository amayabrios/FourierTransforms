from numpy import *
from matplotlib import pyplot as plt

# from mpl_toolkits import mplot3d
# from time import sleep
# from IPython.display import clear_output


class Heat1D:
    def __init__(self, c=1, L=pi, h=pi/16, s=1/4, Tf=1):
        # instance variables are changeable when initialized but most stable at default values
        self.c = c
        self.L = L
        self.h = h
        self.n = int(L/h)    # there are n+1 total SPATIAL points

        self.s = s      # s=1/4 is a stable value for s under these conditions
        self.k = s * h**2 / c**2   # solve for what a stable time is
        self.Tf = Tf
        self.m = int(Tf/self.k)   # there are m+1 total TIME points

        self.x = linspace(0, self.L, self.n+1)
        self.t = linspace(0, self.Tf, self.m+1)
        self.u = zeros((self.n+1, self.m+1))

    def fill_u(self):
        # fill in initial conditions into u:
        self.u[1:self.n, 0] = 100

        # fill in boundary conditions into u:
        self.u[0, :] = 0
        self.u[self.n, :] = 0

        for j in range(self.m):
            for i in range(1, self.n):
                self.u[i, j+1] = (1-2*self.s)*self.u[i, j] + self.s*(self.u[i+1, j] + self.u[i-1, j])

    def exact(self, x, t, N):
        S = zeros(len(x))       # initialization
        for k in range(N):
            n = 2*k+1
            S += (400/pi)*(1/n)*exp(-n**2 * t) * sin(n * x)
            # print(array2string(u, formatter={'float_kind': '{0:.2f}'.format}))
        return S

    def display(self):
        self.fill_u()

        plt.ion()   # default on for PyCharm
        for j in range(0, self.m):
            plt.plot(self.x, self.u[:, j], 'ro')
            plt.plot(self.x, self.exact(self.x, self.t[j], 10), 'b')
            plt.xlim([-.1, self.L+.1])
            plt.ylim([-1, 120])

            plt.draw()
            plt.pause(0.01)

        # keep window open when loop done
        plt.ioff()
        plt.show()

    def display3D(self):
        self.fill_u()

        T, X = meshgrid(self.t, self.x)
        ax = plt.axes(projection='3d')
        ax.plot_surface(T, X, self.u, rstride=1, cstride=1, cmap='cool', edgecolor='none')
        plt.show()

        # conventions
        ax.set_title('u(x, t)')
        ax.set_xlabel('time')
        ax.set_ylabel('x')


def main():
    heat = Heat1D()
    heat.display()
    heat.display3D()


if __name__ == "__main__":
    main()
