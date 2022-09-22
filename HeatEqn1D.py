from numpy import *
from matplotlib import pyplot as plt
from time import sleep
from IPython.display import clear_output


class Heat1D:
    c = 1
    L = pi
    h = pi/16
    n = int(L/h) # there are n+1 total spatial points

    s = 1/4   #this is a stable value for s
    k = s * h**2 / c**2   #solve for what a stable time is
    Tf = 1
    m = int(Tf/k)   #there are m+1 total time points

    x = linspace(0, L, n+1)
    t = linspace(0, Tf, m+1)
    u = zeros((n+1, m+1))

    def fill_u(self):
        # fill in initial conditions into u:
        self.u[1:n, 0] = 100

        # fill in boundary conditions into u:
        self.u[0, :] = 0
        self.u[n, :] = 0

        for j in range(self.m):
          for i in range(1, self.n):
            self.u[i, j+1] = (1-2*self.s)*self.u[i, j] + self.s*(self.u[i+1, j] + self.u[i-1, j])

    def exact(self, x, t, N):
        S = zeros(len(x))   #initialization
        for k in range(N):
          n = 2*k+1
          S += (400/pi)*(1/n)*exp(-n**2 * t) * sin(n * x)
          # print(array2string(u, formatter={'float_kind': '{0:.2f}'.format}))
        return S

    def display(self):
        self.fill_u()
        for j in range(0, self.m):
          plt.plot(self.x, self.u[:, j], 'ro')
          plt.plot(self.x, self.exact(self.x, self.t[j], 10), 'b')
          plt.xlim([-.1, self.L+.1])
          plt.ylim([-1, 120])
          plt.show()
          sleep(0.01)
          clear_output(wait=True)