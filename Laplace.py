from numpy import *
from matplotlib import pyplot as plt
# from mpl_toolkits import mplot3d


class Laplace:
    def __init__(self, m, n):
        self.m = m
        self.n = n

        self.x = linspace(0, 1, m+1)
        self.y = linspace(0, 1, n+1)

        self.h = 1/m
        self.k = 1/n

        self.u = zeros((m+1, n+1))       # will hold all u values including boundary points

    def f2(self, x):
        return sin(pi * x)

    def matrix_A(self):
        g = -2 / self.h ** 2 - 2 / self.k ** 2
        a_1 = diag(g * ones(self.m-1), 0) + diag((1/self.h**2) * ones(self.m-2), 1) + diag((1/self.h**2) * ones(self.m-2), -1)
        a_2 = diag((1/self.k ** 2) * ones(self.m-1), 0)

        I = diag(ones(self.m-1), 0)
        J = diag(ones(self.m-2), 1) + diag(ones(self.m-2), -1)

        A = kron(I, a_1) + kron(J, a_2)
        return A

    def vector_b(self):
        interior_points = (self.m - 1) ** 2
        b = zeros((interior_points, 1))

        # filling in non-zero points of b
        for j in (1, self.m - 1):
            b[-j] = self.f2(self.x[j]) / self.k ** 2

        return b

    def fill_u(self):
        # fill in lower boundary temps
        for j in range(1, self.m):
            self.u[j, 0] = self.f2(self.x[j])

        u_interior = linalg.solve(self.matrix_A(), self.vector_b())     # matrix multiplication
        point = 0
        for j in range(1, self.n):
            for i in range(1, self.m):
                self.u[i, j] = u_interior[point]
                point += 1

    def display3D(self):
        X, Y = meshgrid(self.x, self.y)
        ax = plt.axes(projection='3d')
        ax.plot_surface(X, Y, self.u, cmap='cool')

        # ax.plot_surface(X, Y, self.u.transpose, cmap='warm')

        # conventions
        ax.set_title(f"m = {self.m}, n = {self.n}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

        ax.view_init()  # rotation

        plt.show()

def main():
    lap = Laplace(5, 5)
    lap.fill_u()
    lap.display3D()


if __name__ == "__main__":
    main()
