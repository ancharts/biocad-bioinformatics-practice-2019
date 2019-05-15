from math import sqrt
import numpy as np


x_crd_l = 30
x_crd_r = 38
y_crd_l = 38
y_crd_r = 46
z_crd_l = 46
z_crd_r = 54

class Point:
    def __init__(self, x, y=None, z=None):
        if type(x) is str:
            self.x = float(x[x_crd_l:x_crd_r])
            self.y = float(x[y_crd_l:y_crd_r])
            self.z = float(x[z_crd_l:z_crd_r])
            return
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return '(%f, %f, %f)' % (self.x, self.y, self.z)

    def __add__(self, vector):
        return Point(self.x + vector.x, self.y + vector.y, self.z + vector.z)

    def __sub__(self, vector):
        return Point(self.x - vector.x, self.y - vector.y, self.z - vector.z)

    def between(self, a, b):
        return Vector(a, self)*Vector(a, b) >= 0 and Vector(b, self)*Vector(b, a) >= 0

class Vector:

    def __init__(self, a, b=None):
        if b is None:
            self.x = a.x
            self.y = a.y
            self.z = a.z
        else:
            self.x = b.x - a.x
            self.y = b.y - a.y
            self.z = b.z - a.z

    def __str__(self):
        return '(%f, %f, %f)' % (self.x, self.y, self.z)

    def __add__(self, other):
        return Vector(Point(self.x + other.x, self.y + other.y, self.z + other.z))

    def __sub__(self, other):
        return Vector(Point(self.x - other.x, self.y - other.y, self.z - other.z))

    def __mul__(self, other):
        if type(other) is Vector:
            return self.x * other.x + self.y * other.y + self.z * other.z
        return Vector(Point(self.x * other, self.y * other, self.z * other))

    def __truediv__(self, other):
        return Vector(Point(self.x / other, self.y / other, self.z / other))

    def len(self):
        return sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def to_point(self):
        return Point(self.x, self.y, self.z)

    def normalize(self):
        len = self.len()
        return Vector(Point(self.x / len, self.y / len, self.z / len))

class Plane:

    def __init__(self, vector, point, p3=None):
        if p3 == None:
            # перпендикулярная вектору, проходящая через данную точку
            self.a, self.b, self.c = vector.x, vector.y, vector.z
            self.d = vector.x * point.x + vector.y * point.y + vector.z * point.z
        else:
            A = np.array([[1, vector.y, vector.z],
                          [1, point.y, point.z],
                          [1, p3.y, p3.z]])
            B = np.array([[vector.x, 1, vector.z],
                          [point.x, 1, point.z],
                          [p3.x, 1, p3.z]])
            C = np.array([[vector.x, vector.y, 1],
                          [point.x, point.y, 1],
                          [p3.x, p3.y, 1]])
            D = np.array([[vector.x, vector.y, vector.z],
                          [point.x, point.y, point.z],
                          [p3.x, p3.y, p3.z]])
            self.a = np.linalg.det(A)
            self.b = np.linalg.det(B)
            self.c = np.linalg.det(C)
            self.d = np.linalg.det(D)


    def __str__(self):
        return '%fx + %fy + %fz = %f' % (self.a, self.b, self.c, self.d)

    def __eq__(self, other):
        if (self.a == other.a and self.b == other.b and self.c == other.c and self.d == other.d):
            return True
        if (self.a == -other.a and self.b == -other.b and self.c == -other.c and self.d == -other.d):
            return True
        return False

    def is_valid(self):
        return not (self.a == 0 and self.b == 0 and self.c == 0)


class Line:
    # p1, p2 = Plane(Vector(Point(0, 0, 0)), Point(0, 0, 0)), \
    #          Plane(Vector(Point(0, 0, 0)), Point(0, 0, 0))
    #если что-нибудь не работает, то из-за этого
    def __init__(self, a, b):
        self.p1 = Plane(Vector(Point(0, 0, 0)), Point(0, 0, 0))
        self.p2 = Plane(Vector(Point(0, 0, 0)), Point(0, 0, 0))
        v = Vector(a, b)
        self.p1.a, self.p2.a = v.y, 0
        self.p1.b, self.p2.b = -v.x, v.z
        self.p1.c, self.p2.c = 0, -v.y
        self.p1.d, self.p2.d = v.y * a.x - v.x * a.y, v.z * a.y - v.y * a.z

        # случай с нулями
        if (self.p1 == self.p2 or not self.p1.is_valid()):
            self.p1.a = v.z
            self.p1.b = 0
            self.p1.c = -v.x
            self.p1.d = v.z * a.x - v.x * a.z
        if (not self.p2.is_valid()):
            self.p2.a = v.z
            self.p2.b = 0
            self.p2.c = -v.x
            self.p2.d = v.z * a.x - v.x * a.z

    def __str__(self):
        return str(self.p1) + '\n' + str(self.p2)


def vmul(a, b):
    # векторное умножение
    matrix1 = np.array([[a.y, a.z], [b.y, b.z]])
    matrix2 = np.array([[a.x, a.z], [b.x, b.z]])
    matrix3 = np.array([[a.x, a.y], [b.x, b.y]])
    return Vector(Point(np.linalg.det(matrix1), -np.linalg.det(matrix2), np.linalg.det(matrix3)))


def find_intersection_point(a, b, point):
    line = Line(a, b)
    plane = Plane(Vector(a, b), point)
    matrix0 = np.array([[line.p1.a, line.p1.b, line.p1.c],
                        [line.p2.a, line.p2.b, line.p2.c],
                        [plane.a, plane.b, plane.c]])
    matrix1 = np.array([[line.p1.d, line.p1.b, line.p1.c],
                        [line.p2.d, line.p2.b, line.p2.c],
                        [plane.d, plane.b, plane.c]])
    matrix2 = np.array([[line.p1.a, line.p1.d, line.p1.c],
                        [line.p2.a, line.p2.d, line.p2.c],
                        [plane.a, plane.d, plane.c]])
    matrix3 = np.array([[line.p1.a, line.p1.b, line.p1.d],
                        [line.p2.a, line.p2.b, line.p2.d],
                        [plane.a, plane.b, plane.d]])
    det0 = np.linalg.det(matrix0)
    det1 = np.linalg.det(matrix1)
    det2 = np.linalg.det(matrix2)
    det3 = np.linalg.det(matrix3)
    if det0 == 0:
        raise Exception('You tried to divide by zero.')
    return Point(det1 / det0, det2 / det0, det3 / det0)

def distance_from_point_to_line(a1, b1, b2):
    line = Line(b1, b2)
    perp = Line(Point(0, 0, 0), Point(0, 0, 0))
    perp.p1 = Plane(Vector(b1, b2), a1)
    perp.p2 = Plane(a1, b1, b2)
    point = lines_intersection(line, perp)
    return Vector(a1, point).len()


def point_on_line(x, a1, a2):
    line = Line(a1, a2)
    matrix0 = np.array([[line.p1.b, line.p1.c],
                        [line.p2.b, line.p2.c]])
    matrix1 = np.array([[-line.p1.a * x + line.p1.d, line.p1.c],
                        [-line.p2.a * x + line.p2.d, line.p2.c]])
    matrix2 = np.array([[line.p1.b, -line.p1.a * x + line.p1.d],
                        [line.p2.b, -line.p2.a * x + line.p2.d]])
    det0 = np.linalg.det(matrix0)
    det1 = np.linalg.det(matrix1)
    det2 = np.linalg.det(matrix2)
    return Point(x, det1 / det0, det2 / det0)


def t_search(a1, a2, b1, b2):
    l = min(a1.x, a2.x)
    r = max(a1.x, a2.x)
    EPS = 1e-3

    while r - l >= EPS:
        m1 = l + (r - l) / 3
        m2 = r - (r - l) / 3
        a = point_on_line(m1, a1, a2)
        b = point_on_line(m2, a1, a2)
        d1 = distance_from_point_to_line(a, b1, b2)
        d2 = distance_from_point_to_line(b, b1, b2)
        if (d1 <= 1 or d2 <= 1):
            return True
        if (d1 > d2):
            l = m1
        else:
            r = m2
    return False

def lines_intersection(l1, l2):
    matrix0 = np.array([[l1.p1.a, l1.p1.b, l1.p1.c],
                        [l1.p2.a, l1.p2.b, l1.p2.c],
                        [l2.p1.a, l2.p1.b, l2.p1.c]])
    matrix1 = np.array([[l1.p1.d, l1.p1.b, l1.p1.c],
                        [l1.p2.d, l1.p2.b, l1.p2.c],
                        [l2.p1.d, l2.p1.b, l2.p1.c]])
    matrix2 = np.array([[l1.p1.a, l1.p1.d, l1.p1.c],
                        [l1.p2.a, l1.p2.d, l1.p2.c],
                        [l2.p1.a, l2.p1.d, l2.p1.c]])
    matrix3 = np.array([[l1.p1.a, l1.p1.b, l1.p1.d],
                        [l1.p2.a, l1.p2.b, l1.p2.d],
                        [l2.p1.a, l2.p1.b, l2.p1.d]])
    det0 = np.linalg.det(matrix0)
    det1 = np.linalg.det(matrix1)
    det2 = np.linalg.det(matrix2)
    det3 = np.linalg.det(matrix3)
    if det0 == 0:
        return None
    point = Point(det1 / det0, det2 / det0, det3 / det0)
    if (abs(l2.p2.a * point.x + l2.p2.b * point.y + l2.p2.c * point.z - l2.p2.d) > 0.5):
        return None
    return point

def linear_angle(p1, p2, p3, p4):
    p1 = Plane(p1, p2, p3)
    p2 = Plane(p2, p3, p4)
    return round(abs(p1.a * p2.a + p1.b * p2.b + p1.c * p2.c) / sqrt(p1.a ** 2 + p1.b ** 2 + p1.c ** 2) / sqrt(
        p2.a ** 2 + p2.b ** 2 + p2.c ** 2), 5)

#
# p = []
# for i in range(4):
#     x, y, z = [float(i) for i in input().split()]
#     p.append(Point(x, y, z))
#
# print(linear_angle(p[0], p[1], p[2], p[3]))