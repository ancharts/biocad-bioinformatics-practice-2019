from math import sin, cos, acos, atan2, degrees, sqrt
from collections import namedtuple
from xmlrpc.client import ServerProxy
from random import randint
from scipy.optimize import fsolve
from geometry import Point, Vector, Line, vmul, lines_intersection, t_search, linear_angle
import datetime
import time

time0 = datetime.datetime.now()

pymol = ServerProxy(uri="http://localhost:9123/RPC2")
file_name = '77.pdb'
st, fn = 10, 16
chain = 'H'
iterations = 5
samples = 1

# константы для чтения
# https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
data_type = "ATOM"
serial_number_l = 6
serial_number_r = 11
atom_name_l = 12
atom_name_r = 16
resi_l = 17
resi_r = 20
chain_identifier = 21
resi_seq_num_l = 22
resi_seq_num_r = 26
x_crd_l = 30
x_crd_r = 38
y_crd_l = 38
y_crd_r = 46
z_crd_l = 46
z_crd_r = 54
first_l = 11
first_r = 16
second_l = 16
second_r = 21
third_l = 21
third_r = 26
fourth_l = 26
fourth_r = 31

Angle = namedtuple('Angle', ['sin', 'cos'])


def read():
    with open(file_name, 'r') as fin:

        s = fin.readline()

        while not (s.startswith(data_type) and int(s[resi_seq_num_l:resi_seq_num_r]) == st and s[chain_identifier].strip() == chain):
            s = fin.readline()

        while s.startswith(data_type) and int(s[resi_seq_num_l:resi_seq_num_r]) < fn + 1:
            loop.append(s)
            s = fin.readline()


def renew_crd(t0, t1, t2, axis):
    # координаты

    # первая ось
    axis1 = Vector(crd[index[axis[0]]], crd[index[axis[1]]]).normalize()
    for i in range(index[axis[0]] + 1, index[axis[1]]):

        old = Vector(crd[index[axis[0]]])
        v = Vector(crd[index[axis[0]]], crd[i])
        new = v * t0.cos - vmul(axis1, v) * t0.sin + axis1 * (axis1 * v) * (1 - t0.cos)
        crd[i] = (new + old).to_point()

    # вторая ось
    axis1 = Vector(crd[index[axis[1]]], crd[index[axis[2]]]).normalize()
    for i in range(index[axis[1]] + 1, index[axis[2]]):
        old = Vector(crd[index[axis[1]]])
        v = Vector(crd[index[axis[1]]], crd[i])
        new = v * t1.cos - vmul(axis1, v) * t1.sin + axis1 * (axis1 * v) * (1 - t1.cos)
        crd[i] = (new + old).to_point()
    axis1 = Vector(crd[index[axis[2]]], crd[index[axis[0]]]).normalize()

    # третья ось
    for i in range(index[axis[0]] + 1, index[axis[2]]):
        old = Vector(crd[index[axis[2]]])
        v = Vector(crd[index[axis[2]]], crd[i])
        new = v * (-t2.cos) - vmul(axis1, v) * t2.sin + axis1 * (axis1 * v) * (1 + t2.cos)
        crd[i] = (new + old).to_point()

    renew_len()


def renew_len():
    for i in range(len(crd) - 1):
        for bond in bond_length[i]:
            v = Vector(crd[i], crd[bond[0]])
            k = bond[1] / v.len()
            crd[bond[0]] = crd[i] + (v * k)
            if (bond[0] < i):
                for l in range(i + 1, len(crd)):
                    if not loop[l][resi_l:resi_r].strip() == 'PRO':
                        v = Vector(crd[l - 1], crd[l])
                        crd[l] = crd[l - 1] + v
            else:
                for l in range(i - 1, -1, -1):
                    if not (loop[l][resi_l:resi_r].strip() == 'PRO'):
                        v = Vector(crd[l - 1], crd[l])
                        crd[l] = crd[l - 1] + v


def is_useful(s):
    atom = s[atom_name_l:atom_name_r].strip()
    if atom == 'CA':
        return True
    return False


def fill_crd():
    for i in range(len(loop)):
        s = loop[i]
        crd.append(Point(s))
        if is_useful(s):
            index.append(i)


def save_len():
    first_atom = int(loop[0][serial_number_l:serial_number_r])
    last_atom = int(loop[-1][serial_number_l:serial_number_r])
    with open(file_name, 'r') as fin:

        s = fin.readline()

        while not (s.startswith("CONECT")):
            s = fin.readline()

        while s.startswith("CONECT"):
            atom = int(s[serial_number_l:serial_number_r])
            if (atom >= first_atom and atom <= last_atom):
                atom -= first_atom
                first = s[first_l:first_r]
                second = s[second_l:second_r]
                third = s[third_l:third_r]
                fourth = s[fourth_l:fourth_r]
                if (first.strip() != "" and int(first) >= first_atom and int(first) <= last_atom):
                    first = int(first) - first_atom
                    bond_length[atom].append((first, Vector(crd[atom], crd[first]).len()))

                if (second.strip() != "" and int(second) >= first_atom and int(second) <= last_atom):
                    second = int(second) - first_atom
                    bond_length[atom].append((second, Vector(crd[atom], crd[second]).len()))

                if (third.strip() != "" and int(third) >= first_atom and int(third) <= last_atom):
                    third = int(third) - first_atom
                    bond_length[atom].append((third, Vector(crd[atom], crd[third]).len()))

                if (fourth.strip() != "" and int(fourth) >= first_atom and int(fourth) <= last_atom):
                    fourth = int(fourth) - first_atom
                    bond_length[atom].append((fourth, Vector(crd[atom], crd[fourth]).len()))

            s = fin.readline()


def generate():
    first, second, third = -1, -1, -1
    while not (first < second and second < third and third < number and first != -1 and second != -1 and third != -1):
        first = randint(0, number - 1)
        if loop[index[first]][resi_l:resi_r].strip() == 'PRO':
            first = -1
        if first <= number - 3:
            second = randint(first + 1, number - 1)
            if loop[index[second]][resi_l:resi_r].strip() == 'PRO':
                second = -1
        else:
            second = -1
        if second <= number - 2:
            third = randint(second + 1, number - 1)
            if loop[index[third]][resi_l:resi_r].strip() == 'PRO':
                third = -1
        else:
            third = -1
    return [first, second, third]


def equations(prev_res):
    t0c, t1c, t2c = prev_res

    if not (t0c**2 <= 1 and t1c**2 <=1 and t2c**2 <=1):
        return (100, 100, 100)

    f1 = a[0] + b[0] * (t2c * cos(delta[2]) - sqrt(1-t2c**2) * sin(delta[2])) + c[0] * t0c + d[0] * t0c * (t2c * cos(delta[2]) - sqrt(1-t2c**2) * sin(delta[2])) + e[0] * sqrt(1-t0c**2) * (sqrt(1-t2c**2) * cos(delta[2]) + t2c * sin(delta[2]))
    f2 = a[1] + b[1] * (t0c * cos(delta[0]) - sqrt(1-t0c**2) * sin(delta[0])) + c[1] * t1c + d[1] * t1c * (t0c * cos(delta[0]) - sqrt(1-t0c**2) * sin(delta[0])) + e[1] * sqrt(1-t1c**2) * (sqrt(1-t0c**2) * cos(delta[0]) + t0c * sin(delta[0]))
    f3 = a[2] + b[2] * (t1c * cos(delta[1]) - sqrt(1-t1c**2) * sin(delta[1])) + c[2] * t2c + d[2] * t2c * (t1c * cos(delta[1]) - sqrt(1-t1c**2) * sin(delta[1])) + e[2] * sqrt(1-t2c**2) * (sqrt(1-t1c**2) * cos(delta[1]) + t1c * sin(delta[1]))

    # t0c, t1c, t2c = prev_res
    # f1 = a[0] + b[0] * (t2c * cos(delta[2]) - sqrt(1-t2c**2) * sin(delta[2])) + c[0] * t0c + d[0] * t0c * (t2c * cos(delta[2]) - sqrt(1-t2c**2) * sin(delta[2])) + e[0] * sqrt(1-t0c**2) * (sqrt(1-t2c**2) * cos(delta[2]) + t2c * sin(delta[2]))
    # f2 = a[1] + b[1] * (t0c * cos(delta[0]) - sqrt(1-t0c**2) * sin(delta[0])) + c[1] * t1c + d[1] * t1c * (t0c * cos(delta[0]) - sqrt(1-t0c**2) * sin(delta[0])) + e[1] * sqrt(1-t1c**2) * (sqrt(1-t0c**2) * cos(delta[0]) + t0c * sin(delta[0]))
    # f3 = a[2] + b[2] * (t1c * cos(delta[1]) - sqrt(1-t1c**2) * sin(delta[1])) + c[2] * t2c + d[2] * t2c * (t1c * cos(delta[1]) - sqrt(1-t1c**2) * sin(delta[1])) + e[2] * sqrt(1-t2c**2) * (sqrt(1-t1c**2) * cos(delta[1]) + t1c * sin(delta[1]))

    # f4 = t0s ** 2 + t0c ** 2 - 1
    # f5 = t1s ** 2 + t1c ** 2 - 1
    # f6 = t2s ** 2 + t2c ** 2 - 1
    return (f1, f2, f3)


def run_pymol():
    s = ['delete all',
         'load %s' % file_name,
         'show_as cartoon, all',
         'select loop, /%s//%s/%s-%s' % (file_name[:len(file_name)-4], chain, st, fn),
         'show_as сartoon, loop',
         'color white, loop',
         'color red, elem o',
         'color blue, elem n',
         'orient loop',
         'deselect']
    for i in s:
        pymol.do(i)


def renew_pymol():
    for i in range(len(loop)):
        resi = loop[i][resi_l:resi_r].strip()
        resi_num = loop[i][resi_seq_num_l:resi_seq_num_r].strip()
        atom = loop[i][atom_name_l:atom_name_r].strip()
        s = 'alter_state 1,/%s//%s/%s`%s/%s,(x,y,z)=(%f, %f, %f)' % \
            (file_name[:len(file_name) - 4], chain, resi, resi_num, atom, crd[i].x, crd[i].y, crd[i].z)
        pymol.do(s)
    pymol.do("orient loop")


def check_intersections():

    for i in range(len(bond_length) - 1):
        for j in bond_length[i]:
            for k in range(i + 1, len(bond_length)):
                for l in bond_length[k]:

                    a1 = crd[i]
                    a2 = crd[j[0]]
                    b1 = crd[k]
                    b2 = crd[l[0]]

                    if a1 != b1 and a1 != b2 and a2 != b1 and a2 != b2:
                        p = lines_intersection(Line(a1, a2), Line(b1, b2))
                        if not p is None and p.between(a1, a2) and p.between(b1, b2):
                            print("!!!!!!!!")
                            return True
                        if p is None:
                            return t_search(a1, a2, b1, b2)


for sample in range(samples):
    loop = []  # координаты атомов петли (строки)
    crd = []  # координаты атомов петли
    index = []  # номера CA в общей последовательности петли

    # считываю строки
    read()

    # координаты атомов петли и основной массив
    fill_crd()

    bond_length = [[] for i in range(len(crd))]  # длины связей

   # save_len()

    number = len(index)

    # запускаю PyMOL
    run_pymol()

    ideal = crd.copy()

    iteration = iterations

    while(iteration > 0):
        tmp = crd.copy()
        # случайным образом выбираю три последовательных атома CA в цепи
        axis = generate()

        alpha = []
        eta = []
        xi = []
        theta = []
        delta = []

        z = []
        r_tau = []
        r_sigma = []

        # векторы
        for i in range(3):
            z.append(Vector(crd[index[axis[i]]], crd[index[axis[(i + 1) % 3]]]).normalize())
            r_tau.append(Vector(crd[index[axis[i]]], crd[index[axis[i]] + 1]).normalize())
            r_sigma.append(Vector(crd[index[axis[(i + 1) % 3]]], crd[index[axis[(i + 1) % 3]] - 1]).normalize())

        # углы
        for i in range(3):
            theta.append(111*3.1415926/180)
            alpha.append(acos(z[i] * z[i - 1]))
            eta.append(acos(z[i] * r_tau[i]))
            xi.append(acos(Vector(Point(-z[i].x, -z[i].y, -z[i].z)) * r_sigma[i]))
            delta.append(acos(linear_angle(crd[index[axis[i]] + 1], crd[index[axis[i]]], crd[index[axis[(i + 1) % 3]]],
                                           crd[index[axis[(i + 1) % 3]] - 1])))

        # коэффициенты
        a = [0.0] * 3
        b = [0.0] * 3
        c = [0.0] * 3
        d = [0.0] * 3
        e = [0.0] * 3
        for i in range(3):
            a[i] = -cos(theta[i]) - cos(eta[i]) * cos(xi[i - 1]) * cos(alpha[i])
            b[i] = sin(alpha[i]) * sin(xi[i - 1]) * cos(eta[i])
            c[i] = sin(alpha[i]) * cos(xi[i - 1]) * sin(eta[i])
            d[i] = cos(alpha[i]) * sin(xi[i - 1]) * sin(eta[i])
            e[i] = sin(xi[i - 1]) * sin(eta[i])

        # решения системы
        t0c, t1c, t2c = fsolve(equations, (1, 1, 1), xtol= 1.49012e-16, maxfev = 1000000000)
        print(t0c, t1c, t2c)
        print(equations((t0c, t1c, t2c)))
        t0 = Angle(sqrt(1-t0c**2), t0c)
        t1 = Angle(sqrt(1-t1c**2), t1c)
        t2 = Angle(sqrt(1-t2c**2), t2c)

        # for i in range(len(crd)):
        #     print(loop[i][atom_name_l:atom_name_r], crd[i])
        # print()

        # tau0 = atan2(t0s, t0c)
        # tau1 = atan2(t1s, t1c)
        # tau2 = atan2(t2s, t2c)

        # print(tau0, tau1, tau2)

        # print(t0, t0.cos ** 2 + t0.sin ** 2)

        # поворачиваю на полученные углы
        renew_crd(t0, t1, t2, axis)

        # for i in range(len(crd)):
        #     print(loop[i][atom_name_l:atom_name_r], crd[i])
        # print()

        # if check_intersections():
        #     crd = tmp.copy()
        #     continue

        # dist = 0
        # for i in range(len(crd)):
        #     dist += Vector(crd[i], ideal[i]).len()
        # print(dist / len(crd))
        iteration-=1
        renew_len()
        renew_pymol()
    # pymol.do('load 5mo3_fab.pdb')
    # pymol.do('select l1, /5mo3_fab//H/99-109')
    # pymol.do('color red, l1')
    # # pymol.do('orient l1')
    # pymol.do('align 77, 5mo3_fab, cycles = 0')
    # pymol.do('deselect')
    # pymol.do('save %d.png' % sample)
    # pymol.do("save /home/uue/Documents/5mo3_min/sample_H_cdr3_%d.pdb, all, state = 0" % (sample))
    # print(sample)
    time.sleep(5)

print(datetime.datetime.now() - time0)