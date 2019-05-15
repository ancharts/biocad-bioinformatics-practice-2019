from math import sqrt, sin, cos
from collections import namedtuple
from xmlrpc.client import ServerProxy
from random import uniform, randint
from geometry import Point, Vector, Line, vmul, find_intersection_point, lines_intersection, t_search

pymol = ServerProxy(uri="http://localhost:9123/RPC2")
file_name = 'forccd.pdb'
st, fn = 61, 66
eps = 0.08

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


def parallel_transfer(v):
    # параллельный перенос точек на вектор v
    for i in range(len(crd)):
        crd[i] += v
    for i in range(len(frame)):
        frame[i] += v

def find_angle(atoms, target_coordinates, point1, point2):
    b, c = 0, 0
    for i in range(3):
        o = find_intersection_point(point1, point2, atoms[i])
        r = Vector(o, atoms[i])
        f = Vector(o, target_coordinates[i])

        if (r.len() < eps**10):
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            return Angle(0, 1)

        # если эти 3 точки лежат на одной прямой, что бы я ни делала, это ни на что не повлияет

        axis = Vector(point1, point2) / Vector(point1, point2).len()
        r1 = r / r.len()
        s1 = vmul(r1, axis)
        b += 2 * r.len() * (f * r1)
        c += 2 * r.len() * (f * s1)

    return Angle(c / sqrt(b ** 2 + c ** 2), b / sqrt(b ** 2 + c ** 2))


def read():
    with open(file_name, 'r') as fin:

        s = fin.readline()

        while not (s.startswith(data_type) and int(s[resi_seq_num_l:resi_seq_num_r]) == st):
            s = fin.readline()

        while s.startswith(data_type) and int(s[resi_seq_num_l:resi_seq_num_r]) < fn + 1:
            if (is_useful(s)):
                start.append(s)
            s = fin.readline()

        while not (s.startswith(data_type) and s[chain_identifier] == " "):
            s = fin.readline()

        while s.startswith(data_type):
            to_change.append(s)
            s = fin.readline()


def run_pymol():
    id1 = to_change[0][serial_number_l:serial_number_r].strip()
    id2 = to_change[-1][serial_number_l:serial_number_r].strip()
    s = ['delete all',
         'load %s' % file_name,
         'show_as cartoon, all',
         'select loop, id %s-%s' % (id1, id2),
         'show_as sticks, loop',
         'color argon, all',
         'color brightorange, loop',
         'color red, elem o',
         'color blue, elem n']
    for i in s:
        pymol.do(i)


def fill_target_crd():
    target_crd.append(Point(start[3]))  # n
    target_crd.append(Point(start[4]))  # ca
    target_crd.append(Point(start[5]))  # c


def fill_crd_and_frame():
    for i in range(len(to_change)):
        s = to_change[i]
        crd.append(Point(s))
        if is_useful(s):
            frame.append(Point(s))
            index.append(i)


def renew_crd(theta, id):
    # координаты
    axis1 = Vector(frame[id], frame[id + 1]) / Vector(frame[id], frame[id + 1]).len()

    for i in range(index[id + 2], len(crd), 1):
        old = Vector(crd[index[id]])
        v = Vector(crd[index[id]], crd[i])
        new = v * theta.cos - vmul(axis1, v) * theta.sin + axis1 * (axis1 * v) * (1 - theta.cos)
        crd[i] = (new + old).to_point()

def renew_frame():
    for i in range(len(to_change)):
        if is_useful(to_change[i]):
            frame.append(Point(crd[i].x, crd[i].y, crd[i].z))


def renew_pymol():
    for i in range(len(to_change)):
        num = to_change[i][serial_number_l:serial_number_r].strip()
        s = 'alter_state 1,id %s,(x,y,z)=(%f, %f, %f)' % (num, crd[i].x, crd[i].y, crd[i].z)
        pymol.do(s)
    pymol.do('zoom loop')
    pymol.do('orient loop')
    pymol.do('center loop')


def is_useful(s):
    atom = s[atom_name_l:atom_name_r].strip()
    if atom == 'CA' or atom == 'N' or atom == 'C':
        return True
    return False

# если что-то не работает, убери корень

def dist():
    return Vector(frame[-1], target_crd[-1]).len() ** 2 + \
           Vector(frame[-2], target_crd[-2]).len() ** 2 + \
           Vector(frame[-3], target_crd[-3]).len() ** 2

#
# def save_len():
#     for i in range(len(crd) - 1):
#         bond_length.append(Vector(crd[i], crd[i + 1]).len())

def sampling():
    for i in range(10):
        id = randint(0, len(frame) - 3)
        angle = uniform(0, 2 * 3.1415926)
        theta = Angle(sin(angle), cos(angle))
        renew_crd(theta, id)

def check_intersections():
    cnt = 0
    for i in range(len(crd) - 1):
        for j in range(i + 1, len(crd)):
            if (Vector(crd[i], crd[j]).len() < 1):
                cnt+=1

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
                            cnt+=1
                        if p is None:
                            cnt+= t_search(a1, a2, b1, b2)
    return cnt

def save_len():
    first_atom = int(to_change[0][serial_number_l:serial_number_r])
    last_atom = int(to_change[-1][serial_number_l:serial_number_r])
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

def renew_len():
    for i in range(len(crd) - 1):
        for bond in bond_length[i]:
            if i < bond[0]:
                v = Vector(crd[i], crd[bond[0]])
                k = bond[1] / v.len()
                tmp = crd[bond[0]]
                crd[bond[0]] = crd[i] + (v * k)
                for l in range(bond[0] + 1, len(crd)):
                    v = Vector(tmp, crd[l])
                    tmp = crd[l]
                    crd[l] = crd[l-1] + v
#исправление длин не работает, подумой пожалуста
#тебе надо сохранять старый вектор, а так ты ничего не меняешь


start = []  # заданные координаты начала и конца петли (строки)
to_change = []  # координаты атомов петли (строки)
crd = []  # координаты атомов петли
frame = []  # координаты N, CA, C атомов петли
target_crd = []  # целевые координаты
index = []  # номера N, CA, C в общей последовательности петли

# читаю из файла
read()

# добавляю целевые координаты
fill_target_crd()

# координаты атомов петли и основной массив
fill_crd_and_frame()

# переношу первый атом N на его место
parallel_transfer(Vector(frame[0], Point(start[0])))

# сохраняю длины
bond_length = [[] for i in range(len(crd))]
save_len()

# подключаю PyMOL
run_pymol()
#renew_pymol()

# рандомно гну
sampling()

renew_pymol()

# CCD
count = 0
while(True):
    for i in range(len(frame) - 2):

        # if to_change[index[i]][resi_l:resi_r].strip() == "PRO" or to_change[index[i+1]][resi_l:resi_r].strip() == "PRO":
        #     continue

        print(dist())
        # нахожу угол поворота
        theta = find_angle(frame[(len(frame) - 3):], target_crd, frame[i], frame[i + 1])

        # поворачиваю и обновляю координаты
        renew_crd(theta, i)
        frame = []
        renew_frame()

        if dist() < eps or count >= 5000:
            print(dist())
            print(count)
            renew_pymol()
            exit()

        count+=1
    renew_len()