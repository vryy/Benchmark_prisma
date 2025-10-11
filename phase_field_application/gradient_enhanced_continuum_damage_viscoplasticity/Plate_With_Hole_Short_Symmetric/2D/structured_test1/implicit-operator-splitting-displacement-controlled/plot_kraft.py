import pylab

def add_list(a, b):
    if len(a) != len(b):
        a = b.copy()
        return a
    n = len(b)
    for i in range(0, n):
        a[i] += b[i]
    return a

def div_list(a, m):
    n = len(a)
    for i in range(0, n):
        a[i] /= m
    return a

def extract_data(filename):
    ifile = open(filename, 'r')

    data = {}
    read_flag = False
    for line in ifile.readlines():
        words = line.split()
        if len(words) > 0:
            # print(words[0])

            if words[0] == 'node':
                read_flag = True
                continue
        else:
            read_flag = False

        if read_flag:
            # print(words[0])
            node = int(words[0])
            dx = float(words[1])
            dy = float(words[2])
            dz = float(words[3])
            px = float(words[4])
            py = float(words[5])
            pz = float(words[6])
            if not (node in data):
                data[node] = {}
                data[node]['dx'] = []
                data[node]['dy'] = []
                data[node]['dz'] = []
                data[node]['px'] = []
                data[node]['py'] = []
                data[node]['pz'] = []
            data[node]['dx'].append(dx)
            data[node]['dy'].append(dy)
            data[node]['dz'].append(dz)
            data[node]['px'].append(px)
            data[node]['py'].append(py)
            data[node]['pz'].append(pz)

    ifile.close()

    # node = 14
    # pylab.plot(data[node]['dy'], data[node]['py'])
    # pylab.xlabel('u [mm]')
    # pylab.show()

    sum_p = {}
    sum_p['px'] = []
    sum_p['py'] = []
    sum_p['pz'] = []
    nnode = len(data)
    for node, val in data.items():
        sum_p['px'] = add_list(sum_p['px'], val['px'])
        sum_p['py'] = add_list(sum_p['py'], val['py'])
        sum_p['pz'] = add_list(sum_p['pz'], val['pz'])

    sum_p['px'] = div_list(sum_p['px'], 1.0)
    sum_p['py'] = div_list(sum_p['py'], 1.0)
    sum_p['pz'] = div_list(sum_p['pz'], 1.0)

    return data[list(data.keys())[0]]['dy'], sum_p['py']

def write_data(filename, x, y):
    ifile = open(filename, 'w')

    ifile.write("du\treaction\n")

    for i in range(0, len(x)):
        ifile.write(str(x[i]) + '\t' + str(y[i]) + '\n')

    ifile.close()

x0, y0 = extract_data('../../structured_du=1.0e-3/implicit-displacement-controlled/rate_independent/mesh_200.gid/Disp_Reaction.txt')

x1, y1 = extract_data('dt=1.0e-3/mesh_200.gid/Disp_Reaction.txt')
x2, y2 = extract_data('dt=1.0e-4/mesh_200.gid/Disp_Reaction.txt')
x3, y3 = extract_data('dt=1.0e-5/mesh_200.gid/Disp_Reaction.txt')
x4, y4 = extract_data('dt=1.0e-6/mesh_200.gid/Disp_Reaction.txt')

x1i, y1i = extract_data('../implicit-displacement-controlled/dt=1.0e-3[failed]/mesh_200.gid/Disp_Reaction.txt')
x2i, y2i = extract_data('../implicit-displacement-controlled/dt=1.0e-4/mesh_200.gid/Disp_Reaction.txt')
x3i, y3i = extract_data('../implicit-displacement-controlled/dt=1.0e-5/mesh_200.gid/Disp_Reaction.txt')
x4i, y4i = extract_data('../implicit-displacement-controlled/dt=1.0e-6/mesh_200.gid/Disp_Reaction.txt')

# # x5, y5 = extract_data('dt=2.0e-4/mesh_200.gid/Disp_Reaction.txt')
# x6, y6 = extract_data('dt=1.0e-2/mesh_200.gid/Disp_Reaction.txt')
# x7, y7 = extract_data('dt=2.0e-3/mesh_200.gid/Disp_Reaction.txt')
# x8, y8 = extract_data('dt=1.0e-7/mesh_200.gid/Disp_Reaction.txt')

# write_data('disp_reaction_200_rate_independent.txt', x0, y0)
# write_data('disp_reaction_200_rate=1e-1.txt', x6, y6)
# write_data('disp_reaction_200_rate=1e1.txt', x2, y2)
# write_data('disp_reaction_200_rate=1e2.txt', x3, y3)
# write_data('disp_reaction_200_rate=1e3.txt', x4, y4)

# print(sum_p['py'])

# pylab.plot(x0, y0, label='rate-independent', marker='.')
# # pylab.plot(x6, y6, label='rate=0.1 mm/s', marker='.')
# # pylab.plot(x7, y7, label='rate=0.5 mm/s', marker='.')
# pylab.plot(x1, y1, label='rate=1.0 mm/s', marker='.')
# # pylab.plot(x5, y5, label='rate=5.0 mm/s', marker='.')
# pylab.plot(x2, y2, label='rate=10.0 mm/s', marker='.')
# pylab.plot(x3, y3, label='rate=100.0 mm/s', marker='.')
# pylab.plot(x4, y4, label='rate=1000.0 mm/s', marker='.')

pylab.plot(x0, y0, label='rate-independent', marker='.')
# pylab.plot(x6, y6, label='rate=0.1 mm/s', marker='.')
# pylab.plot(x7, y7, label='rate=0.5 mm/s', marker='.')
pylab.plot(x1, y1, label='rate=1.0 mm/s (ops)')
pylab.plot(x1i, y1i, label='rate=1.0 mm/s (implicit)', marker='x')
# pylab.plot(x5, y5, label='rate=5.0 mm/s', marker='.')
pylab.plot(x2, y2, label='rate=10.0 mm/s (ops)')
pylab.plot(x2i, y2i, label='rate=10.0 mm/s (implicit)', marker='x')
pylab.plot(x3, y3, label='rate=100.0 mm/s (ops)')
pylab.plot(x3i, y3i, label='rate=100.0 mm/s (implicit)', marker='x')
pylab.plot(x4, y4, label='rate=1000.0 mm/s (ops)')
pylab.plot(x4i, y4i, label='rate=1000.0 mm/s (implicit)', marker='x')

# pylab.plot(x8, y8, label='rate=10000.0 mm/s', marker='.')
pylab.legend()
pylab.ylabel('F [N]')
pylab.xlabel('u [mm]')
pylab.show()
