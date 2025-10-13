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

def extract_data(filename, node_id, sd, sp):
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
            data[node]['dy'].append(sd*dy)
            data[node]['dz'].append(dz)
            data[node]['px'].append(px)
            data[node]['py'].append(sp*py)
            data[node]['pz'].append(pz)

    ifile.close()

    for node in data:
        data[node]['py'] = div_list(data[node]['py'], 1.0e3)

    return data[node_id]['dy'], data[node_id]['py']

    # sum_p = {}
    # sum_p['px'] = []
    # sum_p['py'] = []
    # sum_p['pz'] = []
    # nnode = len(data)
    # for node, val in data.items():
    #     sum_p['px'] = add_list(sum_p['px'], val['px'])
    #     sum_p['py'] = add_list(sum_p['py'], val['py'])
    #     sum_p['pz'] = add_list(sum_p['pz'], val['pz'])

    # sum_p['px'] = div_list(sum_p['px'], 1.0e3)
    # sum_p['py'] = div_list(sum_p['py'], 1.0e3)
    # sum_p['pz'] = div_list(sum_p['pz'], 1.0e3)

    # return data[list(data.keys())[0]]['dy'], sum_p['py']

def extract_lambda(filename):
    ifile = open(filename, 'r')

    read = False
    lambda_list = []
    for line in ifile.readlines():
        if read:
            words = line.split()
            lambda_list.append(float(words[1]))

        if not read:
            read = True

    ifile.close()

    return lambda_list

lambda1_q4 = extract_lambda('mesh_1.gid/lambda_arc_length.log')
x1_q4_1, y1_q4_1 = extract_data('mesh_1.gid/Disp_Reaction.txt', 2000, -1.0, -1.0)
x1_q4_2, y1_q4_2 = extract_data('mesh_1.gid/Disp_Reaction.txt', 4068, 1.0, -1.0)
lambda1_q9 = extract_lambda('mesh_1_q9.gid/lambda_arc_length.log')
x1_q9_1, y1_q9_1 = extract_data('mesh_1_q9.gid/Disp_Reaction.txt', 7862, -1.0, -1.0)
x1_q9_2, y1_q9_2 = extract_data('mesh_1_q9.gid/Disp_Reaction.txt', 16020, 1.0, -1.0)

# print(lambda1_q4)
# print(lambda1_q9)
P = 1.0e3

pylab.subplot(211)
pylab.title("Load-displacement curve (single edge notched beam)")
# pylab.plot(x1_q4_1, y1_q4_1, label='3955 (middle)', marker='.')
pylab.plot(x1_q4_1, lambda1_q4, label='3955 (middle)', marker='.')
# pylab.plot(x1_q9_1, y1_q9_1, label='3955 (q9) (middle)', marker='.')
pylab.plot(x1_q9_1, lambda1_q9, label='3955 (q9) (middle)', marker='.')
pylab.axhline(y=40.0, color='r', linestyle='--')
pylab.xlim([0.0, 0.2])
pylab.ylim([0.0, 45.0])
pylab.legend()
pylab.ylabel('F [kN]')
pylab.xlabel('u [mm]')
pylab.subplot(212)
# pylab.plot(x1_q4_2, y1_q4_2, label='3955 (side)', marker='.')
pylab.plot(x1_q4_2, lambda1_q4, label='3955 (side)', marker='.')
# pylab.plot(x1_q9_2, y1_q9_2, label='3955 (q9) (side)', marker='.')
pylab.plot(x1_q9_2, lambda1_q9, label='3955 (q9) (side)', marker='.')
pylab.axhline(y=40.0, color='r', linestyle='--')
pylab.xlim([-0.1, 1.0])
pylab.ylim([0.0, 45.0])
pylab.legend()
pylab.ylabel('F [kN]')
pylab.xlabel('u [mm]')
pylab.show()
