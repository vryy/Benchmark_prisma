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

    return data[list(data.keys())[0]]['dx'], sum_p['px']

x1, y1 = extract_data('implicit-displacement-controlled/quad9.gid/Disp_Reaction.txt')

# print(sum_p['py'])

pylab.plot(x1, y1, label='implicit-dc', marker='.')
pylab.legend()
pylab.ylabel('F [N]')
pylab.xlabel('u [mm]')
pylab.show()
