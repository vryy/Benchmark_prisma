import numpy
import pylab

def read_data(data_file, scale=1.0):
    fid = open(data_file, 'r')
    x = []
    y = []

    lines = fid.readlines()
    first_line = True
    for line in lines:
        if first_line: # skip first line
            first_line = False
            continue
        tmp = line.split()
        if tmp[0] != '#':
            x.append(float(tmp[0]))
            y.append(float(tmp[1])*scale)

    fid.close()
    return x, y


x1, y1 = read_data("displacement.log", 1e3)
pylab.title('Displacement profile')
pylab.xlabel('Time (s)')
pylab.ylabel('uy (mm)')
pylab.plot(x1, y1)

pylab.figure()
x2, y2 = read_data("water_pressure.log")
pylab.title('Water pressure profile')
pylab.xlabel('Time (s)')
pylab.ylabel('wp (Pa)')
pylab.plot(x2, y2)

pylab.show()
