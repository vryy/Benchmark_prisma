import numpy
import pylab

data_files = []

data_files.append('temperature.log')

plot_data = []
plot_legend = []

for data in data_files:
    fid = open(data, 'r')
    x = []
    y = []
    y2 = []
    y3 = []

    lines = fid.readlines()
    first_line = True
    for line in lines:
        if first_line: # skip first line
            first_line = False
            continue
        tmp = line.split()
        if tmp[0] != '#':
            x.append(float(tmp[0]))
            y.append(float(tmp[1]))
            y2.append(float(tmp[2]))
            y3.append(float(tmp[3])*1000)

    fid.close()
    p, = pylab.plot(x, y, label = data)
    p2, = pylab.plot(x, y2, label = data)
    plot_data.append(p)
    plot_legend.append(data.split('.')[0] + " far")
    plot_data.append(p2)
    plot_legend.append(data.split('.')[0] + " near")

pylab.title('Temperature profile')
pylab.xlabel('Time (s)')
pylab.ylabel('Temperature (oC)')
pylab.legend(plot_data, plot_legend, loc="lower right")
pylab.savefig("temperature.png")

pylab.figure()
pylab.title('Displacement profile')
pylab.xlabel('Time (s)')
pylab.ylabel('Displacement (mm)')
pylab.plot(x, y3)
pylab.savefig("displacement.png")

pylab.show()
