import pylab
import numpy as np

def extract_data(filename):
    ifile = open(filename, 'r')

    data = {}
    first_line = True
    coordinates = []
    times = []
    time_values = []
    for line in ifile.readlines():
        words = line.split()

        if first_line:
            # read the coordinates
            for i in range(1, len(words)):
                coordinates.append(float(words[i]))
            first_line = False
        else:
            # read the time values
            time = float(words[0])
            values = []
            for i in range(1, len(words)):
                values.append(float(words[i]))
            times.append(time)
            time_values.append(values)

    ifile.close()

    return coordinates, times, time_values

coordinates, times, time_values = extract_data("air_pressure.log")
# print(times)

time_to_plot = np.arange(2, 7201, 300) # [2.0, 7201.0]
time_to_plot = np.append(time_to_plot, 7201)
# print(time_to_plot)

for t in time_to_plot:
    for i in range(0, len(times)):
        if times[i] == t:
            label = 'time=' + str(t)
            # pylab.plot(coordinates, time_values[i], label=label)
            pylab.plot(time_values[i], coordinates, label=label)

pylab.legend()
pylab.xlabel('air pressure (Pa)')
pylab.ylabel('coordinates (m)')
pylab.show()
