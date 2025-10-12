import sys
import matplotlib.pyplot as plt

# read data
monfile = sys.argv[1]
ifile = open(monfile, "r")

data = {}
data['time'] = []
data['dx'] = []
data['dy'] = []
data['dz'] = []
data['vx'] = []
data['vy'] = []
data['vz'] = []
data['ax'] = []
data['ay'] = []
data['az'] = []
for line in ifile.readlines():
    # print(line)
    words = line.split()
    # print(words)

    if words[0] == "#":
        continue

    if len(words) < 6:
        # data['time'].append(float(words[1]))
        # data['dx'].append(float(words[2]))
        # data['dy'].append(float(words[3]))
        # data['dz'].append(float(words[4]))
        pass
    else:
        time = float(words[1])
        if time >= 0.0:
            data['time'].append(time)
            data['dx'].append(float(words[2]))
            data['dy'].append(float(words[3]))
            data['dz'].append(float(words[4]))
            data['vx'].append(float(words[5]))
            data['vy'].append(float(words[6]))
            data['vz'].append(float(words[7]))
            data['ax'].append(float(words[8]))
            data['ay'].append(float(words[9]))
            data['az'].append(float(words[10]))

ifile.close()

print(len(data['time']))
print(len(data['dy']))
print(len(data['vx']))

# plot data
plt.plot(data['time'], data['dy'])
plt.xlabel('time')
plt.ylabel('dy')

plt.figure()
plt.plot(data['time'], data['vy'])
plt.xlabel('time')
plt.ylabel('vy')

plt.figure()
plt.plot(data['time'], data['ay'])
plt.xlabel('time')
plt.ylabel('ay')

plt.show()
