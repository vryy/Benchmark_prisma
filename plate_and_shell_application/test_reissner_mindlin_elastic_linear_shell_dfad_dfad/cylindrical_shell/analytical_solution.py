import math

class CylindricalShellSolution:

    def __init__(self, filename, wscale = 1.0, ucolumn = 1):

        self.w = []
        self.u = []

        ifile = open(filename, "r")

        cnt = 0
        for line in ifile.readlines():

            if cnt == 0: # skip first line
                cnt += 1
                continue

            words = line.split();
            self.w.append(float(words[0]) / wscale)
            self.u.append(float(words[ucolumn]))

        ifile.close()

    def Interpolate(self, w):

        for i in range(0, len(self.w)-1):
            if (w >= self.w[i]) and (w < self.w[i+1]):
                a = (w - self.w[i]) / (self.w[i+1] - self.w[i])
                v = (1.0 - a) * self.u[i] + a * self.u[i+1]
                return v

        if w < self.w[0]:
            return self.u[0]

        if w > self.w[-1]:
            return self.u[-1]
