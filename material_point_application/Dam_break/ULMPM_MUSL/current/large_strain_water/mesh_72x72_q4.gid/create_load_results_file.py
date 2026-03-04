import os

igid = open("load_results.bch", "w")
igid.write("mescape\nPostprocess\n")
igid.write("Files ReadMultiple {\n")

prefix  = "mesh_72x72_q4_"
postfix = ".post.bin"

itime = open("time.txt", "r")
itime.readline()

first = True
while True:
    line = itime.readline()

    if not line:
        break

    words = line.split()
    time = float(words[1])

    filename = os.getcwd() + "/" + prefix + str(time*1.0e3) + postfix

    igid.write(filename + "\n")

igid.write("}\n")
igid.close()
itime.close()
