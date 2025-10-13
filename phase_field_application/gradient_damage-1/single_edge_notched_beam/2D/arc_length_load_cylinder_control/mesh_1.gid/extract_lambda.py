ifile = open('residuum_arc_length.log', 'r')

read = True
read_1 = False
read_2 = False
time_list = []
lambda_list = []
for line in ifile.readlines():
    words = line.split()
    # print(words)

    if '----' in words[0]:
        lambda_list.append(lambda_)
        read_2 = False

    if words[0] == 'time:':
        read_1 = True
        time_list.append(float(words[1]))
        continue

    if read_1:
        read_1 = False
        read_2 = True
        continue

    if read_2:
        if len(words) > 2:
            lambda_ = float(words[4])

ifile.close()

print(time_list)
print(lambda_list)

ifile2 = open("lambda_arc_length.log", 'w')
ifile2.write("time\tlambda\n")
for i in range(0, len(time_list)):
    ifile2.write(str(time_list[i]) + '\t' + str(lambda_list[i]) + '\n')
ifile2.close()
