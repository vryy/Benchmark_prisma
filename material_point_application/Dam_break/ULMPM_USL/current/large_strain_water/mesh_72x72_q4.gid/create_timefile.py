import math

time = 0.0
hx = 1.0/72
kappa = 2.0e6
rho = 997.5
delta_time = 0.1*hx/math.sqrt(kappa/rho)
print("delta_time: " + str(delta_time))
total_time = 5.0
sample_output = 100

step = 0
itime = open("time.txt", "w")
itime.write("step\ttime\n")
while time < total_time:
    time = time + delta_time
    print("#########################################################")
    print("######### TIME STEP " + str(time) + " BEGIN #############")

    if step % sample_output == 0:

        itime.write(str(step) + "\t" + str(time) + "\n")
        itime.flush()

    step = step + 1
