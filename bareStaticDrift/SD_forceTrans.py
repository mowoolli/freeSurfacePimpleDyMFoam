#!/usr/bin/python

import sys, math

try:
	infilename = sys.argv[1]; outfilename = sys.argv[2]
except:
	print "Usage:", sys.argv[0], "infile outfile"; sys.exit(1)

rampTime = float(input("What is the simulation ramp time? "))

# infilename1 = "%s/forces1/0/forces.dat" % infilename
ifile1= open( infilename, 'r') # open file for reading
ofile = open(outfilename, 'w') # open file for writing

ofile.write("#time     forceX    forceY    momentN  \n")
ofile.write("#(s)      (N)       (N)       (N-m)  \n")

lines1 = ifile1.readlines()[3:]

time = []; forceX = []; forceY = []; momentN = []; 

t = int(0)

#gen = (line for line in lines1[1:len(lines1)] if line > rampTime)
#for line in gen:
for line in lines1[0:len(lines1)]:
	line = line.strip()
        line = line.replace('(',' ')
        line = line.replace(',',' ')
        line = line.replace(')',' ')
        columns = line.split()
        time.append(float(columns[0]))

#	for t in range(1, len(lines1)):
#      	print time
#       print line
	if time[t] < rampTime:
		deltaX = 1.50774*0.5*(time[t] - 9/(3.14159)*math.sin(3.14159*time[t]/(9)))
		deltaY = 1.50774*0.5*(time[t] - 9/(3.14159)*math.sin(3.14159*time[t]/(9)))*math.tan(-10*3.14159/180)
	else:
		deltaX = 1.50774*(time[t] - 9*0.5)
		deltaY = 1.50774*(time[t] - 9*0.5)*math.tan(-10*3.14159/180)
	
#	print "deltaX", deltaX
#	print "deltaY", deltaY
        forceXShip = float(columns[1]) + float(columns[4])
        forceYShip = float(columns[2]) + float(columns[5])
#       print "forceXShip", forceXShip
#       print "forceYShip", forceYShip
        momentNShip = (float(columns[12]) + float(columns[15])) + deltaY*forceXShip - deltaX*forceYShip
#       print "deltaY*forceXShip = ", deltaY*forceXShip
#       print "deltaX*forceYShip = ", deltaX*forceYShip
        forceX.append(forceXShip)
        forceY.append(forceYShip)
        momentN.append(momentNShip)
	t += 1
for i in range(0, len(time), 1):
	ofile.write('%9.7f %9.7f %9.7f %9.7f \n' \
        % (time[i], forceX[i], forceY[i], momentN[i]))

ifile1.close(); ofile.close()

