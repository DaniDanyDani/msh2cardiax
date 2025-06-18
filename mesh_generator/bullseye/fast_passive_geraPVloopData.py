
import numpy as np

ta_file = './active_tension_900ms.txt'
time = 0.9 #s

active_tension = np.loadtxt(ta_file)*2.0 #*0.76 #????
t = np.linspace(0,time,len(active_tension))
pressure = np.zeros(len(t))
id = np.where(active_tension>1)[0][0]

passive_time = np.linspace(0, t[id], 11)

passive_press = (passive_time/t[id])*2000
print passive_press, passive_time
passive_ta = []
for tt  in passive_time:
	passive_ta.append(active_tension[ np.where(t>=tt)[0][0] ])

print passive_ta

passive_r =np.array([passive_time, passive_press, passive_ta]).T

pressure[:id+1] = (t[:id+1]/t[id])*2000
#pressure[:500] = (t[:500]/0.5)*1500
#r = np.array([t[::2], pressure[::2], active_tension[::2]]).T
r = np.array([t, pressure, active_tension]).T
#r = np.array([t[::2], pressure[::2]*0.001, active_tension[::2]*0.001]).T

outputFile = open('pvloop_data.txt','w')

outputFile.write('%d %f\n' % (len(passive_time) + len(r[id+1:]),time))

for i in xrange(len(passive_time)):
	outputFile.write('%f %f %f\n' % (passive_r[i,0],passive_r[i,1],passive_r[i,2]))



for i in range(id+1, len(r)):
	outputFile.write('%f %f %f\n' % (r[i,0],r[i,1],r[i,2]))


outputFile.close()
