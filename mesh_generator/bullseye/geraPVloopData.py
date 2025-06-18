
import numpy as np
import matplotlib.pyplot as plt

#ta_file = './old_active_tension_900ms.txt'
ta_file = 'kerckhoffs_active_stress.txt'
time = 0.9 #s

active_tension = np.loadtxt(ta_file) #1.25 #*0.76 #????
t = np.linspace(0,time,len(active_tension))

plt.plot(t, active_tension, label='Ta')
plt.show()

pressure = np.zeros(len(t))
id = np.where(active_tension>1/50000.)[0][0]

pressure[:id+1] = (t[:id+1]/t[id])*2000
#pressure[:500] = (t[:500]/0.5)*1500
#r = np.array([t[::2], pressure[::2], active_tension[::2]]).T
r = np.array([t, pressure, active_tension]).T
#r = np.array([t[::2], pressure[::2]*0.001, active_tension[::2]*0.001]).T

outputFile = open('pvloop_data.txt','w')

outputFile.write('%d %f\n' % (len(r),time))

for i in range(len(r)):
	outputFile.write('%f %f %f\n' % (r[i,0],r[i,1],r[i,2]))


outputFile.close()

plt.plot(r[:,0], r[:,1], label='P')
plt.plot(r[:,0], r[:,2], label='Ta')
plt.legend()
plt.show()
