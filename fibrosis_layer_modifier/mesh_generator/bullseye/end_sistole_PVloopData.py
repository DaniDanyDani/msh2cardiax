
import numpy as np
import matplotlib.pyplot as plt

ta_file = './kerckhoffs_active_stress.txt'
time = 0.9 #s

active_tension = np.loadtxt(ta_file)*1.5 #1.25 #*0.76 #????
t = np.linspace(0,time,len(active_tension))
pressure = np.zeros(len(t))
id = np.where(active_tension>1)[0][0]
idt = np.where(t>0.4)[0][0]
idt = 200
szt = len(t)-idt
sz = len(active_tension[id:id+szt])

r = np.zeros((len(t), 3))

r[:,0] = t

r[0:idt,1] = np.linspace(0,1600,idt)
print(len(active_tension[id:id+szt]))
print(len(r[idt:id+szt,2]))
r[idt:idt+len(active_tension[id:id+szt]),2] = active_tension[id:id+szt]

r[:,2] = r[:,2]/np.max(r[:,2])

#plt.plot(r[:,0], r[:,1], label='P')
plt.plot(r[:,0], r[:,2], label='Ta')
plt.legend()
plt.show()
#pressure[:id+1] = (t[:id+1]/t[id])*2000
#r = np.array([t, pressure, active_tension]).T

outputFile = open('end_sistole_pvloop_data.txt','w')

outputFile.write('%d %f\n' % (len(r),time))

for i in range(len(r)):
	outputFile.write('%f %f %f\n' % (r[i,0],r[i,1],r[i,2]))


outputFile.close()
