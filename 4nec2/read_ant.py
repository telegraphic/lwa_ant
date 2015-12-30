import numpy as np
import pylab as plt
import h5py

ant_data = np.zeros([15, 72*5+1, 36*5+1])


grid = np.zeros([72*5+1, 36*5+1])
theta = np.arange(-90, 91, 1)
phi   = np.arange(0, 361, 1)
freqs = np.arange(20, 91, 5)

ii = 0
for ii in range(0,15):
    freq = 20 + ii * 5
    print "Opening lwa-3d-%imhz.txt" % freq
    a = np.genfromtxt('lwa-4nec2-sim-output/lwa-3d-%imhz.txt' % freq)

    for line in a:
        x, y, z = line
        
        yidx = np.argmin(np.abs(theta - y))
        xidx = np.argmin(np.abs(phi - x))
        
        #print xidx, yidx
        
        grid[xidx, yidx] = z
    
    ant_data[ii] = grid
    ii += 1

h5 = h5py.File("antenna_data.h5", "w")
h5.create_dataset("data", data=ant_data)
h5.create_dataset("theta", data=theta)
h5.create_dataset("phi", data=phi)
h5.create_dataset("frequency", data=freqs, compression="lzf", shuffle=True)

h5["data"].dims.create_scale(h5["theta"],     "Az angle")
h5["data"].dims.create_scale(h5["phi"],       "Alt angle")
h5["data"].dims.create_scale(h5["frequency"], "Frequency MHz")

h5["data"].dims[0].attach_scale(h5["frequency"])
h5["data"].dims[1].attach_scale(h5["theta"])
h5["data"].dims[2].attach_scale(h5["phi"])

h5.close()    
    