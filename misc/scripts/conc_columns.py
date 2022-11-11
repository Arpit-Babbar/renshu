import numpy as np

data = []
for i in range(6):
    data.append(np.loadtxt('exact'+str(i+1)+'.txt'))
N = data[0].shape[0]
out = np.zeros((N,7))
out[:,0] = data[0][:,0]
for i in range(6):
    out[:,i+1] = data[i][:,1]
np.savetxt('exact.txt',out)
