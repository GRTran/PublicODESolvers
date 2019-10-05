import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('output/bench/benchtest.dat')
legend = ['fwdEul','bkwdEul','RK4','ShampGord','Analytical']

for i in range(1,len(data[1,:])):
    plt.plot(data[:,0], data[:,i], label=legend[i-1], marker='.')
plt.xlabel('Time')
plt.ylabel('X')
plt.legend()
plt.title('Numerical and Analytical solution of dx/dt=ax')
plt.savefig('output/bench/python_output2.png')
#plt.show()
