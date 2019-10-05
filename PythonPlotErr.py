import matplotlib.pyplot as plt
import numpy as np

err = np.loadtxt('output/bench/err_dat.dat')
steps = np.loadtxt('output/bench/step_dat.dat')
time = np.loadtxt('output/bench/time_dat.dat')
legend = ['fwdEul','bkwdEul','RK4','ShampGord']

plt.ylabel('Log Error')
# plt.xlabel('Log CPU Time')
# for i in range(0,len(err[1,:])):
#     plt.plot(steps[:,i], err[:,i], label=legend[i], marker='.')
for i in range(0,len(err[1,:])):
    plt.plot(time[:,i], err[:,i], label=legend[i], marker='.')
plt.xlabel('Log Time Step')
plt.legend()
plt.title('Error Analysis of log error vs log time step')
plt.savefig('output/bench/python_output2.png')
plt.show()
