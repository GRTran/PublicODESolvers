import numpy as np
import matplotlib.pyplot as plt

mean = 82362231336.050201
var = 3.1998599e20
g = 2.902086e-6
result = []
xp = []
for x in range(-9,-1):
    x = x/10.0
    print(x)
    result += [np.log( 1+ x*np.log( (1/(1+x)) * (var / mean + x ) )) - ( x / ((1+x)*mean))* (var/mean - 1)*np.log(1/g)]
    xp += [x]
    #print(x, result)
plt.plot(xp,result)
plt.show()
