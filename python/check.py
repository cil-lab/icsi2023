import numpy as np
import icsi2023
from matplotlib import pyplot as plt





def main():
    a = 50*np.ones((1,50))
    a = np.random.uniform(-10,10,(1,50))
    for i in range(10):
        print(icsi2023.eval(a,i+1))




if __name__ == '__main__':
    main()
    