import os
import argparse
import numpy as np
from evaluation import Evaluation
# from LoTFWA import LoTFWA
import time
import importlib
import icsi2023
import operators as opt


# Run FWA for EDA optimization

os.environ["OMP_NUM_THREADS"] ="1" # export OMP_NUM_THREADS=1
os.environ["OPENBLAS_NUM_THREADS"] ="1" # export OPENBLAS_NUM_THREADS=1
os.environ["MKL_NUM_THREADS"] ="1" # export MKL_NUM_THREADS=1
os.environ["VECLIB_MAXIMUM_THREADS"] ="1" # export VECLIB_MAXIMUM_THREADS=1
os.environ["NUMEXPR_NUM_THREADS"] ="1" # export NUMEXPR_NUM_THREADS=1

def parsing():
    parser=argparse.ArgumentParser()
    parser.add_argument("--algorithm", "-a", help="The algorithm use to optimize")
    parser.add_argument("--dim", "-d", help="The dimension of the function")
    parser.add_argument("--mode", "-m", help="The mode for testing")
    parser.add_argument("--num", "-n", help="The num for testing")
    return parser.parse_args()



def opt(args):
    opt = {}
    print('algorithm name:',args.algorithm)
    print('problem dimension:', args.dim)
    print('mode:', args.mode)
    alg=importlib.import_module('algorithms.' + args.algorithm)
    model = getattr(alg, args.algorithm)
    optimizer = model()
    start_time=time.time()
    if args.mode == 'all':
        for i in range(10):
            opt[i] = []
            for j in range(int(args.num)):
                opt[i].append(optimizer.optimize(Evaluation(i, args.dim)))
    elif int(args.mode) in [i for i in range(10)]:
        opt[int(args.mode)] = []
        opt[int(args.mode)].append(optimizer.optimize(Evaluation(int(args.mode), args.dim)))
    else:
        print('wrong mode is provided')
    end_time=time.time()

    print("final result is ",opt,"\nTime is {}".format(end_time-start_time))


    
if __name__=="__main__":
    args=parsing()
    opt(args)
