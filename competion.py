import numpy as np

def competition(oldPop,A,B,C):

    if isinstance(oldPop, np.ndarray):
        totalPop = float(np.sum(oldPop))    
        a = [A*i for i in list(oldPop)]
        b = (1+B+C*totalPop)
        newPop = [i/b for i in a] 

    elif isinstance(oldPop,float):
        
        totalPop = oldPop
        a = A*oldPop
        b = (1+B+C*totalPop)
        newPop = a/b

    return newPop