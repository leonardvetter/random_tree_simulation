#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 13:58:15 2024

@author: leonardvetter

"""
import numpy as np

def calc_lucas_path(n,p_i):
    
    """
    Simulate a random walk conditioned to hit -1 for the first time at n
    with p_i as step distribution on {-1,0,1,...,k-1}

    The simulation follows the algorithm from [Dev12]

    Input: n, k, p_i

    output: x = (x_1,...,x_n) the jumps/increments
            postion = sum_i x_i
    """

    number_jumps = np.array([0]*len(p_i))
    
    np.random.seed()
    #draw the increments from multinomial untill the sum is n-1
    while np.dot(number_jumps,np.arange(0,len(p_i),1))!= n-1:
        number_jumps = np.random.multinomial(n,p_i)
  
    #now create the vector of jumps
    # Create the jump vector using NumPy for efficiency
    x = np.concatenate([[j] * count for j, count in enumerate(number_jumps)])
    
    #next permute:
    x= np.random.permutation(x)
    #shift by -1
    x = x-1 
    position = np.array([0]*(n+1))
    position[1:] = np.cumsum(x)
    
    #find minima of position and index of first min
    minima = np.min(position)
    
    v = np.min(np.argmin(position))  # Get the index of the first minimum
    #print(v)
    #cyclic shift by v
    x = np.roll(x,-v)
    position[1:] = np.cumsum(x)
    
    return position, x