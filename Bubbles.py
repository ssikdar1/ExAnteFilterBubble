#!/usr/bin/env python
# coding: utf-8

# In[7]:


import random

# number of goods
N = 1000

# number of goods to select in the simulation
# T << N
T = 100

# init the utilities
#u = [1 for i in range(N)]
u = [0 for i in range(N)]

# p ∈ (0, 1) 
p_0 = random.uniform(0,1)


# δ ∈ (0, 1) strength of learning spillovers
delta = .1

# d d(i, j) represents the distance that good j is from the consumed item i.
# so min number of hops on the circle
def d(i,j):
    return min(abs(i-j), abs(i-j-N))

# update
def update(i,j,p):
    rt = u[i]
    delta_d_ij = delta**d(i,j)    
    p[j] = (1-delta_d_ij)* p[j] + delta_d_ij * rt

def simulate(N, T, delta, p_0, u):
    
    p = [ p_0 for i in range(N)]
    
    # To keep track of what goods we chosen so far
    items = [i for i in range(N)]

    selections = []
    for t in range(T):
        #print("iterations: %d" % t)
        if t == 0:
            # for the first iteration choose item N/2 
            idx = N//2
        else:
            # after look for highest probability, skip over previously selected nodes
            val, idx = max((val, idx) for (idx, val) in enumerate(p)if idx not in selections )

        choice = items[idx]
        selections.append(choice)
        #print("chose item: %s " % (choice))

        # update beliefs
        for p_idx in range(len(p)):
            update(idx, p_idx, p)
            
    return selections


# In[8]:


selections = simulate(N, T, delta, p_0, u)
print(sorted(selections))
import json
with open('utility0.json', 'w') as outfile:
    json.dump(sorted(selections), outfile)


# In[17]:





# In[ ]:


d(20,1)
d(20,2)
d(19,20)
d(18,2)
d(8,12)

