```python
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo


def kmer_probability(k): 
    
    #returns the probability of two k-mers being the same assuming an average length of 1000
    
    if k<20:
        raise ValueError("input k must be greater than 20 for probability to be less than 1")
        
    if k - int(k) != 0:
        raise ValueError("Input must be an integer (not type, just k - int_part(k) == 0)")
        
    return ((1001-k)**2)/(2**k)

def np_calculator(n, k):
    
    if k - int(k) != 0 or n - int(n) != 0:
        raise ValueError("Inputs must be integers (not type, just a - int_part(a) == 0)")
    
    return n*kmer_probability(k)

def solve_for_lambda(x, l):
    
    return np.exp(-l*x) + x - 1



def plot_kmers_against_largestCCs(kmer_range):
    
    results_1million = []
    
    results_10million = []
    
    results_100million = []
    
    results_1billion = []
    
    kmer_sizes = np.array([k for k in range(20,20+kmer_range)])
    
    nps_1million = [np_calculator(10**6, k) for k in kmer_sizes]
    nps_10million = [np_calculator(10**7, k) for k in kmer_sizes]
    nps_100million = [np_calculator(10**8, k) for k in kmer_sizes]
    nps_1billion = [np_calculator(10**9, k) for k in kmer_sizes]
    
    
    
    #loop through kmer_sizes and find largest connected component
    for i in range(kmer_range):
        
        if nps_1million[i] == 1 or nps_10million[i] == 1 or nps_100million[i] == 1 or nps_1billion[i] == 1:
            raise ValueError("There's an np == 1 value")
        
        if nps_1million[i] > 1:
            #work out how to get the right % of giant component better
#             X = np.arange(0,1.001,0.001)

#             f = np.array([1-x for x in X])

#             g = np.array([np.exp(-nps_1million[i]*x) for x in X])
            
            
#             idx = np.argwhere(np.diff(np.sign(f - g))).flatten()

#             results_1million.append(X[idx[1]]*10**6)
            
            results_1million.append(10**6 * (spo.newton(solve_for_lambda,1.5,args=[nps_1million[i]])))
            
            
        else:
            results_1million.append(np.log(10**6))
        
        if nps_10million[i] > 1:
            #work out how to get the right % of giant component better
#             X = np.arange(0,1.001,0.001)

#             f = np.array([1-x for x in X])

#             g = np.array([np.exp(-nps_10million[i]*x) for x in X])

#             idx = np.argwhere(np.diff(np.sign(f - g))).flatten()
            

#             results_10million.append(X[idx[1]]*10**7)
            results_10million.append(10**7 * (spo.newton(solve_for_lambda,1.5,args=[nps_10million[i]])))
            
        else:
            results_10million.append(np.log(10**7))
            
        if nps_100million[i] > 1:
            #work out how to get the right % of giant component better
#             X = np.arange(0,1.001,0.001)

#             f = np.array([1-x for x in X])

#             g = np.array([np.exp(-nps_100million[i]*x) for x in X])

#             idx = np.argwhere(np.diff(np.sign(f - g))).flatten()

#             results_100million.append(X[idx[1]]*10**8)
            results_100million.append(10**8 * (spo.newton(solve_for_lambda,1.5,args=[nps_100million[i]])))
            
        else:
            results_100million.append(np.log(10**8))
            
        if nps_1billion[i] > 1:
            #work out how to get the right % of giant component better
#             X = np.arange(0,1.001,0.001)

#             f = np.array([1-x for x in X])

#             g = np.array([np.exp(-nps_1billion[i]*x) for x in X])

#             idx = np.argwhere(np.diff(np.sign(f - g))).flatten()

#             results_1billion.append(X[idx[1]]*10**9)
            results_1billion.append(10**9 * (spo.newton(solve_for_lambda,1.5,args=[nps_1billion[i]])))
            
        else:
            results_1billion.append(np.log(10**9))
            
    
    #make them arrays (should do this earlier)
    results_1million,results_10million,results_100million,results_1billion\
    = np.array(results_1million),np.array(results_10million),np.array(results_100million)\
    ,np.array(results_1billion)
    
#     plt.figure()
#     plt.plot(kmer_sizes, results_1million, label='1 Million sequences')
#     plt.plot(kmer_sizes, results_10million, label='10 Million sequences')
#     plt.plot(kmer_sizes, results_100million, label='100 Million sequences')
#     plt.plot(kmer_sizes, results_1billion, label='1 Billion sequences')
#     plt.hlines(10**6, 20, 60, linestyles='dashed')
#     plt.xlabel('K-mer sizes')
#     plt.ylabel('Largest Connected Component (in millions)')
#     plt.yscale('log')
#     plt.xticks(np.arange(20,60, 5))
#     plt.legend()
    
    return kmer_sizes, results_1million, results_10million, results_100million, results_1billion

kmer_range = int(input("What range of k-mer sizes? "))

results = plot_kmers_against_largestCCs(kmer_range)

plt.figure()
plt.plot(results[0], results[1], label='1 Million sequences')
plt.plot(results[0], results[2], label='10 Million sequences')
plt.plot(results[0], results[3], label='100 Million sequences')
plt.plot(results[0], results[4], label='1 Billion sequences')
plt.hlines(10**6, 20, 20+kmer_range, linestyles='dashed')
plt.xlabel('K-mer sizes')
plt.ylabel('Largest Connected Component (in millions)')
plt.yscale('log')
plt.xticks(np.arange(20,20+kmer_range, 5))
plt.legend()
```

    What range of k-mer sizes? 35





    <matplotlib.legend.Legend at 0x7fec17502b80>




![png](K-Split%20by%20Random%20Graph_files/K-Split%20by%20Random%20Graph_0_2.png)



```python
import numpy as np
import scipy.optimize as spo
X = np.arange(0,1,0.01)

f = np.array([1-x for x in X])

g = np.array([np.exp(-3*x) for x in X])

idx = np.argwhere(np.diff(np.sign(f - g))).flatten()

print(idx, X[idx[1]])


def f(x,k):
    return np.exp(-k*x) + x - 1
    
print(spo.newton(f,1.5,args=[1.6833728296]))
```

    [ 0 94] 0.9400000000000001
    0.6836038066725082



```python
if True:
    print('here1')
    
else:
    print('here2')
    
if False:
    print('h3')
else:
    print('h4')
```

    here1
    h4



```python
np.arange(0,1.001,0.001)[-1]
```




    1.0




```python

```




    3




```python

```
