#################################################################################
#### PLEASE READ ALL COMMENTS BELOW AND MAKE SURE YOU FOLLOW MY INSTRUCTIONS ####
#################################################################################

# This is the skeleton program 'NatAlgReal.py' around which you should build your implementation.
# Please read through this program and follow the instructions given.

# There are no input or output files, with the results printed to the standard output.

# As regards the two values to be entered below
# - make sure that the first two values appear within double quotes
# - make sure that no comments are inserted after you have entered the values.

# Ensure that your implementation works for *arbitrary* hard-coded functions of arbitrary
# dimension and arbitray min- and max-ranges!

##############################
#### ENTER YOUR USER-NAME ####
##############################

username = "frcg69"

###############################################################
#### ENTER THE CODE FOR THE ALGORITHM YOU ARE IMPLEMENTING ####
###############################################################

alg_code = "AB"

################################################################
#### DO NOT TOUCH ANYTHING BELOW UNTIL I TELL YOU TO DO SO! ####
####      THIS INCLUDES IMPORTING ADDITIONAL MODULES!       ####
################################################################

import time
import random
import math
import sys
import os
import datetime

def compute_f(point):
    f = 40 + (point[0]**2 - 10*math.cos(2*math.pi*point[0])) + \
        (point[1]**2 - 10*math.cos(2*math.pi*point[1])) + (point[2]**2 - 10*math.cos(2*math.pi*point[2])) + \
        (point[3]**2 - 10*math.cos(2*math.pi*point[3]))
    return f

n = 4

min_range = [-5.12, -5.12, -5.12, -5.12]
max_range = [5.12, 5.12, 5.12, 5.12]

start_time = time.time()

#########################################################################################
#### YOU SHOULDN'T HAVE TOUCHED *ANYTHING* UP UNTIL NOW APART FROM SUPPLYING VALUES  ####
####                 FOR 'username' and 'alg_code' AS REQUESTED ABOVE.               ####
####                        NOW READ THE FOLLOWING CAREFULLY!                        ####
#########################################################################################

# The function 'f' is 'n'-dimensional and you are attempting to MINIMIZE it.
# To compute the value of 'f' at some point 'point', where 'point' is a list of 'n' integers or floats,
# call the function 'compute_f(point)'.
# The ranges for the values of the components of 'point' are given above. The lists 'min_range' and
# 'max_range' above hold the minimum and maximum values for each component and you should use these
# list variables in your code.

# On termination your algorithm should be such that:
#   - the reserved variable 'min_f' holds the minimum value that you have computed for the
#     function 'f' 
#   - the reserved variable 'minimum' is a list of 'n' entries (integer or float) holding the point at which
#     your value of 'min_f' is attained.

# Note that the variables 'username', 'alg_code', 'f', 'point', 'min_f', 'n', 'min_range', 'max_range' and
# 'minimum' are all reserved.

# FOR THE RESERVED VARIABLES BELOW, YOU MUST ENSURE THAT ON TERMINATION THE TYPE
# OF THE RESPECTIVE VARIABLE IS AS SHOWN.

#  - 'min_f'                int or float
#  - 'minimum'              list of int or float

# You should ensure that your code works on any function hard-coded as above, using the
# same reserved variables and possibly of a dimension different to that given above. I will
# run your code with a different such function/dimension to check that this is the case.

# The various algorithms all have additional parameters (see the lectures). These parameters
# are detailed below and are referred to using the following reserved variables.
#
# AB (Artificial Bee Colony)
#   - 'n' = dimension of the optimization problem       int
#   - 'num_cyc' = number of cycles to iterate           int
#   - 'N' = number of employed bees / food sources      int
#   - 'M' = number of onlooker bees                     int
#   - 'lambbda' = limit threshold                       float or int
#
# FF (Firefly)
#   - 'n' = dimension of the optimization problem       int
#   - 'num_cyc' = number of cycles to iterate           int
#   - 'N' = number of fireflies                         int
#   - 'lambbda' = light absorption coefficient          float or int
#   - 'alpha' = scaling parameter                       float or int
#
# CS (Cuckoo Search)
#   - 'n' = dimension of optimization problem           int
#   - 'num_cyc' = number of cycles to iterate           int
#   - 'N' = number of nests                             int
#   - 'p' = fraction of local flights to undertake      float or int
#   - 'q' = fraction of nests to abandon                float or int
#   - 'alpha' = scaling factor for Levy flights         float or int
#   - 'beta' = parameter for Mantegna's algorithm       float or int
#
# WO (Whale Optimization)
#   - 'n' = dimension of optimization problem           int
#   - 'num_cyc' = number of cycles to iterate           int
#   - 'N' = number of whales                            int
#   - 'b' = spiral constant                             float or int
#
# BA (Bat)
#   - 'n' = dimension of optimization problem           int
#   - 'num_cyc' = number of cycles to iterate           int
#   - 'N' = number of bats                              int
#   - 'sigma' = scaling factor                          float or int
#   - 'f_min' = minimum frequency                       float or int
#   - 'f_max' = maximum frequency                       float or int

# These are reserved variables and need to be treated as such, i.e., use these names for these
# parameters and don't re-use the names. Don't forget to ensure that on termination all the above
# variables have the stated type. In particular, if you use specific numpy types then you'll need
# to ensure that they are changed prior to termination (this is checked).

# INITIALIZE THE ACTUAL PARAMETERS YOU USE FOR YOUR ALGORITHM BELOW. ENSURE THAT YOU INITIALIZE
# *ALL* OF THE PARAMETERS REQUIRED APPROPRIATELY (SEE ABOVE) FOR YOUR CHOSEN ALGORITHM.

# In summary, before you input the bulk of your code, ensure that you:
# - import any (legal) modules you wish to use in the space provided below 
# - initialize your parameters in the space provided below
# - ensure that reserved variables have the correct type on termination.

###########################################
#### NOW YOU CAN ENTER YOUR CODE BELOW ####
###########################################
####################################################
#### FIRST IMPORT ANY MODULES IMMEDIATELY BELOW ####
####################################################

import numpy as np

##########################################################
#### NOW INITIALIZE YOUR PARAMETERS IMMEDIATELY BELOW ####
##########################################################

num_cyc = 800
N = 170
M = 170
lambbda = 5
min_f = math.inf

###########################################
#### NOW INCLUDE THE REST OF YOUR CODE ####
###########################################

basealg=False
#####################################################
### TO ACCESS BASE ALGORITHM UNCOMMENT LINE BELOW ###
#####################################################
# basealg=True


def lhs(N=N, n=n, min_range=min_range, max_range=max_range):
    # generate N intervals to generate points in
    intervals = np.linspace(0, 1, N + 1)

    # generate N points, one in each interval
    sources = np.zeros((N, n))

    for i in range(n):
        sources[:, i] = np.random.uniform(intervals[i], intervals[i + 1], N)

    # convert min_range and max_range to NumPy arrays for easier transforming
    min_range = np.array(min_range)
    max_range = np.array(max_range)

    # transform the samples to the specified range
    sources = min_range + sources * (max_range - min_range)

    sources = sources.tolist()
    fitnesses = [compute_f(x) for x in sources]
    return sources, fitnesses

def initialize(N=N, n=n, min_range=min_range,max_range=max_range):
    sources = []
    fitnesses = []
    for _ in range(N):
        x = []
        for i in range(n):
            x.append(random.uniform(min_range[i],max_range[i]))
        sources.append(x)
        fitnesses.append(compute_f(x))
    return sources, fitnesses


#near neighbor algorithm
def calcnn(xk, xj, b):
    phi = random.uniform(-1,1)
    xk[b] = xk[b] + phi * (xj[b]-xk[b])
    return xk

def roulette(sources, fitnesses):
    # Calculate inverse fitness values
    inverse_fitness_values = [1 / fitness for fitness in fitnesses]

    # Calculate probabilities
    total_inverse_fitness = sum(inverse_fitness_values)
    probabilities = [inverse_fitness / total_inverse_fitness for inverse_fitness in inverse_fitness_values]

    # Randomly select an index based on probabilities
    selected_index = random.choices(range(len(sources)), weights=probabilities)[0]

    return selected_index

#generate the N food sources and store fitnesses
if basealg:
    sources,fitnesses=initialize()
else:
    sources, fitnesses = lhs()

limits = [0] * len(sources)
# cycle counter
t = 1
#main loop
while t <= num_cyc and (time.time() - start_time) < 9.8:
    for i in range(0,N+M) :
        if i < N:
            k=i
        else:
            k = roulette(sources, fitnesses)
        b = random.randint(0,n-1)
        #possible values j can take
        possible = [x for x in range(N) if x != k]
        j = random.choice(possible)
        nn = calcnn(sources[k],sources[j],b)
        for i in range(n):
            if nn[i]<min_range[i] and basealg:
                nn[i] = min_range[i]
            elif nn[i]>max_range[i] and basealg:
                nn[i] = max_range[i]
            elif nn[i]<min_range[i]: #code to "bounce" the bee away from the walls
                #bounciness determines how much the bee bounces as a scaling factor
                bounciness = random.random()
                dist = (min_range[i] - nn[i]) * bounciness
                nn[i] = min_range[i] + dist
            elif nn[i]>max_range[i]:
                bounciness = random.random()
                dist = (nn[i]- max_range[i]) * bounciness
                nn[i] = max_range[i] - dist
        nn_f = compute_f(nn)
        if (nn_f<fitnesses[k]):
            sources[k] = nn
            fitnesses[k] = nn_f
            limits[k] = 0
        else:
            limits[k] = limits[k] + 1
    for i in range(N):
        if limits[i] > lambbda:
            x = []
            for i in range(n):
                x.append(random.uniform(min_range[i],max_range[i]))
            sources[k] = x
            fitnesses[k] = (compute_f(x))
            limits[k] = 0
    t = t+1
    if min_f > min(fitnesses):
        min_f = min(fitnesses)
        minimum = sources[fitnesses.index(min_f)]

    

#########################################################
#### YOU SHOULD HAVE NOW FINISHED ENTERING YOUR CODE ####
####     DO NOT TOUCH ANYTHING BELOW THIS COMMENT    ####
#########################################################

# At this point in the execution, you should have computed your minimum value for the function 'f' in the
# variable 'min_f' and the variable 'minimum' should hold a list containing the values of the point 'point'
# for which function 'f(point)' attains your minimum.

now_time = time.time()
elapsed_time = round(now_time - start_time, 1)

error = []

try:
    n
    try:
        y = n
    except:
        error.append("*** error: 'n' has not been initialized")
        n = -1
except:
    error.append("*** error: the variable 'n' does not exist\n")
    n = -1
try:
    num_cyc
    try:
        y = num_cyc
    except:
        error.append("*** error: 'num_cyc' has not been initialized")
        num_cyc = -1
except:
    error.append("*** error: the variable 'num_cyc' does not exist")
    num_cyc = -1

if alg_code == "AB":
    try:
        N
        try:
           y = N
        except:
            error.append("*** error: 'N' has not been initialized")
            N = -1
    except:
        error.append("*** error: the variable 'N' does not exist")
        N = -1
    try:
        M
        try:
           y = M
        except:
            error.append("*** error: 'M' has not been initialized")
            M = -1
    except:
        error.append("*** error: the variable 'M' does not exist")
        M = -1
    try:
        lambbda
        try:
           y = lambbda
        except:
            error.append("*** error: 'lambbda' has not been initialized")
            lambbda = -1
    except:
        error.append("*** error: the variable 'lambbda' does not exist")
        lambbda = -1
if alg_code == "FF":
    try:
        N
        try:
           y = N
        except:
            error.append("*** error: 'N' has not been initialized")
            N = -1
    except:
        error.append("*** error: the variable 'N' does not exist")
        N = -1
    try:
        alpha
        try:
           y = alpha
        except:
            error.append("*** error: 'alpha' has not been initialized")
            alpha = -1
    except:
        error.append("*** error: the variable 'alpha' does not exist")
        alpha = -1
    try:
        lambbda
        try:
           y = lambbda
        except:
            error.append("*** error: 'lambbda' has not been initialized")
            lambbda = -1
    except:
        error.append("*** error: the variable 'lambbda' does not exist")
        lambbda = -1
if alg_code == "CS":
    try:
        N
        try:
           y = N
        except:
            error.append("*** error: 'N' has not been initialized")
            N = -1
    except:
        error.append("*** error: the variable 'N' does not exist")
        N = -1
    try:
        p
        try:
           y = p
        except:
            error.append("*** error: 'p' has not been initialized")
            p = -1
    except:
        error.append("*** error: the variable 'p' does not exist")
        p = -1
    try:
        q
        try:
           y = q
        except:
            error.append("*** error: 'q' has not been initialized")
            q = -1
    except:
        error.append("*** error: the variable 'q' does not exist")
        q = -1
    try:
        alpha
        try:
           y = alpha
        except:
            error.append("*** error: 'alpha' has not been initialized")
            alpha = -1
    except:
        error.append("*** error: the variable 'alpha' does not exist")
        alpha = -1
    try:
        beta
        try:
           y = beta
        except:
            error.append("*** error: 'beta' has not been initialized")
            beta = -1
    except:
        error.append("*** error: the variable 'beta' does not exist")
        beta = -1
if alg_code == "WO":
    try:
        N
        try:
           y = N
        except:
            error.append("*** error: 'N' has not been initialized")
            N = -1
    except:
        error.append("*** error: the variable 'N' does not exist")
        N = -1
    try:
        b
        try:
           y = b
        except:
            error.append("*** error: 'b' has not been initialized")
            b = -1
    except:
        error.append("*** error: the variable 'b' does not exist")
        b = -1
if alg_code == "BA":
    try:
        sigma
        try:
           y = sigma
        except:
            error.append("*** error: 'sigma' has not been initialized")
            sigma = -1
    except:
        error.append("*** error: the variable 'sigma' does not exist")
        sigma = -1
    try:
        f_max
        try:
           y = f_max
        except:
            error.append("*** error: the variable 'f_max' has not been initialized")
            f_max = -1
    except:
        error.append("*** error: the variable 'f_max' does not exist")
        f_max = -1
    try:
        f_min
        try:
           y = f_min
        except:
            error.append("*** error: 'f_min' has not been initialized")
            f_min = -1
    except:
        error.append("*** error: the variable 'f_min' does not exist")
        f_min = -1

if type(n) != int:
    error.append("*** error: 'n' is not an integer: it is {0} and it has type {1}".format(n, type(n)))
if type(num_cyc) != int:
    error.append("*** error: 'num_cyc' is not an integer: it is {0} and it has type {1}".format(num_cyc, type(num_cyc)))

if alg_code == "AB":
    if type(N) != int:
        error.append("*** error: 'N' is not an integer: it is {0} and it has type {1}".format(N, type(N)))
    if type(M) != int:
        error.append("*** error: 'M' is not an integer: it is {0} and it has type {1}".format(M, type(M)))
    if type(lambbda) != int and type(lambbda) != float:
        error.append("*** error: 'lambbda' is not an integer or a float: it is {0} and it has type {1}".format(lambbda, type(lambbda)))

if alg_code == "FF":
    if type(N) != int:
        error.append("*** error: 'N' is not an integer: it is {0} and it has type {1}".format(N, type(N)))
    if type(lambbda) != int and type(lambbda) != float:
        error.append("*** error: 'lambbda' is not an integer or a float: it is {0} and it has type {1}".format(lambbda, type(lambbda)))
    if type(alpha) != int and type(alpha) != float:
        error.append("*** error: 'alpha' is not an integer or a float: it is {0} and it has type {1}".format(alpha, type(alpha)))

if alg_code == "CS":
    if type(N) != int:
        error.append("*** error: 'N' is not an integer: it is {0} and it has type {1}".format(N, type(N)))
    if type(p) != int and type(p) != float:
        error.append("*** error: 'p' is not an integer or a float: it is {0} and it has type {1}".format(p, type(p)))
    if type(q) != int and type(q) != float:
        error.append("*** error: 'q' is not an integer or a float: it is {0} and it has type {1}".format(q, type(q)))
    if type(alpha) != int and type(alpha) != float:
        error.append("*** error: 'alpha' is not an integer or a float: it is {0} and it has type {1}".format(alpha, type(alpha)))
    if type(beta) != int and type(beta) != float:
        error.append("*** error: 'beta' is not an integer or a float: it is {0} and it has type {1}".format(beta, type(beta)))

if alg_code == "WO":
    if type(N) != int:
        error.append("*** error: 'N' is not an integer: it is {0} and it has type {1}\n".format(N, type(N)))
    if type(b) != int and type(b) != float:
        error.append("*** error: 'b' is not an integer or a float: it is {0} and it has type {1}".format(b, type(b)))

if alg_code == "BA":
    if type(sigma) != int and type(sigma) != float:
        error.append("*** error: 'sigma' is not an integer or a float: it is {0} and it has type {1}".format(sigma, type(sigma)))
    if type(f_min) != int and type(f_min) != float:
        error.append("*** error: 'f_min' is not an integer or a float: it is {0} and it has type {1}".format(f_min, type(f_min)))
    if type(f_max) != int and type(f_max) != float:
        error.append("*** error: 'f_max' is not an integer or a float: it is {0} and it has type {1}".format(f_max, type(f_max)))

if type(min_f) != int and type(min_f) != float:
    error.append("*** error: there is no real-valued variable 'min_f'")
if type(minimum) != list:
    error.append("*** error: there is no tuple 'minimum' giving the minimum point")
elif type(n) == int and len(minimum) != n:
    error.append("*** error: there is no {0}-tuple 'minimum' giving the minimum point; you have a {1}-tuple".format(n, len(minimum)))
elif type(n) == int:
    for i in range(0, n):
        if not "int" in str(type(minimum[i])) and not "float" in str(type(minimum[i])):
            error.append("*** error: the value for component {0} (ranging from 1 to {1}) in the minimum point is not numeric\n".format(i + 1, n))

if error != []:
    print("\n*** ERRORS: there were errors in your execution:")
    length = len(error)
    for i in range(0, length):
        print(error[i])
    print("\n Fix these errors and run your code again.\n")
else:
    print("\nYou have found a minimum value of {0} and a minimum point of {1}.".format(min_f, minimum))
    print("Your elapsed time was {0} seconds.\n".format(elapsed_time))
    
