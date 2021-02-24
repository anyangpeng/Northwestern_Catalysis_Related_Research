'''

	This program describes a simple reaction pathway
	in a batch reactor:

	A  ==> B ==> C (0th order, 1st order or 2nd ONLY)

	The concentration of each species at any given time
	can be calculated.

'''

import math
import matplotlib.pyplot as plt
%matplotlib inline

class Reactions(object):
	"""Define a Class called reaction"""
	def __init__(self, rate, order, time):
		
		self.rate = rate
		self.order = order
		self.time = time

	def __str__(self):
		pass

	"""Method for calculating concentration"""
	def reaction(self, ini_conc):
		if self.order == 0:#Zero oder reaction
			rconc = ini_conc - self.rate * self.time
			pconc = ini_conc - rconc
		elif self.order == 1:#First order reaction
			rconc = ini_conc * (math.exp(-self.rate * self.time))
			pconc = ini_conc - rconc
		elif self.order == 2:#Second order reaction
			rconc = 1/(1/(ini_conc+1e-10) + self.rate * self.time)
			pconc = ini_conc - rconc
		else:
			print('NOT A VALID REACTION ORDER!')
		return rconc, pconc

"""Initializing Parameter"""

rate_ab=1
rate_bc=1
order_ab =1
order_bc =2
time=15
ini_conc =5


"""Define a plotting tool to picture the reaction profile"""

 def plotrxn (time, rate_ab, rate_bc, order_ab, order_bc, ini_conc):
    rconclst = []
    pconclstb = []
    pconclstc = []
    tlab=[]

    for t in range(time):
        a_b = Reactions(rate_ab,order_ab,t)
        b_c = Reactions(rate_bc,order_bc,t)

        rins,pinsb = a_b.reaction(ini_conc)
        rinsc,pinsc = b_c.reaction(pinsb)
        rconclst.append(rins)
        pconclstb.append(rinsc)
        pconclstc.append(pinsc)
        tlab.append(t)


    plt.plot(tlab,rconclst)
    plt.plot(tlab,pconclstb)
    plt.plot(tlab,pconclstc)


plotrxn (time, rate_ab, rate_bc, order_ab, order_bc, ini_conc) #Calling function



		
