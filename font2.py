import sys
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.gridspec as gridspec

from matplotlib import rcParams
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['Tahoma']
#
rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'
#rcParams['font.serif'] = 'Times New Roman'
rcParams['font.size'] = 18

#print dir(rcParams)
#print rcParams

#params = {		'font.family'	: "serif" ,
#			'font.sans-serif' : ['Times New Roman'] ,
#		     	'font.size'	: 10 }
#pylab.rcParams.update(params)

#from matplotlib import rc
##rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
### for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Times New Roman']})
#rc('text', usetex=True)


fig, ax = plt.subplots()
ax.plot([1, 2, 3], label='test'+r'$ABC$')
ax.set_xlabel( 'test'+r'$ABC$' + r'$n_c \times 10^2$' )

ax.legend()
plt.show()
