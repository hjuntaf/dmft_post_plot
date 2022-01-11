#!/opt/python/2.7.15/gcc-4.8.5/bin/python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys, os
import subprocess

from basislabel import *
from hjl_common import *

import argparse

parser = argparse.ArgumentParser(description="Read directories from \"rdminfo.dat\" and Plot it." )  
parser.add_argument('ddir', metavar='D', type=str, help='Specify a data-directory file' ) 
parser.add_argument("-d","--dat", "--data", type=str,    help="Specify name of the data (e.g. LSJd4, SLd4)." )
parser.add_argument("-u", "--use", type=str,    help="Column to be plotted. e.g. 1_2 (default) as x_y " ) 
parser.add_argument("-l", "--label", type=str,  help="Label of each data. e.g. 1's_2's_..." ) 
parser.add_argument("--xlabel", type=str,  help="Set x-label (default=U)" ) 
parser.add_argument("-j", "--jbasis", action="store_true",  help="Plot the j-basis renorm. factor." ) 
parser.add_argument("-t", "--t2gbasis", action="store_true",  help="Plot the t2g-basis renorm. factor." ) 
parser.add_argument("--cptmp", action="store_true",  help="Copy .png to Dropbox." ) 
parser.add_argument("--pdf", action="store_true",  help="Save thd file in .pdf" ) 
parser.add_argument("-x", "--xx", type=str,  help="Set xrange. e.g. _x0_x1" ) 
parser.add_argument("--xcol", type=str,  help="Choose the column of x-axis data." ) 
parser.add_argument("--notitle", action="store_true", help="Remove the title caption." )
args = parser.parse_args()

if args.xlabel : 
	xlab= args.xlabel
else :
	xlab="$i$"
lab=[]
lsall=[]
maglab=[]
mls=[]
msind=[]
if args.label : 
	try : 
		lab.append( args.label.split("_")[j] ) 
	except : 
		lab.append( "" ) 
		lab.append( "" ) 
elif args.t2gbasis : 
	lab.append( "(xy,up)" ) 
	lab.append( "(xy,dn)" ) 
	lab.append( "(yz,up)" ) 
	lab.append( "(yz,dn)" ) 
	lab.append( "(zx,up)" ) 
	lab.append( "(zx,dn)" ) 
	lsall.append( "ro-" ) 
	lsall.append( "bo--" ) 
	lsall.append( "r>-" ) 
	lsall.append( "b>--" ) 
	lsall.append( "r^-" ) 
	lsall.append( "b^--" ) 
	datname = "fillingt2g"
	magind = [[0,1],[2,3],[4,5]]
	maglab.append( "$P_{xy}$" ) 
	maglab.append( "$P_{yz}$" ) 
	maglab.append( "$P_{zx}$" ) 
	mls.append( "ko-" ) 
	mls.append( "k>--" ) 
	mls.append( "k^:" ) 
	msind.append( "o" ) 
	msind.append( ">" ) 
	msind.append( "^" ) 
	basisind = [ 0,0,1,1,2,2 ]
	colorind = [ 'r','b','r','b','r','b' ]
else : 
	lab.append( "(3/2,3/2)" ) 
	lab.append( "(3/2,1/2)" ) 
	lab.append( "(3/2,-1/2)" ) 
	lab.append( "(3/2,-3/2)" ) 
	lab.append( "(1/2,1/2)" ) 
	lab.append( "(1/2,-1/2)" ) 
	lsall.append( "rx-" ) 
	lsall.append( "rD-" ) 
	lsall.append( "bD--" ) 
	lsall.append( "bx--" ) 
	lsall.append( "rs-" ) 
	lsall.append( "bs--" ) 
	datname = "fillingdeg"
	magind = [[0,3],[1,2],[4,5]]
	maglab.append( "$P_{(3/2,3/2)}$" ) 
	maglab.append( "$P_{(3/2,1/2)}$" ) 
	maglab.append( "$P_{(1/2,1/2)}$" ) 
	mls.append( "kx-" ) 
	mls.append( "kD--" ) 
	mls.append( "ks:" ) 
	msind.append( "x" ) 
	msind.append( "D" ) 
	msind.append( "s" ) 
	basisind = [ 0,1 , 1,0, 2,2 ] 
	colorind = [ 'r','r','b','b','r','b' ]

import matplotlib.pylab as pylab
fs = 12
figs=10 ; figsy = figs * 0.5
params = {      'legend.fontsize'       : 0.7*fs ,
		'figure.figsize'        : (figs,figsy) ,
		'axes.labelsize'        : 'x-large' ,
		'axes.titlesize'        :'x-large' ,
		'xtick.labelsize'       :'x-large' ,
		'ytick.labelsize'       :'x-large' ,
		'markers.fillstyle'     : 'none' , # full|left|right|bottom|top|none
		'font.size'             : fs }
		#font = {'family' : 'normal',
			#       'weight' : 'bold',
			#       'size'   : 22}
pylab.rcParams.update(params)

hobj = headobj( args.ddir ) 
hobj.readParameters()
niter = hobj.niter - 1 

toolpath = "/home/jun/tools/dmft/plot/"
rdmName		= args.dat if args.dat else "LSJd4"
if rdmName == "LSJd4" : 
	rdminfoName	= "rdmargLSJd4"
elif rdmName == "SLd4" : 
	rdminfoName	= "rdmargSLd4"
else :
	rdminfoName	= "rdmarg"+rdmName
rdminfoNamefull	= rdminfoName + ".dat"
rdminfo = np.genfromtxt( toolpath + "/"+rdminfoNamefull, dtype=int ) 

fig, ax = plt.subplots()
width = 0.5 ; 
#ax.bar( ii , niorder , width, align='center', alpha=0.9 , color='b' )
offconf=0.3
area = 6
#print "niorder : ", niorder
dimLSJ  = 15
ndat    = dimLSJ
nargLSJ = 4
if rdmName == "LSJd4" : 
	rdmargLabel	= r"($L,S,J,J_z$|$t_{2g}^4$)"
elif rdmName == "SLd4" : 
	rdmargLabel	= r"($S,S_z,L,L_z$|$t_{2g}^4$)"
elif rdmName == "d4SL" : 
	rdmargLabel	= r"($S,S_z,L,L_z$|$t_{2g}^4$)"
else : 
	rdmargLabel	= ""
ax.set_xlabel( rdmargLabel ) 
#ax.text( 0, offconf*(1-0.05), rdmargLabel , size=area*1.4, va="center", ha="center", rotation=0 ) 
ylim	= ax.axis()[:2]
print "axis() : ", ax.axis()
for jj in range(dimLSJ) :
	st = rdminfo[jj]
	sttext = "({:2d},{:2d},{:2d},{:2d})".format( st[0], st[1], st[2], st[3] )
	ax.text( jj, offconf*(ylim[1]-ylim[0]), sttext , size=area*1.4, va="center", ha="center", rotation=90 ) 
rdmargName	= rdminfoName + "div.dat"
print "Reading : \t ", toolpath + "/"+rdmargName 
argdiv = np.genfromtxt( toolpath + "/"+rdmargName , dtype=int ) 
for aa in argdiv :
	ax.axvline(  aa-0.5,  linestyle="-" , color='lightgrey' ) 

if rdmName == "LSJd4" : 
	rdmdatNamefull	= "rdmLSJd4diag.dat"
elif rdmName == "SLd4" : 
	rdmdatNamefull	= "rdmSLd4diag.dat"
else :
	rdmdatNamefull	= "rdm{}diag.dat".format(rdmName)
datname = hobj.pathResult() + "/" + rdmdatNamefull
degen = int( hobj.degen_0 )
print "niter : ", niter
print "degen : ", degen
if degen < 2 :
	if niter < 2 : 
		print "niter<2::"
		rdmdiag = np.genfromtxt( datname, dtype=float )[1:]
	else : 
		print "niter>=2::"
		rdmdiag = np.genfromtxt( datname, dtype=float )[-1][1:]
else :
	rdmdiag = ["",""]
	rdmdiag[0] = np.genfromtxt( datname, dtype=float )[-1][1:]
	rdmdiag[1] = np.genfromtxt( datname, dtype=float )[-2][1:]
	rdmdiag	= np.genfromtxt( datname, dtype=float )[:,1:]
print "rdmdiag.sum : ",  np.sum(rdmdiag) 
print "rdmdiag ({}) : ".format( np.shape(rdmdiag) ) 

boltzdatname	= hobj.pathResult() + "/boltzweight.dat"
boltzweight	= np.genfromtxt( boltzdatname, dtype=float )[1:]
print "boltzweight ({}) : ".format( np.shape(boltzweight) ) , boltzweight 
rdmdiag = np.dot( boltzweight , rdmdiag  ) 
print "rdmdiag.sum (boltzweighted) : ",  np.sum(rdmdiag) 
print "rdmdiag ({}) (boltzweighted) : ".format( np.shape(rdmdiag) ) 

for a in rdmdiag :
	print a

jz = [0]*degen
#if degen > 1 :
#	for a in range(degen) : 
#		for i in range(ndat) : 
#			jz[a] += rdmdiag[a][i]* rdminfo[i][-1]
#			print "multiplying", rdminfo[i][-1]
#		print "jz[", degen, "]", jz[a]

rdmdiagorder    = rdmdiag 
ii = range(ndat)
if degen < 2 :
	for a in ii : 
		rdmdiagorder[a] = rdmdiag[ a ]
else : 
	rdmdiagorder    = [ np.zeros( ndat ),  np.zeros( ndat )]
	#for a in ii : 
	#	rdmdiagorder[0][a] = rdmdiag[0][ a ]
	#	rdmdiagorder[1][a] = rdmdiag[1][ a ]
	rdmdiagsum	= np.sum( rdmdiag, axis=0 )
if degen < 2 :
	ax.bar( ii , rdmdiagorder , width, align='center', alpha=0.9 , color='skyblue', edgecolor='dimgray' )
elif len(rdmdiag.shape) < 2 :
	ydat	= rdmdiag
	ax.bar( ii , ydat , width, align='center', alpha=0.9 , color='skyblue', edgecolor='dimgray' )
elif degen > 1 :
	ydat	= rdmdiagsum
	ax.bar( ii , ydat , width, align='center', alpha=0.9 , color='skyblue', edgecolor='dimgray' )
else : 
	bc = "olive"
	width = width*0.5
	iil = np.array(ii) - 0.15
	iir = np.array(ii) + 0.15
	ax.bar( iil , rdmdiagorder[0] , width, align='center',  alpha=1.0 , color=bc , edgecolor='dimgray' )
	ax.bar( iir , rdmdiagorder[1] , width, align='center',  alpha=0.5 , color=bc , edgecolor='dimgray' )
ax.set_ylabel( r'$p_i$' )
if args.notitle :  pass
else : hobj.headSetTitle( ax ) 

#ax.set_ylim( [0,0.34] ) 
plt.subplots_adjust(left=0.2, bottom=0.2 )

fname = "rdmfillLSJ" 
argsavefigname = sys.argv[1] + "/lattice/" + hobj.fignamepart() + "_"+fname
for ar in sys.argv[2:] :
        print "arg : ", ar
        argsavefigname = argsavefigname + '_' + ar.split("-")[-1]
fsave = argsavefigname
if args.cptmp :
	dtype = "pdf"
	if args.pdf :
		dtype="pdf"
	if args.png :
		dtype="png"
	print "Saved : ", fsave + "."+dtype
	ndpi=400
	plt.savefig( fsave + "."+dtype , dpi=ndpi, transparent=transparent, format=dtype ) 
	cptmpftn( fsave, cmd="scp", destMac=False , dtype=dtype )
	#print "Copying : ", fsave + ".png"
	#print "into : ",  "/home/jun/Dropbox/tmp_linux/" 
	#os.system( "cp " + fsave+".png " + "/home/jun/Dropbox/tmp_linux/" ) 
else : 
	print 'showing'
	plt.show()
