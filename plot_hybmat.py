#!/opt/python/2.7.15/gcc-4.8.5/bin/python

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import sys, os
import argparse

from basislabel import *
from hjl_common import *

parser = argparse.ArgumentParser(description="Plot ldos with j-basis")
parser.add_argument('dirs', metavar='D', type=str, nargs='+',
		                    help='an integer for the accumulator')
parser.add_argument("-t", "--totally", action="store_true", help="Plot total data in one panel" )
parser.add_argument("--oldsoc", action="store_true" , help="SOC strength is in a old-format and multiplied by 2./3." )
parser.add_argument("-y", "--yrange", type=str, help="set yrange e.g. _y0_y1" )
parser.add_argument("--selfyrange", "--sy", type=str, help="set yrange of self-energy e.g. _y0_y1" )
parser.add_argument("--selfxrange", "--sx", type=str, help="set xrange of self-energy e.g. _x0_x1" )
parser.add_argument("--hybyrange", "--hyy", type=str, help="set yrange of hyb-energy e.g. _y0_y1" )
parser.add_argument("--hybxrange", "--hyx", type=str, help="set xrange of hyb-energy e.g. _x0_x1" )
parser.add_argument("-x", "--xxrange", type=str, help="set xrange e.g. _x0_x1" )
parser.add_argument("-b", "--basis", type=str, help="Set basis. e.g. \"-b j | -b t2g\"" ) 
parser.add_argument("-c", "--chdat", type=str, help="Choose a component(s) of basis. e.g. \"-c _0_1_2 (default 0,2,4)\"" ) 
parser.add_argument("-u", "--ups", action="store_true" , help="Plot (pseudo)spin-up components of basis" )
parser.add_argument("-a", "-s", "--add", action="store_true", help="Add plotting of self-energy." ) 
parser.add_argument("--hyb", "--hy",     action="store_true", help="Add plotting of imag. hybridization fucntions." ) 
parser.add_argument("--hybre", "--hyre", action="store_true", help="Add plotting of real hybridization fucntions." ) 
parser.add_argument("--iadd", action="store_true", help="Add plotting of integrated PDOS" ) 
parser.add_argument("--pdf", action="store_true", help="Save the figure in .pdf" ) 
parser.add_argument("--png", action="store_true", help="Save the figure in .png" ) 
parser.add_argument("--cptmp", action="store_true", help="Copy the saved figure into Dropbox/tmp_linux/") 
parser.add_argument("--savetransparent", "--savetp", "--stp", action="store_true", help="Save the fig with transparent background" )
parser.add_argument("--notitle", action="store_true", help="Remove the title caption." )
parser.add_argument("--trans", action="store_true", help="Transpose the plot" ) 
parser.add_argument("--selfm", action="store_true", help="Plot minus self-energy" ) 
parser.add_argument("--dxtick", action="store_true", help="Double the xtick increment." ) 
parser.add_argument("--dytick", action="store_true", help="Double the ytick increment." ) 
parser.add_argument("--total", "-tt", action="store_true", help="Plot sum of densities as a total density." ) 
parser.add_argument("-v", "--vert", action="store_true", help="Plot vertically for ppt-use." ) 
parser.add_argument("--notick", action="store_true", help="Remove ticks." ) 
parser.add_argument("--yoff", action="store_true", help="Remove xticks." ) 
parser.add_argument("--fonts", type=float, help="Set fontsize." ) 
parser.add_argument("--nplot", type=int, help="Set Nplot." ) 
parser.add_argument("--fillb", action="store_true", help="Fill the plotting of blue one." ) 
parser.add_argument("--offd", action="store_true", help="Add plotting the off-diagonal." ) 
parser.add_argument("--chdatoffd", "--co", "--chdatoff", type=str, nargs="+", help="Choose a component(s) of basis. e.g. \"-c _0_1_2_4 (giving (0,1) and (2,4) components.)\"" ) 
parser.add_argument("--grid", "--gr", action='store_true', help="Draw grid-lines.")
args = parser.parse_args()

import matplotlib.pylab as pylab
rat	= 1.618
figsx	= 10
figs	= figsx 
figsrat = 0.618
figsy	= figsx * figsrat
dfigsrat= 1
if args.offd   :  figsy  += figsy
#if args.hyb   :  figsrat += dfigsrat
#if args.hybre :  figsrat += dfigsrat
if args.trans :
	figsydum	= figsx
	figsx	= figsy 
	figsy	= figsxdum
print "figsize = ({}, {})".format( figs, figsy ) 
fs=12 
ms=6
params = {      'legend.fontsize'       : fs*0.7 ,
		'figure.figsize'        : (figs,figsy) ,
		'axes.labelsize'        : fs ,
		'axes.titlesize'        : fs ,
		'xtick.labelsize'       : fs ,
		'ytick.labelsize'       : fs ,
		'lines.markersize'      : ms ,
		#'markers.fillstyle'     : 'none' ,
		'font.size'             : fs }
pylab.rcParams.update(params)


ndirs = len(args.dirs)
hobjarr = [""]*ndirs
nbasis	= 6
if args.basis : 
	if args.basis.find("j")>-1	:
		basis = "j"
		fname = "jldos"
	elif args.basis.find("t")>-1	:
		basis = "t2g"
		fname = "tldos"
	elif args.basis.find("z")>-1	:
		fname = "ldos"
	basis = args.basis
else : 
	#default
	basis = "t2g"
	fname = "tldos"

chdat = range(nbasis)
if args.ups : 
	chdat = [ 0, 2, 4 ]
	if basis.find("j")>-1 : 
		chdat = [ 0, 1, 4 ]
elif args.chdat : 
	chdat = np.array( args.chdat.split("_") , dtype=int ) 

chdatoffd = []
for i in range(nbasis) :
	for j in range(i+1,nbasis) :
		chdatoffd.append( [i,j] ) 
if args.chdatoffd :
	dumchdat = np.array( args.chdatoffd, dtype=int ) 
	chdatoffd = np.reshape( dumchdat , (-1,2) )
	print "chdatoffd : \n", chdatoffd

axarr = []
if args.offd : 
	if args.trans : 
		axr = 1 ; axc = 2
	else : 
		axr = 2 ; axc = 1
	figs, ax = plt.subplots( axr, axc )
	titleax = ax[0]
	axdiag = ax[0]
	axoffd = ax[1]
	axarr	= ax
else : 
	figs, ax = plt.subplots()
	titleax = ax
	axdiag  = ax
	axarr	= [ax]
for dd in range(ndirs) :
	hobj = headobj( args.dirs[dd] ) 
	selfeobj = datobj()
	selfeobj.getplothyb(		axdiag	, args.dirs[dd] , basis, nbasis, chdat, hobj , args ) 
	if args.offd : 
		selfeobj.getplothyboffdiag(	axoffd	, args.dirs[dd] , basis, nbasis, chdatoffd, hobj , args ) 
	#aax.axvline( 0 , color='g' , lw = 0.5 , ls='--' ) 

for axnow in axarr : 
	axnow.grid(linestyle=":")
plt.legend()

if args.notitle :  pass
else : hobj.headSetTitle( titleax ) 

argsavefigname = sys.argv[1] + "/lattice/" + hobj.fignamepart() + "_"+fname
for ar in sys.argv[2:] :
        print "arg : ", ar
        argsavefigname = argsavefigname + '_' + ar.split("-")[-1]
fsave = argsavefigname
#print fsave , "\n", argsavefigname
#plt.savefig( fsave + ".png" ) 
#print "Saved : ", fsave + ".png"
ndpi=400
transparent = False 
if args.savetransparent : transparent = True
if args.cptmp :
	dtype = "pdf"
	if args.pdf :
		dtype="pdf"
	if args.png :
		dtype="png"
	print "Saved : ", fsave + "."+dtype
	plt.savefig( fsave + "."+dtype , dpi=ndpi, transparent=transparent, format=dtype ) 
	cptmpftn( fsave, cmd="scp", destMac=False , dtype=dtype )
else : 
	print 'showing'
	plt.show()
