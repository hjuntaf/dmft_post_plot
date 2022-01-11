#!/opt/python/2.7.15/gcc-4.8.5/bin/python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys, os

import argparse

homepath = "/home/jun/"
sys.path.append(homepath+"/tools/dmft/plot/" ) 
from basislabel import *
from hjl_common import *
from hjl_class  import *

print "ARGV : ", sys.argv
parser = argparse.ArgumentParser(description="Plot ldos with j-basis")
parser.add_argument('ddir', metavar='D', type=str, nargs='+',
		                    help='Data paths' ) 
parser.add_argument("-y", "--yrange", type=str, help="set yrange e.g. _y0_y1" )
parser.add_argument("-x", "--xxrange", type=str, help="set xrange e.g. _x0_x1" )
parser.add_argument("--yyrange", type=str, help="set yrange e.g. _y0_y1" )
parser.add_argument("-b", "--basis", type=str, help="Set basis. e.g. \"-b j | -b t2g\"" ) 
parser.add_argument("-c", "--chdat", type=str, help="Choose a component(s) of basis. e.g. \"-c _0_1_2 (default 0,2,4)\"" ) 
parser.add_argument("-u", "--ups", action="store_true" , help="Plot (pseudo)spin-up components of basis" )
parser.add_argument("--pdf", action="store_true", help="Save the figure in .pdf" ) 
parser.add_argument("--cptmp", action="store_true", help="Copy the saved figure into Dropbox/tmp_linux/") 
parser.add_argument("--notitle", action="store_true", help="Remove the title caption." )
parser.add_argument("--squre", action="store_true", help="Plot in square frame.")
parser.add_argument("-s", "--separate", action="store_true" , help="Plot data in different frame for each orbital." ) 
parser.add_argument("-w", "--w0", action="store_true" , help="Plot data near 0 ." ) 
parser.add_argument("-l", "--llog", action="store_true" , help="Plot data in log-log plot." ) 
parser.add_argument("-m", "--minus", action="store_true" , help="Plot self-energy with the original minus sign." ) 
parser.add_argument("--plusself", action="store_true" , help="Plot self-energy with the original minus sign." ) 
parser.add_argument("--count", type=int, help="Choose the count of which data you want." ) 
parser.add_argument("--ppt", action="store_true" , help="Plot as ppt-use." ) 
parser.add_argument("--paper", action="store_true" , help="Plot as paper-use." ) 
parser.add_argument("--noleg", action="store_true" , help="Remove legend")
parser.add_argument("--tl", action="store_true" , help="Tight-layout") 
parser.add_argument("--label", type=str, help="Set label" ) 
parser.add_argument("--leglabel", "--legl", type=str, help="Set the labels of the legend" ) 
parser.add_argument("--sublabel", "--slab", action="store_true", help="Set the labels of sub-figures with alphabet." ) 
parser.add_argument("-o", "--one", action="store_true" ,  help="Plot all in one frame" ) 
parser.add_argument("--figs", "--figsize", type=str, help="Set the labels of the legend" ) 
parser.add_argument("--line", "--lines", action="store_true" , help="Draw only with lines.") 
parser.add_argument("--insetlog", "--il", action="store_true" , help="Draw an inset of log-plot." ) 
parser.add_argument("--xrangeil", "--xrangeinsetlog", "--xil", type=str, help="Set xrange of the inset-log plot.") 
parser.add_argument("--yrangeil", "--yrangeinsetlog", "--yil", type=str, help="Set yrange of the inset-log plot.") 
parser.add_argument("--ladjust", "--lad", type=str, help="Set   left-adjust value." ) 
parser.add_argument("--radjust", "--rad", type=str, help="Set  right-adjust value." ) 
parser.add_argument("--badjust", "--bad", type=str, help="Set bottom-adjust value." ) 
parser.add_argument("--tadjust", "--tad", type=str, help="Set    top-adjust value." ) 
parser.add_argument("--rowsep", "--rowseparation", action="store_true" , help="Separate the subplots in row." )
parser.add_argument("--plotlast", "--pll", action="store_true", help="Plot the last part of data.")
parser.add_argument("--plotlastoff", "--plo", type=str, help="Specify offset-data in a set of last parts of data.")
parser.add_argument("--plotlabx", "--plx", type=str, help="Specify label of x-axis in plotting." ) 
parser.add_argument("--plotlaby", "--ply", type=str, help="Specify label of y-axis in plotting." ) 
parser.add_argument("--plotabs", "--plabs", action='store_true', help="Plot absolute values")
parser.add_argument("--datafile", "--df", "--datf", 		type=str, help="Specifiy a filename having lists of dataname, otherwise use './tmpplotstatic'	")
parser.add_argument("--ylabfile", "--yf", "--ylabf", "--ylf",	type=str, help="Specifiy a filename having lists of ylabel  , otherwise use './tmpplotstaticylab'	")
parser.add_argument("--nax", type=str, help="Specify the number of frames." )
parser.add_argument("--legall", "--lgl", action='store_true', help="Show all legends")
parser.add_argument("--zerohline", "--zhl", action='store_true', help="Show horizontal lines at zero.")
parser.add_argument("--coloroffset", "--co", type=str, help="Specify the offset for color-regime." ) 
parser.add_argument("--diff", "--diffenergy", action='store_true', help="Plot the difference of energy from ground state.")
parser.add_argument("--nsector", "--nsec", type=str, help="Specify number of sector. (25 w. SOC, 169 w/o SOC" ) 
parser.add_argument("--choosesector", "--csec", type=str, help="Specify a set sectors you want to plot.")
parser.add_argument("--dsector", "--dsec", type=str, help="Specify boundary of sector. (default: --dsec 11_17)" ) 
parser.add_argument("--deg", "--degenergy", action='store_true', help="Plot the energy of degenerate ground states.")
parser.add_argument("--degsector", "--degsec", "--degs", type=str, help="Specify a sector for plotting degenerate ones." ) 
parser.add_argument("--degnvector", "--degnvec", "--degn", type=str, help="Specify # of degenerate ones for plotting." ) 
parser.add_argument("--beta", type=str, help="Specify the beta. (default beta=128)" )
parser.add_argument("--backgroundcolor", "--bgc", action='store_true', help="Change background color." )
parser.add_argument("--all", "-a", action='store_true', help="Plot all the eigenvalues in scatters.")
parser.add_argument("--iimp", "--ii", type=str, help="Specify the index of the impurity.")
parser.add_argument("--noarpack", "--noarp", "--noa", action='store_true', help="Use data not including ARPACK.") 
parser.add_argument("--xtic", "--xtick", nargs='+', type=float, help="Set xtics.")
parser.add_argument("--boltzmann", "--boltz", action='store_true', help="Plot boltzmann factors from energies w.r.t. the specified beta.")
args = parser.parse_args()

import matplotlib.pylab as pylab
fs = 12 
figsy=8 ; figs = figsy * 1.2
if args.figs :
	figsarr = args.figs.split("_")
	figs    = float(figsarr[0])
	figsy   = float(figsarr[1])
print "figsize = ({}, {})".format( figs, figsy ) 
params = {      'legend.fontsize'       : fs ,
		'figure.figsize'        : (figs,figsy) ,
		'axes.labelsize'        : 'x-large' ,
		'axes.titlesize'        :'small' ,
		'xtick.labelsize'       :'x-large' ,
		'ytick.labelsize'       :'x-large' ,
		'xtick.direction'	: 'in' ,
		'ytick.direction'	: 'in' ,
		'xtick.minor.visible'	: True ,
		'ytick.minor.visible'	: True ,
		'lines.markersize'      : fs*0.7 ,
		'markers.fillstyle'     : 'full' , # full|left|right|bottom|top|none
		'font.size'             : fs }
if args.ppt :
	fs=14 
	figs = 3
	params = {      'legend.fontsize'       : fs ,
			'figure.figsize'        : (figs,figs) ,
			'axes.labelsize'        : fs ,
			'axes.titlesize'        : fs ,
			'xtick.labelsize'       : fs ,
			'ytick.labelsize'       : fs ,
			'lines.markersize'      : fs*0.4 ,
			'font.size'             : fs }
if args.paper :
	fs=14 
	figs = 3
	params = {      'figure.figsize'        : (figs,figs) ,
			'axes.labelsize'        : fs ,
			'axes.titlesize'        : fs ,
			'xtick.labelsize'       : fs ,
			'ytick.labelsize'       : fs ,
			'font.size'             : fs ,
			'legend.fontsize'       : fs*0.65 ,
			'lines.markersize'      : fs*0.4 }
#font = {'family' : 'normal',
#       'weight' : 'bold',
#       'size'   : 22}
pylab.rcParams.update(params)
print "\n"
nbasis = 6
basis  = "j"
fname = "Egptlall"
markers = markersdiff

ndir = len(args.ddir) 
if args.one :	nax  = 1
if args.nax :	nax  = int(args.nax)
else :		nax  = 1
f, ax = ( plt.subplots(1, nax, sharey=True, sharex=True) if args.rowsep is False else  plt.subplots(nax, 1, sharey=True, sharex=True) )
if nax ==1 :	axarr = [ax]*ndir
else :		axarr = ax
if ndir ==1 :	colors = returnblankftn
else : 		colors = colordiffftn

hobjarr = []
j=0
for idir in range(ndir) : 
	itemdir = args.ddir[idir]
	print "DIR : ", itemdir
	axnow = axarr[idir] 
	if itemdir.find(".p.o")>-1 :
		hobj = headobj( dirname )
		hobj.readParameters()
	else :
		hobj = headobj(itemdir) ;
		setattr( hobj , "iarg", j ) ; j+=1
		hobj.readParameters()
		hobj.readInitialMak()
		hobj.readPfile()
		beta  =  hobj.beta
	hobj.readStatic()
	if args.beta :	beta = int(args.beta)
	hobj.PDarr = []
	hobjarr.append( hobj ) 
	dname = hobj.tdir+"/next.mak"
	#if args.iimp : dname = hobj.tdir+"/iter/energyall_n{}.txt".format( int(args.iimp) )
	print "Reading :", dname
	with open(dname,'r') as fr :
		frlines	= fr.readlines()
		il	= 2
		ndeg	= int( frlines[il].split("\n")[0] ) 
		energyallArr = np.genfromtxt( frlines[3:3+ndeg] , dtype=float)[:,0]
	print "All : ", np.shape( energyallArr ) 
	print "All : "
	print energyallArr 

	energyallDiffArr = energyallArr - energyallArr[0]
	energyallDiffArr	=  np.sort( energyallDiffArr )
	print "Diff : "
	print energyallDiffArr
	ydat	= np.exp( -beta*energyallDiffArr )
	xdat	= np.arange( len(ydat) )

	m	= "o-"
	axnow.plot( xdat, ydat, m, mec='black', label="%d"%idir )

ylab = "Energy_imp"
if args.diff : ylab = r"$\Delta$" + ylab
xlab = "iteration"
if args.plotlabx : xlab = filterlabel(args.plotlabx)
axSetRange( axnow, args.xxrange , "x" )
axSetRange( axnow, args.yrange , "y" )
axnow.set_xlabel( xlab ) 
axnow.set_ylabel( ylab ) 


#plt.show()
#sys.exit(1)


if nax==1 : 
	xticax	= axnow
	titleax	= axnow
	legax	= axnow
	xlabax	= axnow
	ylabax	= axnow
	axarr	= [axnow]
else : 
	xticax	= ax[0]
	titleax	= ax[0]
	legax	= ax[-1]
	ylabax	= ax[0]
	xlabax	= ax[1]
	if args.rowsep : xlabax = ax[-1] ; ylabax = ax[nax/2] ; 

#for axnow in axarr : 
#	for xticlab in axnow.get_xticklabels() : xticlab.set_rotation(30)
#	for yticlab in axnow.get_yticklabels() : yticlab.set_rotation(30)
#print "XTIC LBA : ", xticlab
if args.xtic : 
	print "xtics : ", args.xtic
	for axnow in axarr : 
		axnow.set_xticks( args.xtic ) 
	
if args.noleg	:  
	for axnow in axarr :
		axnow.get_legend().remove()
else		: legax.legend()
if args.legall	: 
	for axnow in axarr :
		axnow.legend(loc=4)
if args.legall	: 
	for axnow in axarr :
		axnow.axhline( 0, 0,1 , color='g', ls=':', lw=0.5 )
if args.ylabfile	:
	for jj in range(nax) :
		axarr[jj].set_ylabel( '{}'.format(ylabarr[jj]) ) 
hobj = hobjarr[0]
hobj.readParameters()
if args.notitle : pass
else : 
	try :
		titleax.set_title( hobj.title + " n={:.2f}".format( getattr(hobj,"cflattu_%d"%ni) )  ) 
	except : 
		titleax.set_title( hobj.title ) 
if args.tl : plt.tight_layout()

argsavefigname = hobj.tdir+"/"+hobj.fignamepart()+"_boltzmann"
for ar in sys.argv[1:] :
        print "arg : ", ar
	if ar.find("Dir")>-1 : 
		pass
	else : 
        	argsavefigname = argsavefigname + '_' + ar.split("-")[-1]
print "length(argsavefigname)" , len(argsavefigname)
if len(argsavefigname)>350 : argsavefigname = argsavefigname[:350]
fsave = argsavefigname
ladjust=( 0.13 if (args.ladjust is None) else float(args.ladjust) )
radjust=( 0.93 if (args.radjust is None) else float(args.radjust) )
badjust=( 0.10 if (args.badjust is None) else float(args.badjust) )
tadjust=( 0.93 if (args.tadjust is None) else float(args.tadjust) )
plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, wspace=0 , hspace=0 )
pylab.rcParams.update(params)
print "axis : ", plt.axis()

if args.backgroundcolor :
	xy = plt.axis()
	x = np.linspace( xy[0],xy[1], 20 ) 
	y = np.linspace( xy[2],xy[3], 20 ) 
	print "X : ", x
	print "Y : ", y
	cend = 0.2
	ctop = [ cend if xi<0.075 else 1.-cend  for xi in x ] 
	cbot = [ 0.5 for xi in x ] 
	#ctop = np.zeros(20) + 0.5
	#cbot = np.linspace(0,1,20)
	print "X : ", ctop
	print "Y : ", cbot
	z = [[z] * 10 for z in range(10)]
	num_bars = 100  # more bars = smoother gradient
	#plt.contourf(x, y, z, num_bars)

	#plt.imshow([[64, 192], [64, 192]],
	#		cmap = plt.cm.Greens,
	#		interpolation = 'bicubic',
	#		vmin = 0, vmax = 255
	#	     )
	#axarr[0].imshow([[1.,1.], [0.,0.]], 
	#		  cmap = plt.cm.Blues,
	#		  interpolation = 'bicubic',
	#		  extent= (xy[0],0.05,xy[2],xy[3])
	#	     )
	from matplotlib.colors import NoNorm
	axarr[0].imshow([ctop, cbot], 
			  cmap = plt.cm.bwr,
			  interpolation = 'bicubic',
			  extent= (xy[0],xy[1],xy[2],xy[3]) ,
			  norm = NoNorm()
			  #vmin = 0., vmax = 0.5
		     )
	axarr[0].set_aspect('auto')
	print "DAT0 : ", hobjarr[0].PDarr[0].tag
	print "DAT1 : ", hobjarr[0].PDarr[1].tag

	print "XDAT : "
       	for a in xdat : print a,
	print ""

	print "YDAT0: "
       	print [ hobjarr[a].PDarr[0].ydat for a in range(ndat) ]
	print "YDAT0: "
       	print [ hobjarr[a].PDarr[1].ydat for a in range(ndat) ]

#
#	background_color = 'w'
#	plt.fill_between(x, y, y2=max(y), color=background_color)
#
#	plt.show()

ndpi=400
if args.pdf :
	plt.savefig( fsave + ".pdf" ) 
	print "Saved : ", fsave + ".pdf"
if args.cptmp :
	print "Saved : ", fsave + ".png"
	plt.savefig( fsave + ".png" , dpi=ndpi) 
	cptmpftn( fsave, cmd="scp", destMac=True )
	#print "Copying : ", fsave + ".png"
	#print "into : ",  "/home/jun/Dropbox/tmp_linux/" 
	#os.system( "cp " + fsave+".png " + "/home/jun/Dropbox/tmp_linux/" ) 
else : 
	plt.show()
