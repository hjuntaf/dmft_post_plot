#!/opt/python/2.7.15/gcc-4.8.5/bin/python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys, os

import argparse

from basislabel import *
from hjl_common import *

parser = argparse.ArgumentParser(description="Plot Hyb with j-basis")
parser.add_argument('ddir', metavar='D', type=str, nargs='+',
		                    help='Data paths' ) 
parser.add_argument("-y", "--yrange", type=str, help="set yrange e.g. _y0_y1" )
parser.add_argument("-dy", "--datyrange", type=str, help="set datyrange e.g. _y0_y1" )
parser.add_argument("--hybyrange", type=str, help="set yrange of hyb e.g. _y0_y1" )
parser.add_argument("--hybxrange", type=str, help="set xrange of hyb e.g. _x0_x1" )
parser.add_argument("-x", "--xxrange", type=str, help="set xrange e.g. _x0_x1" )
parser.add_argument("-dx", "--datxrange", type=str, help="set datrange e.g. _x0_x1" )
parser.add_argument("-b", "--basis", type=str, help="Set basis. e.g. \"-b j | -b t2g\"" ) 
parser.add_argument("-c", "--chdat", type=str, help="Choose a component(s) of basis. e.g. \"-c _0_1_2 (default 0,2,4)\"" ) 
parser.add_argument("-u", "--ups", action="store_true" , help="Plot (pseudo)spin-up components of basis" )
parser.add_argument("--figs", "--figsize", type=str, help="Set the labels of the legend" ) 
parser.add_argument("--pdf", action="store_true", help="Save the figure in .pdf" ) 
parser.add_argument("--cptmp", action="store_true", help="Copy the saved figure into Dropbox/tmp_linux/") 
parser.add_argument("--notitle", action="store_true", help="Remove the title caption." )
parser.add_argument("--squre", action="store_true", help="Plot in square frame.")
parser.add_argument("-s", "--separate", action="store_true" , help="Plot data in different frame for each orbital." ) 
parser.add_argument("-w", "--w0", action="store_true" , help="Plot data near 0 ." ) 
parser.add_argument("-l", "--llog", action="store_true" , help="Plot data in log-log plot." ) 
parser.add_argument("-m", "--minus", action="store_true" , help="Plot self-energy with the original minus sign." ) 
parser.add_argument("--count", type=int, help="Choose the count of which data you want." ) 
parser.add_argument("--countarr", type=str, help="Choose the count-array of each data you want with delimiter '_'." ) 
parser.add_argument("--ppt", action="store_true" , help="Plot as ppt-use." ) 
parser.add_argument("--paper", action="store_true" , help="Plot as paper-use." ) 
parser.add_argument("--noleg", action="store_true" , help="Remove legend")
parser.add_argument("--tl", action="store_true" , help="Tight-layout") 
parser.add_argument("--label", type=str, help="Set label" ) 
parser.add_argument("--trans", action="store_true" , help="Transpose the plot") 
parser.add_argument("--dxtick", action="store_true", help="Double the xtick increment." ) 
parser.add_argument("--dytick", action="store_true", help="Double the ytick increment." ) 
parser.add_argument("--one", "-o", action="store_true", help="Plot all the data in one frame." ) 
parser.add_argument("--rew", "--wre", action="store_true", help="Plot the data of real-frequency spectrum." ) 
parser.add_argument("--fillb", action="store_true", help="Fill the plotting of blue one." )
parser.add_argument("--realpart", "--realdata", "--rp", "--rd", action="store_true", help="Plot the real-part of the data instead of the imaginary part.")
parser.add_argument("--linecolorbasis", "--lcb", "--lcbasis", action="store_true", help="Plot with different colors for each basis." )
parser.add_argument("--leglabel", "--legl", nargs="+", type=str, help="Set the labels of the legend" ) 
parser.add_argument("--dashedline", "--dashl", action="store_true", help="Plot the data as dashedlines.")
parser.add_argument("--epsilon", "--ep", type=str, help="Specify epsilon of data. (default=0.03)")
parser.add_argument("--nplot", type=str, help="Specify 'Nplot' in DOS-plotting")
parser.add_argument("--nax", type=str, help="Specify the number of ax.")
parser.add_argument("--axcomp", action="store_true",  help="Plot the components in each.")
args = parser.parse_args()

import matplotlib.pylab as pylab
fs = 12 
figsy=5 ; figs = figsy * 1.2
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
chdat = [ 0, 1, 4 ]
if args.ups : 
	chdat = [ 0, 2, 4 ]
	if basis.find("j")>-1 : 
		chdat = [ 0, 1, 4 ]
elif args.chdat : 
	chdat = np.array( args.chdat.split("_") , dtype=int ) 

fname = "self"
if args.basis : 
	if args.basis.find("j")>-1 : 
		basis = args.basis #"j"
		chdat = [ 0, 1, 4 ]
		fname = "self"
	elif args.basis.find("t")>-1 : 
		basis = args.basis #"t2g"
		chdat = [ 0, 2, 4 ]
		fname = "selft2g"
print "basis : ", basis

if args.chdat : 
	chdat = np.array( args.chdat.split("_") , dtype=int ) 

hobjarr = []
ndir = len(args.ddir) 
if args.one :	nax =  1
else : 		nax = ndir 
if args.nax :	nax =  int(args.nax)
f, ax = plt.subplots( nax, 1, sharey=True, sharex=True )

for jj in range(len(args.ddir)) :
	if nax == 1 : axnow = ax 
	elif nax > 1 : axnow = ax[jj]
	axarr = ax if args.axcomp else None
	hobj = headobj( args.ddir[jj] )
	hobj.readParameters()
	hobjarr.append( hobj ) 
	if args.nplot : 
		hobj.Nplot = int(args.nplot)
	if args.countarr :	
		if len(args.countarr.split("_"))>1 : 
			try : 
				count = int( args.countarr.split("_")[0] )
				count = int( args.countarr.split("_")[jj] )
			except  :
				count = int( args.countarr.split("_")[jj+1] )
		else : 
			count = int( args.countarr ) 
	elif args.count :	count = args.count
	else :			count = hobj.count
	if args.leglabel :
		leglabel = "" 
		print "args.leglabel : ",  args.leglabel 
		for it in args.leglabel : 
			leglabel += "${}={}$ ".format(it ,  getattr( hobj, it ) ) 
	else : 			leglabel = False
	dat = args.ddir[jj]
	print "Reading :", dat
	selfeobj = datobj()
	lc = "C{}".format( jj ) 
	if args.linecolorbasis :  lcarr = ["C%d"%a for a in [0,1,4,3,2,5] ] if basis.find('j')>-1 else  ["C%d"%a for a in [0,3,1,4,2,5] ] 
	else : lcarr = None
	if args.dashedline :
		if jj<1 : 
			dashpt = (None,None)
		else : 
			dashpt = [jj+1,1]
	else : dashpt = (None,None)
	if args.rew : 
		selfeobj.getplothyb_lc( axnow , dat , basis, nbasis, chdat, hobj , args , lc , lcarr=lcarr, laboverw=leglabel , dashes=dashpt ,realpart=args.realpart, axarr=axarr ) 
	else : 
		selfeobj.getplothybiw_lc( axnow , dat , basis, nbasis, chdat, hobj , args , lc , laboverw=leglabel , dashes=dashpt , lcarr=lcarr ,realpart=args.realpart , axarr=axarr ) 

	[x0,x1,y0,y1] = axnow.axis()
	if args.w0 : 
		axnow.set_xlim( [0,0.3] )
		if args.llog : 
			axnow.set_xlim( [x0,3e-1] )
			axnow.set_ylim( [4e-3,5e-1] )
	if args.xxrange : 
		x0 = float( args.xxrange.split("_")[-2] ) 
		x1 = float( args.xxrange.split("_")[-1] )
		axnow.set_xlim( [x0,x1] ) 
	if args.yrange : 
		y0 = float( args.yrange.split("_")[-2] ) 
		y1 = float( args.yrange.split("_")[-1] )
		axnow.set_ylim( [y0,y1] ) 
	
	if args.label : 
		axnow.annotate( args.label, (0.05,0.93), xycoords="figure fraction" ,weight='bold') 

#xtic = [pow(10,-1.*a) for a in range(3,-2,-2)]
#ax[0].set_xticks( xtic )
#for axnow in ax : 
	#for xticlab in axnow.get_xticklabels() :
		#xticlab.set_rotation(30)
#print "XTIC LBA : ", xticlab
#ax[0].minorticks_on()
titleax = axnow
hobj = hobjarr[0]
hobj.readParameters()
if args.notitle : pass
else : 
	try :
		titleax.set_title( hobj.title + " n={:.2f}".format( getattr(hobj,"cflattu_%d"%ni) )  ) 
	except : 
		titleax.set_title( hobj.title ) 
if args.noleg :  pass
else : 
	titleax.legend()
if args.tl : plt.tight_layout()

argsavefigname = hobj.tdir+"/"+hobj.fignamepart()+"_hybiw"
for ar in sys.argv[ndir+1:] :
        print "arg : ", ar
        argsavefigname = argsavefigname + '_' + ar.split("-")[-1]
fsave = argsavefigname
plt.subplots_adjust(left=0.2, bottom=0.2, wspace=0 )
pylab.rcParams.update(params)
ndpi=400
#plt.savefig( fsave + ".png" ) 
#print "Saved : ", fsave + ".png"
if args.pdf :
	plt.savefig( fsave + ".pdf" ) 
	print "Saved : ", fsave + ".pdf"
if args.cptmp :
	print "Saved : ", fsave + ".png"
	plt.savefig( fsave + ".png" , dpi=ndpi) 
	cptmpftn( fsave, cmd="scp", destMac=True )
	#plt.savefig( fsave + ".png" , dpi=ndpi) 
	#print "Copying : ", fsave + ".png"
	#print "into : ",  "/home/jun/Dropbox/tmp_linux/" 
	#os.system( "cp " + fsave+".png " + "/home/jun/Dropbox/tmp_linux/" ) 
else : 
	plt.show()
