#!/opt/python/2.7.15/gcc-4.8.5/bin/python
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import numpy.ma as ma
import sys, os
import argparse

from hjl_common import *
from basislabel import *

parser = argparse.ArgumentParser(description="Plot ldos with j-basis")
parser.add_argument('ddir', metavar='D', type=str, nargs='+',
		                    help='Data paths' ) 
parser.add_argument("-s", "--save", action="store_true", help="Only save the fig in .png." ) 
parser.add_argument("--saveeps", action="store_true", help="Save the fig in .eps too" )
parser.add_argument("--pdf", "--savepdf", action="store_true", help="Save the fig in .pdf too" )
parser.add_argument("--png", action="store_true", help="Save the fig in .png too" )
parser.add_argument("--savetransparent", "--savetp", "--stp", action="store_true", help="Save the fig with transparent background" )
parser.add_argument("--paper", action="store_true" , help="Plot as paper-use." ) 
parser.add_argument("--fs", "--fermisurface", action="store_true" , help="Plot the fermi-surface." ) 
parser.add_argument("--oldsoc", action="store_true" , help="SOC strength is in a old-format and multiplied by 2./3." )
parser.add_argument("-l", "--label", type=str, help="set text-label" ) 
parser.add_argument("--labelcolor", "--lc", type=str, help="specify color of the text-label" ) 
parser.add_argument("--labelposition", "--lpos", "--labpos", type=float, nargs=2, help="specify position of the text-label" ) 
parser.add_argument("-y", "--yrange", type=str, help="set yrange e.g. _y0_y1" )
parser.add_argument("-x", "--xxrange", type=str, help="set xrange e.g. _x0_x1" )
parser.add_argument("-z", "--zrange", type=str, help="set zrange e.g. _y0_y1" )
parser.add_argument("--xxrangedos", "--xdos", type=str, help="set xrange e.g. _x0_x1" )
parser.add_argument("--xdatrange", "--xdr", "--xdrange", type=float, nargs=2, help="set xrange and extract e.g. --xdr x0 x1" )
parser.add_argument("--ydatrange", "--ydr", "--ydrange", type=float, nargs=2, help="set xrange and extract e.g. --ydr x0 x1" )
parser.add_argument("-b", "--basis", type=str, help="Set basis. e.g. \"-b j | -b t2g\"" ) 
parser.add_argument("-c", "--chdat", type=str, help="Choose a component(s) of basis. e.g. \"-c _0_1_2 (default _0_2_4)\"" ) 
parser.add_argument("--dns", action="store_true", help="Take dataset as down-spin. (setting chdat=[1,3,5])" ) 
parser.add_argument("-f", "--nofermi", action="store_true", help="Remove the fermi-level line" ) 
parser.add_argument("--bandpath", type=str, help="Specify the path of band-data.") 
parser.add_argument("--band", action="store_true", help="Plot the noninteracting band-data.") 
parser.add_argument("--locmax", action="store_true", help="Plot the local maxima together." ) 
parser.add_argument("--locmaxonly", action="store_true", help="Plot the local maxima only." ) 
parser.add_argument("--exp", type=str, help="Plot the local maxima only with background of exp data specified" ) 
parser.add_argument("--th", type=str, help="Set the threshold of local maxima to plot. " ) 
parser.add_argument("--cb", "--cbar", action="store_true", help="Add the colorbar" ) 
parser.add_argument("--lineplot", "--lp", type=str, help="Plot along a specfied line. e.g. --lineplot x_0.5" ) 
parser.add_argument("--lineplotk", "--lpk", action="store_true", help="Plot along the Fermi level in kline(dispersion).")
parser.add_argument("--lineplotkw", "--lpkw", type=str, help="Set the fraction of energy value in the kline-plot, e.g. y_0.51 ." )
parser.add_argument("--lineplotlocmax", "--lplm", type=str, help="Plot along a specfied line. e.g. --lineplotlocmax x_0.5" ) 
parser.add_argument("--mlineplot", type=str, help="Plot along a specfied line. e.g. --mlineplot x_x0_x1_x2 or y_y1_y2_y3" ) 
parser.add_argument("--mlineplotlocmax", type=str, help="Plot along a specfied line with locmax. e.g. --mlineplotlocmax x_x0_x1_x2 or y_y1_y2_y3" ) 
parser.add_argument("--lineplotonly", "--lponly", "--lpo", action="store_true", help="Show the lineplot only" ) 
parser.add_argument("--lineplotxrange", "--lpx", "--lpxrange", type=str,  help="Set x-range of the lineplot" )
parser.add_argument("--lineplotyrange", "--lpy", "--lpyrange", type=str,  help="Set y-range of the lineplot" )
parser.add_argument("--lineplotfit", "--lpf", "--lpfit", action="store_true", help="Show fitting of lineplot to the Gaussian.")
parser.add_argument("--noleastsq", "--nolsq", action="store_true", help="Fitting without the least-square method algorithm.")
parser.add_argument("--leastsquare", "--lsq", type=str, help="Change the least-square weight scheme. (default=cauchy) (options=huber,soft_l1,linear)")
parser.add_argument("--readfitparam", "--rfit", "--readfp", action="store_true", help="Read the fitting paramteres from inversetau/tmptauparam.")
parser.add_argument("--savefitparam", "--sfit", "--savefp", action="store_true", help="Save the fitting paramteres in inversetau/tmptauparam.")
parser.add_argument("--grid", action="store_true", help="Take gridspec in mlineplot to plot data separately" )
parser.add_argument("--xgrid", "--xg", "--xgr", action="store_true", help="Take gridspec in mlineplot to plot data separately" )
parser.add_argument("--ygrid", "--yg", "--ygr", action="store_true", help="Take gridspec in mlineplot to plot data separately" )
parser.add_argument("--decomp", action="store_true", help="Plot the basis-decomposed one." )
parser.add_argument("--decompall", action="store_true", help="Plot the all-basis-decomposed one." )
parser.add_argument("--decomponly", action="store_true", help="Plot only impurity-decomposed one without total one." )
parser.add_argument("--decompimp", action="store_true", help="Plot the impurity-decomposed one." )
parser.add_argument("--decompsum", action="store_true", help="Plot the impurity-decomposed one." )
parser.add_argument("--decompsumonly", action="store_true", help="Plot only the impurity-decomposed one." )
parser.add_argument("--decompsumadd", type=str, nargs="+",help="Specify and add other components to the sum, and plot.")
parser.add_argument("--notitle", action="store_true", help="Take no title" ) 
parser.add_argument("--levhrat", type=float, help="Set the ratio of plotmax- and spectmax-levels." )
parser.add_argument("--levtot", type=str, nargs=3, help="Specify levtot in the false-color plot.")
parser.add_argument("--cptmp", action="store_true", help="Copy the png file to ~/Dropbox/tmp_linux/" ) 
parser.add_argument("--showaxes", action="store_true", help="Show axes" ) 
parser.add_argument("--epsilon", "--ep", type=str, help="Set Eipsilon value for different data set" )
parser.add_argument("--epsilontag", "--ept", "--eptag", action="store_true",  help="Set Eipsilon tag in file-name." )
parser.add_argument("--epsilontagfloat", "--eptf", "--eptagf", action="store_true",  help="Set Eipsilon tag with float number in file-name." )
parser.add_argument("--broadening", "--br", type=str, help="Set Broadening value for different data set" )
parser.add_argument("--broadeningtag", "--brt", "--brtag", action="store_true", help="Set Broadening tag in file-name." )
parser.add_argument("--nof", "--nofermiline", action="store_true", help="Remove the line of Fermi level." ) 
parser.add_argument("--cmjet",  action="store_true", help="Change the color-scheme into plt.cm.jet." ) 
parser.add_argument("--cm",  "--colormap", "--cmap", type=str, help="Specify the color-scheme, e.g. jet of plt.cm.jet." ) 
parser.add_argument("--cmapnote",  "--colormapnote", "--cmnote", action='store_true', help="Specify the color-scheme, e.g. jet of plt.cm.jet." ) 
parser.add_argument("--proj",  "--project", action="store_true", help="Plot the basis-projected spectrum." ) 
parser.add_argument("--fswn", type=float, help="Set w for FS" )
parser.add_argument("--fswnfloat", "--fswnf", action='store_true', help="Use float-format.")
parser.add_argument("--fsind", type=str, help="Plot a different FS along GX(GV), or with kz-off. e.g. --fsind {vX,vM,or kz0.5}" )
parser.add_argument("--tl",  action="store_true", help="Tight-layout" ) 
parser.add_argument("--count", type=int, help="Set iter")
parser.add_argument("--figs", "--figsize", type=str, help="Set the labels of the legend" ) 
parser.add_argument("--fontsize", "--fonts", type=str, help="Set the size of fonts.")
parser.add_argument("--legendfontsize", "--lfonts", type=str, help="Set the size of fonts.")
parser.add_argument("--ladjust", "--lad", type=str, help="Set   left-adjust value." ) 
parser.add_argument("--radjust", "--rad", type=str, help="Set  right-adjust value." ) 
parser.add_argument("--badjust", "--bad", type=str, help="Set bottom-adjust value." ) 
parser.add_argument("--tadjust", "--tad", type=str, help="Set    top-adjust value." ) 
parser.add_argument("--wspace" , "--wsp", type=str, help="Set wspace-adjust value." ) 
parser.add_argument("--hspace" , "--hsp", type=str, help="Set hspace-adjust value." ) 
parser.add_argument("--noleg", "--nol",  action="store_true", help="Remove the legend.") 
parser.add_argument("--leglabel", "--legl", type=str, nargs="+", help="Set the labels of the legend" ) 
parser.add_argument("--nolab", "--nolabel",  action="store_true", help="Remove the xy-label.") 
parser.add_argument("--notick", action="store_true", help="Remove ticks." ) 
parser.add_argument("--noytick", "--noytic", action="store_true", help="Remove y-ticks." ) 
parser.add_argument("--noxtick", "--noxtic", action="store_true", help="Remove x-ticks." ) 
parser.add_argument("--noylabel", "--noylab", action="store_true", help="Remove y-label." ) 
parser.add_argument("--noxlabel", "--noxlab", action="store_true", help="Remove x-label." ) 
parser.add_argument("--locmaxcolor", "--lmc", "--lmcolor", type=str, help="Set a color of the locmax-plot." ) 
parser.add_argument("--fermivelocity", "--fermiv", action="store_true",  help="Obatina Fermi velocities.")
parser.add_argument("--ldos", "--ld", "--pd", "--pdos", action="store_true",  help="Add PDOS as a subplots")
parser.add_argument("--nplot", type=int, help="Set Nplot." ) 
parser.add_argument("--nplotdos", type=int, help="Set Nplot for dos." ) 
parser.add_argument("--naxis", type=int, help="Set Naxis." ) 
parser.add_argument("--nk", "--Nk", type=int, help="Set Nk." ) 
parser.add_argument("--angle", "--ang", type=str, help="Set kM-kX angle." ) 
parser.add_argument("--trans", "--tr", action="store_true", help="Transpose the plot" ) 
parser.add_argument("--yoff", action="store_true", help="Remove xticks." ) 
parser.add_argument("--dxtick", action="store_true", help="Double the xtick increment." ) 
parser.add_argument("--dytick", action="store_true", help="Double the ytick increment." ) 
parser.add_argument("--tickpad", "--tpad", type=str, help="Set pad of tick labels.")
parser.add_argument("--fillb", action="store_true", help="Fill the plotting of blue one." ) 
parser.add_argument("--sublabel", "--slab", action="store_true", help="Set the labels of sub-figures with alphabet." ) 
parser.add_argument("--sublabelparameter", "--slabpar", "--slabp", type=str, help="Set a parameter of the sub-labels.")
parser.add_argument("--sublabeloff", "--slaboff", "--slabo", type=str, help="Set the off-set of the sub-labels.")
parser.add_argument("--fermisurfacelabel", "--fslabel", "--fslab", action="store_true", help="Set the labels of fermi-surface." ) 
parser.add_argument("--fsonly", "--fso", action="store_true", help="Plot fs only.")
parser.add_argument("--fssoc", "--fss", "--osoc", "--onesoc", "--addsoc", "--asoc", type=str, help="Plot the result of the one-shot SOC with specfited strength.")
parser.add_argument("--vertdos", "--vdos", action="store_true", help="Plot DOS in the dispersion from 'vert'.")
parser.add_argument("--vertself", "--vself", action="store_true", help="Plot self-energy in the dispersion from 'vert'.")
parser.add_argument("--vertselfrealpart", "--vselfre", "--vselfreal", action="store_true", help="Plot self-energy in the dispersion from 'vert'.")
parser.add_argument("--nimp", "--ni", "--nimpurity", type=str, help="Set the number of impurity sites.")
args = parser.parse_args()

import matplotlib.pylab as pylab
figs=8.  
figsrat = 0.8
figsratpar = 0.25
if args.decomp      : figsrat = figsratpar
elif args.decompall : figsrat = figsratpar*2
figsy = figs*figsrat
fsize=10
if args.fontsize : fsize = float( args.fontsize )
fs=fsize
legendfs = fs
if args.legendfontsize : legendfs = float( args.legendfontsize ) 
ms=figs
if args.figs :
	figsarr = args.figs.split("_")
	figs    = float(figsarr[0])
	figsy   = float(figsarr[1])
print "figsize = ({}, {})".format( figs, figsy ) 
def fixfigs() : 
	if args.figs :
		figsarr = args.figs.split("_")
		figs    = float(figsarr[0])
		figsy   = float(figsarr[1])
params = {      'legend.fontsize'       : legendfs ,
		'figure.figsize'        : (figs,figsy) ,
		'axes.labelsize'        : 'x-large' ,
		'axes.titlesize'        :'small' ,
		'xtick.labelsize'       :'x-large' ,
		'ytick.labelsize'       :'x-large' ,
		'lines.markersize'      : ms ,
		#'markers.fillstyle'     : 'none' ,
		'font.size'             : fs }
pylab.rcParams.update(params)
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

ndir = len(args.ddir)
print "ndir : ", ndir
dirarr = args.ddir
i=0 
hobjarr = []
for jj in range(ndir) :
	print "dir%d:"%i, dirarr[jj] ; i+=1
	hobj = headobj( args.ddir[jj] )
	hobj.readParameters()
	hobjarr.append( hobj ) 
	setattr( hobj, "indZ", 0 ) 
	setattr( hobj, "indZdecomp", 0 ) 
	setattr( hobj, "indZfs", 0 ) 
	setattr( hobj, "indZdecompfs", 0 ) 
	dat = args.ddir[jj]

argsavefigname = dat + "/lattice/vdx/" + hobj.fignamepart()
for ar in sys.argv[1:] :
	print "arg : ", ar 
	if ar.find("Dir")>-1 : 
		pass
	else : 
		argsavefigname = argsavefigname + '_' + ar.split("-")[-1]
levhrat = 1
if args.levhrat :
	levhrat = args.levhrat

if args.notitle : hobj.title = ""
basis = "t2g"
chdat = [0,2,4]
chdat2= [1,0,4]
if args.basis : 
	basis = args.basis
	if basis.find("t")>-1 : 
		#basis = "t2g"
		chdat = [0,2,4]
		chdat2= [1,0,4]
	if basis.find("j")>-1 : 
		#basis = "j"
		chdat = [1,0,4]
		chdat2= [0,2,4]
	if basis.find("s")>-1 : 
		#basis = "s"
		chdat = [0,2,3]
		chdat2= [1,0,4]
	if args.basis.find("z")>-1 : 
		#basis = basis +"z"
		pass
print "basis : ", basis 
if args.chdat : 
	chdat = args.chdat.split("_")
elif args.dns :
	chdat = [1,3,5]

epsilon = 0.02
if args.epsilon :
	epsilon = float(args.epsilon)
	for hobj in hobjarr : hobj.title = hobj.title + " $\eta_s$={}".format(epsilon) 
for hobj in hobjarr : setattr( hobj , "epsilon", float(epsilon) ) 
print "epsilon : ", epsilon

broadening = 0.02
if args.broadening :
	broadening = float(args.broadening)
	for hobj in hobjarr : hobj.title = hobj.title + " $\eta_d$={}".format(broadening) 
for hobj in hobjarr : setattr( hobj , "broadening", float(broadening) ) 

epsilontag	= ""
broadeningtag	= ""
if args.epsilontag :
	epsilontag = "ep{:g}_".format(epsilon)
if args.epsilontagfloat :
	epsilontag = "ep{:.2f}_".format(epsilon)
if args.broadeningtag :
	broadeningtag = "br{:g}_".format(broadening)
tagEpBr =  epsilontag + broadeningtag
nplot = 1024
if args.nplot : nplot = int(args.nplot)
nax	= 1
naxis = 4
if args.naxis : naxis = int(args.naxis)
nktot = 341
if args.nk : nktot = int(args.nk)
angle = False ;
if args.angle : angle = int(args.angle)
def fanglelab(angle) : 
	if (angle is False) :
		return ""
	else :
		return "_angle{}".format(angle)
anglelab = fanglelab(angle)
angletick = '${}^\circ$'.format(angle)
if args.fswn is not None : 
	fswn = float( args.fswn ) 
	indfswn = "w{:.2g}".format(fswn)
	if args.fswnfloat :
		indfswn = "w{:.2f}".format(fswn)
else :
	fswn = 0  
	indfswn = ""
if args.fsind : 
	if args.fsind.find("kz")>-1 :
		fsind = "kz{:g}".format( float(args.fsind[2:]) ) 
	elif args.fsind.find("vX")>-1 or args.fsind.find("vM")>-1 :
		fsind = args.fsind 
	else : 
		print "ERROR ::'fsind' is not properly given. (args.fsind : ", args.fsind, ")"
		sys.exit(1)
else : 
	fsind = ""
fsfileind = "fs"+fsind+indfswn
if args.broadeningtag :
	fsfileind = "fs"+fsind+"w{:g}".format(fswn)
dispfilehead = ""
dispfileind = ""
filehead = ""
if args.fssoc :
	filehead	= "addsoc{:g}_".format(float(args.fssoc)) + filehead
	dispfilehead	= "addsoc{:g}_".format(float(args.fssoc)) + dispfilehead

#figs = 12
#fs = 50 ; fsm = 50
cmap = args.cm if args.cm else 'afmhot' 
cmapcolor = plt.cm.afmhot  # cmapcolor = plt.cm.nipy_spectral
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
if cmap == "jet2" :
	ncc = 512
	origcolors = plt.cm.get_cmap( 'jet' , ncc )
	newcolors  = origcolors( np.linspace(0, 1, ncc) )
	#print "ORI CORLOR : ", newcolors[:ncc/3,:]
	#print "ORI CORLOR : ", newcolors[:ncc/8,:]
	#newcolors[:ncc/3,0] = np.sqrt( newcolors[:ncc/3,0] ) 
	#newcolors[:ncc/3,1] = np.sqrt( newcolors[:ncc/3,1] ) 
	#newcolors[:ncc/3,1] = np.sqrt( np.array( np.linspace(0,1,ncc/3) ) )  
	#newcolors[:ncc/6,0] = np.linspace( 0,1, ncc/6 ) 
	#newcolors[:ncc/12,0] = np.zeros( ncc/12 ) +1
	#newcolors[:ncc/12,1] = np.zeros( ncc/12 ) +1
	#newcolors[:ncc/12,2] = np.zeros( ncc/12 ) +1
	newcolors[:ncc/8,:] =  np.power( newcolors[:ncc/8,:] , 2 )
	newcolors[:ncc/10,:] =  np.power( newcolors[:ncc/10,:] , 2 )
	#print "NEW CORLOR : ", newcolors[:ncc/3,:]

	cmapcolor = ListedColormap( newcolors )
	print "Modifying colormap. Done."
elif cmap == "jet3" :
	ncc = 512
	origcolors = plt.cm.get_cmap( 'jet' , ncc )
	newcolors  = origcolors( np.linspace(0, 1, ncc) )
	#print "ORI CORLOR : ", newcolors[:ncc/3,:]
	nccb = 174
	#print "ORI CORLOR : ", newcolors[:ncc/8,:]
	#newcolors[:ncc/3,0] = np.sqrt( newcolors[:ncc/3,0] ) 
	#newcolors[:ncc/3,1] = np.sqrt( newcolors[:ncc/3,1] ) 
	#newcolors[:ncc/3,1] = np.sqrt( np.array( np.linspace(0,1,ncc/3) ) )  
	#newcolors[:ncc/6,0] = np.linspace( 0,1, ncc/6 ) 
	#newcolors[:ncc/12,0] = np.zeros( ncc/12 ) +1
	#newcolors[:ncc/12,1] = np.zeros( ncc/12 ) +1
	#newcolors[:ncc/12,2] = np.zeros( ncc/12 ) +1
	#newcolors[:ncc/8,:] =  np.power( newcolors[:ncc/8,:] , 2 )
	#newcolors[:ncc/10,:] =  np.power( newcolors[:ncc/10,:] , 2 )
	newcolors[:nccb,2] =  np.linspace(0,1,nccb) 
	#print "NEW CORLOR : ", newcolors[:ncc/3,:]

	cmapcolor = ListedColormap( newcolors )
	print "Modifying colormap. Done."
elif cmap == "GreysMid" :
	print cmap
	origcolors	= plt.cm.get_cmap( 'Greys' , 256 )
	newcolors	= origcolors( np.linspace(0.8, 1, 80) )
	newcolorsbottom	= origcolors( np.linspace(0, 0.8, 200) )
	#newcolorsbottom[:,:]	= 1
	print "ORI CORLOR bottom : ", newcolorsbottom
	newcolors	= np.vstack( (newcolorsbottom, newcolors) )
	print "NEW CORLOR : ", newcolors
	#top = cm.get_cmap('Oranges_r', 128)
	#bottom = cm.get_cmap('Blues', 128)
	#newcolors = np.vstack((top(np.linspace(0, 1, 128)), bottom(np.linspace(0, 1, 128))))
	cmapcolor = ListedColormap( newcolors )
	#os.exit(1)
elif args.cmjet : 
	cmapcolor = plt.cm.jet
elif args.cm : 
	cmapcolor = getattr( plt.cm , args.cm ) 
if args.proj : 
	cmapcolor = [ plt.cm.binary , plt.cm.Reds , plt.cm.Greens , plt.cm.Blues ] 

def gaussftn( x, ymax, x0, sigmax ) :
	return ymax * sigmax * sigmax / ( (x-x0)*(x-x0) + sigmax*sigmax )
def func( var, t, y ) :
	return gaussftn( t, var[0], var[1], var[2] ) +  gaussftn( t, var[3], var[4], var[5] ) - y 
def funcmirrorpi( var, t, y ) :
	return gaussftn( t, var[0], var[1], var[2] ) +  gaussftn( t, var[3], var[4], var[5] ) +  gaussftn( -t+200./341, var[0], var[1], var[2] ) +  gaussftn( -t+200./341, var[3], var[4], var[5] ) - y 
def gendata( t, var ) :
	return gaussftn( t, var[0], var[1], var[2] ) +  gaussftn( t, var[3], var[4], var[5] ) 
def gendatamirrorpi( t, var ) :
	return gaussftn( t, var[0], var[1], var[2] ) +  gaussftn( t, var[3], var[4], var[5] ) +  gaussftn( -t+200./341, var[0], var[1], var[2] ) +  gaussftn( -t+200./341, var[3], var[4], var[5] ) 
class kline :
	def __init__ ( self, fig, cmapcolor , ax=False , findlocmax=False, lev=300 )  :
		#count = np.genfromtxt( hobj.ddir + "/parameters_u{:.2f}.dat".format(hobj.UF) , dtype=int )
		#try :
		#	lcount = count.shape[1]
		#	print "csh : ", count.shape
		#	count  = count[-1][0]
		#except IndexError :
		#	count  = count[0]
	
		#fspect = "{}/lattice/vdx/u{:.2f}_ts0.00_d{:.3f}_ep{}_{}th.txt".format( hobj.ddir, hobj.u, hobj.D ,epsilon, hobj.count )
		#faxis  = "{}/lattice/vdx/naxis4_Nplot1024_Nk341_ep{}_{}th.dx".format( hobj.ddir, epsilon, hobj.count )
		
		
		nkgrid=[]
		print "Entering kline().."
		print "1:",   faxis
		print "2:",   fspect
		f = open(faxis, "r")
		for line in f:
			values = line.split()
			if "KNAMES" in line:
				knames = values[1:]
				print knames,
			if "KPOINTS" in line:
				nkgrid = values[1:]
				nkgrid = np.array( nkgrid , dtype=int) 
				print nkgrid,
			if "object 2" in line:
				ny = int(values[-2])
				nx = int(values[-1])
			if "origin" in line:
				ori = float(values[-1])
			if "delta" in line:
				if float(values[-1]) > 0 :
					deltay = float(values[-1])
				else:
					deltax = float(values[-2])
		print ""
		f.close()
		
		xi = []
		yi = []
		xi = np.linspace( 0, nx*deltax, nx )
		ximax = (nx-1)*deltax
		yi = np.linspace( ori, ori+ny*deltay, ny )
		ntot = nx*ny
		
		Z=[]
		f2= open(fspect, "r")
		nbasis = 6
		if fspect.find("band_")>-1 : 
			print "Read", nx, "band data"
			Z = range(nx * nbasis)
		if "u{:.2f}".format(hobj.u) in fspect :
			print "Read", ntot, "density data"
			Z = range(ntot)
		ncount = 0
		ndata  = 0
		for line in f2 :
			values = line.split()
			darr = []
			for j in range(len(values)) :
				Z[ndata] =  float(values[j])
				#darr.append(float(values[j]))
				#Z = np.append( Z, np.array(darr) )
				#Z = np.append( Z, np.array([float(values[0]), float(values[1]), float(values[2]), float(values[3]), float(values[4])]) )
			#for j in range(5) :
				#np.append( Z, float(values[j]) )
				#if( j > 5 ) :
					#print "ERROR"
				ndata += 1
			ncount += 1
			#print "\rnline={:10d}".format(ncount),
			#print "\rnline={:3.1f}%".format(float(ncount)/ntot*100.),
		
		print "nc=", ncount
		print "ndata=", ndata
		
		if "u{:.2f}".format(hobj.u) in fspect :
			X, Y = np.meshgrid(xi, yi)
			Z = np.array(Z)
			Z=Z.reshape(X.shape)
			self.Z = Z
			#if args.xdatrange : 
			#	xmin, xmax	= args.xdatrange 
			#	xind		= np.bitwise_and( xi>xmin, xi<xmax )
			#	Xind		= np.meshgrid( xind, yi ) [0]
			#if args.ydatrange : 
			#	ymin, ymax	= args.ydatrange 
			#	yind		= np.bitwise_and( yi>ymin, yi<ymax )
			#	Yind		= np.meshgrid( xi, yind ) [1]
			#if args.xdatrange or args.ydatrange : 
			#	zind		= 
	
			#Z = np.ma.array( Z, mask=Z < 1 ) 
			if ax : self.cbar = ax.contourf(X, Y, Z, lev, rstride=1, cstride=1, cmap=cmapcolor, linewidth=0, antialiased=False)
			if args.cb :
				if args.fs : pass
				#elif (args.decomp is not None )  : self.cb = plt.colorbar( self.cbar )
			#plt.axis([0,1,-3,6])
			if args.yrange :
				ymin = float( args.yrange.split("_")[-2] ) 
				ymax = float( args.yrange.split("_")[-1] ) 
				if ax : ax.axis( [0,1,ymin,ymax] )
			if ax : axSetRange( ax, args.xxrange , "x" ) 
			textc	= args.labelcolor if args.labelcolor else 'k'
			if args.label :
				if ax : ax.text( 0.8, ymax*1.02, "{}".format( args.label ), color=textc )
			if args.leglabel :	
				textc	= args.labelcolor if args.labelcolor else 'k'
				textpos	= args.labelposition if args.labelposition else [0.1,0.9] 
				leglabel	= ""
				sep		= ", " if len(args.leglabel)>1 else ""
				for leglabelpar	 in args.leglabel :
					leglabelval	= "={}".format( getattr( hobj, leglabelpar ) ) 
					if args.paper : 
						try : 
							leglabelval	= " = {:.3f}".format( float( getattr( hobj, leglabelpar )) ) 
							if leglabelpar == "Treal" :
								leglabelval	= " = {:.0f}".format( float( getattr( hobj, leglabelpar )) ) 
						except : pass
					leglabel += filterlabel(leglabelpar) + leglabelval + sep
				if ax : ax.text( textpos[0], textpos[1], leglabel, color=textc, transform=ax.transAxes )
			if args.bandpath : 
				print "Add band"
				yb = range(nx * nbasis)
				ndata  = 0
				ncount = 0
				with open(args.bandpath) as fband :
					f2 = fband.readlines()
					for line in f2 :
						values = line.split()
						darr = []
						for j in range(len(values)) :
							yb[ndata] =  float(values[j])
							ndata += 1
						ncount += 1
						print "\rnline={:10d} in band-data".format(ncount),
						#print "\rnline={:3.1f}%".format(float(ncount)/ntot*100.),
				yb = np.array(yb)
				yb = yb.reshape( (nx,nbasis) )
				#print yb
				yb = yb.transpose()
				for i in range( nbasis ) :
					if ax : ax.plot( xi, yb[i], "k:", lw = 0.7)
		if fspect.find("band_")>-1 :
			X = xi
			Z = np.array(Z)
			Z = Z.reshape((nx,nbasis))
			Z = Z.transpose() 
			self.Z = Z
			print "len x,z  =", len(X), len(Z)
			for i in range( nbasis ) :
				if ax : ax.plot(X, Z[i], "b-")
			if ax : ax.axis([0,1,-3,2])
			if args.yrange :
				ymin = float( args.yrange.split("_")[-2] ) 
				ymax = float( args.yrange.split("_")[-1] ) 
				if ax : ax.axis( [0,1,ymin,ymax] )
			if args.label :
				textc	= args.labelcolor if args.labelcolor else 'k'
				if ax : ax.text( 0.8, ymax*1.02, "{}".format( args.label ), color=textc )
			if args.leglabel :	
				textc	= args.labelcolor if args.labelcolor else 'k'
				textpos	= args.labelposition if args.labelposition else [0.1,0.9] 
				leglabel	= ""
				sep		= ", " if len(args.leglabel)>1 else ""
				for leglabelpar	 in args.leglabel :
					leglabelval	= "={}".format( getattr( hobj, leglabelpar ) ) 
					if args.paper : 
						try : 
							leglabelval	= " = {:.3f}".format( float( getattr( hobj, leglabelpar )) ) 
							if leglabelpar == "Treal" :
								leglabelval	= " = {:.0f}".format( float( getattr( hobj, leglabelpar )) ) 
						except : pass
					leglabel += filterlabel(leglabelpar) + leglabelval + sep
				if ax : ax.text( textpos[0], textpos[1], leglabel, color=textc, transform=ax.transAxes )
			#plt.plot(X, Z[0], X, Z[1], X, Z[2], X, Z[3], X, Z[4], X, Z[5])
			#plt.plot(X, Z[0], X, Z[1], 300, rstride=1, cstride=1, cmap=cmapcolor, linewidth=0, antialiased=False)
			print "Read", nx, "band data"
		#surf = ax.contourf(X, Y, Z, 300, rstride=1, cstride=1, cmap=cm.nipy_spectral, linewidth=0, antialiased=False)
		#plt.savefig("{}/lattice/vdx/spect_u{}.png".format(hobj.ddir,u) )
		#print "scalex:", scalex
		
		# Fermi level 
		if args.nofermi : None
		else :
			dashpt= [1,1]
			if ax : ax.axhline( 0, 0,1 , color='w', dashes=dashpt , lw=0.5) 
		#plt.plot( [0, 0], [1,0], "k-", lw = 1)
		
		nkgridtot = 0
		xtics = np.zeros(naxis)
		for i in range(len(nkgrid)) :
			nkgridtot += nkgrid[i]
			xtics[i+1] = nkgridtot
		print "KPOINTS, TOTAL : ", xtics, nkgridtot
		xtics = xtics / float(nkgridtot) 
		print "xtics rescaled from KPOINTS/TOTAL : ",
		print xtics
		print "KNAMES converted : ",
		for i in range( len(knames) ) :
			if( knames[i].find('G') > -1 ) :
				knames[i] = '$\Gamma$'
			elif( knames[i].find('a') > -1 ) :
				knames[i] = angletick
		print knames
		
		if ax : ax.set_xticks(xtics)
		if ax : ax.set_xticklabels(knames ) 
		if ax : ax.tick_params( axis='y')
		if ax : ax.tick_params( axis='both', direction="in"  )
		if ax : ax.set_ylabel(r"$E$ [eV]") 
		if ax : axSetRange( ax, args.xxrange , "x" ) 

		if args.noylabel : ax.set_ylabel('') ; ax.set_yticks([])
		
	
		#if args.exp : 
		#	#if args.exp.find("B" )>-1 :
		#	expfig  = "/home/jun/Dropbox/tmp_linux/Ref_Shen2007PRL.png"
		#	img = plt.imread( expfig ) 
		#	xG = float(xtics[0])
		#	xM = float(xtics[1]) #as (Pi,0)
		#	ax.imshow(img, extent=[xM*0.7,xM,-0.1,0] )
		self.X = X
		self.Y = Y
		self.Z = Z
		self.nx = nx
		self.ny = ny
		self.xtics=xtics
		self.knames=knames

		if findlocmax : 
			findLocalMaxSave( Z, "{}/lattice/vdx/".format( hobj.ddir ) , hobj , fhead='disp', XY=[X,Y] )

class klinedecomp :
	def __init__( self, ax , bcomp , lev , cmapcolor , Zdat=None , findlocmax=False ) :
		print "Entering klinedecomp().."
		print "klinedecomp ", bcomp, cmapcolor, lev[0], lev[-1]
		count = np.genfromtxt( hobj.ddir + "/parameters_u{:.2f}.dat".format(hobj.UF) , dtype=int )
		try :
			lcount = count.shape[1]
			count  = count[lcount-1][0]
		except IndexError :
			count  = count[0]
	
		#for fname in os.listdir( os.getcwd()+'/'+hobj.ddir+"/lattice/vdx/" ) :
		#	if ".txt" in fname :
		#		if "decomp" in fname  and "_fs_" not in fname : 
		#			fspectdecomp = "{}/lattice/vdx/{}".format( hobj.ddir, fname ) 
		#		elif "_fs_" not in fname :
		#			fspect = "{}/lattice/vdx/{}".format( hobj.ddir, fname ) 
		#	if ".dx" in fname :
		#		if "decomp" in fname  and "_fs_" not in fname : 
		#			faxisdecomp  = "{}/lattice/vdx/{}".format( hobj.ddir, fname ) 
		#		elif "_fs_" not in fname :
		#			faxis  = "{}/lattice/vdx/{}".format( hobj.ddir, fname ) 
		
		
		nkgrid=[]
		if args.decomp or args.decompall: 
			faxis   = faxisdecomp
			fspect  = fspectdecomp
		faxis   = faxisdecomp
		fspect  = fspectdecomp

		print "1:",   faxis
		print "2:",   fspect
		f = open(faxis, "r")
		for line in f:
			values = line.split()
			if "KNAMES" in line:
				knames = values[1:]
				print knames,
			if "KPOINTS" in line:
				nkgrid = values[1:]
				nkgrid = np.array( nkgrid , dtype=int) 
				print nkgrid,
			if "object 2" in line:
				self.ny = int(values[-2])
				self.nx = int(values[-1])
			if "Nx" in line:
				self.nx = int(values[1])
			if "Ny" in line:
				self.ny = int(values[1])
			if "origin" in line:
				ori = float(values[-1])
			if "delta" in line:
				if float(values[-1]) > 0 :
					deltay = float(values[-1])
				else:
					deltax = float(values[-2])
			if "datafile" in line:
				fspect = hobj.ddir + "/lattice/vdx/" + values[1]
			if "xrange" in line:
				xmin = float( values[1] ) 
				xmax = float( values[2] ) 
			if "yrange" in line:
				ymin = float( values[1] ) 
				ymax = float( values[2] ) 
		f.close()
		print ""
		
		xi = []
		yi = []
		try :
			xi = np.linspace( xmin, xmax , self.nx )
		except :
			xi = np.linspace( 0, (self.nx-1)*deltax , self.nx )
			ximax = (self.nx-1)*deltax
		try :
			yi = np.linspace( ymin, ymax , self.ny )
		except :
			yi = np.linspace( 0, (self.ny-1)*deltay , self.ny )
			yimax = (self.ny-1)*deltay

		if hobj.indZdecomp : 
			Z = hobj.Zdecomp 
		else : 
			Z= np.genfromtxt( fspect, dtype=float ) 
			hobj.Zdecomp  = Z 
			hobj.indZdecomp = 1

		ntot = self.nx*self.ny
		nbasis = 6
		ncount = len(Z) * len(Z[0])
		ndata  = len(Z[0])
		print "nx,ny,ntot,ncount,nZline,ndata : ", self.nx, self.ny, ntot, ncount, len(Z), ndata,
		
		f2= open(fspect, "r")
		nbasis = 6
		ni = 1
		try :	ni = int( hobj.ni )
		except : pass
		if args.nimp :
			ni = int( args.nimp )

		if "u{:.2f}".format(hobj.u) in fspect :
			print "Z[{}][{}]".format( np.shape(Z)[0], np.shape(Z)[1] ),
			#self.Z = [0.]*ntot
		ncount = 0
		ndata  = 0

		NU = ni*nbasis
		if bcomp > -1 :
			if basis.find("j")>-1 or basis.find("z")>-1 :
				self.Z = np.array( Z[bcomp] ,dtype=float )
			elif basis.find("s")>-1 :
				self.Z = np.array( Z[bcomp+NU*2] ,dtype=float )
			else : 
				self.Z = np.array( Z[bcomp+NU] ,dtype=float )
		else : 
			self.Z	= Z
		
		ax.axis( [ xi[0],xi[-1], yi[0],yi[-1] ] ) 
		if "u{:.2f}".format(hobj.u) in fspect :
			X, Y = np.meshgrid(xi, yi)
			self.X = X
			self.Y = Y
			if Zdat is not None : 
				Zdatplot = Zdat
				self.Z	= Zdatplot 
			else :
				self.Z = np.array(self.Z)
				self.Z = self.Z.reshape(X.shape)
				Zdatplot = self.Z
			self.cbar = ax.contourf(X, Y, Zdatplot, lev, rstride=1, cstride=1, cmap=cmapcolor, linewidth=0, antialiased=False)
			#if args.cb : 
			#	cb = plt.colorbar( self.cbar )

			if args.yrange :
				ymin = float( args.yrange.split("_")[-2] ) 
				ymax = float( args.yrange.split("_")[-1] ) 
				ax.axis( [0,1,ymin,ymax] )
			if args.label :
				textc	= args.labelcolor if args.labelcolor else 'k'
				ax.text( 0.8, ymax*1.02, "{}".format( args.label ), color=textc )
			if args.leglabel :	
				textc	= args.labelcolor if args.labelcolor else 'k'
				textpos	= args.labelposition if args.labelposition else [0.1,0.9] 
				leglabel	= ""
				sep		= ", " if len(args.leglabel)>1 else ""
				for leglabelpar	 in args.leglabel :
					leglabelval	= "={}".format( getattr( hobj, leglabelpar ) ) 
					if args.paper : 
						try : 
							leglabelval	= " = {:.3f}".format( float( getattr( hobj, leglabelpar )) ) 
							if leglabelpar == "Treal" :
								leglabelval	= " = {:.0f}".format( float( getattr( hobj, leglabelpar )) ) 
						except : pass
					leglabel += filterlabel(leglabelpar) + leglabelval + sep
				if ax : ax.text( textpos[0], textpos[1], leglabel, color=textc, transform=ax.transAxes )
			#ax.text( 0.8, ymax*1.02, "{}".format( bolabel(basis,bcomp) ) )
			if basis.find("t")>-1 : fsrat = 1
			else : fsrat=1
			if (Zdat is not None) and bcomp==0 :
				bcomp -= 1
			if bcomp<0 :
		       		lab =''	
			else : 
				lab = ''
				#lab = r"${}$".format( bolabel(basis,bcomp) )
			ax.text( 0.5, ymax*0.8, lab , horizontalalignment='center', va='center', color='w' , size=fsize*fsrat )  

			dashpt = [1,1]
			if args.nof : pass
			else : ax.axhline( 0, 0,1, color='w', dashes=dashpt, lw=0.5 )

			print "xtics rescaled : ",
			nkgridtot = 0
			xtics = np.zeros(naxis)
			for i in range(len(nkgrid)) :
				nkgridtot += nkgrid[i]
				xtics[i+1] = nkgridtot
			xtics = xtics / float(nkgridtot) 
			#print nkgrid
			print xtics
			print "knames converted : ",
			for i in range( len(knames) ) :
				if( knames[i].find('G') > -1 ) :
					knames[i] = '$\Gamma$'
				elif( knames[i].find('a') > -1 ) :
					knames[i] = angletick
			print knames
			
			
			ax.set_xticks(xtics)
			ax.set_xticklabels(knames)
			ax.tick_params( axis='y')
			ax.tick_params( axis='both', direction="in"  )
			#ax.set_ylabel(r"$E$ (eV)", fontsize = fs )

		if findlocmax : 
			findLocalMaxSave( Zdatplot, "{}/lattice/vdx/".format( hobj.ddir ) , hobj , fhead='disp', XY=[X,Y] )
	
class fsurface : 
	def __init__( self, ax=False ) :
		#count = np.genfromtxt( hobj.ddir + "/parameters_u{:.2f}.dat".format(hobj.UF) , dtype=int )
		#try :
		#	lcount = count.shape[1]
		#	count  = count[lcount-1][0]
		#except IndexError :
		#	count  = count[0]
		#for fname in os.listdir( os.getcwd()+'/'+hobj.ddir+"/lattice/vdx/" ) :
		#	if "ep{}".format(epsilon) in fname : 
		#		if ".txt" in fname :
		#			if "decomp" in fname : 
		#				fspectdecomp = "{}/lattice/vdx/{}".format( hobj.ddir, fname ) 
		#			else :
		#				fspect = "{}/lattice/vdx/{}".format( hobj.ddir, fname ) 
		#		elif ".dx" in fname :
		#			if "decomp" in fname : 
		#				faxisdecomp  = "{}/lattice/vdx/{}".format( hobj.ddir, fname ) 
		#			else :
		#				faxis  = "{}/lattice/vdx/{}".format( hobj.ddir, fname ) 

		
		print "Entering fsurface().."
		
		nkgrid=[]
		print "1:",   faxis
		print "2:",   fspect
		f = open(faxis, "r")
		for line in f:
			values = line.split()
			if "KNAMES" in line:
				knames = values[1:]
				print knames
			if "KC" in line:
				kc = values[1]
			if "KR" in line:
				kr = values[1]
			if "KPOINTS" in line:
				nkgrid = values[1:]
				nkgrid = np.array( nkgrid , dtype=int) 
				print nkgrid
				self.ny = int(nkgrid[-1])
				self.nx = int(nkgrid[-2])
			if "origin" in line:
				ori = float(values[-1])
			if "delta" in line:
				if float(values[-1]) > 0 :
					deltay = float(values[-1])
				else:
					deltax = float(values[-2])
		f.close()
		
		xi = []
		yi = []
		xi = np.linspace( 0, (self.nx-1)*deltax , self.nx )
		ximax = (self.nx-1)*deltax
		yi = np.linspace( 0, (self.ny-1)*deltay , self.ny )
		ntot = self.nx*self.ny
		
		self.Z=[]
		f2= open(fspect, "r")
		nbasis = 6
		if "u{:.2f}".format(hobj.u) in fspect :
			print "Read", ntot, "density data"
			self.Z = [0.]*ntot
		ncount = 0
		ndata  = 0
		for line in f2 :
			values = line.split()
			darr = []
			for j in range(len(values)) :
				self.Z[ndata] =  float(values[j])
				#darr.append(float(values[j]))
				#Z = np.append( Z, np.array(darr) )
				#Z = np.append( Z, np.array([float(values[0]), float(values[1]), float(values[2]), float(values[3]), float(values[4])]) )
			#for j in range(5) :
				#np.append( Z, float(values[j]) )
				#if( j > 5 ) :
					#print "ERROR"
				ndata += 1
			ncount += 1
			print "\rnline={:10d}".format(ncount),
			#print "\rnline={:3.1f}%".format(float(ncount)/ntot*100.),
		
		print "nc=", ncount
		print "ndata=", ndata
		
		if "u{:.2f}".format(hobj.u) in fspect :
			X, Y = np.meshgrid(xi, yi)
			self.X = X
			self.Y = Y
			self.Z = np.array(self.Z) * nbasis*0.5
			self.Z =self.Z.reshape(X.shape)
			levh   = self.Z.max()
			lev = np.linspace( 0, levh*levhrat, 300 ) 
			if ax : self.cbar = ax.contourf(X, Y, self.Z, lev, rstride=1, cstride=1, cmap=cmapcolor, linewidth=0, antialiased=False)
			if args.cb and (args.decompall is None) and (args.decomp is None) :
				plt.colorbar( self.cbar )

			if ax : ax.axis('scaled')
			if ax : ax.axis([0,1,0,1])
			#plt.tight_layout()
			try :
				if( kc.find('G') > -1 ) : kc = '$\Gamma$'
				if( kr.find('G') > -1 ) : kr = '$\Gamma$'
				if( knames[0].find('G') > -1 ) : knames[0] = '$\Gamma$'
				#fs = 25 
				fc = 'w'
				if args.fermisurfacelabel : 
					if ax : ax.text( 0.5,  0.5, kc		, color=fc , fontsize=fsize , zorder=12 )
					if ax : ax.text( 0.88, 0.5, kr		, color=fc , fontsize=fsize , zorder=12 )
					if args.fs : 
						if ax : ax.text( 0.88, 0.88, knames[0]	, color=fc , fontsize=fsize , zorder=12 )
					else : 
						if ax : ax.text( 1.0, 0.0,   knames[0]	, color=fc , fontsize=fsize , zorder=12 )
			except : pass
			if args.notick :  pass
			else : 
				if ax : ax.set_xticks([])
				if ax : ax.set_xticklabels([])
				if ax : ax.set_yticks([])
				if ax : ax.set_yticklabels([])
			#if args.yrange :
			#	ymin = float( args.yrange.split("_")[-2] ) 
			#	ymax = float( args.yrange.split("_")[-1] ) 
			#	ax.axis( [0,1,ymin,ymax] )
			if args.label :
				textc	= args.labelcolor if args.labelcolor else 'k'
				if ax : ax.text( 0.8, 1.02, "{}".format( args.label ), color=textc )
			if args.leglabel :	
				textc	= args.labelcolor if args.labelcolor else 'k'
				textpos	= args.labelposition if args.labelposition else [0.1,0.9] 
				leglabel	= ""
				sep		= ", " if len(args.leglabel)>1 else ""
				for leglabelpar	 in args.leglabel :
					leglabelval	= "={}".format( getattr( hobj, leglabelpar ) ) 
					if args.paper : 
						try : 
							leglabelval	= " = {:.3f}".format( float( getattr( hobj, leglabelpar )) ) 
							if leglabelpar == "Treal" :
								leglabelval	= " = {:.0f}".format( float( getattr( hobj, leglabelpar )) ) 
						except : pass
					leglabel += filterlabel(leglabelpar) + leglabelval + sep
				if ax : ax.text( textpos[0], textpos[1], leglabel, color=textc, transform=ax.transAxes )
			if fswn != 0. :
				if ax : ax.text( 0, 1.02, r"$\omega=$"+"{}".format( fswn ) )
		#Ztr   = self.Z.transpose()
		#dztr = Ztr - self.Z
		#print "NNZ Z-Z.tr :"
		#for ii in range(self.nx) :
		#	for jj in range(self.nx) :
		#		if dztr[ii][jj] > 1e-6 :
		#			print "({}\t{})".format(ii,jj) 
		findLocalMaxSave( self.Z, "{}/lattice/vdx/".format( hobj.ddir ) , hobj  )
	
class fsurfacedecomp :
	def __init__( self, ax , bcomp , lev , Zdat=None , findlocmax=False) :
		print "Entering fsurfacedecomp().."
		count = np.genfromtxt( hobj.ddir + "/parameters_u{:.2f}.dat".format(hobj.UF) , dtype=int )
		try :
			lcount = count.shape[1]
			count  = count[lcount-1][0]
		except IndexError :
			count  = count[0]
	
		#for fname in os.listdir( os.getcwd()+'/'+hobj.ddir+"/lattice/vdx/" ) :
		#	if ".txt" in fname :
		#		if "decomp" in fname  and "_fs_" in fname : 
		#			fspectdecomp = "{}/lattice/vdx/{}".format( hobj.ddir, fname ) 
		#		elif "_fs_" in fname :
		#			fspect = "{}/lattice/vdx/{}".format( hobj.ddir, fname ) 
		#	if ".dx" in fname :
		#		if "decomp" in fname  and "_fs_" in fname : 
		#			faxisdecomp  = "{}/lattice/vdx/{}".format( hobj.ddir, fname ) 
		#		elif "_fs_" in fname :
		#			faxis  = "{}/lattice/vdx/{}".format( hobj.ddir, fname ) 
		
		
		nkgrid=[]
		if args.decomp or args.decompall: 
			faxis   = faxisdecomp
			fspect  = fspectdecomp
		print ""
		print "1:",   faxis
		print "2:",   fspect
		f = open(faxis, "r")
		for line in f:
			values = line.split()
			if "KNAMES" in line:
				knames = values[1:]
				print knames
			if "KC" in line:
				kc = values[1]
			if "KR" in line:
				kr = values[1]
			if "KPOINTS" in line:
				nkgrid = values[1:]
				nkgrid = np.array( nkgrid , dtype=int) 
				print nkgrid
				self.ny = int(nkgrid[-1])
				self.nx = int(nkgrid[-2])
			if "xrange" in line:
				xxrange = [ float(values[1]), float(values[2]) ] 
			if "yrange" in line:
				yrange = [ float(values[1]), float(values[2]) ]
			if "origin" in line:
				ori = float(values[-1])
			if "Nx" in line:
				nx = float(values[-1])
			if "Ny" in line:
				ny = float(values[-1])
			if "delta" in line:
				if float(values[-1]) > 0 :
					deltay = float(values[-1])
				else:
					deltax = float(values[-2])
		f.close()
		
		xi = []
		yi = []
		try :
			xi = np.linspace( 0, (self.nx-1)*deltax , self.nx )
			ximax = (self.nx-1)*deltax
		except :
			xi = np.linspace( xxrange[0], xxrange[1] , self.nx )
		try :
			yi = np.linspace( 0, (self.ny-1)*deltay , self.ny )
			yimax = (self.ny-1)*deltay
		except :
			yi = np.linspace( yrange[0], yrange[1] , self.ny )
		ntot = self.nx*self.ny
		
		self.Z=[]
		f2= open(fspect, "r")
		nbasis = 6
		ni = 1
		try :	ni = int( hobj.ni )
		except : pass
		if args.nimp :
			ni = int( args.nimp )
		NU = ni*nbasis

		if "u{:.2f}".format(hobj.u) in fspect :
			print "Read", ntot, "density data in", bcomp, "component"
			self.Z = [0.]*ntot
		ncount = 0
		ndata  = 0

		print "hobj.indZdecompfs0 : ", hobj.indZdecompfs 
		if hobj.indZdecompfs : 
			Z = hobj.Zdecompfs 
		else : 
			Z= np.genfromtxt( fspect, dtype=float ) 
			hobj.Zdecompfs  = Z 
			hobj.indZdecompfs = 1
			print "Reading done."
		print "hobj.indZdecompfs : ", hobj.indZdecompfs 
		print "Z0 : ", Z

		if bcomp > -1 :
			if basis.find("j")>-1 or basis.find("z")>-1 :
				self.Z = np.array( Z[bcomp] ,dtype=float )
			elif basis.find("s")>-1 :
				self.Z = np.array( Z[bcomp+NU*2] ,dtype=float )
			else : 
				self.Z = np.array( Z[bcomp+NU] ,dtype=float )
		#for line in f2 :
		#	values = line.split()
		#	darr = []
		#	for j in range(len(values)) :
		#		self.Z[ndata] =  float(values[j])
		#		#darr.append(float(values[j]))
		#		#Z = np.append( Z, np.array(darr) )
		#		#Z = np.append( Z, np.array([float(values[0]), float(values[1]), float(values[2]), float(values[3]), float(values[4])]) )
		#	#for j in range(5) :
		#		#np.append( Z, float(values[j]) )
		#		#if( j > 5 ) :
		#			#print "ERROR"
		#		ndata += 1
		#	ncount += 1
		#	print "\rnline={:10d}".format(ncount),
		#	#print "\rnline={:3.1f}%".format(float(ncount)/ntot*100.),
		
		print "nc=", ncount
		print "ndata=", ndata
		
		if "u{:.2f}".format(hobj.u) in fspect :
			X, Y = np.meshgrid(xi, yi)
			self.X = X
			self.Y = Y
			self.Z = np.array(self.Z)
			print "Z : ", self.Z
			self.Z=self.Z.reshape(X.shape)
			if Zdat is not None : 
				Zdatplot = Zdat
			else :
				Zdatplot = self.Z
			self.cbar = ax.contourf(X, Y, Zdatplot, lev, rstride=1, cstride=1, cmap=cmapcolor, linewidth=0, antialiased=False)
			#if args.cb : 
				#cb = plt.colorbar( self.cbar )
			ax.axis('scaled')
			ax.axis([0,1,0,1])
			#plt.tight_layout()
			try :
				if( kc.find('G') > -1 ) : kc = '$\Gamma$'
				if( kr.find('G') > -1 ) : kr = '$\Gamma$'
				if( knames[0].find('G') > -1 ) : knames[0] = '$\Gamma$'
				fs = fsize*1.8 ; fc = 'g'
				if args.fermisurfacelabel : 
					ax.text( 0.5, 0.5, kc		, fontsize=fs , color=fc )
					ax.text( 1.0, 0.5, kr		, fontsize=fs , color=fc )
					ax.text( 1.0, 1.0, knames[0]	, fontsize=fs , color=fc )
			except : pass
			ax.set_xticks([])
			ax.set_xticklabels([])
			ax.set_yticks([])
			ax.set_yticklabels([])
			[xmin,xmax,ymin,ymax] = ax.axis()
			#if args.yrange :
			#	ymin = float( args.yrange.split("_")[-2] ) 
			#	ymax = float( args.yrange.split("_")[-1] ) 
			#	ax.axis( [0,1,ymin,ymax] )
			if args.label :
				textc	= args.labelcolor if args.labelcolor else 'k'
				ax.text( 1, ymax*1.02, "{}".format( args.label ), color=textc )
			if args.leglabel :	
				textc	= args.labelcolor if args.labelcolor else 'k'
				textpos	= args.labelposition if args.labelposition else [0.1,0.9] 
				leglabel	= ""
				sep		= ", " if len(args.leglabel)>1 else ""
				for leglabelpar	 in args.leglabel :
					leglabelval	= "={}".format( getattr( hobj, leglabelpar ) ) 
					if args.paper : 
						try : 
							leglabelval	= " = {:.3f}".format( float( getattr( hobj, leglabelpar )) ) 
							if leglabelpar == "Treal" :
								leglabelval	= " = {:.0f}".format( float( getattr( hobj, leglabelpar )) ) 
						except : pass
					leglabel += filterlabel(leglabelpar) + leglabelval + sep
				if ax : ax.text( textpos[0], textpos[1], leglabel, color=textc, transform=ax.transAxes )
			if basis.find("t")>-1 : fsrat = 1
			else : fsrat=1
			if (Zdat is not None) and bcomp==0 :
				bcomp -= 1
			if bcomp<0 :
		       		lab =''	
			else : 
				lab = r"${}$".format( bolabel(basis,bcomp) )
			ax.text( 0.5, ymax*0.5, lab  ,
					horizontalalignment='center', va='center' , color='w', size=fsize*fsrat ) 
		#Ztr   = self.Z.transpose()
		#dztr = Ztr - self.Z
		#print "NNZ Z-Z.tr :"
		#for ii in range(self.nx) :
			#for jj in range(self.nx) :
				#if dztr[ii][jj] > 1e-6 :
					#print "({}\t{})".format(ii,jj) 
		if findlocmax : 
			findLocalMaxSave( Zdatplot, "{}/lattice/vdx/".format( hobj.ddir ) , hobj )
	
	
class bandonly : 
	def __init__ ( self, fig, cmapcolor , ax=False )  :
		print "Entering bandonly().."
		nkgrid=[]
		print "1:",   bandfaxis
		print "2:",   bandfspect
		f = open(bandfaxis, "r")
		for line in f:
			values = line.split()
			if "KNAMES" in line:
				knames = values[1:]
				print knames,
			if "KPOINTS" in line:
				nkgrid = values[1:]
				nkgrid = np.array( nkgrid , dtype=int) 
				print nkgrid,
			if "object 2" in line:
				ny = int(values[-2])
				nx = int(values[-1])
			if "origin" in line:
				ori = float(values[-1])
			if "delta" in line:
				if float(values[-1]) > 0 :
					deltay = float(values[-1])
				else:
					deltax = float(values[-2])
		print ""
		f.close()
		
		xi = []
		yi = []
		xi = np.linspace( 0, nx*deltax, nx )
		ximax = (nx-1)*deltax
		yi = np.linspace( ori, ori+ny*deltay, ny )
		ntot = nx*ny
		
		f2= open(bandfspect, "r")
		nbasis = 6
		Z = range(nx * nbasis)
		ncount = 0
		ndata  = 0
		print "LEN : ", len(Z), ", ", nx
		for line in f2 :
			values = line.split()
			darr = []
			for j in range(len(values)) :
				Z[ndata] =  float(values[j])
				ndata += 1
			ncount += 1
			print "\rnline={:10d}".format(ncount),
		print "nc=", ncount
		print "ndata=", ndata
		print "Add band"
		#Z = np.genfromtxt("bandfspect") 
		X = xi
		Z = np.array(Z)
		Z = Z.reshape((nx,nbasis))
		Z = Z.transpose() 
		print "len x,z  =", len(X), len(Z)
		print "Read", nx, "band data"

		for i in range( nbasis ) :
			if ax : ax.plot(X, Z[i], "b-")
		if ax : ax.axis([0,1,-3,2])
		ax.axhline( 0, 0,1 , color='k', dashes=[1,1] , lw=0.5) 

		if args.label :
			textc	= args.labelcolor if args.labelcolor else 'k'
			if ax : ax.text( 0.8, ymax*1.02, "{}".format( args.label ), color=textc )
		if args.leglabel :	
			textc	= args.labelcolor if args.labelcolor else 'k'
			textpos	= args.labelposition if args.labelposition else [0.1,0.9] 
			leglabel	= ""
			sep		= ", " if len(args.leglabel)>1 else ""
			for leglabelpar	 in args.leglabel :
				leglabelval	= "={}".format( getattr( hobj, leglabelpar ) ) 
				if args.paper : 
					try : 
						leglabelval	= " = {:.3f}".format( float( getattr( hobj, leglabelpar )) ) 
						if leglabelpar == "Treal" :
							leglabelval	= " = {:.0f}".format( float( getattr( hobj, leglabelpar )) ) 
					except : pass
				leglabel += filterlabel(leglabelpar) + leglabelval + sep
			if ax : ax.text( textpos[0], textpos[1], leglabel, color=textc, transform=ax.transAxes )

		if args.fermivelocity : 
			for ii in range(nx) :
				if X[ii]>0.2 and X[ii]<0.27 : 
					print "Xi {} : {}\t".format( ii, X[ii] ) , 
					for aa in range( nbasis ) :
						print  "{:.6f}\t".format(Z[aa][ii]) , 
					for bb in range( nbasis/2 ) :
						aa = 2*bb + 1
						dZ3b = -Z[aa][ii+2] +4*Z[aa][ii+1] -3*Z[aa][ii]
						dX3b = X[ii+2] - X[ii]
						dk3b = ( X[ii+2] - X[ii] ) * 341./100 * np.pi
						print  "{:.6f}\t".format(dZ3b/dX3b), 
						print  "| {:.6f}\t".format(dZ3b/dk3b), 
					print ""
			blinetest = range(91-68)
			glinetest = range(91-68)
			xlinetest = np.array(range(91-68))
			for ii in xlinetest :
				blinetest[ii] = 7.08 * ( X[ii+68] - 0.2286 ) 
				glinetest[ii] = 3.5  * ( X[ii+68] - 0.2611 ) 
			ax.plot( X[68+xlinetest], blinetest, "r-")
			ax.plot( X[68+xlinetest], glinetest, "r-")
			

for dd in range(ndir) : 
	hobj = hobjarr[dd] 
	cwd = os.getcwd()+'/'+hobj.ddir+"/lattice/vdx/" 
	if args.th : 
		zthratio = float(args.th) 
	else : 
		zthratio = 0.3
	if args.locmax : 
		findlocmax = True
	else : 
		findlocmax = False

	if args.fs :
		nax	= 1
	
		if fswn == 0 :
			#try : 
			fspect = "{}_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( fsfileind, hobj.UF, hobj.D, tagEpBr, hobj.count ) 
			fspectdecomp = "decomp_{}_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( fsfileind, hobj.UF, hobj.D, tagEpBr, hobj.count ) 
			faxis  = "{}_naxis{}_Nplot1_Nk40401_{}{:d}th.dx".format( fsfileind, naxis, tagEpBr, hobj.count ) 
			faxisdecomp  = "decomp_{}_naxis{}_Nplot1_Nk40401_{}{:d}th.dx".format( fsfileind, naxis, tagEpBr, hobj.count ) 
			#except : 
			#	fspect = "fs_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( hobj.UF, hobj.D, tagEpBr, hobj.count ) 
			#	fspectdecomp = "decomp_fs_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( hobj.UF, hobj.D, tagEpBr, hobj.count ) 
			#	faxis  = "fs_naxis{}_Nplot1_Nk40401_{}{:d}th.dx".format( naxis, tagEpBr, hobj.count ) 
			#	faxisdecomp  = "decomp_fs_naxis{}_Nplot1_Nk40401_{}{:d}th.dx".format( naxis, tagEpBr, hobj.count ) 
		else : 
			fspect = "{}_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( fsfileind, hobj.UF, hobj.D, tagEpBr, hobj.count ) 
			fspectdecomp = "decomp_{}_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( fsfileind, hobj.UF, hobj.D, tagEpBr, hobj.count ) 
			faxis  = "{}_naxis{}_Nplot1_Nk40401_{}{:d}th.dx".format( fsfileind, naxis, tagEpBr, hobj.count ) 
			faxisdecomp  = "decomp_{}_naxis{}_Nplot1_Nk40401_{}{:d}th.dx".format( fsfileind, naxis, tagEpBr, hobj.count ) 
		fspect		= filehead + fspect	 
		fspectdecomp	= filehead + fspectdecomp 
		faxis		= filehead + faxis 	 
		faxisdecomp	= filehead + faxisdecomp  
		fspect = cwd + fspect
		fspectdecomp = cwd + fspectdecomp
		faxis = cwd + faxis
		faxisdecomp = cwd + faxisdecomp
	
		dispfspect	 = "{}u{:.2f}_ts0.00_d{:.3f}_{}{:d}th{}.txt".format( dispfileind, hobj.UF, hobj.D, tagEpBr, hobj.count, anglelab ) 
		dispfspectdecomp = "decomp_{}u{:.2f}_ts0.00_d{:.3f}_{}{:d}th{}.txt".format( dispfileind, hobj.UF, hobj.D, tagEpBr, hobj.count, anglelab ) 
		dispfaxis 	 = "{}naxis{}_Nplot{}_Nk{}_{}{:d}th{}.dx".format( dispfileind, naxis, nplot, nktot, tagEpBr, hobj.count, anglelab) 
		dispfaxisdecomp  = "decomp_{}naxis{}_Nplot{}_Nk{}_{}{:d}th{}.dx".format( dispfileind, naxis, nplot, nktot, tagEpBr, hobj.count, anglelab ) 
		dispfspect	 = dispfilehead + dispfspect	 
		dispfspectdecomp = dispfilehead + dispfspectdecomp 
		dispfaxis 	 = dispfilehead + dispfaxis 	 
		dispfaxisdecomp  = dispfilehead + dispfaxisdecomp  
		dispfspect = cwd + dispfspect
		dispfspectdecomp = cwd + dispfspectdecomp
		dispfaxis = cwd + dispfaxis
		dispfaxisdecomp = cwd + dispfaxisdecomp
	
		if os.system("test -e {}".format(fspect)) < 1 :  pass
		else :
			print "Error :: ", fspect, "doesn't exist." 
			sys.exit(1)
	
		if (args.locmax or args.locmaxonly) and (not args.decomp) : 
			if args.locmaxonly :
				if dd<1 : fig, ax = plt.subplots( figsize=(8,8) ) 
				if args.exp : 
					if args.exp.find("B" )>-1 :
						expfig  = "/home/jun/Dropbox/tmp_linux/Ref_Bergemann2003.png"
					elif args.exp.find("I1" )>-1 :
						expfig  = "/home/jun/Dropbox/tmp_linux/Ref_Iwasawa2010PRL_a_p.png"
					elif args.exp.find("I2" )>-1 :
						expfig  = "/home/jun/Dropbox/tmp_linux/Ref_Iwasawa2010PRL_b_s.png"
					elif args.exp.find("A" )>-1 :
						expfig  = "/home/jun/Dropbox/tmp_linux/Ref_Acharya2017SciRep.png"
					elif args.exp.find("D" )>-1 :
						expfig  = "/home/jun/Dropbox/tmp_linux/Ref_Damascelli2000PRL.png"
					elif args.exp.find("Z" )>-1 :
						expfig  = "/home/jun/Dropbox/tmp_linux/Ref_Zhang2016PRL.png"
					elif args.exp.find("S" )>-1 :
						expfig  = "/home/jun/Dropbox/tmp_linux/Ref_Shen2001PRB.png"
					else : 
						expfig  = "/home/jun/Dropbox/tmp_linux/Ref_Bergemann2003.png"
					if args.exp.find("S" )>-1 :
						img = plt.imread( expfig ) 
						ax.imshow(img, extent=[0.5,1,0.0,0.5] )
					else :
						img = plt.imread( expfig ) 
						ax.imshow(img, extent=[0,1,0,1] )
				plotLocalMax( "{}/lattice/vdx/".format( hobj.ddir ) , ax,    zthratio , hobj , args )
				if args.notick :
					ax.tick_params( axis='both',          # changes apply to the x-axis
							which='both',      # both major and minor ticks are affected
							bottom='off',      # ticks along the bottom edge are off
							left='off',        
							labelleft='off',        
							labelbottom='off')
					ax.axis( [0,1,0,1] )
			else :
				if dd<1 : fig, ax = plt.subplots(1,2, sharex=True, sharey=True ,figsize=(11,6) ) 
				fs = fsurface( ax[0] )
				plotLocalMax( "{}/lattice/vdx/".format( hobj.ddir ) , ax[1], zthratio , hobj , args )
		elif args.lineplot or args.mlineplot or args.lineplotonly : 
			fsx, fsy = (  [11,6]  if args.figs is None else [figs,figsy] ) 
			if args.lineplotonly is False : 
				if dd<1 : fig, ax = plt.subplots(1,2, figsize=(fsx,fsy) , gridspec_kw = {'width_ratios':[3, 1]}) 
				fs = fsurface( ax[0] )
			else :
				if dd<1 : fig, axone = plt.subplots(figsize=(fsx,fsy))
				fs = fsurface( )
				ax = [axone]*5
	
			ax[1].tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
			if args.lineplot :
				linearg = args.lineplot
				lineargs= linearg
				lineaxis = linearg.split("_")[0]
				n0, xdat, ydat, datsection = refineLinePlot( linearg , fs )
				if args.leglabel :	leglabel = "${}={}$".format(args.leglabel ,  getattr( hobj, args.leglabel ) ) 
				else : 			leglabel = linearg
				ax[1].plot( xdat[n0], ydat[n0] , "o-", lw = 1, label=leglabel , markerfacecolor='w')#, markeredgecolor='k')
				saveLinePlot( hobj.ddir , linearg , n0, xdat[n0], ydat[n0] , hobj )
				plotIndLine( ax[0], datsection, lineaxis, (2,1) , 'w' )  if args.lineplotonly is False else '.'
				if   lineaxis=="x" : ax[1].set_xlabel("$A(\omega=0)$") ; ax[1].tick_params( axis="y", which='both', labelleft=False  )#, labelsize = args.fsm ) 
				elif lineaxis=="y" : ax[1].set_ylabel("$A(\omega=0)$") ; ax[1].tick_params( axis="x", which='both', labelbottom=False  )#, labelsize = args.fsm ) 
			elif args.mlineplot :
				j=0
				lineargs = args.mlineplot 
				lineaxis = lineargs.split("_")[0]
				for linearg in refineLineArg( args.mlineplot ) :
					print "linearg : ", linearg 
					n0, xdat, ydat, datsection = refineLinePlot( linearg , fs )
					ax[1].plot( xdat[n0], ydat[n0], lscycler(j), lw = 1, markersize = 2.5 , label=linearg )
					saveLinePlot( hobj.ddir , linearg , n0, xdat[n0], ydat[n0] , hobj )
					j+=1
					plotIndLine( ax[0], datsection, lineaxis, (2,1) , 'w' ) 
			ax[1].legend() if args.noleg is not None else '.'
			axSetRange( ax[1], args.lineplotxrange , "x" ) 
			axSetRange( ax[1], args.lineplotyrange , "y" ) 

			if args.lineplotfit : 
				for ix in range( xdat[n0] ) :
					if xdat[ix]>100./341 : 
						x1 = ix
				print "GM-line data : (", xdat[[0,x1]],")(", ydat[[0,x1]] ,")"
				nxfit = 200
			  	xfit = np.linspace( 0,0.25, nxfit ) 
			  	yfit1= range( nxfit ) 
			  	yfit2= range( nxfit ) 
			  	yfitall= range( nxfit ) 

				param1 = [ 15.87, 0.059, 0.0162/2 ] #[ ymax, xmax, sigmax ]
				param2 = [ 13.04, 0.113, 0.0100/2 ] #[ ymax, xmax, sigmax ]
				text1 = "(1)Ymax,X0,sigmaX : {:.3f}, {:.3f}, {:.3f}".format( param1[0], param1[1], param1[2] ) 
				text2 = "(2)Ymax,X0,sigmaX : {:.3f}, {:.3f}, {:.3f}".format( param2[0], param2[1], param2[2] ) 
				for a in range(nxfit) : 
					yfit1[a] = gaussftn( xfit[a], param1[0], param1[1], param1[2] ) 
					yfit2[a] = gaussftn( xfit[a], param2[0], param2[1], param2[2] )
					yfitall[a] = yfit1[a] + yfit2[a] 
					#print "Gxy : ", xfit[a], "\t", yfit1[a],"\t", yfit2[a]
				ax[1].plot( xfit, yfit1 , "k:" , label = text1 )
				ax[1].plot( xfit, yfit2 , "k:" , label = text2 )
				ax[1].plot( xfit, yfitall , "r-" )
				ax[1].legend()
				vFermibeta  = 2.53
				vFermigamma = 1.27
				invtaugamma = ( param1[2] * 2. / 3.833 ) * vFermigamma
				invtaubeta  = ( param2[2] * 2. / 3.833 ) * vFermibeta
				print "invtau(beta,gamma)[eV] : ", invtaubeta, invtaugamma
				from scipy.optimize import least_squares
				#res_log = least_squares(func, param1, loss='cauchy', f_scale=0.1, args=( xfit, yfit1 ) 
	
		elif args.lineplotlocmax or args.mlineplotlocmax : 
			if args.grid :
				import matplotlib.gridspec as gridspec
				lax = len( args.mlineplotlocmax.split('_') ) -1
				gs = gridspec.GridSpec(2*lax, 2*lax)
				ax = [""]*(lax+1)
				if dd<1 : ax[0] = plt.subplot(gs[:, :lax])
				for xx in range(lax) :
					if dd<1 : ax[xx+1] = plt.subplot(gs[ : , lax+xx:lax+xx+1])
					print "axx : ", ax[xx+1], ";", xx
					xx+=1
				print "axx : ", ax
			else : 
				if dd<1 : fig, ax = plt.subplots(1,2, figsize=(11,6) ) 
				lax = len( args.mlineplotlocmax.split('_') ) -1
				for xx in range(lax) :
					print "xx : ", xx
					np.append( ax, ax[1] )
				print "axx : ", ax
				print "laxx : ", len(ax)
				print "saxx : ", ax.shape
			fs = fsurface( ax[0] )
			ax[0].clear()
			zmax = plotLocalMax( "{}/lattice/vdx/".format( hobj.ddir ) , ax[0],    zthratio , hobj , args )
	
			if args.lineplotlocmax :
				linearg = args.lineplotlocmax
				lineargs= linearg
				lineaxis = linearg.split("_")[0]
				n0, xdat, ydat, datsection = refineLinePlot( linearg , fs )
				ax[1].plot( xdat[n0], ydat[n0] , "r-o", lw = 1, markersize = 1.5 , label=linearg )
				saveLinePlot( hobj.ddir , linearg , n0, xdat[n0], ydat[n0] , hobj )
				plotIndLine( ax[0], datsection, lineaxis, (2,1), 'b' ) 
			elif args.mlineplotlocmax :
				j=0
				lineargs = args.mlineplotlocmax 
				lineaxis = lineargs.split("_")[0]
				for linearg in refineLineArg( args.mlineplotlocmax ) :
					print "linearg : ", linearg 
					n0, xdat, ydat, datsection = refineLinePlot( linearg , fs )
					ax[j+1].plot( xdat[n0], ydat[n0], lscycler(j), lw = 1, markersize = 2.5 , label=linearg )
					setylimInd( ax[j+1], lineaxis, zmax ) 
					saveLinePlot( hobj.ddir , linearg , n0, xdat[n0], ydat[n0] , hobj )
					j+=1
					print "datsection : ", datsection 
					plotIndLine( ax[0], datsection, lineaxis, (2,1), 'b' ) 
					#ax[0].axhline( fs.Y[n0][0] , 0,1 , color='w', dashes=(2,1) , lw=0.5) 
					ax[j].legend(fontsize=7)
	
			#if args.lineplot.split('_')[0].find('x') > -1 :
			#	lineaxis = "x_"
			#	nx0 = round( float(args.lineplot.split('_')[1]) * (fs.nx-1) )
			#	n0 = int(nx0)
			#	print "nx0 : ", nx0
			#	fs.Z = fs.Z.transpose()
			#	fs.X = fs.X.transpose()
			#	fs.Y = fs.Y.transpose()
			#	xdat = fs.Z
			#	ydat = fs.Y
			#	ax[0].axvline( fs.X[n0][0] , 0,1 , color='w', dashes=(2,1) , lw=0.5) 
			#	ax[1].set_ylim( 0,1 ) 
			#elif args.lineplot.split('_')[0].find('y') > -1 :
			#	lineaxis = "y_"
			#	ny0 = round( float(args.lineplot.split('_')[1]) * (fs.ny-1) )
			#	n0 = int(ny0)
			#	print "ny0 : ", ny0
			#	xdat = fs.X
			#	ydat = fs.Z
			#	ax[0].axhline( fs.Y[n0][0] , 0,1 , color='w', dashes=(2,1) , lw=0.5) 
			#	ax[1].set_xlim( 0,1 ) 
			#print "xd : ", fs.X[n0][0], "...",  fs.X[n0][-1], "|",
			#print "yd : ", fs.Y[n0][0], "...",  fs.Y[n0][-1], "|",
			#print "zd : ", fs.Z[n0][0], "...",  fs.Z[n0][-1], "|"
			#plotAndSave( hobj.ddir , lineaxis , n0, fig, ax[1] , xdat[n0], ydat[n0] )
			if args.notitle :  pass
			else :
				ax[1].text( 0.5, 1.15, lineaxis+'='+lineargs.split('_')[1]+'~'+lineargs.split('_')[-1] )
		elif args.decomp: 
			nchdat	= len(chdat) 
			nax	= nchdat # +1
			if args.decomponly : 
				nax = nchdat 
			if args.decompsum : 
				nax += 1 
			fs = [""]*nax
			print "nax : ", nax 
			if dd<1 : fig, ax = plt.subplots( 1, nax , sharex=True, sharey=True  )

			if args.decomponly :
				i=0
				titleax = ax[len(chdat)/2]
			else : 
				i=1
				titleax = ax[i]
			fs[0] = fsurface( ax[0] ) ;
			levtot = fs[0].cbar.levels
			if args.levtot :
				lmin, lmax, npoint = args.levtot
				levtot	= np.linspace( float(lmin), float(lmax), int(npoint), endpoint=True)
			print "levtot : ", type(levtot), np.shape(levtot), levtot
			print "Zmin Zmax [%d]: "%0, fs[0].Z.min(), fs[0].Z.max()
			ax[0].clear()
			if args.decompsum : 
				i=1 

			for ch in chdat :
				"dat : ", basis, ch
				fs[i] = fsurfacedecomp( ax[i] , int(ch) , levtot )
				print "Zmin Zmax [%d]: "%i, fs[i].Z.min(), fs[i].Z.max()
				i += 1
			if args.decompsumonly : 
				nax = 1 
				if args.locmax : 
					nax += 1 
				plt.close()
				if dd<1 : fig, ax = plt.subplots( 1, nax , sharex=True, sharey=True  )
				if nax < 2 : ax=[ax]
			if args.decompsum : 
				Zsum = fs[1].Z * 0.
				for j in range(nchdat) :
					Zsum += fs[j+1].Z
				if args.decompsumadd : 
				 	iax	= 1 
					for j in args.decompsumadd :
						fsdum = fsurfacedecomp( ax[iax] , int(j) , levtot )
						ax[0].clear()
						Zsum	+= fsdum.Z
						iax	+= 1 
					Zsum /= 2.
				print "Zmin Zmax [%d]: "%i, Zsum.min(), Zsum.max()
				fs[0] = fsurfacedecomp( ax[0] , -int(chdat[0]) , levtot, Zdat=Zsum , findlocmax=findlocmax )
				if args.locmax : 
					plotLocalMax( "{}/lattice/vdx/".format( hobj.ddir ) , ax[1], zthratio , hobj , args )

			if args.cb : 
				if args.decompsumonly :
					cax	 = ax
				else : 
					cax	 = ax.ravel().tolist()
				fig.colorbar( fs[0].cbar , ax=cax, extend="max")
			#print "CBlevel : ", fs[0].cbar.levels
			print "nCBlevel : ", len(fs[0].cbar.levels)

			if args.showaxes :  pass
			else :
				try : 
					for aax in ax.flat :
						aax.axis("off") 
				except :
					for aax in ax : 
						aax.axis("off") 
				if args.locmax : 
					ax[1].axis('on')
		elif args.decompall: 
			if dd<1 : fig, ax = plt.subplots( 2,4 )
			fs = [""]*7
			fs[0] = fsurface( ax[0][0] ) ;
			levtot = fs[0].cbar.levels
			i=1
			#if fswn != 0 : 
				#ax[1][0].remove()
			#else : 
			fspectdum = fspect
			faxisdum  = faxis
			fspect = dispfspect
			faxis  = dispfaxis
			if args.fsonly : pass
			else : kline( fig , cmapcolor , ax=ax[1][0]) 
			ax[1][0].axhline( fswn, 0,1 , color='w', dashes=[1,2] , lw=1) 
			ax[1][0].set_ylim( -1,1)
			axSetRange( ax[1][0], args.xxrange,  "x" )
			axSetRange( ax[1][0], args.yrange,   "y" )
			fspect = fspectdum
			faxis  = faxisdum
			for ch in chdat :
				"dat : ", basis, ch
				fs[i] = fsurfacedecomp( ax[0][i] , int(ch) , levtot )
				i += 1
			if basis.find("t")>-1 :   basis2 = "j"
			elif basis.find("j")>-1 : basis2 = "t"
			elif basis.find("s")>-1 : basis2 = "j"
			i=1
			chdat = chdat2
			for ch in chdat :
				basis1= basis
				basis = basis2
				"dat : ", basis, ch
				fs[i] = fsurfacedecomp( ax[1][i] , int(ch) , levtot )
				i += 1
			titleax = ax[0][1]
			if args.cb : 
				xinset = fig.add_axes([0.91,0.18,0.02,0.65])
				cbmax = int(levtot[-1])
				cbmin = int(levtot[0])
				cbmid = float( cbmax-cbmin ) /2.
				cbtick= [ cbmin, cbmid, cbmax ]
				fcb = fig.colorbar( fs[0].cbar , cax=xinset , ticks=cbtick )
				fcb.ax.tick_params(labelsize=fsize) 
			#print "CBlevel ", type(levtot), ": ", fs[0].cbar.levels
			print "nCBlevel : ", len(fs[0].cbar.levels)
			#ax[1][0].axis('scaled')
		else : 
			if dd<1 : fig, ax = plt.subplots( )
			fs = fsurface( ax )
			#ax.axis('off')
	else :
		if args.count : 
			count = int(args.count ) 
		else : 
			count = hobj.count
		cwd = os.getcwd()+'/'+hobj.ddir+"/lattice/vdx/" 
		fspect = "{}u{:.2f}_ts0.00_d{:.3f}_{}{:d}th{}.txt".format( dispfileind, hobj.UF, hobj.D, tagEpBr, count, anglelab ) 
		fspectdecomp = "decomp_{}u{:.2f}_ts0.00_d{:.3f}_{}{:d}th{}.txt".format( dispfileind, hobj.UF, hobj.D, tagEpBr, count, anglelab ) 
		faxis  = "{}naxis{}_Nplot{}_Nk{}_{}{:d}th{}.dx".format( dispfileind, naxis, nplot, nktot, tagEpBr, count, anglelab) 
		faxisdecomp  = "{}decomp_naxis{}_Nplot{}_Nk{}_{}{:d}th{}.dx".format( dispfileind, naxis, nplot, nktot, tagEpBr, count, anglelab) 
		fspect		= dispfilehead + fspect
		fspectdecomp	= dispfilehead + fspectdecomp
		fspect		= cwd + fspect
		fspectdecomp	= cwd + fspectdecomp
		faxis		= dispfilehead + faxis
		faxisdecomp	= dispfilehead + faxisdecomp
		faxis		= cwd + faxis
		faxisdecomp	= cwd + faxisdecomp
		if args.decomp: 
			nchdat	= len(chdat) 
			nax	= nchdat+1
			if args.decomponly : 
				nax = nchdat 
			if args.decompsum : 
				nax += 1 

			naxrow = 1
			if findlocmax : naxrow=2
			if dd<1 : fig, axarr = plt.subplots( naxrow, nax , sharex=True, sharey=True )
			print "subplots : ", naxrow, nax 
			if naxrow > 1 :
				ax	 = axarr[0,:]
				locaxarr = axarr[1,:]
			else :		ax = axarr
			kclass = [""]*nax
			if args.decomponly :
				i=0
				titleax = ax[len(chdat)/2]
			else : 
				i=1
				titleax = ax[i]
	
			kclass[0] = kline( fig , cmapcolor , ax=ax[0]) ;
			levtot = kclass[0].cbar.levels * len(chdat)
			if args.levtot :
				lmin, lmax, npoint = args.levtot
				levtot	= np.linspace( float(lmin), float(lmax), int(npoint), endpoint=True)
			print "levtot : ", type(levtot), np.shape(levtot), levtot
			print "levtot : ", levtot[0], levtot[-1]
			print "Zmax/min 0 : ", np.amax(kclass[0].Z.flatten()), np.amin(kclass[0].Z.flatten())
			ax[0].clear();
	
			i=1
			for ch in chdat :
				print "dat : ", basis, ch
				kclass[i] = klinedecomp( ax[i] , int(ch) , levtot , cmapcolor , findlocmax=findlocmax )
				print "Zmax/min %d : "%i, np.amax(kclass[i].Z.flatten()), np.amin(kclass[i].Z.flatten())
				if args.locmax : 
					plotLocalMax( "{}/lattice/vdx/".format( hobj.ddir ) , axarr[1,i], zthratio , hobj , args , fhead='disp' )
				i += 1

			if args.decompsumonly : 
				nax = 1 
				plt.close()
				if dd<1 : fig, axarr = plt.subplots( naxrow, nax , sharex=True, sharey=True  )
				if naxrow > 1 :
					if nax < 2 :	
						ax	 = [axarr[0]]
						locaxarr = [axarr[1]]
					else :		
						ax	 = [axarr[0,:]]
						locaxarr = axarr[1,:]
				else : 
					if nax < 2 :	ax=[axarr]
					else :		ax=axarr

			if args.decompsum : 
				Zsum = kclass[1].Z * 0.
				for j in range(nchdat) :
					Zsum += kclass[j+1].Z
				if args.decompsumadd : 
				 	iax	= 1 
					for j in args.decompsumadd :
						fsdum = klinedecomp( ax[0] , int(j) , levtot, cmapcolor )
						ax[0].clear()
						Zsum	+= fsdum.Z
						iax	+= 1 
					Zsum /= 2.
				kclass[0] = klinedecomp( ax[0] , -int(chdat[0]) , levtot, cmapcolor, Zdat=Zsum , findlocmax=findlocmax)
				print "Zmax/min %d : "%i, np.amax(kclass[0].Z.flatten()), np.amin(kclass[0].Z.flatten())
				if args.locmax : 
					plotLocalMax( "{}/lattice/vdx/".format( hobj.ddir ) , locaxarr[0], zthratio , hobj , args , fhead='disp' )

			if args.cb : 
				kclasscb	 =  kclass[1].cbar 
				if args.decompsumonly or args.decompsum :
					cax	 = ax
					kclasscb	=  kclass[0].cbar 
				else : 
					cax	 = ax.ravel().tolist()
				fig.colorbar( kclass[1].cbar , ax=cax ) #, extend="max")
			#if args.cb : 
			#	fig.colorbar( kclass[1].cbar , ax=ax.ravel().tolist())
			print "nCBlevel : ", len(kclass[0].cbar.levels)

			# Gamma-(pi,pi) cut
			#for tmpax in ax : 
			#	#axSetRange( tmpax, args.xxrange, "x" )
			#	#axSetRange( tmpax, args.yrange,  "y" )
			#	if args.xgrid :
			#		tmpax.grid( axis='x', linestyle=":")
			#	if args.ygrid :
			#		tmpax.grid( axis='y', linestyle=":")
			#	print "axis  orig : ", tmpax.axis()
			#	xtick = np.array( [ 0.5 , 0.6, 0.7 , 0.8, 1 ] ) 
			#	scaledToOrig = lambda y :  (y-1.)*(1.-0.70805369)/0.5 + 1.
			#	print scaledToOrig( xtick )
			#	tmpax.set_xticks( scaledToOrig( xtick )  )
			#	tmpax.set_xticklabels( xtick )
			#	tmpax.set_xlabel(r"$k / |k_{(\pi,\pi)}|$")
			#	if args.tickpad : 
			#		tickpad = float(args.tickpad)
			#		tmpax.tick_params( axis='both', pad = tickpad )
			# Gamma-(pi,0) cut
			for tmpax in ax : 
				#tmpax.set_xlabel(r"$|k / k_{(0,0),(\pi,0)}|$")
				xtickOrigK	= [0.7919463099999999, 1.]
				xtick		= np.array( [ "($\pi/2$,0)", r"($\pi$,0)" ] )
				tmpax.set_xticks( xtickOrigK )
				tmpax.set_xticklabels( xtick )
				if args.tickpad : 
					tickpad = float(args.tickpad)
					tmpax.tick_params( axis='both', pad = tickpad )
			ax[0].set_ylabel(r"$E$ [eV]") 
			axSetRange( ax[0], args.xxrange, "x" )
			axSetRange( ax[0], args.yrange,  "y" )
		elif args.decompall: 
			if dd<1 : fig, ax = plt.subplots( 2,4 )
			kclass = [""]*4
	
			kclass[0] = kline( fig , cmapcolor , ax=ax[0][0]) ;
			levtot = kclass[0].cbar.levels
	
			i=1
			if len(chdat)==3 :
				for ch in chdat :
					"dat : ", basis, ch
					kclass[i] = klinedecomp( ax[0][i] , int(ch) , levtot , cmapcolor )
					i += 1
			titleax = ax[0][1]
			if args.cb : 
				fig.colorbar( kclass[0].cbar , ax=ax.ravel().tolist())
			print "nCBlevel : ", len(kclass[0].cbar.levels)
	
			if basis.find("t")>-1 :   basis2 = "j"
			elif basis.find("j")>-1 : basis2 = "t"
			basis1 = basis
			basis  = basis2
			chdat1 = chdat
			chdat  = chdat2
			i=1
			if len(chdat)==3 :
				for ch in chdat :
					"dat : ", basis, ch
					kclass[i] = klinedecomp( ax[1][i] , int(ch) , levtot , cmapcolor )
					i += 1
	
			for i in range(1,4) : 
				for kk in range(2) : 
					ax[kk][i].tick_params( axis='both',          # changes apply to the x-axis
							which='both',      # both major and minor ticks are affected
							bottom='off',      # ticks along the bottom edge are off
							left='off',        
							labelleft='off',        
							labelbottom='off')
	
			fswn = 0 
			try : 
				fspect = cwd + "{}_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( fsfileind, hobj.UF, hobj.D, tagEpBr, count ) 
				fspectdecomp = cwd + "decomp_{}_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( fsfileind, hobj.UF, hobj.D, tagEpBr, count ) 
				faxis  = cwd + "{}_naxis{}_Nplot1_Nk40401_{}{:d}th{}.dx".format( fsfileind, naxis, tagEpBr, count, anglelab ) 
				faxisdecomp  = cwd + "decomp_{}_naxis{}_Nplot1_Nk40401_{}{:d}th{}.dx".format( fsfileind, naxis, tagEpBr, count, anglelab ) 
				fsurface( ax[1][0] ) ;
			except : 
				print "Warnning in"
				print fspect 
				print fspectdecomp 
				print faxis  
				print faxisdecomp  
				print "Try another option"
				fspect = cwd + "fs_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( hobj.UF, hobj.D, tagEpBr, count ) 
				fspectdecomp = cwd + "decomp_fs_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( hobj.UF, hobj.D, tagEpBr, count ) 
				faxis  = cwd + "fs_naxis{}_Nplot1_Nk40401_{}{:d}th{}.dx".format( naxis, tagEpBr, count, anglelab ) 
				faxisdecomp  = cwd + "decomp_fs_naxis{}_Nplot1_Nk40401_{}{:d}th{}.dx".format( naxis, tagEpBr, count, anglelab ) 
				fsurface( ax[1][0] ) ;
		elif args.proj: 
			#figs = 12
			#fs = 25 ; fsm = 25
			if dd<1 : fig, ax = plt.subplots( )
			kclass = [""]*4
	
			kclass[0] = kline( fig , cmapcolor[0] , ax=ax) ;
			levtot = kclass[0].cbar.levels
	
			i=1
			if len(chdat)==3 :
				for ch in chdat :
					"dat : ", basis, ch
					kclass[i] = klinedecomp( ax , int(ch) , levtot , cmapcolor[i] )
					i += 1
			titleax = ax
			if args.cb : 
				fig.colorbar( kclass[0].cbar , ax=ax )
			print "nCBlevel : ", len(kclass[0].cbar.levels)
		elif args.band : 
			bandfspect = "band_naxis{}_Nk{}{}.dat".format(naxis,nktot,anglelab)
			bandfaxis  = "naxis{}_Nplot1024_Nk{}_{}{:d}th{}.dx".format( naxis, nktot, tagEpBr, hobj.count, anglelab) 
			bandfspect = cwd + bandfspect
			bandfaxis  = cwd + bandfaxis
	
			if os.system("test -e {}".format(fspect)) < 1 :  pass
			else :
				print "Error :: ", fspect, "doesn't exist." 
				sys.exit(1)
			fig, ax = plt.subplots()
			bo = bandonly( fig , cmapcolor , ax=ax )
			axSetRange( ax, args.xxrange, "x" )
			axSetRange( ax, args.yrange,  "y" )
		elif args.locmax :
			fig, ax = plt.subplots(1,2, sharex=True, sharey=True ,figsize=(11,6) ) 
			kl = kline( fig , cmapcolor , ax=ax[0], findlocmax=findlocmax ) 
			print "nCBlevel : ", len(kl.cbar.levels), "[",kl.cbar.levels[0], ", ..., ", kl.cbar.levels[-1], "]"
			print "zmin,zmax : ", kl.cbar.zmin, kl.cbar.zmax

			plotLocalMax( "{}/lattice/vdx/".format( hobj.ddir ) , ax[1], zthratio , hobj , args, fhead='disp' )
		else :
			#fig = plt.figure(figsize=(figs,figs))
			#ax = plt.subplot2grid((5, 5), (0, 1), colspan=4, rowspan=4 )
			if dd<1 : fig, axone = plt.subplots( ) ; ax = [axone]*5
			if args.ldos : 
				args.trans=True
				basis = "p"
				if args.basis : basis   = args.basis 
				if basis.find("t")>-1  :
					fhead   = "tldos"
					if basis.find("z")>-1  :
						fhead   = "ldos"
				elif basis.find("j")>-1  :
					fhead   = "ldos"
				else : 
					fhead   = "tldos"
				from matplotlib import gridspec


				if args.vertself :
					gs = gridspec.GridSpec(1, 4, width_ratios=[3, 1, 1, 1]) 
				else :
					gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 

				ax0 = plt.subplot(gs[0])
				kl = kline( fig , cmapcolor , ax=ax0) 

				if args.nplotdos : 
					hobj.Nplot = int(args.nplotdos)
				ldos = ldos(args, fname=fhead , basis=basis )
				ax1 = plt.subplot(gs[1])
				ldos.totally( args, fig, ax1 , hobj=hobj, chdat=chdat)

				ax1.set_yticklabels([])
				ax1.set_ylabel("")
				ax1.legend()
				#axSetRange( ax, args.xxrange, "x" )
				axSetRange( ax0, args.yrange,  "y" )
				axSetRange( ax1, args.yrange,  "y" )

				hobj.Nplot = nplot

				if args.vertself :
					ax2 = plt.subplot(gs[2])
					ldos1 = ldos(args , basis=basis , vertdata=True, tagEpBr=tagEp ,vertSelf=True )
					ldos1.totally( args, fig, ax2 , hobj=hobj, chdat=chdat)

					ax3 = plt.subplot(gs[3])
					ldos1 = ldos(args , basis=basis , vertdata=True, tagEpBr=tagEp ,vertSelf=True , vertSelfRealpart=True )
					ldos1.totally( args, fig, ax3 , hobj=hobj, chdat=chdat)

			elif args.vertdos : 
				args.trans=True
				basis = "p"
				if args.basis : basis   = args.basis 
				if basis.find("t")>-1  :
					if basis.find("z")>-1  :
						fhead   = ""
					else : 
						fhead   = "t"
				elif basis.find("j")>-1  :
					fhead   = ""
				else : 
					fhead   = "t"
				tagEp = "ep{:g}_".format(hobj.epsilon)
				if args.broadeningtag :	tagEp = tagEpBr

				from matplotlib import gridspec
				if args.vertself :
					gs = gridspec.GridSpec(1, 4, width_ratios=[3, 1, 1, 1]) 
				else : 
					gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 


				ax0 = plt.subplot(gs[0])
				kl = kline( fig , cmapcolor , ax=ax0) 

				tagEp = "ep{:.3f}_".format(hobj.epsilon)
				if args.broadeningtag :	tagEp = tagEpBr

				axdosarr = []
				print "CHDAT : ", chdat
				ax1 = plt.subplot(gs[1])
				ldos1 = ldos(args , basis=basis , fname=fhead+"ldos", tagEpBr=tagEp )
				ldos1.totally( args, fig, ax1 , hobj=hobj, chdat=chdat)
				axdosarr.append( ax1 )
				axSetRange( ax1, args.xxrangedos, "x" )

				if args. vertself : 
					ax2 = plt.subplot(gs[2])
					ldos1 = ldos(args , basis=basis , fname=fhead+"selfreal", tagEpBr=tagEp , vertSelf=True )
					ldos1.totally( args, fig, ax2 , hobj=hobj, chdat=chdat)

					ax3 = plt.subplot(gs[3])
					ldos1 = ldos(args , basis=basis , fname=fhead+"selfreal", tagEpBr=tagEp, vertSelfRealpart=True )
					ldos1.totally( args, fig, ax3 , hobj=hobj, chdat=chdat)
					axdosarr.append( ax2 )
					axdosarr.append( ax3 )

				titleax = ax0
				for axnow in axdosarr : 
					axnow.set_yticklabels([])
					axnow.set_ylabel("")
					axnow.legend(loc=4)
					axSetRange( axnow, args.yrange,  "y" )
				#axSetRange( ax, args.xxrange, "x" )
				#axSetRange( ax0, args.yrange,  "y" )
				#axSetRange( ax1, args.yrange,  "y" )
			else : 
				lev=300
				if args.zrange :
					zmin = float( args.zrange.split("_")[-2] ) 
					zmax = float( args.zrange.split("_")[-1] ) 
					lev	= np.linspace(zmin,zmax, 300 ) 
				kclass = [ kline( fig , cmapcolor , ax=axone, lev=lev )  ]

				#levtot = kclass[0].cbar.levels
				#if args.cb : 
				#	xinset = fig.add_axes([0.91,0.18,0.02,0.65])
				#	cbmax = int(levtot[-1])
				#	cbmin = int(levtot[0])
				#	cbmid = float( cbmax-cbmin ) /2.
				#	cbtick= [ cbmin, cbmid, cbmax ]
				#	fcb = fig.colorbar( kclass[0].cbar , cax=xinset , ticks=cbtick )
				#	fcb.ax.tick_params(labelsize=fsize) 
				if args.cb : 
					fig.colorbar( kclass[0].cbar , ax=axone )

				print "nCBlevel : ", len(kclass[0].cbar.levels), "[",kclass[0].cbar.levels[0], ", ..., ", kclass[0].cbar.levels[-1], "]"
				print "zmin,zmax : ", kclass[0].cbar.zmin, kclass[0].cbar.zmax
				#print "dir(kclass[0]ine) : ", dir(kclass[0].cbar)
		if args.lineplotk :
			if dd<1 : axLineplot = fig.add_axes([0.3,0.3,0.3,0.3])
			axone	= axLineplot
			kl = kclass[0]
			linearg = "y_0.5"
			if args.lineplotkw : linearg = args.lineplotkw 
			lineaxis = linearg.split("_")[0]
			n0, xdat, ydat, datsection = refineLinePlot( linearg , kl )
			if args.leglabel :	leglabel = "${}={}$".format(args.leglabel ,  getattr( hobj, args.leglabel ) ) 
			else : 			leglabel = linearg
			axone.plot( xdat[n0], ydat[n0] , ".-", lw = 1, label=leglabel ) #, markerfacecolor='w')#, markeredgecolor='k')
			saveLinePlot( hobj.ddir , linearg , n0, xdat[n0], ydat[n0] , hobj, fhead='vert' )
			#plotIndLine( ax[0], datsection, lineaxis, (2,1) , 'w' )  if args.lineplotonly is None else '.'
			if   lineaxis=="x" : axone.set_xlabel("$A(|k|={:.3f})$".format(datsection)) #; axone.tick_params( axis="y", which='both', labelleft=False  )#, labelsize = args.fsm ) 
			elif lineaxis=="y" : axone.set_ylabel("$A(\omega={:.3f})$".format(datsection)) #; axone.tick_params( axis="x", which='both', labelbottom=False  )#, labelsize = args.fsm ) 
			axone.set_xlim(0,1)
			axSetRange( axone, args.lineplotxrange, "x" )
			axSetRange( axone, args.lineplotyrange, "y" )
			ymax = axone.axis()[-1]
			#axone.text( 0.8, ymax*0.95, "$\omega$={:.4f}".format( datsection ) )
			if args.noleg : pass
			else : axone.legend()
			if args.lineplotfit : 
				ik0 =  50./341
				ik1 = 100./341
				for ix in range( len(xdat[n0]) ) :
					if xdat[n0][ix]>ik1 : 
						x1 = ix
						break
				for ix in range( len(xdat[n0]) ) :
					if xdat[n0][ix]>ik0 : 
						x0 = ix
						break
				print "GM-line data ind=[",x0,x1,"] :: x,y=(", xdat[n0][[0,x0,x1]],")(", ydat[n0][[0,x0,x1]] ,")"
				nxfit = 400
			  	xfit = np.linspace( 0.15,ik1, nxfit ) 
			  	yfit1= range( nxfit ) 
			  	yfit2= range( nxfit ) 
			  	yfitall= range( nxfit ) 

				#param1 = [ 15.87, 0.059, 0.0162/2 ] #[ ymax, xmax, sigmax ]
				#param2 = [ 13.04, 0.113, 0.0100/2 ] #[ ymax, xmax, sigmax ]
			  	print "IK1\t: ", ik1
				param1 = [ 23.34, ik1-0.0641, 0.0170/2 ] #[ ymax, xmax, sigmax ]
				param2 = [ 23.18, ik1-0.0718, 0.0150/2 ] #[ ymax, xmax, sigmax ]
				param1 = [ 250.74 ,0.2604 ,0.000518 ]
				param2 = [ 303.8  ,0.2292 ,0.000568 ]
				if args.readfitparam :
					try : 
						paramfile ="inversetau/tmptauparam" 
						paramread = np.genfromtxt( paramfile , dtype=float ) 
						param1 = paramread[0:3]
						param2 = paramread[3:6]
						print "READ : ", paramfile 
					except : 
						print "WARNNING :: reading %s"%paramfile, "is failed."
				text1 = "(1)Ymax,X0,sigmaX : {:.3f}, {:.3f}, {:.3f}".format( param1[0], param1[1], param1[2] ) 
				text2 = "(2)Ymax,X0,sigmaX : {:.3f}, {:.3f}, {:.3f}".format( param2[0], param2[1], param2[2] ) 
				for a in range(nxfit) : 
					yfit1[a] = gaussftn( xfit[a], param1[0], param1[1], param1[2] ) 
					yfit2[a] = gaussftn( xfit[a], param2[0], param2[1], param2[2] )
					yfitall[a] = yfit1[a] + yfit2[a] 
					#print "Gxy : ", xfit[a], "\t", yfit1[a],"\t", yfit2[a]
				#axone.plot( xfit, yfit1 , "k:" , label = text1 )
				#axone.plot( xfit, yfit2 , "k:" , label = text2 )
				#axone.plot( xfit, yfitall , "r-" )
				#axone.legend()
				#axone.set_xticks(np.linspace(0,0.3,3))
				#axone.set_xticklabels(np.linspace(0,0.3,3))
				print "XTICK\t: ", axone.get_xticks()
				vFermigamma = 1.27
				vFermibeta  = 2.53
				invtaugamma = ( param1[2] /ik1 *np.pi/ 3.833 ) * vFermigamma
				invtaubeta  = ( param2[2] /ik1 *np.pi/ 3.833 ) * vFermibeta
				invtaugammaInit = invtaugamma 
				invtaubetaInit  = invtaubeta  
				from scipy.optimize import least_squares
				t_train = xdat[n0][x0:x1+1]
				y_train = ydat[n0][x0:x1+1]*6
				t_trainmirror = xdat[n0][x0:x1+1] * -1. + 200./341 
				y_trainall = np.append( y_train, y_train ) 
				t_trainall = np.append( t_train, t_trainmirror ) 
				paramfit = [ param1[0], param1[1], param1[2], param2[0], param2[1], param2[2] ]
				#res_log = least_squares(func, paramfit, loss='cauchy', f_scale=0.1, args=(t_train, y_train) )
				#t_fit = np.linspace( xdat[n0][x0], xdat[n0][x1], 300 ) 
				#y_fit = gendata(t_fit, res_log.x)
				if args.noleastsq : 
					paramres = paramfit 
				else : 
					lsqloss = 'cauchy'
					if args.leastsquare : lsqloss = args.leastsquare
					res_log = least_squares(funcmirrorpi, paramfit, loss=lsqloss, f_scale=1.1, args=(t_trainall, y_trainall) )
					paramres = res_log.x
				t_fit = np.linspace( xdat[n0][x0], xdat[n0][x1], 300 ) 
				y_fit = gendatamirrorpi(t_fit, paramres)

				textfit = "Ymax,X0,sigmaX :\n{:.3f}, {:.3f}, {:.3f}\n".format( paramres[0], paramres[1], paramres[2] )  +  "{:.3f}, {:.3f}, {:.3f}".format( paramres[3], paramres[4], paramres[5] ) 
				axone.plot( t_fit, y_fit , "r-" , label=textfit )
				#axone.plot( t_trainall, y_trainall , "b-" )
				axone.legend()
				print "Fit params\t: ",
				for a in paramfit : print "{:9.4g}".format(a),
				print "\t(initial condition)"
				print "Fit params\t: ",
				for a in paramres : print "{:9.4g}".format(a),
				print "\t(fitted result)"
				invtaugamma = ( paramres[2] /ik1 *np.pi/ 3.833 ) * vFermigamma
				invtaubeta  = ( paramres[5] /ik1 *np.pi/ 3.833 ) * vFermibeta
				print "invtau(beta,gamma)[eV]\t: ", "{:9.6g}".format(invtaubetaInit),	"{:9.6g}".format(invtaugammaInit)	, "\t(initial condition)"
				print "invtau(beta,gamma)[eV]\t: ", "{:9.6g}".format(invtaubeta), 	"{:9.6g}".format(invtaugamma)	 	, "\t(fitted result)"
				if args.savefitparam :
					fr = open( "inversetau/inversetau.log" , 'a' ) 
					fr.write( "{:10.6f}\t{:10.6f}\t{:10.6f}\t".format( kl.Y[n0][0] , invtaubeta, invtaugamma )  )
					for itemparam in paramres : 
						fr.write( "{:10.6f}\t".format( itemparam ) )
					fr.write( "#{}\n".format( hobj.tdir )  )
					fr.close()
					fr = open( "inversetau/U{}S{}J{}_inversetau_ep{}.dat".format( hobj.UF, hobj.S, hobj.J, epsilon) , 'a' ) 
					fr.write( "{:10.6f}\t{:10.6f}\t{:10.6f}\t".format( kl.Y[n0][0] , invtaubeta, invtaugamma )  )
					for itemparam in paramres : 
						fr.write( "{:10.6f}\t".format( itemparam ) )
					fr.write( "#{}\n".format( hobj.tdir )  )
					fr.close()
					np.savetxt( "{}/inversetau.log".format(hobj.tdir) , paramres , fmt="%12.6f")

				nxfit = 300
			  	xfit = np.linspace( 0.15,ik1, nxfit ) 
			  	yfit1= range( nxfit ) 
			  	yfit2= range( nxfit ) 
			  	yfit1m= range( nxfit ) 
			  	yfit2m= range( nxfit ) 
			  	yfitall= range( nxfit ) 
				if xfit[0]>paramres[1] or xfit[-1]<paramres[1] :  paramres[1] = -paramres[1]
				if xfit[0]>paramres[4] or xfit[-1]<paramres[4] :  paramres[4] = -paramres[4]
				for a in range(nxfit) : 
					yfit1[a] = gaussftn( xfit[a], paramres[0], paramres[1], paramres[2] ) 
					yfit2[a] = gaussftn( xfit[a], paramres[3], paramres[4], paramres[5] )
					yfit1m[a] = gaussftn( xfit[a]*-1.+2*ik1, paramres[0], paramres[1], paramres[2] )
					yfit2m[a] = gaussftn( xfit[a]*-1.+2*ik1, paramres[3], paramres[4], paramres[5] )
				axone.plot( xfit, yfit1 , "k:" , label = text1 )
				axone.plot( xfit, yfit2 , "k--" , label = text2 )
				axone.plot( xfit, yfit1m , "b:" )
				axone.plot( xfit, yfit2m , "b:" )
				ftau = open( "inversetau/tmptauparam" , 'w' )  ; np.savetxt( ftau, np.array(paramres) , delimiter='\t', fmt='%10.6f')  ; ftau.close()
					
			#else : 
			#	axone.set_xticks(kl.xtics)
			#	axone.set_xticklabels(kl.knames) 
			print "YDATpart(frequency) ", n0, "/",np.shape(kl.Y) ," : ", [ kl.Y[a][0] for a in range(n0-3,n0+3) ] , "(dn={}".format(1./1024),")"
	
	if args.oldsoc : hobj.S = hobj.S*2./3.
	if args.notitle :  pass
	else :
		if args.locmax : 
			#ax[0].set_title("U={}, J={}, S={}, D={} (nb={})".format( UF, J, S, D, Nb ) )
			hobj.headSetTitle( ax[0] ) 
		elif args.lineplot : 
			#ax[0].set_title("U={}, J={}, S={}, D={} (nb={})".format( UF, J, S, D, Nb ) )
			hobj.headSetTitle( ax[0] ) 
		elif args.decomp or args.decompall: 
			#titleax.set_title("U={}, J={}, S={}, D={} (nb={})".format( UF, J, S, D, Nb ) )
			hobj.headSetTitle( titleax ) 
		else :
			#plt.title("U={}, J={}, S={}, D={} (nb={})".format( UF, J, S, D, Nb ) )
			try : 
				hobj.headSetTitle( titleax ) 
			except : 
				hobj.headPltTitle( plt ) 

try : 
	for axnow in axarr : 
		pass
	axflatten	= axarr
except :
	if nax==1 :
		axflatten	= [ax]
	else : 
		axflatten	= ax
	#for axnow in ax : 
	#	pass
	#axflatten	= ax
		
try : 
	for axnow in axflatten : 
		if args.noxtick :
			axnow.tick_params( axis="x", left=True,  bottom=False,  labelbottom=False, labelleft=False  )
		if args.noytick :
			axnow.tick_params( axis="y", left=False, bottom=True,   labelbottom=False, labelleft=False  )
		if args.noxlabel :
			axnow.set_xlabel("")
		if args.noylabel :
			axnow.set_ylabel("")
		if args.notick :
			#axnow.tick_params( axis="both", left=True,   labelbottom=False, labelleft=False  )
			axnow.tick_params( axis="both", left=False,   bottom=False, labelbottom=False, labelleft=False  )
except:
	pass

if args.fs :  
	ladjust=( 0.05 if (args.ladjust is None) else float(args.ladjust) )
	radjust=( 0.95 if (args.radjust is None) else float(args.radjust) )
	badjust=( 0.05 if (args.badjust is None) else float(args.badjust) )
	tadjust=( 0.95 if (args.tadjust is None) else float(args.tadjust) )
	wspace =( 0.02 if (args.wspace  is None) else float(args.wspace ) )
	hspace =( 0.00 if (args.hspace  is None) else float(args.hspace ) )
	if args.decompall :
		ladjust=( 0.11 if (args.ladjust is None) else float(args.ladjust) )
		radjust=( 0.99 if (args.radjust is None) else float(args.radjust) )
		badjust=( 0.07 if (args.badjust is None) else float(args.badjust) )
		tadjust=( 0.96 if (args.tadjust is None) else float(args.tadjust) )
		wspace =( 0.02 if (args.wspace  is None) else float(args.wspace ) )
		hspace =( 0.02 if (args.hspace  is None) else float(args.hspace ) )
		ax[1][0].xaxis.labelpad = 5
	plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, wspace=wspace , hspace=hspace )
else : 
	if args.decompall :
		ladjust=( 0.20 if (args.ladjust is None) else float(args.ladjust) )
		hspace =( 0.40 if (args.hspace  is None) else float(args.hspace ) )
		plt.subplots_adjust(left=ladjust , hspace=hspace )
	else : 
		#plt.subplots_adjust(left=0.2)
		ladjust=( 0.13 if (args.ladjust is None) else float(args.ladjust) )
		radjust=( 0.95 if (args.radjust is None) else float(args.radjust) )
		badjust=( 0.10 if (args.badjust is None) else float(args.badjust) )
		tadjust=( 0.93 if (args.tadjust is None) else float(args.tadjust) )
		wspace =( 0.10 if (args.wspace  is None) else float(args.wspace ) )
		hspace =( 0.10 if (args.hspace  is None) else float(args.hspace ) )
		if args.ldos : 
			wspace =( 0.01 if (args.wspace  is None) else float(args.wspace ) )
			hspace =( 0.00 if (args.hspace  is None) else float(args.hspace ) )
		plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, wspace=wspace , hspace=hspace )
if args.tl :  plt.tight_layout()
	
if args.locmaxonly and args.nolab : 
	ax.tick_params( axis='both',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom='off',      # ticks along the bottom edge are off
			left='off',        
			labelleft='off',        
			labelbottom='off')
if args.sublabel :
	if len(ax.flatten())>1  :
		for dd in range(len(ax.flatten())) :
			#print "tF[{}] : ".format(dd), type(ax.flatten())
			axnow = ax.flatten()[dd]
			x0,x1,y0,y1 = axnow.axis()
		        dx=x1-x0;dy=y1-y0#; axnow.text( dx*0.02, dy*0.98 , markerssub[dd] ) 
			slab = markerssub[dd]
			if args.sublabeloff : slab = markerssub[dd+int(args.sublabeloff)]
			if args.sublabelparameter :  slab = slab +" {} = {:.6g}".format( filterlabel(args.sublabelparameter) , getattr(hobj,args.sublabelparameter) )
		        axnow.text(0.01, 0.92, slab,   
	    				horizontalalignment='left', verticalalignment='center', transform=axnow.transAxes,
		                        fontsize=fsize , color='white' , weight='semibold')

if args.xgrid :
	plt.grid( axis='x', linestyle=":")
if args.ygrid :
	plt.grid( axis='y', linestyle=":")

if args.cmapnote : 
	textc	= args.labelcolor if args.labelcolor else 'k'
	textpos	= args.labelposition if args.labelposition else [0.9,1.02] 
	textpos	= [0.9,1.02] 
	leglabel = cmap
	try :
		titleax.text( textpos[0], textpos[1], leglabel, color=textc, transform=titleax.transAxes )
	except :
		ax[0].text( textpos[0], textpos[1], leglabel, color=textc, transform=ax[0].transAxes )
		pass


ndpi = 400 
print "length(argsavefigname)" , len(argsavefigname)
if len(argsavefigname)>300 : argsavefigname = argsavefigname[:300]
fsave = argsavefigname
transparent = False 
if args.savetransparent : transparent = True
if args.save or args.cptmp :
	dtype="pdf"
	if args.pdf : 
		dtype="pdf"
	if args.png : 
		dtype="png"
	if args.saveeps : 
		dtype="eps"
	print "Saved : ", fsave + "." + dtype
	plt.savefig( fsave + "."+dtype , dpi=ndpi, transparent=transparent) 
	cptmpftn( fsave, cmd="scp", destMac=False , dtype=dtype )
	#print "saved in : ", argsavefigname
	#plt.savefig( argsavefigname + ".png" , dpi=ndpi , transparent=transparent)
#	if args.savepdf : 
#		plt.savefig( argsavefigname + ".pdf" , transparent=transparent)
#	if args.saveeps : 
#		plt.savefig( argsavefigname + ".eps" , transparent=transparent)
#	if args.cptmp :
#		os.system( "cp " + argsavefigname + ".png ~/Dropbox/tmp_linux/"  )
else :
	#print "saved in : ", argsavefigname
	#plt.savefig( argsavefigname + ".png" )
	#if args.png : 
	#	plt.savefig( argsavefigname + ".png" , transparent=transparent)
	#if args.pdf : 
	#	plt.savefig( argsavefigname + ".pdf" , transparent=transparent)
	#if args.saveeps : 
	#	plt.savefig( argsavefigname + ".eps" , transparent=transparent)
	#plt.tight_layout()
	plt.show()
