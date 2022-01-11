#!/opt/python/2.7.15/gcc-4.8.5/bin/python
import sys
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.gridspec as gridspec

sys.path.append("/home/jun/tools/dmft/plot/" )
from basislabel import *
from hjl_common import *
from hjl_class  import *

import argparse
parser = argparse.ArgumentParser(description="Plot jzeval" ) 
parser.add_argument("--fname", "-f", type=str, help="Specify dat.")
parser.add_argument("--fnamearr", "--fa", nargs='+', type=str, help="Specify a list of dat.")
parser.add_argument("--fprefix", "--fp", "--fpre", type=str, help="Specify prefix of data file.")
parser.add_argument("--inverse", "--inv", "--iv", action='store_true', help="Use the inverse of the data.")
parser.add_argument("--inverseconstant", "--invconst", "--invc", "--ivc", type=str, help="Substitute the data by a value before taking its inverse.")
parser.add_argument("--Tfit", "--tfit", "--tf", action='store_true', help="Use T_fit instead of beta_fit.")
parser.add_argument("-Z", "--Zfactor", "--Zfac", action='store_true', help="Use Z instead of m^*.")
parser.add_argument("--beta", "-b", type=str, help="Specify beta (default=128)." , default="128" ) 
parser.add_argument("--IT", "-T", type=str, help="Specify beta_real (default=128)." )
parser.add_argument("-J", "--JH", type=str, help="Specify J_H (default=0.4)." ) 
parser.add_argument("--JHarr", "--Jarr", "--jarr", action='store_true', help="Draw with JH-array.")
parser.add_argument("--betaarr", "--barr", action='store_true', help="Draw with beta-array.")
parser.add_argument("--ITarr", "--ITa", action='store_true', help="Draw with IT-array.")
parser.add_argument("--setbetaarr", "--setbarr", nargs='+', type=int, help="Set the beta-array.")
parser.add_argument("-S", "--soc", "--SOC", type=str, help="Specify SOC (default=0.1)." ) 
parser.add_argument("--multiplybeta", "--mulbeta", "--mub", action='store_true', help="Multiply the data by 'beta'.")
parser.add_argument("--multiplycomponent", "--mulcomp", "--mcomp", type=str, help="Multiply the data by the other specified data.")
parser.add_argument("--plot", "-p", type=str, help="Specify a quantity you want to plot among ['vv','lt','ltb','tot'] (default=vv)." )
parser.add_argument("--logy", "--ly", action='store_true', help="Draw y-data in log-scale." )
parser.add_argument("--logx", "--lx", action='store_true', help="Draw x-data in log-scale." )
parser.add_argument("--loglog", "--ll", action='store_true', help="Draw data in log-scale." )
parser.add_argument("--cmap", "--cm", type=str, help="Specify a colormap (default=jet)." , default="jet") 
parser.add_argument("--cptmp", action="store_true", help="Copy the saved figure into Dropbox/tmp_linux/") 
parser.add_argument("--pdf", action="store_true", help="Save the figure in .pdf" ) 
parser.add_argument("--png", action="store_true", help="Save the figure in .png" ) 
parser.add_argument("--eps", action="store_true", help="Save the figure in .eps" ) 
parser.add_argument("--ndpi", type=str, help="Set ndpi, number of dots per inch." ) 
parser.add_argument("--ladjust", "--lad", type=str, help="Set   left-adjust value." )
parser.add_argument("--radjust", "--rad", type=str, help="Set  right-adjust value." )
parser.add_argument("--badjust", "--bad", type=str, help="Set bottom-adjust value." )
parser.add_argument("--tadjust", "--tad", type=str, help="Set    top-adjust value." )
parser.add_argument("--wspace",  "--wsp", type=str, help="Set wspace-adjust value." ) 
parser.add_argument("--hspace",  "--hsp", type=str, help="Set hspace-adjust value." ) 
parser.add_argument("--plotx", "--px", "--xobj", type=str, help="Specify an object for x-axis.")
parser.add_argument("--ploty", "--py", "--yobj", type=str, help="Specify an object for y-axis.")
parser.add_argument("--plotyerr", "--pyerr", "--yerrobj", type=str, nargs='+', help="Specify an error object for y-axis.")
parser.add_argument("--typex", "--tx", type=str, help="Specify the type of object for x-axis.")
parser.add_argument("--typey", "--ty", type=str, help="Specify the type of object for y-axis.")
parser.add_argument("--autolegend", "--alegend", "--aleg", type=str, nargs='+', help="Specify parameters to be used in labelling in legend automatically amogn S, J, and T. ")
parser.add_argument("--nolegend", "--noleg", action='store_true', help="Remove lenged.")
parser.add_argument("--legend", "--leg", "--legpar", type=str, nargs='+', help="Specify parameters and draw them in legend.")
parser.add_argument("--ylabelcomp", "--ylcomp", "--ylabcomp", type=str, help="Specify the y-axis label as some quantity.")
parser.add_argument("--xlabelcomp", "--xlcomp", "--xlabcomp", type=str, help="Specify the x-axis label as some quantity.")
parser.add_argument("--noylabel", "--noylab", "--noyl", action="store_true", help="Remove the y-axis label.")
parser.add_argument("--noxlabel", "--noxlab", "--noxl", action="store_true", help="Remove the x-axis label.")
parser.add_argument("--noytick", "--noytic", action="store_true", help="Remove y-ticks." )
parser.add_argument("--noxtick", "--noxtic", action="store_true", help="Remove x-ticks." )
parser.add_argument("--xticklabel", "--xticklab", "--xtlab", type=float, nargs='+', help="Specify x-ticklabel.")
parser.add_argument("--yticklabel", "--yticklab", "--ytlab", type=float, nargs='+', help="Specify y-ticklabel.")
parser.add_argument("--ylabelcoord", "--ylabcoord", "--ylabc", "--ylc", type=str, help="Set the x-coordinate of y-axis label.")
parser.add_argument("--sublabel", "--slab", "--sublab", action="store_true", help="Use the labels of sub-figures with alphabet." ) 
parser.add_argument("--sublabelarr", "--slabarr", "--sublabarr", type=str, nargs='+', help="Specify the labels of sub-figures with alphabet." ) 
parser.add_argument("--sublabelparameter", "--slabpar", "--slabp", "--sublabpar", nargs='+', type=str, help="Set a parameter of the sub-labels.")
parser.add_argument("--sublabelnoalphabet", "--slabnoalpha", "--sublabnoalpha", action="store_true", help="Set the labels of sub-figures without alphabet." ) 
parser.add_argument("--filelabel", "--flab", "--filelab", type=str, help="Set the filename-label.")
parser.add_argument("-y", "--yrange", type=str, help="set yrange e.g. _y0_y1" )
parser.add_argument("-x", "--xxrange", type=str, help="set xrange e.g. _x0_x1" )
parser.add_argument("--y2range", "--y2", type=str, help="set yrange e.g. _y0_y1" )
parser.add_argument("--xx2range","--x2",  type=str, help="set xrange e.g. _x0_x1" )
parser.add_argument("--figs", "--figsize", type=str, help="Set the labels of the legend" )
parser.add_argument("--fontsize", "--fonts", "--fs", type=str, help="Set font size" )
parser.add_argument("--fontsizelegend", "--fontsleg", "--fsleg", type=float, help="Set font size of lengend.")
parser.add_argument("--fontsizetitle", "--fontstitle", "--fstitle", type=float, help="Set font size of title.")
parser.add_argument("--ngrid", "--ng", "-n", type=str, help="Set the number of grid of data." )
parser.add_argument("--datapoint", "--dp", "--scatter", action='store_true', help="Draw the data points.")
parser.add_argument("--colormesh", "--cmesh", "--pmesh", "--pcolormesh", action='store_true', help="Use colormesh() instead of contourf().")
parser.add_argument("--diff", "--difftwo", action='store_true', help="Use labels of differential or derivatives.")
parser.add_argument("--diffrange", "--diffr", "--drange", type=str, help="Set the boundary in plotting the difference.")
parser.add_argument("--savetransparent", "--savetp", "--tp", action='store_true',  help="Save the figure with transparent background.")
parser.add_argument("--savenottransparent", "--savenottp", "--saventp", "--nottp", action='store_true',  help="Save the figure with transparent background.")
parser.add_argument("--hidelowfreq", "--hidelf", "--hide", "--hidelow", action='store_true',  help="Hide low-frequency part below pi/beta.")
parser.add_argument("--linelowfreq", "--linelf", "--linelow", action='store_true',  help="Draw a line of pi/beta.")
parser.add_argument("--omegafl", "--wfl", "--ofl", action='store_true',  help="Draw omega_FL from devitation-starting points.")
parser.add_argument("--omegafl2", "--wfl2", "--ofl2", action='store_true',  help="Draw omega_FL from extremum values.")
parser.add_argument("--omegaflinset", "--wflinset", "--oflinset", "--wflin", "--oflin", action='store_true',  help="Draw in inset omega_FL from devitation-starting points.")
parser.add_argument("--inset", "--is", action='store_true',  help="Draw an inset.")
parser.add_argument("--zerosoc", "--zsoc", action='store_true',  help="Plot the SOC=0 data.")
parser.add_argument("--operator", "--op", type=str,  help="Set type of operator.")
parser.add_argument("--particlehole", "--ph", "--ptlh", action='store_true', help="Plot the particle-hole partner.")
parser.add_argument("--particleholeonly", "--phonly", "--ptlhonly", action='store_true', help="Plot the particle-hole partner only.")
parser.add_argument("--nparticlehole", "--nph", "--nptlh", type=str, help="Set the maximum number of  particle-hole.")
parser.add_argument("--mergeparticlehole", "--mergeph", "--mph", "--mptlh", action='store_true', help="Merge and plot the particle-hole partner.")
parser.add_argument("--fitlinear", "--fitlin", "--flin", action='store_true', help="Fit the data to the linear-function.")
parser.add_argument("--fitpower", "--fitpow", "--fpow", action='store_true', help="Fit the data to the power-function.")
parser.add_argument("--nfitdat", "--nfit", "--nfd", type=str, help="Set the number of data to be fitted.")
parser.add_argument("--extrapolate", "--extrap", "--epol", action='store_true', help="Exptrapolate the data.")
parser.add_argument("--fitquadratic", "--fitquad", "--fquad", action='store_true', help="Fit the data to the quadratic-function.")
parser.add_argument("--fitinverse", "--fitinv", "--finv", action='store_true', help="Fit the data to a invsere-function.")
parser.add_argument("--fitparam", "--fitpar", "--fpar", action='store_true', help="Show the fitting parameters.")
parser.add_argument("--fitfunction", "--fitfunc", "--ffunc", action='store_true', help="Show the fitting function.")
parser.add_argument("--fitparamtitle", "--fitpartitle", "--fptitle", action='store_true', help="Show the fitting param as title.")
parser.add_argument("--interpolate", "--interp", "--ip", "--interpolate1d", "--ip1d", action='store_true', help="Interpolate data.")
parser.add_argument("--interpolaterange", "--interrange", "--irange", "--ipr", type=float, nargs='+', help="Set the range of values you will pick.")
parser.add_argument("--interpolatemethod", "--interpmethod", "--im", "--ipmethod", "--ipm", type=str, help="Set the direction of axis in interpolation. (default='linear')")
parser.add_argument("--ninterpolate", "--ninterp", "--nip", type=str, help="Set the grid number of interpolation.")
parser.add_argument("--clipon", "--clip", "--cl", "--clon", action='store_true', help="Place the data symbol behind the frame.", default=False)
parser.add_argument("--multiyobj", "--my", "--myobj", "--mpy", type=str, nargs='+', help="Plot a list of y-object.")
parser.add_argument("--multiyobjsplit", "--mys", "--myobjsp", "--mpysp", action='store_true', help="Plot a list of y-object matching with number of files.")
parser.add_argument("--subplotarrobj", "--subobj", "--subpobj", "--spobj", type=str, help="Set which object to be plotted in each subplots. e.g. 'dir' 'y' ")
parser.add_argument("--paper", action="store_true" , help="Plot as paper-use." ) 
parser.add_argument("--solidall", "--solida", action="store_true" , help="Plot as solid-lines." ) 
parser.add_argument("--ncolumnlegend", "--ncolleg", "--ncoll", type=int, help="Set the number of column of legend.")
parser.add_argument("--legendlocation", "--legloc", type=str, help="Set the location of legend.")
parser.add_argument("--legendlabel", "--leglab", type=str, help="Set the label to be used in legend.")
parser.add_argument("--legendlabelarr", "--leglabarr", "--legarr", type=str, nargs='+', help="Set the labels to be used in legend with an array.")
parser.add_argument("--stylearrall", "--starrall", "--styarr", "--stylearrall", "--staa", type=str, nargs='+', help="Set color-marker-line styles all in array.")
parser.add_argument("--kondotick", "--kondot", action='store_true', help="Set other's tick.")
parser.add_argument("--vline", "--vl", type=float, nargs='+', help="Draw vlines")
parser.add_argument("--hline", "--hl", type=float, nargs='+', help="Draw hlines")
parser.add_argument("--dashline", "--dashedline", "--dashl", action="store_true" , help="Plot as dashed-lines." ) 
parser.add_argument("--markerdiffoffset", "--msdiffoffset", "--msoffset", "--msoff", type=str, help="Set offset to use 'markersdiff'.")
parser.add_argument("--fixmarker", "--fixm", type=str, help="Set marker.") 
parser.add_argument("--emptymarker", "--emarker", "--emptym", action='store_true', help="Use open(empty)-marker.") 
parser.add_argument("--nomarker", "--nmarker", "--nom", action='store_true', help="Remove markers.") 
parser.add_argument("--markersize", "--msize", type=str, help="Set the size of marker.") 
parser.add_argument("--mecarr", "--mecarray", type=str, nargs='+', help="Set mec array")
parser.add_argument("--mfcarr", "--mfcarray", type=str, nargs='+', help="Set mfc array")
parser.add_argument("--csarr", "--colorstylearr", "--colorstylearray", type=str, nargs='+', help="Set cs array")
parser.add_argument("--msarr", "--markerarr", "--markerarray", type=str, nargs='+', help="Set ms array")
parser.add_argument("--nax", type=str, help="Set nax.")
parser.add_argument("--derivative", "--deriv", "--dydx", action='store_true', help="Plot derivatives from the data with midpoint.") 
parser.add_argument("--absy", "--absolutey", "--abs",  action='store_true', help="Plot abs(y)" ) 
parser.add_argument("--absx", "--absolutex", action='store_true', help="Plot abs(x)" ) 
args = parser.parse_args()

fs = 12
if args.fontsize : fs = int(args.fontsize)
figsy=3.5 ; figs = figsy * 2.2
figs  = 6
figsy = 3.7
if args.figs :
	figsarr = args.figs.split("_")
	figs    = float(figsarr[0])
	figsy   = float(figsarr[1])
print "figsize = ({}, {})".format( figs, figsy )
fsleg = fs
if args.fontsizelegend : fsleg = float(args.fontsizelegend)
fstitle = 'small'
if args.fontsizetitle : fstitle = float(args.fontsizetitle)
params = {      'legend.fontsize'       : fsleg ,
                'figure.figsize'        : (figs,figsy) ,
                'axes.labelsize'        : 'x-large' ,
                'axes.titlesize'        : fstitle ,
		'axes.linewidth'        : 2. ,
                'xtick.labelsize'       :'x-large' ,
                'ytick.labelsize'       :'x-large' ,
                'xtick.direction'       : 'in' ,
                'ytick.direction'       : 'in' ,
                'xtick.minor.visible'   : False ,
                'ytick.minor.visible'   : False ,
		'xtick.major.pad'	: 7. +(fs-12.)*0.3 ,      # distance to major tick label in points
		'ytick.major.pad'	: 7. +(fs-12.)*0.3 ,      # distance to major tick label in points
		#'axes.labelpad'		: 4* (fs/12) ,
                'lines.markersize'      : fs*0.7 ,
                'markers.fillstyle'     : 'full' , # full|left|right|bottom|top|none
                'font.size'             : fs }
#print  "PADPAD : ", 7. +(12.-fs)*1.0 
#print  "PADPAD : ", 17. +(12.-fs)*1.0 
#exit(1)
pylab.rcParams.update(params)

AdjustSpaceDefault = [ 0.23 , 0.93 , 0.2 , 0.93 , 0    , 0    ]

nax = int(args.nax) if args.nax else 1
if nax == 1 : 
	fig, ax	= plt.subplots()
	axarr	= [ ax ] * 10
else :
	fig, axarr	= plt.subplots(nax, sharex=True)
	ax	= axarr[0]

if nax==1 : 
	ax.set_xlabel( r'$J_H$ [eV]' )
	ax.set_ylabel( r'$\chi^{M}(\beta/2\rightarrow\infty)$') #at $T=0$' ) # [$\mu_B^2$)' )
for axnow in axarr : 
	axnow.tick_params( axis='both' )#,  pad=8)
#ax.tick_params( axis='x', labelbottom=False )
#plt.legend()


beta	= int( args.beta ) 
IT	= int( args.IT )  if args.IT else 128

ms = "C0s"
faffix = ""
if args.zerosoc :	faffix = "_S0"


op = "S"
if args.operator : op = args.operator 

oplab	= r'$^{}$'.format(op)

##--Get Data (pre)-----
def getobjarr( objclass , itemarr )  :
	return [ getattr(objclass,item) for item in itemarr ]
def ftype(tname) :
	if tname.find("float")>-1 :	return float
	if tname.find("str")>-1 :		return str
	if tname.find("complex")>-1 :	return complex
xobj = "nOcculattAra"	; xobjtype = float
yobj = "J"		; yobjtype = float
if args.plotx :	xobj = args.plotx
if args.ploty :	yobj = args.ploty
if args.typex :	xobjtype = ftype(args.typex)
if args.typey :	yobjtype = ftype(args.typey)
itemarr		= [ xobj, yobj ]
if args.multiplycomponent : 
	itemarr.append(	args.multiplycomponent )
itemarr.append(	'beta' )
print "Chosen : ", itemarr, " ({},{})".format( xobjtype, yobjtype )
axarr[-1].set_xlabel( xobj )
for axnow in axarr :
	axnow.set_ylabel( yobj )
xlab = xobj
ylab = yobj
print "xobj : ", xobj
print "yobj : ", yobj
	
try : 
	xlab = filterlabel( xobj , paper=args.paper, unit=True )
except : pass
try : 
	ylab = filterlabel( yobj , paper=args.paper, unit=True )
	if args.derivative : 
		ylab = r'$\partial$'+ylab + " / " + r"$\partial$"+xlab 
except : pass

if args.multiyobj : 
	itemarr		= [ xobj ]
	yobjarr		= args.multiyobj
	for yo in yobjarr : 
		itemarr.append( yo ) 
	if args.multiplycomponent : 
		itemarr.append(	args.multiplycomponent )
	itemarr.append( 'beta' ) 
	print "Chosen (multi-yobj) : ", itemarr, " ({},{})".format( xobjtype, yobjtype )
	ylab = ''
else : 
	yobjarr = [yobj]

quant = args.plot if args.plot else 'power'


def fLabel( name ) :
	if name=="S" : 
		labname = r'$\lambda_{\mathrm{SOC}}$'
	elif name=="J" : 
		labname = r'$J_H$'
	elif name=="power-exponent" : 
		labname = r'$\alpha$'
	if name.find('inv')>-1 : 
		labname = "("+labname+")"+r'$^{-1}$'
	return labname
def fLabelVal( name, value ) :
	labval = "={}".format(value)
	return fLabel(name) + labval
def autoLabelVal( nameArr, S=None, J=None, T=None ) :
	lab	= ""
	for name in nameArr :
		if name=="S" : 
			if len(lab)>0 : lab += ', '
			lab	+= filterlabel('lambdasoc',paper=args.paper,unit=False)+" = {} eV".format(S)
		elif name=="J" : 
			if len(lab)>0 : lab += ', '
			lab	+= filterlabel('J',paper=args.paper,unit=False)+" = {} eV".format(J)
		elif name=="T" : 
			if len(lab)>0 : lab += ', '
			lab	+= filterlabel('T',paper=args.paper,unit=False)+" = {:.2f} K".format(T)
	return lab

for axnow in axarr : 
	axSetRange( axnow , args.xxrange , "x" )
	axSetRange( axnow , args.yrange  , "y" )

betaarr	= [ 128, 256, 512 ]
if args.setbetaarr :
	betaarr = args.setbetaarr 
	args.betaarr	= True
ITarr	= [ 128, 256, 512 ]
JHarr	= [ 0, 0.2, 0.4, 0.5, 0.7 ]
Sarr	= [ float(args.soc) ] if args.soc else [ 0.1, 0 ]
S	= float(args.soc) if args.soc else 0.1 
J	= float( args.JH ) if args.JH else 0.4
if J==0.4 and IT==0.	: 	betaarr.append( 1024 )
if (J==0.4 or J==0.2 or J==0.7 or J==0.5) and IT>0	: 	ITarr.append( 1024 )
if (J==0.4 ) and IT>0	: 	ITarr.append( 0 )
#if (J==0.4 ) and IT>0	: 	ITarr.append( 2048 )
if args.IT : 
	ITarr	= None

arr	= Sarr
if args.JHarr :		arr = JHarr
elif args.betaarr :	arr = betaarr
elif args.ITarr :	arr = ITarr
print "Arr : ", arr
ims=0
icmls	= 0 
fitvalarr = []
fiterrarr = []
yerr	= None
titlefit = ''
if args.fnamearr : arr = args.fnamearr
for ii in range(len(arr)) :
	axii	= axarr[ii]	if args.subplotarrobj=='dir' else ax 
	csii	= ii		if args.subplotarrobj=='dir' else "C%d"%(ii%10)
	print "---", ii, "------------------------"
	if args.JHarr :
		J = arr[ii]
	elif args.betaarr :
		beta = arr[ii]
	elif args.ITarr :
		IT = arr[ii]
	elif args.fnamearr :
		pass
	else : 
		S = Sarr[ii]
	varSJlab	= fLabelVal('S',S)+'eV' + ", " + fLabelVal('J',J)+'eV'
	varSlab	= fLabelVal('S',S)+'eV'
	varJlab	= fLabelVal('J',J) +'eV'
	#datlab	= "{}".format(varSJlab)
	varSJlab	= ""
	datlab		= ""

	##--Get Data (start)-----
	fnameprefix = args.fprefix if args.fprefix else "all"
	if args.ITarr : 
		fname = args.fname+'_IT{}'.format(IT) if args.fname else fnameprefix+"_S{:.1f}_J{:.1f}_IT{}.dat".format(S,J,IT)
		Treal	= round(11605./IT,1) if IT>0 else 0
		datlab	= varSJlab + ", "+ filterlabel("T")+" = {} K".format( Treal )
	else : 
		fname = args.fname if args.fname else fnameprefix+"_S{:.1f}_J{:.1f}_IT_fit_occu.dat".format(S,J)
	if args.IT and float(args.IT)==0. :
		if args.betaarr : 
			fname	= fnameprefix+"_S{:.1f}_J{:.1f}_IT{}_beta{}.dat".format(S,J,IT,beta)
			datlab	= varSJlab + ", "+ filterlabel("betafit", paper=args.paper )+" = {}".format(beta) + r' eV$^{-1}$'
		else :
			fname	= fnameprefix+"_S{:.1f}_J{:.1f}_IT{}_fit_occu.dat".format(S,J,IT)
			datlab	= varSJlab + ", "+ filterlabel("Treal", paper=args.paper )+" = 0 K"
	elif args.IT and int(args.IT)>0 :
		ITint = int(args.IT) 
		fname	= fnameprefix+"_S{:.1f}_J{:.1f}_IT{}.dat".format(S,J,ITint)
		datlab	= varSJlab + ", "+ filterlabel("Treal", paper=args.paper )+" = "+"{:.2f} K".format(11605./float(ITint))
	if args.fnamearr :
		fname	= arr[ii]
		datlab	= fname
	if args.autolegend : 
		Treal	= round(11605./IT,1) if IT>0 else 0
		datlab = autoLabelVal( args.autolegend , S=S,J=J,T=Treal )
	if args.legendlabel : 
		datlab = r"{}".format( args.legendlabel )
	if args.legendlabelarr : 
		datlab = r"{}".format( args.legendlabelarr[ii] )
	print "lab : ", datlab
	print "Reading : ", fname
	#dat = txtclass() ; dat.getData( fname )

	nyobj	= len(yobjarr)
	if args.multiyobjsplit : nyobj	= 1
	for iy in range(nyobj) : 
		axiy	= axarr[iy]	if args.subplotarrobj=='y' else axii 
		csiy	= csii		if args.subplotarrobj=='y' else "C%d"%(iy%10)
		if args.multiyobjsplit :
			yobj	= yobjarr[ii]
		else	:
			yobj	= yobjarr[iy]
		if args.multiyobj : 
			datlab	= filterlabel(yobj, paper=args.paper )
			if args.derivative : 
				dumxlab = filterlabel( xobj , paper=args.paper, unit=True )
				dumylab = filterlabel( yobj , paper=args.paper, unit=True )
				datlab = r'$\partial$'+dumylab + "/" + r"$\partial$"+dumxlab 
			datlab2	= ""
			#if IT==0 : 
			#	datlab2	= " (" + autoLabelVal( ["S"] , S=S,J=J,T=0 ) + ")"
			#else :
			#	datlab2	= " (" + autoLabelVal( ["S"] , S=S,J=J,T=11605./IT ) + ")"
			if args.legendlabel : 
				datlab2 = r" {}".format( args.legendlabel )
			if args.legendlabelarr : 
				datlab2 = r" {}".format( args.legendlabelarr[ii] )
			datlab	= datlab + datlab2
			print "datlab : ", datlab

		dat = np.genfromtxt( fname, names=True, skip_header=0, usecols=(xobj,yobj) , dtype=float )
		xdat	= dat[xobj]
		ydat	= dat[yobj]
		print "xdat : ", xdat
		print "ydat : ", ydat
		if args.derivative : 
			print "Taking derivatives."
			xmid	= ( xdat[:-1] + xdat[1:] ) 
			xdiff	= ( xdat[:-1] - xdat[1:] ) 
			ydiff	= ( ydat[:-1] - ydat[1:] ) 
			xdat	= xmid/2.
			ydat	= ydiff/xdiff
			print "xdat.derivative : ", xdat
			print "ydat.derivative : ", ydat
		if args.plotyerr : 
			yerrobj	= args.plotyerr[iy]
			dat = np.genfromtxt( fname, names=True, skip_header=0, usecols=yerrobj , dtype=float )
			yerr	= dat[yerrobj]
			print "yerr : ", yerr

		msoff	= int(args.markerdiffoffset) if args.markerdiffoffset else 0 
		ms = markersdiff[ims+msoff]
		cs = csiy if args.subplotarrobj else 'C%d'%(ims%10)
		ims+=1
		ls = '-'
		me = 'k'
		mf = None
		msize = fs*0.7
		if args.soc :
			pass
		else : 
			if len(yobjarr)>1 : 
				ms = markersdiff[iy+msoff]
				if S==0 : 
					cs = 'C1'
				else : 
					cs = 'C0'
			if (yobj=="chitotszszStaticwn0zz") or (yobj=="chitotlzlzStaticwn0zz") :
				me	= cs
				mf	= 'w' 
				ls	= ':'
		if ms == 'D' : msize = fs*0.5
		if args.solidall : ls='-'
		dashpt = (None, None)
		if args.dashline : dashpt = dashdiff(ims-1)
		if args.interpolate :  ls =''
		print "cs, ms, ls : ", cs, ms, ls, dashpt
		if args.nomarker :
			ms =''
		if args.fixmarker : 
			ms = args.fixmarker 
		if args.emptymarker : 
			me = cs
			mf = 'w'
		if args.markersize : 
			msize = float( args.markersize ) 
		if args.mecarr : 
			me	= args.mecarr[ii]
		if args.mfcarr : 
			mf	= args.mfcarr[ii]
		if args.csarr :
			cs	= args.csarr[ii]
		if args.msarr :
			print "msarr : ", args.msarr 
			ms	= args.msarr[ii]

		cmlstyle	= cs+ms+ls
		if args.stylearrall :
			print ("stylearrall length : ", len(args.stylearrall))
			cmlstyle	= args.stylearrall[icmls]
			icmls		+= 1
			if cmlstyle.find(":")>-1 :
				dashpt = (1,1)
			elif cmlstyle.find("-.")>-1 :
				dashpt = (3,1,1,1)
			elif cmlstyle.find("--")>-1 :
				dashpt = (3,1)
			elif cmlstyle.find("-")>-1 :
				dashpt = (1,0)
		print "cmlstyle, dashpt : ", cmlstyle, dashpt, type(cmlstyle)
		print "cs, ms, ls : ", cs, ms, ls, dashpt
		if args.inverse :
			ydat = 1./ydat
		if args.absy :
			ydat = np.abs(ydat)
		if args.absx :
			xdat = np.abs(xdat)
		if args.logy : 
			if args.mergeparticlehole : 
				nparticlehole	= int(args.nparticlehole) if args.nparticlehole else 6
				xdatmerge	= np.concatenate( (xdat, nparticlehole-xdat) )
				ydatmerge	= np.concatenate( (ydat, ydat) )
				xdatmergeargsort	= xdatmerge.argsort()
				xdatmerge	= xdatmerge[ xdatmergeargsort ]
				ydatmerge	= ydatmerge[ xdatmergeargsort ]
				axiy.semilogy( xdatmerge , ydatmerge , color=cs, marker=ms, linestyle=ls, mfc=mf, mec=me , label=datlab, clip_on=args.clipon , markersize=msize, dashes=dashpt ) # markersize=msize #markersize=fs*0.3 )
			elif args.particleholeonly : 
				pass
			else : 
				axiy.semilogy( xdat , ydat , color=cs, marker=ms, linestyle=ls, mfc=mf, mec=me , label=datlab, clip_on=args.clipon , markersize=msize, dashes=dashpt ) # markersize=msize #markersize=fs*0.3 )
			if args.particlehole : 
				if args.stylearrall :
					cmlstyle	= args.stylearrall[icmls]
					icmls		+= 1
			if args.particlehole or args.particleholeonly : 
				print ("Particle-hole partner plotted.")
				nparticlehole	= int(args.nparticlehole) if args.nparticlehole else 6
				axiy.semilogy( nparticlehole-xdat , ydat , color=cs, marker=ms, linestyle=ls, mfc=mf, mec=me , label="", clip_on=args.clipon , markersize=msize, dashes=dashpt ) # markersize=msize #markersize=fs*0.3 )
		elif args.logy : 
			axiy.semilogx( xdat , ydat , color=cs, marker=ms, linestyle=ls, mfc=mf, mec=me , label=datlab, clip_on=args.clipon , markersize=msize, dashes=dashpt ) # markersize=msize #markersize=fs*0.3 )
		elif args.loglog : 
			axiy.loglog( xdat , ydat , color=cs, marker=ms, linestyle=ls, mfc=mf, mec=me , label=datlab, clip_on=args.clipon , markersize=msize, dashes=dashpt ) # markersize=msize #markersize=fs*0.3 )
		else : 
			if args.mergeparticlehole : 
				nparticlehole	= int(args.nparticlehole) if args.nparticlehole else 6
				xdatmerge	= np.concatenate( (xdat, nparticlehole-xdat) )
				ydatmerge	= np.concatenate( (ydat, ydat) )
				xdatmergeargsort	= xdatmerge.argsort()
				#print "xdat_merge : ", xdatmerge
				#print "xdat_merge argsort : ", xdatmergeargsort
				xdatmerge	= xdatmerge[ xdatmergeargsort ]
				ydatmerge	= ydatmerge[ xdatmergeargsort ]
				#print "xdat_merge[argsort] : ", xdatmerge
				axiy.plot( xdatmerge , ydatmerge , color=cs, marker=ms, linestyle=ls, mfc=mf, mec=me , label=datlab, clip_on=args.clipon , markersize=msize, dashes=dashpt ) # markersize=msize #markersize=fs*0.3 )
			elif args.particleholeonly : 
				pass
			else :
				if yerr is not None : 
					axiy.errorbar( xdat , ydat , fmt='none', yerr=yerr , ecolor='k', capsize=4, capthick=1  )
				axiy.plot( xdat , ydat , color=cs, marker=ms, linestyle=ls, mfc=mf, mec=me , label=datlab, clip_on=args.clipon , markersize=msize, dashes=dashpt ) # markersize=msize #markersize=fs*0.3 )

			if args.particlehole :
				if args.stylearrall :
					cmlstyle	= args.stylearrall[icmls]
					icmls		+= 1
			if args.particlehole or args.particleholeonly : 
				print ("Particle-hole partner plotted.")
				nparticlehole	= int(args.nparticlehole) if args.nparticlehole else 6
				axiy.plot( nparticlehole-xdat , ydat , color=cs, marker=ms, linestyle=ls, mfc=mf, mec=me , label="", clip_on=args.clipon , markersize=msize, dashes=dashpt ) # markersize=msize #markersize=fs*0.3 )

		if args.interpolate : 
			from scipy.interpolate import interp1d
			ipmethod = args.interpolatemethod if args.interpolatemethod else 'linear'
			nip = int(args.ninterpolate) if args.ninterpolate else 10

			xipmin  = xdat[0]
			xipmax  = xdat[-1]
			if args.interpolaterange :
				xipmin, xipmax = np.array( args.interpolaterange )
			xip     = np.linspace( xipmin, xipmax,   nip)
			finter = interp1d(xdat, ydat, kind=ipmethod )
			
			for ix in range(len(xip)) : 
				ftag	= "{:.3f}".format(xip[ix])+"_{}_{}".format(xobj,yobj)
				if ii<1 : 
					with open("interpolate_{}.out".format(ftag), 'w') as fw :
						fw.write( "#Treal "+xobj+" "+yobj+"\n" )
				with open("interpolate_{}.out".format(ftag), 'a') as fw :
					fw.write( "{} {} {}\n".format(Treal,xip[ix],finter(xip)[ix]) )
			print "xdat-ip : "
			print xip
			print "ydat-ip : "
			print finter(xip)  
			axiy.plot( xip , finter(xip)  , "r:", clip_on=args.clipon, ) # markersize=fs*0.3 )
		if args.fitpower : 
			from scipy.optimize import curve_fit

			nfit = int(args.nfitdat) if args.nfitdat else len(xdat)
			xraw = xdat[-nfit:]
			yraw = ydat[-nfit:]
			print "nfit : ", nfit 
			print "xraw : ", xraw
			print "yraw : ", yraw

			formfunc= 'a x^b'
			def func(x, a, b) :
				return a * np.power(x,b) 

			initvararr = [-0 , -0 ]
			initvararr2= [ 20,  2 ]
			popt, pcov = curve_fit(func, xraw, yraw, bounds=(initvararr,initvararr2) )
			print 'fit({})'.format(formfunc) +' : a=%5.3f, b=%5.3f' % tuple(popt)
			y0 = popt[-1]

			fitvalarr.append(y0)
			axisval = axiy.axis()
			xdatfit = np.linspace(xdat[0],xdat[-1],20)
			axiy.plot(xdatfit, func(xdatfit, *popt), cs+'--', label=r'${}$'.format(formfunc)+'\n'+'fit: a=%5.3f, b=%5.3f' % tuple(popt))
		print "------------------------------"

		if args.fitinverse : 
			from scipy.optimize import curve_fit
			def func(x, a, b, c):
				return 1./(a*x+b)+c
			if yerr is not None : 
				ystdev	= yerr
				popt, pcov = curve_fit(func, xdat, ydat, bounds=(0, [5., 5., 10]), sigma=ystdev, absolute_sigma=True )
			else : 
				popt, pcov = curve_fit(func, xdat, ydat, bounds=(0, [5., 5., 10]) )
			perr	= np.sqrt(np.diag(pcov))
			print 'fit (1/(ax+b) + c) : a=%9.6f, b=%9.6f, c=%9.6f' % tuple(popt) , ', 1/a={:9.6f}, b/a={:9.6f}'.format( 1./popt[0], popt[1]/popt[0] )
			a,b,c	= tuple(popt)
			y0	= 1./b+c
			errbinv	= abs(1./b-1./(b+perr[1])) + abs(1./b-1./(b-perr[1]))
			errbinv	= errbinv/2.
			errc	= abs(perr[2])
			print "perr : ", perr
			print "errbinv errc : ", errbinv, errc
			fitvalarr.append(y0)
			fiterrarr.append(errbinv+errc)
			axisval = axiy.axis()
			xdatfit = np.linspace(axisval[0],10,60)
			xdatfit = np.append( xdatfit, np.linspace(10,axisval[1],60) ) 
			labfit	= ''
			if args.fitfunction :
				labfit	= r'$\frac{1}{ax+b}+c$'
			if args.fitparam :
				if len(labfit)>0 : labfit = labfit +'\n'
				labfit	= labfit + 'fit: a=%5.3f, b=%5.3f, c=%5.3f, b/a=%5.3f' % tuple([a,b,c,b/a]) 
			if args.fitparamtitle : 
				if len(titlefit)>0 : titlefit = titlefit +'\n'
				titlefit = titlefit + 'fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt)
			axiy.plot(xdatfit, func(xdatfit, *popt), cs+'--', label=labfit )
		if args.fitlinear : 
			from scipy.optimize import curve_fit
			def func(x, a, b) :
				return a*x + b
			popt, pcov = curve_fit(func, xdat, ydat)
			print 'fit: a=%5.3f, b=%5.3f' % tuple(popt)
			y0 = popt[-1]
			if args.inverse :
				y0 = 1./y0
				print "1/y : ", y0
			fitvalarr.append(y0)
			axisval = axiy.axis()
			xdatfit = np.linspace(axisval[0],axisval[1],20)
			axiy.plot(xdatfit, func(xdatfit, *popt), cs+'--', label=r'$ax+b$'+'\n'+'fit: a=%5.3f, b=%5.3f' % tuple(popt))
		if args.fitquadratic : 
			from scipy.optimize import curve_fit
			def func(x, a, b, c):
				return a*x*x + b*x + c
			popt, pcov = curve_fit(func, xdat, ydat)
			perr	= np.sqrt(np.diag(pcov))
			print 'fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt)
			y0 = popt[-1]
			if args.inverse :
				y0 = 1./y0
				print "1/y : ", y0
			fitvalarr.append(y0)
			axisval = axiy.axis()
			xdatfit = np.linspace(axisval[0],axisval[1],20)
			labfit	= ''
			if args.fitfunction : 
				labfit	= r'$ax^2+bx+c$'
			if args.fitparam : 
				if len(labfit)>0 : labfit = labfit +'\n'
				labfit	= labfit + 'fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt)
			if args.fitparamtitle : 
				if len(titlefit)>0 : titlefit = titlefit +'\n'
				titlefit = titlefit + 'fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt)
			axiy.plot(xdatfit, func(xdatfit, *popt), cs+'--', label=labfit ) 


if args.fitlinear or args.fitquadratic or args.fitinverse or args.extrapolate or args.fitpower: 
	print "fit-value, err:"
	for ival in range(len(fitvalarr)) : 
		try :
			print fitvalarr[ival], fiterrarr[ival]
		except :
			print fitvalarr[ival]
	#for val in fitvalarr : 
	#	print val 
	print ""

if args.inverse : 
	ylab = "("+ylab+")"+r'$^{-1}$'

if args.kondotick :
	ax.set_xticks( range(0,301,100) )
	ax.set_yticks( range(0,10,3) )

if args.xticklabel : 
	xt =  np.array(args.xticklabel)
	for axnow in axarr : 
		axnow.set_xticks( xt ) 
	print "xticklab : ", xt
if args.yticklabel : 
	yt =  np.array(args.yticklabel)
	for axnow in axarr : 
		axnow.set_yticks( yt ) 
	print "yticklab : ", yt

if args.noxtick :
	for axnow in axarr : 
		axnow.tick_params( axis="x", left=True,  bottom=False,  labelbottom=False, labelleft=False  )
if args.noytick :
	for axnow in axarr : 
		axnow.tick_params( axis="y", left=False, bottom=True,   labelbottom=False, labelleft=False  )
if args.ylabelcomp : 
	ylab = filterlabel( args.ylabelcomp , paper=args.paper, unit=True )
if args.xlabelcomp : 
	xlab = filterlabel( args.xlabelcomp , paper=args.paper, unit=True )
if args.noxlabel : 
	xlab = ''
if args.noylabel : 
	ylab = ''
for axnow in axarr : 
	axnow.set_ylabel( ylab )
axarr[-1].set_xlabel( xlab )

if args.ylabelcoord :
	ylabcoordx = float( args.ylabelcoord )
	for axnow in axarr : 
		axnow.yaxis.set_label_coords( ylabcoordx, 0.5)

if args.ncolumnlegend : 
	ncol = args.ncolumnlegend 
else :
	ncol = 1
if args.nolegend : 
	pass
else : 
	for axnow in axarr : 
		axnow.legend( ncol=ncol ) #, fontsize=20 )
		if args.legendlocation : 
			try : 
				legloc = int(args.legendlocation)
			except :
				legloc = args.legendlocation
			if args.legendlocation == "outside" : 
				axnow.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., ncol=ncol)
			else : 
				axnow.legend( loc=legloc, ncol=ncol ) #, fontsize=20)
		
		#axnow.legend(loc=4)
#legend_elements = [
#			Line2D([0], [0], color='C0', marker='s', mec='k', lw=2, label=filterlabel('chiSStaticTotalzz')+'  '+ r'($\lambda_\mathrm{SOC}=0.1$) '	),
#			Line2D([0], [0], color='C1', marker='s', mec='k', lw=2, label=filterlabel('chiSStaticTotalzz')+'  '+ r'($\lambda_\mathrm{SOC}=0$) '		),
#			Line2D([0], [0], color='C0', marker='v', mfc='w', lw=2, linestyle=':', label=filterlabel('chitotszszStaticwn0zz')+'  '+r'($\lambda_\mathrm{SOC}=0.1$) ' ),
#			Line2D([0], [0], color='C1', marker='v', mfc='w', lw=2, linestyle=':', label=filterlabel('chitotszszStaticwn0zz')+'  '+r'($\lambda_\mathrm{SOC}=0$) '	 ),
#			]
#ax.legend(handles=legend_elements)

if args.sublabel :
	for iax in range(len(axarr)) : 
		axnow	= axarr[iax]
		x0,x1,y0,y1 = axnow.axis()
		#axnow = axarr[dd]
		dx=x1-x0;dy=y1-y0#; axnow.text( dx*0.02, dy*0.98 , markerssub[dd] ) 
		slab = markerssub[iax]
		if args.sublabelnoalphabet : 
			slab = ''
		if args.sublabelarr : 
			slab = slab + args.sublabelarr[iax] 
		#if args.sublabelparameter : 
		#	for param in args.sublabelparameter : 
		#		if param=="nOcculattAraFormer" : 
       		#			slab = slab +" {}={:.2f}".format( filterlabel(param) , float(getattr(hobj,param)) )
		#		elif param.find("tfilling")>-1 : 
       		#			slab = slab +" {}={:.2f}".format( filterlabel(param) , float(getattr(hobj,param)) )
		#		else : 
       		#			slab = slab +" {}={:.6g}".format( filterlabel(param) , float(getattr(hobj,param)) )
        	axnow.text(0.02, 0.76, slab,   
				horizontalalignment='left', verticalalignment='center', transform=axnow.transAxes,
        	                fontsize=fs+4 )# , weight='semibold')

if args.fitparamtitle : 
	plt.title( titlefit ) 

for axnow in axarr :
	axnow.grid(linestyle=":")

axis = ax.axis()
print "axis : ", axis

if args.hline : 
	for hVal in args.hline : 
		for axnow in axarr : 
			axnow.axhline( hVal , 0, 1 ) 
if args.vline : 
	for vVal in args.vline : 
		for axnow in axarr : 
			axnow.axvline( vVal , 0, 1 ) 

#ax.axhline( 3.87, 0, 1 ) 
#ax.axhline( 3.94, 0, 1 ) 
#ax.axhline( 5.5, 0, 1 ) 
#ax.axhline( 3.5, 0, 1 ) 
#ax.axhline( 3.0, 0, 1 ) 

ladjust=( AdjustSpaceDefault[0] if (args.ladjust is None) else float(args.ladjust) )
radjust=( AdjustSpaceDefault[1] if (args.radjust is None) else float(args.radjust) )
badjust=( AdjustSpaceDefault[2] if (args.badjust is None) else float(args.badjust) )
tadjust=( AdjustSpaceDefault[3] if (args.tadjust is None) else float(args.tadjust) )
wspace =( AdjustSpaceDefault[4] if (args.wspace  is None) else float(args.wspace ) )
hspace =( AdjustSpaceDefault[5] if (args.hspace  is None) else float(args.hspace ) )
print "lad, bad, rad, tad, wsp, hsp : " , ladjust, badjust, radjust, tadjust, wspace, hspace
plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, wspace=wspace , hspace=hspace )

transparent = True
if args.savetransparent : transparent = True
if args.savenottransparent : transparent = False
fsave = "plot_{}{}".format(xobj,yobj)
if args.filelabel : 
	fsave += args.filelabel 
ndpi=200
if args.ndpi : ndpi	= int(args.ndpi)
#plt.savefig( fsave+'.pdf', dpi=ndpi , transparent=transparent ) 
if args.cptmp :
	dtype = "pdf"
	if args.pdf :
		dtype="pdf"
	if args.png :
		dtype="png"
	if args.eps :
		dtype="eps"
	print "Saved : ", fsave + "."+dtype
	plt.savefig( fsave + "."+dtype , dpi=ndpi, transparent=transparent) 
	cptmpftn( fsave, cmd="scp", destMac=False , dtype=dtype )
else : 
	plt.show()
