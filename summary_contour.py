#!/opt/python/2.7.15/gcc-4.8.5/bin/python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys, os

import argparse

from basislabel import *
from hjl_common import *

parser = argparse.ArgumentParser(description="Plot ldos with j-basis")
parser.add_argument('ddir', metavar='D', type=str, nargs='+',
		                    help='Data paths' ) 
parser.add_argument("-y", "--yrange", type=str, help="set yrange e.g. _y0_y1" )
parser.add_argument("-x", "--xxrange", type=str, help="set xrange e.g. _x0_x1" )
parser.add_argument("-z", "--zrange", type=str, help="set zrange e.g. _x0_x1" )
parser.add_argument("-dy", "--datyrange", type=str, help="set datyrange e.g. _y0_y1" )
parser.add_argument("-dx", "--datxrange", type=str, help="set datrange e.g. _x0_x1" )
parser.add_argument("-b", "--basis", type=str, help="Set basis. e.g. \"-b j | -b t2g\"" ) 
parser.add_argument("-c", "--chdat", type=str, help="Choose a component(s) of basis. e.g. \"-c _0_1_2 (default 0,2,4)\"" ) 
parser.add_argument("--chdatoffd", "--co", "--chdatoff", type=str, help="Choose a component(s) of basis. e.g. \"-c _0_1_2_4 (giving (0,1) and (2,4) components.)\"" ) 
parser.add_argument("-u", "--ups", action="store_true" , help="Plot (pseudo)spin-up components of basis" )
parser.add_argument("--pdf", action="store_true", help="Save the figure in .pdf" ) 
parser.add_argument("--png", action="store_true", help="Save the figure in .png" ) 
parser.add_argument("--ndpi", type=str, help="Set ndpi, number of dots per inch." ) 
parser.add_argument("--savetransparent", "--savetp", "--stp", action="store_true", help="Save the fig with transparent background" )
parser.add_argument("--cptmp", action="store_true", help="Copy the saved figure into Dropbox/tmp_linux/") 
parser.add_argument("--notitle", action="store_true", help="Remove the title caption." )
parser.add_argument("--titlevar", "--tvar", "--tv", action="store_true", help="Set the title as name of variable.")
parser.add_argument("--squre", action="store_true", help="Plot in square frame.")
parser.add_argument("-s", "--separate", action="store_true" , help="Plot data in different frame for each orbital." ) 
parser.add_argument("-w", "--w0", action="store_true" , help="Plot data near 0 ." ) 
parser.add_argument("-l", "--llog", action="store_true" , help="Plot data in log-log plot." ) 
parser.add_argument("-m", "--minus", action="store_true" , help="Plot self-energy with the original minus sign." ) 
parser.add_argument("--minusx", "--mx", action="store_true" , help="Plot minus x-data." ) 
parser.add_argument("--selfm", "--sm", action="store_true", help="Plot minus self-energy" ) 
parser.add_argument("--selffactor", "--sf", type=str, help="Multiply the self-energy by a factor." ) 
parser.add_argument("--selfoffset", "--so", type=str, help="Add a offset to the self-energy.")
parser.add_argument("--plusself", action="store_true" , help="Plot self-energy with the original minus sign." ) 
parser.add_argument("--count", type=int, help="Choose the count of which data you want." ) 
parser.add_argument("--ppt", action="store_true" , help="Plot as ppt-use." ) 
parser.add_argument("--paper", action="store_true" , help="Plot as paper-use." ) 
parser.add_argument("--noleg", action="store_true" , help="Remove legend")
parser.add_argument("--alllegend", "--allleg", "--aleg", action="store_true" , help="Draw all legends.")
parser.add_argument("--tl", action="store_true" , help="Tight-layout") 
parser.add_argument("--label", type=str, help="Set label" ) 
parser.add_argument("--leglabel", "--legl", type=str, help="Set the labels of the legend" ) 
parser.add_argument("--sublabel", "--slab", action="store_true", help="Set the labels of sub-figures with alphabet." ) 
parser.add_argument("--sublabeloff", "--slaboff", type=str, help="Set the offset of the labels. e.g. '2' will begin with (c) for the first sub-label.")
parser.add_argument("--sublabelparameter", "--slabpar", "--slabp", type=str, help="Set a parameter of the sub-labels.")
parser.add_argument("--setxlabel", "--setxlab", "--xlab", type=str, help="Set x-label.")
parser.add_argument("--setylabel", "--setylab", "--ylab", type=str, help="Set y-label.")
parser.add_argument("--noxticklabel", "--notickxlab", "--noxtlab", action="store_true", help="Unset x-ticklabel.")
parser.add_argument("--noyticklabel", "--notickylab", "--noytlab", action="store_true", help="Unset y-ticklabel.")
parser.add_argument("--noylabel", "--noylab", "--noyl", action="store_true", help="Remove the y-axis label.")
parser.add_argument("--noxlabel", "--noxlab", "--noxl", action="store_true", help="Remove the x-axis label.")
parser.add_argument("--xticklabel", "--xticklab", "--xtlab", type=float, nargs='+', help="Specify x-ticklabel.")
parser.add_argument("--yticklabel", "--yticklab", "--ytlab", type=float, nargs='+', help="Specify y-ticklabel.")
parser.add_argument("--ylabelcoord", "--ylabcoord", "--ylabc", "--ylc", type=str, help="Set the x-coordinate of y-axis label.")
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
parser.add_argument("--wspace",  "--wsp", type=str, help="Set wspace-adjust value." ) 
parser.add_argument("--hspace",  "--hsp", type=str, help="Set hspace-adjust value." ) 
parser.add_argument("--rowsep", "--rowseparation", action="store_true" , help="Separate the subplots in row." )
parser.add_argument("--legendlocation", "--legloc", "--legendloc", type=str, help="Specify location index of legend." ) 
parser.add_argument("--legendlabelspace", "--legls", "--legendlsp", type=str, help="Specify the value of 'labelspacing' of legend." ) 
parser.add_argument("--jnat", "--jn", action="store_true" , help="Plot the natural basis from t2g's." )
parser.add_argument("--fontsize", "--fonts", "--fs", type=str, help="Set the font-size.") 
parser.add_argument("--linearline", "--lline", "--ll", type=str, help="Set m^* and Plot a linear line of a*omega (e.g. '--ll 2' : m^*=2, a=1-m^*=-1)." ) 
parser.add_argument("--fractionline", "--fline", "--fl", type=str, help="Set m^* and Plot a line of a*omega^b (e.g. -'-ll m_b' : m^*=2, a=1-m^*=-1, b=0.5)." ) 
parser.add_argument("--vline", "--vl", type=str, help="Plot a vertical line with on specifited x-value." )
parser.add_argument("--rotatelab", "--rl", "--rotl", action="store_true" , help="Rotate the tick labels.")
parser.add_argument("--backgroundcolor", "--bgc", action='store_true', help="Change background color." )
parser.add_argument("--offdiag", "--offd", action='store_true', help="Plot off-diagonal part.")
parser.add_argument("--rew", "--realfreq", action='store_true', help="Plot the data of real-frequency.")
parser.add_argument("--dashedline", "--dashl", action="store_true", help="Plot the data as dashedlines.")
parser.add_argument("--epsilon", "--ep", type=str, help="Specify epsilon of data. (default=0.03)")
parser.add_argument("--sinv", "--selfinv", action="store_true", help="Plot inverse of self-energies.") 
parser.add_argument("--solidline", "--solidl", action="store_true", help="Plot the data as solidlines.")
parser.add_argument("--trans", "--tr", action="store_true", help="Transpose the plot" ) 
parser.add_argument("--fillb", "--fb", action="store_true", help="Fill the plotting of blue one." ) 
parser.add_argument("--fill", action="store_true", help="Fill the plotting of data." ) 
parser.add_argument("--selfyrange", "--sy", type=str, help="set yrange of self-energy e.g. _y0_y1" )
parser.add_argument("--selfxrange", "--sx", type=str, help="set xrange of self-energy e.g. _x0_x1" )
parser.add_argument("--dxtick", action="store_true", help="Double the xtick increment." ) 
parser.add_argument("--dytick", action="store_true", help="Double the ytick increment." ) 
parser.add_argument("--ztick", "--zt", type=str, help="Set ztick as 'z0_z1_z2'.")
parser.add_argument("--transform", "--transf", type=str, help="Transform the matrix data." )
parser.add_argument("--eigenval", "--eval", "--ev", action="store_true", help="Plot eigenvalues.")
parser.add_argument("--addSOCt2g", "--asoc", "--addsoc", type=str, help="Add the SOC to the matrix-data of t2g.")
parser.add_argument("--addSOCjeff", "--asocj", "--addsocj", type=str, help="Add the SOC to the matrix-data of t2g.")
parser.add_argument("--zerosoc", "--zsoc", action="store_true", help="Turn on the zeroSOC mode (that will change some indices).")
parser.add_argument("--realpart", "--realdata", "--rp", "--rd", action="store_true", help="Plot the real-part of the data instead of the imaginary part.")
parser.add_argument("--write", "--wr", action="store_true", help="Write the last part of data in 'tmpxdat' and 'tmpydat'.")
parser.add_argument("--writegrid", "--wrg", type=str, help="Write the data with max grid.")
parser.add_argument("--linearfit", "--lfit", "--lf", action="store_true", help="Fit the data to a linear-line.")
parser.add_argument("--quadraticfit", "--qfit", "--qf", action="store_true", help="Fit the data to a quadratic-line.")
parser.add_argument("--zerotempquadraticfit", "--zqfit", "--zqf", action="store_true", help="Fit the data to a quadratic-line with T=0.")
parser.add_argument("--fitting", "--fit", type=str, help="Set the type of fitting and do.")
parser.add_argument("--ndatalinearfit", "--nlfit", "--nlf", type=str, help="Specify the number of data to be used in fitting (default=3).")
parser.add_argument("--linearfitfromzero", "--lfitz", "--lfz", action="store_true", help="Do linear-fitting assuming it has (0,0).")
parser.add_argument("--fitoffset", "--fitoff", "--foff", type=str, help="Set the offset in indices of the data.")
parser.add_argument("--xlnxfit", "--xlnx", action="store_true", help="Fit the data to 'x ln(x)'.")
parser.add_argument("--xlnxfitconst", "--xlnxc", type=str, help="Set 'c' of 'x ln(x/c)'.")
parser.add_argument("--removelinear", "--rmlin", "--rlin", action="store_true", help="Fit linear-line and remove the linear-part in the data .")
parser.add_argument("--setlinear", "--setlin", "--slin", type=str, help="Set linear-line steepness and remove the linear-part in the data .")
parser.add_argument("--drawabc", "--dabc", "--abc", type=str, help="Draw a 2nd-order polynomial a*x*x+b*x+c with a_b_c.")
parser.add_argument("--carray", "--carr", type=str, help="Set a constant term of a 2nd-order polynomial a*x*x+b*x+c.")
parser.add_argument("--finegrid", "--fineg", "--fg", action="store_true", help="Plot the data of fine-grid's.")
parser.add_argument("--finegrid2", "--fineg2", "--fg2", action="store_true", help="Plot the data of fine-grid's.")
parser.add_argument("--finegrid3", "--fineg3", "--fg3", action="store_true", help="Plot the data of fine-grid's according to 'beta_real'.")
parser.add_argument("--finegrid4", "--fineg4", "--fg4", action="store_true", help="Plot the data of fine-grid's from 'self4.c'.")
parser.add_argument("--finegridarray", "--fgarr", "--fga", type=str, help="Set the type of finegrid of each data as a array. e.g. --fga fg_fg2_x_fg3")
parser.add_argument("--transformarray", "--trarr", "--tra", type=str, help="Set the type of transform of each data as a array. e.g. --fga t_t_x_j")
parser.add_argument("--markersize", "--msize", type=str, help="Set the marker-size.") 
parser.add_argument("--cutoff", "--coff", type=str, help="Set lower-bound cutoff of x-data to be plotted .") 
parser.add_argument("--lowestmatsu", "--lmatsu", "--wn0", action='store_true', help="Draw a line of lowest-Matsubara frequency pi/beta.")
parser.add_argument("--wnbetaline", "--wbline", "--wbl", action='store_true', help="Draw lines of pi/beta within [1024,512,256,128,64].")
parser.add_argument("--wnbetalinefit", "--wblinefit", "--wblfit", action='store_true', help="Draw lines of pi/beta within beta_fit.")
parser.add_argument("--difftwopoints", "--difftwo", "--difft", "--dtwo", action='store_true', help="Draw the two-point average slope.")
parser.add_argument("--notshow", "--ns", action='store_true', help="Don't plot on the screen.")
parser.add_argument("--nplot", type=str, help="Specify 'Nplot' in DOS-plotting")
parser.add_argument("--diffmarker", "--diffm", action='store_true', help="Use different markers.")
parser.add_argument("--indimp", "--iimp", type=str, help="Specify the index of impurity atom.")
parser.add_argument("--plotz", "--pz", type=str, help="Specify an object.")
parser.add_argument("--plotx", "--px", type=str, help="Specify an object for x-axis.")
parser.add_argument("--ploty", "--py", type=str, help="Specify an object for y-axis.")
parser.add_argument("--typex", "--tx", type=str, help="Specify the type of object for x-axis.")
parser.add_argument("--typey", "--ty", type=str, help="Specify the type of object for y-axis.")
parser.add_argument("--nou", "--noulab", action='store_true', help="Exclude 'U' in labels.")
parser.add_argument("--colormap", "--cm",  type=str, help="Specify the color-scheme, e.g. jet of plt.cm.jet." ) 
parser.add_argument("--zmax", "--zm",  type=str, help="Specify the maximum value of z." ) 
parser.add_argument("--markdegen", "--md",  action='store_true', help="Mark the point of degen=2.")
parser.add_argument("--nocontourline", "--nocline", "--noc",  action='store_true', help="Remove the black lines on contour.")
parser.add_argument("--tflline", "--tfl", action='store_true', help="Draw T_FL line.")
parser.add_argument("--setattrobj", "--seto", type=str, nargs='+', help="Set attributes of n-th's object. e.g. --seto <n> nOcculattAra 3.99 <n'> SzdegMateval 0")
parser.add_argument("--ncolorlevels", "--ncolor","--ncolors",  "--nlevels", "--nlv", "--nlev", type=int, help="Set the number of colors for colorbar.")
parser.add_argument("--customcbar", "--ccbar",  action='store_true', help="Customize the colorbar.")
parser.add_argument("--nocbar", "--ncbar", "--nocb", action='store_true', help="Remove the colorbar.")
parser.add_argument("--cbarlabel", "--cbarlab", "--colorbarlab", "--cblab", type=str,  help="Set label of color-bar.")
parser.add_argument("--cbartick", "--colorbartick", "--cbtick", nargs='+', type=float,  help="Remove the colorbar.")
parser.add_argument("--cbarpos", "--colorbarposition", "--cbpos", nargs=4, type=float,  help="Specify the position of the colorbar.")
parser.add_argument("--clipon", "--clip", action='store_true', help="Remove data-markers outside the frame.", default=False)
parser.add_argument("--interpolate1d", "--interp1e", "--ip1d", action='store_true', help="Interpolate data.")
parser.add_argument("--interpolateaxis", "--interpax", "--iax", "--ipax", "--ipaxis", type=str, help="Set the direction of axis in interpolation. (default='y')")
parser.add_argument("--interpolaterange", "--interrange", "--irange", "--ipr", type=float, nargs='+', help="Set the range of values you will pick.")
parser.add_argument("--interpolatetranspose", "--interpt", "--ipt", "--iptr", action='store_true', help="Transpose data after interpolation.")
parser.add_argument("--interpolatechoosex", "--ipchx", "--ipcx", "--ipchoosex", type=float, nargs='+',  help="Choose x values after interpolating.")
parser.add_argument("--ninterpolate", "--ninterp", "--nip", type=str, help="Set the grid number of interpolation.")
parser.add_argument("--interpolatecontour", "--interpcont", "--ipcont", "--ipc", action='store_true', help="Interpolate data.")
parser.add_argument("--nomarker", "--nomark", "--nom", action='store_true', help="Remove marekers of the data.")
parser.add_argument("--extend", "--ext", type=str, help="Set the extend regime. e.g. max min neither both" )
args = parser.parse_args()

import matplotlib.pylab as pylab
fs = 12 
if args.fontsize : fs = int( args.fontsize ) 
figsy=5 ; figs = figsy * 1.2
if args.figs :
	figsarr = args.figs.split("_")
	figs    = float(figsarr[0])
	figsy   = float(figsarr[1])
print "figsize = ({}, {})".format( figs, figsy ) 
if args.markersize :	msize = float(args.markersize)
else :			msize = fs*0.7
params = {      'legend.fontsize'       : fs ,
		'figure.figsize'        : (figs,figsy) ,
		'axes.labelsize'        : 'x-large' ,
		'axes.titlesize'        :'small' ,
		'axes.linewidth'        : 2 ,
		'xtick.major.pad'	: 7. +(12.-fs)*0.7 ,      # distance to major tick label in points
		'ytick.major.pad'	: 7. +(12.-fs)*0.7 ,      # distance to major tick label in points
		'xtick.labelsize'       :'x-large' ,
		'ytick.labelsize'       :'x-large' ,
		'xtick.direction'	: 'in' ,
		'ytick.direction'	: 'in' ,
		'xtick.minor.visible'	: False ,
		'ytick.minor.visible'	: False ,
		'markers.fillstyle'     : 'full' , # full|left|right|bottom|top|none
		'lines.markersize'      : msize ,
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
#if args.paper :
#	fs=14 
#	figs = 3
#	params = {      'figure.figsize'        : (figs,figs) ,
#			'axes.labelsize'        : fs ,
#			'axes.titlesize'        : fs ,
#			'xtick.labelsize'       : fs ,
#			'ytick.labelsize'       : fs ,
#			'font.size'             : fs ,
#			'legend.fontsize'       : fs*0.65 ,
#			'lines.markersize'      : fs*0.4 }
pylab.rcParams.update(params)

AdjustSpaceDefault = [ 0.18 , 0.91 , 0.14 , 0.82 , 0    , 0    ]
ladjust=( AdjustSpaceDefault[0] if (args.ladjust is None) else float(args.ladjust) )
radjust=( AdjustSpaceDefault[1] if (args.radjust is None) else float(args.radjust) )
badjust=( AdjustSpaceDefault[2] if (args.badjust is None) else float(args.badjust) )
tadjust=( AdjustSpaceDefault[3] if (args.tadjust is None) else float(args.tadjust) )
wspace =( AdjustSpaceDefault[4] if (args.wspace  is None) else float(args.wspace ) )
hspace =( AdjustSpaceDefault[5] if (args.hspace  is None) else float(args.hspace ) )
AdjustSpaceDefault = [ ladjust, radjust, badjust, tadjust, wspace, hspace ]

hobjarr = []
ndir = len(args.ddir) 

nax	= 1
f, ax = ( plt.subplots(1, nax, sharey=True, sharex=True) if args.rowsep is False else  plt.subplots(nax, 1, sharey=True, sharex=True) )

def getobjarr( objclass , itemarr )  :
	return [ getattr(objclass,item) for item in itemarr ]

def ftype(tname) :
	if tname.find("float")>-1 :	return float
	if tname.find("str")>-1 :		return str
	if tname.find("complex")>-1 :	return complex
xobj = "D"	; xobjtype = float
yobj = "J"	; yobjtype = float
zobj = "degen_0"
if args.plotx :	xobj = args.plotx
if args.ploty :	yobj = args.ploty
if args.plotz :	zobj = args.plotz
if args.typex :	xobjtype = ftype(args.typex)
if args.typey :	yobjtype = ftype(args.typey)
itemarr		= [ xobj, yobj, zobj , 'gptl_0' ]
print "itemarr : ", itemarr

for jj in range(ndir) :
	hobj = headobj( args.ddir[jj] )
	hobj.readParameters()
	hobj.readStatic()
	hobj.readTldos()
	hobjarr.append( hobj ) 
	if args.nplot : 
		hobj.Nplot = int(args.nplot)
	dat = args.ddir[jj]
	print "Reading[%d] :"%jj, dat

if args.setattrobj :
	setObjArr = np.reshape( args.setattrobj , (-1,3) )
	print "setObjArr : ", setObjArr 
	for setObj in setObjArr : 
		setattr( hobjarr[int(setObj[0])] , setObj[1], float(setObj[2]) )

datobjarr	= []
for jj in range(ndir) :
	hobj = hobjarr[jj]
	#print "jj : ", jj 
	datobjarr.append(  getobjarr( hobj, itemarr ) )

for jj in range(ndir) :
	hobj = hobjarr[jj]
	print  jj, " : ", datobjarr[jj]

datobjarr = np.array( datobjarr ) 
print np.shape( datobjarr ) 

d	= datobjarr.transpose()
x	= np.array( d[0], dtype=xobjtype )
y	= np.array( d[1], dtype=yobjtype )
#ax.scatter( x, y, marker="o" )

z = np.array( d[2], dtype=float )
extend	= 'neither'
nlevels = int(args.ncolorlevels) if args.ncolorlevels  else 15
zmin, zmax = [z.min(), z.max()]
if args.zrange : 
	zmin, zmax = np.array( [ args.zrange.split("_")[-2], args.zrange.split("_")[-1] ] , dtype=float ) 
levels = np.linspace( zmin, zmax, nlevels )
print "Levels : ", levels
if args.zrange : 
	cbartick = np.linspace( levels[0],levels[-1], 3 )
	ybound = [ y[0], y[-1] ]
#elif zobj.find("mass")>-1 and zobj.find("Rat")<0 : 
#	levels = np.linspace(0,10, nlevels)
#	cbartick = np.linspace( levels[0],levels[-1], 3 )
#	ybound = [ 3, 5 ] 
elif zobj.find("degen_0")>-1 : 
	levels = range(0,7)
	if args.zmax : levels = range( 0,int(args.zmax)+1 )
	levels = np.array(levels)+0.5
	cbartick = range(int(levels[0]),int(levels[-1])+1)
	extend	= 'max'
	#z = np.round(z)
elif zobj.find("gptl_0")>-1 : 
	nlevels = 4
	levels = np.linspace(12,14, nlevels )
	extend	= 'max'
	#z = np.round(z)
elif zobj.find("SzdegMateval")>-1 : 
	levels = np.linspace(-1e-6,0.25, nlevels )
	cbartick = np.linspace( 0,levels[-1], 6 )
elif zobj.find("chiSStaticTotal")>-1 : 
	#levels = np.linspace(0,50,20)
	cbartick = np.linspace( levels[0],levels[-1], 3 )
	pass
elif zobj.find("chitotszszStaticwn0")>-1 : 
	#levels = np.linspace(0,50,20)
	cbartick = np.linspace( levels[0],levels[-1], 3 )
	extend = 'max'
if args.extend : extend = args.extend
#elif zobj.find("powerfit")>-1 : 
#	nlevels = 21
#	nlevels = int(args.ncolorlevels) if args.ncolorlevels  else 21
#	levels = np.linspace(0.4,1, nlevels )
#	cbartick = np.linspace( levels[0],levels[-1], 4 )
else : 
	cbartick = np.linspace( levels[0],levels[-1], 3 )
	ybound = [ y[0], y[-1] ]

if args.cbartick : 
	cbartick = np.array( args.cbartick ) 

#from scipy.interpolate import griddata
#xi = np.linspace(0,1,100)
#yi = np.linspace(3,5,100)
## grid the data.
#zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
## contour the gridded data, plotting dots at the randomly spaced data points.
#levels = np.linspace(0,10,15)
#CS = ax.contour(xi,yi,zi,levels,linewidths=0.5,colors='k')
#CS = ax.contourf(xi,yi,zi,levels,cmap=plt.cm.jet)
#plt.colorbar(CS) # draw colorbar
## plot data points.
#ax.scatter(x,y,marker='o',c='b',s=5)
#ax.set_xlim(0,1)
#ax.set_ylim(3,5)
#npts=len(z)
#ax.set_title( zobj ) 

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

cmap = args.colormap if args.colormap else 'jet'
if cmap.find("bwr_red")>-1 : 
	origcolors = cm.get_cmap( 'bwr' , 256 )
	newcolors  = origcolors( np.linspace(0, 1, 256) )
	#pink = np.array([248./256, 24./256, 148./256, 1])
	#newcolors[:128, :] = pink
	#cmap = ListedColormap( newcolors )
	cmap = ListedColormap( newcolors[128:, :] )
	print "Changed color"
elif cmap == "Blues2" :
	ncc = 256
	cc  = np.linspace(0, 1, ncc )
	cc  = np.power(cc, 0.7 )
	newcolors  = np.array( zip( cc,cc,cc,cc ) )
	newcolors[:,2]=1
	newcolors[:,3]=1
	newcolors = newcolors[-1::-1,:]
	cmap = ListedColormap( newcolors )
elif cmap == "Blues3" :
	origcolors = cm.get_cmap( 'Blues' , 256 )
	newcolors  = origcolors( np.linspace(0, 1, 256) )
	#pink = np.array([248./256, 24./256, 148./256, 1])
	#newcolors[:128, :] = pink
	#print "ORI CORLOR : ", newcolors
	#print "NEW CORLOR : ", newcolors * np.array([0.5,0.5,1.,1.])
	newcolors[:,2]=1
	#sys.exit(1)
	cmap = ListedColormap( newcolors )
	#cmap = ListedColormap( newcolors )
	print "Changed color"
elif cmap == "Reds" :
	origcolors	= cm.get_cmap( 'Reds' , 256 )
	newcolors	= origcolors( np.linspace(0, 1, 256) )
	newcolorsbottom	= origcolors( np.linspace(0, 1, 10) )
	newcolorsbottom[:,:]	= 1
	print "ORI CORLOR bottom : ", newcolorsbottom
	newcolors	= np.vstack( (newcolorsbottom, newcolors) )
	print "NEW CORLOR : ", newcolors
	#top = cm.get_cmap('Oranges_r', 128)
	#bottom = cm.get_cmap('Blues', 128)
	#newcolors = np.vstack((top(np.linspace(0, 1, 128)), bottom(np.linspace(0, 1, 128))))
	cmap = ListedColormap( newcolors )
	#os.exit(1)

import matplotlib.tri as tri

nip = int(args.ninterpolate) if args.ninterpolate else 20
print "# of grid for interpolation (nip) : ", nip
if args.interpolate1d :
        xi = []
        yi = []
        zi = []
        if args.interpolateaxis and args.interpolateaxis=='x' :
                #Tarr = [ 90.6640625, 45.33203125, 22.666015625, 11.3330078125, 0. ]
                Tarr	= 11605. / np.array( [ 128, 256, 512, 1024 ] ) 
		try : 
			jj	= 0
			indJ = np.logical_and( (y-jj)<1e-5 , (y-jj)>-1e-5 )
			if len(x[indJ])>0 :
				Tarr	= np.append( Tarr, 0. ) 
		except : 
			pass
		print "Tarr : ", Tarr 
		nTarr	= len(Tarr)
		xiptrans	= np.array ( [ np.zeros(nTarr) for kk in range(nip) ] )
		yiptrans	= np.array ( [ np.zeros(nTarr) for kk in range(nip) ] )
		ziptrans	= np.array ( [ np.zeros(nTarr) for kk in range(nip) ] )
                for iT in range(nTarr) :
			jj = Tarr[iT]
			print "---------- component : ", jj
                        indJ = np.logical_and( (y-jj)<1e-5 , (y-jj)>-1e-5 )
                        print "indT ", jj, " : " , indJ
                        xfilt = x[indJ]
                        yfilt = y[indJ]
                        zfilt = z[indJ]
                        print "xfilt : " , xfilt
                        print "yfilt : " , yfilt
                        print "zfilt : " , zfilt
                        from scipy.interpolate import interp1d
                        #f = interp1d(y)
                        fi = interp1d(xfilt, zfilt, kind='quadratic')

			xipmin	= 0.
			xipmax	= 0.7
			if args.interpolaterange : 
				xipmin, xipmax = np.array( args.interpolaterange )
                        xip     = np.linspace( xipmin, xipmax,   nip)
                        zip1d   = fi(xip)
                        print jj, "xip : ",     xip
                        print jj, "zip1d : ",   zip1d
                        zi.append( zip1d )
                        xi.append( xip   )
                        yi.append( np.zeros(nip)+jj )

			if args.interpolatetranspose is False : 
				ms = markersdiff[iT] # "s"
				if args.dashedline :
					if iT<1 : 
						dashpt = (None,None)
					else : 
						dashpt = [iT+1,1]
				else : dashpt = (None,None)
				plt.plot( xip,zip1d, "C%d"%iT+"-" , dashes=dashpt ) #, label=filterlabel(yobj,paper=args.paper)+"={:.2f} (interp.)".format(jj) )
				plt.plot( xfilt,zfilt, "C%d"%iT+ms, label=filterlabel(yobj,paper=args.paper)+" = {:.2f}".format(jj) , mec='k' )
			else : 
				for iip in range(nip) :
					xiptrans[iip,iT]	= xip[iip] 
					yiptrans[iip,iT]	= jj 
					ziptrans[iip,iT]	= zip1d[iip]
				print "xiptrans : ", xiptrans
				print "yiptrans : ", yiptrans
				print "ziptrans : ", ziptrans
        else :
		if yobj == "J" : 
        	        JHarr = [0,0.2,0.4,0.5,0.7]
			Yarr	= JHarr
		else : 
			print "ERROR :: invalid in interpolation for yobj rather than 'J'" 
			print "exit."
			sys.exit(1)
                for jj in JHarr :
                        indJ = np.logical_and( (x-jj)<1e-5 , (x-jj)>-1e-5 )
                        print "indJ ", jj, " : " , indJ
                        xfilt = x[indJ]
                        yfilt = y[indJ]
                        zfilt = z[indJ]
                        print "xfilt : " , xfilt
                        print "yfilt : " , yfilt
                        print "zfilt : " , zfilt
                        from scipy.interpolate import interp1d
                        #f = interp1d(y)
                        fi = interp1d(yfilt, zfilt, kind='quadratic')
                        yip     = np.linspace(0, 90.66, nip)
                        zip1d   = fi(yip)
                        print jj, "yip : ",     yip
                        print jj, "zip1d : ",   zip1d
                        zi.append( zip1d )
                        yi.append( yip   )
                        xi.append( np.zeros(nip)+jj )
			ms = markersdiff[jj] # "s"
                        plt.plot( yip,zip1d, label=yobj+"={:.2f}(interp.)".format(jj) )
                        plt.plot( yfilt,zfilt, ms, label=yobj+"={:.2f}".format(jj) , mec='k' )
	#AdjustSpaceDefault = [ 0.18 , 0.93 , 0.21 , 0.95 , 0 , 0 ]
	ax.set_xlabel( filterlabel(xobj,paper=args.paper,unit=True) )
	ax.set_ylabel( filterlabel(zobj,paper=args.paper,unit=True) )
	if args.interpolatetranspose : 
		ax.set_xlabel( filterlabel(yobj,paper=args.paper,unit=True) )
		jj=0
		for iip in range(nip) :
			iplot = False
			if args.interpolatechoosex : 
				for chx in args.interpolatechoosex : 
					if (np.abs(chx-xiptrans[iip][0])<1e-2) : 
						iplot=True
			else :
				iplot = True
			if iplot : 
				if args.dashedline :
					if iip<1 : 
						dashpt = (None,None)
					else : 
						dashpt = [iip+1,1]
				else : dashpt = (None,None)
				ms = markersdiff[iip] # "s"
				cs = "C%d"%jj
				jj+=1
				plt.plot( yiptrans[iip], ziptrans[iip], ms+"-", label=filterlabel(xobj,paper=args.paper)+" = {:.2f} (interp.)".format(xiptrans[iip][0]) , dashes=dashpt, mec=cs, mfc='w' )
	if args.xxrange : 
		axSetRange( ax , args.xxrange, 'x' ) 
	if args.yrange : 
		axSetRange( ax , args.yrange,  'y' ) 
	plt.legend()
	plt.grid(linestyle=":")
	try : 
		argsavefigname = hobj.tdir+"/"+hobj.fignamepart()+"_summarycurve"
	except : 
		argsavefigname = "./plot_summary_ipcurve"
	namesuffix =""
	for ar in sys.argv[1:] :
	        print "arg : ", ar
		if ar.find("Dir")>-1 : 
			pass
		else : 
	        	argsavefigname = argsavefigname + '_' + ar.split("-")[-1]
	print "length(argsavefigname)" , len(argsavefigname)
	if len(argsavefigname)>350 : argsavefigname = argsavefigname[:350]
	fsave = argsavefigname
	ladjust=( AdjustSpaceDefault[0] if (args.ladjust is None) else float(args.ladjust) )
	radjust=( AdjustSpaceDefault[1] if (args.radjust is None) else float(args.radjust) )
	badjust=( AdjustSpaceDefault[2] if (args.badjust is None) else float(args.badjust) )
	tadjust=( AdjustSpaceDefault[3] if (args.tadjust is None) else float(args.tadjust) )
	wspace =( AdjustSpaceDefault[4] if (args.wspace  is None) else float(args.wspace ) )
	hspace =( AdjustSpaceDefault[5] if (args.hspace  is None) else float(args.hspace ) )
	plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, wspace=wspace , hspace=hspace )
	pylab.rcParams.update(params)
	ndpi=400
	if args.ndpi : ndpi	= int(args.ndpi)
	transparent = False 
	if args.savetransparent : transparent = True
	if args.cptmp :
		dtype = "pdf"
		if args.pdf :
			dtype="pdf"
		if args.png :
			dtype="png"
		print "Saved : ", fsave + "."+dtype
		plt.savefig( fsave + "."+dtype , dpi=ndpi, transparent=transparent) 
		cptmpftn( fsave, cmd="scp", destMac=True , dtype=dtype )
	elif args.notshow : 
		pass
	else : 
		plt.show()
        sys.exit(1)
        xi = np.array(xi).flatten()
        yi = np.array(yi).flatten()
        zi = np.array(zi).flatten()
        cntr2 = ax.tricontourf(xi, yi, zi,
                        levels=levels,
                        cmap=cmap,
                        extend=extend
                        )
	if args.tflline :
		#ax.tricontour(xi, yi, zi, levels=[0.8], linewidths=1.5, linestyles='dashed', colors='k')
		#ax.tricontour(xi, yi, zi, levels=[0.735], linewidths=1.5, linestyles='--', colors='k')
		ax.tricontour(xi, yi, zi, levels=[0.6, 0.8], linewidths=0.5, linestyles='-', colors='k')
		#pass
        elif args.nocontourline : pass
        else :
                ax.tricontour(xi, yi, zi, levels=levels, linewidths=0.5, colors='k')
else : 
	if args.nocontourline : pass
	elif args.tflline :
		#ax.tricontour(x, y, z, levels=[0.735], linewidths=1.5, linestyles='--', colors='k')
		ax.tricontour(x, y, z, levels=[0.6, 0.8], linewidths=0.5, linestyles='-', colors='k')
	else : 
		ax.tricontour(x, y, z, levels=levels, linewidths=0.5, colors='k')

	if args.interpolatecontour :
		from scipy.interpolate import griddata
		ipmethod	= "linear" #"cubic" #"linear"
		nip		= 100
		xreg	= np.linspace(0,6,nip)
		yreg	= np.linspace(0,9,nip)
		X, Y = np.meshgrid( xreg, yreg ) 
		Ti = griddata((x, y), z, (X, Y), method=ipmethod)
		cntr2 = ax.contourf(X, Y, Ti,
				levels=levels, 
				cmap=cmap, 
				extend=extend
				)
	else : 
		cntr2 = ax.tricontourf(x, y, z, 
				levels=levels, 
				cmap=cmap, 
				extend=extend
				)

if args.tflline is None :
	xfl = []
	yfl = []
	for jj in range(ndir) :
		hobj = hobjarr[jj]
		if hobj.Treal < 1. :
			xfl.append(  getattr( hobj, xobj ) )
			yfl.append(  getattr( hobj, "Tlinwnmaxxyu" ) )
	xfl = np.array( xfl ) 
	yfl = np.array( yfl ) 
	print "xfl : ", xfl
	print "yfl : ", yfl
	ax.plot( xfl, yfl, 'ko--' ) 

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, FuncFormatter)
if args.customcbar : 
        hoff    = 0.1	#0.02
        voff    = 0.015	#0.02
        vwidth  = 0.015 #0.02
        cbarpos = [AdjustSpaceDefault[0]+hoff,AdjustSpaceDefault[3]+vwidth+voff, AdjustSpaceDefault[1]-AdjustSpaceDefault[0]-hoff*2, vwidth]
	if args.cbarpos :
		cbarposval	 = np.array( args.cbarpos, dtype=float ) 
        	cbarpos = [cbarposval[0]+hoff,cbarposval[3]+vwidth+voff, cbarposval[1]-cbarposval[0]-hoff*2, vwidth]
        cax = f.add_axes(cbarpos) #vertleft
	#if abs(levels[-1])<0.02 : 
	#	cbar = f.colorbar(cntr2, cax=cax, extend='neither'  ,orientation='horizontal' , format=FuncFormatter(fmt10) )
	#else : 	
	#	cbar = f.colorbar(cntr2, cax=cax, extend='neither'  ,orientation='horizontal' )
	cbar = f.colorbar(cntr2, cax=cax, extend='neither'  ,orientation='horizontal' )
	cax.xaxis.set_label_position('top')
	cax.xaxis.set_ticks_position('top')
	if zobj=="powerfittxyu"	:
		cbar.ax.set_xlabel(r'$\alpha^{(xy)}$')
	elif zobj=="powerfittyzu"	:
		cbar.ax.set_xlabel(r'$\alpha^{(yz)}$')
	else : 
		cbar.ax.set_xlabel( filterlabel(zobj,args.paper,unit=True) ) 
	if args.cbarlabel : 
		cbarlab = args.cbarlabel
		if cbarlab=='blank' : cbarlab=''
		cbar.ax.set_xlabel( filterlabel(cbarlab,args.paper,unit=True) ) 
		print "cbar label : ", args.cbarlabel
	cbar.ax.tick_params( axis='x', pad = 3 )
	cbar.ax.xaxis.labelpad = 11
elif args.nocbar : 
	pass
else : 
	def fmt(x, pos):
		a, b = '{:.1e}'.format(x).split('e')
		b = int(b)
		return r'${} \times 10^{{{}}}$'.format(a, b)
	def fmt_once(x, pos):
		a, b = '{:.1e}'.format(x).split('e')
		b = int(b)
		return r'${}$'.format(a)
	cbar = f.colorbar(cntr2, ax=ax, extend='neither' )

if args.markdegen : 
	for jj in range(ndir) :
		hobj = hobjarr[jj]
		xi = float( getattr( hobj, xobj )  )
		#print  "GETT : ", yobj, getattr( hobj, yobj ) ,  type(  getattr( hobj, yobj ) )
		yi = float( getattr( hobj, yobj )  )
		deg = int(hobj.degen_0)
		if deg>1 : 
			cs = 'yellow'
			ma = 's'
		elif deg<2 : 
			cs = 'w'
			ma = 'o'
		#if hobj.IT<1 :
		#	print "DEGEN : ", jj, deg, type(deg), deg>1, deg<2, cs, ma
		if 1 : 
			gptl = int(hobj.gptl_0)
			colorsarr	= [ 'w', 'yellow' ] + [ "C%d"%i for i in range(10) ] * 3
			markersarr	= [ 'o', 's' ] + markersdiff[::-1]
			cs = colorsarr[gptl%12]
			ma = markersarr[gptl%12]
			if (gptl/12<1) : 
				cs = colorsarr[gptl%12+1]
				ma = markersarr[gptl%12+1]
			print "GPTL : ", jj, gptl, type(gptl), gptl%12, cs, ma
		ax.plot(xi, yi, ma, ms=5, mec='k', color=cs , clip_on=args.clipon, zorder=9 )
elif args.nomarker : 
	pass
else : 
	ax.plot(x, y, 'wo', mec='k',ms=fs*0.5, clip_on=args.clipon, zorder=9 ) #ms=3, 
#ax.set_xlim(0,1)
#ax.set_ylim(ybound)

ax.grid(linestyle=":")
valax	= plt.axis()
print "axis : ", valax

try : 
	print "cbar label : ", cbartick
	cbar.set_ticks( cbartick )
except : 
	pass

if xobj=="J" and yobj=="U" :
	x = np.linspace( 0,valax[1] , 10 )
	Jlimit = lambda a, J : a*J
	ax.plot( x, Jlimit(3,x) 	, "-" , label=r"$J_H/U=1/{}$".format(3) ) 
	ax.plot( x, Jlimit(6,x)		, ":" , label=r"$J_H/U=1/{}$".format(6) ) 
	ax.legend()

if args.ztick : 
	try :
		zt =  np.array(args.ztick.split("_") , dtype=float)
	except : 
		zt =  np.array(args.ztick.split("_")[1:] , dtype=float)
	cbar.set_ticks( zt ) 
if args.xticklabel : 
	xt =  np.array(args.xticklabel)
	ax.set_xticks( xt ) 
if args.yticklabel : 
	yt =  np.array(args.yticklabel)
	ax.set_yticks( yt ) 


ax.set_xlabel( xobj )
ax.set_ylabel( yobj )

xyobjCompArr = [	
		[ xobj, 'x'] ,
	 	[ yobj, 'y'] ,	
		]
def setAxisLabel( aax, lab, axis )  :
	if axis.find("x")>-1 : 
		aax.set_xlabel( lab )
	elif axis.find("y")>-1 : 
		aax.set_ylabel( lab )
	else : 
		print "ERROR in 'setAxisLabel'."
		sys.exit(1)
		
for xyobjComp in xyobjCompArr : 
	xyobj	= xyobjComp[0]
	axis	= xyobjComp[1]
	if xyobj=="J" : setAxisLabel( ax, r"$J_H$ [eV]",	axis )
	elif xyobj=="U" : setAxisLabel( ax, r"$U$ [eV]",	axis )
	elif xyobj=="nOcculattAra" or xyobj=="nOcculattAraFormer" : 
		if args.paper : 
			setAxisLabel( ax, filterlabel(xyobj,paper=args.paper) ,	axis )
		else : 
			setAxisLabel( ax, r"$n_{\mathrm{latt}}$",	axis )
	elif xyobj=="Tfit"	: setAxisLabel( ax, r"$T_{\mathrm{fit}}$ [K]",	axis )
	elif xyobj=="Treal"	:
		if args.paper : 
			setAxisLabel( ax, r"$T$ [K]",	axis )
		else : 
			setAxisLabel( ax, r"$T_{\mathrm{real}}$ [K]",	axis )
	elif xyobj=="cfilling_0"	:
		setAxisLabel( ax, r"$n_{\mathrm{impurity}}$",	axis )
	else :
		setAxisLabel( ax, filterlabel(xyobj,paper=args.paper), axis ) 

if args.setxlabel is not None :
	print "xlab : ",  args.setxlabel
	setAxisLabel( ax, args.setxlabel,	"x" )
if args.setylabel is not None :
	print "ylab : ",  args.setylabel
	setAxisLabel( ax, args.setylabel,	"y" )
if args.noxlabel : 
	xlab = ''
	setAxisLabel( ax, '',	"x" )
if args.noylabel : 
	ylab = ''
	setAxisLabel( ax, '',	"y" )

if args.ylabelcoord :
	ylabcoordx = float( args.ylabelcoord )
	ax.yaxis.set_label_coords( ylabcoordx, 0.5)


try : 
	print "ax.axis : ", ax.axis()
except : 
	print "plt.axis : ", plt.axis()
if args.xxrange : 
	axSetRange( ax , args.xxrange, 'x' ) 
if args.yrange : 
	axSetRange( ax , args.yrange,  'y' ) 
if args.noxticklabel :
	ax.set_xticklabels([])
if args.noyticklabel :
	ax.set_yticklabels([])
try : 
	print "ax.axis (mod) : ", ax.axis()
except : 
	print "plt.axis (mod) : ", plt.axis()

titleax = ax
if args.notitle : pass
elif args.titlevar :
	ax.set_title( zobj ) 
else : 
	try :
		titleax.set_title( hobj.title + " n={:.2f}".format( getattr(hobj,"cflattu_%d"%ni) )  ) 
	except : 
		titleax.set_title( hobj.title ) 


try : 
	argsavefigname = hobj.tdir+"/"+hobj.fignamepart()+"_summary"
except : 
	argsavefigname = "./plot_summary_contour"
namesuffix =""
for ar in sys.argv[1:] :
        print "arg : ", ar
	if ar.find("Dir")>-1 : 
		pass
	else : 
        	argsavefigname = argsavefigname + '_' + ar.split("-")[-1]
print "length(argsavefigname)" , len(argsavefigname)
if len(argsavefigname)>350 : argsavefigname = argsavefigname[:350]
fsave = argsavefigname
ladjust=( AdjustSpaceDefault[0] if (args.ladjust is None) else float(args.ladjust) )
radjust=( AdjustSpaceDefault[1] if (args.radjust is None) else float(args.radjust) )
badjust=( AdjustSpaceDefault[2] if (args.badjust is None) else float(args.badjust) )
tadjust=( AdjustSpaceDefault[3] if (args.tadjust is None) else float(args.tadjust) )
wspace =( AdjustSpaceDefault[4] if (args.wspace  is None) else float(args.wspace ) )
hspace =( AdjustSpaceDefault[5] if (args.hspace  is None) else float(args.hspace ) )
plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, wspace=wspace , hspace=hspace )
print "lad, bad, rad, tad, wsp, hsp : " , ladjust, badjust, radjust, tadjust, wspace, hspace 
pylab.rcParams.update(params)
ndpi=400
if args.ndpi : ndpi	= int(args.ndpi)
#ax.clear()
#ax.axis('off')

#import matplotlib.axis as axis
#aa=  axis.Axis(ax)
#print "LABEL PAD : ",  aa.get_tick_padding()
transparent = False 
if args.savetransparent : transparent = True
if args.cptmp :
	dtype = "pdf"
	if args.pdf :
		dtype="pdf"
	if args.png :
		dtype="png"
	print "Saved : ", fsave + "."+dtype
	plt.savefig( fsave + "."+dtype , dpi=ndpi, transparent=transparent) 
	cptmpftn( fsave, cmd="scp", destMac=False , dtype=dtype )
elif args.notshow : 
	pass
else : 
	plt.show()

