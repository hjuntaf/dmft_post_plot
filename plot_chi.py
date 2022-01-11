#!/opt/python/2.7.15/gcc-4.8.5/bin/python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys, os

import argparse

from basislabel import *
from hjl_common import *

parser = argparse.ArgumentParser(description="Plot chi.") 
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
parser.add_argument("--pdf", action="store_true", help="Save the figure in .pdf" ) 
parser.add_argument("--cptmp", action="store_true", help="Copy the saved figure into Dropbox/tmp_linux/") 
parser.add_argument("--notitle", action="store_true", help="Remove the title caption." )
parser.add_argument("--squre", action="store_true", help="Plot in square frame.")
parser.add_argument("-s", "--separate", action="store_true" , help="Plot data in different frame for each orbital." ) 
parser.add_argument("-w", "--w0", action="store_true" , help="Plot data near 0 ." ) 
parser.add_argument("-l", "--llog", action="store_true" , help="Plot data in log-log plot." ) 
parser.add_argument("-m", "--minus", action="store_true" , help="Plot self-energy with the original minus sign." ) 
parser.add_argument("--minusarray", "--minusarr", "--marr", type=str, help="Specify if you will use a minus sign of the data.")
parser.add_argument("--minusx", "--mx", action="store_true" , help="Plot minus x-data.")
parser.add_argument("--minuswithoutsign", "--mwos", action="store_true" , help="Remove the minus sign on label.")
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
parser.add_argument("-r", "--rew", "--wre", action="store_true", help="Plot the data of real-frequency spectrum." ) 
parser.add_argument("--fillb", action="store_true", help="Fill the plotting of blue one." )
parser.add_argument("--offsetecluster", "--offecl", "--offsetecl", action="store_true", help="Set off-set values of Ecluster for diagonal.")
parser.add_argument("--operator", "--op", type=str, help="Specify the type of operator of chi to plot")
parser.add_argument("--offdiag", "--offd", action="store_true", help="Plot the off-diag.") 
parser.add_argument("-d", "--diag", action="store_true", help="Plot the diagonal only") 
parser.add_argument("--sum", action="store_true", help="Plot the summation of the data in 'chdat'.") 
parser.add_argument("--dimdat", "--dimd", type=str, help="Specify the dimension of the data matrix") 
parser.add_argument("--realdata", "--realdat", "--rd", action="store_true", help="Plotting real-data" )
parser.add_argument("--upper", type=str, help="Specify the # of coeff of cont. fract.")
parser.add_argument("--upperarray", "--upperarr", "--uparr", type=str, help="Specify the # of coeff of cont. fract.")
parser.add_argument("--leglabel", "--legl", type=str, help="Set the labels of the legend" ) 
parser.add_argument("--ladjust", "--lad", type=str, help="Set   left-adjust value." ) 
parser.add_argument("--radjust", "--rad", type=str, help="Set  right-adjust value." ) 
parser.add_argument("--badjust", "--bad", type=str, help="Set bottom-adjust value." ) 
parser.add_argument("--tadjust", "--tad", type=str, help="Set    top-adjust value." ) 
parser.add_argument("--guideline", "--guide", action="store_true",  help="Plot the guide-line.")
parser.add_argument("--figs", "--figsize", type=str, help="Set the labels of the legend" ) 
parser.add_argument("--vline", "--vl", type=str, help="Plot a vertical line of the specified value on x-axis.")
parser.add_argument("--showfermidat", "--sf", action="store_true",  help="Show  data at the Fermi level.")
parser.add_argument("--halfnratio", "--hnr", "--halfratio", "--halfaniso", action="store_true",  help="Plot ratio of (x_{00}/2.)/x_{22}.")
parser.add_argument("--halffirstcomp", "--hfc", "--halffirst", "--hfirst", "--hf", action="store_true",  help="Denominate the value of 0-th component by 2.")
parser.add_argument("--nomarker", "--nm", "--nom", action="store_true",  help="Plot with lines, hiding markers.")
parser.add_argument("--diffmarker", "--dm", "--diffm", action="store_true",  help="Plot each data with different markers.")
parser.add_argument("--markersize", "--msize", type=str,  help="Set the size of markers.")
parser.add_argument("--solidline", "--solidl", "--sl", action="store_true",  help="Plot with sold lines.")
parser.add_argument("--dashedline", "--dashl", action="store_true", help="Plot the data as dashedlines.")
parser.add_argument("--filesuffix", "--fsuffix", "--fsuff", type=str,  help="Specify a suffix-folder name of the data.")
parser.add_argument("--filesuffixarray", "--fsuffixarr", "--fsuffarr", type=str,  help="Specify a suffix-folder name of the data.")
parser.add_argument("--write", "--wr", action="store_true", help="Write the data.")
parser.add_argument("--noshow", "--ns", action="store_true", help="Do not plot the data.")
parser.add_argument("--epsilon", "--ep", type=str,  help="Set Epsilon.")
args = parser.parse_args()

import matplotlib.pylab as pylab
fs = 12
if args.diag : 
	figrat = 1.2
else : 
	figrat = 2.2 
figsy=5 ; figs = figsy * figrat
if args.figs :
	figsarr = args.figs.split("_")
	figs    = float(figsarr[0])
	figsy   = float(figsarr[1])
print "figsize = ({}, {})".format( figs, figsy ) 
if args.markersize :	msize = float(args.markersize)
else :			msize = fs*0.4
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
			'lines.markersize'      : msize,
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
			'lines.markersize'      : msize }
#font = {'family' : 'normal',
#       'weight' : 'bold',
#       'size'   : 22}
pylab.rcParams.update(params)
print "\n"
nbasis = 6
basis  = "j"
chdat = [ 0, 1, 2 ]
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

if args.chdat : 
	chdat = np.array( args.chdat.split("_") , dtype=int ) 
if args.epsilon : 
	epsilon = float( args.epsilon )
else :
	epsilon = False

hobjarr = []
ndir = len(args.ddir) 
if args.one :	nax =  1
else : 		nax = ndir 

if args.diag : pass
else : 
	nax *= 2 
if args.dimdat :
	dimdat = int(args.dimdat)
else : dimdat  = 3
chdatoffd = []
chdatoffd2= []
for i in range(dimdat) :
	for j in range(i+1,dimdat) :
		chdatoffd.append(  [i,j] ) 
		chdatoffd2.append( [j,i] ) 

f, ax = plt.subplots( 1, nax, sharey=True, sharex=True )

if args.operator : 
	operator = args.operator
	operatorarr = [ operator ]
	if args.operator.find("_")>-1 : 
		operatorarr = args.operator.split("_")
		operator = operatorarr[0]
else :
	operator = "jzjz"
	operatorarr = [ operator ]
if args.realdata :	realdat = True 
else :			realdat = False

for jj in range(len(args.ddir)) :
	if nax == 1 : axnow = ax 
	elif nax > 1 : axnow = ax[jj] 
        if nax > 1 and args.diag is False :	axnowoffd = ax[jj+1]
	hobj = headobj( args.ddir[jj] )
	hobj.readParameters()
	hobjarr.append( hobj ) 
	if args.upper :	setattr( hobj, "upper", int(args.upper) )
	else :		setattr( hobj, "upper", 100 ) 
	if args.upperarray :
		upperarray = args.upperarray.split("_") 
		upper = int(upperarray[jj])
		setattr( hobj, "upper", upper ) 
	if args.filesuffix :	filesuffix = args.filesuffix
	else :			filesuffix = None
	if args.filesuffixarray :
		filesuffixarray	= args.filesuffixarray.split("_") 
		filesuffix 	= filesuffixarray[jj]
		print "filesuffix : ", filesuffix
	if args.minusarray :
		minusarray = args.minusarray.split("_") 
		minus = True if int(minusarray[jj]) else False
		setattr( args, "minus", minus ) 
	if args.count :	count = args.count
	else :		count = hobj.count
	if args.countarr :	
		if len(args.countarr.split("_"))>1 : 
			try : 
				count = int( args.countarr.split("_")[0] )
				count = int( args.countarr.split("_")[jj] )
			except  :
				count = int( args.countarr.split("_")[jj+1] )
		else : 
			count = int( args.countarr ) 
	dat = args.ddir[jj]
	print "Reading :", dat
	selfeobj = datobj()
	lc = ""
	lc2= ""
	#lc = "C{}".format( jj ) 
	#lc2= "C{}".format( jj+1 ) 
	#lc3= "C{}".format( jj+2 )
	labelarrsquare = ['xx+yy','','zz']
	if args.halffirstcomp : labelarrsquare[0] = '(xx+yy)/2'
	if args.sum :
		if args.leglabel : 
			#leglabel = "${}={}$".format(args.leglabel ,  getattr( hobj, args.leglabel ) ) 
			labsuffix = filterlabel(args.leglabel) + "={}".format( getattr( hobj, args.leglabel ) )
			if args.leglabel.find("filling")>-1 : 
				labsuffix = filterlabel(args.leglabel) + "={:.2f}".format( float(getattr( hobj, args.leglabel )) )
			#leglabel = [labsuffix]*nbasis
			leglabel = labsuffix
		else : 
			leglabel = "Sum"
		labelarrsquare = [ leglabel ] * 3 
	if args.realdata :	dtypelab = "Re"
	else :			dtypelab = "Im"
	if args.rew :	freqtypelab = "\omega"
	else :		freqtypelab = "i\omega_n"
	nop = len(operatorarr)
	for iop in range(nop) : 
		operator = operatorarr[iop]
		functlab = "{}$\chi_{{{:s}{:s}}}( {})$".format(dtypelab, operator[0].upper(), operator[0].upper() , freqtypelab )
		if nop>1 : 
			functlab = "{}$\chi_{{{:s}{:s}}}( {})$".format(dtypelab, "","" , freqtypelab )
			if args.sum : 
				labelarrsquare = [ "{}$\chi_{{{:s}{:s}}}( {})$".format(dtypelab, operator[0].upper(), operator[0].upper() , freqtypelab ) ]*3
			lc =""
		ms = False
		if args.diffmarker : 
			if ndir > 1 :	ms = markersdiff[jj]
		if args.minus : 
			if args.minuswithoutsign :  pass
			else :
				functlab = "-"+functlab
		if args.dashedline :
			if jj<1 : 
				dashpt = (None,None)
			else : 
				dashpt = [jj+1,1]
			if nop>1 : 
				dashpt = [jj*nop+iop+1,1]
		else : dashpt = (None,None)
		if args.rew : 
			selfeobj.getplotdataiw_lc( axnow     , dat , "alldiff", nbasis, chdat,     hobj , args , lc , "chiw" +operator , functlab=functlab,   varlab="$\omega$",   nondiag=False, realdata=realdat, labelarr = labelarrsquare, logplot=args.llog, marker=ms , dashes=dashpt, epsilon=epsilon ,filesuffix=filesuffix ) 
			if args.diag is False : 
				selfeobj.getplotdataiw_lc( axnowoffd , dat , "alldiff", nbasis, chdatoffd, hobj , args , lc2, "chiw" +operator , functlab="", varlab="$\omega$",   nondiag=True, realdata=realdat, labelarr = labelarrsquare, logplot=args.llog, marker=ms , dashes=dashpt, epsilon=epsilon , filesuffix=filesuffix) 
			#selfeobj.getplotdataiw_lc( axnowoffd , dat , "alldiff", nbasis, chdatoffd2,hobj , args , lc3, "chiw" +operator , functlab=""						, varlab="$\omega$",   nondiag=True) 
		else : 
			selfeobj.getplotdataiw_lc( axnow     , dat , "alldiff", nbasis, chdat,     hobj , args , lc , "chiiw"+operator , functlab=functlab,   varlab="$\omega_n$", nondiag=False, realdata=realdat, labelarr = labelarrsquare, logplot=args.llog, marker=ms , dashes=dashpt, epsilon=epsilon , filesuffix=filesuffix) 
			if args.diag is False : 
				selfeobj.getplotdataiw_lc( axnowoffd , dat , "alldiff", nbasis, chdatoffd, hobj , args , lc2, "chiiw"+operator , functlab="", varlab="$\omega_n$", nondiag=True, realdata=realdat, labelarr = labelarrsquare, logplot=args.llog, marker=ms , dashes=dashpt, epsilon=epsilon , filesuffix=filesuffix) 
	
	if args.vline : 
		redline = float(args.vline) # 0.048594999999998834
		axnow.axvline(  redline , color='red' , lw = 1 , ls='--' ) 
		axnow.axvline( -redline , color='red' , lw = 1 , ls='--' ) 

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

print "AXIS : ", plt.axis()
if args.guideline : 
	axisarr = plt.axis()
	plt.axis()
	def guideftn( x ) :
		y = 50 * x 
		return y
	xguide  = np.linspace( axisarr[0], axisarr[1], 30 )  
	axnow.loglog( xguide , guideftn(xguide) , "k-.", label=r"$x^{1.0}$" )

	def guideftn( x ) :
		y = np.power(x,-2)
		return y
	xguide  = np.linspace( axisarr[0], axisarr[1], 30 )  
	axnow.loglog( xguide , guideftn(xguide) , "k--", label=r"$x^{-2.0}$" )
if args.halfnratio : 
	ylab = r"$\chi^{}$".format(operator.upper()[0]) + "$_{xy}$"
	ylab = ylab + r" / $\chi^{}$".format(operator.upper()[0]) + "$_{z}$"
	axnow.set_ylabel( ylab )
	axnow.tick_params( axis='x',  pad=8)
#xtic = [pow(10,-1.*a) for a in range(3,-2,-2)]
#ax[0].set_xticks( xtic )
#for axnow in ax : 
	#for xticlab in axnow.get_xticklabels() :
		#xticlab.set_rotation(30)
#print "XTIC LBA : ", xticlab
#ax[0].minorticks_on()
titleax = axnow
if args.diag is None : titleax2= axnowoffd
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
	if args.diag is None : titleax2.legend()
if args.tl : plt.tight_layout()

dataname = "chi" 
if args.rew :	dataname=dataname+"w" +operator
else :		dataname=dataname+"iw"+operator
argsavefigname = hobj.tdir+"/"+hobj.fignamepart()+"_"+dataname
for ar in sys.argv[ndir+1:] :
        print "arg : ", ar
        argsavefigname = argsavefigname + '_' + ar.split("-")[-1]
fsave = argsavefigname

ladjust=( 0.18 if (args.ladjust is None) else float(args.ladjust) )
radjust=( 0.93 if (args.radjust is None) else float(args.radjust) )
badjust=( 0.22 if (args.badjust is None) else float(args.badjust) )
tadjust=( 0.96 if (args.tadjust is None) else float(args.tadjust) )
plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, wspace=0.08 , hspace=0 )

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
elif args.noshow : 
	pass
else : 
	plt.show()
