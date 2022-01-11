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
parser.add_argument("-dy", "--datyrange", type=str, help="set datyrange e.g. _y0_y1" )
parser.add_argument("-dx", "--datxrange", type=str, help="set datrange e.g. _x0_x1" )
parser.add_argument("-b", "--basis", type=str, help="Set basis. e.g. \"-b j | -b t2g\"" ) 
parser.add_argument("-c", "--chdat", type=str, help="Choose a component(s) of basis. e.g. \"-c _0_1_2 (default 0,2,4)\"" ) 
parser.add_argument("--chdatoffd", "--co", "--chdatoff", type=str, help="Choose a component(s) of basis. e.g. \"-c _0_1_2_4 (giving (0,1) and (2,4) components.)\"" ) 
parser.add_argument("-u", "--ups", action="store_true" , help="Plot (pseudo)spin-up components of basis" )
parser.add_argument("--pdf", action="store_true", help="Save the figure in .pdf" ) 
parser.add_argument("--cptmp", action="store_true", help="Copy the saved figure into Dropbox/tmp_linux/") 
parser.add_argument("--notitle", action="store_true", help="Remove the title caption." )
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
parser.add_argument("--offdiagratio", "--offdrat", "--offdr", action='store_true', help="Show ratios of off-diagonal to diagonal part of imaginary.")
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
parser.add_argument("--transform", "--transf", type=str, help="Transform the matrix data." )
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
parser.add_argument("--removelinear", "--rmlin", "--rlin", action="store_true", help="Fit linear-line and remove the linear-part in the data .")
parser.add_argument("--setlinear", "--setlin", "--slin", type=str, help="Set linear-line steepness and remove the linear-part in the data .")
parser.add_argument("--drawabc", "--dabc", "--abc", type=str, help="Draw a 2nd-order polynomial a*x*x+b*x+c with a_b_c.")
parser.add_argument("--carray", "--carr", type=str, help="Set a constant term of a 2nd-order polynomial a*x*x+b*x+c.")
parser.add_argument("--finegrid", "--fineg", "--fg", action="store_true", help="Plot the data of fine-grid's.")
parser.add_argument("--finegrid2", "--fineg2", "--fg2", action="store_true", help="Plot the data of fine-grid's.")
parser.add_argument("--finegrid3", "--fineg3", "--fg3", action="store_true", help="Plot the data of fine-grid's according to 'beta_real'.")
parser.add_argument("--finegrid4", "--fineg4", "--fg4", action="store_true", help="Plot the data of fine-grid's from 'self4.c'.")
parser.add_argument("--markersize", "--msize", type=str, help="Set the marker-size.") 
parser.add_argument("--cutoff", "--coff", type=str, help="Set lower-bound cutoff of x-data to be plotted .") 
parser.add_argument("--lowestmatsu", "--lmatsu", "--wn0", action='store_true', help="Draw a line of lowest-Matsubara frequency pi/beta.")
parser.add_argument("--wnbetaline", "--wbline", "--wbl", action='store_true', help="Draw lines of pi/beta within [1024,512,256,128,64].")
parser.add_argument("--difftwopoints", "--difftwo", "--difft", "--dtwo", action='store_true', help="Draw the two-point average slope.")
parser.add_argument("--object", "--dobj", "--obj", type=str, help="Specify the type of data you want to plot. (e.g. 'selfiw','glocaliw')")
parser.add_argument("--indimp", "--iimp", type=str, help="Specify the index of impurity atom.")
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
else :			msize = fs*0.4
params = {      'legend.fontsize'       : fs ,
		'figure.figsize'        : (figs,figsy) ,
		'axes.labelsize'        : 'x-large' ,
		'axes.titlesize'        :'small' ,
		'axes.linewidth'        : 1.5 ,
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
fname = "self"
markers = markersj
if args.basis : 
	if args.basis.find("t")>-1 : 
		basis = "t2g"
		chdat = [ 0, 2, 4 ]
		fname = "selft2g"
		markers = markerst
if args.transform : 
	if args.transform.find("t")>-1 : 
		chdat = [ 0, 2, 4 ]
	elif args.transform.find("j")>-1 : 
		chdat = [ 0, 1, 4 ]
	elif args.transform.find("l")>-1 : 
		chdat = [ 0, 2, 4 ]
chdatoffd = []
if args.offdiag or args.offdiagratio : 
	for i in range(nbasis) :
		for j in range(i+1,nbasis) :
			chdatoffd.append( [i,j] ) 
if args.chdatoffd :
	dumchdat = np.array( args.chdatoffd.split("_") , dtype=int )
	chdatoffd = np.reshape( dumchdat , (-1,2) )
obj = 'glocaliw'
if args.object : 
	obj = args.object
if obj.find("selfiw")>-1 : 
	ylaborig = r"$\Sigma$"
elif obj.find("glociw")>-1 : 
	ylaborig = r"$G_{loc}$"
elif obj.find("glocw")>-1 : 
	ylaborig = r"$G_{loc}$"
else : 
	print "Warnning :: Setting 'ylaborig' in obj."
	ylaborig = filterlabel( obj ) 
print "OBJ : ", obj
	
#basis = "j"
#chdat = [ 0, 1, 4 ]
#fname = "self"
#markers = markersj

def simplelogplot( xdat, ydat, labelcomp, ax , laboverw=False, msls=False , msize=None )  :
	if laboverw :	labelstr = laboverw
	else :		labelstr = bolabel(basis,labelcomp)
	mslsinit = markers[labelcomp]+linestyles[labelcomp]
	if msls : pass
	else : msls = mslsinit
	try : 
		ax.loglog( xdat, ydat , msls, label=labelstr , mec=None , markersize=msize ) #(None if (labelcomp==0 and basis.find("j")>-1) else 'gray') ) #, mfc='w', mew=2) #
	except :
		print "not-plotted"
		pass

sgnself = 1. ;
minusind = ""
minusindx= ""
if args.chdat : 
	chdat = np.array( args.chdat.split("_") , dtype=int ) 
if args.plusself : 
	sgnself = 1. ;
	minusind = ""
if args.minus : 
	sgnself = -1. ;
	minusind = "-"
if args.minusx : 
	minusindx= "-"


hobjarr = []
ndir = len(args.ddir) 
if args.one :	nax  = 1
else :		nax  = ndir
f, ax = ( plt.subplots(1, nax, sharey=True, sharex=True) if args.rowsep is False else  plt.subplots(nax, 1, sharey=True, sharex=True) )

for jj in range(len(args.ddir)) :
	hobj = headobj( args.ddir[jj] )
	hobj.readParameters()
	hobjarr.append( hobj ) 
	dat = args.ddir[jj]
	print "Reading :", dat

	selfeobj = datobj()
	#if args.jnat : 
	#	gdatReDiag , gdatImDiag = matrixReturnDiagTransfjnat( gdat , nbasis, hobj.matTjnat ) 
	#else : 
	#	gdatReDiag , gdatImDiag = matrixDiagReturnReIm( gdat , nbasis ) 

	if nax<2 :	axnow = ax 
	else : 		axnow = ax[jj]
	if args.leglabel :	leglabel = "${}={}$".format(args.leglabel ,  getattr( hobj, args.leglabel ) ) 
	else : 			leglabel = False
	if args.line :	msls = markersn[jj]
	else :		msls = False
	if args.dashedline :
		if jj<1 : 
			dashpt = (None,None)
		else : 
			dashpt = [jj+1,1]
	else : dashpt = (None,None)
	if args.epsilon : setattr( hobj, "epsilon", float(args.epsilon) ) 
	else : setattr( hobj, "epsilon", 0.03 )
	if args.llog :	logplot=True
	else :		logplot=False
	if args.transform :	
		if args.transform.find("t")>-1 :	T=Taraph
		elif args.transform.find("j")>-1 :	T=np.conjugate(Taraph).transpose()
		elif args.transform.find("l")>-1 :	T=np.dot( Taraph , np.conjugate(Tlefft2g).transpose() ) 
	else :						T=False

	selfeobj.getGFgeneral( dat, nbasis, hobj, args, obj=obj ) 
	gf	= getattr( selfeobj, obj ) 
	#print gf
	print "gf.shape : ", np.shape( gf ) 
	
	gftest	=  gf[:,:,0] 
	print "gftest : "
	printArrMatAuto( gftest )

	#Rorb	= [ [ i,i,	 1 ]  for i in range( 12 ) ]		# Identity 
	#Rorb	= [ [ i,i,	-1 ]  for i in range( 12 ) ]		# Inversion
	NU	= 6
	#Rorb	= [
	#		[ 0, NU,	 1 ] ,
	#		[ 2, NU+2, 	-1 ] ,
	#		[ 4, NU+4,  	-1 ] ,
	#		[ NU, 0,	 1 ] ,
	#		[ NU+2, 2,	-1 ] ,
	#		[ NU+4, 4,  	-1 ] ,

	#		[ 1, NU+1,	 1 ] ,
	#		[ 3, NU+3, 	-1 ] ,
	#		[ 5, NU+5,  	-1 ] ,
	#		[ NU+1, 1,	 1 ] ,
	#		[ NU+3, 3,	-1 ] ,
	#		[ NU+5, 5,  	-1 ] 
	#		]
	#Rorb	= [
	#		[ 0, NU,	 1 ] ,
	#		[ 2, NU+2, 	-1 ] ,
	#		[ 4, NU+4,  	-1 ] ,
	#		[ NU, 0,	 1 ] ,
	#		[ NU+2, 2,	-1 ] ,
	#		[ NU+4, 4,  	-1 ] ,

	#		[ 1, NU+1,	 1 ] ,
	#		[ 3, NU+3, 	-1 ] ,
	#		[ 5, NU+5,  	-1 ] ,
	#		[ NU+1, 1,	 1 ] ,
	#		[ NU+3, 3,	-1 ] ,
	#		[ NU+5, 5,  	-1 ] 
	#		]
	RorbMirrorxy	= [
			[ 0, 0+NU,	 1 ] ,
			[ 2, 4+NU, 	 1 ] ,
			[ 4, 2+NU,  	 1 ] ,
			[ 0+NU, 0,	 1 ] ,
			[ 2+NU, 4, 	 1 ] ,
			[ 4+NU, 2,  	 1 ] ,

			#[ 1, 1,		-1 ] ,
			#[ 3, 5, 	 1 ] ,
			#[ 5, 3,  	-1 ] ,
			#[ 1+NU, 1+NU,	-1 ] ,
			#[ 3+NU, 5+NU, 	 1 ] ,
			#[ 5+NU, 3+NU,  	-1 ] 
			]
	RorbC4z	= [
			[ 0, 0,		-1 ] ,
			[ 2, 4, 	 1 ] ,
			[ 4, 2,  	-1 ] ,
			[ 0+NU, 0+NU,	-1 ] ,
			[ 2+NU, 4+NU, 	 1 ] ,
			[ 4+NU, 2+NU,  	-1 ] ,

			#[ 1, 1,		-1 ] ,
			#[ 3, 5, 	 1 ] ,
			#[ 5, 3,  	-1 ] ,
			#[ 1+NU, 1+NU,	-1 ] ,
			#[ 3+NU, 5+NU, 	 1 ] ,
			#[ 5+NU, 3+NU,  	-1 ] 
			]
	RspinC4z	= [
			[ 0, 0,	 	-1j ] ,
			[ 1, 1,	 	 1j ] ,
			[ 2, 2,	 	-1j ] ,
			[ 3, 3,	 	 1j ] ,
			[ 4, 4,	 	-1j ] ,
			[ 5, 5,	 	 1j ] ,

			[ 0+NU, 0+NU,	 	-1j ] ,
			[ 1+NU, 1+NU,	 	 1j ] ,
			[ 2+NU, 2+NU,	 	-1j ] ,
			[ 3+NU, 3+NU,	 	 1j ] ,
			[ 4+NU, 4+NU,	 	-1j ] ,
			[ 5+NU, 5+NU,	 	 1j ] ,
			]
	RspinId		= [ [ i,i, 1 ] for  i in range(12) ]
	#R	= Rorb
	#R	= mul_sparseAB_complex( RorbC4z, RspinC4z, NU*2 )
	R	= mul_sparseAB_complex( RorbMirrorxy, RspinId, NU*2 )
	print "R : ", R
	print "R.R.inv : "
	#printArrMatAuto(  transf_ArrMat_sparseTwf_complex( np.identity(12) , R )  )
	print np.dot( np.conjugate(R.transpose()) , R ) 

	#gftestT	=  transf_ArrMat_sparseTwf_complex( gftest, R ) 
	gftestT	=  np.dot( np.conjugate(R) , np.dot( gftest , R.transpose() )  )
	print "gftestT : "
	printArrMatAuto( gftestT )

	print "Spin-up block : "
	print "gftest : "
	printArrMatAuto(gftest[::2,::2])
	print "gftestT : "
	printArrMatAuto(gftestT[::2,::2])
	print "gftest - gftestT : "
	gres	= gftest - gftestT 
	printArrMatAuto( gres[::2,::2] )


	ecl	= hobj.readEcluster()
	eclT	=  np.dot( np.conjugate(R) , np.dot( ecl , R.transpose() )  )
	print ""
	print "ecl : "
	printArrMatAuto( ecl[::2,::2] )
	print "eclT : "
	printArrMatAuto( eclT[::2,::2] )
	print ""
	print "ecl - eclT : "
	printArrMatAuto( ecl[::2,::2] - eclT[::2,::2] )

	xdat	= selfeobj.wn
	ydat	= gf.imag

	mscsls	= ".-"
	
	for j in chdat : 
		pass
		#axnow.plot( xdat, ydat[j][j] , mscsls ) 

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
	if args.sublabel :
		dx=x1-x0;dy=y1-y0#; axnow.text( dx*0.02, dy*0.98 , markerssub[jj] ) 
		slab = markerssub[jj]
		if args.sublabeloff : slab = markerssub[jj+int(args.sublabeloff)]
		if args.sublabelparameter :  slab = slab +" {}={:.6g}".format( filterlabel(args.sublabelparameter) , getattr(hobj,args.sublabelparameter) )
		axnow.text(0.015, 0.94, slab,
				horizontalalignment='left', verticalalignment='center', transform=axnow.transAxes,
				fontsize=fs+4 )#, weight='semibold')

if nax<2 : 
	xticax	= axnow
	titleax	= axnow
	legax	= axnow
	xlabax	= axnow
	ylabax	= axnow
	axarr	= [ax] * ndir 
else : 
	xticax	= ax[0]
	titleax	= ax[0]
	legax	= ax[0]
	ylabax	= ax[0]
	xlabax	= ax[1]
	axarr	= ax
	if args.rowsep : xlabax = ax[-1] ; ylabax = ax[nax/2] 

if args.llog :
	pass
	#xtic = [pow(10,-1.*a) for a in range(3,-2,-2)]
	#xticax.set_xticks( xtic )
if args.rotatelab :
	for axnow in axarr : 
		for xticlab in axnow.get_xticklabels() : xticlab.set_rotation(30)
		for yticlab in axnow.get_yticklabels() : yticlab.set_rotation(30)
#print "XTIC LBA : ", xticlab
print "axis : ", plt.axis()

if args.lowestmatsu : 
	for jj in range(ndir) :
		if hobjarr[jj].IT > 0  : 
			val = np.pi / hobjarr[jj].IT
		else : 
			val = np.pi / hobjarr[jj].beta
		axarr[jj].axvline( val , linestyle="--", color='k' ) 
if args.wnbetaline : 
	betaarr = [ 1024, 512, 256, 128, 64 ] 
	for jj in range(ndir) :
		for dumbeta in betaarr : 
			val = np.pi / dumbeta
			axarr[jj].axvline( val , linestyle="--", color='k' ) 

if args.backgroundcolor :
	for iax in range(len(ax)) : 
		xy = ax[iax].axis()
		x = np.linspace( xy[0],xy[1], 20 ) 
		y = np.linspace( xy[2],xy[3], 20 ) 
		print "X : ", x
		print "Y : ", y
		cend = 0.4
		cendarr = [cend , 1.-cend, 1.-cend ]
		ctop = [ cendarr[iax] for xi in x ]
		cbot = [ 0.5 for xi in x ] 
		print "X : ", ctop
		print "Y : ", cbot
		z = [[z] * 10 for z in range(10)]
		num_bars = 100  # more bars = smoother gradient
		from matplotlib.colors import NoNorm
		axarr[iax].imshow([ctop, cbot], 
				  cmap = plt.cm.bwr,
				  interpolation = 'bicubic',
				  extent= (xy[0],xy[1],xy[2],xy[3]) ,
				  norm = NoNorm()
				  #vmin = 0., vmax = 0.5
			     )
		axarr[iax].set_aspect('auto')

if args.linearline or args.fractionline : 
	#[x0,x1,y0,y1] = axnow.axis()
	x0 = wn[0]
	xll = np.linspace( x0/10.,x0*10., num=5 ) 
	if args.fractionline : 
		mstar = float( args.fractionline.split("_")[0] )
		slope = 1. - mstar
		fract = float( args.fractionline.split("_")[1] )
		yll = np.array( -slope * np.power(xll,fract) )
		leglab = lambda m,b : r"$m^*$={} ; $\omega^{{{}}}$".format(m,b)
		if args.minus : yll = -yll
		axnow.plot( xll, yll , "k-", label=leglab(mstar,fract) ) 
	if args.linearline : 
		if len(args.linearline.split("_"))>1 : 
			nmstar = len(args.linearline.split("_"))
		else : 
			mstar = float( args.linearline ) 
			slope = 1. - mstar
			yll = np.array( -slope * xll )
		fract = 1
		leglab = lambda m,b : r"$m^*$={}".format(m)
		if args.minus : yll = -yll

		if len(args.linearline.split("_"))>1 : 
			for mstar in args.linearline.split("_") :
				slope = 1. - float(mstar)
				yll = np.array( -slope * xll )
				axnow.plot( xll, yll , "-", label=leglab(mstar,fract) ) 
		else : 
			axnow.plot( xll, yll , "c-", label=leglab(mstar,fract) ) 

Comp="Im"
if args.realpart : Comp = "Re"
xticax.minorticks_on()
if args.rew : 
	xlab = r"$\omega$"
	ylab = Comp+ylaborig +"("+xlab+")"
else :
	xlab = r"$\omega_n$"
	ylab = Comp+ylaborig +"("+xlab+")"
if args.selffactor :
	ylab = args.selffactor+ylab
if args.selfoffset :
	ylab = ylab +"+"+ args.selfoffset
ylab = minusind  + ylab
xlab = minusindx + xlab
labbasis = basis
if args.transform : labbasis = args.transform
if len(chdat)<2 : ylab = ylab + "$^{{{}}}$".format(bolabel(labbasis,j))
print "YLAB : ", ylab
ylabax.set_ylabel( ylab )
xlabax.set_xlabel( xlab )
if args.noleg :  pass
else : 
	labsp = 0.5
	if args.legendlabelspace :  labsp = float( args.legendlabelspace ) 
	if args.legendlocation : 
		legax.legend( loc=int(args.legendlocation), labelspacing=labsp  )
	else : 
		legax.legend( labelspacing=labsp )
if args.alllegend  :
	for jj in range(ndir) :
		axarr[jj].legend()
hobj = hobjarr[0]
hobj.readParameters()
if args.notitle : pass
else : 
	try :
		titleax.set_title( hobj.title + " n={:.2f}".format( getattr(hobj,"cflattu_%d"%ni) )  ) 
	except : 
		titleax.set_title( hobj.title ) 
if args.tl : plt.tight_layout()

for axnow in axarr : 
	axnow.tick_params( axis='x',  pad=8)
#alph = [ "(a)" , "(d)", "(e)" ]
#plt.text(-0.41, 3,  alph[0], horizontalalignment='left', verticalalignment='center', transform=axnow.transAxes, fontsize=fs+4 )# , weight='semibold')
#axarr[0].set_ylabel("")
#axarr[2].set_ylabel("")

argsavefigname = hobj.tdir+"/"+hobj.fignamepart()+"_"+obj
for ar in sys.argv[ndir+1:] :
        print "arg : ", ar
        argsavefigname = argsavefigname + '_' + ar.split("-")[-1]
fsave = argsavefigname
ladjust=( 0.18 if (args.ladjust is None) else float(args.ladjust) )
radjust=( 0.93 if (args.radjust is None) else float(args.radjust) )
badjust=( 0.20 if (args.badjust is None) else float(args.badjust) )
tadjust=( 0.93 if (args.tadjust is None) else float(args.tadjust) )
wspace =( 0    if (args.wspace  is None) else float(args.wspace ) )
hspace =( 0    if (args.hspace  is None) else float(args.hspace ) )
plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, wspace=wspace , hspace=hspace )
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
	pass
	#plt.show()
