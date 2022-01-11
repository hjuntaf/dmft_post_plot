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

def findgs( gptlenergyarr ) :
	ngptl = len(gptlenergyarr)
	niter = len(gptlenergyarr[0])
	e0arr	 = [999. for a in range(niter) ]
	e0secarr = [-1 for a in range(niter) ]
	for iiter in range(niter) :
		for igptl in range(ngptl) :
			try : 
				if e0arr[iiter] > float(gptlenergyarr[igptl][iiter]) : 
					e0arr[iiter]    = float(gptlenergyarr[igptl][iiter])
					e0secarr[iiter] = igptl
			except : 
				pass
	return np.array(e0arr, dtype=float), e0secarr

def findgsgptl( eAll, count, method ) :
	e=9999.
	for line in eAll :
		lcount = int(line[0])
		lmethod= line[2]
		emin   = float(line[3])
		if lcount==count and lmethod==method :
			if emin<e : 
				e = emin
		#print lcount,count , lmethod,method 
		#print lcount==count , lmethod==method 
		#print line
		#sys.exit(1)
	return e
def returnDecomp( energyallLine ) :
	count  = int( energyallLine[0] ) 
	gindex = int( energyallLine[1] ) 
	method =      energyallLine[2]
	egptlarr  = np.array(energyallLine[3:],dtype=float)
	return count, gindex, method, egptlarr
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
	dname = hobj.tdir+"/iter/energyall.txt"  
	if args.iimp : dname = hobj.tdir+"/iter/energyall_n{}.txt".format( int(args.iimp) )
	print "Reading :", dname
	energyallArr = np.genfromtxt( dname , dtype=str)
	#print "All : "
	#print energyallArr 

	#Find the number of iteractions 
	count0 = int(energyallArr[ 0][0])
	countE = int(energyallArr[-1][0])
	ncount = countE-count0+1
	print "count0, countE, ncount: ", count0,countE,ncount
	
	#Find the numbers of gindex_sector and of vectors in Davidson/Arpack method.
	nsector = 25-1 
	if args.nsector :
		nsector = int(args.nsector)
	if energyallArr.shape[0]/ncount > 150 :
		nsector=169
	if nsector > 150 : 
		gptlinfo = np.genfromtxt(homepath+"/tools/dmft/plot/gptlinfo_zeroSOC", dtype=str )
	else : 
		gptlinfo = [["",""]]*nsector
	print energyallArr[-1][2], "D" , energyallArr[-1][1], nsector 
	if args.noarpack is None : 
		if energyallArr[-1][2]=="D" and (int(energyallArr[-1][1])!=nsector if nsector==24 else nsector-1) : ncount -=1
	nvec = len(energyallArr[0][3:])
	print "nvec,nsector : ", nvec,nsector

	#Find ground state energy in each iteration
	eminarr = np.zeros( ncount )
	for ic in range(ncount) :
		count = ic + count0
		eminarr[ic] = findgsgptl( energyallArr, count, "D" )
		if args.noarpack : 
			pass
		else : 
			eminarr[ic] = findgsgptl( energyallArr, count, "A" )
	print "eminarr :\n", eminarr

	egindexdataDA = np.zeros( (ncount,nsector,nvec) )
	egindexdataA = [ [] for a in range(ncount)]
	for line in energyallArr : 
		count  = int( line[0] ) 
		gindex = int( line[1] ) 
		method =      line[2]
		egptlarr  = np.array(line[3:],dtype=float)
		if count < ncount+count0 :
			egindexdataDA[count-count0][gindex] = egptlarr
			if method=='A' :
				egindexdataA[count-count0].append( line )
		else : pass
	print "shape(egindexdataDA) :\n", np.shape(egindexdataDA)
	print "shape(egindexdataA ) :\n", np.shape(egindexdataA)

	hobj.ncount = ncount
	hobj.count0 = count0
	hobj.eminarr = eminarr
	hobj.egindexdataDA = egindexdataDA 
	hobj.egindexdataA  = egindexdataA 
	hobj.nvec = nvec

def plotalliter( hobj , axnow , xval=False , label=False ) : 
	ncount = hobj.ncount
	count0 = hobj.count0  
	eminarr= hobj.eminarr 
	nvec   = hobj.nvec
	if args.plotlast : countarr = [-1]
	else : countarr = range(ncount)
	for ic in countarr :
		count = ic+count0 if ic>-1 else  ncount-1+count0
		eicarr = hobj.egindexdataDA[ic]
		for ig in range(nsector) :
			gindex = ig
			egptlarr = eicarr[ig]
			if args.diff : egptlarr = egptlarr - hobj.eminarr[ic]
			area= 8 #if method=='D' else 12
			c   = 'C{}'.format(gindex%10)
			colors = colordiffftn(gindex%10)
			m   = markersdiff[gindex%10] #if method=='D' else "_"
			if label : lab = "gindex={}".format(gindex)+gptlinfo[ig][-2] if ic<1 else ""
			else : lab=''

			xdat = [ count for a in range(nvec) ] if xval is False else [ xval for a in range(nvec) ]
			ydat = egptlarr
			axnow.plot( xdat, ydat, m, alpha=0.4, c=colors, ms=area , mec='black', label=lab )

	for ig in range(nsector) :
		gindex = ig
		egptlarr0 = hobj.egindexdataDA[:,ig,0]
		xdat = np.arange(ncount)+count0 if xval is False else [ xval for a in range(nvec) ]
		if args.plotlast :
			egptlarr0 = hobj.egindexdataDA[-1,ig,0]
			xdat      = count0 if xval is False else xval
		if args.diff : egptlarr0 = egptlarr0 - hobj.eminarr
		if args.plotlast :
			pass
		else :
			colors = colordiffftn(gindex%10)
			ydat = egptlarr0 
			try :
				for xpt,ypt in zip(xdat,ydat) :
					axnow.text( xpt , ypt , "{}".format( gindex )+gptlinfo[ig][-2] , alpha=0.5 ) #, color='grey')
			except : 
				axnow.text( xdat , ydat , "{}".format( gindex )+gptlinfo[ig][-2] , alpha=0.5 ) #, color='grey')
			axnow.plot( xdat , ydat , '-', alpha=0.4, c=colors )
			#axnow.plot(xpt, ypt, marker=m, alpha=0.4, c=colors, ms=area , mec='black', label='gptl={}'.format(gindex) )#, c=colors)

	axnow.legend()
	# reverse the order
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles[::-1], labels[::-1])

if args.plotlast :
	edata0 = np.zeros( (nsector,ndir) ) 
	for ig in range(nsector) : 
		for idir in range(ndir) : 
			edata0[ig][idir] =  hobjarr[idir].egindexdataDA[-1,ig,0]
			if args.diff : edata0[ig][idir] -=  hobjarr[idir].eminarr[-1]
	if args.plotlabx : 
	  	xdat = []
		for idir in range(ndir) : 
			hobj=hobjarr[idir]
			xdat.append( getattr(hobj, args.plotlabx) )
	else : 
		xdat = range(ndir)
	xdat = np.array( xdat, dtype=float ) 
	for ig in range(nsector) : 
		ydat = edata0[ig]
		colors = colordiffftn(ig%10)
		axarr[0].plot( xdat , ydat , '-', alpha=0.4, c=colors )
		print ig, ydat
	labonarr = np.zeros(ndir) ; labonarr[0]=1
	if args.noleg : labonarr[0]=0
	for idir in range(ndir) : 
		axnow = axarr[idir] 
		hobj  = hobjarr[idir]
		plotalliter( hobj, axnow , xval=xdat[idir], label=labonarr[idir] )
else  : 
	for idir in range(ndir) : 
		axnow = axarr[idir] 
		hobj  = hobjarr[idir]
		plotalliter( hobj, axnow )
	

ylab = "Energy_imp"
if args.diff : ylab = r"$\Delta$" + ylab
xlab = "iteration"
if args.plotlabx : xlab = filterlabel(args.plotlabx)
axSetRange( axnow, args.xxrange , "x" )
axSetRange( axnow, args.yrange , "y" )
axnow.set_xlabel( xlab ) 
axnow.set_ylabel( ylab ) 
axnow.grid( color='grey' , linestyle=':' , linewidth=0.8 )


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

argsavefigname = hobj.tdir+"/"+hobj.fignamepart()+"_gptlE"
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
