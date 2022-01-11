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
parser.add_argument('ddir', metavar='D', type=str, nargs='+',
		                    help='Data paths' ) 
parser.add_argument("-t", "--totally", action="store_true", help="Plot total data in one panel" )
parser.add_argument("--oldsoc", action="store_true" , help="SOC strength is in a old-format and multiplied by 2./3." )
parser.add_argument("--selfyrange", "--sy", type=str, help="set yrange of self-energy e.g. _y0_y1" )
parser.add_argument("--selfxrange", "--sx", type=str, help="set xrange of self-energy e.g. _x0_x1" )
parser.add_argument("--hybyrange", "--hyy", type=str, help="set yrange of hyb-energy e.g. _y0_y1" )
parser.add_argument("--hybxrange", "--hyx", type=str, help="set xrange of hyb-energy e.g. _x0_x1" )
parser.add_argument("-y", "--yrange", type=str, help="set yrange e.g. _y0_y1" )
parser.add_argument("-x", "--xxrange", type=str, help="set xrange e.g. _x0_x1" )
parser.add_argument("-b", "--basis", type=str, help="Set basis. e.g. \"-b j | -b t2g\"" ) 
parser.add_argument("--basisarr", "--barr", "--ba", type=str, help="Set basis. e.g. \"-b j | -b t2g\"" ) 
parser.add_argument("-c", "--chdat", type=str, help="Choose a component(s) of basis. e.g. \"-c _0_1_2 (default 0,2,4)\"" ) 
parser.add_argument("-u", "--ups", action="store_true" , help="Plot (pseudo)spin-up components of basis" )
parser.add_argument("-a", "-s", "--add", action="store_true", help="Add plotting of self-energy." ) 
parser.add_argument("--realpart", "--realdata", "--rp", "--rd", action="store_true", help="Plot the real-part of the data instead of the imaginary part of self-energy.")
parser.add_argument("-o", "--one", action="store_true", help="Plot in one frame.." ) 
parser.add_argument("--savetransparent", "--savetp", "--stp", action="store_true", help="Save the fig with transparent background" )
parser.add_argument("--selfonly", "--so", action="store_true", help="Plot self-energy only." ) 
parser.add_argument("--self", "--self",     action="store_true", help="Add plotting of imag. self-energy fucntions." ) 
parser.add_argument("--hyb", "--hy",     action="store_true", help="Add plotting of imag. hybridization fucntions." ) 
parser.add_argument("--hybre", "--hyre", action="store_true", help="Add plotting of real hybridization fucntions." ) 
parser.add_argument("--iadd", action="store_true", help="Add plotting of integrated PDOS" ) 
parser.add_argument("--pdf", action="store_true", help="Save the figure in .pdf" ) 
parser.add_argument("--png", action="store_true", help="Save the figure in .png" ) 
parser.add_argument("--eps", action="store_true", help="Save the figure in .eps" ) 
parser.add_argument("--cptmp", action="store_true", help="Copy the saved figure into Dropbox/tmp_linux/") 
parser.add_argument("--notitle", action="store_true", help="Remove the title caption." )
parser.add_argument("--paper", action="store_true" , help="Plot as paper-use." ) 
parser.add_argument("--trans", "--tr", action="store_true", help="Transpose the plot" ) 
parser.add_argument("--rowsep", "--rs", action="store_true", help="Transpose the plot" ) 
parser.add_argument("--selfm", "--sm", action="store_true", help="Plot minus self-energy" ) 
parser.add_argument("--minusx", "--mx", action="store_true", help="Plot minus x-data.")
parser.add_argument("--dxtick", action="store_true", help="Double the xtick increment." ) 
parser.add_argument("--dytick", action="store_true", help="Double the ytick increment." ) 
parser.add_argument("--xtick", "--xt", type=str, help="Set xtick as 'x0_x1_x2'.")
parser.add_argument("--ytick", "--yt", type=str, help="Set ytick as 'x0_x1_x2'.")
parser.add_argument("--total", "--tt", action="store_true", help="Plot sum of densities as a total density." ) 
parser.add_argument("--totalcolor", "--ttc", type=str, help="Specify color of the total density." ) 
parser.add_argument("-v", "--vert", action="store_true", help="Plot vertically for ppt-use." ) 
parser.add_argument("--notick", action="store_true", help="Remove ticks." ) 
parser.add_argument("--noytick", "--noytic", action="store_true", help="Remove y-ticks." ) 
parser.add_argument("--noxtick", "--noxtic", action="store_true", help="Remove x-ticks." ) 
parser.add_argument("--noylabel", "--noylab", action="store_true", help="Remove y-label." ) 
parser.add_argument("--noxlabel", "--noxlab", action="store_true", help="Remove x-label." ) 
parser.add_argument("--yoff", action="store_true", help="Remove xticks." ) 
parser.add_argument("--fontsize", "--fonts", "--fs", type=float, help="Set fontsize." ) 
parser.add_argument("--nplot", type=int, help="Set Nplot." ) 
parser.add_argument("--fillb", "--fb", action="store_true", help="Fill the plotting of blue one." ) 
parser.add_argument("--fill", "--fl", action="store_true", help="Fill the plotting of data." ) 
parser.add_argument("--fillcomp", "--flc", type=str, help="Fill the plotting of data according to the specified component(s)." ) 
parser.add_argument("--fillx", "--flx", action="store_true", help="Use fill_betweenx()")
parser.add_argument("--figs", "--figsize", type=str, help="Set the size of figure." ) 
parser.add_argument("--ladjust", "--lad", type=str, help="Set   left-adjust value." ) 
parser.add_argument("--radjust", "--rad", type=str, help="Set  right-adjust value." ) 
parser.add_argument("--badjust", "--bad", type=str, help="Set bottom-adjust value." ) 
parser.add_argument("--tadjust", "--tad", type=str, help="Set    top-adjust value." ) 
parser.add_argument("--wspace",  "--wsp", type=str, help="Set wspace-adjust value." ) 
parser.add_argument("--hspace",  "--hsp", type=str, help="Set hspace-adjust value." ) 
parser.add_argument("--noleg", "--nol",  action="store_true", help="Remove the legend.") 
parser.add_argument("--leglabel", "--legl", "--leglabl", "--leglab", type=str, help="Set the labels of the legend" ) 
parser.add_argument("--leglabelcomp", "--leglc", "--lglc", action="store_true", help="Add the labels of component" ) 
parser.add_argument("--leglabelarr", "--leglabarr", "--legla", type=str, nargs='+', help="Set the labels of the legend" ) 
parser.add_argument("--leglabeltail", "--legltail", "--leglt", type=str, help="Add tails to the labels of the legend as parameters." ) 
parser.add_argument("--legncol", "--legnc", "--legendncol", type=str, help="Set the number of columns in the legend.")
parser.add_argument("--alllegend", "--allleg", "--aleg", "--legall", action="store_true" , help="Draw all legends.")
parser.add_argument("--sublabel", "--slab", "--sublab", action="store_true", help="Set the labels of sub-figures with alphabet." ) 
parser.add_argument("--sublabelparameter", "--slabpar", "--slabp", "--sublabpar", nargs='+', type=str, help="Set a parameter of the sub-labels.")
parser.add_argument("--sublabelnoalphabet", "--slabnoalpha", "--sublabnoalpha", action="store_true", help="Set the labels of sub-figures without alphabet." ) 
parser.add_argument("--insetdos", "--insetd", action="store_true", help="Draw an inset of the PDOS." ) 
parser.add_argument("--xrangein", "--xrangeinsetlog", "--xin", type=str, help="Set xrange of the inset-dos plot.") 
parser.add_argument("--yrangein", "--yrangeinsetlog", "--yin", type=str, help="Set yrange of the inset-dos plot.") 
parser.add_argument("--sinv", "--selfinv", action="store_true", help="Plot inverse of self-energies.") 
parser.add_argument("--legloc", type=str, help="Specify the location of logened.")
parser.add_argument("--line", "--lines", action="store_true" , help="Draw only with lines.") 
parser.add_argument("--linestylearr", "--lstylearr", "--lstylea", "--lsa", "--lstya", nargs='+', type=str, help="Specify linestyles.")
parser.add_argument("--linestylecompoarr", "--lstylecomparr", "--lstylecompa", "--lsca", nargs='+', type=str, help="Specify linestyles.")
parser.add_argument("--solidline", "--solidl", action="store_true", help="Plot the data as solidlines.")
parser.add_argument("--dashedline", "--dashl", action="store_true", help="Plot the data as dashedlines.")
parser.add_argument("--diffcolordir", "--diffcdir", "--diffc", action="store_true", help="Plot the data with different color for different directories")
parser.add_argument("--diffcolororb", "--diffcorb", "--diffco", action="store_true", help="Plot the data with different color for different orbitals")
parser.add_argument("--diffcolororbcomp", "--diffcorbcomp", "--diffcocomp", nargs="+", type=str, help="Plot the data with different color for different orbitals as specified.")
parser.add_argument("--epsilon", "--ep", type=str, help="Specify epsilon of data. (default=0.03)")
parser.add_argument("--epsilontag", "--ept", "--eptag", action="store_true",  help="Set Eipsilon tag in file-name." )
parser.add_argument("--epsilontagfloat", "--eptf", "--eptagf", action="store_true",  help="Set Eipsilon tag in file-name." )
parser.add_argument("--broadening", "--br", type=str, help="Set Broadening value for different data set" )
parser.add_argument("--broadeningtag", "--brt", "--brtag", action="store_true", help="Set Broadening tag in file-name." )
parser.add_argument("--setcount", "--count", "--sc", type=str, help="Set the number of iteration-index ('count') of the 'Dir' file." )
parser.add_argument("--ylabelcoord", "--ylabcoord", "--ylabc", "--ylc", type=str, help="Set the x-coordinate of y-axis label.")
parser.add_argument("--ylabel", "--yl", "--ylab", type=str, help="Set the y-axis label.")
parser.add_argument("--vline", "--vl", type=str, help="Plot a vertical line of the specified value on x-axis.")
parser.add_argument("--vline2", "--vl2", type=str, help="Plot a vertical line of the specified value on x-axis.")
parser.add_argument("--fermiline", "--fline", action="store_true", help="Draw a line in the Fermi level.")
parser.add_argument("--cutoff", "--coff", type=str, help="Set lower-bound cutoff of x-data to be plotted .") 
parser.add_argument("-l", "--llog", action="store_true" , help="Plot data in log-log plot." ) 
parser.add_argument("--nimp", "--ni", "--nimpurity", type=str, help="Set the number of impurity sites.")
parser.add_argument("--findpeak", "--fp", "--fpeak", "--findp", nargs="+", type=float, help="Find a peak with quadratic-fit with specified range. e.g. --fp a b ")
parser.add_argument("--ndatalinearfit", "--nlfit", "--nlf", type=str, help="Specify the number of data to be used in fitting (default=3).")
parser.add_argument("--fitoffset", "--fitoff", "--foff", type=str, help="Set the offset in indices of the data.")
parser.add_argument("--grid", "--gr", action='store_true', help="Draw grid-lines.")
parser.add_argument("--legequation", "--legeq", "--legendeq", "--legendequation", action='store_true', help="Use legend with equation mode.")
parser.add_argument("--chooselegax", "--clegax", "--chooselax", "--clax", type=str, help="Choose one of ax's indices you want to put the legend into.")
parser.add_argument("--write", "--wr", action="store_true", help="Write the last part of data in 'tmpxdat' and 'tmpydat'.")
args = parser.parse_args()

from matplotlib.patches import Ellipse
import matplotlib.pylab as pylab
figs=5.5	#6
figsrat = 0.5	#0.416667
dfigsrat= 0.5	#0.416667
if args.add   :  figsrat += dfigsrat
if args.hyb   :  figsrat += dfigsrat
if args.hybre :  figsrat += dfigsrat
figsy = figs*figsrat
if args.trans :
	figsy=figs ; figs=figsy*figsrat
elif args.trans :
	figs=2*figs ; figsy = figs*figsrat
fs=12 
if args.fontsize : fs = int( args.fontsize ) 
ms=6
if args.figs :
	figsarr = args.figs.split("_")
	figs    = float(figsarr[0])
	figsy   = float(figsarr[1])
print "figsize = ({}, {})".format( figs, figsy ) 
params = {      'legend.fontsize'       : fs ,
		'figure.figsize'        : (figs,figsy) ,
		'axes.labelsize'        : 'x-large' ,
		'axes.titlesize'        : 'small' ,
		'axes.linewidth'        : 1.5 ,
		'xtick.labelsize'       : 'x-large' ,
		'ytick.labelsize'       : 'x-large' ,
		'lines.markersize'      : ms ,
		#'markers.fillstyle'     : 'none' ,
		'font.size'             : fs }
pylab.rcParams.update(params)

AdjustSpaceDefault = [ 0.15 , 0.93 , 0.15 , 0.95 , 0    , 0    ]

print "\n"
nbasis = 6
if args.basis : 
	if args.basis.find("j")>-1	:
		basis = args.basis
		fname = "ldos"
	elif args.basis.find("t")>-1	:
		basis = args.basis
		fname = "tldos"
	elif args.basis.find("e")>-1	:
		basis = args.basis
		nbasis = 4 
		fname = "ldos"
	if args.basis.find("z")>-1	:
		fname = "ldos"
else : 
	#default
	basis = "t2g"
	fname = "tldos"
print "nbasis : ", nbasis

chdat = range(nbasis)
if args.ups : 
	chdat = [ 0, 2, 4 ]
	if basis.find("j")>-1 : 
		chdat = [ 0, 1, 4 ]
	if basis.find("e")>-1 : 
		chdat = [ 0, 2 ]
if args.chdat : 
	chdat = np.array( args.chdat.split("_") , dtype=int ) 

logplot=False
if args.llog :	logplot=True

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
	dat = args.ddir[jj]
	if args.epsilon : setattr( hobj, "epsilon", float(args.epsilon) ) 
	else : setattr( hobj, "epsilon", 0.03 )

	epsilon = hobj.epsilon
	if args.epsilon :
		epsilon = float(args.epsilon)
		if jj<1 :
			for hobj in hobjarr : hobj.title = "$\eta_s$={} ".format(epsilon)  + hobj.title 
	for hobj in hobjarr : setattr( hobj , "epsilon", float(epsilon) ) 
	print "epsilon ", jj, " : ", epsilon


	broadening = 0.02
	if args.broadening :
		broadening = float(args.broadening)
		if jj<1 :
			for hobj in hobjarr : hobj.title = " $\eta_d$={} ".format(broadening)  + hobj.title  

	epsilontag	= ""
	broadeningtag	= ""
	if args.epsilontag :
		epsilontag = "ep{:g}_".format(epsilon)
	if args.epsilontagfloat :
		epsilontag = "ep{:.3f}_".format(epsilon)
	if args.broadeningtag :
		broadeningtag = "br{:g}_".format(broadening)
	tagEpBr =  epsilontag + broadeningtag
	

	if args.basisarr  :
		basis = args.basisarr.split("_")[jj]
		chdat = chdatFtn(basis)
		fname = fnameDOSFtn(basis)
		print "BASIS : ", basis
		if args.chdat : 
			chdat = np.array( args.chdat.split("_") , dtype=int ) 


	wn   	= [""]*  hobj.Nplot
	coldata = [""]*( hobj.Nplot )
	if args.nplot : hobj.Nplot = int(args.nplot) 
	#hobj.title = ""
	for i in [1] :
		#if args.basis and args.basis.find("z")>-1 : 
		#	datf= "{}/lattice/u{:.3f}/{}_ep{:.3f}_{}th.dat".format(dat, hobj.UF, fname, hobj.epsilon, hobj.count )
		#	print i, " :", dat
		#	print i, "f:", datf,
		#	datarr = np.genfromtxt( datf, dtype=float ) 
		#	print "reading-successed."
		#	wn      = datarr.transpose()[0]
		#	coldata = datarr.transpose()[4:]
		#	coldata = coldata.transpose()
		#else : 
		try : 
			icount = int(args.setcount) if args.setcount else hobj.count
			datf= "{}/lattice/u{:.3f}/{}_Nplot{}_{}{}th.dat".format(dat, hobj.UF, fname, hobj.Nplot, tagEpBr, icount  )
			if args.epsilontag or args.broadeningtag : 
				datf= "{}/lattice/u{:.3f}/{}_Nplot{}_{}{}th.dat".format(dat, hobj.UF, fname, hobj.Nplot, tagEpBr, icount  )
			print i, "f:", datf,
			datarr = np.genfromtxt( datf, dtype=float ) 
			print "reading-successed."
		except :
			print "reading-faild."
			datf= "{}/lattice/{}_u{:.2f}_Nplot{}.dat".format(dat, fname, hobj.UF, hobj.Nplot )
			print i, " :", dat
			print i, "f:", datf,
			datarr = np.genfromtxt( datf, dtype=float ) 
			print "reading-successed."
		wn      = datarr.transpose()[0]
		coldata = datarr.transpose()[1:]
		coldata = coldata.transpose()

	wn   = np.array( wn   )
	print len(wn ), "points (wn)"
	coldata = np.array( coldata )
	coldata = coldata.transpose()
	print len(coldata), "components (basis vectors + sum )"
	hobj.wn = wn
	hobj.coldata = coldata

if args.oldsoc :
	S = 2./3 * S
#figs = 8
#figsy= 8
#fs = figs*1.4 ; fsm = figs*1.4
#if args.fonts :  
#	fs = args.fonts
#	fsm= args.fonts*0.5
#setattr( args , "figs", figs )
#setattr( args , "figsy", figsy )
#setattr( args , "fs", fs )
#setattr( args , "fsm", fsm )
def totally( ffig, aax , laboverw=False , labeladdtail=False , dashes=(None,None) , hobj=False, cs='', csarr=False, lstyleArr=None  ) :
	jj=0 
	sumdos = np.array( coldata[0] ) 
	sumdos = sumdos * 0 
	for ich in range(len(chdat)) :
		j	= chdat[ich]
		ldos = coldata[j]
		ldos = np.array( ldos )
		sumdos += ldos 
		
		#lc = [ "-" ]*2 + [ ":" ] *2 + [ "-" ] *2 
		if laboverw :	labelstr = laboverw[j]
		else :		labelstr = bolabel(basis,j)
		if args.legequation : 
			labelstr = r'${}$'.format(labelstr)
		if labeladdtail : labelstr += labeladdtail
		print "LAB ", j , " : ", labelstr
		xdat = wn
		if args.write : 
			dum = np.array( [ xdat, ldos] ).transpose()
			np.savetxt("tmpdata_%d"%j, dum ) 
		print "basis : ", basis
		if lstyleArr :
			try : 
				print "Use lstyleArr : ", lstyleArr
				print "shape : ",  np.shape(lstyleArr) 
				if np.shape(lstyleArr)[0]>1 : 
					lstyle = lstyleArr[ich]
			except : 
				print "Not-Use lstyleArr : "
				lstyle = lstyleArr
				pass
		else:
			print "blt-Use lstyleArr : "
			lstyle	= blt(basis,j) 
		if csarr is not False	:
			#print "csarr : ", csarr
			cs	= csarr[int( (j%hobj.tNC) )]
		print "cs ls : ", cs, lstyle
		if args.minusx :
			xdat = -1.*wn
		if args.trans : 
			if args.one :
				aax.plot( ldos, xdat, cs+lstyle , label=labelstr , lw=1 , dashes=dashes )
			else :
				aax.plot( ldos, xdat, cs+lstyle , label=labelstr , lw=1 )
			aax.set_xlabel(r'PDOS' )#, fontsize=fs )
			aax.set_ylabel(r'$E$ (eV)')#, fontsize=fs )
			if args.fillb :
				if j==1 :
					if len(cs) < 1 : 
						cs = 'C%d'%dd
					aax.fill( ldos, xdat, cs+lstyle ,alpha=0.3 ) 
				if j==0 :
					aax.fill_betweenx( xdat, 0, ldos, facecolor='C%d'%j ,alpha=0.3 ) 
		else : 
			if args.llog : 
				print "BLT ", j , basis, cs+lstyle
				if args.one : 
					aax.loglog( xdat, ldos, cs+lstyle , label=labelstr , dashes=dashes )
				else :
					aax.loglog( xdat, ldos, cs+lstyle , label=labelstr )
			else : 
				if args.one : 
					aax.plot( xdat, ldos, cs+lstyle , label=labelstr , dashes=dashes )
				else :
					aax.plot( xdat, ldos, cs+lstyle , label=labelstr )
					#aax.plot( xdat, ldos, "C0-" , linewidth=2.5 ,  label=labelstr )
			#aax.set_ylabel(r'PDOS' )#, fontsize=fs )
			aax.set_xlabel(r'$\omega$ [eV]')#, fontsize=fs )
			if args.fill :
				aax.fill( xdat, ldos, cs ,alpha=0.2 ) 
				#if (basis.find("j")>-1 and j==1) or (basis.find("t")>-1 and (j==0 or j==1 or j==3 or j==5) ) : 
				#	aax.fill( xdat, ldos, cs ,alpha=0.2 ) 
			if j==0 :
				if args.fillb : aax.fill_between( xdat, ldos, facecolor='C%d'%dd, linestyle=lstyle ,alpha=0.3 ) 
			if args.fillcomp : 
				if args.fillcomp.find("{}".format(j))>-1  :
				 	carr = { 0 : 'xkcd:lightblue',  2 : 'xkcd:orange' }
				 	carr = { 0 : 'xkcd:lightblue',  2 : 'tab:green' }
				 	aarr = { 0 : 0.2 , 2 : 0.2 }
					alph = aarr[j]
					color= carr[j]
					aax.fill( xdat, ldos, color ,alpha=alph ) 
		#aax.axvline( 0 , color='gray', linestyle="--", lw=0.5 ) 
		#aax.axhline( 0 , color='gray', linestyle="--", lw=0.5 ) 
		#aax.legend()# fontsize=fsm*0.8 ) 
		if args.findpeak : 
			w0, w1 = args.findpeak 
			print "finding a peak in [", w0,',',w1,"]"
			ind1	= xdat>w0	; xdat1 = xdat[ind1] 	; ldos1 = ldos[ind1]	
			ind2	= xdat1<w1	; xdat2 = xdat1[ind2]	; ldos2 = ldos1[ind2]	

			from scipy.optimize import curve_fit
			def funclin( x, a, b ) : return a*x 
			typefit = 'quad'
			nfit = len(xdat2)
			nresult = 50
			print "Entering linearfit ..."
			if args.ndatalinearfit :
				nfit = int(args.ndatalinearfit)
				if nfit > nresult  : 
					nresult = nfit*2
				print "nfit : ", nfit 
			iinit = 0 
			if args.fitoffset : iinit += int(args.fitoffset)
			xraw = xdat2[iinit:iinit+nfit] 
       			yraw = ldos2[iinit:iinit+nfit] 
			print "FIT_xraw : " , xraw[ [0,1,-2,-1] ]
			print "FIT_yraw : " , yraw[ [0,1,-2,-1] ]

			formfunc= 'ax^2+bx+c'
			nvarfunc= 3
			def funcpoly2( x, a, b, c ) : return a*x*x + b*x + c
			fitfunc = funcpoly2

			initvararr = [99999]*nvarfunc
			popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=(-99999 , initvararr) )
			if args.fitoffset : 
				xfit = np.linspace( xdat[0],  xraw[-1] , nresult ) 
			else : 
				xfit = np.linspace( xraw[0],  xraw[-1] , nresult ) 
			yfit = fitfunc(xfit,*popt) 
			print "XFIT : ", xfit[ [0,1,-2-1] ]
			print "YFIT : ", yfit[ [0,1,-2-1] ]
			print "FIT({};{}) : ".format( formfunc, typefit ),
			for jj in range(nvarfunc) :
				print "{}".format(popt[jj]),
			print ""
			fitls = "k--"
			aax.plot( xfit, yfit, fitls , lw=1 )

			ldosmax		= yfit.max()
			indmax		= yfit == ldosmax 
			xldosmax	= xfit[indmax]
			xldosmaxfit	= -popt[1]/2./popt[0] 	# peak  at x=-b/2a
			print "xldosmaxfit, xldosmax, ldosmax  : ", xldosmaxfit, xldosmax, ldosmax 
			if basis.find("t")>-1 and j==0 and hobj :
				fname 	= hobj.tdir+"/result/u{:.3f}/vHSpeak_fit_t2g.dat".format(hobj.UF)
			elif basis.find("j")>-1 and j==1 and hobj :
				fname 	= hobj.tdir+"/result/u{:.3f}/vHSpeak_fit.dat".format(hobj.UF)
			if (basis.find("j")>-1 and j==1 and hobj) or ( basis.find("t")>-1 and j==0 and hobj ) :
				with open( fname, 'w' ) as fwfit :
					icount = int(args.setcount) if args.setcount else hobj.count
					fwfit.write( "{:3d} ".format(icount) )
					fwfit.write( "{:12f} ".format(xldosmaxfit) ) 
					fwfit.write( "\n" )
				
		jj+=1
	#plt.title("U={} J={} S={:.2f} D={} ; nb={}".format(UF,J,S,D,Nb))
	if args.total :
		if args.trans : 
			if args.one :
				aax.plot( sumdos, xdat, "--" , label=bolabel("total",0) , lw=1 , dashes=dashes )
			else :
				aax.plot( sumdos, xdat, "--" , label=bolabel("total",0) , lw=1 )
			aax.set_xlabel(r'PDOS' )#, fontsize=fs )
			aax.set_ylabel(r'$E$ (eV)')#, fontsize=fs )
		else : 
			if args.one : 
				aax.plot( xdat, sumdos, "--" , label=bolabel("total",0) , dashes=dashes )
			else :
				aax.plot( xdat, sumdos, "--" , label=bolabel("total",0) )
			aax.set_xlabel(r'$\omega$')#, fontsize=fs )
			if args.fillcomp : 
				if args.fillcomp.find("t")>-1  :
					aax.fill( xdat, sumdos, 'C%d'%jj ,alpha=0.2 ) 
	aax.tick_params( axis="both", direction="in" )#, labelsize=fs , length=10, width=2 ) 
	aax.tick_params( axis="x", bottom=True, labelbottom=True )
	#if args.trans : 
	#	aax.set_xlim( 0, 1.2 ) 
	#else :
	#	aax.set_ylim( -0.05,1.05 ) 
	if args.xxrange : 
		axSetRange( aax , args.xxrange, 'x' ) 
	if args.yrange : 
		axSetRange( aax , args.yrange,  'y' ) 
	if args.dxtick :
		aax.set_xticks(aax.get_xticks()[0::2])
	if args.dytick :
		aax.set_yticks(aax.get_yticks()[1::2])
def totallyi( ffig, aax , laboverw=False , labeladdtail=False , dashes=(None,None) , hobj=False, cs='', csarr=False, lstyleArr=None, factor=None  ) :
	dwn = float( wn[1] )  - float( wn[0]  )
	sldosall = [0]*len(coldata[0]) 
	xdat = wn
	if args.minusx :
		xdat = -wn 
	for j in chdat :
		ldos = coldata[j]
		ldos = np.array( ldos , dtype=float )
		ildos = np.zeros( len(ldos) , dtype=float )
		sldos = 0. 
		for jj in range( len(ldos) )  :
			sldos += ldos[jj]
			ildos[jj] = sldos * dwn
			sldosall[jj] += sldos * dwn
		#plt.figure()
		lc = [ "-" ]*2 + [ ":" ] *2 + [ "-" ] *2 
		if factor is not None :  ildos = ildos*factor
		aax.plot( xdat, ildos, blt(basis,j) , label=blabel(basis , j) )
		aax.axvline( 0 , color='gray', linestyle="--", lw=0.5 ) 
		aax.axhline( 0 , color='gray', linestyle="--", lw=0.5 ) 
	#plt.title("U={} J={} S={:.2f} D={} ; nb={}".format(UF,J,S,D,Nb))
	sldosall = np.array( sldosall, dtype=float) * 2
	intfill = [0] * (nbasis+1)
	j=0
	for jj in range( len(sldosall)-1 )  :
		if sldosall[jj+1] > j and sldosall[jj] < j :
			intfill[j] = jj 
			aax.axvline( xdat[jj] , color='gray' , lw = 0.5 , ls='--' ) 
			print j, "; ", jj , " : ", sldosall[jj], sldosall[jj+1]
			j += 1
				
	aax.plot( xdat, -sldosall,		'k-' , label="-total" )
	aax.plot( xdat, -sldosall/nbasis,	'k:' , label="-total/{}".format(nbasis) )
	aax.axvline( 0 , color='gray', linestyle="--", lw=0.5 ) 
	aax.axhline( 0 , color='gray', linestyle="--", lw=0.5 ) 
	#aax.set_ylabel(r'sum(PDOS)' ) 
	aax.set_xlabel(r'$\omega$')
	#aax.legend()
	#aax.set_ylim( -1.05,1.05 ) 
	#aax.set_yticks( range(7) )
	axSetRange( aax , args.xxrange, 'x' ) 
	axSetRange( aax , args.yrange,  'y' ) 
def separate(dd, ax, ini=0 ) :
	xdat = wn
	if args.minusx :
		xdat = -wn 
	print ini, "-th impurity ::"
	for j in range(nbasis/2) :

		lab  = updnplotlabel(j,basis)
		ldos = hobj.coldata[ lab[0] +ini*nbasis ]
		ldos2= hobj.coldata[ lab[1] +ini*nbasis ]
		ldos = np.array( ldos )
		ldos2= np.array( ldos2)
		if args.write : 
			dum = np.array( [ xdat, ldos, ldos2 ] ).transpose()
			np.savetxt("tmpdata_%d"%j, dum ) 

		ax[j].plot( xdat, ldos , "C0-" ,  label=blabel(basis,lab[0])  )
		ax[j].plot( xdat, ldos2, "C3-" , label=blabel(basis,lab[1])  )
		#ax[j].axvline( 0 , color='gray', linestyle="--", lw=0.5 ) 
		#ax[j].axhline( 0 , color='gray', linestyle="--", lw=0.5 ) 
		#ax[j].set_ylabel('PDOS' + lab[2] )
		axSetRange( ax[j] , args.xxrange, 'x' ) 
		axSetRange( ax[j] , args.yrange,  'y' ) 
		kwargs_leg  =   dict()
		#if args.legncol : 
		#	leg_ncol   = int(args.legncol)
		#	kwargs_leg.update( {"ncol":leg_ncol} )
		if args.legloc : 
			loc = int(args.legloc)
			kwargs_leg.update( {"loc":loc} )
		ax[j].legend(**kwargs_leg)
		if args.fill :
			#ax[j].fill( xdat, ldos,  "C0" ,alpha=0.2 ) 
			#ax[j].fill( xdat, ldos2, "C3" ,alpha=0.2 ) 
			ax[j].fill_between( xdat, ldos,  xdat*0, color="C0" ,alpha=0.2 ) 
			ax[j].fill_between( xdat, ldos2, xdat*0, color="C3" ,alpha=0.2 ) 
	#ax[0].set_title("U={} J={} S={:.2f} D={} ; nb={}".format(UF,J,S,D,Nb))
	titleax = ax[1]
	plt.xlabel(r'$\omega$')
	#titleax.set_ylabel('PDOS' + lab[2] )
	titleax.set_ylabel('PDOS')
	x1,x2,y1,y2 = ax[-1].axis()
	print "ax.axis() : ", x1, ',', x2, ',', y1, ',', y2
	if args.yrange :
		print args.yrange.split('_')
		ax[-1].set_ylim( float( args.yrange.split('_')[-2] ) , float( args.yrange.split('_')[-1] ) ) 
def total(ax, dashes=(None,None), label=None ) :
	xdat = wn
	if args.minusx :
		xdat = -wn 
	#if dd<1 : fig, ax = plt.subplots()
	ldossum= np.zeros( len(xdat) ) 
	for j in range(nbasis/2) :

		lab  = updnplotlabel(j,basis)
		ldos = coldata[ lab[0] ]
		ldos2= coldata[ lab[1] ]
		ldos = np.array( ldos , dtype=float )
		ldos2= np.array( ldos2, dtype=float )

		ldossum = ldossum + ldos
		ldossum = ldossum + ldos2

	ax.plot( xdat, ldossum , "-" ,  label=label ,dashes=dashes  )
	#ax.set_ylabel('DOS' + "$_{tot}$" )
	ax.axvline( 0 , color='gray', linestyle="--", lw=0.5 ) 
	ax.axhline( 0 , color='gray', linestyle="--", lw=0.5 ) 
	axSetRange( ax , args.xxrange, 'x' ) 
	axSetRange( ax , args.yrange,  'y' ) 
	kwargs_leg  =   dict()
	#if args.legncol : 
	#	leg_ncol   = int(args.legncol)
	#	kwargs_leg.update( {"ncol":leg_ncol} )
	#if args.legloc : 
	#	loc = int(args.legloc)
	#	kwargs_leg.update( {"loc":loc} )
	ax.legend(**kwargs_leg) 
	#ax[0].set_title("U={} J={} S={:.2f} D={} ; nb={}".format(UF,J,S,D,Nb))
	#titleax = ax 
	plt.xlabel(r'$\omega$')
	x1,x2,y1,y2 = ax.axis()
	print "ax.axis() : ", x1, ',', x2, ',', y1, ',', y2
	if args.yrange :
		print args.yrange.split('_')
		ax.set_ylim( float( args.yrange.split('_')[-2] ) , float( args.yrange.split('_')[-1] ) ) 
	#return [fig, ax]

print "ndir : ", ndir
nax	= ndir
for dd in range(ndir)  :
	hobj = hobjarr[dd] 
	wn = hobj.wn
	coldata = hobj.coldata
	if args.leglabel or args.leglabelarr :	
		basisCalc = basis
		#if args.transform :	
		#	if args.transform.find("t")>-1 :	basisCalc='t'
		#	elif args.transform.find("j")>-1 :	basisCalc='j'
		nchdat=len(chdat)
		legparArr	= args.leglabelarr if args.leglabelarr else [args.leglabel]
		labsuffix	= ''
		labsuffixSep	= '' if len(legparArr)==1 else ', '
		for legpar in legparArr :
			if legpar.find("Treal")>-1 : 
				if hobj.Treal<1e-5 : 
					labsuffix += r"$T = $"+"{:.0f} K ".format(hobj.Treal)
				else : 
					labsuffix += r"$T = $"+"{:.0f} K ".format(hobj.Treal)
				leglabel = [labsuffix]*nbasis
			elif legpar.find("nOcculattAra")>-1 or legpar.find("tfilling")>-1 :
				if args.paper : 
					labsuffix += filterlabel(legpar,paper=args.paper) + " = {:.2f}".format( float(getattr( hobj, legpar )) )
				else : 
					labsuffix += r"$n$"+" = {:.2f}".format( float(getattr(hobj,legpar)) )
				leglabel = [labsuffix]*nbasis
			elif args.paper : 
				labsuffix += filterlabel(legpar) + " = {}".format( getattr( hobj, legpar ) ) 
				leglabel = [labsuffix]*nbasis
			else : 
				#labsuffix = "${}={}$".format(legpar ,  getattr( hobj, legpar ) ) 
				labsuffix += filterlabel(legpar) + " = {}".format( getattr( hobj, legpar ) , paper=args.paper )
				if legpar.find("filling")>-1 : 
					labsuffix += filterlabel(legpar) + " = {:.2f}".format( float(getattr( hobj, legpar )) )
				leglabel = [labsuffix]*nbasis
			labsuffix += labsuffixSep
		if args.leglabelcomp : 
			leglabel = [ leglabel[j] +"; "+ bolabelArr(basisCalc)[j] for j in range(nbasis) ]
	else : 		leglabel = False
	if args.leglabel and args.leglabeltail :
		leglabeltail	= leglabel 
		leglabel	= False
	else :
		leglabeltail	= False
	print "leglabel :: ", leglabel 

	if args.basisarr  :
		basis = args.basisarr.split("_")[dd]
		chdat = chdatFtn(basis)
		fname = fnameDOSFtn(basis)
		print "BASIS : ", basis
		if args.chdat : 
			chdat = np.array( args.chdat.split("_") , dtype=int ) 
	if args.totally :
		print "Plotting with args.totally()."
		if args.add :
			naxdef	= 2 
			nax	= 2 
			if args.hyb : nax += 1 ; iaxhybim = nax-1
			if args.hybre : nax += 1 ; iaxhybre = nax-1
			naxrow	= nax 
			naxcol	= 1
			if args.trans : 
				naxrow	= 1
				naxcol	= nax
				if dd<1 : fig, ax = plt.subplots(naxrow, naxcol, sharey=True )
			else :
				if dd<1 : fig, ax = plt.subplots(naxrow, naxcol, sharex=True )
			if nax>1 : axarr =ax 
			else :	   axarr =[ax]*ndir
			#if dd==ndir-1 : args.fillb=None #True
			#else : args.fillb=None
			if args.dashedline :
				print "dashed-ON"
				if dd<1 : 
					dashpt = (None,None)
				else : 
					dashpt = [dd+1,1]
			else : dashpt = (None,None)
			axdos , axself = [ ax[0], ax[1] ]
			titleax = axdos
			legax   = axself
			totally( fig, axdos, laboverw=leglabel, dashes=dashpt )
			selfeobj = datobj()
			selfeobj.getplotselfw( axself , hobj.ddir , basis, nbasis, chdat, hobj , args , idir=dd , laboverw=leglabel , dashes=dashpt , logplot=logplot ) 
			if args.hyb : 
				selfeobj.getplothyb(   ax[iaxhybim]   , hobj.ddir , basis, nbasis, chdat, hobj , args ) 
			if args.hybre : 
				selfeobj.getplothybre( ax[iaxhybre] , hobj.ddir , basis, nbasis, chdat, hobj , args ) 
			#if args.selfyrange :
				#print args.selfyrange.split('_')
				#ax[1].set_ylim( float( args.selfyrange.split('_')[-2] ) , float( args.selfyrange.split('_')[-1] ) ) 
			if args.selfxrange : axSetRange( axself, args.selfxrange, "x" )
			if args.selfyrange : axSetRange( axself, args.selfyrange, "y" )
			#axdos.axvline(  0 , color='gray' , lw = 0.5 , ls='--' ) 
			#axself.axvline( 0 , color='gray' , lw = 0.5 , ls='--' ) 
			if args.vline : 
				redline = float(args.vline) # 0.048594999999998834
				axdos.axvline(  redline , color='gray' , lw = 1 , ls='-.' ) 
				axself.axvline( redline , color='gray' , lw = 1 , ls='-.' ) 
				axdos.axvline( -redline , color='gray' , lw = 1 , ls='-.' ) 
				axself.axvline(-redline , color='gray' , lw = 1 , ls='-.' ) 
			if args.vline2 : 
				redlineL, redlineR = np.array( args.vline2.split("_") , dtype=float )
				axdos.axvline(  redlineR, color='gray' , lw = 1 , ls='-.' ) 
				axself.axvline( redlineR, color='gray' , lw = 1 , ls='-.' ) 
				axdos.axvline( -redlineL, color='gray' , lw = 1 , ls='-.' ) 
				axself.axvline(-redlineL, color='gray' , lw = 1 , ls='-.' ) 
			#ax[1].autoscale(enable=True, axis='x', tight=False)
			if args.trans : 
				axdos.set_ylabel(r'PDOS')#, fontsize=fs )
				axself.set_ylabel(r'Im$\Sigma(\omega)$')#, fontsize=fs )
				axdos.set_xlabel(r'$\omega$ (eV)')
				axself.set_xlabel(r'$\omega$ (eV)')
				axdos.yaxis.set_label_position("right")

			else : 
				axdos.tick_params( axis="x", labelbottom=False )#, labelsize = args.fsm ) 
				axself.tick_params( axis="x", labelbottom=True )#, labelsize = args.fsm ) 
				#axself.tick_params( axis="y", labelright='on' )#, labelsize = args.fsm ) 
				#axdos.tick_params( axis="y", labelright='on' )#, labelsize = args.fsm ) 
				#axself.tick_params( axis="y", labelleft='off' )#, labelsize = args.fsm ) 
				#axdos.tick_params( axis="y", labelleft='off' )#, labelsize = args.fsm ) 
				pass
				#axdos.set_xlabel(r'PDOS')#, fontsize=fs )
				#axself.set_xlabel(r'Im$\Sigma(\omega)$')#, fontsize=fs )
				#axdos.set_ylabel('')
				#axself.set_ylabel(r'$\omega$ (eV)')

			if args.insetdos :
			        axinset = fig.add_axes([0.76,0.27,0.18,0.18])
			        axinset.xaxis.set_tick_params(labelsize=fs)
			        axinset.yaxis.set_tick_params(labelsize=fs)
			        axinset.tick_params( axis="both", which='both', direction="in" )#, length=10, width=2 ) 
			        for j in chdat :
					totally( fig, axinset , laboverw='' )
			        if args.xrangein : axSetRange( axinset, args.xrangein, "x" )
			        if args.yrangein : axSetRange( axinset, args.yrangein, "y" )
				axinset.set_xlabel("")
				axinset.set_ylabel("")
		elif args.iadd :
			if nax>1 :	nax +=1 
			else :		nax = 3
			if dd<1 : fig, ax = plt.subplots(nax,1, sharex=True )
			if nax > 1 :
				axarr = ax
			elif nax==1 : 
				axarr = [ax]
			print "nax : ", nax
			titleax = ax[0]
			totally(  fig, ax[0] , laboverw=leglabel) 
			totallyi( fig, ax[1] , laboverw=leglabel, factor=6 ) 
			selfeobj = datobj()
			selfeobj.getplotselfw( ax[2] , hobj.ddir , basis, nbasis, chdat, hobj , args )
			ax[0].axvline( 0 , color='gray' , lw = 0.5 , ls='--' ) 
			ax[1].axvline( 0 , color='gray' , lw = 0.5 , ls='--' ) 
			ax[2].axvline( 0 , color='gray' , lw = 0.5 , ls='--' ) 
		else :
			nax = ndir 
			if args.one : nax = 1
			naxrow	= nax
			naxcol	= 1
			if args.self	:	naxcol += 1 ; iaxselfim	= naxcol-1
			if args.hyb	:	naxcol += 1 ; iaxhybim	= naxcol-1
			if args.hybre	:	naxcol += 1 ; iaxhybre	= naxcol-1

			#if args.trans : 
			#	naxrow	= 1
			#	naxcol	= nax
			if dd<1 :
				print "naxrow, naxcol : ", naxrow, naxcol
				fig, ax = plt.subplots(naxrow, naxcol, sharex=True)
				print "ax shape : ", np.shape(ax)
				if nax > 1 : 
					if naxcol>1 : 
						axarr	= ax[:,0]
					elif naxcol==1 :
						axarr	= ax 
				if nax==1 : 
					if naxcol>1 : 
						axarr	= ax
					elif naxcol==1 :
						axarr	= [ax]

			if args.dashedline :
				if dd<1 : 
					dashpt = (None,None)
				else : 
					dashpt = [dd+1,1]
			else : dashpt = (None,None)

			cs	= ''
			lstyle	= args.linestylearr[dd] if args.linestylearr else None 
			if args.linestylecompoarr : 
				lstyle	= args.linestylecompoarr 
			if ndir > 1 : 
				if args.diffcolordir :  
					cs	= 'C%d'%dd
			if args.diffcolororb :  
				csarr	= [ ["C{}".format(ic)]*2 for ic in range(nbasis/2) ] * hobj.Ni
				csarr	= np.array( csarr ).flatten()
			elif args.diffcolororbcomp :  
				csarr	= [ ic for ic in args.diffcolororbcomp ] * hobj.Ni
			else :
				csarr	= False
			if csarr is not False : 
				print "csarr : ", csarr
			if nax > 1 : 
				totally( fig, axarr[dd]  , laboverw=leglabel, labeladdtail=leglabeltail, dashes=dashpt , hobj=hobj , cs=cs, lstyleArr=lstyle, csarr=csarr )
				titleax = axarr[0]
			else :
				totally( fig, ax , laboverw=leglabel, labeladdtail=leglabeltail, dashes=dashpt , hobj=hobj , cs=cs, lstyleArr=lstyle, csarr=csarr )
				titleax = ax 

			selfeobj = datobj()
			iaxmid  = ndir/2
			if args.self : 
				axself	= ax[dd,iaxselfim]
				selfeobj.getplotselfw( axself	, hobj.ddir , basis, nbasis, chdat, hobj , args , idir=dd , laboverw=leglabel , dashes=dashpt , logplot=logplot ) 
				ylabself	= axself.get_ylabel()
				axself.set_ylabel("")
				ax[iaxmid,iaxselfim].set_ylabel(ylabself)
			if args.hyb : 
				axhyb	= ax[dd,iaxhybim]
				selfeobj.getplothyb(   axhyb	, hobj.ddir , basis, nbasis, chdat, hobj , args ) 
				ylabhyb	= axhyb.get_ylabel()
				axhyb.set_ylabel("")
				ax[iaxmid,iaxhybim].set_ylabel(ylabhyb)
			if args.hybre : 
				axhybre	= ax[dd,iaxhybre]
				selfeobj.getplothybre( axhybre	, hobj.ddir , basis, nbasis, chdat, hobj , args ) 
				ylabhybre	= axhybre.get_ylabel()
				axhybre.set_ylabel("")
				ax[iaxmid,iaxhybre].set_ylabel(ylabhybre)
				
			if nax > 1 :
				legax   = ax[-2]
				axdos = axarr[iaxmid]
			else : 
				legax   = ax
				axdos	= ax
			legax = axarr[0]
			ylab = r'PDOS'
			#if len(chdat)<2 : 
			#	ylab = ylab + updnplotlabel(chdat[0],basis)[2]
			axdos.set_ylabel( ylab )

			if nax>1 and dd<ndir-1 :
				axarr[dd].tick_params( axis="x", labelbottom=False )#, labelsize = args.fsm ) 
			if nax>1 :
				for ia in range(iaxmid)+range(iaxmid+1,ndir) :
					axarr[ia].set_ylabel("")
			if args.vline : 
				gaparr = [-0.004469000000000278, -0.005885999999996727, -0.0003070000000029438, 0.01162000000000063, 0.017754999999997523, 0.03054599999999752, 0.042317000000000604, 0.05029700000000048]
				gaparr = [0, 0.052, 0.09 ]
				#gaparr = [0.019147000000000247, 0.03054599999999752, 0.042317000000000604, 0.05029700000000048] #J_H = [0.2,0.3,0.4,0.5]
				#redline = float(args.vline) # 0.048594999999998834
				gapval  = gaparr[dd]
				#gapval  = float(args.vline)
				if gapval > 0 : 
					axarr[dd].axvline( gapval, color='gray' , lw = 1 , ls='-.' ) 
					axarr[dd].axvline(-gapval, color='gray' , lw = 1 , ls='-.' ) 
				if gapval < 0 : 
					pass
					#axarr[dd].axvline( gapval, color='blue' , lw = 1 , ls='--' ) 
					#axarr[dd].axvline(-gapval, color='blue' , lw = 1 , ls='--' ) 
			if args.vline2 : 
				gapval  = float(args.vline2)
				axarr[dd].axvline( gapval, color='gray' , lw = 1 , ls='-.' ) 
				axarr[dd].axvline(-gapval, color='gray' , lw = 1 , ls='-.' ) 
		try : axnow = axarr[0] ; x0,x1,y0,y1 = axnow.axis()
		except : axnow = axarr ; x0,x1,y0,y1 = axnow.axis()
		#axnow.set_yticks([])
        	if args.sublabel :
			axnow = axarr[dd]
        	        dx=x1-x0;dy=y1-y0#; axnow.text( dx*0.02, dy*0.98 , markerssub[dd] ) 
			slab = markerssub[dd]
			if args.sublabelnoalphabet : 
				slab = ''
			if args.sublabelparameter : 
				for param in args.sublabelparameter : 
					if param=="nOcculattAraFormer" : 
		       				slab = slab +" {} = {:.2f}".format( filterlabel(param,paper=args.paper) , float(getattr(hobj,param)) )
					elif param.find("tfilling")>-1 : 
		       				slab = slab +" {} = {:.2f}".format( filterlabel(param,paper=args.paper) , float(getattr(hobj,param)) )
					elif param.find("Treal")>-1 : 
		       				slab = slab +" {} = {:.0f} K".format( filterlabel(param,paper=args.paper) , float(getattr(hobj,param)) )
					else : 
		       				slab = slab +" {} = {:.6g}".format( filterlabel(param,paper=args.paper) , float(getattr(hobj,param)) )
        	        axnow.text(0.02, 0.76, slab,   
	    				horizontalalignment='left', verticalalignment='center', transform=axnow.transAxes,
        	                        fontsize=fs+4 )# , weight='semibold')
	elif args.total :
		print "Plotting with args.total ."
		if dd<1 : fig, ax = plt.subplots()
		nax	= 1
		axarr	= [ax]
		titleax	= ax
		if args.dashedline :
			print "dashed-ON"
			if dd<1 : 
				dashpt = (None,None)
			else : 
				dashpt = [dd+1,1]
		if args.leglabel or args.leglabelarr :	
			basisCalc = basis
			#if args.transform :	
			#	if args.transform.find("t")>-1 :	basisCalc='t'
			#	elif args.transform.find("j")>-1 :	basisCalc='j'
			nchdat=len(chdat)
			legparArr	= args.leglabelarr if args.leglabelarr else [args.leglabel]
			labsuffix	= ''
			labsuffixSep	= '' if len(legparArr)==1 else ', '
			for legpar in legparArr :
				if legpar.find("Treal")>-1 : 
					if hobj.Treal<1e-5 : 
						labsuffix += r"$T = $"+"{:.0f} K ".format(hobj.Treal)
					else : 
						labsuffix += r"$T = $"+"{:.0f} K ".format(hobj.Treal)
					leglabel = [labsuffix]*nbasis
				elif legpar.find("nOcculattAra")>-1 or legpar.find("tfilling")>-1 :
					if args.paper : 
						labsuffix += filterlabel(legpar,paper=args.paper) + " = {:.2f}".format( float(getattr( hobj, legpar )) )
					else : 
						labsuffix += r"$n$"+" = {:.2f}".format( float(getattr(hobj,legpar)) )
					leglabel = [labsuffix]*nbasis
				elif args.paper : 
					labsuffix += filterlabel(legpar) + " = {}".format( getattr( hobj, legpar ) ) 
					leglabel = [labsuffix]*nbasis
				else : 
					#labsuffix = "${}={}$".format(legpar ,  getattr( hobj, legpar ) ) 
					labsuffix += filterlabel(legpar) + " = {}".format( getattr( hobj, legpar ) , paper=args.paper )
					if legpar.find("filling")>-1 : 
						labsuffix += filterlabel(legpar) + " = {:.2f}".format( float(getattr( hobj, legpar )) )
					leglabel = [labsuffix]*nbasis
				labsuffix += labsuffixSep
			if args.leglabelcomp : 
				leglabel = [ leglabel[j] +"; "+ bolabelArr(basisCalc)[j] for j in range(nbasis) ]
		else : 		leglabel = None
		#if args.leglabel :	
		#	if args.leglabel.find("nOcculattAra")>-1 : 
		#		leglabel = "$n$"+" = {:.2f}".format(float(getattr( hobj, args.leglabel )) ) 
		#	else : 
		#		labsuffix = filterlabel(args.leglabel,paper=args.paper) + " = {}".format( getattr( hobj, args.leglabel ) )
		#		if args.leglabel.find("filling")>-1 : 
		#			labsuffix = filterlabel(args.leglabel,paper=args.paper) + " = {:.2f}".format( float(getattr( hobj, args.leglabel )) )
		#		leglabel = labsuffix
		#	#else :
		#	#	leglabel = "${}={}$".format(args.leglabel ,  getattr( hobj, args.leglabel ) ) 
		#else : 		leglabel = None

		total(ax,dashes=dashpt,label=leglabel)
		#titleax = ax[0]
	elif args.selfonly :
		if args.one : 
			naxdef = 1 
			nax = 1 
		else :
			naxdef = ndir 
			nax = ndir 
		if args.hyb : nax += 1 ; iaxhybim = nax-1
		if args.hybre : nax += 1 ; iaxhybre = nax-1
		print "nax : ", nax
		if args.trans : 
			if dd<1 : fig, ax = plt.subplots(1,nax, sharey=True )
		else :
			if dd<1 : fig, ax = plt.subplots(nax,1 )
		if args.dashedline :
			if dd<1 : 
				dashpt = (None,None)
			else : 
				dashpt = [dd+1,1]
		else : dashpt = (None,None)
		if nax > 1 :
			axarr = ax
		else : 
			axarr = [ax]*ndir
		axself = axarr[dd]
		titleax = axself
		legax   = axself

		selfeobj = datobj()
		selfeobj.getplotselfw( axself , hobj.ddir , basis, nbasis, chdat, hobj , args , idir=dd , laboverw=leglabel , dashes=dashpt ) 
		if args.hyb : 
			selfeobj.getplothyb(   ax[iaxhybim]   , hobj.ddir , basis, nbasis, chdat, hobj , args ) 
		if args.hybre : 
			selfeobj.getplothybre( ax[iaxhybre] , hobj.ddir , basis, nbasis, chdat, hobj , args ) 
		if args.selfxrange : axSetRange( axself, args.selfxrange, "x" )
		if args.selfyrange : axSetRange( axself, args.selfyrange, "y" )
		if args.vline : 
			redline = float(args.vline) # 0.048594999999998834
			axself.axvline( redline , color='gray' , lw = 1 , ls='--' ) 
			axself.axvline(-redline , color='gray' , lw = 1 , ls='--' ) 
		if args.trans : 
			axself.set_ylabel(r'Im$\Sigma(\omega)$')#, fontsize=fs )
			axself.set_xlabel(r'$\omega$ (eV)')

		else : 
			axself.tick_params( axis="x", labelbottom=True )#, labelsize = args.fsm ) 
	else :
		print "Plotting with separate() ."
		ni = 1 
		try :
			ni = int( hobj.ni ) 
		except : pass
		if args.nimp :
	       		ni = int(args.nimp)	
		print "ni : ", ni

		if dd<1 :
			if ni<2 : 
				nax = 3
				fig, ax = plt.subplots(3,1, sharex=True, sharey=True )
				separate(dd,ax)
				titleax = ax[0]
				axarr = ax
			else : 
				nax = 3*ni
				fig, axmul = plt.subplots(3,ni, sharex=True, sharey=True )
				for ii in range(ni) :
					separate(dd, axmul.transpose()[ii], ini=ii )
				for ii in range(ni) :
					axnow = axmul[0][ii]
					nilab = "atom%d"%(ii)
					axnow.text(0.20, 0.85, nilab, transform=axnow.transAxes, bbox=dict(facecolor='white', edgecolor='none', pad=0), fontsize=fs*0.7 )
				titleax = axmul[0][0]
				ax = axmul[-1]
				axarr = axmul.flatten()
		else : 
			#print "Check the number of data-set :: ", ndir 
			if ni<2 : 
				separate(dd,ax)
				titleax = ax[0]
			else : 
				for ii in range(ni) :
					separate(dd, axmul.transpose()[ii], ini=ii )
				for ii in range(ni) :
					axnow = axmul[0][ii]
					nilab = "atom%d"%(ii)
					axnow.text(0.20, 0.85, nilab, transform=axnow.transAxes, bbox=dict(facecolor='white', edgecolor='none', pad=0), fontsize=fs*0.7 )
				titleax = axmul[0][0]
				ax = axmul[-1]
		for aax in axarr : 
			aax.tick_params( axis="both", direction="in" )#, labelsize=fs , length=10, width=2 ) 

#ax[-1].annotate('',
#            xy=(-0.017525, 0.4), xycoords='data',
#            xytext=(-15, 0), textcoords='offset points',
#            arrowprops=dict(arrowstyle="->") )
#ax[-1].annotate('',
#            xy=(0.017525, 0.4), xycoords='data',
#            xytext=(15, 0), textcoords='offset points',
#            arrowprops=dict(arrowstyle="->") )
#ax[-1].text( 0.1, 0.37, r'2|$\Delta_{sd}$|' ) 

try : 
	if nax > 1	: axarr = ax
	else		: axarr = [ax]
except :
	axarr = [ax]
try : 
	for axnow in axarr : 
		if args.ylabelcoord :
			ylabcoordx = float( args.ylabelcoord ) 
			axnow.yaxis.set_label_coords( ylabcoordx, 0.5)
		#axnow.tick_params( axis="y", labelleft=False )#, length=10, width=2 )
except : 
	pass
if args.ylabel : 
	axdos.set_ylabel( args.ylabel )

if args.minusx : 
	for axnow in axarr : 
		if axnow.get_xlabel()>0 : 
			axnow.set_xlabel( "-"+axnow.get_xlabel() ) 
for axnow in axarr : 
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
	#elif args.yoff :
	#	axnow.tick_params( axis="both", left=True,   labelleft=False  )
		
	
hobj=hobjarr[0]
if ndir > 1 and nax > 1  :
	print "shape of axarr : ", np.shape(axarr) 
	try : 
		for axnow in axarr[:-1,0] :
			axnow.set_xlabel("")
		for axnow in axarr[:-1,1] :
			axnow.set_xlabel("")
	except : 
		for axnow in axarr[:-1] :
			axnow.set_xlabel("")
if args.noylabel :
	plt.ylabel( "" ) 
if args.noxlabel :
	plt.xlabel( "" ) 
if args.chooselegax : 
	legax	 = ax[ int(args.chooselegax) ]
kwargs_leg  =   dict()
if args.legncol : 
	leg_ncol   = int(args.legncol)
	kwargs_leg.update( {"ncol":leg_ncol} )
	kwargs_leg.update( {"columnspacing":0.2} )
	#kwargs_leg.update( {"labelspacing":0.05} )
if args.legloc : 
	loc = int(args.legloc)
	kwargs_leg.update( {"loc":loc} )
if args.noleg : pass
elif args.alllegend :
	for axnow in axarr : 
		l = axnow.legend( **kwargs_leg ) ; l.set_zorder( 10 ) 
else : 
	try : 
		l = legax.legend(**kwargs_leg) ; l.set_zorder( 10 ) 
	except : 
		l = plt.legend(**kwargs_leg) ; l.set_zorder( 10 ) 
print "kwargs_leg: ", kwargs_leg

#axarr[1].set_xticks( [-1,-0.5,0,0.5,1] )
if args.xtick :
	try :
		xt =  np.array(args.xtick.split("_") , dtype=float)
	except : 
		xt =  np.array(args.xtick.split("_")[1:] , dtype=float)
	axarr[-1].set_xticks( xt )
if args.ytick :
	try :
		yt =  np.array(args.ytick.split("_") , dtype=float)
	except : 
		yt =  np.array(args.ytick.split("_")[1:] , dtype=float)
	axarr[-1].set_yticks( yt )
if args.fermiline : 
	for axnow in axarr : 
		axnow.axvline( 0 , color='k', linestyle="--", lw=1 ) 

try : 
	for axnow in axarr : 
		axnow.tick_params( axis='x',  pad=8)
except : 
	pass

ladjust=( AdjustSpaceDefault[0] if (args.ladjust is None) else float(args.ladjust) )
radjust=( AdjustSpaceDefault[1] if (args.radjust is None) else float(args.radjust) )
badjust=( AdjustSpaceDefault[2] if (args.badjust is None) else float(args.badjust) )
tadjust=( AdjustSpaceDefault[3] if (args.tadjust is None) else float(args.tadjust) )
wspace =( AdjustSpaceDefault[4] if (args.wspace  is None) else float(args.wspace ) )
hspace =( AdjustSpaceDefault[5] if (args.hspace  is None) else float(args.hspace ) )
print "lad, bad, rad, tad, wsp, hsp : " , ladjust, badjust, radjust, tadjust, wspace, hspace
plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, wspace=wspace , hspace=hspace )
if args.notitle :  pass
else : hobj.headSetTitle( titleax ) 

if args.grid :
	try : 
		n1, n2 = np.shape(ax)
		print "Shape of ax : ", n1, n2  
		if (n1*n2) > 1  :
			for axnow in ax.flatten() : 
				axnow.grid(linestyle=":")
		else :
			plt.grid(linestyle=":")
	except : 
		print "Shape of ax is ", nax
		if nax > 1 : 
			for axnow in axarr : 
				axnow.grid(linestyle=":")
		else : 
			plt.grid(linestyle=":")

print "Axis : ", plt.axis()

#axnow.axvline( 0 , color='k', linestyle=":" )
#axnow.axvline( 0.1 , color='k', linestyle=":" )
#axnow.set_facecolor('w')

#alph = [ "(b)" , "(d)", "(e)" ]
#plt.text(-0.18, 3,  alph[0], horizontalalignment='left', verticalalignment='center', transform=axnow.transAxes, fontsize=fs+4 )# , weight='semibold')
#
#for i in range(len(axarr)) :
#alph = [ "(a)" , "(b)", "(c)" ]  
#axarr[0].text(0.02, 1.96,  alph[0], horizontalalignment='left', verticalalignment='center', transform=axnow.transAxes, fontsize=fs+4 )# , weight='semibold')
#axarr[1].text(0.02, 0.88,  alph[1], horizontalalignment='left', verticalalignment='center', transform=axnow.transAxes, fontsize=fs+4 )# , weight='semibold')

argsavefigname = hobj.ddir + "/lattice/" + hobj.fignamepart() + "_"+fname
for ar in sys.argv[ndir+1:] :
	if ar.find("/")>-1 : 
		pass
	else : 
        	print "arg : ", ar
        	argsavefigname = argsavefigname + '_' + ar.split("-")[-1]
print "length(argsavefigname)" , len(argsavefigname)
if len(argsavefigname)>320 : argsavefigname = argsavefigname[:320]
fsave = argsavefigname
#print fsave , "\n", argsavefigname
#plt.savefig( fsave + ".png" ) 
#print "Saved : ", fsave + ".png"
ndpi=400
transparent = False 
if args.savetransparent : transparent = True
#if args.pdf :
#	plt.savefig( fsave + ".pdf" ) 
#	print "Saved : ", fsave + ".pdf"
if args.cptmp :
	dtype = "pdf"
	if args.pdf :
		dtype="pdf"
	if args.png :
		dtype="png"
	if args.eps :
		dtype="eps"
	print "Saved : ", fsave + "."+dtype
	plt.savefig( fsave + "."+dtype , dpi=ndpi, transparent=transparent, format=dtype ) 
	cptmpftn( fsave, cmd="scp", destMac=False , dtype=dtype )
	#plt.savefig( fsave + ".png" , dpi=ndpi ) 
	#print "Copying : ", fsave + ".png"
	#print "into : ",  "/home/jun/Dropbox/tmp_linux/" 
	#os.system( "cp " + fsave+".png " + "/home/jun/Dropbox/tmp_linux/" ) 
else : 
	print 'showing'
	plt.show()
