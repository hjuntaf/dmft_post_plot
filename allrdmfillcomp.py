#!/usr/bin/python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys, os
import subprocess

from basislabel import *
from hjl_common import *

import argparse

parser = argparse.ArgumentParser(description="Read directories from \"rdminfo.dat\" and Plot it." )  
parser.add_argument('dat', metavar='D', type=str, help='Specify a data-directory file' ) 
parser.add_argument('--objarr', "--obja", "--oarr", type=str, nargs='+', help='Specify name(s) of the data you want.' ) 
parser.add_argument("-u", "--use", type=str,    help="Column to be plotted. e.g. 1_2 (default) as x_y " ) 
parser.add_argument("-l", "--label", type=str,  help="Label of each data. e.g. 1's_2's_..." ) 
parser.add_argument("--ylabel", type=str,  help="Set y-label" ) 
parser.add_argument("-j", "--jbasis", action="store_true",  help="Plot the j-basis renorm. factor." ) 
parser.add_argument("-t", "--t2gbasis", action="store_true",  help="Plot the t2g-basis renorm. factor." ) 
parser.add_argument("--cptmp", action="store_true",  help="Copy .png to Dropbox." ) 
parser.add_argument("--pdf", action="store_true",  help="Save thd file in .pdf" ) 
parser.add_argument("--png", action="store_true", help="Save the figure in .png" ) 
parser.add_argument("--savetransparent", "--savetp", "--stp", action="store_true", help="Save the fig with transparent background" )
parser.add_argument("--figs", "--figsize", type=str, help="Set the size of figure." ) 
parser.add_argument("--fontsize", "--fonts", "--fs", type=float, help="Set fontsize." ) 
parser.add_argument("--ladjust", "--lad", type=str, help="Set   left-adjust value." ) 
parser.add_argument("--radjust", "--rad", type=str, help="Set  right-adjust value." ) 
parser.add_argument("--badjust", "--bad", type=str, help="Set bottom-adjust value." ) 
parser.add_argument("--tadjust", "--tad", type=str, help="Set    top-adjust value." ) 
parser.add_argument("--wspace",  "--wsp", type=str, help="Set wspace-adjust value." ) 
parser.add_argument("--hspace",  "--hsp", type=str, help="Set hspace-adjust value." ) 
parser.add_argument("-x", "--xxrange", "--xx", type=str,  help="Set xrange. e.g. _x0_x1" ) 
parser.add_argument("--tx", "--txxrange", "--txx", type=str,  help="Set xrange of jeff-states. e.g. _x0_x1" ) 
parser.add_argument("--xcol", type=str,  help="Choose the column of x-axis data." ) 
parser.add_argument("--notitle", action="store_true", help="Remove the title caption." )
parser.add_argument("--ni", type=int, help="Specify the number of impurity site (Ni)." ) 
parser.add_argument("--niter", type=int, help="Specify the number of DFMT-iterations." ) 
parser.add_argument("-d", "--degen", type=int, help="Specify the number of degeneracy." ) 
parser.add_argument("-y", "--yrange", type=str, help="Specify the range of y-value." ) 
parser.add_argument("--trans", action="store_true", help="Transverse the plot" )
parser.add_argument("-o", "--ord", action="store_true", help="Sorting the data in an ordinary order." )
parser.add_argument("--iimp", "--indimp", "--indimpurity", type=str, help="Specify the impurity you want to display and add the prefix s.t. 'imp0_'.")
parser.add_argument("--showrawdat", "--showdat", "--sdat", action='store_true', help="Show the raw data of the plot.")
args = parser.parse_args()


from matplotlib.patches import Ellipse
import matplotlib.pylab as pylab
figs=12
figsrat = 1
figsy = figs*figsrat
if args.trans :
	figsy=figs ; figs=figsy*figsrat
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


AdjustSpaceDefault = [ 0.10 , 0.93 , 0.10 , 0.95 , 0.1    , 0    ]

lab=[]
lsall=[]
maglab=[]
mls=[]
msind=[]
if args.label : 
	try : 
		lab.append( args.label.split("_")[j] ) 
	except : 
		lab.append( "" ) 
		lab.append( "" ) 
elif args.t2gbasis : 
	lab.append( "(xy,up)" ) 
	lab.append( "(xy,dn)" ) 
	lab.append( "(yz,up)" ) 
	lab.append( "(yz,dn)" ) 
	lab.append( "(zx,up)" ) 
	lab.append( "(zx,dn)" ) 
	lsall.append( "ro-" ) 
	lsall.append( "bo--" ) 
	lsall.append( "r>-" ) 
	lsall.append( "b>--" ) 
	lsall.append( "r^-" ) 
	lsall.append( "b^--" ) 
	datname = "fillingt2g"
	magind = [[0,1],[2,3],[4,5]]
	maglab.append( "$P_{xy}$" ) 
	maglab.append( "$P_{yz}$" ) 
	maglab.append( "$P_{zx}$" ) 
	mls.append( "ko-" ) 
	mls.append( "k>--" ) 
	mls.append( "k^:" ) 
	msind.append( "o" ) 
	msind.append( "o" ) 
	msind.append( ">" ) 
	msind.append( ">" ) 
	msind.append( "^" ) 
	msind.append( "^" ) 
	basisind = [ 0,0,1,1,2,2 ]
	colorind = [ 'r','b','r','b','r','b' ]
	labbasis    = r'$t_{{2g}}$'
else : 
	lab.append( "(3/2,3/2)" ) 
	lab.append( "(3/2,1/2)" ) 
	lab.append( "(3/2,-1/2)" ) 
	lab.append( "(3/2,-3/2)" ) 
	lab.append( "(1/2,1/2)" ) 
	lab.append( "(1/2,-1/2)" ) 
	lsall.append( "rx-" ) 
	lsall.append( "rD-" ) 
	lsall.append( "bD--" ) 
	lsall.append( "bx--" ) 
	lsall.append( "rs-" ) 
	lsall.append( "bs--" ) 
	datname = "fillingdeg"
	magind = [[0,3],[1,2],[4,5]]
	maglab.append( "$P_{(3/2,3/2)}$" ) 
	maglab.append( "$P_{(3/2,1/2)}$" ) 
	maglab.append( "$P_{(1/2,1/2)}$" ) 
	mls.append( "kx-" ) 
	mls.append( "kD--" ) 
	mls.append( "ks:" ) 
	msind.append( "x" ) 
	msind.append( "D" ) 
	msind.append( "D" ) 
	msind.append( "x" ) 
	msind.append( "s" ) 
	msind.append( "s" ) 
	basisind = [ 0,1 , 1,0, 2,2 ] 
	colorind = [ 'r','r','b','b','r','b' ]
	labbasis    = r'$j_{{eff}}$'

if args.ord : 
	ordCFjjJd4 = range(15) 
	ordCFLSJd4 = range(15)
else : 
	ordCFjjJd4 = [ 0, 3, 1, 2, 8, 6, 7, 4, 5, 9, 10, 13, 11, 12, 14 ]
	ordCFLSJd4 = [ 0, 3, 1, 2, 4, 5, 6, 7, 8, 13, 11, 12, 9, 10, 14 ]
CFordd4 = {	'jjJd4' : ordCFjjJd4 ,
		'LSJd4' : ordCFLSJd4 , 
		'LSJcubicd4' : ordCFLSJd4 } 

hobj = headobj( args.dat ) 
hobj.readParameters()
niter = hobj.niter - 1 
dimjjJd4  = 15
dimLSJd4  = 15
dimjeff   = 64
dimt2g    = 64
nargjjJ = 4


toolpath = "/home/jun/tools/dmft/plot/"

if args.degen :  degen = int( args.degen ) 
else : degen = int( hobj.degen_0 )
if args.ni : ni = int(args.ni)
else : ni = hobj.Ni 
if args.niter : niter = int(args.niter)
print "nimp  \t\t: ", ni
print "niter \t\t: ", niter
print "degen \t\t: ", degen

def plotrdm( ax , objname, kwargs=None ) :
	if objname.find("jjJd4")>-1	: ndat = dimjjJd4
	elif objname.find("LSJd4")>-1 or objname.find("LSJcubicd4")>-1	: ndat = dimLSJd4
	else 				: ndat = dimjeff
	width = 0.5 ; 
	datname		= "rdm"+objname
	datnamefull	= hobj.pathResult() + "/rdm"+objname +"diag.dat"
	iimp  = 0
	if args.iimp : 
		iimp  = int(args.iimp)
		prefixiimp  = "imp%d_"%int(args.iimp)
		datnamefull	= hobj.pathResult() + "/"+prefixiimp +"rdm"+objname +"diag.dat"
	print "Reading \t\t\t: ", datnamefull
	print "degen, ni : ", degen , ni 
	if degen < 2 and ni < 2 :
		if niter < 2 : 
			print "niter<2::"
			rdmdiag = np.genfromtxt( datnamefull, dtype=float )[1:]
		else : 
			print "niter>=2::"
			rdmdiag = np.genfromtxt( datnamefull, dtype=float )[-1][1:]
	elif degen<2 :
		rdmdiag	= np.genfromtxt( datnamefull, dtype=float )[1:]
	else :
		#rdmdiag = [""]*degen
		#fulldata = np.genfromtxt( datnamefull, dtype=float )
		#for j in range(degen) :
		#	rdmdiag[j] = fulldata[-j-1][1:]
		#rdmdiag[0] = np.genfromtxt( datnamefull, dtype=float )[-1][1:]
		#rdmdiag[1] = np.genfromtxt( datnamefull, dtype=float )[-2][1:]
		rdmdiag	= np.genfromtxt( datnamefull, dtype=float )[:,1:]
	print "rdmdiag.sum \t\t\t: ", np.sum(rdmdiag) 
	print "rdmdiag.max \t\t\t: ", np.max(rdmdiag) 
	print datname, "shape \t\t :({}) ".format( np.shape(rdmdiag) ) 

	boltzdatname	= hobj.pathResult() + "/boltzweight.dat"
	if ni < 2 :
		boltzweight	= np.genfromtxt( boltzdatname, dtype=float )[1:]
	else : 
		boltzweight	= np.genfromtxt( boltzdatname, dtype=float )[iimp,1:]
	print "boltzweight ({}) \t\t\t: ".format( np.shape(boltzweight) ) , boltzweight 
	if degen < 2 :
		rdmdiag = boltzweight * rdmdiag 
	else :
		rdmdiag = np.dot( boltzweight , rdmdiag  ) 
	print "rdmdiag.sum (boltzweighted) \t\t: ",  np.sum(rdmdiag) 
	print "rdmdiag ({}) (boltzweighted) \t: ".format( np.shape(rdmdiag) ) 
	#for a in rdmdiag :
	#	print a
	jz = [0]*degen
	#rdminfo = np.genfromtxt( toolpath + "/rdmarg{}.dat".format(datname) , dtype=str ) 
	#if degen > 1 :
	#	for a in range(degen) : 
	#		for i in range(ndat) : 
	#			jz[a] += rdmdiag[a][i]* rdminfo[i][-1]
	#			print "multiplying", rdminfo[i][-1]
	#		print "jz[", degen, "]", jz[a]
	ii = range(ndat)
	niorderargprev = np.genfromtxt( toolpath+"rdmniorderarg.dat", dtype=int ).transpose()[1] 
	if objname.find("jjJd4")>-1	: iiorder = ordCFjjJd4 ; print "jjJd4order",iiorder
	elif objname.find("LSJd4")>-1	: iiorder = ordCFLSJd4 ; print "LSJd4order",iiorder
	elif objname.find("LSJcubicd4")>-1	: iiorder = ordCFLSJd4 ; print "LSJd4order",iiorder
	elif objname == "" or (objname == "t2g") 		: iiorder = niorderargprev ; print "t2g",iiorder
	else 				: iiorder = range(ndat) ; print "integer-order",iiorder
	print "iiorder.shape : ", np.shape(iiorder) 
	if degen < 2 and ni < 2 :
		#rdmdiagorder    = np.zeros( ndat ) 
		#for a in ii : 
		#	rdmdiagorder[a] = rdmdiag[ iiorder[a] ]
		rdmdiagorder	= rdmdiag[ iiorder ]
	else : 
		#rdmdiagorder    = [ "" ] *degen 
		#for j in range(degen) : 
		#	rdmdiagorder[j] = np.zeros( ndat )
		#	for a in ii : 
		#		rdmdiagorder[j][a] = rdmdiag[j][ iiorder[a] ]
		#for a in ii : 
		#	rdmdiagorder[0][a] = rdmdiag[0][ iiorder[a] ]
		#	rdmdiagorder[1][a] = rdmdiag[1][ iiorder[a] ]
		rdmdiagorder	= rdmdiag[iiorder]
		rdmdiagordersum	= np.sum( rdmdiagorder, axis=0 )

	#bc  = "#151B8D" #  Denim Dark Blue
	#bc  = "#82CAFF" #  Day Sky Blue
	#bc  = "#C6DEFF" #  Heavenly Blue
	#bc  = "#FAF0DD" #  Soft Ivory
	#bc  = "#FFFFF0" #  Ivory (W3C)
	bc  = "#F8B88B" #  Pastel Orange
	if degen < 2 and ni < 2 :
		xdat	= ii
		ydat	= rdmdiagorder
		if args.showrawdat : 
			print "XDAT: "
			print xdat
			print "YDAT: "
			print ydat
		ax.barh( ii , ydat	 , width, align='center', alpha=0.9 , color=bc , edgecolor='gray' )
	elif len(rdmdiagorder.shape) < 2 :
		xdat	= ii
		ydat	= rdmdiagorder
		if args.showrawdat : 
			print "XDAT: "
			print xdat
			print "YDAT: "
			print ydat
		ax.barh( ii , ydat	 , width, align='center', alpha=0.9 , color=bc , edgecolor='gray' )
	else : 
		bc = "olive"
		if ni > 1 : bc="red"
		if degen > 3 : dii = 0.2
		else : dii = 0.15
		width = float(width)/degen
		iil = np.array(ii) - dii
		iir = np.array(ii) + dii
		iiarr		= np.linspace( -dii, dii, degen ) 
		alphaarr	= np.linspace( 0., 1., degen+1 ) 
		for ideg in range(degen ) :
			xdat	= np.array(ii)+iiarr[ideg]
			ydat	= rdmdiagorder[ideg]
			if args.showrawdat : 
				print "XDAT[degen=%d]: "%ideg
				print xdat
				print "YDAT[degen=%d]: "%ideg
				print ydat
			ax.barh( xdat , ydat , width, align='center',  alpha=alphaarr[ideg+1] , color=bc , edgecolor='gray' )
		#ax.barh( iil , rdmdiagorder[0] , width, align='center',  alpha=1.0 , color=bc )
		#ax.barh( iir , rdmdiagorder[1] , width, align='center',  alpha=0.5 , color=bc )
	ax.set_xlabel( r'$p_i$' )

	[x0,x1,y0,y1] = ax.axis()
	print "xylim \t\t\t\t: ", [x0,x1,y0,y1]
	#y1 = float( xylim[3] )
	#x1 = float( xylim[1] )
	if args.yrange :
		y = args.yrange.split("_")
		y0 = float(y[-2])
		y1 = float(y[-1])
		ax.set_ylim( y0, y1 )
	if objname.find("LSJ")>-1 or objname.find("jjJ")>-1 : 
		if args.xxrange :
			x = args.xxrange.split("_")
			x0 = float(x[-2])
			x1 = float(x[-1])
			ax.set_xlim( x0, x1 )
	elif args.tx : 
			x = args.tx.split("_")
			x0 = float(x[-2])
			x1 = float(x[-1])
			ax.set_xlim( x0, x1 )
	ax.set_yticks( np.arange(0,ndat,1) , minor=True ) 
	if args.ylabel : 
		if objname.find("LSJ")>-1 : 
			ax.set_ylabel( r'$i$-th of $|(S,L,J,J_z);$'+ labbasis + r'$^4>$' ) 
		elif objname.find("jjJ")>-1 : 
			ax.set_ylabel( r'$i$-th of $|(j_1,j_2,J,J_z);$'+ labbasis + r'$^4>$' ) 
		else : 
			ax.set_ylabel( r'$i$-th of $\{|$' + labbasis+ r'$>\}$' ) 
	ax.tick_params( axis='both', direction="in" , which='both',
			labelleft=True )
	print ""

def plotrdmarg( ax, datname )  :
	[x0,x1,y0,y1] = ax.axis()
	offconf=0.5*x1
	area = 6
	if datname == "" or (datname == "t2g") : 
	 	fnameinfo	= "/rdminfo.dat"
		rdminfo = np.genfromtxt( toolpath+fnameinfo, dtype=int ) 
		nc = 3 ; tnc  = 2*nc
		nb = 9 ; tnb  = 2*nb
		ndat = 64 
		ni = rdminfo.transpose()[-1]
		ii = range(ndat)
		niorderarg     = np.argsort( ni )
		niorderargprev = np.genfromtxt( toolpath+"rdmniorderarg.dat", dtype=int ).transpose()[1] 
		print "niorderarg test : "
		for j in range( len(niorderarg) ) :
			if( niorderarg[j] != niorderargprev[j]  )  : print "orderargdiff : ", j, niorderarg[j], niorderargprev[j] ; sys.exit(1)
		niorderdiv = np.genfromtxt( toolpath+"rdmniorderargdiv.dat" )
		for jj in range(ndat) :
			j = niorderarg[jj]
			for mu in range(nc) :
				st = rdminfo[j][mu]
				if ((st>>1)&01) : 
					#print 'u',
					tmu = 2*mu 
					ax.scatter( offconf*(0.9+float(tmu)*0.06), jj , s=area, marker=msind[tmu], c=colorind[tmu] , zorder=10 )
				else : 
					pass
					#print '0',
				if (st&01) : 
					#print 'd',
					tmu = 2*mu+1
					ax.scatter( offconf*(0.9+float(tmu)*0.06), jj , s=area, marker=msind[tmu], c=colorind[tmu] , zorder=10 )
				else : 
					pass
					#print '0',
				#print ' ',
			#print ''
		print "niorderdiv : ", niorderdiv
		for j in range(len(niorderdiv)) :
			ax.axhline(  niorderdiv[j],  linestyle=":" , linewidth=0.25 , color='black' ) 
			if j<len(niorderdiv)-1 :
				interval = (niorderdiv[j+1]+niorderdiv[j])/2. 
				#print "interval : ", interval 
				ax.text( offconf*(1-0.3), interval, "("+labbasis+")"+r'$^{:d}$'.format(j) , size=area*2., va="center", ha="center", rotation=0 ) 
	else : 
		fnameinfo	= "/rdmarg{}.dat".format(datname) 
		print "Reading fnameinfo \t\t: ", toolpath+fnameinfo
		rdminfo = np.genfromtxt( toolpath + fnameinfo , dtype=str ) 
		#print "niorder : ", niorder
		dimDum	= dimjjJd4
		if datname.find("d4")>-1 : 
			if datname.find("LSJ")>-1 : 
				basistext= r"($S,L;J,J_z$)" + ";("+labbasis+")"+r'$^{:d}$'.format(4) 
			elif datname.find("jjJ")>-1 : 
				basistext= r"($j_1,j_2;J,J_z$)" + ";("+labbasis+")"+r'$^{:d}$'.format(4) 
			else : 
				basistext= r"\{|j_{eff}>\}"  + ";("+labbasis+")"+r'$^{:d}$'.format(4) 
			stateTextFtn	= lambda tag	: "({:s},{:s};{:>2s},{:>2s})".format( *tag )
			dimDum	= dimjjJd4
			for jj in range(dimDum) :
				st = rdminfo[ CFordd4[datname][jj] ]
				#sttext = "({:s},{:s};{:>2s},{:>2s})".format( st[0], st[1], st[2], st[3] )
				sttext = stateTextFtn( st )
				ax.text( offconf*(1-0.1), jj, sttext , size=area*2.3, va="center", ha="center", rotation=0 ) 

			argdiv = np.genfromtxt( toolpath + "/rdmarg{}div.dat".format(datname) , dtype=int ) 
			for j in range(len(argdiv)) :
				ax.axhline(  argdiv[j]-0.5,  linestyle=":" , linewidth=0.25 , color='black' ) 
		else : 
			if datname == "SL" :
				basistext= r"($S,S_z,L,L_z$);$n_{el}$" 
			dimDum	= dimt2g 
			stateTextFtn	= lambda tag	: "({:s},{:s};{:>2s},{:>2s}|{:s})".format( *tag )

			#print "rdminfo.shape : ", np.shape(rdminfo)
			for jj in range(dimDum) :
				#print "jj : ", jj 
				st = rdminfo[jj]
				sttext = stateTextFtn( st )
				ax.text( offconf*(1-0.1), jj, sttext , size=area*1.4, va="center", ha="center", rotation=0 ) 

			argdiv = np.genfromtxt( toolpath + "/rdmarg{}div.dat".format(datname) , dtype=int ) 
			for j in range(len(argdiv)) :
				ax.axhline(  argdiv[j]-0.5,  linestyle=":" , linewidth=0.25 , color='black' ) 
				if j<len(argdiv)-1 :
					interval = (float(argdiv[j+1])+argdiv[j])/2. - 0.5
					#print "interval : ", interval 
					ax.text( offconf*(1-0.3), interval, "("+labbasis+")"+r'$^{:d}$'.format(j) , size=area*1.4, va="center", ha="center", rotation=0 ) 

		ax.text( offconf*(1-0.1), y1+0.5, basistext , size=area*1.4, va="center", ha="center", rotation=0 ) 
		#for aa in argdiv :
		#	ax.axhline(  aa-0.5,  linestyle=":" , linewidth=0.25 , color='black' ) 
	print ""
	
datarr  = args.objarr if args.objarr else ["", "jjJd4", "LSJd4"] #, "LSJcubicd4"]
print "objarr : ", datarr 
nax = len(datarr) 
fig, ax = plt.subplots(1,nax)
if nax < 2 : 
	axarr	= [ax]
	titleax	= ax
else : 
	axarr	= ax
	titleax	= ax[0]
	
for aax in axarr : 
	if args.xxrange : 
		axSetRange( aax , args.xxrange, 'x' ) 
	if args.yrange : 
		axSetRange( aax , args.yrange,  'y' ) 

for j in range( nax ) :
	objname = datarr[j]
	plotrdm(    axarr[j], objname )
	plotrdmarg( axarr[j], objname )

#titleax = axarr[0]
if args.notitle :  pass
else :
	hobj.readParameters()
	hobj.headSetTitle( titleax ) 

ladjust=( AdjustSpaceDefault[0] if (args.ladjust is None) else float(args.ladjust) )
radjust=( AdjustSpaceDefault[1] if (args.radjust is None) else float(args.radjust) )
badjust=( AdjustSpaceDefault[2] if (args.badjust is None) else float(args.badjust) )
tadjust=( AdjustSpaceDefault[3] if (args.tadjust is None) else float(args.tadjust) )
wspace =( AdjustSpaceDefault[4] if (args.wspace  is None) else float(args.wspace ) )
hspace =( AdjustSpaceDefault[5] if (args.hspace  is None) else float(args.hspace ) )
print "lad, bad, rad, tad, wsp, hsp : " , ladjust, badjust, radjust, tadjust, wspace, hspace
plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, wspace=wspace , hspace=hspace )


fname = "rdmfillcomp" 
argsavefigname = hobj.tdir + "/"+hobj.fignamepart() + "_"+fname
for ar in sys.argv[2:] :
        print "arg : ", ar
        argsavefigname = argsavefigname + '_' + ar.split("-")[-1]
fsave = argsavefigname
if args.cptmp :
	dtype = "pdf"
	if args.pdf :
		dtype="pdf"
	if args.png :
		dtype="png"
	print "Saved : ", fsave + "."+dtype
	ndpi=400
	transparent = False 
	if args.savetransparent : transparent = True
	plt.savefig( fsave + "."+dtype , dpi=ndpi, transparent=transparent, format=dtype ) 
	cptmpftn( fsave, cmd="scp", destMac=False , dtype=dtype )
else : 
	plt.show()
