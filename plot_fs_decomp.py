#!/usr/bin/python
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
parser.add_argument("--savepdf", action="store_true", help="Save the fig in .pdf too" )
parser.add_argument("--savetransparent", "--savetp", "--stp", action="store_true", help="Save the fig with transparent background" )
parser.add_argument("--fs", action="store_true" , help="Plot the fermi-surface." ) 
parser.add_argument("--oldsoc", action="store_true" , help="SOC strength is in a old-format and multiplied by 2./3." )
parser.add_argument("-l", "--label", type=str, help="set text-label" ) 
parser.add_argument("-y", "--yrange", type=str, help="set yrange e.g. _y0_y1" )
parser.add_argument("-x", "--xxrange", type=str, help="set xrange e.g. _x0_x1" )
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
parser.add_argument("--cb", action="store_true", help="Add the colorbar" ) 
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
parser.add_argument("--decomp", action="store_true", help="Plot the basis-decomposed one." )
parser.add_argument("--decompall", action="store_true", help="Plot the all-basis-decomposed one." )
parser.add_argument("--notitle", action="store_true", help="Take no title" ) 
parser.add_argument("--levhrat", type=float, help="Set the ratio of plotmax- and spectmax-levels." )
parser.add_argument("--cptmp", action="store_true", help="Copy the png file to ~/Dropbox/tmp_linux/" ) 
parser.add_argument("--showaxes", action="store_true", help="Show axes" ) 
parser.add_argument("--epsilon", "--ep", type=str, help="Set Eipsilon value for different data set" )
parser.add_argument("--epsilontag", "--ept", "--eptag", action="store_true",  help="Set Eipsilon tag in file-name." )
parser.add_argument("--epsilontagfloat", "--eptf", "--eptagf", action="store_true",  help="Set Eipsilon tag with float number in file-name." )
parser.add_argument("--broadening", "--br", type=str, help="Set Broadening value for different data set" )
parser.add_argument("--broadeningtag", "--brt", "--brtag", action="store_true", help="Set Broadening tag in file-name." )
parser.add_argument("--notick", action="store_true", help="Remove the ticks" ) 
parser.add_argument("--nof", "--nofermiline", action="store_true", help="Remove the line of Fermi level." ) 
parser.add_argument("--cmjet",  action="store_true", help="Change the color-scheme into plt.cm.jet." ) 
parser.add_argument("--cm",  type=str, help="Specify the color-scheme, e.g. jet of plt.cm.jet." ) 
parser.add_argument("--proj",  "--project", action="store_true", help="Plot the basis-projected spectrum." ) 
parser.add_argument("--fswn", type=float, help="Set w for FS" )
parser.add_argument("--fsind", type=str, help="Plot a different FS along GX(GV), or with kz-off. e.g. --fsind {vX,vM,or kz0.5}" )
parser.add_argument("--tl",  action="store_true", help="Tight-layout" ) 
parser.add_argument("--count", type=int, help="Set iter")
parser.add_argument("--figs", "--figsize", type=str, help="Set the labels of the legend" ) 
parser.add_argument("--fontsize", "--fonts", type=str, help="Set the size of fonts.")
parser.add_argument("--ladjust", "--lad", type=str, help="Set   left-adjust value." ) 
parser.add_argument("--radjust", "--rad", type=str, help="Set  right-adjust value." ) 
parser.add_argument("--badjust", "--bad", type=str, help="Set bottom-adjust value." ) 
parser.add_argument("--tadjust", "--tad", type=str, help="Set    top-adjust value." ) 
parser.add_argument("--noleg", "--nol",  action="store_true", help="Remove the legend.") 
parser.add_argument("--leglabel", "--legl", type=str, help="Set the labels of the legend" ) 
parser.add_argument("--noylab", "--noylabel",  action="store_true", help="Remove the y-label.") 
parser.add_argument("--nolab", "--nolabel",  action="store_true", help="Remove the xy-label.") 
parser.add_argument("--locmaxcolor", "--lmc", "--lmcolor", type=str, help="Set a color of the locmax-plot." ) 
parser.add_argument("--fermivelocity", "--fermiv", action="store_true",  help="Obatina Fermi velocities.")
parser.add_argument("--ldos", "--ld", "--pd", "--pdos", action="store_true",  help="Add PDOS as a subplots")
parser.add_argument("--nplot", type=int, help="Set Nplot." ) 
parser.add_argument("--naxis", type=int, help="Set Naxis." ) 
parser.add_argument("--angle", "--ang", type=str, help="Set kM-kX angle." ) 
parser.add_argument("--trans", "--tr", action="store_true", help="Transpose the plot" ) 
parser.add_argument("--yoff", action="store_true", help="Remove xticks." ) 
parser.add_argument("--dxtick", action="store_true", help="Double the xtick increment." ) 
parser.add_argument("--dytick", action="store_true", help="Double the ytick increment." ) 
parser.add_argument("--fillb", action="store_true", help="Fill the plotting of blue one." ) 
parser.add_argument("--sublabel", "--slab", action="store_true", help="Set the labels of sub-figures with alphabet." ) 
parser.add_argument("--sublabelparameter", "--slabpar", "--slabp", type=str, help="Set a parameter of the sub-labels.")
parser.add_argument("--sublabeloff", "--slaboff", "--slabo", type=str, help="Set the off-set of the sub-labels.")
parser.add_argument("--fsonly", "--fso", action="store_true", help="Plot fs only.")
parser.add_argument("--sympointwhite", "--symw", "--symptw", action="store_true", help="Plot points of high-symmetry as white.")
args = parser.parse_args()

import matplotlib.pylab as pylab
figs=8.  
figsrat = 0.8
figsratpar = 0.25
if args.decomp      : figsrat = figsratpar
elif args.decompall : figsrat = figsratpar*2
figsy = figs*figsrat
fsize=10
if args.fontsize : fsize = int( args.fontsize )
fs=fsize
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
params = {      'legend.fontsize'       : fs ,
		'figure.figsize'        : (figs,figsy) ,
		'axes.labelsize'        : 'x-large' ,
		'axes.titlesize'        :'small' ,
		'xtick.labelsize'       :'x-large' ,
		'ytick.labelsize'       :'x-large' ,
		'lines.markersize'      : ms ,
		#'markers.fillstyle'     : 'none' ,
		'font.size'             : fs }
pylab.rcParams.update(params)

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

argsavefigname = dat + "/lattice/vdx/" + hobj.fignamepart()
for ar in sys.argv[ndir+1:] :
	print "arg : ", ar 
	argsavefigname = argsavefigname + ar
levhrat = 1
if args.levhrat :
	levhrat = args.levhrat

if args.notitle : hobj.title = ""
if args.basis : 
	basis = args.basis
	if basis.find("t")>-1 : 
		basis = "t2g"
		chdat = [0,2,4]
		chdat2= [1,0,4]
	if basis.find("j")>-1 : 
		basis = "j"
		chdat = [1,0,4]
		chdat2= [0,2,4]
	if basis.find("s")>-1 : 
		basis = "s"
		chdat = [0,2,3]
		chdat2= [0,1,4]
else : 
	basis = "t2g"
	chdat = [0,2,4]
	chdat2= [1,0,4]
if args.chdat : 
	chdat = args.chdat.split("_")
elif args.dns :
	chdat = [1,3,5]

epsilon = 0.02
if args.epsilon :
	epsilon = float(args.epsilon)
	for hobj in hobjarr : hobj.title = hobj.title + " ep={}".format(epsilon) 
for hobj in hobjarr : setattr( hobj , "epsilon", float(epsilon) ) 
print "epsilon : ", epsilon
nplot = 1024
if args.nplot : nplot = int(args.nplot)
naxis = 4
if args.naxis : naxis = int(args.naxis)
angle = False ;
if args.angle : angle = int(args.angle)
def fanglelab(angle) : 
	if (angle is False) :
		return ""
	else :
		return "_angle{}".format(angle)
anglelab = fanglelab(angle)
angletick = '${}^\circ$'.format(angle)
if args.fswn : 
	fswn = float( args.fswn ) 
else :
	fswn = 0  
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
fileind = "fs"+fsind+"w{:.2f}".format(fswn)
if args.broadeningtag :
	fileind = "fs"+fsind+"w{:g}".format(fswn)


broadening = 0.02
if args.broadening :
	broadening = float(args.broadening)
	for hobj in hobjarr : hobj.title = hobj.title + " $\eta_d$={}".format(broadening) 
epsilontag	= ""
broadeningtag	= ""
if args.epsilontag :
	epsilontag = "ep{:g}_".format(epsilon)
if args.epsilontagfloat :
	epsilontag = "ep{:.2f}_".format(epsilon)
if args.broadeningtag :
	broadeningtag = "br{:g}_".format(broadening)
tagEpBr =  epsilontag + broadeningtag

#figs = 12
#fs = 50 ; fsm = 50
cmapcolor = plt.cm.nipy_spectral
if args.cmjet : 
	cmapcolor = plt.cm.jet
if args.cm : 
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
	def __init__ ( self, fig, cmapcolor , ax=False )  :
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
		if "band" in fspect :
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
			print "\rnline={:10d}".format(ncount),
			#print "\rnline={:3.1f}%".format(float(ncount)/ntot*100.),
		
		print "nc=", ncount
		print "ndata=", ndata
		
		if "u{:.2f}".format(hobj.u) in fspect :
			X, Y = np.meshgrid(xi, yi)
			Z = np.array(Z)
			Z=Z.reshape(X.shape)
	
			#Z = np.ma.array( Z, mask=Z < 1 ) 
			if ax : self.cbar = ax.contourf(X, Y, Z, 300, rstride=1, cstride=1, cmap=cmapcolor, linewidth=0, antialiased=False)
			if args.cb :
				if args.fs : pass
				else : self.cb = plt.colorbar( self.cbar )
			#plt.axis([0,1,-3,6])
			if args.yrange :
				ymin = float( args.yrange.split("_")[-2] ) 
				ymax = float( args.yrange.split("_")[-1] ) 
				if ax : ax.axis( [0,1,ymin,ymax] )
			if ax : axSetRange( ax, args.xxrange , "x" ) 
			if args.label :
				if ax : ax.set_text( 0.8, ymax*1.02, "{}".format( sys.argv[i+1] ) )
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
		if "band" in fspect :
			X = xi
			Z = np.array(Z)
			Z = Z.reshape((nx,nbasis))
			Z = Z.transpose() 
			print "len x,z  =", len(X), len(Z)
			for i in range( nbasis ) :
				if ax : ax.plot(X, Z[i], "b-")
			if ax : ax.axis([0,1,-3,2])
			if args.yrange :
				ymin = float( args.yrange.split("_")[-2] ) 
				ymax = float( args.yrange.split("_")[-1] ) 
				if ax : ax.axis( [0,1,ymin,ymax] )
			if args.label :
				if ax : ax.text( 0.8, ymax*1.02, "{}".format( args.label ) )
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
		
		print "xtics rescaled : ",
		nkgridtot = 0
		xtics = np.array([0., 0., 0., 0.])
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
		
		if ax : ax.set_xticks(xtics)
		if ax : ax.set_xticklabels(knames ) 
		if ax : ax.tick_params( axis='y')
		if ax : ax.tick_params( axis='both', direction="in"  )
		if ax : ax.set_ylabel(r"$E$ (eV)") 
		if ax : axSetRange( ax, args.xxrange , "x" ) 

		if args.noylab : ax.set_ylabel('') ; ax.set_yticks([])
		
	
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

class klinedecomp :
	def __init__( self, ax , bcomp , lev , cmapcolor ) :
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

		Z= np.genfromtxt( fspect, dtype=float ) 
		ntot = self.nx*self.ny
		nbasis = 6
		ncount = len(Z) * len(Z[0])
		ndata  = len(Z[0])
		print "nx,ny,ntot,ncount,nZline,ndata : ", self.nx, self.ny, ntot, ncount, len(Z), ndata,
		
		f2= open(fspect, "r")
		nbasis = 6
		if "u{:.2f}".format(hobj.u) in fspect :
			print "[", bcomp, "component]"
			#self.Z = [0.]*ntot
		ncount = 0
		ndata  = 0

		if basis.find("j")>-1 :
			self.Z = np.array( Z[bcomp] ,dtype=float )
		else : 
			self.Z = np.array( Z[bcomp+6] ,dtype=float )
		
		ax.axis( [ xi[0],xi[-1], yi[0],yi[-1] ] ) 
		if "u{:.2f}".format(hobj.u) in fspect :
			X, Y = np.meshgrid(xi, yi)
			self.X = X
			self.Y = Y
			self.Z = np.array(self.Z)
			self.Z = self.Z.reshape(X.shape)
			self.cbar = ax.contourf(X, Y, self.Z, lev, rstride=1, cstride=1, cmap=cmapcolor, linewidth=0, antialiased=False)
			if args.cb : 
				cb = plt.colorbar( self.cbar )

			if args.yrange :
				ymin = float( args.yrange.split("_")[-2] ) 
				ymax = float( args.yrange.split("_")[-1] ) 
				ax.axis( [0,1,ymin,ymax] )
			if args.label :
				ax.text( 0.8, ymax*1.02, "{}".format( args.label ) )
			#ax.text( 0.8, ymax*1.02, "{}".format( bolabel(basis,bcomp) ) )
			if basis.find("t")>-1 : fsrat = 1
			else : fsrat=1
			ax.text( 0.5, ymax*0.8, r"${}$".format( bolabel(basis,bcomp) ) , horizontalalignment='center', va='center', color='w' , size=fsize*fsrat )  

			dashpt = [1,1]
			if args.nof : pass
			else : ax.axhline( 0, 0,1, color='w', dashes=dashpt, lw=0.5 )

			print "xtics rescaled : ",
			nkgridtot = 0
			xtics = np.array([0., 0., 0., 0.])
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
			if args.cb and (args.decompall is False) :
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
				if args.sympointwhite : ax.scatter( 0.5,  0.5, 'o',  color='w' , zorder=12 )
#				if ax : ax.text( 0.5,  0.5, kc		, color=fc , fontsize=fsize , zorder=12 )
#				if ax : ax.text( 0.5,  0.5, kc		, color=fc , fontsize=fsize , zorder=12 )
#				if ax : ax.text( 0.88, 0.5, kr		, color=fc , fontsize=fsize , zorder=12 )
#				if args.fs : 
#					if ax : ax.text( 0.88, 0.88, knames[0]	, color=fc , fontsize=fsize , zorder=12 )
#				else : 
#					if ax : ax.text( 1.0, 0.0,   knames[0]	, color=fc , fontsize=fsize , zorder=12 )
			except : pass
			if args.sympointwhite : ax.scatter( 0.5,  0.5, marker='o',  color='w' , zorder=12 , s=5)
			if args.sympointwhite : ax.scatter( 0. ,  1. , marker='o',  color='w' , zorder=12 , s=5)
			if args.sympointwhite : ax.scatter( 0. ,  0.5, marker='o',  color='w' , zorder=12 , s=5)
#if ax : ax.text( 0.5,  0.5, kc		, color=fc , fontsize=fsize , zorder=12 )
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
				if ax : ax.text( 0.8, 1.02, "{}".format( args.label ) )
			if fswn != 0 :
				if ax : ax.text( 0.1, 1.02, "w={}".format( fswn ) )
		#Ztr   = self.Z.transpose()
		#dztr = Ztr - self.Z
		#print "NNZ Z-Z.tr :"
		#for ii in range(self.nx) :
		#	for jj in range(self.nx) :
		#		if dztr[ii][jj] > 1e-6 :
		#			print "({}\t{})".format(ii,jj) 
		findLocalMaxSave( self.Z, "{}/lattice/vdx/".format( hobj.ddir ) , hobj  )
	
class fsurfacedecomp :
	def __init__( self, ax , bcomp , lev ) :
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
		if "u{:.2f}".format(hobj.u) in fspect :
			print "Read", ntot, "density data in", bcomp, "component"
			self.Z = [0.]*ntot
		ncount = 0
		ndata  = 0
		Z = np.genfromtxt( fspect , dtype=float ) 
		if basis.find("j")>-1 :
			self.Z = np.array( Z[bcomp] ,dtype=float )
		elif basis.find("s")>-1 :
			self.Z = np.array( Z[bcomp+12] ,dtype=float )
		else : 
			self.Z = np.array( Z[bcomp+6] ,dtype=float )
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
			self.Z=self.Z.reshape(X.shape)
			self.cbar = ax.contourf(X, Y, self.Z, lev, rstride=1, cstride=1, cmap=cmapcolor, linewidth=0, antialiased=False)
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
#				ax.text( 0.5, 0.5, kc		, fontsize=fs , color=fc )
#				ax.text( 1.0, 0.5, kr		, fontsize=fs , color=fc )
#				ax.text( 1.0, 1.0, knames[0]	, fontsize=fs , color=fc )
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
				ax.text( 1, ymax*1.02, "{}".format( args.label ) )
			if basis.find("t")>-1 : fsrat = 1
			else : fsrat=1
#			ax.text( 0.5, ymax*0.5, r"${}$".format( bolabel(basis,bcomp) ) ,
#					horizontalalignment='center', va='center' , color='w', size=fsize*fsrat ) 
		#Ztr   = self.Z.transpose()
		#dztr = Ztr - self.Z
		#print "NNZ Z-Z.tr :"
		#for ii in range(self.nx) :
			#for jj in range(self.nx) :
				#if dztr[ii][jj] > 1e-6 :
					#print "({}\t{})".format(ii,jj) 
		findLocalMaxSave( self.Z, "{}/lattice/vdx/".format( hobj.ddir ) , hobj )
	
	
class bandonly : 
	def __init__ ( self, fig, cmapcolor , ax=False )  :
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
			if ax : ax.text( 0.8, ymax*1.02, "{}".format( args.label ) )

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
	if args.fs :
		if args.th : 
			zthratio = float(args.th) 
		else : 
			zthratio = 0.3
	
		if fswn == 0 :
			try : 
				fspect = "{}_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( fileind, hobj.UF, hobj.D, tagEpBr, hobj.count ) 
				fspectdecomp = "decomp_{}_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( fileind, hobj.UF, hobj.D, tagEpBr, hobj.count ) 
				faxis  = "{}_naxis{}_Nplot1_Nk40401_{}{:d}th.dx".format( fileind, naxis, tagEpBr, hobj.count ) 
				faxisdecomp  = "decomp_{}_naxis{}_Nplot1_Nk40401_{}{:d}th.dx".format( fileind, naxis, tagEpBr, hobj.count ) 
			except : 
				fspect = "fs_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( hobj.UF, hobj.D, tagEpBr, hobj.count ) 
				fspectdecomp = "decomp_fs_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( hobj.UF, hobj.D, tagEpBr, hobj.count ) 
				faxis  = "fs_naxis{}_Nplot1_Nk40401_{}{:d}th.dx".format( naxis, tagEpBr, hobj.count ) 
				faxisdecomp  = "decomp_fs_naxis{}_Nplot1_Nk40401_{}{:d}th.dx".format( naxis, tagEpBr, hobj.count ) 
		else : 
			fspect = "{}_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( fileind, hobj.UF, hobj.D, tagEpBr, hobj.count ) 
			fspectdecomp = "decomp_{}_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( fileind, hobj.UF, hobj.D, tagEpBr, hobj.count ) 
			faxis  = "{}_naxis{}_Nplot1_Nk40401_{}{:d}th.dx".format( fileind, naxis, tagEpBr, hobj.count ) 
			faxisdecomp  = "decomp_{}_naxis{}_Nplot1_Nk40401_{}{:d}th.dx".format( fileind, naxis, tagEpBr, hobj.count ) 
		fspect = cwd + fspect
		fspectdecomp = cwd + fspectdecomp
		faxis = cwd + faxis
		faxisdecomp = cwd + faxisdecomp
	
		dispfspect = "u{:.2f}_ts0.00_d{:.3f}_{}{:d}th{}.txt".format( hobj.UF, hobj.D, tagEpBr, hobj.count, anglelab ) 
		dispfspectdecomp = "decomp_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th{}.txt".format( hobj.UF, hobj.D, tagEpBr, hobj.count, anglelab ) 
		dispfaxis  = "naxis{}_Nplot{}_Nk341_{}{:d}th{}.dx".format( naxis, nplot, tagEpBr, hobj.count, anglelab) 
		dispfaxisdecomp  = "decomp_naxis{}_Nplot{}_Nk341_{}{:d}th{}.dx".format( naxis, nplot, tagEpBr, hobj.count, anglelab ) 
		dispfspect = cwd + dispfspect
		dispfspectdecomp = cwd + dispfspectdecomp
		dispfaxis = cwd + dispfaxis
		dispfaxisdecomp = cwd + dispfaxisdecomp
	
		if os.system("test -e {}".format(fspect)) < 1 :  pass
		else :
			print "Error :: ", fspect, "doesn't exist." 
			sys.exit(1)
	
		if args.locmax or args.locmaxonly : 
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
			axSetRange( ax[0], args.lineplotxrange , "x" ) 
			axSetRange( ax[0], args.lineplotyrange , "y" ) 

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
			if len(chdat)==3 : 
				if dd<1 : fig, ax = plt.subplots( 1,4 )
				fs = [""]*4
				fs[0] = fsurface( ax[0] ) ;
				levtot = fs[0].cbar.levels
				i=1
				for ch in chdat :
					"dat : ", basis, ch
					fs[i] = fsurfacedecomp( ax[i] , int(ch) , levtot )
					i += 1
				titleax = ax[1]
				if args.cb : 
					fig.colorbar( fs[0].cbar , ax=ax.ravel().tolist())
				#print "CBlevel : ", fs[0].cbar.levels
				print "nCBlevel : ", len(fs[0].cbar.levels)
			else :
				if dd<1 : fig, ax = plt.subplots( 1,len(chdat)+1 )
				i = 0 
				fs = fsurface( ax[i] ) ; i+=1
				levtot = fs.cbar.levels
				for ch in chdat :
					"dat : ", basis, ch
					fs = fsurfacedecomp( ax[i] , int(ch) , levtot )
					i += 1
				titleax = ax[len(chdat)/2]
			if args.showaxes :  pass
			else :
				for aax in ax.flat :
					aax.axis("off") 
		elif args.decompall: 
			if dd<1 : fig, axarr = plt.subplots( 2,3 )
			ax = axarr.flatten()
			fs = [""]*7
			fs[0] = fsurface( ax[0] ) ;
			levtot = fs[0].cbar.levels
			i=0
			#if fswn != 0 : 
				#ax[1][0].remove()
			#else : 
			fspectdum = fspect
			faxisdum  = faxis
#			fspect = dispfspect
#			faxis  = dispfaxis
#			if args.fsonly : pass
#			else : kline( fig , cmapcolor , ax=ax[1][0]) 
#			ax[1][0].axhline( fswn, 0,1 , color='w', dashes=[1,2] , lw=1) 
#			ax[1][0].set_ylim( -1,1)
			fspect = fspectdum
			faxis  = faxisdum
			for ch in chdat :
				"dat : ", basis, ch
				fs[i] = fsurfacedecomp( ax[i] , int(ch) , levtot )
				i += 1
			if basis.find("t")>-1 :   basis2 = "j"
			elif basis.find("j")>-1 : basis2 = "t"
			elif basis.find("s")>-1 : basis2 = "j"
			i=0
			chdat = chdat2
			for ch in chdat :
				basis1= basis
				basis = basis2
				"dat : ", basis, ch
				fs[i] = fsurfacedecomp( ax[i+3] , int(ch) , levtot )
				i += 1
			titleax = ax[0]
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
		fspect = "u{:.2f}_ts0.00_d{:.3f}_{}{:d}th{}.txt".format( hobj.UF, hobj.D, tagEpBr, count, anglelab ) 
		fspectdecomp = "decomp_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th{}.txt".format( hobj.UF, hobj.D, tagEpBr, count, anglelab ) 
		faxis  = "naxis{}_Nplot{}_Nk341_{}{:d}th{}.dx".format( naxis, nplot, tagEpBr, count, anglelab) 
		faxisdecomp  = "decomp_naxis{}_Nplot{}_Nk341_{}{:d}th{}.dx".format( naxis, nplot, tagEpBr, count, anglelab) 
		fspect = cwd + fspect
		fspectdecomp = cwd + fspectdecomp
		faxis = cwd + faxis
		faxisdecomp = cwd + faxisdecomp
		if args.decomp: 
			#figs = 12
			#fs = 25 ; fsm = 25
			if dd<1 : fig, ax = plt.subplots( 2,2 )
			kclass = [""]*4
	
			kclass[0] = kline( fig , cmapcolor , ax=ax[0][0]) ;
			levtot = kclass[0].cbar.levels
	
			i=1
			if len(chdat)==3 :
				for ch in chdat :
					"dat : ", basis, ch
					kclass[i] = klinedecomp( ax[i/2][i%2] , int(ch) , levtot , cmapcolor )
					i += 1
			titleax = ax[0][0]
			if args.cb : 
				fig.colorbar( kclass[0].cbar , ax=ax.ravel().tolist())
			print "nCBlevel : ", len(kclass[0].cbar.levels)
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
				fspect = cwd + "{}_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( fileind, hobj.UF, hobj.D, tagEpBr, count ) 
				fspectdecomp = cwd + "decomp_{}_u{:.2f}_ts0.00_d{:.3f}_{}{:d}th.txt".format( fileind, hobj.UF, hobj.D, tagEpBr, count ) 
				faxis  = cwd + "{}_naxis{}_Nplot1_Nk40401_{}{:d}th{}.dx".format( naxis, fileind, tagEpBr, count, anglelab ) 
				faxisdecomp  = cwd + "decomp_{}_naxis{}_Nplot1_Nk40401_{}{:d}th{}.dx".format( naxis, fileind, tagEpBr, count, anglelab ) 
				fsurface( ax[1][0] ) ;
			except : 
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
			bandfspect = "band_naxis{}_Nk341{}.dat".format(naxis,anglelab)
			bandfaxis  = "naxis{}_Nplot1024_Nk341_{}{:d}th{}.dx".format( naxis, tagEpBr, hobj.count, anglelab) 
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
		else :
			#fig = plt.figure(figsize=(figs,figs))
			#ax = plt.subplot2grid((5, 5), (0, 1), colspan=4, rowspan=4 )
			if dd<1 : fig, axone = plt.subplots( ) ; ax = [axone]*5
			if args.lineplotk :
				kl = kline( fig , cmapcolor )
				linearg = "y_0.5"
				if args.lineplotkw : linearg = args.lineplotkw 
				lineaxis = linearg.split("_")[0]
				n0, xdat, ydat, datsection = refineLinePlot( linearg , kl )
				if args.leglabel :	leglabel = "${}={}$".format(args.leglabel ,  getattr( hobj, args.leglabel ) ) 
				else : 			leglabel = linearg
				axone.plot( xdat[n0], ydat[n0]*nbasis , "o-", lw = 1, label=leglabel , markerfacecolor='w')#, markeredgecolor='k')
				saveLinePlot( hobj.ddir , linearg , n0, xdat[n0], ydat[n0] , hobj )
				#plotIndLine( ax[0], datsection, lineaxis, (2,1) , 'w' )  if args.lineplotonly is None else '.'
				if   lineaxis=="x" : axone.set_xlabel("$A(\omega={:.3f})$".format(datsection)) #; axone.tick_params( axis="y", which='both', labelleft=False  )#, labelsize = args.fsm ) 
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
						
				else : 
					axone.set_xticks(kl.xtics)
					axone.set_xticklabels(kl.knames) 
				print "YDATpart(frequency) ", n0, "/",np.shape(kl.Y) ," : ", [ kl.Y[a][0] for a in range(n0-3,n0+3) ] , "(dn={}".format(1./1024),")"
			else : 
				if args.ldos : 
					args.trans=True
					ldos = ldos(args, fname="tldos" , basis="p" )
					from matplotlib import gridspec
					gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 

					ax0 = plt.subplot(gs[0])
					kl = kline( fig , cmapcolor , ax=ax0) 

					ax1 = plt.subplot(gs[1])
					ldos.totally( args, fig, ax1 , hobj=hobj, chdat=chdat)
					ax1.set_yticklabels([])
					ax1.set_ylabel("")
					ax1.legend()
				else : 
					kl = kline( fig , cmapcolor , ax=axone) 
					print "nCBlevel : ", len(kl.cbar.levels), "[",kl.cbar.levels[0], ", ..., ", kl.cbar.levels[-1], "]"
					print "zmin,zmax : ", kl.cbar.zmin, kl.cbar.zmax
					#print "dir(kline) : ", dir(kl.cbar)
	
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
			hobj.headPltTitle( plt ) 

if args.fs :  
	ladjust=( 0.05 if (args.ladjust is None) else float(args.ladjust) )
	radjust=( 0.95 if (args.radjust is None) else float(args.radjust) )
	badjust=( 0.05 if (args.badjust is None) else float(args.badjust) )
	tadjust=( 0.95 if (args.tadjust is None) else float(args.tadjust) )
	if args.decompall :
		plt.subplots_adjust(left=0.11 , bottom=0.07, right=0.99, top=0.96, hspace=0.02 , wspace=0.02 )
#ax[1][0].xaxis.labelpad = 5
		plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, hspace=0.02 , wspace=0.02 )
	else :
		plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, wspace=0.02 , hspace=0 )
else : 
	if args.decompall :
		plt.subplots_adjust(left=0.2 , hspace=0.4 )
	else : 
		#plt.subplots_adjust(left=0.2)
		ladjust=( 0.18 if (args.ladjust is None) else float(args.ladjust) )
		radjust=( 0.93 if (args.radjust is None) else float(args.radjust) )
		badjust=( 0.20 if (args.badjust is None) else float(args.badjust) )
		tadjust=( 0.93 if (args.tadjust is None) else float(args.tadjust) )
		if args.ldos : 
			plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, wspace=0.2 , hspace=0 )
		else : 
			plt.subplots_adjust(left=ladjust, bottom=badjust, right=radjust, top=tadjust, wspace=0 , hspace=0 )
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
			if args.sublabelparameter :  slab = slab +" {}={:.6g}".format( filterlabel(args.sublabelparameter) , getattr(hobj,args.sublabelparameter) )
		        axnow.text(0.01, 0.92, slab,   
	    				horizontalalignment='left', verticalalignment='center', transform=axnow.transAxes,
		                        fontsize=fsize , color='white' , weight='semibold')



ndpi = 400 
fsave = argsavefigname
transparent = False 
if args.savetransparent : transparent = True
if args.save or args.cptmp :
	print "Saved : ", fsave + ".png"
	plt.savefig( fsave + ".png" , dpi=ndpi, transparent=transparent) 
	cptmpftn( fsave, cmd="scp", destMac=True )
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
	if args.savepdf : 
		plt.savefig( argsavefigname + ".pdf" , transparent=transparent)
	if args.saveeps : 
		plt.savefig( argsavefigname + ".eps" , transparent=transparent)
	#plt.tight_layout()
	plt.show()
