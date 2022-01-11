#!/opt/python/2.7.15/gcc-4.8.5/bin/python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

from hjl_common import *
import basislabel

import argparse

parser = argparse.ArgumentParser(description="Plot ldos with j-basis")
parser.add_argument( "mak", type=str, nargs='+', help="Data path or filename" )
parser.add_argument("-t", "--transpose", action="store_true", help="Plot the transposed data" )
parser.add_argument("--oldsoc", action="store_true" , help="SOC strength is in a old-format and multiplied by 2./3." )
parser.add_argument("-y", "--yrange", type=str, help="set yrange e.g. _y0_y1" )
parser.add_argument("-x", "--xxrange", type=str, help="set xrange e.g. _x0_x1" )
parser.add_argument("--decomp", action="store_true" , help="Plot mak with decomposition into t2g (default) or j." ) 
parser.add_argument("--decompbath", action="store_true" , help="Plot mak with decomposition into t2g (default) or j in each bath." ) 
parser.add_argument("--t2g", action="store_true" , help="Set the basis as t2g.")
parser.add_argument("--toj", "--convj", action="store_true" , help="Transform to j-basis" ) 
parser.add_argument("--nosum", action="store_true" , help="Remove the plot of sum." ) 
parser.add_argument("-n", "--noshow",  action="store_true" , help="Do not show the plot." ) 
parser.add_argument("--diffcolor", "--diffc", "--dc",  action="store_true" , help="Use different color between orbital-components and same one between trev-pair.")
parser.add_argument("--hatch", "--hat", action="store_true" , help="Use hatches.")
parser.add_argument("--ni", "--Ni", "--nimp", type=str, help="Set Ni.")
parser.add_argument("--tNC", "--tnc", "--tNc", type=str, help="Set tNC.")
parser.add_argument("--tNb", "--tnb", type=str, help="Set tNb.")
parser.add_argument("--Nb", "--nb", type=str, help="Set Nb.")
parser.add_argument("-w", "--write",  action="store_true" , help="Write new maks." ) 
parser.add_argument("-a", "--analyze",  action="store_true" , help="Analyze file.")
parser.add_argument("-i", "--iter", "--it",  action="store_true" , help="Plot .mak from all iterations.")
parser.add_argument("-e", "--energy", "--energyonly",  action="store_true" , help="Plot bath energies from all iterations.")
parser.add_argument("--iimp", "--ii", "--indimp", "--indi",  type=str , help="Specify impurity site.")
parser.add_argument("--iimpall", "--iia", "--indimpall", "--indiall",  type=str , help="Plot all impurity sites'.")
args = parser.parse_args()
#ndat = len(sys.argv) - 1
#x = [""] * ndat
#y = [""] * ndat
print "\n"
Ni	= int(args.ni)  if args.ni else 1
tNC	= int(args.tNC) if args.tNC else 6
tNb	= int(args.tNb) if args.tNb else 18
Nb	= int(args.Nb) if args.Nb else tNb/2
nbasis = tNC
basis  = "j"
chdatups = [0,1,4]
chdatdns = [4,2,5]
hatchftn = [ '' ]*nbasis
colorftn = [ 'b' ]*nbasis
if args.diffcolor  :
	colorftn = [ 'C0', 'C1', 'C1', 'C0', 'C2', 'C2' ]
if args.hatch  :
	hatchftn = [ 'o', '.', '.', 'o', '//', '//' ]
if args.t2g :
	basis = "t2g"
	chdatups = [0,2,4]
	chdatdns = [1,3,5]
	if args.diffcolor :
		colorftn = [ 'C0', 'C0', 'C1', 'C1', 'C2', 'C2' ]
	if args.hatch :
		hatchftn = [ '', '', '-', '-', '.', '.' ]

matT = np.matrix( [
		[ 0.,			0.,		-np.sqrt(1./2.),	0.,			-np.sqrt(1./2.)*1j ,	0.			], 
		[ np.sqrt(2./3.),	0.,		0.,		-1./np.sqrt(6.),	0.,			-1/np.sqrt(6.)*1j	],
		[ 0.,			np.sqrt(2./3.),	1./np.sqrt(6.),	0.,			-1./np.sqrt(6.)*1j,	0.			],
		[ 0.,			0.,		0.,		1./np.sqrt(2.),		0.,			-1./np.sqrt(2.)*1j	],
		[-1./np.sqrt(3.),	0.,		0.,		-1./np.sqrt(3.),		0.,		-1./np.sqrt(3.)*1j	],
		[ 0.,			1./np.sqrt(3.),	-1./np.sqrt(3.),	0.,			1./np.sqrt(3.)*1j,	0.		]] ) 
TRANS=False
if args.transpose :
	TRANS = True
if args.oldsoc : 
 	S = 2./ 3 * S

ndat	= len(args.mak)
hobjarr	= []
for makfile in args.mak : 
	if( makfile.find('Dir')<0 )  :
		hobj = autoheadobj( makfile , ["nb", "U", "S", "J", "D"] ) 
		setattr( hobj, "titlefoot", makfile.split('/')[-1] )
		setattr( hobj, "tNb", tNb ) 
		try : 
			setattr( hobj, "UF", float(hobj.U) ) 
		except : pass
	else :
	      	imak = makfile.find('.mak')>-1
	      	islash = makfile.rfind('/')
		if( imak )  :
			hobj = headobj( makfile[:islash] )
		else : 
			hobj = headobj( makfile )
		hobj.readMakClass()
	hobjarr.append( hobj )
titlefoot=hobj.titlefoot

mclassArr = []
for i in range( ndat ) :
	dat = args.mak[i]
	print "arg ", i, ":", dat
	if( dat.find('.mak') > -1 )  :
		datfile = dat
	else :
	   	datfile = "{}/next.mak".format(dat)
	   	print "d : ", datfile
	mclass = makclass( datfile , Ni=Ni, tNC=tNC, tNb=tNb, basis=basis ) 
	[ eldat , valdat ] = mclass.returnarr()
	mclassArr.append( mclass )


if  (hobj.tdir.find("mak")<0) and args.write :
	fwmak = open( hobj.tdir + "/next_j.mak" , 'w' ) 
	with open( datfile ) as fo :
		lines = fo.readlines()
		ndeg = int( lines[2] ) 
		for j in range( 3+ndeg ) :
			fwmak.write( lines[j] )
		for eli in el : 
			fwmak.write( "{:20.17f}\t\t\t\t".format( eli )  ) 
		fwmak.write("\n\n")
		for vlmu in vl : 
			for vlmui in vlmu : 
				fwmak.write( "{:20.17f} {:20.17f}\t".format( vlmui.real, vlmui.imag )  ) 
			fwmak.write("\n")
		fwmak.write("\n\n")
	fwmak.close()
	
	fwmakj		= open( hobj.tdir + "/next_j_sort.mak" , 'w' ) 
	fwmakjinfo	= open( hobj.tdir + "/next_j_sort.mak.info" , 'w' ) 
	fwmakt		= open( hobj.tdir + "/next_t2g_sort.mak" , 'w' ) 
	fwmaktinfo	= open( hobj.tdir + "/next_t2g_sort.mak.info" , 'w' ) 
	with open( datfile ) as fo :
		lines = fo.readlines()
		ndeg = int( lines[2] ) 
		for j in range( 3+ndeg ) :
			fwmakj.write(  lines[j] )
			fwmakt.write( lines[j] )
		sort_perm = el.argsort()
		vlsort    = vl[:, sort_perm ]
		elsort = np.copy(el)
		elsort.sort()
		vlt2gsort = vlt2g.transpose()[:, sort_perm ]
		for eli in elsort : 
			fwmakj.write(  "{:20.17f}\t\t\t\t".format( eli )  ) 
			fwmakt.write(  "{:20.17f}\t\t\t\t".format( eli )  ) 
		fwmakj.write("\n\n")
		fwmakt.write("\n\n")
	
		for vlmu in vlsort : 
			for vlmui in vlmu : 
				fwmakj.write( "{:20.17f} {:20.17f}\t".format( vlmui.real, vlmui.imag )  ) 
			fwmakj.write("\n")
		fwmakj.write("\n\n")
	
		for vlmu in vlt2gsort : 
			for vlmui in vlmu : 
				fwmakt.write( "{:20.17f} {:20.17f}\t".format( vlmui.real, vlmui.imag )  ) 
			fwmakt.write("\n")
		fwmakt.write("\n\n")
	
		for j in range( len(elsort) ) : 
			fwmakjinfo.write(  "{:9.6f}\t".format( elsort[j] )  ) 
			fwmaktinfo.write(  "{:9.6f}\t".format( elsort[j] )  ) 
			for mu in range(nbasis) :
				fwmakjinfo.write( "{:9.6f} ".format( np.absolute(   vlsort[mu][j]) )  )
				fwmaktinfo.write( "{:9.6f} ".format( np.absolute(vlt2gsort[mu][j]) )  )
			fwmakjinfo.write( "\n" ) 
			fwmaktinfo.write( "\n" ) 

	fwmakj.close()
	fwmakt.close()
	fwmakjinfo.close()
	fwmaktinfo.close()

if args.toj and args.write :
	mclass.writemak( "vjl" , "_toj" )	


if args.analyze :
	bathsym = bathsym()
	bathcubict2gpartner = []
	bathtimerevt2gpartner = []
	bathCTt2gpartner = []
	TOL = 1e-4
	for ii in range(hobj.tNb) :
		for jj in range(ii+1,hobj.tNb) :
			ndijc = 0
			ndijt = 0
			ndijct = 0
			for mu in range(nbasis) :
				nu  = bathsym.cubict2g.partner[mu][-2]
				alp = bathsym.timerevt2g.partner[mu][-2]
				bet = bathsym.timerevt2g.partner[nu][-2]
				dijcubic   = abs(vlt2gabs[ii][mu] - vlt2gabs[jj][nu])
				dijtimerev = abs(vlt2gabs[ii][mu] - vlt2gabs[jj][alp])
				dijCT      = abs(vlt2gabs[ii][mu] - vlt2gabs[jj][bet])
				if dijcubic < TOL :
					ndijc += 1 
				if dijtimerev < TOL :
					ndijt += 1 
				if dijCT < TOL :
					ndijct += 1 
				print ii, jj, ";", mu, nu, alp, bet, ";", ndijc, ndijt, ndijct, dijcubic, dijtimerev, dijCT
			if ndijc > 5 :
				print ii, vlt2gabs[ii]
				print jj, vlt2gabs[jj]
				print "ndijc : ", ndijc
				#sys.exit(1)
				bathcubict2gpartner.append( [ii,jj] ) 
			if ndijt > 5 :
				print ii, vlt2gabs[ii]
				print jj, vlt2gabs[jj]
				print "ndijt : ", ndijt 
				#sys.exit(1)
				bathtimerevt2gpartner.append( [ii,jj] ) 
			if ndijct > 5 :
				print ii, vlt2gabs[ii]
				print jj, vlt2gabs[jj]
				print "ndijct : ", ndijct 
				#sys.exit(1)
				bathCTt2gpartner.append( [ii,jj] ) 
	print "bathcubict2gpartner : ",   bathcubict2gpartner
	print "bathtimerevt2gpartner : ", bathtimerevt2gpartner
	print "bathCTt2gpartner : ",      bathCTt2gpartner
#sys.exit(1)

if args.decomp :
	fig, ax = plt.subplots(1,nbasis+1, figsize=(11,6) , sharex=True, sharey=True) 
	width = 0.09
	if args.nosum : 
	   	titleax = ax[1]
	   	labelax = ax[1]
	else : 
		vsumbottom = np.zeros( (len(el)) )
		for i in range(nbasis) :
        		ax[0].barh( el, vlt2gabs.transpose()[i], width, left=vsumbottom, alpha=(10-1.8*i)/10. , color='b' , label=blabel(basis,i) )
			for j in range(len(el)) :
				vsumbottom[j] += vlt2gabs.transpose()[i][j]
		ax[0].plot( vsumbottom, el, 'ko' , ms=2 )
		titleax = ax[0]
		labelax = ax[0]
		x1,x2,y1,y2 = titleax.axis()
		#ax.axis((-0.1,x2,y1,y2))
		a=0
		for i,j in zip(el,[0]*len(el)):
			corr = 0.02*(x2-x1) # adds a little correction to put annotation in marker's centrum
			#ax.annotate(str(a),  xy=(i + 0.01, j + corr*a))
			titleax.annotate(str(a),  xy=( (x2-x1)/2. + j + corr*a, i + 0.01 ), fontsize=5 , rotation=-45)
			a+=1
	for i in chdatups :  #0,2,4
		#print vlt2gabs.transpose()[i]
        	ax[i/2+1].barh( el, vlt2gabs.transpose()[i], width, alpha=(10-1.8*i)/10. , color='b' , label=blabel(basis,i) )
        	ax[i/2+1].legend(fontsize=8)
	for i in chdatdns :  #1,3,5
		#print vlt2gabs.transpose()[i]
        	ax[i/2+4].barh( el, vlt2gabs.transpose()[i], width, alpha=(10-1.8*i)/10. , color='g' , label=blabel(basis,i) )
        	ax[i/2+4].legend(fontsize=8)
elif args.decompbath :
	fig, ax = plt.subplots(hobj.tNb,1, figsize=(11,6) , sharex=True, sharey=True) 
	titleax = ax[0]
	labelax = ax[hobj.tNb-1]
	width = 0.09
	for bi in range(hobj.tNb) :
		for i in chdatups :  #0,2,4
        		ax[bi].bar( i/2,   vlt2gabs[bi][i], width, alpha=(10-1.8*i)/10. , color='b' , label=blabel(basis,i) )
        		#ax[bi].legend(fontsize=8)
		for i in chdatdns :  #1,3,5
        		ax[bi].bar( i/2+3, vlt2gabs[bi][i], width, alpha=(10-1.8*i)/10. , color='g' , label=blabel(basis,i) )
        		#ax[bi].legend(fontsize=8)
else :
	nraw	= 1
	ncol	= ndat
	if args.iter :
		ncol	= hobjarr[0].ncount
		nraw	= 1
	if args.energy : 
		ncol	= 1
		nraw	= 1

	print "subplot dimension : ", nraw, ncol
	fig, axarr = plt.subplots( nraw, ncol, sharex=True, sharey=True )

	if nraw*ncol == 1 : 
		axarr = [axarr]


	if args.iimpall : 
		pass
	else : 
		iimp	= int(args.iimp) if args.iimp else 0 
		if args.energy :
       			if args.iter : 
				for idat in range( 1 )  :
					hobj	= hobjarr[idat]
					EbathIter	= hobj.returnEbathIterImp(iimp=iimp)
					EbathIterT	= EbathIter.transpose()
					print "EbathIter : ", EbathIter
					print "EbathIter dim : ", np.shape( EbathIter )
					for ib in range(tNb) :
						axarr[idat].plot( range(hobj.countstart,hobj.count+1), EbathIterT[ib], 'o-', mfc='w' )
			else : 
				EbathLastArr	= np.array([ getattr( hobjarr[idat], 'makclass_i-1' ).eildat[iimp] for idat in range(ndat) ])
				EbathLastArrT	= EbathLastArr.transpose()
				print "EbathLastArr : ", EbathLastArr
				print "EbathLastArr dim : ", np.shape( EbathLastArr )
				for ib in range(tNb) :
					axarr[0].plot( range(ndat), EbathLastArrT[ib], 'o-', mfc='w' )
		elif args.iter : 
			for idat in range( 1 )  :
				hobj = hobjarr[idat]
				nmak = hobj.nmak 
				istart	 = hobj.countstart
				iend	 = hobj.count+1
				for imak in range(istart,iend) :
					getattr( hobjarr[idat], 'makclass_i{}'.format(imak) ).plotAbsImp( axarr[imak-istart] , iimp=iimp ) 
		else : 
			for idat in range( ndat )  :
				hobj = hobjarr[idat]
				nmak = hobj.nmak 
				istart	 = hobj.countstart
				iend	 = hobj.count+1
				getattr( hobjarr[idat], 'makclass_i-1' ).plotAbsImp( axarr[idat] , iimp=iimp ) 

	ax	= axarr[0]
	ax.legend()
	#mclass.plotAbs( ax )

	titleax = ax
	labelax = ax

for axnow in axarr : axnow.grid(linestyle=":")

if( len(titlefoot)>10 ) :
	titlefoot = titlefoot[-10:]
if args.mak[0].find(".mak")<0 : 
	hobj.headSetTitle( titleax ) 
else :
	titleax.set_title( titlefoot )
if TRANS is False :
	labelax.set_xlabel(r'$\sum_{\mu}$ V$_{\mu l}$',		fontsize='large' )
	labelax.set_ylabel(r'$\epsilon_l$',			fontsize='large' )
else :
	labelax.set_ylabel(r'$\sum_{\mu}$ V$_{\mu l}$',		fontsize='large' )
	labelax.set_xlabel(r'$\epsilon_l$',			fontsize='large' )

if args.yrange :
	y0 = float( args.yrange.split('_')[-2] )
	y1 = float( args.yrange.split('_')[-1] )
	titleax.set_ylim(y0,y1)
if args.xxrange :
	x0 = float( args.xxrange.split('_')[-2] )
	x1 = float( args.xxrange.split('_')[-1] )
	titleax.set_xlim(x0,x1)

print "axis : ", plt.axis()

print datfile.split('/')
fsave = ""
for fs in datfile.split('/') :
	if fs.find(".mak")>-1 : 
		fsave = fsave + "plot_" + fs 
	else :
		fsave = fsave + fs +'/'
if sys.argv[1].find("Dir")>-1 :
	pass
	#print "Saved : ", fsave
	#plt.savefig( fsave + '.png' ) 
#plt.savefig( fsave + '.pdf' ) 

if args.noshow : pass
else :
	plt.show()

#print x
#print y
