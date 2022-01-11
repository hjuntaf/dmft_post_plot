#!/opt/python/2.7.15/gcc-4.8.5/bin/python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys, os

from hjl_common import *
import basislabel

import argparse

parser = argparse.ArgumentParser(description="Plot ldos with j-basis")
parser.add_argument('ddir', metavar='D', type=str, nargs=1,
		                    help='Data path or .mak-file' ) 
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
parser.add_argument("--ni", "--nimp", type=int , help="Specify Ni.")
parser.add_argument("--imak", "-i", type=str , help="Specify i-th mak.")
parser.add_argument("--next", "--nxt", action="store_true" , help="Plot next.mak")
parser.add_argument("--sharex", "--shx", type=str, help="Set sharex in plot.")
parser.add_argument("--sharey", "--shy", type=str, help="Set sharey in plot.")
args = parser.parse_args()

ndat = 1
x   = [""] * ndat
y   = [""] * ndat
tag = [""] * ndat
ntagc = 2	# 2 components : UF, D
print "\n"

nbasis = 6

dirname = args.ddir[0]
UF=None
J=None
D=None
JT=None

i=0
print sys.argv
if( dirname.find('.mak') > -1 )  :
	hobj = autoheadobj( dirname , ["nb", "U", "S", "J", "D", "JTD"] ) 
	fname	= "%s"%dirname
	dirname = dirname[ :dirname.rfind("/") ]
	setattr( hobj, "titlefoot", "" ) 
	setattr( hobj, "tNb", int(hobj.nb)*2 ) 
	try : 
		setattr( hobj, "UF", float(hobj.U) ) 
	except : pass
	Ni = args.ni if args.ni else 1
else :
	hobj = headobj( dirname )
	Ni = hobj.Ni

UF = hobj.UF
#JT = hobj.JT
J  = hobj.J 
S  = hobj.S 
D  = hobj.D 
print "Nimp : ", Ni

fcheck = True
nmak=1
inmak=0
if inmak is 0 :
	fr = open( dirname + "/output_u{:.2f}.txt".format(UF) )
	ldata = fr.readlines()
	firstdat = ldata[2].split()
	for data in firstdat :
		if data.find("th")>-1 :
			inmak = int( data.split("th")[0] ) 
			print ":: Starting with {}-th read ::".format( inmak )  
			inmak -= 1
	
nmak += inmak
while fcheck :
	fname = "u{:.2f}_{}th.mak".format(UF,nmak) 
	if( os.system("test -f {}/{}".format(dirname,fname)) ) :
		if( nmak < 2 ) :
			if( inmak > 0 ) :
				fcheck  = False 
				print "End.\n"
			else : 
				print "{} does not exist.".format(fname)
				inmak += 1
				nmak  += 1
		else :
			fcheck  = False 
			print "End.\n"
	else :
		print "{} exists.".format(fname)
		nmak += 1


nmak -= 1+inmak
print "nmak : ", nmak

if (args.imak) or args.next : 
	nmak = 1
	if args.imak :
		inmak = int(args.imak)

sharex	= True
sharey	= True
if args.sharex : 
	if args.sharex.find("rue")>-1 : 	sharex	= True
	else :					sharex	= False
if args.sharey : 
	if args.sharey.find("rue")>-1 : 	sharey	= True
	else :					sharey	= False
fig, ax = plt.subplots(Ni,nmak, sharex=sharex, sharey=sharey)
if Ni < 2 : 
	ax = [ax]
if nmak < 2 : 
	ax = [ax]
	if Ni > 1 : 
		ax = np.array(ax).transpose()
for k in range(nmak) :
	dat = "{}/u{:.2f}_{}th.mak".format(dirname,UF,k+inmak+1)
	if args.next : 
		dat = "{}/next.mak".format(dirname)
	print k, ":", dat
	with open( dat ) as fo : 
		lines = fo.readlines()

		degenarr = np.zeros(Ni,dtype=int)
		istart = 2
		for iimp in range(Ni) :
			degendat	= int( lines[istart].split()[0] ) 
			degenarr[iimp]	= degendat
			istart		+= degendat
			eldat		= lines[istart+1].split()
			valdat = [""]*nbasis
			for i in range(nbasis) :
				valdat[i] = lines[istart+3].split()

			el = np.array( eldat , dtype=float )
			vl = np.zeros( (nbasis, len(el)) , dtype=complex)
			npel=0;nnel=0
			for i in range( len(el) ) :
				if el[i] > 0 : 
					npel += 1
				elif el[i] < 0 : 
					nnel += 1
			print "negative/positive %dth-imp : "%iimp, nnel, '/', npel 

			for i in range(nbasis) :
				for j in range(len(el)) :
					vl[i][j] = float(valdat[i][2*j]) + float(valdat[i][2*j+1])*1j
			vlabs = np.zeros( (nbasis, len(el)) , dtype=float)
			for i in range(nbasis) :
				for j in range(len(el)) :
					vlabs[i][j] = np.real( vl[i][j] * np.conjugate( vl[i][j] ) )
			
			istart		+= 11
			#print "istart : ", istart
			#sys.exit(1)

			width = 0.1
			vsumbottom = np.zeros( (len(el)) , dtype=float)
			axnow	= ax[iimp][k]
			if args.yrange :
				y0 = float( args.yrange.split('_')[-2] )
				y1 = float( args.yrange.split('_')[-1] )
				axnow.set_ylim(y0,y1)
			if args.xxrange :
				x0 = float( args.xxrange.split('_')[-2] )
				x1 = float( args.xxrange.split('_')[-1] )
				axnow.set_xlim(x0,x1)

			for i in range(nbasis) :
        			axnow.barh( el, vlabs[i], width, align='center', left=vsumbottom, alpha=(10-1.8*i)/10. , color='b' , label=r'$\mu={}$'.format(i) )
				for j in range(len(el)) :
					vsumbottom[j] += vlabs[i][j]
			axnow.plot( vsumbottom, el, 'ko' , ms=2 )
			plt.setp( axnow.xaxis.get_majorticklabels(), rotation=-45 )
			axnow.annotate( "{}".format(k+inmak+1),  xy=(0.5,0.98), xycoords=axnow.transAxes )
			x1,x2,y1,y2 = axnow.axis()
			if( k > 0 ) :
				for i in range(nbasis) :
					for j in range(len(el)) :
						ymax = 0 
						if ( ymax < vlabs[i][j] ) :
							ymax = vlabs[i][j]
				#ax[k].set_xlim( 0,vlabs[-1][-1] ) 
			a = 0
			for i,j in zip(el,[0]*len(el)):
				corr = 0.03*(x2-x1) # adds a little correction to put annotation in marker's centrum
				axnow.annotate(str(a),  xy=( (x2-x1)*0.4 + j + corr*a, i + 0.01 ), fontsize=7 , rotation=-45)
				a += 1
			#ax[k].tick_params( labelsize = 1 )
			axnow.xaxis.set_tick_params( labelsize = 8 )

for axnow in ax.flatten() : 
	axnow.grid(linestyle=":")
#ax[1].set_xlim( 0,0.5 )
titlefoot=""
fig.subplots_adjust(wspace=0.1)
ax[0][nmak/2].set_title('(U={}, J={}, S={}, D={}) {}'.format( UF, J, S, D, titlefoot ) )
ax[0][nmak/2].set_xlabel(r'$\sum_{\mu}$ V$_{\mu l}$',	fontsize='large' )
ax[0][0     ].set_ylabel(r'$\epsilon_l$',			fontsize='large' )

i=1 ;

plt.show()
