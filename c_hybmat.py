#!/usr/bin/python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
from numpy import linalg as LA

from hjl_common import *
import basislabel

import argparse

parser = argparse.ArgumentParser(description="Plot ldos with j-basis")
parser.add_argument( "mak", type=str, help="Data path or filename" )
#parser.add_argument( "Ecluster", type=str, help="Data path or filename" )
parser.add_argument("-t", "--transpose", action="store_true", help="Plot the transposed data" )
parser.add_argument("--oldsoc", action="store_true" , help="SOC strength is in a old-format and multiplied by 2./3." )
parser.add_argument("-y", "--yrange", type=str, help="set yrange e.g. _y0_y1" )
parser.add_argument("-x", "--xxrange", type=str, help="set xrange e.g. _x0_x1" )
parser.add_argument("--decomp", action="store_true" , help="Plot mak with decomposition into t2g (default) or j." ) 
parser.add_argument("--decompbath", action="store_true" , help="Plot mak with decomposition into t2g (default) or j in each bath." ) 
parser.add_argument("--basis", type=str , help="Set the basis.")
parser.add_argument("--t2g", action="store_true" , help="Set the basis as t2g.")
parser.add_argument("--toj", "--convj", action="store_true" , help="Transform to jeff-basis" ) 
parser.add_argument("--toleff", "--convleff", action="store_true" , help="Transform to leff-basis" ) 
parser.add_argument("--tot2g", "--convt", action="store_true" , help="Transform to t2g-basis" ) 
parser.add_argument("--ton", "--convn", action="store_true" , help="Transform to natural-basis" ) 
parser.add_argument("--tojfromn", "--convjfromn", action="store_true" , help="Transform to jeff-basis from natural-basis" ) 
parser.add_argument("--nosum", action="store_true" , help="Remove the plot of sum." ) 
parser.add_argument("-n", "--noshow",  action="store_true" , help="Do not show the plot." ) 
parser.add_argument("-w", "--write",  action="store_true" , help="Write the new mak." ) 
parser.add_argument("--shownumer", "--shown",  action="store_true" , help="Show the Deltanumersum." )
parser.add_argument("--showbathftn", "--showb",  action="store_true" , help="Show the Bathftn." )
parser.add_argument("--diagbathftn", "--diagb",  action="store_true" , help="Diagonalize the Bathftn." )
parser.add_argument("--comp", "--complex",  action="store_true" , help="Show the complex-components." )
parser.add_argument("--c4rot", "--c4",  action="store_true" , help="Perform the C4-rotation on the Vml." )
parser.add_argument("--c4rotwrite", "--c4write",  action="store_true" , help="Write the new mak of c4rot." ) 
parser.add_argument("-M", "--Matrix",  type=str,  help="Show the hybridization function at w(i)." ) 
parser.add_argument("--iw", "--matsu",  action="store_true" , help="Using imaginary freq.") 
parser.add_argument("--nmax", type=int , help="Set Nmax for iwn." )
parser.add_argument("--beta", type=int , help="Set beta for iwn." )
parser.add_argument("--cnnz", "--checknnzHyb",  action="store_true" , help="Check Hyb function for rea/imag frequency." )
parser.add_argument("--printnnz", "--pnz", action="store_true", help="Print the nnz-information for each freq." )
parser.add_argument("--printnnzmat", "--pnzm", "--pnzmat", action="store_true", help="Print the nnz-matrix for each freq." )
parser.add_argument("--showe", "--showecluster", action="store_true", help="Print the ecluster." )
parser.add_argument("--count", type=int , help="Choose the data of the i-th iteration." )
parser.add_argument("--sphase", "--checkphasesel", action="store_true", help="Check the relative phase of symmetry partner of t2g in the selective way.") 
parser.add_argument("--cphase", "--checkphase", "--cp", action="store_true", help="Check the relative phase of symmetry partner in t2g.") 
parser.add_argument("--cphaset2g", "--checkphaset2g", action="store_true", help="Check the relative phase of symmetry partner in t2g.") 
parser.add_argument("--cphasenosoc", "--checkphaset2gnosoc", "--cpns", action="store_true", help="Check the relative phase of symmetry partner in t2g-case without SOC.") 
parser.add_argument("--makesym", "--maksym", action="store_true", help="Make .mak file invariant under time-rev. and C4-rotation based on the input .mak file." ) 
parser.add_argument("--makesymnew", "--maksymnew", action="store_true", help="Make a new .mak file invariant under time-rev. and C4-rotation." ) 
parser.add_argument("--makephase", "--makphase", "--makph", action="store_true", help="Make .mak file from changing '-1'*[0,4,5]-basis-phases of the input .mak file." ) 
parser.add_argument("--makeconverge", "--makconv", "--makc", action="store_true", help="Save Make .mak file from converged data with sorting ." ) 
parser.add_argument("--unsort", "--usort", action="store_true", help="Show the unsorted original data." )
parser.add_argument("--nojeff", "--noj", action="store_true", help="Do not display the sorted data." )
parser.add_argument("--t2gph", "--tp", "--tph", action="store_true" , help="Change T to Taramod.")
parser.add_argument("--factor", "--fact", type=str, help="Multiply a factor to e_l and V_{mu,l}.") 
parser.add_argument("--extractandsave", "--extsave", "--extractsave", action="store_true" , help="Extract and save the raw data.")
args = parser.parse_args()

#ndat = len(sys.argv) - 1
#x = [""] * ndat
#y = [""] * ndat
print "\n"
nbasis = 6
basis  = "j"
if args.basis : basis = args.basis
if basis == "e" : nbasis = 4
if args.t2g : basis = "t2g"
if args.t2gph : 
	matT = np.matrix( [
		[ 0.,			0.,		-np.sqrt(1./2.),	0.,			-np.sqrt(1./2.)*1j ,	0.			], 
		[ np.sqrt(2./3.),	0.,		0.,		-1./np.sqrt(6.),	0.,			-1/np.sqrt(6.)*1j	],
		[ 0.,			np.sqrt(2./3.),	1./np.sqrt(6.),	0.,			-1./np.sqrt(6.)*1j,	0.			],
		[ 0.,			0.,		0.,		1./np.sqrt(2.),		0.,			-1./np.sqrt(2.)*1j	],
		[-1./np.sqrt(3.),	0.,		0.,		-1./np.sqrt(3.),		0.,		-1./np.sqrt(3.)*1j	],
		[ 0.,			1./np.sqrt(3.),	-1./np.sqrt(3.),	0.,			1./np.sqrt(3.)*1j,	0.		]] ) 
	matTlefft2g = np.matrix( [
		[ np.sqrt(2),			0.,		0.,		0.,			0.,			0.			], 
		[ 0.,			np.sqrt(2),		0.,		0.,			0.,			0.			], 
		[ 0.,				0.,		-1.,		0.,			-1j,			0.			], 
		[ 0.,				0.,		0.,		-1.,			0.,			-1j			], 
		[ 0.,				0.,		1.,		0.,			-1j,			0.			], 
		[ 0.,				0.,		0.,		1.,			0.,			-1j			]] ) / np.sqrt(2.)
else :
	matT = np.matrix( [
		[ 0.,			0.,		np.sqrt(1./2.),	0.,			np.sqrt(1./2.)*1j ,	0.			], 
		[ np.sqrt(2./3.),	0.,		0.,		-1./np.sqrt(6.),	0.,			-1/np.sqrt(6.)*1j	],
		[ 0.,			np.sqrt(2./3.),	1./np.sqrt(6.),	0.,			-1./np.sqrt(6.)*1j,	0.			],
		[ 0.,			0.,		0.,		1./np.sqrt(2.),		0.,			-1./np.sqrt(2.)*1j	],
		[ 1./np.sqrt(3.),	0.,		0.,		1./np.sqrt(3.),		0.,			1./np.sqrt(3.)*1j	],
		[ 0.,			-1./np.sqrt(3.),	1./np.sqrt(3.),	0.,			-1./np.sqrt(3.)*1j,	0.		]] ) 
	matTlefft2g = np.matrix( [
		[ 1.,			0.,		0.,		0.,			0.,			0.			], 
		[ 0.,			1.,		0.,		0.,			0.,			0.			], 
		[ 0.,			0.,		1.,		0.,			1j,			0.			], 
		[ 0.,			0.,		0.,		1.,			0.,			1j			], 
		[ 0.,			0.,		1.,		0.,			-1j,			0.			], 
		[ 0.,			0.,		0.,		1.,			0.,			-1j			]] ) / np.sqrt(2.)
matTjeffleff = np.dot( np.conjugate(matT),  matTlefft2g.transpose() )
a = 0.980031
b = 0.198847
matTjnat = np.matrix( [
		[ 1, 0, 0, 0, 0, 0 ],
		[ 0, 0, 0, 1, 0, 0 ],
		[ 0,-a, 0, 0,-b, 0 ],
		[ 0, 0, a, 0, 0,-b ],
		[ 0,-b, 0, 0, a, 0 ],
		[ 0, 0, b, 0, 0, a ]
		] ) 

#for i in range( len(dataset) ) :
dat = args.mak
print "arg : ", dat
if( dat.find('.mak') > -1 )  :
	datfile = dat
   	print "d : ", datfile
	if dat.find("Dir") > -1 : 
		hobj = headobj( dat.split("/")[-2] )
		try : 
			hobj.readParameters()
		except :
			pass
		try : 
			hobj.readEcluster()
			matTjnat = hobj.returnNatTransf()
		except : pass
	mclass = makclass( datfile , basis=basis ) 
	if args.t2gph :  mclass.setTt2gjeff( matT ) 
	[ eldat , valdat ] = mclass.returnarr()
else :
   	datfile = "{}/next.mak".format(dat)
	hobj = headobj( dat )
	try : 
		hobj.readParameters()
	except :
		pass
	try : 
		hobj.readEcluster()
		matTjnat = hobj.returnNatTransf()
	except : pass
	if args.count : 
		count = int(args.count)
		datfile = "{}/u{:.2f}_{}th.mak".format(dat, hobj.UF, count )
   	print "d : ", datfile
	mclass = makclass( datfile,basis=basis ) 
	if args.t2gph :  mclass.setTt2gjeff( matT ) 
	[ eldat , valdat ] = mclass.returnarr()

mclass.convnparray()
mclass.showNelpn()
vml = mclass.npvmldat ; vl=vml
el = mclass.npeldat
[eln,elp] = [ mclass.neln, mclass.nelp ]
if basis =="e" :
	pass
else :
	deltanumersum = mclass.returnDeltanumersum()

print "vl ", np.shape(vl), " :"
#print vl
print "el ", np.shape(el), " :"
#print el
if args.extractandsave : 
	print "vl ", np.shape(vl), " :"
	print vl
	print "el ", np.shape(el), " :"
	print el
	fname_param_save  = hobj.tdir + "/D{:g}_vml.dat".format(hobj.D) ; print "Saving : ", fname_param_save 
	np.savetxt( fname_param_save, vl ) 
	dumdat  = np.loadtxt( fname_param_save, dtype=complex ) 
	print "cheking :: ", np.abs( vl - dumdat ) 
	fname_param_save  = hobj.tdir + "/D{:g}_el.dat".format(hobj.D) ; print "Saving : ", fname_param_save 
	np.savetxt( fname_param_save, el ) 
	dumdat  = np.loadtxt( fname_param_save )
	print "cheking :: ", np.abs( el - dumdat ) 
	exit(1)

if args.shownumer : 
	#print "Deltanumersum in jeff's (before rotation) :  "
	#print deltanumersum
	#print "Deltanumersum in t2g's (after rotation) :  "
	#print deltanumersumt2g
	deltanumersumt2g = transfMjefftot2g( deltanumersum, matT ) 
	deltanumersumnat = transfMt2gtojeff( deltanumersum, matTjnat ) 
	print "deltanumersum   " ; printArrMat( deltanumersum    ) 
	print "deltanumersumt2g" ; printArrMat( deltanumersumt2g ) 
	print "deltanumersumnat" ; printArrMat( deltanumersumnat ) 
#def printArrMat( arrMat )  :
#	for i in range(nbasis) :
#		for j in range(nbasis) :
#			print "{:17.5f} ".format(arrMat[i][j]),
#		print ""
#def printArrMat( arrMat )  :
#	for i in range(nbasis) :
#		for j in range(nbasis) :
#			print "{:17.5f} ".format(arrMat[i][j]),
#		print ""
#deltawn =  mclass.returnCalcDelta( -6 ) 
#print "deltawn" ; printArrMat( deltawn ) 

def printCompVmlel( vl, el )  :
	for j in range(len(vl[0])) :
		for i in range(nbasis) :
			print "{:19.6f} ".format( vl[i][j] ),
		print "| {:9.6f}".format( el[j] )

def printAbsVmlel( vl, el )  :
	for j in range(len(vl[0])) :
		for i in range(nbasis) :
			print "{:9.6f} ".format( np.abs(vl[i][j]) ),
		print "| {:9.6f}".format( el[j] )

def printAbsVmlelnpmat( vl, el )  :
	for j in range(np.shape(vl)[1]) :
		for i in range(nbasis) :
			print "{:9.6f} ".format( np.abs(vl[i,j]) ),
		print "| {:9.6f}".format( el[j] )

def printAbsSumVmlel( vl, el )  :
	for j in range(len(vl[0])) :
		for i in range(nbasis) :
			print "{:9.6f} ".format( np.abs(vl[i][j]) ),
		print "| {:9.6f} | ".format( el[j] ),
		dumsum = 0 
		for i in range(nbasis/2) :
			print "{:9.6f} ".format( np.abs(vl[2*i][j]) + np.abs(vl[2*i+1][j]) ),
			dumsum += np.abs(vl[2*i][j]) + np.abs(vl[2*i+1][j])
		print "| ",
		print "{:9.6f} ".format( dumsum )

if args.unsort : 
	print "Vml | el (unsorted,read) : " ; printAbsVmlel(vml,el)
	if args.comp : 
		printCompVmlel( mclass.npvmldat ,el )
wrtag = ""
if args.factor : 
	fact = float(args.factor)
	mclass.npvmldat *= fact
	vml = mclass.npvmldat ; vl=vml
	mclass.npeldat *= fact
	el = mclass.npeldat
	wrtag = "_fact{}".format(args.factor)
#print "Vjl | el : " ; printAbsVmlel(vl,el)
#vlsort , elsort = sortVmlel( vl, el ) 
vlsort , elsort = mclass.sortVmlel( )
if args.nojeff or args.tot2g :  pass 
else : print "Vml | el (sorted,read) : " ; printAbsVmlel(vlsort,elsort)
if args.comp : 
	printCompVmlel( mclass.npvmlsort ,elsort)
#print "Vml | el (sorted,read) : " ; printAbsSumVmlel(vlsort,elsort)
#
if args.tot2g :
	mclass.toVtl()
	mclass.toVtlusort()
	if args.unsort : 
		print "Vtl | el (unsorted): " ; printAbsVmlel( mclass.npvtldatusort ,el )
	else : 
		print "Vtl | el (sorted): " ; printAbsVmlel( mclass.npvtldat ,elsort)
	if args.comp : 
		if args.unsort : 
			printCompVmlel( mclass.npvtldatusort ,el )
		else :
			printCompVmlel( mclass.npvtldat ,elsort)
		a = np.abs( mclass.npvtldat[0][0] )
		b = np.abs( mclass.npvtldat[1][0] )
		c = mclass.npvtldat[1][0] /  mclass.npvtldat[1][1] 
		#print "c : ", c 
		#print "0 : ", ["{:19.6f}".format( mclass.npvtldat[0][0] ) , "{:19.6f}".format( mclass.npvtldat[1][0] ) ]
		#print "1 : ", ["{:19.6f}".format( mclass.npvtldat[0][1] ) , "{:19.6f}".format( mclass.npvtldat[1][1] ) ]
		#print "0' : ", ["{:19.6f}".format( mclass.npvtldat[0][0]          - mclass.npvtldat[0][1]*c ) , "{:19.6f}".format( mclass.npvtldat[1][0]          - mclass.npvtldat[1][1]*c ) ]
		#print "1' : ", ["{:19.6f}".format( mclass.npvtldat[0][0]*np.conj(c) + mclass.npvtldat[0][1]   ) , "{:19.6f}".format( mclass.npvtldat[1][0]*np.conj(c) + mclass.npvtldat[1][1] ) ]
	if args.write : 
		mclass.writemaksort( "npvtldat" , "_t2g"+wrtag )	
if args.toj :
	print "T : "
	printNpMat( mclass.Tt2gjeff )
	print ""
	mclass.toVjl()
	print "Vjl | el (sorted,read) : " ; printAbsVmlel(mclass.npvjldat,elsort)
	if args.comp : 
		print "Vjl | el (sorted,read,complex) : " ; printCompVmlel(mclass.npvjldat,el)
	if args.write : 
		mclass.writemaksort( "npvjldat" , "_jeff"+wrtag )	
if args.ton :
	mclass.toVnl()
	#print "Vnl | el : " ; printAbsVmlel(mclass.npvnldat,el)
	mclass.npvnldat , mclass.npeldat = sortVmlel( mclass.npvnldat, el ) 
	print "Vnl | el (sorted): " ; printAbsVmlel( mclass.npvnldat ,elsort)
	if args.comp : 
		printCompVmlel( mclass.npvnldat ,elsort)
	if args.write : 
		mclass.writemaksort( "npvnldat" , "_nat"+wrtag )	
if args.tojfromn :
	mclass.toVjlfromnat()
	if args.comp : 
		printCompVmlel( mclass.npvjldat ,elsort)
	if args.write : 
		mclass.writemaksort( "npvjldat" , "_jefffromnat"+wrtag )	
if args.toleff :
	mclass.toVfromjeff( matTjeffleff )
	if args.comp : 
		printCompVmlel( mclass.npvfromjeff ,el )
	if args.write : 
		mclass.writemak( "npvleff" , "_lefffromjeff"+wrtag )	
	print "Vleffl | el (unsorted): " ; printAbsVmlel( mclass.npvfromjeff ,el )
if args.write : 
	mclass.sortVmlel()
	mclass.writemaksort( "npvmlsort" , "_sort"+wrtag )	

def dotArrNpmat( Arr, N , dim ) :
	resArr = np.zeros( (dim,dim) , dtype=complex)
	for mu in range(dim) :
		for nu in range(dim) :
			for lu in range(dim) :
				resArr[mu][nu] += Arr[mu][lu] * N[lu,nu] 
	return resArr
def dotArrArr( Arr, Arr2 , dim ) :
	resArr = np.zeros( (dim,dim) , dtype=complex)
	for mu in range(dim) :
		for nu in range(dim) :
			for lu in range(dim) :
				resArr[mu][nu] += Arr[mu][lu] * Arr2[lu][nu] 
	return resArr

if args.mak.find("Dir")>-1 : 
	try : 
		Eclusterjeff= hobj.readEclustersoc()
		Eclustert2g = hobj.readEclustersoct2g()
	except : pass
	if args.showe :
		print "Ecluster : " ;  printArrMat( Eclusterjeff ) 
		if args.tot2g : 
			print "Eclustert2g : " ;  printArrMat( Eclustert2g ) 
		if args.toj : 
			print "Eclusterjeff : " ;  printArrMat( Eclusterjeff ) 
		if args.ton : 
			Eclusternat = hobj.readEclustersocnat()
			print "Eclusternat : " ;  printArrMat( Eclusternat ) 
		if args.toleff : 
			Eclusterleff = dotArrArr( matTjeffleff.transpose() , dotArrArr(Eclusternat,np.conjugate(matTjeffleff)) )
			print "Eclusterleff : " ;  printArrMat( Eclusterleff ) 

bathftn = mclass.returnBathftn()
if args.showbathftn : 
	print "Bathftn : "
	printAbsArrMatdim( bathftn , len(bathftn) ) 
if args.diagbathftn : 
	w, v  = LA.eig( np.matrix(bathftn) )
	print "w : \n", w.real
	vmldiagb =  np.dot( np.matrix(vlsort) , v ) 
	#testdiagb =  np.dot( v.transpose().conjugate() , np.dot( np.matrix(bathftn) , v ) ) ; print np.diag( testdiagb ).real 
	elsortmat = np.diag( np.array(elsort) ) 
	elsortdiagb =  np.dot( v.transpose().conjugate(), np.dot( elsortmat , v ) )
	print "elsortdiagb : \n", elsortdiagb[0]

if args.c4rot : 
	c4VmlMat = c4rotVml( mclass.npvmlsort, "j" ) 
	#print "c4VmlMat : ", c4VmlMat
	print "c4VmlMat : " 
	printAbsVmlelnpmat( c4VmlMat,elsort)
	if args.comp :
		printCompVmlelnp( c4VmlMat,  elsort ) 
	if args.tot2g : 
		c4VmlMatt2g = c4rotVml( mclass.npvtldat, "t" ) 
		print "c4VmlMatt2g : " 
		printAbsVmlelnpmat( c4VmlMatt2g,elsort)
		if args.comp :
			printCompVmlelnp( c4VmlMatt2g,  elsort ) 
	if args.c4rotwrite : 
		setattr( mclass , "c4VmlMat" , c4VmlMat )
		mclass.writenpmaksort( "c4VmlMat" , "_c4rot"+wrtag )	

if args.mak.find("Dir")>-1  : hobj = headobj( args.mak.split("/")[0] ) ; hobj.readParameters() ; nmax = hobj.Nmax ; beta = hobj.beta 
else : nmax = 512 ; beta=128
if args.nmax : nmax = args.nmax
if args.beta : beta = args.beta
print "beta : ", beta
print "nmax : ", nmax
if args.Matrix :
        #printCompVmlel( mclass.npvmldat ,mclass.npeldat)
        if args.iw :
                a = int( args.Matrix )
                wi = 1./beta * (2*a+1) * np.pi
                print "wi ; ", wi
                print "Deltaiw {}".format( int(args.Matrix) ) ; printArrMat( mclass.returnDeltaiw(wi)    )
        else :
                a = int( args.Matrix )
		wn = np.linspace(-6,6,1025)
                print "wi ; ", wn[a]
                print "Deltaw {}".format( int(args.Matrix) ) ; printArrMat( mclass.returnDeltaw(wn[a])    )

checkcomparr = [ 			#for purely 2-4 block-diagonal in jeff
	[0,1], [0,2], [0,4], [0,5],
	[1,3], 
	[2,3], 
	[3,4], [3,5],
	]
def checknnz( Mcomp, nzind, nnz, nzMat , a )  :
	for ij in checkcomparr :
		i = ij[0]
		j = ij[1]
		nnz = 0  
		d  = Mcomp[i][j] 
		t  = Mcomp[j][i] 
		nzind = ""
		if abs(d.real)> 1e-4 :
			nzMat[i][j] += d.real
			nzind += "[r{}{}] ".format(i,j)
			nnz+=1 
		if abs(d.imag)> 1e-4 :
			nzMat[i][j] += d.imag
			nzind += "[i{}{}] ".format(i,j)
			nnz+=1 
		if abs(t.real)> 1e-4 :
			nzMat[j][i] += t.real
			nzind += "[r{}{}] ".format(j,i)
			nnz+=1 
		if abs(t.imag)> 1e-4 :
			nzMat[j][i] += t.imag
			nzind += "[i{}{}] ".format(j,i)
			nnz+=1 
		if nnz > 0 : 
			nzind2 = ""
			nzind2 = "i{}; ".format(a) + nzind2 
			nzind2 += "absij={}; ".format(abs(d))
			nzind2 += "absji={};".format(abs(t))
			nzind2 += "nnz={};\t".format(nnz)
			nzind  = nzind2+nzind
			if args.printnnz : 
				print nzind
			if args.printnnzmat : 
				printArrMatdim( Mcomp, nbasis ) 
if args.cnnz :
	nnzsum = 0
	nzMat  = np.array( [ [ 0 for i in range(6) ] for j in range(6) ] , dtype=complex) 
	hybw = []
        if args.iw :
		nplot = nmax
		ii = range(nplot) 
		wi = range(nplot) 
		for a in ii : 
                	wi[a]  = 1./beta * (2*a+1) * np.pi
                	Diw =  mclass.returnDeltaiw(wi[a]) 
			nzind = "" ; nnz=0 ; 
			checknnz( Diw, nzind, nnz, nzMat , a)  ; nnzsum += nnz
                print "iwn : {} to {} ".format( wi[0] , wi[-1] )
		wn = wi 
        else :
		nplot = 1025
		ii = range(nplot) 
		wn = np.linspace(-6,6,nplot)
		for a in ii : 
                	Dw  =  mclass.returnDeltaw(wn[a])  
			nzind = "" ; nnz=0 ; 
			checknnz( Dw, nzind, nnz, nzMat , a)  ; nnzsum += nnz
			#hybw.append( Dw[0][0].real )
                print "wn : {} to {} ".format( wn[0], wn[-1] )
	
	#print "hybw : ", len(hybw)
	#fig, ax = plt.subplots()
	#ax.plot( hybw, wn, 'k-o' , ms=2 )
	#ax.set_ylim( -2,2)
	#plt.show()



#print "(nzMat)" ; printArrMat( nzMat ) 

def pp(aa) : 
	if aa==np.nan : 
		print "{:30}".format(aa),
	elif np.abs(aa.real)<1e-15 and np.abs(aa.imag)<1e-15 : 
		print "{:>30}".format("."),
	else : 
		print "{:>30.8g}".format(aa),
def ppr(aa) : 
	if aa==np.nan : 
		print "{:22}".format(aa),
	else : 
		print "{:19.8f}".format(aa),

def convRadianinpi( z ) :
	try : 
		return np.angle(z) / np.pi 
	except : 
		return np.nan
def relativeangle( v, mu1, l1, mu2, l2 , p="") :
	try :
		v2=v[mu2][l2]
		if np.abs(v2)<1e-6 : 
			dd=np.nan
		else : 
			dd= convRadianinpi( v[mu1][l1]/v[mu2][l2] )
		if p.find("print") >-1 : 
			pp(dd) 
		return dd
	except : 
		if p.find("print") >-1 : 
			print np.nan,
		return np.nan
def printrelativeangle( v, mu1, l1, mu2, l2 ) :
	try :
		v2=v[mu2][l2]
		if np.abs(v2)<1e-6 : 
			dd=np.nan
		else : 
			dd = convRadianinpi( v[mu1][l1] / v[mu2][l2] )
		pp(dd)
		return dd
	except : 
		return np.nan
def relativeconjangle( v, mu1, l1, mu2, l2 , p="" ) :
	try : 
		v2=v[mu2][l2]
		if np.abs(v2)<1e-6 : 
			dd=np.nan
		else : 
			dd = convRadianinpi( v[mu1][l1] / (v[mu2][l2].conjugate()) )
		if p.find("print") >-1 : 
			pp(dd)
		return dd
	except : 
		if p.find("print") >-1 : 
			print np.nan,
		return np.nan
def printrelativeconjangle( v, mu1, l1, mu2, l2 ) :
	try : 
		v2=v[mu2][l2]
		if np.abs(v2)<1e-6 : 
			dd=np.nan
		else : 
			dd = convRadianinpi( v[mu1][l1] / (v[mu2][l2].conjugate()) )
		pp(dd)
		return dd
	except : 
		return np.nan
def printdiffangleC4( v, mu1, l1, mu2, l2 ) :
	try : 
		th1 =  convRadianinpi( v[mu1][l1] / (v[mu2][l2]) )
		th2 =  convRadianinpi( v[mu2][l1] / (v[mu1][l2]) )
		pp(th1-th2)
		return th1-th2
	except : 
		return np.nan
def printdiffconjangleTrev( v, mu1, l1, mu2, l2 ) :
	try : 
		th1 =  convRadianinpi( v[mu1][l1] / (v[mu2][l2].conjugate()) )
		th2 =  convRadianinpi( v[mu2][l1] / (v[mu1][l2].conjugate()) )
		pp(th1-th2)
		return th1-th2
	except : 
		return np.nan
def ratio( v, mu1, l1, mu2, l2 , p="" ) :
	try :
		if np.abs( v[mu1][l1] ) < 1e-3 or np.abs( v[mu2][l2] ) < 1e-3 : 
			dd = np.nan 
		else :
			dd = v[mu1][l1] / v[mu2][l2]
		if p.find("print")>-1 :
			ppr(dd)
		return dd
	except : 
		return np.nan
def conjratio( v, mu1, l1, mu2, l2 , p="" ) :
	try :
		if np.abs( v[mu1][l1] ) < 1e-3 or np.abs( v[mu2][l2] ) < 1e-3 : 
			dd = np.nan 
		else :
			dd = v[mu1][l1] / (v[mu2][l2].conjugate())
		if p.find("print")>-1 :
			ppr(dd)
		return dd
	except : 
		return np.nan
def printratio( v, mu1, l1, mu2, l2 ) :
	try :
		if np.abs( v[mu1][l1] ) < 1e-3 or np.abs( v[mu2][l2] ) < 1e-3 : 
			dd = np.nan 
		else :
			dd = v[mu1][l1] / v[mu2][l2]
		ppr(dd)
		return dd
	except : 
		return np.nan
def printratioabs( v, mu1, l1, mu2, l2 ) :
	try :
		if np.abs( v[mu1][l1] ) < 1e-3 or np.abs( v[mu2][l2] ) < 1e-3 : 
			dd = np.nan 
		else :
			dd = np.abs( v[mu1][l1] / v[mu2][l2] )
		pp(dd)
		return dd
	except : 
		return np.nan

def tf( val1, val2 ) :
	dd =  np.abs(val1-val2)
	if dd <1e-4 : 
		print "{:10}".format("True"),
	else : 
		print "{:10}".format("False"),

def pairsum(v, mu, nu, l ) :
	return v[mu][l] * v[nu][l].conjugate() +  v[mu][l+1] * v[nu][l+1].conjugate()

def printpairsum(v, mu, nu, l ) :
	dd = v[mu][l] * (v[nu][l].conjugate()) +  v[mu][l+1] * (v[nu][l+1].conjugate())
	ppr(dd) 
	return dd

def lsum( v, mu,nu, larr ) : 
	lsum = 0.
	for ll in larr : 
		lsum += v[mu][ll] * v[nu][ll].conjugate()
	return lsum

def cphaseftnselectivet2g( v , elsort )  :
	print "Time-reversal condition :: cross-difference of relative-conjugate angle (pi pi pi)"
	for l in range(nbath/2) :
		print "{:9.6f} {:5d} : ".format(elsort[2*l],2*l) ,
		trat0 = [ printrelativeconjangle( v, 0, 2*l, 1, 2*l+1 ) , printrelativeconjangle( v, 0, 2*l+1, 1, 2*l )  ] 
		trat1 = [ printrelativeconjangle( v, 2, 2*l, 3, 2*l+1 ) , printrelativeconjangle( v, 2, 2*l+1, 3, 2*l )  ] 
		trat2 = [ printrelativeconjangle( v, 4, 2*l, 5, 2*l+1 ) , printrelativeconjangle( v, 4, 2*l+1, 5, 2*l )  ] 
		print "|",
		pp( trat0[0] - trat0[1] ) 
		pp( trat1[0] - trat1[1] ) 
		pp( trat2[0] - trat2[1] ) 
		print ""
	print "Time-reversal condition :: ratio of bath up-dn (-D* D D -D* D -D*)" 
	for l in range(nbath/2) :
		print "{:9.6f} {:5d} : ".format(elsort[2*l],2*l) ,
		a0=ratio( v, 0, 2*l,   0, 2*l+1 ,p="print") 
		b0=ratio( v, 1, 2*l+1, 1, 2*l   ,p="print") 
		a1=ratio( v, 2, 2*l+1, 2, 2*l   ,p="print") 
		b1=ratio( v, 3, 2*l,   3, 2*l+1 ,p="print") 
		a2=ratio( v, 4, 2*l+1, 4, 2*l   ,p="print") 
		b2=ratio( v, 5, 2*l,   5, 2*l+1 ,p="print") 
		if np.abs(a0)>1e-4 :   D = -a0.conjugate() 
		elif np.abs(a1)>1e-4 : D = a1 
		tf( D, -a0.conjugate() ) 
		tf( D, b0 )
		tf( D, a1 ) 
		tf( D, -b1.conjugate() ) 
		tf( D, a2 ) 
		tf( D, -b2.conjugate() ) 
		print ""
	print "Time-reversal condition :: ratio of orbital up-dn (-D* D* D*)" 
	for l in range(nbath) :
		print "{:9.6f} {:5d} : ".format(elsort[l],l) ,
		a0=conjratio( v, 0, l, 1, l ,p="print") 
		a1=conjratio( v, 3, l, 2, l ,p="print") 
		a2=conjratio( v, 5, l, 4, l ,p="print") 
		if np.abs(a0)>1e-4 :   D = a0
		elif np.abs(a1)>1e-4 : D = -a1
		tf( D, a0 ) 
		tf( D, -a1 ) 
		tf( D, -a2 ) 
		print ""
	print "C4-rotation condition :: angle of (yz sigma / zx sigma) "
	for l in range(nbath) :
		print "{:9.6f} {:5d} : ".format(elsort[l],l) ,
		printrelativeangle( v, 2, l, 4, l ) 
		printrelativeangle( v, 3, l, 5, l ) 
		print ""
	print "C4-rotation condition :: sum_l(pair) of Vmul x Vnul^* (0 for all) :: [0,1], [2,3], [4,5], [0,2], [2,5], [3,4]"
	for l in range(nbath/2) :
		print "{:9.6f} {:5d} : ".format(elsort[2*l],2*l) ,
		sumpair = [ [0,1], [2,3], [4,5], [0,2], [2,5], [3,4] ]
		[
		printpairsum( v, munu[0], munu[1], 2*l )
		for munu in sumpair
		]
		print ""
	print "C4-rotation condition :: exchanging ratio and argument between bath-orbs :: [mu1,ll1;mu2;ll2] , [+1,0;+1,0], [0,+1;0,+1], [+1,+1;+1,+1]"
	Nsym = 3 
	tNC  = 6
	for l in range(Nsym) :
		mu1 = 2
		mu2 = 4 
		ll1 = tNC*l + mu1 ;
		ll2 = tNC*l + mu2 ;
		print "{:9.6f} {:5d} : ".format(elsort[ll1],ll1) ,
		printratio( v, mu1  , ll1  , mu2  , ll2   ) 
		printratio( v, mu1+1, ll1  , mu2+1, ll2   ) 
		printratio( v, mu1  , ll1+1, mu2  , ll2+1 ) 
		printratio( v, mu1+1, ll1+1, mu2+1, ll2+1 ) 
		printrelativeangle( v, mu1  , ll1  , mu2  , ll2   ) 
		printrelativeangle( v, mu1+1, ll1  , mu2+1, ll2   ) 
		printrelativeangle( v, mu1  , ll1+1, mu2  , ll2+1 ) 
		printrelativeangle( v, mu1+1, ll1+1, mu2+1, ll2+1 ) 
		print ""
	print "C4-rotation condition :: orbital up-dn ratio and argument in parallel bath-orbs :: [mu1,ll1;mu2;ll1] , [+1,0;+1,0], [0,+1;0,+1], [+1,+1;+1,+1]"
	for l in range(Nsym) :
		mu1 = 2
		mu2 = 4 
		ll1 = tNC*l + mu1 ;
		ll2 = tNC*l + mu2 ;
		pairarr = [	[mu1,	ll1,	mu2,	ll1	],
				[mu1+1,	ll1,	mu2+1,	ll1	], 
				[mu1  ,	ll2,	mu2,	ll2	], 
				[mu1+1,	ll2,	mu2+1,	ll2	] ]
		print "{:9.6f} {:5d} : ".format(elsort[ll1],ll1) ,
		[ printratio(		v, pa[0], pa[1] , pa[2], pa[3] ) for pa in pairarr ]
		[ printrelativeangle(	v, pa[0], pa[1] , pa[2], pa[3] ) for pa in pairarr ]
		print ""
	print "C4-rotation condition :: ratio and argument between parallel bath-orbs in each orbital :: [mu1,ll1;mu1;ll2] , [+1,0;+1,0], [0,+1;0,+1], [+1,+1;+1,+1]"
	for l in range(Nsym) :
		mu1 = 2
		mu2 = 4 
		ll1 = tNC*l + mu1 ;
		ll2 = tNC*l + mu2 ;
		pairarr = [	[mu1,	ll1,	mu1,	ll2	],
				[mu1+1,	ll1,	mu1+1,	ll2	], 
				[mu2  ,	ll1,	mu2,	ll2	], 
				[mu2+1,	ll1,	mu2+1,	ll2	] ]
		print "{:9.6f} {:5d} : ".format(elsort[ll1],ll1) ,
		[ printratio(		v, pa[0], pa[1] , pa[2], pa[3] ) for pa in pairarr ]
		[ printrelativeangle(	v, pa[0], pa[1] , pa[2], pa[3] ) for pa in pairarr ]
		print ""
	print "C4-rotation condition :: 2nd part of exchanging ratio and argument between bath-orbs :: [mu1,ll2;mu1;ll1] , [+1,0;+1,0], [0,+1;0,+1], [+1,+1;+1,+1]"
	Nsym = 3 
	tNC  = 6
	for l in range(Nsym) :
		mu1 = 2
		mu2 = 4 
		ll1 = tNC*l + mu1 ;
		ll2 = tNC*l + mu2 ;
		print "{:9.6f} {:5d} : ".format(elsort[ll1],ll1) ,
		printratio( v, mu1  , ll2  , mu2  , ll1   ) 
		printratio( v, mu1+1, ll2  , mu2+1, ll1   ) 
		printratio( v, mu1  , ll2+1, mu2  , ll1+1 ) 
		printratio( v, mu1+1, ll2+1, mu2+1, ll1+1 ) 
		printrelativeangle( v, mu1  , ll2  , mu2  , ll1   ) 
		printrelativeangle( v, mu1+1, ll2  , mu2+1, ll1   ) 
		printrelativeangle( v, mu1  , ll2+1, mu2  , ll1+1 ) 
		printrelativeangle( v, mu1+1, ll2+1, mu2+1, ll1+1 ) 
		print ""

def cphaseftnselectivejeff( v , elsort )  :
	print "Time-reversal condition :: cross-difference of relative-conjugate angle (pi pi pi)"
	for l in range(nbath/2) :
		print "{:9.6f} {:5d} : ".format(elsort[2*l],2*l) ,
		trevpair = [ [0,3],[1,2],[4,5] ]
		[
		pp( relativeconjangle( v, munu[0], 2*l, munu[1], 2*l+1 ) - relativeconjangle( v, munu[0], 2*l+1, munu[1], 2*l )  ) 
		for munu in trevpair 
		]
		print ""
	print "C4-rotation condition :: ratio of orbital up-dn (S* -S* S*)" 
	for l in range(nbath) :
		print "{:9.6f} {:5d} : ".format(elsort[l],l) ,
		a0=conjratio( v, 3, l, 0, l ,p="print") 
		a1=conjratio( v, 1, l, 2, l ,p="print") 
		a2=conjratio( v, 4, l, 5, l ,p="print") 
		if np.abs(a0)>1e-4 :   S = a0.conjugate()
		elif np.abs(a1)>1e-4 : S = -a1.conjugate()
		tf( S,  a0.conjugate() ) 
		tf( S, -a1.conjugate() ) 
		tf( S,  a2.conjugate() ) 
		print ""
	print "C4-rotation condition :: ratio of bath up-dn (S -S* S -S* -S* S)" 
	for l in range(nbath/2) :
		print "{:9.6f} {:5d} : ".format(elsort[2*l],2*l) ,
		a0=ratio( v, 0, 2*l+1, 0, 2*l   ,p="print") 
		a1=ratio( v, 1, 2*l  , 1, 2*l+1 ,p="print") 
		b0=ratio( v, 2, 2*l+1, 2, 2*l   ,p="print") 
		b1=ratio( v, 3, 2*l  , 3, 2*l+1 ,p="print") 
		a2=ratio( v, 4, 2*l  , 4, 2*l+1 ,p="print") 
		b2=ratio( v, 5, 2*l+1, 5, 2*l   ,p="print") 
		if np.abs(a0)>1e-4 :   S =  a0
		elif np.abs(a1)>1e-4 : S = -a1.conjugate()
		tf( S,  a0 ) 
		tf( S, -a1.conjugate() ) 
		tf( S,  b0 ) 
		tf( S, -b1.conjugate() ) 
		tf( S, -a2.conjugate() ) 
		tf( S,  b2 ) 
		print ""
	print "C4-rotation condition :: sum_l(pair) of Vmul x Vnul^* (0 for all)"
	for l in range(nbath/2) :
		print "{:9.6f} {:5d} : ".format(elsort[2*l],2*l) ,
		sumpair = [ [0,3],[1,2],[4,5], [1,5], [0,1],[0,4], [3,5],[2,3] ]
		[
		printpairsum( v, munu[0], munu[1], 2*l )
		for munu in sumpair
		]
		print ""

def cphaseftnNoSOC( v , elsort )  :
	print "\tr(xyu,l/xyd,l+1), r(yzu,l/yzd,l+1), r(zxu,l/zxd,l+1)"
	for l in range(nbath/2) :
		print "{:9.6f} {:5d} : ".format(elsort[2*l],2*l) ,
		printrelativeangle( v, 0, 2*l, 1, 2*l+1 ) 
		printrelativeangle( v, 2, 2*l, 3, 2*l+1 ) 
		printrelativeangle( v, 4, 2*l, 5, 2*l+1 ) 
		print ""
	print "\tr(yzu,l/zxu,l), r(yzu,l/yzu,l+2), r(yzu,l/zxu,l+2)"
	for l in range(nbath/6) :
		ll1 = l*6+2
		ll2 = l*6+4
		print "{:9.6f} {:5d} : ".format(elsort[ll1],ll1)
		print "{:9.6f} {:5d} : ".format(elsort[ll2],ll2) ,
		printrelativeangle( v, 2, ll1, 4, ll1 ) 
		printrelativeangle( v, 2, ll1, 2, ll2 ) 
		printrelativeangle( v, 2, ll1, 4, ll2 ) 
		print ""
	print "\t sum_l V(mu,l)V(nu,l)^*"
	for l in range(nbath/6) :
		ll0 = l*6
		ll1 = l*6+2
		ll2 = l*6+4
		print "{:9.6f} {:5d} : ".format(elsort[ll0],ll0)
		for mu in range(2) :
			print "\t",
			for nu in range(2) :
				pp( lsum( v, mu, nu, [ll0,ll0+1] ) )
			print ""
		print "{:9.6f} {:5d} : ".format(elsort[ll1],ll1)
		print "{:9.6f} {:5d} : ".format(elsort[ll2],ll2) 
		for mu in range(4) :
			print "\t",
			for nu in range(4) :
				pp( lsum( v, mu+2, nu+2, [ll1,ll1+1,ll2,ll2+1] ) )
			print ""

def cphaseftn( v , elsort )  :
	print "\tr(xyu/xyd), r(yzu/yzd), r(yzu/zxu), r(yzu/zxd), r(yzd/zxd), r(zxu/zxd) | anglediff24-35"
	for l in range(nbath) :
		print "{:9.6f} {:5d} : ".format(elsort[l],l) ,
		printrelativeangle( v, 0, l, 1, l ) 
		printrelativeangle( v, 2, l, 3, l ) 
		printrelativeangle( v, 2, l, 4, l ) 
		printrelativeangle( v, 2, l, 5, l ) 
		printrelativeangle( v, 3, l, 5, l ) 
		printrelativeangle( v, 4, l, 5, l ) 
		print "|",
		print "{:10.4f}".format( relativeangle( v, 2, l, 4, l ) - relativeangle( v, 3, l, 5, l ) ),
		print ""
	print "\tr(xyu,l/xyu,l+1), (xyd,l/xyd,l+1), (xyu,l/xyd,l+1 *), (xyd,l/xyu,l+1 *) | crossconjdiff01, 23, 45"
	for l in range(nbath/2) :
		print "{:9.6f} {:5d} : ".format(elsort[2*l],2*l) ,
		printrelativeangle( v, 0, 2*l, 0, 2*l+1 ) 
		printrelativeangle( v, 1, 2*l, 1, 2*l+1 ) 
		printrelativeconjangle( v, 0, 2*l, 1, 2*l+1 ) 
		printrelativeconjangle( v, 1, 2*l, 0, 2*l+1 ) 
		print "|",
		printdiffconjangleTrev( v, 0, 2*l, 1, 2*l+1 ) 
		printdiffconjangleTrev( v, 2, 2*l, 3, 2*l+1 ) 
		printdiffconjangleTrev( v, 4, 2*l, 5, 2*l+1 ) 
		print ""
	print "\trat(xyu,l/xyu,l+1), (yzu,l+1/yzu,l), (zxu,l+1/zxu,l), (xyd,l/xyd,l+1), (yzd,l+1/yzd,l), (zxd,l+1/zxd,l), (yzu,l/yzu,l+1), (yzd,l/yzd,l+1) "
	for l in range(nbath/2) :
		print "{:9.6f} {:5d} : ".format(elsort[2*l],2*l) ,
		printratio( v, 0, 2*l, 0, 2*l+1 ) 
		printratio( v, 2, 2*l+1, 2, 2*l ) 
		printratio( v, 4, 2*l+1, 4, 2*l ) 
		printratio( v, 1, 2*l, 1, 2*l+1 ) 
		printratio( v, 3, 2*l+1, 3, 2*l ) 
		printratio( v, 5, 2*l+1, 5, 2*l ) 
		printratio( v, 2, 2*l, 2, 2*l+1 ) 
		printratio( v, 3, 2*l, 3, 2*l+1 ) 
		print ""
	print "\trat(xyu,xyd), (yzu,yzd), (zxu,zxd), (xyu/yzu), (yzu/zxu), (xyd/yzd), (yzd/zxd)"
	for l in range(nbath) :
		print "{:9.6f} {:5d} : ".format(elsort[l],l) ,
		printratio( v, 0, l, 1, l ) 
		printratio( v, 2, l, 3, l ) 
		printratio( v, 4, l, 5, l ) 
		printratio( v, 0, l, 2, l ) 
		printratio( v, 2, l, 4, l ) 
		printratio( v, 1, l, 3, l ) 
		printratio( v, 3, l, 5, l ) 
		print ""
	print "e1-e2"
	print "  xyu0/xyd0  xyu1/xyd1  xyu0/xyu1  xyd0/xyd1  xyu0/xyd1  xyd0/xyu1"
	aa = [
	printrelativeangle( v, 0, 0, 1, 0 ) ,
	printrelativeangle( v, 0, 1, 1, 1 ) ,
	printrelativeangle( v, 0, 0, 0, 1 ) ,
	printrelativeangle( v, 1, 0, 1, 1 ) ,
	printrelativeangle( v, 0, 0, 1, 1 ) ,
	printrelativeangle( v, 1, 0, 0, 1 ) ]
	print "|",
	[ pp( aa[i]-aa[i+1] ) for i in [0,2,4] ]
	[ pp( aa[i]+aa[i+1] ) for i in [0,2,4] ]
       	print ""
	aa = [
	printrelativeconjangle( v, 0, 0, 1, 0 ) ,
	printrelativeconjangle( v, 0, 1, 1, 1 ) ,
	printrelativeconjangle( v, 0, 0, 0, 1 ) ,
	printrelativeconjangle( v, 1, 0, 1, 1 ) ,
	printrelativeconjangle( v, 0, 0, 1, 1 ) ,
	printrelativeconjangle( v, 1, 0, 0, 1 ) ]
	print "|",
	[ pp( aa[i]-aa[i+1] ) for i in [0,2,4] ]
	[ pp( aa[i]+aa[i+1] ) for i in [0,2,4] ]
	print ""
	print "  yzu2/yzd2  yzu3/yzd3  yzu2/yzu3  yzd2/yzd3  yzu2/yzd3  yzd2/yzu3"
	aa = [
	printrelativeangle( v, 2, 2, 3, 2 ) ,
	printrelativeangle( v, 2, 3, 3, 3 ) ,
	printrelativeangle( v, 2, 2, 2, 3 ) ,
	printrelativeangle( v, 3, 2, 3, 3 ) ,
	printrelativeangle( v, 2, 2, 3, 3 ) ,
	printrelativeangle( v, 3, 2, 2, 3 ) ]
	print "|",
	[ pp( aa[i]-aa[i+1] ) for i in [0,2,4] ]
	[ pp( aa[i]+aa[i+1] ) for i in [0,2,4] ]
       	print ""
	aa = [
	printrelativeconjangle( v, 2, 2, 3, 2 ) ,
	printrelativeconjangle( v, 2, 3, 3, 3 ) ,
	printrelativeconjangle( v, 2, 2, 2, 3 ) ,
	printrelativeconjangle( v, 3, 2, 3, 3 ) ,
	printrelativeconjangle( v, 2, 2, 3, 3 ) ,
	printrelativeconjangle( v, 3, 2, 2, 3 ) ]
	print "|",
	[ pp( aa[i]-aa[i+1] ) for i in [0,2,4] ]
	[ pp( aa[i]+aa[i+1] ) for i in [0,2,4] ]
	print ""
	print "  zxu4/zxd4  zxu5/zxd5  zxu4/zxu5  zxd4/zxd5  zxu4/zxd5  zxd4/zxu5"
	aa = [
	printrelativeangle( v, 4, 4, 5, 4 ) ,
	printrelativeangle( v, 4, 5, 5, 5 ) ,
	printrelativeangle( v, 4, 4, 4, 5 ) ,
	printrelativeangle( v, 5, 4, 5, 5 ) ,
	printrelativeangle( v, 4, 4, 5, 5 ) ,
	printrelativeangle( v, 5, 4, 4, 5 ) ]
	print "|",
	[ pp( aa[i]-aa[i+1] ) for i in [0,2,4] ]
	[ pp( aa[i]+aa[i+1] ) for i in [0,2,4] ]
       	print ""
	aa = [
	printrelativeconjangle( v, 4, 4, 5, 4 ) ,
	printrelativeconjangle( v, 4, 5, 5, 5 ) ,
	printrelativeconjangle( v, 4, 4, 4, 5 ) ,
	printrelativeconjangle( v, 5, 4, 5, 5 ) ,
	printrelativeconjangle( v, 4, 4, 5, 5 ) ,
	printrelativeconjangle( v, 5, 4, 4, 5 ) ]
	print "|",
	[ pp( aa[i]-aa[i+1] ) for i in [0,2,4] ]
	[ pp( aa[i]+aa[i+1] ) for i in [0,2,4] ]
	print ""

def cphaseftnjeff( v , elsort ) :
#	print "T(qq+/qq-), T(qd+/qd-), T(d+/d-)"
#	for l in range(nbath/2) :
#		print "{:9.6f} {:5d} : ".format(elsort[2*l],2*l) ,
#		printdiffconjangleTrev( v, 0, 2*l, 3, 2*l+1 ) 
#		printdiffconjangleTrev( v, 1, 2*l, 2, 2*l+1 ) 
#		printdiffconjangleTrev( v, 4, 2*l, 5, 2*l+1 ) 
#		print ""
#	print "\trat(qq+,l/qq+,l+1), (qd+,l+1/qd+,l), (dd+,l+1/dd+,l), (qq-,l/qq-,l+1), (qd-,l+1/qd-,l), (dd-,l+1/dd-,l)"
#	for l in range(nbath/2) :
#		print "{:9.6f} {:5d} : ".format(elsort[2*l],2*l) ,
#		printratio( v, 0, 2*l, 0, 2*l+1 ) 
#		printratio( v, 1, 2*l, 1, 2*l+1 ) 
#		printratio( v, 4, 2*l, 4, 2*l+1 ) 
#		printratio( v, 3, 2*l+1, 3, 2*l ) 
#		printratio( v, 2, 2*l+1, 2, 2*l ) 
#		printratio( v, 5, 2*l+1, 5, 2*l ) 
#		print ""
	print "\t sum_l V(mu,l)V(nu,l)^* (mu,nu of {alpha=q3/2,gamma=q1/2,beta=d1/2})"
	for l in range(nbath/6) :
		ll0 = l*6
		ll1 = l*6+2
		ll2 = l*6+4
		print "{:9.6f} {:5d} : ".format(elsort[ll0],ll0)
		indalphagammabeta = [0,3,1,2,4,5] 
		for mu in indalphagammabeta :
			print "\t",
			for nu in indalphagammabeta :
				pp( lsum( v, mu, nu, [ll0,ll0+1] ) )
			print ""
		print "{:9.6f} {:5d} : ".format(elsort[ll1],ll1)
		print "{:9.6f} {:5d} : ".format(elsort[ll2],ll2) 
		for mu in indalphagammabeta :
			print "\t",
			for nu in indalphagammabeta :
				pp( lsum( v, mu, nu, [ll1,ll1+1,ll2,ll2+1] ) )
			print ""

if (args.cphase and args.tot2g) :
	print type(mclass.npvtldat) , np.shape(mclass.npvtldat)
	nbath = np.shape(mclass.npvtldat)[1]
	if args.unsort :
		v = mclass.npvtldatusort
		dumel = el 
	else : 
		v = mclass.npvtldat
		dumel = elsort
	cphaseftn( v, dumel ) 
elif args.cphasenosoc and args.tot2g :
	print type(mclass.npvtldat) , np.shape(mclass.npvtldat)
	nbath = np.shape(mclass.npvtldat)[1]
	if args.unsort :
		v = mclass.npvtldatusort
		dumel = el 
	else : 
		v = mclass.npvtldat
		dumel = elsort
	cphaseftnNoSOC( v, dumel ) 
elif (args.cphasenosoc and args.toleff) :
	print type(mclass.npvfromjeff) , np.shape(mclass.npvfromjeff)
	nbath = np.shape(mclass.npvfromjeff)[1]
	v = mclass.npvfromjeff
	dumel = el
	cphaseftnNoSOC( v, dumel ) 
elif args.cphaset2g : 
	print type(mclass.npvmlsort) , np.shape(mclass.npvmlsort)
	nbath = np.shape(mclass.npvmlsort)[1]
	v = mclass.npvmlsort
	cphaseftn( v, elsort ) 
elif args.cphase :
	print type(mclass.npvmlsort) , np.shape(mclass.npvmlsort)
	nbath = np.shape(mclass.npvmlsort)[1]
	v = mclass.npvmlsort
	if args.unsort :
		v = mclass.npvmldat
		dumel = el 
	else : 
		v = mclass.npvmlsort
		dumel = elsort
	cphaseftnjeff( v, dumel ) 
elif args.sphase and args.tot2g :
	print type(mclass.npvtldat) , np.shape(mclass.npvtldat)
	nbath = np.shape(mclass.npvtldat)[1]
	if args.unsort :
		v = mclass.npvtldatusort
		dumel = el 
	else : 
		v = mclass.npvtldat
		dumel = elsort
	cphaseftnselectivet2g( v, dumel ) 
elif args.sphase :
	print type(mclass.npvmlsort) , np.shape(mclass.npvmlsort)
	nbath = np.shape(mclass.npvmlsort)[1]
	if args.unsort :
		v = mclass.npvmldat
		dumel = el
	else : 
		v = mclass.npvmlsort
		dumel = elsort
	cphaseftnselectivejeff( v, dumel ) 


import random
def randscale( a, b ) :
	return random.random() *(b-a) + a
def writenewmak( v, elarr, fhead='symmak', nb=9 , nc=3 , ncfill=2 , elrange=[-3,3] , fname=None ) :
	if fname :  pass
	else : 
		fname = fhead+'_nb{}_nc{}_ncf{}_l{}_r{}.mak'.format( nb, nc, ncfill , elrange[0], elrange[1] ) 
	with open( fname , 'w' ) as fw :
		fw.write( "0\t0\t0\t0\n" )
		fw.write( "128\t512\n" )
		fw.write( "1\n" )
		fw.write( "-99\t7036530\t2704156\n" )
		for i in range(nb*2) :
			fw.write( "{:19.16f}\t\t\t\t".format( elarr[i] )  )
		fw.write('\n\n')
		for i in range(nc*2) :
				for j in range(nb*2) :
					fw.write( "{:20.17f}  ".format( v[i][j].real )  )
					fw.write( "{:20.16f}\t".format( v[i][j].imag )  )
				fw.write( '\n' ) 
		fw.write( '\n\n' ) 
if args.makesymnew : 
	nbasis = 6 
	nbath  = 18 
	bathperbasis  = nbath/nbasis
	if args.tot2g : 
		fill = 2
		structarr = [ [0,1],[2,3,4,5] ]
		#
		tsymarr = [	[0,0,1,1, np.exp( 1j*np.pi)	] ,
				[1,0,0,1, 1			] ,
				[2,2,3,3, np.exp( 1j*np.pi)	] ,
				[3,2,2,3, 1			] ,
				[2,4,3,5, np.exp( 1j*np.pi)	] ,
				[3,4,2,5, 1			] ,
				[4,2,5,3, np.exp( 1j*np.pi)	] ,
				[5,2,4,3, 1			] ,
				[4,4,5,5, np.exp( 1j*np.pi)	] ,
				[5,4,4,5, 1			] ]
		#
		c4symarr = [	[2,2,4,2, np.exp(-1j*np.pi/2)	] ,
				[3,2,5,2, np.exp( 1j*np.pi/2)	] ,
				[2,3,4,3, np.exp(-1j*np.pi/2)	] ,
				[3,3,5,3, np.exp( 1j*np.pi/2)	] ,
				[2,4,4,4, np.exp(-1j*np.pi/2)	] ,
				[3,4,5,4, np.exp( 1j*np.pi/2)	] ,
				[2,5,4,5, np.exp(-1j*np.pi/2)	] ,
				[3,5,5,5, np.exp( 1j*np.pi/2)	] ]
	else :
		structarr = [ 2, 4]
		tsymarr = [ [0,3], [1,2,4,5] ]
		print "Available only in t2g. Exit(1)."
		sys.exit(1)

	elrange = [-3,3]
	elarr = []
	vnl  = np.array( [ [ 0 for i in range(nbath) ] for j in range(nbasis) ] , dtype=complex ) 
	for ii in range(bathperbasis) :
		for jj in structarr   :
			lenjj = len(jj)
			if ii<fill :  
				elrand  = randscale( elrange[0], 0 )
			else :
				elrand  = randscale( 0, elrange[1] )
			for k in range(lenjj) :
				elarr.append( elrand )
			for mu in jj :
				for nu in jj :
					ll = nbasis*ii + nu
					vnl[mu][ll] = randscale( -1., 1. ) + 1j* randscale( -1., 1. )
	print "elarr : ", elarr  
	print "shape : ", np.shape(elarr)
	print "shape : ", np.shape(vnl)
	print "Setting randomly in the ", structarr, " structure."
	
	for ii in range(bathperbasis) :
		mu=0
		for jj in tsymarr  :
			mu1, l1, mu2, l2  = [ jj[0], nbasis*ii+jj[1], jj[2], nbasis*ii+jj[3] ]
			#print "l1l2 jj : ", l1, l2 , jj[:4], "{:9.5f}\t {}\t{}".format(jj[4] ,  vnl[mu2][l2] , vnl[mu1][l1].conjugate() * jj[4] )
			vnl[mu2][l2]  = vnl[mu1][l1].conjugate() * jj[4]
			#print mu1, mu2, l1, l2
			mu += 1 
	for ii in range(bathperbasis) :
		for jj in c4symarr  :
			mu1, l1, mu2, l2  = [ jj[0], nbasis*ii+jj[1], jj[2], nbasis*ii+jj[3] ]
			vnl[mu2][l2]  = vnl[mu1][l1] * jj[4]
			#print "l1l2 jj : ", l1, l2 , jj[:4], "{:9.5f}\t {}\t{}".format(jj[4] ,  vnl[mu2][l2] , vnl[mu1][l1] * jj[4] )

	print "vnl : "
	printAbsVmlel( vnl , elarr ) 
	vjl = returnVjlfromVtl( vnl ) 
	print "vjl : "
	printAbsVmlel( vjl , elarr ) 
	writenewmak( vjl, elarr , fhead="random" ) 
	
	
def returnPartnerel ( el, tnb ) : 
	partnerel = [ [] for a in range(tnb) ]
	for l1 in range( tnb ) :
		dumel = el[l1] 
		for l2 in range( tnb ) :
			if abs(dumel-el[l2])<5e-4 :
				partnerel[l1].append( l2 ) 
	#print "partnerel : "
	#print partnerel

	partnerelstruct = []
	dum = ""
	for a in partnerel : 
		if dum==a 	: pass
		else 		: partnerelstruct.append( a ) ; dum=a 
	#print "partnerelstruct : "
	#print partnerelstruct
	return partnerel, partnerelstruct

def averageel ( el, tnb ) : 
	partnerel, partnerelstruct = returnPartnerel( el, tnb ) 
	elavg = np.zeros( len(el) )
	for elindarr in partnerelstruct :
		dumavg = np.sum( el[elindarr] ) / len(elindarr)
		elavg[elindarr] = np.array( [ dumavg for aa in elindarr ] ) 
	return elavg
			

if args.makesym : 
	if args.tot2g :
		structarr = [ [0,1],[2,3,4,5] ]
		c4symarr = [	[2,2,4,2, np.exp(-1j*np.pi/2)	] ,
				[3,2,5,2, np.exp( 1j*np.pi/2)	] ,
				[2,3,4,3, np.exp(-1j*np.pi/2)	] ,
				[3,3,5,3, np.exp( 1j*np.pi/2)	] ,
				[2,4,4,4, np.exp(-1j*np.pi/2)	] ,
				[3,4,5,4, np.exp( 1j*np.pi/2)	] ,
				[2,5,4,5, np.exp(-1j*np.pi/2)	] ,
				[3,5,5,5, np.exp( 1j*np.pi/2)	] ]
	#find activated basis in block 
	nbath = np.shape(mclass.npvmlsort)[1]
	bathperbasis  = nbath/nbasis
	v = mclass.npvmlsort
	vnlsym  = np.array( [ [ 0 for i in range(nbath) ] for j in range(nbasis) ] , dtype=complex ) 
	#find el-parner
	partnerel, partnerelstruct = returnPartnerel( elsort , nbath )

	orbitalinel = [ [] for a in range(nbath) ]
	for elindarr in partnerelstruct :
		for ll  in elindarr :
			for mu in range(nbasis) :
				if np.abs( v[mu][ll] ) < 2e-3 : pass
				else : vnlsym[mu][ll] = v[mu][ll] ; orbitalinel[ll].append(mu) 
	print "orbitalinel : "
	print orbitalinel
	print "vnlsym : "
	printAbsVmlel( vnlsym , elsort ) 
	vtlsym = returnVtlfromVjl( vnlsym ) 
	print "vtlsym : "
	printAbsVmlel( vtlsym , elsort ) 
	print "elsort : "
	print elsort
	elavg = averageel( elsort , nbath ) 
	print "elavg : "
	print elavg
	writenewmak( vnlsym, elavg , fhead="cleaned" ) 

import subprocess
if args.makephase : 
	vl = mclass.npvmlsort
	nbath = np.shape( vl )[-1]
	vlphase  = np.array( [ [ 0 for i in range(nbath) ] for j in range(nbasis) ] , dtype=complex ) 
	for ll  in range(len(elsort)) :
		for mu in range(nbasis) :
			if mu==0 or mu==4 or mu==5 :
				vlphase[mu][ll] = vl[mu][ll] * -1.
			else : 
				vlphase[mu][ll] = vl[mu][ll] 
	thisfoldertag = subprocess.check_output( 'cat this.tag', shell=True ).split()[0]
	copymak = hobj.maknameonly( thisfoldertag  )
	writenewmak( vlphase, elsort , fname="./makchphase/"+copymak ) 

if args.makeconverge : 
	vl = mclass.npvmlsort
	nbath = np.shape( vl )[-1]
	thisfoldertag = subprocess.check_output( 'cat this.tag', shell=True ).split()[0]
	copymak = hobj.maknameonly( thisfoldertag  )
	writenewmak( vl, elsort , fname="./makconv/"+copymak ) 
	


sys.exit(1)
