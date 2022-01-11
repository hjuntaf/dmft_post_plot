import numpy as np
import sys, os
import re
import subprocess

nbasis = 6
Tara =  np.matrix( [
                [ 0.,                   0.,             np.sqrt(1./2.), 0.,                     np.sqrt(1./2.)*1j ,     0.                      ],  
                [ np.sqrt(2./3.),       0.,             0.,             -1./np.sqrt(6.),        0.,                     -1/np.sqrt(6.)*1j       ],  
                [ 0.,                   np.sqrt(2./3.), 1./np.sqrt(6.), 0.,                     -1./np.sqrt(6.)*1j,     0.                      ],  
                [ 0.,                   0.,             0.,             1./np.sqrt(2.),         0.,                     -1./np.sqrt(2.)*1j      ],  
                [ 1./np.sqrt(3.),       0.,             0.,             1./np.sqrt(3.),         0.,                     1./np.sqrt(3.)*1j       ],  
                [ 0.,                   -1./np.sqrt(3.),        1./np.sqrt(3.), 0.,                     -1./np.sqrt(3.)*1j,     0.              ]] )
a = 0.980031
b = 0.198847
matTjnat = np.matrix( [
		[ 1, 0, 0, 0, 0, 0 ],
		[ 0, 0, 0, 1, 0, 0 ],
		[ 0,-a, 0, 0,-b, 0 ],
		[ 0, 0, a, 0, 0,-b ],
		[ 0,-b, 0, 0, a, 0 ],
		[ 0, 0, b, 0, 0, a ]
] , dtype=complex)

Taraph = np.matrix( [	# = Tt2gjeff
	[ 0.,			0.,		-np.sqrt(1./2.),	0.,			-np.sqrt(1./2.)*1j ,	0.			], 
	[ np.sqrt(2./3.),	0.,		0.,		-1./np.sqrt(6.),	0.,			-1/np.sqrt(6.)*1j	],
	[ 0.,			np.sqrt(2./3.),	1./np.sqrt(6.),	0.,			-1./np.sqrt(6.)*1j,	0.			],
	[ 0.,			0.,		0.,		1./np.sqrt(2.),		0.,			-1./np.sqrt(2.)*1j	],
	[-1./np.sqrt(3.),	0.,		0.,		-1./np.sqrt(3.),		0.,		-1./np.sqrt(3.)*1j	],
	[ 0.,			1./np.sqrt(3.),	-1./np.sqrt(3.),	0.,			1./np.sqrt(3.)*1j,	0.		]] ) 
Tlefft2g = np.matrix( [
	[ np.sqrt(2),			0.,		0.,		0.,			0.,			0.			], 
	[ 0.,			np.sqrt(2),		0.,		0.,			0.,			0.			], 
	[ 0.,				0.,		-1.,		0.,			-1j,			0.			], 
	[ 0.,				0.,		0.,		-1.,			0.,			-1j			], 
	[ 0.,				0.,		1.,		0.,			-1j,			0.			], 
	[ 0.,				0.,		0.,		1.,			0.,			-1j			]] ) / np.sqrt(2.)

lst2g = np.matrix( [
	[ 0.,				0.,		0.,		.5,			0.,			-1j*.5			], 
	[ 0.,				0.,		-.5,		0.,			-1j*.5,			0.			], 
	[ 0.,				-.5,		0.,		0.,			 1j*.5,			0.			], 
	[ .5,				0.,		0.,		0.,			0.,			-1j*.5			], 
	[ 0.,				1j*.5,		-1j*.5,		0.,			0.,			0.			], 
	[ 1j*.5,			0.,		0.,		1j*.5,			0.,			0.			]] ) 
Tleffjeff = np.dot( Taraph , np.conjugate(Tlefft2g).transpose() ) 
Tjeffleff = np.conjugate(Tleffjeff).transpose()
lsjeff = np.diag([ -0.5, -0.5, -0.5, -0.5, 1., 1. ])
	
Trevt2g = np.matrix( [
                [ 0,-1, 0, 0, 0, 0 ],
                [ 1, 0, 0, 0, 0, 0 ],
                [ 0, 0, 0,-1, 0, 0 ],
                [ 0, 0, 1, 0, 0, 0 ],
                [ 0, 0, 0, 0, 0,-1 ],
                [ 0, 0, 0, 0, 1, 0 ]
], dtype=complex )
Trevjeff = np.matrix( [
                [ 0, 0, 0, 1, 0, 0 ],
                [ 0, 0,-1, 0, 0, 0 ],
                [ 0, 1, 0, 0, 0, 0 ],
                [-1, 0, 0, 0, 0, 0 ],
                [ 0, 0, 0, 0, 0, 1 ],
                [ 0, 0, 0, 0,-1, 0 ]
], dtype=complex )

def filtPrintformat( item , val)  :
	item4	= ( (item.find('gptl')>-1) or (item.find('degen')>-1) or (item.find('npos')>-1) )
	item5	= ( (item=='U') or (item=='J') or (item=='S') or (item=='D') or (item=='JT') or (item=='JTD') )
	item6	= ( (item.find('iter')>-1) or (item=='beta') or (item=='Nmax') or (item=='IT') )
	if item4 :
		out	= "{:<4} ".format(val) 
	elif item5 :
		out	= "{:<5} ".format(val) 
	elif item6 :
		out	= "{:<6} ".format(val) 
	else : 
		out	= "{:<17} ".format(val) 
	return out

def cast_type( a ) :
	try :
		return int(a)
	except ValueError :
		pass
	try :
		return float(a)
	except ValueError :
		return a
class dstruct  :
	def __init__(self , ddir ) :
		self.d = ddir 
		self.dic = get_dic_from_self()

	def get_dic_from_str(self, string):
		dic = {}
		pattern = re.compile(r"[-]?\d*\.\d+|\d+")
		sp = re.split(r'[/_]', string)
		for x in sp:
			numbers = [ cast_type(y) for y in re.findall(pattern, x)]
			strings = filter(None, re.split(pattern, x))
			dic.update(dict(zip(strings, numbers)))
	
		#print dic
		return dic
	def get_dic_from_self(self) :
		string = self.d
		dic = {}
		pattern = re.compile(r"[-]?\d*\.\d+|\d+")
		sp = re.split(r'[/_]', string)
		for x in sp:
			numbers = [ cast_type(y) for y in re.findall(pattern, x)]
			strings = filter(None, re.split(pattern, x))
			dic.update(dict(zip(strings, numbers)))
	
		#print dic
		return dic

#	def rmak(self ) :
#		print "\n#initial_mak"
#		fargs[j] = args.dirs[j].split('_')
#		for farg in fargs[j] :
#			if farg.find("UF")>-1 :
#				UF = float( farg[farg.find("UF")+2:] )
#		fw = open( args.dirs[j] + "/output_u{:.2f}.txt".format( UF )  ) 
#		fpars[j] = fw.readlines()
#		if len( fpars[j] ) > 0 :
#			print j, "\t", fpars[j][0].split()[-1]
#		else :
#			print j, "\tNone"
#		fw.close()

class cdstruct :  
	def __init__( self, dirs ) :
		dsts = [ dstruct(string) for string in  dirs ] 
		print "dsts : ", dsts

class datclass :
	def __init__( self, dirs ) :
		self.d = dirs
		self.tdir = self.d
	def get_params( self ) :
		for t in self.d.split("_") :
			if( t.find("UF") > -1 ) : 
				UF = float(t[2:])
			if( t.find("J") > -1 ) : 
				if( t.find("JT") > -1 ) : 
					JT = float(t[2:])
				else :
					J = float(t[1:])
			if( t.find("S") > -1 ) : 
				S = float(t[1:])
			if( t.find("D") > -1 ) : 
				if( t.find("Dir") > -1 ) :  pass
				else :
					D = float(t[1:])
			if( t.find("nb") > -1 ) : 
				Nb = int( t[ t.find("nb")+2: ] )

def mToD( m ) :  #month to day
	if m==1 or m==3 or m==5 or m==7 or m==8 or m==10 or m==12 :
		return 31 
	elif m==4 or m==6 or m==9 or m==11 :
		return 30
	elif m==2 :
		return 28

daydic = { "Jan" : "1" , "Feb" : "2" , "Mar" : "3" ,
	"Apr" : "4" , "May" : "5" , "Jun" : "6" ,
	"Jul" : "7" , "Aug" : "8" , "Sep" : "9" ,
	"Oct" : "10" , "Nov" : "11" , "Dec" : "12" ,
}

def transfMjefftot2g( origMat, transfMat ) :
	#print "transfMjefftot2g done."
	resultMat = [ range(nbasis) for j in range(nbasis) ]
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			resultMat[mu][nu] = origMat[mu][nu]*0.
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					resultMat[mu][nu] += transfMat[j1,mu] * origMat[j1][j2] * np.conj( transfMat[j2,nu] )
	return resultMat 

def eigMat( origMat ) :
	res = np.reshape( origMat, (nbasis*nbasis,-1) )
	res = res.transpose()
	for ii in range( len(res) )  :
		mati =  np.matrix( np.reshape(res[ii], (nbasis,nbasis) ) )
		matisym = mati + np.conjugate(mati).transpose()
		matisym = matisym / 2.
		w,v = np.linalg.eig( matisym )
		#print "Eig Shape : ", np.shape( a[ii] ) 
		#print "Eig w : ", w
		#print "Eig w.shape : ", np.shape(w)
		#print "Eig Shape : ", np.diag(w)
		#w = np.sort(w)
		#print "w0 : "
		#print w 
		#print "v0 : "
		#for j1 in range(nbasis) :
		#	for j2 in range(nbasis) :
		#		print "{:30.8f}\t".format( v[j1,j2] ),
		#	print ""
		w = np.dot( v.transpose() , np.dot( mati, v.conjugate() ) )
		idx = np.diag(-w.imag).argsort()[::-1]
		v = v[:,idx]
		w = np.dot( v.transpose() , np.dot( mati, v.conjugate() ) )
		#print "w : ", type(w), np.shape(w)
		#for j1 in range(nbasis) :
		#	for j2 in range(nbasis) :
		#		print "{:30.8f}\t".format( w[j1,j2] ),
		#	print ""
		#sys.exit(1)
		res[ii] = w.flatten() 
		#sys.exit(1)
	#print "Eig type : ", type( a ) 
	#print "Eig Shape : ", np.shape( a ) 
	res = np.reshape( res.transpose() , (nbasis,nbasis,-1) )
	#print "Eig Shape : ", np.shape( a ) 
	return res

def addSOCMatt2g( origMat, S=False, socMat=False ) :
	if socMat is False :	
		socMat = lst2g
	if S      is False :	S = 1.
	socMat = S*socMat
	resultMat = [ range(nbasis) for j in range(nbasis) ]
	for j1 in range(nbasis) :
		for j2 in range(nbasis) :
			resultMat[j1][j2] = origMat[j1][j2]*0.
	for j1 in range(nbasis) :
		for j2 in range(nbasis) :
			resultMat[j1][j2] = origMat[j1][j2] + socMat[j1,j2]
	return resultMat 
def addSOCMat( origMat, basis="t2g", S=False, socMat=False ) :
	if socMat is False :	
		if basis.find("t")>-1 : socMat = lst2g
		if basis.find("j")>-1 : socMat = lsjeff
		else :			print "ERROR in 'addSOCMat()'."
	if S      is False :	S = 1.
	socMat = S*socMat
	resultMat = [ range(nbasis) for j in range(nbasis) ]
	for j1 in range(nbasis) :
		for j2 in range(nbasis) :
			resultMat[j1][j2] = origMat[j1][j2]*0.
	for j1 in range(nbasis) :
		for j2 in range(nbasis) :
			resultMat[j1][j2] = origMat[j1][j2] + socMat[j1,j2]
	return resultMat 

def transfMtrconj( origMat, transfMat ) :
	#print "transfMtrconj done."
	resultMat = [ range(nbasis) for j in range(nbasis) ]
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			resultMat[mu][nu] = origMat[mu][nu]*0.
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					resultMat[mu][nu] += transfMat[j1,mu] * origMat[j1][j2] * np.conj( transfMat[j2,nu] )
	return resultMat 

def transfMconjtr( origMat, transfMat ) :
	#print "transfMconjtr done."
	resultMat = [ range(nbasis) for j in range(nbasis) ]
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			resultMat[mu][nu] = origMat[mu][nu]*0.
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					resultMat[j1][j2] += np.conj( transfMat[j1,mu] )  * origMat[mu][nu] * transfMat[j2,nu]
	return resultMat 
def transfMwnconjtr( origMat, transfMat ) :
	resultMat = [ range(nbasis) for j in range(nbasis) ]
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			resultMat[mu][nu] = origMat[mu][nu]*0.
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					resultMat[j1][j2] += np.conj( transfMat[j1,mu] )  * origMat[mu][nu] * transfMat[j2,nu]
	return resultMat 

def transfMt2gtojeff( origMat, transfMat ) :
	#print "transfMt2gtojeff done."
	resultMat = [ range(nbasis) for j in range(nbasis) ]
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			resultMat[mu][nu] = origMat[mu][nu]*0.
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					resultMat[j1][j2] += np.conj( transfMat[j1,mu] )  * origMat[mu][nu] * transfMat[j2,nu]
	return resultMat 

def transfM_Trev( origMat, transfMatTrev ) :
	#print "transfMconjtr done."
	resultMat = [ range(nbasis) for j in range(nbasis) ]
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			resultMat[mu][nu] = origMat[mu][nu]*0.
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					resultMat[mu][nu] += transfMatTrev[j1,nu] * -origMat[j1][j2] * transfMatTrev[mu,j2]
	return resultMat 

def transfM_non_dagger( origMat, transfMat ) :
	#print "transfMconjtr done."
	resultMat = [ range(nbasis) for j in range(nbasis) ]
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			resultMat[mu][nu] = origMat[mu][nu]*0.
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					resultMat[j1][j2] += transfMat[j1,mu] * origMat[mu][nu] * np.conj( transfMat[j2,nu] )
	return resultMat 

def mul_sparseAB_complex( sA, sB , dim )  :		# sT.conjugateTranspose * M * sT
	dumMat	= np.zeros( dim*dim , dtype=complex ).reshape( dim, dim )
	for ist0 in sA : 
		mu0, nu0, dat0	= ist0
		for ist1 in sB : 
			mu1, nu1, dat1	= ist1
			# [ nu0, mu0 ] . [mu0, mu1] . [mu1, nu1]
			if nu0 == mu1 : 
				dumMat[mu0][nu1]	+=  dat0 * dat1
	return dumMat

def transf_ArrMat_sparseTwf_complex( arrMat , sT )  :		# sT.conjugateTranspose * M * sT
	dim	= np.shape(arrMat)[0]
	dumMat	= np.zeros( dim*dim , dtype=complex ).reshape( dim, dim )
	for ist0 in sT : 
		print "ist0 : ",  ist0
		mu0, nu0, dat0	= ist0
		for ist1 in sT : 
			mu1, nu1, dat1	= ist1
			# [ nu0, mu0 ] . [mu0, mu1] . [mu1, nu1]
			dumMat[nu0][nu1]	+= np.conjugate( dat0 ) * arrMat[mu0][mu1] * dat1
	return dumMat

def printArrMatAuto( arrMat )  :
	dim	= np.shape(arrMat)[0]
	for i in range(dim) :
		for j in range(dim) :
			print "{:17.5f} ".format(arrMat[i][j]),
		print ""

def printNpMat( npMat )  :
	for i in range(nbasis) :
		for j in range(nbasis) :
			print "{:17.5f} ".format(npMat[i,j]),
		print ""

def printArrMat( arrMat )  :
	for i in range(nbasis) :
		for j in range(nbasis) :
			print "{:17.5f} ".format(arrMat[i][j]),
		print ""

def printArrMatdim( arrMat , dim)  :
	for i in range(dim) :
		for j in range(dim) :
			print "{:17.5f} ".format(arrMat[i][j]),
		print ""

def printAbsArrMatdim( arrMat , dim)  :
	for i in range(dim) :
		for j in range(dim) :
			print "{:7.5f} ".format( np.abs(arrMat[i][j]) ),
		print ""

def transfconjTOtransT( origMat, transfMat ) :
	print "transfconjTOtransT done."
	resultMat = np.zeros( (nbasis, nbasis) , dtype=complex)
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			resultMat[mu][nu] = origMat[mu][nu]*0.
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					resultMat[j1][j2] += np.conj( transfMat[j1,mu] )  * origMat[mu][nu] * transfMat[j2,nu]
	return resultMat 
def transftransTOconjT( origMat, transfMat ) :
	print "transftransTOconjT done." 
	resultMat = np.zeros( (nbasis, nbasis) , dtype=complex)
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			resultMat[mu][nu] = origMat[mu][nu]*0.
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					resultMat[mu][nu] += transfMat[j1,mu] * origMat[j1][j2] * np.conj( transfMat[j2,nu] )
	return resultMat 


def printCompVmlelnp( vl, el )  :
	for j in range(np.shape(vl)[1]) :
		for i in range(np.shape(vl)[0]) :
			print "{:19.6f} ".format( vl[i,j] ),
		print "| {:9.6f}".format( el[j] )

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

def sortVmlel( npvml, npel )  : 
	sort_perm = npel.argsort()
	npvmlsort    = npvml[:, sort_perm ]
	npelsort     =  npel[sort_perm ]
	return npvmlsort, npelsort

def returnVjlfromVtl( vtl, tnb=18 ) :
	vjl = np.zeros( (nbasis, tnb ) , dtype=complex)
	for mu in range( nbasis ) :
        	for aa in range( tnb ) :
			for nu in range( nbasis ) :
				vjl[mu][aa] += vtl[nu][aa] * Tara[mu,nu]
	return vjl 
def returnVtlfromVjl( vjl, tnb=18 ) :
	vtl = np.zeros( (nbasis,tnb) , dtype=complex)
	for mu in range( nbasis ) :
        	for aa in range(tnb) :
			for nu in range( nbasis ) :
				vtl[mu][aa] += vjl[nu][aa] * np.conj( Tara[nu,mu] )
	return vtl

def c4rotVml( vmlArrMat , basis )  : 
	phaseh		= np.exp( -1j * 0.25 * np.pi )	#spinhalfphase
	phasehinv	= np.exp(  1j * 0.25 * np.pi )	#spinhalfphase
	phaseq		= np.exp( -1j * 0.75 * np.pi )	#spinquadphase	
	phaseqinv	= np.exp(  1j * 0.75 * np.pi )	#spinquadphase	
	if basis.find("t")>-1 : 
		c4Mat	= np.matrix( [
			       	[ phaseh	, 0 			, 0, 0, 0, 0 ] , 
				[ 0		, phasehinv		, 0, 0, 0, 0 ] , 
				[ 0, 0, 0,	0, 		phaseh,	0 ] , 
				[ 0, 0, 0,	0,		0,	phasehinv ] , 
				[ 0, 0, -phaseh,0, 		0, 		0 ] , 
				[ 0, 0, 0,	-phasehinv,	0, 		0 ]
				] , dtype=complex 
				)
	elif basis.find("j")>-1 : 
		c4Mat	= np.matrix( 
				np.diag(
					[ phaseq, phaseh, phasehinv, phaseqinv, phaseh, phasehinv]
				       )
				)
	rotvmlMat	= np.dot( c4Mat, np.matrix(vmlArrMat)  ) 
	#print "printNpMat( rotvmlMat )  : " 
	#printCompVmlel( rotvmlMat ) 
	return rotvmlMat

def returnEclusterarr( path , basis)  : 
	if basis.find("t")>-1 : 
		fname =  path + "/Ecluster_before_rotation.dat" 
	elif basis.find("j")>-1 : 
		fname =  path + "/Ecluster.dat" 
	rawEcluster = np.genfromtxt( fname )
	Ecluster = np.zeros( (nbasis, nbasis) , dtype=complex)
	for arrdat in rawEcluster :
		Ecluster[ int(arrdat[0]) ][ int(arrdat[1]) ] = arrdat[2] + 1j*arrdat[3] 
	return Ecluster

def pA( M ) :
        for i in range(len(M)) :
                print M[i]
def pAF( M ) :
        for i in range(len(M)) :
        	for j in range(len(M[i])) :
                	print "{:%9.6f}".format(M[i][j]),
		print ""
def pAr( M ) :
        for i in range(len(M)) :
        	for j in range(len(M[i])) :
                	print "{}".format(M[i][j]),
		print ""
def pMF( M ) :
        for i in range(len(M)) :
        	for j in range(len(M[i])) :
                	print "{:%9.6f}".format(M[i,j]),
		print ""
def pMr( M ) :
        for i in range(len(M)) :
        	for j in range(len(M[i])) :
                	print "{}".format(M[i,j]),
		print ""
def Rspinhalf( n, phi ) :
        R       = np.zeros(4,dtype=complex).reshape(2,2)
        R[0][0] = np.cos(phi/2.) - 1j*n[2]*np.sin(phi/2.)
        R[1][1] = np.cos(phi/2.) + 1j*n[2]*np.sin(phi/2.)
        R[0][1] = (-1j*n[0]-n[1])*np.sin(phi/2.)
        R[1][0] = (-1j*n[0]+n[1])*np.sin(phi/2.)
        return R

class dateclass : 
	def __init__( self ) :
		self.date = subprocess.check_output( 'date +"%x %T"', shell=True )
		[self.day , self.time ] = self.date.split()
	def readdate( self , dateout ) :
		self.date = dateout
		[self.day , self.time ] = self.date.split()
	def call( self ) :
		print "Date.call() : ", self.date,
		print "Date.day,time : ", self.day, self.time
	def show( self ) :
		print self.day, self.time
	def elapsedd ( self , dat ) :
		d1 = np.array( self.day.split('/') ,dtype=int ) #MM/DD/YYYY
       		d2 = np.array( dat.split('/')      ,dtype=int ) #MM/DD/YYYY
		if d1[2]< 2000 : d1[2] += 2000
		el = d1 - d2 
		if el[2]!=0 :
			return el[1] + mToD( d2[0] ) 
		else :
			return el[1] + el[0] * mToD( d2[0] ) 
	def elapsedt ( self , dat, hday ) :
		d1 = np.array( self.time.split(':'),dtype=int )
       		d2 = np.array( dat.split(':')      ,dtype=int ) 
		el = d1 - d2 
		if el[0]<0 : el[0]+=24 ; hday -= 1
		if el[1]<0 : el[1]+=60 ; el[0]-= 1
		if el[2]<0 : el[2]+=60 ; el[1]-= 1
		el[0] += hday * 24
		return "{}:{}:{}".format( el[0], el[1], el[2] ) 
	def diffdate ( self , dclass ) :
		return self.elapsedt( dclass.time, self.elapsedd( dclass.day ) )
	def convnormal ( self , normaldate ) :
		arr = normaldate.split()
		try : 
			self.day  = "{}/{}/{}".format( daydic[ arr[1] ], arr[2], arr[-1] ) 
			self.time = arr[-3]
		except : 
			self.day  = arr[0]
			self.time = arr[1]
class qstatclass : 
	def __init__( self ) :
		self.qstat = subprocess.check_output( 'qstat -r | egrep "jun|jobname"', shell=True )
		qstatlines = self.qstat.split('\n')	#split line
		datlen = len(qstatlines)/2
		self.qstatarrind = [ "id","status","day","time", "node", "name", "workdir" ]
		self.qstatarr = [""]* datlen
		for j in range( 0, 2*datlen , 2 ) :
			items  = qstatlines[j].split()
			items2 = qstatlines[j+1].split()
			self.qstatarr[j/2] = [ items[0], items[4], items[5], items[6], items[7] , items2[-1] ]
			dirarr =  subprocess.check_output( 'qstat -j {} | egrep workdir'.format( self.qstatarr[j/2][0] ), shell=True  ) 
			self.qstatarr[j/2].append( dirarr.split()[-1] ) 
	def call( self ) :
		print "qstat.call() : ", self.qstat,
	def show( self , lastdat ) :
		datenow = dateclass()
		pfileind = ["JTD", "S", "J", "U", "D", "beta", "IT" ]
		npar = len(pfileind)
		for q in self.qstatarrind[:6] :
			print "{:5s}".format(q),
		print ""
		queuenamearr = "short_24,short_28,long_24,long_28,single".split(',')
		nq	     = len(queuenamearr)
		njobsqueue   = np.zeros( nq ,dtype=int ) # short_24, short_28, long_24, long_28, single
		for q in self.qstatarr :
			print q[0], "{:3s}".format(q[1]), q[2], q[3], "{:21s}".format(q[4]), #q[5] 

			ns =  q[5].split( "_" )
			beta = q[5][ q[5].find("beta"): ].split("_")[0]
			ns.append( beta ) # =ns[6]
			infostr = ""
			try : 
				it = 0
				infostr = infostr+ "{:8s}".format( ns[it] ) ; it+=1
				infostr = infostr+ "{:6s}".format( ns[it] ) ; it+=1
				infostr = infostr+ "{:7s}".format( ns[it] ) ; it+=1
				infostr = infostr+ "{:6s}".format( ns[it] ) ; it+=1
				infostr = infostr+ "{:6s}".format( ns[it] ) ; it+=1
				if ns[it].find("Hz")>-1 or ns[it].find("IT")>-1 :
					infostr = infostr+ "{:6s}".format( ns[it] ) ; it+=1
				for jj in range(3) :
					if ns[it+jj].find("beta")>-1 :
						infostr = infostr+ "{:8s}".format( ns[it+jj] ) ; it+=1
					if ns[it+jj].find("IT")>-1 :
						infostr = infostr+ "{:6s}".format( ns[it+jj] ) ; it+=1
					elif ns[it+jj].find("next")>-1 :
						break
					else : 
						break
				#infostr = infostr+ "{:8s}".format( ns[6] ) 
			except : pass
			print infostr,
			ns =  q[5].split( "_" , it+jj )
			try : 
				if lastdat=="p" :
					print "{:52s}".format( ns[it+jj] ),
				elif lastdat=="workdir" :
					stdchar = "distortion"
					charcursor =  q[6].find( stdchar ) 
					print "{:52s}".format( q[6][charcursor+len(stdchar)+1:] ),
			except : 
				print "",

			print datenow.elapsedt( q[3], datenow.elapsedd(q[2]) )

			#for qq in q :
			#	print qq,
			#print ""
			for jj in range(nq) :
				if q[4].find( queuenamearr[jj] )>-1 and q[1].find('r')>-1 :
					njobsqueue[jj] += 1 
		for jj in range(nq) :
			print "{}=({}), ".format( queuenamearr[jj] , njobsqueue[jj] ),
		print ""
		print "{:>21s} {:6s}".format( datenow.day , datenow.time )


class pfileclass : 
	def __init__( self ) :
		pass
	def readinfo( self, mak ) :
		self.mak = mak 
		makarr = mak.split("_")
		self.IT  = int(0)
		for j in range(len(makarr)) : 
			if makarr[j].find("J")>-1 : 
				if makarr[j].find("JTD")>-1  :
					self.JTD = float( makarr[j][3:] )
					self.JT  = self.JTD
		       		elif	makarr[j].find("JT")>-1 : 
					self.JT  = float( makarr[j][2:] )
					self.JTD = self.JT
				else :
					self.J   = float( makarr[j][1:] )
			elif makarr[j].find("S")>-1 : 
				self.S   = float( makarr[j][1:] )
			elif makarr[j].find("U")>-1 : 
				self.U   = float( makarr[j][1:] )
			elif makarr[j].find("J")>-1 : 
				self.J   = float( makarr[j][1:] )
			elif makarr[j].find("D")>-1 : 
				self.D   = float( makarr[j][1:] )
			elif makarr[j].find("beta")>-1 : 
				self.beta= int( makarr[j][4:] )
				self.nmax= 4*self.beta
			elif makarr[j].find("IT")>-1 : 
				try : 
					self.IT  = int( makarr[j][2:] )
				except : 
					pass
	def pfilename( self ) :
		self.pfilename = "JTD{}_S{}_J{}_U{}_D{}_beta{}_IT{}_{}.p".format( self.JTD, self.S, self.J, self.U, self.D, self.beta, self.IT, self.mak.split('/')[-1] )
		return self.pfilename

class makclass : 
	def __init__( self , datfile , tNC=6, tNb=18, Ni=1, basis='j' , Tt2gjeff=None ) :
		self.tNC	= tNC
		self.tNb	= tNb
		self.Ni		= Ni
		self.basis	= basis
		self.datfile = datfile 
		if Tt2gjeff is not None : 
			self.Tt2gjeff = Tt2gjeff
		else : 
			self.Tt2gjeff = Taraph
		titlefoot = datfile.split('/')[-1]
		valdatimp	= []
		eldatimp	= []
                with open( datfile ) as fo : 
                        lines = fo.readlines()
			iofflineimp	= 0
			for iimp in range(Ni) :
				ioffline	= 2+iofflineimp
                        	degendat =  int( lines[ioffline].split()[0] ) 
				ioffline	+= 1
                        	eldat =  lines[ioffline+degendat].split()
				ioffline	+= degendat+2
                        	#valdat = [""]*tNC
                        	#for i in range(tNC) :
                        	#        valdat[i] = lines[ioffline+i].split()
                        	valdat	= np.genfromtxt( lines[ioffline:ioffline+tNC], dtype=float ).view(complex)
				#print np.shape(valdat)
				#print valdat
				#sys.exit(1)
				valdatimp.append( valdat )
				eldatimp.append( np.array(eldat,dtype=float) )

				iofflineimp	= ioffline+tNC
				#print "iofflineimp  : ", iofflineimp
				#print "el : " , eldat
				#print "Vml : " , valdat

			#self.irrellines = lines[:2]
			if ( self.tNb != len(eldat) ) : 
				print "ERROR :: tNb is not the same as in .mak file.."
			self.tnb	= len(eldat)
			self.nb		= len(eldat) / 2
		self.ndeg  = degendat
		self.eldat = eldatimp[0]
		self.el    = eldatimp[0]
		self.vmldat= valdatimp[0]
		self.vmlabs= valdatimp[0] * np.conjugate( valdatimp[0] ) ;  self.vmlabs= self.vmlabs.real
		self.eildat	= np.array( eldatimp )
		self.viml	= np.array( valdatimp )
		self.titlefoot = titlefoot

		self.npeldat 	= self.eldat
		self.npvmldat	= self.vmldat

		if basis.find("j")>-1 : 
			self.vjl	= self.vmldat
			self.vjlabs	= self.vmlabs
			self.convertt2g()
		elif basis.find("t")>-1 : 
			self.vtl	= self.vmldat
			self.vtlabs	= self.vmlabs
			self.convertjeff()
	def setTt2gjeff( self , T) :
		self.Tt2gjeff = T
	def returnarr( self ) :
		return [ self.eldat, self.vmldat ]
	def convnparray( self ) :
		pass
		#self.npeldat = np.array( self.eldat , dtype=float ) 
		#self.npvmldat= np.zeros( (nbasis, self.tnb) , dtype=complex)
		#for i in range(nbasis) :
 	        #	for j in range(self.tnb) :
 	        #        	self.npvmldat[i][j] = float(self.vmldat[i][2*j]) + float(self.vmldat[i][2*j+1])*1j
	def returnNpArr( self ) :
		return [ self.npeldat, self.npvmldat ]
	def showNelpn( self ) :
		eln = 0 ; elp = 0
		for i in self.npeldat :
		        if i > 0 : elp+=1
		        else : eln+=1 
		print "negative/positive el : ", eln, '/', elp
		self.nelp = elp 
		self.neln = eln
	def convertt2g( self , T=None ) :
		if T is not None :
			self.Tt2gjeff = T
		else :
			T = self.Tt2gjeff
		if self.basis.find("t")>-1 : 
			pass
		elif self.basis.find("j")>-1 : 
			self.toVtl()
	def convertjeff( self , T=None ) :
		if T is not None :
			self.Tt2gjeff = T
		else :
			T = self.Tt2gjeff
		if self.basis.find("t")>-1 : 
			self.toVjl()
		elif self.basis.find("j")>-1 : 
			pass
	def toVjl( self ) :
		self.vjl	= np.zeros( (nbasis, self.tnb) , dtype=complex)
		for mu in range( nbasis ) :
 	        	for aa in range(self.tnb) :
				for nu in range( nbasis ) :
					self.vjl[mu][aa] += self.Tt2gjeff[mu,nu].conjugate() * self.vtl[nu][aa]
		self.vjlabs	= self.vjl * np.conjugate( self.vjl ) ; self.vjlabs	= self.vjlabs.real
		return self.vjl
	def toVjlSort( self ) :
		self.npvjldat= np.zeros( (nbasis, self.tnb) , dtype=complex)
		for mu in range( nbasis ) :
 	        	for aa in range(self.tnb) :
				for nu in range( nbasis ) :
					self.npvjldat[mu][aa] += self.Tt2gjeff[mu,nu].conjugate() * self.npvmldat[nu][aa]
		self.vjl = self.npvjldat
		return self.npvjldat
	def toVtl( self ) :
		self.vtl	= np.zeros( (nbasis, self.tnb) , dtype=complex)
		for mu in range( nbasis ) :
 	        	for aa in range(self.tnb) :
				for nu in range( nbasis ) :
					self.vtl[mu][aa] += self.vjl[nu][aa] * self.Tt2gjeff[nu,mu]
		self.vtlabs	= self.vtl * np.conjugate( self.vtl ) ; self.vtlabs	= self.vtlabs.real
		return self.vtl
	def toVtlSort( self ) :
		self.npvtldat= np.zeros( (nbasis, self.tnb) , dtype=complex)
		for mu in range( nbasis ) :
 	        	for aa in range(self.tnb) :
				for nu in range( nbasis ) :
					self.npvtldat[mu][aa] += self.npvmlsort[nu][aa] * self.Tt2gjeff[nu,mu]
		self.vtlsort =  self.npvtldat
		return self.npvtldat
	def toVtlusort( self ) :
		self.npvtldatusort= np.zeros( (nbasis, self.tnb) , dtype=complex)
		for mu in range( nbasis ) :
 	        	for aa in range(self.tnb) :
				for nu in range( nbasis ) :
					self.npvtldatusort[mu][aa] += self.npvmldat[nu][aa] * np.conj( self.Tt2gjeff[nu,mu] )
		self.vtl =  self.npvtldatusort
		return self.npvtldatusort
	def toVnl( self ) :
		self.npvnldat= np.zeros( (nbasis, self.tnb) , dtype=complex)
		for mu in range( nbasis ) :
 	        	for aa in range(self.tnb) :
				for nu in range( nbasis ) :
					self.npvnldat[mu][aa] += self.npvmlsort[nu][aa] * matTjnat[mu,nu] 
		return self.npvnldat
	def toVjlfromnat( self ) :
		self.npvjldat= np.zeros( (nbasis, self.tnb) , dtype=complex)
		for mu in range( nbasis ) :
 	        	for aa in range(self.tnb) :
				for nu in range( nbasis ) :
					self.npvjldat[mu][aa] += self.npvmlsort[nu][aa] * matTjnat[nu,mu]
		return self.npvjldat
	def toVfromjeff( self , T ) :
		self.npvfromjeff = np.zeros( (nbasis, self.tnb) , dtype=complex)
		for mu in range( nbasis ) :
 	        	for aa in range(self.tnb) :
				for nu in range( nbasis ) :
					self.npvfromjeff[mu][aa] += self.npvmldat[nu][aa] * np.conj( T[nu,mu] )
		return self.npvfromjeff
	def sortVmlel( self ) :
		sort_perm = self.npeldat.argsort()
		self.npvmlsort    = self.npvmldat[:, sort_perm ]
		self.npelsort     = self.npeldat[sort_perm ]
		return self.npvmlsort, self.npelsort
	def returnDeltanumersum( self ) :
		self.Deltanumersum = np.zeros( (nbasis, nbasis) , dtype=complex)
		for i in range(nbasis) :
			for j in range(nbasis) :
				for l in range(self.tnb) :
					self.Deltanumersum[i][j] += self.npvmldat[i][l] * np.conj( self.npvmldat[j][l] )
		return self.Deltanumersum
	def returnDeltaw( self , w) :
		self.Deltaw = np.zeros( (nbasis, nbasis) , dtype=complex)
		for i in range(nbasis) :
			for j in range(nbasis) :
				for l in range(self.tnb) :
					self.Deltaw[i][j] += self.npvmldat[i][l] * np.conj( self.npvmldat[j][l] ) / ( w - self.npeldat[l] )
		return self.Deltaw
	def returnDeltaiw( self , w) :
		self.Deltaw = np.zeros( (nbasis, nbasis) , dtype=complex)
		for i in range(nbasis) :
			for j in range(nbasis) :
				for l in range(self.tnb) :
					self.Deltaw[i][j] += self.npvmldat[i][l] * np.conj( self.npvmldat[j][l] ) / ( w*1j - self.npeldat[l] )
					#print "{} {} {:3d}".format(i,j,l), 
					#print "{:19.6f} {:19.6f} {:19.6f} {:19.6f}".format( self.npvmldat[i][l] * np.conj( self.npvmldat[j][l] ) , ( w*1j - self.npeldat[l] ) , 1./ ( w*1j - self.npeldat[l] ) , self.Deltaw[i][j] )
		return self.Deltaw
	def returnBathftn( self ) :
		self.Bathftn = np.zeros( (self.tnb, self.tnb) , dtype=complex)
		for i in range(self.tnb) :
			for j in range(self.tnb) :
				for l in range(nbasis) :
					self.Bathftn[i][j] += self.npvmlsort[l][i] * np.conj( self.npvmlsort[l][j] ) 
		return self.Bathftn
	def returnCalcDelta( self , wn ) :
		Deltawn = np.zeros( (nbasis, nbasis) , dtype=complex)
		for i in range(nbasis) :
			for j in range(nbasis) :
				for l in range(nbasis) :
					Deltawn[i][j] += np.conj( self.npvmldat[i][l] ) * self.npvmldat[j][l] / (wn-self.npeldat[j])
		return Deltawn
	def writemak( self , datname, footname ) :
		fnameind = self.datfile.find(".mak")
		if fnameind > -1 : 
			fname    = self.datfile[:fnameind] + footname + ".mak" 
		else : 
			print "Need to specify mak-file explicitly."
			exit(1)
		print "Write : ", fname
		fwmak = open( fname, 'w' )
		with open( self.datfile ) as fo :
		        lines = fo.readlines()
		        ndeg = int( lines[2] )
		        for j in range( 3+ndeg ) :
		                fwmak.write( lines[j] )
		        for eli in self.npeldat :
		                fwmak.write( "{:20.17f}\t\t\t\t".format( eli )  )
		        fwmak.write("\n\n")
		        for vlmu in getattr(self,datname) :
		                for vlmui in vlmu :
		                        fwmak.write( "{:20.17f} {:20.17f}\t".format( vlmui.real, vlmui.imag )  )
		                fwmak.write("\n")
		        fwmak.write("\n\n")
		fwmak.close()
	def writemaksort( self , datname, footname ) :
		fnameind = self.datfile.find(".mak")
		if fnameind > -1 : 
			fname    = self.datfile[:fnameind] + footname + ".mak" 
		else : 
			print "Need to specify mak-file explicitly."
			exit(1)
		print "Write : ", fname
		fwmak = open( fname, 'w' )
		with open( self.datfile ) as fo :
		        lines = fo.readlines()
		        ndeg = int( lines[2] )
		        for j in range( 3+ndeg ) :
		                fwmak.write( lines[j] )
		        for eli in self.npelsort :
		                fwmak.write( "{:20.17f}\t\t\t\t".format( eli )  )
		        fwmak.write("\n\n")
			i=0
		        for vlmu in getattr(self,datname) :
				print i, vlmu ; i+=1
		                for vlmui in vlmu :
		                        fwmak.write( "{:20.17f} {:20.17f}\t".format( vlmui.real, vlmui.imag )  )
		                fwmak.write("\n")
		        fwmak.write("\n\n")
		fwmak.close()
	def writenpmaksort( self , datname, footname ) :
		fnameind = self.datfile.find(".mak")
		if fnameind > -1 : 
			fname    = self.datfile[:fnameind] + footname + ".mak" 
		else : 
			print "Need to specify mak-file explicitly."
			exit(1)
		print "Write : ", fname
		fwmak = open( fname, 'w' )
		with open( self.datfile ) as fo :
		        lines = fo.readlines()
		        ndeg = int( lines[2] )
		        for j in range( 3+ndeg ) :
		                fwmak.write( lines[j] )
		        for eli in self.npelsort :
		                fwmak.write( "{:20.17f}\t\t\t\t".format( eli )  )
		        fwmak.write("\n\n")
		        vlmu = getattr(self,datname) 
		        for mu in range( np.shape(vlmu)[0] ) :
		        	for i in range( np.shape(vlmu)[1] ) :
		                        fwmak.write( "{:20.17f} {:20.17f}\t".format( vlmu[mu,i].real, vlmu[mu,i].imag )  )
		                fwmak.write("\n")
		        fwmak.write("\n\n")
		fwmak.close()
	def plotAbs( self , ax, TRANS=False ) :
		self.plotAbsImp( ax=ax, TRANS=TRANS )
	def plotAbsImp( self , ax, TRANS=False , iimp=0 ) :
		width = 0.09
		hatchftn = [ '' ]*nbasis
		colorftn = [ 'b' ]*nbasis
		el	= self.eildat[iimp]
		vlabs	= self.viml[iimp]
		vlabs	= vlabs * np.conjugate(vlabs)
		vsumbottom = np.zeros( (len(el)) )
		if TRANS is False :
			for i in range(nbasis) :
				cs = colorftn[i]
		        	ax.barh( el, vlabs[i], width, left=vsumbottom, alpha=(9.95-1.8*i)/10. , color=cs , label=r'$\mu={}$'.format(i) , hatch=hatchftn[i] )
				for j in range(len(el)) :
					vsumbottom[j] += vlabs[i][j]
			ax.plot( vsumbottom, el, 'ko' , ms=2 )
		else :
			for i in range(nbasis) :
				cs = colorftn[i]
		        	ax.bar( el, vlabs[i], width, bottom=vsumbottom, alpha=(9.95-1.8*i)/10. , color=cs , label=r'$\mu={}$'.format(i) , hatch=hatchftn[i] )
				for j in range(len(el)) :
					vsumbottom[j] += vlabs[i][j]
			ax.plot( el, vsumbottom, 'ko' )
		titleax = ax
		labelax = ax
		#ax.legend()
		
		x1,x2,y1,y2 = titleax.axis()
		#ax.axis((-0.1,x2,y1,y2))
		a=0
		for i,j in zip(el,[0]*len(el)):
			corr = 0.02*(x2-x1) # adds a little correction to put annotation in marker's centrum
			#ax.annotate(str(a),  xy=(i + 0.01, j + corr*a))
			titleax.annotate(str(a),  xy=( (x2-x1)/2. + j + corr*a, i + 0.01 ), fontsize=5 , rotation=-45)
			a+=1
		#print "negative/positive el : ", eln, '/', elp
		self.showNelpn()

def pltSetRange( plt, xystr, xyind ) :
	if xystr : 
		if xyind.find("y") > -1 :
			plt.ylim( float( xystr.split("_")[-2]) , float( xystr.split("_")[-1])  ) 
		if xyind.find("x") > -1 :
			plt.xlim( float( xystr.split("_")[-2]) , float( xystr.split("_")[-1])  ) 


class plotdat  :
	def __init__(self) :  pass
	def set( self, ax=None, xdat=None, ydat=None, lt=None, lc=None, ls=None, mec=None, dashes=None , tag=None ) :
		self.ax	= ax
		self.xdat	= xdat
		self.ydat	= ydat
		self.lt	= lt
		self.lc	= lc
		self.ls	= ls
		self.mec	= mec
		self.dashes	= dashes
		self.tag	= tag

class Crystalclass  :
	def __init__(self) :  
		self.ni			= 1
		self.nbasis		= 6
		self.nbasish		= 3
	def readW90win( self, path , ni=1 , nbasis=6 ) :
	        self.ni		= ni
	        self.nbasis	= nbasis
	        nbasish	= nbasis/2
	        self.nbasish	= nbasish
	        if path.find("wannier90.win")>-1 :
	                fr = open( path , 'r' )
	                datlines = fr.readlines()
	                #print datlines
	
	        nlines = len(datlines)
	        self.avec = [""]*3
	        for i in range(nlines):
	                ii = i
	                line = datlines[ii]
	                if line.find( "begin unit_cell_cart" )>-1  :
	                        #print datlines[ii-1].splitlines()[0].split()
	                        #print datlines[ii-2]
	                        #print datlines[ii-3]
				print "READING avec :"
	                        self.avec[0] = np.array(  datlines[ii+1].splitlines()[0].split() , dtype=float )
	                        self.avec[1] = np.array(  datlines[ii+2].splitlines()[0].split() , dtype=float )
	                        self.avec[2] = np.array(  datlines[ii+3].splitlines()[0].split() , dtype=float )
	                        print self.avec[0]
	                        print self.avec[1]
	                        print self.avec[2]
				break
	        self.atomposFrac = [""]*ni
	        self.localaxz = [""]*ni
	        self.localaxx = [""]*ni
	        self.localaxy = [""]*ni
	        for i in range(nlines):
	                ii = i
	                line = datlines[ii]
	                if line.find( "begin projections" )>-1  :
	                        #print datlines[ii-1].splitlines()[0].split()
	                        #print datlines[ii-2]
	                        #print datlines[ii-3]
			
				print "READING atomposFrac :"
				for a in range(ni) :
					line = datlines[ii+1+a*nbasish]
					istart = line.find("=")+1
					iend = line.find(":")
	                        	self.atomposFrac[a] = np.array(  line[istart:iend].split(",") , dtype=float )
	                        	print self.atomposFrac[a]

				print "READING local axis :"
				
				for a in range(ni) :
					line = datlines[ii+1+a*nbasish]
					linesp = line.split(":")
					for il in range(len(linesp)) :
						ifind = linesp[il].find("z=")
						if ifind>-1 :
							self.localaxz[a] = np.array( linesp[il][ifind+2:].split(",")[:3] , dtype=float ) 
						ifind = linesp[il].find("x=")
						if ifind>-1 :
							self.localaxx[a] = np.array( linesp[il][ifind+2:].split(",")[:3] , dtype=float ) 
					try : 
	                        		self.localaxy[a]  = np.cross( self.localaxz[a] , self.localaxx[a] )
					except : 
	      					print "May NOT using local axes in {}-th's.".format(a)
	                        		self.localaxx[a] = np.array( [1,0,0] ) 
	                        		self.localaxy[a] = np.array( [0,1,0] ) 
	                        		self.localaxz[a] = np.array( [0,0,1] ) 
	                        	print "local axis-[x',y',z'] as f(x,y,z) : ", self.localaxx[a], self.localaxy[a], self.localaxz[a] 
	                        	x = np.linalg.norm( self.localaxx[a] ) 
	                        	y = np.linalg.norm( self.localaxy[a] ) 
	                        	z = np.linalg.norm( self.localaxz[a] ) 
	                        	self.localaxx[a] = self.localaxx[a] /x
	                        	self.localaxy[a] = self.localaxy[a] /y
	                        	self.localaxz[a] = self.localaxz[a] /z
	      			 
	      			#print "Not using local axes."
				#for a in range(ni) :
	                        #	self.localaxx[a] = np.array( [1,0,0] ) 
	                        #	self.localaxy[a] = np.array( [0,1,0] ) 
	                        #	self.localaxz[a] = np.array( [0,0,1] ) 
				break
		print "local axis-[x',y',z'] as f(x,y,z) (normalized) : "
		for a in range(ni) :
			print self.localaxx[a], self.localaxy[a], self.localaxz[a]
		print "Transforming atomposFrac to atomposCart :"
	        self.atomposCart = [""]*ni
		for a in range(ni) :
			self.atomposCart[a] = np.dot( self.atomposFrac[a] , self.avec ) 
			print self.atomposCart[a] 

		dumiOxy=0
		dumnOxy=0
	        for i in range(nlines):
	                line = datlines[i]
			j=0
			dumiOxy = 0
			#print i , line , line.find("begin atoms_frac")>-1 
			if line.find("begin atoms_frac")>-1 : 
				while( dumiOxy<1 ) : 
					if datlines[i+1+j].find("O ")>-1 : 
						dumiOxy = j
					else : 
						j+=1
				dumnOxy = 0
				while( datlines[i+1+j].find("end atoms_frac")<0 ) : 
					if datlines[i+1+j].find("O ")>-1 : 
						dumnOxy += 1
					else : 
						pass
					j+=1
				break
		print "iOxy : ", dumiOxy
		print "nOxy : ", dumnOxy
		self.iOxy = dumiOxy
		self.nOxy = dumnOxy

	def readPOSCAR( self, path ) :
	        ni	= self.ni
	        nbasis	= self.nbasis
	        nbasish	= self.nbasish
	        if path.find("POSCAR")>-1 :
	                fr = open( path , 'r' )
	                datlines = fr.readlines()
	                #print datlines
	
	        nlines = len(datlines)
		try : 
			if self.nOxy>0 :
				nOxy = self.nOxy  
				iOxy = self.iOxy  
			else : 
				nOxy = 16
				iOxy = 12
		except : 
			nOxy = 16
			iOxy = 12
		print "iOxy : ", iOxy
		print "nOxy : ", nOxy
	        self.atomposFracO = [""]*nOxy
	        for i in range(nlines):
	                ii = i
	                line = datlines[ii]
	                if line.find( "irect" )>-1  :
	                        #print datlines[ii-1].splitlines()[0].split()
	                        #print datlines[ii-2]
	                        #print datlines[ii-3]
			
				print "READING atomposFracO :"
				for a in range(nOxy) :
					line = datlines[ii+1+iOxy+a]
	                        	self.atomposFracO[a] = np.array(  line.split()[:3] , dtype=float )
	                        	print self.atomposFracO[a]
				break
		print "Transforming atomposFracO to atomposCartO :"
	        self.atomposCartO = [""]*nOxy
		for a in range(nOxy) :
			self.atomposCartO[a] = np.dot( self.atomposFracO[a] , self.avec ) 
			print self.atomposCartO[a] 
		
	def readMoment( self, hobj, op="J" ) :
		print "Reading op"+op+" : "
		localaxz= self.localaxz
		localaxx= self.localaxx
		localaxy= self.localaxy
		ni	= self.ni
	        vecFrac = [""]*ni
		for a in range(ni) :
			vecFrac[a] = getattr( hobj, "op"+op+"%d"%a ) 
			print op+"-vecFrac %d : "%a, vecFrac[a] 

	        vecCart = [""]*ni
		for a in range(ni) :
			#print [ localaxx[a], localaxy[a], localaxz[a] ], vecFrac[a]
			#vecCart[a] = np.dot( [ localaxx[a], localaxy[a], localaxz[a] ], vecFrac[a] )
			vecCart[a] = localaxx[a] * vecFrac[a][0]  + localaxy[a] * vecFrac[a][1] + localaxz[a] * vecFrac[a][2] 
			print op+"-vecCart,direct,norm %d : "%a, vecCart[a] , vecCart[a] / np.linalg.norm( vecCart[a] ), np.linalg.norm( vecCart[a] )

	        self.vecFrac = vecFrac 
	        self.vecCart = vecCart 

        def writeRLocalt2g( self , fname='TRANSFLOCALt2g') :
                localaxz= self.localaxz
                localaxx= self.localaxx
                localaxy= self.localaxy
                ni      = self.ni
                dict2gToCart = {        'xy' : 2 ,
                                        'yz' : 0 ,
                                        'zx' : 1        }
                RLocalt2g = [""]*ni
                f = open( fname , 'w' )
                orderCart        = [ dict2gToCart[ii] for ii in ['xy','yz','zx'] ]
                for i in range(ni) :
                        RLocalt2g[i]     = np.array([ localaxz[i][orderCart], localaxx[i][orderCart], localaxy[i][orderCart] ]).transpose()
                        np.savetxt( f, RLocalt2g[i] ,fmt='%20.16f')
                        f.write("\n")
                f.close()
                self.RLocalt2g = RLocalt2g
        def writeRLocalspinhalf( self , fname="TRANSFLOCALspinhalf" ,) :
                localaxz= self.localaxz
                localaxx= self.localaxx
                localaxy= self.localaxy
                ni      = self.ni
                RLocalspinhalf = [""]*ni
                from scipy.spatial.transform import Rotation as R
                f = open( fname , 'w' )
                for i in range(ni) :
                        ntmp    = np.array( [ localaxx[i],localaxy[i],localaxz[i] ] )
                        r = R.from_dcm( ntmp )
                        v = r.as_rotvec() ; vn = np.linalg.norm(v)
                        print "rotvec : ", vn , v/vn
                        print "RLocalspinhTest : ", i
                        pA(  Rspinhalf( v/vn, -vn ) )
                        RLocalspinhalf[i] = Rspinhalf( v/vn, -vn )
			if np.abs(vn)<1e-4 :
                        	RLocalspinhalf[i] = np.array( [ [1.+0.*1j,0.],[0.,1.+0.*1j] ] )
                        np.savetxt( f, RLocalspinhalf[i].view(float) ,fmt='%20.16f')
                        f.write("\n")
                f.close()
                self.RLocalspinhalf = RLocalspinhalf
        def writeRLocalspinhalfBlock( self , fname="TRANSFBLOCKLOCALspinhalf" ) :
		ni	= self.ni
		nbasis	= self.nbasis
		nbasish	= self.nbasish
		ntot	= ni*nbasis

                RLocalspinhalfimp		= [""]*ni
                for i in range(ni) :
			idm	= np.eye( nbasish ) 
                        RLocalspinhalfimp[i] = np.kron( idm, self.RLocalspinhalf[i] )
                        print "DP : ", i
			pMr( RLocalspinhalfimp[i] )
                self.RLocalspinhalfimp		= np.array( RLocalspinhalfimp )
                RLocalspinhalfimpBlock		= np.zeros( (ntot,ntot) , dtype=complex)
                for i in range(ni) :
                	RLocalspinhalfimpBlock[i*nbasis:(i+1)*nbasis, i*nbasis:(i+1)*nbasis] = RLocalspinhalfimp[i] +  1j*np.eye(nbasis)*0
                self.RLocalspinhalfimpBlock		= np.array( RLocalspinhalfimpBlock )

                print "Writing RLocalspinhalfimpBlock ", np.shape( self.RLocalspinhalfimpBlock )
                f = open( fname , 'w' )
                np.savetxt( f, self.RLocalspinhalfimpBlock.view(float) ,fmt='%20.16f')
                #for i in range(ntot) :
                #	print self.RLocalspinhalfimpBlock[i].view(float) 
                #	np.savetxt( f, self.RLocalspinhalfimpBlock[i].view(float) ,fmt='%20.16f')
                #	f.write("\n")
                f.close()
                pMr(self.RLocalspinhalfimpBlock)
        def writeRLocalspinhalfSiteBlock( self , fname="TRANSFSITEBLOCKLOCALspinhalf" ) :
		try : 
			print "Checking RLocalspinhalfimp... ", np.shape( self.RLocalspinhalfimp ), " (done)"
		except : 
        		self.writeRLocalspinhalfBlock( fname="tmp" )
			print "Checking RLocalspinhalfimp... ", np.shape( self.RLocalspinhalfimp ), " (done)"

		# M is in local-axes, Mout would be in global-axes
		ni			= self.ni
		RLocalspinhalfimp	= self.RLocalspinhalfimp
		print "Checking dimension :  ", np.shape(RLocalspinhalfimp)
                f = open( fname , 'w' )
                for i in range(ni) :
                        np.savetxt( f, RLocalspinhalfimp[i].view(float) ,fmt='%20.16f')
                        f.write("\n")
                f.close()
        def writeProjectRGlobalspinhalfSiteBlock( self , M, fname="MATBlockGlobalspinhalf" ) :
		try : 
			print "Checking RLocalspinhalfimp... ", np.shape( self.RLocalspinhalfimp ), " (done)"
		except : 
        		self.writeRLocalspinhalfBlock( fname="tmp" )
			print "Checking RLocalspinhalfimp... ", np.shape( self.RLocalspinhalfimp ), " (done)"

		# M is in local-axes, Mout would be in global-axes
		print "Checking dimension :  ", np.shape(M)
		ni			= self.ni
		RLocalspinhalfimp	= self.RLocalspinhalfimp
                MoutBlock		= [""]*ni
                for i in range(ni) :
			#print "Checking unitarity ..."
			#pA(np.dot( np.conjugate(RLocalspinhalfimp[i].transpose()),  RLocalspinhalfimp[i] ) )
			MoutBlock[i]	= np.dot( np.conjugate(RLocalspinhalfimp[i].transpose()) ,  np.dot( M, RLocalspinhalfimp[i] ) )
			#print "Checking hermiticity ..."
			#pA( np.conjugate(MoutBlock[i].transpose()) - MoutBlock[i] )

                f = open( fname , 'w' )
                for i in range(ni) :
                        np.savetxt( f, MoutBlock[i].view(float) ,fmt='%20.16f')
                        f.write("\n")
                f.close()
        def writeRLocalcart( self , fname='TRANSFLOCALcart') :
                localaxz= self.localaxz
                localaxx= self.localaxx
                localaxy= self.localaxy
                ni      = self.ni
                dicCart = {        'x' : 0 ,
                                        'y' : 1 ,
                                        'z' : 2        }
                RLocalcart = [""]*ni
                f = open( fname , 'w' )
                orderCart        = [ dicCart[ii] for ii in ['x','y','z'] ]
                for i in range(ni) :
                        RLocalcart[i]     = np.array([ localaxx[i][orderCart], localaxy[i][orderCart] , localaxz[i][orderCart] ])
                        np.savetxt( f, RLocalcart[i] ,fmt='%20.16f')
                        f.write("\n")
                f.close()




class datafilt :
	def __init__(self , fname ) :
		self.fname = fname
	def filt( self, filtarr=False, comparr= [0,5]  , filttype=str ) :
		self.filttype=filttype
		self.dat = np.genfromtxt( self.fname , dtype=float ).transpose()
		self.dats= np.genfromtxt( self.fname , dtype=filttype ).transpose()
		if filtarr :
			print "filt-type : ", filtarr
			colFilt , filtname = filtarr
			self.colArr = self.colArrFilt( colFilt, filtname )
			return self.dat[ comparr ].transpose()[self.colArr].transpose()
		else : 
			return self.dat[ comparr ]
	def colArrFilt( self , col , filtname ) :
		nd	= len(self.dats[col])
		self.colArr	= range(nd)
		for i in range(nd) :
			if self.filttype == str : 
				if self.dats[col][i].find(filtname)>-1 : 
					self.colArr[i] = True
				else :
					self.colArr[i] = False
			elif self.filttype == float : 
				if self.dats[col][i] == filtname : 
					self.colArr[i] = True
				else :
					self.colArr[i] = False
		return self.colArr

class tmpclass :
	def __init__(self) :
		pass
class txtclass :
	def __init__(self) :
		self.objarr = []
		self.nobj = 0
	def getData( self , txtname ) :
		with open(txtname,'r') as fr :
			line0 = fr.readlines()[0]
		#if line0[0]=='#' :
		# 	line0 = line0[1:]
		itemArr = line0.split()
		dataArr = np.genfromtxt(txtname)
		nitem	= len(itemArr)
		nobj	= len(dataArr)

		print "Item : ", 
		for jj in range(nitem) :
			print itemArr[jj],
		print ""

		self.nobj = nobj
		self.objarr = [""] * nobj
		print nobj, "-data loaded"
		for ii in range(nobj) :
			self.objarr[ii] = tmpclass()
			obj = self.objarr[ii] 
			for jj in range(nitem) :
				setattr( obj, itemArr[jj], dataArr[ii][jj] ) 

