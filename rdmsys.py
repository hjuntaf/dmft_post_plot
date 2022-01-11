#!/opt/python/2.7.15/gcc-4.8.5/bin/python

import numpy as np
import sys
import argparse

sys.path.append("/home/jun/tools/dmft/plot/" ) 
from hjl_class import *
from hjl_common import *

parser = argparse.ArgumentParser(description="Rearranging and showing data")
parser.add_argument('dirs', metavar='vars', type=str, nargs='+', help="Collection of args" ) 
parser.add_argument("-a", "--all",  action="store_true", help="Show the original complex-data." )
parser.add_argument("-c", "--comp",  action="store_true", help="Show the original complex-data only." )
parser.add_argument("-d", "--div",  action="store_true", help="Show the partition of d4jjJ or d4LSJ." )
parser.add_argument("-l", "--line",  type=str,  help="Show the specified line." )
parser.add_argument("--dtype", type=str,  help="Specify the type of the data, e.g. float or complex.") 
parser.add_argument("-t", "--type", type=str,  help="Specify data that you want, e.g. rdmevald4, rdmevecd4jjJ, rdmevecd4_rdmevald4, or etc." )
parser.add_argument("-p", "--plot", type=str,  help="Plot data for rdmeval." ) 
args = parser.parse_args()

if	args.dirs[0].find('d4jjJ')>-1 : argdiv = [ 1 , 4 , 9 , 14, 15 ]
elif	args.dirs[0].find('d4LSJ')>-1 : argdiv = [ 1 , 4 , 9 , 14, 15 ]

if args.line : indline = int(args.line) ; print "argsline : ", args.line, indline
else :		indline = -2 

datnamearr	= args.type.split("_") 
ndat		= len(datnamearr)
datoutarr	= range(ndat)
lengtharr	= range(ndat)
lengthorigarr	= range(ndat)
dimarr		= range(ndat)
if args.dirs[0].find("Dir")>-1 :
	hobj = headobj( args.dirs[0].split("/")[0] )
	hobj.readParameters()
	#print "dir(hobj) : ", dir(hobj)
	print "hobj.degen_0 : ", hobj.degen_0

	fnamearr  = [ hobj.pathResult()+"/"+name+".dat"  for name in datnamearr ]
	datarr = [ np.genfromtxt( hobj.pathResult()+"/"+name+".dat" , dtype=float ) for name in datnamearr ]

	for j in range(ndat) :
		print "fname {} : {}".format( j, fnamearr[j] )  
		datnow = datarr[j]
		print "shape(dat_i) : ",  np.shape(datnow) 
		shapedatnow = np.shape(datnow) 
		try : 
			length = shapedatnow[1] - 1 
			if int(hobj.degen_0) < 2 : 
				datout = [ datnow[indline] ]  
			else : 
				datout = [ datnow[indline],datnow[indline+1] ] 
		except : 
			length = shapedatnow[0] - 1 
			datout = [ datnow ]
		lengthorig = length 
		length /= 2 
		if args.dtype : 
			if args.dtype.find("f")>-1 or  args.dtype.find("float")>-1 : 
				length *= 2
		dim = int(np.sqrt(length)) #dim of each-iteration data
		lengtharr[j]	= length
		lengthorigarr[j]= lengthorig
		dimarr[j]	= dim
		datoutarr[j]	= datout
else : 
	print "ERROR :: Invalid filename."


#print "length	: ", length
#print "dirarr	: ", datarr
#print "dirarr[0] : ", datarr[0]
#print "lengthorigarr	: ", lengthorigarr
#print "lengtharr	: ", lengtharr
#print "dimarr	: ", dimarr
for dd in range(ndat) : 
	datout		= datoutarr[dd]
	dim		= dimarr[dd]
	length		= lengtharr[dd]
	dname		= datnamearr[dd]
	print datnamearr[dd]
	if length == dim*dim :
		print "Square-matrix of ::", dim, length
		if args.all or args.comp : 
			for datline in datout :
				print "GS [%d]: "%indline
				for i in range( dim )  :
					for j in range( dim )  :
						print "{:6.3f}+i{:6.3f}\t".format(datline[1+2*dim*i+2*j],datline[2+2*dim*i+2*j]),
					print ""
				print ""
		if args.comp : pass 
		else : 
			for i in range( dim )  : print "{:8d} ".format(  i   ),
			print ""
			for datline in datout :
				print "GSabssq [%d]: "%indline
				for i in range( dim )  :
					inddiv = 0 
					for j in range( dim )  :
						if args.div : 
							if j<argdiv[inddiv] : pass
							else : print "|", ; inddiv +=1
						print "{:8.5f} ".format(  np.absolute( datline[1+2*dim*i+2*j]+ 1j* datline[2+2*dim*i+2*j] )**2  ),
					print ""
				print ""
	else :
		length = lengthorigarr[dd]
		print "Non-square-matrix data ::", dim, length
		for i in range( length )  : print "{:8d} ".format(  i   ),
		print ""
		for datline in datout :
			print "GSeval [%d]: "%indline
			for i in range( length )  :
				print "{:8.5f} ".format(  datline[1+i] ),
			print ""

	nbasis = 6 
	if dim==nbasis :
		datMat66 = np.zeros( (nbasis, nbasis) , dtype=complex)
		for i in range( dim )  :
			for j in range( dim )  :
				datMat66[i][j] = dat[indline][1+2*dim*i+2*j] + 1j* dat[indline][2+2*dim*i+2*j]
		print "datMat66 : "  
		printArrMat( datMat66 ) 
		datMat66nat = transfconjTOtransT( datMat66, matTjnat ) 
		print "datMat66nat : "  
		printArrMat( datMat66nat ) 

