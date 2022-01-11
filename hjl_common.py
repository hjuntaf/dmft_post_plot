import numpy as np
import sys, os
import subprocess

from basislabel import *
from hjl_class   import *

from matplotlib.lines import Line2D
import numpy as np

t = np.arange(0.0, 1.0, 0.1)
s = np.sin(2*np.pi*t)
linestyles = [ '-', '--', ':', '-.']*2
markers = []
#for m in Line2D.markers:
#    try:
#        if len(m) == 1 and m != ' ':
#            markers.append(m)
#    except TypeError:
#        pass
markers = ["o", "v", "^", "<", ">", "s" ]*2

#styles = markers + [
#    r'$\lambda$',
#    r'$\bowtie$',
#    r'$\circlearrowleft$',
#    r'$\clubsuit$',
#    r'$\checkmark$']
def lscycler(x) : return linestyles[x%3]
def ptcycler(x) : return markers[x%6]

class dumclass :
	def __init__(self) :
		pass

def returnDatafname( dataname , dat, hobj, count, fsuffix="sus1.20" , epsilon=None ) : 
	responsefolder = fsuffix
	tail_epsilon = ""
	if epsilon : 
		tail_epsilon = "_ep%.3lf"%epsilon
		print "TAIL : " , tail_epsilon
	if dataname.find("hybiw")>-1 :
		hybiwfolder = "ongoing"
		datafname = "{}/{}/bathinverse{}_u{:.2f}_{}th".format(dat, hybiwfolder, hobj.ni-1, hobj.UF, count )
	elif dataname.find("chiiwmzmz")>-1 :
		try : 
			datafname = "{}/{}/chitotmzmz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , hobj.upper )
		except :
			datafname = "{}/{}/chitotmzmz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , 100 )
	elif dataname.find("chiwmzmz")>-1 :
		try : 
			datafname = "{}/{}/real_chitotmzmz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , hobj.upper )
		except :
			datafname = "{}/{}/real_chitotmzmz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , 100 )
	elif dataname.find("chiiwjzjz")>-1 :
		try : 
			datafname = "{}/{}/chitotjzjz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , hobj.upper )
		except :
			datafname = "{}/{}/chitotjzjz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , 100 )
	elif dataname.find("chiwjzjz")>-1 :
		try : 
			datafname = "{}/{}/real_chitotjzjz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , hobj.upper )
		except :
			datafname = "{}/{}/real_chitotjzjz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , 100 )
	elif dataname.find("chiiwszsz")>-1 :
		try : 
			datafname = "{}/{}/chitotszsz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , hobj.upper )
		except :
			datafname = "{}/{}/chitotszsz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , 100 )
	elif dataname.find("chiwszsz")>-1 :
		try : 
			datafname = "{}/{}/real_chitotszsz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , hobj.upper )
		except :
			datafname = "{}/{}/real_chitotszsz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , 100 )
	elif dataname.find("chiiwlzlz")>-1 :
		try : 
			datafname = "{}/{}/chitotlzlz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , hobj.upper )
		except :
			datafname = "{}/{}/chitotlzlz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , 100 )
	elif dataname.find("chiwlzlz")>-1 :
		try : 
			datafname = "{}/{}/real_chitotlzlz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , hobj.upper )
		except :
			datafname = "{}/{}/real_chitotlzlz_u{:.2f}_{}th_n{}_upper{}".format(dat, responsefolder, hobj.UF, count, hobj.ni-1 , 100 )
	datafname = datafname + tail_epsilon + ".txt"
	return datafname
class autoheadobj  :
	def __init__(self, ddir , varlist ) :
		self.ddir = ddir 
		self.tdir = ddir 
		ddirsp = ddir.split("_")
		for item in ddirsp :
			for v in varlist :
				var    =  v.split("_")[0] 
				posind =  item.find( var ) 
				if posind > -1 : 
					try :
						setattr( self, v , item[posind+len(var):] ) 
					except : pass
	def headSetTitle(self, titleax ) :
		titleax.set_title( "" ) 
	def readMakClass(self) :
		pass

class headobj  :
	def __init__(self, ddir ) :
		if ddir[-1]=="/" : 
			self.tdir = ddir[:-1]
		else : 
			self.tdir = ddir 
		self.ddir = self.tdir
		for path in self.tdir.split("/") :
			if path.find("Dir")>-1 : 
				self.Dirname = path
		self.path = self.ddir[:self.ddir.find("Dir")-1]
		if (self.Dirname.find("nDir")>-1) or (self.Dirname.find("oDir")>-1) :
			self.basis	 = 't'
		elif (self.Dirname.find("jDir")>-1) or (self.Dirname.find("iDir")>-1) :
			self.basis	 = 'j'
		else : 
			print "ERROR :: invalid header in your directory."
			sys.exit(1)

		try : 
			if self.path.find("wnalpha")>-1 and self.path.find("_")<1 :
				pathalpha	=  self.path[ self.path.find("wnalpha") : ]
				self.alpha =  float( pathalpha[ pathalpha.find("wnalpha")+7 : pathalpha.find("/") ] ) 
			if self.path.find("2dNNN")>-1 and self.path.find("_")<1 :
				pathpre	=  self.path[ self.path.find("2dNNN") : ]
				self.tNNN =  float( pathpre[ pathpre.find("2dNNN")+5 : pathpre.find("/") ] ) 
			elif self.path.find("tNNN")>-1 and self.path.find("_")<1 :
				pathpre	=  self.path[ self.path.find("tNNN") : ]
				self.tNNN =  float( pathpre[ pathpre.find("tNNN")+4 : pathpre.find("/") ] ) 
			else : self.tNNN  = 0.1
		except :
			self.alpha  = 0.

		try : 
			ddirsp = self.Dirname.split("_")
		except : 
		       print "ERROR ::"
		       print "check : ", self.tdir
		       print "check : ", ddir
		       print "check : ", self.tdir.split("/") 
		       sys.exit(1)
		if ddirsp[-1].find("R") > -1 :
			self.tdirind = ""
		else :
			self.tdirind = ddirsp[-1].split('/')[0]
		self.IT=int(0)
		for t in ddirsp :
			if( t.find("U") > -1 ) : 
				if( t.find("UF") > -1 ) : 
					self.UF = float(t[2:])
					self.u = self.UF
					self.U = self.UF
				elif( t.find("UT") > -1 ) : 
					self.UT = float(t[2:])
					self.u  = self.UT
					self.UF = self.UT
				else :
					self.U = float(t[1:])
					self.u  = self.U
					self.UF = self.U
			if( t.find("S") > -1 ) : 
				self.S = float(t[1:])
			if( t.find("J") > -1 ) : 
				if( t.find("JT") > -1 ) : 
					self.JT  = float(t[2:])
					self.JTD = self.JT
				else :
					self.J = float(t[1:])
					self.JH = float(t[1:])
					self.J_H = float(t[1:])
			if( t.find("D") > -1 ) : 
				if( t.find("Dir") > -1 ) :  pass
				else : 
					self.D  = float(t[1:])
			if( t.find("nb") > -1 ) : 
				self.nb = int(t[t.find("nb")+2:])
				self.Nb = self.nb
				self.tNb = 2*self.Nb
				self.nc = int(t.split("nb")[0][t.find("nc")+2:])
				self.Nc = self.nc
				self.NC = self.nc
				self.tNC = 2*self.Nc
			if( t.find("Nmax") > -1 ) : 
				self.Nmax = int(t[4:])
			if( t.find("beta") > -1 ) : 
				self.beta = int(t[4:])
				self.Tfict = 11605./self.beta
			if( t.find("Hz") > -1 ) : 
				self.Hz = float(t[2:])
			if( t.find("IT") > -1 ) : 
				try : 
					self.IT = int(t[2:])
				except : 
					pass
		if self.IT<1 :	self.Treal = 0.
		else :		self.Treal = 11605./self.IT
		try :
			self.Tfit = 11605./self.beta
		except :
		 	self.Tfit = '.'

		try :
			if self.Hz  :  pass
		except :
		 	self.Hz = None 
		try : 
			if abs(self.JT)<1e-5 : 
				self.title = r"U={},J={},S={:.2f},D={};$\beta^*={}$,$\beta={}$;({})".format( self.UF , self.J  , self.S  , self.D  , self.beta, self.IT, self.tdirind ) 
			else :
				self.title = r"U={},J={},S={:.2f},D={},JT={};$\beta^*={}$,$\beta={}$;({})".format( self.UF , self.J  , self.S  , self.D , self.JT , self.beta, self.IT, self.tdirind ) 
		except :
			self.title = ddir.split('/')[-1]
		self.titlefoot  = ddir.split("/")[-1]
		try :
			with open( "local.c" ) as floc :
				lines  = floc.readlines()
				for  ll in range(45) :
					if (lines[ll].find("Nplot") > - 1) and (lines[ll].find("Nplot_") < 0) : 
						self.Nplot = int( lines[ll].split()[2] ) 
						break 
		except : 
			self.Nplot = 1025 
		try : 
			fpar = open( ddir + "/parameters_u{:.2f}.dat".format(self.UF) )
			parlines = fpar.readlines()
			self.parind = parlines[0]
			self.parindlist = self.parind.split()
			self.niter =  len( parlines )
			if( len ( parlines[-1] ) > 0 ) :
				self.count	= int( parlines[-1].split()[0] )
				self.countstart	= int( parlines[1].split()[0] )
				self.ncount	= self.count - self.countstart + 1
			else :
				print "The parameters-file reading is incomplete."
		except :
			print "The parameters-file reading is incomplete."
			pass
		try : 
			self.readParameters()
		except :
			pass
		self.readStatic()
	def headPrint(self ) :
		print "hP UF : ", self.UF
	def headSetTitle(self , ax ) :
		try : 
			ax.set_title( self.title+",n={:.3f}".format(self.tfilling_0) ) 
		except : 
			ax.set_title( self.title ) 
	def headPltTitle(self , plt) :
		plt.title( self.title ) 
	def maknameonly(self , lab ) : #to save
		if float(self.JT) == 0 :
			return "/next_nb{}_U{}_S{}_J{}_D{}_beta{}_IT{}_{}d{}gptl{}{}.mak".format(self.Nb, self.UF, self.S, self.J, self.D, self.beta, self.IT, self.tdirind, self.degen_0, self.gptl_0,  lab ) 
		else :
			return "/next_nb{}_JT{}_U{}_S{}_J{}_D{}_beta{}_IT{}_{}d{}gptl{}{}.mak".format(self.Nb, self.JT, self.UF, self.S, self.J, self.D, self.beta, self.IT, self.tdirind, self.degen_0, self.gptl_0,  lab ) 
	def makname(self , lab ) : #to save
		return self.tdir + self.maknameonly(lab)
	def fignamepart(self ) :
		try : 
			if float(self.JT) == 0 :
				return "plot_nb{}_U{}_S{}_J{}_D{}_beta{}_IT{}_{}d{}gptl{}".format(self.Nb, self.UF, self.S, self.J, self.D, self.beta, self.IT, self.tdirind, self.degen_0, self.gptl_0  )
			else : 
				return "plot_nb{}_JT{}_U{}_S{}_J{}_D{}_beta{}_IT{}_{}d{}gptl{}".format(self.Nb, self.JT, self.UF, self.S, self.J, self.D, self.beta, self.IT, self.tdirind , self.degen_0, self.gptl_0)
		except :
			return "plot_nb{}_JT{}_U{}_S{}_J{}_D{}_beta{}_IT{}_{}".format(self.Nb, self.JT, self.UF, self.S, self.J, self.D, self.beta, self.IT, self.tdirind)
	def readInitialMak(self ) :
		try :
		        foutput = open( self.tdir + "/output_u{:.2f}.txt".format( self.UF )  )
		        foutputlines = foutput.readlines()
		        if len( foutputlines ) > 0 :
		                self.initmakfull =  foutputlines[0].split()[-1]
				self.initmak =  self.initmakfull.split("/")[-1]
		                #print self.iarg, "\t", foutputlines[0].split()[-1]
			else : 
				self.initmak 	 =  "."
				self.initmakfull =  "."
		        foutput.close()
		except :
		        print self.iarg, "\tNone"
			self.initmak =  "None"
			self.initmakfull =  "None"
	def printInitialMak(self ) :
		print self.iarg, "\t", self.initmakfull
	def readParameters(self ) :
		try :
		        fpars = open( self.tdir + "/parameters_u{:.2f}.dat".format( self.UF )  )
		except :
			fpars = None 
			pass
		if fpars : 
			fparslines = fpars.readlines()
                        parindlist = fparslines[0].split()
                        nparindimp = 0
                        niind = 0
                        for parind in parindlist :
                                if parind.find("_0") > -1 :
                                        nparindimp += 1
                        nimax = 30
                        for parind in parindlist :
                                for ni in range(nimax) :
                                        if parind.find("_{}".format(ni))>-1 :
                                                niind += 1
                        #print "niind : ", niind
                        #print "nparind : ", nparindimp
                        niind /= nparindimp
			self.ni = niind
			self.Ni = niind
			self.nparindimp = nparindimp

                        pardatlist = fparslines[-1].split("\n")[0].split() 
			for kk in range(len(parindlist)) : 
				setattr( self , parindlist[kk], pardatlist[kk] ) 

			tfilling_all = 0.
			for parind in  parindlist : 
				try : 
					for item in ["filling", "flatt"] : 
						for j in range( self.ni ) : 
							if parind == (item+"_%d"%j) : 
								setattr( self, "c"+parind, float(getattr(self,parind))*self.nc ) 
				except : pass
				try : 
					for item in ["flattu"] : 
						for j in range( self.ni ) : 
							if parind == (item+"_%d"%j) : 
								setattr( self, "c"+parind, float(getattr(self,parind))*self.nc ) 
				except : pass 
				try : 
					for item in ["tfilling"] : 
						for j in range( self.ni ) : 
							if parind == (item+"_%d"%j) : 
								setattr( self, parind[1:], float(getattr(self,parind))/self.nc ) 
				except : pass 
				try : 
					for item in ["tfilling"] : 
						for j in range( self.ni ) : 
							if parind.find(item+"_%d"%j)>-1 : 
								tfilling_all += float(getattr(self,parind))
				except : pass 
			setattr( self, item+"_all", tfilling_all )
			setattr( self, item+"_avg", tfilling_all/self.ni )
			fpars.close()
		try : 
			degen = int(self.degen_0)
			setattr( self, "degen" , int(self.degen_0) ) 
		except : 
			pass
	def printParametersInd(self , indlist ) :
		itemlist =  indlist
		try : 
			ni=1
			for item in itemlist : 
				dat = item
                	        if item.find("chi2")>-1 :
                	                if len(dat)<16 :
                	                        bl = 16 - len(dat)
                	                        dat = dat + " "*bl
                	                        print dat, "\t",
                	                else :
                	                        print dat, "\t",
					ni+=1
                	        else :
                	                print dat, "\t",
                	        if item.find("chi2_{}".format(ni))>-1 :
					print "\t", ni-1, "\n"
			print "\t", "iargs"
		except :
			print "None"
			pass
	def printParametersIndOut(self , indlist , fileout ) :
		itemlist =  indlist
		try : 
			ni=1
			for item in itemlist : 
				dat = item
                	        if item.find("chi2")>-1 :
                	                if len(dat)<16 :
                	                        bl = 16 - len(dat)
                	                        dat = dat + " "*bl
					ni+=1
				fileout.write( "{}\t".format(dat) )
                	        if item.find("chi2_{}".format(ni))>-1 :
					fileout.write( "\t{}\n".format(ni-1) )
			fileout.write( "\tiargs\n" )
		except :
			fileout.write( "None\n" )
			pass
	def printParameters(self , itemlist) :
		try : 
			ni=0
			for item in itemlist : 
                                dat = getattr( self, item )
                                if item.find("chi2")>-1 and ni > 0 :
                                        print "\n{:5s}\t".format( "" ),
                                if item.find("chi2")>-1 :
                                        ni+=1
                                        if len(dat)<15 :
                                                bl = 15 - len(dat)
                                                dat = dat + " "*bl
                                                print dat, "\t",
                                        else :
                                                print dat, "\t",
                                elif item.find("gm")>-1 :
                                        print "{:>16}\t".format(dat),
                                elif item.find("count")>-1 :
                                        print "{:5s}\t".format( getattr(self,item) ),
                                elif item.find("fill")>-1 :
                                        print "{:10.8f}\t".format(float(dat)),
                                elif item.find("mag")>-1 :
                                        print "{:19.16f}\t".format(float(dat)),
                                else :
                                        print dat, "\t",
			print "\t", self.iarg
		except :
			print "None", self.iarg
			pass
	def printParametersOut(self , itemlist, fileout) :
		try : 
			ni=0
			for item in itemlist : 
                                dat = getattr( self, item )
                                if item.find("chi2")>-1 and ni > 0 :
                                        fileout.write( "\n{:5s}\t".format( "" ) )
                                if item.find("chi2")>-1 :
                                        ni+=1
                                        if len(dat)<15 :
                                                bl = 15 - len(dat)
                                                dat = dat + " "*bl
					fileout.write( "{}\t".format(dat) )
                                elif item.find("gm")>-1 :
                                        fileout.write( "{:>16}\t".format(dat) )
                                elif item.find("count")>-1 :
                                        fileout.write( "{:5s}\t".format(getattr(self,item)) )
                                elif item.find("fill")>-1 :
                                        fileout.write( "{:10.8f}\t".format(float(dat)) ) 
                                elif item.find("mag")>-1 :
                                        fileout.write( "{:19.16f}\t".format(float(dat)) )
                                else :
					fileout.write( "{}\t".format(dat) )
			fileout.write( "\t{}\n".format(self.iarg) )
		except :
			fileout.write( "None {}\n".format(self.iarg) )
			pass
	def readStatic(self ) :
		try : 
			nimp = int(self.Ni)
		except : 
			nimp = 0
		self.tldosind    = ["xyu","xyd","yzu","yzd","zxu","zxd"]
		indOrb = ["xy","yz","zx" ]
		try : 
			ffill = open( self.tdir + "/result/u{:.3f}/filling.dat".format( self.UF )  )	# 0.total 1.total_imp0 2-7. 8.total_imp1 9-14. 15.total_imp2 ...
			ffilllines = ffill.readlines()
			self.filling = ffilllines[-1].split("\n")[0].split()[1]
			self.jind = ["q3u","q1u","q1d","q3d","d1u","d1d","tot"]
			indOrb = ["q3","q1","d1" ]
			basisind	= self.jind
			if self.basis	== 't' : 
				basisind	= self.tldosind
				indOrb 		= ["xy","yz","zx" ]
			fillitem = ffilllines[-1].split("\n")[0].split()
			for nj in range( self.tNC ) :
				setattr( self, "fillimp"+basisind[nj] , float(fillitem[3+nj]) ) 
				setattr( self, "fillimp0"+basisind[nj] , float(fillitem[3+nj]) ) 
			for nj in range( self.NC ) :
				setattr( self, "fillimp"+basisind[2*nj][:2] , float(fillitem[3+2*nj]) + float(fillitem[3+2*nj+1]) ) 
				setattr( self, "fillimp0"+basisind[2*nj][:2] , float(fillitem[3+2*nj]) + float(fillitem[3+2*nj+1]) ) 
			if nimp>1 :
				for iimp in range( nimp ) :
					idat	 = 2 +  iimp*(self.tNC+1) 
					setattr( self, "fillimp%d"%iimp , float(fillitem[idat]) ) 
					for nj in range( self.NC ) :
						setattr( self, "fillimp%d"%iimp+indOrb[nj] , 0. ) 
					for nj in range( self.tNC ) :
						idat	 = 2 +  iimp*(self.tNC+1) + 1 + nj
						setattr( self, "fillimp%d"%iimp+basisind[nj] , float(fillitem[idat]) ) 
						for mj in range( self.NC ) :
							if basisind[nj].find(indOrb[mj])>-1 : 
								dum = getattr( self, "fillimp%d"%iimp+indOrb[mj] )
								setattr( self, "fillimp%d"%iimp+indOrb[mj] , dum+float(fillitem[idat]) ) 
			setattr( self, "fillimptot",    float( self.filling )    ) 
			setattr( self, "fillimptotavg", float( self.filling )/6. ) 
			ffill.close()
		except : 
			self.filling = "."
		try : 
			fmag  = open( self.tdir + "/result/u{:.3f}/mag.dat".format( self.UF )  )
			fmaglines = fmag.readlines()
			self.mag  = fmaglines[-1].split("\n")[0].split()[1]
			fmag.close()
		except : 
			self.mag = "."
		#try : 
		#	fmag  = open( self.tdir + "/result/u{:.3f}/double.dat".format( self.UF )  )
		#	fmaglines = fmag.readlines()
		#	self.doccu  = fmaglines[-1].split("\n")[0].split()
		#	for i in range(nbasis/2) :
		#		setattr( self, "doccu{}".format(i),    float( self.doccu[3+i] )    ) 
		#	fmag.close()
		#except : 
		#	self.doccu = "."
		#try : 
		#	fmag  = open( self.tdir + "/result/u{:.3f}/double_t2g.dat".format( self.UF )  )
		#	fmaglines = fmag.readlines()
		#	self.doccut2g  = fmaglines[-1].split("\n")[0].split()
		#	for i in range(nbasis/2) :
		#		setattr( self, "doccut2g{}".format(i),    float( self.doccut2g[3+i] )    ) 
		#	fmag.close()
		#except : 
		#	self.doccut2g = "."
		try : 
			fmag  = open( self.tdir + "/result/u{:.3f}/t2g_filling.dat".format( self.UF )  )
			fmaglines = fmag.readlines()
			self.fillt2g  = fmaglines[-1].split("\n")[0].split()
			for i in range(nbasis) :
				setattr( self, "fillt2g{}".format(i),    float( self.fillt2g[3+i] )    ) 
			if nimp>1 :
				for iimp in range( nimp ) :
					for nj in range( self.tNC ) :
						idat	 = 2 + iimp*(self.tNC+1) + 1 + nj
						setattr( self, "fillt2g%d"%iimp+self.tldosind[nj] , float(self.fillt2g[idat]) ) 
			fmag.close()
		except : 
			self.fillt2g = "."
		try : 
			fmag  = open( self.tdir + "/result/u{:.3f}/filling_t2g.dat".format( self.UF )  )
			fmaglines = fmag.readlines()
			self.fillt2g  = fmaglines[-1].split("\n")[0].split()
			for i in range(nbasis) :
				setattr( self, "fillt2g{}".format(i),    float( self.fillt2g[3+i] )    ) 
			fmag.close()
		except : 
			self.fillt2g = "."
		try : 
			fmag  = open( self.tdir + "/result/u{:.3f}/filling.dat".format( self.UF )  )
			fmaglines = fmag.readlines()
			self.fill  = fmaglines[-1].split("\n")[0].split()
			for i in range(nbasis) :
				setattr( self, "fill{}".format(i),    float( self.fill[3+i] )    ) 
			fmag.close()
		except : 
			self.fill = "."
		try : 
			ffill = self.tdir + "/result/u{:.3f}/fillingdeg.dat".format( self.UF )
			self.filldeg = np.genfromtxt( ffill, dtype=float ) 
		except : 
			self.filldeg = ["."]

		
		#try :
		#	ffexc = self.tdir + "/result/u{:.3f}/excitonic.dat".format( self.UF )
		#	self.exc = np.genfromtxt( ffexc, dtype=float ) 
		#	ffshape = np.shape( self.exc ) 
		#	if len(ffshape) > 1 :	# multi-line data
		#		lastdat = self.exc[-1][1:]
		#	else : 			# sinlgle-line data
		#		lastdat = self.exc[1:]
		#	exc = complexMatDatReturn( lastdat, nbasis ) 
		#	for ii in range( nbasis ) :
		#		for jj in range( nbasis ) :
		#			setattr( self, "exc{}{}re".format(ii,jj) , exc[ii][jj].real ) 
		#			setattr( self, "exc{}{}im".format(ii,jj) , exc[ii][jj].imag ) 
		#	for ii in range( len(exc[0]) ) :
		#		setattr( self, "fillimp"+self.tldosind[ii] , exc[ii][ii].real ) 
		#except :
		#	self.exc = '.'

		#try :
		#	ffexct2g = self.tdir + "/result/u{:.3f}/t2g_excitonic.dat".format( self.UF )
		#	self.exct2g = np.genfromtxt( ffexct2g, dtype=float ) 
		#	ffshape = np.shape( self.exct2g ) 
		#	if len(ffshape) > 1 :	# multi-line data
		#		lastdat = self.exct2g[-1][1:]
		#	else : 			# sinlgle-line data
		#		lastdat = self.exct2g[1:]
		#	exct2g = complexMatDatReturn( lastdat, nbasis ) 
		#	for ii in range( nbasis ) :
		#		for jj in range( nbasis ) :
		#			setattr( self, "exct2g{}{}re".format(ii,jj) , exct2g[ii][jj].real ) 
		#			setattr( self, "exct2g{}{}im".format(ii,jj) , exct2g[ii][jj].imag ) 
		#		for jj in range( nbasis ) :
		#			setattr( self, "exct2g{}{}reFrac".format(ii,jj) , exct2g[ii][jj].real / exct2g[ii][ii].real ) 
		#			setattr( self, "exct2g{}{}imFrac".format(ii,jj) , exct2g[ii][jj].imag / exct2g[ii][ii].real ) 
		#	for ii in range( len(exct2g[0]) ) :
		#		setattr( self, "fillimp"+self.tldosind[ii] , exct2g[ii][ii].real ) 
		#except :
		#	self.exct2g = '.'

		
		try :
			fraw  = self.tdir + "/result/u{:.3f}/EEall.dat".format( self.UF )
			frawdat = np.genfromtxt( fraw )
			ffshape	= np.shape(frawdat)
			if len(ffshape) > 1 :	# multi-line data
				lastdat = frawdat[-1][1:]
			else : 			# sinlgle-line data
				lastdat = frawdat[1:]
			nameop	= "EE"
			#print "lastdat : ", lastdat 
			eetot	= lastdat[0]
			setattr( self, nameop+"tot", lastdat[0] )
			namearrlist = [ nameop		, [ 'orb1', 'orb2' ] ]
			eeorb1 = np.zeros( self.Nc )
			eeorb2 = np.zeros( self.Nc )
			for imu in range(self.Nc) :
				orb1tail	= "_%d"%imu
				orb2tail	= ""
				for inu in range(self.Nc) :
					if imu!=inu :
						orb2tail	+= "_%d"%inu
				#print "orb1tail %d: "%imu , orb1tail	
				#print "orb2tail %d: "%imu , orb2tail	
				setattr( self, namearrlist[0] + namearrlist[1][0] + orb1tail, lastdat[1+imu*2  ] )
				setattr( self, namearrlist[0] + namearrlist[1][1] + orb2tail, lastdat[1+imu*2+1] )
				eeorb1[imu]	= lastdat[1+imu*2  ]
				eeorb2[imu]	= lastdat[1+imu*2+1]
			nameop	= "MInfo"
			namearrlist = [ nameop		, [ 'orb1', 'orb2' ] ]
			#print "eeorb1 : ", eeorb1
			#print "eeorb2 : ", eeorb2
			for imu in range(self.Nc) :
				orb1tail	= "_%d"%imu
				orb2tail	= ""
				MInfo1	= eeorb1[imu] + eeorb2[imu] - eetot
				MInfo2	= -eetot
				Minfotail = ""
				for inu in range(self.Nc) :
					if imu!=inu :
						orb2tail	+= "_%d"%inu
						MInfo2		+= eeorb2[inu]
						Minfotail	+= "+orb2.%d"%inu
					else : 
						MInfo2		-= eeorb1[inu]
						Minfotail	+= "-orb1.%d"%inu
				#print "Minfotail %d: "%imu, Minfotail
				#print "orb1tail %d: "%imu , orb1tail	
				#print "orb2tail %d: "%imu , orb2tail	
				setattr( self, namearrlist[0] + namearrlist[1][0] + orb1tail, MInfo1 )
				setattr( self, namearrlist[0] + namearrlist[1][1] + orb2tail, MInfo2 )
		except : 
			pass
		

		try :
			fjsq = self.tdir + "/result/u{:.3f}/jsq.dat".format( self.UF )
			self.jsqdat = np.genfromtxt( fjsq, dtype=float ) 
			ffshape = np.shape( self.jsqdat ) 
			if len(ffshape) > 1 :	# multi-line data
				lastdat = self.jsqdat[-1][1:]
			else : 			# sinlgle-line data
				lastdat = self.jsqdat[1:]
			self.jsq = float(lastdat[0]) + float(lastdat[1])*1j
		except :
			self.jsq = '.'

		try :
			fjsq = self.tdir + "/result/u{:.3f}/jsqdeg.dat".format( self.UF )
			self.jsqdegdat = np.genfromtxt( fjsq, dtype=float ) 
			ffshape = np.shape( self.jsqdegdat ) 
		except :
			self.jsqdegdat = '.'

		namearrlist = []
		namearrlist.append(  [ "JJdeg"			, [ 'ic', 'sq', 'zz', 'xyxy' ]	] )
		#namearrlist.append(  [ "JaJadeg"		, [ 'ic', 'sq', 'zz', 'xyxy' ]	] )
		namearrlist.append(  [ "LLdeg"			, [ 'ic', 'sq', 'zz', 'xyxy' ] 	] )
		namearrlist.append(  [ "SSdeg"			, [ 'ic', 'sq', 'zz', 'xyxy' ] 	] )
		namearrlist.append(  [ "MMdeg"			, [ 'ic', 'sq', 'zz', 'xyxy' ]	] )
		namearrlist.append(  [ "Jzxydeg"		, [ 'ic', 'z', 'x', 'y' ]	] )
		#namearrlist.append(  [ "Jazxydeg"		, [ 'ic', 'z', 'x', 'y' ]	] )
		namearrlist.append(  [ "Lzxydeg"		, [ 'ic', 'z', 'x', 'y' ]	] )
		namearrlist.append(  [ "Szxydeg"		, [ 'ic', 'z', 'x', 'y' ]	] )
		namearrlist.append(  [ "Mzxydeg"		, [ 'ic', 'z', 'x', 'y' ]	] )
		for namearr in namearrlist : 
			setattrwithlist( self, namearr[0], namearr[1] )
		namearrlist = []
		namearrlist.append(  [ "JzdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "JxdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "JydegMateval"		, [ 'ic', '' ]			] )
		#namearrlist.append(  [ "JazdegMateval"		, [ 'ic', '' ]			] )
		#namearrlist.append(  [ "JaxdegMateval"		, [ 'ic', '' ]			] )
		#namearrlist.append(  [ "JaydegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "LzdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "LxdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "LydegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "SzdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "SxdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "SydegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "MzdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "MxdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "MydegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "JJdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "LLdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "SSdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "MMdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "JJzzdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "LLzzdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "SSzzdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "MMzzdegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "JJxydegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "LLxydegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "SSxydegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "MMxydegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QSx2y2Mateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QLx2y2Mateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QJx2y2Mateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QMx2y2Mateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QSz2r2Mateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QLz2r2Mateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QJz2r2Mateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QMz2r2Mateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QSxyMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QLxyMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QJxyMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QMxyMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QSyzMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QLyzMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QJyzMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QMyzMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QSzxMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QLzxMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QJzxMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "QMzxMateval"		, [ 'ic', '' ]			] )
		for namearr in namearrlist : 
			setattrwithlistdeg( self, namearr[0], namearr[1] , real=True )
		namearrlist = []
		namearrlist.append(  [ "Jzxydeg"		, [ 'ic', 'z', 'x', 'y' ] ] )
		#namearrlist.append(  [ "Jazxydeg"		, [ 'ic', 'z', 'x', 'y' ] ] )
		namearrlist.append(  [ "Lzxydeg"		, [ 'ic', 'z', 'x', 'y' ] ] )
		namearrlist.append(  [ "Szxydeg"		, [ 'ic', 'z', 'x', 'y' ] ] )
		namearrlist.append(  [ "Mzxydeg"		, [ 'ic', 'z', 'x', 'y' ] ] )
		for namearr in namearrlist : 
			setattrTranscomponent( self, namearr[0], namearr[1] )
		namearrlist = []
		namearrlist.append(  [ "JJxydegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "LLxydegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "SSxydegMateval"		, [ 'ic', '' ]			] )
		namearrlist.append(  [ "MMxydegMateval"		, [ 'ic', '' ]			] )
		for namearr in namearrlist : 
			try : 
				dumval = float(getattr( self, namearr[0] ))
				setattr( self, namearr[0]+"half", dumval*0.5 ) 
			except :
				setattr( self, namearr[0]+"half", '.' ) 
		namearrlist = []
		namearrlist.append(  [ "JzdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "JxdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "JydegMateval"		, [ 'ic', 'Tr' ]			] )
		#namearrlist.append(  [ "JazdegMateval"		, [ 'ic', 'Tr' ]			] )
		#namearrlist.append(  [ "JaxdegMateval"		, [ 'ic', 'Tr' ]			] )
		#namearrlist.append(  [ "JaydegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "LzdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "LxdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "LydegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "SzdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "SxdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "SydegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "MzdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "MxdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "MydegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "JJdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "LLdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "SSdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "MMdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "JJzzdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "LLzzdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "SSzzdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "MMzzdegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "JJxydegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "LLxydegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "SSxydegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "MMxydegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QSx2y2Mateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QLx2y2Mateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QJx2y2Mateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QMx2y2Mateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QSz2r2Mateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QLz2r2Mateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QJz2r2Mateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QMz2r2Mateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QSxyMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QLxyMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QJxyMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QMxyMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QSyzMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QLyzMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QJyzMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QMyzMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QSzxMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QLzxMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QJzxMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "QMzxMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "JJxydegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "LLxydegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "SSxydegMateval"		, [ 'ic', 'Tr' ]			] )
		namearrlist.append(  [ "MMxydegMateval"		, [ 'ic', 'Tr' ]			] )
		for namearr in namearrlist : 
			setattrwithlistdeg( self, namearr[0], namearr[1] , real=True , iftrace=True)
			try : 
				dumval = float(getattr( self, namearr[0] ))
				#setattr( self, namearr[0][:-4]+"Tr", dumval.5 ) 
			except :
				pass
				#setattr( self, namearr[0]+"half", '.' ) 
		namearrlist = []
		namearrlist.append(  [ "JzdegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		namearrlist.append(  [ "JxdegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		namearrlist.append(  [ "JydegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		#namearrlist.append(  [ "JazdegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		#namearrlist.append(  [ "JaxdegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		#namearrlist.append(  [ "JaydegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		namearrlist.append(  [ "LzdegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		namearrlist.append(  [ "LxdegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		namearrlist.append(  [ "LydegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		namearrlist.append(  [ "SzdegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		namearrlist.append(  [ "SxdegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		namearrlist.append(  [ "SydegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		namearrlist.append(  [ "MzdegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		namearrlist.append(  [ "MxdegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		namearrlist.append(  [ "MydegMatLtcorrT"		, [ 'ic', '', 'Diag', 'Split' ]			] )
		for namearr in namearrlist : 
			try : 
				setattrwithlistdeg( self, namearr[0], namearr[1] , real=True , iftrace=False)
				setattr( self, namearr[0]+"DiagBeta",		self.IT*getattr(self,namearr[0]+"Diag")  )
				setattr( self, namearr[0]+"Beta",		self.IT*getattr(self,namearr[0])  )
				setattr( self, namearr[0]+"DiagBetafit",	self.beta*getattr(self,namearr[0]+"Diag")  )
				setattr( self, namearr[0]+"Betafit",		self.beta*getattr(self,namearr[0])  )
				dumval = float(getattr( self, namearr[0] ))
				#setattr( self, namearr[0][:-4]+"Tr", dumval.5 ) 
			except :
				pass
				#setattr( self, namearr[0]+"half", '.' ) 
		namearrlist = []
		namearrlist.append(  [ "MzdegMat"		, [ 'ic', '0','1' ]			] )
		namearrlist.append(  [ "MxdegMat"		, [ 'ic', '0','1' ]			] )
		for namearr in namearrlist : 
			setattrwithlistdeg( self, namearr[0], namearr[1] , real=False , iftrace=False, diag=True )
			try : 
				dumval = float(getattr( self, namearr[0] ))
				#setattr( self, namearr[0][:-4]+"Tr", dumval.5 ) 
			except :
				pass
				#setattr( self, namearr[0]+"half", '.' ) 
		namearrlist = []
		namearrlist.append(  [ "JzdegMatevalT"		, [ 'ic', 'TrAbs' ]			] )
		for namearr in namearrlist : 
			setattrwithlistdeg( self, namearr[0], namearr[1] , real=True , iftrace=True, ifabs=True )
			try : 
				dumval = float(getattr( self, namearr[0] ))
				#setattr( self, namearr[0], dumval*0.5 ) 
			except :
				setattr( self, namearr[0], '.' ) 
			setattrwithlistdeg( self, namearr[0], ['ic','Tr'] , real=True , iftrace=True, ifabs=False )

		namearrlist = []
		namearrlist.append(  [ "JzdegMateval"		, [ 'Sq' ]			] )
		namearrlist.append(  [ "JxdegMateval"		, [ 'Sq' ]			] )
		namearrlist.append(  [ "JydegMateval"		, [ 'Sq' ]			] )
		namearrlist.append(  [ "LzdegMateval"		, [ 'Sq' ]			] )
		namearrlist.append(  [ "LxdegMateval"		, [ 'Sq' ]			] )
		namearrlist.append(  [ "LydegMateval"		, [ 'Sq' ]			] )
		namearrlist.append(  [ "SzdegMateval"		, [ 'Sq' ]			] )
		namearrlist.append(  [ "SxdegMateval"		, [ 'Sq' ]			] )
		namearrlist.append(  [ "SydegMateval"		, [ 'Sq' ]			] )
		namearrlist.append(  [ "MzdegMateval"		, [ 'Sq' ]			] )
		namearrlist.append(  [ "MxdegMateval"		, [ 'Sq' ]			] )
		namearrlist.append(  [ "MydegMateval"		, [ 'Sq' ]			] )
		for namearr in namearrlist : 
			try : 
				dumval = float(getattr( self, namearr[0] ))
				setattr( self, namearr[0]+namearr[1][0], dumval*dumval ) 
			except :
				setattr( self, namearr[0]+namearr[1][0], '.' ) 
				pass
				#setattr( self, namearr[0]+"half", '.' ) 

		try :
			oplist		= ['S','L','J','M']
			for op in oplist : 
				fname		= self.tdir + "/result/u{:.3f}/{}vec.dat".format( self.UF, op )
				dumdat		= np.genfromtxt( fname, dtype=float ) 
				dumdatshape	= np.shape( dumdat ) 
				ioffset		= 1
				if len(dumdatshape) > 1 :	# multi-line data
					lastdat = dumdat[-1][ioffset:]
				else : 			# sinlgle-line data
					lastdat = dumdat[ioffset:]
				namearrlist = [ op		, [ 'x', 'y', 'z' ] ]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], lastdat[ii] )
				for iimp in range(nimp) : 
					for ii in range(len(namearrlist[1])) : 
						setattr( self, namearrlist[0]+namearrlist[1][ii]+'%d'%iimp, lastdat[ii + 3*iimp] )
		except :
			pass

		
		try : 
			RLocalaxCart = np.genfromtxt("inputs/TRANSFLOCALcart")
			if nimp>1 : 
				RLocalaxCart.reshape( nimp, 3, 3 ) 
			oplist		= ['S','L','J','M']
			for op in oplist : 
				#opvec	= [ getattr( self, op+mu ) for mu in ['x','y','z'] ]
				#print "OPVEC ", op, " : ", opvec
				for iimp in range(nimp) : 
					opvec	= [ getattr( self, op+mu+'%d'%iimp ) for mu in ['x','y','z'] ]
					#print "OPVEC ", op, iimp, " : ", opvec
					opvecglob	= np.dot( RLocalaxCart[iimp].transpose(), np.dot(np.array(opvec), RLocalaxCart[iimp]) )
					for imu in range(3) :
						mu = ['x','y','z'][imu]
						setattr( self, op+mu+'%d'%iimp+'glob' , opvecglob[imu] )
		except : 
			pass

		try : 
			fdatArr = [  self.tdir + "/lattice/occupancy_Nplot3072_wh{:g}_{}th.dat".format(wh,self.count) for wh in [10000, 100000, 1000000]  ] 
		except : 
			fdatArr = []
		for fdat in fdatArr : 
			try : 
				#fdat = self.tdir + "/lattice/occupancy_Nplot3072_wh10000_{}th.dat".format(self.count)
				fgoccudat = np.genfromtxt( fdat, dtype=float ) 
				dimdat = np.shape( fgoccudat ) 
				#print "FDAT : ", fdat, dimdat
				try :
					if dimdat[1]>1 :
						self.nOcculattAraFormer = float(fgoccudat[-1][0])
						break
				except : 
					self.nOcculattAraFormer = float(fgoccudat[0])
					break
			except :
				self.nOcculattAraFormer = '.'
		
		try : 
			fdatArr = [  self.tdir + "/lattice/filling_Nplot3072_wh{:g}_{}th.dat".format(wh,self.count) for wh in [10000, 100000, 1000000]  ] 
		except : 
			fdatArr = []
		for fdat in fdatArr : 
			try : 
				fgoccudat	= np.genfromtxt( fdat, dtype=float ) 
				dimdat		= np.shape( fgoccudat ) 
				NU	= self.Ni * self.tNC
				for inu in range(NU) : 
					setattr( self, "nOcculattAraFormerOrb{}".format(inu), float(fgoccudat[inu,1]) )
				for ini in range(self.ni) : 
					setattr( self, "nOcculattAraFormerImp{}".format(ini), np.sum(fgoccudat[ini*self.ni:ini*self.ni+self.tNC,1]) )
			except :
				pass
		
		try : 
			fdat = self.tdir + "/lattice/occupancy_Nplot1024_wh10000_{}th.dat".format(self.count)
			fgoccudat = np.genfromtxt( fdat, dtype=float ) 
			dimdat = np.shape( fgoccudat ) 
			try :
				if dimdat[1]>1 :
					self.nOcculattAra = float(fgoccudat[-1][0])
			except : 
				self.nOcculattAra = float(fgoccudat[0])
		except :
			self.nOcculattAra = '.'
		
		try : 
			fdat = self.tdir + "/result/u{:.3f}/EEdeg.dat".format( self.UF )
			rawdat = np.genfromtxt( fdat, dtype=float ) 
			namearrlist = [ "EEdeg"		, [ 'ic', '' ] ]
			setattrwithlistdeg( self, namearrlist[0], namearrlist[1] , real=True )
		except :
			self.EEdeg = '.'

		try :	
			fdat = self.tdir + "/renorm_im_j.dat"
			self.massimj = np.genfromtxt( fdat, dtype=float ) 
			namearrlist = [ "massimj"		, [ 'q3u', 'q1u', 'q1d', 'q3d', 'd1u', 'd1d' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.massimj[ii] )
				setattr( self, 'Zimj'+namearrlist[1][ii], 1./self.massimj[ii] )
		except :
			self.massimj = '.'

		try :	
			fdat = self.tdir + "/result/u{:.3f}/mass_fitIm_eg.dat".format( self.UF )
			self.massfitimt = np.genfromtxt( fdat, dtype=float ) 
			ioffset = 3
			namearrlist = [ "massfitimt"		, [ 'xyu', 'yzu' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimt[ii+ioffset] )
				setattr( self, 'Zfitimt'+namearrlist[1][ii], 1./self.massfitimt[ii+ioffset] )
		except :
			self.massfitimt = '.'

		try :	
			fdat = self.tdir + "/result/u{:.3f}/mass_fitIm_t2g.dat".format( self.UF )
			self.massfitimt = np.genfromtxt( fdat, dtype=float ) 
			ioffset = 4
			namearrlist = [ "massfitimt"		, [ 'xyu', 'yzu', 'zxu' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimt[ii+ioffset] )
				setattr( self, 'Zfitimt'+namearrlist[1][ii], 1./self.massfitimt[ii+ioffset] )
		except :
			self.massfitimt = '.'

		try :	
			fdat = self.tdir + "/result/u{:.3f}/mass_fitIm.dat".format( self.UF )
			self.massfitimj = np.genfromtxt( fdat, dtype=float ) 
			ioffset = 4
			namearrlist = [ "massfitimj"		, [ 'q3u', 'q1u', 'd1u' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimj[ii+ioffset] )
				setattr( self, 'Zfitimj'+namearrlist[1][ii], 1./self.massfitimj[ii+ioffset] )
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii]+"Ratq1", self.massfitimj[ii+ioffset] / self.massfitimj[1+ioffset] )
		except :
			self.massfitimj = '.'

		try :	
			fdat = self.tdir + "/result/u{:.3f}/poly2_fitSelfIm.dat".format( self.UF )
			self.massfitimt = np.genfromtxt( fdat, dtype=float ) 
			ioffset = 4
			namearrlist = [ "masspoly2fit"		, [ 'xyu', 'yzu', 'zxu' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimt[ii+ioffset] )
				setattr( self, 'Zpoly2fit'+namearrlist[1][ii], 1./self.massfitimt[ii+ioffset] )
		except :
			self.massfitimt = '.'

		try :	
			fdat = self.tdir + "/result/u{:.3f}/poly2_fitSelfIm".format( self.UF )
			self.massfitimj = np.genfromtxt( fdat, dtype=float ) 
			ioffset = 4
			namearrlist = [ "masspoly2fit"		, [ 'q3u', 'q1u', 'd1u' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimj[ii+ioffset] )
				setattr( self, 'Zpoly2fit'+namearrlist[1][ii], 1./self.massfitimj[ii+ioffset] )
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii]+"Ratq1", self.massfitimj[ii+ioffset] / self.massfitimj[1+ioffset] )
		except :
			self.massfitimj = '.'

		try :	
			fdat = self.tdir + "/result/u{:.3f}/poly3_fitSelfIm_t2g.dat".format( self.UF )
			self.massfitimt = np.genfromtxt( fdat, dtype=float ) 
			ioffset = 4
			namearrlist = [ "masspoly3fit"		, [ 'xyu', 'yzu', 'zxu' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimt[ii+ioffset] )
				setattr( self, 'Zpoly3fit'+namearrlist[1][ii], 1./self.massfitimt[ii+ioffset] )

				setattr( self, 'err'+namearrlist[0]+namearrlist[1][ii], self.massfitimt[ii+ioffset+3] )
				m	 =  self.massfitimt[ii+ioffset]
				errm	 =  self.massfitimt[ii+ioffset+3]
				errZ	 =  ( abs( 1./m - 1./(m+errm) )  + abs( 1./m - 1./(m-errm) )  )/2.
				setattr( self, 'errZpoly3fit'+namearrlist[1][ii], errZ )
		except :
			self.massfitimt = '.'

		try :	
			fdat = self.tdir + "/result/u{:.3f}/poly3_fitSelfIm.dat".format( self.UF )
			self.massfitimj = np.genfromtxt( fdat, dtype=float ) 
			ioffset = 4
			namearrlist = [ "masspoly3fit"		, [ 'q3u', 'q1u', 'd1u' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimj[ii+ioffset] )
				setattr( self, 'Zpoly3fit'+namearrlist[1][ii], 1./self.massfitimj[ii+ioffset] )
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii]+"Ratq1", self.massfitimj[ii+ioffset] / self.massfitimj[1+ioffset] )
		except :
			self.massfitimj = '.'

		try :	
			fdat = self.tdir + "/result/u{:.3f}/poly4_fitSelfIm_t2g.dat".format( self.UF )
			self.massfitimt = np.genfromtxt( fdat, dtype=float ) 
			ioffset = 4
			namearrlist = [ "masspoly4fit"		, [ 'xyu', 'yzu', 'zxu' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimt[ii+ioffset] )
				setattr( self, 'Zpoly4fit'+namearrlist[1][ii], 1./self.massfitimt[ii+ioffset] )

				setattr( self, 'err'+namearrlist[0]+namearrlist[1][ii], self.massfitimt[ii+ioffset+3] )
				m	 =  self.massfitimt[ii+ioffset]
				errm	 =  self.massfitimt[ii+ioffset+3]
				errZ	 =  ( abs( 1./m - 1./(m+errm) )  + abs( 1./m - 1./(m-errm) )  )/2.
				setattr( self, 'errZpoly4fit'+namearrlist[1][ii], errZ )
			ioffset = 7
			namearrlist = [ "poly4fitconst"		, [ 'xyu', 'yzu', 'zxu' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimt[ii+ioffset] )
		except :
			self.massfitimt = '.'

		try :	
			fdat = self.tdir + "/result/u{:.3f}/poly4_fitSelfIm.dat".format( self.UF )
			self.massfitimj = np.genfromtxt( fdat, dtype=float ) 
			ioffset = 4
			namearrlist = [ "masspoly4fit"		, [ 'q3u', 'q1u', 'd1u' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimj[ii+ioffset] )
				setattr( self, 'Zpoly4fit'+namearrlist[1][ii], 1./self.massfitimj[ii+ioffset] )
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii]+"Ratq1", self.massfitimj[ii+ioffset] / self.massfitimj[1+ioffset] )
			ioffset = 7
			namearrlist = [ "poly4fitconst"		, [ 'q3u', 'q1u', 'd1u' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimj[ii+ioffset] )
		except :
			self.massfitimj = '.'

		fheadArr = ['power', 'powerplusconst' ]
		for fhead in fheadArr : 
			try :	
				fdat = self.tdir + "/result/u{:.3f}/{}_fitSelfIm_t2g.dat".format( self.UF, fhead )
				self.powerfitt = np.genfromtxt( fdat, dtype=float ) 
				ioffset = 4
				namearrlist = [ fhead+"fitt"		, [ 'xyu', 'yzu', 'zxu' ] ]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.powerfitt[ii+ioffset] )
				ioffset = 7
				namearrlist = [ fhead+"fitconst"		, [ 'xyu', 'yzu', 'zxu' ] ]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.powerfitt[ii+ioffset] )
			except :
				self.powerfitt = '.'
	
			try :	
				fdat = self.tdir + "/result/u{:.3f}/{}_fitSelfIm_eg.dat".format( self.UF, fhead )
				self.powerfitt = np.genfromtxt( fdat, dtype=float ) 
				ioffset = 3
				namearrlist = [ fhead+"fitt"		, [ 'xyu', 'yzu' ] ]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.powerfitt[ii+ioffset] )
				ioffset = 5
				namearrlist = [ fhead+"fitconst"		, [ 'xyu', 'yzu' ] ]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.powerfitt[ii+ioffset] )
			except :
				self.powerfitt = '.'
	
			try :	
				fdat = self.tdir + "/result/u{:.3f}/{}_fitSelfIm.dat".format( self.UF, fhead )
				self.powerfitj = np.genfromtxt( fdat, dtype=float ) 
				ioffset = 4
				namearrlist = [ fhead+"fitj"		, [ 'q3u', 'q1u', 'd1u' ] ]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.powerfitj[ii+ioffset] )
				ioffset = 7
				namearrlist = [ fhead+"fitconst"		, [ 'q3u', 'q1u', 'd1u' ] ]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.powerfitj[ii+ioffset] )
			except :
				self.powerfitj = '.'
			#try :	
			#	namearrlist = [ fhead+"fitt"		, [ 'xyu', 'yzu', 'zxu' ] ]
			#	fillpowerfractot = 0.
			#	for item in namearrlist[1] :
			#		dumpower = getattr( self, "powerfitt"+item )
			#		dumfill  = getattr( self, "fillimp"  +item )
			#		setattr( self, "fill"+fhead+item , dumpower*dumfill )
			#		setattr( self, "fill"+fhead+"frac"+item , dumpower*dumfill/self.fillimptot*2 )
			#		fillpowerfractot += dumpower*dumfill/self.fillimptot*2 
			#	setattr( self, "fill"+fhead+"ractot" , fillpowerfractot )
	
			#	namearrlist = [ fhead+"fitj"		, [ 'q3u', 'q1u', 'd1u' ] ]
			#	fillpowerfractot = 0.
			#	for item in namearrlist[1] :
			#		dumpower = getattr( self, fhead+"fitj"+item )
			#		dumfill  = getattr( self, "fillimp"  +item )
			#		setattr( self, "fill"+fhead+item , dumpower*dumfill )
			#		setattr( self, "fill"+fhead+"frac"+item , dumpower*dumfill/self.fillimptot*2 )
			#		fillpowerfractot += dumpower*dumfill/self.fillimptot*2 
			#	setattr( self, "fill"+fhead+"fractot" , fillpowerfractot )
			#except :
			#	pass

		fileIndArr	= [ "" ]
		for ix	 in range(0,nbasis+1) :
			fileIndArr.append( "d%d"%(ix) )
		#print "filndIndArr : ", fileIndArr 
		for fileInd	in fileIndArr :
			try :	
				fdat = self.tdir + "/result/u{:.3f}/rdm{}diag.dat".format( self.UF, fileInd )
				rdmdxdiag = np.genfromtxt( fdat, dtype=float ) 
				ioffset = 1
				if( len(np.shape(rdmdxdiag)) == 1 ) :
					rdmdxdiagsum	= np.sum( rdmdxdiag[ioffset:] ) 
				else :
					ndeg		= np.shape(rdmdxdiag)[0]
					rdmdxdiagsum	= np.sum( rdmdxdiag[:,ioffset:] )  / ndeg
				nameVal	= "rdm{}diagsum".format( fileInd )
				setattr( self, nameVal, rdmdxdiagsum )
			except :
				nameVal	= "rdm{}diagsum".format( fileInd )
				setattr( self, nameVal, '.' )

		for iimp in range(nimp) : 
			try :	
				fdat = self.tdir + "/result/u{:.3f}/mass_fitIm_eg_{}.dat".format( self.UF , iimp )
				self.massfitimt = np.genfromtxt( fdat, dtype=float ) 
				ioffset = 3
				namearrlist  = [ "massfitimt{}".format(iimp)		, [ 'xyu', 'yzu' ] ]
				namearrlistZ = [ "Zfitimt{}".format(iimp)		]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimt[ii+ioffset] )
					setattr( self, namearrlistZ[0]+namearrlist[1][ii], 1./self.massfitimt[ii+ioffset] )
			except :
				self.massfitimt = '.'
	
			try :	
				fdat = self.tdir + "/result/u{:.3f}/mass_fitIm_t2g_{}.dat".format( self.UF , iimp )
				self.massfitimt = np.genfromtxt( fdat, dtype=float ) 
				ioffset = 4
				namearrlist  = [ "massfitimt{}".format(iimp)		, [ 'xyu', 'yzu', 'zxu' ] ]
				namearrlistZ = [ "Zfitimt{}".format(iimp)		]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimt[ii+ioffset] )
					setattr( self, namearrlistZ[0]+namearrlist[1][ii], 1./self.massfitimt[ii+ioffset] )
			except :
				self.massfitimt = '.'
	
			try :	
				fdat = self.tdir + "/result/u{:.3f}/mass_fitIm_{}.dat".format( self.UF , iimp )
				self.massfitimj = np.genfromtxt( fdat, dtype=float ) 
				ioffset = 4
				namearrlist  = [ "massfitimj{}".format(iimp)		, [ 'q3u', 'q1u', 'd1u' ] ]
				namearrlistZ = [ "Zfitimj{}".format(iimp)		]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimj[ii+ioffset] )
					setattr( self, namearrlistZ[0]+namearrlist[1][ii], 1./self.massfitimj[ii+ioffset] )
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii]+"Ratq1", self.massfitimj[ii+ioffset] / self.massfitimj[1+ioffset] )
			except :
				self.massfitimj = '.'
	
			try :	
				fdat = self.tdir + "/result/u{:.3f}/poly4_fitSelfIm_eg_{}.dat".format( self.UF, iimp )
				self.massfitimt = np.genfromtxt( fdat, dtype=float ) 
				ioffset = 3
				namearrlist = [ "masspoly4fit{}".format(iimp)		, [ 'xyu', 'yzu' ] ]
				namearrlistZ = [ "Zpoly4fit{}".format(iimp)		]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimt[ii+ioffset] )
					setattr( self, namearrlistZ[0]+namearrlist[1][ii], 1./self.massfitimt[ii+ioffset] )
			except :
				self.massfitimt = '.'

			try :	
				fdat = self.tdir + "/result/u{:.3f}/poly4_fitSelfIm_t2g_{}.dat".format( self.UF, iimp )
				self.massfitimt = np.genfromtxt( fdat, dtype=float ) 
				ioffset = 4
				namearrlist = [ "masspoly4fit{}".format(iimp)		, [ 'xyu', 'yzu', 'zxu' ] ]
				namearrlistZ = [ "Zpoly4fit{}".format(iimp)		]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimt[ii+ioffset] )
					setattr( self, namearrlistZ[0]+namearrlist[1][ii], 1./self.massfitimt[ii+ioffset] )
			except :
				self.massfitimt = '.'

			try :	
				fdat = self.tdir + "/result/u{:.3f}/poly4_fitSelfIm_{}.dat".format( self.UF, iimp )
				self.massfitimj = np.genfromtxt( fdat, dtype=float ) 
				ioffset = 4
				namearrlist = [ "masspoly4fit{}".format(iimp)		, [ 'q3u', 'q1u', 'd1u' ] ]
				namearrlistZ = [ "Zpoly4fit{}".format(iimp)		]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimj[ii+ioffset] )
					setattr( self, namearrlistZ[0]+namearrlist[1][ii], 1./self.massfitimj[ii+ioffset] )
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii]+"Ratq1", self.massfitimj[ii+ioffset] / self.massfitimj[1+ioffset] )
			except :
				self.massfitimj = '.'

			try :	
				fdat = self.tdir + "/result/u{:.3f}/poly3_fitSelfIm_t2g_{}.dat".format( self.UF, iimp )
				self.massfitimt = np.genfromtxt( fdat, dtype=float ) 
				ioffset = 4
				namearrlist = [ "masspoly3fit{}".format(iimp)		, [ 'xyu', 'yzu', 'zxu' ] ]
				namearrlistZ = [ "Zpoly3fit{}".format(iimp)		]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimt[ii+ioffset] )
					setattr( self, namearrlistZ[0]+namearrlist[1][ii], 1./self.massfitimt[ii+ioffset] )
			except :
				self.massfitimt = '.'

			try :	
				fdat = self.tdir + "/result/u{:.3f}/poly3_fitSelfIm_{}.dat".format( self.UF, iimp )
				self.massfitimj = np.genfromtxt( fdat, dtype=float ) 
				ioffset = 4
				namearrlist = [ "masspoly3fit{}".format(iimp)		, [ 'q3u', 'q1u', 'd1u' ] ]
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii], self.massfitimj[ii+ioffset] )
					setattr( self, 'Zpoly3fit'+namearrlist[1][ii], 1./self.massfitimj[ii+ioffset] )
				for ii in range(len(namearrlist[1])) : 
					setattr( self, namearrlist[0]+namearrlist[1][ii]+"Ratq1", self.massfitimj[ii+ioffset] / self.massfitimj[1+ioffset] )
			except :
				self.massfitimj = '.'

			fheadArr = ['power', 'powerplusconst' ]
			for fhead in fheadArr : 
				try :	
					fdat = self.tdir + "/result/u{:.3f}/{}_fitSelfIm_eg_{}.dat".format( self.UF , fhead, iimp )
					self.powerfitt = np.genfromtxt( fdat, dtype=float ) 
					ioffset = 3
					namearrlist = [ fhead+"fitt{}".format(iimp)		, [ 'xyu', 'yzu' ] ]
					for ii in range(len(namearrlist[1])) : 
						setattr( self, namearrlist[0]+namearrlist[1][ii], self.powerfitt[ii+ioffset] )
					ioffset = 5
					namearrlist = [ fhead+"fitconst{}".format(iimp)		, [ 'xyu', 'yzu' ] ]
					for ii in range(len(namearrlist[1])) : 
						setattr( self, namearrlist[0]+namearrlist[1][ii], self.powerfitt[ii+ioffset] )
				except :
					self.powerfitt = '.'
	
				try :	
					fdat = self.tdir + "/result/u{:.3f}/{}_fitSelfIm_t2g_{}.dat".format( self.UF , fhead, iimp )
					self.powerfitt = np.genfromtxt( fdat, dtype=float ) 
					ioffset = 4
					namearrlist = [ fhead+"fitt{}".format(iimp)		, [ 'xyu', 'yzu', 'zxu' ] ]
					for ii in range(len(namearrlist[1])) : 
						setattr( self, namearrlist[0]+namearrlist[1][ii], self.powerfitt[ii+ioffset] )
					ioffset = 7
					namearrlist = [ fhead+"fitconst{}".format(iimp)		, [ 'xyu', 'yzu', 'zxu' ] ]
					for ii in range(len(namearrlist[1])) : 
						setattr( self, namearrlist[0]+namearrlist[1][ii], self.powerfitt[ii+ioffset] )
				except :
					self.powerfitt = '.'
	
				try :	
					fdat = self.tdir + "/result/u{:.3f}/{}_fitSelfIm_{}.dat".format( self.UF , fhead, iimp )
					self.powerfitj = np.genfromtxt( fdat, dtype=float ) 
					ioffset = 4
					namearrlist = [ fhead+"fitj{}".format(iimp)		, [ 'q3u', 'q1u', 'd1u' ] ]
					for ii in range(len(namearrlist[1])) : 
						setattr( self, namearrlist[0]+namearrlist[1][ii], self.powerfitj[ii+ioffset] )
					ioffset = 7
					namearrlist = [ fhead+"fitconst{}".format(iimp)		, [ 'q3u', 'q1u', 'd1u' ] ]
					for ii in range(len(namearrlist[1])) : 
						setattr( self, namearrlist[0]+namearrlist[1][ii], self.powerfitj[ii+ioffset] )
				except :
					self.powerfitj = '.'
			try :	
				namearrlist = [ "powerfitt{}".format(iimp)		, [ 'xyu', 'yzu', 'zxu' ] ]
				fillpowerfractot = 0.
				for item in namearrlist[1] :
					dumpower = getattr( self, "powerfitt"+item )
					dumfill  = getattr( self, "fillimp"  +item )
					setattr( self, "fillpower"+item , dumpower*dumfill )
					setattr( self, "fillpowerfrac"+item , dumpower*dumfill/self.fillimptot*2 )
					fillpowerfractot += dumpower*dumfill/self.fillimptot*2 
				setattr( self, "fillpowerfractot" , fillpowerfractot )
	
				namearrlist = [ "powerfitj{}".format(iimp)		, [ 'q3u', 'q1u', 'd1u' ] ]
				fillpowerfractot = 0.
				for item in namearrlist[1] :
					dumpower = getattr( self, "powerfitj"+item )
					dumfill  = getattr( self, "fillimp"  +item )
					setattr( self, "fillpower"+item , dumpower*dumfill )
					setattr( self, "fillpowerfrac"+item , dumpower*dumfill/self.fillimptot*2 )
					fillpowerfractot += dumpower*dumfill/self.fillimptot*2 
				setattr( self, "fillpowerfractot" , fillpowerfractot )
			except :
				pass


		try :	
			fdat = self.tdir + "/result/u{:.3f}/linear_wn_max.dat".format( self.UF )
			linwnmax = np.genfromtxt( fdat, dtype=float ) 
			ioffset = 1
			namearrlist = [ "linwnmax"		, [ 'xyu', 'yzu', 'zxu' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], linwnmax[ii+ioffset] )
				setattr( self, 'Tlinwnmax'+namearrlist[1][ii], linwnmax[ii+ioffset]*11605./np.pi )
			#for ii in range(len(namearrlist[1])) : 
			#	setattr( self, namearrlist[0]+namearrlist[1][ii]+"Ratq1", self.powerfitj[ii+ioffset] / self.powerfitj[1+ioffset] )
		except :
			pass

		try :	
			self.NmaxOpt	= 512
			self.NplotOpt	= 128
			self.NintOpt	= 64
			self.EpsilonOpt	= 0.06
			self.dwOpt	= 0.03
			fdat = self.tdir + "/opt/opt_M{}_Np{}_N{}_dp{:.2f}_dw{:.2f}.dat".format( self.NmaxOpt, self.NplotOpt, self.NintOpt, self.EpsilonOpt, self.dwOpt )
			#print "FDAT OPT : ", fdat
			self.optdat = np.genfromtxt( fdat, dtype=float ) 
			if len( np.shape(self.optdat) )>1 : self.optdat = self.optdat[-1]
			ioffset = 1
			namearrlist = [ "opt"		, [ '' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.optdat[ii+ioffset] )
				setattr( self, "inverse"+namearrlist[0]+namearrlist[1][ii], 1./self.optdat[ii+ioffset] )
		except :
			self.opt = '.'
			self.inverseopt = '.'

		try :	
			self.NmaxOpt	= 512
			self.NplotOpt	= 128
			self.NintOpt	= 64
			self.EpsilonOpt	= 0.06
			self.dwOpt	= 0.03
			fdat = self.tdir + "/opt/x_u{:.2f}_M{}_Np{}_N{}_ep{:.2f}_dw{:.2f}.dat".format( self.UF, self.NmaxOpt, self.NplotOpt, self.NintOpt, self.EpsilonOpt, self.dwOpt )
			self.optdat = np.genfromtxt( fdat, dtype=float )[0]
			ioffset = 1
			namearrlist = [ "optw0"		, [ '' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.optdat[ii+ioffset] )
				setattr( self, "inverse"+namearrlist[0]+namearrlist[1][ii], 1./self.optdat[ii+ioffset] )
		except :
			self.optw0 = '.'
			self.inverseoptw0 = '.'

		try :	
			fdat = self.tdir + "/result/u{:.3f}/vHSpeak_fit_t2g.dat".format( self.UF )
			self.vHSpeakfitt = np.genfromtxt( fdat, dtype=float ) 
			ioffset = 1
			namearrlist = [ "vHSpeakfitt"		, [ 'xyu' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.vHSpeakfitt[ii+ioffset] )
				setattr( self, 'Zfitimt'+namearrlist[1][ii], 1./self.vHSpeakfitt[ii+ioffset] )
		except :
			self.vHSpeakfitt = '.'

		try :	
			fdat = self.tdir + "/result/u{:.3f}/vHSpeak_fit.dat".format( self.UF )
			self.vHSpeakfitj = np.genfromtxt( fdat, dtype=float ) 
			ioffset = 1
			namearrlist = [ "vHSpeakfitj"		, [ 'q1u' ] ]
			for ii in range(len(namearrlist[1])) : 
				setattr( self, namearrlist[0]+namearrlist[1][ii], self.vHSpeakfitj[ii+ioffset] )
				setattr( self, 'Zfitimj'+namearrlist[1][ii], 1./self.vHSpeakfitj[ii+ioffset] )
		except :
			self.vHSpeakfitj = '.'

		try :
			fdat = self.tdir + "/ongoing/energyLatt_u{:.2f}_{}th.txt".format( self.UF, self.count )
			self.energyLatt = np.genfromtxt( fdat, dtype=float ) 
		except :
			self.energyLatt = []
		try  :
			if len(self.energyLatt)>0 : 
				namearr = [ "Elatt", "KElatt", "IElatt", "ElattmuN", "KElattmuN" ]
				for a in range(len(self.energyLatt)) :
					setattr( self, namearr[a], self.energyLatt[a] )
		except :pass

		try :
			fdat = self.tdir + "/ongoing/occupationLatt_u{:.2f}_{}th.txt".format( self.UF, self.count )
			self.occupationLatt = np.genfromtxt( fdat, dtype=float ) 
			namearr = [ "nOcculatt" ]
			setattr( self, namearr[0], self.occupationLatt )
		except :
			self.occupationLatt = None

		oparr = [ "szsz", "lzlz", "jzjz", "mzmz", "szszStatic", "lzlzStatic", "jzjzStatic", "mzmzStatic" ] 
		try :		
			self.upper  = int( np.genfromtxt("inputs/CHIUPPER") )
		except :	self.upper=600
		try :		self.epsilon>0 
		except :	self.epsilon=5e-3
		for oo in range(len(oparr)) :
			try : 
				op = oparr[oo]
				name = "chitot"+op
				fdat = self.tdir + "/sus/chitot{}_u{:.2f}_{}th_n0_upper{}_ep{}.txt".format( op, self.UF, self.count, self.upper, self.epsilon )
				dat = np.genfromtxt( fdat, dtype=float ) 
				dat = dat[0] if len(np.shape(dat))>1 else dat
				setattr( self, name, dat ) 
				ioffset = [1,17]
				ifactor = [.5 , 1. ]
				namearrlist = [ name+"wn0"		, [ 'xy', 'zz' ] ]
				dumsum = 0
				for ii in range(len(namearrlist[1])) : 
					valii = getattr(self,name)[ioffset[ii]] 
					setattr( self, namearrlist[0]+namearrlist[1][ii], valii * ifactor[ii] )
					dumsum += valii
				setattr( self, name+"wn0", dumsum ) 
			except : 
				setattr( self, name, None ) 

		oparr = [ "szsz", "lzlz", "jzjz", "mzmz" ]
		for oo in range(len(oparr)) :
			try : 
				op = oparr[oo]
				opUpper = op[:1].upper()
				name = "chitot"+op

				chiStaticOffdiag	= getattr( self, name+"Staticwn0" )
				chiStaticDiagZ		= getattr( self, opUpper+"zdegMatLtcorrTDiagBeta" )
				chiStaticDiagX		= getattr( self, opUpper+"xdegMatLtcorrTDiagBeta" )
				chiStaticTotalxy	= getattr( self, name+"Staticwn0xy" ) + chiStaticDiagX
				chiStaticTotalzz	= getattr( self, name+"Staticwn0zz" ) + chiStaticDiagZ
				setattr( self, "chi"+opUpper+"StaticTotal", chiStaticOffdiag+chiStaticDiagZ+chiStaticDiagX*2 )
				setattr( self, "chi"+opUpper+"StaticTotalxy", chiStaticTotalxy )
				setattr( self, "chi"+opUpper+"StaticTotalzz", chiStaticTotalzz )
				setattr( self, "chi"+opUpper+"StaticTotalxyOverzz", chiStaticTotalxy/chiStaticTotalzz )

				chiLtBetaxy		= getattr( self, opUpper+"xdegMatLtcorrT" ) * self.IT
				chiLtBetazz		= getattr( self, opUpper+"zdegMatLtcorrT" ) * self.IT
				chiStaticDeltaxy	= chiStaticTotalxy - chiLtBetaxy
				chiStaticDeltazz	= chiStaticTotalzz - chiLtBetazz
				setattr( self, "chi"+opUpper+"StaticDeltaxy", chiStaticDeltaxy )
				setattr( self, "chi"+opUpper+"StaticDeltazz", chiStaticDeltazz )
			except : 
				pass
			op = oparr[oo]

		#try :
		#	fdat  = open( self.tdir + "/result/u{:.3f}/J.dat".format( self.UF )  )
		#	fdatlines = fdat.readlines()
		#	self.opJ  = fdatlines[-1].split("\n")[0].split()
		#	self.opJ  = np.array( self.opJ , dtype=float ) 
		#	for i in range( self.ni ) :
		#		ii = 3*i+1
		#		setattr( self, "opJ"+"%d"%i, np.array(self.opJ)[ [ii, ii+1, ii+2] ] )
		#	fdat.close()
		#except :
		#	self.opJ  = None
		for op in [ "J", "S", "L", "M" ] :
			try :
				count = -1
				fdat  = open( self.tdir + "/result/u{:.3f}/{}.dat".format( self.UF, op )  )
				fdatlines = fdat.readlines()
				if count == -1 : 
					ic = count
				else : 
					ic = count - self.count -1
				if len(fdatlines) > 1 : 
					setattr( self, "op"+op  , np.array( fdatlines[ic].split("\n")[0].split() , dtype=float ) )
				else : 
					setattr( self, "op"+op  , np.array( fdatlines[ic].split("\n")[0].split() , dtype=float ) )
				#print "reading count=", getattr( self, "op"+op )[0]
				for i in range( self.ni ) :
					ii = 3*i+1
					opdat = getattr( self, "op"+op  )
					setattr( self, "op"+op+"%d"%i, opdat[ [ii, ii+1, ii+2] ] )
				fdat.close()
			except :
				setattr( self, "op"+op  , None )

		
		try : 
			iimp = self.ni-1
			#if args.indimp : 
			#	iimp = int(args.indimp)
			datfile = "{}/ongoing/self{}_u{:.2f}_{}th.txt".format( self.tdir, iimp, self.UF, self.count )
			dat = open(datfile,'r')
			j = 0
			for line in dat : 
				if j==1 :
					datline = line.split("\n")[0]
					break 
				j += 1
			#print "size  : ", np.shape( datline.split()[1:] ) 
			selfmat		= np.array( datline.split()[1:] , dtype=float ).reshape( self.tNC, self.tNC*2 ) 
			selfmatcomp	= selfmat[:,::2] + 1j* selfmat[:,1::2]
			#if self.Treal == 0. :
			#	selfmatcomp = selfmatcomp * 0.
			if self.tdir.find("iDir")>-1 : 
				self.jind = ["q3u","q1u","q1d","q3d","d1u","d1d","tot"]
				for i  in  range(self.tNC) : 
					#print "imsigmawn0"+self.jind[i] , selfmatcomp[i][i].imag
					setattr( self, "imsigmawn0"+self.jind[i] , -selfmatcomp[i][i].imag )
				self.tind    = ["xyu","xyd","yzu","yzd","zxu","zxd"]
				selfmatcompt2g = transfMjefftot2g( selfmatcomp, Taraph )
				for i  in  range(self.tNC) : 
					#print "imsigmawn0"+self.tind[i] , -selfmatcompt2g[i][i].imag
					setattr( self, "imsigmawn0"+self.tind[i] , -selfmatcompt2g[i][i].imag )
			elif self.tdir.find("nDir")>-1 : 
				self.tind    = ["xyu","xyd","yzu","yzd","zxu","zxd"]
				for i  in  range(self.tNC) : 
					#print "imsigmawn0"+self.tind[i] , selfmatcomp[i][i].imag
					setattr( self, "imsigmawn0"+self.tind[i] , -selfmatcomp[i][i].imag )
				self.jind = ["q3u","q1u","q1d","q3d","d1u","d1d","tot"]
				selfmatcompjeff = transfMjefftot2g( selfmatcomp, np.conjugate(Taraph).transpose() )
				for i  in  range(self.tNC) : 
					#print "imsigmawn0"+self.jind[i] , -selfmatcompjeff[i][i].imag
					setattr( self, "imsigmawn0"+self.jind[i] , -selfmatcompjeff[i][i].imag )
			#print "SELFMATCOMP ::"
			#for i in selfmatcomp : 
			#	for j in i : 
			#		print "{:9.6f}  ".format( j ),
			#	print ""
			#sys.exit(1)
		except : pass

	def readMoment(self , op="J", count=int(-1) , args=None ) :
		try :
			fname = self.tdir + "/result/u{:.3f}/{}vec.dat".format( self.UF, op )
			print "Reading : ", fname
			fdat  = open( fname  )
			fdatlines = fdat.readlines()
			if count == -1 : 
				ic = count
			else : 
				ic = count - self.count -1
			setattr( self, "op"+op  , np.array( fdatlines[ic].split("\n")[0].split() , dtype=float ) )
			print "reading count=", getattr( self, "op"+op )[0]
			for i in range( self.ni ) :
				ii = 3*i+1
				opdat = getattr( self, "op"+op  )
				setattr( self, "op"+op+"%d"%i, opdat[ [ii, ii+1, ii+2] ] )
			fdat.close()
		except :
			setattr( self, "op"+op  , None )

	def readMomentLocal(self , op="J", count=int(-1) ) :
		self.readMomentLocal(op=op, count=count )
		try :
			fname = 'inputs/TRANSFLOCALcart'
			print "Reading : ", fname
			fdat  = np.genfromtxt( fname )

			for i in range( self.ni ) :
				ii = 3*i+1
				opdat = getattr( self, "op"+op  )
				setattr( self, "op"+op+"%d"%i, opdat[ [ii, ii+1, ii+2] ] )
			fdat.close()
		except :
			setattr( self, "op"+op  , None )
	def readTldos(self ) :
		self.tldosind    = ["xyu","xyd","yzu","yzd","zxu","zxd"]
		self.tldosmagind = [ 1.  , -1. , 1.  , -1. , 1.  , -1. ]
		self.maglatt     = 0
		try : 
			ftldos = open( self.tdir +"/lattice/tldos_u{:.2f}_Nplot{}.dat".format( self.UF, self.Nplot ) ) 
			ftldoslines = ftldos.readlines()
			efermi =  float( ftldoslines[self.Nplot/2].split()[0] ) 
			self.efermi	= efermi
			if abs(efermi)>1e-5 : 
				self.fillinglatt = "invalid:fermilevel:{}".format(efermi)
			else : 
				self.fillinglatt =  float( ftldoslines[self.Nplot/2].split()[-1] ) 
			datarr = ftldoslines[self.Nplot/2].split()[1:]
			fermidostot = 0.
			for ii in range(len(self.tldosind)) : 
				setattr( self, "fermidos"+self.tldosind[ii], float(datarr[ii]) ) 
				fermidostot += float(datarr[ii]) 
			for ii in range(len(self.tldosind)) : 
				setattr( self, "fermidosnorm"+self.tldosind[ii], float(datarr[ii])/fermidostot ) 
			self.fermidostot =  fermidostot
		except : 
			self.fillinglatt = "."
		try : 
			tldosupind    = ["xyu","yzu","zxu"]
			dospowertot = 0.
			for ii in range(len(tldosupind)) : 
				dospowertot += getattr( self, "fermidos"+tldosupind[ii] ) * getattr( self, "powerfitt"+tldosupind[ii] )  * 2
			self.dospowertot = dospowertot / self.fermidostot
		except : 
			pass
		try : 
			tldosdat = np.genfromtxt( self.tdir +"/lattice/tldos_u{:.2f}_Nplot{}.dat".format( self.UF, self.Nplot ) , dtype=float ) 
			tldos = np.zeros( self.tNC+1 ) 
			Etldos = 0.
			nldosdat = (len(tldosdat)+1)/2
			for j in range(nldosdat) : 
				for mu in range(self.tNC) : 
					tldos[mu] = tldos[mu] + tldosdat[j][mu+1]
					Etldos   += tldosdat[j][mu+1]*tldosdat[j][0]
			dw = tldosdat[1][0] - tldosdat[0][0] 
			tldos  = tldos * dw 
			Etldos*= dw
			#print "tldossum :", tldos[-1] , self.fillinglatt, dw
			for mu in range(self.tNC) :
				tldos[-1]    += tldos[mu]
				self.maglatt += tldos[mu] * self.tldosmagind[mu]
				setattr( self, self.tldosind[mu], tldos[mu] ) 
			setattr( self, "Etldos", Etldos ) 
			setattr( self, "filltldos", tldos[-1] ) 
			setattr( self, "tldossum", tldos[-1] ) 
			self.maglatt = "%9.6f"%self.maglatt
			#print "tldossum :", tldos[-1] , self.fillinglatt
		except :
			for mu in range(self.tNC) :
				setattr( self, self.tldosind[mu], ". " ) 
			self.maglatt   = ". "

		oparr = [ "S", "J", "M", "L" ]
		for op in oparr :
			try : 
				valsq = float( getattr( self, op+op+"degMateval" ) )
				valz  = float( getattr( self,   op+"zdegMateval" ) ) 
				valx  = float( getattr( self,   op+"xdegMateval" ) ) 
				valy  = float( getattr( self,   op+"ydegMateval" ) ) 
				dval = valsq - valz*valz - valx*valx - valy*valy 
				setattr( self, "Delta"+op , dval ) 
			except:
				setattr( self, "Delta"+op , '.' ) 
	def readLdos(self ) :
		self.ldosind    = ["q3u","q1u","q1d","q3d","d1u","d1d"]
		self.ldosmagind = [3./2 , 1./2,-1./2,-3./2, 1./2,-1./2]
		self.maglattj   = 0
		try : 
			#fldos = open( self.tdir +"/lattice/ldos_u{:.2f}_Nplot{}.dat".format( self.UF, self.Nplot ) ) 
			#fldoslines = fldos.readlines()
			ldosdat = np.genfromtxt( self.tdir +"/lattice/ldos_u{:.2f}_Nplot{}.dat".format( self.UF, self.Nplot ) ) 
			ndatldos = len(ldosdat) 
			efermi =  float( ldosdat[(ndatldos)/2][0] )
			self.efermi	= efermi 
			if abs(efermi)>1e-5 : 
				self.fillinglattj = "invalid:fermilevel"
			else : 
				self.fillinglattj =  float( ldosdat[(ndatldos)/2][-1] )
			datarr = ldosdat[ndatldos/2][1:]
			fermidostot = 0.
			for ii in range(len(self.ldosind)) : 
				setattr( self, "fermidos{}".format(ii), float(datarr[ii]) ) 
				fermidostot += float(datarr[ii]) 
			for ii in range(len(self.ldosind)) : 
				setattr( self, "fermidosnorm{}".format(ii), float(datarr[ii])/fermidostot ) 
			self.fermidostot =  fermidostot
		except : 
			self.fillinglattj = "."
		
		try :
			ldosdat = np.genfromtxt( self.tdir +"/lattice/ldos_u{:.2f}_Nplot{}.dat".format( self.UF, self.Nplot ) , dtype=float ) 
			fillldos = np.zeros( self.tNC+1 ) 
			nldosdat = (len(ldosdat)+1)/2
			for j in range(nldosdat) : 
				for mu in range(self.tNC) : 
					fillldos[mu] = fillldos[mu] + ldosdat[j][mu+1]
			dw = ldosdat[1][0] - ldosdat[0][0] 
			fillldos = fillldos * dw 
			#print "ldossum :", ldos[-1] , self.fillinglattj, dw
			for mu in range(self.tNC) :
				fillldos[-1]  += fillldos[mu]
				self.maglattj += fillldos[mu] * self.ldosmagind[mu]
				setattr( self, self.ldosind[mu], fillldos[mu] ) 
			self.maglattj = "%9.6f"%self.maglattj
			self.fillldos = fillldos[-1]
			setattr( self, "ldossum", fillldos[-1] ) 
			#print "ldossum :", ldos[-1] , self.fillinglatt
		except : 
			for mu in range(self.tNC) :
				setattr( self, self.ldosind[mu], ". " )
			self.maglattj = ". "
			self.fillldos = ". "
	def readPfile(self ) :
		try : 
			if self.Hz is None : 
				pfilelist = "JTD{}_S{}_J{}_U{}_D{}_{}.p.o* ".format( self.JT, self.S, self.J, self.U, self.D ,self.initmak ) 
				pfilelist +="JTD{}_S{}_J{}_U{}_D{}_beta{}_{}.p.o* ".format( self.JT, self.S, self.J, self.U, self.D, self.beta, self.initmak ) 
				if abs(self.JT)<1e-5 : 
					pfilelist +="JTD0_S{}_J{}_U{}_D{}_{}.p.o* ".format( self.S, self.J, self.U, self.D ,self.initmak ) 
					pfilelist +="JTD0_S{}_J{}_U{}_D{}_beta{}_{}.p.o* ".format( self.S, self.J, self.U, self.D, self.beta, self.initmak ) 
			else : 
				pfilelist = "JTD{}_S{}_J{}_U{}_D{}_Hz{}_{}.p.o* ".format( self.JT, self.S, self.J, self.U, self.D, self.Hz ,self.initmak ) 
				pfilelist +="JTD{}_S{}_J{}_U{}_D{}_Hz{}_IT{}_{}.p.o* ".format( self.JT, self.S, self.J, self.U, self.D, self.Hz, self.IT ,self.initmak ) 
				pfilelist +="JTD{}_S{}_J{}_U{}_D{}_Hz{}_beta{}_{}.p.o* ".format( self.JT, self.S, self.J, self.U, self.D, self.Hz, self.beta, self.initmak ) 
				pfilelist +="JTD{}_S{}_J{}_U{}_D{}_Hz{}_beta{}_IT{}_{}.p.o* ".format( self.JT, self.S, self.J, self.U, self.D, self.Hz, self.beta, self.IT, self.initmak ) 
				if abs(self.JT)<1e-5 : 
					pfilelist +="JTD0_S{}_J{}_U{}_D{}_Hz{}_{}.p.o* ".format( self.S, self.J, self.U, self.D, self.Hz ,self.initmak ) 
					pfilelist +="JTD0_S{}_J{}_U{}_D{}_Hz{}_IT{}_{}.p.o* ".format( self.S, self.J, self.U, self.D, self.Hz , self.IT ,self.initmak ) 
					pfilelist +="JTD0_S{}_J{}_U{}_D{}_Hz{}_beta{}_{}.p.o* ".format( self.S, self.J, self.U, self.D, self.Hz, self.beta, self.initmak ) 
					pfilelist +="JTD0_S{}_J{}_U{}_D{}_Hz{}_beta{}_IT{}_{}.p.o* ".format( self.S, self.J, self.U, self.D, self.Hz, self.beta, self.IT , self.initmak ) 
			dumlist = ""
			for a in pfilelist.split() : dumlist = dumlist + "{}/{} ".format( self.tdir,a )
			for a in pfilelist.split() : dumlist = "{} ".format( dumlist ) + a 
			pfilelist = dumlist
			proc = subprocess.Popen(['ls -r '+pfilelist], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
			(out, err) = proc.communicate()
			searchpfilelist = out.split("\n")
			self.pfileoutyn = "x"
			for pfile in searchpfilelist : 
				if len(pfile)>0 : 
					proceg = subprocess.Popen(['egrep '+self.Dirname+"/"+" "+pfile], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
					(outeg, erreg) = proceg.communicate()
					for outegline in outeg.split("\n") :
						if outegline.find("A new directory")>-1 : 
							indstr = outegline.find("A new directory")
							self.pfileout = pfile
							self.pfileoutyn = "y"
							self.pfileoutid = self.pfileout.split("o")[-1]
							self.pfile	= self.pfileout.split(".p.o")[0] + ".p"
				if self.pfileoutyn=="y" : break
		except :
			self.pfile = "None"
		try : 
			self.pfileoutyn.find("y")
		except : 
			self.pfileout   = "None"
			self.pfileoutyn = "x"
			self.pfileoutid = "None"
		if self.pfileoutyn=="y" :
			fpfile = open( self.pfileout ) 
			fpfilelines = fpfile.readlines()
			if fpfilelines[-1].find("time-end")>-1 :
				self.pfileend = "e"
			else : 
				self.pfileend = "n"
		try : 
			if self.pfileout != "None" :
				if self.pfileout.find(self.tdir)>-1 : 
					rawpfileout = self.pfileout[ self.pfileout.find(self.tdir)+len(self.tdir) :]
					dirpfileout = self.pfileout
					if rawpfileout[0].find("/")>-1 : rawpfileout = rawpfileout[1:]
				else :
					rawpfileout = self.pfileout
					dirpfileout = self.tdir + "/" + self.pfileout
				cmdcheck =  "cmp -s "+rawpfileout+" "+dirpfileout
				if os.system(cmdcheck) and os.system("test -e "+rawpfileout)<1 :
					os.system( "cp -u "+rawpfileout+" "+self.tdir+"/" ) 
					os.system( "cp -n "+self.pfile   +" "+self.tdir+"/" ) 
			elif self.pfileout == "None" :
				self.pfileoutyn = "n"
				self.pfileend   = "n"
		except : 
			self.pfileoutyn = "x"
			self.pfileend   = "x"
		try : 
			if self.pfileoutid.find(".")>-1 : 
				idgrep = self.pfileoutid[:self.pfileoutid.find(".")]
			else :
				idgrep = self.pfileoutid
			proc = subprocess.Popen(["qstat -j "+idgrep+" | grep "+idgrep], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
			(out, err) = proc.communicate()
			if len(out) > 0 :
			  	#print "IDG : ", idgrep
				proc = subprocess.Popen(["qstat | grep "+idgrep], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
				(out, err) = proc.communicate()
				outlist = out.split()
				self.qstat = outlist[4]
			else :
				self.qstat = "n"
		except : 
			self.qstat = "x"
		self.status = self.pfileoutyn + self.pfileend + self.qstat
		#read the time-elapsed
		try :
			fpout = open( self.pfileout ) 
			fpoutlines = fpout.readlines()
			dtstart = dateclass()
			for line in fpoutlines[:20] : 
				if line.find("time")>-1 and line.find("start")>-1  :
					dtstart.convnormal( line[line.find("start")+6:] )
			dtend   = dateclass()
			for line in fpoutlines[-2:] : 
				if line.find("time")>-1 and line.find("end")>-1  :
					dtend.convnormal( line[line.find("end")+5:] )
			self.elapsedt = dtend.diffdate( dtstart ) 
		except :
			self.elapsedt = "."
			pass
	def readEcluster(self ) :
		eclusterdat = np.genfromtxt(self.tdir+"/Ecluster.dat")
		dimdim = len(eclusterdat) 
		self.nbasis = int(np.sqrt(dimdim))
		nbasis = self.nbasis
		#if np.sqrt(dimdim) != 6	: print "ERROR::reading Ecluster.dat, size of matrix is incompatible." ; sys.exit(1)
		#else			: nbasis = 6 
		self.ecluster = np.array( [np.zeros(nbasis)]*nbasis , dtype=complex ) 
		for edat in eclusterdat : 
			i = int(edat[0])
			j = int(edat[1])
			d = edat[2] + 1j*edat[3]
			self.ecluster[i][j] = d 
		return self.ecluster 
	def returnNatTransf( self ) :
		try : 
			a = self.eclustersoc[1][1] 
			b = self.eclustersoc[1][4] 
			c = self.eclustersoc[4][4] 
			d = self.eclustersoc[5][5] 
		except : 
			a = self.ecluster[1][1] 
			b = self.ecluster[1][4] 
			c = self.ecluster[4][4]  + self.S*1.5
			d = self.ecluster[5][5]  + self.S*1.5
		if (abs( a - self.ecluster[2][2] ) < 1e-5 ) and (abs( b + self.ecluster[2][5] ) < 1e-5 ) and (abs( b - self.ecluster[4][1].conjugate() ) < 1e-5 ) and (abs( c - d ) < 1e-5 ) :  pass
		else :
		      print (abs( a - self.ecluster[2][2] ) < 1e-5 ) ,(abs( b + self.ecluster[2][5] ) < 1e-5 ) ,(abs( b - self.ecluster[4][1].conjugate() ) < 1e-5 ) ,(abs( c - d ) < 1e-5 ) 
		      print "ERROR :: Ecluster's symmetry is not the same as one expected." 
		try : dum = a*a + 4*b*b - 2*a*c + c*c ; dum = np.sqrt(dum)
		except : print "ERROR :: sqrt({})".format(dum) 
		v1 =  np.array( [ (a-c-dum)/2./b , 1. ] ) ; v1 = v1/np.sqrt(np.abs(( v1[0]*v1[0]+v1[1]*v1[1] )))
		v2 =  np.array( [ (a-c+dum)/2./b , 1. ] ) ; v2 = v2/np.sqrt(np.abs(( v2[0]*v2[0]+v2[1]*v2[1] )))
		#print "v1 : ", v1 , v1[0]*v1[0] , v1[1]*v1[1] 
		#print "v2 : ", v2,  v2[0]*v2[0] , v2[1]*v2[1] 
		self.matTjnat = np.matrix( [
				[ 1, 0, 0, 0, 0, 0 ],
				[ 0, 0, 0, 1, 0, 0 ],
				[ 0,v1[0], 0, 0,v1[1], 0 ],
				[ 0, 0, -v1[0], 0, 0, v1[1] ],
				[ 0,v2[0], 0, 0, v2[1], 0 ],
				[ 0, 0, -v2[0], 0, 0, v2[1] ]
		] , dtype=complex)
		return self.matTjnat
	def readEclustersoc(self ) :
		self.eclustersoc = self.readEcluster()
		try : 
			self.eclustersoc[4][4] += 1.5 *self.S
			self.eclustersoc[5][5] += 1.5 *self.S
		except:
			pass 
		return self.eclustersoc
	def readEclustersocnat(self ) :
		self.readEclustersoc() 
		Tnat = self.returnNatTransf()
		self.eclustersocnat = transfMconjtr( self.eclustersoc , Tnat )
		return self.eclustersocnat
	def readEclustersoct2g(self ) :
		self.readEclustersoc() 
		self.eclustersoct2g = transfMtrconj( self.eclustersoc , Tara )
		return self.eclustersoct2g
	def returnLabel( self , var ) :
		if var == "beta" :
			varname = r"\beta" 
		elif var == "J" :
			varname = "J_H" 
		elif var == "S" :
			varname = r"\lambda_{SOC}" 
		elif var == "UF" :
			varname = "U"
		else :
			varname = var 
		dum = r"${}$={}".format( varname , getattr( self, var ) )
		return dum 
	def readDensity(self ) :
		fden = open( self.tdir + "/ongoing/flattj_u{:.2f}_{}th.txt".format(self.UF, self.count) ) 
		fdlines = fden.readlines()
		self.jind = ["q3u","q1u","q1d","q3d","d1u","d1d","tot"]
		for ll in range( len(fdlines) ) :
			fdlinelist = fdlines[ll].split()
			for kk in range( len(fdlinelist) ) :
				setattr( self, "flattj"+self.jind[kk]+"%d"%(ll), float(fdlinelist[kk])  ) 
				setattr( self, "flattj"+self.jind[kk],		 float(fdlinelist[kk])  ) 

		fden = open( self.tdir + "/ongoing/flattut2g_u{:.2f}_{}th.txt".format(self.UF, self.count) ) 
		fdlines = fden.readlines()
		self.tldosind    = ["xyu","xyd","yzu","yzd","zxu","zxd","tot"]
		for ll in range( len(fdlines) ) :
			fdlinelist = fdlines[ll].split()
			for kk in range( len(fdlinelist) ) :
				setattr( self, "flattut2g"+self.tldosind[kk]+"%d"%(ll), float(fdlinelist[kk])  ) 
				setattr( self, "flattut2g"+self.tldosind[kk],		float(fdlinelist[kk])  ) 
	def readGptlEnergy(self,args, deg=None) :
		try :
			pfile = self.pfileout
		except : 
			pfile = None
		if pfile : 
		        nsector = 25
		        if args.nsector : nsector = int(args.nsector)
			stdout, stderr = stdoutgrep( "egrep \"rank\" {} ".format( pfile )  )
			Egptlarr = [ [] for a in range(nsector) ]
			for a in  stdout.split("\n") :
				print a.split()
				if len(a)>0 :
					if a.split()[0].find("rank")>-1 :
						sector , Esector = a.split()[2] , a.split()[5]
						intsector = int(sector.split("th")[0])
						Egptlarr[ intsector ].append(  float(Esector) )
	        	print "shape(EGPTLARR)\t: ", np.shape(Egptlarr)
			E0arr, sectorE0arr = findgs( Egptlarr )
			print "e0arr : ", E0arr
		        print "e0secarr : ",  sectorE0arr
        	if deg :
        	        sectorDeg = int( hobj.gptl_0 )
        	        if args.degsector :
        	                sectorDeg = int( args.degsector )
        	        stdout, stderr = stdoutgrep( "egrep \"{}:\" {} ".format( sectorDeg, pfile )  )
        	        E0degarrAll = []
        	        nline = 0
        	        for a in  stdout.split("\n") :
        	                if len(a)>0 and a.find("ndeg")>-1 :
        	                        nline+=1
        	                        if nline%2 < 1 :
        	                                print "O : ", a.split()
        	                                e0dum = [ a.split()[1][:-1] ]
        	                                for item in a.split()[2:-6:4] :
        	                                        e0dum.append( item )
        	                                E0degarr = np.array( e0dum , dtype=float )
        	                                print "E : ", E0degarr
        	                                E0degarrAll.append( E0degarr )
        	        E0degarrAll = np.array( E0degarrAll )
        	        print "sectorDeg : ", sectorDeg
        	        print "shape(E0degarrAll) : ", np.shape(E0degarrAll)
        	        ndat, nvec = np.shape(E0degarrAll)
        	        if args.degnvector :
        	                nvec = int( args.degnvector )
        	        xdat = np.arange( ndat )
        	        print "E0deg[vec0] : ", E0degarrAll.transpose()[0]
        	        for ivec in range(nvec) :
        	                ydat = E0degarrAll.transpose()[ivec]
        	                print "YDAT [vec{}]: ".format(ivec), ydat
        	                axnow.plot( xdat , ydat , "o-" , label="{}th-vec".format(ivec) , dashes=[ivec+1,ivec+1], mec='gray' )
        	        ylab = "Energy of $n$th lowest state in gptl={}".format( sectorDeg )
        	else :
        	        for isec in range(nsector) :
        	                if args.diff :
        	                        dat = np.array(Egptlarr[isec],dtype=float) - E0arr
        	                else :
        	                        dat = Egptlarr[isec]
        	                xdat= range( len(dat) )
        	                print "xDAT,yDAT %4d: "%a, xdat, " ", dat
        	                if args.dsector :
        	                        sec0 = int( args.dsector.split("_")[0] )
        	                        sec1 = int( args.dsector.split("_")[1] )
        	                else :
        	                        sec0 = 12
        	                        sec1 = 17
        	                if isec > sec0 and isec<sec1  :
        	                        leglabel = "{}".format(a)
        	                        axnow.plot( xdat , dat , colors(idir)+markersdiff[isec%10]+"-" , label=leglabel, mec='gray' )
        	                        axnow.text( xdat[0] , dat[0] , "{}".format( isec ) )
        	                else                    :
        	                        leglabel = None
        	        ylab = "Energy difference in gptl={}".format( sectorDeg )

	def printParametersAll(self ) :
		print self.tdir +'/parameters_u{:.2f}.dat'.format(self.UF), "\t", self.iarg
		os.system( "cat " + self.tdir +'/parameters_u{:.2f}.dat'.format(self.UF) )
	def printMagAll(self ) :
		fout = self.tdir +"/result/u{:.3f}/mag.dat".format(self.UF)
		print fout, "\t", self.iarg
		os.system( "cat " + fout ) 
	def printFillAll(self ) :
		fout = self.tdir +"/result/u{:.3f}/filling.dat".format(self.UF)
		print fout, "\t", self.iarg
		os.system( "cat " + fout ) 
	def printFillt2gAll(self ) :
		fout = self.tdir +"/result/u{:.3f}/fillingt2g.dat".format(self.UF)
		print fout, "\t", self.iarg
		os.system( "cat " + fout ) 
	def printMassAll(self , datname=None ) :
		if datname is None : 
			fout = self.tdir +"/renorm_*" 
			print fout, "\t", self.iarg
			os.system( "cat " + fout ) 
		elif datname.find("renorm")>-1 :
			fout = self.tdir +"/"+datname+".dat"
			proc = subprocess.Popen(["cat " + fout], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
			(out, err) = proc.communicate()
			if out=='' : 
				print '.'
			else : 
				print out,
	def printDensityAll(self , datname=None ) :
		if datname is None : 
			fout = self.tdir + "/ongoing/flattj_u{:.2f}_{}th.txt".format(self.UF, self.count)
			print fout, "\t", self.iarg
			os.system( "cat " + fout ) 
			fout = self.tdir + "/ongoing/flattut2g_u{:.2f}_{}th.txt".format(self.UF, self.count)
			print fout, "\t", self.iarg
			os.system( "cat " + fout ) 
		elif datname.find("j_")>-1 : 
			fout = self.tdir + "/ongoing/flattj_u{:.2f}_{}th.txt".format(self.UF, self.count)
			proc = subprocess.Popen(["cat " + fout], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
			(out, err) = proc.communicate()
			print out,
		elif datname.find("t2g_")>-1 : 
			fout = self.tdir + "/ongoing/flattut2g_u{:.2f}_{}th.txt".format(self.UF, self.count)
			proc = subprocess.Popen(["cat " + fout], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
			(out, err) = proc.communicate()
			print out,
		#fout = open( self.tdir + "/ongoing/felatt_u{:.2f}_{}th.txt".format(self.UF, self.count) ) 
		#for it in fout.readlines()[0].split() : 
		#	print self.nc*float(it), "\t\t",
		#print ""
	def printPfilename(self ) :
		try : 
			print self.pfile, "\t", self.pfileoutid, "\t", self.iarg
		except : 
			print ".\t",  self.iarg
	def printEtldos(self ) :
		try : 
			print self.Etldos, "\t", self.iarg
		except : 
			print ".\t",  self.iarg
	def calcJtot(self ) :
		jjp1 = [ 15./4 , 15./4 , 15./4 , 15./4 , 3./4 , 3./4 ] 
		#self.jind = ["q3u","q1u","q1d","q3d","d1u","d1d","tot"]
		self.jjp1 = 0.
		for jj in range(self.tNC) :
			self.jjp1 += jjp1[jj] * getattr( self, "fillimp"+self.jind[jj] )
		self.jjp1unit = self.jjp1/float(self.filling)
		self.jjp1cell = self.jjp1/float(self.filling)*round(float(self.filling)) 
		self.jjp1hole = self.jjp1/float(self.filling)*( self.tNC - round(float(self.filling)) ) 

		self.jtot     = np.sqrt( self.jjp1 +0.25 ) - 0.5  
		self.jtotunit = np.sqrt( self.jjp1unit +0.25 ) - 0.5  
		self.jtotcell = np.sqrt( self.jjp1cell +0.25 ) - 0.5  
		self.jtothole = np.sqrt( self.jjp1hole +0.25 ) - 0.5  
		self.jtotprintind = [ "jjp1", "jjp1unit", "jtot", "jtotunit", "jtotcell" , "jtothole" ]
	def printJtot(self ) :
		try : 
		 	for item in self.jtotprintind : 
		 		print getattr( self, item ) , "\t",
			print self.iarg
		except : 
			print ".\t",  self.iarg

	def readMassLambda(self ) :
		try : 
			fmass = open( self.tdir +"/renorm_im.dat" )
			lines = fmass.readlines()
			mimt2g = lines[0].split()
			fmass.close()
			self.mimt2g = np.array( mimt2g , dtype=float ) 
		except : 
			self.mimt2g = ["."]*self.nb
		try : 
			fmass = open( self.tdir +"/renorm_im_j.dat" )
			lines = fmass.readlines()
			mimj = lines[0].split()
			self.mimj = np.array( mimj , dtype=float ) 
			fmass.close()
		except : 
			self.mimj = ["."]*self.nb
		try : 
			fmass = open( self.tdir +"/renorm_re.dat" )
			lines = fmass.readlines()
			lambdaxy = lines[0].split()[2:]
			lambdaz  = lines[1].split()[2:]
			if abs( float(lambdaxy[0]) - float(lambdaxy[1]) ) > 1e-5 : 
				print "Discrepancy in lambdaxy : ", lambdaxy[0] , lambdaxy[1] 
			if abs( float(lambdaz[0]) - float(lambdaz[1]) ) > 1e-5 : 
				print "Discrepancy in lambdaz : ", -float(lambdaz[0]) , -float(lambdaz[1]) 
			mret2g = lines[2].split()
			self.mret2g = np.array( mret2g , dtype=float ) 
			self.lambdaxy = float(lambdaxy[0])
			self.lambdaz  = -float(lambdaz[0])
			fmass.close()
		except : 
			self.mret2g = ["."]*self.nb
			self.lambdaxy = "."
			self.lambdaz  = "."

	def printMassLambda(self ) :
		print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
				self.UF, self.J, self.S, self.D, self.lambdaxy, self.lambdaz,
				self.mret2g[0], self.mret2g[2],
				self.mimt2g[0], self.mimt2g[2],
				self.mimj[0], self.mimj[1], self.mimj[4] ) 
	def pathStaticRdm(self ) :
		rdmeval   = self.tdir +"/result/u{:.3f}/rdm.dat".format(self.UF)
		ncrdmevec = self.tdir +"/result/u{:.3f}/ncrdmevec.dat".format(self.UF)
		return  rdmeval, ncrdmevec 
	def pathResult(self ) :
		return self.tdir +"/result/u{:.3f}".format(self.UF)
	def pathParameters(self ) :
		return self.tdir +'/parameters_u{:.2f}.dat'.format(self.UF)
	def returnattrarr(self , attrname , args=None ) :
		plist = self.parindlist
		print "PLIST : ", plist 
		j = 0 
		j0 = 0 
		for j in range(len(plist)) :
			if plist[j].find(attrname)>-1 : 
				j0 = j 
		if args and args.impindex :
			for j in range(len(plist)) :
				if plist[j].find(attrname+"_%s"%args.impindex)>-1 : 
					j0 = j 
		print "Returning :: ", plist[j0]
		pardatarr = np.genfromtxt( self.pathParameters() ) ;
		if attrname.find(".dat")>-1 : 
			print "Returning :: ", attrname
			return self.returndatfile(attrname)
		try :
			return pardatarr.transpose()[j0]
		except : 
			print  "ERROR in 'returnattrarr()'."
			return "nan"
	def returnattrall(self , attrname ) :
		try :
			return getattr( self, attrname ) 
		except :
			return float('nan')
	def returndatfile(self , attrname ) :
		datname = attrname.split(".dat")[0]
		datind  = int(attrname.split(".dat")[-1])
		dat = np.genfromtxt( self.pathResult()+"/"+datname+".dat" ) ;
		try :
			return dat.transpose()[datind]
		except :
			return float('nan')
	def readMakClass(self, nmak=None, imakstart=None, imakend=None ) :
		makclassArr	= []
		makcheck	= True
		if nmak is None		: nmak		= self.ncount 
		if imakstart is None	: imakstart	= self.countstart
		if imakend is None	: imakend	= self.count+1
		self.nmak	= nmak
		for imak in range(imakstart,imakend) :
			makname = self.tdir+"/u{:.2f}_{}th.mak".format(self.UF,imak)
			if( os.system("test -f {}".format(makname)) ) :
				print "WARNNING :: file missed. ", makname 
			else :
				print "Loading: {} ".format(makname)
				setattr( self, 'makclass_i{}'.format(imak), makclass(makname, Ni=self.Ni, tNC=self.tNC, tNb=self.tNb, basis=self.basis) )
			imak	+= 1

		makname = self.tdir+"/next.mak".format(self.UF,imak)
		if( os.system("test -f {}".format(makname)) ) :
			print "WARNNING :: file missed. ", makname 
		else :
			print "Loading: {} ".format(makname)
			setattr( self, 'makclass_i-1', makclass(makname, Ni=self.Ni, tNC=self.tNC, tNb=self.tNb, basis=self.basis) )
	def returnEbathIter(self, nmak=None, imakstart=None, imakend=None ) :
		returnEbathIterImp(nmak=nmak, imakstart=imakstart, imakend=imakend)
	def returnEbathIterImp(self, nmak=None, imakstart=None, imakend=None, iimp=0 ) :
		try : 
			self.nmak > 0 
		except : 
			self.readMakClass()
		makclassArr	= []
		makcheck	= True
		if nmak is None		: nmak		= self.ncount 
		if imakstart is None	: imakstart	= self.countstart
		if imakend is None	: imakend	= self.count+1
		self.EbathIter	= [ getattr( self, 'makclass_i{}'.format(imak) ).eildat[iimp] for imak in range(imakstart,imakend) ]
		self.EbathIter	= np.array( self.EbathIter )
		return self.EbathIter
	def setValue(self, item, valArr, iimp=0 ) :
		if (item == "powerfitt") or (item == "powerfitj") or (item == "powerplusconstfitt") or (item == "powerplusconstfitj") or (item == "poly4fitt") or (item == "poly4fitj") or (item == "linfitt") or (item == "linfitj") :
			if len(valArr) != 3 :
				print "Invalid number of values. exit."
			else : 
				try :	
					fbody	= 'fitSelfIm'
					if (item == "powerfitt") :
						fhead	= "power_"
						ftail	= "_t2g"
					elif (item == "powerfitj") :
						fhead	= "power_"
						ftail	= ""
					if (item == "powerplusconstfitt") :
						fhead	= "powerplusconst_"
						ftail	= "_t2g"
					elif (item == "powerplusconstfitj") :
						fhead	= "powerplusconst_"
						ftail	= ""
					elif (item == "poly4fitt") :
						fhead	= "poly4_"
						ftail	= "_t2g"
					elif (item == "poly4fitj") :
						fhead	= "poly4_"
						ftail	= ""
					elif (item == "linfitt") :
						fhead	= "mass_"
						ftail	= "_t2g"
						fbody	= 'fitIm'
					elif (item == "linfitj") :
						fhead	= "mass_"
						ftail	= ""
						fbody	= 'fitIm'
					fdat = self.tdir + "/result/u{:.3f}/{}{}{}.dat".format( self.UF , fhead, fbody, ftail )
					print "Modifying : ", fdat

					ftest = open( fdat, 'r' ) 
					ftest.readlines()
					ftest.close()

					dat = np.genfromtxt( fdat, dtype=float ) 
				except :
					dat = np.zeros(18)
					print "File is newly created : ", fdat

				with open( fdat, 'w' ) as fwfit :
					i=0;
					#print "DAT : ", dat
					fwfit.write( "{:3d} ".format( int(dat[i])) ) ; i+=1
					for j in range(3) : 
						fwfit.write( "{} ".format(dat[i]) ) ; i+=1
					for j in range(3) : 
						fwfit.write( "{:12f} ".format( valArr[j] ) ) ; i+=1
					for j in range(3) : 
						fwfit.write( "{:12f} ".format( dat[i] ) ) ; i+=1
					fwfit.write( "\t" ) 
					for it in dat[i:] : 
						fwfit.write( "{}\t".format(it) )
					fwfit.write( "\n" ) 
		elif (item == "powerfitconstt") or (item == "powerfitconstj") or (item == "powerplusconstfitconstt") or (item == "powerplusconstfitconstj") or (item == "poly4fitconstt") or (item == "poly4fitconstj") :
			if len(valArr) != 3 :
				print "Invalid number of values. exit."
			else : 
				try :	
					if (item == "powerfitconstt") :
						fhead	= "power_"
						ftail	= "_t2g"
					elif (item == "powerfitconstj") :
						fhead	= "power_"
						ftail	= ""
					elif (item == "powerplusconstfitconstt") :
						fhead	= "powerplusconst_"
						ftail	= "_t2g"
					elif (item == "powerplusconstfitconstj") :
						fhead	= "powerplusconst_"
						ftail	= ""
					elif (item == "poly4fitconstt") :
						fhead	= "poly4_"
						ftail	= "_t2g"
					elif (item == "poly4fitconstj") :
						fhead	= "poly4_"
						ftail	= ""
					fdat = self.tdir + "/result/u{:.3f}/{}fitSelfIm{}.dat".format( self.UF , fhead, ftail )

					ftest = open( fdat, 'r' ) 
					ftest.readlines()
					ftest.close()

					dat = np.genfromtxt( fdat, dtype=float ) 
				except :
					dat = np.zeros(18)
					print "File is newly created : ", fdat

				with open( fdat, 'w' ) as fwfit :
					i=0;
					#print "DAT : ", dat
					fwfit.write( "{:3d} ".format( int(dat[i])) ) ; i+=1
					for j in range(3) : 
						fwfit.write( "{} ".format(dat[i]) ) ; i+=1
					for j in range(3) : 
						fwfit.write( "{:12f} ".format( dat[i] ) ) ; i+=1
					for j in range(3) : 
						fwfit.write( "{:12f} ".format( valArr[j] ) ) ; i+=1
					fwfit.write( "\t" ) 
					for it in dat[i:] : 
						fwfit.write( "{}\t".format(it) )
					fwfit.write( "\n" ) 
		elif (item[:-1] == "chistaticVV") :
			if len(valArr) != 2 :
				print "Invalid number of values in '{}'. exit.".format(item)
			else : 
				op	= item[-1]
				oplower	= op.lower()
				opupper	= op.upper()
				fhead	= "chitot{}z{}zStatic".format( oplower, oplower ) 
				ftail	= ""

				fdat = self.tdir + "/sus/{}_u{:.2f}_{}th_n0_upper600_ep0.005.txt".format( fhead, self.UF, self.count ) 

				try :	
					ftest = open( fdat, 'r' ) 
					fdum = ftest.readlines()
					fdum[0]
					ftest.close()

					dat = np.genfromtxt( fdat, dtype=float )
				except :
					dat = np.zeros( (5,19) )
					os.system( "mkdir "+self.tdir+"/sus" )
					print "File is newly created : ", fdat

				with open( fdat, 'w' ) as fwfit :
					for irow in range( dat.shape[0] ) :
						datrow	 = dat[irow,:]
						i=0;
						#print "DAT : ", datrow
						fwfit.write( "{:.3f} \t".format( datrow[i] ) ) ; i+=1
						fwfit.write( "{:12f} {:12f}\t".format( 2.*valArr[0],0. ) ) ; i+=2
						for j in range(7) : 
							fwfit.write( "{:12f} {:12f}\t".format( 0.,0. ) ) ; i+=2
						fwfit.write( "{:12f} {:12f}\t".format( valArr[1],0. ) ) ; i+=2
						fwfit.write( "\n" ) 
		elif (item[:-1] == "chistaticCurie") or (item[:-1] == "chistaticLT") :
			if len(valArr) != 2 :
				print "Invalid number of values in '{}'. exit.".format(item)
			else : 
				directArr = [ "x", "z" ]
				for idirect in range( len(directArr) ) :
					direct	= directArr[idirect]
					op	= item[-1]
					oplower	= op.lower()
					opupper	= op.upper()
					fhead	= opupper + direct
					ftail	= ""

					fdat = self.tdir + "/result/u{:.3f}/{}degMatLtcorrT.dat".format( self.UF, fhead )

					
					try : 
						ftest = open( fdat, 'r' ) 
						fdum = ftest.readlines()
						ftest.close()
	
						dat = np.genfromtxt( fdat, dtype=float )
					except : 
						dat = np.zeros( 8 )
						dat[0] = self.count
						print "File is newly created : ", fdat
	
					with open( fdat, 'w' ) as fwfit :
						i=0;
						#print "DAT : ", dat
						fwfit.write( "{:3d} ".format( int(dat[i])) )	; i+=1
						if (item[:-1] == "chistaticCurie") :
							fwfit.write( "{:12f} ".format( dat[i] ) )		; i+=1
							fwfit.write( "{:12f} ".format( valArr[idirect] ) )	; i+=1
						elif (item[:-1] == "chistaticLT") :
							fwfit.write( "{:12f} ".format( valArr[idirect] ) )	; i+=1
							fwfit.write( "{:12f} ".format( dat[i] ) )		; i+=1
						fwfit.write( "{:12f} ".format( dat[i] ) )	; i+=1
						for j in range(2) : 
							fwfit.write( "{:12f} {:12f}\t".format( dat[i],dat[i+1] ) ) ; i+=2
						fwfit.write( "\n" ) 


def axSetRange( ax, xystr, xyind ) :
	if xystr : 
		if xyind.find("y") > -1 :
			ax.set_ylim( float( xystr.split("_")[-2]) , float( xystr.split("_")[-1])  ) 
		if xyind.find("x") > -1 :
			ax.set_xlim( float( xystr.split("_")[-2]) , float( xystr.split("_")[-1])  ) 
		if xyind.find("z") > -1 :
			ax.set_zlim( float( xystr.split("_")[-2]) , float( xystr.split("_")[-1])  ) 
def axSetAxis( ax, xstr, ystr ) :
	ax.axis( [ float( xstr.split("_")[-2]) , float( xstr.split("_")[-1]) , float( ystr.split("_")[-2]) , float( ystr.split("_")[-1]) ] ) 
		

class datobj  :
	def __init__(self ) : pass
	def getGFgeneral(self, dat , nbasis, hobj , args , obj="glocaliw" ) :
		if( dat.find(".dat") > -1 ) :
			datfile = dat
		else : 
			if args.count :	count = args.count
			else :		count = hobj.count
			if args.finegrid : 
				ep=0.03
				nplot=1024
				datfile = "{}/lattice/u{:.3f}/greenimag_Nplot{}_ep{:.3f}_{}th.dat".format(dat, hobj.UF, nplot, ep , count )
				if hobj.ni>1 : 
					print "Need to modify this code to handle multi-impurity data."
					sys.exit(1)
			elif args.finegrid2 or args.finegrid3 or args.finegrid4 : 
				nplot=2048
				if args.finegrid3 : 
					nplot=2049
				elif args.finegrid4 : 
					nplot=2050
				datfile = "{}/ongoing/greenlociw_Nplot{}_{}th.dat".format(dat, nplot, count )
				if hobj.ni>1 : 
					print "Need to modify this code to handle multi-impurity data."
					sys.exit(1)
			else : 
				fhead	= "Glocal"
				datfile = "{}/ongoing/Glocal_u{:.2f}_{}th.txt".format(dat, hobj.UF, count )
				if obj == "greennew" :
					iimp	= 0
					datfile = "{}/ongoing/greennew{:d}_u{:.2f}_{}th.txt".format(dat, iimp, hobj.UF, count )
				elif obj == "glocalinverseiw" :
					fhead	= "Glocalinverse"
					datfile = "{}/ongoing/Glocalinverse_u{:.2f}_{}th.txt".format(dat, hobj.UF, count )
		print "Reading :", datfile
		rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		self.wn   = rawdatarr.transpose()[0]
		gfdat		= rawdatarr.transpose()[1:]
		gfRe		= gfdat[::2]
		gfIm		= gfdat[1::2]
		gfComp		= gfRe + gfIm * 1j
		print	"gfRe.shape : ", np.shape( gfRe )
		gfComp		= gfComp.reshape( hobj.Ni * hobj.tNC , hobj.Ni * hobj.tNC, -1 ) 
		print	"gfComp.shape : ", np.shape( gfComp )
		setattr( self, obj, gfComp ) 
	def getglocaliw(self, dat , nbasis, hobj , args , obj="glocaliw" ) :
		hobj.readEclustersoc()		#saved in hobj.ecluster[i][j]
		if( dat.find(".dat") > -1 ) :
			datfile = dat
		else : 
			if args.count :	count = args.count
			else :		count = hobj.count
			if args.finegrid : 
				ep=0.03
				nplot=1024
				datfile = "{}/lattice/u{:.3f}/greenimag_Nplot{}_ep{:.3f}_{}th.dat".format(dat, hobj.UF, nplot, ep , count )
				if hobj.ni>1 : 
					print "Need to modify this code to handle multi-impurity data."
					sys.exit(1)
			elif args.finegrid2 or args.finegrid3 or args.finegrid4 : 
				nplot=2048
				if args.finegrid3 : 
					nplot=2049
				elif args.finegrid4 : 
					nplot=2050
				datfile = "{}/ongoing/greenlociw_Nplot{}_{}th.dat".format(dat, nplot, count )
				if hobj.ni>1 : 
					print "Need to modify this code to handle multi-impurity data."
					sys.exit(1)
			else : 
				fhead	= "Glocal"
				datfile = "{}/ongoing/Glocal_u{:.2f}_{}th.txt".format(dat, hobj.UF, count )
				if obj == "greennew" :
					iimp	= 0
					datfile = "{}/ongoing/greennew{:d}_u{:.2f}_{}th.txt".format(dat, iimp, hobj.UF, count )
		print "Reading :", datfile
		rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		self.wn   = rawdatarr.transpose()[0]
		self.arrdat = rawdatarr.transpose()[1:]
		self.compMat = complexArrDatReturn( self.arrdat, nbasis ) 
	def getplotiw_lc(self, ax, dat , basis, nbasis, chdat ,hobj , args , lc , laboverw=False , dashes=(None,None) , logplot=False, T=False , lstylearr=False, labsuffix=False , obj='selfiw' ) :
		if obj.find("selfiw")>-1 :
			self.getselfiw( dat, nbasis, hobj, args );
			ylaborig = "Im$\Sigma(i\omega_n)$"
		elif obj.find("glocaliw")>-1 or obj.find("Glocaliw")>-1 or obj.find("glociw")>-1 or obj.find("Glociw")>-1 or obj.find("greennew")>-1 :
			self.getglocaliw( dat, nbasis, hobj, args , obj=obj );
			ylaborig = "Im$G_{loc}(i\omega_n)$"
		else : 
			print "ERROR :: Wrong argument for object with '{}'.".format(obj)
			sys.exit(1)
		self.gloc = self.arrdat
		nplot = len( self.wn )
		ncomp = len( self.gloc )
		print nplot, " points(lines)"
		print ncomp, " components(of matrix)"
		self.glocReDiag = range(nbasis)
		self.glocImDiag = range(nbasis)
		self.glocMat = self.compMat
		if args.zerosoc :
			basis = "t"
		print "0-th matrix :"
		for j1 in range(nbasis) :
			for j2 in range(nbasis) :
				print "{:20.8f}\t".format( self.glocMat[j1][j2][0] ),
			print ""
		if T is not False : 
			self.glocMat = transfMjefftot2g( self.glocMat, T ) 
			basis = args.transform
			print "0-th matrix (transformed) :"
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					print "{:20.8f}\t".format( self.glocMat[j1][j2][0] ),
				print ""

		Zarr = np.zeros( nbasis ) 
		Zavgtargetarr = np.zeros( nbasis ) 
		Gammatargetarr = np.zeros( nbasis ) 

		print "Entering getplotiw() ..."
		for j in chdat :
			self.glocReDiag[j] = self.glocMat[j][j].real
			self.glocImDiag[j] = self.glocMat[j][j].imag
			re = np.array( self.glocReDiag[j] )
			im = np.array( self.glocImDiag[j] )
			xdat = self.wn
			ydat = im
			if args.cutoff : 
				cutoff = float( args.cutoff ) 
				print "Cutoff : ", cutoff 
				ydat = ydat[ xdat > cutoff ] 
				xdat = xdat[ xdat > cutoff ] 
			if args.realpart : ydat = re
		
			if laboverw :	labelstr = laboverw
			else :		labelstr = bolabel(basis,j)
			if labsuffix :	labelstr = labsuffix+labelstr
			lt = [ "-" ]   + [ "-" ] *2 + [ "-" ]   +  [ ":" ] *2
			if basis.find("t")>-1 : 
				lt = [ "-" ]*2 + [ "-" ] *2 + [ ":" ] *2
			elif args.transform and args.transform("t")>-1 : 
				lt = [ "-" ]*2 + [ "-" ] *2 + [ ":" ] *2
			lstyle = lc+markers[j]+linestyles[j]
			if args.line : lstyle = lc+lt[j]
			if lstylearr : lstyle=lstylearr[j]
			ylab = ylaborig
			xlab = "$\omega_n$"
			if args.minus :
				ydat = -ydat
				ylab = "-"+ylab
			if args.minusx :
				xdat = -xdat
			if args.selffactor :
				ydat = ydat*float(args.selffactor)
				ylab = args.selffactor+"*"+ylab
			if args.selfoffset :
				print "OFFSET : ", float(args.selfoffset)
				ydat = ydat+float(args.selfoffset)
				ylab = ylab +"+"+ args.selfoffset
			if args.removelinear : 
				from scipy.optimize import curve_fit
				def funclin( x, a, b ) : return a*x + b
				typefit = "lin"
				formfunc= 'ax+b'
				nvarfunc= 2
				nfit = 3 
				if args.ndatalinearfit :
					nfit = int(args.ndatalinearfit)
				iinit = 0 
				if args.fitoffset : iinit += int(args.fitoffset)
				xraw = xdat[iinit:iinit+nfit] 
       				yraw = ydat[iinit:iinit+nfit] 
				if args.linearfitfromzero : 
					xraw = np.insert( xraw, 0, 0 )
					yraw = np.insert( yraw, 0, 0 )
				popt, pcov = curve_fit(funclin, xraw, yraw, bounds=([-99999,-1.], [999000., 1]))
				print "removing-FIT({};{}) : a={} b={}".format( formfunc, typefit, *popt )
				a = popt[0]
				devidat = xdat*a - ydat
				devidat = np.array( [xdat, devidat, ydat, xdat*a, devidat/xdat/a , xdat*0.+a ] ).transpose()
				np.savetxt(  "tmprmlinear_J{}_{}{}".format(hobj.J,basis,chdat[0]) , devidat ) 
			if args.difftwopoints : 
				xdum1= xdat[:-1]
				xdum2= xdat[1:]
				ydum1= ydat[:-1]
				ydum2= ydat[1:]
				xdat = xdat[:-1] 
				ydat = (ydum1-ydum2)/(xdum1-xdum2)
			if args.trans :
				if logplot :
					ax.loglog( ydat, xdat, lstyle , label=labelstr , linewidth=1.5 , dashes=dashes )
				else :
					ax.plot(   ydat, xdat, lstyle , label=labelstr , linewidth=1.5 , dashes=dashes )
				if j==0 and args.fillb :
					ax.fill( ydat, xdat, lstyle ,alpha=0.2 )
				ylab = ""
			else :
				if args.dashedline :
					if logplot :
						ax.loglog( xdat, ydat, lstyle , label=labelstr , dashes=dashes )
					else :
						ax.plot(   xdat, ydat, lstyle , label=labelstr , dashes=dashes )
				else : 
					if logplot :
						ax.loglog( xdat, ydat, lstyle , label=labelstr )
					else :
						ax.plot(   xdat, ydat, lstyle , label=labelstr )
			if len(chdat)<2 : ylab = ylab + "$^{{{}}}$".format(blabel(basis,j))
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
			if args.write : 
				dum = np.array( [xdat,ydat] ) 
				np.savetxt( "tmpplotdat", dum.transpose() ) 
			if args.drawabc : 
				def funcabc( x, a, b ,c ) : return a*x*x + b*x + c 
				nfit = 3 
				nresult = 100
				if args.ndatalinearfit :
					nfit = int(args.ndatalinearfit)
					if nfit > nresult  : 
						nresult = nfit*2
				iinit = 0 
				if args.fitoffset : iinit += int(args.fitoffset)
				xraw = xdat[iinit:iinit+nfit] 
       				yraw = ydat[iinit:iinit+nfit] 
				if logplot :
					xfit = np.logspace( np.log10(xraw[0]*0.1),  np.log10(xdat[-10]) , nresult ) 
				else : 
					xfit = np.linspace( xraw[0]*0.1,  xraw[-1]*2. , nresult ) 
				cba = args.drawabc.split("_")[-1:-4:-1]
				c,b,a = np.array( cba, dtype=float ) 
				carr = [c]
				print "ABC : ", a,b,c
				if args.carray : carr = np.array( args.carray.split("_"), dtype=float ) 
				for c in carr : 
					print "ABC : ", a,b,c
					yfit = funcabc(xfit, a,b,c ) 
					fitls = "b:"
					fitlab = None
					if args.carray : 
						fitls = "--"
						fitlab = '$c$={:.3g}'.format(c)
					if logplot :
						ax.loglog( xfit, yfit, fitls , lw=1 , label = fitlab )
					else :
						ax.plot( xfit, yfit, fitls , lw=1 , label = fitlab )
				if args.carray : 
				 	ax.legend()
			if args.linearfit or args.fitting or args.quadraticfit: 
				from scipy.optimize import curve_fit
				def funclin( x, a, b ) : return a*x 
				fitfunc = funclin
				typefit = "lin"
				formfunc= 'ax+b'
				nvarfunc= 2
				if args.fitting : 
					typefit = args.fitting
				elif args.quadraticfit : 
					typefit = 'quad'
				nfit = 3 
				nresult = 30
				print "Entering linearfit ..."
				if args.ndatalinearfit :
					nfit = int(args.ndatalinearfit)
					if nfit > nresult  : 
						nresult = nfit*2
					print "nfit : ", nfit 
				iinit = 0 
				if args.fitoffset : iinit += int(args.fitoffset)
				xraw = xdat[iinit:iinit+nfit] 
       				yraw = ydat[iinit:iinit+nfit] 
				if args.linearfitfromzero : 
					xraw = np.insert( xraw, 0, 0 )
					yraw = np.insert( yraw, 0, 0 )
				print "FIT_xraw : " , xraw[ [0,1,-2,-1] ]
				print "FIT_yraw : " , yraw[ [0,1,-2,-1] ]

				if typefit.find("poly3")>-1 :
					formfunc= 'ax^2+bx+c+d*x^3'
					nvarfunc= 4
					def funcpoly3( x, a, b, c, d ) : return a*x*x + b*x + c + d*x*x*x
					fitfunc = funcpoly3
				elif typefit.find("poly2")>-1 :
					formfunc= 'ax^2+bx+c'
					nvarfunc= 3
					def funcpoly2( x, a, b, c ) : return a*x*x + b*x + c
					fitfunc = funcpoly2
				elif typefit.find("fixlin")>-1 :
					formfunc= 'ax^2+b'
					nvarfunc= 2
					def funcfixlin( x, a, b ) : return a*x*x + b
					fitfunc = funcfixlin
				elif typefit.find("fixtemp")>-1 :
					formfunc= 'a(x^2-pi^2*T^2)'
					nvarfunc= 1
				elif typefit == "const" :
					formfunc= 'a'
					nvarfunc= 1
					def funcconst( x, a ) : return x*0.+a
					fitfunc = funcconst

				if typefit.find("fixlin")>-1 :
					def funclin( x, a, b ) : return a*x
					fitfunc = funclin
					nfitlin = 10
					xraw = xdat[0:nfitlin] 
       					yraw = ydat[0:nfitlin] 
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=(0., [999000., 1]))
					print "FIT({};{}) : ".format( formfunc, 'linear'+typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					xfit = np.linspace( xraw[0]*0.1,  xraw[-1]*2. , nresult ) 
					yfit = fitfunc(xfit,*popt) 
					fitls = "k:"
					if args.setlinear : 
						b = float( args.setlinear ) 
					else :
						if logplot :
							ax.loglog( xfit, yfit, fitls , lw=1 )
						else :
							ax.plot( xfit, yfit, fitls , lw=1 )
						b = popt[0]

					xraw = xdat[iinit:iinit+nfit] 
       					yraw = ydat[iinit:iinit+nfit] 
					def funcfixlin( x, a, c ) : return a*x*x + b*x + c
					fitfunc = funcfixlin
					initvarlowerarr = [-99999, 0 ] 
					initvarupperarr = [99999 , 1 ]
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=( initvarlowerarr , initvarupperarr ) )
					if args.fitoffset : 
						if logplot : 
							xfit = np.logspace( np.log10(xdat[0]),  np.log10(xraw[-1]) , num=nresult ) 
						else : 
							xfit = np.linspace( xdat[0],  xraw[-1] , nresult ) 
					else : 
						if logplot : 
							xfit = np.logspace( np.log10(xraw[0]),  np.log10(xraw[-1]) , num=nresult ) 
							xfit3= np.logspace( np.log10(xraw[0]/1e3),  np.log10(xraw[0]) , num=nresult ) 
						else : 
							xfit = np.linspace( xraw[0],  xraw[-1] , nresult ) 
					yfit = fitfunc(xfit,*popt) 
					yfit3= fitfunc(xfit3,*popt) 
					print "XFIT : ", xfit[ [0,1,-2-1] ]
					print "YFIT : ", yfit[ [0,1,-2-1] ]
					print "FIT({}+{}x;{}) : ".format( formfunc, b, 'quad'+typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					if args.minus :
						for jj in range(nvarfunc) :
							popt[jj] = -popt[jj]
					def returnT( a,c ) :
						T = np.sqrt( - c / a / np.pi/np.pi ) 
						return T, T*11605., 1./T
					print "\tT,T[K],beta : ", "{}".format( returnT(*popt) )
					fitls = "k--"
					fitls3= "r--"
					if logplot :
						ax.loglog( xfit, yfit, fitls , lw=1 )
						ax.loglog( xfit3, yfit3, fitls3 , lw=1 )
					else :
						ax.plot( xfit, yfit, fitls , lw=1 )
					with open( "tmpfitfix{}".format(j), 'a' ) as fwfit :
						outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT ] 
						for it in outarr : 
							fwfit.write( "{}\t".format(it) )
						fwfit.write( "{} {} {} {} {}\n".format( nfit, popt[0], b , popt[1], returnT(*popt) ) )
				elif typefit.find("fixtemp")>-1 :
					def funclin( x, a ) : return a*x
					fitfunc = funclin
					nfitlin = 10
					xraw = xdat[0:nfitlin] 
       					yraw = ydat[0:nfitlin] 
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=([0.], [999000.]) )
					print "FIT({};{}) : ".format( formfunc, 'linear'+typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					if hobj.IT :
						wtemp = np.pi/hobj.IT
					else :
						wtemp = np.pi/hobj.beta
					tempcutind = 0
					for ii in range(len(xdat)) : 
						if wtemp < xdat[ii] :
							tempcutind = ii
							break
					print "Ind_cuttemp : ", tempcutind , wtemp, xdat[tempcutind], xdat[tempcutind-1]
					xfitmax = xraw[-1]*2.
					if args.fitoffset : 
						xfitmax = xdat[tempcutind-1]
					xfit = np.linspace( xraw[0]*0.1, xfitmax , nresult ) 
					yfit = fitfunc(xfit,*popt) 
					fitls = "k:"
					if args.setlinear : 
						b = float( args.setlinear ) 
					else :
						if logplot :
							ax.loglog( xfit, yfit, fitls , lw=1 )
						else :
							ax.plot( xfit, yfit, fitls , lw=1 )
						b = popt[0]

					xraw = xdat[iinit:iinit+nfit] 
       					yraw = ydat[iinit:iinit+nfit] 
					if hobj.IT :
						temp = 1./hobj.IT
					else :
						temp = 0.
					piTsq = np.pi*np.pi*temp*temp
					def funcfixtemp( x, a ) : return a*( x*x - piTsq) + b*x 
					fitfunc = funcfixtemp
					initvarlowerarr = [-99999]
					initvarupperarr = [99999 ]
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=( initvarlowerarr , initvarupperarr ) )
					a = popt[0]
					if args.fitoffset : 
						if logplot : 
							xfit2= np.logspace( np.log10(xdat[0]),  np.log10(xraw[0]) , num=nresult ) 
							xfit = np.logspace( np.log10(xraw[0]),  np.log10(xraw[-1]) , num=nresult ) 
							xfit3= np.logspace( np.log10(xraw[0]/1e3),  np.log10(xraw[0]) , num=nresult ) 
						else : 
							xfit2= np.linspace( xdat[0],  xraw[0] , nresult ) 
							xfit = np.linspace( xraw[0],  xraw[-1] , nresult ) 
					else : 
						if logplot : 
							xfit = np.logspace( np.log10(xraw[0]),  np.log10(xraw[-1]) , num=nresult ) 
							xfit3= np.logspace( np.log10(xraw[0]/1e3),  np.log10(xraw[0]) , num=nresult ) 
						else : 
							xfit = np.linspace( xraw[0],  xraw[-1] , nresult ) 
							xfit3= np.linspace( xraw[0]/1e3,  xraw[0] , nresult ) 
					if args.fitoffset : 
						yfit2= fitfunc(xfit2,popt) 
						print "XFIT2: ", xfit2[ [0,1,-2-1] ]
						print "YFIT2: ", yfit2[ [0,1,-2-1] ]
					yfit = fitfunc(xfit,popt) 
					yfit3= fitfunc(xfit3,popt) 
					print "XFIT : ", xfit[ [0,1,-2-1] ]
					print "YFIT : ", yfit[ [0,1,-2-1] ]
					print "FIT({}+{}x;{}) : ".format( formfunc, b, 'quad'+typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					print "(pi*T)^2, -a*piTsq, -a*pi*pi*T : ", piTsq, -a*piTsq , -a*np.pi*np.pi*temp
					fitls = "k--"
					fitls2= "r--"
					fitls3= "b--"
					if logplot :
						if args.fitoffset : 
							ax.loglog( xfit,  yfit,  fitls  , lw=1 )
							ax.loglog( xfit2, yfit2, fitls2 , lw=1 )
						else : 	
							ax.loglog( xfit, yfit, fitls , lw=1 )
							ax.loglog( xfit3, yfit3, fitls3, lw=1 )
					else :
						if args.fitoffset : 
							ax.plot( xfit,  yfit,  fitls , lw=1 )
							ax.plot( xfit2, yfit2, fitls2, lw=1 )
						else : 
							ax.plot( xfit, yfit, fitls , lw=1 )
							ax.plot( xfit3, yfit3, fitls3, lw=1 )
					with open( "tmpfitfixtemp{}".format(j), 'a' ) as fwfit :
						outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT ] 
						for it in outarr : 
							fwfit.write( "{}\t".format(it) )
						fwfit.write( "{} {} {}\n".format( nfit, popt[0], b  ) )
				elif (typefit.find("poly2")>-1) or (typefit.find("poly3")>-1) :
					if nvarfunc == 3 :
						initvararr = [-20]*(nvarfunc-1) + [-0.5]
						initvararr2= [ 20]*(nvarfunc-1) + [ 0.5]
					elif nvarfunc == 4 :
						initvararr = [-20]*(nvarfunc-2) + [-0.5] +[-20]
						initvararr2= [ 20]*(nvarfunc-2) + [ 0.5] +[ 20]
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=( initvararr , initvararr2 ) )
					if args.fitoffset : 
						if logplot : 
							xfit = np.logspace( np.log10(xdat[0]),  np.log10(xraw[-1]) , num=nresult ) 
						else : 
							xfit = np.linspace( xdat[0],  xraw[-1] , nresult ) 
					else : 
						if logplot : 
							xfit = np.logspace( np.log10(xraw[0]),  np.log10(xraw[-1]) , num=nresult ) 
						else : 
							xfit = np.linspace( xraw[0],  xraw[-1] , nresult ) 
					yfit = fitfunc(xfit,*popt) 
					print "XFIT : ", xfit[ [0,1,-2-1] ]
					print "YFIT : ", yfit[ [0,1,-2-1] ]
					print "FIT({};{}) : ".format( formfunc, typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					if args.minus :
						for jj in range(nvarfunc) :
							popt[jj] = -popt[jj]
					if nvarfunc==3 : 
						def mZCT( a,b,c ) :
							m = 1.-b
							Z = 1./m
							C = a /m
							T = np.sqrt( - c / a / np.pi/np.pi ) 
							return m, Z, C, T, T*11605., 1./T
						print "\tm,Z,C,T,T[K],beta : ", "{}".format( mZCT(*popt) )
					else : 
						def mZCT( a,b,c,d ) :
							m = 1.-b
							Z = 1./m
							C = a /m
							T = np.sqrt( - c / a / np.pi/np.pi ) 
							return m, Z, C, T, T*11605., 1./T
						print "\tm,Z,C,T,T[K],beta : ", "{}".format( mZCT(*popt) )
					fitls = "k--"
					if logplot :
						ax.loglog( xfit, yfit, fitls , lw=1 )
					else :
						ax.plot( xfit, yfit, fitls , lw=1 )
				elif typefit == "const" :
					initvararr = [99999]*nvarfunc
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=(-99999 , initvararr) )
					if args.fitoffset : 
						if logplot : 
							xfit = np.logspace( np.log10(xdat[0]),  np.log10(xraw[-1]) , num=nresult ) 
							xfit2= np.logspace( np.log10(xraw[-1]),  np.log10(xdat[-1]) , num=nresult ) 
						else : 
							xfit = np.linspace( xdat[0],  xraw[-1] , nresult ) 
					else : 
						if logplot : 
							xfit = np.logspace( np.log10(xraw[0]),  np.log10(xraw[-1]) , num=nresult ) 
							xfit2= np.logspace( np.log10(xraw[-1]),  np.log10(xdat[-1]) , num=nresult ) 
						else : 
							xfit = np.linspace( xraw[0],  xraw[-1] , nresult ) 
					yfit = fitfunc(xfit,*popt) 
					yfit2= fitfunc(xfit2,*popt) 
					print "XFIT : ", xfit[ [0,1,-2-1] ]
					print "YFIT : ", yfit[ [0,1,-2-1] ]
					print "FIT({};{}) : ".format( formfunc, typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					if args.minus :
						for jj in range(nvarfunc) :
							popt[jj] = -popt[jj]
					fitls = "k--"
					if logplot :
						ax.loglog( xfit, yfit, fitls , lw=1 )
					else :
						ax.plot( xfit, yfit, fitls , lw=1 )
					fitls2= "r--"
					if logplot :
						ax.loglog( xfit,  yfit,  fitls  , lw=1 )
						ax.loglog( xfit2, yfit2, fitls2 , lw=1 )
					else :
						ax.plot( xfit,  yfit,  fitls , lw=1 )
						ax.plot( xfit2, yfit2, fitls2, lw=1 )
				elif typefit.find("lin")>-1 :
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=(0., [999000., 1]))
					xfit = np.linspace( xdat[:nfit][0]*0.1,  1 , nresult ) 
					yfit = fitfunc(xfit,*popt) 
					print "FIT(a*x) : a={}".format( popt[0] )
					Zftn = lambda a,b : 1./(1.+a)
					print "Z : {} {} {}".format( j, hobj.J, Zftn(*popt) ) 
					print "m* = 1/Z : {} {} {} {}".format( j, hobj.J, 1./Zftn(*popt), hobj.D ) , "\n"
					fitls = "k--"
					Zarr[j] = Zftn(*popt)
					if logplot :
						ax.loglog( xfit, yfit, fitls , lw=1 )
					else :
						ax.plot( xfit, yfit, fitls , lw=1 )

					with open( "tmplinearfitavg{}".format(j), 'w' ) as fwfitavg :
						for ii in range(30) : 
							ydum = ydat[ii*5]
							xdum = xdat[ii*5]
							steepness = (0.-ydum) / (0.-xdum)
							Zavg = 1./(1.+steepness)
							Gamma= Zavg * ydum	#Inverse quasiparticle lifetime (Mravlje201!)
							outarr = [ xdum, ydum, steepness, 1./Zavg, Zavg, xdum*11605. , Gamma ] 
							for it in outarr : 
								fwfitavg.write( "{:8f}\t".format(it) )
							fwfitavg.write( "\n" )

					ydum = ydat[39]
					xdum = xdat[39]
					tempdum = xdum*11605. 
					steepness = (0.-ydum) / (0.-xdum)
					Zavgtargetarr[j] = 1./(1.+steepness)
					Gammatargetarr[j] = ydum * Zavgtargetarr[j]

		if args.linearfit : 
			with open( "tmplinearfitdat", 'a' ) as fwfit :
				for j in chdat : 
					fwfit.write( "{} ".format(j) )
				for j in chdat : 
					fwfit.write( "{:12f} ".format(1./Zarr[j]) )
				fwfit.write( "\t" ) 
				outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT ] 
				for it in outarr : 
					fwfit.write( "{}\t".format(it) )
				fwfit.write( "\n" ) 
			with open( "tmplinearfitavg_b512", 'a' ) as fwfitavg :
				for j in chdat : 
					fwfitavg.write( "{} ".format(j) )
				for j in chdat : 
					fwfitavg.write( "{:12f} ".format(1./Zavgtargetarr[j]) )
				outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra, hobj.IT , tempdum ] 
				for it in outarr : 
					try : 
						fwfitavg.write( "{:8f}\t".format(it) )
					except : 
						fwfitavg.write( "{}\t".format(it) )
				fwfitavg.write( "\n" )
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		#ax.axhline( 0, 0,1 , color='g', dashes=dashpt , lw=0.5 )
		if args.datyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.datyrange ) 
			else :
				axSetRange( ax , args.datyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 
		if args.datxrange : 
			axSetRange( ax , args.datxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])
		if args.write : 
			dum = np.array( [ xdat, ydat] ).transpose()
			np.savetxt("tmpdata", dum ) 
		if args.writegrid : 
			dumn = int(args.writegrid) 
			dum = np.array( [ xdat[:dumn], ydat[:dumn]] ).transpose()
			np.savetxt("tmpdatagrid_n{}_J{:.2f}".format(dumn,hobj.J), dum ) 
	def getplotiwoffdiag_lc(self, ax, dat , basis, nbasis, chdat ,hobj , args , lc , T=False , obj='selfiw' ) :
		if obj.find("selfiw")>-1 :
			self.getselfiw( dat, nbasis, hobj, args );
			ylaborig = "Im$\Sigma(i\omega_n)$"
		elif obj.find("glocaliw")>-1 or obj.find("Glocaliw")>-1 or obj.find("glociw")>-1 or obj.find("Glociw")>-1 :
			self.getglocaliw( dat, nbasis, hobj, args );
			ylaborig = "Im$G_{loc}(i\omega_n)$"
		else : 
			print "ERROR :: Wrong argument for object with '{}'.".format(obj)
			sys.exit(1)
		self.gloc = self.arrdat
		nplot = len( self.wn )
		ncomp = len( self.gloc )
		print nplot, " points(lines)"
		print ncomp, " components(of matrix)"
		self.glocReOffDiag = [ range(nbasis) for j in range(nbasis) ]
		self.glocImOffDiag = [ range(nbasis) for j in range(nbasis) ]
		self.glocMat = self.compMat
		print "0-th matrix :"
		for j1 in range(nbasis) :
			for j2 in range(nbasis) :
				print "{:20.8f}\t".format( self.glocMat[j1][j2][0] ),
			print ""
		if T is not False : 
			self.glocMat = transfMjefftot2g( self.glocMat, T ) 
			basis = args.transform
			print "0-th matrix (transformed) :"
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					print "{:20.8f}\t".format( self.glocMat[j1][j2][0] ),
				print ""
		
		for k in chdat :
			i , j = [ k[0], k[1] ]
			self.glocReOffDiag[i][j] = self.glocMat[i][j].real
			self.glocImOffDiag[i][j] = self.glocMat[i][j].imag
			re = np.array( self.glocReOffDiag[i][j] ) 
			im = np.array( self.glocImOffDiag[i][j] ) 
			xdat = self.wn
			ydat = im
			if args.realpart : ydat = re

			ls = [ "-" ] + [ "-." ] + [ "--" ] *3 +  [ "-" ] + [ "-." ] + [ "--" ] *3  
			ms = [ "x" ] + [ "3" ] + [ ">" ] + [ "^" ] + [ "x" ] + [ "+" ] 
			ylab = "Im$\Delta(i\omega_n)$"
			xlab = "$\omega_n$"
			if args.minus :
				ydat = -ydat
				ylab = "-"+ylab
			if args.minusx :
				xdat = -xdat
			if args.selffactor :
				ydat = ydat*float(args.selffactor)
				ylab = args.selffactor+"*"+ylab
			if args.selfoffset :
				ydat = ydat+float(args.selfoffset)
				ylab = ylab +"+"+ args.selfoffset
			if args.trans :
				ax.plot( ydat, xdat, ms[i]+ls[j], label="{}{}".format(i,j)  , linewidth=1.5 )
				if i==0 and args.fillb :
					ax.fill( ydat, xdat, ms[i]+ls[j] ,alpha=0.2 )
				tmplab = xlab
				xlab   = ylab
				ylab   = tmplab
			else :
				ax.plot( xdat, ydat, ms[i]+ls[j], label="{}{}".format(i,j) )
				#ax.legend()
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		#ax.axhline( 0, 0,1 , color='g', dashes=dashpt , lw=0.5 )
		if args.datyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.datyrange ) 
			else :
				axSetRange( ax , args.datyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 
		if args.datxrange : 
			axSetRange( ax , args.datxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])
	def getplotiwoffdiagratio_lc(self, ax, dat , basis, nbasis, chdat ,hobj , args , lc , T=False , obj='selfiw' ) :
		if obj.find("selfiw")>-1 :
			self.getselfiw( dat, nbasis, hobj, args );
			ylaborig = "Im$\Sigma(i\omega_n)$"
		elif obj.find("glocaliw")>-1 or obj.find("Glocaliw")>-1 or obj.find("glociw")>-1 or obj.find("Glociw")>-1 :
			self.getglocaliw( dat, nbasis, hobj, args );
			ylaborig = "Im$G_{loc}(i\omega_n)$"
		else : 
			print "ERROR :: Wrong argument for object with '{}'.".format(obj)
			sys.exit(1)
		self.gloc = self.arrdat
		nplot = len( self.wn )
		ncomp = len( self.gloc )
		print nplot, " points(lines)"
		print ncomp, " components(of matrix)"
		self.glocReOffDiag = [ range(nbasis) for j in range(nbasis) ]
		self.glocImOffDiag = [ range(nbasis) for j in range(nbasis) ]
		self.glocMat = self.compMat
		print "0-th matrix :"
		for j1 in range(nbasis) :
			for j2 in range(nbasis) :
				print "{:20.8f}\t".format( self.glocMat[j1][j2][0] ),
			print ""
		if T is not False : 
			self.glocMat = transfMjefftot2g( self.glocMat, T ) 
			basis = args.transform
			print "0-th matrix (transformed) :"
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					print "{:20.8f}\t".format( self.glocMat[j1][j2][0] ),
				print ""
		
		for k in chdat :
			i , j = [ k[0], k[1] ]
			self.glocReOffDiag[i][j] = self.glocMat[i][j].real
			self.glocImOffDiag[i][j] = self.glocMat[i][j].imag
			re = np.array( self.glocReOffDiag[i][j] ) 
			im = np.array( self.glocImOffDiag[i][j] ) 
			rerat = re /  self.glocMat[i][i].real
			imrat = im /  self.glocMat[i][i].imag
			xdat = self.wn
			ydat = imrat
			if args.realpart : ydat = rerat

			ls = [ "-" ] + [ "-." ] + [ "--" ] *3 +  [ "-" ] + [ "-." ] + [ "--" ] *3  
			ms = [ "x" ] + [ "3" ] + [ ">" ] + [ "^" ] + [ "x" ] + [ "+" ] 
			ylab = "Im$\Delta(i\omega_n)$"
			xlab = "$\omega_n$"
			if args.minus :
				ydat = -ydat
				ylab = "-"+ylab
			if args.minusx :
				xdat = -xdat
			if args.selffactor :
				ydat = ydat*float(args.selffactor)
				ylab = args.selffactor+"*"+ylab
			if args.selfoffset :
				ydat = ydat+float(args.selfoffset)
				ylab = ylab +"+"+ args.selfoffset
			if args.trans :
				ax.plot( ydat, xdat, ms[i]+ls[j], label="{}{}/{}{}".format(i,j,i,i)  , linewidth=1.5 )
				if i==0 and args.fillb :
					ax.fill( ydat, xdat, ms[i]+ls[j] ,alpha=0.2 )
				tmplab = xlab
				xlab   = ylab
				ylab   = tmplab
			else :
				ax.plot( xdat, ydat, ms[i]+ls[j], label="{}{}/{}{}".format(i,j,i,i) )
				#ax.legend()
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		#ax.axhline( 0, 0,1 , color='g', dashes=dashpt , lw=0.5 )
		if args.datyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.datyrange ) 
			else :
				axSetRange( ax , args.datyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 
		if args.datxrange : 
			axSetRange( ax , args.datxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])
	def getselfw(self, dat , basis, nbasis, hobj , args ) :
		if args.count :	count = int(args.count)
		else :		count = int(hobj.count)
		if args.epsilon :	epsilon = float(args.epsilon)
		else :			epsilon = float(hobj.epsilon)
		if( dat.find(".dat") > -1 ) :
			datfile = dat
		else : 
		 	if basis.find("t") > -1  : 
		 		if basis.find("z")>-1 : 
					datfile = "{}/lattice/u{:.3f}/selfreal_Nplot{}_ep{:.3f}_{}th.dat".format(dat,hobj.UF,hobj.Nplot,epsilon,count)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
				else : 
					datfile = "{}/lattice/tselfreal_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
				#try : 
				#except :
				#	datfile = "{}/lattice/u{:.3f}/tselfreal_ep{:.3f}_{}th.dat".format(dat,hobj.UF,epsilon,count)
				#	rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		 	elif basis.find("j") > -1  : 
				try : 
					datfile = "{}/lattice/selfreal_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
				except :
					datfile = "{}/lattice/u{:.3f}/selfreal_Nplot{}_ep{:.3f}_{}th.dat".format(dat,hobj.UF,hobj.Nplot,epsilon,count)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		print "Reading :", datfile
		self.wn   = rawdatarr.transpose()[0]
		self.arrdat = rawdatarr.transpose()[1:]
		self.compMat = complexArrDatReturn( self.arrdat, nbasis ) 
	def getglocw(self, dat , basis, nbasis, hobj , args ) :
		if args.count :	count = int(args.count)
		else :		count = int(hobj.count)
		if args.epsilon :	epsilon = float(args.epsilon)
		else :			epsilon = float(hobj.epsilon)
		datname = "ldos"
		if basis.find("t") > -1 :
			datname = "tldos"
		try : 
			datfile = "{}/lattice/{}_u{:.2f}_Nplot{}.dat".format(dat, datname, hobj.UF, hobj.Nplot )
			rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		except :
			datfile = "{}/lattice/u{:.3f}/{}_Nplot{}_ep{:.3f}_{}th.dat".format(dat,hobj.UF,datname,hobj.Nplot,epsilon,count)
			rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		print "Reading :", datfile
		self.wn   = rawdatarr.transpose()[0]
		self.arrdat = rawdatarr.transpose()[1:]
		self.compMat = complexArrDatReturn( self.arrdat, nbasis ) 
	def getplotw(self, ax, dat , basis, nbasis, chdat ,hobj , args , idir=0 , laboverw=None, dashes=(None,None) , T=False , realpart=False , obj="selfw" ) :
		if obj.find("selfw")>-1 :
			self.getselfw(dat , basis, nbasis, hobj , args )
			ylaborig = "Im$\Sigma(i\omega_n)$"
		elif obj.find("glocalw")>-1 or obj.find("Glocalw")>-1 or obj.find("glocw")>-1 or obj.find("Glocw")>-1 :
			self.getglocw(dat , basis, nbasis, hobj , args )
			ylaborig = "Im$G_{loc}(i\omega_n)$"
		else : 
			print "ERROR :: Wrong argument for object with '{}'.".format(obj)
			sys.exit(1)
		self.gloc = self.arrdat
		self.glocReDiag , self.glocImDiag = matrixDiagReturnReIm( self.gloc , nbasis )
		print len( self.wn ), "points"
		print len( self.gloc ), "components"
		try : 
			if args.sum : 
				xDatSum = self.wn
				yDatSum = np.zeros( len(self.wn) ) 
				labsum  = "Sum["
		except : 
			setattr( args, "sum" , None )  
		Zarr = np.zeros(nbasis)
		effmasserrarr = np.zeros(nbasis)
		GammaperZarr = np.zeros(nbasis)
		for j in chdat :
			self.glocReDiag[j] = self.gloc[2*j   +2*nbasis*j]
			self.glocImDiag[j] = self.gloc[2*j+1 +2*nbasis*j]
			re = np.array( self.glocReDiag[j] ,dtype=float )
			im = np.array( self.glocImDiag[j] ,dtype=float )
			if args.sinv : 
				re = 1./re
				im = 1./im
		
			lc = [ "-" ]*2 + [ "-" ] *2 + [ ":" ] *2
			xdat = self.wn
			ydat = im
			xlab = "$\omega$"
			ylab = "Im$\Sigma(\omega)$"
			if args.realpart :
				ydat = re
				ylab = "Re$\Sigma(\omega)$"
			if args.minus :
				ydat = -ydat
				ylab = "-"+ylab
			if args.selffactor :
				ydat = ydat*float(args.selffactor)
				ylab = args.selffactor+"*"+ylab
			if args.selfoffset :
				ydat = ydat+float(args.selfoffset)
				ylab = ylab +"+"+ args.selfoffset
			if args.sum : 
				yDatSum = yDatSum + ydat 
			if laboverw :	labelstr = laboverw
			else :		labelstr = bolabel(basis,j)
			if args.solidline :	plotbasis = "jsolid"+basis
			else :			plotbasis = basis
			if args.trans :
				xlab = ylab
				ylab = ""
				ydatdum = ydat
				xdat = ydatdum
				ydat = xdat
			if args.selfm :
				ydat = np.array(-ydat , dtype=float) 
				ylab = "-"+ylab
			if args.fillb :
				if args.trans : 
					if j==0 :
						ax.fill_betweenx( self.wn, 0, ydat, facecolor='C%d'%j ,alpha=0.3 )
					if j==1 :
						ax.fill_betweenx( self.wn, 0, ydat, facecolor='C%d'%idir ,alpha=0.3 )
				else :
					if j==0 :
						ax.fill_between( self.wn, ydat, facecolor='C%d'%j ,alpha=0.3 )
			if args.fill :
				if j==1 :
					ax.fill( self.wn, ydat, 'C%d'%j ,alpha=0.3 ) 
			if args.sum : 
				labsum = labsum + labelstr + ','
			else : 
				ax.plot( xdat, ydat, blt(plotbasis,j) , label=labelstr , linewidth=1.5 , dashes=dashes )
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )

			if args.linearfit or args.quadraticfit or args.fitting : 
				from scipy.optimize import curve_fit
				def funclin( x, a, b ) : return a*x + b
				fitfunc = funclin
				typefit = "lin"
				formfunc= 'ax+b'
				nvarfunc= 2
				if args.quadraticfit : 
					typefit="quad"
					formfunc= 'ax^2+b'
					nvarfunc= 2
				if args.fitting :
					typefit=args.fitting
				if typefit.find("quad")>-1 : 
					if args.zerotempquadraticfit : 
						formfunc= 'ax^2'
						nvarfunc= 2
						def funclin( x, a, b ) : return a*x*x 
					else : 
						def funclin( x, a, b ) : return a*x*x + b
				elif typefit.find("poly2")>-1 : 
					formfunc= 'ax^2+bx+c'
					nvarfunc= 2
					def funclin( x, a, b, c ) : return a*x*x + b*x + c
				nfit = 10 
				nmiddle = hobj.Nplot/2
				if args.fitoffset : nmiddle += int(args.fitoffset)
				print "FIT from xdat, idat : ", xdat[nmiddle], nmiddle
				if args.ndatalinearfit : nfit = int(args.ndatalinearfit)
				xraw = xdat[nmiddle-nfit:nmiddle+nfit] 
       				yraw = ydat[nmiddle-nfit:nmiddle+nfit] 
				initvar = [999]*nvarfunc
				popt, pcov = curve_fit(funclin, xraw, yraw, bounds=(-999., initvar ))
				#zfit = np.polyfit( xraw , yraw , 1 )
				#ffit = np.poly1d( zfit ) 
				print "FIT_xraw : " , xraw
				print "FIT_yraw : " , yraw
				xfit = np.linspace( xraw[0],xraw[-1] , 30 ) 
				yfit = funclin(xfit,*popt) 
				fitls = "k--"
				ax.plot( xfit, yfit, fitls , lw=1 )
				if typefit.find("quad")>-1 :
					print "FIT({};{}) : a={}, b={}".format( formfunc, typefit,  *popt )
					print "b ~ T^2 : {} {} {} {}".format( j, hobj.J, popt[1], hobj.D ) , "\n"
					GammaperZarr[j] = popt[1]
				elif typefit.find("lin")>-1 :	
					print "FIT(a*x+b;{}) : a={}, b={}".format( typefit, *popt )
					Zftn = lambda a,b : 1./(1.-a)
					print "Z : {} {} {}".format( j, hobj.J, Zftn(*popt) ) 
					print "m* = 1/Z : {} {} {} {}".format( j, hobj.J, 1./Zftn(*popt), hobj.D ) , "\n"
					Zarr[j] = Zftn(*popt)
					effmasserrarr[j] = np.sqrt( np.diag(pcov)[0] )

					with open( "tmplinearfitRe{}_{}_b{}".format(j,basis,hobj.beta), 'w' ) as fwfitavg :
						for ii in range(30) : 
							ydum = ydat[ii*5]
							xdum = xdat[ii*5]
							steepness = (0.-ydum) / (0.-xdum)
							Zavg = 1./(1.+steepness)
							Gamma= Zavg * ydum	#Inverse quasiparticle lifetime (Mravlje201!)
							outarr = [ xdum, ydum, steepness, 1./Zavg, Zavg, xdum*11605. , Gamma ] 
							for it in outarr : 
								fwfitavg.write( "{:8f}\t".format(it) )
							fwfitavg.write( "\n" )
				else : 
					print "ERROR :: something wrong in fitting , " + typefit 
					sys.exit(1)
		if args.sum : 
			xdat = xDatSum
			ydat = yDatSum
			if args.trans :
				ydatdum = ydat
				xdat = ydatdum
				ydat = xdat
			ax.plot( xdat, ydat, "ro-" , label=labsum+"]" , linewidth=1.5 , dashes=dashes )
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		#ax.axhline( 0, 0,1 , color='g', dashes=dashpt , lw=0.5 )
		if args.realpart : 
			pass
		else : 
			ax.axhline( 0 , color='gray', linestyle="--", lw=0.5 ) 
		if args.selfyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.selfyrange ) 
			else :
				axSetRange( ax , args.selfyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 
		if args.selfxrange : 
			axSetRange( ax , args.selfxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])

		if args.linearfit or args.quadraticfit : 
			if args.quadraticfit : 
				foutname = "tmpquadfitdat"
			else :
				foutname = "tmplinearfitRedat_{}_J{:.2f}_b{}".format(basis,hobj.J,hobj.beta)
			with open( foutname, 'a' ) as fwfit :
				for j in chdat : 
					fwfit.write( "{} ".format(j) )
				for j in chdat : 
					mass = 1./Zarr[j] 
					if args.quadraticfit : mass = GammaperZarr[j]
					fwfit.write( "{:12f} ".format(mass) )
				for j in chdat : 
					fwfit.write( "{:12f} ".format(effmasserrarr[j]) )
				fwfit.write( "\t" ) 
				outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT , nfit , xraw[0],xraw[-1] ] 
				for it in outarr : 
					fwfit.write( "{}\t".format(it) )
				fwfit.write( "\n" ) 

	def getplotselfw2(self, ax, dat , basis, nbasis, chdat ,hobj , args , idir=0 , laboverw=None, dashes=(None,None) , T=False , realpart=False ) :
		self.getselfw(dat , basis, nbasis, hobj , args )
		self.gloc = self.arrdat
		self.glocReDiag , self.glocImDiag = matrixDiagReturnReIm( self.gloc , nbasis )
		print len( self.wn ), "points"
		print len( self.gloc ), "components"
		try : 
			if args.sum : 
				xDatSum = self.wn
				yDatSum = np.zeros( len(self.wn) ) 
				labsum  = "Sum["
		except : 
			setattr( args, "sum" , None )  
		Zarr = np.zeros(nbasis)
		effmasserrarr = np.zeros(nbasis)
		GammaperZarr = np.zeros(nbasis)
		for j in chdat :
			self.glocReDiag[j] = self.gloc[2*j   +2*nbasis*j]
			self.glocImDiag[j] = self.gloc[2*j+1 +2*nbasis*j]
			re = np.array( self.glocReDiag[j] ,dtype=float )
			im = np.array( self.glocImDiag[j] ,dtype=float )
			if args.sinv : 
				re = 1./re
				im = 1./im
		
			lc = [ "-" ]*2 + [ "-" ] *2 + [ ":" ] *2
			xdat = self.wn
			ydat = im
			xlab = "$\omega$"
			ylab = "Im$\Sigma(\omega)$"
			if args.realpart :
				ydat = re
				ylab = "Re$\Sigma(\omega)$"
			if args.minus :
				ydat = -ydat
				ylab = "-"+ylab
			if args.selffactor :
				ydat = ydat*float(args.selffactor)
				ylab = args.selffactor+"*"+ylab
			if args.selfoffset :
				ydat = ydat+float(args.selfoffset)
				ylab = ylab +"+"+ args.selfoffset
			if args.sum : 
				yDatSum = yDatSum + ydat 
			if laboverw :	labelstr = laboverw
			else :		labelstr = bolabelArr(basis)
			if args.solidline :	plotbasis = "jsolid"+basis
			else :			plotbasis = basis
			if args.trans :
				xlab = ylab
				ylab = ""
				ydatdum = ydat
				xdat = ydatdum
				ydat = xdat
			if args.selfm :
				ydat = np.array(-ydat , dtype=float) 
				ylab = "-"+ylab
			if args.fillb :
				if args.trans : 
					if j==0 :
						ax.fill_betweenx( self.wn, 0, ydat, facecolor='C%d'%j ,alpha=0.3 )
					if j==1 :
						ax.fill_betweenx( self.wn, 0, ydat, facecolor='C%d'%idir ,alpha=0.3 )
				else :
					if j==0 :
						ax.fill_between( self.wn, ydat, facecolor='C%d'%j ,alpha=0.3 )
			if args.fill :
				if j==1 :
					ax.fill( self.wn, ydat, 'C%d'%j ,alpha=0.3 ) 
			if args.sum : 
				labsum = labsum + labelstr[j] + ','
			else : 
				print "PBASIS : ", blt(plotbasis,j), plotbasis,j 
				if dashes[0] is None : 
					ax.plot( xdat, ydat, blt(plotbasis,j) , label=labelstr[j] , linewidth=1.5 )
				else : 
					ax.plot( xdat, ydat, blt(plotbasis,j) , label=labelstr[j] , linewidth=1.5 , dashes=dashes )
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )

			if args.linearfit or args.quadraticfit or args.fitting : 
				from scipy.optimize import curve_fit
				def funclin( x, a, b ) : return a*x + b
				fitfunc = funclin
				typefit = "lin"
				formfunc= 'ax+b'
				nvarfunc= 2
				if args.quadraticfit : 
					typefit="quad"
					formfunc= 'ax^2+b'
					nvarfunc= 2
				if args.fitting :
					typefit=args.fitting
				if typefit.find("quad")>-1 : 
					if args.zerotempquadraticfit : 
						formfunc= 'ax^2'
						nvarfunc= 2
						def funclin( x, a, b ) : return a*x*x 
					else : 
						def funclin( x, a, b ) : return a*x*x + b
				elif typefit.find("poly2")>-1 : 
					formfunc= 'ax^2+bx+c'
					nvarfunc= 2
					def funclin( x, a, b, c ) : return a*x*x + b*x + c
				nfit = 10 
				nmiddle = hobj.Nplot/2
				if args.fitoffset : nmiddle += int(args.fitoffset)
				print "FIT from xdat, idat : ", xdat[nmiddle], nmiddle
				if args.ndatalinearfit : nfit = int(args.ndatalinearfit)
				xraw = xdat[nmiddle-nfit:nmiddle+nfit] 
       				yraw = ydat[nmiddle-nfit:nmiddle+nfit] 
				initvar = [999]*nvarfunc
				popt, pcov = curve_fit(funclin, xraw, yraw, bounds=(-999., initvar ))
				#zfit = np.polyfit( xraw , yraw , 1 )
				#ffit = np.poly1d( zfit ) 
				print "FIT_xraw : " , xraw
				print "FIT_yraw : " , yraw
				xfit = np.linspace( xraw[0],xraw[-1] , 30 ) 
				yfit = funclin(xfit,*popt) 
				fitls = "k--"
				ax.plot( xfit, yfit, fitls , lw=1 )
				if typefit.find("quad")>-1 :
					print "FIT({};{}) : a={}, b={}".format( formfunc, typefit,  *popt )
					print "b ~ T^2 : {} {} {} {}".format( j, hobj.J, popt[1], hobj.D ) , "\n"
					GammaperZarr[j] = popt[1]
				elif typefit.find("lin")>-1 :	
					print "FIT(a*x+b;{}) : a={}, b={}".format( typefit, *popt )
					Zftn = lambda a,b : 1./(1.-a)
					print "Z : {} {} {}".format( j, hobj.J, Zftn(*popt) ) 
					print "m* = 1/Z : {} {} {} {}".format( j, hobj.J, 1./Zftn(*popt), hobj.D ) , "\n"
					Zarr[j] = Zftn(*popt)
					effmasserrarr[j] = np.sqrt( np.diag(pcov)[0] )

					with open( "tmplinearfitRe{}_{}_b{}".format(j,basis,hobj.beta), 'w' ) as fwfitavg :
						for ii in range(30) : 
							ydum = ydat[ii*5]
							xdum = xdat[ii*5]
							steepness = (0.-ydum) / (0.-xdum)
							Zavg = 1./(1.+steepness)
							Gamma= Zavg * ydum	#Inverse quasiparticle lifetime (Mravlje201!)
							outarr = [ xdum, ydum, steepness, 1./Zavg, Zavg, xdum*11605. , Gamma ] 
							for it in outarr : 
								fwfitavg.write( "{:8f}\t".format(it) )
							fwfitavg.write( "\n" )
				else : 
					print "ERROR :: something wrong in fitting , " + typefit 
					sys.exit(1)
		if args.sum : 
			xdat = xDatSum
			ydat = yDatSum
			if args.trans :
				ydatdum = ydat
				xdat = ydatdum
				ydat = xdat
			ax.plot( xdat, ydat, "ro-" , label=labsum+"]" , linewidth=1.5 , dashes=dashes )
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		#ax.axhline( 0, 0,1 , color='g', dashes=dashpt , lw=0.5 )
		if args.realpart : 
			pass
		else : 
			ax.axhline( 0 , color='gray', linestyle="--", lw=0.5 ) 
		if args.selfyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.selfyrange ) 
			else :
				axSetRange( ax , args.selfyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 
		if args.selfxrange : 
			axSetRange( ax , args.selfxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])

		if args.linearfit or args.quadraticfit : 
			if args.quadraticfit : 
				foutname = "tmpquadfitdat"
			else :
				foutname = "tmplinearfitRedat_{}_J{:.2f}_b{}".format(basis,hobj.J,hobj.beta)
			with open( foutname, 'a' ) as fwfit :
				for j in chdat : 
					fwfit.write( "{} ".format(j) )
				for j in chdat : 
					mass = 1./Zarr[j] 
					if args.quadraticfit : mass = GammaperZarr[j]
					fwfit.write( "{:12f} ".format(mass) )
				for j in chdat : 
					fwfit.write( "{:12f} ".format(effmasserrarr[j]) )
				fwfit.write( "\t" ) 
				outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT , nfit , xraw[0],xraw[-1] ] 
				for it in outarr : 
					fwfit.write( "{}\t".format(it) )
				fwfit.write( "\n" ) 

	def getplotselfwoffdiag(self, ax, dat , basis, nbasis, chdat ,hobj , args , idir=0 , laboverw=None, dashes=(None,None), logplot=False, T=False ) :
		self.getselfw(dat , basis, nbasis, hobj , args )
		self.gloc = self.arrdat
		print len( self.wn ), "points"
		print len( self.gloc ), "components"
		self.glocMat = self.compMat
		print "0-th matrix :"
		for j1 in range(nbasis) :
			for j2 in range(nbasis) :
				print "{:20.8f}\t".format( self.glocMat[j1][j2][0] ),
			print ""
		if T is not False : 
			self.glocMat = transfMjefftot2g( self.glocMat, T ) 
			basis = args.transform
			print "0-th matrix (transformed) :"
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					print "{:20.8f}\t".format( self.glocMat[j1][j2][0] ),
				print ""
		try : 
			if args.sum : 
				xDatSum = self.wn
				yDatSum = np.zeros( len(self.wn) ) 
				labsum  = "Sum["
		except : 
			setattr( args, "sum" , None )  
		self.glocReOffDiag = [ range(nbasis) for j in range(nbasis) ]
		self.glocImOffDiag = [ range(nbasis) for j in range(nbasis) ]
		for k in chdat :
			i , j = [ k[0], k[1] ]
			self.glocReOffDiag[i][j] = self.glocMat[i][j].real
			self.glocImOffDiag[i][j] = self.glocMat[i][j].imag
			re = np.array( self.glocReOffDiag[i][j] ,dtype=float )
			im = np.array( self.glocImOffDiag[i][j] ,dtype=float )
			if args.sinv : 
				re = 1./re
				im = 1./im
		
			ls = [ "-" ] + [ "-." ] + [ "--" ] *3 +  [ "-" ] + [ "-." ] + [ "--" ] *3  
			ms = [ "x" ] + [ "3" ] + [ ">" ] + [ "^" ] + [ "x" ] + [ "+" ] 
			xdat = self.wn
			ydat = im
			xlab = "$\omega$"
			ylab = "Im$\Sigma(\omega)$"
			if args.minus :
				ydat = -ydat
				ylab = "-"+ylab
			if args.minusx :
				xdat = -xdat
			if args.selffactor :
				ydat = ydat*float(args.selffactor)
				ylab = args.selffactor+"*"+ylab
			if args.selfoffset :
				ydat = ydat+float(args.selfoffset)
				ylab = ylab +"+"+ args.selfoffset
			if args.sum : 
				yDatSum = yDatSum + ydat 
			if laboverw :	labelstr = laboverw
			else :		labelstr = bolabel(basis,j)
			if args.solidline :	plotbasis = "jsolid"+basis
			else :			plotbasis = basis
			if args.trans :
				xlab = "Im$\Sigma(\omega)$"
				ylab = ""
				ydatdum = ydat
				xdat = ydatdum
				ydat = xdat
			if args.selfm :
				ydat = np.array(-ydat , dtype=float) 
				ylab = "-"+ylab
			if args.fillb :
				if args.trans : 
					if j==0 :
						ax.fill_betweenx( self.wn, 0, ydat, facecolor='C%d'%j ,alpha=0.3 )
					if j==1 :
						ax.fill_betweenx( self.wn, 0, ydat, facecolor='C%d'%idir ,alpha=0.3 )
				else :
					if j==0 :
						ax.fill_between( self.wn, ydat, facecolor='C%d'%j ,alpha=0.3 )
			if args.fill :
				if j==1 :
					ax.fill( self.wn, ydat, 'C%d'%j ,alpha=0.3 ) 
			if args.sum : 
				labsum = labsum + labelstr + ','
			else : 
				ax.plot( xdat, ydat, ms[i]+ls[j], label="{}{}".format(i,j)   , linewidth=1.5 , dashes=dashes )
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		if args.sum : 
			xdat = xDatSum
			ydat = yDatSum
			if args.trans :
				ydatdum = ydat
				xdat = ydatdum
				ydat = xdat
			ax.plot( xdat, ydat, "ro-" , label=labsum+"]" , linewidth=1.5 , dashes=dashes )
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		#ax.axhline( 0, 0,1 , color='g', dashes=dashpt , lw=0.5 )
		ax.axhline( 0 , color='gray', linestyle="--", lw=0.5 ) 
		if args.selfyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.selfyrange ) 
			else :
				axSetRange( ax , args.selfyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 
		if args.selfxrange : 
			axSetRange( ax , args.selfxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])
	def getplotselfw(self, ax, dat , basis, nbasis, chdat ,hobj , args , idir=0 , laboverw=None, dashes=(None,None) , logplot=False ) :
		print "Reading :", dat
		datfile=None
		if( dat.find(".dat") > -1 ) :
			datfile = dat
		else : 
		 	if basis.find("t") > -1  : 
		 		fhead = '' if basis.find("z")>-1 else 't'
				try : 
					datfile = "{}/lattice/{}selfreal_u{:.2f}_Nplot{}.dat".format(dat,fhead,hobj.UF,hobj.Nplot)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
				except : 
					broadeningtag = ""
					if args.broadeningtag :
						broadeningtag = "_br{:g}".format(float(args.broadening))
					if args.epsilontagfloat :
						datfile = "{}/lattice/u{:.3f}/{}selfreal_Nplot{}_ep{:.3f}{}_{}th.dat".format(dat,hobj.UF,fhead,hobj.Nplot,hobj.epsilon,broadeningtag,hobj.count)
					else :
						datfile = "{}/lattice/u{:.3f}/{}selfreal_Nplot{}_ep{:.3g}{}_{}th.dat".format(dat,hobj.UF,fhead,hobj.Nplot,hobj.epsilon,broadeningtag,hobj.count)
		 	elif basis.find("j") > -1  or basis.find("e") > -1  : 
				try : 
					datfile = "{}/lattice/selfreal_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
					print "Reading : ", datfile,
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
					print "successed."
				except :
					print "failed."
					broadeningtag = ""
					if args.broadeningtag :
						broadeningtag = "_br{:g}".format(float(args.broadening))
					if args.epsilontagfloat :
						datfile = "{}/lattice/u{:.3f}/selfreal_Nplot{}_ep{:.3f}{}_{}th.dat".format(dat,hobj.UF,hobj.Nplot,hobj.epsilon,broadeningtag,hobj.count)
					else :
						datfile = "{}/lattice/u{:.3f}/selfreal_Nplot{}_ep{:.3g}{}_{}th.dat".format(dat,hobj.UF,hobj.Nplot,hobj.epsilon,broadeningtag,hobj.count)
					print "Reading : ", datfile,
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
					print "successed."
		print "Reading : ", datfile
		rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		self.wn   = rawdatarr.transpose()[0]
		print len( self.wn ), "points (lines/rows)"
		self.gloc = rawdatarr.transpose()[1:]
		print len( self.gloc ), "components (columns)"
		glocRe		= self.gloc[::2]
		glocIm		= self.gloc[1::2]
		glocComp	= glocRe + glocIm * 1j
		dim		= int(np.sqrt( len(glocComp)/hobj.Ni ))
		if hobj.Ni*dim*dim != len(glocComp) :
			print "ERROR: hobj.Ni*dim*dim != len(glocComp) ", hobj.Ni*dim*dim , np.shape(glocComp)
			sys.exit(1)
		glocComp	= glocRe + glocIm * 1j
		glocComp	= glocComp.reshape( hobj.Ni, dim,dim,len(self.wn) ) 
		#print "TNC : ", hobj.tNC
		#sys.exit(1)

		#self.glocReDiag = range(nbasis)
		#self.glocImDiag = range(nbasis)
		a = range( len(self.wn) ) 
		try : 
			if args.sum : 
				xDatSum = self.wn
				yDatSum = np.zeros( len(self.wn) ) 
				labsum  = "Sum["
		except : 
			setattr( args, "sum" , None )  
		for j in chdat :
			#self.glocReDiag[j] = self.gloc[2*j   +2*nbasis*j]
			#self.glocImDiag[j] = self.gloc[2*j+1 +2*nbasis*j]
			#re = np.array( self.glocReDiag[j] ,dtype=float )
			#im = np.array( self.glocImDiag[j] ,dtype=float )
			tnc	= hobj.tNC
			re	= glocComp[int(j/tnc),j%tnc,j%tnc].real
			im	= glocComp[int(j/tnc),j%tnc,j%tnc].imag
			if args.sinv : 
				re = 1./re
				im = 1./im
		
			lc = [ "-" ]*2 + [ "-" ] *2 + [ ":" ] *2
			xdat = self.wn
			ydat = im
			xlab = "$\omega$"
			ylab = "Im$\Sigma(\omega)$"
			if args.realpart :
				ydat = re
				ylab = "Re$\Sigma(\omega)$"
			if args.sum : 
				yDatSum = yDatSum + ydat 
			if laboverw :	labelstr = laboverw
			else :		labelstr = bolabelArr(basis)
			if args.solidline :	plotbasis = "jsolid"+basis
			else :			plotbasis = basis
			if args.minusx :
				xdat = -xdat
			if args.trans :
				ylabdum	= "{}".format(ylab)
				ylab = xlab
				ylab = ylabdum
				ydatdum = ydat
				xdat = ydatdum
				ydat = xdat
			if args.selfm :
				ydat = np.array(-ydat , dtype=float) 
				ylab = "-"+ylab
			if args.fillb :
				if args.trans : 
					if j==0 :
						ax.fill_betweenx( self.wn, 0, ydat, facecolor='C%d'%j ,alpha=0.3 )
					if j==1 :
						ax.fill_betweenx( self.wn, 0, ydat, facecolor='C%d'%idir ,alpha=0.3 )
				else :
					if j==0 :
						ax.fill_between( self.wn, ydat, facecolor='C%d'%j ,alpha=0.3 )
			if args.fill :
				if j==1 :
					#ax.fill( self.wn, ydat, 'C%d'%j ,alpha=0.3 ) 
					ax.fill_between( self.wn, ydat, facecolor='C%d'%j ,alpha=0.3 ) 
			if args.sum : 
				labsum = labsum + labelstr[j] + ','
			else : 
				if logplot :
					if args.one : 
						ax.loglog( xdat, ydat, blt(plotbasis,j) , label=labelstr[j] , linewidth=1.5 , dashes=dashes )
					else : 
						ax.loglog( xdat, ydat, blt(plotbasis,j) , label=labelstr[j] , linewidth=1.5 )
				else : 
					ax.plot( xdat, ydat, blt(plotbasis,j) , label=labelstr[j] , linewidth=1.5 , dashes=dashes )
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		if args.sum : 
			xdat = xDatSum
			ydat = yDatSum
			if args.trans :
				ydatdum = ydat
				xdat = ydatdum
				ydat = xdat
			if logplot :
				ax.loglog( xdat, ydat, "ro-" , label=labsum+"]" , linewidth=1.5 , dashes=dashes )
			else :
				ax.plot( xdat, ydat, "ro-" , label=labsum+"]" , linewidth=1.5 , dashes=dashes )
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		#ax.axhline( 0, 0,1 , color='g', dashes=dashpt , lw=0.5 )
		ax.axhline( 0 , color='gray', linestyle="--", lw=0.5 ) 
		if args.selfyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.selfyrange ) 
			else :
				axSetRange( ax , args.selfyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 
		if args.selfxrange : 
			axSetRange( ax , args.selfxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])
	def getselfiw(self, dat , nbasis, hobj , args ) :
		hobj.readEclustersoc()		#saved in hobj.ecluster[i][j]
		if( dat.find(".dat") > -1 ) :
			datfile = dat
		else : 
			if args.count :	count = args.count
			else :		count = hobj.count
			if args.finegrid : 
				ep=0.03
				nplot=1024
				if args.nplot : nplot = int(args.nplot)
				datfile = "{}/lattice/u{:.3f}/selfimag_Nplot{}_ep{:.3f}_{}th.dat".format(dat, hobj.UF, nplot, ep , count )
				if hobj.ni>1 : 
					print "Need to modify this code to handle multi-impurity data."
					sys.exit(1)
			elif args.finegrid2 or args.finegrid3 or args.finegrid4 : 
				nplot=2048
				if args.finegrid3 : 
					nplot=2049
				if args.finegrid4 : 
					nplot=2050
				datfile = "{}/ongoing/selfiw_Nplot{}_{}th.dat".format(dat, nplot, count )
				if hobj.ni>1 : 
					print "Need to modify this code to handle multi-impurity data."
					sys.exit(1)
			else : 
				iimp = hobj.ni-1
				if args.indimp : 
					iimp = int(args.indimp)
				datfile = "{}/ongoing/self{}_u{:.2f}_{}th.txt".format(dat, iimp, hobj.UF, count )
		print "Reading :", datfile
		rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		self.wn   = rawdatarr.transpose()[0]
		self.arrdat = rawdatarr.transpose()[1:]
		self.compMat = complexArrDatReturn( self.arrdat, nbasis ) 
	def getplotselfiw(self, ax, dat , basis, nbasis, chdat ,hobj , args, T=False ) :
		self.getplotselfiw_lc(ax, dat , basis, nbasis, chdat ,hobj , args , "", T=T ) 
	def getplotselfiw_lc(self, ax, dat , basis, nbasis, chdat ,hobj , args , lc , laboverw=False , dashes=(None,None) , logplot=False, T=False , lstylearr=False, labsuffix=False , eig=False, mec=None, mfc=None ) :
		print "Calling 'getplotselfiw_lc()' .."
		self.getselfiw( dat, nbasis, hobj, args );
		self.gloc = self.arrdat
		nplot = len( self.wn )
		ncomp = len( self.gloc )
		print nplot, " points(lines)"
		print ncomp, " components(of matrix)"
		self.glocReDiag = range(nbasis)
		self.glocImDiag = range(nbasis)
		self.glocMat = self.compMat
		if args.addSOCt2g :
			print "ONESOC : ", float(args.addSOCt2g)
			self.glocMat = addSOCMatt2g( self.glocMat, S=float(args.addSOCt2g) ) 
		if args.addSOCjeff :
			print "ONESOC : ", float(args.addSOCjeff)
			self.glocMat = addSOCMat( self.glocMat, basis="jeff", S=float(args.addSOCjeff) ) 
		if args.zerosoc :
			basis = "t"
		print "0-th matrix :"
		for j1 in range(nbasis) :
			for j2 in range(nbasis) :
				print "{:20.8f}\t".format( self.glocMat[j1][j2][0] ),
			print ""
		if T is not False : 
			self.glocMat = transfMjefftot2g( self.glocMat, T ) 
			basis = args.transform
			print "0-th matrix (transformed) :"
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					print "{:20.8f}\t".format( self.glocMat[j1][j2][0] ),
				print ""
		#print "SHape : ", np.shape( self.glocMat ) 
		#print "SHape : ", np.shape( np.reshape( self.glocMat,  (36, -1) ) ) 

		if args.eigenval : 
			self.glocMat = eigMat( self.glocMat ) 
			print "Diagonalized."
			print "0-th matrix (diagonalized) :"
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					print "{:20.8f}\t".format( self.glocMat[j1][j2][0] ),
				print ""

		Zarr = np.zeros( nbasis ) 
		Zavgtargetarr = np.zeros( nbasis ) 
		Gammatargetarr = np.zeros( nbasis ) 
		effmasserrarr = np.zeros( nbasis ) 
		quantclass = dumclass()
		setattr( quantclass , 'power',		np.zeros( nbasis )  ) 
		setattr( quantclass , 'effmass',		np.zeros( nbasis )  ) 
		setattr( quantclass , 'powerfitconst',	np.zeros( nbasis )  ) 
		setattr( quantclass , 'polyfitconst',	np.zeros( nbasis )  ) 
		setattr( quantclass , 'scattrate',	np.zeros( nbasis )  ) 
		setattr( quantclass , 'erreffmass',	np.zeros( nbasis )  ) 
		for j in chdat :
			print "shapes : ", np.shape( self.glocReDiag ), np.shape( self.glocMat )
			self.glocReDiag[j] = self.glocMat[j][j].real
			self.glocImDiag[j] = self.glocMat[j][j].imag
			re = np.array( self.glocReDiag[j] )
			im = np.array( self.glocImDiag[j] )
			xdat = self.wn
			ydat = im
			if args.cutoff : 
				cutoff = float( args.cutoff ) 
				print "Cutoff : ", cutoff 
				ydat = ydat[ xdat > cutoff ] 
				xdat = xdat[ xdat > cutoff ] 
			if args.realpart : ydat = re
		
			if laboverw :	labelstr = laboverw
			else :		labelstr = bolabelArr(basis)
			if labsuffix :	labelstr = [labsuffix+labelstr[j2] for j2 in range(nbasis) ]
			lt = [ "-" ]   + [ "-" ] *2 + [ "-" ]   +  [ ":" ] *2
			if basis.find("t")>-1 : 
				lt = [ "-" ]*2 + [ "-" ] *2 + [ ":" ] *2
			elif args.transform and args.transform.find("t")>-1 : 
				lt = [ "-" ]*2 + [ "-" ] *2 + [ ":" ] *2
			lstyle = markers[j]+lt[j] #lstyle = lc+markers[j]+linestyles[j]
			print "MFC : ", mfc
			print "lstyle, mfc, mec : ", lstyle, mfc, mec
			if args.emptymarker :
				mfc = 'w'
			if args.line : lstyle = lt[j]
			if lstylearr : lstyle=lstylearr[j]
			ylab = "Im$\Sigma(i\omega_n)$"
			xlab = "$\omega_n$"
			if args.minus :
				ydat = -ydat
				ylab = "-"+ylab
			if args.minusx :
				xdat = -xdat
			if args.selffactor :
				ydat = ydat*float(args.selffactor)
				ylab = args.selffactor+"*"+ylab
			if args.selfoffset :
				print "OFFSET : ", float(args.selfoffset)
				ydat = ydat+float(args.selfoffset)
				ylab = ylab +"+"+ args.selfoffset
			if args.removelinear : 
				from scipy.optimize import curve_fit
				def funclin( x, a, b ) : return a*x + b
				typefit = "lin"
				formfunc= 'ax+b'
				nvarfunc= 2
				nfit = 3 
				if args.ndatalinearfit :
					nfit = int(args.ndatalinearfit)
				iinit = 0 
				if args.fitoffset : iinit += int(args.fitoffset)
				xraw = xdat[iinit:iinit+nfit] 
       				yraw = ydat[iinit:iinit+nfit] 
				if args.linearfitfromzero : 
					xraw = np.insert( xraw, 0, 0 )
					yraw = np.insert( yraw, 0, 0 )
				popt, pcov = curve_fit(funclin, xraw, yraw, bounds=([-99999,-1.], [999000., 1]))
				print "removing-FIT({};{}) : a={} b={}".format( formfunc, typefit, *popt )
				a = popt[0]
				devidat = xdat*a - ydat
				devidat = np.array( [xdat, devidat, ydat, xdat*a, devidat/xdat/a , xdat*0.+a ] ).transpose()
				np.savetxt(  "tmprmlinear_J{}_{}{}".format(hobj.J,basis,chdat[0]) , devidat ) 
			if args.difftwopoints : 
				xdum1= xdat[:-1]
				xdum2= xdat[1:]
				ydum1= ydat[:-1]
				ydum2= ydat[1:]
				xdat = xdat[:-1] 
				ydat = (ydum1-ydum2)/(xdum1-xdum2)
			if args.trans :
				if logplot :
					ax.loglog( ydat, xdat, lstyle , label=labelstr[j] , dashes=dashes, color=lc )
				else :
					ax.plot(   ydat, xdat, lstyle , label=labelstr[j] , dashes=dashes, color=lc )
				if j==0 and args.fillb :
					ax.fill( ydat, xdat, lstyle ,alpha=0.2 )
				ylab = ""
			else :
				if args.dashedline :
					if logplot :
						ax.loglog( xdat, ydat, lstyle , label=labelstr[j] , dashes=dashes, color=lc , mfc=mfc, mec=mec )
					else :
						ax.plot(   xdat, ydat, lstyle , label=labelstr[j] , dashes=dashes, color=lc , mfc=mfc, mec=mec )
				else : 
					if logplot :
						ax.loglog( xdat, ydat, lstyle , label=labelstr[j], color=lc , mfc=mfc, mec=mec )
					else :
						ax.plot(   xdat, ydat, lstyle , label=labelstr[j], color=lc , mfc=mfc, mec=mec )
			if len(chdat)<2 : ylab = ylab + "$^{{{}}}$".format(blabel(basis,j))
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
			if args.write : 
				dum = np.array( [xdat,ydat] ) 
				np.savetxt( "tmpplotdat", dum.transpose() ) 
			if args.drawabc : 
				def funcabc( x, a, b ,c ) : return a*x*x + b*x + c 
				nfit = 3 
				nresult = 100
				if args.ndatalinearfit :
					nfit = int(args.ndatalinearfit)
					if nfit > nresult  : 
						nresult = nfit*2
				iinit = 0 
				if args.fitoffset : iinit += int(args.fitoffset)
				xraw = xdat[iinit:iinit+nfit] 
       				yraw = ydat[iinit:iinit+nfit] 
				if logplot :
					xfit = np.logspace( np.log10(xraw[0]*0.1),  np.log10(xdat[-10]) , nresult ) 
				else : 
					xfit = np.linspace( xraw[0]*0.1,  xraw[-1]*2. , nresult ) 
				cba = args.drawabc.split("_")[-1:-4:-1]
				c,b,a = np.array( cba, dtype=float ) 
				carr = [c]
				print "ABC : ", a,b,c
				if args.carray : carr = np.array( args.carray.split("_"), dtype=float ) 
				for c in carr : 
					print "ABC : ", a,b,c
					yfit = funcabc(xfit, a,b,c ) 
					fitls = "b:"
					fitlab = None
					if args.carray : 
						fitls = "--"
						fitlab = '$c$={:.3g}'.format(c)
					if logplot :
						ax.loglog( xfit, yfit, fitls , lw=1 , label = fitlab )
					else :
						ax.plot( xfit, yfit, fitls , lw=1 , label = fitlab )
				if args.carray : 
				 	ax.legend()
			if args.linearfit or args.fitting or args.quadraticfit or args.xlnxfit : 
				from scipy.optimize import curve_fit
				def funclin( x, a, b ) : return a*x 
				fitfunc = funclin
				typefit = "lin"
				formfunc= 'ax+b'
				nvarfunc= 2
				if args.fitting : 
					typefit = args.fitting
				elif args.quadraticfit : 
					typefit = 'quad'
				elif args.xlnxfit : 
					typefit = 'xlnx'
				nfit = 5 
				nresult = 30
				print "Entering linearfit() ..."
				if args.ndatalinearfit :
					nfit = int(args.ndatalinearfit)
					if nfit > nresult  : 
						nresult = nfit*2
				nfitresult = 5*nresult
				#if hobj.J<0.1 : 
				#	nfit = 20 
				print "nfit : ", nfit
				iinit = 0 
				if args.fitoffset : iinit += int(args.fitoffset)
				xraw = xdat[iinit:iinit+nfit] 
       				yraw = ydat[iinit:iinit+nfit] 
				if args.linearfitfromzero : 
					xraw = np.insert( xraw, 0, 0 )
					yraw = np.insert( yraw, 0, 0 )
				if nfit > 4 : 
					print "FIT_xraw : " , xraw[ [0,1,-2,-1] ]
					print "FIT_yraw : " , yraw[ [0,1,-2,-1] ]
				else : 
					print "FIT_xraw : " , xraw
					print "FIT_yraw : " , yraw

				if typefit.find("poly4")>-1 :
					formfunc= 'ax^2 + bx + c + d*x^3 + e*x^4'
					nvarfunc= 5
					def funcpoly4( x, a, b, c ,d, e) : return a*x*x + b*x + c + d*x*x*x + e*x*x*x*x
					fitfunc = funcpoly4
					if args.fitnoconst :
						nvarfunc= 4
						def funcpoly4( x, a, b, d, e) : return a*x*x + b*x + d*x*x*x + e*x*x*x*x
						fitfunc = funcpoly4
				if typefit.find("poly3")>-1 :
					formfunc= 'ax^2+bx+c+d*x^3'
					nvarfunc= 4
					def funcpoly3( x, a, b, c ,d ) : return a*x*x + b*x + c + d*x*x*x
					fitfunc = funcpoly3
					if args.fitnoconst :
						nvarfunc= 3
						def funcpoly3( x, a, b, d ) : return a*x*x + b*x + d*x*x*x
						fitfunc = funcpoly3
				if typefit.find("poly2")>-1 :
					formfunc= 'ax^2+bx+c'
					nvarfunc= 3
					def funcpoly2( x, a, b, c ) : return a*x*x + b*x + c
					fitfunc = funcpoly2
				elif typefit.find("fixlin")>-1 :
					formfunc= 'ax^2+b'
					nvarfunc= 2
					def funcfixlin( x, a, b ) : return a*x*x + b
					fitfunc = funcfixlin
				elif typefit.find("fixtemp")>-1 :
					formfunc= 'a(x^2-pi^2*T^2)'
					nvarfunc= 1
				elif typefit =="const" :
					formfunc= 'a'
					nvarfunc= 1
					def funcconst( x, a ) : return x*0.+a
					fitfunc = funcconst
				elif typefit.find("xlnx")>-1 :
					formfunc= 'a x ln(x/c)'
					nvarfunc= 1
					const = 1.
					if args.xlnxfitconst : const = float(args.xlnxfitconst)
					def funcconst( x, a ) : return a * x * np.log(x/const)
					fitfunc = funcconst
				elif typefit =="power" : 
					formfunc= 'a x^b'
					nvarfunc= 2
					def funcpower( x, a , b ) : return a * np.power(x,b)
					fitfunc = funcpower
				elif typefit =="powerplusconst" : 
					formfunc= 'a x^b + c '
					nvarfunc= 3
					def funcpower( x, a , b, c ) : return a * np.power(x,b) + c
					fitfunc = funcpower

				if typefit.find("fixlin")>-1 :
					def funclin( x, a, b ) : return a*x
					fitfunc = funclin
					nfitlin = 10
					xraw = xdat[0:nfitlin] 
       					yraw = ydat[0:nfitlin] 
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=(0., [999000., 1]))
					print "FIT({};{}) : ".format( formfunc, 'linear'+typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					xfit = np.linspace( xraw[0]*0.1,  xraw[-1]*2. , nfitresult ) 
					yfit = fitfunc(xfit,*popt) 
					fitls = "k:"
					if args.setlinear : 
						b = float( args.setlinear ) 
					else :
						if logplot :
							ax.loglog( xfit, yfit, fitls , lw=1 )
						else :
							ax.plot( xfit, yfit, fitls , lw=1 )
						b = popt[0]

					xraw = xdat[iinit:iinit+nfit] 
       					yraw = ydat[iinit:iinit+nfit] 
					def funcfixlin( x, a, c ) : return a*x*x + b*x + c
					fitfunc = funcfixlin
					initvarlowerarr = [-99999, 0 ] 
					initvarupperarr = [99999 , 1 ]
					if nfit >1  :
						popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=( initvarlowerarr , initvarupperarr ) )
					else :
						popt = [ xraw / yraw ] 
					if args.fitoffset : 
						if logplot : 
							xfit = np.logspace( np.log10(xdat[0]),  np.log10(xraw[-1]) , num=nfitresult ) 
						else : 
							xfit = np.linspace( xdat[0],  xraw[-1] , nfitresult ) 
					else : 
						if logplot : 
							xfit = np.logspace( np.log10(xraw[0]),  np.log10(xraw[-1]) , num=nfitresult ) 
							xfit3= np.logspace( np.log10(xraw[0]/1e3),  np.log10(xraw[0]) , num=nfitresult ) 
						else : 
							xfit = np.linspace( xraw[0],  xraw[-1] , nfitresult ) 
					yfit = fitfunc(xfit,*popt) 
					yfit3= fitfunc(xfit3,*popt) 
					print "XFIT : ", xfit[ [0,1,-2-1] ]
					print "YFIT : ", yfit[ [0,1,-2-1] ]
					print "FIT({}+{}x;{}) : ".format( formfunc, b, 'quad'+typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					if args.minus :
						for jj in range(nvarfunc) :
							popt[jj] = -popt[jj]
					def returnT( a,c ) :
						T = np.sqrt( - c / a / np.pi/np.pi ) 
						return T, T*11605., 1./T
					print "\tT,T[K],beta : ", "{}".format( returnT(*popt) )
					fitls = "k--"
					fitls3= "r--"
					if logplot :
						ax.loglog( xfit, yfit, fitls , lw=1 )
						ax.loglog( xfit3, yfit3, fitls3 , lw=1 )
					else :
						ax.plot( xfit, yfit, fitls , lw=1 )
					with open( "tmpfitfix{}".format(j), 'a' ) as fwfit :
						outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT ] 
						for it in outarr : 
							fwfit.write( "{}\t".format(it) )
						fwfit.write( "{} {} {} {} {}\n".format( nfit, popt[0], b , popt[1], returnT(*popt) ) )
				elif typefit.find("fixtemp")>-1 :
					def funclin( x, a ) : return a*x
					fitfunc = funclin
					nfitlin = 10
					xraw = xdat[0:nfitlin] 
       					yraw = ydat[0:nfitlin] 
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=([0.], [999000.]) )
					print "FIT({};{}) : ".format( formfunc, 'linear'+typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					if hobj.IT :
						wtemp = np.pi/hobj.IT
					else :
						wtemp = np.pi/hobj.beta
					tempcutind = 0
					for ii in range(len(xdat)) : 
						if wtemp < xdat[ii] :
							tempcutind = ii
							break
					print "Ind_cuttemp : ", tempcutind , wtemp, xdat[tempcutind], xdat[tempcutind-1]
					xfitmax = xraw[-1]*2.
					if args.fitoffset : 
						xfitmax = xdat[tempcutind-1]
					xfit = np.linspace( xraw[0]*0.1, xfitmax , nfitresult ) 
					yfit = fitfunc(xfit,*popt) 
					fitls = "k:"
					if args.setlinear : 
						b = float( args.setlinear ) 
					else :
						if logplot :
							ax.loglog( xfit, yfit, fitls , lw=1 )
						else :
							ax.plot( xfit, yfit, fitls , lw=1 )
						b = popt[0]

					xraw = xdat[iinit:iinit+nfit] 
       					yraw = ydat[iinit:iinit+nfit] 
					if hobj.IT :
						temp = 1./hobj.IT
					else :
						temp = 0.
					piTsq = np.pi*np.pi*temp*temp
					def funcfixtemp( x, a ) : return a*( x*x - piTsq) + b*x 
					fitfunc = funcfixtemp
					initvarlowerarr = [-99999]
					initvarupperarr = [99999 ]
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=( initvarlowerarr , initvarupperarr ) )
					a = popt[0]
					if args.fitoffset : 
						if logplot : 
							xfit2= np.logspace( np.log10(xdat[0]),  np.log10(xraw[0]) , num=nfitresult ) 
							xfit = np.logspace( np.log10(xraw[0]),  np.log10(xraw[-1]) , num=nfitresult ) 
							xfit3= np.logspace( np.log10(xraw[0]/1e3),  np.log10(xraw[0]) , num=nfitresult ) 
						else : 
							xfit2= np.linspace( xdat[0],  xraw[0] , nfitresult ) 
							xfit = np.linspace( xraw[0],  xraw[-1] , nfitresult ) 
					else : 
						if logplot : 
							xfit = np.logspace( np.log10(xraw[0]),  np.log10(xraw[-1]) , num=nfitresult ) 
							xfit3= np.logspace( np.log10(xraw[0]/1e3),  np.log10(xraw[0]) , num=nfitresult ) 
						else : 
							xfit = np.linspace( xraw[0],  xraw[-1] , nfitresult ) 
							xfit3= np.linspace( xraw[0]/1e3,  xraw[0] , nfitresult ) 
					if args.fitoffset : 
						yfit2= fitfunc(xfit2,popt) 
						print "XFIT2: ", xfit2[ [0,1,-2-1] ]
						print "YFIT2: ", yfit2[ [0,1,-2-1] ]
					yfit = fitfunc(xfit,popt) 
					yfit3= fitfunc(xfit3,popt) 
					print "XFIT : ", xfit[ [0,1,-2-1] ]
					print "YFIT : ", yfit[ [0,1,-2-1] ]
					print "FIT({}+{}x;{}) : ".format( formfunc, b, 'quad'+typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					print "(pi*T)^2, -a*piTsq, -a*pi*pi*T : ", piTsq, -a*piTsq , -a*np.pi*np.pi*temp
					fitls = "k--"
					fitls2= "r--"
					fitls3= "b--"
					if logplot :
						if args.fitoffset : 
							ax.loglog( xfit,  yfit,  fitls  , lw=1 )
							ax.loglog( xfit2, yfit2, fitls2 , lw=1 )
						else : 	
							ax.loglog( xfit, yfit, fitls , lw=1 )
							ax.loglog( xfit3, yfit3, fitls3, lw=1 )
					else :
						if args.fitoffset : 
							ax.plot( xfit,  yfit,  fitls , lw=1 )
							ax.plot( xfit2, yfit2, fitls2, lw=1 )
						else : 
							ax.plot( xfit, yfit, fitls , lw=1 )
							ax.plot( xfit3, yfit3, fitls3, lw=1 )
					with open( "tmpfitfixtemp{}".format(j), 'a' ) as fwfit :
						outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT ] 
						for it in outarr : 
							fwfit.write( "{}\t".format(it) )
						fwfit.write( "{} {} {}\n".format( nfit, popt[0], b  ) )
				elif typefit.find("poly2")>-1  or (typefit.find("poly3")>-1) or (typefit.find("poly4")>-1) :
					if nvarfunc == 3 :
						initvararr = [-20]*(nvarfunc-1) + [-0.5]
						initvararr2= [ 20]*(nvarfunc-1) + [ 0.5]
						if args.fitnoconst : 
							initvararr = [-100, -00 , -100 ]
							initvararr2= [ 000,  30 ,  000 ]
					elif nvarfunc == 4 :
						#initvararr = [-20]*(nvarfunc-2) + [-0.5] +[-20]
						#initvararr2= [ 20]*(nvarfunc-2) + [ 0.5] +[ 20]
						initvararr = [-100, -00 , -0.1, -100 ]
						initvararr2= [ 000,  30 ,  0.1,  000 ]
						if args.fitnoconst : 
							#initvararr = [-100, -00 , -100, -100 ]
							#initvararr2= [ 000,  30 ,  000,  100 ]
							initvararr = [-100, -00 , -100, -100 ]
							initvararr2= [ 000,  30 ,  000,  1000 ]
					elif nvarfunc == 5 :
						#initvararr = [-100, -00 , -0.1, -100, -100 ]
						#initvararr2= [ 000,  30 ,  0.1,  000,  100 ]
						#initvararr = [-1000, -000 , -10, -1000, -100 ]
						#initvararr2= [ 0000,  300 ,  10,  0000,  1000000 ]
						initvararr = [-1000, -00 , -0, -1000, -100 ]
						initvararr2= [ 0000,  30 ,  1,  0000,  1000000 ]
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=( initvararr , initvararr2) )
					perr	= np.sqrt(np.diag(pcov))
					if args.fitoffset : 
						if logplot : 
							xfit = np.logspace( np.log10(xdat[0]),  np.log10(xraw[-1]) , num=nfitresult ) 
							xfit2= np.logspace( np.log10(xraw[-1]),  np.log10(xdat[2*nresult]) , num=nfitresult ) 
						else : 
							xfit = np.linspace( xdat[0],  xraw[-1] , nfitresult ) 
							xfit2= np.linspace( 0	   ,  xdat[2*nresult] , nfitresult ) 
					else : 
						if logplot : 
							xfit = np.logspace( np.log10(xraw[0]),  np.log10(xraw[-1]) , num=nfitresult ) 
							xfit2= np.logspace( np.log10(xraw[-1]),  np.log10(xdat[2*nresult]) , num=nfitresult ) 
						else : 
							xfit = np.linspace( xraw[0],  xraw[-1] , nfitresult ) 
							xfit2= np.linspace( 0	   ,  xdat[2*nresult] , nfitresult ) 
					yfit = fitfunc(xfit,*popt) 
					yfit2 = fitfunc(xfit2,*popt) 
					print "XFIT : ", xfit[ [0,1,-2-1] ]
					print "YFIT : ", yfit[ [0,1,-2-1] ]
					print "FIT({};{}) : ".format( formfunc, typefit )
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					print "FITerr : "
					for jj in range(nvarfunc) :
						print "{}".format(perr[jj]),
					print ""
					if args.minus :
						for jj in range(nvarfunc) :
							popt[jj] = -popt[jj]
					Zftn = lambda b : 1./(1.-b)
					Zarr[j] = Zftn( popt[1] ) 
					(quantclass.effmass)[j]		= 1./Zftn( popt[1] )
					(quantclass.erreffmass)[j]	= perr[1]
					if typefit=='poly2' or typefit=='poly3' or typefit=='poly4' : 
						(quantclass.polyfitconst)[j]	= -popt[2]
						if args. fitnoconst : 
							(quantclass.polyfitconst)[j]	= 0.
					else : 
						(quantclass.polyfitconst)[j]	= -999.
					if nvarfunc == 3 :
						def mZCT( a,b,c ) :
							m = 1.-b
							Z = 1./m
							C = a /m
							T = np.sqrt( - c / a / np.pi/np.pi ) 
							return m, Z, C, T, T*11605., 1./T
						print "\tm,Z,C,T,T[K],beta : ", "{}".format( mZCT(*popt) )
					elif nvarfunc == 4 :
						def mZCT( a,b,c,d ) :
							m = 1.-b
							Z = 1./m
							C = a /m
							T = np.sqrt( - c / a / np.pi/np.pi ) 
							return m, Z, C, T, T*11605., 1./T
						print "\tm,Z,C,T,T[K],beta : ", "{}".format( mZCT(*popt) )
					elif nvarfunc == 5 :
						def mZCT( a,b,c,d,e ) :
							m = 1.-b
							Z = 1./m
							C = a /m
							T = np.sqrt( - c / a / np.pi/np.pi ) 
							return m, Z, C, T, T*11605., 1./T
						print "\tm,Z,C,T,T[K],beta : ", "{}".format( mZCT(*popt) )
					fitls = "k--"
					fitls2= "r--"
					if logplot :
						ax.loglog( xfit , yfit , fitls  , lw=1 )
						ax.loglog( xfit2, yfit2, fitls2 , lw=1 )
					else :
						ax.plot( xfit , yfit , fitls  , lw=1 )
						ax.plot( xfit2, yfit2, fitls2 , lw=1 )
				elif typefit == "const" : 
					initvararr = [99999]*nvarfunc
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=(-99999 , initvararr) )
					if args.fitoffset : 
						if logplot : 
							xfit = np.logspace( np.log10(xdat[0]),  np.log10(xraw[-1]) , num=nfitresult ) 
							xfit2= np.logspace( np.log10(xraw[-1]),  np.log10(xdat[-1]) , num=nfitresult ) 
						else : 
							xfit = np.linspace( xdat[0],  xraw[-1] , nfitresult ) 
					else : 
						if logplot : 
							xfit = np.logspace( np.log10(xraw[0]),  np.log10(xraw[-1]) , num=nfitresult ) 
							xfit2= np.logspace( np.log10(xraw[-1]),  np.log10(xdat[-1]) , num=nfitresult ) 
						else : 
							xfit = np.linspace( xraw[0],  xraw[-1] , nfitresult ) 
					yfit = fitfunc(xfit,*popt) 
					yfit2= fitfunc(xfit2,*popt) 
					print "XFIT : ", xfit[ [0,1,-2-1] ]
					print "YFIT : ", yfit[ [0,1,-2-1] ]
					print "FIT({};{}) : ".format( formfunc, typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					if args.minus :
						for jj in range(nvarfunc) :
							popt[jj] = -popt[jj]
					fitls = "k--"
					if logplot :
						ax.loglog( xfit, yfit, fitls , lw=1 )
					else :
						ax.plot( xfit, yfit, fitls , lw=1 )
					fitls2= "r--"
					if logplot :
						ax.loglog( xfit,  yfit,  fitls  , lw=1 )
						ax.loglog( xfit2, yfit2, fitls2 , lw=1 )
					else :
						ax.plot( xfit,  yfit,  fitls , lw=1 )
						ax.plot( xfit2, yfit2, fitls2, lw=1 )
				elif typefit.find("xlnx")>-1 :
					initvararr = [99999]*nvarfunc
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=(-99999 , initvararr) )
					if args.fitoffset : 
						if logplot : 
							xfit = np.logspace( np.log10(xdat[0]),  np.log10(xraw[-1]) , num=nfitresult ) 
							xfit2= np.logspace( np.log10(xraw[-1]),  np.log10(xdat[-1]) , num=nfitresult ) 
						else : 
							xfit = np.linspace( xdat[0],  xraw[-1] , nfitresult ) 
					else : 
						if logplot : 
							xfit = np.logspace( np.log10(xraw[0]),  np.log10(xraw[-1]) , num=nfitresult ) 
							xfit2= np.logspace( np.log10(xraw[-1]),  np.log10(xdat[-1]) , num=nfitresult ) 
						else : 
							xfit = np.linspace( xraw[0],  xraw[-1] , nfitresult ) 
					yfit = fitfunc(xfit,*popt) 
					yfit2= fitfunc(xfit2,*popt) 
					print "XFIT : ", xfit[ [0,1,-2-1] ]
					print "YFIT : ", yfit[ [0,1,-2-1] ]
					print "FIT({};{}) : ".format( formfunc, typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					if args.minus :
						for jj in range(nvarfunc) :
							popt[jj] = -popt[jj]
					fitls = "k--"
					if logplot :
						ax.loglog( xfit, yfit, fitls , lw=1 )
					else :
						ax.plot( xfit, yfit, fitls , lw=1 )
					fitls2= "r--"
					if logplot :
						ax.loglog( xfit,  yfit,  fitls  , lw=1 )
						ax.loglog( xfit2, yfit2, fitls2 , lw=1 )
					else :
						ax.plot( xfit,  yfit,  fitls , lw=1 )
						ax.plot( xfit2, yfit2, fitls2, lw=1 )
				elif typefit.find("lin")>-1 :
					initvararr = [-0000, -00 ]
					initvararr2= [ 1000,  30 ]
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=(initvararr , initvararr2) )
					xfit = np.linspace( xdat[:nfit][0]*0.1,  5 , nfitresult ) 
					yfit = fitfunc(xfit,*popt) 
					print "FIT({};{}) : ".format( formfunc, typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					Zftn = lambda a,b : 1./(1.+a)
					print "Z : {} {} {}".format( j, hobj.J, Zftn(*popt) ) 
					print "m* = 1/Z : {} {} {} {}".format( j, hobj.J, 1./Zftn(*popt), hobj.D ) , "\n"
					fitls = "k-."
					Zarr[j] = Zftn(*popt)
					effmasserrarr[j] = np.sqrt( np.diag(pcov)[0] )
					if logplot :
						#if hobj.IT<1 : 
						ax.loglog( xfit, yfit, fitls , lw=1 )
					else :
						ax.plot( xfit, yfit, fitls , lw=1 )

					with open( "tmplinearfitavg{}".format(j), 'w' ) as fwfitavg :
						for ii in range(30) : 
							ydum = ydat[ii*5]
							xdum = xdat[ii*5]
							steepness = (0.-ydum) / (0.-xdum)
							Zavg = 1./(1.+steepness)
							Gamma= Zavg * ydum	#Inverse quasiparticle lifetime (Mravlje201!)
							outarr = [ xdum, ydum, steepness, 1./Zavg, Zavg, xdum*11605. , Gamma ] 
							for it in outarr : 
								fwfitavg.write( "{:8f}\t".format(it) )
							fwfitavg.write( "\n" )

					ydum = ydat[39]
					xdum = xdat[39]
					tempdum = xdum*11605. 
					steepness = (0.-ydum) / (0.-xdum)
					Zavgtargetarr[j] = 1./(1.+steepness)
					Gammatargetarr[j] = ydum * Zavgtargetarr[j]
				elif typefit.find("power")>-1 :
					if typefit=='power' :
						initvararr 	= [-5 , -2 ]
						initvararr2	= [ 50,  1 ]
						sigma		= None 
					elif typefit=='powerplusconst' :
						initvararr 	= [-0 , -0 , -0.000 ]
						initvararr2	= [ 50,  1 , 10.0 ]
						if args.fitweight0 : 
							sigma		= [1./nfit]+[1]*(nfit-1)
						else :
							sigma		= None
					popt, pcov = curve_fit(fitfunc, xraw, yraw, bounds=( initvararr , initvararr2 ) , sigma=sigma )
					if args.fitoffset : 
						if logplot : 
							xfit = np.logspace( np.log10(xdat[0]),  np.log10(xraw[-1]) , num=nfitresult ) 
							xfit2= np.logspace( np.log10(xraw[-1]),  np.log10(xdat[-1]) , num=nfitresult ) 
						else : 
							xfit = np.linspace( xdat[0],  xraw[-1] , nfitresult ) 
							xfit2= np.linspace( xraw[-1], xdat[-1] , nfitresult ) 
					else : 
						if logplot : 
							xfit = np.logspace( np.log10(xraw[0]),  np.log10(xraw[-1]) , num=nfitresult ) 
							xfit2= np.logspace( np.log10(xraw[-1]),  np.log10(xdat[-1]) , num=nfitresult ) 
						else : 
							xfit = np.linspace( xraw[0],  xraw[-1] , nfitresult ) 
							xfit2= np.linspace( xraw[-1], xdat[-1] , nfitresult ) 
							if args.xfitdata0 : 
								xfit = np.linspace( 0,  xraw[-1] , nfitresult ) 
					yfit = fitfunc(xfit,*popt) 
					yfit2= fitfunc(xfit2,*popt) 
					print "XFIT : ", xfit[ [0,1,-2-1] ]
					print "YFIT : ", yfit[ [0,1,-2-1] ]
					print "FIT({};{}) : ".format( formfunc, typefit ),
					for jj in range(nvarfunc) :
						print "{}".format(popt[jj]),
					print ""
					Zarr[j] = popt[1]
					(quantclass.power)[j]		= popt[1]
					if typefit=='powerplusconst' : 
						(quantclass.powerfitconst)[j]	= popt[2]
					elif typefit=='power' : 
						(quantclass.powerfitconst)[j]	= -999.
					fitls = "k:"
					if logplot :
						ax.loglog( xfit, yfit, fitls , lw=1 )
					else :
						ax.plot( xfit, yfit, fitls , lw=1 )
					fitls2= "k:"
					if logplot :
						ax.loglog( xfit,  yfit,  fitls  , lw=1 )
						ax.loglog( xfit2, yfit2, fitls2 , lw=1 )
					else :
						ax.plot( xfit,  yfit,  fitls , lw=1 )
						ax.plot( xfit2, yfit2, fitls2, lw=1 )

		if args.linearfit : 
			with open( "tmplinearfitdat", 'a' ) as fwfit :
				fwfit.write( "{:3d} ".format(hobj.count) )
				for j in chdat : 
					fwfit.write( "{} ".format(j) )
				for j in chdat : 
					fwfit.write( "{:12f} ".format(1./Zarr[j]) )
				for j in chdat : 
					fwfit.write( "{:12f} ".format( effmasserrarr[j] ) )
				fwfit.write( "\t" ) 
				try : 
					outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT ] 
				except : 
					outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.filling_0 , hobj.nOcculattAra ,hobj.IT ] 
				for it in outarr : 
					fwfit.write( "{}\t".format(it) )
				fwfit.write( "\n" ) 
			with open( "effmass_linearfitdat", 'a' ) as fwfit :
				fwfit.write( "{:3d} ".format(hobj.count) )
				for j in chdat : 
					fwfit.write( "{} ".format(j) )
				for j in chdat : 
					fwfit.write( "{:12f} ".format(1./Zarr[j]) )
				fwfit.write( "\t" ) 
				try : 
					outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT, hobj.Treal ] 
				except : 
					outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.filling_0 , hobj.nOcculattAra ,hobj.IT, hobj.Treal ] 
				for it in outarr : 
					fwfit.write( "{}\t".format(it) )
				fwfit.write( "\n" ) 
			if (len(chdat)==3) or (len(chdat)==2) : 
				if os.system("test -e "+hobj.tdir+"/result/u{:.3f}".format(hobj.UF) ) : 
					os.system( "mkdir "+hobj.tdir+"/result/u{:.3f} -p".format(hobj.UF) )
				fname 				= hobj.tdir+"/result/u{:.3f}/mass_fitIm".format(hobj.UF)
				if basis.find("t")>-1 :	fname	= fname+"_t2g"
				elif basis.find("e")>-1 :	fname	= fname+"_eg"
				if args.indimp :	fname	= fname+"_{}".format(int(args.indimp))
				fname	= fname+".dat"
				with open( fname, 'w' ) as fwfit :
					fwfit.write( "{:3d} ".format(hobj.count) )
					for j in chdat : 
						fwfit.write( "{} ".format(j) )
					for j in chdat : 
						fwfit.write( "{:12f} ".format(1./Zarr[j]) )
					for j in chdat : 
						fwfit.write( "{:12f} ".format( effmasserrarr[j] ) )
					fwfit.write( "\t" ) 
					try : 
						outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT ] 
					except : 
						outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.filling_0 , hobj.nOcculattAra ,hobj.IT ] 
					for it in outarr : 
						fwfit.write( "{}\t".format(it) )
					fwfit.write( "\n" ) 
			with open( "tmplinearfitavg_b512", 'a' ) as fwfitavg :
				for j in chdat : 
					fwfitavg.write( "{} ".format(j) )
				for j in chdat : 
					fwfitavg.write( "{:12f} ".format(1./Zavgtargetarr[j]) )
				try : 
					outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra, hobj.IT , tempdum ] 
				except : 
					outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.filling_0 , hobj.nOcculattAra ,hobj.IT , tempdum ] 
				for it in outarr : 
					try : 
						fwfitavg.write( "{:8f}\t".format(it) )
					except : 
						fwfitavg.write( "{}\t".format(it) )
				fwfitavg.write( "\n" )
		if args.fitting and typefit.find("power")>-1 : 
			with open( "tmppowerfitdat", 'a' ) as fwfit :
				fwfit.write( "{:3d} ".format(hobj.count) )
				for j in chdat : 
					fwfit.write( "{} ".format(j) )
				for j in chdat : 
					fwfit.write( "{:12f} ".format(Zarr[j]) )
				for j in chdat : 
					fwfit.write( "{:12f} ".format( (quantclass.powerfitconst)[j] ) )
				fwfit.write( "\t" ) 
				try : 
					outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT ] 
				except : 
					outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.filling_0 , hobj.nOcculattAra ,hobj.IT ] 
				for it in outarr : 
					fwfit.write( "{}\t".format(it) )
				fwfit.write( "\n" ) 
			if (len(chdat)==3) or (len(chdat)==2) : 
				if os.system("test -e "+hobj.tdir+"/result/u{:.3f}".format(hobj.UF) ) : 
					os.system( "mkdir "+hobj.tdir+"/result/u{:.3f} -p".format(hobj.UF) )
				fhead	= "power"
				if args.fitting=="powerplusconst" : fhead = "powerplusconst"
				fname 				= hobj.tdir+"/result/u{:.3f}/{}_fitSelfIm".format(hobj.UF,fhead)
				if basis.find("t")>-1 :	fname	= fname+"_t2g"
				elif basis.find("e")>-1 :	fname	= fname+"_eg"
				if args.indimp :	fname	= fname+"_{}".format(int(args.indimp))
				fname	= fname+".dat"
				with open( fname, 'w' ) as fwfit :
					fwfit.write( "{:3d} ".format(hobj.count) )
					for j in chdat : 
						fwfit.write( "{} ".format(j) )
					for j in chdat : 
						fwfit.write( "{:12f} ".format( Zarr[j]) )
					for j in chdat : 
						fwfit.write( "{:12f} ".format( (quantclass.powerfitconst)[j] ) )
					fwfit.write( "\t" ) 
					try : 
						outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT ] 
					except : 
						outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.filling_0 , hobj.nOcculattAra ,hobj.IT ] 
					for it in outarr : 
						fwfit.write( "{}\t".format(it) )
					fwfit.write( "\n" ) 

		if args.fitting and ( (typefit=='poly2') or (typefit=='poly3') or (typefit=='poly4') ):
			with open( "tmppolyfitdat", 'a' ) as fwfit :
				fwfit.write( "{:3d} ".format(hobj.count) )
				for j in chdat : 
					fwfit.write( "{} ".format(j) )
				for j in chdat : 
					fwfit.write( "{:12f} ".format( (quantclass.effmass)[j] ) )
				for j in chdat : 
					fwfit.write( "{:12f} ".format( (quantclass.erreffmass)[j] ) )
				fwfit.write( "\t" ) 
				try : 
					outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT ] 
				except : 
					outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.filling_0 , hobj.nOcculattAra ,hobj.IT ] 
				for it in outarr : 
					fwfit.write( "{}\t".format(it) )
				fwfit.write( "\n" ) 
			if (len(chdat)==3) or (len(chdat)==2) : 
				if os.system("test -e "+hobj.tdir+"/result/u{:.3f}".format(hobj.UF) ) : 
					os.system( "mkdir "+hobj.tdir+"/result/u{:.3f} -p".format(hobj.UF) )
				fname 				= hobj.tdir+"/result/u{:.3f}/{}_fitSelfIm".format(hobj.UF,typefit)
				if basis.find("t")>-1 :	fname	= fname+"_t2g"
				elif basis.find("e")>-1 :	fname	= fname+"_eg"
				if args.indimp :	fname	= fname+"_{}".format(int(args.indimp))
				fname	= fname+".dat"
				with open( fname, 'w' ) as fwfit :
					fwfit.write( "{:3d} ".format(hobj.count) )
					for j in chdat : 
						fwfit.write( "{} ".format(j) )
					for j in chdat : 
						fwfit.write( "{:12f} ".format( (quantclass.effmass)[j] ) )
					for j in chdat : 
						fwfit.write( "{:12f} ".format( (quantclass.erreffmass)[j] ) )
					fwfit.write( "\t" ) 
					try : 
						outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT ] 
					except : 
						outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.filling_0 , hobj.nOcculattAra ,hobj.IT ] 
					for it in outarr : 
						fwfit.write( "{}\t".format(it) )
					fwfit.write( "\n" ) 
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		#ax.axhline( 0, 0,1 , color='g', dashes=dashpt , lw=0.5 )
		if args.datyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.datyrange ) 
			else :
				axSetRange( ax , args.datyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 
		if args.datxrange : 
			axSetRange( ax , args.datxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])
		if args.write : 
			dum = np.array( [ xdat, ydat] ).transpose()
			np.savetxt("tmpdata", dum ) 
		if args.writegrid : 
			dumn = int(args.writegrid) 
			dum = np.array( [ xdat[:dumn], ydat[:dumn]] ).transpose()
			np.savetxt("tmpdatagrid_n{}_D{:.2f}".format(dumn,hobj.D), dum ) 
	def getplotselfiwoffdiag(self, ax, dat , basis, nbasis, chdat ,hobj , args, T=False ) :
		self.getplotselfiwoffdiag_lc(ax, dat , basis, nbasis, chdat ,hobj , args , "", T=T ) 
	def getplotselfiwoffdiag_lc(self, ax, dat , basis, nbasis, chdat ,hobj , args , lc , T=False ) :
		self.getselfiw( dat, nbasis, hobj, args );
		self.gloc = self.arrdat
		nplot = len( self.wn )
		ncomp = len( self.gloc )
		print nplot, " points(lines)"
		print ncomp, " components(of matrix)"
		self.glocReOffDiag = [ range(nbasis) for j in range(nbasis) ]
		self.glocImOffDiag = [ range(nbasis) for j in range(nbasis) ]
		self.glocMat = self.compMat
		print "0-th matrix :"
		for j1 in range(nbasis) :
			for j2 in range(nbasis) :
				print "{:20.8f}\t".format( self.glocMat[j1][j2][0] ),
			print ""
		if T is not False : 
			self.glocMat = transfMjefftot2g( self.glocMat, T ) 
			basis = args.transform
			print "0-th matrix (transformed) :"
			for j1 in range(nbasis) :
				for j2 in range(nbasis) :
					print "{:20.8f}\t".format( self.glocMat[j1][j2][0] ),
				print ""
		
		for k in chdat :
			i , j = [ k[0], k[1] ]
			self.glocReOffDiag[i][j] = self.glocMat[i][j].real
			self.glocImOffDiag[i][j] = self.glocMat[i][j].imag
			re = np.array( self.glocReOffDiag[i][j] ) 
			im = np.array( self.glocImOffDiag[i][j] ) 
			xdat = self.wn
			ydat = im
			if args.realpart : ydat = re

			ls = [ "-" ] + [ "-." ] + [ "--" ] *3 +  [ "-" ] + [ "-." ] + [ "--" ] *3  
			ms = [ "x" ] + [ "3" ] + [ ">" ] + [ "^" ] + [ "x" ] + [ "+" ] 
			ylab = "Im$\Delta(i\omega_n)$"
			xlab = "$\omega_n$"
			if args.minus :
				ydat = -ydat
				ylab = "-"+ylab
			if args.minusx :
				xdat = -xdat
			if args.selffactor :
				ydat = ydat*float(args.selffactor)
				ylab = args.selffactor+"*"+ylab
			if args.selfoffset :
				ydat = ydat+float(args.selfoffset)
				ylab = ylab +"+"+ args.selfoffset
			if args.trans :
				ax.plot( ydat, xdat, ms[i]+ls[j], label="{}{}".format(i,j)  , linewidth=1.5 )
				if i==0 and args.fillb :
					ax.fill( ydat, xdat, ms[i]+ls[j] ,alpha=0.2 )
				tmplab = xlab
				xlab   = ylab
				ylab   = tmplab
			else :
				ax.plot( xdat, ydat, ms[i]+ls[j], label="{}{}".format(i,j) )
				#ax.legend()
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		#ax.axhline( 0, 0,1 , color='g', dashes=dashpt , lw=0.5 )
		if args.datyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.datyrange ) 
			else :
				axSetRange( ax , args.datyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 
		if args.datxrange : 
			axSetRange( ax , args.datxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])
	def gethybw(self, dat , basis, nbasis, hobj , args ) :
		hobj.readEclustersoc()		#saved in hobj.ecluster[i][j]
		if( dat.find(".dat") > -1 ) :
			datfile = dat
		else : 
		 	if basis.find("t2g") > -1  : 
				datfile = "{}/lattice/tdelta_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
				rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		 	elif basis.find("j") > -1  : 
				try : 
					datfile = "{}/lattice/u{:.3f}/delta_Nplot{}_ep{:.3f}_{}th.dat".format(dat,hobj.UF,hobj.Nplot,hobj.epsilon,hobj.count)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
				except :
					datfile = "{}/lattice/delta_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		print "Reading :", datfile
		self.wn   = rawdatarr.transpose()[0]
		self.arrdat = rawdatarr.transpose()[1:]
		self.compMat = complexArrDatReturn( self.arrdat, nbasis ) 
	def gethybiw(self, dat , nbasis, hobj , args ) :
		hobj.readEclustersoc()		#saved in hobj.ecluster[i][j]
		if( dat.find(".dat") > -1 ) :
			datfile = dat
		else : 
			if args.count :	count = args.count
			else :		count = hobj.count
			datfile = "{}/ongoing/bathinverse{}_u{:.2f}_{}th.txt".format(dat, hobj.ni-1, hobj.UF, count )
		print "Reading :", datfile
		rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		self.wn   = rawdatarr.transpose()[0]
		self.arrdat = rawdatarr.transpose()[1:]
		self.compMat = complexArrDatReturn( self.arrdat, nbasis ) 
		for i in range(nbasis) : 		#getting the pure Hybridization function from the impurity green function_0.
			for j in range(nbasis) :
				self.compMat[i][j] = - self.compMat[i][j] - hobj.eclustersoc[i][j]
			self.compMat[i][i] += 1j * self.wn + hobj.D + hobj.UF/2 + hobj.JT
		for i in range(4,6) : 		#getting the pure Hybridization function from the impurity green function_0.
			self.compMat[i][i] -= 1.5* hobj.S
	def getplothybiwoffdiag(self, ax, dat , basis, nbasis, chdat ,hobj , args ) :
		self.getplothybiwoffdiag_lc(ax, dat , basis, nbasis, chdat ,hobj , args , "" ) 
	def getplothybiwoffdiag_lc(self, ax, dat , basis, nbasis, chdat ,hobj , args , lc ) :
		hobj.readEclustersoc()
		if args.count :	count = args.count
		else :		count = hobj.count
		if( dat.find(".dat") > -1 ) :
			datfile = dat
		else : 
			datfile = "{}/ongoing/bathinverse{}_u{:.2f}_{}th.txt".format(dat, hobj.ni-1, hobj.UF, count )
		print "Reading :", datfile
		rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		self.wn   = rawdatarr.transpose()[0]
		#self.gloc = rawdatarr.transpose()[1:]
		gdatReOffDiag , gdatImOffDiag = matrixOffDiagReturnReIm( rawdatarr.transpose()[1:] , nbasis )
		nplot = len( self.wn )
		ncomp = len( self.gloc )
		print nplot, " points(lines)"
		print ncomp, " components(of matrix)"
		self.glocReOffDiag = [ range(nbasis) for j in range(nbasis) ]
		self.glocImOffDiag = [ range(nbasis) for j in range(nbasis) ]
		a = range( nplot ) 
		for k in chdat :
			i , j = [ k[0], k[1] ]
			self.glocReOffDiag[i][j] = gdatReOffDiag[i][j]
			self.glocImOffDiag[i][j] = gdatImOffDiag[i][j]
			re = np.array( self.glocReOffDiag[i][j] ) * -1 - hobj.eclustersoc[i][j].real 
			im = np.array( self.glocImOffDiag[i][j] ) * -1 - hobj.eclustersoc[i][j].imag + self.wn *int(i==j)
			try : 
				if args.hybnat : 
					mm=range(nbasis)
					Tnat = hobj.returnNatTransf()
					glocnat = [ [0 for ii in mm] for jj in mm ]
					for kk in mm : 
						for ll in mm : 
							for ii in mm : 
								for jj in mm : 
									glocnat[kk][ll] = glocnat[kk][ll] + Tnat[kk,ii]*(self.gloc[2*ii+2*nbasis*jj]+1j*self.gloc[2*ii+1+2*nbasis*jj]) * Tnat[ll,jj]
					print "glocnat shpe : ", np.shape(glocnat)
					eclsocnat = hobj.readEclustersocnat()
					re = np.array( glocnat[i][j].real ,dtype=float ) * -1 - eclsocnat[i][j].real
					im = np.array( glocnat[i][j].imag ,dtype=float ) * -1 - eclsocnat[i][j].imag + self.wn *int(i==j)
			except : pass
		
			lc = [ "-" ] + [ "-." ] + [ "--" ] *3
			lt = [ "x" ] + [ "3" ] + [ ">" ] + [ "^" ] + [ "x" ] + [ "+" ] 
			if args.trans :
				if args.minus :
					im = np.array( im , dtype=float) 
					ax.plot( -im, self.wn, lc[i]+lt[j], label="{}{}".format(i,j)  , linewidth=1.5 )
					xlab = "-Im$\Delta(i\omega_n)$"
				else :
					ax.plot( im, self.wn, lc[i]+lt[j], label="{}{}".format(i,j)  , linewidth=1.5 )
					xlab = "Im$\Delta(i\omega_n)$"
					if i==0 and args.fillb :
						ax.fill( im, self.wn, lc[i]+lt[j] ,alpha=0.2 )
				ylab = ""
			else :
				ax.plot( self.wn, im, lc[i]+lt[j], label="{}{}".format(i,j) )
				#ax.legend()
				ylab = "Im$\Delta(i\omega_n)$"
				xlab = "$\omega_n$"
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		#ax.axhline( 0, 0,1 , color='g', dashes=dashpt , lw=0.5 )
		if args.datyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.datyrange ) 
			else :
				axSetRange( ax , args.datyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 
		if args.datxrange : 
			axSetRange( ax , args.datxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])
	def getplothybiw_lc(self, ax, dat , basis, nbasis, chdat ,hobj , args , lc , laboverw=False , dashes=(None,None) , lcarr=None , realpart=False ) :
		hobj.readEclustersoc()
		if args.count :	count = args.count
		else :		count = hobj.count
		if( dat.find(".dat") > -1 ) :
			datfile = dat
		else : 
			datfile = "{}/ongoing/bathinverse{}_u{:.2f}_{}th.txt".format(dat, hobj.ni-1, hobj.UF, count )
		print "Reading :", datfile
		rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		self.wn   = rawdatarr.transpose()[0]
		self.gloc = rawdatarr.transpose()[1:]
		nplot = len( self.wn )
		ncomp = len( self.gloc )
		print nplot, " points(lines)"
		print ncomp, " components(of matrix)"
		self.glocReDiag = range(nbasis)
		self.glocImDiag = range(nbasis)
		a = range( nplot ) 
		for j in chdat :
			self.glocReDiag[j] = self.gloc[2*j   +2*nbasis*j]
			self.glocImDiag[j] = self.gloc[2*j+1 +2*nbasis*j]
			re = np.array( self.glocReDiag[j] ,dtype=float ) * -1 - hobj.eclustersoc[j][j].real
			im = np.array( self.glocImDiag[j] ,dtype=float ) * -1 - hobj.eclustersoc[j][j].imag + self.wn 
			
			try : 
				if args.hybnat : 
					mm=range(nbasis)
					Tnat = hobj.returnNatTransf()
					glocnat = [ [0 for ii in mm] for jj in mm ]
					for kk in mm : 
						for ll in mm : 
							for ii in mm : 
								for jj in mm : 
									glocnat[kk][ll] = glocnat[kk][ll] + Tnat[kk,ii]*(self.gloc[2*ii+2*nbasis*jj]+1j*self.gloc[2*ii+1+2*nbasis*jj]) * Tnat[ll,jj]
					print "glocnat shpe : ", np.shape(glocnat)
					eclsocnat = hobj.readEclustersocnat()
					re = np.array( glocnat[j][j].real ,dtype=float ) * -1 - eclsocnat[j][j].real
					im = np.array( glocnat[j][j].imag ,dtype=float ) * -1 - eclsocnat[j][j].imag + self.wn 
			except : pass
		
			if laboverw :	labelstr = laboverw
			else :		labelstr = blabel(basis,j)
			lt = [ "-" ]*2 + [ "-" ] *2 + [ ":" ] *2
			if lcarr : lc = lcarr[j]
			indreim = "Im"
			ydat = np.array( im , dtype=float) 
			if realpart :
				ydat = np.array( re , dtype=float) 
				indreim = "Re"
			if args.trans :
				if args.minus :
					ax.plot( -ydat, self.wn, lc+blt(basis,j) , label=labelstr , linewidth=1.5 , dashes=dashes )
					xlab = '-'+indreim+"$\Delta(i\omega_n)$"
				else :
					ax.plot( ydat, self.wn, lc+blt(basis,j) , label=labelstr , linewidth=1.5 , dashes=dashes )
					xlab = indreim+"$\Delta(i\omega_n)$"
					if j==0 and args.fillb :
						ax.fill( im, self.wn, lc+blt(basis,j) ,alpha=0.2 )
				ylab = ""
			else :
				print "LC : ", lc
				ax.plot( self.wn, ydat, lc+blt(basis,j) , label=labelstr , dashes=dashes)
				#ax.legend()
				ylab = indreim+"$\Delta(i\omega_n)$"
				xlab = "$\omega_n$"
			if len(chdat)<2 : ylab = ylab + "$^{{{}}}$".format(blabel(basis,j))
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		#ax.axhline( 0, 0,1 , color='g', dashes=dashpt , lw=0.5 )
		if args.datyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.datyrange ) 
			else :
				axSetRange( ax , args.datyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 
		if args.datxrange : 
			axSetRange( ax , args.datxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])
	def getplothybiw(self, ax, dat , basis, nbasis, chdat ,hobj , args ) :
		self.getplothybiw_lc(ax, dat , basis, nbasis, chdat ,hobj , args , "" ) 
	def getplothyb_lc(self, ax, dat , basis, nbasis, chdat ,hobj , args, lc , lcarr=None , laboverw=False , dashes=(None,None) , realpart=False, axarr=None ) :
		self.wn   = [""]*  hobj.Nplot
		if( dat.find(".dat") > -1 ) :
			datfile = dat
		else : 
		 	if basis.find("t") > -1  : 
		 		if basis.find("z") > -1  : 
					datfile = "{}/lattice/delta_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
				else : 
					datfile = "{}/lattice/tdelta_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		 	elif basis.find("j") > -1 or basis.find("e") > -1  : 
				with open( dat +'/parameters_u{:.2f}.dat'.format(hobj.UF) ) as fo : 
					lines = fo.readlines()
					maxiter = int( lines[-1].split()[0] )
				try : 
					datfile = "{}/lattice/u{:.3f}/delta_Nplot{}_ep{:.3f}_{}th.dat".format(dat,hobj.UF,hobj.Nplot,hobj.epsilon,maxiter)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
				except :
					datfile = "{}/lattice/delta_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		print "Reading :", datfile

		rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		self.wn   = rawdatarr.transpose()[0]
		print len( self.wn ), "points (lines/rows)"
		self.gloc = rawdatarr.transpose()[1:]
		print len( self.gloc ), "components (columns)"
		glocRe		= self.gloc[::2]
		glocIm		= self.gloc[1::2]
		glocComp	= glocRe + glocIm * 1j
		dim		= int(np.sqrt( len(glocComp) ))		# delta_{}.txt has NU by NU matrix
		if dim*dim != len(glocComp) :
			print "ERROR: dim*dim != len(glocComp) ", hobj.Ni*dim*dim , np.shape(glocComp)
			print "dim, shape(glocComp) ", dim , np.shape(glocComp)
			sys.exit(1)
		glocComp	= glocRe + glocIm * 1j
		glocComp	= glocComp.reshape( dim,dim,len(self.wn) ) 

		for j in chdat :
			if axarr is not None : ax = axarr[j]
			#self.glocReDiag[j] = gdatReDiag[j]
			#self.glocImDiag[j] = gdatImDiag[j]
			#re = np.array( self.glocReDiag[j] )
			#im = np.array( self.glocImDiag[j] )

			#tnc	= hobj.tNC
			re	= glocComp[j,j].real
			im	= glocComp[j,j].imag
		
			if laboverw :	labelstr = laboverw
			else :		labelstr = blabel(basis,j)
			lt = [ "-" ]*2 + [ "-" ] *2 + [ ":" ] *2
			lcj= lc
			if lcarr : lcj = lcarr[j] 
			print "lc %d :"%j, lcj
			xdat = self.wn
			ydat = np.array( im , dtype=float) 
			indreim = "Im"
			if realpart :
				ydat = np.array( re , dtype=float) 
				indreim = "Re"
			if args.trans :
				if args.selfm :
					ydat = -ydat 
					xlab = "-"+indreim+"$\Delta_{hyb}(\omega)$"
				else :
					xlab = indreim+"$\Delta_{hyb}(\omega)$"
				ylab = ""
				tmp  = ydat
				ydat = xdat 
				xdat = tmp
			else :
				#ax.legend()
				ylab = indreim+"$\Delta_{hyb}(\omega)$"
				xlab = "$\omega$"
			ax.plot( xdat, ydat, lcj+blt(basis,j) , label=labelstr , linewidth=1.5 , dashes=dashes )
			if j==0 and args.fillb : ax.fill( xdat, ydat, blt(basis,j) , alpha=0.2 )
			if len(chdat)<2 : ylab = ylab + "$^{{{}}}$".format(blabel(basis,j))
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		ax.axhline( 0, 0,1 , color='gray', dashes=dashpt , lw=0.5 )
		ax.axvline( 0, 0,1 , color='gray', dashes=dashpt , lw=0.5 )
		if args.hybxrange : 
			axSetRange( ax , args.hybxrange , 'x' ) 
		if args.hybyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.hybyrange ) 
			else :
				axSetRange( ax , args.hybyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])
	def getplothyb(self, ax, dat , basis, nbasis, chdat ,hobj , args) :
		self.getplothyb_lc(ax, dat , basis, nbasis, chdat ,hobj , args, "" )
	def getplothybre(self, ax, dat , basis, nbasis, chdat ,hobj , args) :
		self.wn   = [""]*  hobj.Nplot
		self.gloc = [""]*( hobj.Nplot )
		if( dat.find(".dat") > -1 ) :
			datfile = dat
		else : 
		 	if basis.find("t") > -1  : 
		 		if basis.find("z") > -1  : 
					datfile = "{}/lattice/delta_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
				else : 
					datfile = "{}/lattice/tdelta_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		 	elif basis.find("j") > -1 or basis.find("e") > -1  : 
				with open( dat +'/parameters_u{:.2f}.dat'.format(hobj.UF) ) as fo : 
					lines = fo.readlines()
					maxiter = int( lines[-1].split()[0] )
				try : 
					datfile = "{}/lattice/u{:.3f}/delta_Nplot{}_ep{:.3f}_{}th.dat".format(dat,hobj.UF,hobj.Nplot,hobj.epsilon,maxiter)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
				except :
					datfile = "{}/lattice/delta_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		print "Reading :", datfile
		gdatReDiag , gdatImDiag = matrixDiagReturnReIm( rawdatarr.transpose()[1:] , nbasis )
		self.wn   = rawdatarr.transpose()[0]
		print len( self.wn ), "points"
		self.glocReDiag = range(nbasis)
		self.glocImDiag = range(nbasis)
		print len( gdatImDiag ), "components"
		for j in chdat :
			self.glocReDiag[j] = gdatReDiag[j]
			self.glocImDiag[j] = gdatImDiag[j]
			re = np.array( self.glocReDiag[j] )
			im = np.array( self.glocImDiag[j] )
			ydat = re
		
			lc = [ "-" ]*2 + [ "-" ] *2 + [ ":" ] *2
			if args.trans :
				if args.selfm :
					ydat = np.array( ydat , dtype=float) 
					ax.plot( -ydat, self.wn, lc[j] , label=bolabel(basis,j) , linewidth=1.5 )
					xlab = "-Re$\Delta_{hyb}(\omega)$"
				else :
					ax.plot( ydat, self.wn, lc[j] , label=bolabel(basis,j) , linewidth=1.5 )
					xlab = "Re$\Delta_{hyb}(\omega)$"
				ylab = ""
			else :
				ax.plot( self.wn, ydat, lc[j] , label=blabel(basis,j) )
				#ax.legend()
				ylab = "Re$\Delta_{hyb}(\omega)$"
				xlab = "$\omega$"
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		ax.axhline( 0, 0,1 , color='gray', dashes=dashpt , lw=0.5 )
		if args.hybyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.hybyrange ) 
			else :
				axSetRange( ax , args.hybyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])
	def getplothyboffdiag(self, ax, dat , basis, nbasis, chdat ,hobj , args) :
		self.wn   = [""]*  hobj.Nplot
		self.gloc = [""]*( hobj.Nplot )
		print "Reading :", dat
		datfile = ""
		if( dat.find(".dat") > -1 ) :
			datfile = dat
		else : 
		 	if basis.find("t") > -1  : 
		 		if basis.find("z") > -1  : 
					datfile = "{}/lattice/delta_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
				else : 
					datfile = "{}/lattice/tdelta_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		 	elif basis.find("j") > -1  : 
				with open( dat +'/parameters_u{:.2f}.dat'.format(hobj.UF) ) as fo : 
					lines = fo.readlines()
					maxiter = int( lines[-1].split()[0] )
				try : 
					datfile = "{}/lattice/u{:.3f}/delta_Nplot{}_ep{:.3f}_{}th.dat".format(dat,hobj.UF,hobj.Nplot,hobj.epsilon,maxiter)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
				except :
					datfile = "{}/lattice/delta_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
					rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		 	#if basis.find("t") > -1  : 
			#	datfile = "{}/lattice/tdelta_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
		 	#elif basis.find("j") > -1  : 
			#	with open( dat +'/parameters_u{:.2f}.dat'.format(hobj.UF) ) as fo : 
			#		lines = fo.readlines()
			#		maxiter = int( lines[-1].split()[0] )
			#	try : 
			#		datfile = "{}/lattice/u{:.3f}/delta_ep{:.3f}_{}th.dat".format(dat,hobj.UF,hobj.epsilon,maxiter)
			#	except :
			#		datfile = "{}/lattice/delta_u{:.2f}_Nplot{}.dat".format(dat,hobj.UF,hobj.Nplot)
		print "Reading : ", datfile
		rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		self.wn   = rawdatarr.transpose()[0]
		print len( self.wn ), "points (lines/rows)"
		self.gloc = rawdatarr.transpose()[1:]
		print len( self.gloc ), "components (columns)"
		glocRe		= self.gloc[::2]
		glocIm		= self.gloc[1::2]
		glocComp	= glocRe + glocIm * 1j
		dim		= int(np.sqrt( len(glocComp) ))		# delta_{}.txt has NU by NU matrix
		if dim*dim != len(glocComp) :
			print "ERROR: dim*dim != len(glocComp) ", hobj.Ni*dim*dim , np.shape(glocComp)
			print "dim, shape(glocComp) ", dim , np.shape(glocComp)
			sys.exit(1)
		glocComp	= glocRe + glocIm * 1j
		glocComp	= glocComp.reshape( dim,dim,len(self.wn) ) 

		for k in chdat :
			i , j = [ k[0], k[1] ]
			#self.glocReOffDiag[i][j] = gdatReOffDiag[i][j]
			#self.glocImOffDiag[i][j] = gdatImOffDiag[i][j]
			#re = np.array( self.glocReOffDiag[i][j] )
			#im = np.array( self.glocImOffDiag[i][j] )
			re	= glocComp[i,j].real
			im	= glocComp[i,j].imag

			lc = [ "-" ] + [ "-." ] + [ "--" ] *3
			xdat = self.wn
			ydat = np.array( im , dtype=float) 
		
			if args.trans :
				if args.selfm :
					ydat = -ydat 
					xlab = "-Im$\Delta_{hyb}(\omega)$"
				else :
					xlab = "Im$\Delta_{hyb}(\omega)$"
				ylab = ""
				tmp  = ydat
				ydat = xdat 
				xdat = tmp
			else :
				#ax.legend()
				ylab = "Im$\Delta_{hyb}(\omega)$"
				xlab = "$\omega$"
			ax.plot( xdat, ydat, lc[i], label="{}{}".format(i,j) , linewidth=1.5 )
			#if j==0 and args.fillb : ax.fill( xdat, ydat, lc[j] , alpha=0.2 )
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		ax.axhline( 0, 0,1 , color='gray', dashes=dashpt , lw=0.5 )
		ax.axvline( 0, 0,1 , color='gray', dashes=dashpt , lw=0.5 )
		if args.hybxrange : 
			axSetRange( ax , args.hybxrange , 'x' ) 
		if args.hybyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.hybyrange ) 
			else :
				axSetRange( ax , args.hybyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])
		print "axaxis : ", ax.axis()
	def getplotdataiw_lc(self, ax, dat , basis, nbasis, chdat ,hobj , args , lc , dataname , functlab=False, varlab=False , nondiag=False , realdata=False , labelarr=False, logplot=False , marker=False , dashes=(None,None) ,
		       epsilon=False	, filesuffix=None ) :
		hobj.readEclustersoc()
		if args.count :	count = args.count
		else :		count = hobj.count
		if( dat.find(".dat") > -1 ) :
			datfile = dat
		else : 
		 	if filesuffix : 
				datfile = returnDatafname( dataname, dat, hobj, count , fsuffix= filesuffix , epsilon=epsilon )
		 	else : 
				datfile = returnDatafname( dataname, dat, hobj, count , epsilon=epsilon )
		print "Reading :", datfile
		rawdatarr = np.genfromtxt( datfile, dtype=float ) 
		self.wn   = rawdatarr.transpose()[0]
		self.gloc = rawdatarr.transpose()[1:]
		nplot = len( self.wn )
		ncomp = len( self.gloc )
		print nplot, " points(lines),\t",
		print ncomp, " components(of matrix)"
		datnbasis = int(np.sqrt(ncomp/2))
		if (datnbasis*datnbasis) != (ncomp/2) : print "ERROR :: Not square-matrix data." ; sys.exit(1)
		self.glocReDiag = range(datnbasis)
		self.glocImDiag = range(datnbasis)
		a = range( nplot ) 
		if realdata :	dtypelab = "Re"
		else :		dtypelab = "Im"
		if functlab==False : 
			functlab = "{}{}$(i\omega_n)$".format(dtypelab, dataname[:3])
		if varlab==False : 
			varlab   = "$\omega_n$"
		labftn   = nlabel 
		if labelarr :
			def labftn(b,j) :
				return labelarr[j]
		if args.leglabel : 
			labsuffix = filterlabel(args.leglabel) + "={}".format( getattr( hobj, args.leglabel ) )
			if args.leglabel.find("filling")>-1 :
				labsuffix = filterlabel(args.leglabel) + "={:.2f}".format( float(getattr( hobj, args.leglabel )) )
			elif args.leglabel.find("nOcculatt")>-1 :
				labsuffix = filterlabel(args.leglabel) + "={:.2f}".format( float(getattr( hobj, args.leglabel )) )
			def labftn(b,j) :
				return labsuffix
		k=0
		try : 
			if args.sum : 
				xDatSum = self.wn
				yDatSum = np.zeros( len(self.wn) ) 
				labsum  = "Sum["
		except : 
			setattr( args, "sum" , None )  
		w0datarr = range(nbasis)
		for j in chdat :
			if nondiag==False : 
				self.glocReDiag[j] = self.gloc[2*j   +2*datnbasis*j]
				self.glocImDiag[j] = self.gloc[2*j+1 +2*datnbasis*j]
				if args.offsetecluster : 
					offreal = hobj.eclustersoc[j][j].real
					offimag = hobj.eclustersoc[j][j].imag - self.wn
				else : 
					offreal = 0.
					offimag = 0.
				re = np.array( self.glocReDiag[j] ,dtype=float ) + offreal
				im = np.array( self.glocImDiag[j] ,dtype=float ) + offimag
				indj = j 
			else : 
				mu1, mu2 = [int(j[0]), int(j[1]) ] 
				print "cMM : ", chdat
				print "MM : ", 2*mu1   +2*datnbasis*mu2
				self.glocReDiag[k] = self.gloc[2*mu1   +2*datnbasis*mu2]
				self.glocImDiag[k] = self.gloc[2*mu1+1 +2*datnbasis*mu2]
				if args.offsetecluster : 
					offreal = hobj.eclustersoc[j][j].real
					offimag = hobj.eclustersoc[j][j].imag - self.wn
				else : 
					offreal = 0.
					offimag = 0.
				re = np.array( self.glocReDiag[k] ,dtype=float ) + offreal
				im = np.array( self.glocImDiag[k] ,dtype=float ) + offimag
				k+=1
				labftn = offdlabel
				indj = k 
			try : 
				if args.hybnat : 
					mm=range(datnbasis)
					Tnat = hobj.returnNatTransf()
					glocnat = [ [0 for ii in mm] for jj in mm ]
					for kk in mm : 
						for ll in mm : 
							for ii in mm : 
								for jj in mm : 
									glocnat[kk][ll] = glocnat[kk][ll] + Tnat[kk,ii]*(self.gloc[2*ii+2*datnbasis*jj]+1j*self.gloc[2*ii+1+2*datnbasis*jj]) * Tnat[ll,jj]
					print "glocnat shpe : ", np.shape(glocnat)
					eclsocnat = hobj.readEclustersocnat()
					re = np.array( glocnat[j][j].real ,dtype=float ) * -1 - eclsocnat[j][j].real
					im = np.array( glocnat[j][j].imag ,dtype=float ) * -1 - eclsocnat[j][j].imag + self.wn 
			except : pass
		
			lt = [ "-" ]*2 + [ "-" ] *2 + [ ":" ] *2
			if len(args.ddir) > 1  : 
				lcj =""
			else : 
				try : 
					lcj = "C{}".format( int(j) ) 
				except : 
					lcj = "C{}".format( int(j[0]) ) 
			dat = im
			if realdata : dat = re
			dat = np.array( dat , dtype=float) 
			ydat = dat
			if args.halffirstcomp :
				if nondiag==False and j==0 : ydat = ydat/2.
			xdat = self.wn
			if args.nomarker :	ms = ""
			elif args.diffmarker :	ms = markersdiff[indj]
			else :			ms = "."
			if marker : 		ms = marker
			if args.solidline :	ls = "-"
			else :			ls = blt(basis,indj)
			if args.sum : 
				yDatSum = yDatSum + ydat 
			if args.minus :
				ydat = -ydat
			if args.minusx :
				xdat = -xdat
			if args.trans :
				xdat = ydat
				ydat = self.wn
				xlab = functlab
				ylab = ""
				if indj==0 and args.fillb :
					ax.fill( dat, self.wn, lcj+ms+ls ,alpha=0.2 )
			else :
				ylab = functlab
				xlab = varlab
			if args.halfnratio :
				pass
			elif args.sum : 
				labsum = labsum + labftn(basis,j) + ','
			elif logplot : 
				ax.loglog(	xdat, ydat, lcj+ms+ls , label=labftn(basis,j) , dashes=dashes ) #, linewidth=1.5 )
			else : 
				ax.plot(	xdat, ydat, lcj+ms+ls , label=labftn(basis,j) , dashes=dashes ) #, linewidth=1.5 )
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
			if args.diag : 
				w0datarr[j] = ydat[0]
			if args.showfermidat : 
				print "Near Fermi-level [{}]: ".format(j)
				ndata = len(xdat)
				fw = open( 'datafermilevel', 'a' ) 
				if( abs(xdat[ndata/2-1]) < 1e-1 ) :
					print "{}\t{}\t{:19.16f} {:19.16f}".format( hobj.J, j, xdat[ndata/2-1], ydat[ndata/2-1] )
					print "{}\t{}\t{:19.16f} {:19.16f}".format( hobj.J, j, xdat[ndata/2  ], ydat[ndata/2  ] )
					print "{}\t{}\t{:19.16f} {:19.16f}".format( hobj.J, j, xdat[ndata/2+1], ydat[ndata/2+1] )
				else : 
					print "{}\t{}\t{:19.16f} {:19.16f} {:19.16f}".format( hobj.J, j, xdat[0], ydat[0] , (ydat[1]-ydat[0])/(xdat[1]-xdat[0]) )
					print "{}\t{}\t{:19.16f} {:19.16f} {:19.16f}".format( hobj.J, j, xdat[1], ydat[1] , (ydat[2]-ydat[1])/(xdat[2]-xdat[1]) )
					print "{}\t{}\t{:19.16f} {:19.16f} {:19.16f}".format( hobj.J, j, xdat[2], ydat[2] , (ydat[3]-ydat[2])/(xdat[3]-xdat[2]) )
					for aa in range(3) :
						fw.write( "{}\t{}\t{:19.16f} {:19.16f} {:19.16f}\n".format( hobj.J, j, xdat[aa], ydat[aa] , (ydat[aa+1]-ydat[aa])/(xdat[aa+1]-xdat[aa]) ) )
				fw.close()
		if args.showfermidat : 
			if nondiag==False : 
				if realdata : 
					gxdat = self.glocReDiag[0]
					gzdat = self.glocReDiag[2]
				else :
					gxdat = self.glocImDiag[0]
					gzdat = self.glocImDiag[2]
				if( abs(xdat[ndata/2-1]) < 1e-1 ) :
					#print "{}\t{}\t{:19.16f} {:19.16f} {:19.16f}".format( hobj.J, j, xdat[ndata/2-2], gzdat[ndata/2-2], gxdat[ndata/2-2]/2. )
					#print "{}\t{}\t{:19.16f} {:19.16f} {:19.16f}".format( hobj.J, j, xdat[ndata/2-2], gzdat[ndata/2-2], gxdat[ndata/2-2]/2. )
					#print "{}\t{}\t{:19.16f} {:19.16f} {:19.16f}".format( hobj.J, j, xdat[ndata/2-1], gzdat[ndata/2-1], gxdat[ndata/2-1]/2. )
					#print "{}\t{}\t{:19.16f} {:19.16f} {:19.16f}".format( hobj.J, j, xdat[ndata/2  ], gzdat[ndata/2  ], gxdat[ndata/2  ]/2. )
					#print "{}\t{}\t{:19.16f} {:19.16f} {:19.16f}".format( hobj.J, j, xdat[ndata/2+1], gzdat[ndata/2+1], gxdat[ndata/2+1]/2. )
					navg = 6 # 94 (omega=0.18) # 190 (omega=0.35) # 150 (omega=0.3) # 230(omega=0.45) # 26 (omega=0.05) for orginal figure 190617
					wavg =   xdat[ndata/2-navg/2:ndata/2+navg/2]
					zavg =  gzdat[ndata/2-navg/2:ndata/2+navg/2]
					xavg =  gxdat[ndata/2-navg/2:ndata/2+navg/2]/2.
					print "data-length : ", len(wavg)
					print "wavg, zavg, xavg : "
					print wavg , wavg.mean()
					print zavg , zavg.mean()
					print xavg , xavg.mean()
					zstd = zavg.std()
					xstd = xavg.std()
					ravg = xavg/zavg
					rstd = ravg.std()
					ravg = ravg.mean()
					wavg = wavg.mean()
					zavg = zavg.mean()
					xavg = xavg.mean()
					for nn in range(navg) :
						n0 = ndata/2+nn-navg/2
						print "{}\t{}\t{:19.16f} {:19.16f} {:19.16f} {:19.16f}".format( hobj.J, j, xdat[n0], gzdat[n0], gxdat[n0]/2. , (gxdat[n0]/2.)/gzdat[n0] )
					fw = open( 'datafermilevelavg_b{}'.format(hobj.beta), 'a' ) 
					print dataname
					fw.write( "{}\t{}\t{:19.16f} {:19.16f} {:19.16f} {:19.16f} {:19.16f} {:19.16f} {:19.16f} {:.2f}\n".format( hobj.J, dataname, wavg, zavg, xavg, ravg , zstd, xstd, rstd , hobj.Treal) )
					fw.close()
				else :
					print "{}\t{}\t{:19.16f} {:19.16f} {:19.16f}".format( hobj.J, j, xdat[0], gzdat[0], gxdat[0]/2. )
					print "{}\t{}\t{:19.16f} {:19.16f} {:19.16f}".format( hobj.J, j, xdat[1], gzdat[1], gxdat[1]/2. )
					print "{}\t{}\t{:19.16f} {:19.16f} {:19.16f}".format( hobj.J, j, xdat[2], gzdat[2], gxdat[2]/2. )
					fw = open( 'datafermilevel_J{}_b{}'.format(hobj.J,hobj.beta), 'a' ) 
					for aa in range(3) :
						fw.write( "{}\t{}\t{:19.16f} {:19.16f} {:19.16f} {:19.16f}\n".format( hobj.J, j, xdat[aa], gzdat[aa], gxdat[aa]/2., gxdat[aa]/2./gzdat[aa] ) )
					fw.close()
		if args.sum : 
			xdat = xDatSum
			ydat = yDatSum
			if args.minus :
				ydat = -ydat
			if args.minusx :
				xdat = -xdat
			if args.trans :
				ydatdum = ydat
				xdat = ydatdum
				ydat = xdat
			if labelarr : 
				labsum = labelarr[0]
			else : 
				labsum = labsum[:-1]+"]"
			if logplot : 
				ax.loglog( xdat, ydat, lc+"-"+ms , label=labsum , dashes=dashes ) #, linewidth=1.5 )
			else : 
				ax.plot( xdat, ydat, lc+"-"+ms , label=labsum , dashes=dashes ) #, linewidth=1.5 )
			ax.set_ylabel(r'{}'.format(ylab) )#, fontsize = args.fs )
			ax.set_xlabel(r'{}'.format(xlab) )#, fontsize = args.fs )
		if args.halfnratio : 
			gxdat = self.glocReDiag[0]
			gzdat = self.glocReDiag[2]
			xdat  = self.wn
			ydat  = gxdat/2./gzdat
			leglab = None
			if args.leglabel : leglab =  "${}={}$".format(args.leglabel ,  getattr( hobj, args.leglabel ) ) 
			if args.nomarker : 	ms = ""
			elif marker : 		ms = marker
			else : 			ms = "."
			if logplot : 
				ax.loglog(	xdat, ydat, ms+"-", label=leglab , dashes=dashes )
			else : 
				ax.plot(	xdat, ydat, ms+"-", label=leglab , dashes=dashes )
			ax.axhline(	1, 0, self.wn[-1],  ls=":", color='gray' ) 
		else : 
			ax.axhline( 0, 0,1 , color='gray')
		if args.write :
			print "W0DATARR : ", w0datarr 
			with open("tmpwn0data_"+dataname,'a') as fw :
				for j in chdat : 
					fw.write( "{} ".format(j) )
				for j in chdat : 
					fw.write( "{:12f} ".format(float(w0datarr[j])) )
				for j in chdat : 
					fw.write( "{:12f} ".format(float(w0datarr[j])/float(w0datarr[chdat[-1]])) )
				fw.write( "\t" ) 
				outarr = [ hobj.UF, hobj.S, hobj.J, hobj.D, hobj.beta , hobj.fillimptot , hobj.nOcculattAra ,hobj.IT , hobj.degen_0 ] 
				for it in outarr : 
					fw.write( "{}\t".format(it) )
				fw.write( "\n" ) 
		ax.tick_params( axis="both", direction="in" )#, labelsize = args.fsm ) 
		#ax.tick_params( axis="y", labelleft="off" )
		dashpt = [1,1]
		if args.datyrange : 
			if args.xxrange :
				axSetAxis( ax , args.xxrange , args.datyrange ) 
			else :
				axSetRange( ax , args.datyrange , 'y' ) 
		elif args.xxrange :
			axSetRange( ax , args.xxrange , 'x' ) 
		if args.datxrange : 
			axSetRange( ax , args.datxrange , 'x' ) 

		if args.dxtick : 
			ax.set_xticks(ax.get_xticks()[1::2])
		if args.dytick : 
			ax.set_yticks(ax.get_yticks()[1::2])

def findLocalMaxSave( Z , fdir , hobj, fhead='fs' , XY=None ) :
	xdim, ydim = Z.shape
	zmax = 0. 
	zminpeak = -99. 
	fmax = open( fdir + "/{}_localmax_ep{:g}_br{:g}.dat".format(fhead, hobj.epsilon, hobj.broadening) , 'w' ) 
	print "Writing : ", fdir + "/{}_localmax_ep{:g}_br{:g}.dat".format(fhead, hobj.epsilon, hobj.broadening)
	if XY is None :
		xflat	= np.array(range(xdim),dtype=float) / (xdim-1)
		yflat	= np.array(range(ydim),dtype=float) / (ydim-1)
		X, Y = np.meshgrid(xflat, yflat)
	else : 
		X, Y = XY
	for x in range(xdim) :
		dspacing	= 1
		da		= dspacing
		r = x + da
		l = x - da
		if  x < da :
			r = x+da 
			l = xdim-1 -(da-1)
		elif x > xdim-1-da :
			r = 0 +(da-1)
			l = x-da
		for y in range(ydim) :
			if Z[x][y] > zmax  :
				zmax = Z[x][y]
			u = y + da
			d = y - da
			if  y < da :
				u = y+da 
				d = ydim-1 -(da-1)
			elif y > ydim-1-da :
				u = 0 +(da-1)
				d = y-da
			dz = 0 
			if Z[x][y] > Z[r][y] :
				if Z[x][y] > Z[l][y] :
					dz += 1
			if Z[x][y] > Z[x][u] :
				if Z[x][y] > Z[x][d] :
					dz += 1
			if Z[x][y] > Z[r][u] :
				if Z[x][y] > Z[l][d] :
					dz += 1
			if Z[x][y] > Z[l][u] :
				if Z[x][y] > Z[r][d] :
					dz += 1
			if dz > 0 : 
				#if Z[x][y] > 2.1 :
				fmax.write( "{}\t{}\t{}\n".format( X[x][y], Y[x][y] , Z[x][y] )  )
				if zminpeak < Z[x][y] :
					zminpeak = Z[x][y]
	fmax.close()
	fmaxinfo = open( fdir + "/{}_info_localmax_ep{:g}_br{:g}.dat".format(fhead, hobj.epsilon, hobj.broadening) , 'w' ) 
	fmaxinfo.write( "maximum_value\t{}\n".format(zmax) ) 
	fmaxinfo.write( "minimum_value_of_peak\t{}\n".format(zminpeak) ) 
	print "maximum_value\t{}\n".format(zmax) 
	print "minimum_value_of_peak\t{}\n".format(zminpeak) 
	fmaxinfo.close()

def plotLocalMax( fdir , ax , zthratio , hobj , args , fhead='fs' ) :
	fmaxinfo = open( fdir + "/{}_info_localmax_ep{:g}_br{:g}.dat".format(fhead, hobj.epsilon, hobj.broadening) )
	for lines in fmaxinfo.readlines() :
		if lines.find("maximum_value") > -1 :
			zmax = float( lines.split()[1] ) 
			zth = zthratio * zmax 
	fmaxinfo.close()

	fmaxdata = np.genfromtxt( fdir + "/{}_localmax_ep{:g}_br{:g}.dat".format(fhead, hobj.epsilon, hobj.broadening) , dtype=float )
	fmaxdata = fmaxdata.transpose()
	filtarr = fmaxdata[2]>zth
	x = fmaxdata[0][filtarr]
	y = fmaxdata[1][filtarr]

	if args.locmaxcolor :	lmcolor = args.locmaxcolor 
	else : 			lmcolor = 'k'
	ax.plot( x, y , ".", color=lmcolor, lw = 1, ms=1 )
	#ax.plot( fmaxdata[0], fmaxdata[1] , "k.", lw = 1)
	#ax.axis('equal')
	if args.label :
		ax.text( 0.8, 1.02, "{}".format(args.label) )
	if args.notitle : pass
	else : 
		ax.set_title( "threshold={:.3f} (={}*max)".format(zth,zthratio) , fontsize=8 )
	if args.fs : ax.axis('scaled')
	if args.fs : ax.axis([0,1,0,1])
	return zmax 

def saveLinePlot( fdir, lineind, lineval , x, y , hobj, fhead='fs' ) :
	fdname = fdir + "/lattice/vdx/{}_ep{:g}_br{:g}_along_".format(fhead, hobj.epsilon, hobj.broadening)+lineind+"_{}.dat".format(lineval)
	print "saved in : ", fdname
	fw = open( fdname , 'w' )
	#ax.plot( x, y , "r-o", lw = 1, markersize = 1.5 )
	for i in range( len(x) ) :
		fw.write( "{}\t{}\n".format( x[i], y[i] )  ) 
	fw.close()

def refineLineArg( mlineargs ) :
	mlinearg = mlineargs.split('_')
	lineargs = [] 
	for arg in mlinearg[1:] :
		lineargs.append( mlinearg[0]+"_"+arg ) 
	return  lineargs

def refineLinePlot( linearg , fs ) : #, axl ) :
	if linearg.find("x")>-1 : 
	        lineaxis = "x_"
	        nx0 = round( float(linearg.split('_')[1]) * (fs.nx-1) )
	        n0 = int(nx0)
	        print "nx0 : ", nx0 
	        xdat = fs.Z.transpose()
	        ydat = fs.Y.transpose()
	        rdat = fs.X.transpose()
	        datsection = rdat[n0][0]
		#axl.axvline( rdat[n0][0] , 0,1 , color='w', dashes=(2,1) , lw=0.5) 
	elif linearg.split('_')[0].find('y') > -1 :
	        lineaxis = "y_"
	        ny0 = round( float(linearg.split('_')[1]) * (fs.ny-1) )
	        n0 = int(ny0)
	        print "ny0 : ", ny0 
	        xdat = fs.X
	        ydat = fs.Z
	        rdat = fs.Y
		datsection = rdat[n0][0]
	        #axl.axhline( fs.Y[n0][0] , 0,1 , color='w', dashes=(2,1) , lw=0.5)
	print "xd : ", xdat[n0][0], "...",  xdat[n0][-1], "|",
	print "yd : ", ydat[n0][0], "...",  ydat[n0][-1], "|",
	print "rd : ", rdat[n0][0], "...",  rdat[n0][-1], "|"
	return n0, xdat, ydat, datsection

def plotIndLine( ax, datsection, axis, dash , c ) :
	if axis.find("y")>-1  :
		ax.axhline( datsection, 0,1 , color=c, lw=0.5 ) 
	if axis.find("x")>-1  :
		ax.axvline( datsection, 0,1 , color=c, lw=0.5 ) 
def setylimInd( ax, axis ,zmax ) :
	if axis.find("y")>-1  :
		ax.set_ylim( 0,zmax ) 
	if axis.find("x")>-1  :
		ax.set_xlim( 0,zmax ) 

class sympartner : 
	def __init__(self, basis ) :
		self.basis   = basis
		self.partner = []
		self.npartner = 0
	def add ( self, pair ) :
		self.partner.append( pair ) 
		self.npartner += 1

	
class bathsym : 
	def __init__( self ) :
		self.timerevt2g	=  sympartner( "t2g" ) 
		self.timerevt2g.add( [0,1,1] ) 
		self.timerevt2g.add( [1,0,1] ) 
		self.timerevt2g.add( [2,3,1] ) 
		self.timerevt2g.add( [3,2,1] ) 
		self.timerevt2g.add( [4,5,1] ) 
		self.timerevt2g.add( [5,4,1] ) 
		self.timerevj	=  sympartner( "j" ) 
		self.timerevj.add( [0,3,1] ) 
		self.timerevj.add( [1,2,1] ) 
		self.timerevj.add( [2,1,1] ) 
		self.timerevj.add( [3,0,1] ) 
		self.timerevj.add( [4,5,1] ) 
		self.timerevj.add( [5,4,1] ) 
		self.cubict2g	=  sympartner( "t2g" ) 
		self.cubict2g.add( [0,0,-1] ) 
		self.cubict2g.add( [1,1,-1] ) 
		self.cubict2g.add( [2,4,-1] ) 
		self.cubict2g.add( [3,5,-1] ) 
		self.cubict2g.add( [4,2,1] ) 
		self.cubict2g.add( [5,3,1] ) 

class transfMat : 
	def __init__(self ) :
		self.matT = np.matrix( [
		                [ 0.,                   0.,             np.sqrt(1./2.), 0.,                     np.sqrt(1./2.)*1j ,     0.                      ],  
		                [ np.sqrt(2./3.),       0.,             0.,             -1./np.sqrt(6.),        0.,                     -1/np.sqrt(6.)*1j       ],  
		                [ 0.,                   np.sqrt(2./3.), 1./np.sqrt(6.), 0.,                     -1./np.sqrt(6.)*1j,     0.                      ],  
		                [ 0.,                   0.,             0.,             1./np.sqrt(2.),         0.,                     -1./np.sqrt(2.)*1j      ],  
		                [ 1./np.sqrt(3.),       0.,             0.,             1./np.sqrt(3.),         0.,                     1./np.sqrt(3.)*1j       ],  
		                [ 0.,                   -1./np.sqrt(3.),        1./np.sqrt(3.), 0.,                     -1./np.sqrt(3.)*1j,     0.              ]] ) 
		self.matTi = np.linalg.inv( self.matT ) 
		self.matTt  = self.matT.transpose()
		self.matTti = np.linalg.inv( self.matTt ) 
		self.matTc  = self.matT.conjugate()
	def vecjtot2gNp (self, V ) :
		resV = np.dot( V , self.matT ) 
		return resV
	def matjtot2gNp (self, M ) :
		resM = np.dot( M , self.matTc ) 
		return np.dot( self.matTt , resM ) 

def complexMatDatFilt( M , nbasis )  :
	Mre = np.zeros( ( nbasis, nbasis ) )
	Mim = np.zeros( ( nbasis, nbasis ) )
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			Mre[mu][nu] = M[ 2*nu + 2*nbasis*mu ]
			Mim[mu][nu] = M[ 2*nu + 2*nbasis*mu + 1 ]
	return Mre, Mim

def complexMatDatReturn( M , nbasis)  :
	Mcomp = [ range(nbasis) for j in range(nbasis) ]
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			Mcomp[mu][nu] = M[ 2*nu + 2*nbasis*mu ] +  M[ 2*nu + 2*nbasis*mu +1 ]* 1j 
	return Mcomp

def complexArrDatReturn( M , nbasis)  : # M is an array having Re,Im separately.
	Mcomp = [ range(nbasis) for j in range(nbasis) ]
	for mu in range(nbasis) :
		for nu in range(nbasis) :
			try : 
				Mcomp[mu][nu] = M[ 2*nu + 2*nbasis*mu ] +  M[ 2*nu + 2*nbasis*mu +1 ]* 1j 
			except :
				print "ERROR:: nbasis : ", nbasis 
				exit(1)
	return Mcomp

def complexWholetxtReturnRemove1line( textarr , nbasis)  :
	dat = np.genfromtxt( textarr[1:] , dtype=float ).transpose() 
	wn  = np.array( dat[0] ) 
	gdat= np.array( dat[1:] )
	return wn, gdat
	
def matrixDiagReturnReIm( gdat , nbasis)  :
	datReDiag = range(nbasis)
	datImDiag = range(nbasis)
	for j in range(nbasis) :
		datReDiag[j] = gdat[2*j   +2*nbasis*j]
		datImDiag[j] = gdat[2*j+1 +2*nbasis*j]
	datReDiag = np.array( datReDiag )
	datImDiag = np.array( datImDiag )
	return datReDiag, datImDiag
def matrixOffDiagReturnReIm( gdat , nbasis)  :
	datReOffDiag = [ range(nbasis) for j in range(nbasis) ]
	datImOffDiag = [ range(nbasis) for j in range(nbasis) ]
	for i in range(nbasis) :
		for j in range(i+1,nbasis) :
			datReOffDiag[i][j] = np.array( gdat[2*j   +2*nbasis*i] )
			datImOffDiag[i][j] = np.array( gdat[2*j+1 +2*nbasis*i] )
	return datReOffDiag, datImOffDiag
def matrixReturnReIm( gdat , nbasis)  :
	datRe = [ range(nbasis) for j in range(nbasis) ]
	datIm = [ range(nbasis) for j in range(nbasis) ]
	for i in range(nbasis) :
		for j in range(nbasis) :
			datRe[i][j] = np.array( gdat[2*j   +2*nbasis*i] )
			datIm[i][j] = np.array( gdat[2*j+1 +2*nbasis*i] )
	return datRe, datIm
def matrixReturnDiagTransffromj( gdat , nbasis, Tjsth)  : # Normal(T).V.Dagger(T)
	datReDiag, datImDiag  = matrixReturnDiagTransfjnat( gdat , nbasis, np.conjugate(Tjsth).transpose() )
	return datReDiag, datImDiag
def matrixReturnDiagTransfjnat( gdat , nbasis, Tjnat)  :
	datRe   = [ range(nbasis) for j in range(nbasis) ]
	datIm   = [ range(nbasis) for j in range(nbasis) ]
	datComp = [ range(nbasis) for j in range(nbasis) ]
	datCompNat = [ range(nbasis) for j in range(nbasis) ]
	for i in range(nbasis) :
		for j in range(nbasis) :
			datRe[i][j]   = np.array( gdat[2*j   +2*nbasis*i] )
			datIm[i][j]   = np.array( gdat[2*j+1 +2*nbasis*i] )
			datComp[i][j] = datRe[i][j] + 1j * datIm[i][j]
	datCompNat = transfMwnconjtr( datComp , Tjnat )
	for i in range(nbasis) :
		for j in range(nbasis) :
			datRe[i][j]   = datCompNat[i][j].real
			datIm[i][j]   = datCompNat[i][j].imag
	print "DAT COMP : "
	for i in range(nbasis) :
		for j in range(nbasis) :
			print "{:19.6f}\t".format( datComp[i][j][0] ),
		print ""
	print "DAT COMP NAT : "
	for i in range(nbasis) :
		for j in range(nbasis) :
			print "{:19.6f}\t".format( datCompNat[i][j][0] ),
		print ""
	datReDiag = range(nbasis)
	datImDiag = range(nbasis)
	for i in range(nbasis) :
		datReDiag[i]   = datCompNat[i][i].real
		datImDiag[i]   = datCompNat[i][i].imag
	return datReDiag, datImDiag
def matrixDiagonality( gdat , nbasis ) :
	datRe   = [ range(nbasis) for j in range(nbasis) ]
	datIm   = [ range(nbasis) for j in range(nbasis) ]
	datComp = [ range(nbasis) for j in range(nbasis) ]
	for i in range(nbasis) :
		for j in range(nbasis) :
			datRe[i][j]   = np.array( gdat[2*j   +2*nbasis*i] )
			datIm[i][j]   = np.array( gdat[2*j+1 +2*nbasis*i] )
			datComp[i][j] = datRe[i][j] + 1j * datIm[i][j]
	datDiagonality = range(nbasis)
	for i in range(nbasis) :
		dum = 0.*1j
		for j in range(nbasis) :
			dum += np.array(datComp[i][j]) 
		dum -= np.array(datComp[i][i]) 
		datDiagonality[i] = [dum.real/np.array(datComp[i][i]).real, dum.imag/np.array(datComp[i][i]).imag ]
	return datDiagonality
def matrixDiagonalityTransffromj( gdat , nbasis, Tjsth)  :
	datRe   = [ range(nbasis) for j in range(nbasis) ]
	datIm   = [ range(nbasis) for j in range(nbasis) ]
	datComp = [ range(nbasis) for j in range(nbasis) ]
	datCompNat = [ range(nbasis) for j in range(nbasis) ]
	for i in range(nbasis) :
		for j in range(nbasis) :
			datRe[i][j]   = np.array( gdat[2*j   +2*nbasis*i] )
			datIm[i][j]   = np.array( gdat[2*j+1 +2*nbasis*i] )
			datComp[i][j] = datRe[i][j] + 1j * datIm[i][j]
	datCompNat = transfMwnconjtr( datComp , np.conjugate(Tjsth).transpose() )
	print "DAT COMP fromj : "
	for i in range(nbasis) :
		for j in range(nbasis) :
			print "{:19.6f}\t".format( datCompNat[i][j][0] ),
		print ""
	datDiagonality = range(nbasis)
	for i in range(nbasis) :
		dum = 0.*1j
		for j in range(nbasis) :
			dum += np.abs(np.array(datCompNat[i][j]).real) + np.abs(np.array(datCompNat[j][i]).real) + 1j*np.abs(np.array(datCompNat[i][j]).imag) + 1j*np.abs(np.array(datCompNat[j][i]).imag)
		dum -= np.abs(np.array(datCompNat[i][i]).real) + np.abs(np.array(datCompNat[i][i]).real) + 1j*np.abs(np.array(datCompNat[i][i]).imag) + 1j*np.abs(np.array(datCompNat[i][i]).imag)
		datDiagonality[i] = [dum.real/np.array(datCompNat[i][i]).real, dum.imag/np.array(datCompNat[i][i]).imag ]
	#print "DatDiagRe : ", [ datDiagonality[j][0][0] for j in range(nbasis) ]
	#print "DatDiagIm : ", [ datDiagonality[j][1][0] for j in range(nbasis) ]
	#print "DDiagR0: ",  datDiagonality[0][0]
	#print "DDiagI0: ",  datDiagonality[0][1]
	#print "DDiagR2: ",  datDiagonality[2][0]
	#print "DDiagI2: ",  datDiagonality[2][1]
	#print "DDiagR4: ",  datDiagonality[4][0]
	#print "DDiagI4: ",  datDiagonality[4][1]
	return datDiagonality

def splitFiltArg( filtArg , indlist , delimit="_" ) :
	filtArr = []
	k=-1
	rawArr = filtArg.split(delimit)  
	for i in range(0,len(rawArr),2) :
		arg = rawArr[i:i+2]
		k+=1
		#filtArr.append( [] ) 
		filtArr.append( arg ) 
	return filtArr

def splitFiltArgImp( filtArg , indlist , ni, delimit="_" ) :
	filtArr = []
	k=-1
	argarr = filtArg.split(delimit)[::2] 
	valarr = filtArg.split(delimit)[1::2] 
	lenarr = len(argarr)
	for iarg in range(lenarr) :
		arg = argarr[iarg]
		val = valarr[iarg]
		for iimp in range(ni) : 
			argimp = arg+delimit+"{}".format(iimp)
			for ind in indlist : 
				if ind == argimp  :
					k+=1
					filtArr.append( [] ) 
		filtArr[k] =  [argimp,val] 
	return filtArr

def simpledatacol( data ) :
	datarr = np.genfromtxt( data )
	return datarr.transpose()


def setattrwithlist( hobj, name, dname , real=False ) :
	try :
		fname = name ; dnamearr = dname
		f = hobj.tdir + "/result/u{:.3f}/{}.dat".format( hobj.UF , fname )
		rawdat = np.genfromtxt( f, dtype=float )
		setattr( hobj,  fname ,  rawdat  )
		ffshape = np.shape( rawdat ) 
		if len(ffshape) > 1 :	# multi-line data
			lastdat = rawdat[-1][1:] 
		else : 			# sinlgle-line data
			lastdat = rawdat[1:]
		for j in range(len(lastdat)/2) : 
			setattr( hobj,  fname+dnamearr[1+j] ,  lastdat[2*j]  )
			if real :
				if  abs(lastdat[2*j+1]) > 1e-5  : print "ERROR :: imaginary moment found on {}. (Stopped)".format(fname)  ; sys.exit(1)
	except :
		setattr( hobj,  fname ,  '.'  )
	try : 
		degen_0 = int(hobj.degen_0)
	except : 
		degen_0 = 1
	if degen_0 > 1 :
		for idegen in range(degen_0) :
			try :
				fnamedeg = fname+"%d"%idegen
				setattr( hobj,  fname ,  rawdat  )
				ffshape = np.shape( rawdat ) 
				if len(ffshape) > 1 :	# multi-line data
					lastdat = rawdat[-1-idegen][1:]
				else : 			# sinlgle-line data
					print "ERROR :: Importing static data of degenerate states : incompatible file-format with degeneracy. "
					sys.exit(1)
				for j in range(len(lastdat)/2) : 
					setattr( hobj,  fnamedeg+dnamearr[1+j] ,  lastdat[2*j]  )
					if  abs(lastdat[2*j+1]) > 1e-5  : print "ERROR :: imaginary moment found on {}. (Stopped)".format(fname)  ; sys.exit(1)
			except :
				setattr( hobj,  fnamedeg ,  '.'  )

def setattrwithlistdeg( hobj, fname, dnamearr , real=False , mustreal=False , iftrace=False , ifabs=False , diag=False ) :
	try :
		f = hobj.tdir + "/result/u{:.3f}/{}.dat".format( hobj.UF , fname )
		rawdat = np.genfromtxt( f, dtype=float )
		ffshape = np.shape( rawdat ) 
		if len(ffshape) > 1 :	# multi-line data
			lastdat = rawdat[-1][1:]
		else : 			# sinlgle-line data
			lastdat = rawdat[1:]
		if iftrace : 
			lenlastdat = len(lastdat)
			trace = 0.
			for idat in range(lenlastdat) :
				if ifabs : 
					trace = trace + np.abs(lastdat[idat]) + 1j*0
				else :
					trace = trace + lastdat[idat] + 1j*0
			lastdat = [trace.real , trace.imag]
		else : 
			setattr( hobj,  fname ,  rawdat  )
		if diag : 
			for j in range(len(dnamearr)-1) : 
				setattr( hobj,  fname+"diag"+dnamearr[1+j] ,  lastdat[2*j+2*2*j]  )
		else :
			for j in range(len(dnamearr)-1) : 
				if real :
					setattr( hobj,  fname+dnamearr[1+j] ,  lastdat[j]  )
				else :	
					setattr( hobj,  fname+dnamearr[1+j] ,  lastdat[2*j]  )
					if  mustreal and abs(lastdat[2*j+1]) > 1e-5  : print "ERROR :: imaginary moment found on {}. (Stopped)".format(fname)  ; sys.exit(1)
	except :
		setattr( hobj,  fname ,  '.'  )

def setattrTranscomponent( hobj, name, dname , real=False ) :
	try :
		fname = name ; dnamearr = dname
		dx = getattr( hobj,  fname+'x' )
		dy = getattr( hobj,  fname+'y' )
		setattr( hobj,  fname+'t' , np.sqrt(np.square(dx)+np.square(dy)) )
	except : 
		setattr( hobj,  fname ,  '.'  )
	try : 
		degen_0 = int(hobj.degen_0)
	except : 
		degen_0 = 1
	if degen_0 > 1 :
		for idegen in range(degen_0) :
			try :
				fnamedeg = fname+"%d"%idegen
				dx = getattr( hobj,  fnamedeg+'x' )
				dy = getattr( hobj,  fnamedeg+'y' )
				setattr( hobj,  fnamedeg+'t' , np.sqrt(np.square(dx)+np.square(dy))  )
			except : 
				setattr( hobj,  fnamedeg ,  '.'  )

def stdoutgrep( cmd ) :
	proc = subprocess.Popen( [ cmd ] , stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	return [ out, err ] 

def findgs( gptlenergyarr ) :
        ngptl = len(gptlenergyarr[0])
        e0arr = [999. for a in range(ngptl) ]
        for energyarr in gptlenergyarr  :
                for a in range(ngptl) :
                        try :
                                if e0arr[a] > float(energyarr[a]) :
                                        e0arr[a] = float(energyarr[a])
                                        sectorE0 = a
                        except :
                                pass
        return np.array(e0arr, dtype=float), sectorE0

class ldos  :
	def __init__(self, args , fname=None , basis="j", vertdata=False , vertSelf=False, vertSelfRealpart=False, tagEpBr="" ) :
		ndir = len(args.ddir)
		print "ndir : ", ndir
		dirarr = args.ddir
		i=0 
		self.basis = basis
		hobjarr = []
		for jj in range(ndir) :
			print "dir%d:"%i, dirarr[jj] ; i+=1
			hobj = headobj( args.ddir[jj] )
			hobj.readParameters()
			hobjarr.append( hobj ) 
			dat = args.ddir[jj]
		
			wn   	= [""]*  hobj.Nplot
			coldata = [""]*( hobj.Nplot )
			if args.nplot : hobj.Nplot = int(args.nplot) 
			if args.nplotdos : hobj.Nplot = int(args.nplotdos) 
			#hobj.title = ""
			for i in [1] :
				if vertdata : 
					fname = "ldos"
					if vertSelf or vertSelfRealpart : fname = "selfreal"  
					datf= "{}/lattice/u{:.3f}/{}_{}{}th.dat".format(dat, hobj.UF, fname, tagEpBr, hobj.count )
				elif fname.find("selfreal")>-1 : 
					datf= "{}/lattice/u{:.3f}/{}_Nplot{}_{}{}th.dat".format(dat, hobj.UF, fname, hobj.Nplot, tagEpBr, hobj.count )
				else : 
					datf= "{}/lattice/{}_u{:.2f}_Nplot{}.dat".format(dat, fname, hobj.UF, hobj.Nplot )
				print i, " :", dat
				print i, "f:", datf
				datarr = np.genfromtxt( datf, dtype=float ) 
				wn      = datarr.transpose()[0]
				if vertdata : 
					coldata = datarr.transpose()[4:]
				else : 
					coldata = datarr.transpose()[1:]
				if vertSelf or vertSelfRealpart : 
					istart = 2
					if vertSelfRealpart : istart=1
					diagarr = range( istart, 73 , 14 ) 
					print "vertSelf datarr : ", diagarr
					coldata = datarr.transpose()[ diagarr ]
				coldata = coldata.transpose()
		
			wn   = np.array( wn   )
			print len(wn ), "points (wn)"
			coldata = np.array( coldata )
			coldata = coldata.transpose()
			print len(coldata), "components (basis vectors + sum )"
			self.wn = wn
			self.coldata = coldata
			self.datlab = r"PDOS"
			if fname.find("self")>-1 : self.datlab = "Im$\Sigma$"
			if vertSelfRealpart : self.datlab = "Re$\Sigma$"
	def totally( self, args, ffig, aax , xdat=None, ydat=None, hobj=None, chdat=None, laboverw=False , dashes=(None,None) ) :
		jj=0 
		coldata	= self.coldata
		wn	= self.wn
		basis	= self.basis
		for j in chdat :
			ldos = coldata[j]
			ldos = np.array( ldos )
			lc = [ "-" ]*2 + [ "-" ] *2 + [ ":" ] *2 
			if laboverw :	labelstr = laboverw
			else :		labelstr = bolabel(basis,j)
			if args.trans : 
				aax.plot( ldos, wn, blt(basis,j) , label=labelstr ) #, dashes=dashes )
				aax.set_xlabel( self.datlab )#, fontsize=fs )
				aax.set_ylabel(r'$E$ (eV)')#, fontsize=fs )
				if j==1 :
					if args.fillb : aax.fill( ldos, wn, 'C%d'%dd+blt(basis,j) ,alpha=0.3 ) 
				if j==0 :
					if args.fillb : aax.fill_betweenx( wn, 0, ldos, facecolor='C%d'%j ,alpha=0.3 ) 
				print "if LEG : ", basis, j, bolabel(basis,j), blt(basis,j)
			else : 
				aax.plot( wn, ldos, blt(basis,j) , label=labelstr ) #, dashes=dashes )
				aax.set_ylabel( self.datlab )#, fontsize=fs )
				aax.set_xlabel(r'$\omega$')#, fontsize=fs )
				print "else LEG : ", basis, j, bolabel(basis,j), blt(basis,j)
			#aax.axvline( 0 , color='gray', linestyle="--", lw=0.5 ) 
			aax.axhline( 0 , color='gray', linestyle="--", lw=0.5 ) 
			#aax.legend()# fontsize=fsm*0.8 ) 
			jj+=1
		#plt.title("U={} J={} S={:.2f} D={} ; nb={}".format(UF,J,S,D,Nb))
		aax.tick_params( axis="both", direction="in" )#, labelsize=fs , length=10, width=2 ) 
		aax.tick_params( axis="x", bottom=True, labelbottom=True )
		if args.notick :
			aax.tick_params( axis="both", left=True,   labelbottom=False, labelleft=False  )
		elif args.yoff :
			aax.tick_params( axis="both", left=True,   labelleft=False  )
		if args.trans : 
			pass
			#aax.set_xlim( 0, 1.2 ) 
		else :
			aax.set_ylim( -0.05,1.05 ) 
		if args.xxrange : 
			axSetRange( aax , args.xxrange, 'x' ) 
		if args.yrange : 
			axSetRange( aax , args.yrange,  'y' ) 
		if args.dxtick :
			aax.set_xticks(aax.get_xticks()[1::2])
		if args.dytick :
			aax.set_yticks(aax.get_yticks()[1::2])

def cptmpftn( fname, cmd="cp", destMac=None, dtype="png" , affix=None, ip="10.10.23.25" ) :
	if destMac is True :
		destFolder = "Realitystone@{}:~/Dropbox/tmp_linux/".format(ip)
	elif (destMac is None) or (destMac is False)  :
		destFolder = "/home/jun/Dropbox/tmp_linux/" 
	print "Copying	: ",	fname+"."+dtype
	print "into	: ",	destFolder
	cmdline = cmd+" "+fname+"."+dtype+" "+destFolder
	os.system( cmdline )
	cmdline = 'mv'+" "+fname+"."+dtype+" "+'/home/jun/recycle/'
	os.system( cmdline )

def fmt10(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	return r'${} \times 10^{{{}}}$'.format(a, b)
