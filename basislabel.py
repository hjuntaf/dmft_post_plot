mulabel={
	0 : "\mu=0",
	1 : "\mu=1",
	2 : "\mu=2",
	3 : "\mu=3",
	4 : "\mu=4",
	5 : "\mu=5" }
jlabel={
	0 : "3/2,  3/2",
	1 : "3/2,  1/2",
	2 : "3/2, -1/2",
	3 : "3/2, -3/2",
	4 : "1/2,  1/2",
	5 : "1/2, -1/2" }
slabel={
	0 : "xy, up",
	1 : "xy, dn",
	2 : "yz, up",
	3 : "yz, dn",
	4 : "zx, up",
	5 : "zx, dn" }
splabel={
	0 : "xy, up",
	1 : "xy, dn",
	2 : "+1, up",
	3 : "+1, dn",
	4 : "-1, up",
	5 : "-1, dn" }
plabel={
	0 : "xy",
	1 : "xy",
	2 : "yz",
	3 : "yz",
	4 : "zx",
	5 : "zx" }
plabel={
	0 : "$d_{xy}$",
	1 : "$d_{xy}$",
	2 : "$d_{yz}$",
	3 : "$d_{yz}$",
	4 : "$d_{zx}$",
	5 : "$d_{zx}$" }
plabelspin = slabel
plabelcompactdorb={
	0 : "$d_{xy}$",
	1 : "$d_{xy}$",
	2 : "$d_{yz/zx}$",
	3 : "$d_{yz/zx}$",
	4 : "$d_{yz/zx}$",
	5 : "$d_{yz/zx}$" }
plabelcompact={
	0 : "xy",
	1 : "xy",
	2 : "yz/zx",
	3 : "yz/zx",
	4 : "yz/zx",
	5 : "yz/zx" }
plabelcompactspin={
	0 : "xy, up",
	1 : "xy, dn",
	2 : "yz/zx, up",
	3 : "yz/zx, dn",
	4 : "yz/zx, up",
	5 : "yz/zx, dn" }
pjlabel={
	0 : "(3/2, 3/2)",
	1 : "(3/2, 1/2)",
	2 : "(3/2,-1/2)",
	3 : "(3/2,-3/2)",
	4 : "(1/2, 1/2)",
	5 : "(1/2,-1/2)" }
def nlabel(b,j) :
	return "{}".format(j)
def offdlabel(b,j) :
	return "{},{}".format(j[0],j[1])
def blabel(b,j, nbasis=6) :
	lab =''
	iimp = j/nbasis
	j = j%nbasis
	if( b.find("t2g")>-1 ) :
		lab	= slabel[j]
	elif( b.find("t")>-1 ) :
		lab	= slabel[j]
	elif( b.find("sbasis")>-1 ) :
		lab	= slabel[j]
	elif( b.find("jbasis")>-1 ) :
		lab	= jlabel[j]
	elif( b.find("j")>-1 ) :
		lab	= jlabel[j]
	elif( b.find("p")>-1 ) :
		lab	= plabel[j]
	else :
		lab	= mulabel[j]
	print "BOLABEL : ",  lab
	if iimp>1 : 
		return "atom%d"%iimp + lab
	else :
		return lab
def bolabel(b,j, nbasis=6) :
	iimp = j/nbasis
	if iimp<1 : 
		labimp = ''
	else : 
		labimp = '{}'.format(iimp)+'^{\mathrm{th}};'
		labimp = ''

	if j<0 : 
		if labimp==''  : 
			if iimp < 0 : 
				iimp = (-j)/nbasis
			labimp = '{}'.format(iimp)+'^{\mathrm{th}};'
		return labimp+'\mathrm{total}'

	j = j%nbasis
	if( b.find("total")>-1 ) :
		lab = "total"
	elif( b.find("t2g")>-1 and b.find("compact")>-1 ) :
		if( b.find("spin")>-1 ) :
			lab = plabelcompactspin[j]
		elif( b.find("dorb")>-1 ) :
			lab = plabelcompactdorb[j]
		else  : 
			lab = plabelcompact[j]
	elif( b.find("t2g")>-1 ) :
		if( b.find("spin")>-1 ) :
			lab = plabelspin[j]
		else : 
			lab = plabel[j]
	elif( b.find("t")>-1 ) :
		if( b.find("spin")>-1 ) :
			lab = plabelspin[j]
		else : 
			lab = plabel[j]
	elif( b.find("s")>-1 ) :
		lab = splabel[j]
	elif( b.find("sbasis")>-1 ) :
		lab = slabel[j]
	elif( b.find("jbasis")>-1 ) :
		lab = pjlabel[j]
	elif( b.find("j")>-1 ) :
		lab = pjlabel[j]
	elif( b.find("p")>-1 ) :
		lab = plabel[j]
	else :
		lab = mulabel[j]
	lab = labimp + lab
	#print "label ", iimp, j, lab
	return lab
def bolabelArr(b) :
	if( b.find("total")>-1 ) :
		return "total"
	elif( b.find("t2g")>-1 and b.find("compact")>-1 ) :
		if( b.find("spin")>-1 ) :
			return plabelcompactspin
		elif( b.find("dorb")>-1 ) :
			return plabelcompactdorb
		else  : 
			return plabelcompact
	elif( b.find("t2g")>-1 ) :
		if( b.find("spin")>-1 ) :
			return plabelspin
		else : 
			return plabel
	elif( b.find("t")>-1 ) :
		if( b.find("spin")>-1 ) :
			return plabelspin
		else : 
			return plabel
	elif( b.find("s")>-1 ) :
		return splabel
	elif( b.find("sbasis")>-1 ) :
		return slabel
	elif( b.find("jbasis")>-1 ) :
		return pjlabel
	elif( b.find("j")>-1 ) :
		return pjlabel
	elif( b.find("p")>-1 ) :
		return plabel
	else :
		return mulabel
def updnplotlabel(j, order) :
	if order.find("j")>-1 :
		if   j==0 :  return [0,3,'$_{j=3/2}$']
		elif j==1 :  return [1,2,'$_{j=3/2}$']
		elif j==2 :  return [4,5,'$_{j=1/2}$']
		else : print "ERROR in updnplot.0 range"
	else :
		if   j==0 :  return [0,1,'$_{xy}$']
		elif j==1 :  return [2,3,'$_{yz}$']
		elif j==2 :  return [4,5,'$_{zx}$']
		else : print "ERROR in updnplot.1 range"
ltt2g={
	0 : "-" ,
	1 : "-" ,
	2 : "-" ,
	3 : "-" ,
	4 : ":" ,
	5 : ":" }
ltjeff={
	0 : "-" ,
	1 : "-" ,
	2 : "-" ,
	3 : "-" ,
	4 : ":" ,
	5 : ":" }
ltt2gdash={
	0 : "-" ,
	1 : "-" ,
	2 : "--" ,
	3 : "--" ,
	4 : ":" ,
	5 : ":" }
ltjeffdash={
	0 : "-" ,
	1 : "--" ,
	2 : "--" ,
	3 : "-" ,
	4 : ":" ,
	5 : ":" }
allsolid={
	0 : "-" ,
	1 : "-" ,
	2 : "-" ,
	3 : "-" ,
	4 : "-" ,
	5 : "-" ,
	6 : "-" ,
	7 : "-" }
alldiff={
	0 : ":" ,
	1 : "-" ,
	2 : "--" ,
	3 : "-." ,
	4 : ":" ,
	5 : "--" ,
	6 : "-" ,
	7 : "-." }
def dashdiff(j) :
	if j<1 : 
		return (None,None)
	else : 
		return [j+1,1]
def blt(b,j, nbasis=6) :
	j = j%nbasis
	if( ( b.find("alldiff")>-1 ) ) :
		return alldiff[j]
	elif( ( b.find("solid")>-1 )or( b.find("line")>-1 ) ) :
		return allsolid[j]
	elif( b.find("t2gdash")>-1 ) :
		return ltt2gdash[j]
	elif( b.find("t2g")>-1 ) :
		return ltt2g[j]
	elif( b.find("t")>-1 ) :
		return ltt2g[j]
	elif( b.find("sbasis")>-1 ) :
		return ltt2g[j]
	elif( b.find("jdash")>-1 ) :
		return ltjeffdash[j]
	elif( b.find("jbasis")>-1 ) :
		return ltjeff[j]
	elif( b.find("j")>-1 ) :
		return ltjeff[j]
	elif( b.find("p")>-1 ) :
		return ltjeff[j]
	else :
		return ltt2g[j]

def orbfilt( orb, nbasis=6) :
	if( (orb == 'q3') or (orb=='q3u') ) :
		return '3/2, 3/2'
	if( (orb == 'q1') or (orb=='q1u') ) :
		return '3/2, 1/2'
	if( (orb == 'd1') or (orb=='d1u') ) :
		return '1/2, 1/2'
	else :
		return orb

msindj= [ "x","D","D","x","s","s" ]
markersjeff	= msindj
markersj		= msindj
msindt= [ "o","o",">",">","^","^" ] 
markerst2g 	= msindt
markerst	= msindt
msindnone = [ '-', '--', ':', '-.']*4
markersnone	= msindnone
markersn	= msindnone
markersdiff = ["o", "v", "^", "<", ">", "s" , "D", "P", "X" ]*2
markersd	= markersdiff

markerssub = ["({})".format(chr(a)) for a in range(97,123) ] + ["(1{})".format(chr(a)) for a in range(97,123) ]
markerssubfig	= markerssub

def bmarkers(b,j) :
	if( b.find("t2g")>-1 ) :	return msindt[j]
	elif( b.find("t")>-1 ) :	return msindt[j]
	elif( b.find("s")>-1 ) :	return msindt[j]
	elif( b.find("sbasis")>-1 ) :	return msindt[j]
	elif( b.find("jbasis")>-1 ) :	return msindj[j]
	elif( b.find("j")>-1 ) :	return msindj[j]
	else :				return markersdiff[j]

def itemlab(ind) : 
	if ind.find('fillimpq3u')>-1  : return "(3/2, 3/2)"
	elif ind.find('fillimpq1u')>-1  : return "(3/2, 1/2)"
	elif ind.find('fillimpq1d')>-1  : return "(3/2,-1/2)"
	elif ind.find('fillimpq3d')>-1  : return "(3/2,-3/2)"
	elif ind.find('fillimpd1u')>-1  : return "(1/2, 1/2)"
	elif ind.find('fillimpd1d')>-1  : return "(1/2,-1/2)"
	elif ind.find('fillimpxyu')>-1  : return "xy, up"
	elif ind.find('fillimpxyd')>-1  : return "xy, dn"
	elif ind.find('fillimpyzu')>-1  : return "yz, up"
	elif ind.find('fillimpyzd')>-1  : return "yz, dn"
	elif ind.find('fillimpzxu')>-1  : return "zx, up"
	elif ind.find('fillimpzxd')>-1  : return "zx, dn"
	elif ind.find('fillimptot')>-1 	: return "average"
	elif ind.find('JJdegsq')>-1 	: return r"$JJ$"
	elif ind.find('JJdegzz')>-1 	: return r"$J_zJ_z$"
	elif ind.find('JJdegxyxy')>-1 	: return r"$J_xJ_x+J_yJ_y$"
	elif ind.find('JaJadegsq')>-1 	: return r"$KK$"
	elif ind.find('JaJadegzz')>-1 	: return r"$K_zK_z$"
	elif ind.find('JaJadegxyxy')>-1 : return r"$K_xK_x+K_yK_y$"
	elif ind.find('LLdegsq')>-1 	: return r"$LL$"
	elif ind.find('LLdegzz')>-1 	: return r"$L_zL_z$"
	elif ind.find('LLdegxyxy')>-1 	: return r"$L_xL_x+L_yL_y$"
	elif ind.find('SSdegsq')>-1 	: return r"$SS$"
	elif ind.find('SSdegzz')>-1 	: return r"$S_zS_z$"
	elif ind.find('SSdegxyxy')>-1 	: return r"$S_xS_x+S_yS_y$"
	elif ind.find('MMdegsq')>-1 	: return r"$MM$"
	elif ind.find('MMdegzz')>-1 	: return r"$M_zM_z$"
	elif ind.find('MMdegxyxy')>-1 	: return r"$M_xM_x+M_yM_y$"
	elif ind.find('Jzxydegz')>-1 	: return r"$J_z$"
	elif ind.find('Jzxydegx')>-1 	: return r"$J_x$"
	elif ind.find('Jzxydegy')>-1 	: return r"$J_y$"
	elif ind.find('Jazxydegz')>-1 	: return r"$K_z$"
	elif ind.find('Jazxydegx')>-1 	: return r"$K_x$"
	elif ind.find('Jazxydegy')>-1 	: return r"$K_y$"
	elif ind.find('Lzxydegz')>-1 	: return r"$L_z$"
	elif ind.find('Lzxydegx')>-1 	: return r"$L_x$"
	elif ind.find('Lzxydegy')>-1 	: return r"$L_y$"
	elif ind.find('Szxydegz')>-1 	: return r"$S_z$"
	elif ind.find('Szxydegx')>-1 	: return r"$S_x$"
	elif ind.find('Szxydegy')>-1 	: return r"$S_y$"
	elif ind.find('Mzxydegz')>-1 	: return r"$M_z$"
	elif ind.find('Mzxydegx')>-1 	: return r"$M_x$"
	elif ind.find('Mzxydegy')>-1 	: return r"$M_y$"
	elif ind.find( 'Jzxydegt')>-1 	: return r"$J_\perp=J_x+J_y$"
	elif ind.find('Jazxydegt')>-1 	: return r"$K_\perp=K_x+K_y$"
	elif ind.find( 'Lzxydegt')>-1 	: return r"$L_\perp=L_x+L_y$"
	elif ind.find( 'Szxydegt')>-1 	: return r"$S_\perp=S_x+S_y$"
	elif ind.find( 'Mzxydegt')>-1 	: return r"$M_\perp=M_x+M_y$"
	elif ind.find('Jzdeg')>-1 	: return r"$\hat{J}_z$"
	elif ind.find('Jxdeg')>-1 	: return r"$\hat{J}_x$"
	elif ind.find('Jydeg')>-1 	: return r"$\hat{J}_y$"
	elif ind.find('Lzdeg')>-1 	: return r"$\hat{L}_z$"
	elif ind.find('Lxdeg')>-1 	: return r"$\hat{L}_x$"
	elif ind.find('Lydeg')>-1 	: return r"$\hat{L}_y$"
	elif ind.find('Szdeg')>-1 	: return r"$\hat{S}_z$"
	elif ind.find('Sxdeg')>-1 	: return r"$\hat{S}_x$"
	elif ind.find('Sydeg')>-1 	: return r"$\hat{S}_y$"
	elif ind.find('Mzdeg')>-1 	: return r"$\hat{M}_z$"
	elif ind.find('Mxdeg')>-1 	: return r"$\hat{M}_x$"
	elif ind.find('Mydeg')>-1 	: return r"$\hat{M}_y$"
	elif ind.find('JJxydegMatevalhalf')>-1 	: return r"$(\hat{J}_x\hat{J}_x+\hat{J}_y\hat{J}_y)/2$"
	elif ind.find('LLxydegMatevalhalf')>-1 	: return r"$(\hat{L}_x\hat{L}_x+\hat{L}_y\hat{L}_y)/2$"
	elif ind.find('SSxydegMatevalhalf')>-1 	: return r"$(\hat{S}_x\hat{S}_x+\hat{S}_y\hat{S}_y)/2$"
	elif ind.find('MMxydegMatevalhalf')>-1 	: return r"$(\hat{M}_x\hat{M}_x+\hat{M}_y\hat{M}_y)/2$"
	elif ind.find('JJdeg')>-1 	: return r"$\hat{J}\hat{J}$"
	elif ind.find('LLdeg')>-1 	: return r"$\hat{L}\hat{L}$"
	elif ind.find('SSdeg')>-1 	: return r"$\hat{S}\hat{S}$"
	elif ind.find('MMdeg')>-1 	: return r"$\hat{M}\hat{M}$"
	elif ind.find('JJzzdeg')>-1 	: return r"$\hat{J}_z\hat{J}_z$"
	elif ind.find('LLzzdeg')>-1 	: return r"$\hat{L}_z\hat{L}_z$"
	elif ind.find('SSzzdeg')>-1 	: return r"$\hat{S}_z\hat{S}_z$"
	elif ind.find('MMzzdeg')>-1 	: return r"$\hat{M}_z\hat{M}_z$"
	elif ind.find('JJxydeg')>-1 	: return r"$\hat{J}_x\hat{J}_x+\hat{J}_y\hat{J}_y$"
	elif ind.find('LLxydeg')>-1 	: return r"$\hat{L}_x\hat{L}_x+\hat{L}_y\hat{L}_y$"
	elif ind.find('SSxydeg')>-1 	: return r"$\hat{S}_x\hat{S}_x+\hat{S}_y\hat{S}_y$"
	elif ind.find('MMxydeg')>-1 	: return r"$\hat{M}_x\hat{M}_x+\hat{M}_y\hat{M}_y$"
	elif ind.find('x2y2')>-1 	: return r"$\hat{Q}$" + r"$_{x2-y2}$"	 + "({})".format(ind[1])  
	elif ind.find('z2r2')>-1 	: return r"$\hat{Q}$" + r"$_{z2-r2}$"	 + "({})".format(ind[1]) 
	elif ind.find('xy')>-1 		: return r"$\hat{Q}$" + r"$_{xy}$"	 + "({})".format(ind[1]) 
	elif ind.find('yz')>-1 		: return r"$\hat{Q}$" + r"$_{yz}$"	 + "({})".format(ind[1]) 
	elif ind.find('zx')>-1 		: return r"$\hat{Q}$" + r"$_{zx}$"	 + "({})".format(ind[1]) 
	else : return ind

def itemmarker(ind) :
	if ind.find('fillimpq3u')>-1 : return markersjeff[0]
	elif ind.find('fillimpq1u')>-1 : return markersjeff[1]
	elif ind.find('fillimpq1d')>-1 : return markersjeff[2]
	elif ind.find('fillimpq3d')>-1 : return markersjeff[3]
	elif ind.find('fillimpd1u')>-1 : return markersjeff[4]
	elif ind.find('fillimpd1d')>-1 : return markersjeff[0]
	elif ind.find('fillimpxyu')>-1 : return markerst2g[0] 
	elif ind.find('fillimpxyd')>-1 : return markerst2g[1] 
	elif ind.find('fillimpyzu')>-1 : return markerst2g[2] 
	elif ind.find('fillimpyzd')>-1 : return markerst2g[3] 
	elif ind.find('fillimpzxu')>-1 : return markerst2g[4] 
	elif ind.find('fillimpzxd')>-1 : return markerst2g[5] 
	elif ind.find('fillimptot')>-1 : return 'o'
	elif ind.find('JJdegsq')>-1 	: return 'o'
	elif ind.find('JJdegzz')>-1 	: return 'o'
	elif ind.find('JJdegxyxy')>-1 	: return 'o'
	elif ind.find('JaJadegsq')>-1 	: return 'P'
	elif ind.find('JaJadegzz')>-1 	: return 'P'
	elif ind.find('JaJadegxyxy')>-1 : return 'P'
	elif ind.find('LLdegsq')>-1 	: return '^'
	elif ind.find('LLdegzz')>-1 	: return '^'
	elif ind.find('LLdegxyxy')>-1 	: return '^'
	elif ind.find('SSdegsq')>-1 	: return 's'
	elif ind.find('SSdegzz')>-1 	: return 's'
	elif ind.find('SSdegxyxy')>-1 	: return 's'
	elif ind.find('MMdegsq')>-1 	: return '*'
	elif ind.find('MMdegzz')>-1 	: return '*'
	elif ind.find('MMdegxyxy')>-1 	: return '*'
	elif ind.find('Jzxydegz')>-1 	: return 'o'
	elif ind.find('Jzxydegx')>-1 	: return 'o'
	elif ind.find('Jzxydegy')>-1 	: return 'o'
	elif ind.find('Jazxydegz')>-1 	: return 'P'
	elif ind.find('Jazxydegx')>-1 	: return 'P'
	elif ind.find('Jazxydegy')>-1 	: return 'P'
	elif ind.find('Lzxydegz')>-1 	: return '^'
	elif ind.find('Lzxydegx')>-1 	: return '^'
	elif ind.find('Lzxydegy')>-1 	: return '^'
	elif ind.find('Szxydegz')>-1 	: return 's'
	elif ind.find('Szxydegx')>-1 	: return 's'
	elif ind.find('Szxydegy')>-1 	: return 's'
	elif ind.find('Mzxydegz')>-1 	: return '*'
	elif ind.find('Mzxydegx')>-1 	: return '*'
	elif ind.find('Mzxydegy')>-1 	: return '*'
	elif ind.find( 'Jzxydegt')>-1 	: return 'o'
	elif ind.find('Jazxydegt')>-1 	: return 'P'
	elif ind.find( 'Lzxydegt')>-1 	: return '^'
	elif ind.find( 'Szxydegt')>-1 	: return 's'
	elif ind.find( 'Mzxydegt')>-1 	: return '*'
	#elif ind.find( 'JJdegMateval')>-1 	: return 'o'
	#elif ind.find( 'LLdegMateval')>-1 	: return '^'
	#elif ind.find( 'SSdegMateval')>-1 	: return 's'
	#elif ind.find( 'MMdegMateval')>-1 	: return '*'
	#elif ind.find( 'JJzzdegMateval')>-1 	: return 'd'
	#elif ind.find( 'LLzzdegMateval')>-1 	: return 'd'
	#elif ind.find( 'SSzzdegMateval')>-1 	: return 'd'
	#elif ind.find( 'MMzzdegMateval')>-1 	: return 'd'
	#elif ind.find( 'JJxydegMateval')>-1 	: return '*'
	#elif ind.find( 'LLxydegMateval')>-1 	: return '*'
	#elif ind.find( 'SSxydegMateval')>-1 	: return '*'
	#elif ind.find( 'MMxydegMateval')>-1 	: return '*'
	#elif ind.find( 'Jzdeg')>-1 	: return 'o'
	#elif ind.find( 'Lzdeg')>-1 	: return '^'
	#elif ind.find( 'Szdeg')>-1 	: return 's'
	#elif ind.find( 'Mzdeg')>-1 	: return '*'
	#elif ind.find( 'Jxdeg')>-1 	: return 'o'
	#elif ind.find( 'Lxdeg')>-1 	: return '^'
	#elif ind.find( 'Sxdeg')>-1 	: return 's'
	#elif ind.find( 'Mxdeg')>-1 	: return '*'
	#elif ind.find( 'Jydeg')>-1 	: return 'o'
	#elif ind.find( 'Lydeg')>-1 	: return '^'
	#elif ind.find( 'Sydeg')>-1 	: return 's'
	#elif ind.find( 'Mydeg')>-1 	: return '*'
	elif ind.find( 'Jzdeg')>-1 	: return 'd'
	elif ind.find( 'Lzdeg')>-1 	: return 'd'
	elif ind.find( 'Szdeg')>-1 	: return 'd'
	elif ind.find( 'Mzdeg')>-1 	: return 'd'
	elif ind.find( 'Jxdeg')>-1 	: return '*'
	elif ind.find( 'Lxdeg')>-1 	: return '*'
	elif ind.find( 'Sxdeg')>-1 	: return '*'
	elif ind.find( 'Mxdeg')>-1 	: return '*'
	elif ind.find( 'Jydeg')>-1 	: return 'o'
	elif ind.find( 'Lydeg')>-1 	: return '^'
	elif ind.find( 'Sydeg')>-1 	: return 's'
	elif ind.find( 'Mydeg')>-1 	: return '*'
	elif ind.find( 'JJ')>-1 	: return 'o'
	elif ind.find( 'LL')>-1 	: return '^'
	elif ind.find( 'SS')>-1 	: return 's'
	elif ind.find( 'MM')>-1 	: return '*'
	elif ind.find( 'QJ')>-1 	: return 'o'
	elif ind.find( 'QL')>-1 	: return '^'
	elif ind.find( 'QS')>-1 	: return 's'
	elif ind.find( 'QM')>-1 	: return '*'
	else : return '.'
def colorftn(ind) :
	if   ind.find( 'Mz' )>-1 				: return 'C1'
	elif ind.find( 'zz' )>-1 				: return 'C1'
	elif ind.find( 'xy' )>-1 				: return 'C4'
	elif ind.find( 'MM' )>-1 				: return 'C1'
	elif ind.find( 'Mx' )>-1 				: return 'C1'
	elif ind.find( 'My' )>-1 				: return 'C1'
	elif ind.find( 'My' )>-1 				: return 'C1'
	elif ind.find( 'J' )>-1 				: return 'C2'
	elif ind.find( 'L' )>-1 				: return 'C3'
	elif ind.find( 'S' )>-1 				: return 'C4'
	else : return ''
def colordiffftn(ind) :
	return "C{}".format(ind)
def returnblankftn(ind) :
	return ""
def filterlabel( lab , paper=False, unit=False , unitlab='' ) :
	labconverted = lab
	if (lab == "nOcculatt") or (lab == "nOcculattAra") or (lab == "nOcculattAraFormer") :
		if paper : 
			labconverted = r"$n_{\mathrm{el}}$"
		else : 
			labconverted = r"$n_{latt}$"
		if lab.find("Imp")>-1 : 
			iimp = lab[lab.find("Imp")+3:]
			#labconverted = labconverted + r'$^{\mathrm{{{}}}}$'.format( 'site'+iimp )
			labconverted += r"$^{{({})}}$".format( 'site'+iimp )
	elif lab == "fillimp0xy" :
		labconverted = r"$n^{(xy)}_{imp}$"
		if paper : 
			labconverted = r"$n^{(xy)}_{\mathrm{el}}$"
	elif lab == "fillimp0yz" :
		labconverted = r"$n^{(yz)}_{imp}$"
		if paper : 
			labconverted = r"$n^{(yz)}_{\mathrm{el}}$"
	elif lab == "fillimp0zx" :
		labconverted = r"$n^{(zx)}_{imp}$"
		if paper : 
			labconverted = r"$n^{(zx)}_{\mathrm{el}}$"
	elif lab == "tfilling_0" :
		labconverted = r"$n_{imp}$"
		if paper : 
			labconverted = r"$n_{\mathrm{el}}$"
	elif lab =="rho" :
		labconverted = r"$\rho$"
	elif lab =="omega" :
		labconverted = r"$\omega$"
	elif lab =="D" :
		labconverted = r"$\Delta\mu_{\mathrm{chem}}$"
	elif lab =="U" :
		labconverted = r"$U$"
	elif lab =="J" or lab.find("J_H")>-1 : 
		labconverted = r"$J_H$"
	elif lab =="soc" or lab.find("lambdasoc")>-1 : 
		labconverted = r"$\lambda_{\mathrm{SOC}}$"
	elif (lab.find("powerplusconstfitt")>-1) or (lab.find("powerplusconstfitj")>-1) : 
		orblab = orbfilt( lab[-3:-1] )
		labconverted = r"$\tilde{{\alpha}}^{{({})}}_{{{}}}$".format(orblab,'')
	elif (lab.find("powerfitt")>-1) or (lab.find("powerfitj")>-1) : 
		orblab = orbfilt( lab[-3:-1] )
		labconverted = r"$\alpha^{{({})}}_{{{}}}$".format(orblab,'')
	elif lab.find("powerfitconst")>-1 : 
		orblab = orbfilt( lab[-3:-1] )
		labconverted = r"$C^{{({})}}_{{{}}}$".format(orblab,'')
	elif lab.find("powerplusconstfitconst")>-1 : 
		orblab = orbfilt( lab[-3:-1] )
		labconverted = r"$\tilde{{C}}^{{({})}}_{{{}}}$".format(orblab,'')
	elif lab =="powerfitconstyzu" : 
		labconverted = r"$C^{(yz)}$"
	elif lab =="powerfitconstxyu" : 
		labconverted = r"$C^{(xy)}$"
	elif lab =="powerfittyzu" : 
		labconverted = r"$\alpha^{(yz)}$"
	elif lab =="powerfittxyu" : 
		labconverted = r"$\alpha^{(xy)}$"
	elif lab =="powerfittq3u" : 
		labconverted = r"$\alpha^{(3/2,3/2)}$"
	elif lab =="powerfittq1u" : 
		labconverted = r"$\alpha^{(3/2,1/2)}$"
	elif lab =="powerfittq1u" : 
		labconverted = r"$\alpha^{(1/2,1/2)}$"
	elif lab.find("massfitimt")>-1 : 
		orblab = lab[-3:-1]
		labconverted = r"$m^{{*({})}}$".format(orblab)
	elif lab.find("masspoly4fit")>-1 : 
		orblab = orbfilt( lab[-3:-1] )
		labconverted = r"$m^{{*({})}}_{{{}}}$".format(orblab,'poly4')
		if paper : 
			labconverted = r"$m^{{*({})}}$".format(orblab)
	elif lab.find("masspoly3fit")>-1 : 
		orblab = orbfilt( lab[-3:-1] )
		labconverted = r"$m^{{*({})}}_{{{}}}$".format(orblab,'poly3')
	elif lab.find("masspoly2fit")>-1 : 
		orblab = orbfilt( lab[-3:-1] )
		labconverted = r"$m^{{*({})}}_{{{}}}$".format(orblab,'poly2')
	elif lab.find("massfitimj")>-1 : 
		orb = lab[-3:]
		if orb=='q3u' :
			orblab = "3/2,3/2"
		elif orb=='q1u' :
			orblab = "3/2,1/2"
		elif orb=='d1u' :
			orblab = "1/2,1/2"
		labconverted = r"$m^{{*({})}}$".format(orblab)
	elif lab.find("Zpoly4fit")>-1 : 
		orblab = orbfilt( lab[-3:-1] )
		labconverted = r"$Z^{{({})}}_{{{}}}$".format(orblab,'poly4')
		if paper : 
			labconverted = r"$Z^{{({})}}$".format(orblab)
	elif lab.find("Zpoly3fit")>-1 : 
		orblab = orbfilt( lab[-3:-1] )
		labconverted = r"$Z^{{({})}}_{{{}}}$".format(orblab,'poly3')
	elif lab.find("Zpoly2fit")>-1 : 
		orblab = orbfilt( lab[-3:-1] )
		labconverted = r"$Z^{{({})}}_{{{}}}$".format(orblab,'poly2')
	elif lab.find("Zfitimt")>-1 : 
		orblab = lab[-3:-1]
		labconverted = r"$Z^{{({})}}$".format(orblab)
	elif lab.find("Zfitimj")>-1 : 
		orb = lab[-3:]
		if orb=='q3u' :
			orblab = "3/2,3/2"
		elif orb=='q1u' :
			orblab = "3/2,1/2"
		elif orb=='d1u' :
			orblab = "1/2,1/2"
		labconverted = r"$Z^{{({})}}$".format(orblab)
	elif lab =="alpha" :
		labconverted = r"$\alpha$"
	elif lab.find("fermidos")>-1 :
		labcomp		= lab[8:] 
		if labcomp[-1]=="u" : labcomp=labcomp[:-1] 
		labconverted	= r"$D(\epsilon_{\mathrm{F}})$"+r'$_{{{}}}$'.format(labcomp)
	elif lab =="T" :
		labconverted = r"$T$"
	elif lab =="Tfit" :
		unitlab	= ' '+r'[K]'
		labconverted = r"$T_{\mathrm{fit}}$"
	elif lab =="betafit" :
		unitlab	= ' '+r'[K]'
		labconverted = r"$\beta_{\mathrm{fit}}$"
	elif lab =="Treal" :
		unitlab	= ' '+r'[K]'
		if paper : 
			labconverted = r"$T$"
		else : 
			labconverted = r"$T_{\mathrm{real}}$"
	elif (lab[1:] =="zdegMatLtcorrTDiag") or (lab[1:] =="xdegMatLtcorrTDiag") :
		oplab	= lab[0]
		oplab	= oplab + ';diag.'
		direct	= lab[1]+lab[1]
		if direct == "xx" : direct = "xy"
		unitlab	= ' '+r'[$\mu_B^2$]'
		labconverted = r'${\chi}$'+r'$^{{{}}}_{{{}}}$'.format(oplab, direct)+r'$(\beta/2)$'	#+' '+r'[$\mu_B^2$]'
	elif (lab[1:] =="zdegMatLtcorrTDiagBeta") or (lab[1:] =="xdegMatLtcorrTDiagBeta") :
		oplab	= lab[0]
		oplab	= oplab + ';Curie'
		direct	= lab[1]+lab[1]
		if direct == "xx" : direct = "xy"
		unitlab	= ' '+r'[$\mu_B^2$/eV]'
		labconverted = r'${\chi}$'+r'$^{{{}}}_{{{}}}$'.format(oplab, direct) #+' '+r'[$\mu_B^2$]'
	elif (lab[1:] =="zdegMatLtcorrT") or (lab[1:] =="xdegMatLtcorrT") :
		oplab	= lab[0]
		direct	= lab[1]+lab[1]
		if direct == "xx" : direct = "xy"
		unitlab	= ' '+r'[$\mu_B^2$]'
		labconverted = r'${\chi}$'+r'$^{{{}}}_{{{}}}$'.format(oplab, direct)+r'$(\beta/2)$'	#+' '+r'[$\mu_B^2$]'
	elif (lab[1:] =="zdegMatLtcorrTBeta") or (lab[1:] =="xdegMatLtcorrTBeta") :
		oplab	= lab[0]
		direct	= lab[1]+lab[1]
		if direct == "xx" : direct = "xy"
		unitlab	= ' '+r'[$\mu_B^2$]'
		labconverted = r'$\beta {\chi}$'+r'$^{{{}}}_{{{}}}$'.format(oplab, direct)+r'$(\beta/2)$'	#+' '+r'[$\mu_B^2$]'
	elif lab.find("chi")>-1 and lab.find("StaticTotal")>-1 :
		oplab	= lab[3]
		direct	= lab[ lab.find("Total")+5: ]
		unitlab	= ' '+r'[$\mu_B^2$/eV]'
		labconverted = r'$\chi$'+r'$^{{{}}}_{{{}}}$'.format(oplab, direct) #+' '+r'[$\mu_B^2$/eV]'
	elif lab.find("chitot")>-1 and lab.find("Staticwn0")>-1 :
		oplab	= lab[6].upper()
		direct	= lab[ lab.find("wn0")+3: ]
		unitlab	= ' '+r'[$\mu_B^2$/eV]'
		oplab	= oplab + ';VV'
		labconverted = r'$\chi$'+r'$^{{{}}}_{{{}}}$'.format(oplab, direct)#+r'$^{\mathrm{;VV}}$'		#+' '+r'[$\mu_B^2$/eV]'
	elif lab=="chistatic" or lab=="chi" : 
		unitlab	= ' '+r'[$\mu_B^2$/eV]'
		labconverted = r'$\chi$'
	elif lab=="alpha" :
		labconverted = r'$\alpha$'
	elif lab=="chitau" : 
		unitlab	= ' '+r'[$\mu_B^2$]'
		labconverted = r'${\chi}(\tau=\beta/2)$'
	elif lab=="JT" or lab=="JTD" : 
		unitlab	= ' '+r'[eV]'
		labconverted = r'$\Delta_{xy}$'
	elif lab=="tNNN" :
		unitlab	= ' '+r'[eV]'
		labconverted = r'$t_{\mathrm{NNN}}$'
	elif lab=="EEtot" :
		unitlab	= ''
		labconverted = r'$S_{E}$'
	elif lab.find("EEorb1_")>-1 :
		unitlab	= ''
		mulab	= lab.split("_")[-1]
		labconverted = r'$S_{E}^{(1)}$' + r'$[{}]$'.format(mulab)
	elif lab.find("EEorb2_")>-1 :
		unitlab	= ''
		mu1lab	= lab.split("_")[-2]
		mu2lab	= lab.split("_")[-1]
		labconverted = r'$S_{E}^{(2)}$' + r'$[{},{}]$'.format(mu1lab,mu2lab)
	elif lab.find("MInfoorb1_")>-1 :
		unitlab	= ''
		mulab	= lab.split("_")[-1]
		labconverted = r'$I$' + r'$[{}]$'.format(mulab)
	elif lab.find("MInfoorb2_")>-1 :
		unitlab	= ''
		mu1lab	= lab.split("_")[-2]
		mu2lab	= lab.split("_")[-1]
		labconverted = r'$I$' + r'$[{},{}]$'.format(mu1lab,mu2lab)
	elif lab.find("double_0")>-1 :
		unitlab	= ''
		labconverted = r'$\langle n_{\uparrow} n_{\downarrow} \rangle$' 
	if unit : 
		#if len(unitlab) > 0 : 
		#	unitlab = " " + unitlab
		labconverted = labconverted  + unitlab
	return labconverted

def markerdefault( ind ) :
	#defarr = [".", ",", "o", "v", "^", "<", ">", "1", "2", "3", "4", "8", "s", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d"]
	defarr = ["o", "s", "^", "v", "<", ">", "1", "2", "3", "4", "8", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d"]
	return defarr[ind]

def chdatFtn( basis ) :
	if basis.find("t")>-1 : 
		return [ 0, 2, 4 ]
	if basis.find("j")>-1 : 
		return [ 0, 1, 4 ]
def fnameDOSFtn( basis ) :
	if basis.find("t")>-1 : 
		if basis.find("z")>-1 : 
			return "ldos"
		else : 
			return "tldos"
	if basis.find("j")>-1 : 
		return "ldos"
