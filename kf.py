# A vertex-centric representation for kite fractals
# based on dk.py and the paper
# Fractal tilings based on kite- and dart-shaped prototiles by Fathauer (2001)
# https://doi.org/10.1016/S0097-8493(00)00134-5

from __future__ import print_function
from sys import argv

# depth of refinement
R= int(argv[1]) if len(argv)>1 else 4

# make array
def array(n):
	return list(range(n))

# normalize 3-adic lattice coordinates [a,b,m]=(a+b*w^2)/(3^m)
def normalize(a,b,m):
	a,b=int(a),int(b)
	while m>0 and a%3==0 and b%3==0:
		a=a/3
		b=b/3
		m=m-1
	return [a,b,m]

# add 3-adic lattice points with scaling
def add(v1,v2,n):
	a1,b1,m1=v1
	a2,b2,m2=v2
	m2=m2+n//2
	d1=3**m1
	d2=3**m2
	a=d2*a1+d1*a2
	b=d2*b1+d1*b2
	m=m1+m2
	return normalize(a,b,m)

# multiply 3-adic lattice points
# PolynomialRemainder[(a_1 + b_1 z)(a_2 + b_2 z),z^2-z+1,z]
#	z (a_2 b_1 + a_1 b_2 + b_2 b_1) + a_1 a_2 - b_1 b_2
def multiply(v1,v2):
	a1,b1,m1=v1
	a2,b2,m2=v2
	a=a1*a2-b1*b2
	b=a2*b1+a1*b2+b2*b1
	m=m1+m2
	return normalize(a,b,m)

# rotate star
def rotate(a,z):
	n=len(a)
	b=array(n)
	for i in range(n):
		b[i]=multiply(a[i],z)
	return b

# 3-adic lattice coordinates of basic directions; odd directions are shorter
W=array(12)
W[0]=[1,0,0]	# 1
W[1]=[1,1,1]	# (1+w^2)/3
W[2]=[0,1,0]	# w^2
for i in range(3,12):
	W[i]=multiply(W[i-2],W[2])

# stars
TYPE=[20,31,32,41,42,43,50,60]
STAR={}
for s in TYPE:
	d=s//10
	STAR[s]=array(12)
	STAR[s][0]=array(d)
# 20
STAR[20][0]=[W[4],W[8]]
# 31
STAR[31][0]=[W[3],W[6],W[9]]
# 41
STAR[41][0]=[W[3],W[6],W[9],W[0]]
STAR[41][0][3]=[1,0,1]
# 32
STAR[32][0]=[W[2],W[6],W[10]]
# 42
STAR[42][0]=[W[2],W[6],W[10],W[0]]
# 43
STAR[43][0]=[W[3],W[6],W[9],W[0]]
# 50
STAR[50][0]=[W[4],W[8],W[10],W[0],W[2]]
# 60
STAR[60][0]=[W[0],W[2],W[4],W[6],W[8],W[10]]
# rotate basic stars
for s in TYPE:
	STAR[s][1]=rotate(STAR[s][0],W[1])
	for k in range(2,12):
		STAR[s][k]=rotate(STAR[s][k-2],W[2])

# opposites
OPP={}
for s in TYPE:
	d=s//10
	if d<4: d=d-1
	OPP[s]=array(12)
	OPP[s][0]=array(d)
# 20
OPP[20][0][0]=[-2,0,0]
# 31
OPP[31][0][0]=[-1,1,0]
OPP[31][0][1]=[0,-1,0]
# 41
OPP[41][0][0]=[-1,1,0]
OPP[41][0][1]=[0,-1,0]
OPP[41][0][2]=[2,-1,1]
OPP[41][0][3]=[1,1,1]
# 32
OPP[32][0][0]=[-2,2,0]
OPP[32][0][1]=[0,-2,0]
# 42
OPP[42][0][0]=[-2,2,0]
OPP[42][0][1]=[0,-2,0]
OPP[42][0][2]=[4,-2,1]
OPP[42][0][3]=[2,2,1]
# 43
OPP[43][0][0]=[-1,1,0]
OPP[43][0][1]=[0,-1,0]
OPP[43][0][2]=[1,-1,0]
OPP[43][0][3]=[0,1,0]
# 50
OPP[50][0][0]=[-2,0,0]
OPP[50][0][3]=[2,2,1]
OPP[50][0][4]=multiply(OPP[50][0][3],W[2])
OPP[50][0][2]=multiply(OPP[50][0][3],W[10])
OPP[50][0][1]=multiply(OPP[50][0][3],W[8])
# 60
OPP[60][0][0]=[2,2,1]
for i in range(1,6):
	OPP[60][0][i]=multiply(OPP[60][0][i-1],W[2])
# rotate basic opposites
for s in TYPE:
	OPP[s][1]=rotate(OPP[s][0],W[1])
	for k in range(2,12):
		OPP[s][k]=rotate(OPP[s][k-2],W[2])

# Cartesian coordinates
w2=complex(0.5,0.866025403784438646763723170752936183471402626905190314027)

# vertex cloud
V={}

# topological data
U={}
E={}
F={}
nU=0

# use dot notation for dict -- https://stackoverflow.com/a/74214556/107090
class DOTTED(dict):
	__getattr__ = dict.get
	__setattr__ = dict.__setitem__
	__delattr__ = dict.__delitem__

# add vertex to cloud
def addvertex(a,b,m,s,k,n):
	global nU
	a,b,m=normalize(a,b,m)
	if not (a,b,m) in V:
		z=(a+b*w2)/(3**m)
		t=(a,b,m)
		d=s//10
		U[nU]=V[a,b,m]=DOTTED({'t':t,'d':d,'k':k,'n':n,'z':z,'id':nU,'s':s})
		nU=nU+1
	return V[a,b,m]

# base mesh
def basemesh():
	v=addvertex(0,0,0,60,0,0)
	s=thestar(v)
	for i in range(6):
		w=s[i]
		w.s=31
		w.d=3
		w.k=2*i
		w.n=0
	s=theopposites(v)
	for i in range(6):
		w=s[i]
		w.s=20
		w.d=2
		w.k=2*i+1
		w.n=1

def loadmesh(filename):
	N=0
	for line in open(filename):
		if N>0:
			a,b,m,d,k,n=list(map(int,line.split(",")))
			addvertex(a,b,m,d,k,n)
		N=N+1

# stars
def adda(v,i,A):
	w=A[v.s][v.k][i]
	return add(v.t,w,v.n)

def adj(v,i):
	return adda(v,i,STAR)

def adjacent(v,i):
	a,b,m=adj(v,i)
	return V[a,b,m]

def star(v):
	s=array(v.d)
	for i in range(v.d):
		s[i]=adjacent(v,i)
	return s

def thestar(v):
	s=array(v.d)
	for i in range(v.d):
		a,b,m=adj(v,i)
		s[i]=addvertex(a,b,m,0,0,0)
	return s

# opposites
def opp(v,i):
	return adda(v,i,OPP)

def opposite(v,i):
	a,b,m=opp(v,i)
	return V[a,b,m]

def opposites(v):
	s=array(v.d)
	for i in range(v.d):
		s[i]=opposite(v,i)
	return s

def theopposites(v):
	s=array(v.d)
	for i in range(v.d):
		a,b,m=opp(v,i)
		s[i]=addvertex(a,b,m,0,0,0)
	return s

# faces and edges
def addedge(v1,v2):
	a=[v1,v2]
	j=a.index(min(a))
	a=(a[j],a[(j+1)%2])
	E[a]=True

def addface(v1,v2,v3,v4):
	a=[v1.id,v2.id,v3.id,v4.id]
	j=a.index(min(a))
	a=(a[j],a[(j+1)%4],a[(j+2)%4],a[(j+3)%4])
	F[a]=True
	addedge(a[0],a[1])
	addedge(a[1],a[2])
	addedge(a[2],a[3])
	addedge(a[3],a[0])

def addfaces(v):
	s=star(v)
	o=opposites(v)
	for i in range(v.d):
		addface(v,s[i],o[i],s[(i+1)%v.d])

def facetype(k):
	v1,v2,v3,v4=k
	n=max(U[v1].n,U[v2].n,U[v3].n,U[v4].n)
	return str(n)+" qt"

# refinement
def refine(v):
	if v.s==0:
		return
	elif v.s==20:
		v.s=50
		v.d=5
	elif v.s==31:
		v.s=41
		v.d=4
	elif v.s==32:
		v.s=42
		v.d=4
	thestar(v)
	theopposites(v)
	if v.s==41:
		w=adjacent(v,3)
		w.s=32
		w.d=3
		w.k=v.k
		w.n=v.n+2
	elif v.s==42:
		w=adjacent(v,3)
		w.s=31
		w.d=3
		w.k=v.k
		w.n=v.n
	elif v.s==50:
		for j in [2,3,4]:
			w=adjacent(v,j)
			old= w.s!=0
			w.s=31
			w.d=3
			w.k=(v.k+2*j-6)%12
			w.n=v.n
			if old:
				w.s=43
				w.d=4
		for j in [2,3]:
			w=opposite(v,j)
			k1=w.k
			old= w.s!=0
			w.s=20
			w.d=2
			w.k=(v.k+2*j-5)%12
			w.n=v.n+1
			if old:
				k2=w.k
				if (k1+4)%12==k2:
					w.k=(k1+2)%12
				else:
					w.k=(k2+2)%12
				w.s=32
				w.d=3

# main -----------------------------------------------------------------------

# initial mesh
basemesh()

# refine mesh
for r in range(R):
	#print "%=U",r,R,len(U)
	for k in range(len(U)):
		v=U[k]
		if v.s<40:
			refine(v)
	#print "%=U",r,R,len(U)

# reconstruct mesh
for k in V:
	v=V[k]
	if v.s>40:
		addfaces(v)

# output PostScript
print('''\
%!PS-Adobe-2.0 EPSF-2.0
%%BoundingBox: 0 0 500 500
%%BoundingBox: 400 200 460 300
250 dup translate
90 dup scale
0 setlinewidth
1 setlinejoin
1 setlinecap
/c { 0.02 0 360 arc fill } bind def
/p { 0.005 0 360 arc fill } bind def
/p { 0.04 0 360 arc fill } bind def
/p { 0.01 0 360 arc fill } bind def
/l { moveto lineto stroke } bind def
/q { moveto lineto lineto lineto closepath fill } bind def
/p0 { 0.5 0.5 0.5 setrgbcolor p } bind def
/p20 { 0.5 0.5 1 setrgbcolor p } bind def
/p31 { 0 1 0 setrgbcolor p } bind def
/p32 { 0 1 1 setrgbcolor p } bind def
/p41 { 1 0 0 setrgbcolor p } bind def
/p42 { 1 0 1 setrgbcolor p } bind def
/p43 { 1.0 0.65 0 setrgbcolor p } bind def
/p50 { 1 1 0 setrgbcolor p } bind def
/p60 { 1 1 1 setrgbcolor p } bind def
/qt { 0.2 mul 0.1 add 1 sub neg dup 1 setrgbcolor q } bind def
/qt { 0.8 6 div mul 0.1 add 1 sub neg dup 1 setrgbcolor q } bind def
''')
print("%=ARG","R",R)
print("%=VEF",R,len(V),len(U),len(E),len(F))

print("")
print("% faces")
print("1 0.8 0.8 setrgbcolor")
for k in F:
	z1=U[k[0]].z
	z2=U[k[1]].z
	z3=U[k[2]].z
	z4=U[k[3]].z
	print(z1.real,z1.imag, z2.real,z2.imag, z3.real,z3.imag, z4.real,z4.imag,facetype(k))

print("")
print("% edges")
print("0 0 0 setrgbcolor")
print("0.01 setlinewidth")
for k in E:
	z1=U[k[0]].z
	z2=U[k[1]].z
	print(z1.real,z1.imag, z2.real,z2.imag, "l")

"""
print ""
print "% vertices"
print "0 0 0 setrgbcolor"
for k in U:
	v=U[k]
	if v.s<90:
		print v.z.real, v.z.imag, "p" + str(v.s)
"""

print("")
print("showpage")
print("%%EOF")

exit()	# if extra files not needed

# output CSV
print("")
PREFIX="%=CSV"
print(PREFIX, "a,b,m,s,k,n")
for k in range(len(U)):
	v=U[k]
	print(PREFIX, ','.join(map(str,[v.t[0],v.t[1],v.t[2],v.s,v.k,v.n])))

# done
print("")
print("%=EOF")

