# A vertex-centric representation for adaptive diamond-kite meshes
#
# usage: python R INPUT -- all arguments optional

from __future__ import print_function
from sys import argv
from random import random

# depth of refinement
R= int(argv[1]) if len(argv)>1 else 6

# file containing initial mesh
INPUT= argv[2] if len(argv)>2 else "-"

# size of base mesh
N=6

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
STAR=array(6+1)
for d in range(3,6+1):
	STAR[d]=array(12)
# basic stars from Eppstein (2014) fig 4
STAR[3][0]=[W[0],W[4],W[8]]
STAR[4][0]=[W[0],W[4],W[7],W[9]]
STAR[5][0]=[W[0],W[3],W[5],W[7],W[9]]
STAR[6][0]=[W[0],W[2],W[4],W[6],W[8],W[10]]
# rotate basic stars
for d in range(3,6+1):
	STAR[d][1]=rotate(STAR[d][0],W[1])
	for k in range(2,12):
		STAR[d][k]=rotate(STAR[d][k-2],W[2])

# opposites
OPP=array(6+1)
for d in range(3,6+1):
	OPP[d]=array(12)
	OPP[d][0]=array(2*d)
# basic opposites for degree 3
OPP[3][0][0]=[0,1,0]
OPP[3][0][1]=[0,2,0]
OPP[3][0][2]=multiply(OPP[3][0][0],W[4])
OPP[3][0][3]=multiply(OPP[3][0][1],W[4])
OPP[3][0][4]=multiply(OPP[3][0][2],W[4])
OPP[3][0][5]=multiply(OPP[3][0][3],W[4])
# basic opposites for degree 4
OPP[4][0][0]=[0,1,0]
OPP[4][0][1]=[0,2,0]
OPP[4][0][2]=multiply([0,1,0],W[4])
OPP[4][0][3]=[0,0,0]
OPP[4][0][4]=multiply([2,2,1],W[7])
OPP[4][0][5]=multiply([1,1,0],W[7])
OPP[4][0][6]=[1,-1,0]
OPP[4][0][7]=[0,0,0]
# basic opposites for degree 5
OPP[5][0][0]=[0,1,0]
OPP[5][0][1]=[0,0,0]
OPP[5][0][2]=multiply([2,2,1],W[3])
OPP[5][0][3]=multiply([1,1,0],W[3])
OPP[5][0][4]=multiply([2,2,1],W[5])
OPP[5][0][5]=multiply([1,1,0],W[5])
OPP[5][0][6]=multiply([2,2,1],W[7])
OPP[5][0][7]=multiply([1,1,0],W[7])
OPP[5][0][8]=[1,-1,0]
OPP[5][0][9]=[0,0,0]
# basic opposites for degree 6
OPP[6][0][0]=[2,2,1]
OPP[6][0][1]=[1,1,0]
for i in range(2,2*6):
	OPP[6][0][i]=multiply(OPP[6][0][i-2],W[2])
# rotate basic opposites
for d in range(3,6+1):
	OPP[d][1]=rotate(OPP[d][0],W[1])
	for k in range(2,12):
		OPP[d][k]=rotate(OPP[d][k-2],W[2])

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
def addvertex(a,b,m,d,k,n):
	global nU
	a,b,m=normalize(a,b,m)
	if not (a,b,m) in V:
		z=(a+b*w2)/(3**m)
		t=(a,b,m)
		U[nU]=V[a,b,m]=DOTTED({'t':t,'d':d,'k':k,'n':n,'z':z,'w':False,'id':nU})
		nU=nU+1
	return V[a,b,m]

# add vertex in base mesh
def basevertex(w,d,k):
	a,b=w.real,w.imag
	v=addvertex(a,b,0,d,k,0)
	v.w=w

# base mesh: hexagonal grid
def basemesh(N):
	# use complex numbers for easy addition
	for i in range(12):
		a,b,m=W[i]
		W[i]=complex(a,b)
	for i in range(N):
		for j in range(N):
			c=0+i*(W[2]+W[4])+j*3*W[0]
			basevertex(c,3,0)
			for k in range(0,12,4):
				basevertex(c+W[k],6,0)
			for k in range(2,12,4):
				basevertex(c+W[k],3,2)
			if i<N-1 and j<N-1:
				basevertex(c+W[0]+W[2],3,0)
	boundary()

# base mesh: mark boundary vertices
def boundary():
	for k in V:
		v=V[k]
		if isboundary(v):
			v.d=0
			v.k=0

def isboundary(v):
	for k in range(0,12,2):
		w=v.w+W[k]
		a,b=w.real,w.imag
		a,b,m=normalize(a,b,0)
		if not (a,b,m) in V:
			return True
	return False

# load mesh from csv file
def loadmesh(filename):
	N=0
	for line in open(filename):
		if N>0:
			a,b,m,d,k,n=list(map(int,line.split(",")))
			addvertex(a,b,m,d,k,n)
		N=N+1
	addredundant()

# add vertices absent in file
# vertices of degree 3 that are adjacent to a vertex of degree 6
def addredundant():
	for k in range(len(U)):
		v=U[k]
		if v.d==6:
			for i in range(v.d):
				a,b,m=adj(v,i)
				if not (a,b,m) in V:
					k=(6+2*i+v.k)%12
					addvertex(a,b,m,3,k,v.n)

# stars
def adda(v,i,A):
	w=A[v.d][v.k][i]
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

def makestar(v):
	s=array(v.d)
	for i in range(v.d):
		a,b,m=adj(v,i)
		assert not (a,b,m) in V,(a,b,m)
		s[i]=addvertex(a,b,m,0,0,0)
	return s

# local subdivision step
def subdivide(v):
	assert v.d==6,v.d
	k=v.k
	n=v.n
	s0=star(v)
	v.k=(k+1)%12
	v.n=n+1
	s1=makestar(v)
	for j in range(v.d):
		w=s1[j]
		w.d=3
		w.k=(6+2*j+k+1)%12
		w.n=n+1
	for j in range(v.d):
		w=s0[j]
		kk=(6+2*j+k)%12
		if w.d==0:
			pass
		elif w.d==3:
			w.d=4
			w.k=(kk+4)%12
		elif w.d==4:
			w.d=5
			if w.k==kk:
				w.k=(kk+4)%12
			else:
				w.k=(kk-4)%12
		elif w.d==5:
			w.d=6
			w.k=(kk-1)%12
			w.n=n+1
		else:
			assert False,w.d
	return s0

# opposites
def opp(v,i):
	return adda(v,i,OPP)

def opposites(v):
	s=array(v.d)
	for i in range(0,2*v.d,2):
		a,b,m=opp(v,i)
		if not (a,b,m) in V:
			a,b,m=opp(v,i+1)
		s[i//2]=V[a,b,m]
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

# diamond or kite?
def facetype(k):
	v1,v2,v3,v4=k
	a1,b1,m1=U[v1].t
	a2,b2,m2=U[v2].t
	a3,b3,m3=U[v3].t
	a4,b4,m4=U[v4].t
	m=max(m1,m2,m3,m4)
	d1=3**(m-m1); a1=a1*d1; b1=b1*d1
	d2=3**(m-m2); a2=a2*d2; b2=b2*d2
	d3=3**(m-m3); a3=a3*d3; b3=b3*d3
	d4=3**(m-m4); a4=a4*d4; b4=b4*d4
	if (a1+a3)==(a2+a4) and (b1+b3)==(b2+b4):
		return "qd"
	else:
		return "qk"

# color face according to level
def facetype(k):
	v1,v2,v3,v4=k
	n=max(U[v1].n,U[v2].n,U[v3].n,U[v4].n)
	return str(n)+" qt"

# refinement
queue=set()

# Eppstein's prerequisite structure of replacement operations
# must select all adjacent vertices beforehand: refine(w) can change v.d
# v.d!=6 can happen if object is near the boundary when base mesh is too small
def refine(v):
	if v.d==4:
		w0=adjacent(v,0)
		w1=adjacent(v,1)
		refine(w0)
		refine(w1)
	elif v.d==5:
		w0=adjacent(v,0)
		refine(w0)
	#assert v.d==6,v
	if v.d!=6:
		return
	if v.n>=R:
		return
	s=subdivide(v)
	for w in s:
		if w.d==4:
			queue.add(w.t)
	queue.add(v.t)

# refinement criterion: uniform
def needsrefinement(v):
	return v.n < R

# refinement criterion: random
def needsrefinement(v):
	return v.n < R and random() < 0.45

# refinement criterion: implicit curve
def needsrefinement(v):
	z = min([f(w)*f(v) for w in star(v)])
	return v.n < R and z<=0

# implicit curve (Taubin, 1994)
def f(v):
	z=v.z-complex(8,4)
	x,y=z.real,z.imag
	z=0.004+0.110*x-0.177*y-0.174*x*x+0.224*x*y-0.303*y*y-0.168*x*x*x+0.327*x*x*y-0.087*x*y*y-0.013*y*y*y+0.235*x*x*x*x-0.667*x*x*x*y+0.745*x*x*y*y-0.029*x*y*y*y+0.072*y*y*y*y;
	return z

# reduction: vertices of degree 3 that are adjacent to a vertex of degree 6
#    same as vertices of degree 3 that are adjacent to an internal vertex
def redundant(v):
	if v.d!=3:
		return False
	for w in star(v):
		#if w.d!=0:
		if w.d==6:
			return True
	return False

# main -----------------------------------------------------------------------

# initial mesh
if INPUT=="-":
	basemesh(N)
else:
	loadmesh(INPUT)

# refine mesh
queue = {k for k in V if V[k].d==6}
while queue:
	k=next(iter(queue))
	queue.remove(k)
	v=V[k]
	if needsrefinement(v):
		refine(v)

# reconstruct mesh
for k in V:
	if V[k].d>0:
		addfaces(V[k])

# output PostScript
print('''\
%!PS-Adobe-2.0 EPSF-2.0
%%BoundingBox: 0 0 500 300
%%BoundingBox: 20 25 455 295
%%BoundingBox: 135 80 390 240
50 50 translate
25 dup scale
0 setlinewidth
1 setlinejoin
1 setlinecap
/c { 0.02 0 360 arc fill } bind def
/p { 0.04 0 360 arc fill } bind def
/l { moveto lineto stroke } bind def
/q { moveto lineto lineto lineto closepath fill } bind def
/p0 { 0.5 0.5 0.5 setrgbcolor p } bind def
/p3 { 0 0 1 setrgbcolor p } bind def
/p4 { 0 1 0 setrgbcolor p } bind def
/p5 { 1 0 1 setrgbcolor p } bind def
/p6 { 1 0 0 setrgbcolor p } bind def
/qd { 1 0.8 0.8 setrgbcolor q } bind def
/qk { 0.8 0.8 1 setrgbcolor q } bind def
/qt { 0.9 6 div mul 0.1 add 1 sub neg dup 1 setrgbcolor q } bind def
''')
print("%=ARG",R,INPUT)
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

print("")
print("% curve")
print("1 0 0 setrgbcolor")
for k in E:
	v0=U[k[0]]; z0=v0.z; f0=f(v0)
	v1=U[k[1]]; z1=v1.z; f1=f(v1)
	if f0*f1<=0:
		t=(f0-0)/(f0-f1)
		print((1-t)*z0.real+t*z1.real,(1-t)*z0.imag+t*z1.imag, "c")

"""
print ""
print "% vertices"
print "0 0 0 setrgbcolor"
for k in U:
	v=U[k]
	print v.z.real, v.z.imag, "p" + str(v.d) # + "0" # + str(v.k)
"""

print("")
print("showpage")
print("%%EOF")

#exit()	# if extra files not needed

# output CSV
print("")
PREFIX="%=CSV"
print(PREFIX, "a,b,m,d,k,n")
for k in range(len(U)):
	v=U[k]
	print(PREFIX, ','.join(map(str,[v.t[0],v.t[1],v.t[2],v.d,v.k,v.n])))

# output OBJ
print("")
PREFIX="%=OBJ"
print(PREFIX, "OBJ")
print(PREFIX, len(U),len(F),0)
for k in range(len(U)):
	v=U[k]
	print(PREFIX, "v", v.z.real, v.z.imag, 0)
for k in F:
	print(PREFIX, "f",k[0]+1,k[1]+1,k[2]+1,k[3]+1)

# output OFF
print("")
PREFIX="%=OFF"
print(PREFIX, "OFF")
print(PREFIX, len(U),len(F),0)
for k in range(len(U)):
	v=U[k]
	print(PREFIX, v.z.real, v.z.imag, 0)
for k in F:
	print(PREFIX, "4",k[0],k[1],k[2],k[3])

# done
print("")
print("%=EOF")

