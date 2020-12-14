TITLE 
'H9.2 hussam43'
SELECT
errlim=1e-8
ngrid=12
spectral_colors
COORDINATES
cartesian3

VARIABLES
u	!Displacement in x
v !Displacement in y
w !Displacement in z

DEFINITIONS
Lz=5
diam=1 !diameter
nu=0.320
E=4.1e9
G=E/(2*(1+nu))

r = sqrt(x^2+y^2)
phi=atan2(y,x)
!Calculate the twist angle in flexPDE directly:
thetatest = atan2(y+v,x+u)-phi !twist angle, but possibly outside the range -Pi to Pi
theta = if(thetatest<-pi) then thetatest+2*pi else if (thetatest>pi) then thetatest-2*pi else thetatest !twist angle


mag = 500 !magnification

C11 =E*(1-nu)/((1+nu)*(1-2*nu))
C22 = C11
C33 = C11

C12 = E*nu/((1+nu)*(1-2*nu))
C13 = C12
C21 = C12
C23 = C12
C31 = C12
C32 = C12

!! Strain
!Axial Strain
ex=dx(u)
ey=dy(v)
ez=dz(w)
!Engineering Shear Strain
gxy=(dx(v)+dy(u))
gyz=(dy(w)+dz(v))
gxz=(dz(u)+dx(w))

!!Stress via Hooke's law
!Axial Stress
sx = C11*ex+C12*ey+C13*ez
sy = C21*ex+C22*ey+C23*ez
sz = C31*ex+C32*ey+C33*ez
!Shear stress
sxy=G*gxy
sxz=G*gxz
syz=G*gyz


EQUATIONS
!FNet = 0
u:		dx(sx)+dy(sxy)+dz(sxz)=0
v: 	dx(sxy)+dy(sy)+dz(syz)=0
w:	dx(sxz)+dy(syz)+dz(sz)=0

EXTRUSION
surface 'bottom' z=0
surface 'top' z=Lz

BOUNDARIES
surface 'bottom'
value(u)=0
value(v)=0
value(w)=0
surface 'top'
load(u)=-(y)*2e6
load(v)=(x)*2e6
load(w)=0

  	REGION 1
    	START(-diam/2,0)
		load(u)=0
		load(v)=0
		load(w)=0
		arc(center=0,0) angle=360

  		TO CLOSE


PLOTS
	grid(x+u, y+v, z+w)
	grid(x+u*mag, y+v*mag, z+w*mag)
	contour(sxz) painted on surface z=Lz
	contour(syz) painted on surface z=Lz
	contour(sz) painted on surface z=Lz
	contour(u) painted on surface z=Lz
	contour(theta) painted on surface z=Lz
	elevation(u,v,w) from (0,0,0) to (0,0,Lz)
	elevation(theta) from (diam/4,0,0) to (diam/4,0,Lz)
	
	summary
		report val(u,diam/2,0,Lz)
		report val(v,diam/2,0,Lz)
		report val(w,diam/2,0,Lz)
		report val(theta,diam/2,0,Lz)
end

