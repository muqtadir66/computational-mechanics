TITLE 
'H7.2 hussam43'
SELECT
errlim=1e-4
ngrid=6
spectral_colors
COORDINATES
cartesian3

VARIABLES
u	!Displacement in x
v !Displacement in y
w !Displacement in z

DEFINITIONS
Lx=7
Ly=3
Lz=3
nu=0.42
E=400e6
G=E/(2*(1+nu))

s_applied = 150e6 !applied stress


C11 =E*(1-nu)/(1+nu)/(1-2*nu)
C22 = C11
C33 = C11

C12 = E*nu/(1+nu)/(1-2*nu)
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
u:	dx(sx)+dy(sxy)+dz(sxz)=0
v:	dx(sxy)+dy(sy)+dz(syz)=0
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
load(u)=s_applied 
load(v)=0
load(w)=0

  	REGION 1
    	START(0,0) !y=0 surface:
		load(u)=0
		load(v)=0
		load(w)=0
    	LINE TO (Lx,0) !x=Lx surface
		load(u)=0
		load(v)=0
		load(w)=s_applied
		LINE TO (Lx,Ly) !y=Ly surface
		load(u)=0
		load(v)=0
		load(w)=0
		LINE TO (0,Ly) !x=0 surface
		load(u)=0
		load(v)=0
		load(w)=-s_applied
		LINE TO CLOSE


MONITORS
contour(u) painted on x=0
contour(v) painted on x=0
contour(w) painted on x=0
grid(x+u,y+v,z+w)

PLOTS
	grid(x+u, y+v, z+w)
 	contour(u) on surface z=0
	elevation(sx,sy,sz) from (0,0,0) to (0,0,Lz)
	
	summary
		report val(u,Lx,Ly,Lz)
		report val(v,Lx,Ly,Lz)
		report val(w,Lx,Ly,Lz)
end
