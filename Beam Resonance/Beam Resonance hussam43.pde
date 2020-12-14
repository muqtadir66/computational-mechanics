TITLE 
'H10.2 hussam43'
SELECT
errlim=1e-4
ngrid=15
spectral_colors
stages = 10

COORDINATES
cartesian3

VARIABLES
u	!Displacement in x
v !Displacement in y
w !Displacement in z

DEFINITIONS
mag=1000

Lx=2
Ly=.4
Lz=.2
E=if z<0.3333*Lz then 120e9 else if z<0.6666*Lz then 74e9 else 47e9
nu=if z<0.3333*Lz then 0.316 else if z<0.6666*Lz then 0.415 else 0.36
rho=if z<0.3333*Lz then 7640 else if z<0.6666*Lz then 19320 else 7310
G=E/(2*(1+nu))

omega = 565+0.1*stage

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
gyz = dz(v) + dy(w)
gxz = dz(u) + dx(w)
gxy = dy(u) + dx(v)

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
u:	dx(sx)+dy(sxy)+dz(sxz)= -rho*omega^2*u
v:	dx(sxy)+dy(sy)+dz(syz) =-rho*omega^2*v
w:	dx(sxz)+dy(syz)+dz(sz) -rho*9.81=-rho*omega^2*w !gravity in thickness direction is really just to get it started

EXTRUSION
surface 'bottom' z=0
surface 'top' z=Lz

BOUNDARIES
surface 'bottom'
		load(u)=0
		load(v)=0 
		load(w)=0
surface 'top'
		load(u)=0
		load(v)=0 
		load(w)=0

  	REGION 1
    	START(0,0) !y=0 surface:
		load(u)=0
		load(v)=0 
		load(w)=0
    	LINE TO (Lx,0) !x=Lx surface
		value(u)=0
		load(v)=0
		value(w)=0
		LINE TO (Lx,Ly) !y=Ly surface
		load(u)=0
		load(v)=0
		load(w)=0
		LINE TO (0,Ly) !x=0 surface
		load(u)=0 
		load(v)=0
		value(w)=0
		LINE TO CLOSE


MONITORS
!contour(u) painted on x=0
!contour(v) painted on x=0
!contour(w) painted on x=0
!grid(x+u,y+v,z+w)

PLOTS
	grid(x+u*mag, y+v*mag, z+w*mag)
 	contour(sxz) on surface z=Lz
	elevation(sx,sy,sz,syz,sxz,sxy) from (0,0,0) to (0,0,Lz)
	elevation(w) from (0,Ly/2,Lz/2) to(Lx,Ly/2,Lz/2)
	history(w) at (Lx/2,Ly/2,Lz/2) !the middle will show a lot of displacement
	
	summary
		report val(u,Lx,Ly,Lz)
		report val(v,Lx,Ly,Lz)
		report val(w,Lx,Ly,Lz)
		report val(w,Lx/2,Ly/2,Lz)
end
