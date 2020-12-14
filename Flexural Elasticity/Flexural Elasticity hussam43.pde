TITLE 'H8.2 hussam43'
COORDINATES cartesian3
VARIABLES        { system variables }
	u   
	v
	w
SELECT         { method controls }
ngrid=5
DEFINITIONS    { parameter definitions }
mag=1e4

grav=9.81

wM = .3
hM = .75
wG = .5
hG = .3
Lz = 2.5

E= if y<0 then 329e9 else 74e9
nu= if y<0 then 0.29 else 0.415
rho= if y<0 then 10220 else 19320

G=E/(2*(1+nu))

C11 =E*(1-nu)/((1+nu)*(1-2*nu))
C22 = C11
C33 = C11

C12 = E*nu/((1+nu)*(1-2*nu))
C13 = C12
C21 = C12
C23 = C12
C31 = C12
C32 = C12

sigzTip =328964.2765*y + 54279.10562

!!Strain
!Axial Strain
ex=dx(u)
ey=dy(v)
ez=dz(w)
!Engineering Shear Strain
gxy=dx(v)+dy(u)
gyz=dy(w)+dz(v)
gxz=dz(u)+dx(w)

!Stress via Hooke's Laws
!Axial Stress
sx = C11*ex + C12*ey + C13*ez
sy = C21*ex + C22*ey + C23*ez
sz = C31*ex + C32*ey + C33*ez
!Shear Stress
syz= G*gyz
sxz= G*gxz 
sxy= G*gxy

EQUATIONS        { PDE's, one for each variable }
!Fnet = 0
u:		dx(sx)+dy(sxy)+dz(sxz)=0
v:		dx(sxy)+dy(sy)+dz(syz)-rho*grav=0
w:	dx(sxz)+dy(syz)+dz(sz)=0

EXTRUSION
surface 'bottom' z=0
surface 'top' z=Lz

BOUNDARIES       { The domain definition }
surface 'bottom'
	value(u)=0
	value(v)=0
	value(w)=0
surface 'top'
	load(u)=0
	load(v)=0
	load(w)=sigzTip

  REGION 1       { For each material region }
    START(-wM/2,-hM)
	load(u)=0
	load(v)=0
	load(w)=0
	LINE TO(wM/2,-hM) TO(wM/2,0) TO(wG/2,0) TO(wG/2,hG)TO(-wG/2,hG)TO(-wG/2,0)TO(-wM/2,0)TO CLOSE
PLOTS            
	grid(x+mag*u,y+mag*v,z+mag*w)
	contour(sz) painted on z=Lz/2
	contour(ez) painted on z=Lz/2
	elevation(ex,ey,ez) from (0,0,0) to (0,0,Lz)
SUMMARY
	report val(v,0,0,Lz)
END

