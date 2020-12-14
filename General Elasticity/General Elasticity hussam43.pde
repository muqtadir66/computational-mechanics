TITLE 'H11.2 hussam43' { the problem identification }
COORDINATES cartesian3 { coordinate system, 1D,2D,3D, etc }
VARIABLES { system variables }
 u
v
w
SELECT { method controls }
ngrid = 5
DEFINITIONS { parameter definitions }
mag = 200
Lx = .1
Ly = .1
Lz = .1
sapplied=-5e6
!Matrix components
C11 =20e9
C12 = -9e9
C13 = 19e9
C14 =1e9
C15 = -10e9
C16 = 12e9
C22 = 52e9
C23 = 0
C24 = -9e9
C25 = 7e9
C26 = 0
C33 = 91e9
C34 = 0
C35 = 0
C36 = 2e9
C44 = 17e9
C45 = 0
C46 = -12e9
C55 = 21e9
C56 = 5e9
C66 = 48e9
C21=C12 !matrix is symmetric
C31=C13
C32=C23
C41=C14
C42=C24
C43=C34
C51=C15
C52=C25
C53=C35
C61=C16
C62=C26
C63=C36
C54=C45
C64=C46
C65=C56
!Strain definitions from displacements
ex = dx(u)
ey = dy(v)
ez = dz(w)
gyz = dy(w) + dz(v)
gxz = dx(w) + dz(u)
gxy = dx(v) + dy(u)
!Hookes Law
sx = C11*ex + C12*ey + C13*ez + C14*gyz + C15*gxz + C16*gxy
sy = C21*ex + C22*ey + C23*ez + C24*gyz + C25*gxz + C26*gxy
sz = C31*ex + C32*ey + C33*ez + C34*gyz + C35*gxz + C36*gxy
syz = C41*ex + C42*ey + C43*ez + C44*gyz + C45*gxz + C46*gxy
sxz = C51*ex + C52*ey + C53*ez + C54*gyz + C55*gxz + C56*gxy
sxy = C61*ex + C62*ey + C63*ez + C64*gyz + C65*gxz + C66*gxy
EQUATIONS { PDE's, one for each variable }
u: dx(sx) + dy(sxy) + dz(sxz) = 0
v: dx(sxy) + dy(sy) + dz(syz) = 0
w: dx(sxz) + dy(syz) + dz(sz) = 0
EXTRUSION
surface 'bottom' z = 0
surface 'top' z = Lz
BOUNDARIES { The domain definition }
surface 'bottom'
load(u) =0
value(v) = 0
load(w) = 0
surface 'top'
value(v) = 0
 REGION 1 { For each material region }
 START(0,0) !y=0
load(u) =0
value(v) = 0
load(w) = 0
 LINE TO (Lx,0) !x = Lx
load(u) =sapplied
LINE TO (Lx,Ly) !y = Ly
load(u) =0
LINE TO (0,Ly) !x = 0
load(u) =-sapplied
!value(u) = 0
LINE TO CLOSE
PLOTS { save result displays }
grid(x+mag*u, y+mag*v, z+mag*w)
CONTOUR(sxz) on x = Lx/3
CONTOUR(v) on x = Lx/3
CONTOUR(u) on x = Lx/3
CONTOUR(ey) on x = Lx/3
SUMMARY
report val(sy, Lx/2,Ly/2,Lz/2)
report val(ex, Lx/2,Ly/2,Lz/2)
report val(ey, Lx/2,Ly/2,Lz/2)
report val(ez, Lx/2,Ly/2,Lz/2)
report val(gyz, Lx/2,Ly/2,Lz/2)
report val(gxz, Lx/2,Ly/2,Lz/2)
report val(gxy, Lx/2,Ly/2,Lz/2)
report val(sy, Lx*.8, Ly*.8, Lz*.8)
report val(ex, Lx*.8,Ly*.8,Lz*.8)
report val(ey, Lx*.8,Ly*.8,Lz*.8)
report val(ez, Lx*.8,Ly*.8,Lz*.8)
report val(gyz, Lx*.8,Ly*.8,Lz*.8)
report val(gxz, Lx*.8,Ly*.8,Lz*.8)
report val(gxy, Lx*.8,Ly*.8,Lz*.8)
report val(sapplied/ex, Lx,Ly,Lz) as 'Effective Stiffness in x'
report val(sapplied/ex, Lx/2,Ly/2,Lz/2) as 'Effective Stiffness'
report val(sapplied/ex, Lx/2,.8*Ly,.8*Lz) as 'Effective Stiffness'
report val(sapplied/ex, Lx/2,.8*Ly,.4*Lz/2) as 'Effective Stiffness'
END