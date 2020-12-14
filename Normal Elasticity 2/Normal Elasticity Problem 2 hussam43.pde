TITLE 'H6.2 Assignment hussam43'
COORDINATES cartesian3
 VARIABLES
u
v
w
DEFINITIONS
E=0.46e9
nu=0.325
mag=0.5*globalmax(magnitude(x,y,z))/globalmax(magnitude(u,v,w))
Lx = 0.3
Ly = 3
Lz = 0.2
Fy =414.7784425 
Ay =Lz*Lx
sigmaX = 0
sigmaY = Fy/Ay
sigmaZ = 0
epsilonx = dx(u)
epsilony = dy(v)
epsilonz = dz(w)
C11=E/((1+nu)*(1-2*nu))*(1-nu)
C12=E/((1+nu)*(1-2*nu))*nu
C13=C12
C21=C12
C22=C11
C23=C12
C31=C12
C32=C12
C33=C11
Sx = C11*epsilonx+C12*epsilony+C13*epsilonz
Sy = C21*epsilonx+C22*epsilony+C23*epsilonz
Sz = C31*epsilonx+C32*epsilony+C33*epsilonz
LxNew = val(Lx+u,Lx,Ly,Lz)
LyNew = val(Ly+v,Lx,Ly,Lz)
LzNew = val(Lz+w,Lx,Ly,Lz)
VolOrig=Lx*Ly*Lz
VolNew = LxNew*LyNew*LzNew
VolChange=VolNew-VolOrig
EQUATIONS
u: dx(Sx)=0
v: dy(Sy)=0
w: dz(Sz)=0
EXTRUSION
surface 'bottom' z=0
surface 'top' z=Lz
BOUNDARIES
surface 'bottom'
load(u)=0
load(v)=0
value(w)=0
surface 'top'
load(u)=0
load(v)=0
load(w)=0
REGION 1
     START(0,0)
load(u)=0
value(v)=0
load(w)=0
LINE TO (Lx,0)
 load(u)=0
load(v)=0
load(w)=0
LINE TO (Lx,Ly)
 load(u)=0
load(v)= sigmaY
load(w)=0
LINE TO (0,Ly)
 value(u)=0
load(v)=0
load(w)=0
LINE TO CLOSE
PLOTS
 grid(x,y,z)
grid(x+u*mag,y+v*mag,z+w*mag)
SUMMARY
    report LxNew as 'New x length (m)'
report LyNew as 'New y length (m)'
report LzNew as 'New z length (m)'
report VolNew as 'New Volume (m^3):'
report VolChange as 'Change in volume (m^3)'
END