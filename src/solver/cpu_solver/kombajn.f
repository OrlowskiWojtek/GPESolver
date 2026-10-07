program gaussian
implicit double precision(a-b,d-h,o-z)
implicit double complex(c)
parameter(N=75,ny=75,nz=16)
c     parameter(N=34,ny=21,nz=12)
c     allocatable cpsi(:,:,:)
allocatable cpsi1(:,:,:)
allocatable cpsi2(:,:,:)
c     allocatable cpsin(:,:,:)
allocatable cpsin1(:,:,:)
allocatable cpsin2(:,:,:)
c     allocatable cpsii(:,:,:)
allocatable cpsii1(:,:,:)
allocatable cpsii2(:,:,:)
allocatable dmin(:,:)
allocatable fi3do1(:,:,:)
allocatable fi3d1(:,:,:)
allocatable fi3d2(:,:,:)
allocatable fi3do2(:,:,:)
allocatable flhy1(:,:,:)
allocatable flhy2(:,:,:)
allocatable flhyo1(:,:,:)
allocatable flhyo2(:,:,:)
dimension cds(4),f(4)
dimension wg(16),ug(16)
dimension xnorma(2),xm(2),xncz(2),xnormak(2),ene(2),wd(2)
dimension vx(-n:n),vy(-ny:ny),vz(-nz:nz)
allocatable rdy(:,:,:)
common/xnorma/xnorma,wrl,wzl,dx,xncz,ggp11,dz,ggp12,wd,y0,omega
>,ggp22,ggp21
common/loc/xnormak,ene,eod11,eod22,eod12,ix,iy,iz
common/add/add,cdd,gamma1,gamma2,cdd2,cdd12,edd,edd2,edd12
common/xl/xl
common/skladowe/e1,e2,e3,e4,e5,e6,e7,e8,e9
common/gauss/wg,ug
ug(1)  = 0.005299532504175031d0
ug(2)  = 0.027712488463383700d0
ug(3)  = 0.067184398806084122d0
ug(4)  = 0.122297795822498501d0
ug(5)  = 0.191061877798678115d0
ug(6)  = 0.270991611171386315d0
ug(7)  = 0.359198224610370542d0
ug(8)  = 0.452493745081181287d0
ug(9)  = 0.547506254918818769d0
ug(10) = 0.640801775389629458d0
ug(11) = 0.729008388828613629d0
ug(12) = 0.808938122201321885d0
ug(13) = 0.877702204177501555d0
ug(14) = 0.932815601193915933d0
ug(15) = 0.972287511536616300d0
ug(16) = 0.994700467495824969d0

wg(1)  = 0.013576229705877088d0
wg(2)  = 0.031126761969323728d0
wg(3)  = 0.047579255841246303d0
wg(4)  = 0.062314485627767036d0
wg(5)  = 0.074797994408288354d0
wg(6)  = 0.084578259697501323d0
wg(7)  = 0.091301707522461820d0
wg(8)  = 0.094725305227534320d0
wg(9)  = 0.094725305227534320d0
wg(10) = 0.091301707522461820d0
wg(11) = 0.084578259697501323d0
wg(12) = 0.074797994408288354d0
wg(13) = 0.062314485627767036d0
wg(14) = 0.047579255841246303d0
wg(15) = 0.031126761969323728d0
wg(16) = 0.013576229705877088d0

ug(1) = 0.019855071751231884d0
ug(2) = 0.101666761293186630d0
ug(3) = 0.237233795041835507d0
ug(4) = 0.408282678752175098d0
ug(5) = 0.591717321247824902d0
ug(6) = 0.762766204958164493d0
ug(7) = 0.898333238706813370d0
ug(8) = 0.980144928248768116d0

wg(1) = 0.050614268145188129d0
wg(2) = 0.111190517226687235d0
wg(3) = 0.156853322938943644d0
wg(4) = 0.181341891689180991d0
wg(5) = 0.181341891689180991d0
wg(6) = 0.156853322938943644d0
wg(7) = 0.111190517226687235d0
wg(8) = 0.050614268145188129d0

xlp=11000/.05292 
ci=(0.d0,1.d0)
xncz(1)=4e4
xncz(2)=4e4
allocate (cpsi1(-N:N,-Ny:Ny,-Nz:Nz))
allocate (cpsin1(-N:N,-Ny:Ny,-Nz:Nz))
allocate (cpsii1(-N:N,-Ny:Ny,-Nz:Nz))
allocate (cpsi2(-N:N,-Ny:Ny,-Nz:Nz))
allocate (cpsin2(-N:N,-Ny:Ny,-Nz:Nz))
allocate (cpsii2(-N:N,-Ny:Ny,-Nz:Nz))
allocate (fi3do1(-N:N,-Ny:Ny,-Nz:Nz))
allocate (fi3d1(-N:N,-Ny:Ny,-Nz:Nz))
allocate (fi3d2(-N:N,-Ny:Ny,-Nz:Nz))
allocate (fi3do2(-N:N,-Ny:Ny,-Nz:Nz))
allocate (flhy1(-N:N,-Ny:Ny,-Nz:Nz))
allocate (flhy2(-N:N,-Ny:Ny,-Nz:Nz))
allocate (flhyo1(-N:N,-Ny:Ny,-Nz:Nz))
allocate (flhyo2(-N:N,-Ny:Ny,-Nz:Nz))
allocate (rdy(-N*2:N*2,-Ny*2:Ny*2,-Nz*2:Nz*2))
c*3.4
pi=4*atan(1.0)
wzl=120*4.1356e-12/27211.6
wrl=60*4.1356e-12/27211.6
y0=13500/.05292*0
omega0=4.1357e-6/27.2116
nxs=n
nys=n
nzs=n
dx=xlp/(nxs-1)
dx=125/.05292
dx=dx
dz=350/.05292
v=-200/27211.6
xm(1)=164/5.486e-4
xm(2)=162/5.486e-4
c     xm(2)=163.929/5.486e-4
eha=27211.6




zred1=xm(1)**2/(xm(1))
zred2=xm(2)**2/(xm(2))
add=131
cdd=12*pi*add/xm(1)
cdd2=cdd
cdd12=cdd

c     do b=0,200,.1
    b=25
    delta1=1.9
    delta2=0.14
    b01=21.91
    b02=26.902
    a162=220*(1-delta1/(b-b01)-delta2/(b-b02))
    c%https://arxiv.org/pdf/1803.10676
    delta1=30.9
    b01=76.9
    delta2=23.6
    b02=178.8
    a164=91*(1-delta1/(b-b01)-delta2/(b-b02))
    c https://arxiv.org/pdf/1506.01875
    delta12=1
    abg6264=105
    b01=10.8
    c rozsadne miedzy 80 a 130
    a6264=abg6264*(1-delta12/(b-b01))
    ggp11=4*pi*a164/xm(1)
    ggp22=4*pi*a162/xm(2)
    ggp12=4*pi*a6264/(xm(1)*xm(2)*2/(xm(1)+xm(2)))
    ggp21=ggp12
    edd2=add/a162
    edd=add/a164
    edd12=add/a6264
    write(9,989) b,a162,a164,a6264,ggp11,ggp22,ggp12,edd,edd2,edd12
    c     enddo
c     stop
write(*,*) a,a162,a12,ggp11,ggp22,ggp12
write(*,*) cdd2,cdd12


do ix=-2*N,2*N
    do iy=-2*Ny,2*Ny
        do iz=-2*Nz,2*Nz
            x=ix*dx
            y=iy*dx
            z=iz*dz
            if(ix**2+iy**2+iz**2.gt.0) then
                rdijk=1/sqrt((x-xp)**2+(y-yp)**2+(z-zp)**2)/4/pi
                c ! uwaga: rdy sluzy do liczenia warunku brzegowego
                c calkowania gestosci, ktora na brzegi jest 0 tak czy inaczej
            endif
            rdy(ix,iy,iz)=rdijk
        enddo
    enddo
enddo

iomega=0
omega=iomega*omega0

psi=0
xnorma=0
iczytaj=1
if(iczytaj.eq.0) then 
    rrr=n*dx
    do ix=-n+1,n-1
        do iy=-ny+1,ny-1
            do iz=-nz+1,nz-1
                x=ix*dx
                y=iy*dx
                z=iz*dz
                cpsi1(ix,iy,iz)=cos(Pi*x/rrr)*cos(Pi*y/rrr)*cos(Pi*z/rrr)
                >*(1+x*x+x*y)
                cpsi2(ix,iy,iz)=cos(Pi*x/rrr)*cos(Pi*y/rrr)*cos(Pi*z/rrr)
                >*(1+y*y+x*y)
                xnorma(1)=xnorma(1)+cdabs(cpsi1(ix,iy,iz))**2*dx**2*dz
                xnorma(2)=xnorma(2)+cdabs(cpsi2(ix,iy,iz))**2*dx**2*dz
            enddo
        enddo
    enddo
    cpsin1=cpsi1
    cpsin2=cpsi2
else
    cpsi1=0
    cpsi2=0 
    read(123,*) i1,i2,i3
    nl=(2*i1+1)*(2*i2+1)*(2*i3+1) 
    allocate (dmin(nl,7))
    do ix=1,nl
        read(123,*) x,y,z,fr,fi,gr,gi
        dmin(ix,1)=x
        dmin(ix,2)=y
        dmin(ix,3)=z
        dmin(ix,4)=fr
        dmin(ix,5)=fi
        dmin(ix,6)=gr
        dmin(ix,7)=gi
    enddo
    do ix=-n+1,n-1
        do iy=-ny+1,ny-1
            do iz=-nz+1,nz-1
                xs=ix*dx
                ys=iy*dx
                zs=iz*dz
                rmin=1e29
                do i=1,nl
                    rc=(xs-dmin(i,1))**2+(ys-dmin(i,2))**2+(zs-dmin(i,3))**2
                    if(rc.lt.rmin) then
                        rmin=rc
                        imin=i
                        if(rc.lt.dx**2/2) goto 343
                    endif 
                enddo 
                343   continue
                cpsi1(ix,iy,iz)=dmin(imin,4)+ci*dmin(imin,5)
                cpsi2(ix,iy,iz)=dmin(imin,6)+ci*dmin(imin,7)
                xnorma(1)=xnorma(1)+cdabs(cpsi1(ix,iy,iz))**2*dx**2*dz
                xnorma(2)=xnorma(2)+cdabs(cpsi2(ix,iy,iz))**2*dx**2*dz
                c     write(125,989) xs,ys,zs,cpsi1(ix,iy,iz),cpsi2(ix,iy,iz),xnorma
            enddo 
        enddo 
    enddo
    cpsi1=cpsi1/dsqrt(xnorma(1))
    cpsi2=cpsi2/dsqrt(xnorma(2))
    cpsin1=cpsi1
    cpsin2=cpsi2
    deallocate(dmin)
endif

xl=5000/.05292 
imax=10000
eold=1e6
dt=2e10
t=0
bb=0.5*xm(1)*wrl**2
dd=2500/.05292
aa=xm(1)*wrl**2/4/dd**2
do ix=-n,n
    x=ix*dx
    vx(ix)=bb*(x-0*dx)**2
    write(11,*) x*.05292,vx(ix)*27211.6
enddo
do iy=-ny,ny
    y=iy*dx
    vy(iy)=0.5*xm(1)*(y-0*dx)**2*wrl**2
enddo
do iz=-nz,nz
    z=iz*dz
    vz(iz)=0.5*xm(1)*z**2*wzl**2
enddo
i192=192
i498=1498
istart=-10100
write(*,*) xncz(1),xncz(2)
c     do iter=0,istart
do iter=0,imax
if(iter.lt.istart) then
callcnd
>(cpsii1,cpsin1,cpsi1,cpsii2,cpsin2,cpsi2,
>n,ny,nz,vx,vy,vz,xlp,xm,dt,fi3d1,fi3do1,fi3d2,fi3do2,
>flhy1,flhyo1,flhy2,flhyo2,rdy)
else

i192=197
i498=498
write(9,989) b,a164,a162,a6264,ggp11,ggp22,ggp12,edd,edd2,edd12

callcndt
>(cpsii1,cpsin1,cpsi1,cpsii2,cpsin2,cpsi2,
>n,ny,nz,vx,vy,vz,xlp,xm,dt,fi3d1,fi3do1,fi3d2,fi3do2,
>flhy1,flhyo1,flhy2,flhyo2,rdy)
endif
t=t+dt
call norm(cpsin1,cpsin2,n,ny,nz)
cpsi1=cpsin1
cpsi2=cpsin2
w1=xncz(1)/xnorma(1)
w2=xncz(2)/xnorma(2)
if(mod(iter,10).eq.0) then
    wredna=energiacnd(cpsi1,cpsi2
    >,n,ny,nz,vx,vy,vz,xlp,xm,dt,fi3d1,fi3d2,
    >flhy1,flhy2,rdy)
    w1=xncz(1)/xnorma(1)
    w2=xncz(2)/xnorma(2)
    eold=wredna
    do ix=-n,n
        x=ix*dx
        do iy=-ny,ny
            x=ix*dx
            y=iy*dx
            write(i498,988)x*.05292,y*.05292,
            >cdabs(cpsi1(ix,iy,0))**2*1e10,
            >cdabs(cpsi2(ix,iy,0))**2*1e10,
            >fi3d1(ix,iy,0)*27211.6,
            >fi3d2(ix,iy,0)*27211.6,
            >flhy1(ix,iy,0)*27211.6,
            >flhy2(ix,iy,0)*27211.6
        enddo
    enddo
    c     stop
    write(i192,88) iter*1.,wredna*eha,e1*eha,e2*eha,
    >e3*eha,e4*eha,e5*eha,e6*eha,e7*eha,e8*eha,e9*eha,dt,
    >t,
    >b,edd,edd2,edd12,ggp11,ggp22,ggp12,
    >a164,a162,a6264
    write(i498,*)
    write(501,*)
    write(499,*)
endif
if(mod(iter,10).eq.0) then
    if(iter.le.istart) then
        open(124,file='ff.dat',action='write')
    else
        open(124,file='fr.dat',action='write')
    endif
    write(124,*) n,ny,nz
    do i=-n,n
        do j=-ny,ny
            do k=-nz,nz
                x=i*dx
                y=j*dx
                z=k*dz
                write(124,989) x,y,z,cpsi1(i,j,k),cpsi2(i,j,k)
            enddo
        enddo
    enddo
    close(124)
endif
enddo
88    format(30f50.32)
988   format(30f30.12)
989   format(30g30.18)
89    format(5g17.8)
end

subroutine norm(cpsi1,cpsi2,n,ny,nz)
implicit double precision(a,b,d-h,o-z)
implicit double complex(c)
dimension cpsi1(-N:N,-Ny:Ny,-Nz:Nz)
dimension cpsi2(-N:N,-Ny:Ny,-Nz:Nz)
dimension xnorma(2),xm(2),ene(2),v(2),xncz(2),xnormak(2),wd(2)
common/loc/xnormak,ene,eod11,eod22,eod12,iix,iiy,iiz
common/xnorma/xnorma,wrl,wzl,dx,xncz,ggp11,dz,ggp12,wd,y0,omega,
>ggp22,ggp21
c     g=9.109e-31*9.8*(0.05292*1e-9)*6.242e18/27.2116
g=.1083e-21*0
xnorma=0

ene=0
do ix=-n+1,n-1
    do iy=-ny+1,ny-1
        do iz=-nz+1,nz-1
            xnorma(1)=xnorma(1)+cdabs(cpsi1(ix,iy,iz))**2*dx**2*dz
            xnorma(2)=xnorma(2)+cdabs(cpsi2(ix,iy,iz))**2*dx**2*dz
        enddo
    enddo
enddo

do ix=-n+1,n-1
    do iy=-ny+1,ny-1
        do iz=-nz+1,nz-1
            cpsi1(ix,iy,iz)=cpsi1(ix,iy,iz)/sqrt(xnorma(1))
            cpsi2(ix,iy,iz)=cpsi2(ix,iy,iz)/sqrt(xnorma(2))
        enddo
    enddo
enddo
end      


subroutine cnd(cpsii1,cpsin1,cpsi1,cpsii2,cpsin2,cpsi2,
>n,ny,nz,vx,vy,vz
>,xlp,xm,dt,fi3d1,fi3do1,fi3d2,fi3do2,flhy1,flhyo1,flhy2,flhyo2,
>rdy)
implicit double precision(a,b,d-h,o-z)
implicit double complex(c)
dimension cpsi1(-N:N,-Ny:Ny,-Nz:Nz)
dimension cpsin1(-N:N,-Ny:Ny,-Nz:Nz)
dimension cpsii1(-N:N,-Ny:Ny,-Nz:Nz)
dimension cpsi2(-N:N,-Ny:Ny,-Nz:Nz)
dimension cpsin2(-N:N,-Ny:Ny,-Nz:Nz)
dimension cpsii2(-N:N,-Ny:Ny,-Nz:Nz)
dimension vdd(-N:N,-Ny:Ny,-Nz:Nz)
dimension fi3d1(-N:N,-Ny:Ny,-Nz:Nz)
dimension fi3d2(-N:N,-Ny:Ny,-Nz:Nz)
dimension flhy1(-N:N,-Ny:Ny,-Nz:Nz)
dimension flhy2(-N:N,-Ny:Ny,-Nz:Nz)
dimension flhyo1(-N:N,-Ny:Ny,-Nz:Nz)
dimension flhyo2(-N:N,-Ny:Ny,-Nz:Nz)
dimension fi3do1(-N:N,-Ny:Ny,-Nz:Nz)
dimension fi3do2(-N:N,-Ny:Ny,-Nz:Nz)
dimension rdy(-2*N:2*N,-2*Ny:2*Ny,-2*Nz:2*Nz)
dimension xnorma(2),xm(2),ene(2),v(2),xncz(2),xnormak(2),wd(2)
dimension vx(-n:n),vy(-ny:ny),vz(-nz:nz)
common/loc/xnormak,ene,eod11,eod22,eod12,iix,iiy,iiz
common/xnorma/xnorma,wrl,wzl,dx,xncz,ggp11,dz,ggp12,wd,y0,omega,
>ggp22,ggp21
common/add/add,cdd,gamma1,gamma2,cdd2,cdd12,edd,edd2,edd12
xnorma=0
ci=(0.d0,1.d0)
ene=0
cpsin1=cpsi1
cpsii1=cpsi1
cpsin2=cpsi2
cpsii2=cpsi2
cdt=dt*(-ci)

xnorma=0
do 131 ix=-n+1,n-1
do 131 iy=-ny+1,ny-1
do 131 iz=-nz+1,nz-1
xnorma(1)=xnorma(1)+cdabs(cpsi1(ix,iy,iz))**2*dx**2*dz
xnorma(2)=xnorma(2)+cdabs(cpsi2(ix,iy,iz))**2*dx**2*dz
131   continue

call lifi3(fi3d1,fi3d2,flhy1,flhy2,
>cpsi1,cpsi2,n,ny,nz,dx,dz,rdy,xm,ggp11,ggp22,ggp12,ggp21)
fi3do1=fi3d1
fi3do2=fi3d2
flhyo1=flhy1
flhyo2=flhy2

do 10 icn=1,1
ene=0
w1=(xncz(1))/xnorma(1)
w2=(xncz(2))/xnorma(2)
do ix=-n+1,n-1
do iy=-ny+1,ny-1
do iz=-nz+1,nz-1
    x=ix*dx
    y=iy*dx-y0
    z=iz*dz
    vvv=vx(ix)+vy(iy)+vz(iz)
    c1=-0.5/xm(1)/dx**2*(
    >cpsi1(ix-1,iy,iz)+cpsi1(ix+1,iy,iz)+
    >cpsi1(ix,iy-1,iz)+cpsi1(ix,iy+1,iz)-4*cpsi1(ix,iy,iz))
    >   -0.5/xm(1)/dz**2*(
    >cpsi1(ix,iy,iz-1)+cpsi1(ix,iy,iz+1)
    >-2*cpsi1(ix,iy,iz))
    >+cpsi1(ix,iy,iz)*(vvv+fi3do1(ix,iy,iz)+flhyo1(ix,iy,iz)
    > +ggp11*cdabs(cpsi1(ix,iy,iz))**2*w1
    >+ ggp12*cdabs(cpsi2(ix,iy,iz))**2*w2)
    cpsii1(ix,iy,iz)=cpsi1(ix,iy,iz)
    >+cdt/ci*c1
    c2=-0.5/xm(2)/dx**2*(
    >cpsi2(ix-1,iy,iz)+cpsi2(ix+1,iy,iz)+
    >cpsi2(ix,iy-1,iz)+cpsi2(ix,iy+1,iz)-4*cpsi2(ix,iy,iz))
    >   -0.5/xm(2)/dz**2*(
    >cpsi2(ix,iy,iz-1)+cpsi2(ix,iy,iz+1)
    >-2*cpsi2(ix,iy,iz))
    >+cpsi2(ix,iy,iz)*(vvv+fi3do2(ix,iy,iz)+flhyo2(ix,iy,iz)
    > +ggp22*cdabs(cpsi2(ix,iy,iz))**2*w2
    >+ ggp21*cdabs(cpsi1(ix,iy,iz))**2*w1)
    cpsii2(ix,iy,iz)=cpsi2(ix,iy,iz)
    >+cdt/ci*c2
    enddo
    enddo
enddo
cpsin1=cpsii1
cpsin2=cpsii2
10    continue
88    format(30g30.12)      
end      


subroutine cndt(cpsii1,cpsin1,cpsi1,cpsii2,cpsin2,cpsi2,
>n,ny,nz,vx,vy,vz
>,xlp,xm,dt,fi3d1,fi3do1,fi3d2,fi3do2,flhy1,flhyo1,flhy2,flhyo2,
>rdy)
implicit double precision(a,b,d-h,o-z)
implicit double complex(c)
dimension cpsi1(-N:N,-Ny:Ny,-Nz:Nz)
dimension cpsin1(-N:N,-Ny:Ny,-Nz:Nz)
dimension cpsii1(-N:N,-Ny:Ny,-Nz:Nz)
dimension cpsi2(-N:N,-Ny:Ny,-Nz:Nz)
dimension cpsin2(-N:N,-Ny:Ny,-Nz:Nz)
dimension cpsii2(-N:N,-Ny:Ny,-Nz:Nz)
dimension vdd(-N:N,-Ny:Ny,-Nz:Nz)
dimension fi3d1(-N:N,-Ny:Ny,-Nz:Nz)
dimension fi3d2(-N:N,-Ny:Ny,-Nz:Nz)
dimension flhy1(-N:N,-Ny:Ny,-Nz:Nz)
dimension flhy2(-N:N,-Ny:Ny,-Nz:Nz)
dimension flhyo1(-N:N,-Ny:Ny,-Nz:Nz)
dimension flhyo2(-N:N,-Ny:Ny,-Nz:Nz)
dimension fi3do1(-N:N,-Ny:Ny,-Nz:Nz)
dimension fi3do2(-N:N,-Ny:Ny,-Nz:Nz)
dimension rdy(-2*N:2*N,-2*Ny:2*Ny,-2*Nz:2*Nz)
dimension xnorma(2),xm(2),ene(2),v(2),xncz(2),xnormak(2),wd(2)
dimension vx(-n:n),vy(-ny:ny),vz(-nz:nz)
common/loc/xnormak,ene,eod11,eod22,eod12,iix,iiy,iiz
common/xnorma/xnorma,wrl,wzl,dx,xncz,ggp11,dz,ggp12,wd,y0,omega,
>ggp22,ggp21
common/add/add,cdd,gamma1,gamma2,cdd2,cdd12,edd,edd2,edd12
xnorma=0
ci=(0.d0,1.d0)
ene=0
cpsin1=cpsi1
cpsii1=cpsi1
cpsin2=cpsi2
cpsii2=cpsi2
cdt=dt

xnorma=0
do 131 ix=-n+1,n-1
do 131 iy=-ny+1,ny-1
do 131 iz=-nz+1,nz-1
xnorma(1)=xnorma(1)+cdabs(cpsi1(ix,iy,iz))**2*dx**2*dz
xnorma(2)=xnorma(2)+cdabs(cpsi2(ix,iy,iz))**2*dx**2*dz
131   continue

call lifi3(fi3d1,fi3d2,flhy1,flhy2,
>cpsi1,cpsi2,n,ny,nz,dx,dz,rdy,xm,ggp11,ggp22,ggp12,ggp21)
fi3do1=fi3d1
fi3do2=fi3d2
flhyo1=flhy1
flhyo2=flhy2

do 10 icn=1,3
ene=0
w1=(xncz(1))/xnorma(1)
w2=(xncz(2))/xnorma(2)
do ix=-n+1,n-1
do iy=-ny+1,ny-1
do iz=-nz+1,nz-1
    x=ix*dx
    y=iy*dx-y0
    z=iz*dz
    vvv=vx(ix)+vy(iy)+vz(iz)
    c1=-0.5/xm(1)/dx**2*(
    >cpsi1(ix-1,iy,iz)+cpsi1(ix+1,iy,iz)+
    >cpsi1(ix,iy-1,iz)+cpsi1(ix,iy+1,iz)-4*cpsi1(ix,iy,iz))
    >   -0.5/xm(1)/dz**2*(
    >cpsi1(ix,iy,iz-1)+cpsi1(ix,iy,iz+1)
    >-2*cpsi1(ix,iy,iz))
    >+cpsi1(ix,iy,iz)*(vvv+fi3do1(ix,iy,iz)+flhyo1(ix,iy,iz)
    > +ggp11*cdabs(cpsi1(ix,iy,iz))**2*w1
    >+ ggp12*cdabs(cpsi2(ix,iy,iz))**2*w2)
    c1n=-0.5/xm(1)/dx**2*(
    >cpsin1(ix-1,iy,iz)+cpsin1(ix+1,iy,iz)+
    >cpsin1(ix,iy-1,iz)+cpsin1(ix,iy+1,iz)-4*cpsin1(ix,iy,iz))
    >   -0.5/xm(1)/dz**2*(
    >cpsin1(ix,iy,iz-1)+cpsin1(ix,iy,iz+1)
    >-2*cpsin1(ix,iy,iz))
    >+cpsin1(ix,iy,iz)*(vvv+fi3d1(ix,iy,iz)+flhy1(ix,iy,iz)
    > +ggp11*cdabs(cpsin1(ix,iy,iz))**2*w1
    >+ ggp12*cdabs(cpsin2(ix,iy,iz))**2*w2)
    cpsii1(ix,iy,iz)=cpsi1(ix,iy,iz)
    >+cdt/ci*(c1+c1n)/2
    c2=-0.5/xm(2)/dx**2*(
    >cpsi2(ix-1,iy,iz)+cpsi2(ix+1,iy,iz)+
    >cpsi2(ix,iy-1,iz)+cpsi2(ix,iy+1,iz)-4*cpsi2(ix,iy,iz))
    >   -0.5/xm(2)/dz**2*(
    >cpsi2(ix,iy,iz-1)+cpsi2(ix,iy,iz+1)
    >-2*cpsi2(ix,iy,iz))
    >+cpsi2(ix,iy,iz)*(vvv+fi3do2(ix,iy,iz)+flhyo2(ix,iy,iz)
    > +ggp22*cdabs(cpsi2(ix,iy,iz))**2*w2
    >+ ggp21*cdabs(cpsi1(ix,iy,iz))**2*w1)
    c2n=-0.5/xm(2)/dx**2*(
    >cpsin2(ix-1,iy,iz)+cpsin2(ix+1,iy,iz)+
    >cpsin2(ix,iy-1,iz)+cpsin2(ix,iy+1,iz)-4*cpsin2(ix,iy,iz))
    >   -0.5/xm(2)/dz**2*(
    >cpsin2(ix,iy,iz-1)+cpsin2(ix,iy,iz+1)
    >-2*cpsin2(ix,iy,iz))
    >+cpsin2(ix,iy,iz)*(vvv+fi3d2(ix,iy,iz)+flhy2(ix,iy,iz)
    > +ggp22*cdabs(cpsin2(ix,iy,iz))**2*w2
    >+ ggp21*cdabs(cpsin1(ix,iy,iz))**2*w1)
    cpsii2(ix,iy,iz)=cpsi2(ix,iy,iz)
    >+cdt/ci*(c2+c2n)/2
    enddo
    enddo
enddo
cpsin1=cpsii1
cpsin2=cpsii2
call lifi3(fi3d1,fi3d2,flhy1,flhy2,
>cpsii1,cpsii2,n,ny,nz,dx,dz,rdy,xm,ggp11,ggp22,ggp12,ggp21)
10    continue
88    format(30g30.12)      
end      


subroutine lifi3(fi3d1,fi3d2,
>flhy1,flhy2,cpsi1,cpsi2,n,ny,nz,dx,dz,
>rdy,xm,ggp11,ggp22,ggp12,ggp21)
implicit double precision(a,b,d-h,o-z)
implicit double complex (c)
dimension fi3d1(-N:N,-Ny:Ny,-Nz:Nz),xm(2)
dimension fi3d2(-N:N,-Ny:Ny,-Nz:Nz)
dimension flhy1(-N:N,-Ny:Ny,-Nz:Nz)
dimension flhy2(-N:N,-Ny:Ny,-Nz:Nz)
dimension rdy(-2*N:2*N,-2*Ny:2*Ny,-2*Nz:2*Nz)
dimension cpsi1(-N:N,-Ny:Ny,-Nz:Nz),xnorma(2),xncz(2)
dimension cpsi2(-N:N,-Ny:Ny,-Nz:Nz)
dimension wg(16),ug(16)
allocatable fi1(:,:,:)
allocatable fi2(:,:,:)
allocatable psi1(:,:,:)
allocatable psi2(:,:,:)
common/add/add,cdd,gamma1,gamma2,cdd2,cdd12,edd,edd2,edd12
common/xnorma/xnorma,wr,wz,ddx,xncz,gp11,ddz,
>gp12,wwd,yy0,oomega,gp22,gp21
common/gauss/wg,ug
pi=4*atan(1.0)
allocate(fi1(-N:N,-Ny:Ny,-Nz:Nz))
allocate(fi2(-N:N,-Ny:Ny,-Nz:Nz))
allocate(psi1(-N:N,-Ny:Ny,-Nz:Nz))
allocate(psi2(-N:N,-Ny:Ny,-Nz:Nz))
c     fi3d=0
c     return
c     xfct=dx**2*dz*xncz/xnorma
xfct=dx**2*dz
w1=xncz(1)/xnorma(1)
w2=xncz(2)/xnorma(2)
c ############### wyliczamy najpierw fi3d1 -- na gaz 1, a potem fi3d2 -- na gaz2
do ix=-N,N
do iy=-Ny,Ny
do iz=-Nz,Nz
    psi1(ix,iy,iz)=xncz(1)/xnorma(1)*cdabs(cpsi1(ix,iy,iz))**2*cdd
    >+xncz(2)/xnorma(2)*cdabs(cpsi2(ix,iy,iz))**2*cdd12
    psi2(ix,iy,iz)=xncz(1)/xnorma(1)*cdabs(cpsi1(ix,iy,iz))**2*cdd12
    >+xncz(2)/xnorma(2)*cdabs(cpsi2(ix,iy,iz))**2*cdd2
    enddo
    enddo
enddo

do ix=-N,N,2*N
    do iy=-Ny,Ny
        do iz=-Nz,Nz
            fi1(ix,iy,iz)=0
            fi2(ix,iy,iz)=0
        enddo
    enddo
enddo

do ix=-N,N
    do iy=-Ny,Ny,2*Ny
        do iz=-Nz,Nz
            fi1(ix,iy,iz)=0
            fi2(ix,iy,iz)=0
        enddo
    enddo
enddo

do ix=-N,N
    do iy=-Ny,Ny
        do iz=-Nz,Nz,2*nz
            fi1(ix,iy,iz)=0
            fi2(ix,iy,iz)=0
        enddo
    enddo
enddo

do ix=-N,N,2*N
    x=ix*dx
    do iy=-Ny,Ny
        do iz=-Nz,Nz
            y=iy*dx
            z=iz*dz
            do ixp=-N+1,N-1
                do iyp=-Ny+1,Ny-1
                    do izp=-Nz+1,Nz-1
                        xp=ixp*dx
                        yp=iyp*dx
                        zp=izp*dz  
                        rdijk=rdy(ix-ixp,iy-iyp,iz-izp)
                        fi1(ix,iy,iz)=fi1(ix,iy,iz)+
                        >rdijk*psi1(ixp,iyp,izp)*xfct
                        fi2(ix,iy,iz)=fi2(ix,iy,iz)+
                        >rdijk*psi2(ixp,iyp,izp)*xfct
                    enddo
                enddo
            enddo
        enddo
    enddo
enddo

do iy=-Ny,Ny,2*Ny
    y=iy*dx
    do ix=-N,N
        do iz=-Nz,Nz
            x=ix*dx
            z=iz*dz
            do ixp=-N+1,N-1
                do iyp=-Ny+1,Ny-1
                    do izp=-Nz+1,Nz-1
                        xp=ixp*dx
                        yp=iyp*dx
                        zp=izp*dz  
                        rdijk=rdy(ix-ixp,iy-iyp,iz-izp)
                        fi1(ix,iy,iz)=fi1(ix,iy,iz)+
                        >rdijk*psi1(ixp,iyp,izp)*xfct
                        fi2(ix,iy,iz)=fi2(ix,iy,iz)+
                        >rdijk*psi2(ixp,iyp,izp)*xfct
                    enddo
                enddo
            enddo

        enddo
    enddo
enddo

do iz=-Nz,Nz,2*Nz
    z=iz*dz
    do ix=-N,N
        do iy=-Ny,Ny
            x=ix*dx
            y=iy*dx
            do ixp=-N+1,N-1
                do iyp=-Ny+1,Ny-1
                    do izp=-Nz+1,Nz-1
                        xp=ixp*dx
                        yp=iyp*dx
                        zp=izp*dz  
                        rdijk=rdy(ix-ixp,iy-iyp,iz-izp)
                        fi1(ix,iy,iz)=fi1(ix,iy,iz)+
                        >rdijk*psi1(ixp,iyp,izp)*xfct
                        fi2(ix,iy,iz)=fi2(ix,iy,iz)+
                        >rdijk*psi2(ixp,iyp,izp)*xfct
                    enddo
                enddo
            enddo

        enddo
    enddo
enddo

c rownanie poissona: nabla^2\phi=-rho
c dla rho=delta, phi=1/r/4/pi

c     st=system("date")

omega=1.92
do iter=1,1200
do ix=-N+1,N-1
do iy=-Ny+1,Ny-1
    do iz=-Nz+1,Nz-1
        fi1(ix,iy,iz)=fi1(ix,iy,iz)*(1-omega)+
        >omega*
        >(
        >(fi1(ix,iy+1,iz)+fi1(ix,iy-1,iz)
        >+fi1(ix+1,iy,iz)+fi1(ix-1,iy,iz))*dz**2
        >+(fi1(ix,iy,iz+1)+fi1(ix,iy,iz-1))*dx**2
        >+psi1(ix,iy,iz)*dx**2*dz**2
        >)
        >/(4*dz**2+2*dx**2)
        fi2(ix,iy,iz)=fi2(ix,iy,iz)*(1-omega)+
        >omega*
        >(
        >(fi2(ix,iy+1,iz)+fi2(ix,iy-1,iz)
        >+fi2(ix+1,iy,iz)+fi2(ix-1,iy,iz))*dz**2
        >+(fi2(ix,iy,iz+1)+fi2(ix,iy,iz-1))*dx**2
        >+psi2(ix,iy,iz)*dx**2*dz**2
        >)
        >/(4*dz**2+2*dx**2)
        enddo
        enddo
    enddo

enddo

fi3d1=0
fi3d2=0
do ix=-N+1,N-1
    do iy=-Ny+1,Ny-1
        do iz=-Nz+1,Nz-1
            fi3d1(ix,iy,iz)=
            >-
            >(fi1(ix,iy,iz+1)+fi1(ix,iy,iz-1)-2*fi1(ix,iy,iz))/dz**2
            >-1./3*psi1(ix,iy,iz)
            fi3d2(ix,iy,iz)=
            >-
            >(fi2(ix,iy,iz+1)+fi2(ix,iy,iz-1)-2*fi2(ix,iy,iz))/dz**2
            >-1./3*psi2(ix,iy,iz)
            rn1=cdabs(cpsi1(ix,iy,iz))**2*w1
            xmno1=4./3/pi**2*xm(1)**1.5
            xmno2=4./3/pi**2*xm(2)**1.5
            rn2=cdabs(cpsi2(ix,iy,iz))**2*w2
            wlhy1=0
            wlhy2=0
            do ig=1,8
                aggp11=ggp11+cdd/3*(3*ug(ig)**2-1)
                aggp22=ggp22+cdd2/3*(3*ug(ig)**2-1)
                aggp12=ggp12+cdd12/3*(3*ug(ig)**2-1)
                aggp21=aggp12
                D=sqrt((aggp11*rn1-aggp22*rn2)**2
                >+4*aggp12**2*rn1*rn2)
                if(D.gt.1d-41) then
                    xlp=0.5*
                    >(aggp11*rn1+aggp22*rn2+D)
                    cxlm=0.5*
                    >(aggp11*rn1+aggp22*rn2-D)
                    xlppn1=0.5*
                    >(aggp11+(aggp11*(aggp11*rn1-aggp22*rn2)+2*aggp12**2*rn2)/D)
                    xlmpn1=0.5*
                    >(aggp11-(aggp11*(aggp11*rn1-aggp22*rn2)+2*aggp12**2*rn2)/D)
                    xlppn2=0.5*
                    >(aggp22+(-aggp22*(aggp11*rn1-aggp22*rn2)+2*aggp12**2*rn1)/D)
                    xlmpn2=0.5*
                    >(aggp22-(-aggp22*(aggp11*rn1-aggp22*rn2)+2*aggp12**2*rn1)/D)
                    slp=xlp**1.5
                    slm=cxlm**1.5
                    wlhy1=wlhy1+xmno1*real(slp*xlppn1+slm*xlmpn1)*wg(ig)
                    wlhy2=wlhy2+xmno2*real(slp*xlppn2+slm*xlmpn2)*wg(ig)
                endif
            enddo
            flhy1(ix,iy,iz)=wlhy1
            flhy2(ix,iy,iz)=wlhy2
            c     if((.not.wlhy1.gt.0.and..not.wlhy1.le.0).or.
                c    >(.not.wlhy2.gt.0.and..not.wlhy1.le.0))
                c    > then
                c     write(*,*) rn1,rn2,wlhy1,wlhy2,xlp,xlppn1,cxml,xlmpn1,xlppn2,
                c    >xlmpn2,D,
                c    >((aggp11*rn1-aggp22*rn2)**2+4*aggp12**2*rn1*rn2)
                c     write(*,*) 'rn1,rn2,wlhy1,wlhy2,xlp,xlppn1,cxml,xlmpn1,
                c    >xlppn2,xlmpn2,
                c    >D,arg'
                c     stop
            c     endif
        enddo
    enddo
enddo

end


function energiacnd
>(cpsi1,cpsi2,n,ny,nz,vx,vy,vz,xlp,xm,dt,
>fi3d1,fi3d2,flhy1,flhy2,rdy)
implicit double precision(a,b,d-h,o-z)
implicit double complex(c)
dimension cpsi1 (-N:N,-Ny:Ny,-Nz:Nz)
dimension cpsi2 (-N:N,-Ny:Ny,-Nz:Nz)
dimension rdy(-2*N:2*N,-2*Ny:2*Ny,-2*Nz:2*Nz)
dimension vdd(-N:N,-Ny:Ny,-Nz:Nz)
dimension fi3d1(-N:N,-Ny:Ny,-Nz:Nz)
dimension fi3d2(-N:N,-Ny:Ny,-Nz:Nz)
dimension flhy1(-N:N,-Ny:Ny,-Nz:Nz)
dimension flhy2(-N:N,-Ny:Ny,-Nz:Nz),wg(16),ug(16)
dimension xnorma(2),xm(2),ene(2),v(2),xncz(2),xnormak(2),wd(2)
dimension vx(-n:n),vy(-ny:ny),vz(-nz:nz)
common/loc/xnormak,ene,eod11,eod22,eod12,iix,iiy,iiz
common/xnorma/xnorma,wrl,wzl,dx,xncz,ggp11,dz,ggp12,wd,y0,omega,
>ggp22,ggp21
common/add/add,cdd,gamma1,gamma2,cdd2,cdd12,edd,edd2,edd12
common/skladowe/e1,e2,e3,e4,e5,e6,e7,e8,e9
common/gauss/wg,ug
xnorma=0
pi=4*atan(1.0)
ci=(0.d0,1.d0)
ene=0
cdt=dt*(-ci)

xnorma=0
do 131 ix=-n+1,n-1
do 131 iy=-ny+1,ny-1
do 131 iz=-nz+1,nz-1
xnorma(1)=xnorma(1)+cdabs(cpsi1(ix,iy,iz))**2*dx**2*dz
xnorma(2)=xnorma(2)+cdabs(cpsi2(ix,iy,iz))**2*dx**2*dz
131   continue

e1=0
e2=0
e3=0
e4=0
e5=0
e6=0
e7=0
e8=0
e9=0
w1=xncz(1)/xnorma(1)
w2=xncz(2)/xnorma(2)
do ix=-n+1,n-1
do iy=-ny+1,ny-1
do iz=-nz+1,nz-1
    rn1=cdabs(cpsi1(ix,iy,iz))**2*w1
    rn2=cdabs(cpsi2(ix,iy,iz))**2*w2

    elhy1=0
    elhy2=0
    do ig=1,8
        aggp11=ggp11+cdd/3*(3*ug(ig)**2-1)
        aggp22=ggp22+cdd2/3*(3*ug(ig)**2-1)
        aggp12=ggp12+cdd12/3*(3*ug(ig)**2-1)
        aggp21=aggp12
        D=sqrt((aggp11*rn1-aggp22*rn2)**2
        >+4*aggp12**2*rn1*rn2)
        xlp=0.5*
        >(aggp11*rn1+aggp22*rn2+D)
        cxlm=0.5*
        >(aggp11*rn1+aggp22*rn2-D)
        slm=cxml**2.5
        slp=xlp**2.5
        elhy1=8./15/pi**2*xm(1)**1.5*slm
        >*wg(ig)+elhy1
        elhy2=8./15/pi**2*xm(2)**1.5*slp
        >*wg(ig)+elhy2
        enddo
            vvv=vx(ix)+vy(iy)+vz(iz)
            c1=-0.5/xm(1)/dx**2*(
            >cpsi1(ix-1,iy,iz)+cpsi1(ix+1,iy,iz)+
            >cpsi1(ix,iy-1,iz)+cpsi1(ix,iy+1,iz)-4*cpsi1(ix,iy,iz))
            >   -0.5/xm(1)/dz**2*(
            >cpsi1(ix,iy,iz-1)+cpsi1(ix,iy,iz+1)
            >-2*cpsi1(ix,iy,iz))
            c2=
            >+cpsi1(ix,iy,iz)*(vvv)
            c3=
            >+cpsi1(ix,iy,iz)*fi3d1(ix,iy,iz)/2
            c4=ggp11*w1*cdabs(cpsi1(ix,iy,iz))**2/2
            >+ggp12*w2*cdabs(cpsi2(ix,iy,iz))**2/2
            c5=elhy1+elhy2
            e1=e1+c1*dconjg(cpsi1(ix,iy,iz))*dx**2*xncz(1)/xnorma(1)*dz
            e2=e2+c2*dconjg(cpsi1(ix,iy,iz))*dx**2*xncz(1)/xnorma(1)*dz
            e3=e3+c3*dconjg(cpsi1(ix,iy,iz))*dx**2*xncz(1)/xnorma(1)*dz
            e4=e4+c4*cdabs(cpsi1(ix,iy,iz))**2*dx**2*xncz(1)/xnorma(1)*dz
            e9=e9+c5*dx**2*dz

            c1=-0.5/xm(2)/dx**2*(
            >cpsi2(ix-1,iy,iz)+cpsi2(ix+1,iy,iz)+
            >cpsi2(ix,iy-1,iz)+cpsi2(ix,iy+1,iz)-4*cpsi2(ix,iy,iz))
            >   -0.5/xm(2)/dz**2*(
            >cpsi2(ix,iy,iz-1)+cpsi2(ix,iy,iz+1)
            >-2*cpsi2(ix,iy,iz))
            c2=
            >+cpsi2(ix,iy,iz)*(vvv)
            c3=
            >+cpsi2(ix,iy,iz)*fi3d2(ix,iy,iz)/2
            c4=ggp22*w2*cdabs(cpsi2(ix,iy,iz))**2/2
            >+ggp21*w1*cdabs(cpsi1(ix,iy,iz))**2/2
            e5=e5+c1*dconjg(cpsi2(ix,iy,iz))*dx**2*xncz(2)/xnorma(2)*dz
            e6=e6+c2*dconjg(cpsi2(ix,iy,iz))*dx**2*xncz(2)/xnorma(2)*dz
            e7=e7+c3*dconjg(cpsi2(ix,iy,iz))*dx**2*xncz(2)/xnorma(2)*dz
            e8=e8+c4*cdabs(cpsi2(ix,iy,iz))**2*dx**2*xncz(2)/xnorma(2)*dz
        enddo
    enddo
enddo

eod22=ene(1)+eod11
energiacnd=e1+e2+e3+e4+e5+e6+e7+e8+e9
88    format(30g30.12)      
end      

