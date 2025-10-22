program gaussian
    implicit double precision(a-b,d-h,o-z)
    implicit double complex(c)
    ! declare variables !
    parameter(N=10,ny=10,nz=3)
    allocatable cpsi(:,:,:)
    allocatable cpsin(:,:,:)
    allocatable cpsii(:,:,:)
    allocatable dmin(:,:)
    allocatable fi3do(:,:,:)
    allocatable fi3d(:,:,:)
    dimension cds(4),f(4)
    dimension xnorma(2),xm(2),xncz(2),xnormak(2),ene(2),wd(2)
    dimension vx(-n:n),vy(-ny:ny),vz(-nz:nz)
    allocatable rdy(:,:,:)
    common/xnorma/xnorma,wrl,wzl,dx,xncz,ggp11,dz,ggp12,wd,y0,omega
    common/loc/xnormak,ene,eod11,eod22,eod12,ix,iy,iz
    common/add/add,cdd,gamma
    common/xl/xl
    common/skladowe/e1,e2,e3,e4,e5

    ! initialize variables !
    ci=(0.d0,1.d0) ! just imaginary double
    xncz(1)=.5e4
    allocate (cpsi(-N:N,-Ny:Ny,-Nz:Nz))   ! wavefunction
    allocate (cpsin(-N:N,-Ny:Ny,-Nz:Nz))  !
    allocate (cpsii(-N:N,-Ny:Ny,-Nz:Nz))  !
    allocate (fi3do(-N:N,-Ny:Ny,-Nz:Nz))
    allocate (fi3d(-N:N,-Ny:Ny,-Nz:Nz))
    allocate (rdy(-N*2:N*2,-Ny*2:Ny*2,-Nz*2:Nz*2))
    pi=4*atan(1.0)
    wzl=120*4.1356e-12/27211.6
    wrl=60*4.1356e-12/27211.6
    y0=13500/.05292*0
    omega0=4.1357e-6/27.2116
    nxs=n
    nys=n
    nzs=n
    dx=60/.05292
    dz=350/.05292
    v=-200/27211.6
    xm(1)=163.929/5.486e-4
    eha=27211.6

    zred1=xm(1)**2/(xm(1))
    add=131
    cdd=12*pi*add/xm(1)
    edd=1.5
    a=add/edd
    ggp11=4*pi*a/zred1
    gamma=128*sqrt(pi)*a**2.5/3/xm(1)*(1+1.5*edd**2)

    do ix=-2*N,2*N
        do iy=-2*Ny,2*Ny
            do iz=-2*Nz,2*Nz
                x=ix*dx
                y=iy*dx
                z=iz*dz
                if(ix**2+iy**2+iz**2.gt.0) then
                    rdijk=1/sqrt((x-xp)**2+(y-yp)**2+(z-zp)**2)/4/pi
! uwaga: rdy sluzy do liczenia warunku brzegowego
! calkowania gestosci, ktora na brzegi jest 0 tak czy inaczej
                endif
                rdy(ix,iy,iz)=rdijk
            enddo
        enddo
    enddo

    iomega=0
    omega=iomega*omega0

    psi=0
    xnorma=0
    iczytaj=0
    if(iczytaj.eq.0) then
        rrr=n*dx
        do ix=-n+1,n-1
            do iy=-ny+1,ny-1
                do iz=-nz+1,nz-1
                    x=ix*dx
                    y=iy*dx
                    z=iz*dz
                    cpsi(ix,iy,iz)=cos(2*Pi*x/rrr)*cos(2*Pi*y/rrr)*cos(Pi*z/rrr)
                    xnorma(1)=xnorma(1)+cdabs(cpsi(ix,iy,iz))**2*dx**2*dz
                enddo
            enddo
        enddo
    else
        read(123,*) i1,i2,i3
        nl=(2*i1+1)*(2*i2+1)*(2*i3+1)
        allocate (dmin(nl,5))
        do ix=1,nl
            read(123,*) x,y,z,fr,fi
            dmin(ix,1)=x
            dmin(ix,2)=y
            dmin(ix,3)=z
            dmin(ix,4)=fr
            dmin(ix,5)=fi
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
343                 continue
                    cpsi(ix,iy,iz)=dmin(imin,4)+ci*dmin(imin,5)
                    xnorma=xnorma+cdabs(cpsi(ix,iym,iz))**2*dx**2*dz
                enddo
            enddo
        enddo
        deallocate(dmin)
    endif

    xl=5000/.05292
    imax=6000*200
    eold=1e6
    dt=1.25e11
    t=0
    b=0.5*xm(1)*wrl**2
    dd=1500/.05292
    aa=xm(1)*wrl**2/4/dd**2
    do ix=-n,n
        x=ix*dx
        vx(ix)=-b*x**2+aa*x**4
        write(11,*) x*.05292,vx(ix)*27211.6
    enddo
    do iy=-ny,ny
        y=iy*dx
        vy(iy)=0.5*xm(1)*y**2*wrl**2
    enddo
    do iz=-nz,nz
        z=iz*dz
        vz(iz)=0.5*xm(1)*z**2*wzl**2
    enddo
    i192=192
    i498=1498
    istart=8000
    do iter=istart+1,imax

        if(iter.le.istart) then
            call cnd(cpsii,cpsin,cpsi,n,ny,nz,vx,vy,vz,xm,dt,fi3d,fi3do,rdy)
        else
            t=t+dt
            i192=197
            i498=498
            t=0

            xl=0
            dt=1.0e10
            b=0
            do ix=-n,n
                x=ix*dx
                vx(ix)=-b*x**2+aa*x**4
            enddo
            call cndt(cpsii,cpsin,cpsi,n,ny,nz,vx,vy,vz,xm,dt,fi3d,fi3do,rdy)
        endif

        call norm(cpsin,n,ny,nz)
        cpsi=cpsin
        if(mod(iter,300).eq.0) then
            wredna=energiacnd(cpsi,n,ny,nz,vx,vy,vz,xm,dt,fi3d,rdy)
            w=xncz(1)/xnorma(1)
            eold=wredna
            do ix=-n,n
                x=ix*dx
                iy=0
                write(499,989)x*.05292,cdabs(cpsi(ix,0,0))**2,                      &
                    dreal(cdd*fi3d(ix,iy,0)-cdd/3*cdabs(cpsi(ix,iy,0))**2*w)*eha,   &
                    (ggp11*cdabs(cpsi(ix,iy,0))**2*w)*eha,                          &
                    (gamma*cdabs(cpsi(ix,iy,0))**3*w**1.5)*eha

                do iy=-ny,ny
                    x=ix*dx
                    y=iy*dx
                    write(i498,989)x*.05292,y*.05292,cdabs(cpsi(ix,iy,0))**2,           &
                        dreal(cdd*fi3d(ix,iy,0)-cdd/3*cdabs(cpsi(ix,iy,0))**2*w)*eha,   &
                        (ggp11*cdabs(cpsi(ix,iy,0))**2*w)*eha,                          &
                        (gamma*cdabs(cpsi(ix,iy,0))**3*w**1.5)*eha
                enddo
            enddo
            do ix=-nz,nz
                y=ix*dz
                write(501,88)y*.05292,cdabs(cpsin(0,0,ix))**2
            enddo
            write(i192,88) iter*1.,wredna*eha,e1*eha,e2*eha,e3*eha,e4*eha,e5*eha,dt,t
            write(i498,*)
            write(501,*)
            write(499,*)
        endif
        if(mod(iter,300).eq.0) then
            open(124,file='ff.dat',action='write')
            write(124,*) n,ny,nz
            do i=-n,n
                do j=-ny,ny
                    do k=-nz,nz
                        x=i*dx
                        y=j*dx
                        z=k*dz
                        write(124,89) x,y,z,cpsi(i,j,k)
                    enddo
                enddo
            enddo
            close(124)
        endif
    enddo
88  format(30f50.32)
988 format(30f30.12)
989 format(30g30.12)
89  format(5g17.8)
end

subroutine norm(cpsi,n,ny,nz)
    implicit double precision(a,b,d-h,o-z)
    implicit double complex(c)
    dimension cpsi(-N:N,-Ny:Ny,-Nz:Nz)
    dimension xnorma(2),xm(2),ene(2),v(2),xncz(2),xnormak(2),wd(2)
    common/loc/xnormak,ene,eod11,eod22,eod12,iix,iiy,iiz
    common/xnorma/xnorma,wrl,wzl,dx,xncz,ggp11,dz,ggp12,wd,y0,omega
    g=.1083e-21*0
    xnorma=0

    ene=0
    do ix=-n+1,n-1
        do iy=-ny+1,ny-1
            do iz=-nz+1,nz-1
                do is=1,1
                    xnorma(is)=xnorma(is)+cdabs(cpsi(ix,iy,iz))**2*dx**2*dz
                enddo
            enddo
        enddo
    enddo

    do ix=-n+1,n-1
        do iy=-ny+1,ny-1
            do iz=-nz+1,nz-1
                cpsi(ix,iy,iz)=cpsi(ix,iy,iz)/sqrt(xnorma(1))
            enddo
        enddo
    enddo
end


subroutine cnd(cpsii,cpsin,cpsi,n,ny,nz,vx,vy,vz,xm,dt,fi3d,fi3do,rdy)
    implicit double precision(a,b,d-h,o-z)
    implicit double complex(c)
    dimension cpsi (-N:N,-Ny:Ny,-Nz:Nz)
    dimension cpsin(-N:N,-Ny:Ny,-Nz:Nz)
    dimension cpsii(-N:N,-Ny:Ny,-Nz:Nz)
    dimension vdd(-N:N,-Ny:Ny,-Nz:Nz)
    dimension fi3d(-N:N,-Ny:Ny,-Nz:Nz)
    dimension fi3do(-N:N,-Ny:Ny,-Nz:Nz)
    dimension rdy(-2*N:2*N,-2*Ny:2*Ny,-2*Nz:2*Nz)
    dimension xnorma(2),xm(2),ene(2),v(2),xncz(2),xnormak(2),wd(2)
    dimension vx(-n:n),vy(-ny:ny),vz(-nz:nz)
    common/loc/xnormak,ene,eod11,eod22,eod12,iix,iiy,iiz
    common/xnorma/xnorma,wrl,wzl,dx,xncz,ggp11,dz,ggp12,wd,y0,omega
    common/add/add,cdd,gamma
    g=.1083e-21*0
    xnorma=0
    ci=(0.d0,1.d0)
    ene=0
    cpsin=cpsi
    cpsii=cpsi
    cdt=dt*(-ci)

    xnorma=0
    do ix=-n+1,n-1
        do iy=-ny+1,ny-1
            do iz=-nz+1,nz-1
                xnorma(1)=xnorma(1)+cdabs(cpsi(ix,iy,iz))**2*dx**2*dz
            enddo
        enddo
    enddo

    call lifi3(fi3d,cpsi,n,ny,nz,dx,dz,xnorma(1),xncz(1),rdy)
    fi3do=fi3d
    fi3d=fi3do

    do 10 icn=1,1
        ene=0
        w=sqrt(xncz(1))
        do ix=-n+1,n-1
            do iy=-ny+1,ny-1
                do iz=-nz+1,nz-1
                    x=ix*dx
                    y=iy*dx-y0
                    z=iz*dz
                    vvv=vx(ix)+vy(iy)+vz(iz)
                    c1=-0.5/xm(1)/dx**2*(                                       &
                        cpsi(ix-1,iy,iz)+cpsi(ix+1,iy,iz)+                      &
                        cpsi(ix,iy-1,iz)+cpsi(ix,iy+1,iz)-4*cpsi(ix,iy,iz))     &
                        -0.5/xm(1)/dz**2*(                                      &
                        cpsi(ix,iy,iz-1)+cpsi(ix,iy,iz+1)                       &
                        -2*cpsi(ix,iy,iz))                                      &
                        +cpsi(ix,iy,iz)*(vvv+cdd*fi3do(ix,iy,iz))
                    cpsii(ix,iy,iz)=cpsi(ix,iy,iz)+cdt/ci*c1
                enddo
            enddo
        enddo
        xnormak=xnorma
        eod11=0
        eod12=0
        eod22=0
        w=xncz(1)/xnorma(1)
        do ix=-n+1,n-1
            do iy=-ny+1,ny-1
                do iz=-nz+1,nz-1
                    cpsii(ix,iy,iz)=cpsii(ix,iy,iz)                             &
                        +cdt/ci*((ggp11-cdd/3)                                      &
                        *cdabs(cpsi(ix,iy,iz))**2*cpsi(ix,iy,iz)*(w)                &
                        +gamma*cdabs(cpsi(ix,iy,iz))**3*cpsi(ix,iy,iz)*w**1.5)
                enddo
            enddo
        enddo
        cpsin=cpsii

        eod22=ene(1)+eod11
10  continue
88  format(30g30.12)
end

function energiacnd(cpsi,n,ny,nz,vx,vy,vz,xm,dt,fi3d,rdy)
    implicit double precision(a,b,d-h,o-z)
    implicit double complex(c)
    dimension cpsi (-N:N,-Ny:Ny,-Nz:Nz)
    dimension rdy(-2*N:2*N,-2*Ny:2*Ny,-2*Nz:2*Nz)
    dimension vdd(-N:N,-Ny:Ny,-Nz:Nz)
    dimension fi3d(-N:N,-Ny:Ny,-Nz:Nz)
    dimension xnorma(2),xm(2),ene(2),v(2),xncz(2),xnormak(2),wd(2)
    dimension vx(-n:n),vy(-ny:ny),vz(-nz:nz)
    common/loc/xnormak,ene,eod11,eod22,eod12,iix,iiy,iiz
    common/xnorma/xnorma,wrl,wzl,dx,xncz,ggp11,dz,ggp12,wd,y0,omega
    common/add/add,cdd,gamma
    common/skladowe/e1,e2,e3,e4,e5
    g=.1083e-21*0
    xnorma=0
    ci=(0.d0,1.d0)
    ene=0
    cdt=dt*(-ci)

    xnorma=0
    do ix=-n+1,n-1
        do iy=-ny+1,ny-1
            do iz=-nz+1,nz-1
                xnorma(1)=xnorma(1)+cdabs(cpsi(ix,iy,iz))**2*dx**2*dz
            enddo
        enddo
    enddo

    e1=0
    e2=0
    e3=0
    e4=0
    e5=0
    w=sqrt(xncz(1))
    do ix=-n+1,n-1
        do iy=-ny+1,ny-1
            do iz=-nz+1,nz-1
                vvv=vx(ix)+vy(iy)+vz(iz)
                c1=-0.5/xm(1)/dx**2*(                                   &
                    cpsi(ix-1,iy,iz)+cpsi(ix+1,iy,iz)+                      &
                    cpsi(ix,iy-1,iz)+cpsi(ix,iy+1,iz)-4*cpsi(ix,iy,iz))     &
                    -0.5/xm(1)/dz**2*(                                      &
                    cpsi(ix,iy,iz-1)+cpsi(ix,iy,iz+1)                       &
                    -2*cpsi(ix,iy,iz))
                c2=+cpsi(ix,iy,iz)*(vvv+0*cdd*fi3d(ix,iy,iz)/2)
                c3=+cpsi(ix,iy,iz)*(0*vvv+cdd*fi3d(ix,iy,iz)/2)
                e1=e1+c1*dconjg(cpsi(ix,iy,iz))*dx**2*xncz(1)/xnorma(1)*dz
                e2=e2+c2*dconjg(cpsi(ix,iy,iz))*dx**2*xncz(1)/xnorma(1)*dz
                e3=e3+c3*dconjg(cpsi(ix,iy,iz))*dx**2*xncz(1)/xnorma(1)*dz
            enddo
        enddo
    enddo
    xnormak=xnorma
    eod11=0
    eod12=0
    eod22=0
    w=xncz(1)/xnorma(1)
    do ix=-n+1,n-1
        do iy=-ny+1,ny-1
            do iz=-nz+1,nz-1
                e3=e3-cdd/3*cdabs(cpsi(ix,iy,iz))**4*dx**2/2*xncz(1)**2/xnorma(1)**2*dz
                e4=e4+(ggp11-0*cdd/3)*cdabs(cpsi(ix,iy,iz))**4*dx**2/2*xncz(1)**2/xnorma(1)**2*dz
                e5=e5+gamma*cdabs(cpsi(ix,iy,iz))**5*w**2.5*.4*dx**2*dz
            enddo
        enddo
    enddo
    cpsin=cpsii

    eod22=ene(1)+eod11
    energiacnd=e1+e2+e3+e4+e5
88  format(30g30.12)
end

subroutine lifi3(fi3d,cpsi,n,ny,nz,dx,dz,xnorma,xncz,rdy)
    implicit double precision(a,b,d-h,o-z)
    implicit double complex (c)
    dimension fi3d(-N:N,-Ny:Ny,-Nz:Nz)
    dimension rdy(-2*N:2*N,-2*Ny:2*Ny,-2*Nz:2*Nz)
    dimension cpsi(-N:N,-Ny:Ny,-Nz:Nz)
    allocatable fi(:,:,:)
    allocatable psi(:,:,:)
    common/add/add,cdd,gamma
    pi=4*atan(1.0)
    allocate(fi(-N:N,-Ny:Ny,-Nz:Nz))
    allocate(psi(-N:N,-Ny:Ny,-Nz:Nz))

    xfct=dx**2*dz*xncz/xnorma
    do ix=-N,N
        do iy=-Ny,Ny
            do iz=-Nz,Nz
                psi(ix,iy,iz)=cdabs(cpsi(ix,iy,iz))**2
            enddo
        enddo
    enddo
    do ix=-N,N,2*N
        do iy=-Ny,Ny
            do iz=-Nz,Nz
                fi(ix,iy,iz)=0
            enddo
        enddo
    enddo

    do ix=-N,N
        do iy=-Ny,Ny,2*Ny
            do iz=-Nz,Nz
                fi(ix,iy,iz)=0
            enddo
        enddo
    enddo

    do ix=-N,N
        do iy=-Ny,Ny
            do iz=-Nz,Nz,2*nz
                fi(ix,iy,iz)=0
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
                            fi(ix,iy,iz)=fi(ix,iy,iz)+rdijk*psi(ixp,iyp,izp)*xfct
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
                            fi(ix,iy,iz)=fi(ix,iy,iz)+rdijk*psi(ixp,iyp,izp)*xfct
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
                            fi(ix,iy,iz)=fi(ix,iy,iz)+rdijk*psi(ixp,iyp,izp)*xfct
                        enddo
                    enddo
                enddo

            enddo
        enddo
    enddo
    ! rownanie poissona: nabla^2\phi=-rho
    ! dla rho=delta, phi=1/r/4/pi
    omega=1.92
    do iter=1,520
        do ix=-N+1,N-1
            do iy=-Ny+1,Ny-1
                do iz=-Nz+1,Nz-1
                    fi(ix,iy,iz)=fi(ix,iy,iz)*(1-omega)+        &
                        omega*                                      &
                        (                                           &
                        (fi(ix,iy+1,iz)+fi(ix,iy-1,iz)              &
                        +fi(ix+1,iy,iz)+fi(ix-1,iy,iz))*dz**2       &
                        +(fi(ix,iy,iz+1)+fi(ix,iy,iz-1))*dx**2      &
                        +psi(ix,iy,iz)*dx**2*dz**2*xncz/xnorma)/(4*dz**2+2*dx**2)
                enddo
            enddo
        enddo

    enddo

    fi3d=0
    do ix=-N+1,N-1
        do iy=-Ny+1,Ny-1
            do iz=-Nz+1,Nz-1
                fi3d(ix,iy,iz)=-(fi(ix,iy,iz+1)+fi(ix,iy,iz-1)-2*fi(ix,iy,iz))/dz**2
            enddo
        enddo
    enddo
end

subroutine cndt(cpsii,cpsin,cpsi,n,ny,nz,vx,vy,vz,xm,dt,fi3d,fi3do,rdy)
    implicit double precision(a,b,d-h,o-z)
    implicit double complex(c)
    dimension cpsi (-N:N,-Ny:Ny,-Nz:Nz)
    dimension rdy(-2*N:2*N,-2*Ny:2*Ny,-2*Nz:2*Nz)
    dimension cpsin(-N:N,-Ny:Ny,-Nz:Nz)
    dimension cpsii(-N:N,-Ny:Ny,-Nz:Nz)
    dimension vdd(-N:N,-Ny:Ny,-Nz:Nz)
    dimension fi3d(-N:N,-Ny:Ny,-Nz:Nz)
    dimension fi3do(-N:N,-Ny:Ny,-Nz:Nz)
    dimension xnorma(2),xm(2),ene(2),v(2),xncz(2),xnormak(2),wd(2)
    dimension vx(-n:n),vy(-ny:ny),vz(-nz:nz)
    common/loc/xnormak,ene,eod11,eod22,eod12,iix,iiy,iiz
    common/xnorma/xnorma,wrl,wzl,dx,xncz,ggp11,dz,ggp12,wd,y0,omega
    common/add/add,cdd,gamma
    g=.1083e-21*0
    xnorma=0
    ci=(0.d0,1.d0)
    ene=0
    cpsin=cpsi
    cpsii=cpsi
    cdt=dt

    xnorma=0
    do ix=-n+1,n-1
        do iy=-ny+1,ny-1
            do iz=-nz+1,nz-1
                xnorma(1)=xnorma(1)+cdabs(cpsi(ix,iy,iz))**2*dx**2*dz
            end do
        end do
    end do

    call lifi3(fi3d,cpsi,n,ny,nz,dx,dz,xnorma(1),xncz(1),rdy)
    fi3do=fi3d
    fi3d=fi3do

    do icn = 1, 2
        ene=0
        w=sqrt(xncz(1))
        do ix=-n+1,n-1
            do iy=-ny+1,ny-1
                do iz=-nz+1,nz-1
                    do is=1,1
                        vvv=vx(ix)+vy(iy)+vz(iz)
                        c1=-0.5/xm(1)/dx**2*(                                   &
                            cpsi(ix-1,iy,iz)+cpsi(ix+1,iy,iz)+                      &
                            cpsi(ix,iy-1,iz)+cpsi(ix,iy+1,iz)-4*cpsi(ix,iy,iz))     &
                            -0.5/xm(1)/dz**2*(                                      &
                            cpsi(ix,iy,iz-1)+cpsi(ix,iy,iz+1)                       &
                            -2*cpsi(ix,iy,iz))                                      &
                            +cpsi(ix,iy,iz)*(vvv+cdd*fi3do(ix,iy,iz))

                        c2=-0.5/xm(1)/dx**2*(                                   &
                            cpsin(ix-1,iy,iz)+cpsin(ix+1,iy,iz)+                    &
                            cpsin(ix,iy-1,iz)+cpsin(ix,iy+1,iz)-4*cpsin(ix,iy,iz))  &
                            -0.5/xm(1)/dz**2*(                                      &
                            cpsin(ix,iy,iz-1)+cpsin(ix,iy,iz+1)                     &
                            -2*cpsin(ix,iy,iz))                                     &
                            +cpsin(ix,iy,iz)*(vvv+cdd*fi3d(ix,iy,iz))               

                            cpsii(ix,iy,iz)=cpsi(ix,iy,iz)+cdt/ci*(c1+c2)/2
                    end do
                end do
            end do
        end do
        xnormak=xnorma
        eod11=0
        eod12=0
        eod22=0
        w=xncz(1)/xnorma(1)
        do ix=-n+1,n-1
            do iy=-ny+1,ny-1
                do iz=-nz+1,nz-1
                    cpsii(ix,iy,iz)=cpsii(ix,iy,iz)                             &
                        +cdt/ci*((ggp11-cdd/3)                                      &
                        *cdabs(cpsi(ix,iy,iz))**2*cpsi(ix,iy,iz)*(w)                &
                        +gamma*cdabs(cpsi(ix,iy,iz))**3*cpsi(ix,iy,iz)*w**1.5)/2    &
                        +cdt/ci*((ggp11-cdd/3)                                      &
                        *cdabs(cpsin(ix,iy,iz))**2*cpsin(ix,iy,iz)*(w)              &
                        +gamma*cdabs(cpsin(ix,iy,iz))**3*cpsin(ix,iy,iz)*w**1.5)/2  
                end do
            end do
        end do
        cpsin=cpsii

        eod22=ene(1)+eod11
        if(mod(iter-1,2).eq.0) then
            call lifi3(fi3d,cpsin,n,ny,nz,dx,dz,xnorma(1),xncz(1),rdy)
            fi3d=fi3d
        endif
    end do
88  format(30g30.12)
end
