C
C Compensator plugin for NEK5000
C
C Nicolo Fabbiane (nicolo@mech.kth.se), January 2015
C Pierluigi Morra (pmorra@mech.kth.se), August  2017
C      
      module cmp

         implicit none
         
         include 'CTRL'

         character*80 infile
         parameter (infile='cmp.i')
         
         logical             :: restart,reference,fxlms,savek,meanflag

         real                :: t,dt,tsave

         integer             :: nyy,nuu
         real                :: yy(nyiov),uu(nuiov)         

         real                :: dtint,dtsub,dterr
         real                :: yyint(nyiov)
                  
         integer             :: nyybuf
         real   ,allocatable :: yybuf(:,:)

         integer             :: ny,nu,nz
         integer,allocatable :: iu(:),iy(:),iz(:)
         real   ,allocatable :: y(:),u(:),z(:),f(:)

         integer             :: nybuf,nfbuf
         real   ,allocatable :: ybuf(:,:),fbuf(:,:)
      
         character*80        :: kfile,pfile
         integer             :: ntk(2),nyk(2),ntp(2),nzp(2)
         real                :: dtk,dtp
         real   ,allocatable :: K(:,:),P(:,:)

         real                :: mumax,Lp2Deq,tonfxlms,zref         
         integer             :: ksave,iocmp
         character*11        :: ioformat


      end module cmp

c-----------------------------------------------------------------------

      subroutine cmp_init(nyysim,nuusim,fid)

      use cmp
      implicit none

      integer,intent(in) :: nyysim,nuusim,fid

      real    tau2,sigma2
      logical cmpstop
      integer i,j
      character*80 name

c
c     CHECK FOR INPUT FILE
c
      inquire(file=infile,exist=cmpstop); cmpstop=.not.cmpstop
      if (cmpstop) then

         write(*,*) 'Warning: ',trim(infile),' does not exist:',
     &              ' simulation will stop.'

         stop

      end if
c
c     Read infile
c
      nyy = nyysim
      nuu = nuusim

      open(unit=fid,file=infile,form='formatted')
      
      yy(:) = 0.0d0; yyint(:) = 0.0d0
      uu(:) = 0.0d0

      read(fid,*) ny
      read(fid,*) nu

      allocate(iy(ny),y(ny))
      allocate(iu(nu),u(nu))
      allocate(f(ny))

      read(fid,*) iy
      read(fid,*) iu

      read(fid,*) dt
      dterr = 1d-6 * dt

      read(fid,*) kfile

      read(fid,*) fxlms
      
      meanflag = .false.
      reference = .false.
      if (fxlms) then

         read(fid,*) tonfxlms
         read(fid,*) mumax

         read(fid,*) nz
         allocate(iz(nz),z(nz))
         read(fid,*) iz

         read(fid,*) pfile
        
         read(fid,*) meanflag

         read(fid,*) reference
         if (reference) then
            read(fid,*) zref
         end if

      else
         
         nz = 0;
         allocate(iz(nz),z(nz))
         
      end if

      read(fid,*) savek
      if (savek) then

         ksave = 1
         read(fid,*) tsave

      end if

      read(fid,*) restart

      close(fid)
c
c     Read kernel files
c
      if (restart) then
         open(unit=fid,file='K.restart',form='formatted')
      else
         open(unit=fid,file=kfile,form='formatted')
      end if

      read(fid,*) dtk; if (dt.ne.dtk) stop
      read(fid,*) ntk
      read(fid,*) nyk
      allocate(K(nyk(2)-nyk(1)+1,ntk(2)-ntk(1)+1))
      do i=1,ntk(2)-ntk(1)+1
         read(fid,*) K(:,i)
      end do

      close(fid)

      if (fxlms) then

         open(unit=fid,file=pfile,form='formatted')
      
         read(fid,*) dtp; if (dt.ne.dtp) stop
         read(fid,*) ntp
         read(fid,*) nzp
         allocate(P(nzp(2)-nzp(1)+1,ntp(2)-ntp(1)+1))
         do i=1,ntp(2)-ntp(1)+1
            read(fid,*) P(:,i)
         end do

         close(fid)

      else
         ntp = 0; nzp = nyk
      end if
c      
c     mumax scaling
c
c     For reference see: Theoretical convergence analysis of
c     FxLMS algorithm (I. Tabatabaei Ardekani, W.H.Abdulla,
c     2010).
c
      if (fxlms) then
      ! tau2, sigma2
      tau2   = 0.0
      sigma2 = 0.0
      do j=1,nzp(2)-nzp(1)+1
        do i=1,ntp(2)-ntp(1)+1
          tau2 = tau2 + i*P(j,i)**2
          sigma2 = sigma2 + P(j,i)**2
        end do
      end do

      !L+2*delta_eq definition
      Lp2Deq = ntk(2)+2*tau2/sigma2
c
c     mumax rescaling w/o 1/Pxf (see mentioned paper)
c     Pxf is computed in lms_step
c     NOTE: the 0.5 factor is due to the problem definition.
c           See Nicolò Fabbiane PhD Thesis (Paper 1 or 5).
c

      mumax = 0.5*mumax/Lp2Deq
      end if

c
c     Initialise buffers
c
      nybuf = max(ntp(2),ntk(2))
      allocate(ybuf(ny,nybuf)); ybuf = 0.0d0
     
      nyybuf = nybuf
      allocate(yybuf(nyiov,nyybuf)); yybuf = 0.0d0

      if (restart) then
      name='ybuf.restart'; call read_buffer(name,t,ybuf,ny,nybuf,fid)
      name='yybuf.restart';call read_buffer(name,t,yybuf,nyy,nyybuf,fid)
      end if

      if (fxlms) then

         nfbuf = ntk(2)
         allocate(fbuf(ny,nfbuf)); fbuf = 0.0d0

         if (restart) then
         name='fbuf.restart'; call read_buffer(name,t,fbuf,ny,nfbuf,fid)  
         end if

      end if
c
c     Initialize output
c
      iocmp=fid
      open(unit=iocmp,file='cmp.o',status='unknown')
      write(ioformat,'(a,I3.0,a)') '(',(1+nyy+nuu),'E24.16)'
c
c     Initlialize time
c
      if (.not.restart) t=0.0d0

      end subroutine cmp_init

c-----------------------------------------------------------------------

      subroutine cmp_step(uusim,yysim,dtsim,fid)

      use cmp, only: nuiov,nyiov,
     &               t,dt,
     &               dtint,dtsub,dterr,yyint,
     &               meanflag,reference,zref,
     &               yy,nyy,yybuf,nyybuf,
     &               uu,nuu,
     &               y,ny,iy,ybuf,nybuf,
     &               z,nz,iz,
     &               u,nu,iu,
     &               f,fbuf,nfbuf,
     &               K,nyk,ntk,
     &               fxlms,mumax,tonfxlms,
     &               P,nzp,ntp,
     &               savek,ksave,tsave,
     &               iocmp,ioformat

      implicit none

      real,   intent(out):: uusim(nuiov)
      real,   intent(in) :: yysim(nyiov),dtsim
      integer,intent(in) :: fid

      integer i,j,nmean
      character*80 name

c
c     Downsampling
c
      dtsub = min(dtsim,dt-dtint)

      dtint = dtint + dtsub
      yyint = yyint + dtsub*yysim

      if (dtint.ge.(dt-dterr)) then

         do i=1,nyy
            yy(i) = yyint(i)/dt
         end do

         dtint = dtsim - dtsub
         yyint = dtint * yysim
c
c     Subtract the mean to measurement
c
         do i=nyybuf,2,-1
            yybuf(:,i) = yybuf(:,i-1)
         end do
         yybuf(:,1) = yy(:)

         if (meanflag) then

           nmean = min(nyybuf,int(t/dt)+1)
           do j=1,nmean
              do i=1,ny
                 yy(iy(i)) = yy(iy(i)) - yybuf(iy(i),j)/nmean
              end do
           end do

           if (.not.reference) then
              do j=1,nmean
                 do i=1,nz
                    yy(iz(i)) = yy(iz(i)) - yybuf(iz(i),j)/nmean
                 end do
              end do
           end if
         end if

         if (reference) then
            do i=1,nz
               yy(iz(i)) = yy(iz(i)) - zref
            end do
         end if
c
c     Measurement decomposition
c
         do i=1,ny
            y(i) = yy(iy(i))
         end do

         if (fxlms) then
            do i=1,nz
               z(i) = yy(iz(i))
            end do
         end if

         do i=nybuf,2,-1
            ybuf(:,i) = ybuf(:,i-1)
         end do
         ybuf(:,1) = y(:)
c
c     Computing the inputs
c        
         write(*,*)
         write(*,*) 'COMPUTING THE INPUTS'
         write(*,*) 
         write(*,*) ' FOLLOWING... KAPPA: ',K(1,1),K(1,2)
         call calc_fir(u,K,ybuf,nyk,ntk,nu,ny,nybuf)
         uu = 0.0d0
         do i=1,nu
            uu(iu(i)) = u(i)
         end do
         write(*,*) 'THE uu,u : ',uu(2),u(1)
c
c     Adaptation step
c
         if (fxlms) then

            call calc_fir(f,P,ybuf,nzp,ntp,ny,ny,nybuf)

            do i=nfbuf,2,-1
               fbuf(:,i) = fbuf(:,i-1)
            end do
            fbuf(:,1) = f(:)

            if (t.ge.tonfxlms) then
               call lms_step(mumax,z,K,fbuf,nyk,ntk,nz,ny,nfbuf)
            end if

         end if
c
c     Save Kernel & buffers
c
         if (savek) then
         if (t.ge.(ksave*tsave - dt/2)) then

            ksave = ksave+1
            
            write(name,'(a,E9.3)') 'K.restart'
            call save_fir(name,K,dt,nyk,ntk,fid)
            
            write(name,'(a,E9.3)') 'ybuf.restart'
            call save_buffer(name,t,ybuf,ny,nybuf,fid)
            
            write(name,'(a,E9.3)') 'yybuf.restart'
            call save_buffer(name,t,yybuf,nyy,nyybuf,fid)
            
            if (fxlms) then
               write(name,'(a,E9.3)') 'fbuf.restart'
               call save_buffer(name,t,fbuf,ny,nfbuf,fid)
            end if
            
         end if
         end if
c     
c     Finalise compensator step
c     
         write(iocmp,ioformat) t,
     &        (yy(j),j=1,nyy),(uu(j),j=1,nuu)
         flush(iocmp)

         t=t+dt

      end if
c
c     Send control signal
c
      do i=1,nuu
         uusim(i) = uu(i)
      end do
      write(*,*) ' STE CAZZO di uusmi: ',uusim(:) 
      end subroutine cmp_step

c-----------------------------------------------------------------------

      subroutine cmp_stop(fid)

      use cmp, only: t,dt,
     &               yybuf,nyy,nyybuf,
     &               ybuf,ny,nybuf,
     &               K,nyk,ntk,
     &               fxlms,fbuf,nfbuf,
     &               savek,ksave,tsave,
     &               iocmp
      implicit none

      integer fid
      character*80 name

c
c     Finilize output
c
      close(iocmp)

c
c     Save Kernel & buffers
c
      if (savek) then
         if (t.ge.(ksave*tsave - dt/2)) then
            
            ksave = ksave+1
            
            write(name,'(a,E9.3)') 'K.restart'
            call save_fir(name,K,dt,nyk,ntk,fid)
            
            write(name,'(a,E9.3)') 'ybuf.restart'
            call save_buffer(name,t,ybuf,ny,nybuf,fid)
            
            write(name,'(a,E9.3)') 'yybuf.restart'
            call save_buffer(name,t,yybuf,nyy,nyybuf,fid)
            
            if (fxlms) then
               write(name,'(a,E9.3)') 'fbuf.restart'
               call save_buffer(name,t,fbuf,ny,nfbuf,fid)
            end if
            
         end if
      end if
      
      end subroutine cmp_stop

c-----------------------------------------------------------------------

      subroutine calc_fir(out,H,inbuf,ninh,nth,nout,nin,nbuf)
      
      implicit none
      
      integer,intent(in) :: ninh(2),nth(2),nout,nin,nbuf
      
      real,intent(out)   :: out(nout)
      real,intent(in)    :: inbuf(nin,nbuf),
     &                      H(ninh(2)-ninh(1)+1,nth(2)-nth(1)+1)

      integer j,l,m,ml,l0

      out = 0.0d0

      do m = 1,nout
         do j = nth(1),nth(2)
            do l = ninh(1),ninh(2)
               ml = mod(m+l-1 + nin,nin) + 1
c                             + nin: to fix a mod(a,b) problem when a crosses zero
               out(m) = out(m) + H(l-ninh(1)+1,j-nth(1)+1) * inbuf(ml,j)
            end do
         end do
      end do

      end subroutine calc_fir

c-----------------------------------------------------------------------

      subroutine lms_step(mumax,e,H,ubuf,nyh,nth,ne,nu,nbuf)

      implicit none
      
      integer, intent(in) :: nyh(2),nth(2),ne,nu,nbuf
      
      real,intent(in)     :: mumax,e(ne),ubuf(nu,nbuf)
      real,intent(inout)  :: H(nyh(2)-nyh(1)+1,nth(2)-nth(1)+1)

      integer j,l,m,ml
      real lambda(nyh(2)-nyh(1)+1,nth(2)-nth(1)+1),mu,Pxf

c
c     Update direction
c
      lambda = 0.0d0
      do l = 1,ne
         do j = nth(1),nth(2)
            do m = nyh(1),nyh(2) 
               ml = mod(m+l-1 + nu,nu) + 1
c                             + nu: to fix a mod(x,y) problem when x crosses zero
               lambda(m-nyh(1)+1,j-nth(1)+1) = 
     &              lambda(m-nyh(1)+1,j-nth(1)+1) - 2*e(l)*ubuf(ml,j)
            end do
         end do
      end do
c
c     Step length
c
      mu = 0.0d0
c
c     Calculate Pxf
c     (For reference: see mentioned paper on mumax rescaling)
c
      Pxf = 0.0d0
      do j = 1,nbuf
         do m = 1,nu
            Pxf = Pxf + ubuf(m,j)**2/nbuf/nu
         end do
      end do
      
      mu = mumax/Pxf

c
c     Kernel update
c
      do j = 1,nth(2)-nth(1)+1
         do l = 1,nyh(2)-nyh(1)+1
            H(l,j) = H(l,j) + mu*lambda(l,j)
         end do
      end do
      
      end subroutine lms_step

c-----------------------------------------------------------------------

      subroutine read_buffer(name,t,buf,nu,nbuf,fid)
      
      implicit none

      character*80,intent(in) :: name
      integer,intent(in)      :: nu,nbuf,fid
      real,intent(out)        :: t,buf(nu,nbuf)

      integer j

      open(unit=fid,file=trim(name),form='formatted')

      read(fid,*) t
      do j=1,nbuf
         read(fid,*) buf(:,j)
      end do

      close(fid)

      end subroutine read_buffer

c-----------------------------------------------------------------------

      subroutine save_buffer(name,t,buf,nu,nbuf,fid)
      
      implicit none

      character*80,intent(in) :: name
      integer,intent(in)      :: nu,nbuf,fid
      real,intent(in)         :: t,buf(nu,nbuf)

      integer i,j

      open(unit=fid,file=trim(name))

      write(fid,'(E24.15)') t
      do j=1,nbuf
         write(fid,'(10000E24.15)') (buf(i,j),i=1,nu)
      end do

      close(fid)

      end subroutine save_buffer

c-----------------------------------------------------------------------

      subroutine save_fir(name,H,dt,nyh,nth,fid)

      implicit none

      character*80,intent(in) :: name
      integer,intent(in) :: nyh(2),nth(2),fid

      real,intent(in)    :: H(nyh(2)-nyh(1)+1,nth(2)-nth(1)+1),dt

      integer i,j
            
      open(unit=fid,file=trim(name))
      
      write(fid,'(E24.15)') dt
      write(fid,'(2I10)') nth(1),nth(2)
      write(fid,'(2I10)') nyh(1),nyh(2)
      
      do i=1,nth(2)-nth(1)+1
         write(fid,'(20E24.15)') (H(j,i),j=1,nyh(2)-nyh(1)+1)
      end do
      
      close(fid)
      
      end subroutine save_fir
