C
C	Control plugin for NEK5000
C
C	Nicolo Fabbiane (nicolo@mech.kth.se), October 2015
C      
      module ctrl

         implicit none

         include 'SIZE'
         include 'CTRL'

         character*80 infile
         parameter (infile='ctrl.i')
         
         logical ctrlflag,cmpflag
         integer nshio,nyio,nuio,ndio,ngio
         
         integer ishyio(nyiov),ishuio(nuiov)
         real shiox(lx1,ly1,lz1,lelv,nshiov)
         real shioy(lx1,ly1,lz1,lelv,nshiov)
         real shioz(lx1,ly1,lz1,lelv,nshiov)

         integer diotype(ndiov),idio(ndiov)
         integer giotype(ngiov),igio(ngiov)
         real ampdio(ndiov),ampgio(ngiov)
      
         real yio(nyiov),uio(nuiov),dio(ndiov),gio(ngiov)
         
         real buf(lx1,ly1,lz1,lelv)

         integer ioctrl(2),fid
         character*11 ioformat

         integer iobj_sen(nyiov) ! JF

      end module ctrl

c-----------------------------------------------------------------------
      
      subroutine ctrl_init
      
      use ctrl
      implicit none

      character*80 fname
      integer i
      integer ishio ! JF
      
      nshio=0; nyio=0; nuio=0; ndio=0; ngio=0;
      ishyio=0; ishuio=0;
      shiox=0.0d0; shioy=0.0d0; shioz=0.0d0
      yio=0.0d0; uio=0.0d0 
      diotype=0; idio=0; ampdio=0.0d0;
      giotype=0; igio=0; ampgio=0.0d0;  

      ! check for control.i to exist
      inquire(file=infile,exist=ctrlflag)
         
      if (ctrlflag) then

         ! open file
         fid = 41
         open(unit=fid,file=infile,form='formatted')
         
         if (NID.eq.0) then

            write(*,*) 'Reading control.i:'

         end if

         ! number of i/o shapes
         read(fid,*) nshio

         if (NID.eq.0) then

            write(*,*) '  nshio       : ',nshio

         end if

         ! read i/o shapes
         do i=1,nshio

            read(fid,*) fname
            call readnek(fname,
     &               shiox(1,1,1,1,i),shioy(1,1,1,1,i),shioz(1,1,1,1,i))

         end do

         ! number of outputs
         read(fid,*) nyio

         if (NID.eq.0) then

            write(*,*) '  nyio        : ',nyio

            if (nyio.gt.nyiov) then

               write(*,*) 'ERROR: nyio > nyiov (',nyio,' > ',nyiov,')'
               stop

            end if

         end if

         ! read outputs (shape number)
         do i=1,nyio

            read(fid,*) ishyio(i)

         end do

         ! number of inputs
         read(fid,*) nuio

         if (NID.eq.0) then

            write(*,*) '  nuio        : ',nuio

            if (nuio.gt.nuiov) then

               write(*,*) 'ERROR: nuio > nuiov (',nuio,' > ',nuiov,')'
               stop

            end if

         end if

         ! read inputs (shape number)
         do i=1,nuio

            read(fid,*) ishuio(i)

         end do

         ! compensator flag
         read(fid,*) cmpflag

         ! number of input disturbances
         read(fid,*) ndio

         if (NID.eq.0) then

            write(*,*) '  ndio        : ',ndio

            if (ndio.gt.ndiov) then

               write(*,*) 'ERROR: ndio > ndiov (',ndio,' > ',ndiov,')'
               stop

            end if

         end if

         ! read input disturbances
         do i=1,ndio

            read(fid,*) diotype(i)
            read(fid,*) idio(i)
            read(fid,*) ampdio(i)
            
            if (NID.eq.0) then

               select case (diotype(i))

               case (1)

                  write(*,*) '    diotype   : 1. Impulse'

               case (2)
                  write(*,*) '    diotype   : 2. Random noise'

               case (3)
                  write(*,*) '    diotype   : 3. Constant'

               end select

               write(*,*) '    idio      : ',idio(i)
               write(*,*) '    ampdio    : ',ampdio(i)

            end if

         end do

         ! number of output disturbances
         read(fid,*) ngio

         if (NID.eq.0) then

            write(*,*) '  ngio        : ',ngio

            if (ngio.gt.ngiov) then

               write(*,*) 'ERROR: ngio > ngiov (',ngio,' > ',ngiov,')'
               stop

            end if

         end if

         ! read input disturbances
         do i=1,ngio

            read(fid,*) giotype(i)
            read(fid,*) igio(i)
            read(fid,*) ampgio(i)
            
            if (NID.eq.0) then

               select case (giotype(i))

               case (1)
                  write(*,*) '    giotype   : 1. Impulse'

               case (2)
                  write(*,*) '    giotype   : 2. Random noise'

               case (3)
                  write(*,*) '    giotype   : 3. Constant'

               end select

               write(*,*) '    igio      : ',igio(i)
               write(*,*) '    ampgio    : ',ampgio(i)

            end if

         end do

         ! close control.i
         close(fid)

         ! initialize output
         if (NID.eq.0) then

            ioctrl(1) = 42
            ioctrl(2) = 43

            open(unit=ioctrl(1),file='ctrl.o',status='unknown')
C            open(unit=ioctrl(2),file='energy.o',status='unknown')
            write(ioformat,'(a,I3.0,a)') '(',(1+nyio+nuio
     &                                         +ndio+ngio),'E24.16)'

         end if

         !initialize compensator
         if (NID.eq.0) then

            if (cmpflag) call cmp_init(nyio,nuio,44) ! TODO: replace 44 with free-fid function

         end if
c---------------- Added by Jose
         ! Initalize objects for sensor
         do i=1,nyio
           ishio = ishyio(i)
           call set_sensor(shiox(1,1,1,1,ishio),iobj_sen(i))
         enddo
         
c----------------------
      else
         
         if (NID.eq.0) write(*,*) 'ctrl.i not found.'
         
      end if

      end subroutine ctrl_init

c-----------------------------------------------------------------------

      subroutine ctrl_get_dist
      
      use ctrl, only: ctrlflag,
     &                ndio,diotype,idio,ampdio,dio,
     &                ngio,giotype,igio,ampgio,gio
      implicit none

      include 'SIZE'
      include 'TSTEP'

      integer i,j

      dio = 0.0d0;
      gio = 0.0d0;
      
      if (ctrlflag) then
         
         if (NID.eq.0) then
         
            ! compute input disturbances
            do i=1,ndio

               select case (diotype(i))

               case (1) ! Impulse
                  if (istep.lt.1) dio(i) = ampdio(i)/DT;

               case (2) ! Random
                  call random_number(dio(i))
                  dio(i) = 2.*(dio(i)-.5d0)*ampdio(i)

               case (3) ! Constant
                  dio(i) = ampdio(i)

               end select

            end do
         
            ! compute output disturbances
            do i=1,ngio

               select case (giotype(i))

               case (1) ! Impulse
                  if (TIME.lt.DT) gio(i) = ampgio(i)/DT;

               case (2) ! Random
                  call random_number(gio(i))
                  gio(i) = 2.*(gio(i)-.5d0)*ampgio(i)

               case (3) ! Constant
                  gio(i) = ampgio(i)

               end select
            end do
            
         end if

         ! broadcast disturbances
         call bcast(dio(1),8*ndio)
         call bcast(gio(1),8*ngio)
      
      end if
      
      end subroutine ctrl_get_dist

c-----------------------------------------------------------------------
      subroutine ctrl_get_output2()
c    This function makes use of objects to define sensors,so make sure
c    of being consistent if you create new objects or they might be
c    considered as sensors.
      use ctrl, only: ctrlflag,
     &                nyio,ishyio,yio,
     &                ngio,igio,gio,
     &                shiox,shioy,shioz,buf,
     &                iobj_sen
      
      INCLUDE 'SIZE'  
      INCLUDE 'TOTAL' 

      real*8  duT(lx1,ly1,lz1,lelv)
      real yiol,yiog
      integer ntot, ish
      real Stot, dsl,dss
      real to_file(nyio+1)
c
      ntot = lx1*ly1*lz1*lelv
c	
c     
      yio(:)=0.0d0

      if (ctrlflag) then      ! compute output
         
         if (npert.eq.0) then ! non linear simulation

            do i=1,nyio
             ish = ishyio(i)      ! get shape        
             yiol = 0
             dsl = 0
             iobj = iobj_sen(i)
c             do iobj = 1,nobj  !nobj should match with nyio
                memtot = nmember(iobj)
             do mem  = 1,memtot
                ieg   = object(iobj,mem,1)
                ifc   = object(iobj,mem,2)
                if (gllnid(ieg).eq.nid) then ! this processor has a contribution
                   ie = gllel(ieg)
                   call surface_int_J(yiog,dss,
     $                  shiox(1,1,1,ie,ish),ie,ifc)
                   yiol = yiol + yiog 
                   dsl = dsl + dss  
                endif
             enddo
c             enddo
             Stot = glsum(dsl,1) !Total area
             yio(i) = glsum(yiol,1) !Total integral 
             yio(i) = yio(i)/Stot
c             if(nid.eq.0) write(*,*) 'yio jose ',yio(i)
             enddo !nyio do
         else
             write(*,*) "IMPLEMENT PERTURBATIOOON!!"
c       TODO  I think the only change should be when computing the
c         the tangent velocity in surface_int_J    
           endif
        
        
        ! add disturbances
         do i=1,ngio
            yio(igio(i)) = yio(igio(i)) + gio(i);
         end do
         
         if(nid.eq.0) then
            to_file(1) = time
            do i=2,nyio+1
              to_file(i) = yio(i-1)
            enddo
c            write(*,*) "TEST", nyio
c            do i=1,nyio+1
c              write(*,*) to_file(i)
c            enddo
            if(istep.eq.0) then 
              open(888,file='output.dat',status='new',
     $              action='write')
            else 
              open(888,file="output.dat",status="old",
     $          position="append",action="write")
            endif
c            do i=1,nyio
              write(888,111) (to_file(i),i=1,nyio+1)
  111  format(600(F16.10))      
c            enddo
            close(888)
         endif
        
        endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ctrl_get_output3()
c    This function makes use of objects to define sensors,so make sure
c    of being consistent if you create new objects or they might be
c    considered as sensors.
      use ctrl, only: ctrlflag,
     &                nyio,ishyio,yio,
     &                ngio,igio,gio,
     &                shiox,shioy,shioz,buf,
     &                iobj_sen
      
      INCLUDE 'SIZE'  
      INCLUDE 'TOTAL' 

      real yiol,yiog
      integer nlxyz, ish
      real*8  uxp(lx1,ly1,lz1),uyp(lx1,ly1,lz1),
     $     uzp(lx1,ly1,lz1)
      real*8  derT(lx1,ly1,lz1)
      real to_file(nyio+1)
      
      real*8  ux_bc(lx1,ly1,lz1,lelv),uy_bc(lx1,ly1,lz1,lelv),
     $     uz_bc(lx1,ly1,lz1,lelv)
      COMMON / usrboundc / ux_bc,uy_bc,uz_bc
c
      nlxyz = lx1*ly1*lz1
c	
c     
      yio(:)=0.0d0

      if (ctrlflag) then      ! compute output
         
         if (npert.eq.0) then ! non linear simulation

            do i=1,nyio
             ish = ishyio(i)      ! get shape        
             yiol = 0
             iobj = iobj_sen(i)
c             do iobj = 1,nobj  !nobj should match with nyio
                memtot = nmember(iobj)
             do mem  = 1,memtot
                ieg   = object(iobj,mem,1)
                ifc   = object(iobj,mem,2)
                if (gllnid(ieg).eq.nid) then ! this processor has a contribution
                ie = gllel(ieg)
c                write(*,*) mem, nid, ie, ieg
                call sub3(uxp,vx(1,1,1,ie),
     $             ux_bc(1,1,1,ie),nlxyz)
                call sub3(uyp,vy(1,1,1,ie),
     $             uy_bc(1,1,1,ie),nlxyz)
                call sub3(uzp,vz(1,1,1,ie),
     $              uz_bc(1,1,1,ie),nlxyz)
                call rzero(derT,nlxyz) ! initialize buffer

                call Xaddcol3(derT,shiox(1,1,1,ie,ish),
     $                  uxp,nlxyz) ! cumulative point-wise product
                call Xaddcol3(derT,shioy(1,1,1,ie,ish),
     $                  uyp,nlxyz) ! cumulative point-wise product
                call Xaddcol3(derT,shioz(1,1,1,ie,ish),
     $                  uzp,nlxyz) ! cumulative point-wise product

                yiog = vlsum(derT,nlxyz) ! sum element contribution
                   
                yiol = yiol + yiog 
c                write(*,*) yiog,yiol
                endif
             enddo
c             enddo
             yio(i) = glsum(yiol,1) !Total integral 
c             if(nid.eq.0) write(*,*) 'yio jose 3 ',yio(i)
             enddo !nyio do
         else
             write(*,*) "IMPLEMENT PERTURBATIOOON!!"
c       TODO     
           endif
        
        
        ! add disturbances
         do i=1,ngio
            yio(igio(i)) = yio(igio(i)) + gio(i);
         end do
         
         if(nid.eq.0) then
            to_file(1) = time
            do i=2,nyio+1
              to_file(i) = yio(i-1)
            enddo
            if(istep.eq.0) then 
              open(888,file='output2.dat',status='new',
     $              action='write')
            else 
              open(888,file="output2.dat",status="old",
     $          position="append",action="write")
            endif
c            do i=1,nyio
              write(888,111) (to_file(i),i=1,nyio+1)
  111  format(600(F16.10))      
c            enddo
            close(888)
         endif
         
        
        endif
c
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine ctrl_get_output

      use ctrl, only: ctrlflag,
     &                nyio,ishyio,yio,
     &                ngio,igio,gio,
     &                shiox,shioy,shioz,buf
      implicit none
      
      include 'SIZE'
      include 'SOLN'

      integer i,ish,ntot
      
      real glsum
      
      real*8  uxp(lx1,ly1,lz1,lelv),uyp(lx1,ly1,lz1,lelv),
     $     uzp(lx1,ly1,lz1,lelv)
      
      real*8  ux_bc(lx1,ly1,lz1,lelv),uy_bc(lx1,ly1,lz1,lelv),
     $     uz_bc(lx1,ly1,lz1,lelv)
      COMMON / usrboundc / ux_bc,uy_bc,uz_bc
      

      ntot = lx1*ly1*lz1*lelv
      yio(:)=0.0d0
c     Perturbation velocity
      call sub3(uxp,vx,ux_bc,ntot)
      call sub3(uyp,vy,uy_bc,ntot)
      call sub3(uzp,vz,uz_bc,ntot)
      
c      call copy(uxp,vx,ntot)
c      call copy(uyp,vy,ntot)
c      call copy(uzp,vz,ntot)


      if (ctrlflag) then      ! compute output
         
         if (npert.eq.0) then ! non linear simulation

            do i=1,nyio

               ish = ishyio(i)      ! get shape        

               call rzero(buf,ntot) ! initialize buffer

               call Xaddcol3(buf,shiox(1,1,1,1,ish),uxp,ntot) ! cumulative point-wise product
               call Xaddcol3(buf,shioy(1,1,1,1,ish),uyp,ntot) ! cumulative point-wise product
               call Xaddcol3(buf,shioz(1,1,1,1,ish),uzp,ntot) ! cumulative point-wise product

               yio(i) = glsum(buf(1,1,1,1),ntot) ! sum over the domain
               if(nid.eq.0) write(*,*) 'yio',yio(i)

            end do
            
         else ! perturbation

            do i=1,nyio

               ish = ishyio(i)      ! get shape        

               call rzero(buf,ntot) ! initialize buffer

               call Xaddcol3(buf,shiox(1,1,1,1,ish),vxp(1,1),ntot) ! cumulative point-wise product
               call Xaddcol3(buf,shioy(1,1,1,1,ish),vyp(1,1),ntot) ! cumulative point-wise product
               call Xaddcol3(buf,shioz(1,1,1,1,ish),vzp(1,1),ntot) ! cumulative point-wise product

               yio(i) = glsum(buf(1,1,1,1),ntot) ! sum over the domain

            end do         
            
         end if

         ! add disturbances
         do i=1,ngio
            yio(igio(i)) = yio(igio(i)) + gio(i);
         end do

         
      end if

      end subroutine ctrl_get_output

c-----------------------------------------------------------------------

      subroutine ctrl_get_input

      use ctrl, only: ctrlflag,cmpflag,
     &                nyio,yio,
     &                nuio,uio,
     &                ndio,idio,dio,
     &                fid
      implicit none
      
      include 'SIZE'
      include 'TSTEP'

      integer i
      
      uio=0.0d0
      
      if (ctrlflag) then
         
         if (NID.eq.0) then
            
            if (cmpflag) then
               
               call cmp_step(uio,yio,DT,fid)
            
            end if

            ! add disturbances
            do i=1,ndio
               uio(idio(i)) = uio(idio(i)) + dio(i)
            end do

         end if
         
         ! broadcast inputs
         call bcast(uio(1),8*nuio)

      end if
      
      end subroutine ctrl_get_input

c-----------------------------------------------------------------------

      subroutine ctrl_dump

      use ctrl, only: ctrlflag,
     &                nyio,yio,
     &                nuio,uio,
     &                ndio,dio,
     &                ngio,gio,
     &                ioctrl,ioformat
      implicit none

      include 'SIZE'
      include 'TSTEP'
      
      if (ctrlflag) then

         ! write i/o file
         if (NID.eq.0) then
            write(*,*) 'uio: ',uio(1:nuio)
            write(*,*) 'uio(:): ',uio
            write(ioctrl(1),ioformat) time,yio(1:nyio),uio(1:nuio),
     &                                dio(1:ndio),gio(1:ngio)
            flush(ioctrl(1))

         end if
      
      end if
      
      end subroutine ctrl_dump

c-----------------------------------------------------------------------

      subroutine ctrl_stop

      use ctrl, only: ctrlflag,cmpflag,
     &                ioctrl,fid
      implicit none

      include 'SIZE'
      
      if (ctrlflag) then
         
         if (NID.eq.0) then

            ! close outputs
            close(ioctrl(1))
c           close(iocltr(2))

            ! finalise compensator
            if (cmpflag) call cmp_stop(fid)
            
         end if

      end if
      
      end subroutine ctrl_stop

c-----------------------------------------------------------------------

      subroutine ctrl_forcing (fx,fy,fz,ix,iy,iz,iel)

      use ctrl, only: ctrlflag,
     &                nuio,ishuio,uio,
     &                shiox,shioy,shioz
      
      implicit none

      real,   intent(out) :: fx,fy,fz
      integer,intent(in)  :: ix,iy,iz,iel

      integer i

      fx = 0.0d0; fy = 0.0d0; fz = 0.0d0;

      if (ctrlflag) then

         do i=1,nuio
            fx = fx + uio(i) * shiox(ix,iy,iz,iel,ishuio(i))
            fy = fy + uio(i) * shioy(ix,iy,iz,iel,ishuio(i))
            fz = fz + uio(i) * shioz(ix,iy,iz,iel,ishuio(i))
         end do
      end if

      end subroutine ctrl_forcing

c-----------------------------------------------------------------------

      subroutine readnek (fname,rdx,rdy,rdz)

      implicit none
      
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      
      character*80,intent(in) :: fname
      real,intent(out)      :: rdx(lx1,ly1,lz1,lelv),
     &                         rdy(lx1,ly1,lz1,lelv),
     &                         rdz(lx1,ly1,lz1,lelv)

      real  dumx(lx1,ly1,lz1,lelv),
     &      dumy(lx1,ly1,lz1,lelv),
     &      dumz(lx1,ly1,lz1,lelv),
     &      dumt(lx1,ly1,lz1,lelv),
     &      dumxp(lx1*ly1*lz1*lelv,lpert),
     &      dumyp(lx1*ly1*lz1*lelv,lpert),
     &      dumzp(lx1*ly1*lz1*lelv,lpert),
     &      dumtp(lx1*ly1*lz1*lelv,ldimt,lpert)
      
      integer i,j

      call opcopy(dumx,dumy,dumz,vx,vy,vz)
      call copy(dumt,t,lx1*ly1*lz1*lelv)

      do j = 1,lpert

         call opcopy(dumxp(1,j),dumyp(1,j),dumzp(1,j),
     &               vxp(1,j),vyp(1,j),vzp(1,j))

         do i = 1,ldimt

            call copy(dumtp(1,i,j),tp(1,i,j),lx1*ly1*lz1*lelv)

         end do

      end do

      initc(1) = fname
      call setics
      
      call opcopy(rdx,rdy,rdz,vx,vy,vz)
      
      call opcopy(vx,vy,vz,dumx,dumy,dumz)
      call copy(t,dumt,lx1*ly1*lz1*lelv)

      do j = 1,lpert

         call opcopy(vxp(1,j),vyp(1,j),vzp(1,j),
     &               dumxp(1,j),dumyp(1,j),dumzp(1,j))

         do i = 1,ldimt

            call copy(tp(1,i,j),dumtp(1,i,j),lx1*ly1*lz1*lelv)

         end do

      end do
      
      end subroutine readnek



c-----------------------------------------------------------------------
      subroutine set_sensor(ishiox,iobj_s)  ! define objects for surface integrals
c
      include 'SIZE'
      include 'TOTAL'
c
      real*8 ishiox(lx1,ly1,lz1,lelv)
      integer e,f
      integer nxz
      real sumio
c
c     Define new objects for sensors
c
      nobj = nobj + 1 
      iobj = nobj
      
      if (maxobj.lt.nobj) write(6,*) 'increase maxobj in SIZEu. rm *.o'
      if (maxobj.lt.nobj) call exitt
c
      nxyz = nx1*ny1*nz1
      nxz = lx1*lz1
      do e=1,nelv
      do f=1,2*ndim
         if (cbc(f,e,1).eq.'W  ') then
         sumio = 0
         do j=1,nxyz ! check if the face has sensor support
          sumio = sumio + abs(ishiox(j,1,1,e))
         enddo
         if (sumio.gt.0) then
               nmember(iobj) = nmember(iobj) + 1
               mem = nmember(iobj)
               ieg = lglel(e)
               object(iobj,mem,1) = ieg
               object(iobj,mem,2) = f
              write(*,*) iobj,mem,f,ieg,e,nid,' SEN'
c              write(*,*) xm1(1,1,1,e), zm1(1,1,1,e)
c    1          format(6i9,a4)
         endif
c
         endif
      enddo
      enddo
      iobj_s = iobj
c     write(6,*) 'number',(nmember(k),k=1,4)
c
      return
      end
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine surface_int_J(sint,sarea,ishiox,e,f)
C
c     Copy of surface_int in navier5
c
      include 'SIZE'
      include 'GEOM'
      include 'PARALLEL'
      include 'TOPOL'
      include 'SOLN'
c      real a(lx1,ly1,lz1,1)
      real ishiox(lx1,ly1,lz1)
      real*8 duT(lx1,ly1,lz1), ut(lx1,ly1,lz1)
      real*8 dudx(lx1*ly1*lz1,3) ! field derivative
      
      real*8  uxp(lx1,ly1,lz1),uyp(lx1,ly1,lz1),
     $     uzp(lx1,ly1,lz1)
      
      real*8  ux_bc(lx1,ly1,lz1,lelv),uy_bc(lx1,ly1,lz1,lelv),
     $     uz_bc(lx1,ly1,lz1,lelv)
      COMMON / usrboundc / ux_bc,uy_bc,uz_bc
      
      real*8 sin_a(lx1,ly1,lz1), cos_a(lx1,ly1,lz1),
     $        msin_a(lx1,ly1,lz1)

      integer ix,iy,iz,e,f,nzx

c     Not quite sure about the index ordering here. TODO 
c     But far from the leading edge, the normal over the element
c     surface shouldnt change significantly.
      do iz=1,lz1
      do ix=1,lx1
        do iy=1,ly1
          sin_a(ix,iy,iz)= unx(ix,iz,f,e)
          cos_a(ix,iy,iz)= -uny(ix,iz,f,e)
          msin_a(ix,iy,iz)= -unx(ix,iz,f,e)
        enddo
      enddo
      enddo

c     Perturbation velocity
      call sub3(uxp,vx(1,1,1,e),ux_bc(1,1,1,e),lx1*ly1*lz1)
      call sub3(uyp,vy(1,1,1,e),uy_bc(1,1,1,e),lx1*ly1*lz1)
      call sub3(uzp,vz(1,1,1,e),uz_bc(1,1,1,e),lx1*ly1*lz1)
c      call copy(uxp,vx(1,1,1,e),lx1*ly1*lz1)
c      call copy(uyp,vy(1,1,1,e),lx1*ly1*lz1)
c      call copy(uzp,vz(1,1,1,e),lx1*ly1*lz1)


c      Compute Tangent velocity      
      call vdot2(ut,uxp,uyp,cos_a,sin_a,lx1*ly1*lz1)
c      Compute normal derivative
      call gradm1_elem(dudx(1,1),dudx(1,2),dudx(1,3),ut,e)
      call vdot2(duT,dudx(1,1),dudx(1,2),
     $  msin_a,cos_a,lx1*ly1*lz1)


      call dsset(lx1,ly1,lz1)

      iface  = eface1(f)
      js1    = skpdat(1,iface)
      jf1    = skpdat(2,iface)
      jskip1 = skpdat(3,iface)
      js2    = skpdat(4,iface)
      jf2    = skpdat(5,iface)
      jskip2 = skpdat(6,iface)
c      write(*,*) iface
      sarea = 0.
      sint  = 0.
      i     = 0

      do 100 j2=js2,jf2,jskip2
      do 100 j1=js1,jf1,jskip1
         i = i+1
         if(abs(ishiox(j1,j2,1)).gt.0) then
c         write(*,*) 'SIN: ', sin_a(j1,j2,1)
c         write(*,*) 'COS: ', cos_a(j1,j2,1)
         sarea = sarea+area(i,1,f,e)
         sint  = sint +area(i,1,f,e)*duT(j1,j2,1)
c         write(*,*) "AREA:   ",e ,area(i,1,f,e)
         endif
  100 continue

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine gradm1_elem(ux,uy,uz,u,e)
c     Copy of gradm1 in navier5 
c     Compute gradient of T -- mesh 1 to mesh 1 (vel. to vel.)
c
      include 'SIZE'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
c
      parameter (lxyz=lx1*ly1*lz1)
      real ux(lxyz),uy(lxyz),uz(lxyz),u(lxyz)

      common /ctmp1/ ur(lxyz),us(lxyz),ut(lxyz)

      integer e

      nxyz = lx1*ly1*lz1
c      ntot = nxyz*nelt

      N = lx1-1
c      do e=1,nelt
         if (if3d) then
            call elem_grad3(ur,us,ut,u,N,dxm1,dxtm1)
            do i=1,lxyz
               ux(i) = jacmi(i,e)*(ur(i)*rxm1(i,1,1,e)
     $                             + us(i)*sxm1(i,1,1,e)
     $                             + ut(i)*txm1(i,1,1,e) )
               uy(i) = jacmi(i,e)*(ur(i)*rym1(i,1,1,e)
     $                             + us(i)*sym1(i,1,1,e)
     $                             + ut(i)*tym1(i,1,1,e) )
               uz(i) = jacmi(i,e)*(ur(i)*rzm1(i,1,1,e)
     $                             + us(i)*szm1(i,1,1,e)
     $                             + ut(i)*tzm1(i,1,1,e) )
            enddo
         else
            if (ifaxis) call setaxdy (ifrzer(e))
            call elem_grad2(ur,us,u,N,dxm1,dytm1)
            do i=1,lxyz
               ux(i) =jacmi(i,e)*(ur(i)*rxm1(i,1,1,e)
     $                            + us(i)*sxm1(i,1,1,e) )
               uy(i) =jacmi(i,e)*(ur(i)*rym1(i,1,1,e)
     $                            + us(i)*sym1(i,1,1,e) )
            enddo
         endif
c      enddo
c
      return
      end
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine elem_grad3(ur,us,ut,u,N,D,Dt)
c     Copy of local_grad3 in navier5 
c     Output: ur,us,ut         Input:u,N,e,D,Dt
      real ur(0:N,0:N,0:N),us(0:N,0:N,0:N),ut(0:N,0:N,0:N)
      real u (0:N,0:N,0:N)
      real D (0:N,0:N),Dt(0:N,0:N)
      integer e
c
      m1 = N+1
      m2 = m1*m1
c
      call mxm(D ,m1,u(0,0,0),m1,ur,m2)
      do k=0,N
         call mxm(u(0,0,k),m1,Dt,m1,us(0,0,k),m1)
      enddo
      call mxm(u(0,0,0),m2,Dt,m1,ut,m1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine elem_grad2(ur,us,u,N,D,Dt)
c     Copy of local_grad2 in navier5 
c     Output: ur,us         Input:u,N,e,D,Dt
      real ur(0:N,0:N),us(0:N,0:N)
      real u (0:N,0:N)
      real D (0:N,0:N),Dt(0:N,0:N)
      integer e
c
      m1 = N+1
c
      call mxm(D ,m1,u(0,0),m1,ur,m1)
      call mxm(u(0,0),m1,Dt,m1,us,m1)
c
      return
      end
c-----------------------------------------------------------------------

