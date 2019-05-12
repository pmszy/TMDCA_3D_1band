c**************************************************************************** 
        subroutine readin
c**************************************************************************** 
	use Global
	implicit none
c**************************************************************************** 
        integer ios
        character*8 waslin
c****************************************************************************
c	Open the data to be outputed at the end of the simulation 
c**************************************************************************** 
c       open the output file
        open(unit=40,file='sigma.dat',status='unknown')
        open(unit=41,file='dos.dat',status='unknown')
        open(unit=42,file='run.dat',status='unknown')
        open(unit=43,file='chit.dat',status='unknown')
        open(unit=44,file='chieta.dat',status='unknown')
        open(unit=47,file='gamma_wt.dat',status='unknown')
        open(unit=48,file='sigma_only.dat',status='unknown')
c**************************************************************************** 
c       Open the input file
        open(unit=5,file='in1.dat',status='old',action='read')
c****************************************************************************
c       The cluster and Nc
c****************************************************************************
        read(5,"(1x,a8)") waslin
	read(5,*) cluster
	write(6,*) 'PROGRAM CTMT FOR CLUSTER ', cluster
c       Extract nc from the string cluster

	read(cluster,"(i5)",iostat=ios) Nc     ! should barf with characters
c       ibetts=0
        if (ios.gt.0) then
c         ibetts=1
          read(cluster,"(i4)",iostat=ios) Nc   ! ignore trailing charac.
          if (ios.gt.0) then
            read(cluster,"(i3)",iostat=ios) Nc ! ignore trailing charac.
            if (ios.gt.0) then
              read(cluster,"(i2)",iostat=ios) Nc ! ignore trailing charac.
              if (ios.gt.0) then
        	read(cluster,"(i1)",iostat=ios) Nc ! ignore trailing charac.
        	if (ios.gt.0) then
        	  write(lud,*) "Cannot read Nc from string ",cluster
        	endif
              endif
            endif
          endif
        endif
	write(6,*) "Nc=", Nc
c****************************************************************************
c       The other variables
c****************************************************************************
        write(6,*) 'enter meas,run,niter,nover'
        read(5,"(1x,a8)") waslin
        read(5,*) meas,run,niter,nover
        write(6,"(' meas=',i3, ' run=',i3,' niter= ',i3,' nover= ',i3)") meas,run,niter,nover
        if(meas.lt.1000) write(6,*) 'increase meas'
        if(run.lt.10) write(6,*) 'increase run'

        write(6,*) 'enter rand #,iso,dfac,ifac,isg'
        read(5,"(1x,a8)") waslin
        read(5,*) iseed,iso,dfac,ifac,isg
        write(6,"('iseed= ',i10,' iso=',i3,' dfac=',f7.4,' ifac= ',f7.4,' isg=',i3)") iseed,iso,dfac,ifac,isg
        isi=0
	if(ifac.gt.1.00001) isi=1
        isd=0
	if(abs(dfac-1.0).gt.0.01) isd=1

        write(6,*) 'enter ed,tprime,V'
        read(5,"(1x,a8)") waslin
        read(5,*) ed,tprime,V0
        write(6,"('ed = ',f9.6,' tprime=',f9.6,' V=',f9.6)") ed,tprime,V0
 	if(V0.ge.0) then
	  isv=1	! include diagonal disorder of strength V0.
	else
	  isv=-1	! The system is uniform (box) disorder
	  V0=-V0
	end if

	write(6,*) 'enter delta,nwn'
	read(5,"(1x,a8)") waslin
	read(5,*) delta,nwn
	write(6,"('delta= ',f9.6, ' nwn=',i6)") delta,nwn


	write(6,*) 'enter ca'
	read(5,"(1x,a8)") waslin
	read(5,*) ca
	write(6,"('ca= ',f9.6)") ca
        write(6,*)'************************************************'
        write(6,*)'Data Read Successfully from in1.dat'
        close(unit=5)

        if (ndim.eq.1) then
          N_sp = Nc/2 
          N_sp = 2 + N_sp !  for test!!!
        else if (ndim.eq.2) then
          N_sp=sqrt(float(Nc))
          N_sp=N_sp+2
        else if (ndim.eq.3.) then
          N_sp=sqrt(float(Nc))
          N_sp = 4 + N_sp !  for test!!!
        end if

!	read(5,*) V0
!       Read in some parameters
!	if(ndim==1) then
!	N_sp=Nc/2.
!	N_sp = N_sp+2.
!	elseif(ndim==2) then
!	N_sp=sqrt(float(Nc))
!	N_sp=N_sp+2
!	elseif(ndim==3) then
!	N_sp=sqrt(float(Nc))
!	N_sp=N_sp+2
!	end if

c 	if(delta.lt.1.0e-5) then
c	  isw=1
c	else
c	  isw=0
c	end if

c   set up the random number 
c       call ranset(iseed)	!CRAY
c	call srand(iseed)	!HP
c	r1=rand(iseed)	!g77
c	r1=ran(iseed)	!compaq alpha/linux

c	set some parameters

c	iter=0

          write(lud,*)   '******************************************************************************'
          write(lud,*)   '* SET SWITCHES FOR EASY IDENTIFICATION:  				        *'
	  call userinfo(lud)
          write(lud,*)    'The Dimension of the calculation is', ndim
	  select case(disorder_type)
	    case(1)	! Binary distrubution
	  write(lud,*) '################################################################################'
	  write(lud,*) 'Binary distribution is used'
	  write(lud,*) '################################################################################'
	  case(2)	!Box distribution
	  write(lud,*) '################################################################################'
	  write(lud,*) 'Box distribution is used'
	  write(lud,*) '################################################################################'
	case(3) ! Gaussian 
	  write(lud,*) '################################################################################'
	  write(lud,*) 'Gaussian(Normal) distribution is used'
	  write(lud,*) '################################################################################'
	case(4)	! Cauchy or Lorentzian distribution
	  write(lud,*) '################################################################################'
	  write(lud,*) 'Cauchy distribution is used'
	  write(lud,*) '################################################################################'
	case(5)	! Lognormal distribution
	  write(lud,*) '################################################################################'
	  write(lud,*) 'Lognormal distribution is used'
	  write(lud,*) '################################################################################'
    	     case default
             print*, "Please, select appropriate distribution function"
             stop
              end select
	   select case(tables_type)
          case(1)
            write(lud,*) '* Bettis Lattice Used                   				      *'
          case(2)
            write(lud,*) '* Conventional Lattice Used               				       *'
          case(3)
            write(lud,*) '* Cubic Lattice Used    	             				       *'
    	     case default
             print*, "Please, choose between case 1-3. See module_global.for for details"
             stop
              end select   
          if(iso.eq.0) then
            write(lud,*) '* Initialized with Guessed Hybridization Rate 				    *'
          else if(iso.eq.2) then
            write(lud,*) '* Initialized with Supplied Potential         				*'
          end if
	  if(iso.ge.1.) then
            write(lud,*) '* Initialized with Supplied Sigma             				*'
	end if
          if(tmt)  then
            write(lud,*) '*CTMT used: rho_typ(K) = exp(<Ln(rho(K)>)	 				*'
       else
          if(tmt_new)  then
            if(tmt1)  then
            write(lud,*) '*TMDCA used: exp(<Ln(rho_ii)>)exp(<Ln(rho(K)/rho_ii)>)	 		*'
	 end if
            if(tmt2)  then
            write(lud,*) '*TMDCA used: rho_typ(K) = exp(<Ln(rho_ii)>)<(rho(K)/rho_ii)>	             *'
	  end if
            if(tmt3)  then
            write(lud,*) '*TMDCA used: rho_typ(K) = (1-a)*exp(<Ln(rho_ii) + a*<rho(K)	 		*'
	end if
            if(tmt4)  then
            write(lud,*) '*TMDCA used: rho_typ(K) = exp(<Ln(rho_ii)>)	 				*'
	  end if
	else
            write(lud,*) '*Code using the DCA: rho_typ(K) = rho(K)	 				*'
	end if !tmt new
	end if ! if tmt

c	Print date and time to keep notice when calculation is done
          write(lud,*)   '*-----------------------------------------------------------------------------*'
          write(lud,*)   'Calculation Started at the Time                                         *'
              call date_and_time(date,time,zone,values)
              call date_and_time(DATE=date,ZONE=zone)
              call date_and_time(TIME=time)
              call date_and_time(VALUES=values)
	      write(lud,*) '*DATE(Y-M-D)  TIME(H-M-S.MS)   TIME-ZONE      			       	   *'                    
              write(lud, '(1x,a,3x,a,3x,a)') date, time, zone 
	      write(lud,*) '*Y  M   D T-DIFF(UTC) H MIN SEC MILLIS         				*'                    
              write(lud,'(8i5)') values     
          write(lud,*)   '******************************************************************************'

        return
        end subroutine readin

        subroutine init_parameters
c************************************************************
c       Initialize IPT in preparation for the first iteration
c************************************************************
        use Global
        implicit none
c************************************************************
        integer n,m,i,j,jc,ic
        real(kind) ::  r1,r2,r3
	complex*16 :: z1
        character*72 linemc
c************************************************************
				
	  if(iso.ge.1.) then ! readin sigma from a file	
									
 2002       read(40,"(a72)") linemc				      
            if(index(linemc(1:72),'ed,V').ne.0 ) then	      
	      read(40,*) ed					      
            else						      
              goto 2002 					      
            endif						      
	        						      	        						      
 2004       read(40,"(a72)") linemc				      
            if(index(linemc(1:72),'sigma').ne.0 ) then  	      
	      do ic=1,Nc 					      
	      do n=-nwn,nwn				      
	  	read(40,*,end=997) m,m,r1,r2,r3
		sigma(ic,n)=cmplx(r2,r3)
              end do						      
              end do
 	      write(6,*) '****************************'						      
	      write(6,*) 'code initialized with supplied sigma'
 	      write(42,*) '****************************'      
	      write(42,*) 'code init. with supplied sigma'
	
          else if(iso.ge.2.and.				      
     &          index(linemc(1:72),'potentials').ne.0 ) then	      
	     do i=1,Nc 					      
	       read(40,*) j,V(i)				      
	     end do	
 	     write(6,*) '***********************************'						      
	     write(6,*) 'code initialized with supplied potentials'	     
 	     write(42,*) '***********************************'
	     write(42,*) 'code init. with supplied potentials'	      
	     isf=1						      
	     goto 2004 					      
          else if(index(linemc(1:72),'G(1,w)').ne.0 ) then  	      
	     do n=-nwn,nwn
	       read(40,*) r1,r2,r3
	       gft(n)=cmplx(r2,r3)
	     end do
	     goto 2004 					      
          else if(iso.eq.3.and.				      
     &        index(linemc(1:72),'Chiloc').ne.0 ) then  	      
	    do n=1,neta
	      read(40,*) r1,r2
	      chit(n)=r2
	    end do
	    goto 2004 					      
	  else						      
            goto 2004 					      
          endif						      
 997	  continue
		call init_coarsegraining
            write(6,*) '**********************************'
            write(6,*) 'Code initialized with supplied Self-Energy'
            write(lud,*) '*********************************'
            write(lud,*) 'Code initialized with supplied Self-Energy'	
  
	  else if(iso.eq.0.) then ! set sigma to a constant
	call init_coarsegraining   


cc	do n=-nwn+1,nwn-1
cc	   do m=-nwn,n-1
cc		G_typ_real(:,n) = G_typ_real(:,n) +  p_rho_GcKf(:,m)*dwn(m)/(wn(n)-wn(m))
cc	   enddo
cc	   do m=n+1,nwn
cc		G_typ_real(:,n) = G_typ_real(:,n) +  p_rho_GcKf(:,m)*dwn(m)/(wn(n)-wn(m))
cc	   enddo
cc	enddo
cc	Now load the real and imaginary part into G_typ
cc	do n=-nwn,nwn
cc	  GcKf(:,n)= G_typ_real(:,n) - ii*pi*(p_rho_GcKf(:,n))
cc	      GcKfsc(:,n)=GcKf(:,n)
cc	end do
c	    write(6,*) 'Code initialized with sigma=-0.01*i'
c	    write(42,*) 'Code initialized with sigma=-0.01*i'


 	    write(6,*) '**********************************'						      
	    write(6,*) 'code initialized with guess Hybridization Rate'	     
 	    write(42,*) '*********************************'						      
	    write(42,*) 'code initialized with guess Hybridization Rate'	 

	    end if
	return
	end subroutine init_parameters
c**************************************************************************** 
        subroutine Bin_Header
c****************************************************************************
	use Global
	implicit none
c**************************************************************************** 
       integer rmeas,n
	  write(42,*) 'cluster type '
	  write(42,*) cluster
          write(42,*) 'Results from CTMT for nwn,Nc='
          write(42,*) nwn,Nc
          write(42,*) 'parameters follow'
          write(42,*) 'meas,run,niter,dfac,ifac'
          write(42,"(i5,2x,i5,2x,i5,1x,f5.2,2x,f5.2)") meas,run,niter,dfac,ifac
	  rmeas=meas
	  do n=1,niter-1
	    rmeas=int(ifac*rmeas)
 	  end do
	  write(42,*) 'random number, iso'
	  write(42,*) iseed,iso
          write(42,*) '  ed    ,  tprime,    V    '
          write(42,"(1x,f9.4,1x,f9.4,1x,f9.4,1x,f9.4,1x,f9.4)") ed,tprime,V0
          write(42,*) ' '
	  return
	end subroutine Bin_Header

c**************************************************************************** 
        subroutine init_coarsegraining
c****************************************************************************
	use Global
	implicit none
c****************************************************************************
        integer n,m,i,j,jc,ic
        real(kind) ::  r1,r2,r3,z1r,z1i,z2r,z2i
        complex*16 :: z1,z2
c****************************************************************************

	  GcKf = zeroc
	  GcKfsc = zeroc
	  GcKfsc_av = zeroc
	  Gamma_old = zeroc
	  Gamma_new = zeroc
	  Gammarun_old = zeroc
	  Gammarun_new = zeroc
c       Now on every processor, pass the information around
	if(iso.eq.0) then
	    do jc=1,Nc 
	      do n=-nwn,nwn
	      sigma(jc,n)=dcmplx(0.d0,-0.04D0) ! Use large self-energy to initialize instead of eta
	      end do
	    end do
	else
	end if

        r3=0.d0
	do n=-nwn,nwn
	  do ic=1,Nc
	    z1 = wn(n)-ed-sigma(ic,n)-ed
	    z1r= dreal(z1) 
	    z1i= dimag(z1) 
	    do i=0,e_div
	      r1=estart+(dfloat(i)-0.50D0)*estep 
	      r2=estart+(dfloat(i)+0.50D0)*estep
	      z2r=-0.50D0*dlog(((z1r-r2)**2 + (z1i**2))/
     .     	          ((z1r-r1)**2 + z1i**2))
	      z2i=-(datan((r2-z1r)/dabs(z1i))-datan((r1-z1r)/dabs(z1i)))
              GcKf(ic,n)=GcKf(ic,n)+ PDOS(i,ic)*(z2r+ii*z2i)
	    enddo
	    GcKf(ic,n) = GcKf(ic,n)*dfloat(Nc)
!            r3=r3+(-dimag(GcKf(ic,n))/pi)*dwn(n)
	  end do
	end do
!        write(6,*) 'Norm in readin of GcKf=',r3

        r3=0.d0
	do jc=1,Nc 
	  do n=-nwn,nwn
                    GcKfsc(jc,n)=1.0D0/(sigma(jc,n)+ 1.0D0/GcKf(jc,n)) 
	  end do  
	end do


	  do jc=1,Nc 
	    do n=-nwn,nwn
	    Gamma_old(jc,n) = wn(n)-Epsbar(jc)-1.0D0/GcKfsc(jc,n)
	    Gammarun_old(jc,n) = wn(n)-Epsbar(jc)-1.0D0/GcKfsc(jc,n)
	  end do  
	end do
	end subroutine init_coarsegraining
c	This writes out user info and time
c****************************************************************************
	subroutine userinfo(fh_info)
	use Global
	implicit none
c****************************************************************************

	  integer fh_info
	  integer hostnm,ihost
	  character*24 fdate
	  character*24 getlog
	  character*30 HOSTNAME,hostnam
	  real*8 dtime,tarray(2)
	  real*8 exec_time,user_time
!	 ihost = hostnm(HOSTNAME)
!	CALL getenv("HOST",HOSTNAME)
!	   CALL GET_ENVIRONMENT_VARIABLE("hostname",HOSTNAME) 
	  write(fh_info,*)
	  write(fh_info,*) "USR_INFO: Date : ",fdate()
	  write(fh_info,*) "          User : ",getlog()
	  write(fh_info,*) "          Host : ",HOSTNAM(hostname)


	  exec_time = dtime(tarray)

	  return

	  entry exitinfo(fh_info)

	  exec_time = dtime(tarray)
	  user_time = tarray(1)

	  write(fh_info,*)
	  write(fh_info,*) "USR_INFO: Time of computation"
	  write(fh_info,'(11x,"Total run-time : ",f10.3,1x,"seconds")') exec_time
	  write(fh_info,'(11x,"Total cpu-time : ",f10.3,1x,"seconds")') user_time

	  return
	end subroutine userinfo
