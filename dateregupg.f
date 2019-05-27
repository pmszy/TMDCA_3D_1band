c**************************************************************************** 
        subroutine update
c**************************************************************************** 
c       Update generates a new random configuration of the disorder potentials
c       and calculates the corressponding cluster green function.
 
c**************************************************************************** 
	use Global
	use omp_lib
	implicit none
c**************************************************************************** 
 
        integer jc
        real(kind) :: Vt,start,end1,ran1,ran2
c**************************************************************************** 
	time1= secnds(0.0)
        call profstart('dateregupg')
        select case(disorder_type) 
	  case(1)	! Binary distrubution P(rand()) = ca*delta(rand())+(1-ca)*delta(rand())
	    do jc=1,Nc
    !           Now update the jc th potential V(jc)
    !           if(rand(0).gt.0.5)then 	!g77
    !           Propose a change
		if(isv == 1) then
		    V(jc)=V0*sign(1.0,ca-ran(iseed)) 
	    endif
	    enddo
	  case(2)	!Box distribution P(rand()) = 1/2*W(step_function(W-rand())
	    do jc=1,Nc
		    Vt= 2.d0*V0*(ran(iseed)-0.50D0)! 
		if(ran(iseed) > 0.5 .OR. isv == -1)then
		    V(jc)= Vt
		end if
	    enddo
	  case(3) ! Gaussian or Normal. 
	V0=dabs(V0)
	    do jc=1,Nc
	      V(jc) = V0*dsqrt( -log(ran(iseed))/6.0D0)*dcos( 2.0D0*pi*ran(iseed) )
	    enddo

	
	  case(4)	! Cauchy or Lorentzian distribution
	    do jc=1,Nc
	       V(jc)=V0*dtan(pi*(ran(iseed)-0.50D0))
	     
	    enddo
	  case(5)	! Lognormal distribution
	    V0=dabs(V0)
	    do jc=1,Nc  
	      V(jc)=dexp(-dsqrt(2.0D0)*V0*(1.0D0/derfc(2.0D0*ran(iseed)-0.50D0)/6. ))! Not yet working
	    enddo

    	  case default
            print*, "Please, select appropriate distribution function"
            stop
        end select
 
        call reg
        call profend('dateregupg')
 	time2= secnds(time1)
        update_time=time2
        return
        end
c**************************************************************************** 
        subroutine reg
c       This subroutine reconfigures G given only the potential configuration
c       and the initial green's functions Gcsfsc (which corresponds to
c	all V=0).
 
c**************************************************************************** 
	use Global
	use omp_lib
	implicit none
c**************************************************************************** 
 
        integer n,ic,jc,info
 
c       Gcsfsc is the greens function for all V=0, hence we will
c       sweep through the lattice and reconfigure Gcsfsc to coorespond
c       to the lattice configuration.
c**************************************************************************** 

	V_temp = zeroc
        do jc=1,Nc
	    V_temp(jc,jc) = cmplx(V(jc),0.d0)
	end do

c*************************************************************************************************************************
ccc	In this section, we will use lapack subroutines to invert a matrix, we will do this via LU decomposition 
ccc	using zgetri which should be called after calling the factorization subroutine zgetrf. Note, these subroutines 
ccc	are for double precision for single precision, change kind in module to 4 and the subroutines as cgetri and cgetrf, 
ccc	respectively. C. E. Ekuma: 09-02-12.
c*************************************************************************************************************************
ccc	Perform the operation 1/G = 1/Gs - V
	Ginv = Gcinv 
!$OMP PARALLEL
!$OMP DO SCHEDULE(STATIC,10)
	do n=-nwn,nwn
	do jc =1,Nc
	Ginv(jc,jc,n) = Gcinv(jc,jc,n)-V_temp(jc,jc)
	end do
	end do
    !$OMP END DO NOWAIT  
!$OMP END PARALLEL     	



!$OMP PARALLEL
!$OMP DO ORDERED SCHEDULE(GUIDED,4)
	do n=-nwn,nwn  ! Frequency loop starts
cc	Factorization of a complex general matrix (Needed before inversion) (Lapack) 
!$omp ordered
	  call zgetrf(Nc,Nc,Ginv(1,1,n),Nc,ipvt,info)
cc 	Matrix inversion
       	  if (info.eq. 0.d0)then

	     call zgetri(Nc,Ginv(1,1,n),Nc,ipvt,work,lwork,info)
      else
           print*, 'info=/0 Matrix Inversion Failed in update: Info=',info
	  stop
         end if
    !$OMP END ORDERED 
c       move G_temp to G
	    do ic=1,Nc
	     do jc=1,Nc
		G(ic,jc,n) = Ginv(ic,jc,n)
	       end do
	      end do
	end do ! Frequency loop ends
    !$OMP END DO NOWAIT
!$OMP END PARALLEL

        return
        end
	

