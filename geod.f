c***************************************************************************************
        subroutine geodgen				
c       This block makes the initial greens function.  This
c	function corresponds to the field configuration where all 
c 	s(i)=0.0.  It does this by summing the Fourier transform. 
c	The coarse-graining is done in parallel.  
 
c***************************************************************************************
	use Global
	use omp_lib
	implicit none
c***************************************************************************************
c	integer,parameter:: nt_broyden = (2*5+1)*20
        integer i,n,ic,jc,info
        real(kind) :: r1,r2
        complex(kind) :: csum,G_av_temp(Nc,Nc,-nwn:nwn) 
       call profstart('geodgen')
	time1= secnds(0.0)
c	wtime= OMP_get_wtime()
c	time1= secnds(0.0)
c***************************************************************************************
c	First, we need a map from the process rank to the cluster site
c	indices.  Each process gets its K-value, where its is determined by
c	calculate the initial reciprocal-space green function Gckf and GcKfsc.
c       Symmetrize the GcKfsc and  x <---> +/-x       
c       hence the self energy by      \ /	      
c       imposing the point group       X			    
c       symmetries.		      / \			    
c       			   y <---> +/-y 
c***************************************************************************************

	call Symm_PG_Complex(GcKfsc)  
 

                       
c	do ic=1,Nc
c	write(808,*) real(GcKfsc(ic,0)),aimag(GcKfsc(ic,0))
c	end do
c	stop
c***************************************************************************************                     
c	Now transform the initial greens function to cluster space  
c	for Gcsfsc.
c***************************************************************************************

	do ic=1,Nc
           do n=-nwn,nwn
	    csum=zeroc
	     do jc=1,Nc
 	    csum=csum+GcKfsc(jc,n)*FTCoefs_K_to_R(jc,ic)			  
 	   end do
          Gcsfsc(ic,n)=csum
	  end do
	end do

c***************************************************************************************
c	form G_script inverse
c***************************************************************************************
        do n=-nwn,nwn
        do ic=1,Nc
        do jc=1,Nc
          Gsc_temp(ic,jc,n)=Gcsfsc(icrdiff(ic,jc),n)
	end do
	end do
	end do

c***************************************************************************************
cc	For G_script to be used in dateregup.f to calculate 1/G = 1/Gs + V
cc	Factorization of a complex general matrix (Needed before inversion) (Lapack)
c***************************************************************************************
c	if(n.eq.0) then
c	do ic=1,Nc
c	do jc=1,Nc
c	write(808,*) real(Gsc_temp(ic,jc)), aimag(Gsc_temp(ic,jc))
c	end do
c	end do
c	stop
c	end if
	do n=-nwn,nwn  ! Frequency loop starts
	call zgetrf(Nc,Nc,Gsc_temp(1,1,n),Nc,ipvt,info)
cc 	Matrix inversion
	if (info.eq. 0.d0)then
	call zgetri(Nc,Gsc_temp(1,1,n),Nc,ipvt,work,lwork,info)
      else
           print*, 'info=/0 Matrix Inversion Failed in geod: Info=',info
	  stop
         end if

	do ic =1,Nc
	   do jc = 1,Nc
	    Gcinv(ic,jc,n) = Gsc_temp(ic,jc,n)
	   end do
	 end do
	end do ! Frequency loop ends


c***************************************************************************************
c	Initialize some variables and accumulators
c***************************************************************************************
	  gf=zeroc
	  GcKf=zeroc
	  Gdff=zeroc
	  G_typ=zeroc
	  GcKf_av=zeroc
	  rho=zeror; rho_typ=zeror
	  rhorun=zeror; logrhorun=zeror; logrholoc_run=zeror
	  rho_logf=zeror; rho_logfs=zeror
	  rho_local=zeror; rho_typ_local = zeror
          rhof=zeror; rhofs=zeror
          rholoc_f=zeror; rholoc_fs=zeror
	  logrholoc_f=zeror
	  rholoc_run=zeror
	  rholoc_av=zeror
	  chimeas = zeror
	  chit = zeror
	  gammaf=zeror;gammafs=zeror
	  Gtyprun=zeroc
	  Gbarrun=zeroc
	  Gscriptrun=zeroc
	  chicharge_meas=zeroc
	  G_chg=zeroc

          rho_meas=zeror
          logrho_meas=zeror
          logrholoc_meas=zeror
	  rholoc_meas=zeror






c	time2= secnds(time1) 
	wtime= OMP_get_wtime()-wtime
c	write(6,*) 'time in geod= ',time2

       call profend('geodgen')
 	time2= secnds(time1)
	write(6,*) 'time in geod= ',time2
        geod_time=geod_time+time2
        return
        end


!**********************************************************************
	    subroutine Symm_PG_Real(Gc)
!**********************************************************************
!       Average the green funtion Gc over all allowed point-group
!       operations.
!                         x <---> +/-x
!                            \ /
!                             X
!                            / \
!                         y <---> +/-y
!**********************************************************************
	    use Global
	    implicit none
!**************************************************************
	    integer :: n,ic
	    real(kind) :: Gc(Nc,-nwn:nwn)
!**************************************************************
	    do n=-nwn,nwn
		do ic=1,Nc
		    data_rho_typ(ic,n)=zeror
		    do igroup=1,ngroup
		        data_rho_typ(ic,n)=data_rho_typ(ic,n) + 
     &		        (1.d0/real(ngroup))*Gc(ickequ(ic,igroup),n)
		    end do
		    Gc(ic,n)=data_rho_typ(ic,n)
		end do
	    end do
		                                                       

	    return
	    end subroutine Symm_PG_Real
		    
		    

	    subroutine Symm_PG_Complex(Gc)
!**********************************************************************
!       Average the green funtion Gc over all allowed point-group
!       operations.
!                         x <---> +/-x
!                            \ /
!                             X
!                            / \
!                         y <---> +/-y
!**********************************************************************
	    use Global
	    implicit none
!**************************************************************
	    integer :: n,ic
	    complex(kind) :: Gc(Nc,-nwn:nwn)
!**************************************************************
	    do n=-nwn,nwn
		do ic=1,Nc
		    chi(ic,n)=zeroc
		    do igroup=1,ngroup
		        chi(ic,n)=chi(ic,n) + 
     &		        (1.d0/real(ngroup))*Gc(ickequ(ic,igroup),n)
		    end do
		    chi(ic,n)=Gc(ic,n)
		end do
	    end do
		   
	Gc(1:Nc,-nwn:nwn)=chi(1:Nc,-nwn:nwn)
	                                                    

	    return
	    end subroutine Symm_PG_Complex
      
!**********************************************************************
	    subroutine Symm_PH_Real(Gc)
!**********************************************************************
!       Impose particle-hole symmetry

!       G (iwn) = -G   (-iwn)
!        k          Q-k

!       H is particle-hole symmetric whenever ed=0 and t'=0.  For a
!       two-dimensional system, this means that

!               +             +    i
!       C  --> C exp(iQ.i) = C (-1)
!        i      i             i
!                                                +              +
!       C  = SUM exp(ik.i) C  --> SUM exp(ik.i) C exp(-iQ.i) = C
!        k    i             i                    i              Q-k

!       or if ed=0, G (tau) = -G   (-tau) = G   (beta-tau)
!                    k          Q-k          Q-k
!       and
!                   G (iwn) = - G   (-iwn)  real parts opposite sign
!                    k           Q-k        imaginary parts equal sign

!**********************************************************************
	    use Global
	    implicit none
!**************************************************************
	    integer :: n,ic
	    real(kind) :: Gc(Nc,-nwn:nwn)
!**************************************************************
!       Impose particle-hole symmetry

!       G (iwn) = -G   (-iwn)
!        k          Q-k



	    data_rho_typ(:,:)=zeror
	    do n=-nwn,nwn
		do ic=1,Nc
		    if(gsflag == 0) then
			if(Nc>1) then
		        data_rho_typ(ic,n)=data_rho_typ(ic,n)+halfr*(Gc(ic,n)+ 
     &	 	        Gc(ickdiff(icKpi,ic),-n))
		    else
		        data_rho_typ(ic,n)=data_rho_typ(ic,n)+halfr*(Gc(ic,n)+Gc(ic,-n))
		     end if
		    elseif(gsflag == 1) then
		        data_rho_typ(ic,n)=Gc(ic,n)
		    end if
		end do
	    end do
	    Gc(1:Nc,-nwn:nwn)=data_rho_typ(1:Nc,-nwn:nwn)
	    return
	    end subroutine Symm_PH_Real


!**********************************************************************
	    subroutine Symm_PH_Local(Gc)
!**********************************************************************
!       Impose particle-hole symmetry on local rhos

!**********************************************************************
	    use Global
	    implicit none
!**************************************************************
	    integer :: n,ic
	    real(kind) :: Gc(-nwn:nwn)
!**************************************************************
!       Impose particle-hole symmetry

!       G (iwn) = -G   (-iwn)

	    rho_geo(:)=zeror
	    do n=-nwn,nwn
		if(gsflag == 0) then
		    rho_geo(n)=halfr*(Gc(n)+ Gc(-n))
		elseif(gsflag == 1) then
		    rho_geo(n)=Gc(n)
		end if
	    end do
	   Gc(-nwn:nwn)=rho_geo(-nwn:nwn)
	    return
	    end subroutine Symm_PH_Local


!**********************************************************************
	    subroutine Symm_PH_Complex(Gc)
!**********************************************************************
!       Impose particle-hole symmetry

!       G (iwn) = -G   (-iwn)
!        k          Q-k

!       H is particle-hole symmetric whenever ed=0 and t'=0.  For a
!       two-dimensional system, this means that

!               +             +    i
!       C  --> C exp(iQ.i) = C (-1)
!        i      i             i
!                                                +              +
!       C  = SUM exp(ik.i) C  --> SUM exp(ik.i) C exp(-iQ.i) = C
!        k    i             i                    i              Q-k

!       or if ed=0, G (tau) = -G   (-tau) = G   (beta-tau)
!                    k          Q-k          Q-k
!       and
!                   G (iwn) = - G   (-iwn)  real parts opposite sign
!                    k           Q-k        imaginary parts equal sign

!**********************************************************************
	    use Global
	    implicit none
!**************************************************************
	    integer :: n,ic
	    complex(kind) :: Gc(Nc,-nwn:nwn)
!**************************************************************
!       Impose particle-hole symmetry

!       G (iwn) = -G   (-iwn)
!        k          Q-k


	    chi(:,:)=zeroc
	    do n=-nwn,nwn
		do ic=1,Nc
		    if(gsflag == 0) then
		        chi(ic,n)=chi(ic,n)+halfr*(Gc(ic,n)+ 
     &		        Gc(ickdiff(icKpi,ic),-n))
		    elseif(gsflag == 1) then
		        chi(ic,n)=Gc(ic,n)
		    end if
		end do
	    end do
		                    
	    Gc(1:Nc,-nwn:nwn)=chi(1:Nc,-nwn:nwn)
	    return
	    end subroutine Symm_PH_Complex


!**************************************************************
	    subroutine HilbertT(Gt,Gtr,Rhot)
!**********************************************************************
!c	Do Hilbert transformation to calculated the Typical Green function (TGF).
!       First the real part.  Note, this is a principal value integral
!       call profstart('HilbertT_Sumup') ! To calculate the time spent in HT.

!		      oo
!		       /
!	G_typ(k,w)= P | dw' rho_typ(k,w') -i*pi*rho_typ(k,w)
!		       /     -----------
!	             -oo       w - w'
!**********************************************************************
	    use Global
	    implicit none
!**************************************************************
	    integer :: n,m,ic
	    real(kind) :: Gtr(Nc,-nwn:nwn),Rhot(Nc,-nwn:nwn)
	    complex(kind) :: Gt(Nc,-nwn:nwn),sigmacg(Nc,-nwn:nwn)
!**************************************************************

	    Gtr = zeror
	    do ic=1,Nc
		do n=-nwn,nwn
		    do m=-nwn,n-1
		        Gtr(ic,n) = Gtr(ic,n) + 
     &		        Rhot(ic,m)*dwn(m)/(wn(n)-wn(m))
		    enddo
		    do m=n+1,nwn
		        Gtr(ic,n) = Gtr(ic,n) + 
     &		        Rhot(ic,m)*dwn(m)/(wn(n)-wn(m))
		    enddo
		enddo
	    enddo

!c	Now load the real and imaginary part into G_typ
	    Gt = zeroc
	    do ic=1,Nc
		do n=-nwn,nwn
		    Gt(ic,n)= Gtr(ic,n) - ii*pi*(Rhot(ic,n))
		end do
	    enddo

	    end subroutine HilbertT
