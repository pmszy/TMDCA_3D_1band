c***************************************************************************
        subroutine analysis
c***************************************************************************
c       Compute vertex function and two particle quantities
c***************************************************************************
	use Global
	use omp_lib
c***************************************************************************
 	implicit none
        integer i,n,m,p,ic,jc,jcr,jck
	complex*16 :: c1, csum(Nc,-nwn:nwn),Gdkk_av(Nc,-nwn:nwn)
	real*8 :: r1,r2,r3,r4,histogram(Nc,-nwn:nwn)
c***************************************************************************
	time1= secnds(0.0)
        call profstart('analysis')
	wtime= OMP_get_wtime()
c***************************************************************************	
c  	  The particle-hole response needed to calculate
c  	  the conductivity was computed in measure.  We will now compute
c         the bare bubble, chi0, and extract the vertex function from the
c         equation below.
c
c  	  K,w	    K',w	K,w	     K,w	     P,w      K',w
c  	    ----<----	     ----<----      ----<----|-------|----<----
c  	      / / /				     |       |  / / /
c  	      /Chi/	 =   delta	 +	     | Gamma |  /Chi/
c  	      / / /		  K,K'  	     |       |  / / /
c  	    ---->----	     ---->----      ---->----|-------|---->----
c  	  K,w+w'    K',w+w'	K,w+w'       K,w+w'    P,w+w'	      K',w+w'
c
c
c  	  Only K,K',P are matrix labels, all of the diagrams, at the
c  	  two-particle level, carry the same frequency labels, so for each
c  	  w and w', we solve the matrix equation for Gamma
c
c  	  Chi(K,K',w,w') = G(K,w)G(K,w+w') +
c  			   G(K,w)G(K,w+w')Gamma(K,P,w,w')Chi(P,K',w,w')
c
c
c
c
c  	  Chi0(K,K',w,w') = \delta     Gc(K,K,w) Gc(K,K,w+w') 
c                                 K,K'
c	
c         To calculate the DC limit (w'-->0) we will set w'=wn(1)

c       
c***************************************************************************************
c        Read-in the chicharge_meas and Gsc_temp that were computed in
c        meas.f

c***************************************************************************************
c         Now calculate
c
c  	  Chi0(K,K,w,w') = Gc(K,K,w) Gc(K,K,w+w')
c
c              K      w     K'     ick   <-> K'
c         or   -------<-------     jck   <-> K
c              ------->-------     wn(1) <-> w'
c              K     w+w'   K'
c
c         We will use the time-ordered Fourier transformed cluster
c         Green's function (named Gsc_temp) computed in meas.f

          chicharge_meas0=dcmplx(0.d0,0.d0)
          do ick=1,Nc
          do n=-nwn,nwn-1
	    chicharge_meas0(ick,ick,n) = Gsc_temp(ick,ick,n)*Gsc_temp(ick,ick,n+1)
	  end do
	  end do
c***************************************************************************************
c       Now invert the chi(K,K',w,w') and chi0(K,K',w,w') for each (w); (w'=w+nu)
c       and nu is the smallest non-zero frequency, wn(1)
c       Then get the vertex function

        do iw=-nwn,nwn-1

          do ick=1,Nc
            do jck=1,Nc
              chi_inv0(ick,jck)=dcmplx(0.d0,0.d0)
              if(ick.eq.jck) chi_inv0(ick,jck)=1.d0/chicharge_meas0(ick,ick,iw)
              chi_inv(ick,jck)=chicharge_meas(ick,jck,iw)
            end do
          end do


c	  Factorization of a complex general matrix (Needed before inversion) (Lapack) 
	  call zgetrf(Nc,Nc,chi_inv(1,1),Nc,ipvt,info)
c 	  Matrix inversion
       	  if (info.eq. 0.d0)then
	     call zgetri(Nc,chi_inv(1,1,n),Nc,ipvt,work,lwork,info)
          else
             print*, 'info=/0 Matrix Inversion Failed in update: Info=',info
	     stop
          end if

          do ick=1,Nc
            do jck=1,Nc
              vertexf(ick,jck,iw)=chi_inv0(ick,jck) - chi_inv(ick,jck)
            end do
          end do


        end do


c***************************************************************************************
c       Now that we have the vertex function, we can go ahead and
c       compute the conductivity.
c       
c            
c                             _________________________      
c       Chi_cl_CG(K,K',w,w') = V_kx Chi_cl(k,k',w,w') V_k'x  
c
c                             __________________________
c       Chi_cl0_CG(K,K',w,w')= V_kx Chi_cl0(k,k',w,w') V_k'x 
c
c                                  _____________________            
c       Chi_cl0_CGL(K,K',w,w') =   V_kx Chi_cl0(k,k',w,w')
c
c                                 _____________________
c       Chi_cl0_CGR(K',K,w,w') =  Chi_cl0(k',k,w,w') V_kx
c
c        
c       Chi_cl_CG(K,K',w,w') = Chi_cl0_CG(K,K',w,w') 
c
c        _______                  +
c        \
c         \    
c         /      Chi_cl0_CGL(K,K',w,w') Fl(K',K",w,w') Chi_cl0_CGR(K",K,w,w') 
c        /
c        _______
c         K',K"
c       
c                                           Vertexf(K,K',w,w')
c        Fl(K,K',w,w')=   ---------------------------------------------------           
c                         1    -    Vertexf(K,K',w,w') Chi_cl0_CG(K,K',w,w') 
c
c
c       The lattice Green's functions obtained in spectra_standalone.f
c       will be time-ordered and then used to calculate the
c       Chi_cl0(k,k',w,w'). Then the three Coarse Grained quantities,
c       Chi_cl0_CG, Chi_cl0_CGL, Chi_cl0_CGR will be obtained.
c**************************************************************************************
         
        call profend('analysis')
 	time2= secnds(time1)
        analysis_time=analysis_time+time2
        return
        end
	
        Subroutine Invert_Nc(mat)
        use Global
        implicit none
        complex(kind) :: mat(Nc,Nc), CLAPW(3*Nc)
        integer :: i,n,m,j, ipiv(Nc), LCLAPW
        
        LCLAPW=3Nc
        
c         ZGETRF  -  compute an LU factorization of a general M-by-N
c         matrix A using partial pivoting with row interchanges.
c         ZGETRI - compute the inverse of a matrix using the LU
c         factorization computed by ZGETRF
          call zgetrf(Nc, Nc, mat, Nc, ipiv, info )
          if(info.ne.0) then
            write(6,*) 'zgetrf error: info = ',info
            stop
          endif
          call zgetri(Nc,mat,Nc,ipiv,CLAPW,LCLAPW,info)
          if(info.ne.0) then
            write(6,*) 'zgetri error: info = ',info
            stop
          endif
        
        return
        end subroutine Invert_Nc
