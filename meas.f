c***************************************************************************************
        subroutine measure
c***************************************************************************************
c       Make a series of measurements.  
c***************************************************************************************
	use Global
	use omp_lib
c***************************************************************************************
 	implicit none
        integer i,n,m,p,ic,jc,jcr,jck,count
	complex*16 :: c1, csum(Nc,-nwn:nwn),Gdkk_av(Nc,-nwn:nwn),bin_histo(Nc,Nc,-nwn:nwn,3000)
	real*8 :: r1,r2,r3,r4,histogram(Nc,-nwn:nwn)
c***************************************************************************************
	time1= secnds(0.0)
        call profstart('meas')
	wtime= OMP_get_wtime()
c***************************************************************************************	
c       In addition to the average green function, we need to calculate the angle (K)
c       resolved typical density of states.  It has the form
c       
c                                                     /   rho(K,w)      \
c       rho_typ(K,w) = exp[ < sum_i ln rho_ii(w) > ] <------------------->    (1)
c                                                     \ sum_i rho_ii(w) /
c
c       which will have to be built from multiple parts below.
c***************************************************************************************
c       Form the cluster averaged real-space Green function.  Each entry in the G(i,j,w) 
c       matrix contributes to the average of G(i-j,w).
c***************************************************************************************
	Gdff = zeroc
        do n=-nwn,nwn
	  do ic=1,Nc
	  do jc=1,Nc	
	    Gdff(icrdiff(ic,jc),n)= Gdff(icrdiff(ic,jc),n)+G(ic,jc,n)/float(Nc)
          end do
	  end do
        end do


c***************************************************************************************
cc      Fourier transform Gdff to K-space (Note, we dont need 1/Nc here) forming 
c       GdKK(ick,n).  Note that  rho(K,w) = -1/pi Im G(K,w).
c***************************************************************************************
	do n=-nwn,nwn		
	  do ick=1,Nc						      
            csum(ick,n)=zeroc
!	    GdKK(ick,n)=zeroc
            do icr=1,Nc				      
              csum(ick,n)=csum(ick,n)+
     &                  FTCoefs_R_to_K(ick,icr)*Gdff(icr,n)	        
	    end do
c for some Nc, it happens that at the boundary at +/- nw_max, causality is broken with 	aimag(csum).ge. 0.d0
c here we avoid this numerical problem by making it -depsi 
            if((dimag(csum(ick,n)) >= 0.d0).OR.(dimag(csum(ick,n)) /= dimag(csum(ick,n)))) then
                csum(ick,n) = dcmplx(dreal(csum(ick,n)),-depsi**2)
            end if
            GdKK(ick,n)=csum(ick,n)
	  end do  
	end do 

	
c*********************************************************************************************************
c 	prepare local terms, first the average local dos to be used in the denominator of the linear 
c       average in Eq. 1, 
c
c       rho_ii(w)=(1/Nc) sum_i G_ii (w)
c*********************************************************************************************************
! 	
	do n=-nwn,nwn
          rho_ii(n)=0.d0
          do ic=1,Nc          
            rho_ii(n)=rho_ii(n)+(-(1./pi)*aimag(G(ic,ic,n)))
          end do
          rho_ii(n)=rho_ii(n)/Nc
        end do    
	
c*********************************************************************************************************
cc	Calculate the Algebraic K-resolved density of states which is the numerator of the linear 
c       average in Eq. 1.
c*********************************************************************************************************
        rho(:,:)=(-1.0/pi)*aimag(GdKK(:,:)) 

c   -   Accumulate the agebraically averaged rho(K,w)			      
	rho_meas(:,:)=rho_meas(:,:)+rho(:,:)
!	Prepare the part with the same boundary as in the DCA for spectra
!	rho_av_meas(:,:) = rho_av_meas(:,:)-1.0D0/pi*aimag(GdKK_av(:,:))

  
c*********************************************************************************************************
c 	Now construct rho_typ(K,w) for the different TMT schemes	  
c*********************************************************************************************************
	  
	if (tmt_new) then ! if true use one of new TMT schemes
           if (tmt1) then       !ver1 with log   
	     do n=-nwn,nwn    
             do ick=1,Nc						      
	      logrho_meas(ick,n) = logrho_meas(ick,n) + 
     &          log(abs(rho(ick,n)/rho_ii(n))) !ver 1 exp(Ln(rho_ii))*exp(<Ln(rho(K)/rho_ii)>)	    
	      end do  						      
	      end do 
           end if !if ver1
         
          if (tmt2) then    !ver2 with no log-good for multiorbitals 
	    do n=-nwn,nwn        
            do ick=1,Nc						      
	      logrho_meas(ick,n) = logrho_meas(ick,n) +
     &          rho(ick,n)/rho_ii(n) !ver 2 exp(Ln(rho_ii))*<rho(K)/rho_ii>	    
	    end do  						      
	    end do 
          end if !if ver2

          if (tmt3) then    !ver3 with no log-good for multiorbitals 
            do n=-nwn,nwn
            do i=1,Nc
              rholoc_meas(n) = rholoc_meas(n)-aimag(G(i,i,n))/pi 
            end do                                      
            end do 
          end if !if ver3
c***************************************************************************************				      
c --- local typ dos
c***************************************************************************************				      
	    do n=-nwn,nwn        
            do i=1,Nc
              logrholoc_meas(n) = logrholoc_meas(n) +
     &          log(abs(-aimag(G(i,i,n))/pi))                                      
            end do  
            end do          
        else 
	  	 	  	
c use tmt -original scheme: CTMT
          if (tmt) then
	    do n=-nwn,nwn    
	    do ick=1,Nc						      
c	      logrho_meas(ick,n) = logrho_meas(ick,n) + log(abs(rho(ick,n)))	
	      logrho_meas(ick,n) = logrho_meas(ick,n) + 
     &                          log(rho(ick,n))	        
	    end do  						      
	    end do 
       
	   
c --- local typ dos
	    do n=-nwn,nwn        
            do i=1,Nc
              logrholoc_meas(n) = logrholoc_meas(n) +
     &          log(abs(-aimag(G(i,i,n))/pi))                                        
            end do 
            end do 
     
                  
          else
cc          Calculate Algebraic local-density of states
c ---       local dos for algebaic aver case  =<sum_ii G(i,i,w)>alg    W/out taking log 
            do n=-nwn,nwn
            do i=1,Nc
              logrholoc_meas(n) = logrholoc_meas(n)-aimag(G(i,i,n))/pi 
c     logrholoc_meas(n) = logrholoc_meas(n)+log(abs(-aimag(G(i,i,n)))/pi )
            end do                                      
            end do      
          end if !if tmt     
                  
        end if   !new tmt procedure
                   

  	if(iter==niter) then
        open(unit=111,file='binned_histogram.dat',status='unknown')
	write(111,'(10f15.6)') (-dimag(G(jc,jc,0))/(pi*Nc),jc=1,Ncw)
!	close(111)
  	endif

c***************************************************************************************				      
c	Prepare the data for the Histogram. 
c***************************************************************************************

	if(histo_enable) then ! Bin the data needed for calculating histogram in put.f
	  if(iter.eq.niter-1) then
	    do ic=1,Nc
	       histo_meas(ic,nmeas)=dabs(-dimag(Gdkk(ic,Epsbar_histo(ic)))/(pi*Nc)) ! K-rho Rho(K,w)
	     end do
	  end if
	end if

c***************************************************************************************				      
c	  Measure the probabilty that an electron remains at x
c	  for all time (Thouless, Physics Reports, 1974, p121).
c         The wavefunction evolves according to the Green function:
c
c                       / /  3
c         Psi(x',t') =  | | d x dt G(x',x,t'-t) Psi(x,t)
c         	        / /
c
c         Thus, if at time t=0 there is an electron at x=0 so that
c         Psi(x,t) = delta(x) delta(t), then the probability that the
c	  electron remains at x=0 (the return probability) is given by
c
c                         2        b  / oo               2
c          lim  |G(x,x,t)| = lim  --- |  de |G(x,x,e+ib)|    b is real ->0
c         t->oo              b->0  pi /-oo
c
c                                / oo        A(w) A(w')
c                        =   2bi |   dw dw' -----------
c                                /-oo        w-w'-2ib
c
c                                / oo         A(w) A(w')
c                        =  4b^2 |   dw dw' ---------------
c                           ---   /-oo       (w-w')^2+ 4b^2
c			    pi^2
c***************************************************************************************
	if(iter.eq.niter-1.and.ichim.eq.1) then
	  r4=1.d0/pi**2
	  do p=1,neta
	    r3=(2.d0*weta(p))**2
	    do n=-nwn,nwn
	    do m=-nwn,nwn
	      r2= r3/((wn(n)-wn(m))**2 + r3)
	      r1=zeror
	      do ic=1,Nc
	  	r1=r1+dimag(G(ic,ic,n))*dimag(G(ic,ic,m))
	      end do
	      chimeas(p)= chimeas(p) + r1*r2*r4*dwn(n)*dwn(m)
	    end do
	    end do
	  end do
	end if
c***************************************************************************************				      
c         If the integrand of the first expression has a finite limit
c         in any energy range, the electrons at this energy are 
c         localized.  Thus, a more informative quantity might be
c
c               b   | / oo     A(w)   | 2
c         lim  ---  | |   dw -------- |
c         b->0  pi  | /-oo    w-e-ib  |
c
c         If we set e=0, then we measure the probablity that electrons
c	  at the fermi energy are localized.
c
c	  Yet another useful thing to measure is the probability that 
c	  the electron remains at the site as a function of time:
c
c                           2                                          2
c                 |G(x,x,t)|   	   1       - | / oo    - iwt          |
c          P(t) = ----------  = -------   \  | |   dw e      G(x,x,w) |
c                    pi*pi      pi^2 Nc   /  | /-oo                   |
c                                          -
c                                          x
c
c	  which is what we will measure (order the next 2 loops to optimize).
c	  Convergence is improved by adding and subtracting an exact result
c	             1
c         G0(w) = -------        G0(t) = -i exp(-dt)
c                  w + id
c***************************************************************************************
	if(iter.eq.niter-1.and.ichit.eq.1) then
	  r1=1./(2.*pi)
	  do p=1,nt
	  do ic=1,Nc
	    c1=zeroc
	    do n=-nwn,nwn
	      c1=c1+dwn(n)*(G(ic,ic,n)-G0w(n))*four(n,p)
	    end do
	    chit(p) = chit(p) + (abs(c1*r1 + G0t(p)))**2
	  end do
	  end do

	end if

	if(isg.eq.1.and.iter.eq.niter-1) then

	
c
c  	  Now measure the particle-hole response needed to calculate
c  	  the conductivity.  This is simplified by the fact that the
c  	  disorder potential is non-dynamical (i.e. since it is
c  	  quadratic, no HHS field needs to be introduced).  Thus,
c  	  frequency is conserved in the scattering processes,
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
c  	  In this subroutine, we must measure Chi(K,K',w,w'), in the
c  	  subroutine put, we will use the impurity potential averaged
c  	  Chi(K,K',w,w') and the equation above, to extract the Nc by Nc
c  	  matrix Gamma for each w and w'.
c
c  	  In calculating Chi, it is important to remember that causal
c  	  green functions must be used; whereas, in this code, the green
c  	  functions are all retarded!  The causal averaged green function
c  	  has an imaginary part that is not positive definite, in fact
c
c  	  Im[Gc(K,w)] >0  for  w<ed
c
c  	  Im[Gc(K,w)] <0  for  w>ed
c
c  	  The relationship between the causal green function Gc and the
c  	  retarded one used in this code is
c
c  	  Gc(K,K',w) =  G(K,K',w)  for w>ed
c
c  	  Gc(K,K',w) =  G*(K,K',w) for w<ed
c
c  	  Chi, defined above, has the form
c
c  	  Chi(K,K',w,w') = < Gc(K,K',w) Gc(K',K,w+w') >_V
c
c  	  where the brackets <>_V indicate an average over the disorder
c  	  configurations.
c	
c         To calculate the DC limit (w'-->0) we will set w'=wn(1)

c***************************************************************************************
c         First, calculate G(K,K',w) using a two step process which is more efficient: 
c         G(X,X',w) --> G(K,X',w) --> G(K,K',w).  We will use G_chg, Gsc_temp and Ginv 
c         as temporary arrays.

c         G(X,X',w) --> G(K,X',w) 
	  Gsc_temp=zeroc ! I changed Gcinv to G_chg. It is over-writing the Gcinv coming from geod.f needed in update.
          G_chg=zeroc
          do ic=1,Nc
          do jc=1,Nc
          do ick=1,Nc
          do n=-nwn,nwn
            G_chg(ick,jc,n) = G_chg(ick,jc,n) + FTCoefs_R_to_K(ick,ic)*G(ic,jc,n)
          end do
          end do
          end do
          end do
c         G(K,X',w) --> G(K,K',w)
          Ginv=zeroc
          do ick=1,Nc
          do jc=1,Nc
          do jck=1,Nc
          do n=-nwn,nwn
            Ginv(ick,jck,n) = Ginv(ick,jck,n) + G_chg(ick,jc,n)*dconjg(FTCoefs_R_to_K(jck,jc))
          end do
          end do
          end do
          end do

c         Now calculate 
c  	  Gc(K,K',w) =  G(K,K',w)  for w<ed
c
c  	  Gc(K,K',w) =  G*(K,K',w) for w>ed
          do n=-nwn,ned
            Gsc_temp(:,:,n) = Ginv(:,:,n)
          end do
          do n=ned+1,nwn
            Gsc_temp(:,:,n) = dconjg(Ginv(:,:,n))
          end do

c         Now calculate
c
c  	  Chi(K,K',w,w') = < Gc(K,K',w) Gc(K',K,w+w') >_V
c
c             / K      w     K'\     ick   <-> K'
c         or /  -------<------- \    jck   <-> K
c            \  ------->------- /    wn(1) <-> w'
c             \ K     w+w'   K'/
c
          do ick=1,Nc
          do jck=1,Nc
          do n=-nwn,nwn-1
	    chicharge_meas(ick,jck,n) = chicharge_meas(ick,jck,n) 
     &	    + Gsc_temp(jck,ick,n)*Gsc_temp(ick,jck,n+1)
	  end do
	  end do
	  end do
         
	end if
        call profend('meas')
 	time2= secnds(time1)
        meas_time=meas_time+time2
        return
        end
	
