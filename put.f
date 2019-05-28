c**************************************************************************** 
        subroutine output(istate)
c**************************************************************************** 

c	This block is designed to process the binned data
c	contained in gf, etc. (where the binning was
c	done sequentially).
 
cc***************************************************************************************
	use Global
	use omp_lib
	implicit none
c***************************************************************************************
       integer,parameter:: nt_broyden = (2*100+1)*20
       integer, parameter :: histo_nbins = 100
       integer i,jc,jck,istate,ic,n,rmeas,m,r1_pole,k_histo,histogram(Nc,histo_nbins+1)
       real*8 :: sigg,sigg1,sigg2,sigg3,r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,histo_index(histo_nbins),
     &  sum1(Nc),d_step,min_w,m_r1(Nc,-nwn:nwn),
     &  rhoaveloc_f(-nwn:nwn),rholoc_typ(-nwn:nwn),drholoc_typ(-nwn:nwn),
     &  chir(Nc,-nwn:nwn),data_typ(Nc,-nwn:nwn),locl1(Nc,-nwn:nwn),locl2(Nc,-nwn:nwn)
       real*8 :: xmin,r11(-nwn:nwn)
       real*8 ::  rloc(-nwn:nwn),Gamma_plus_G_Real(1:Nc,-nwn:nwn),Gamma_plus_G_Imag(1:Nc,-nwn:nwn)
	complex(kind):: locl(Nc,-nwn:nwn),csum(Nc,-nwn:nwn)
	complex*16:: GcKfsc_noeta(1:Nc,-nwn:nwn),z1,z2
        real*8:: z1r,z1i,z2r,z2i

c***************************************************************************************
	time1= secnds(0.0)
	if(istate.eq.0) then
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
	  write(42,*) 'concentration, ca'
	  write(42,*) ca
          write(42,*) ' '
	  return
	end if
c	call profstart('put')
c	time1= secnds(0.0)
cc	wtime= OMP_get_wtime()				      
c***************************************************************************************
c  1.--- Impose point group symmetry on ADOS and LogDOS
c***************************************************************************************
!	call Symm_PG_Real(rhof)
!	call Symm_PG_Real(rho_logf)
	
c***************************************************************************************
c  2. ---  normalize quantities e.g., <ln(rho(K))>		
c***************************************************************************************

	
c  --- normazile DOSs
	
	 rhof=rhof/(run)	            ! rho_aver=<-1/pi(G(K,w))>
	 rho_logf = rho_logf/(run)	    ! rho_typ=<Ln(-1/pi(G(K,w)))>
	 logrholoc_f=logrholoc_f/(run)  ! rho_loc=<Ln(-1/pi(G(ii,w)))>
	 rhof_rs_loc= rhof_rs_loc/float(run)  ! loc _rii-separate piece
	 rholoc_av=rholoc_av/float(run)

	 
c ---  for error bars	 
         rhofs=rhofs/float(run)
         rho_logfs=rho_logfs/real(run)
         sigmaf=sigmaf/real(run)
	 gammaf=gammaf/real(run);gammafs=gammafs/real(run)
         rholoc_f=rholoc_f/real(run);rholoc_fs=rholoc_fs/real(run)
         
         
c***************************************************************************************
c  3. ---  Calculate the typical and average density of states (TDOS)
c***************************************************************************************
c       calculate the local average dos for this run
        rhoaveloc_f=zeror
	rho_typ=zeror
        do i=1,Nc
          rhoaveloc_f(:) = rhoaveloc_f(:) + rhof(i,:)
        end do
        rhoaveloc_f=rhoaveloc_f/real(Nc)



c  --- calculate local _typ/avg dos ! The TDOS(R=0) is calculated for all methods	  
	     if(tmt.or.tmt_new) then
	       do n=-nwn,nwn
               rholoc_typ(n)=exp(logrholoc_f(n))
	       end do
	     else
	       do n=-nwn,nwn
               rholoc_typ(n)=logrholoc_f(n)
c              rholoc_typ(n)=exp(logrholoc_f(n)) !for  tmt-dca
	       end do
	     end if  


c  --- calculate typ/avg dos	  	   
c  --- calcualte rho_typ 
         rho_typ=zeror	    
       if(tmt)  then
 	  write(6,*) '***********************************************************'
 	  write(6,*) 'Code using the original CTMT: rho_typ(K) = exp(<Ln(rho(K)>)'		      		      
 	  write(42,*) '***********************************************************'
 	  write(42,*) 'Code using the original CTMT: rho_typ(K) = exp(<Ln(rho(K)>)'		      		      		      	
	  do n=-nwn,nwn        
	    do ic=1,Nc
	      rho_typ(ic,n)=exp(rho_logf(ic,n)) !orig-KTMT procedure	 
	      end do
	     end do	     
c ************************************************     
       else
         if(tmt_new) then
              if (tmt1) then       !ver1 with log  
 	  write(6,*) '***********************************************************'
 	  write(6,*) 'Code using the CTMT: rho_typ(K) = exp(<Ln(rho_ii)>)exp(<Ln(rho(K)/rho_ii)>)'		      		      
 	  write(42,*) '***********************************************************'
 	  write(42,*) 'Code using the CTMT: rho_typ(K) = exp(<Ln(rho_ii)>)exp(<Ln(rho(K)/rho_ii)>)'	         
!$OMP PARALLEL DO
	     do n=-nwn,nwn	    
	      do ic=1,Nc
               rho_typ(ic,n)=exp(logrholoc_f(n))*exp(rho_logf(ic,n)) !ver1 exp^<ln_rhoii>exp(<Ln(rho(K)/rho_ii)>)
	       end do
	      end do
!$OMP END PARALLEL DO
	     end if !if ver1 with log
	     
	   if (tmt2) then       !ver2 -good for multiorbitals
 	  write(6,*) '***********************************************************'
 	  write(6,*) 'Code using the CTMT: rho_typ(K) = exp(<Ln(rho_ii)>)<(rho(K)/rho_ii)>'		      		      
 	  write(42,*) '***********************************************************'
 	  write(42,*) 'Code using the CTMT: rho_typ(K) = exp(<Ln(rho_ii)>)<(rho(K)/rho_ii)>'
!$OMP PARALLEL DO
	    do n=-nwn,nwn	    	    
	     do ic=1,Nc
              rho_typ(ic,n)=exp(logrholoc_f(n))*rho_logf(ic,n) !ver2 - exp^<ln_rhoii><rho(K)/rho_ii> 
	      end do
	     end do
!$OMP END PARALLEL DO	     	
	   end if  !if ver2-w/out orbitals 

	   if(tmt3) then
 	  write(6,*) '***********************************************************'
 	  write(6,*) 'Code using the CTMT: rho_typ(K) = (1-a)*exp(<Ln(rho_ii) + a*<rho(K)'		      		      
 	  write(42,*) '***********************************************************'
 	  write(42,*) 'Code using the CTMT: rho_typ(K) = (1-a)*exp(<Ln(rho_ii) + a*<rho(K)'

c	Calculate the ratio exp(Ln(rho_ii))/rho_ii). Needed for tmt3
	rloc=zeror
	do n=-nwn,nwn
	rloc(n)=rholoc_typ(n)/rholoc_av(n)
	end do

	      do n=-nwn,nwn
	       do ic=1,Nc
	         rho_typ(ic,n)=(1.d0-rloc(n))*exp(logrholoc_f(n))+rloc(n)*rhof(ic,n) ! ver3 (1-a)*exp(<Ln(rho_ii) + a*<rho(K); a=rloc(n)
	        end do
	      end do
	  end if ! ver3.
	
	   if(tmt4) then
 	  write(6,*) '***********************************************************'
 	  write(6,*) 'Code using the CTMT: rho_typ(K) = exp(<Ln(rho_ii)>)'		      		      
 	  write(42,*) '***********************************************************'
 	  write(42,*) 'Code using the CTMT: rho_typ(K) = exp(<Ln(rho_ii)>)'

           do n=-nwn,nwn 
	   rloc(n)=zeror
          do ic=1,Nc
            rloc(n)=rloc(n)+exp(logrholoc_f(n))
           end do
           rloc(n)=rloc(n)/Nc
           end do 
 
           do n=-nwn,nwn 
             do ic=1,Nc
            rho_typ(ic,n)=rloc(n)
           end do
          end do 
	  end if ! ver4.
 

	 else
c   do alg averaging
 	  write(6,*) '***********************************************************'
 	  write(6,*) 'Code using the DCA: rho_typ(K) = rho(K)'		      		      
 	  write(42,*) '***********************************************************'
 	  write(42,*) 'Code using the DCA: rho_typ(K) = rho(K)'
	       do n=-nwn,nwn	  	  
	        do ic=1,Nc
	         rho_typ(ic,n)=rhof(ic,n)
	         end do
	        end do	    
	 end if !tmt new  
	    	    	     	     	     
       end if   ! if tmt  
	
	 
c --- impose p-h symmetry on rho_typ
c --- impose p-h symmetry on rho_typ 
!	call Symm_PH_Real(rho_typ)
 
	
	 
c	Impose point group on TDOS
!	call Symm_PG_Real(rho_typ)
      
  
	    
c***************************************************************************************
c  4.---- Do Hilbert transformation to calculated the Typical Green function (TGF).
c       First the real part.  Note, this is a principal value integral
c***************************************************************************************
	call HilbertT(G_typ,G_typ_real,rho_typ) ! Both these HT give the same result

c***************************************************************************************
c	calcuate local G_typ=(1/Nc)sum_K(G_typ)
c***************************************************************************************
	gf=G_typ	
	chi=zeroc
	if(gsflag==0) then
	do n=-nwn,nwn
	  do ic=1,Nc
            if(ed/=0.0D0) then
	   chi(1,n) = gf(ic,n)
	else
	   chi(1,n) = chi(1,n) + halfr*(gf(ic,n)-conjg(gf(ic,-n)))
	end if
	  end do
	 chi(1,n) = chi(1,n)/(float(Nc))
	end do
	  do n=-nwn,nwn
	    gf(1,n)=chi(1,n)
	  end do
	else
	do n=-nwn,nwn
	  do ic=1,Nc
	   chi(1,n) = chi(1,n) + gf(ic,n)
	  end do
	 chi(1,n) = chi(1,n)/(float(Nc))	
	  end do
	  do n=-nwn,nwn
	    gf(1,n)=chi(1,n)
	  end do
	end if
	


c***************************************************************************************
c 5.---	Coarse-grain with Gamma using the partial density of states calculated using look up
c	table. This is an integration over energy only and as such cheaper than the k_tilde
c	summation.I.e., Calculate GcKf with energy integration via partial 
c	density of states
c***************************************************************************************

!	For the real and imaginary part of the Integral Int_{de/(Z-e_i)}
	r1=0.0D0
	r2=0.0D0
        r3=0.d0
	GcKf = zeroc
	do n=-nwn,nwn
	  do ic=1,Nc
            !if(dimag(G_typ(ic,n)).gt.0.d0) 
     .      !      write(6,*) 'Acausal G_typ',dimag(G_typ(ic,n))
	    z1 = 1.0D0/(G_typ(ic,n))+Gamma_old(ic,n)-ed+Epsbar(ic)
	    z1r= dreal(z1) 
	    z1i= dimag(z1) 
            !if(z1i.le.0.d0) then
            !   write(6,*) 'Acausal z1i',z1i
            !   z1i=dabs(z1i)
            !end if
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
	  enddo
	enddo
!        write(6,*) 'Norm=',r3


c***************************************************************************************
c 6. ----  update Gamma 
c***************** converge the Gamma with broyden method *******
c***************************************************************************************
	 Gamma_new(:,:)=zeroc		

	if(is_broyden) then ! Gamma starts					        	

	Gamma_new = Gamma_old + (1.d0-dfac)*(onec/G_typ - onec/GcKf)			

cc	damp the self-energy with broyden method
	  rms = zeror
	do ic=1,Nc
	     do n=-5,5
	      vector((n+5)+ic*(2*(5+1)),1) = real(Gamma_old(ic,n))
	      vector(nt_broyden+(n+5)+ic*(2*(5+1)),1) = aimag(Gamma_old(ic,n))
	      vector((n+5)+ic*(2*(5+1)),2) = real(Gamma_new(ic,n))
	      vector(nt_broyden+(n+5)+ic*(2*(5+1)),2) = aimag(Gamma_new(ic,n))
		rms=rms+abs(Gamma_new(ic,n)-Gamma_old(ic,n))**2
	    end do
	end do
	  vlen = nt_broyden
	  rms = sqrt(rms/dble(vlen))
	  alpha = 0.4d0
	  broylen = 4
      
	  call broyden(vector, vlen, alpha, rms, iter, broylen,
     &		ub,vt_broyden,fb,df,vold,a_broyden,b_broyden,
     &		d_broyden,cm,wb,ipiv,mbroylen,mvlen)


	do ic=1,Nc
	  do n=-5,5
	      Gamma_new(ic,n) = cmplx(vector(n+5+ic*(2*(5+1)),2),
     &		vector(nt_broyden+n+5+ic*(2*(5+1)),2))
	    end do
	end do 

	else ! Not using broyden starts
c***************************************************************************************
cc             converge the Gamma with linear mixing 
c***************************************************************************************    
         do ic=1,Nc
            do n=-nwn,nwn
	      Gamma_new(ic,n) =  Gamma_old(ic,n) +(1.d0-dfac)*(onec/G_typ(ic,n) - onec/GcKf(ic,n)) 
            end do
         end do

	end if ! Gamma ends

c***************************************************************************************				      				
ctest for convergence by inspecting the difference between New and Old Hybridization rate
c***************************************************************************************				      	      
        testl=0.d0      
        test=0.d0						      
        do ic=1,Nc
	  do n=-nwn,nwn
	    testl = testl + abs(aimag(Gamma_old(ic,n)))
            test = test + abs(aimag(Gamma_new(ic,n))) 
	   end do
        end do 
	
        if(iter.ge.2) then					      
          r1=abs(testl-test)/(test*float(Nc))
          write(6,*) ' test= ',r1 
 	  write(6,*) '**********************************************************'		      
          write(42,*) ' test= ',r1 
 	  write(42,*) '**********************************************************'		      
        end if 
	Gamma_old(:,:) = Gamma_new(:,:)
	
c***************************************************************************************
c. ----   update Gscript 
c ***************************************************************************************
           GcKfsc(:,:)=zeroc  ! This is the conventional way to evaluate Gscript. Okay away from criticality
	  if(iter.eq.niter-1) then
          do ic=1,Nc
           do n=-nwn,nwn          
	   GcKfsc(ic,n) = 1.0D0/(wn(n)-Gamma_new(ic,n)-Epsbar(ic)+eta)
                GcKfsc_noeta(ic,n) = 1.0D0/(wn(n)-Gamma_new(ic,n)-Epsbar(ic))
	   end do  
	   end do
	else
          do ic=1,Nc
           do n=-nwn,nwn          
	   GcKfsc(ic,n) = 1.0D0/(wn(n)-Gamma_new(ic,n)-Epsbar(ic))
	   end do  
	   end do
	endif

c else using poles
c***************************************************************************************
	if(gscript_pole) then	! This method uses delta-function to calculate Gscript close to criticality
	
	do ic =1,Nc
	 sum1(ic) = zeror
	 do n=-nwn,nwn
         sum1(ic)=sum1(ic)+dwn(n)*(-1.d0/pi)*aimag(GcKfsc(ic,n))
	 end do
c	 write(*,*) 'Norm_GcKfsc-before taking pole=',sum1(ic), 'ic=',ic
	end do
	
	if(iter.eq.niter-1) then ! Note eta is added at the last iteration only for Average quantities
	  do ic=1,Nc
	    do n=-nwn,nwn
              m_r1(ic,n)=real(wn(n)+eta-Epsbar(ic)-Gamma_new(ic,n))
c	      write(667,*) wn(n), abs(m_r1(ic,n))
    	    end do
c	    write(667,*) '    '
    	  end do
	else
	  do ic=1,Nc
	    do n=-nwn,nwn
              m_r1(ic,n)=real(wn(n)-Epsbar(ic)-Gamma_new(ic,n))
c	      write(667,*) wn(n), abs(m_r1(ic,n))
    	    end do
c	    write(667,*) '    '
    	  end do
	end if

    	  do ic=1,Nc    ! ic loop starts
    	    min_w=abs(m_r1(ic,-nwn))
	    x_pole = -nwn
    	    do n=-nwn,nwn
    	      if(min_w.gt.abs(m_r1(ic,n))) then ! finding minimum starts
    	        min_w=abs(m_r1(ic,n))
    	        x_pole=n
    	      end if	! finding minimum end        
            end do ! n
 
    	if((-1.d0/pi)*aimag(Gamma_new(ic,x_pole)) .ge. afact*dwn(x_pole)) then 
c        if( abs(1.-sum1(ic)) .lt. 0.2) then !checking norm 
	  if(iter.eq.niter-1) then
	    do n=-nwn,nwn
	    GcKfsc(ic,n) = onec/(wn(n)-Gamma_new(ic,n)-Epsbar(ic)+eta)
	    end do 
	else
	    do n=-nwn,nwn
	    GcKfsc(ic,n) = onec/(wn(n)-Gamma_new(ic,n)-Epsbar(ic))
	    end do   
	end if  
	else 
	    write(6,*) 'Pole of cell', ic, 'occurs at', x_pole, wn(x_pole),abs(aimag(Gamma_new(ic,x_pole)))
	    write(42,*) 'Pole of cell', ic, 'occurs at', x_pole, wn(x_pole),abs(aimag(Gamma_new(ic,x_pole)))
	    	    
	    do n=-nwn,nwn
	    if(n .eq. x_pole) then ! Calculating Gscript when pole occurs. Approximate using delta function
             GcKfsc(ic,n) = cmplx(0.d0,-pi/dwn(x_pole))
                	     
	    else if (n .ne. x_pole) then		
c              GcKfsc(ic,n) = cmplx(1.d0/(wn(n)-real(wn(x_pole))),0.d0)                  
             GcKfsc(ic,n) = cmplx(1.d0/(wn(n)-Epsbar(ic)-real(Gamma_new(ic,n))),0.d0) !  This is p.v.: at w=wpole, skip, else keep the 1/w-(epsk-Gamma) the same as in Hilbert transform above              
	    end if	
    	   end do
    	   
    	 end if
    	 end do    ! ic loop ends
 
	end if	! gscript_pole ends			   
	   		   
	do n=-nwn,nwn
          if(isd.eq.1.and.iter.gt.1.or.isd.eq.1.and.iso.eq.1) then    
c         add some damping  to the calculation of gft
            gft(n)=dfac*gf(1,n) + (1.d0-dfac)*gft(n)
          else
            gft(n)=gf(1,n)
          end if                                       
	end do						      

c***************************************************************************************				      
c	Write some quantity for self-consistency check on the screen	
c***************************************************************************************

	r0=zeror;r1=zeror;r4=zeror
	do n=-nwn,nwn
	r4=r4+rhoaveloc_f(n)*dwn(n)
	do ic=1,Nc
	r0=r0+rhof(ic,n)*dwn(n)
	r1=r1+(-1.d0/pi)*aimag(GcKfsc(ic,n))*dwn(n)
	end do
	end do
	r0=r0/float(Nc)
	r1=r1/float(Nc)
c	write(777,*) iter,r1,r4

 	  write(42,*) '**************************************************'
 	  write(42,*) 'Self-Check of the calculations, sum_rules and others'
 	  write(42,*) '**************************************************'
	
	do ic =1,Nc
	sum1(ic) = zeror
	do n=-nwn,nwn
        sum1(ic)=sum1(ic)+dwn(n)*(-1.d0/pi)*aimag(GcKfsc(ic,n))
	end do
	write(*,*) 'Norm_GcKfsc-after=',sum1(ic), 'ic=',ic
	write(42,*) 'Norm_GcKfsc-after=',sum1(ic), 'ic=',ic
	end do
	
c	Calculate gamma at the band center
          r3=zeror
           do ic=1,Nc
           r3=r3+aimag(-Gamma_new(ic,0))
           end do
c           write(*,*) 'Gamma-w0-oc',r3,r3/Nc
c***************************************************************************************
c	Write out some observables vs iteration to file 
c***************************************************************************************
      open(210,file='observ-vs-iter.dat',status='unknown')
	if (iter.eq.1) write(210,"('# iter,TDOS(R=0),TDOS(K,w),Gamma(K,w),ADOS_GcKfsc,ADOS_rhoaveloc_f(n)') ") 
	write(210,"(i3,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6)") iter,rholoc_typ(0),-aimag(gft(0))/pi,r3/Nc,r1,r4

	write(6,*) 'tdos(0)=',-aimag(gft(0))/pi, 'gamma(0)=', r3/Nc
	write(6,*) 'GdKK_sum_rule=',r0
	write(6,*) 'G(i,i,w)_sum_rule=',r4
	write(6,*) 'GcKfsc_sum_rule=',r1
	write(6,*) 'nwsc=',nwsc	
c	write(6,*) 'Check for singularity in Gscript. x_pole=',x_pole			
	write(42,*) 'tdos(0)=',-aimag(gft(0))/pi, 'gamma(0)=', r3/Nc
	write(42,*) 'nwsc=',nwsc		
	write(42,*) 'GdKK_sum_rule=',r0	
	write(42,*) 'G(i,i,w)_sum_rule=',r4	
	write(42,*) 'GcKfsc_sum_rule=',r1
c	write(42,*) 'Check for singularity in Gscript. x_pole=',x_pole
c	If m_r1 .ne. -10000, then singularity occured in Gscript									      
c***************************************************************************************
c       Calculate the error bars for the measurements. For e.g., error bar for 
c	the linearly averaged K-resolved rho(K,w) denoted as drhof, the error 
c	bar for the typicaly averaged K-resolved rho_typ(K,w) denoted as drho_log
c***************************************************************************************				      
cc        write(6,*) ' '  					      
cc        write(6,*) ' ic    sigg      sigg1      sigg2     sigg3'		
cc        write(42,*) ' '  					      
cc        write(42,*) ' ic   sigg      sigg1      sigg2     sigg3'		

        drhof=sqrt((abs(rhofs-rhof**2))/float(run-1))
	drho_log=sqrt((abs(rho_logfs-rho_logf**2))/float(run-1))
	d_gamma=sqrt((abs(gammafs-gammaf**2))/float(run-1))
        drholoc_typ=sqrt((abs(rholoc_fs-rholoc_f**2))/float(run-1))
cccc

 	do n =-nwn,nwn	
	    sigg3=sqrt((abs(rholoc_fs(n)-rholoc_f(n)**2))/float(run-1))				      
	do ic=1,Nc						      
 	  if(run.gt.1) then
            sigg=sqrt((abs(rhofs(ic,n)-rhof(ic,n))**2)/float(run-1))
	    sigg1=(sqrt((abs(rho_logfs(ic,n)-rho_logf(ic,n)**2)))*-aimag(gft(n))/float(run-1))
	    sigg2=sqrt((abs(gammafs(ic,n)-gammaf(ic,n)**2))/float(run-1))
	  end if
cc         if(n.eq.0) write(6,"(i2,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6)")  ic, sigg,      sigg1,      sigg2,     sigg3			   
cc          if(n.eq.0) write(42,"(i2,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6)") ic, sigg,      sigg1,      sigg2,     sigg3
	end do
	end do

c***************************************************************************************
cc	Write outputs at the last iteration
c***************************************************************************************
       if(iter.eq.niter) then 
!#######################################################################
!	Calculate the average G and then self-energy needed for spectra
!#######################################################################
	call HilbertT(GcKf_av,G_typ_real,rhof) ! Calculate the average cluster G only at the last iteration
!	Add a small imaginary part to reduce spikes



! 1. The self-energy
          sigma_loc(1:Nc,-nwn:nwn)=1.0D0/GcKfsc(1:Nc,-nwn:nwn)-1.d0/GcKf_av(1:Nc,-nwn:nwn) 

!	Print Average DOS that is smooth
      open(488,file='ADoS.dat',status='unknown')
	do n=-nwn,nwn
	r1=0.0D0
	do ic=1,Nc
	r1=r1-dimag(GcKf_av(ic,n))/pi
	enddo
	write(488,*) wn(n), r1/Nc
	enddo
	close(488)


!###################################################################
!	Print the self-energy needed for spectra
!####################################################################
!1. Print the self-energy to file
          write(40,*) ' n   ic     w        sigma_loc(ic,w)  '  		      
          do ic=1,Nc						      
          do n=-nwn,nwn						      
            write(40,"(i5,1x,i2,1x,f11.7,1x,e15.8,1x,e15.8)") n,ic,wn(n),sigma_loc(ic,n)-eta			      
          end do						      
          end do 


          write(40,*) ' n   ic     w        GcKf_av(ic,w)  '  		      
          do ic=1,Nc						      
          do n=-nwn,nwn						      
            write(40,"(i5,1x,i2,1x,f11.7,1x,e15.8,1x,e15.8)") n,ic,wn(n),GcKf_av(ic,n)			      
          end do						      
          end do
	end if !Iter==niter

	if(iter.eq.niter-1) then !We have to print out typical observables before adding eta in the calculation of Average quantities above
c  ----- Calculate the typical self-energy 
          sigma(:,:)=1.0D0/GcKfsc_noeta(:,:)-1.d0/GcKf(:,:) ! To make sure there is no eta in sigma typical

!***************************************************************************************
!	Calculate the scattering rate 1/tau = -2*Imag(SelfEnergy(\omega=0)) 
!	= -2*Imag(Go(\omega=0)*V_bar^2). Note, this is directly related transport
!***************************************************************************************
	open(799,file='tau.dat',status='unknown')
!	do n=-nwn,nwn
	r1=0.0D0
	do ic=1,Nc
	r1=r1-2.0D0*dimag(sigma(ic,0))
	enddo
	write(799,*) V0,r1/dfloat(Nc)
	close(799)

c***************************************************************************************				      
c	The localization probablity is only calcualted on the last
c	iteration.
c***************************************************************************************
        if(ichit.eq.1) then
	do n=1,nt
	  chit(n)=chit(n)/float(meas*Nc*run)
	  end do
	end if 

        if(ichim.eq.1) then
	do n=1,neta
	  chimeas(n)=chimeas(n)/float(meas*Nc*run)
	  end do
	end if	 
	if(histo_enable) then	! This calculates the histogram.


c ---step1: finding min and max for histo_meas at each ic-------------------:
        do ic=1,Ncw	  
	  	histo_min(ic) = histo_meas(ic,1) ! initialization of min
	        histo_max(ic) = histo_meas(ic,1) ! initialization of max
	        
    	     do i=1,meas      
    	      if(histo_meas(ic,i).gt.histo_max(ic)) then
    	         histo_max(ic)= histo_meas(ic,i) 
    	      end if
    	      if(histo_meas(ic,i).lt.histo_min(ic)) then
    	        histo_min(ic)= histo_meas(ic,i) 
    	      end if
             end do !end meas
	end do  !end ic
c	write(68,"(200(f12.6))") (histo_min(ic),ic=1,Nc)	! Remove after debugging
c	write(69,"(200(f12.6))") (histo_max(ic),ic=1,Nc)



c ----step2: binning the data------------
          histogram(:,:)=0	
	   do ic=1,Ncw
	    do i=1,meas
	     k_histo=int(1.+(histo_meas(ic,i)-histo_min(ic))*histo_nbins/(histo_max(ic)-histo_min(ic)) )
	     histogram(ic,k_histo)=histogram(ic,k_histo)+1
	 end do
	end do

      open(121,file='histogram.dat',status='unknown')
c--- step3: -----writing PDF data to file ------
	do ic=1,Ncw
	 do i=1,histo_nbins
         write(121,*) histo_min(ic)+(i-1)*(histo_max(ic)-histo_min(ic))/histo_nbins,histogram(ic,i)
	 end do
	 write(121,*) ' ' 
	end do
	close(121)


	end if ! End the calculation of the histogram


c  -- writing Gamma to file to see how it develops for diff K using wedge symmetry	
	
      open(999,file='gamma_K.dat',status='unknown')
            do ic=1,Ncw
              do n=-nwn,nwn     
	      write(999,*) wn(n), -aimag(Gamma_new(ic,n)),ic
             end do
              write(999,*) '  '   
         end do
	close(999)   

c  -- writing Gsript to file to see how it develops for diff frequency	
      open(124,file='gscript_w.dat',status='unknown')
           do ic=1,Ncw
       	   do n=-nwn,nwn
	   write(124,*) wn(n), -aimag(GcKfsc(ic,n))
           end do
           write(124,*) '   '
         end do
	close(124)                  

c  -- writing Rho_typical to file to see how it develops for diff frequency	
      open(125,file='tdos_w.dat',status='unknown')
             do ic=1,Ncw
              do n=-nwn,nwn
               write(125,*) wn(n),real(rho_typ(ic,n))             
              end do
               write(125,*) '  '
             end do
	close(125)  
       if(isg.eq.1) then ! Chinedu removed iter==niter since it is inside the loop already
c	 Symmetrize and write the two-particle chicharge_meas to a file
	 open(49,file='chicharge.dat',status='unknown')
c	 normalize the data
         chicharge_meas=chicharge_meas/dfloat(run*meas)
c	 symmetrize the data
	  G_chg=zeroc
	 do n=-nwn,nwn !The last number is zero. This removes it
	   do ick=1,Nc
	   do jck=1,Nc
           G_chg(ick,jck,n) = zeroc
             do igroup=1,ngroup
	       G_chg(ick,jck,n)=G_chg(ick,jck,n) +
     &  		(1.d0/dfloat(ngroup))*chicharge_meas(ickequ(ick,igroup),ickequ(jck,igroup),n)
	     end do
             chicharge_meas(ick,jck,n)=G_chg(ick,jck,n)
	     write(49,"(1x,i2,1x,i2,2x,i5,1x,2(1x,e15.8))") ick,jck,n,chicharge_meas(ick,jck,n)
	   end do
	   end do
	 end do
	 close(49)
       end if	 
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Writing observables at the band center for easy access
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c......1. Gamma for each K. Mark said we need this to look at Sector selected Anderson Transition
		
        open(unit=38,file='observ_w0.dat',status='unknown')
c	Calculate important observables at the band center
        r1=zeror;r2=zeror;r3=zeror;r4=zeror
        do ic=1,Nc
          r1=r1+aimag(-Gamma_new(ic,0))
          r2=r2+drhof(ic,0)
	  r3=r3+drho_log(ic,0)
	  r4=r4+d_gamma(ic,0)
        end do

	write(38,"('#     This is the values of Gamma for each K at wn(0)')") 
        write(38,"(90(f10.6))") (aimag(-Gamma_new(ic,0)),ic=1,Nc)

c......2.  Other observables just as it is in the dos.dat output
	write(38,"('#     TDOS(K,w)    Gamma(0,w)           ADOS(K,w)       ADOS_Error '
     &		      '       TDOS_Error    Gamma_Error    TDOS_Loc(R=0)    TDOS_Loc_Error')")
        write(38,"(8(f15.6))") -aimag(gft(0))/pi,r1/Nc,rhoaveloc_f(0),r2/Nc,-aimag(gft(0))*r3/Nc,r4/Nc,rholoc_typ(0),
     &            drholoc_typ(0)

	close(38)
          rewind(40)						      
          write(40,*) 'Results for nwn,Nc,nover,cluster type='	      
          write(40,*) nwn,Nc,nover,cluster 		      
          write(40,*) 'ed,V,tprime,iter,run,meas'  	      
          write(40,"(1x,f9.6,1x,f9.6,1x,f9.6,1x,i4,1x,i4,1x,i4) ") ed,V0,tprime,niter,run,meas


          if(ichit.eq.1) then
 	   write(43,"('#    t     Chit(t)')") 
	    do n=1,nt
	      write(43,*) t(n),chit(n)
	    end do
	  end if

          if(ichim.eq.1) then						      
 	    write(44,"('#   eta   Chit(eta)')") 
	    do n=1,neta
	      write(44,*) weta(n),chimeas(n)
	    end do
	  end if						      
          write(40,*) '   w    Re[G(1,w)]     Im[G(1,w)]   Gamma(0,w) '
	  write(41,"('#     w           TDOS(K,w)           Gamma(0,w)           ADOS(K,w)       ADOS_Error '
     &		      '       TDOS_Error    Gamma_Error    TDOS_Loc(R=0)       TDOS_Loc_Error')") 
	  write(47,"('#     w        drhof       drho_log         d_gamma    ')") 
	  write(48,"('#     w    Real_sigma(ic,w)  -Imag_sigma(ic,w)')  ") 

          do n=-nwn,nwn
	    r1=0.d0
	    r2=0.d0
	    r3=0.d0
	    r4=0.d0
	    r5=0.d0
	    r6=0.d0
	    r7=0.d0
	    r8=0.d0
	    do ic=1,Nc
	      r1=r1+aimag(-Gamma_new(ic,n))
c	      r2=r2+rhof(ic,n)
c	      r2=r2+rho_local(ic,n)
	      r3=r3+drhof(ic,n)
	      r5=r5+real(sigma(ic,n))
	      r6=r6-aimag(sigma(ic,n))
	      r7=r7+drho_log(ic,n)
	      r8=r8+d_gamma(ic,n)
	    end do
	    r1=r1/Nc
c	    r2=r2/float(Nc)
	    r2=rhoaveloc_f(n)
	    r3=r3/float(Nc)
c	    r4=r4/float(Nc)
	    r5=r5/float(Nc)
	    r6=r6/float(Nc)
	    r7=r7/float(Nc)
	    r8=r8/float(Nc)
!	    r9=-1.0D0/pi*aimag(GcKf_av(1,n))
!	    r6=-aimag(gft(n))/pi
            write(40,"(1x,f10.6,10(2x,e15.8))") wn(n),real(gft(n)),aimag(gft(n)),r1
            write(41,"(1x,f15.6,20(2x,e15.8))") wn(n),-aimag(gft(n))/pi,r1,r2,r3,-aimag(gft(n))*r7,r8,rholoc_typ(n),
     &            drholoc_typ(n)
            write(48,"(1x,f10.6,2x,e15.8,2x,e15.8)") wn(n),r5,r6	
	  end do
c	Print constant Gamma by integrating out frequency(W) and momentum (K).
        open(unit=39,file='Gamma_const.dat',status='unknown')
	r6=0.d0	
	do n=-nwn,nwn
            do ic=1,Nc
              r6=r6+dwn(n)*aimag(-Gamma_new(ic,n))
            end do
c	r6=r6/Nc
	end do 
	write(39,*) V0,r6/Nc
	close(39) 
						
          write(40,*) '  n    V(n)  the potentials    ' 	      
          do n=1,Nc						      
            write(40,*) n,V(n)  				      
          end do						      
          write(40,*) ' n   ic     w        sigma(ic,w)  '  		      
          do ic=1,Nc						      
          do n=-nwn,nwn						      
            write(40,"(i5,1x,i2,1x,f11.7,1x,e15.8,1x,e15.8)") n,ic,wn(n),sigma(ic,n)			      
          end do						      
          end do


        end if        ! iter.eq.niter				      
c	write(200,*) iter, -aimag(gft(0))/pi
c	write(201,*) iter, -aimag(GcKfsc(1,0)),-aimag(GcKfsc(4,0))
	wtime= OMP_get_wtime()-wtime
 	time2= secnds(time1)
c        write(6,*) 'time in put= ',time2
        put_time=put_time+time2
        write(42,*) 'time in put= ',time2

c        call profend('put')
        return
        end
        
c        include  'quadpack.f90 '
