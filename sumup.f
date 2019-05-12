c***************************************************************************************
	subroutine sumup
c***************************************************************************************
c
c***************************************************************************************
	use Global
	use omp_lib
	implicit none
c***************************************************************************************
	integer i,ic,n,m
        real(kind) :: x_ratio(-nwn:nwn),rhoaveloc_run(-nwn:nwn),
     &          d_step,m_r1(Nc,-nwn:nwn),min_w,r1,
     &          data_rho(-nwn:nwn)
c***************************************************************************************
	time1= secnds(0.0)
	call profstart('sumup')
c       Now normalize rho_meas and rho_Gdkk 
	rhorun = rho_meas/real(meas)
	logrhorun = logrho_meas/real(meas)
        logrholoc_run = logrholoc_meas/real(meas*Nc)
	rholoc_run=rholoc_meas/real(meas*Nc)
        
c***************************************************************************************
c  Impose point group symmetry on ADOS and LogDOS
c***************************************************************************************

	call Symm_PG_Real(logrhorun)
	call Symm_PG_Real(rhorun)
	call Symm_PG_Real(rho_avrun)
!	call Symm_PG_Complex(G_av_run)

c***************************************************************************************
c       Impose PH symmetry on the logrholoc
c***************************************************************************************
	call Symm_PH_Local(logrholoc_run)
	call Symm_PH_Local(rholoc_run)
	call Symm_PH_Real(rhorun)
	call Symm_PH_Real(rho_avrun)
!	call Symm_PH_Complex(G_av_run)


c***************************************************************************************
c --- Here we do some calcualtions needed for ERROR estimates
c       calculate the local average dos for this run
c***************************************************************************************

        rhoaveloc_run=zeror
        do i=1,Nc
          rhoaveloc_run(:) = rhoaveloc_run(:) + rhorun(i,:)
        end do
        rhoaveloc_run=rhoaveloc_run/real(Nc)

          if(ityp) then ! To enable calculation of error bar for all TMT methods
	     rhotyprun(:,:) = exp(logrhorun(:,:)) !exp((logrhorun(:,:)))  
 	else
            rhotyprun(:,:) = rhorun(:,:) 
	end if
		

c***************************************************************************************
cc	Do Hilbert transformation to calculated the Typical Green function (TGF).
c       First the real part.  Note, this is a principal value integral	
c       call profstart('HilbertT_Sumup') ! To calculate the time spent in HT.	
c***************************************************************************************
	if(errorbar) then	      				      		      				      
	call HilbertT(Gtyprun,Gtyprun_real,rhotyprun) 						      		      				      
		      				          	
cc****************** coarse graining ***************************************************				      						      
c****************************************************************************************
c	Coarse-graining the G_bar using the partial density of 
c	states calculated using look up table
c	This is an integration over energy only and as such cheaper than 
c	the k_tilde summation.
cc	Calculate Gbarrun and Gscriptrun with energy integration via partial density of states
c****************************************************************************************
c       call profstart('EnergyInegralSumup')

	Gbarrun = zeroc
	!$OMP PARALLEL DO 
	do n =-nwn,nwn
	do ic=1,Nc
            do i=0,e_div
              Gbarrun(ic,n)=Gbarrun(ic,n)+ 
     &         PDOS(i,ic)*estep/(1.d0/Gtyprun(ic,n)+
     &          Gammarun_old(ic,n)- 
     &		(estart+float(i)*estep)-ed + Epsbar(ic))
            	  end do
		Gbarrun(ic,n) = Gbarrun(ic,n)*real(Nc)
	  end do
	end do
	!$OMP END PARALLEL DO


c	Update Gamma
	!$OMP PARALLEL DO 
	do n =-nwn,nwn
	do ic=1,Nc
	Gammarun_new(ic,n) = Gammarun_old(ic,n) + 
     &	(1.d0/Gtyprun(ic,n) - 1.d0/Gbarrun(ic,n))
	  end do
	end do
	!$OMP END PARALLEL DO
ccccccc
c	Update Gscript
	if(gscript_pole) then	! This method uses delta-function to calculate Gscript close to criticality
	do ic=1,Nc
	do n=-nwn,nwn
        m_r1(ic,n)=real(wn(n)-Epsbar(ic)-Gammarun_new(ic,n))
    	end do
    	end do

	!$OMP PARALLEL DO 
    	do ic=1,Nc    ! ic loop starts
    	     min_w=abs(m_r1(ic,-nwn))
	      x_pole = -nwn
    	do n=-nwn,nwn
    	if(min_w.gt.abs(m_r1(ic,n))) then
    	min_w=abs(m_r1(ic,n))
    	x_pole=n
    	end if
    	end do
    	if((-1.d0/pi)*aimag(Gammarun_new(ic,x_pole)).ge.
     &          afact*dwn(x_pole)) then
	     do n=-nwn,nwn
	        Gscriptrun(ic,n) = 1.d0/(wn(n)-
     &          Gammarun_new(ic,n)-Epsbar(ic)+eta)
	     end do   	    
	else 
	    do n=-nwn,nwn
	     if(n .eq. x_pole) then ! Calculating Gscript when pole occurs. Approximate using delta function
    	       Gscriptrun(ic,n) = 
     &          cmplx(0.d0,-pi/dwn(x_pole))
	     else if (n .ne. x_pole) then		! wn(n) .ne. wn(n`) => Real Gscript.
               Gscriptrun(ic,n) = cmplx(1.d0/(wn(n)-
     &          real(wn(x_pole))),0.d0)                 
c               Gscriptrun(ic,n) = cmplx(1.d0/(wn(n)-
c     &          real(wn(x_pole)-Epsbar(ic)-Gammarun_new(ic,x_pole))),0.d0) 
             !real(wn(x_pole)-Epsbar(ic)-Gamma_new(ic,x_pole))
	     end if	
    	    end do
    	   end if  !end if check pole is there
    	   end do    ! ic loop ends
	!$OMP END PARALLEL DO 

	Gammarun_new(:,:)= onec/Gscriptrun(:,:) 
	end if	! gscript_pole ends	


	Gammarun_old(:,:) = Gammarun_new(:,:)
	end if

c       call profend('EnergyInegralSumup')
c***************************************************************************************
c	Load the measurements for this run into those for the iteration
c***************************************************************************************
         do ic=1,Nc
	  do n =-nwn,nwn
	   rhof(ic,n) = rhof(ic,n) + rhorun(ic,n)
	   rho_logf(ic,n) = rho_logf(ic,n)+logrhorun(ic,n)  	
         
           end do
         end do  
         
         do n =-nwn,nwn
           logrholoc_f(n)=logrholoc_f(n) + logrholoc_run(n)
           rholoc_av(n)=rholoc_av(n) + rholoc_run(n)
         end do  
        
	if(errorbar) then
c  ---  for error bar calcualtions
        rholoc_f = rholoc_f + exp(logrholoc_run(:))
        rholoc_fs = rholoc_fs + (exp(logrholoc_run(:)))**2
cc	sigmaf = sigmaf + sigmarun
	G_typf = G_typf + Gtyprun
	gammaf = gammaf+aimag(Gammarun_new)
	rhofs = rhofs + rhorun**2
	rho_logfs = rho_logfs + logrhorun**2
	gammafs = gammafs + (aimag(Gammarun_new))**2
	end if
c***************************************************************************************
c	Now zero the accumulator for a single run
c***************************************************************************************

        rho_meas=zeror
	logrho_meas = zeror
        logrholoc_meas=zeror
	rholoc_meas=zeror
	Gammarun = zeror
	Gtyprun = zeroc
	Gammarun_new=zeroc
	G_av_meas=zeroc
	rho_av_meas=zeror
!	data=zeroc
	
	rhorun=zeror
	logrhorun=zeror
	logrholoc_run=zeror
	rholoc_run=zeror
	rho_avrun=zeror
	G_av_run=zeroc

        call profend('sumup')
 	time2= secnds(time1)
        sumup_time=sumup_time+time2
	return
	end

