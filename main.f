cc**************************************************************************** 
        program fmain
	use Global
	use omp_lib
	implicit none
c**************************************************************************** 
 
c       This program simulates the 3-dimensional disorded Anderson-hubbard 
c	model in the Cluster Typical Medium Theory Approximation
c
c                        H=H1+H2
c                             +
c                H1= -t SUM ((C (i,s)C(j,s)) + h.c. )
c                       <ij>sproperty of
c                               +  
c                    -t' SUM ((C (i,s)C(j,s)) + h.c. )
c                       [ij]s                           i,j n.n.n.
c
c                H2=   SUM (ed + V ) * n(i,s)
c                       i         i    
c
c	V is a random binary variable with <V>=0.  H is particle-hole 
c	symmetric whenever ed=0 and t'=0.  For a two-dimensional system, 
c	this means that 
c
c               +             +    i
c       C  --> C exp(iQ.i) = C (-1)
c        i      i             i
c                                                +              +
c       C  = SUM exp(ik.i) C  --> SUM exp(ik.i) C exp(-iQ.i) = C
c        k    i             i                    i              Q-k
c
c       or if ed=0, G (tau) = -G   (-tau) = G   (beta-tau)
c                    k          Q-k          Q-k
c       and
c                   G (iwn) = - G   (-iwn)  real parts opposite sign
c                    k           Q-k        imaginary parts equal sign
c
c	in real-frequencies, this means that 
c
c                    R          A            *
c                   G (w) = - G   (-w)  =  -G   (-w)
c                    k         Q-k           Q-k
c

c**************************************************************************** 
c	Original serial code by Mark Jarrell 1999-2000
c
c	The modification to implement Typical Medium Theory within a cluster (CTMT) is
c	by Chinedu Ekuma, Zi Yang Meng(1d & 2d), Hanna Terletska and Mark Jarrell,  2012--2014
c
c       Checks:
c       ( ) cat *.f > onefile ; mv onefile onefile.f
c           ftnchek -usage=303 -portability=all onefile.f
c
c
c**************************************************************************** 

        integer n,m,i,j,ic
        real(kind) :: time_main1,time_main2,r1,r2,r3
        character*72 linemc
c**************************************************************************** 
        call profinit()
        call profstart('main')
	time_main1= secnds(0.0)
c	time1= secnds(0.0)
c	wtime= OMP_get_wtime()
c	Write out the user info
	call userinfo(lud)
c       open the output file
        call readin
c	Allocate arrays
	call allocate_arrays
c       we must now make the lookup tables
        call tables
 	write(6,*) '***********************************'					      
	write(6,*) "table completed successfully."
	call init_parameters
	call Bin_Header
	isf=0
	test=zeror	!convergence criteria
	testl=zeror
        histo_meas=zeror

	do 999 iter=1,niter					
	  if(isi.ge.1.and.iter.gt.1) meas=int(meas*ifac)		
c	  Added to enable the sampling of Large meas for histogram
	  if(histo_enable) then
	    if(isi.ge.1.and.iter.gt.1.and.iter.ne.niter) meas=int(meas*ifac)		
	    if(isi.ge.1.and.iter.eq.niter) meas=int(meas*10)	
	  end if ! End histogram 	
	  write(6,*) 'iteration',iter,' runs=',run,' meas=',meas
	  write(42,*) 'iteration',iter,' runs=',run,' meas=',meas

c         Make Gcsfsc, the greens function for all V(l)=0 and initialize all the accumulators
          call geodgen
          call ginit        ! Initialize the greens functions and the fields.
 
c         The following is the main part of the main program.  It directs
c         measurement, and updating the lattice and the d greens functions
c         between the measurements.
          do 30 nrun=1,run
            do 40 nmeas=1,meas
              call update
              call measure  
 40         continue
c           Now load the data into bins and accumulators.
            call sumup
 30       continue

          call output(1)


	if(iter.eq.niter) then
 	write(6,*) "Writing Data to Output"
	end if
        time_main2= secnds(time_main1)						      

	  if(iter.eq.niter-1) then
 	    write(6,*) "Last Iteration. Major output will be written into dos.dat"
 	    write(6,*) "The TotalDOS and PDOS will be written into p_rho.dat"
 	    write(6,*) "To understand other outputs, see README, output section"
 	    write(6,*) '*********************************************************'						      
 	    write(42,*) "Last Iteration. Major output will be written into dos.dat"
 	    write(42,*) "The TotalDOS and PDOS will be written into p_rho.dat"
 	    write(42,*) "To understand the many outputs, see README output section"
 	    write(42,*) '*********************************************************'		      
	  end if


	if(iter.eq.niter) then
	if(interpolation) then
	call spectra
	end if
c	Write out the date and time information for when the iteration finished
          write(lud,*)   '******************************************************************************'
          write(lud,*)   'Time for when the Calculation Ended                                         *'
              call date_and_time(date,time,zone,values)
              call date_and_time(DATE=date,ZONE=zone)
              call date_and_time(TIME=time)
              call date_and_time(VALUES=values)
	      write(lud,*) '*DATE(Y-M-D)  TIME(H-M-S.MS)   TIME-ZONE      			       	   *'                    
              write(lud, '(1x,a,3x,a,3x,a)') date, time, zone 
	      write(lud,*) '*Y  M   D T-DIFF(UTC) H MIN SEC MILLIS         				*'                    
              write(lud,'(8i5)') values     
          write(lud,*)   '******************************************************************************'


 	  write(42,*) '*******************************************************'
 	  write(6,*) '********************************************************'		      		      
	  write(lud,*) 'done iteration',iter,'at time',time_main2
	  write(6,*) 'done iteration',iter,'at time',time_main2
          write(6,*) 'cpu  ,time=  ',time_main2                    
          write(lud,*) 'cpu  ,time=  ',time_main2                   
          write(lud,*) 'geod ,time=  ',geod_time                     
          write(lud,*) 'update ,time=  ',update_time                            
          write(lud,*) 'meas ,time=  ',meas_time                        	
          write(lud,*) 'put  ,time=  ',put_time                         	
          write(lud,*) 'sumup,time=  ',sumup_time                       	
 	  write(42,*) '******************************************************'		      
 	  write(6,*) '*******************************************************'		      

        end if

 999	continue	! done iteration loop

	call deallocate_arrays
c        time2= secnds(time1)
c	wtime= OMP_get_wtime()-wtime
c	write(6,*) 'cpu ,time=  ',time2
c        write(42,*) 'cpu ,time=  ',time2

c	write(6,*) 'cpu ,time=  ',time2
c	write(42,*) 'cpu ,time=  ',time2

        call profend('main')
        call profstat()
        stop
        end
