    program spectra

!       to compile type make spectra (or make all) using the
!       associated Makefile.  This code was designed to read the
!       magority of its inputs from a qmc sigma file.
!	This code uses the trilinear and star interpolation method.
!	It is for 1,2,and 3 D

    use Global
    implicit none
    integer :: filenum, istat, i, j,l, ic,jc,n,ickx,icky,ickz,ikz, &
    nkpoints, info, infot, ikx, iky, deg, itunnel,i_energy, & 
    ISPLINE,i1,j1,i2,j2,l1,l2,ic1,ic2,ic3,ic4,m,p, &
    ic5,ic6,ic7,ic8,jck
    real(kind) :: r1, r2, r3, k,kx, ky, kz,fill, Epsk, kx1, ky1,kx2,k1,k2, &
    ky2,offset, res, hall, vx, vyy, dfermi,kz1,kz2,Energy, &
    dx,dy,dz,Kx3,Ky3,Kz3,r7,r8,ksym(3,8)
    complex(kind) :: c1, c2, blint,z1,z2,z3,z4,z5,z6,z7,z8,c_lat,c_lat_EBS,c_Gonly,blint_spline
    complex(kind), allocatable :: c1_sigma(:,:),c_Gl1(:,:),c_Gl2(:,:),c_Gext(:,:), &
    c2_sigma(:,:),c_Gl(:,:),c_latmod(:,:),c2_sigma_EBS(:,:), &
    c_Gl_EBS(:,:),c1_sigma_EBS(:,:)
    		    
    character(72) :: linemc,fmt1,string1
    CHARACTER F*128, FN*128, FILENAME*128
    external blint,blint_spline

!###########################################
!	For Bandstructure
!##########################################

    integer, parameter::maxkp=3000
    integer:: nkp, mend, ndiv,u,npi
    real(kind), dimension(:, :), allocatable::kp
    real(kind), dimension(:), allocatable::kr
    real(kind)::c, s, r, weit
    real(kind)::kx_b, ky_b, kz_b, dist
    real(kind)::endkp(maxkp, 3), dkp(3), b3(3, 3)
    allocate(kp(maxkp, 3), kr(maxkp))
    write(6,*) '***************************************************************'
    write(6,*) 'Please, enter the needed input. Basically self-energy'
    write(6,*) 'Self-energy file is read directly from sigma.dat'
    write(6,*) 'You can rename it and change the name accordingly'
    write(6,*) 'Please, read off the W_e from the Mobility Edge'
    pause 'PRESS ANY KEY TO CONTINUE'
    open (unit=40, file="sigma.dat", status='old', &
    form='formatted', action='read' )
!    open(unit=10, file="bandkpt.out", action="read") ! High symmetry files
!    open(unit=50, file="chicharge.dat", action="read") 
    1000 format(a72)
    401 read(40,1000) linemc
    if(index(linemc(1:72),'Results') /= 0 ) then
        read(40,*) nwn,Nc,nover,cluster
    else
        goto 401
    end if
!################################
!	Set the flag for spline only for case(3) in module_Global
	ispline=0
       acoef = 0.5
	icum=0.0 ! Be careful with the cummulant
!################################
    404 read(40,1000) linemc
    if(index(linemc(1:72),'ed,V,tprime') /= 0 ) then
        read(40,*) ed,V0,tprime,niter,run,meas !Only ed may be needed here but the out has others
    else
        goto 404
    end if

!    405 read(40,1000) linemc
!    if(index(linemc(1:72),'ca Hubbard-U') /= 0 ) then
!        read(40,*) ca,Uint,beta,ph_asym,gsflag,ndim !Only ed may be needed here but the out has others
!    else
!        goto 405
!    end if


!###############################################################################
! This part is for bandstructure grid is needed in many high-symmetry points along
!	the defined high symmetry points
!###############################################################################

!    407 read(10,1000) linemc
!    if(index(linemc(1:72),'total') /= 0 ) then
!        read(10,*) nkp
!    else
!        goto 407
!    end if

!    2006 read(10,1000) linemc
!    do i=1,nkp ! Specify the no of k-points from the bandkpt file
!        if(index(linemc(1:72),'Bandkpoints') /= 0 ) then
!            read(10,*,IOSTAT=istat) m,kp(i, 1), kp(i, 2), kp(i, 3),weit
!            if(istat /= 0) goto 999
!        else
!            goto 2006
!        endif
!        999 continue
!	enddo

!       Read in some parameters
	if(ndim==1) then
	N_sp=Nc/2.
	N_sp = N_sp+4.
	elseif(ndim==2) then
	N_sp=sqrt(float(Nc))
	N_sp=N_sp+4.
	elseif(ndim==3) then
	N_sp=sqrt(float(Nc))
	N_sp=N_sp+10
	end if
	
!	    Energy= 3.50D0
    nkpoints = 1000! This is the maximum. The input nkpoints overides this
!	    nf = 300
!	Give the frequency to nf
    nf = nwn

  			       				 	        	
    call allocate_arrays
    call tables

!#################################################
!Read in the Cluster Self-energy
!#################################################
    2004 read(40,1000) linemc
    do ic=1,Nc
        if(index(linemc(1:72),'sigma_loc') /= 0 ) then
            do i=-nwn,nwn
                read(40,*,IOSTAT=istat) m,m,wn(i),r2,r3
                if(istat /= 0) goto 997
            !	  	read(40,*,end=997) wn(n),m,r1,r2,r3
                sigma_loc(ic,i)=cmplx(r2,r3)
            end do
        else
            goto 2004
        endif
        997 continue

        if(ic == 1) then
            nf=i-1
        else if(nf /= i-1) then       ! trap for errors
            write(6,*) 'Sigma files are not compatible'
            stop
        end if
    end do

!#################################################
!Read in the Average Green function
!#################################################
    2009 read(40,1000) linemc
    do ic=1,Nc
        if(index(linemc(1:72),'GcKf_av') /= 0 ) then
            do i=-nwn,nwn
                read(40,*,IOSTAT=istat) m,m,wn(i),r2,r3
                if(istat /= 0) goto 998
                GcKf_av(ic,i)=cmplx(r2,r3)
            end do
        else
            goto 2009
        endif
        998 continue

        if(ic == 1) then
            nf=i-1
        else if(nf /= i-1) then       ! trap for errors
            write(6,*) 'GcKf_av files are not compatible'
            stop
        end if
    end do



!#################################################
!Read in the chicharge_meas
!#################################################
!    2007 read(50,1000) linemc
!    do ick=1,Nc
!	do jck=1,Nc
!!        if(index(linemc(1:72),'ick') /= 0 ) then
!            do i=-nwn,nwn
!                read(50,*,IOSTAT=istat) m,m,m,r2,r3
!                if(istat /= 0) goto 1001
!                chicharge_meas(ick,jck,i)=cmplx(r2,r3)
!            end do
!!        else
!            goto 2007
!        endif
!        1001 continue
!    end do
!    enddo
!!#################################################
!       Make sure that nf is le the limit set above
!#################################################
    if(nf > 1000) then ! Force it to reduce computational cost
        write(6,*) 'nf in the data file is too large. I have adjusted it'
        nf = 1000
    end if
    	
!############################################################################################
!       Now do the interpolation, and calculate the spectrum at
!       a sequence of locations.
! Note the difference between the spectra data and energy_surface. The energy surface is
! at constant energy.
!############################################################################################
    if(ndim==1) then
        write(6,"('Enter k1,k2,offset,nkpoints,w_e:  ',$)")
        read(5,*) k1,k2,offset,nkpoints,w_e
    elseif(ndim==2) then
        write(6,"('Enter kx1,ky1,kx2,ky2,offset,nkpoints,w_e:  ',$)")
        read(5,*) kx1,ky1,kx2,ky2,offset,nkpoints,w_e
    else if(ndim==3) then
        !write(6,"('Enter kx1,ky1,kz1,kx2,ky2,kz2,offset,nkpoints,w_e:  ',$)")
        !read(5,*) kx1,ky1,kz1,kx2,ky2,kz2,offset,nkpoints,w_e
        write(6,*) 'Enter the mobility edge'
        read(5,*) w_e
    end if

!	Now that we have inputed all the data, allocate the internal arrays as to be the same dimension as input
    allocate(c1_sigma(1:nkpoints,-nf:nf),c_Gl1(1:nkpoints,-nf:nf),c_Gl2(1:nkpoints,-nf:nf),c_Gext(1:nkpoints,-nf:nf), &
    c2_sigma(1:nkpoints,-nf:nf),c_Gl(1:nkpoints,-nf:nf),c_latmod(1:nkpoints,-nf:nf), &
    c2_sigma_EBS(1:nkp,-nf:nf),c1_sigma_EBS(1:nkp,-nf:nf),c_Gl_EBS(1:nkp,-nf:nf))

!	Fin the W_e on the nf grid
!	Both is the same
!        do n=0,nwn
!	nwn_c = n
!	if((wn(n)-w_e)>0) exit
!       end do

    do n=0,nwn
        if((wn(n)-w_e)>0) exit
    enddo
    nwn_c = n

!       Call two-particle code here for debugging

    call analysis_twoparticle
    stop

!         set symmetrize table for kx,ky,kz. For Betts lattice mainly for 3D
          ksym(1,1)= 1.0; ksym(2,1)= 1.0; ksym(3,1)= 1.0
          ksym(1,2)=-1.0; ksym(2,2)= 1.0; ksym(3,2)= 1.0
          ksym(1,3)=-1.0; ksym(2,3)=-1.0; ksym(3,3)= 1.0
          ksym(1,4)=-1.0; ksym(2,4)=-1.0; ksym(3,4)=-1.0
          ksym(1,5)=-1.0; ksym(2,5)= 1.0; ksym(3,5)=-1.0
          ksym(1,6)= 1.0; ksym(2,6)= 1.0; ksym(3,6)=-1.0
          ksym(1,7)= 1.0; ksym(2,7)=-1.0; ksym(3,7)=-1.0
          ksym(1,8)= 1.0; ksym(2,8)=-1.0; ksym(3,8)= 1.0




!###################################################
!	Prepare the self-energy
!###################################################

    sigma_temp(:,:)=sigma_loc(:,:)

!###################################################
!	Prepare the G_local
!###################################################
    G_loc=zeroc
    do n =-nwn,-nwn_c
        do ic=1,Nc
            G_loc(ic,n) = 1.0D0/(wn(n)-Epsbar(ic)-sigma_temp(ic,n))
        enddo
    enddo


    do n = nwn_c+1,nwn
        do ic=1,Nc
            G_loc(ic,n) = 1.0D0/(wn(n)-Epsbar(ic)-sigma_temp(ic,n))
        enddo
    enddo
	


    open(unit=21,file='Akw_spectra.dat',action='write',status='unknown')
    open(unit=22,file='False_color.dat',action='write',status='unknown')
    open(unit=23,file='False_BES.dat',action='write',status='unknown')
    write(6,*) '******************************************************'
    write(6,*) '*** A(k,w) is in Akw_spectra.dat             *********'
    write(21,"('#       w      A(k,w)=-Im G(k,w)/pi ')")
    write(22,"('#      k      w        A(k,w)=-Im G(k,w)/pi ')")
    write(23,"('#      k      w        A(k,w)=-Im G(k,w)/pi ')")


!	Plot the spectra using xmgrace
    open (unit=8,status='unknown',file='Awk_spectra.agr')
    write(8,"('@    default font 0')")
    write(8,90)
    90 format('@    yaxis  label "A(K,\xw\N)"')
    write(8,91)
    91 format('@    xaxis  label "\xw\N"')

    if(ndim==1) then
    !	Interpolate the self-energy first on all frequencies.
	if(ispline==1) then
	call spline_1_3D(sigma_temp)
        do ick=1,nkpoints
            k=k1 + real(ick-1)*(k2-k1)/max(nkpoints-1,1)
            do i=-nf,nf
                c2_sigma(ick,i)=blint(0,k,0,0,i)  
            enddo
	enddo
!	Interpolate the local Green function
	call spline_1_3D(G_loc)
        do ick=1,nkpoints
            k=k1 + real(ick-1)*(k2-k1)/max(nkpoints-1,1)
            do i=-nf,nf
                c_Gl(ick,i)=blint(0,k,0,0,i)
            enddo
	enddo

	else
        do ick=1,nkpoints
            k=k1 + real(ick-1)*(k2-k1)/max(nkpoints-1,1)
            do i=-nf,nf
                c2_sigma(ick,i)=blint(sigma_temp,k,0,0,i)
                c_Gl(ick,i)=blint(G_loc,k,0,0,i)

	enddo
	enddo
	end if

        do ick=1,nkpoints
            k=k1 + real(ick-1)*(k2-k1)/max(nkpoints-1,1)
        !	Restrict G_ext to the -nwn_c and nwn_c to calculate the G_ext
            do i=-nwn_c,nwn_c
                c1_sigma(ick,i)=1.d0/(wn(i)-Epsk(k,0,0)-c2_sigma(ick,i))
            enddo
        !	Next, combine G_ext and G_loc to get G_lattice
            do i=-nf,nf
                c_lat = c1_sigma(ick,i)+c_Gl(ick,i)
                write (8,"(1x,f15.6,', ',f15.6)") wn(i), &
                -dimag(c_lat)/pi+(ick-1)*offset
                write(21,"(5(f15.6,1x))") wn(i), &
                -dimag(c_lat)/pi+(ick-1)*offset
                write(22,"(i4,1x,2(f15.6,1x))") ick,wn(i),-dimag(c_lat)/pi
            end do

            write(8,*) ' & '
            write(21,*) ' & '
            write(22,*) '  '
        end do
    !##############################################################################
    ! This part uses a special grid that gives the band structure as specified above
    !##############################################################################
	if(ispline==1) then
!	Interpolate the self-energy on all frequencies
	call spline_1_3D(sigma_temp)
        do ick=1,nkp
            kx=kp(ick,1);ky=kp(ick,2);kz=kp(ick,3)
            do i=-nf,nf
                c2_sigma_EBS(ick,i)=blint(0,kx, 0,0,i)
            enddo
	enddo
!	Interpolate the local Green function (G_loc)
	call spline_1_3D(G_loc)
        do ick=1,nkp
            kx=kp(ick,1);ky=kp(ick,2);kz=kp(ick,3)
            do i=-nf,nf
                c_Gl_EBS(ick,i)=blint(0,kx, 0,0,i)	
            enddo
	enddo
	else
        do ick=1,nkp
            kx_b=kp(ick,1);ky_b=kp(ick,2);kz_b=kp(ick,3)
            do i=-nf,nf
                c2_sigma_EBS(ick,i)=blint(sigma_temp,kx_b,0,0,i)
                c_Gl_EBS(ick,i)=blint(G_loc,kx_b,0,0,i)
            enddo
	enddo

	end if
!	Calculate the Extended Green function(G_ext)
        do ick=1,nkp
            kx=kp(ick,1);ky=kp(ick,2);kz=kp(ick,3)
            do i=-nwn_c,nwn_c
                c1_sigma_EBS(ick,i)=1.d0/(wn(i)-Epsk(kx,0,0)-c2_sigma_EBS(ick,i))
            enddo
        enddo

        !	Next, combine G_ext and G_loc to get G_lattice
        do ick=1,nkp
            kx=kp(ick,1);ky=kp(ick,2);kz=kp(ick,3)
            do i=-nf,nf
                c_lat_EBS= c1_sigma_EBS(ick,i)+c_Gl_EBS(ick,i)
                write(23,"(i4,1x,2(f15.6,1x))") ick,wn(i),-dimag(c_lat_EBS)/pi
            end do


            write(23,*) ' '
        enddo
    elseif(ndim==2) then
    !	Interpolate the self-energy first on all frequencies.
	if(ispline==1) then
	call spline_1_3D(sigma_temp)
        do ick=1,nkpoints
            kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
            ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)
            do i=-nf,nf
                c2_sigma(ick,i)=blint(0,kx,ky,0,i)  
            enddo
	enddo
!	Interpolate the local Green function
	call spline_1_3D(G_loc)
        do ick=1,nkpoints
            kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
            ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)
            do i=-nf,nf
                c_Gl(ick,i)=blint(0,kx,ky,0,i)
            enddo
	enddo

	else
        do ick=1,nkpoints
            kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
            ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)
            do i=-nf,nf
                c2_sigma(ick,i)=blint(sigma_temp,kx,ky,0,i)
                c_Gl(ick,i)=blint(G_loc,kx,ky,0,i)

	enddo
	enddo
	end if

        do ick=1,nkpoints
            kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
            ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)
        !	Restrict G_ext to the -nwn_c and nwn_c to calculate the G_ext
            do i=-nwn_c,nwn_c
                c1_sigma(ick,i)=1.d0/(wn(i)-Epsk(kx,ky,0)-c2_sigma(ick,i))
            enddo
        !	Next, combine G_ext and G_loc to get G_lattice
            do i=-nf,nf
                c_lat = c1_sigma(ick,i)+c_Gl(ick,i)
                write (8,"(1x,f15.6,', ',f15.6)") wn(i), &
                -dimag(c_lat)/pi+(ick-1)*offset
                write(21,"(5(f15.6,1x))") wn(i), &
                -dimag(c_lat)/pi+(ick-1)*offset
                write(22,"(i4,1x,2(f15.6,1x))") ick,wn(i),-dimag(c_lat)/pi
            end do

            write(8,*) ' & '
            write(21,*) ' & '
            write(22,*) '  '
        end do
    !##############################################################################
    ! This part uses a special grid that gives the band structure as specified above
    !##############################################################################
	if(ispline==1) then
!	Interpolate the self-energy on all frequencies
	call spline_1_3D(sigma_temp)
        do ick=1,nkp
            kx=kp(ick,1);ky=kp(ick,2);kz=kp(ick,3)
            do i=-nf,nf
                c2_sigma_EBS(ick,i)=blint(0,kx, ky, 0,i)
            enddo
	enddo
!	Interpolate the local Green function (G_loc)
	call spline_1_3D(G_loc)
        do ick=1,nkp
            kx=kp(ick,1);ky=kp(ick,2);kz=kp(ick,3)
            do i=-nf,nf
                c_Gl_EBS(ick,i)=blint(0,kx, ky, 0,i)	
            enddo
	enddo
	else
        do ick=1,nkp
            kx_b=kp(ick,1);ky_b=kp(ick,2);kz_b=kp(ick,3)
            do i=-nf,nf
                c2_sigma_EBS(ick,i)=blint(sigma_temp,kx_b,ky_b,0,i)
                c_Gl_EBS(ick,i)=blint(G_loc,kx_b, ky_b,0,i)
            enddo
	enddo

	end if
!	Calculate the Extended Green function(G_ext)
        do ick=1,nkp
            kx=kp(ick,1);ky=kp(ick,2);kz=kp(ick,3)
            do i=-nwn_c,nwn_c
                c1_sigma_EBS(ick,i)=1.d0/(wn(i)-Epsk(kx,ky,0)-c2_sigma_EBS(ick,i))
            enddo
        enddo

        !	Next, combine G_ext and G_loc to get G_lattice
        do ick=1,nkp
            kx=kp(ick,1);ky=kp(ick,2);kz=kp(ick,3)
            do i=-nf,nf
                c_lat_EBS= c1_sigma_EBS(ick,i)+c_Gl_EBS(ick,i)
                write(23,"(i4,1x,2(f15.6,1x))") ick,wn(i),-dimag(c_lat_EBS)/pi
            end do


            write(23,*) ' '
        enddo
    elseif(ndim==3) then
    !	Interpolate the self-energy first on all frequencies.
	if(ispline==1) then
	call spline_1_3D(sigma_temp)
        do ick=1,nkpoints
            kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
            ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)
            kz=kz1 + real(ick-1)*(kz2-kz1)/max(nkpoints-1,1)
            do i=-nf,nf
             do j=1,8
                c2_sigma(ick,i)=blint(0,kx*ksym(1,j),ky*ksym(2,j),kz*ksym(3,j),i)  
            enddo
	enddo
	enddo
!	Interpolate the local Green function
	call spline_1_3D(G_loc)
        do ick=1,nkpoints
            kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
            ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)
            kz=kz1 + real(ick-1)*(kz2-kz1)/max(nkpoints-1,1)
            do i=-nf,nf
	do j=1,8
                c_Gl(ick,i)=blint(0,kx*ksym(1,j),ky*ksym(2,j),kz*ksym(3,j),i)
            enddo
	enddo
	enddo
	else
        do ick=1,nkpoints
            kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
            ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)
            kz=kz1 + real(ick-1)*(kz2-kz1)/max(nkpoints-1,1)
            do i=-nf,nf
	do j=1,8
                c2_sigma(ick,i)=blint(sigma_temp,kx*ksym(1,j),ky*ksym(2,j),kz*ksym(3,j),i)
                c_Gl(ick,i)=blint(G_loc,kx*ksym(1,j),ky*ksym(2,j),kz*ksym(3,j),i)
	enddo

	enddo
	enddo
	end if

        do ick=1,nkpoints
            kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
            ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)
            kz=kz1 + real(ick-1)*(kz2-kz1)/max(nkpoints-1,1)
        !	Restrict G_ext to the -nwn_c and nwn_c to calculate the G_ext
            do i=-nwn_c,nwn_c
	    do j=1,8
                c1_sigma(ick,i)=1.d0/(wn(i)-Epsk(kx*ksym(1,j),ky*ksym(2,j),kz*ksym(3,j))-c2_sigma(ick,i))
            enddo
	enddo
        !	Next, combine G_ext and G_loc to get G_lattice
            do i=-nf,nf
                c_lat = c1_sigma(ick,i)+c_Gl(ick,i)
                write (8,"(1x,f15.6,', ',f15.6)") wn(i), &
                -dimag(c_lat)/pi+(ick-1)*offset
                write(21,"(5(f15.6,1x))") wn(i), &
                -dimag(c_lat)/pi+(ick-1)*offset
                write(22,"(i4,1x,2(f15.6,1x))") ick,wn(i),-dimag(c_lat)/pi
            end do

            write(8,*) ' & '
            write(21,*) ' & '
            write(22,*) '  '
        end do
    !##############################################################################
    ! This part uses a special grid that gives the band structure as specified above
    !##############################################################################
	if(ispline==1) then
!	Interpolate the self-energy on all frequencies
	call spline_1_3D(sigma_temp)
        do ick=1,nkp
            kx=kp(ick,1);ky=kp(ick,2);kz=kp(ick,3)
            do i=-nf,nf
	do j=1,8
                c2_sigma_EBS(ick,i)=blint(0,kx*ksym(1,j),ky*ksym(2,j),kz*ksym(3,j),i)
            enddo
	enddo
	enddo
!	Interpolate the local Green function (G_loc)
	call spline_1_3D(G_loc)
        do ick=1,nkp
            kx=kp(ick,1);ky=kp(ick,2);kz=kp(ick,3)
            do i=-nf,nf
	do j=1,8
                c_Gl_EBS(ick,i)=blint(0,kx*ksym(1,j),ky*ksym(2,j),kz*ksym(3,j),i)	
            enddo
	enddo
	enddo
	else
        do ick=1,nkp
            kx=kp(ick,1);ky=kp(ick,2);kz=kp(ick,3)
            do i=-nf,nf
	do j=1,8
                c2_sigma_EBS(ick,i)=blint(sigma_temp,kx*ksym(1,j),ky*ksym(2,j),kz*ksym(3,j),i)
                c_Gl_EBS(ick,i)=blint(G_loc,kx*ksym(1,j),ky*ksym(2,j),kz*ksym(3,j),i)
            enddo
	enddo
	enddo
	end if
!	Calculate the Extended Green function(G_ext)
        do ick=1,nkp
            kx=kp(ick,1);ky=kp(ick,2);kz=kp(ick,3)
            do i=-nwn_c,nwn_c
	do j=1,8
                c1_sigma_EBS(ick,i)=1.d0/(wn(i)-Epsk(kx*ksym(1,j),ky*ksym(2,j),kz*ksym(3,j))-c2_sigma_EBS(ick,i))
            enddo
        enddo

	enddo
        !	Next, combine G_ext and G_loc to get G_lattice
        do ick=1,nkp
            kx=kp(ick,1);ky=kp(ick,2);kz=kp(ick,3)
            do i=-nf,nf
                c_lat_EBS= c1_sigma_EBS(ick,i)+c_Gl_EBS(ick,i)
                write(23,"(i4,1x,2(f15.6,1x))") ick,wn(i),-dimag(c_lat_EBS)/pi
            end do


            write(23,*) ' '
        enddo
        	

    endif
    close(8)
    close(21)
    close(22)
    close(23)


    		   
    stop
    END PROGRAM
    		            

    complex*16 function blint(GcgSigm,kx,ky,kz,n)

!       This subroutine performs a trilinear interpolation of
!       the self energy at wn, to calculate
!       the self energy at any location (kx,ky,kz) in the zone.
!	It is for 1,2,and 3D.
!	Various interpolation schemes are included. They are:
!	Star interpolation, trilinear interpolation, and spline interpolation
!	Note, there is no spline interpolation for 3D.

    use Global
    implicit none
    integer :: filenum, istat, i,j,ic,n,ickx,icky,ickz,ikz,ik, &
    nkpoints, info, infot, ikx, iky, deg, itunnel,i_energy, &  
    ISPLINE,i1,j1,i2,j2,l1,l2,ic1,ic2,ic3,ic4, &
    ic5,ic6,ic7,ic8
    real(kind) :: r1, r2, r3, k,kx, ky, kz,fill, Epsk, kx1, ky1,kx2,k1,k2, &
    ky2,offset,vx, vyy,kz1,kz2,Energy,dx,dy,dz,Kx3,Ky3,Kz3,r7,r8,g1,g2,g3
    real(8) :: delta_g, delta_g1, delta_g2, delta_g3,resapprox,imsapprox
    		    
    complex(kind) :: c1, c2,z1,z2,z3,z4,z5,z6,z7,z8,Sapprox(1:Nc),GcgSigm(1:Nc,-nwn:nwn),ztemp
    complex(kind) :: m1,m2,m3,m4,m5,m6,m7
    select case(interpolation_type)
    case(1) !if(trilinear_interpolation) then
!	    if(trilinear_interpolation) then
    if(ndim==1) then
    !         Use a one-dimensional linear interpolation within the square

    !         1*        *2


    !         First find the coordinates of the cell which contains (kx)
        i1=int(kx/gvector(1,1)+Nc)-Nc
        i2=i1+1
    !         now get the K-point on the left.
        Kx3=i1*gvector(1,1)
        z1=GcgSigm(i1,n)
        z2=GcgSigm(i2,n)
    !         calculate the slopes
        dx=(kx-Kx3)/gvector(1,1)
    !         interpolate.
        blint=z1+(z2-z1)*dx
	r1=dimag(blint)	!Add damping
	r2=real(blint)
          if(r1.gt.-damp) r1=-damp    
	blint= dcmplx(r2,r1)
    elseif(ndim==2) then
    !         Use a two-dimensional bilinear interpolation within the square


    !         1*        *2




    !         3*        *4

    !         First find the coordinates of the cell which contains (kx,ky)
        i1=int((kx*gvector(1,1)+ky*gvector(1,2))/(gvector(1,1)*gvector(1,1)+gvector(1,2)*gvector(1,2))+Nc) &
        -Nc
        i2=i1+1
        j1=int((kx*gvector(2,1)+ky*gvector(2,2))/(gvector(2,1)*gvector(2,1)+gvector(2,2)*gvector(2,2))+Nc) &
        -Nc
        j2=j1+1
    !         now get the K-point on the lower left.
        Kx3=i1*gvector(1,1)+j1*gvector(2,1)
        Ky3=i1*gvector(1,2)+j1*gvector(2,2)
    !         index the four corners
        ic1=ict(i1,j2,0)
        ic2=ict(i2,j2,0)
        ic3=ict(i1,j1,0)
        ic4=ict(i2,j1,0)
        z1=GcgSigm(ic1,n)
        z2=GcgSigm(ic2,n)
        z3=GcgSigm(ic3,n)
        z4=GcgSigm(ic4,n)
    !         calculate the slopes
        dx=((kx-Kx3)*gvector(1,1)+(ky-Ky3)*gvector(1,2)) &
        /(gvector(1,1)*gvector(1,1)+gvector(1,2)*gvector(1,2))
        dy=((kx-Kx3)*gvector(2,1)+(ky-Ky3)*gvector(2,2)) &
        /(gvector(2,1)*gvector(2,1)+gvector(2,2)*gvector(2,2))
    !         interpolate.
        blint=z3+(z4-z3)*dx+(z1-z3)*dy+(z2-z1-z4+z3)*dx*dy

	r1=dimag(blint)	!Add damping
	r2=real(blint)
          if(r1.gt.-damp) r1=-damp    
	blint= dcmplx(r2,r1)
    elseif(ndim==3) then
    !         Use a three-dimensional trilinear interpolation within a box
    !             8______7
    !            /|     /|
    !           / |    / |
    !          4-------3 |
    !          |  5___ |_|6
    !          | /     | /
    !          |/      |/
    !          1-------2


    !         First find the coordinates of the cell which contains (kx,ky,kz)
        i1=int((kx*gvector(1,1)+ky*gvector(1,2)+kz*gvector(1,3))/ &
        (gvector(1,1)*gvector(1,1)+gvector(1,2)*gvector(1,2)+gvector(1,3)*gvector(1,3))+Nc)-Nc
        i2=i1+1
        j1=int((kx*gvector(2,1)+ky*gvector(2,2)+kz*gvector(2,3))/ &
        (gvector(2,1)*gvector(2,1)+gvector(2,2)*gvector(2,2)+gvector(2,3)*gvector(2,3))+Nc)-Nc
        j2=j1+1
        l1=int((kx*gvector(3,1)+ky*gvector(3,2)+kz*gvector(3,3))/ &
        (gvector(3,1)*gvector(3,1)+gvector(3,2)*gvector(3,2)+gvector(3,3)*gvector(3,3))+Nc)-Nc
        l2=l1+1
    !         now get the K-point on the lower left.
        Kx3=i1*gvector(1,1)+j1*gvector(2,1)+l1*gvector(3,1)
        Ky3=i1*gvector(1,2)+j1*gvector(2,2)+l1*gvector(3,2)
        Kz3=i1*gvector(1,3)+j1*gvector(2,3)+l1*gvector(3,3)
    !         index the four corners
        ic1=ict(i1,j1,l1)
        ic2=ict(i2,j1,l1)
        ic3=ict(i1,j2,l1)
        ic4=ict(i1,j1,l2)
        ic5=ict(i2,j1,l2)
        ic6=ict(i1,j2,l2)
        ic7=ict(i2,j2,l1)
        ic8=ict(i2,j2,l2)
        z1=GcgSigm(ic1,n)
        z2=GcgSigm(ic2,n)
        z3=GcgSigm(ic3,n)
        z4=GcgSigm(ic4,n)
        z5=GcgSigm(ic5,n)
        z6=GcgSigm(ic6,n)
        z7=GcgSigm(ic7,n)
        z8=GcgSigm(ic8,n)

    !         calculate the slopes
        dx=((kx-Kx3)*gvector(1,1)+(ky-Ky3)*gvector(1,2)+(kz-Kz3)*gvector(1,3))/ &
        (gvector(1,1)*gvector(1,1)+gvector(1,2)*gvector(1,2)+gvector(1,3)*gvector(1,3))
        dy=((kx-Kx3)*gvector(2,1)+(ky-Ky3)*gvector(2,2)+(kz-Kz3)*gvector(2,3))/ &
        (gvector(2,1)*gvector(2,1)+gvector(2,2)*gvector(2,2)+gvector(2,3)*gvector(2,3))
        dz=((kx-Kx3)*gvector(3,1)+(ky-Ky3)*gvector(3,2)+(kz-Kz3)*gvector(3,3))/ &
        (gvector(3,1)*gvector(3,1)+gvector(3,2)*gvector(3,2)+gvector(3,3)*gvector(3,3))
    !         interpolate.

	m1=z1*(1.0D0-dx)+z2*dx ! See wikipedia http://en.wikipedia.org/wiki/Trilinear_interpolation
	m2=z3*(1.0D0-dx)+z7*dx
	m3=z4*(1.0D0-dx)+z5*dx
	m4=z6*(1.0D0-dx)+z8*dx

	m5=m1*(1.0D0-dy)+m2*dy
	m6=m3*(1.0D0-dy)+m4*dy
	
        blint=m5*(1-dz)+m6*dz
!        blint=z1*(1.0d0-dx)*(1.0d0-dy)*(1.0d0-dz)+ &
!        z2*       dx *(1.0d0-dy)*(1.0d0-dz)+ 	   &
!        z3*       (1.0D0-dx) *dy* (1.0D0-dz)+ 	   &
!        z4*(1.0d0-dx)*(1.0d0-dy)* 	dz+ 	   &
!        z5*       dx*(1.0D0-dy) * dz+ 		   &
!        z6*(1.0D0-dx) *       dy *dz+ 		   &
!        z7*       dx *       dy *     (1.0D0-dz) + &
!        z8*	  dx*       dy *       dz
	r1=dimag(blint)	!Add cut off
	r2=real(blint)
          if(r1.gt.-damp) r1=-damp    
	blint= dcmplx(r2,r1)
    endif
    case(2) !else if(star_interpolation) then
!	    else if(star_interpolation) then
!        Using star interpolation method
!      interpolation by introducing a factor exp(i kt X) in sigma
    Sapprox=zeroc
    do ic=1,Nc
        do ik=1,Nc
!	z8=1.0/(wn(n)+0.001-GcgSigm(ik,n))    
	z8=GcgSigm(ik,n)            
            Sapprox(ic)= Sapprox(ic)+FTCoefs_K_to_R(ik,ic)*z8
        enddo
 enddo
    blint=zeroc
        if(ndim==1) then
    do ic=1,Nc
!            z1=exp(+ii*(k*Rc(1,ic)))
!            blint=blint+z1*Sapprox(ic)
	z2=cos(k*Rc(1,ic))
      	z3=cos(k*Rc(1,ic))
        blint=blint+z2*real(Sapprox(ic))+ii*z3*dimag(Sapprox(ic))
	enddo
	r1=dimag(blint)
	r2=real(blint)
          if(r1.gt.-damp) r1=-damp    
	blint= dcmplx(r2,r1) 
        elseif(ndim==2) then
    do ic=1,Nc
!            z1=exp(+ii*(kx*Rc(1,ic)+ky*Rc(2,ic)))
!            blint=blint+z1*Sapprox(ic)
	z2=cos(kx*Rc(1,ic)+ky*Rc(2,ic))
      	z3=cos(kx*Rc(1,ic)+ky*Rc(2,ic))
        blint=blint+z2*real(Sapprox(ic))+ii*z3*dimag(Sapprox(ic))
	enddo
	r1=dimag(blint)
	r2=real(blint)
          if(r1.gt.-damp) r1=-damp    
	blint= dcmplx(r2,r1) 
        elseif(ndim==3) then
    do ic=1,Nc
!            z1=exp(+ii*(kx*Rc(1,ic)+ky*Rc(2,ic)+kz*Rc(3,ic)))
	z2=cos(kx*Rc(1,ic)+ky*Rc(2,ic)+kz*Rc(3,ic))
      	z3=cos(kx*Rc(1,ic)+ky*Rc(2,ic)+kz*Rc(3,ic))
        blint=blint+z2*real(Sapprox(ic))+ii*z3*dimag(Sapprox(ic))
!            blint=blint+z1*Sapprox(ic)
	enddo
	r1=dimag(blint)
	r2=real(blint)
        if(r1.gt.-damp) r1=-damp   ! Add damping to the interpolation
	blint= dcmplx(r2,r1)    
        endif

    case(3) !else if(spline_interpolation) then
!	    else if(spline_interpolation) then
    if(ndim==1) then
        g1=kx/gvector(1,1)
    !         Interpolate for the real part of GcgSigm(kx,n)
        call splint_1D(g1_sp(-N_sp),Sig_sp_r(-N_sp,0,n), &
        Derivs_sp_r(-N_sp,0,n),2*N_sp,g1,r7)
    !         Interpolate for the imaginary part of GcgSigm(kx,n)
        call splint_1D(g1_sp(-N_sp),Sig_sp_i(-N_sp,0,n), &
        Derivs_sp_i(-N_sp,0,n),2*N_sp,g1,r8)
	
          if(r8.gt.-damp) r8=-damp         
        blint=dcmplx(r7,r8)
 !       blint=(1-icum)*blint + icum*(wn(n)-oner/blint)          ! use Sigma or cumulant, depending on icum

    elseif(ndim==2) then
    !         Use a bicubic spline.
    !         Find the location on the grid of principle translation
    !         vectors in K-space.
        g1=(kx*gvector(2,2)-ky*gvector(2,1))/(gvector(1,1)*gvector(2,2)-gvector(1,2)*gvector(2,1))
        g2=(kx*gvector(1,2)-ky*gvector(1,1))/(gvector(2,1)*gvector(1,2)-gvector(2,2)*gvector(1,1))
    !         Interpolate for the real part of GcgSigm(kx,ky,n)
        call splin2(g1_sp(-N_sp),g2_sp(-N_sp),Sig_sp_r(-N_sp,-N_sp,n), &
        Derivs_sp_r(-N_sp,-N_sp,n),2*N_sp+1,2*N_sp+1,g1,g2,r7)
    !         Interpolate for the imaginary part of GcgSigm(kx,ky,n)
        call splin2(g1_sp(-N_sp),g2_sp(-N_sp),Sig_sp_i(-N_sp,-N_sp,n), &
        Derivs_sp_i(-N_sp,-N_sp,n),2*N_sp+1,2*N_sp+1,g1,g2,r8)
        	

          if(r8.gt.-damp) r8=-damp         
        blint=dcmplx(r7,r8)
!        blint=(1-icum)*blint + icum*(wn(n)-oner/blint)          ! use Sigma or cumulant, depending on icum	         

    elseif(ndim==3) then

          delta_g  = gvector(1,1)*gvector(2,2)*gvector(3,3) + gvector(1,2)*gvector(2,3)*gvector(3,1)& 
                  + gvector(1,3)*gvector(2,1)*gvector(3,2)&
                  - gvector(1,1)*gvector(2,3)*gvector(3,2) - gvector(1,2)*gvector(2,1)*gvector(3,3)&
                  - gvector(1,3)*gvector(2,2)*gvector(3,1)
          delta_g1 = kx    *gvector(2,2)*gvector(3,3) + gvector(1,2)*gvector(2,3)*kz &    
                  + gvector(1,3)*ky    *gvector(3,2)&
                  - kx    *gvector(2,3)*gvector(3,2) - gvector(1,2)*ky    *gvector(3,3) &
                  - gvector(1,3)*gvector(2,2)*kz    
          delta_g2 = gvector(1,1)*ky    *gvector(3,3) + kx    *gvector(2,3)*gvector(3,1) &
                  + gvector(1,3)*gvector(2,1)*kz&    
                  - gvector(1,1)*gvector(2,3)*kz     - kx    *gvector(2,1)*gvector(3,3)&
                  - gvector(1,3)*ky    *gvector(3,1)
          delta_g3 = gvector(1,1)*gvector(2,2)*kz     + gvector(1,2)*ky    *gvector(3,1)&
                  + kx    *gvector(2,1)*gvector(3,2)&
                  - gvector(1,1)*ky    *gvector(3,2) - gvector(1,2)*gvector(2,1)*kz &    
                  - kx    *gvector(2,2)*gvector(3,1)
          g1 = delta_g1/delta_g
          g2 = delta_g2/delta_g
          g3 = delta_g3/delta_g

!c         Interpolate for the real part of Sigma of Green function
          call splint3(g1_sp(-N_sp),g2_sp(-N_sp),g3_sp(-N_sp),&
     &              Sig_sp_r3(-N_sp,-N_sp,-N_sp,n),&
     &              Derivs_sp_r3(-N_sp,-N_sp,-N_sp,n),&
     &              2*N_sp+1,2*N_sp+1,2*N_sp+1,&
     &              g1,g2,g3,r7)
!c         Interpolate for the imaginary part of Sigma of Green function
          call splint3(g1_sp(-N_sp),g2_sp(-N_sp),g3_sp(-N_sp),&
     &              Sig_sp_i3(-N_sp,-N_sp,-N_sp,n),&
     &              Derivs_sp_i3(-N_sp,-N_sp,-N_sp,n),&
     &              2*N_sp+1,2*N_sp+1,2*N_sp+1,&
     &              g1,g2,g3,r8)
     
          if(r8.gt.-damp) r8=-damp         
        blint=dcmplx(r7,r8)
!        r1=dimag(blint) ! For some K-direction, uncausality occurs
!        r2=real(blint)
!        if(r1>epsilon) r1=-epsilon
!        blint=dcmplx(r2,r1)
!          if(icum.gt.0) blint = onec/blint + acoef*ii
!          blint=(1-icum)*blint + icum*(wn(n)-oner/blint)          ! use Sigma or cumulant, depending on icum
    end if

    case default
    print*, "Please, choose between case 1-3. See module_global.for for details"
    stop
    end select
    return
    end function blint

    		    
    real*8 function Epsk(kx,ky,kz)
    use Global
    implicit none
    real(kind) :: kx, ky ,kz

    if(ndim==1) then
        Epsk=-0.5_kind*(cos(kx)) - &
        tprime*(cos(kx)-1.0_kind)
    elseif(ndim==2) then
        Epsk=-0.5_kind*(cos(kx)+cos(ky)) - &
        tprime*(cos(kx)*cos(ky)-1.0_kind)

    elseif(ndim==3) then
        Epsk=-0.5_kind*(cos(kx)+cos(ky)+cos(kz)) &
        - tprime*(cos(kx)*cos(ky) &
        +cos(ky)*cos(kz)+cos(kx)*cos(kz)-1.0_kind)
    end if
    		    
    end function Epsk
    		    
    SUBROUTINE gaussj(p,n,np,b,m,mp)
    use Global
    INTEGER :: m,mp,n,np,NMAX
    REAL(kind) p(np,np),b(np,mp)
    PARAMETER (NMAX=50)
! inear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a(1:n,1:n)
! s an input matrix stored in an array of physical dimensions np by np. b(1:n,1:m) is an input
! atrix containing the m right-hand side vectors, stored in an array of physical dimensions
! p by mp. On output, a(1:n,1:n) is replaced by its matrix inverse, and b(1:n,1:m) is
! eplaced by the corresponding set of solution vectors.
! arameter: NMAX is the largest anticipated value of n.
    INTEGER :: i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiva(NMAX)
    do  j=1,n
        ipiva(j)=0
    enddo
    do  i=1,n
        big=0.
        do j=1,n
            if(ipiva(j) /= 1)then
                do  k=1,n
                    if (ipiva(k) == 0) then
                        if (abs(a(j,k)) >= big)then
                            big=abs(p(j,k))
                            irow=j
                            icol=k
                        endif
                    endif
                enddo
            endif
        enddo
        ipiva(icol)=ipiva(icol)+1

        if (irow /= icol) then
            do  l=1,n
                dum=a(irow,l)
                p(irow,l)=a(icol,l)
                p(icol,l)=dum
            enddo
            do  l=1,m
                dum=b(irow,l)
                b(irow,l)=b(icol,l)
                b(icol,l)=dum
            enddo
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (p(icol,icol) == 0.) pause 'singular matrix in gaussj'
        pivinv=1./p(icol,icol)
        a(icol,icol)=1.
        do l=1,n
            p(icol,l)=p(icol,l)*pivinv
        enddo
        do  l=1,m
            b(icol,l)=b(icol,l)*pivinv
        enddo
        do  ll=1,n
            if(ll /= icol)then
                dum=a(ll,icol)
                a(ll,icol)=0.
                do  l=1,n
                    p(ll,l)=a(ll,l)-p(icol,l)*dum
                enddo
                do  l=1,m
                    b(ll,l)=b(ll,l)-b(icol,l)*dum
                enddo
            endif
        enddo
    enddo
    do  l=n,1,-1
        if(indxr(l) /= indxc(l))then
            do  k=1,n
                dum=p(k,indxr(l))
                p(k,indxr(l))=p(k,indxc(l))
                p(k,indxc(l))=dum
            enddo
        endif
    enddo
    return
    END subroutine gaussj
    		    

!       NUMERICAL RECIPES 2D SPLINE

    SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A)
    PARAMETER (NN=100)
    IMPLICIT real(8) (A - H, O - Z)
    real(8) X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN)
    DO J=1,M
        DO K=1,N
            YTMP(K)=YA(J,K)
        END DO
        CALL SPLINE2D(X2A,YTMP,N,1.d30,1.d30,Y2TMP)
        DO K=1,N
            Y2A(J,K)=Y2TMP(K)
        END DO
    END DO
    RETURN
    END SUBROUTINE SPLIE2


    SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)
    PARAMETER (NN=100)
    IMPLICIT real(8) (A - H, O - Z)
    real(8) X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN),YYTMP(NN)
    DO J=1,M
        DO K=1,N
            YTMP(K)=YA(J,K)
            Y2TMP(K)=Y2A(J,K)
        END DO
        CALL SPLINT(X2A,YTMP,Y2TMP,N,X2,YYTMP(J))
    END DO
    CALL SPLINE2D(X1A,YYTMP,M,1.d30,1.d30,Y2TMP)
    CALL SPLINT(X1A,YYTMP,Y2TMP,M,X1,Y)
    RETURN
    END SUBROUTINE SPLIN2

    SUBROUTINE SPLINE2D(X,Y,N,YP1,YPN,Y2)
    IMPLICIT real(8) (A - H, O - Z)
    PARAMETER (NMAX=100)
    real(8) X(N),Y(N),Y2(N),U(NMAX)
    IF (YP1 > .99E30) THEN
        Y2(1)=0.
        U(1)=0.
    ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
    ENDIF
    DO I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
        /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
    END DO
    IF (YPN > .99E30) THEN
        QN=0.
        UN=0.
    ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
    ENDIF
    Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
    DO K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
    END DO
    RETURN
    END subroutine spline2D


    SUBROUTINE spline1D(x,y,n,c)
    use Global
    INTEGER :: n
    REAL(kind) x(0:n),y(0:n), aa(0:n), M(0:n-1,0:n-1), gg(n),c(0:n)
    aa=y

    M(:,:)=0
    do i=0, n-2
        M(i,i)=2*(x(i+2)-x(i))
        M(i+1,i)=(x(i+2)-x(i+1))
        M(i,i+1)=(x(i+2)-x(i+1))
    enddo
    M(n-1,n-1)=2*(x(n)-x(n-1)+x(1)-x(0))
    M(0,n-1)=x(1)-x(0)
    M(n-1,0)=x(1)-x(0)

    do i=1, n-1
        gg(i)=3/(x(i+1)-x(i))*(aa(i+1)-aa(i))-3/(x(i)-x(i-1)) &
        *(aa(i)-aa(i-1))
    enddo
    gg(n)=3/(x(1)-x(0))*(aa(1)-aa(n))-3/(x(n)-x(n-1))*(aa(n)-aa(n-1))


    call gaussj(M,n,n,gg,1,1) !Solve the matrix eqn. M*aa=g
! igenvalues will be returned in g
    do i=1, n
        c(i)=gg(i)
    enddo
    c(0)=c(n)
    RETURN
    END SUBROUTINE spline1D

    subroutine spline_1_3D(GcgSigm)
    use Global
    implicit none
    integer :: filenum, istat, i,j,l,ic,n,ickx,icky,ickz,ikz,ik, &
    nkpoints, info, infot, ikx, iky, deg, itunnel,i_energy, &  
    ISPLINE,i1,j1,i2,j2,l1,l2,ic1,ic2,ic3,ic4, &
    ic5,ic6,ic7,ic8
    real(kind) :: r1, r2, r3, k,kx, ky, kz,fill, Epsk, kx1, ky1,kx2,k1,k2, &
    ky2,offset,vx, vyy,kz1,kz2,Energy,dx,dy,dz,Kx3,Ky3,Kz3,r7,r8,g1,g2,g3
    real(8) :: delta_g, delta_g1, delta_g2, delta_g3
    		    
    complex(kind) :: c1, c2,z1,z2,z3,z4,z5,z6,z7,z8,Sapprox(Nc),GcgSigm(1:Nc,-nwn:nwn),ztemp
	
	if(ndim==1) then

    !         Form the 1D interpolation arrays needed to interpolate
    !         the self energy with a bicubic spline.
    !         1. form the grid in the space of g1
        do i=-N_sp,N_sp
            g1_sp(i)=i
        end do

        do n=-nf,nf
            do i=-N_sp,N_sp
                ick=ict(i,0,0)
                Sig_sp_r(i,0,n)=real(GcgSigm(ick,n))
                Sig_sp_i(i,0,n)=aimag(GcgSigm(ick,n))
            end do
            call spline1D(g1_sp(-N_sp), &
            Sig_sp_r(-N_sp,0,n),2*N_sp, &
            Derivs_sp_r(-N_sp,0,n))
            call spline1D(g1_sp(-N_sp), &
            Sig_sp_i(-N_sp,0,n),2*N_sp, &
            Derivs_sp_i(-N_sp,0,n))

        end do
	elseif(ndim==2) then
    !         Form the 2D interpolation arrays needed to interpolate
    !         the self energy with a bicubic spline.
    !         1. form the grid in the space of g1 and g2
        do i=-N_sp,N_sp
            g1_sp(i)=i
            g2_sp(i)=i
        end do

        do n=-nf,nf
            do i=-N_sp,N_sp
                do j=-N_sp,N_sp
                    ick=ict(i,j,0)
                    Sig_sp_r(i,j,n)=real(GcgSigm(ick,n))
                    Sig_sp_i(i,j,n)=aimag(GcgSigm(ick,n))
                end do
            end do
            call splie2(g1_sp(-N_sp),g2_sp(-N_sp), &
            Sig_sp_r(-N_sp,-N_sp,n),2*N_sp+1,2*N_sp+1, &
            Derivs_sp_r(-N_sp,-N_sp,n))
            call splie2(g1_sp(-N_sp),g2_sp(-N_sp), &
            Sig_sp_i(-N_sp,-N_sp,n),2*N_sp+1,2*N_sp+1, &
            Derivs_sp_i(-N_sp,-N_sp,n))
        end do
	elseif(ndim==3) then
!       interpolate \GcgSigm(K,\omega) to a 3 digital dense frequency grid
!  set up frequency grid, wh: 3 digital homogenous grid.  w3i: 3 digital round off inhomogenous grid



      do i=-N_sp,N_sp
         g1_sp(i)=i
         g2_sp(i)=i
         g3_sp(i)=i
      end do
      do n=-nf,nf
         do i=-N_sp,N_sp
         do j=-N_sp,N_sp
         do l=-N_sp,N_sp
            ick=ict(i,j,l)
!	if(icum>0) then
!           ztemp = (1.0D0-icum)*GcgSigm(ick,n)&
!                + icum/(GcgSigm(ick,n)-acoef*ii)
!	else
	    ztemp=GcgSigm(ick,n)
!	endif
            Sig_sp_r3(i,j,l,n)=real(ztemp)
            Sig_sp_i3(i,j,l,n)=dimag(ztemp)
         end do
         end do
         end do

         call spline3(g1_sp(-N_sp),g2_sp(-N_sp),g3_sp(-N_sp),&
                    Sig_sp_r3(-N_sp,-N_sp,-N_sp,n),&
                    2*N_sp+1,2*N_sp+1,2*N_sp+1,&
                    Derivs_sp_r3(-N_sp,-N_sp,-N_sp,n))
         call spline3(g1_sp(-N_sp),g2_sp(-N_sp),g3_sp(-N_sp),&
                    Sig_sp_i3(-N_sp,-N_sp,-N_sp,n),&
                    2*N_sp+1,2*N_sp+1,2*N_sp, &
                    Derivs_sp_i3(-N_sp,-N_sp,-N_sp,n))
      end do
	else
	end if

!        do ick=1,nkpoints
!            kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
!            ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)
!            kz=kz1 + real(ick-1)*(kz2-kz1)/max(nkpoints-1,1)
!      do n=-nf,nf
!         c2_sigma(ick,n) = blint_spline(kx,ky,kx,n)
!      end do
!      end do
	return
	end subroutine spline_1_3D



    SUBROUTINE SPLINT_1D(XA,YA,Y2A,N,X,Y)
    IMPLICIT real(8) (A - H, O - Z)
    real(8) XA(N),YA(N),Y2A(N)
    KLO=1
    KHI=N
    1 IF (KHI-KLO > 1) THEN
        K=(KHI+KLO)/2
        IF(XA(K) > X)THEN
            KHI=K
        ELSE
            KLO=K
        ENDIF
        GOTO 1
    ENDIF
    H=XA(KHI)-XA(KLO)
    IF (H == 0.) PAUSE 'Bad XA input in SPLINT_1D.'
    if(ndim==2) then
        A=(XA(KHI)-X)/H
        B=(X-XA(KLO))/H
        Y=A*YA(KLO)+B*YA(KHI)+ &
        ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
    elseif(ndim==1) then
        p=KLO
        H=XA(p+1)-XA(p)
        B=1/H*(YA(p+1)-YA(p))-h/3*(Y2A(p+1)+2*Y2A(p))
        D=1/(3*H)*(Y2A(p+1)-Y2A(p))
        Y=YA(p)+B*(X-XA(p))+Y2A(p)*(X-XA(p))**2+D*(X-XA(p))**3
    end if
    RETURN
    END SUBROUTINE SPLINT_1D

      ! calculate the interpolated value using the input and f''
      subroutine splint(xa,ya,y2a,n,x,y) 
      implicit real(8) (a - h, o - z) 
      real(8) xa(n),ya(n),y2a(n) 
      klo=1 
      khi=n 
1     if (khi-klo.gt.1) then 
        k=(khi+klo)/2 
        if(xa(k).gt.x)then 
          khi=k 
        else 
          klo=k 
        endif 
      goto 1 
      endif 
      h=xa(khi)-xa(klo) 
      if (h.eq.0.) pause 'bad xa input in splint.' 
      a=(xa(khi)-x)/h 
      b=(x-xa(klo))/h 
      y=a*ya(klo)+b*ya(khi)+ &
           ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6. 
      return 
      end subroutine splint

      subroutine spline3(x1a,x2a,x3a,ya,m,n,nz,y2a) 
      parameter (nn=100) 
      implicit real(8) (a - h, o - z) 
      real(8) x1a(m),x2a(n),x3a(nz)
      real(8) ya(m,n,nz),y2a(m,n,nz)
      do 13 jz=1,nz 
        call spline2(x1a,x2a,ya(1,1,jz),m,n,y2a(1,1,jz))
13    continue 
      return 
      end subroutine spline3


      subroutine splint3(x1a,x2a,x3a,ya,y2a,m,n,nz,x1,x2,x3,y) 
      parameter (nn=100) 
      implicit real(8) (a - h, o - z) 
      real(8) x1a(m),x2a(n),x3a(nz)
      real(8) ya(m,n,nz),y2a(m,n,nz)
      real(8) ytmp(nn),y2tmp(nn),yytmp(nn)
      do 12 jz=1,nz 
        call splint2(x1a,x2a,ya(1,1,jz),y2a(1,1,jz),m,n,x1,x2,yytmp(jz))
12    continue 
      call SPLINE2D(x3a,yytmp,nz,1.d30,1.d30,y2tmp) 
      call splint(x3a,yytmp,y2tmp,nz,x3,y) 
      return 
      end subroutine splint3

!c       NUMERICAL RECIPES 2D SPLINE

      !  2d f'' calculation using 1D spline f'' subroutine
      subroutine spline2(x1a,x2a,ya,m,n,y2a) 
      parameter (nn=100) 
      implicit real(8) (a - h, o - z) 
      real(8) x1a(m),x2a(n),ya(m,n),y2a(m,n),ytmp(nn), y2tmp(nn) 
      do 13 j=1,m 
        do 11 k=1,n 
          ytmp(k)=ya(j,k) 
11      continue 
        call SPLINE2D(x2a,ytmp,n,1.d30,1.d30,y2tmp) 
        do 12 k=1,n 
          y2a(j,k)=y2tmp(k) 
12      continue 
13    continue 
      return 
      end subroutine spline2

      !  calculate the second-order derivative
      subroutine spline(x,y,n,yp1,ypn,y2) 
      implicit real(8) (a - h, o - z) 
      parameter (nmax=100) 
      real(8) x(n),y(n),y2(n),u(nmax) 
      if (yp1.gt..99e30) then 
        y2(1)=0. 
        u(1)=0. 
      else 
        y2(1)=-0.5 
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) 
      endif 
      do 11 i=2,n-1 
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1)) 
        p=sig*y2(i-1)+2. 
        y2(i)=(sig-1.)/p 
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))& 
           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p 
11    continue 
      if (ypn.gt..99e30) then 
        qn=0. 
        un=0. 
      else 
        qn=0.5 
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1))) 
      endif 
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.) 
      do 12 k=n-1,1,-1 
        y2(k)=y2(k)*y2(k+1)+u(k) 
12    continue 
      return 
      end subroutine spline

      !  2d interpolation evaluation using 1d spline
      subroutine splint2(x1a,x2a,ya,y2a,m,n,x1,x2,y) 
      parameter (nn=100) 
      implicit real(8) (a - h, o - z) 
      real(8) x1a(m),x2a(n),ya(m,n),y2a(m,n),ytmp(nn),y2tmp(nn), yytmp(nn) 
      do 12 j=1,m 
        do 11 k=1,n 
          ytmp(k)=ya(j,k) 
          y2tmp(k)=y2a(j,k) 
11      continue 
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j)) 
12    continue 
      call SPLINE2D(x1a,yytmp,m,1.d30,1.d30,y2tmp) 
      call splint(x1a,yytmp,y2tmp,m,x1,y) 
      return 
      end subroutine splint2

