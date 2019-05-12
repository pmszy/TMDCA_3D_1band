	    subroutine spectra

!       to compile type make spectra (or make all) using the
!       associated Makefile.  This code was designed to read the
!       magority of its inputs from a qmc sigma file.
!	This code uses the trilinear and star interpolation method.
!	It is for 1,2,and 3 D

	    use Global
	    implicit none
	    integer :: filenum, istat, i, j, ic,jc,n,ickx,icky,ickz,ikz, 
     &	    nkpoints, info, infot, ikx, iky, deg, itunnel,i_energy,m,p,
     &	    ISPLINE,i1,j1,i2,j2,l1,l2,ic1,ic2,ic3,ic4, 
     &	    ic5,ic6,ic7,ic8
	    real(kind) :: r1, r2, r3, kx, ky, kz,fill, Epsk, kx1, ky1,kx2,k1,k,k2, 
     &	    ky2,offset, res, hall, vx, vyy, dfermi,kz1,kz2,Energy, 
     &	    dx,dy,dz,Kx3,Ky3,Kz3,r7,r8
	    complex(kind) :: c1, c2, blint,z1,z2,z3,z4,z5,z6,z7,z8,c_lat
	complex(kind), allocatable :: c1_sigma(:,:),c_Gl1(:,:),c_Gl2(:,:),c_Gext(:,:),
     &			c2_sigma(:,:)
		    
	    character(72) :: linemc,fmt1,string1
	    CHARACTER F*128, FN*128, FILENAME*128
	    external blint


	!       Read in some parameters
	    N_sp=sqrt(float(Nc))
		nf = nwn
	    !icum = 0
!	    Energy= 3.50D0
	    nkpoints = 1000! This is the maximum. The input nkpoints overides this

	!       Make sure that nf is le the limit set above
	    if(nf > 1000) then
		write(6,*) 'nf in the data file is too large. I have adjusted it'
		nf = 1000
	    end if

!       Now do the interpolation, and calculate the spectrum at
!       a sequence of locations.
	    if(ndim==1) then
		write(6,"('Enter k1,k2,offset,nkpoints,w_c:  ',$)")
		read(5,*) k1,k2,offset,nkpoints,w_c
	    elseif(ndim==2) then
		write(6,"('Enter kx1,ky1,kx2,ky2,offset,nkpoints,w_c:  ',$)")
		read(5,*) kx1,ky1,kx2,ky2,offset,nkpoints,w_c
	    else if(ndim==3) then
		write(6,"('Enter kx1,ky1,kz1,kx2,ky2,kz2,offset,nkpoints,w_c:  ',$)")
		read(5,*) kx1,ky1,kz1,kx2,ky2,kz2,offset,nkpoints,w_c
	    end if

		    
!	Now that we have inputed all the data, allocate the internal arrays as to be the same dimension as input
	allocate(c1_sigma(1:nkpoints,-nf:nf),c_Gl1(1:nkpoints,-nf:nf),c_Gl2(1:nkpoints,-nf:nf),c_Gext(1:nkpoints,-nf:nf),
     &			c2_sigma(1:nkpoints,-nf:nf))

!	Fin the W_c on the nf grid
	do n=0,nwn
	if((wn(n)-w_c)>0) exit 
	enddo
	nwn_c = n

! You have to know w_c aprior.

!###################################################
!	Prepare the self-energy
!###################################################

	sigma_temp(:,:)=sigma_loc(:,:)  

!###################################################
!	Prepare the G_local
!###################################################
	G_loc1=zeroc
	G_loc2=zeroc
	do n =-nwn,-nwn_c
	do ic=1,Nc
	G_loc1(ic,n) = 1.0D0/(wn(n)-Epsbar(ic)-sigma_temp(ic,n)) 
	enddo
	enddo


	do n = nwn_c+1,nwn
	do ic=1,Nc
	G_loc2(ic,n) = 1.0D0/(wn(n)-Epsbar(ic)-sigma_temp(ic,n)) 
	enddo
	enddo
		    

	    open(unit=21,file='Akw_spectra.dat',action='write',status='unknown')
	    open(unit=22,file='False_color.dat',action='write',status='unknown')
	    write(6,*) '******************************************************'
	    write(6,*) '*** A(k,w) is in Akw_spectra.dat             *********'
            write(21,"('#       w      A(k,w)=-Im G(k,w)/pi ')")
            write(22,"('#      k      w        A(k,w)=-Im G(k,w)/pi ')")


!	Plot the spectra using xmgrace
          open (unit=8,status='unknown',file='Awk_spectra.agr')
         write(8,"('@    default font 0')") 
         write(8,90) 
 90      format('@    yaxis  label "A(K,\xw\N)"')
         write(8,91) 
 91      format('@    xaxis  label "\xw\N"')

	    if(ndim==1) then
!	Interpolate the self-energy first on all frequencies. Also, interpolate G_loc
	 do ick=1,nkpoints
		      k=k1 + real(ick-1)*(k2-k1)/max(nkpoints-1,1) !To ensure it runs even for 1 k-point
		       do i=-nf,nf
		           c2_sigma(ick,i)=blint(sigma_temp,k,0,0,i)
		      enddo

!	Interpolate the -ve side of the G_loc
		do n=-nwn,-nwn_c
		  c_Gl1(ick,i)=blint(G_loc1,k,0,0,i)
	       enddo	

!	Interpolate the +ve side of the G_loc
		do n=nwn_c+1,nwn
		  c_Gl2(ick,i)=blint(G_loc2,k,0,0,i)
	       enddo		

!	Restrict G_ext to the -nwn_c and nwn_c
		       do i=-nwn_c,nwn_c
		            c1_sigma(ick,i)=1.d0/(wn(i)-Epsk(k,0,0)-c2_sigma(ick,i))
		  enddo
!	Next, combine G_ext and G_locs to get G_lattice
		       do i=-nf,nf
			c_lat = c1_sigma(ick,i)+c_Gl1(ick,i)+c_Gl2(ick,i)
           		     write (8,"(1x,f15.6,', ',f15.6)") wn(i), 
     &		            -dimag(c_lat)/pi+(ick-1)*offset
		            write(21,"(5(f15.6,1x))") wn(i), 
     &		            -dimag(c_lat)/pi+(ick-1)*offset
		            write(22,"(i4,1x,2(f15.6,1x))") ick,wn(i),-dimag(c_lat)/pi
		end do

		        write(8,*) ' & '
		        write(21,*) ' & '
		        write(22,*) '  '
		    end do

	    elseif(ndim==2) then
!	Interpolate the self-energy first on all frequencies. Also, interpolate G_loc
		    do ick=1,nkpoints
		        kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
		        ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)
		       do i=-nf,nf
		            c2_sigma(ick,i)=blint(sigma_temp,kx,ky,0,i)
		   enddo

!	Interpolate the -ve side of the G_loc
		do n=-nwn,-nwn_c
		  c_Gl1(ick,i)=blint(G_loc1,kx,ky,0,i)
	       enddo	

!	Interpolate the +ve side of the G_loc
		do n=nwn_c+1,nwn
		  c_Gl2(ick,i)=blint(G_loc2,kx,ky,0,i)
	       enddo		

!	Restrict G_ext to the -nwn_c and nwn_c
		       do i=-nwn_c,nwn_c
		            c1_sigma(ick,i)=1.d0/(wn(i)-Epsk(kx,ky,0)-c2_sigma(ick,i))
		  enddo
!	Next, combine G_ext and G_locs to get G_lattice
		       do i=-nf,nf
			c_lat = c1_sigma(ick,i)+c_Gl1(ick,i)+c_Gl2(ick,i)
           		     write (8,"(1x,f15.6,', ',f15.6)") wn(i), 
     &		            -dimag(c_lat)/pi+(ick-1)*offset
		            write(21,"(5(f15.6,1x))") wn(i), 
     &		            -dimag(c_lat)/pi+(ick-1)*offset
		            write(22,"(i4,1x,2(f15.6,1x))") ick,wn(i),-dimag(c_lat)/pi
		end do

		        write(8,*) ' & '
		        write(21,*) ' & '
		        write(22,*) '  '
		    end do
	    elseif(ndim==3) then
!	Interpolate the self-energy first on all frequencies.
	do ick=1,nkpoints 
		        kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
		        ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)
		        kz=kz1 + real(ick-1)*(kz2-kz1)/max(nkpoints-1,1)
		       do i=-nf,nf
		            c2_sigma(ick,i)=blint(sigma_temp,kx,ky,kz,i)
		   enddo

!	Interpolate the -ve side of the G_loc
		do n=-nwn,-nwn_c
		  c_Gl1(ick,i)=blint(G_loc1,kx,ky,kz,i)
	       enddo	

!	Interpolate the +ve side of the G_loc
		do n=nwn_c+1,nwn
		  c_Gl2(ick,i)=blint(G_loc2,kx,ky,kz,i)
	       enddo		

!	Restrict G_ext to the -nwn_c and nwn_c
		       do i=-nwn_c,nwn_c
		            c1_sigma(ick,i)=1.d0/(wn(i)-Epsk(kx,ky,kz)-c2_sigma(ick,i))
		  enddo
!	Next, combine G_ext and G_locs to get G_lattice
		       do i=-nf,nf
			c_lat = c1_sigma(ick,i)+c_Gl1(ick,i)+c_Gl2(ick,i)
           		     write (8,"(1x,f15.6,', ',f15.6)") wn(i), 
     &		            -dimag(c_lat)/pi+(ick-1)*offset
		            write(21,"(5(f15.6,1x))") wn(i), 
     &		            -dimag(c_lat)/pi+(ick-1)*offset
		            write(22,"(i4,1x,2(f15.6,1x))") ick,wn(i),-dimag(c_lat)/pi
		end do

		        write(8,*) ' & '
		        write(21,*) ' & '
		        write(22,*) '  '
		    end do


	    endif
			close(8)
			close(21)
			close(22)
		
            stop
	    end
		            

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
	    integer :: filenum, istat, i, j, ic,n,ickx,icky,ickz,ikz,ik, 
     &	    nkpoints, info, infot, ikx, iky, deg, itunnel,i_energy,
     &	    ISPLINE,i1,j1,i2,j2,l1,l2,ic1,ic2,ic3,ic4, 
     &	    ic5,ic6,ic7,ic8
	    real(kind) :: r1, r2, r3, kx, ky, kz,fill, Epsk, kx1, ky1,kx2,k,k1,k2, 
     &	    ky2,offset, res, hall, vx, vyy, dfermi,kz1,kz2,Energy, 
     &	    dx,dy,dz,Kx3,Ky3,Kz3,r7,r8,g1,g2
		    
	    complex(kind) :: c1, c2,z1,z2,z3,z4,z5,z6,z7,z8,Sapprox(Nc),GcgSigm(1:Nc,-nwn:nwn)
		    
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
		    dx=(k-Kx3)/gvector(1,1)
		!         interpolate.
		    blint=z1+(z2-z1)*dx

		elseif(ndim==2) then
!         Use a two-dimensional bilinear interpolation within the square


!         1*        *2




!         3*        *4

!         First find the coordinates of the cell which contains (kx,ky)
		    i1=int((kx*gvector(1,1)+ky*gvector(1,2))/(gvector(1,1)*gvector(1,1)+gvector(1,2)*gvector(1,2))+Nc) 
     &		    -Nc
		    i2=i1+1
		    j1=int((kx*gvector(2,1)+ky*gvector(2,2))/(gvector(2,1)*gvector(2,1)+gvector(2,2)*gvector(2,2))+Nc) 
     &		    -Nc
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
		    z2=GcgSigm(n,ic2)
		    z3=GcgSigm(ic3,n)
		    z4=GcgSigm(ic4,n)
		!         calculate the slopes
		    dx=((kx-Kx3)*gvector(1,1)+(ky-Ky3)*gvector(1,2)) 
     &		    /(gvector(1,1)*gvector(1,1)+gvector(1,2)*gvector(1,2))
		    dy=((kx-Kx3)*gvector(2,1)+(ky-Ky3)*gvector(2,2)) 
     &		    /(gvector(2,1)*gvector(2,1)+gvector(2,2)*gvector(2,2))
		!         interpolate.
		    blint=z3+(z4-z3)*dx+(z1-z3)*dy+(z2-z1-z4+z3)*dx*dy


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
		    i1=int((kx*gvector(1,1)+ky*gvector(1,2)+kz*gvector(1,3))/ 
     &		    (gvector(1,1)*gvector(1,1)+gvector(1,2)*gvector(1,2)+gvector(1,3)*gvector(1,3))+Nc)-Nc
		    i2=i1+1
		    j1=int((kx*gvector(2,1)+ky*gvector(2,2)+kz*gvector(2,3))/ 
     &		    (gvector(2,1)*gvector(2,1)+gvector(2,2)*gvector(2,2)+gvector(2,3)*gvector(2,3))+Nc)-Nc
		    j2=j1+1
		    l1=int((kx*gvector(3,1)+ky*gvector(3,2)+kz*gvector(3,3))/ 
     &		    (gvector(3,1)*gvector(3,1)+gvector(3,2)*gvector(3,2)+gvector(3,3)*gvector(3,3))+Nc)-Nc
		    l2=l1+1
		!         now get the K-point on the lower left.
		    Kx3=i1*gvector(1,1)+j1*gvector(2,1)+l1*gvector(3,1)
		    Ky3=i1*gvector(1,2)+j1*gvector(2,2)+l1*gvector(3,2)
		    Kz3=i1*gvector(1,3)+j1+gvector(2,3)+l1*gvector(3,3)
		!         index the four corners
		    ic1=ict(i1,j1,l1)
		    ic2=ict(i2,j1,l1)
		    ic3=ict(i2,j1,l2)
		    ic4=ict(i1,j1,l2)
		    ic5=ict(i1,j2,l1)
		    ic6=ict(i2,j2,l1)
		    ic7=ict(i2,j2,l2)
		    ic8=ict(i1,j2,l2)
		    z1=GcgSigm(ic1,n)
		    z2=GcgSigm(ic2,n)
		    z3=GcgSigm(ic3,n)
		    z4=GcgSigm(ic4,n)
		    z5=GcgSigm(ic5,n)
		    z6=GcgSigm(ic6,n)
		    z7=GcgSigm(ic7,n)
		    z8=GcgSigm(ic8,n)
		!         calculate the slopes
		    dx=((kx-Kx3)*gvector(1,1)+(ky-Ky3)*gvector(1,2)+(kz-Kz3)*gvector(1,3))/ 
     &		    (gvector(1,1)*gvector(1,1)+gvector(1,2)*gvector(1,2)+gvector(1,3)*gvector(1,3))
		    dy=((kx-Kx3)*gvector(2,1)+(ky-Ky3)*gvector(2,2)+(kz-Kz3)*gvector(2,3))/ 
     &		    (gvector(2,1)*gvector(2,1)+gvector(2,2)*gvector(2,2)+gvector(2,3)*gvector(2,3))
		    dz=((kx-Kx3)*gvector(3,1)+(ky-Ky3)*gvector(3,2)+(kz-Kz3)*gvector(3,3))/ 
     &		    (gvector(3,1)*gvector(3,1)+gvector(3,2)*gvector(3,2)+gvector(3,3)*gvector(3,3))
		!         interpolate.
		    blint=z1*(1.0d0-dx)*(1.0d0-dy)*(1.0d0-dz)+ 
     &		    z2*       dx *(1.0d0-dy)*(1.0d0-dz)+ 
     &		    z3*       dx *(1.0d0-dy)* dz       + 
     &		    z4*(1.0d0-dx)*(1.0d0-dy)* dz       + 
     &		    z5*(1.0d0-dx)*       dy *(1.0d0-dz)+ 
     &		    z6*       dx *       dy *(1.0d0-dz)+ 
     &		    z7*       dx *       dy *       dz + 
     &		    z8*(1.0d0-dx)*       dy *       dz
		endif
	    case(2) !else if(star_interpolation) then         
!	    else if(star_interpolation) then
	    !        Using star interpolation method
		Sapprox=zeroc
		do ic=1,Nc
		    do ik=1,Nc
		        Sapprox(ic)= Sapprox(ic) 
     &		        +FTCoefs_K_to_R(ik,ic)*GcgSigm(ik,n)
		    enddo
		enddo

		blint=(0.d0,0.d0)
		do ic=1,Nc
		    if(ndim==1) then
		        z1=exp(+ii*(k*Rc(1,ic)))
		        blint=blint+z1*Sapprox(ic)
		    elseif(ndim==2) then
		        z1=exp(+ii*(kx*Rc(1,ic)+ky*Rc(2,ic)))
		        blint=blint+z1*Sapprox(ic)
		    elseif(ndim==3) then
		        z1=exp(+ii*(kx*Rc(1,ic)+ky*Rc(2,ic)+kz*Rc(3,ic)))
		        blint=blint+z1*Sapprox(ic)
		    endif
		enddo
            case(3) !else if(spline_interpolation) then
!	    else if(spline_interpolation) then
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
		        call spline1D(g1_sp(-N_sp), 
     &		        Sig_sp_r(-N_sp,0,n),2*N_sp, 
     &		        Derivs_sp_r(-N_sp,0,n))
		        call spline1D(g1_sp(-N_sp), 
     &		        Sig_sp_i(-N_sp,0,n),2*N_sp, 
     &		        Derivs_sp_i(-N_sp,0,n))

		    end do


		    g1=k/gvector(1,1)
		!         Interpolate for the real part of GcgSigm(kx,n)
		    call splint(g1_sp(-N_sp),Sig_sp_r(-N_sp,0,n), 
     &		    Derivs_sp_r(-N_sp,0,n),2*N_sp,g1,r7)
		!         Interpolate for the imaginary part of GcgSigm(kx,n)
		    call splint(g1_sp(-N_sp),Sig_sp_i(-N_sp,0,n), 
     &		    Derivs_sp_i(-N_sp,0,n),2*N_sp,g1,r8)
		         
		    blint=cmplx(r7,r8)


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
		        call splie2(g1_sp(-N_sp),g2_sp(-N_sp), 
     &		        Sig_sp_r(-N_sp,-N_sp,n),2*N_sp+1,2*N_sp+1, 
     &		        Derivs_sp_r(-N_sp,-N_sp,n))
		        call splie2(g1_sp(-N_sp),g2_sp(-N_sp), 
     &		        Sig_sp_i(-N_sp,-N_sp,n),2*N_sp+1,2*N_sp+1, 
     &		        Derivs_sp_i(-N_sp,-N_sp,n))
		    end do

		!         Use a bicubic spline.
		!         Find the location on the grid of principle translation
		!         vectors in K-space.
		    g1=(kx*gvector(2,2)-ky*gvector(2,1))/(gvector(1,1)*gvector(2,2)-gvector(1,2)*gvector(2,1))
		    g2=(kx*gvector(1,2)-ky*gvector(1,1))/(gvector(2,1)*gvector(1,2)-gvector(2,2)*gvector(1,1))
		!         Interpolate for the real part of GcgSigm(kx,ky,n)
		    call splin2(g1_sp(-N_sp),g2_sp(-N_sp),Sig_sp_r(-N_sp,-N_sp,n), 
     &		    Derivs_sp_r(-N_sp,-N_sp,n),2*N_sp+1,2*N_sp+1,g1,g2,r7)
		!         Interpolate for the imaginary part of GcgSigm(kx,ky,n)
		    call splin2(g1_sp(-N_sp),g2_sp(-N_sp),Sig_sp_i(-N_sp,-N_sp,n), 
     &		    Derivs_sp_i(-N_sp,-N_sp,n),2*N_sp+1,2*N_sp+1,g1,g2,r8)
		         
		    blint=cmplx(r7,r8)

		elseif(ndim==3) then
		    stop
		    write(6,*) 'No Spline for 3D. Use star or trilinear interpolation '
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
        real(kind) :: kx, ky ,kz,k 

	    if(ndim==1) then
		Epsk=-0.5_kind*(cos(k)) - 
     &		tprime*(cos(k)-1.0_kind)
	    elseif(ndim==2) then
		Epsk=-0.5_kind*(cos(kx)+cos(ky)) - 
     &		tprime*(cos(kx)*cos(ky)-1.0_kind)

	    elseif(ndim==3) then
		Epsk=-0.5_kind*(cos(kx)+cos(ky)+cos(kz)) 
     &		- tprime*(cos(kx)*cos(ky) 
     &		+cos(ky)*cos(kz)+cos(kx)*cos(kz)-1.0_kind)
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
		U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) 
     &		/(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
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
		gg(i)=3/(x(i+1)-x(i))*(aa(i+1)-aa(i))-3/(x(i)-x(i-1)) 
     &		*(aa(i)-aa(i-1))
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



	    SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
	    IMPLICIT real(8) (A - H, O - Z)
	    real(8) XA(N),YA(N),Y2A(N)
	    KLO=1
	    KHI=N
1 	    IF (KHI-KLO > 1) THEN
		K=(KHI+KLO)/2
		IF(XA(K) > X)THEN
		    KHI=K
		ELSE
		    KLO=K
		ENDIF
		GOTO 1
	    ENDIF
	    H=XA(KHI)-XA(KLO)
	    IF (H == 0.) PAUSE 'Bad XA input.'
	    if(ndim==2) then
		A=(XA(KHI)-X)/H
		B=(X-XA(KLO))/H
		Y=A*YA(KLO)+B*YA(KHI)+ 
     &		((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
	    elseif(ndim==1) then
		p=KLO
		H=XA(p+1)-XA(p)
		B=1/H*(YA(p+1)-YA(p))-h/3*(Y2A(p+1)+2*Y2A(p))
		D=1/(3*H)*(Y2A(p+1)-Y2A(p))
		Y=YA(p)+B*(X-XA(p))+Y2A(p)*(X-XA(p))**2+D*(X-XA(p))**3
	    end if
	    RETURN
	    END SUBROUTINE SPLINT


