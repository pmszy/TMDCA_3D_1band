!***************************************************************************
subroutine analysis_twoparticle
!***************************************************************************
!Compute vertex function and two particle quantities
!***************************************************************************
use Global
!***************************************************************************
implicit none
integer :: i,n,m,p,ic,jc,lc,iw,icx,icy,icz,info
complex*16,allocatable :: chi_inv(:,:),chi_inv0(:),Fl(:,:),tmp_fullQ(:,:),&
Glatt(:),Qcl(:,:,:),Qcl0(:,:),Qcl0_LR(:,:),Vertexf(:,:,:),chicl0(:,:),&
chicharge_meas0(:,:),Gsc_tmp(:,:)
complex*16 :: tmp_z1
real*8 :: sum_v,Kx,Ky,Kz,epsilonk,vkx,t_hop,deltaw
real*8 :: r1,r2,r3,r4,conduc,analysis_time
integer :: i1,i2,i3,i4

t_hop=0.25d0

!***************************************************************************
	time1= secnds(0.0)
        call profstart('analysis_twoparticle')
!***************************************************************************	



!***************************************************************************	
!       Allocate arrays

        allocate(chi_inv(1:Nc,1:Nc),stat=info)
        allocate(chi_inv0(1:Nc),stat=info)
        allocate(Fl(1:Nc,1:Nc),stat=info)
        allocate(tmp_fullQ(1:Nc,1:Nc),stat=info)
        allocate(Glatt(-nwn:nwn),stat=info)
        allocate(Qcl(1:Nc,1:Nc,-nwn:nwn),stat=info)
        allocate(Qcl0(1:Nc,-nwn:nwn),stat=info)
        allocate(Qcl0_LR(1:Nc,-nwn:nwn),stat=info)
        allocate(Vertexf(1:Nc,1:Nc,-nwn:nwn),stat=info)
        allocate(chicl0(1:Nc,-nwn:nwn),stat=info)
        allocate(chicharge_meas0(1:Nc,-nwn:nwn),stat=info)
        allocate(Gsc_tmp(1:Nc,-nwn:nwn),stat=info)

!***************************************************************************	
!  	  The particle-hole response needed to calculate
!  	  the conductivity was computed in measure.  We will now compute
!         the bare bubble, chi0, and extract the vertex function from the
!         equation below.
!
!  	  K,w	    K',w	K,w	     K,w	     P,w      K',w
!  	    ----<----	     ----<----      ----<----|-------|----<----
!  	      / / /				     |       |  / / /
!  	      /Chi/	 =   delta	 +	     | Gamma |  /Chi/
!  	      / / /		  K,K'  	     |       |  / / /
!  	    ---->----	     ---->----      ---->----|-------|---->----
!  	  K,w+w'    K',w+w'	K,w+w'       K,w+w'    P,w+w'	      K',w+w'
!
!
!  	  Only K,K',P are matrix labels, all of the diagrams, at the
!  	  two-particle level, carry the same frequency labels, so for each
!  	  w and w', we solve the matrix equation for Gamma
!
!  	  Chi(K,K',w,w') = G(K,w)G(K,w+w') +
!  			   G(K,w)G(K,w+w')Gamma(K,P,w,w')Chi(P,K',w,w')
!
!
!
!
!  	  Chi0(K,K',w,w') = \delta     Gc(K,K,w) Gc(K,K,w+w') 
!                                 K,K'
!	
!         To calculate the DC limit (w'-->0) we will set w'=wn(1)

!       
!***************************************************************************************
!        Read-in the chicharge_meas that was computed in meas.f
!
!***************************************************************************************
!         Now calculate
!
!  	  Chi0(K,K,w,w') = Gc(K,K,w) Gc(K,K,w+w')
!
!              K      w     K'     ick   <-> K'
!         or   -------<-------     jck   <-> K
!              ------->-------     wn(1) <-> w'
!              K     w+w'   K'
!
!         Time-order the cluster average Green's function 
!         (named GcKf_av in put.f; This is also read in
!         the spectra code)

          do ic=1,Nc
            do n=-nwn,-1
              Gsc_tmp(ic,n)=GcKf_av(ic,n)
            end do
            do n=0,nwn
              Gsc_tmp(ic,n)=dconjg(GcKf_av(ic,n))
            end do
          end do
              

          chicharge_meas0=dcmplx(0.d0,0.d0)
          do ic=1,Nc
            do n=-nwn,nwn-1
	    chicharge_meas0(ic,n) = Gsc_tmp(ic,n)*Gsc_tmp(ic,n+1)
	    end do
	  end do

!***************************************************************************************
!       Now invert the chi(K,K',w,w') and chi0(K,K',w,w') for each (w); (w'=w+nu)
!       and nu is the smallest non-zero frequency, wn(1)
!       Then get the vertex function
!***************************************************************************************

       open(unit=162,file='chicharge.dat',status='unknown')
        do iw=-nwn,nwn
         do ic=1,Nc
          do jc=1,Nc
            read(162,*) i1,i2,i3,r1,r2
            chicharge_meas(ic,jc,iw)=dcmplx(r1,r2)
          end do
         end do
        end do
        close(162)


        do iw=-nwn,nwn-1

          do ic=1,Nc
            chi_inv0(ic)=1.d0/chicharge_meas0(ic,iw)
            do jc=1,Nc
              chi_inv(ic,jc)=chicharge_meas(ic,jc,iw)
            end do
          end do

          call Invert_Nc(chi_inv)

          do ic=1,Nc
            do jc=1,Nc
              if(ic.eq.jc) then
                vertexf(ic,jc,iw)=chi_inv0(ic) - chi_inv(ic,jc)
              else
                vertexf(ic,jc,iw)= - chi_inv(ic,jc)
              end if
            end do
          end do


        end do


!***************************************************************************************
!       Now that we have the vertex function, we can go ahead and
!       compute the conductivity.
!       
!            
!                             _________________________      
!       Chi_cl_CG(K,K',w,w') = V_kx Chi_cl(k,k',w,w') V_k'x  
!
!                             __________________________
!       Chi_cl0_CG(K,K',w,w')= V_kx Chi_cl0(k,k',w,w') V_k'x 
!
!                                  _____________________            
!       Chi_cl0_CGL(K,K',w,w') =   V_kx Chi_cl0(k,k',w,w')
!
!                                 _____________________
!       Chi_cl0_CGR(K',K,w,w') =  Chi_cl0(k',k,w,w') V_kx
!
!        
!       Chi_cl_CG(K,K',w,w') = Chi_cl0_CG(K,K',w,w') 
!
!        _______                  +
!        \
!         \    
!         /      Chi_cl0_CGL(K,K',w,w') Fl(K',K",w,w') Chi_cl0_CGR(K",K,w,w') 
!        /
!        _______
!         K',K"
!       
!        Fl(K,K',w,w')=   Vertexf * (1    -    Vertexf * Chi_cl0_CG)^{-1} 

!       Note that the above equation has matrix multiplications and inversions 
!
!
!**************************************************************************************
!**************************************************************************************
!       The lattice Green's functions will be obtained, and
!       time-ordered and then used to calculate the
!       Chi_cl0(k,k',w,w'). Then the three Coarse Grained quantities,
!       Chi_cl0_CG, Chi_cl0_CGL, Chi_cl0_CGR will be obtained.

!**************************************************************************************
!**************************************************************************************

!       Calculate the Coarse-grained susceptiblities

        do ic=1,Nc
          do iw=-nwn,nwn
            Qcl0(ic,iw)=zeroc
            Qcl0_LR(ic,iw)=zeroc
            chicl0(ic,iw)=zeroc
          end do
        end do

        write(6,*) 'nwn_c=',nwn_c

        do ic=1,Nc

          ! Construct the localized part of the time-ordered lattice Green's function

          do iw=-nwn,-nwn_c-1
             Glatt(iw)=1.d0/(wn(iw)-Epsbar(ic)-sigma_temp(ic,iw))
          end do
          do iw=nwn_c+1,nwn
             Glatt(iw)=dconjg(1.d0/(wn(iw)-Epsbar(ic)-sigma_temp(ic,iw)))
          end do

          ! The extended part will be done inside the kt loop

          write(6,*) 'CG beginning.'


          r1=0.5d0/dfloat(nover)
          sum_v=0.d0
          do icx=-nover+1,nover
          do icy=-nover+1,nover
          do icz=-nover+1,nover
            Kx=Kc(1,ic)+ r1*((dfloat(icx)-0.5d0)*gvector(1,1)+&
                            (dfloat(icy)-0.5d0)*gvector(2,1)+&
                            (dfloat(icz)-0.5d0)*gvector(3,1))
            Ky=Kc(2,ic)+ r1*((dfloat(icx)-0.5d0)*gvector(1,2)+&
                            (dfloat(icy)-0.5d0)*gvector(2,2)+&
                            (dfloat(icz)-0.5d0)*gvector(3,2))
            Kz=Kc(3,ic)+ r1*((dfloat(icx)-0.5d0)*gvector(1,3)+&
                            (dfloat(icy)-0.5d0)*gvector(2,3)+&
                            (dfloat(icz)-0.5d0)*gvector(3,3))
            
            epsilonk=-2.d0*t_hop*(dcos(Kx)+dcos(Ky)+dcos(Kz))
            vkx=2.d0*t_hop*(dsin(Kx)+dsin(Ky)+dsin(Kz))

            sum_v=sum_v+vkx


            ! Extended part of the lattice Green's function

            do iw=-nwn_c,nwn_c      ! WITHIN MOBILITY EDGE
              tmp_z1=sigma_temp(ic,iw)
              Glatt(iw)=1.d0/(wn(iw)-epsilonk-tmp_z1)
              if(iw.ge.0) Glatt(iw)=dconjg(Glatt(iw))
            end do


            ! Bare 2-p quantities within the mobility edge

            do iw=-nwn_c,nwn_c
              Qcl0(ic,iw)=Qcl0(ic,iw)+vkx**2*Glatt(iw)*Glatt(iw+1)
              Qcl0_LR(ic,iw)=Qcl0_LR(ic,iw)+vkx*Glatt(iw)*Glatt(iw+1)
              chicl0(ic,iw)=chicl0(ic,iw)+Glatt(iw)*Glatt(iw+1)
            end do


          end do  ! icx
          end do  ! icy
          end do  ! icz

          ! Normalize the Coarse-grained quantities

          do iw=-nwn_c,nwn_c
            Qcl0(ic,iw)=Qcl0(ic,iw)/dfloat((2*nover)**3)
            Qcl0_LR(ic,iw)=Qcl0_LR(ic,iw)/dfloat((2*nover)**3)
            chicl0(ic,iw)=chicl0(ic,iw)/dfloat((2*nover)**3)
            sum_v=sum_v/dfloat((2*nover)**3)
          end do
          write(6,*) 'sum_v=',sum_v

          ! Bare 2-p quantities outside the mobility edge
          do iw=-nwn,-nwn_c-1
            Qcl0(ic,iw)=(sum_v**2)*Glatt(iw)*Glatt(iw+1)
            Qcl0_LR(ic,iw)=sum_v*Glatt(iw)*Glatt(iw+1)
            chicl0(ic,iw)=Glatt(iw)*Glatt(iw+1)
          end do
          do iw=nwn_c+1,nwn-1
            Qcl0(ic,iw)=(sum_v**2)*Glatt(iw)*Glatt(iw+1)
            Qcl0_LR(ic,iw)=sum_v*Glatt(iw)*Glatt(iw+1)
            chicl0(ic,iw)=Glatt(iw)*Glatt(iw+1)
          end do

        end do  ! ic loop
        


        do iw=-nwn,nwn-1
          do ic=1,Nc
            do jc=1,Nc
              
              Qcl(ic,jc,iw)=zeroc       ! Initializing the full charge
                                        !susceptiblity

              ! Find Fl(K,K') for each w
              Fl(ic,jc)=zeroc
              if(ic.eq.jc) Fl(ic,jc)=dcmplx(1.d0,0.d0)

              Fl(ic,jc)=Fl(ic,jc)-chicl0(ic,iw)*vertexf(ic,jc,iw)
            end do      ! ic
          end do        ! jc

          call Invert_Nc(Fl)

          do ic=1,Nc
            do jc=1,Nc

              ! Find Vertexf*Fl
              tmp_fullQ(ic,jc)=zeroc
              do lc=1,Nc
                tmp_fullQ(ic,jc)=tmp_fullQ(ic,jc)+Vertexf(ic,lc,iw)*Fl(lc,jc)
              end do    ! lc

              tmp_fullQ(ic,jc)=tmp_fullQ(ic,jc)/dfloat(Nc)      ! Normalization

              Qcl(ic,jc,iw)=Qcl0_LR(ic,iw)*tmp_fullQ(ic,jc)*Qcl0_LR(jc,iw)

              if(ic.eq.jc) Qcl(ic,jc,iw)=Qcl0(ic,iw)+Qcl(ic,jc,iw)


            end do      ! ic
          end do        ! jc

        end do  ! iw loop

        ! Finally the Conductivity 
        tmp_z1=zeroc
        deltaw=wn(1)-wn(0)
        do iw=-nwn,nwn-1

          do ic=1,Nc
            do jc=1,Nc

              tmp_z1=tmp_z1+Qcl(ic,jc,iw)*deltaw
            end do      ! ic
          end do        ! jc

        end do  ! iw loop
        tmp_z1=tmp_z1/dfloat(Nc**2)     ! Normalization

        write(6,*) 'Conductivity =',Qcl(1,1,1)/(ii*wn(1))
        write(6,*) 'Imag part=',dimag(tmp_z1/(ii*wn(1)))
        conduc=dreal(tmp_z1/(ii*wn(1)))

        write(6,*) 'Conductivity=',conduc

        deallocate(chi_inv)
        deallocate(chi_inv0)
        deallocate(Fl)
        deallocate(tmp_fullQ)
        deallocate(Glatt)
        deallocate(Qcl)
        deallocate(Qcl0)
        deallocate(Qcl0_LR)
        deallocate(Vertexf)
        deallocate(chicl0)
        deallocate(chicharge_meas0)
        deallocate(Gsc_tmp)

        call profend('analysis_twoparticle')
 	time2= secnds(time1)
        analysis_time=analysis_time+time2
return
end
	
!*************************************************************************************

        Subroutine Invert_Nc(mat)
        use Global
        implicit none
        complex(kind) :: mat(Nc,Nc), CLAPW(3*Nc)
        integer :: i,n,m,j, ipiv1(Nc), LCLAPW,info
        
        LCLAPW=3*Nc
        
!         ZGETRF  -  compute an LU factorization of a general M-by-N
!         matrix A using partial pivoting with row interchanges.
!         ZGETRI - compute the inverse of a matrix using the LU
!         factorization computed by ZGETRF
          call zgetrf(Nc, Nc, mat, Nc, ipiv1, info )
          if(info.ne.0) then
            write(6,*) 'zgetrf error: info = ',info
            stop
          endif
          call zgetri(Nc,mat,Nc,ipiv1,CLAPW,LCLAPW,info)
          if(info.ne.0) then
            write(6,*) 'zgetri error: info = ',info
            stop
          endif
        
        return
        end subroutine Invert_Nc

!*************************************************************************************
