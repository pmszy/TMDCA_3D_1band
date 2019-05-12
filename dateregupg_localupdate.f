        subroutine update(ntimes)

 
c       Update updates the potentials by sweeping sequentially through
c       them and proposing changes.  If changes are made
c       it updates the greens functions appropriately.
 
c************************************************************
	use Global
c************************************************************
 
        integer i,j,ic,jc,ntimes,ntime,n
        real delV,Vt
	complex accept,b(Ncm),d(Ncm)
 
        do ntime=1,ntimes
          nslr=nslr+1
          do jc=1,Nc
c           Now update the jc'th potential V(jc)
c           if(rand(0).gt.0.5)then 	!g77
c           Propose a change
            if(isv.eq.1) then
              Vt = -V(jc)
	    else
	      Vt= 2.0*V0*(ran(iseed)-0.5)
	    end if
            if(ran(iseed).gt.0.5.or.isv.eq.-1)then 
c	      accept the change
	      delV = Vt-V(jc)
	      V(jc)= Vt
	      do n=-nwn,nwn  
                accept=delV/(1.0-delV*G(jc,jc,n))
cc		if (n .eq. 0) write(6,*) accept
cc		if (n .eq. 0) write(6,*) G(1,1,0)
                call upg(jc,n,accept,G(1,1,n))
cc		if (n .eq. 0) write(6,*) G(1,1,0)
              end do
            end if
          end do

c         Now we must check to see if it is time to reconfigure the
c         Gd matrix.  Reconfiguration is necessary to deal with roundoff error.
 
          if(nslr.ge.ntr)then
            call reg
            nslr=0
          end if
	end do

        return
        end

        subroutine reg
c       This subroutine reconfigures G given only the potential configuration
c       and the initial green's functions Gcsfsc (which corresponds to
c	all V=0).
 
c************************************************************
	use Global
c************************************************************
 
        complex accept,d(Ncm),b(Ncm)
        integer n,ic,jc,i,j
 
c       Gcsfsc is the greens function for all V=0, hence we will
c       sweep through the lattice and reconfigure Gcsfsc to coorespond
c       to the lattice configuration.
 
c       move Gcsfsc to G
        do n=-nwn,nwn
        do ic=1,Nc
        do jc=1,Nc
          G(ic,jc,n)=Gcsfsc(icrdiff(ic,jc),n)
	end do
	end do
	end do

c       sweep over all the time slices and sites
        do jc=1,Nc
	  do n=-nwn,nwn
            accept=V(jc)/(1.0-V(jc)*G(jc,jc,n))
            call upg(jc,n,accept,G(1,1,n))
	  end do
        end do
 
        return
        end

        subroutine upg(ic,n,accept,gp)
 
c************************************************************
	use Global
c************************************************************
 
        integer ic,i,j,n
        complex accept,d(Ncm),b(Ncm),gp(Ncm,Ncm)
 
c       subroutine upg updates gd after a change in one of the field
c       variables has been accepted.
 
c       FORM THE VECTORS
	do i=1,Nc
          b(i)=gp(i,ic)*accept
          d(i)=gp(ic,i)
        end do
 
c       FORM THE OUTER PRODUCT, AND THE NEW G
c	5 WAYS TO DO IT (see which is fastest on your machine)
        do i=1,Nc
        do j=1,Nc
          gp(j,i)=gp(j,i)+b(j)*d(i)
        end do
        end do
c
c	This is equivalent to the level 3  blas calls (usually best)
c
c        call cgeru(Nc,Nc,1.0,b,1,d,1,gp,Nc)	! 32 bit
c        call zgeru(Nc,Nc,(1.0,0.0),b,1,d,1,gp,Ncm)	! 64 bit (even with -r8)
c
c	and
c
c	call cgemm('N','T',Nc,Nc,1,1.0,b,Nc,d,Nc,1.0,gp,Nc)   !UNTESTED
c
c	or to the loop using level 1 blas call
c	
c	do i=1,Nc
c	  call zaxpy(Nc,d(i),b,1,gp(1,i),1)      ! UNTESTED
c	end do
c
c	and to the equivalent fortran90 call
c       do i=1,Nc                                ! UNTESTED
c 	 gp(:,i)=gp(:,i) +b(:)*d(i)
c       end do

c	Subroutines that may be needed
c	call dgetrf(Nblocka,Nblocka,Mu,Nblock,indxu,info)




        return
        end
