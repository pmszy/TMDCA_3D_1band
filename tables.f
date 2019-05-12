        subroutine tables                    
c***************************************************************************************
c       This subroutine makes the lookup tables used throughout the
c       code.
c     This subroutine makes the lookup tables used throughout the
c     full THREE DIMENSIONAL code.  It employs the definitions and 
c     cluster geometries developed by D.D. Betts and G.E. Stewart,
c     Canadian Journal of Physics, 75, 47-66.
c***************************************************************************************
c	Create a frequency grid for the Green functions. 
c***************************************************************************************
c***************************************************************************************
        use Global
        implicit none
c***************************************************************************************
      integer i,j,n,m,k,m1,n1,m2,n2,iequiv,j1
      real (kind) :: yy,del,deltat,tc
c***************************************************************************************

           
     	  do n=0,nwn
	if(isw) then
c	  write(6,*) 'tables-spectra-1'
     	    yy=float(n)/float(2*(nwn+1))
     	    wn(n)=delta*tan(pi*yy)
     	    wn(-n)=-wn(n)
     	    dwn(n)=pi*(delta**2+(wn(n))**2)/(delta*2*nwn)
     	    dwn(-n)=dwn(n)
c     	    if(wn(n).ge.wc) go to 401
c    	  end do
c 401      nwn=n-1
c     	  write(6,*) 'nwn=',nwn
	else !if(isw.eq.1) then
c	  write(6,*) 'tables-spectra-2'
c    	  do n=0,nwn
    	    wn(n)=delta*float(n)/float(nwn) ! Use large enough delta for homogeneous distribution
     	    wn(-n)=-wn(n)
     	    dwn(n)=delta/float(nwn)
     	    dwn(-n)=dwn(n)
	end if	
     	  end do
c       determine the index ned of ed
        do n=-nwn,nwn
          ned=n
          if(wn(n).gt.ed) exit
        end do

c	Frequency grid for the measurement of p(eta) (the return probability) (see cond-mat/0006431).
c	Typically, only a few low frequecies are needed
	do n=1,neta
	if(isw) then
	  weta(n)=wn(n)
	else
	  weta(n)=halfr*(n-1)/float(neta)+0.001
	end if
	end do  
c	write(6,*) 'tables-spectra-3'
	
c	To calculate the localization time. This is the time much greater
c	than hbar/W (W the bandwidth 2) are necessary.  Here 
c	tc is the longest time considered, nt the number of
c	time samples, and deltat determines the (lorentzian) spacing 
c	of the time samples 

	deltat=8.d0
	tc=100.d0
	do n=1,ntm
	if(isw) then
     	  yy=float(n)/float(2*ntm)
     	  t(n)=deltat*tan(pi*yy)
	else
	  t(n)=tmax*(n-1)/float(ntm)+0.001
	end if
     	  if(t(n).ge.tc) go to 402
	end do
 402	nt=n-1
c      write(6,*) 'tables-spectra-3'
	
c	G0 is a bare green function used for conditioning some 
c	integrals.  Since it involves an inverse, it is better
c	to calculate it once here and use the table later.

	del=0.02
        do n=-nwn,nwn
   	  G0w(n)=1.d0/(wn(n)+ii*del)
        end do
	do n=1,nt
	  G0t(n)=-ii*exp(-del*t(n))
	end do
c	Now make a table for the Fourier factors 
	do i=1,nt
	do n=-nwn,nwn
	  four(n,i)= exp(-ii*wn(n)*t(i))
	end do
	end do
		


        open(unit=53,file='p_rho.dat',status='unknown') ! This is the output for the bare partial DOS and Total DOS

c**************************************************************
c Decide whether to use Betts lattice or conventional lattice
c**************************************************************
!        write(6,*) 'tables-spectra-4',ibetts
	    select case(tables_type)
	    case(1) !Betts
        !if(ibetts.eq.0) then
          if (ndim.eq.3) then
            call tables_3D_Betts
 	write(6,*) '***************************************************************'					      
	write(6,*) 'Betts 3D Lattice is utilized'
          else
            if(myrank.eq.0) write(lud,*) 'bad dimensionality'
          end if
	case(2)
       ! else if(ibetts.eq.1) then
         ! write(6,*) 'table-betts'
          if (ndim.eq.3) then
             call tables_3D
 	write(6,*) '***************************************************************'					      
	write(6,*) 'Conventional Lattice is utilized'
          else
            if(myrank.eq.0) write(lud,*) 'bad dimensionality'
          end if
	case(3)
          if (ndim.eq.3) then
             call  tables_3D_cubic
 	write(6,*) '***************************************************************'					      
	write(6,*) 'Perfect Cubic Lattice is utilized'
          else
            if(myrank.eq.0) write(lud,*) 'bad dimensionality'
          end if 
    	     case default
             print*, "Please, choose between case 1-3. See module_global.for for details"
             stop
              end select
                 
      	  if(myrank.eq.0.and.iprint.ge.1) write(lud,*)  'done tables';flush(lud)     
        return
        end    


      subroutine tables_3D_Betts      ! THREE DIMENSIONAL BETTS TABLES
c***************************************************************************************
c     This subroutine makes the lookup tables used throughout the
c     full THREE DIMENSIONAL code.  It employs the definitions and 
c     cluster geometries developed by D.D. Betts and G.E. Stewart,
c     Canadian Journal of Physics, 75, 47-66.


c***************************************************************************************
        use Global
        implicit none
c***************************************************************************************
      logical, parameter :: use_original = .false.
      real (kind) :: yy,del,deltat,tc
c***************************************************************************************
      interface transform3D
      function gtransform3D(a,igroup,ixy,ixz,iyz,isx,isy,isz)
c      use global
      integer :: a(3,3),gtransform3D(3,3),igroup,ixy(48),ixz(48),
     &           iyz(48),isx(48),isy(48),isz(48)
      end
      
      function Rtransform3D(Rc1,igroup,ixy,ixz,iyz,isx,isy,isz)
      integer :: Rc1(3),Rtransform3D(3),igroup,ixy(48),ixz(48),iyz(48),
     &     isx(48),isy(48),isz(48)
      end
      
      function Ktransform3D(Kc1,igroup,ixy,ixz,iyz,isx,isy,isz)
c      use global
      integer, parameter :: kind=8
      integer    :: ixy(48),ixz(48),iyz(48),
     &     isx(48),isy(48),isz(48)
      real(kind) :: Kc1(3),Ktransform3D(3)
      end
      end interface 
      
      interface Inverse
      function Inverse_3x3(a)
c      use global
      integer, parameter :: kind=8
      integer :: a(3,3)
      real(kind) :: Inverse_3x3(3,3)
      end function
      end interface
      
      interface determinant
      function Ideterminant(b)
c      use global
      integer, parameter :: kind=8
      integer :: b(3,3),Ideterminant
      end function
      function Rdeterminant(c)
      use global
      real(kind) :: c(3,3),Rdeterminant
      end function
      end interface

c***************************************************************************
c     Locally defined
c***************************************************************************
      logical    :: allowed,denied,denied_lin_alg
      integer, parameter :: maxnc=300
      integer, parameter :: nshellmax=100
      integer :: i,i1,i2,i3,it,ic,ic1,ic2,idim,I_F,Nm,I_FShell,
     &           ik,itest,j,jc,k,kroniker,
     &           l,m,n,p,Lc,iv1(3),iv2(3),
     &           ig(maxnc),jg(maxnc),kg(maxnc),a1(3,3),
     &           ixy(48),ixz(48),iyz(48),isx(48),isy(48),isz(48),
     &           irwgroup(maxnc),group(48),
     &		 p_r,nneigh(nshellmax),nneighl(nshellmax),e_index,
     &     nneighp(nshellmax),ishellmax,nshell,nneigh1(108,6)
      real(kind) :: r1,r2,r3,r4,Dkt,c1(3,3),ener(Nc),
     &              rv1(3),rv2(3),rv3(3),n_fact(Nc),
     &              corner(3,8),Bettsd(4),Bettsf(6),Bettsl(3),
     &              meand,meanf,meanl,cubicity,temp_min,temp_max
c For lin algebra test
      real(kind) :: ainv(3,3),a1ainv(3,3)
      real(kind) :: funct(Nc,-nwn:nwn),min_funct(Nc)
c***************************************************************************************
c Select use of "brute force" vs linear algebra method for 
c checking lattice coincidence. Both methods should agree; lin alg
c method is much faster.
c***************************************************************************************
      logical, parameter :: use_brute_force = .false.
      logical, parameter :: use_lin_alg = .true.

        if(myrank.eq.0.and.iprint.gt.19) then
          write(lud,*) 'start Table'
          flush(lud)
        end if
c***************************************************************************
c       Local arrays
c       ixy(igroup)      D=3 cubic point group table(see below)
c       ixz(igroup)      D=3 cubic point group table(see below)
c       iyz(igroup)      D=3 cubic point group table(see below)
c       isx(igroup)      D=3 cubic point group table(see below)
c       isy(igroup)      D=3 cubic point group table(see below)
c       isz(igroup)      D=3 cubic point group table(see below)
c       group(i)         Indices of allowed group operations
c       irwgroup(ic)     the group operation maps K(ic) into IRW
c       ig(ic)           i in the g1 and j in the g2 direction
c       jg(jc)           k in the g3 direction
c       kg(ic)           corresponding to K point ic
c       ict(i,j,k)       maps i*g1+j*g2+k*g3 to ic
c       nneigh(i)        number of neighbors on shell i within the cluster
c       nneighl(i)       number of neighbors on shell i on the lattice
c       nneighp(i)       number of neighbors on shell i on a perfect cluster
c***************************************************************************
        
c       CLUSTER MOMENTUM and DUAL SPACE for the 3D cubic lattice
c                                                             _   _   _
c       Different Real-space Cluster Geometries.  The vectors a1, a2, a3
c       define the parallelpipeds which define the periodically replicated 
c       clusters.  Here, we only use the clusters which Betts considers "good".
c               _             _            _
c  cluster      a1            a2           a3       group  Ng  I_F  bipartide 
c    1A     ( 1, 0, 0)    ( 0, 1, 0)   ( 0, 0, 1)     Oh   48   0     no
c
c    8A     ( 2, 0, 0)    ( 0, 2, 0)   ( 0, 0, 2)     Oh   48   4     yes
c    8B     ( 2, 1, 0)    ( 0, 1, 2)   ( 1,-1, 1)     C2h  4    0     no  (=C2 x Ci)
c
c   10A     ( 1, 0, 2)    ( 1, 2,-1)   ( 1,-2, 0)     Ci   2    0     no
c
c   12A     ( 2, 1, 1)    ( 1,-2, 0)   ( 1, 0,-2)     C2h  4    0     no
c   12B     ( 0, 2, 0)    ( 0, 0, 3)   ( 2, 0, 1)     C2h  4    2     no
c
c   14A     ( 2, 2, 2)    ( 2,-1, 1)   ( 1, 2,-1)     C3i  6    1    yes
c   14B     ( 2, 2, 1)    ( 0,-1, 2)   ( 2,-1, 0)     Ci   2    0     no
c   14C     ( 2, 2, 0)    (-1, 2, 1)   ( 1, 0, 2)     Ci   2    0     no
c
c   16A     ( 2, 2, 0)    ( 2, 0, 2)   ( 0, 2, 2)     Oh  48    2    yes*
c   16B     ( 2, 2, 0)    ( 1,-1,-2)   (-1, 1,-2)     D2h  8    2    yes* (=D2 x Ci)
c   16C     ( 4, 2, 0)    ( 1, 2,-1)   (-1, 2, 1)     D2h  8    2    yes*
c   16D     ( 2,-2, 0)    ( 0, 2,-2)   ( 1, 2, 1)     C2h  4    2    yes
c   16E     ( 2, 2, 0)    ( 2,-2, 0)   ( 1, 0,-2)     C2h  4    0     no
c   16F     ( 1, 2, 0)    (-3, 1,-1)   ( 0, 2, 2)     C2h  4    0     no
c   16Z     ( 2, 0, 0)    ( 0, 2, 0)   ( 0, 0, 4) [Made up BAD cluster]
c
c   18A     ( 2, 2, 2)    ( 2,-1,-1)   ( 1, 1,-2)     D3h 12    3    yes*
c   18B     ( 3, 2, 1)    ( 2,-1,-1)   ( 1, 1,-2)     C3i  6    3    yes*
c   18C     ( 2, 2, 2)    ( 2, 1,-1)   (-1, 2,-1)      Ci  2    3    yes
c   18D     ( 0, 0, 3)    ( 2,-2,-1)   ( 2, 1,-1)     C2h  4    1     no
c   18E     ( 3, 3, 0)    ( 2, 0, 2)   ( 2, 0,-1)     C2h  4    1     no
c   18F     ( 3, 2, 1)    ( 0, 2,-2)   ( 0, 2, 1)     C2h  4    2     no
c
c   20A     ( 3, 2, 1)    ( 2,-1,-1)   (-1, 1,-2)     Ci   2    4     yes
c   20B     ( 4, 1, 1)    ( 0, 2,-2)   ( 1, 2, 1)     Ci   2    4     yes
c   20C     ( 2, 2, 2)    ( 2, 0,-2)   (-1, 1,-2)     Ci   2    4     yes
c   20D     ( 3, 1, 0)    ( 0, 2, 2)   ( 1,-2, 1)     Ci   2    4     yes
c   20E     ( 4, 1, 1)    ( 2,-2, 0)   ( 1,-1,-2)     C2h  4    4     yes
c   20F     ( 2, 2, 1)    ( 2,-2, 1)   ( 1, 0,-2)     C2h  4    0     yes
c   20G     ( 2, 2, 1)    ( 2,-2, 0)   ( 1, 1,-2)     C2h  4    2     no
c
c   22A     ( 4, 2, 0)    ( 0, 3,-1)   (-1, 2, 1)     Ci   2    5     yes
c   22B     ( 2, 2, 2)    ( 3, 0,-1)   (-1, 1,-2)     Ci   2    5     yes
c   22C     ( 3, 2, 0)    ( 0, 2, 2)   ( 1,-1, 2)     Ci   2    4     no
c   22D     ( 1, 2, 0)    ( 2, 1, 2)   ( 2,-2, 1)     Ci   2          no  (BAD)
c   22E     ( 4, 1, 1)    ( 2, 2,-3)   ( 0, 2,-1)     Ci   2    6     no
c
c   24A     ( 3, 1, 0)    ( 0, 2,-2)   ( 0, 2, 2)     C2h  4    6     yes
c   24B     ( 2, 2, 2)    ( 0, 2,-2)   ( 2, 0,-2)     D3h 12    6     yes
c   24C     ( 2, 2, 2)    ( 3, 0,-1)   (-1, 2,-1)     Ci   2    6     yes
c   24D     ( 0, 0, 3)    ( 1, 3, 0)   (-2, 2, 0)     D2h  8    6     no
c   24E     ( 3, 2, 0)    ( 3,-2, 0)   ( 1, 0,-2)     C2h  4    4     no
c                                                           
c   26A     ( 3, 2, 1)    ( 2,-2,-2)   ( 1, 1,-2)     Ci   2    6     yes
c   26B     ( 1, 3, 2)    ( 3,-2, 1)   ( 3, 1, 0)     C3i  6    6     yes
c   26C     ( 2, 1, 3)    ( 2,-1,-2)   ( 2, 2,-1)     Ci   2    2     no
c
c   26Z     ( 1, 2, 3)    ( 3, 3, -2)  (-3,-2,-3)         20          yes
c
c  ----------------------------------------------------------------------------
c Some information on the perfection of the cluster and the completion of 
c the shells. Please, take a look to make sure that lattice and cluster is 
c complete for each shell for a given cluster size. This determines the 
c nearest neighbor, next, and next-next neighbor, etc, respectively.
c  ----------------------------------------------------------------------------
c  Neighbors of Cluster"2A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           1
c  Neighbors of Cluster"2B" Perf= 0 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           1
c  Neighbors of Cluster"2C" Perf= 0 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           1
c  Neighbors of Cluster"2D" Perf= 0 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           1
c  Neighbors of Cluster"4A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           3
c  Neighbors of Cluster"4B" Perf= 1 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           2
c           2          18           1
c  Neighbors of Cluster"4C" Perf= 1 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           2
c           2          18           1
c  Neighbors of Cluster"4D" Perf= 1 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           2
c           2          18           1
c  Neighbors of Cluster"6A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           5
c  Neighbors of Cluster"6B" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           3
c           2          18           2
c  Neighbors of Cluster"6C" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           3
c           2          18           2
c  Neighbors of Cluster"6D" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           3
c           2          18           2
c  Neighbors of Cluster"8A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           1
c  Neighbors of Cluster"8B" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           4
c           2          18           3
c  Neighbors of Cluster"8C" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           4
c           2          18           3
c  Neighbors of Cluster"8D" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           4
c           2          18           3
c  Neighbors of Cluster"10A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           3
c  Neighbors of Cluster"10B" Perf= 1 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           5
c           2          18           4
c  Neighbors of Cluster"10C" Perf= 1 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           5
c           2          18           4
c  Neighbors of Cluster"10D" Perf= 1 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           5
c           2          18           4
c  Neighbors of Cluster"12A" Perf= 0 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           5
c  Neighbors of Cluster"12B" Perf= 0 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           5
c  Neighbors of Cluster"12C" Perf= 0 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           5
c  Neighbors of Cluster"12D" Perf= 0 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           5
c  Neighbors of Cluster"14A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           7
c  Neighbors of Cluster"14B" Perf= 1 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           6
c           3          38           1
c  Neighbors of Cluster"14C" Perf= 1 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           6
c           3          38           1
c  Neighbors of Cluster"14D" Perf= 1 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           6
c           3          38           1
c  Neighbors of Cluster"16A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           9
c  Neighbors of Cluster"16B" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           7
c           3          38           2
c  Neighbors of Cluster"16C" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           7
c           3          38           2
c  Neighbors of Cluster"16D" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           7
c           3          38           2
c  Neighbors of Cluster"18A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          11
c  Neighbors of Cluster"18B" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           8
c           3          38           3
c  Neighbors of Cluster"18C" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           8
c           3          38           3
c  Neighbors of Cluster"18D" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           8
c           3          38           3
c  Neighbors of Cluster"20A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          13
c  Neighbors of Cluster"20B" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           9
c           3          38           4
c  Neighbors of Cluster"20C" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           9
c           3          38           4
c  Neighbors of Cluster"20D" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18           9
c           3          38           4
c  Neighbors of Cluster"22A" Perf= 1 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          14
c           3          38           1
c  Neighbors of Cluster"22B" Perf= 5 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          10
c           3          38           5
c  Neighbors of Cluster"22C" Perf= 5 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          10
c           3          38           5
c  Neighbors of Cluster"22D" Perf= 5 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          10
c           3          38           5
c  Neighbors of Cluster"24A" Perf= 1 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          16
c           3          38           1
c  Neighbors of Cluster"24B" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          11
c           3          38           6
c  Neighbors of Cluster"24C" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          11
c           3          38           6
c  Neighbors of Cluster"24D" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          11
c           3          38           6
c  Neighbors of Cluster"26A" Perf= 2 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          16
c           3          38           3
c  Neighbors of Cluster"26B" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          12
c           3          38           7
c  Neighbors of Cluster"26C" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          12
c           3          38           7
c  Neighbors of Cluster"26D" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          12
c           3          38           7
c  Neighbors of Cluster"28A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38           3
c  Neighbors of Cluster"28B" Perf= 5 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          13
c           3          38           8
c  Neighbors of Cluster"28C" Perf= 5 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          13
c           3          38           8
c  Neighbors of Cluster"28D" Perf= 5 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          13
c           3          38           8
c  Neighbors of Cluster"30A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38           5
c  Neighbors of Cluster"30B" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          14
c           3          38           9
c  Neighbors of Cluster"30C" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          14
c           3          38           9
c  Neighbors of Cluster"30D" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          14
c           3          38           9
c  Neighbors of Cluster"32A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38           7
c  Neighbors of Cluster"32B" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          15
c           3          38          10
c  Neighbors of Cluster"32C" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          15
c           3          38          10
c  Neighbors of Cluster"32D" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          15
c           3          38          10
c  Neighbors of Cluster"34A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38           9
c  Neighbors of Cluster"34B" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          16
c           3          38          11
c  Neighbors of Cluster"34C" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          16
c           3          38          11
c  Neighbors of Cluster"34D" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          16
c           3          38          11
c  Neighbors of Cluster"36A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          11
c  Neighbors of Cluster"36B" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          16
c           3          38          12
c           4          66           1
c  Neighbors of Cluster"36C" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          16
c           3          38          12
c           4          66           1
c  Neighbors of Cluster"36D" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          16
c           3          38          12
c           4          66           1
c  Neighbors of Cluster"38A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          13
c  Neighbors of Cluster"38B" Perf= 0 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          13
c  Neighbors of Cluster"38C" Perf= 0 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          13
c  Neighbors of Cluster"38D" Perf= 0 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          13
c  Neighbors of Cluster"40A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          15
c  Neighbors of Cluster"40B" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          17
c           3          38          14
c           4          66           2
c  Neighbors of Cluster"40C" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          17
c           3          38          14
c           4          66           2
c  Neighbors of Cluster"40D" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          17
c           3          38          14
c           4          66           2
c  Neighbors of Cluster"42A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          17
c  Neighbors of Cluster"42B" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          15
c           4          66           2
c  Neighbors of Cluster"42C" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          15
c           4          66           2
c  Neighbors of Cluster"42D" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          15
c           4          66           2
c  Neighbors of Cluster"44A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          19
c  Neighbors of Cluster"44B" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          16
c           4          66           3
c  Neighbors of Cluster"44C" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          16
c           4          66           3
c  Neighbors of Cluster"44D" Perf= 3 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          16
c           4          66           3
c  Neighbors of Cluster"46A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          21
c  Neighbors of Cluster"46B" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          17
c           4          66           4
c  Neighbors of Cluster"46C" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          17
c           4          66           4
c  Neighbors of Cluster"46D" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          17
c           4          66           4
c  Neighbors of Cluster"48A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          23
c  Neighbors of Cluster"48B" Perf= 5 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          18
c           4          66           5
c  Neighbors of Cluster"48C" Perf= 5 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          18
c           4          66           5
c  Neighbors of Cluster"48D" Perf= 5 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          18
c           4          66           5
c  Neighbors of Cluster"50A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c          2          18          18
c           3          38          25
c  Neighbors of Cluster"50B" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          19
c           4          66           6
c  Neighbors of Cluster"50C" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          19
c           4          66           6
c  Neighbors of Cluster"50D" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          19
c           4          66           6
c  Neighbors of Cluster"52A" Perf= 1 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          26
c           4          66           1
c  Neighbors of Cluster"52B" Perf= 7 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          20
c           4          66           7
c  Neighbors of Cluster"52C" Perf= 7 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          20
c           4          66           7
c  Neighbors of Cluster"52D" Perf= 7 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          20
c           4          66           7
c  Neighbors of Cluster"54A" Perf= 2 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          27
c           4          66           2
c  Neighbors of Cluster"54B" Perf= 8 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          21
c           4          66           8
c  Neighbors of Cluster"54C" Perf= 8 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          21
c           4          66           8
c  Neighbors of Cluster"54D" Perf= 8 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          21
c           4          66           8
c  Neighbors of Cluster"56A" Perf= 2 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          29
c           4          66           2
c  Neighbors of Cluster"56B" Perf= 9 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          22
c           4          66           9
c  Neighbors of Cluster"56C" Perf= 9 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          22
c           4          66           9
c  Neighbors of Cluster"56D" Perf= 9 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          22
c           4          66           9
c  Neighbors of Cluster"58A" Perf= 3 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          30
c           4          66           3
c  Neighbors of Cluster"58B" Perf=10 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          23
c           4          66          10
c  Neighbors of Cluster"58C" Perf=10 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          23
c           4          66          10
c  Neighbors of Cluster"58D" Perf=10 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          23
c           4          66          10
c  Neighbors of Cluster"60A" Perf= 4 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          31
c           4          66           4
c  Neighbors of Cluster"60B" Perf=11 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          24
c           4          66          11
c  Neighbors of Cluster"60C" Perf=11 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          24
c           4          66          11
c  Neighbors of Cluster"60D" Perf=11 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          24
c           4          66          11
c  Neighbors of Cluster"62A" Perf= 5 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          32
c           4          66           5
c  Neighbors of Cluster"62B" Perf=12 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          25
c           4          66          12
c  Neighbors of Cluster"62C" Perf=12 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          25
c           4          66          12
c  Neighbors of Cluster"62D" Perf=12 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          25
c           4          66          12
c  Neighbors of Cluster"64A" Perf= 5 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          17
c           3          38          34
c           4          66           6
c  Neighbors of Cluster"64B" Perf=12 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          26
c           4          66          13
c  Neighbors of Cluster"64C" Perf=12 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          26
c           4          66          13
c  Neighbors of Cluster"64D" Perf=12 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          26
c           4          66          13
c  Neighbors of Cluster"66A" Perf= 4 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          34
c           4          66           7
c  Neighbors of Cluster"66B" Perf=11 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          27
c           4          66          14
c  Neighbors of Cluster"66C" Perf=11 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          27
c           4          66          14
c  Neighbors of Cluster"66D" Perf=11 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          27
c           4          66          14
c  Neighbors of Cluster"68A" Perf= 3 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          35
c           4          66           8
c  Neighbors of Cluster"68B" Perf=10 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          28
c           4          66          15
c  Neighbors of Cluster"68C" Perf=10 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          28
c           4          66          15
c  Neighbors of Cluster"68D" Perf=10 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          28
c           4          66          15
c  Neighbors of Cluster"70A" Perf= 1 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          37
c           4          66           8
c  Neighbors of Cluster"70B" Perf= 9 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          29
c           4          66          16
c  Neighbors of Cluster"70C" Perf= 9 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          29
c           4          66          16
c  Neighbors of Cluster"70D" Perf= 9 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          29
c           4          66          16
c  Neighbors of Cluster"72A" Perf= 1 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          37
c           4          66          10
c  Neighbors of Cluster"72B" Perf= 8 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          30
c           4          66          17
c  Neighbors of Cluster"72C" Perf= 8 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          30
c           4          66          17
c  Neighbors of Cluster"72D" Perf= 8 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          30
c           4          66          17
c  Neighbors of Cluster"74A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          11
c  Neighbors of Cluster"74B" Perf= 7 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          31
c           4          66          18
c  Neighbors of Cluster"74C" Perf= 9 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          30
c           4          66          18
c           5         102           1
c  Neighbors of Cluster"74D" Perf=10 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          30
c           4          66          17
c           5         102           1
c           6         146           1
c  Neighbors of Cluster"76A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          13
c  Neighbors of Cluster"76B" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          32
c           4          66          19
c  Neighbors of Cluster"76C" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          32
c           4          66          19
c  Neighbors of Cluster"76D" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          32
c           4          66          19
c  Neighbors of Cluster"78A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          15
c  Neighbors of Cluster"78B" Perf= 5 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          33
c           4          66          20
c  Neighbors of Cluster"78C" Perf= 9 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          31
c           4          66          20
c           5         102           2
c  Neighbors of Cluster"78D" Perf= 9 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          31
c           4          66          20
c           5         102           2
c  Neighbors of Cluster"80A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          17
c  Neighbors of Cluster"80B" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          34
c           4          66          21
c  Neighbors of Cluster"80C" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          34
c           4          66          21
c  Neighbors of Cluster"80D" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          34
c           4          66          21
c  Neighbors of Cluster"82A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          19
c  Neighbors of Cluster"82B" Perf= 7 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          33
c           4          66          22
c           5         102           2
c  Neighbors of Cluster"82C" Perf= 7 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          33
c           4          66          22
c           5         102           2
c  Neighbors of Cluster"82D" Perf= 8 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          33
c           4          66          21
c           5         102           2
c           6         146           1
c  Neighbors of Cluster"84A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          21
c  Neighbors of Cluster"84B" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          34
c           4          66          23
c           5         102           2
c  Neighbors of Cluster"84C" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          34
c           4          66          23
c           5         102           2
c  Neighbors of Cluster"84D" Perf= 6 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          34
c           4          66          23
c           5         102           2
c  Neighbors of Cluster"86A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          23
c  Neighbors of Cluster"86B" Perf= 5 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          35
c           4          66          24
c           5         102           2
c  Neighbors of Cluster"86C" Perf= 5 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          35
c           4          66          24
c           5         102           2
c  Neighbors of Cluster"86D" Perf= 5 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          35
c           4          66          24
c           5         102           2
c  Neighbors of Cluster"88A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          25
c  Neighbors of Cluster"88B" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          36
c           4          66          25
c           5         102           2
c  Neighbors of Cluster"88C" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          36
c           4          66          25
c           5         102           2
c  Neighbors of Cluster"88D" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          36
c           4          66          25
c           5         102           2
c  Neighbors of Cluster"90A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          27
c  Neighbors of Cluster"90B" Perf= 7 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          35
c           4          66          26
c           5         102           4
c  Neighbors of Cluster"90C" Perf= 7 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          35
c           4          66          26
c           5         102           4
c  Neighbors of Cluster"90D" Perf= 7 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          35
c           4          66          26
c           5         102           4
c  Neighbors of Cluster"92A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          29
c  Neighbors of Cluster"92B" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          27
c           5         102           2
c  Neighbors of Cluster"92C" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          27
c           5         102           2
c  Neighbors of Cluster"92D" Perf= 2 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          27
c           5         102           2
c  Neighbors of Cluster"94A" Perf= 1 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          37
c           4          66          32
c  Neighbors of Cluster"94B" Perf= 7 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          36
c           4          66          28
c           5         102           5
c  Neighbors of Cluster"94C" Perf= 9 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          35
c           4          66          28
c           5         102           6
c  Neighbors of Cluster"94D" Perf= 9 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          35
c           4          66          28
c           5         102           6
c  Neighbors of Cluster"96A" Perf= 0 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          33
c  Neighbors of Cluster"96B" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          29
c           5         102           4
c  Neighbors of Cluster"96C" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          29
c           5         102           4
c  Neighbors of Cluster"96D" Perf= 4 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          29
c           5         102           4
c  Neighbors of Cluster"98A" Perf= 1 Bipartide=F
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          38
c           4          66          34
c           5         102           1
c  Neighbors of Cluster"98B" Perf= 7 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          37
c           4          66          30
c           5         102           6
c  Neighbors of Cluster"98C" Perf= 8 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          37
c           4          66          29
c           5         102           6
c           6         146           1
c  Neighbors of Cluster"98D" Perf= 9 Bipartide=T
c         Shell      Lattice      Cluster
c           1           6           6
c           2          18          18
c           3          38          36
c           4          66          30
c          5         102           7
c
c  The bipartide lattices had a epsilon_0^HA in table 8 (and common sense), 
c  I_F (the feerro imperfection) was calculated below and *not* taken from Betts
c  A * means equivlent to others with same Nc.
c  I.e., lattices 16A=16B=16C, 18A=18B,  20B=20C=20D.
c
c  Larger lattices were calculated with gen_cluster.f Each "A" cluster is the
c  most perfect cluster with the smallest cubicity for a given size. Each "B"
c  cluster is the most perfect cluster for a given size that is also bipartide.
c
c  ----------------------------------------------------------------------------
c                       GROUP OPERATIONS
c  ----------------------------------------------------------------------------
c  Only 16A has all 48 operations and has the point group Oh.  The remaining clusters 
c  have lower symmetry than the cubic lattice, with symmetries designated by the
c  Schoenflies symbols (From Ibach and Luth):
c  Cj (j=2,3,4, 6) j-fold rotation axis 
c  Sj j-fold rotation-inversion axis 
c  Dj j 2-fold rotation axes perpendicular to a j-fold principle rotation axis 
c  T  4 three-and 3 two-fold rotation axes, as in a tetrahedron 
c  O  4 three-and 3 four-fold rotation axes, as in a octahedron 
c  Ci a center of inversion 
c  Cs a mirror plane
c  In addition, their are sufixes for mirror planes 
c     h: horizontal=perpendicular to the c  rotation axis, 
c     v: vertical=parallel to the main rotation axis in the plane, 
c     d: diagonal=parallel to the main rotation axis in the plane  
c        bisecting the two-fold rotation axes.
c 
c  The Clusters identified by eetts, have several different symmetries
c
c  C2h   A two fold axis, plus a single mirror plane that is perpendicular to the axis
c  Ci    Contains only the identity and inversion
c  C3i   Triagonal with a 3-fold axis and inversion (GUESS)??
c          * rotations by 2pi/3 about the axix x=y=z
c          * inversion wrt the origin.
c  Oh    The full symmetry of the cubic lattice 
c  D2h   2 2-fold axes perp. to the 2-fold rotation axis, plus a horizonal mirror plane
c  D3h   3 2-fold axes perp. to the 3-fold rotation axis, plus a horizonal mirror plane
c
c  We need to automate the process of finding which of the 48 point group operations are
c  retained by the different clusters.  To do this, we must perform all 48 cubic group
c  operations on the {a} and see if they form an equivalent parallelpiped.  
c
c  Group Operations for a cubic lattice:
c
c  x --> +/- x                  There are 2^D D! operations in D dimensions
c    \/                         or 48 operations in 3D.
c    /\
c  y --> +/- y   
c    \/
c    /\
c  z --> +/- z   
c        .
c  a'=transform3D(a,igroup,ixy,ixz,iyz,isx,isy,isz)
c  ixy=0(1) dont (do) exchange x and y
c  ixz=0(1) dont (do) exchange x and z
c  iyz=0(1) dont (do) exchange y and z
c  isx      sign of x 
c  isy      sign of y
c  isz      sign of z
c
c  x <-- isx*[(1-ixy)*x + ixy*y]  <-- [(1-ixz)*x + ixz*z]<----------------------
c  y <-- isy*[(1-ixy)*y + ixy*x]  <----------------------<-- [(1-iyz)*y + iyz*z]
c  z <-- isz*<------------------- <-- [(1-ixz)*z + ixz*x]<-- [(1-iyz)*z + iyz*y]
c

      ixy( 1)=0;ixz( 1)=0; iyz( 1)=0; isx( 1)=+1; isy( 1)=+1; isz( 1)=+1 !identity  xyz e
      ixy( 2)=1;ixz( 2)=0; iyz( 2)=0; isx( 2)=+1; isy( 2)=+1; isz( 2)=+1 !ref. x=y  yxz o
      ixy( 3)=0;ixz( 3)=1; iyz( 3)=0; isx( 3)=+1; isy( 3)=+1; isz( 3)=+1 !ref. x=z  zyx o
      ixy( 4)=0;ixz( 4)=0; iyz( 4)=1; isx( 4)=+1; isy( 4)=+1; isz( 4)=+1 !ref. y=z  xzy o
      ixy( 5)=1;ixz( 5)=1; iyz( 5)=0; isx( 5)=+1; isy( 5)=+1; isz( 5)=+1 !          yzx e
      ixy( 6)=1;ixz( 6)=0; iyz( 6)=1; isx( 6)=+1; isy( 6)=+1; isz( 6)=+1 !          zxy e
      
      ixy( 7)=0;ixz( 7)=0; iyz( 7)=0; isx( 7)=-1; isy( 7)=+1; isz( 7)=+1 !
      ixy( 8)=1;ixz( 8)=0; iyz( 8)=0; isx( 8)=-1; isy( 8)=+1; isz( 8)=+1 !
      ixy( 9)=0;ixz( 9)=1; iyz( 9)=0; isx( 9)=-1; isy( 9)=+1; isz( 9)=+1 !
      ixy(10)=0;ixz(10)=0; iyz(10)=1; isx(10)=-1; isy(10)=+1; isz(10)=+1 !
      ixy(11)=1;ixz(11)=1; iyz(11)=0; isx(11)=-1; isy(11)=+1; isz(11)=+1 !
      ixy(12)=1;ixz(12)=0; iyz(12)=1; isx(12)=-1; isy(12)=+1; isz(12)=+1 !
     
      ixy(13)=0;ixz(13)=0; iyz(13)=0; isx(13)=+1; isy(13)=-1; isz(13)=+1 !
      ixy(14)=1;ixz(14)=0; iyz(14)=0; isx(14)=+1; isy(14)=-1; isz(14)=+1 !
      ixy(15)=0;ixz(15)=1; iyz(15)=0; isx(15)=+1; isy(15)=-1; isz(15)=+1 !
      ixy(16)=0;ixz(16)=0; iyz(16)=1; isx(16)=+1; isy(16)=-1; isz(16)=+1 !
      ixy(17)=1;ixz(17)=1; iyz(17)=0; isx(17)=+1; isy(17)=-1; isz(17)=+1 !
      ixy(18)=1;ixz(18)=0; iyz(18)=1; isx(18)=+1; isy(18)=-1; isz(18)=+1 !
      
      ixy(19)=0;ixz(19)=0; iyz(19)=0; isx(19)=-1; isy(19)=-1; isz(19)=+1 !
      ixy(20)=1;ixz(20)=0; iyz(20)=0; isx(20)=-1; isy(20)=-1; isz(20)=+1 !
      ixy(21)=0;ixz(21)=1; iyz(21)=0; isx(21)=-1; isy(21)=-1; isz(21)=+1 !
      ixy(22)=0;ixz(22)=0; iyz(22)=1; isx(22)=-1; isy(22)=-1; isz(22)=+1 !
      ixy(23)=1;ixz(23)=1; iyz(23)=0; isx(23)=-1; isy(23)=-1; isz(23)=+1 !
      ixy(24)=1;ixz(24)=0; iyz(24)=1; isx(24)=-1; isy(24)=-1; isz(24)=+1 !
      
      ixy(25)=0;ixz(25)=0; iyz(25)=0; isx(25)=+1; isy(25)=+1; isz(25)=-1 !
      ixy(26)=1;ixz(26)=0; iyz(26)=0; isx(26)=+1; isy(26)=+1; isz(26)=-1 !
      ixy(27)=0;ixz(27)=1; iyz(27)=0; isx(27)=+1; isy(27)=+1; isz(27)=-1 !
      ixy(28)=0;ixz(28)=0; iyz(28)=1; isx(28)=+1; isy(28)=+1; isz(28)=-1 !
      ixy(29)=1;ixz(29)=1; iyz(29)=0; isx(29)=+1; isy(29)=+1; isz(29)=-1 !
      ixy(30)=1;ixz(30)=0; iyz(30)=1; isx(30)=+1; isy(30)=+1; isz(30)=-1 !
      
      ixy(31)=0;ixz(31)=0; iyz(31)=0; isx(31)=-1; isy(31)=+1; isz(31)=-1 !
      ixy(32)=1;ixz(32)=0; iyz(32)=0; isx(32)=-1; isy(32)=+1; isz(32)=-1 !
      ixy(33)=0;ixz(33)=1; iyz(33)=0; isx(33)=-1; isy(33)=+1; isz(33)=-1 !
      ixy(34)=0;ixz(34)=0; iyz(34)=1; isx(34)=-1; isy(34)=+1; isz(34)=-1 !
      ixy(35)=1;ixz(35)=1; iyz(35)=0; isx(35)=-1; isy(35)=+1; isz(35)=-1 !
      ixy(36)=1;ixz(36)=0; iyz(36)=1; isx(36)=-1; isy(36)=+1; isz(36)=-1 !
      
      ixy(37)=0;ixz(37)=0; iyz(37)=0; isx(37)=+1; isy(37)=-1; isz(37)=-1 !
      ixy(38)=1;ixz(38)=0; iyz(38)=0; isx(38)=+1; isy(38)=-1; isz(38)=-1 !
      ixy(39)=0;ixz(39)=1; iyz(39)=0; isx(39)=+1; isy(39)=-1; isz(39)=-1 !
      ixy(40)=0;ixz(40)=0; iyz(40)=1; isx(40)=+1; isy(40)=-1; isz(40)=-1 !
      ixy(41)=1;ixz(41)=1; iyz(41)=0; isx(41)=+1; isy(41)=-1; isz(41)=-1 !
      ixy(42)=1;ixz(42)=0; iyz(42)=1; isx(42)=+1; isy(42)=-1; isz(42)=-1 !
      
      ixy(43)=0;ixz(43)=0; iyz(43)=0; isx(43)=-1; isy(43)=-1; isz(43)=-1 !
      ixy(44)=1;ixz(44)=0; iyz(44)=0; isx(44)=-1; isy(44)=-1; isz(44)=-1 !
      ixy(45)=0;ixz(45)=1; iyz(45)=0; isx(45)=-1; isy(45)=-1; isz(45)=-1 !
      ixy(46)=0;ixz(46)=0; iyz(46)=1; isx(46)=-1; isy(46)=-1; isz(46)=-1 !
      ixy(47)=1;ixz(47)=1; iyz(47)=0; isx(47)=-1; isy(47)=-1; isz(47)=-1 !
      ixy(48)=1;ixz(48)=0; iyz(48)=1; isx(48)=-1; isy(48)=-1; isz(48)=-1 !

c***********************************************************************        
c       determine the principle translation vectors a1,a2, and a3 
c***********************************************************************
      if(iprint.ge.19.and.myrank.eq.0)	
     &   write(lud,*) ' cluster= ',cluster !0

      If(cluster.eq.'1A') then                              
         a(1,:)=(/1, 0, 0/); a(2,:)=(/0, 1, 0/); a(3,:)=(/0, 0, 1/)
      else if(cluster.eq."2A") then
! Cubicity=  1.00208840 Perfection=  0 Betts Shell=  1 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 0, 2/); a(2,:)=(/ 1, 1,-2/); a(3,:)=(/-4,-1,-2/)
      else if(cluster.eq."2B") then
! Cubicity=  1.00793005 Perfection=  0 Betts Shell=  1 Bipartide=T ID#=  5
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4,-1,-1/); a(3,:)=(/ 2, 0, 0/)
      else if(cluster.eq."2C") then
! Cubicity=  1.02906655 Perfection=  0 Betts Shell=  1 Bipartide=T ID#= 54
         a(1,:)=(/ 1, 0, 1/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 2, 1, 1/)
      else if(cluster.eq."2D") then
! Cubicity=  1.02970629 Perfection=  0 Betts Shell=  1 Bipartide=T ID#= 57
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 3, 3,-4/); a(3,:)=(/-1,-1, 2/)
      else if(cluster.eq."4A") then
! Cubicity=  1.00626104 Perfection=  0 Betts Shell=  1 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 0, 2/); a(2,:)=(/ 4, 1, 3/); a(3,:)=(/-1,-1,-1/)
      else if(cluster.eq."4B") then
! Cubicity=  1.00425949 Perfection=  1 Betts Shell=  1 Bipartide=T ID#=***
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 1, 1/); a(3,:)=(/ 2,-2,-4/)
      else if(cluster.eq."4C") then
! Cubicity=  1.01405382 Perfection=  1 Betts Shell=  1 Bipartide=T ID#=***
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 3,-3/); a(3,:)=(/ 3, 3, 2/)
      else if(cluster.eq."4D") then
! Cubicity=  1.01448493 Perfection=  1 Betts Shell=  1 Bipartide=T ID#=***
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 3,-1/); a(3,:)=(/ 0, 0, 4/)
      else if(cluster.eq."6A") then
! Cubicity=  1.00173407 Perfection=  0 Betts Shell=  1 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 1, 1, 1/); a(3,:)=(/ 1,-1,-1/)
      else if(cluster.eq."6B") then
! Cubicity=  1.01424527 Perfection=  2 Betts Shell=  1 Bipartide=T ID#=***
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 4, 1,-1/); a(3,:)=(/ 2, 0, 0/)
      else if(cluster.eq."6C") then
! Cubicity=  1.01570937 Perfection=  2 Betts Shell=  1 Bipartide=T ID#=***
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 1,-1/); a(3,:)=(/ 2, 0, 0/)
      else if(cluster.eq."6D") then
! Cubicity=  1.01869738 Perfection=  2 Betts Shell=  1 Bipartide=T ID#=***
         a(1,:)=(/ 1, 0, 1/); a(2,:)=(/ 2, 1,-1/); a(3,:)=(/ 1, 2, 1/)
      else if(cluster.eq."8A") then
! Cubicity=  1.00503637 Perfection=  0 Betts Shell=  2 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 0, 2/); a(2,:)=(/ 1, 2, 0/); a(3,:)=(/ 1,-1,-1/)
      else if(cluster.eq."8B") then
! Cubicity=  1.00918638 Perfection=  2 Betts Shell=  2 Bipartide=T ID#=718
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 3, 1, 4/); a(3,:)=(/ 3,-1,-2/)
      else if(cluster.eq."8C") then
! Cubicity=  1.01251883 Perfection=  2 Betts Shell=  2 Bipartide=T ID#=720
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4,-3,-1/); a(3,:)=(/-2, 1, 1/)
      else if(cluster.eq."8D") then
! Cubicity=  1.01345813 Perfection=  2 Betts Shell=  2 Bipartide=T ID#=721
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 2, 0,-2/); a(3,:)=(/-3,-1,-4/)
      else if(cluster.eq."10A") then
! Cubicity=  1.00291827 Perfection=  0 Betts Shell=  2 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 0, 2/); a(2,:)=(/ 4, 0,-2/); a(3,:)=(/ 1, 1,-2/)
      else if(cluster.eq."10B") then
! Cubicity=  1.00504956 Perfection=  1 Betts Shell=  2 Bipartide=T ID#=285
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 4, 1,-3/); a(3,:)=(/-2,-1,-1/)
      else if(cluster.eq."10C") then
! Cubicity=  1.01272001 Perfection=  1 Betts Shell=  2 Bipartide=T ID#=292
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 3, 3, 4/); a(3,:)=(/ 1, 1,-2/)
      else if(cluster.eq."10D") then
! Cubicity=  1.01754753 Perfection=  1 Betts Shell=  2 Bipartide=T ID#=297
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 3, 2, 1/); a(3,:)=(/ 3,-1,-4/)
      else if(cluster.eq."12A") then
! Cubicity=  1.00397119 Perfection=  0 Betts Shell=  2 Bipartide=T ID#=  1
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 1,-3/); a(3,:)=(/ 2, 2, 0/)
      else if(cluster.eq."12B") then
! Cubicity=  1.00999985 Perfection=  0 Betts Shell=  2 Bipartide=T ID#=  7
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 1,-3/); a(3,:)=(/ 3, 3, 2/)
      else if(cluster.eq."12C") then
! Cubicity=  1.01745058 Perfection=  0 Betts Shell=  2 Bipartide=T ID#= 18
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 2,-4, 2/); a(3,:)=(/-2, 1,-1/)
      else if(cluster.eq."12D") then
! Cubicity=  1.01774583 Perfection=  0 Betts Shell=  2 Bipartide=T ID#= 19
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4,-2, 2/); a(3,:)=(/ 1,-2,-3/)
      else if(cluster.eq."14A") then
! Cubicity=  1.00195100 Perfection=  0 Betts Shell=  2 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 0, 2/); a(2,:)=(/ 3, 1, 0/); a(3,:)=(/ 1,-2, 0/)
      else if(cluster.eq."14B") then
! Cubicity=  1.00759080 Perfection=  1 Betts Shell=  2 Bipartide=T ID#=272
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 2, 1, 1/); a(3,:)=(/ 1, 4,-3/)
      else if(cluster.eq."14C") then
! Cubicity=  1.01131134 Perfection=  1 Betts Shell=  2 Bipartide=T ID#=278
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4,-1,-3/); a(3,:)=(/ 3, 2, 1/)
      else if(cluster.eq."14D") then
! Cubicity=  1.01838258 Perfection=  1 Betts Shell=  2 Bipartide=T ID#=289
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 2, 1,-1/); a(3,:)=(/ 1,-2, 1/)
      else if(cluster.eq."16A") then
! Cubicity=  1.00364153 Perfection=  0 Betts Shell=  2 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 0, 2/); a(2,:)=(/ 0, 2, 4/); a(3,:)=(/-2,-2, 0/)
      else if(cluster.eq."16B") then
! Cubicity=  1.01104614 Perfection=  2 Betts Shell=  2 Bipartide=T ID#=869
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 1, 3/); a(3,:)=(/ 0, 2,-2/)
      else if(cluster.eq."16C") then
! Cubicity=  1.01192146 Perfection=  2 Betts Shell=  2 Bipartide=T ID#=871
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 2,-2, 0/); a(3,:)=(/ 1, 1,-2/)
      else if(cluster.eq."16D") then
! Cubicity=  1.01210811 Perfection=  2 Betts Shell=  2 Bipartide=T ID#=872
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 1, 3/); a(3,:)=(/ 2, 1,-3/)
      else if(cluster.eq."18A") then
! Cubicity=  1.00289626 Perfection=  0 Betts Shell=  2 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 0, 2/); a(2,:)=(/ 2, 2,-1/); a(3,:)=(/ 1,-2,-2/)
      else if(cluster.eq."18B") then
! Cubicity=  1.00465339 Perfection=  3 Betts Shell=  2 Bipartide=T ID#=***
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 2,-1,-1/)
      else if(cluster.eq."18C") then
! Cubicity=  1.00553296 Perfection=  3 Betts Shell=  2 Bipartide=T ID#=***
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 2, 1,-1/); a(3,:)=(/ 1,-2,-1/)
      else if(cluster.eq."18D") then
! Cubicity=  1.00578376 Perfection=  3 Betts Shell=  2 Bipartide=T ID#=***
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 3, 1,-4/); a(3,:)=(/ 2,-1, 1/)
      else if(cluster.eq."20A") then
! Cubicity=  1.01120197 Perfection=  0 Betts Shell=  2 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 0, 2/); a(2,:)=(/ 3, 2, 1/); a(3,:)=(/ 2,-2,-1/)
      else if(cluster.eq."20B") then
! Cubicity=  1.00650759 Perfection=  4 Betts Shell=  2 Bipartide=T ID#=***
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 1, 1/); a(3,:)=(/ 0, 2,-2/)
      else if(cluster.eq."20C") then
! Cubicity=  1.00719961 Perfection=  4 Betts Shell=  2 Bipartide=T ID#=***
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 2, 4,-2/); a(3,:)=(/ 2, 1,-1/)
      else if(cluster.eq."20D") then
! Cubicity=  1.00723268 Perfection=  4 Betts Shell=  2 Bipartide=T ID#=***
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 2, 4, 0/); a(3,:)=(/ 1,-2, 1/)
      else if(cluster.eq."22A") then
! Cubicity=  1.01355078 Perfection=  1 Betts Shell=  2 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 3, 1,-1/); a(3,:)=(/ 2,-2, 1/)
      else if(cluster.eq."22B") then
! Cubicity=  1.00750218 Perfection=  5 Betts Shell=  2 Bipartide=T ID#=782
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 2, 4, 0/); a(3,:)=(/ 2, 1,-1/)
      else if(cluster.eq."22C") then
! Cubicity=  1.01279707 Perfection=  5 Betts Shell=  2 Bipartide=T ID#=791
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 4, 1,-1/); a(3,:)=(/ 1, 2,-1/)
      else if(cluster.eq."22D") then
! Cubicity=  1.01327852 Perfection=  5 Betts Shell=  2 Bipartide=T ID#=793
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 0, 2/); a(3,:)=(/ 2, 1,-3/)
      else if(cluster.eq."24A") then
! Cubicity=  1.07490229 Perfection=  1 Betts Shell=  2 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 3, 0, 0/); a(3,:)=(/ 1, 2,-2/)
      else if(cluster.eq."24B") then
! Cubicity=  1.00539291 Perfection=  6 Betts Shell=  2 Bipartide=T ID#=***
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 3, 3, 0/); a(3,:)=(/ 2,-2,-2/)
      else if(cluster.eq."24C") then
! Cubicity=  1.00663346 Perfection=  6 Betts Shell=  2 Bipartide=T ID#=***
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 1, 4,-1/); a(3,:)=(/ 1,-2,-1/)
      else if(cluster.eq."24D") then
! Cubicity=  1.00813937 Perfection=  6 Betts Shell=  2 Bipartide=T ID#=***
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 0, 2/); a(3,:)=(/ 0, 4, 0/)
      else if(cluster.eq."26A") then
! Cubicity=  1.05390076 Perfection=  2 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 2/); a(2,:)=(/ 2, 2,-1/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."26B") then
! Cubicity=  1.00526694 Perfection=  6 Betts Shell=  3 Bipartide=T ID#=553
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."26C") then
! Cubicity=  1.00587487 Perfection=  6 Betts Shell=  3 Bipartide=T ID#=554
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 2, 1,-1/); a(3,:)=(/ 1,-3,-2/)
      else if(cluster.eq."26D") then
! Cubicity=  1.00703271 Perfection=  6 Betts Shell=  3 Bipartide=T ID#=556
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 2, 3, 1/); a(3,:)=(/ 2,-1,-1/)
      else if(cluster.eq."28A") then
! Cubicity=  1.06296811 Perfection=  0 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 3,-1, 1/); a(3,:)=(/ 1, 2,-2/)
      else if(cluster.eq."28B") then
! Cubicity=  1.01547442 Perfection=  5 Betts Shell=  3 Bipartide=T ID#=695
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 0, 0/); a(3,:)=(/ 1, 2,-3/)
      else if(cluster.eq."28C") then
! Cubicity=  1.01563841 Perfection=  5 Betts Shell=  3 Bipartide=T ID#=696
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4,-1, 1/); a(3,:)=(/ 2, 1,-3/)
      else if(cluster.eq."28D") then
! Cubicity=  1.01651932 Perfection=  5 Betts Shell=  3 Bipartide=T ID#=701
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 0, 0/); a(3,:)=(/ 1,-4,-1/)
      else if(cluster.eq."30A") then
! Cubicity=  1.00723632 Perfection=  0 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 2/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 2,-2, 1/)
      else if(cluster.eq."30B") then
! Cubicity=  1.00966474 Perfection=  4 Betts Shell=  3 Bipartide=T ID#=390
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 3,-2,-1/); a(3,:)=(/ 1, 1,-4/)
      else if(cluster.eq."30C") then
! Cubicity=  1.01095376 Perfection=  4 Betts Shell=  3 Bipartide=T ID#=391
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 4, 2, 0/); a(3,:)=(/ 3,-2, 1/)
      else if(cluster.eq."30D") then
! Cubicity=  1.01243421 Perfection=  4 Betts Shell=  3 Bipartide=T ID#=392
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 3, 1,-2/); a(3,:)=(/ 3,-2, 1/)
      else if(cluster.eq."32A") then
! Cubicity=  1.02208161 Perfection=  0 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 2,-2, 1/)
      else if(cluster.eq."32B") then
! Cubicity=  1.02748952 Perfection=  3 Betts Shell=  3 Bipartide=T ID#=319
         a(1,:)=(/ 2, 0, 2/); a(2,:)=(/ 3, 4, 3/); a(3,:)=(/ 2,-2,-2/)
      else if(cluster.eq."32C") then
! Cubicity=  1.02814321 Perfection=  3 Betts Shell=  3 Bipartide=T ID#=320
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 2, 0,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."32D") then
! Cubicity=  1.03975266 Perfection=  3 Betts Shell=  3 Bipartide=T ID#=330
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 4,-2, 0/); a(3,:)=(/ 2, 0,-2/)
      else if(cluster.eq."34A") then
! Cubicity=  1.00907609 Perfection=  0 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 3,-2, 0/); a(3,:)=(/ 1, 2,-2/)
      else if(cluster.eq."34B") then
! Cubicity=  1.05741559 Perfection=  2 Betts Shell=  3 Bipartide=T ID#=105
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 1,-3,-2/)
      else if(cluster.eq."34C") then
! Cubicity=  1.06605340 Perfection=  2 Betts Shell=  3 Bipartide=T ID#=117
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 2, 3, 1/); a(3,:)=(/ 2,-2,-2/)
      else if(cluster.eq."34D") then
! Cubicity=  1.06656422 Perfection=  2 Betts Shell=  3 Bipartide=T ID#=118
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 3, 2, 1/); a(3,:)=(/ 3,-1,-4/)
      else if(cluster.eq."36A") then
! Cubicity=  1.00419451 Perfection=  0 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 2/); a(2,:)=(/ 3, 0,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."36B") then
! Cubicity=  1.02544507 Perfection=  3 Betts Shell=  3 Bipartide=T ID#=461
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 2, 4, 4/); a(3,:)=(/ 2, 2,-4/)
      else if(cluster.eq."36C") then
! Cubicity=  1.04017240 Perfection=  3 Betts Shell=  3 Bipartide=T ID#=495
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 2,-2,-2/)
      else if(cluster.eq."36D") then
! Cubicity=  1.05175695 Perfection=  3 Betts Shell=  3 Bipartide=T ID#=520
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 4, 2, 2/); a(3,:)=(/ 2,-2,-2/)
      else if(cluster.eq."38A") then
! Cubicity=  1.00210705 Perfection=  0 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 3, 1,-3/); a(3,:)=(/ 2,-2, 1/)
      else if(cluster.eq."38B") then
! Cubicity=  1.08723039 Perfection=  0 Betts Shell=  3 Bipartide=T ID#= 29
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-1,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."38C") then
! Cubicity=  1.11696799 Perfection=  0 Betts Shell=  3 Bipartide=T ID#= 37
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-1,-2/); a(3,:)=(/ 2, 3,-1/)
      else if(cluster.eq."38D") then
! Cubicity=  1.14352525 Perfection=  0 Betts Shell=  3 Bipartide=T ID#= 43
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 3, 2, 1/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."40A") then
! Cubicity=  1.00283386 Perfection=  0 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 2,-2,-1/)
      else if(cluster.eq."40B") then
! Cubicity=  1.04098349 Perfection=  3 Betts Shell=  3 Bipartide=T ID#=547
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."40C") then
! Cubicity=  1.07455219 Perfection=  3 Betts Shell=  3 Bipartide=T ID#=610
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 2,-2,-2/)
      else if(cluster.eq."40D") then
! Cubicity=  1.07559994 Perfection=  3 Betts Shell=  3 Bipartide=T ID#=611
         a(1,:)=(/ 2, 2, 2/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."42A") then
! Cubicity=  1.00339418 Perfection=  0 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 2/); a(2,:)=(/ 4, 1,-2/); a(3,:)=(/ 2,-3, 0/)
      else if(cluster.eq."42B") then
! Cubicity=  1.05586176 Perfection=  2 Betts Shell=  3 Bipartide=T ID#=224
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-1, 2/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."42C") then
! Cubicity=  1.05685802 Perfection=  2 Betts Shell=  3 Bipartide=T ID#=226
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 2, 1,-3/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."42D") then
! Cubicity=  1.08112408 Perfection=  2 Betts Shell=  3 Bipartide=T ID#=247
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 2,-2,-2/)
      else if(cluster.eq."44A") then
! Cubicity=  1.00292586 Perfection=  0 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 2/); a(2,:)=(/ 3, 1,-2/); a(3,:)=(/ 2,-4, 0/)
      else if(cluster.eq."44B") then
! Cubicity=  1.03550915 Perfection=  3 Betts Shell=  3 Bipartide=T ID#=261
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."44C") then
! Cubicity=  1.07179973 Perfection=  3 Betts Shell=  3 Bipartide=T ID#=298
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 3, 1,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."44D") then
! Cubicity=  1.07552907 Perfection=  3 Betts Shell=  3 Bipartide=T ID#=300
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 2,-2, 2/); a(3,:)=(/ 1, 4,-3/)
      else if(cluster.eq."46A") then
! Cubicity=  1.01111785 Perfection=  0 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 2/); a(2,:)=(/ 4, 2, 0/); a(3,:)=(/ 2,-4, 1/)
      else if(cluster.eq."46B") then
! Cubicity=  1.01749200 Perfection=  4 Betts Shell=  3 Bipartide=T ID#=271
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 1,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."46C") then
! Cubicity=  1.04917476 Perfection=  4 Betts Shell=  3 Bipartide=T ID#=296
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 1,-3, 2/)
      else if(cluster.eq."46D") then
! Cubicity=  1.05356388 Perfection=  4 Betts Shell=  3 Bipartide=T ID#=299
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."48A") then
! Cubicity=  1.00077008 Perfection=  0 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 4,-2,-1/)
      else if(cluster.eq."48B") then
! Cubicity=  1.00194693 Perfection=  5 Betts Shell=  3 Bipartide=T ID#=643
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-2, 1/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."48C") then
! Cubicity=  1.03330670 Perfection=  5 Betts Shell=  3 Bipartide=T ID#=752
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 1,-2/); a(3,:)=(/ 3,-2, 1/)
      else if(cluster.eq."48D") then
! Cubicity=  1.03490335 Perfection=  5 Betts Shell=  3 Bipartide=T ID#=758
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 3,-2, 1/)
      else if(cluster.eq."50A") then
! Cubicity=  1.08506285 Perfection=  0 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 2,-3,-2/); a(3,:)=(/-3,-2, 0/)
      else if(cluster.eq."50B") then
! Cubicity=  1.01815330 Perfection=  6 Betts Shell=  3 Bipartide=T ID#=555
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 2,-3, 1/)
      else if(cluster.eq."50C") then
! Cubicity=  1.02361407 Perfection=  6 Betts Shell=  3 Bipartide=T ID#=566
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 3,-2,-1/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."50D") then
! Cubicity=  1.02591793 Perfection=  6 Betts Shell=  3 Bipartide=T ID#=571
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 3,-2, 1/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."52A") then
! Cubicity=  1.01420868 Perfection=  1 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 0, 4/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 1,-3,-1/)
      else if(cluster.eq."52B") then
! Cubicity=  1.00288782 Perfection=  7 Betts Shell=  3 Bipartide=T ID#=559
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 1,-2/); a(3,:)=(/ 2,-3, 1/)
      else if(cluster.eq."52C") then
! Cubicity=  1.01621458 Perfection=  7 Betts Shell=  3 Bipartide=T ID#=579
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 4, 0,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."52D") then
! Cubicity=  1.02363889 Perfection=  7 Betts Shell=  3 Bipartide=T ID#=593
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 2,-4, 0/)
      else if(cluster.eq."54A") then
! Cubicity=  1.00777859 Perfection=  2 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 2/); a(2,:)=(/ 3, 1,-3/); a(3,:)=(/ 2,-2, 4/)
      else if(cluster.eq."54B") then
! Cubicity=  1.00461045 Perfection=  8 Betts Shell=  3 Bipartide=T ID#=650
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-3, 0/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."54C") then
! Cubicity=  1.00958360 Perfection=  8 Betts Shell=  3 Bipartide=T ID#=659
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."54D") then
! Cubicity=  1.02803190 Perfection=  8 Betts Shell=  3 Bipartide=T ID#=713
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 3,-3, 0/)
      else if(cluster.eq."56A") then
! Cubicity=  1.06932775 Perfection=  2 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 4, 0,-2/); a(3,:)=(/-3, 2,-2/)
      else if(cluster.eq."56B") then
! Cubicity=  1.00603823 Perfection=  9 Betts Shell=  3 Bipartide=T ID#=755
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 2,-4, 2/)
      else if(cluster.eq."56C") then
! Cubicity=  1.00721056 Perfection=  9 Betts Shell=  3 Bipartide=T ID#=757
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 4,-2, 0/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."56D") then
! Cubicity=  1.01528712 Perfection=  9 Betts Shell=  3 Bipartide=T ID#=773
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-4,-3/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."58A") then
! Cubicity=  1.01370257 Perfection=  3 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 4,-2, 4/); a(3,:)=(/ 3, 2,-2/)
      else if(cluster.eq."58B") then
! Cubicity=  1.00431818 Perfection= 10 Betts Shell=  3 Bipartide=T ID#=446
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 4, 3,-1/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."58C") then
! Cubicity=  1.00538357 Perfection= 10 Betts Shell=  3 Bipartide=T ID#=447
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 1,-4/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."58D") then
! Cubicity=  1.01101713 Perfection= 10 Betts Shell=  3 Bipartide=T ID#=453
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-3, 2/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."60A") then
! Cubicity=  1.00088509 Perfection=  4 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 2, 0, 3/); a(2,:)=(/ 2, 3,-2/); a(3,:)=(/ 2,-3,-2/)
      else if(cluster.eq."60B") then
! Cubicity=  1.00967632 Perfection= 11 Betts Shell=  3 Bipartide=T ID#=822
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 4,-2, 0/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."60C") then
! Cubicity=  1.01047548 Perfection= 11 Betts Shell=  3 Bipartide=T ID#=825
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-4,-1/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."60D") then
! Cubicity=  1.01114956 Perfection= 11 Betts Shell=  3 Bipartide=T ID#=830
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-3, 2/); a(3,:)=(/ 2, 1,-3/)
      else if(cluster.eq."62A") then
! Cubicity=  1.01784219 Perfection=  5 Betts Shell=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 4/); a(2,:)=(/ 3, 0,-2/); a(3,:)=(/ 2,-3, 2/)
      else if(cluster.eq."62B") then
! Cubicity=  1.00308829 Perfection= 12 Betts Shell=  3 Bipartide=T ID#=337
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."62C") then
! Cubicity=  1.02249726 Perfection= 12 Betts Shell=  3 Bipartide=T ID#=355
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-3, 4/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."62D") then
! Cubicity=  1.02482374 Perfection= 12 Betts Shell=  3 Bipartide=T ID#=360
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 4,-3,-1/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."64A") then
! Cubicity=  1.11459751 Perfection=  5 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 2, 4,-2/); a(3,:)=(/ 2,-3, 2/)
      else if(cluster.eq."64B") then
! Cubicity=  1.00186622 Perfection= 12 Betts Shell=  4 Bipartide=T ID#=522
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 4, 2,-2/); a(3,:)=(/ 3,-2, 1/)
      else if(cluster.eq."64C") then
! Cubicity=  1.00282252 Perfection= 12 Betts Shell=  4 Bipartide=T ID#=523
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 4,-1/); a(3,:)=(/ 3,-2, 1/)
      else if(cluster.eq."64D") then
! Cubicity=  1.00344091 Perfection= 12 Betts Shell=  4 Bipartide=T ID#=525
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-2, 1/); a(3,:)=(/ 2, 4,-2/)
      else if(cluster.eq."66A") then
! Cubicity=  1.10811031 Perfection=  4 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 4/); a(2,:)=(/ 3,-2, 2/); a(3,:)=(/ 2, 1,-4/)
      else if(cluster.eq."66B") then
! Cubicity=  1.01066253 Perfection= 11 Betts Shell=  4 Bipartide=T ID#=407
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 3, 0,-3/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."66C") then
! Cubicity=  1.02569249 Perfection= 11 Betts Shell=  4 Bipartide=T ID#=437
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 0,-3/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."66D") then
! Cubicity=  1.03330335 Perfection= 11 Betts Shell=  4 Bipartide=T ID#=453
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 3, 2,-3/); a(3,:)=(/ 2,-3,-1/)
      else if(cluster.eq."68A") then
! Cubicity=  1.05145635 Perfection=  3 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 4/); a(2,:)=(/ 2, 2,-3/); a(3,:)=(/-3, 2,-2/)
      else if(cluster.eq."68B") then
! Cubicity=  1.01821557 Perfection= 10 Betts Shell=  4 Bipartide=T ID#=286
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 2,-4, 0/)
      else if(cluster.eq."68C") then
! Cubicity=  1.02042012 Perfection= 10 Betts Shell=  4 Bipartide=T ID#=288
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 4, 2, 0/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."68D") then
! Cubicity=  1.02151537 Perfection= 10 Betts Shell=  4 Bipartide=T ID#=290
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 2, 4,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."70A") then
! Cubicity=  1.05956510 Perfection=  1 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 4/); a(2,:)=(/ 3, 1,-3/); a(3,:)=(/ 3,-2, 2/)
      else if(cluster.eq."70B") then
! Cubicity=  1.02798300 Perfection=  9 Betts Shell=  4 Bipartide=T ID#=315
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/ 2,-4, 0/)
      else if(cluster.eq."70C") then
! Cubicity=  1.03037072 Perfection=  9 Betts Shell=  4 Bipartide=T ID#=320
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 4, 0,-2/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."70D") then
! Cubicity=  1.03408521 Perfection=  9 Betts Shell=  4 Bipartide=T ID#=326
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/ 3,-2, 3/)
      else if(cluster.eq."72A") then
! Cubicity=  1.01027869 Perfection=  1 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 3,-2, 2/)
      else if(cluster.eq."72B") then
! Cubicity=  1.04698290 Perfection=  8 Betts Shell=  4 Bipartide=T ID#=331
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-3, 2/); a(3,:)=(/ 2, 4,-2/)
      else if(cluster.eq."72C") then
! Cubicity=  1.04995627 Perfection=  8 Betts Shell=  4 Bipartide=T ID#=337
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 4,-2,-2/); a(3,:)=(/ 3, 4,-1/)
      else if(cluster.eq."72D") then
! Cubicity=  1.05040701 Perfection=  8 Betts Shell=  4 Bipartide=T ID#=340
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 4,-1,-3/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."74A") then
! Cubicity=  1.05679459 Perfection=  0 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 2,-2, 4/)
      else if(cluster.eq."74B") then
! Cubicity=  1.02163739 Perfection=  7 Betts Shell=  4 Bipartide=T ID#=119
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 4,-1,-3/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."74C") then
! Cubicity=  1.02702810 Perfection=  9 Betts Shell=  4 Bipartide=T ID#=260
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-3, 2/); a(3,:)=(/ 2, 3,-3/)
      else if(cluster.eq."74D") then
! Cubicity=  1.04744508 Perfection= 10 Betts Shell=  4 Bipartide=T ID#=383
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 4,-2,-4/); a(3,:)=(/-3, 3,-2/)
      else if(cluster.eq."76A") then
! Cubicity=  1.01047869 Perfection=  0 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 3,-3, 1/)
      else if(cluster.eq."76B") then
! Cubicity=  1.02568825 Perfection=  6 Betts Shell=  4 Bipartide=T ID#=108
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 4,-1/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."76C") then
! Cubicity=  1.02666879 Perfection=  6 Betts Shell=  4 Bipartide=T ID#=110
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-3, 2/); a(3,:)=(/ 2, 2,-4/)
      else if(cluster.eq."76D") then
! Cubicity=  1.08009543 Perfection=  6 Betts Shell=  4 Bipartide=T ID#=131
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 4,-2, 2/); a(3,:)=(/ 3, 2,-3/)
      else if(cluster.eq."78A") then
! Cubicity=  1.04329249 Perfection=  0 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 4/); a(2,:)=(/ 3,-2, 2/); a(3,:)=(/ 2, 3,-3/)
      else if(cluster.eq."78B") then
! Cubicity=  1.10279190 Perfection=  5 Betts Shell=  4 Bipartide=T ID#=115
         a(1,:)=(/ 2, 3, 3/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/ 3,-3,-2/)
      else if(cluster.eq."78C") then
! Cubicity=  1.00243983 Perfection=  9 Betts Shell=  4 Bipartide=T ID#=450
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."78D") then
! Cubicity=  1.01559986 Perfection=  9 Betts Shell=  4 Bipartide=T ID#=467
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 4,-1, 3/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."80A") then
! Cubicity=  1.03158985 Perfection=  0 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 4/); a(2,:)=(/ 3, 1,-3/); a(3,:)=(/ 3,-3, 1/)
      else if(cluster.eq."80B") then
! Cubicity=  1.05410609 Perfection=  4 Betts Shell=  4 Bipartide=T ID#= 68
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 3, 2,-3/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."80C") then
! Cubicity=  1.06772187 Perfection=  4 Betts Shell=  4 Bipartide=T ID#= 73
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 4,-2, 2/); a(3,:)=(/ 2, 3,-3/)
      else if(cluster.eq."80D") then
! Cubicity=  1.07851804 Perfection=  4 Betts Shell=  4 Bipartide=T ID#= 77
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 4, 3, 1/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."82A") then
! Cubicity=  1.00158475 Perfection=  0 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 4/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 3,-3, 1/)
      else if(cluster.eq."82B") then
! Cubicity=  1.07422003 Perfection=  7 Betts Shell=  4 Bipartide=T ID#=156
         a(1,:)=(/ 2, 0, 4/); a(2,:)=(/ 3, 2,-3/); a(3,:)=(/ 2,-3,-3/)
      else if(cluster.eq."82C") then
! Cubicity=  1.11679270 Perfection=  7 Betts Shell=  4 Bipartide=T ID#=175
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 3,-3,-2/); a(3,:)=(/-3,-2, 3/)
      else if(cluster.eq."82D") then
! Cubicity=  1.11679270 Perfection=  8 Betts Shell=  4 Bipartide=T ID#=245
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 3, 2,-3/); a(3,:)=(/ 3,-3,-2/)
      else if(cluster.eq."84A") then
! Cubicity=  1.00700892 Perfection=  0 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 4/); a(2,:)=(/ 4, 1,-2/); a(3,:)=(/ 2,-3, 2/)
      else if(cluster.eq."84B") then
! Cubicity=  1.05048444 Perfection=  6 Betts Shell=  4 Bipartide=T ID#=252
         a(1,:)=(/ 2, 2, 4/); a(2,:)=(/ 3, 0,-3/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."84C") then
! Cubicity=  1.05418398 Perfection=  6 Betts Shell=  4 Bipartide=T ID#=256
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 4,-2,-2/); a(3,:)=(/ 2, 3,-3/)
      else if(cluster.eq."84D") then
! Cubicity=  1.05447790 Perfection=  6 Betts Shell=  4 Bipartide=T ID#=257
         a(1,:)=(/ 1, 1, 4/); a(2,:)=(/ 4, 1,-3/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."86A") then
! Cubicity=  1.00713644 Perfection=  0 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 4/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."86B") then
! Cubicity=  1.07279883 Perfection=  5 Betts Shell=  4 Bipartide=T ID#= 87
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 2, 3,-3/); a(3,:)=(/ 2,-4, 0/)
      else if(cluster.eq."86C") then
! Cubicity=  1.09705570 Perfection=  5 Betts Shell=  4 Bipartide=T ID#= 98
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 4, 2, 0/); a(3,:)=(/ 1, 4,-3/)
      else if(cluster.eq."86D") then
! Cubicity=  1.13387392 Perfection=  5 Betts Shell=  4 Bipartide=T ID#=107
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 3,-1, 4/); a(3,:)=(/ 2, 3,-3/)
      else if(cluster.eq."88A") then
! Cubicity=  1.01571974 Perfection=  0 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 4, 2,-2/); a(3,:)=(/-3, 3,-1/)
      else if(cluster.eq."88B") then
! Cubicity=  1.08854203 Perfection=  4 Betts Shell=  4 Bipartide=T ID#= 87
         a(1,:)=(/ 2, 2, 4/); a(2,:)=(/ 4, 2,-2/); a(3,:)=(/-3, 3, 2/)
      else if(cluster.eq."88C") then
! Cubicity=  1.10221909 Perfection=  4 Betts Shell=  4 Bipartide=T ID#= 89
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 4,-2,-2/); a(3,:)=(/ 3, 3,-2/)
      else if(cluster.eq."88D") then
! Cubicity=  1.10424407 Perfection=  4 Betts Shell=  4 Bipartide=T ID#= 90
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/-2, 2,-4/)
      else if(cluster.eq."90A") then
! Cubicity=  1.01430822 Perfection=  0 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 4/); a(2,:)=(/ 3, 3,-3/); a(3,:)=(/-3, 2,-2/)
      else if(cluster.eq."90B") then
! Cubicity=  1.03337789 Perfection=  7 Betts Shell=  4 Bipartide=T ID#=283
         a(1,:)=(/ 2, 0, 4/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/ 2,-3,-3/)
      else if(cluster.eq."90C") then
! Cubicity=  1.05141314 Perfection=  7 Betts Shell=  4 Bipartide=T ID#=299
         a(1,:)=(/ 2, 3, 3/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."90D") then
! Cubicity=  1.05901809 Perfection=  7 Betts Shell=  4 Bipartide=T ID#=311
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 4, 0,-2/); a(3,:)=(/ 3,-2, 3/)
      else if(cluster.eq."92A") then
! Cubicity=  1.01189995 Perfection=  0 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 4/); a(2,:)=(/ 4, 2,-3/); a(3,:)=(/ 3,-2, 2/)
      else if(cluster.eq."92B") then
! Cubicity=  1.08525112 Perfection=  2 Betts Shell=  4 Bipartide=T ID#= 15
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 3,-2, 3/); a(3,:)=(/ 2, 4,-2/)
      else if(cluster.eq."92C") then
! Cubicity=  1.10173577 Perfection=  2 Betts Shell=  4 Bipartide=T ID#= 16
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 4,-1,-3/); a(3,:)=(/ 3,-2, 3/)
      else if(cluster.eq."92D") then
! Cubicity=  1.11859778 Perfection=  2 Betts Shell=  4 Bipartide=T ID#= 17
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 4,-1,-3/); a(3,:)=(/ 2, 4,-2/)
      else if(cluster.eq."94A") then
! Cubicity=  1.02841762 Perfection=  1 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 4/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 3,-3, 4/)
      else if(cluster.eq."94B") then
! Cubicity=  1.06338945 Perfection=  7 Betts Shell=  4 Bipartide=T ID#=123
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 3,-2, 3/); a(3,:)=(/ 2, 3,-3/)
      else if(cluster.eq."94C") then
! Cubicity=  1.08824835 Perfection=  9 Betts Shell=  4 Bipartide=T ID#=224
         a(1,:)=(/ 2, 3, 3/); a(2,:)=(/ 4, 4,-2/); a(3,:)=(/ 3,-3,-2/)
      else if(cluster.eq."94D") then
! Cubicity=  1.10270358 Perfection=  9 Betts Shell=  4 Bipartide=T ID#=231
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 3, 2,-3/); a(3,:)=(/ 3,-1, 4/)
      else if(cluster.eq."96A") then
! Cubicity=  1.03544043 Perfection=  0 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 2, 2, 4/); a(2,:)=(/ 4,-2, 1/); a(3,:)=(/-3,-3, 2/)
      else if(cluster.eq."96B") then
! Cubicity=  1.06142184 Perfection=  4 Betts Shell=  4 Bipartide=T ID#= 65
         a(1,:)=(/ 2, 2, 4/); a(2,:)=(/ 4,-2, 2/); a(3,:)=(/ 3, 2,-3/)
      else if(cluster.eq."96C") then
! Cubicity=  1.06257992 Perfection=  4 Betts Shell=  4 Bipartide=T ID#= 68
         a(1,:)=(/ 2, 2, 4/); a(2,:)=(/ 4,-2, 2/); a(3,:)=(/ 2, 3,-3/)
      else if(cluster.eq."96D") then
! Cubicity=  1.06462040 Perfection=  4 Betts Shell=  4 Bipartide=T ID#= 69
         a(1,:)=(/ 2, 2, 4/); a(2,:)=(/ 3,-3,-2/); a(3,:)=(/-2,-2, 4/)
      else if(cluster.eq."98A") then
! Cubicity=  1.02323401 Perfection=  1 Betts Shell=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 2, 2, 3/); a(2,:)=(/ 3,-4, 2/); a(3,:)=(/ 2, 2,-4/)
      else if(cluster.eq."98B") then
! Cubicity=  1.08051541 Perfection=  7 Betts Shell=  4 Bipartide=T ID#=126
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/ 4,-1,-3/); a(3,:)=(/-3, 3,-2/)
      else if(cluster.eq."98C") then
! Cubicity=  1.08051541 Perfection=  8 Betts Shell=  4 Bipartide=T ID#=175
         a(1,:)=(/ 1, 3, 4/); a(2,:)=(/-4, 1, 3/); a(3,:)=(/-3, 3,-2/)
      else if(cluster.eq."98D") then
! Cubicity=  1.07583908 Perfection=  9 Betts Shell=  4 Bipartide=T ID#=228
         a(1,:)=(/ 2, 3, 3/); a(2,:)=(/ 3,-2,-3/); a(3,:)=(/-4, 4,-2/)
      else if(cluster.eq."100A") then
! Cubicity=  1.01194070 Perfection=  8 Bipartide=T ID#=  1
         a(1,:)=(/ 2, 3, 3/); a(2,:)=(/ 3, 2,-3/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."102A") then
! Cubicity=  1.00756315 Perfection=  5 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 3, 3,-3/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."108A") then
! Cubicity=  1.09583955 Perfection=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 3, 3, 3/); a(2,:)=(/ 3,-3,-3/); a(3,:)=(/-3, 3,-3/)
      else
         if(myrank.eq.0) write(lud,*) 'Bad cluster. Go directly to',
     &        ' jail. Do not pass go, do not collect $200'
         stop
      end if  
      if(iprint.ge.19.and.myrank.eq.0)       
     &  write(lud,"('a=',3('(',i2,1x,i2,1x,i2,') '))") 
     &		a(1,:),a(2,:),a(3,:)  

c     Test the a vectors      
      if(abs(determinant(a)).ne.Nc) then
        if(myrank.eq.0) write(lud,*) 'This cell volume isnt Nc',
     &                                abs(determinant(a))
        stop
      end if

c     Check fixed dimensions (PK: Could be allocated)
      if (nc.gt.maxnc) then
        if(myrank.eq.0) write(lud,*) 'Nc>maxnc parameter',nc,maxnc,
     &        'Need to recompile'
        stop
      end if
c     Calculate the cubicity of the clusters      
c
c       c8-------c7    body diagonals   face diagonals
c      /|       /|     --------------   --------------
c     / |      / |       d1 (c1,c7)       f1 (c1,c3)
c    /  |     /  |       d2 (c2,c8)       f2 (c2,c4)
c   c5-------c6  |       d3 (c3,c5)       f3 (c1,c6)
c   |   c4---|---c3      d4 (c4,c6)       f4 (c2,c5)
c a3|  /     |  /                         f5 (c1,c9)
c   | / a2   | /                          f6 (c4,c5)
c   |/       |/
c   c1-------c2
c        a1
c
      corner(:,1) = zeror
      corner(:,2) = a(1,:)
      corner(:,3) = a(1,:)+a(2,:)
      corner(:,4) = a(2,:)
      corner(:,5) = a(3,:)
      corner(:,6) = a(3,:)+a(1,:)
      corner(:,7) = a(3,:)+a(1,:)+a(2,:)
      corner(:,8) = a(3,:)+a(2,:)
      
      Bettsd(1)=sqrt(sum((corner(:,1)-corner(:,7))**2))
      Bettsd(2)=sqrt(sum((corner(:,2)-corner(:,8))**2))
      Bettsd(3)=sqrt(sum((corner(:,3)-corner(:,5))**2))
      Bettsd(4)=sqrt(sum((corner(:,4)-corner(:,6))**2))
      
      Bettsf(1)=sqrt(sum((corner(:,1)-corner(:,3))**2))
      Bettsf(2)=sqrt(sum((corner(:,2)-corner(:,4))**2))
      Bettsf(3)=sqrt(sum((corner(:,1)-corner(:,6))**2))
      Bettsf(4)=sqrt(sum((corner(:,2)-corner(:,5))**2))
      Bettsf(5)=sqrt(sum((corner(:,1)-corner(:,8))**2))
      Bettsf(6)=sqrt(sum((corner(:,4)-corner(:,5))**2))
      
      Bettsl(1)=sqrt(dot_product(corner(:,2),corner(:,2)))
      Bettsl(2)=sqrt(dot_product(corner(:,4),corner(:,4)))
      Bettsl(3)=sqrt(dot_product(corner(:,5),corner(:,5)))
      
      meand=product(Bettsd)**(0.25000000000_kind)
      meanf=product(Bettsf)**(0.16666666667_kind)
      meanl=product(Bettsl)**(0.33333333333_kind)
      r1=1.732050808_kind*meanl/meand
      r2=1.414213562_kind*meanl/meanf
      cubicity=max(r1,oner/r1)*max(r2,oner/r2)
      if(myrank.eq.0) 
     &  write(lud,*) 'cubicity=',cubicity
      

c*********************************************************************************
c     Determine which group operations are allowed    
c*********************************************************************************
      
      ngroup=0; Lc=2*(Nc**.3333+1)
      if(myrank.eq.0.and.iprint.ge.19) then
        write(lud,*) '--------------------------------------------',
     &               '------------------' 
        write(lud,*) '             trans(a)        igroup L  ixy ix',
     &               'z iyz isx isy isz'
        write(lud,*) '--------------------------------------------',
     &               '------------------'
      end if
      do igroup=1,48
         a1=transform3D(a,igroup,ixy,ixz,iyz,isx,isy,isz)

         if (use_brute_force) then

c        look at the lattice points generated by a1 and see if they
c        correspond to those produced by a.
c        First, generate the a1 lattice points
         do i=-Lc,Lc
         do j=-Lc,Lc
         do k=-Lc,Lc
           denied=.true.
           iv1=(/i,j,k/)
           do l=1,3         ! x,ymz components of this lattice point
             rv1(l)=sum(iv1(:)*a1(:,l))
           end do
c          Now generate the a lattice points.  If a and a1 are the 
c          same lattice one of these points MUST correspond to ij
           do n=-4*Lc,4*Lc
           do m=-4*Lc,4*Lc
           do p=-4*Lc,4*Lc
             iv2=(/n,m,p/)
             do l=1,3         ! x,y components of this lattice point
               rv2(l)=sum(iv2(:)*a(:,l))
             end do
             if(sum(abs(rv1-rv2)).lt.epsilon) denied=.false.
           end do
           end do
           end do
           if(denied) exit !no corrs. point was found
         end do
         if(denied) exit
         end do
         if(denied) exit
         end do
         end if
         if (use_lin_alg) then
c PK If lattice points generated by "a1" vectors can be found in 
c lattice points generated by "a" vectors then (a1)x(a^-1) can only 
c contain integers.
            ainv=Inverse_3x3(a)
            a1ainv=matmul(a1,ainv)
            if (all(abs(a1ainv-nint(a1ainv)).le.epsilon)) then 
               denied_lin_alg=.false. ! matrix only contains integers
            else
               denied_lin_alg=.true.
            end if
         end if

c If both algorithms are enabled, check they agree
         if (use_brute_force.and.use_lin_alg) then
            if (denied.neqv.denied_lin_alg) then
               write (*,*) '*** BUG : Denied algorithm disagree'
               write (*,*) 'a11',a1(:,1)
               write (*,*) 'a12',a1(:,2)
               write (*,*) 'a13',a1(:,3)
               
               write (*,*) 'a1',a(:,1)
               write (*,*) 'a2',a(:,2)
               write (*,*) 'a3',a(:,3)
               write (*,*)
               write (*,*) 'ainv',ainv(:,1)
               write (*,*) 'ainv',ainv(:,2)
               write (*,*) 'ainv',ainv(:,3)
               write (*,*)
               
               write (*,*) 'a1ainv',a1ainv(:,1)
               write (*,*) 'a1ainv',a1ainv(:,2)
               write (*,*) 'a1ainv',a1ainv(:,3)
               write (*,*) 'Denied_lin_alg=',denied_lin_alg
               write (*,*) 'Denied_brute_force=',denied
               write (*,*)
               stop
            end if
         end if
         if (use_lin_alg) denied=denied_lin_alg

         if(denied.eqv..false.) then ! a corresponding point was found
            ngroup=ngroup+1
            group(ngroup)=igroup
         end if    

         if(myrank.eq.0.and.iprint.ge.19) then
           allowed=.true.;if (denied) allowed=.false.
           write(lud,"(3('(',i2,1x,i2,1x,i2,')'),2x,i2,2x,l,6(2x,i2))")
     &     a1(1,:),a1(2,:),a1(3,:),igroup,allowed,ixy(igroup),
     &     ixz(igroup),iyz(igroup),isx(igroup),isy(igroup),isz(igroup)
         end if
      end do

     
      if(myrank.eq.0.and.iprint.ge.19) 
     &  write(lud,"('cluster=',a3,' ngroup=',i2,' groups=',48(i2,1x))") 
     &  cluster,ngroup,group(1:ngroup)

c*********************************************************************************
c       calculate the table of cluster point locations Rc(i,ic), i=1,2,3 for x,y,z
c*********************************************************************************       
c 
c  ^ y
c  |
cc-|-c---c-*-c---c---c---c               A point P is within the cluster 
c| | |   |   |   |   |   |               if a vector v from the box origin, 
cc-|-c---c---c---c---c---c               O, to the point satisfies
ca2^ * * * * * * * * |   |               
cc-|-c---c---c---c-*-c---c                        
c| | |   P   |   | * |   |                            _   _
cc-|-c--/c---c---c-*-c---c                   -epsilon  <  v__ai  < 1-epsilon
c| | | /v|   |   | * |   |                        
cc-|-c/--c---c---c-*-c---c 
c| | /   |   |   | * |   |               for all i=1,3 where v__ai is the
cc-|/0---c---c---c-*-c---c               projection of the vector v along
c| O--------------->---------> x         ai
cc---c---c---c---ca1-c---c                              -1
c|   |   |   |   |   |   |               v = A n;  n = A  v ; where A=transpose(a)
cc---c---c---c---c---c---c
c                                        v = i*a1 + j*a2 + k*a3
c
c                                            [ a1 ]
c                                        a = [ a2 ]
c                                            [ a3 ]
c                                        
cc      do i=1,3
cc        a1(:,i)=a(i,:)
cc      end do
	a1=transpose(a)
      c1=Inverse(a1)
      if(iprint.ge.19.and.myrank.eq.0) write(lud,*) '  ' !0
      if(iprint.ge.19.and.myrank.eq.0) 
     &     write(lud,*) 'ic,Rcx(ic),Rcy(ic),Rcz(ic)' !0
      ic=0
      do i=-Nc,Nc
      do j=-Nc,Nc
      do k=-Nc,Nc
         rv1=(/i,j,k/)          ! The coordinates of the point.
         rv2=matmul(c1,rv1)
         if(rv2(1).lt.oner-epsilon.and.rv2(1).gt.zeror-epsilon.and.
     &      rv2(2).lt.oner-epsilon.and.rv2(2).gt.zeror-epsilon.and.
     &      rv2(3).lt.oner-epsilon.and.rv2(3).gt.zeror-epsilon ) then
            ic=ic+1
            Rc(:,ic)=(/i,j,k/) ! Rc, like i.e. a(1,1) is integer
            if(iprint.ge.19.and.myrank.eq.0)                             !0
     &           write(lud,"(1x,i2,3x,3(i2,6x))") ic,Rc(1,ic),Rc(2,ic), !0
     &                                            Rc(3,ic)              !0
            if(sum(abs(Rc(:,ic))).eq.0) icR0=ic
         end if
      end do
      end do
      end do
      if(ic.ne.Nc) then
         if(myrank.eq.0) write(lud,*) 'Bug in Rc(ic) loop' !0
         stop
      end if
      if(iprint.ge.19.and.myrank.eq.0) write(lud,*) 'icR0=',icR0
 

c***********************************************************************        
c     create a table mapping points outside the cluster back into it.
c***********************************************************************        
c     Note we use c1 from above!!        

      if(iprint.ge.19.and.myrank.eq.0) then                  !0
         write(lud,*) '  '                                  !0
         write(lud,*) 'Equivalency table r'                 !0
         write(lud,*) 'ic   igroup   icrequ(ic,ieqiv)'      !0
      end if                                                !0
      do ic=1,Nc
        do igroup=1,ngroup
c       Generate the equivalent points using the ngroup operations
          rv3=transform3D(Rc(:,ic),group(igroup),ixy,ixz,iyz,
     &                    isx,isy,isz)

! Map rv3 back into cluster
          rv2=matmul(c1,rv3)    !project rv3 onto a1,a2,a3
          do while (any(rv2.lt.zeror-epsilon).or.
     &		any(rv2.gt.oner-epsilon))
             do l=1,3
                if(rv2(l).lt.zeror-epsilon) then
                   rv3(:)=rv3(:)+a(l,:)
                else if(rv2(l).gt.oner-epsilon) then
                   rv3(:)=rv3(:)-a(l,:)
                end if
             end do
             rv2=matmul(c1,rv3) !project rv3 onto a1,a2,a3
          end do
c         Figure out which point rv3 is
          itest=0
          do jc=1,Nc
            if(sum(abs(rv3(:)-Rc(:,jc))).lt.epsilon) then ! rv3=Rc(:,ic)
               icrequ(ic,igroup)=jc
                itest=1
               exit
            end if
          end do
          if(itest.eq.0) then
            write(lud,*) 'No equivalent found in icrequ',ic,igroup
            stop
          end if 
          if(iprint.ge.19.and.myrank.eq.0) 
     &      write(lud,"(2x,i2,4x,i2,8x,i2)") ic,igroup,icrequ(ic,igroup) !0
        end do
      end do                               
      

c***********************************************************************        
c     Now create the tables for the differences of R
c***********************************************************************        
c     Note we use c1 from above!!   
      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) '  ' !0
        write(lud,*) 'differences in r' !0
        write(lud,*) '  ic   jc  icrdiff(ic,jc)' !0
      end if
      do i=1,Nc
      do j=1,Nc
         iv1(:)=Rc(:,i)-Rc(:,j)
         itest=0
         do n=-1,1        ! Look in the surrounding 
         do m=-1,1        ! 27 clusters for R by adding and
         do p=-1,1        ! subtracting a1, a2 and a3
            iv2=iv1 + n*a(1,:)+m*a(2,:)+p*a(3,:)
            rv2=matmul(c1,iv2)
            if(rv2(1).lt.oner-epsilon.and.rv2(1).gt.zeror-epsilon.and.
     &         rv2(2).lt.oner-epsilon.and.rv2(2).gt.zeror-epsilon.and.
     &         rv2(3).lt.oner-epsilon.and.rv2(3).gt.zeror-epsilon ) then
c              This point is within the cluster
c              Now we must find which point it is.
               do ic=1,Nc
                  if(sum(abs(iv2(:)-Rc(:,ic))).eq.0) then
                     icrdiff(i,j)=ic
                     itest=itest+1   ! complain if multiple solutions are found
                  end if
               end do
               if(iprint.ge.5.and.myrank.eq.0) 
     &              write(lud,"(2x,i2,4x,i2,6x,i2)") i,j,icrdiff(i,j) !0
            end if
         end do
         end do
         end do
         if(itest.ne.1) then
            write(lud,*) 'no icrdiff(i,j) found',i,j
            stop
         end if
      end do
      end do  
        
c***********************************************************************        
c     Now create the tables for the neighbors of cluster points
c***********************************************************************               
c     Note we use c1 from above, don't overwrite it!!  
c

      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) '  ' !0
        write(lud,*) 'neighbors of ic' !0
        write(lud,*) '  ic   neigh  nneigh(ic,jc)' !0
      end if
      neighbor(1,:)=(/ 1, 0, 0/) ! neighbor in the +x direction
      neighbor(2,:)=(/-1, 0, 0/) ! neighbor in the -x direction
      neighbor(3,:)=(/ 0, 1, 0/) ! neighbor in the +y direction
      neighbor(4,:)=(/ 0,-1, 0/) ! neighbor in the yy direction
      neighbor(5,:)=(/ 0, 0, 1/) ! neighbor in the +z direction
      neighbor(6,:)=(/ 0, 0,-1/) ! neighbor in the -z direction
      do i=1,Nc
      do j=1,2*ndim
         iv1(:)=Rc(:,i)-neighbor(j,:)
         itest=0
         do n=-1,1        ! Look in the surrounding 
         do m=-1,1        ! 27 clusters for R by adding and
         do p=-1,1        ! subtracting a1, a2 and a3
            iv2=iv1 + n*a(1,:)+m*a(2,:)+p*a(3,:)
            rv2=matmul(c1,iv2)
            if(rv2(1).lt.oner-epsilon.and.rv2(1).gt.zeror-epsilon.and.
     &         rv2(2).lt.oner-epsilon.and.rv2(2).gt.zeror-epsilon.and.
     &         rv2(3).lt.oner-epsilon.and.rv2(3).gt.zeror-epsilon ) then
c              This point is within the cluster
c              Now we must find which point it is.
               do ic=1,Nc
                  if(sum(abs(iv2(:)-Rc(:,ic))).eq.0) then
                     nneigh1(i,j)=ic
                     itest=itest+1   ! complain if multiple solutions are found
                  end if
               end do
               if(iprint.ge.19.and.myrank.eq.0) 
     &              write(lud,"(2x,i2,4x,i2,6x,i2)") i,j,nneigh1(i,j) !0
            end if
         end do
         end do
         end do

      end do
      end do  
      
!      if(iprint.ge.19.and.myrank.eq.0) then
!        write(lud,*) 'nearest neighbors of icR0'
!        do j=1,2*ndim
 !         write(lud,"(i2)") nneigh1(icR0,j)
 !       end do
  !    end if

ccccccccccccccccccccccccccccccccccccccccccccccc
c	Test the imperfection

      if(iprint.ge.19.and.myrank.eq.0) then
        write(42,*) '  ' !0
        write(42,*) 'neighbors of ic' !0
        write(42,*) 
     &    'shell  nneigh(shell)  nneighl(shell) nneighp(shell)' !0
      end if

      Nm=abs(determinant(a))

      i1 = 2*(Nm+1)**(1d0/ndim); rv1=0d0; nneigh=0; ishellmax=0
      nneighl(:)=0 
      nneighp(:)=0
      
      do i=-i1,i1
      do j=-i1,i1
      do k=-i1,i1
        iv1=(/ i, j, k/)
        i2=sum(abs(iv1)) ! the lattice neighbor index (according to Betts)
        if(i2.gt.0) then
          if(i2.le.nshellmax) then
            nneighl(i2)=nneighl(i2)+1
          end if
          do n=1,ndim
            rv1(n)=iv1(n) + 0.5d0*sum(a(:,n))
          end do
          rv2=matmul(c1,rv1)
          if(rv2(1).lt.1d0-epsilon.and.rv2(1).gt.0d0-epsilon.and.
     &       rv2(2).lt.1d0-epsilon.and.rv2(2).gt.0d0-epsilon.and.
     &       rv2(3).lt.1d0-epsilon.and.rv2(3).gt.0d0-epsilon ) then
c            This point is within the cluster. Now we need to find
c            the shortest distance between this point and icR0, and
c            the corresponding cluster neighbor index i3.
             i3=i2
             do n=-1,1
             do m=-1,1
             do p=-1,1
               iv2=iv1 + n*a(1,:) + m*a(2,:) + p*a(3,:)
               if(sum(abs(iv2)).lt.i2) i3=sum(abs(iv2))
             end do !p
             end do !m
             end do !n
             nneigh(i3)=nneigh(i3)+1
             if(i3.gt.ishellmax) ishellmax=i3
          end if
        end if
      end do !k
      end do !j
      end do !i
      
      i2=0
      do while (sum(nneighp).lt.Nc-1)
        i2=i2+1
        nneighp(i2)=nneighl(i2)
      nneighp(i2)=Nc-1-sum(nneighp(1:i2-1))
      end do
      
      if(sum(nneigh)+1.ne.Nc.or.sum(nneighp)+1.ne.Nc) then
        if(myrank.eq.0) write(lud,*) 'error in nneigh'
        if(myrank.eq.0) write(lud,*) 'sum(nneigh)=',sum(nneigh)
        if(myrank.eq.0) write(lud,*) 'sum(nneighp)=',sum(nneighp)
        stop
      endif
        
      if(iprint.ge.19.and.myrank.eq.0) then
        do i1=1,ishellmax
          write(lud,"(2x,i3,8x,i3,10x,i3,10x,i3)") 
     &      i1,nneigh(i1),nneighl(i1),nneighp(i1)
        enddo
        write(lud,*) 'ishellmax=',ishellmax
      endif
c     Now determine the perfection.  It seems that Betts defines the
c     ferromagnetic imperfection I_F as the numer of additional spins in
c     the shells beyond j plus the number missing the shells before j,
c     and not counting those in the jth shell.  The integer j seems to
c     be the number which yeilds the smallest I_F

      do j=1,ishellmax
        i3=0                                 ! count of imperfections
        do i1=1,j-1
          i3=i3+abs(nneigh(i1)-nneighl(i1))  ! number in shell i1 for cluster-lattice
        end do
        do i1=j+1,ishellmax                  ! cluster occupancy beyond j is an imperfection
          i3=i3+nneigh(i1)
        end do
        if(j.gt.1) then
          if(i3.lt.I_F) then
             I_F=i3 
             I_FShell=j
           end if
        else  ! j=1
          i_F=i3
          i_FShell=j
        end if
      end do

!	if(iprint.ge.19.and.myrank.eq.0)
 !    &  write(lud,*) 'Ferromagnetic imperfection I_F=',I_F
     
c     Betts proposes a simpler algorithm to calculate the perfection,
c     involving the difference of two quantities.
c
c     F_f = SUM k nneigh(k)
c            k
c
c     F_p = SUM k nneighp(k)
c            k

      i2=0; i3=0
      do i1=1,ishellmax
        i2=i2+i1*nneigh(i1)
        i3=i3+i1*nneighp(i1)
      end do
	if(iprint.ge.19.and.myrank.eq.0)
     &  write(lud,*) 'Betts imperfection I_F=',i2-i3
     
     
      

c***********************************************************************        
c       Now set up the K-space tables
c***********************************************************************        

c       Now calculate the principle translations in K-space
c   
c  
c                              _    _
c    _     _    _        2 pi aj x ak          
c    a --> g    g_i =  ------------------      plus cyclic permutations
c                      | ai . (aj x ak)|       
c                                           
c       Principle translation vectors for lattice tiling       
c                                        
c
      r1=determinant(a)
      
      gvector(1,1)=2*pi*( a(2,2)*a(3,3)-a(3,2)*a(2,3) )/r1  !a2y*a3z-a3y*a2z
      gvector(1,2)=2*pi*( a(2,3)*a(3,1)-a(3,3)*a(2,1) )/r1  !a2z*a3x-a3z*a2x
      gvector(1,3)=2*pi*( a(2,1)*a(3,2)-a(3,1)*a(2,2) )/r1  !a2x*a3y-a3x*a2y
      
      gvector(2,1)=2*pi*( a(3,2)*a(1,3)-a(1,2)*a(3,3) )/r1  !a3y*a1z-a1y*a3z
      gvector(2,2)=2*pi*( a(3,3)*a(1,1)-a(1,3)*a(3,1) )/r1  !a3z*a1x-a1z*a3x
      gvector(2,3)=2*pi*( a(3,1)*a(1,2)-a(1,1)*a(3,2) )/r1  !a3x*a1y-a1x*a3y
        
      gvector(3,1)=2*pi*( a(1,2)*a(2,3)-a(2,2)*a(1,3) )/r1  !a1y*a2z-a2y*a1z
      gvector(3,2)=2*pi*( a(1,3)*a(2,1)-a(2,3)*a(1,1) )/r1  !a1z*a2x-a2z*a1x
      gvector(3,3)=2*pi*( a(1,1)*a(2,2)-a(2,1)*a(1,2) )/r1  !a1x*a2y-a2x*a1y
      
      r1=abs(determinant(gvector))
cc      r1=determinant(gvector)

      if(abs((2*pi)**3/r1-Nc).gt.epsilon) then
        if(myrank.eq.0) write(lud,*) 'This cell volume aint Nc',r1
        stop
      end if
      do i=1,3
      do j=1,3
        r1=sum(a(i,:)*gvector(j,:))/(2*pi)
        if(i.eq.j.and.abs(r1-oner).gt.epsilon) then
          if(myrank.eq.0) write(lud,*) 'bug in g',i,j,r1
          stop
        end if
        if(i.ne.j.and.abs(r1).gt.epsilon) then
          if(myrank.eq.0) write(lud,*) 'bug in g',i,j,r1
          stop
        end if
      end do
      end do
      
      if(iprint.ge.19.and.myrank.eq.0) 
     &  write(lud,"(' gvector=',3(' (',f5.2,1x,f5.2,1x,f5.2,')'))") 
     &  (gvector(l,:),l=1,3)
     

c***********************************************************************        
c     Create the cluster momentum tables Kc(i,ic), i=x,y,z
c***********************************************************************        

c     First find the points in the irreducible wedge.  They will be
c     indexed from 1 to Ncw
      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) '  Wedge K points  '    !0
        write(lud,*) 'ic  i  j  k  Kcx(ic)  Kcy(ic)  Kcz(ic)  ' !0
      end if
      ic=0; Ncw=0; icK0=0; icKpi=0
      do i=-Nc,Nc    !times a1
      do j=-Nc,Nc    !times a2
      do k=-Nc,Nc    !times a3
        iv1=(/i,j,k/)
        do l=1,3         !x,y,z components 
          rv1(l)=sum(iv1(:)*gvector(:,l))
        end do
        if(rv1(1).gt.-pi+epsilon.and.rv1(1).lt.pi+epsilon.and.
     &     rv1(2).gt.-pi+epsilon.and.rv1(2).lt.pi+epsilon.and.
     &     rv1(3).gt.-pi+epsilon.and.rv1(3).lt.pi+epsilon) then
c          The point is in the first SC Brillouin zone.  Now see if
c          the point falls within the IRW.  First see if it is 
c          equivalent to any known point in the IRW (if so it cannot be in IRW)
          Kcops: do igroup=2,ngroup ! exclude the identity
            rv2=transform3D(rv1,group(igroup),ixy,ixz,iyz,isx,isy,isz)
            itest=0
            do jc=1,Ncw
              if(sum(abs(Kc(:,jc)-rv2(:))).lt.epsilon) then ! rv2 = Kc(jc)
                 exit Kcops
              end if
            end do
          end do  Kcops
          if(igroup.gt.ngroup) then ! rv1 is in 1BZ and not equivalent to any point  
             Ncw=Ncw+1        ! in the known IRW.  It must be in the IRW
             ic=ic+1
             Kc(:,ic)=rv1(:)
             ig(ic)=i
             jg(ic)=j
             kg(ic)=k
             if(sum(abs(rv1-pi)).lt.epsilon) icKpi=ic
             if(sum(abs(rv1-zeror)).lt.epsilon) icK0=ic
             if(iprint.ge.19.and.myrank.eq.0) 
     &            write(lud,5) ic,i,j,k,Kc(1,ic),Kc(2,ic),Kc(3,ic) !0
          end if
        end if 
      end do 
      end do 
      end do 
	
!      if(icK0.eq.0) then  ! we failed to find (0,0,0)
!        if(myrank.eq.0) 
!     &    write(lud,*) 'failed to find (0,0) in the IRW'    !0
!        stop
!      else if(icKpi.eq.0) then  ! we failed to find (pi,pi,pi)
!        icKpi=icK0
!        if(myrank.eq.0 .and. Nc .ne. 1) then ! To make sure only bipartite cells are used. Exception for Nc = 1
!          write(6,*) '*****************************************'    !0
!          write(6,*) '* failed to find (pi,pi,pi) in the IRW  *'    !0
!          write(6,*) '* apparently the lattice isnt bipartide *'    !0
!	  write(6,*) '* You must change to a Biparite lattice *'
!	  write(6,*) '*******************************************'
!	stop
!       end if
!     end if

c     Now find the points in the 1BZ but outside the irreducible wedge.  
c     They will be indexed from Ncw+1 to Nc
      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) ' Non Wedge K points '    !0
      end if
      do i=-Nc,Nc              
      do j=-Nc,Nc
      do k=-Nc,Nc
         iv1=(/i,j,k/)
         do l=1,3         !x,y,z components 
            rv1(l)=sum(iv1(:)*gvector(:,l)) ! K point in xyz
         end do
         if(rv1(1).gt.-pi+epsilon.and.rv1(1).lt.pi+epsilon.and.
     &      rv1(2).gt.-pi+epsilon.and.rv1(2).lt.pi+epsilon.and.
     &      rv1(3).gt.-pi+epsilon.and.rv1(3).lt.pi+epsilon) then
c           The point is in the first SC Brillouin zone.  Now see if
c           the point falls within the IRW.  
            do jc=1,Ncw
               if(sqrt(sum((Kc(:,jc)-rv1(:))**2)).lt.epsilon) then
c                 rv1 and Kc(ic) are the same point
                  exit
               end if
            end do
            if(jc.gt.Ncw) then !rv1 is in 1BZ and not in the IRW.
               ic=ic+1
               Kc(:,ic)=rv1(:)
               ig(ic)=i
               jg(ic)=j
               kg(ic)=k
               if(iprint.ge.19.and.myrank.eq.0) 
     &              write(lud,5) ic,i,j,k,Kc(1,ic),Kc(2,ic),Kc(3,ic) !0
            end if
         end if 
      end do 
      end do 
      end do 
        
 5        format(1x,i2,1x,i2,1x,i2,1x,i2,2x,f6.3,3x,f6.3,3x,f6.3)
      
      ntw=nl*Ncw      
      
      if(ic.ne.Nc) then
         if(myrank.eq.0) write(lud,*) 'Bug in Kc(ic) loop' !0
         stop
      end if
      if(iprint.ge.19.and.myrank.eq.0) 
     &  write(lud,"('icK0=',i2,' icKpi=',i2,' Ncw=',i2,' Nc=',i2)") 
     &  icK0,icKpi,Ncw,ic

c***********************************************************************
c     Create the indirect addressing table ict(i,j,k) which maps a 
c     point indexed by the principle translations in k-space to ic.
c***********************************************************************        
      N_sp=sqrt(float(Nc))
      do i=-N_sp,N_sp
      do j=-N_sp,N_sp
      do k=-N_sp,N_sp
         iv1=(/i,j,k/)
         do l=1,ndim        !x,y,z components 
            rv1(l)=sum(iv1(:)*gvector(:,l)) ! K point in xyz
         end do
! Map rv1 into -pi->+pi BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+epsilon) rv1=rv1-twor*pi
         where (rv1.lt.-pi+epsilon) rv1=rv1+twor*pi

         do jc=1,Nc
            if(sum(abs(rv1(:)-Kc(:,jc))).lt.epsilon) then    
               ict(i,j,k)=jc
               exit
            end if
         end do
         if(jc.gt.Nc) then
            write(lud,*) 'ict failed for i,j,k=',i,j,k
            write(lud,*) rv1
            stop
         end if
      end do
      end do
      end do

c***********************************************************************        
c     Create the table ickdiff(ic1,ic2), which takes the difference
c     between K-points
c***********************************************************************        
      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) '   '                           !0
        write(lud,*) 'ic1   ic2  ickdiff(ic1,ic2)  ' !0
      end if
      do ic1=1,Nc
      do ic2=1,Nc
         rv1=Kc(:,ic1)-Kc(:,ic2)
! Map rv1 into -pi->+pi BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+epsilon) rv1=rv1-twor*pi
         where (rv1.lt.-pi+epsilon) rv1=rv1+twor*pi
c        Now we need to figure out where this point is!
         do ic=1,Nc
            if(sum(abs(rv1(:)-Kc(:,ic))).lt.epsilon) then    
               ickdiff(ic1,ic2)=ic
               exit
            end if
         end do
         if(ic.gt.Nc) then
            write(lud,5) 'ickdiff failed for ic1,ic2=',ic1,ic2
            stop
         end if
         if(iprint.ge.19.and.myrank.eq.0) 
     &        write(lud,"(1x,i2,4x,i2,6x,i2)") ic1,ic2,ickdiff(ic1,ic2) !0
      end do
      end do
          
c***********************************************************************        
c     Create the table ickplus(ic1,ic2), which takes the sum
c     between to K-points
c***********************************************************************        
      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) '   '                           !0
        write(lud,*) 'ic1   ic2  ickplus(ic1,ic2)  ' !0
      end if
      do ic1=1,Nc
      do ic2=1,Nc
         rv1=Kc(:,ic1)+Kc(:,ic2)
! Map rv1 into -pi->+pi BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+epsilon) rv1=rv1-twor*pi
         where (rv1.lt.-pi+epsilon) rv1=rv1+twor*pi
c        Now we need to figure out where this point is!
         do ic=1,Nc
            if(sum(abs(rv1(:)-Kc(:,ic))).lt.epsilon) then    
               ickplus(ic1,ic2)=ic
               exit
            end if
         end do
         if(ic.gt.Nc) then
            write(lud,5) 'ickplus failed for ic1,ic2=',ic1,ic2
            stop
         end if
         if(iprint.ge.19.and.myrank.eq.0) 
     &        write(lud,"(1x,i2,4x,i2,6x,i2)") ic1,ic2,ickplus(ic1,ic2) !0
      end do
      end do

c***********************************************************************
c     Now find the points equivalent to ic, ickequ(ic,igroup) where 
c     igroup is the group operation.
c***********************************************************************        
c
c
       if(iprint.ge.19.and.myrank.eq.0) write(lud,*) '   '               !0
        if(iprint.ge.19.and.myrank.eq.0) write(lud,*) 'equivalency in k'  !0
        if(iprint.ge.19.and.myrank.eq.0) 
     &              write(lud,*) '   ic     igroup   ickequ(ic,igroup)'  !0

      do ic=1,Nc
         i=ig(ic)
         j=jg(ic)
         k=kg(ic)
         iv1=(/i,j,k/)
         do l=1,ndim
           rv2(l)=sum(iv1(:)*gvector(:,l))
         end do
c        Kc(ic) falls within the zone, by construction, now look 
c        for equivalent points
c     
         do igroup=1,ngroup
            rv1=transform3D(rv2,group(igroup),ixy,ixz,iyz,isx,isy,isz) 
c           These are the equivalent points.  Now map to 1st BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+epsilon) rv1=rv1-twor*pi
         where (rv1.lt.-pi+epsilon) rv1=rv1+twor*pi
c           Now we need to figure out where this point is.
            do jc=1,Nc
               if(sqrt(sum((rv1(:)-Kc(:,jc))**2)).lt.epsilon) then
                  ickequ(ic,igroup)=jc
                  exit
               end if
            end do
            if(iprint.ge.5.and.myrank.eq.0) 
     &        write(lud,"(2x,i2,4x,i2,6x,i2)") 
     &        ic,igroup,ickequ(ic,igroup)                  !0
            if(jc.gt.Nc) then
               if(myrank.eq.0) then
                 write(lud,*)'no equivalent K-point found' !0
                 write(lud,"('ic,igroup,group(igroup)',3(2x,i2))")
     &             ic,igroup,group(igroup)         !0
                 write(lud,"('K =',3(f6.3,1x))") rv2        !0
                 write(lud,"('K1=',3(f6.3,1x))") rv1       !0
               end if
               stop
            end if
         end do
      end do


c***********************************************************************        
c       Now generate a table which maps any K to an equivalent point
c       in the irreducible wedge (IW), and a table for the degeneracy 
c       of each point in the IW. 
c 
c       The use of the wedge can greatly reduce the storage required.
c       for example for Nc=8, only 4 of the 8 points are in the irreducible 
c       wedge.  Thus, all storage arrays may be reduced in size by a 
c       factor of 2.  As Nc-->oo, the saving factor approaches 48 for
c       clusters of high symmetry.  
c
c***********************************************************************        
c       
      
      ickdeg=0; irwgroup=0
      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) '   '                            !0
        write(lud,*) ' degeneracy in k'              !0
        write(lud,*) ' ic  ickmap(ic)  ickdeg(ic) irwgroup(ic)' !0
      end if
c     scan over all of the K points in the 1BZ
      do ic=1,Nc
         i=ig(ic)
         j=jg(ic)
         k=kg(ic)
         iv1=(/i,j,k/)
         do l=1,ndim
            rv2(l)=sum(iv1(:)*gvector(:,l))
         end do
c        Scan over all allowed point group operations, to find the
c        one that maps ic into the IRW
         groupops: do igroup=1,ngroup
           rv1=transform3D(rv2,group(igroup),ixy,ixz,iyz,isx,isy,isz) 
c          These are the equivalent points.  Now map to 1st BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+epsilon) rv1=rv1-twor*pi
         where (rv1.lt.-pi+epsilon) rv1=rv1+twor*pi
c          See if this point is in the IRW!
           do jc=1,Ncw
              if(sum(abs(rv1(:)-Kc(:,jc))).lt.epsilon) then
                 ickmap(ic)=jc
                 ickdeg(jc)=ickdeg(jc)+1
                 irwgroup(ic)=igroup
                 exit groupops
              end if
           end do
         end do groupops
         if(igroup.gt.ngroup) then
            if(myrank.eq.0) write(lud,*) 'mapping failed for' !0
            if(myrank.eq.0) write(lud,*) 'ic',ic !0
            if(myrank.eq.0) write(lud,"('K=',3(f6.3,1x))") rv2 !0
            stop
         end if
         
      end do !ic

      if(iprint.ge.19.and.myrank.eq.0) then 
        do ic=1,Nc    
          write(lud,"(2x,i2,6x,i2,9x,i2,9x,i2)") 
     &    ic,ickmap(ic),ickdeg(ic),irwgroup(ic)
        end do
      end if

c     Test the degeneracy table.  Since ickdeg holds the degeneracy
c     of points in the IRW, its sum must equal Nc.
      if(sum(ickdeg(1:Ncw)).ne.Nc) then
         if(myrank.eq.0) then
           write(lud,*) 'bug in ickdeg',sum(ickdeg(1:Ncw))
           write(lud,*) 'ickdeg=',ickdeg(1:Ncw)
         end if
         stop
      end if
        
c     Test the group table.  irwgroup(ic within wedge)=1, the identity
      do ic=1,Ncw
         if(irwgroup(ic).ne.1) then
            if(myrank.eq.0) then
              write(lud,*) 'irwgroup failed for wedge pt.'     !0
              write(lud,*) 'ic',ic,'irwgroup(ic)',irwgroup(ic) !0
              write(lud,"('K=(',3(f6.3,1x),')')") Kc(:,ic)      !0
            end if
            stop
         end if
      end do
     
c***********************************************************************
c     generate the kt-points in the Brillouin-zone of the 
c     super-lattice, i.e. in the Wigner-Seitz cell of the  
c     reciprocal space defined by the cluster K-points
c***********************************************************************

      ntot = nover
      Dkt=oner/real(ntot,kind)
      nwsc=0 ! number of k-points in Wigner-Seitz cell 
      if (ntot.gt.0) then
      do i=-ntot+1,ntot
        do j=-ntot+1,ntot
          do k=1,ntot
           iv1=(/i,j,k/)
           rv1(:)=Dkt*(iv1(:)-halfr)*pi
c  Check if this k-point is in Wigner-Seitz cell, i.e. if the closest
c  K-point is K=0 
            r1=sum(rv1(:)*rv1(:)) ! distance from K=0
            itest=0
            do jc=1,Nc
            if (jc.eq.icK0) cycle
             r2=rv1(1)-Kc(1,jc)
             r3=rv1(2)-Kc(2,jc)
             r4=rv1(3)-Kc(3,jc)
             if((r2-pi).gt.zeror) then; r2=r2-twor*pi; endif
             if((r2+pi).lt.zeror) then; r2=r2+twor*pi; endif
             if((r3-pi).gt.zeror) then; r3=r3-twor*pi; endif
             if((r3+pi).lt.zeror) then; r3=r3+twor*pi; endif
             if((r4-pi).gt.zeror) then; r4=r4-twor*pi; endif
             if((r4+pi).lt.zeror) then; r4=r4+twor*pi; endif
             r2=(r2*r2+r3*r3+r4*r4)
             if(r2.lt.r1) then
               itest=1
               exit
             endif
c
            enddo ! jc
            if(itest.eq.0) then 
             nwsc=nwsc+1
             kt(1:3,nwsc) = rv1(:)
             nwsc=nwsc+1 ! -k must also be part of the WS cell
             kt(1:3,nwsc) = -rv1(:)
            endif
          enddo ! k
        enddo ! j
       enddo ! i
       else ! needed for FSS calculation
          nwsc=1
          kt(:,1)=zeror
       endif
c
       if(myrank.eq.0.and.iprint.ge.19) then
          write(lud,*) 'Number of points in WS, nwsc', nwsc
          write(lud,*) 'Nc times nwsc', Nc*nwsc
          write(lud,*) 'Total number of points in BZ', (2*ntot)**3 
       endif

c***********************************************************************
c     Form Epsd(K) for calculating the partial density of states
c***********************************************************************     
      call allocate_PDOS
      estep = (eend-estart)/dble(e_div)
      PDOS(:,:) = 0.d0
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) '    '
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) 
     &  ' ic      Epsbar(ic)        temp_min      temp_max    min_Eps(ic)    max_Eps(ic)'
	Epsbar=zeror
      do ic=1,Nc
	temp_max = - halfr*(cos(Kc(1,ic)+kt(1,1))+
     &			cos(Kc(2,ic)+kt(2,1))+cos(Kc(3,ic)+kt(3,1))) 
     &             - tprime*(cos(Kc(1,ic)+kt(1,1))*cos(Kc(2,ic)+kt(2,1))
     &                      +cos(Kc(2,ic)+kt(2,1))*cos(Kc(3,ic)+kt(3,1))
     &                      +cos(Kc(1,ic)+kt(1,1))*cos(Kc(3,ic)+kt(3,1))-oner)
	temp_min = temp_max
         do i=1,nwsc
            rv1(1) = Kc(1,ic)+kt(1,i)
            rv1(2) = Kc(2,ic)+kt(2,i)
            rv1(3) = Kc(3,ic)+kt(3,i)
            Epsbar(ic)= Epsbar(ic)
c     &            -   halfr*sum(cos(rv1)) 
     &      -   halfr*(cos(rv1(1))+cos(rv1(2))+cos(rv1(3))) 
     &            - tprime*(cos(rv1(1))*cos(rv1(2))
     &                     +cos(rv1(2))*cos(rv1(3))
     &                     +cos(rv1(1))*cos(rv1(3))-oner)
c            Epsd = - halfr*sum(cos(rv1))
             Epsd = - halfr*(cos(rv1(1))+cos(rv1(2))+cos(rv1(3)))
     &            - tprime*(cos(rv1(1))*cos(rv1(2))
     &                     +cos(rv1(2))*cos(rv1(3))
     &                     +cos(rv1(1))*cos(rv1(3))-oner)
		e_index = NINT((Epsd-estart)/estep)
		if (e_index >= 0 .and. e_index <= e_div) then
		    PDOS(e_index,ic) = PDOS(e_index,ic) + 1.d0
		else
		    write(*,*) "WARNING: ENERGY EXCEEDS : ", Epsd, estart, eend
		endif
	  if (Epsd.gt.temp_max) temp_max = Epsd
	  if (Epsd.lt.temp_min) temp_min = Epsd   
         end do ! i
        Epsbar(ic)=Epsbar(ic)/real(nwsc,kind)
	dEps(ic) = (temp_max - temp_min)/real(nover,kind) ! the partial DOS also has nover ticks
	min_Eps(ic) = temp_min
	max_Eps(ic) = temp_max
       if(myrank.eq.0.and.iprint.ge.19) 
     &	write(lud,"(i6,4x,f12.8,4x,f12.8,4x,f12.8,4x,f12.8,4x,f12.8)") 
     &     ic, Epsbar(ic), temp_min, temp_max,min_Eps(ic), max_Eps(ic)
      	end do !ic

c     Check 
c      r1=sum(Epsbar); r2=sum(Epsbar*ickdeg)
c      if(abs(tprime).lt.epsilon.and.abs(r1).gt.epsilon
c     &   .or.
c     &	 abs(tprime).lt.epsilon.and.abs(r2).gt.epsilon) then
c        if(myrank.eq.0.and.iprint.ge.19) then
c          write(lud,*) 'sum Epsbar wrong.',r1,r2
c	end if
c        stop
c      end if
c***********************************************************************
c     Calculate the partial density of states
c
c                     Nc            
c       p_rho(K,w) =  --   SUM_kt delta( eps(K+kt) - w )
c                     N
c                       
c***********************************************************************  

      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) '    '
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) 
     &  '  ic         ener(ic)         sum(ener) n_fact' !n_fact is the number of points within each cell ic
	
	n_fact = zeror
	do ic=1,Nc
	do i=0,e_div
	n_fact(ic) = n_fact(ic) + PDOS(i,ic)*estep
	end do
	end do

	do ic=1,Nc
	do i=0,e_div
	PDOS(i,ic) = PDOS(i,ic)/(n_fact(ic)*Nc)
	write(53,*) estart+i*estep, PDOS(i,ic) 
	end do
	write(53,*) '   '
	end do

	do i=0,e_div
        write(53,*) estart+i*estep, sum(PDOS(i,:))
        end do
	close(53)

 	  write(4005,403) 
403      format('#    count         pdos')
	ener = zeror
	do ic =1,Nc
	do i=0,e_div
	write(4005,*) estart+i*estep,PDOS(i,ic)
		      ener(ic) = ener(ic) + PDOS(i,ic)*estep
	end do
c	write(4006,*) ic, ener(ic), sum(ener)
       if(myrank.eq.0.and.iprint.ge.19) 
     &	write(lud,"(i6,4x,f12.8,4x,f12.8,4x,f18.8)") 
     &   ic, ener(ic),  sum(ener), n_fact(ic)
	end do

      p_rho_GcKf = zeror
      do ic=1,Nc
        do n=-nwn+1,nwn-1
	  do i=0,e_div-1
	    if ((wn(n) - (estart+float(i)*estep))*(wn(n) - (estart+float(i-1)*estep)).le.0) then
		p_rho_GcKf(ic,n) = halfr*(PDOS(i,ic) + PDOS(i+1,ic))
	    end if
	  end do
	end do
      end do


c***********************************************************************
c     Form Epsbar_histo(K) for calculating the Histogram
c***********************************************************************     
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) '    '
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) 
     &  ' ic,Epsbar_histo(ic),Epsbar(ic),wn(Epsbar_histo(ic))'
c	Epsbar=zeror
c	do ic=1,Nc
c         do i=1,nwsc
c            rv1(1) = Kc(1,ic)+kt(1,i)
c            rv1(2) = Kc(2,ic)+kt(2,i)
c            rv1(3) = Kc(3,ic)+kt(3,i)
c            Epsbar(ic)= Epsbar(ic)
c     &      -   halfr*(cos(rv1(1))+cos(rv1(2))+cos(rv1(3))) 
c     &            - tprime*(cos(rv1(1))*cos(rv1(2))
c     &                     +cos(rv1(2))*cos(rv1(3))
c     &                     +cos(rv1(1))*cos(rv1(3))-oner)
c                         end do ! i
c                      Epsbar(ic)=Epsbar(ic)/real(nwsc,kind)
c                      do n=-nwn,nwn ! Frequency loop starts
c	           
c	       if (wn(n) - Epsbar(ic).lt. epsilon) Epsbar_histo(ic) = n
c          end do !n
c        if(myrank.eq.0.and.iprint.ge.19) 
c     &	write(lud,"(i6,4x,i8,4x,f14.8,4x,f14.8)") 
c     &     ic,Epsbar_histo(ic),Epsbar(ic),wn(Epsbar_histo(ic))
c	end do !ic

!   finding when wn=eps_bar(ic)
!  1.  define function to find min of:
          funct(:,:)=zeror
      	  do ic=1,Nc
	    do n=-nwn,nwn
              funct(ic,n)=wn(n)-Epsbar(ic)
    	    end do
    	  end do
!  2.  find min of this function:  
          
          do ic=1,Nc   
    	     min_funct(ic)=abs(funct(ic,-nwn))
    	     Epsbar_histo(ic)=-nwn 
    	      do n=-nwn,nwn
    	       if(min_funct(ic) .gt. abs(funct(ic,n))) then 
    	        min_funct(ic)=abs(funct(ic,n))
    	        Epsbar_histo(ic)=n      ! index at which wn=eps_bar
    	       end if	
              end do ! n
        if(myrank.eq.0.and.iprint.ge.19) 
     &	write(lud,"(i6,4x,i8,4x,f14.8,4x,f14.8)") 
     &     ic,Epsbar_histo(ic),Epsbar(ic),wn(Epsbar_histo(ic))
           end do 
  	  
    	  



c***********************************************************************        
c       Form the Fourier transform tables.  
c       R-->K,  FTCoefs_K_to_R
c       K-->R,  FTCoefs_K_to_R
c
c                    __          the factor of 1/N is required
c                 1  \           so that both G(R)~1/w and G(K)~1/w
c       G(R=0) = --- /  G(K)     for large w
c                N   --
c                    K
c
c***********************************************************************        
      do ic=1,Nc
         do ik=1,Nc
            r1=zeror
            do idim=1,ndim
               r1=r1+Kc(idim,ik)*Rc(idim,ic)
            end do
            FTCoefs_K_to_R(ik,ic)=exp(-ii*r1)/real(Nc,kind)
            FTCoefs_R_to_K(ik,ic)=exp(+ii*r1)
         end do
      end do

c	do ic = 1,Nc
c	do i = 1,ndim
c	write(235,*) ic,i,Kc(i,ic),Rc(i,ic)
c	end do
c	end do

c     Now test for consistency of Kc and Rc
      do ic=1,Nc
        r1= real(sum(FTCoefs_K_to_R(1:Nc,ic)))-kroniker(ic,icR0)
        r2=aimag(sum(FTCoefs_K_to_R(1:Nc,ic)))
        r3= real(sum(FTCoefs_K_to_R(ic,1:Nc)))-kroniker(ic,icK0)
        r4=aimag(sum(FTCoefs_K_to_R(ic,1:Nc)))
        
        if(abs(r1).gt.epsilon.or.abs(r2).gt.epsilon) then
          write(lud,*) 'Kc and Rc are inconsistent'
          write(lud,*) 'icr,r1,r2=',ic,r1,r2
          stop
        end if
        if(abs(r3).gt.epsilon.or.abs(r4).gt.epsilon) then
          write(lud,*) 'Kc and Rc are inconsistent, K=0'
          write(lud,*) 'ick,r3,r4=',ic,r3,r4
          stop
        end if
      end do
      end

c***********************************************************************        
c***********************************************************************        
c         3D Betts subroutines and functions
c***********************************************************************        
c***********************************************************************

        
      function gtransform3D(a,igroup,ixy,ixz,iyz,isx,isy,isz)
c
c  x <-- isx*[(1-ixy)*x + ixy*y]  <-- [(1-ixz)*x + ixz*z]<----------------------
c  y <-- isy*[(1-ixy)*y + ixy*x]  <----------------------<-- [(1-iyz)*y + iyz*z]
c  z <-- isz*<------------------- <-- [(1-ixz)*z + ixz*x]<-- [(1-iyz)*z + iyz*y]
c
      implicit none
      integer :: i,a(3,3),g(3,3),gtransform3D(3,3),igroup,
     &           ixy(48),ixz(48),iyz(48),isx(48),isy(48),isz(48)
      
      do i=1,3  
         gtransform3D(i,1)=a(i,1)
         gtransform3D(i,2)=(1-iyz(igroup))*a(i,2) + iyz(igroup)*a(i,3)
         gtransform3D(i,3)=(1-iyz(igroup))*a(i,3) + iyz(igroup)*a(i,2)
      end do
               
      do i=1,3  
         g(i,1)= (1-ixz(igroup))*gtransform3D(i,1)
     &         + ixz(igroup)*gtransform3D(i,3)
         g(i,2)=                gtransform3D(i,2)
         g(i,3)= (1-ixz(igroup))*gtransform3D(i,3) 
     &         + ixz(igroup)*gtransform3D(i,1)
      end do

      do i=1,3  
         gtransform3D(i,1)=isx(igroup)*((1-ixy(igroup))*g(i,1) 
     &                  + ixy(igroup)*g(i,2))
         gtransform3D(i,2)=isy(igroup)*((1-ixy(igroup))*g(i,2) 
     &                  + ixy(igroup)*g(i,1))
         gtransform3D(i,3)=isz(igroup)*                 g(i,3)         
      end do
      
      return
      end
        
      function Rtransform3D(Rc1,igroup,ixy,ixz,iyz,isx,isy,isz)
c
c  x <-- isx*[(1-ixy)*x + ixy*y]  <-- [(1-ixz)*x + ixz*z]<----------------------
c  y <-- isy*[(1-ixy)*y + ixy*x]  <----------------------<-- [(1-iyz)*y + iyz*z]
c  z <-- isz*<------------------- <-- [(1-ixz)*z + ixz*x]<-- [(1-iyz)*z + iyz*y]
c
      implicit none
      integer :: Rc1(3),g(3),Rtransform3D(3),igroup, 
     &           ixy(48),ixz(48),iyz(48),isx(48),isy(48),isz(48)
      
      Rtransform3D(1)=            Rc1(1)
      Rtransform3D(2)=(1-iyz(igroup))*Rc1(2) + iyz(igroup)*Rc1(3)
      Rtransform3D(3)=(1-iyz(igroup))*Rc1(3) + iyz(igroup)*Rc1(2)
      
      g(1)=(1-ixz(igroup))*Rtransform3D(1) + ixz(igroup)*Rtransform3D(3)
      g(2)=                Rtransform3D(2)
      g(3)=(1-ixz(igroup))*Rtransform3D(3) + ixz(igroup)*Rtransform3D(1)
      
      Rtransform3D(1)=isx(igroup)*((1-ixy(igroup))*g(1) 
     &             + ixy(igroup)*g(2))
      Rtransform3D(2)=isy(igroup)*((1-ixy(igroup))*g(2) 
     &             + ixy(igroup)*g(1))
      Rtransform3D(3)=isz(igroup)*                 g(3) 
      
      return
      end
      
      function Ktransform3D(Kc1,igroup,ixy,ixz,iyz,isx,isy,isz)
c     
c  x <-- isx*[(1-ixy)*x + ixy*y]  <-- [(1-ixz)*x + ixz*z]<----------------------
c  y <-- isy*[(1-ixy)*y + ixy*x]  <----------------------<-- [(1-iyz)*y + iyz*z]
c  z <-- isz*<------------------- <-- [(1-ixz)*z + ixz*x]<-- [(1-iyz)*z + iyz*y]
c
      implicit none
      integer :: igroup,
     &           ixy(48),ixz(48),iyz(48),isx(48),isy(48),isz(48)
      real(8) :: Kc1(3),g(3),Ktransform3D(3)
      
      Ktransform3D(1)=                Kc1(1)
      Ktransform3D(2)=(1-iyz(igroup))*Kc1(2) + iyz(igroup)*Kc1(3)
      Ktransform3D(3)=(1-iyz(igroup))*Kc1(3) + iyz(igroup)*Kc1(2)
      
      g(1)=(1-ixz(igroup))*Ktransform3D(1) + ixz(igroup)*Ktransform3D(3)
      g(2)=                Ktransform3D(2)
      g(3)=(1-ixz(igroup))*Ktransform3D(3) + ixz(igroup)*Ktransform3D(1)
      
      Ktransform3D(1)=isx(igroup)*((1-ixy(igroup))*g(1) 
     &             + ixy(igroup)*g(2))
      Ktransform3D(2)=isy(igroup)*((1-ixy(igroup))*g(2) 
     &             + ixy(igroup)*g(1))
      Ktransform3D(3)=isz(igroup)*                 g(3) 
      
      return
      end
      
      function Inverse_3x3(a)   
      integer :: a(3,3)
      real(8) :: Inverse_3x3(3,3),r1
c
c       [ a(1,1) a(1,2) a(1,3) ]
c     A=[ a(2,1) a(2,2) a(2,3) ]
c       [ a(3,1) a(3,2) a(3,3) ]
c
c      -1        i+j 
c     A_ji = (-1)   C_ij  /|A|    C_ij is the ij cofactor of A
c     
      
      r1=   a(1,1)*a(2,2)*a(3,3)  ! the determinent of a
     &    + a(1,2)*a(2,3)*a(3,1)
     &    + a(1,3)*a(2,1)*a(3,2)
     &    - a(3,1)*a(2,2)*a(1,3)
     &    - a(3,2)*a(2,3)*a(1,1)
     &    - a(3,3)*a(2,1)*a(1,2) 
      
      Inverse_3x3(1,1) = +( a(2,2)*a(3,3)-a(3,2)*a(2,3) )/r1
      Inverse_3x3(1,2) = -( a(1,2)*a(3,3)-a(3,2)*a(1,3) )/r1
      Inverse_3x3(1,3) = +( a(1,2)*a(2,3)-a(2,2)*a(1,3) )/r1
      
      Inverse_3x3(2,1) = -( a(2,1)*a(3,3)-a(3,1)*a(2,3) )/r1
      Inverse_3x3(2,2) = +( a(1,1)*a(3,3)-a(3,1)*a(1,3) )/r1
      Inverse_3x3(2,3) = -( a(1,1)*a(2,3)-a(2,1)*a(1,3) )/r1
      
      Inverse_3x3(3,1) = +( a(2,1)*a(3,2)-a(3,1)*a(2,2) )/r1
      Inverse_3x3(3,2) = -( a(1,1)*a(3,2)-a(3,1)*a(1,2) )/r1
      Inverse_3x3(3,3) = +( a(1,1)*a(2,2)-a(2,1)*a(1,2) )/r1
      
      
      return
      end
      
      function Ideterminant(b)
      use global
      integer::  Ideterminant,b(3,3)

      Ideterminant=
     &         b(1,1)*b(2,2)*b(3,3)
     &       + b(1,2)*b(2,3)*b(3,1)
     &       + b(1,3)*b(2,1)*b(3,2)
     &       - b(3,1)*b(2,2)*b(1,3)
     &       - b(3,2)*b(2,3)*b(1,1)
     &       - b(3,3)*b(2,1)*b(1,2)
      return
      end
      
      function Rdeterminant(c)
      use global
      real(kind)::  c(3,3),Rdeterminant

      Rdeterminant=
     &         c(1,1)*c(2,2)*c(3,3)
     &       + c(1,2)*c(2,3)*c(3,1)
     &       + c(1,3)*c(2,1)*c(3,2)
     &       - c(3,1)*c(2,2)*c(1,3)
     &       - c(3,2)*c(2,3)*c(1,1)
     &       - c(3,3)*c(2,1)*c(1,2)
      return
      end

      function kroniker(i,j)
      integer i,j,kroniker
      
      if(i.eq.j) then
        kroniker=1
      else
        kroniker=0   
      end if

        end 
c***************************************************************************
        subroutine tables_3D                  !3 DIMENSIONAL VERSION
c***************************************************************************
c       This subroutine makes the lookup tables used throughout the
c       3 DIMENSIONAL code.
c       Note: For now we only implement Lc x Lc x 2 clusters with tetragonal 
c             symmetry. E.g. Nc=32 corresponds to a 4 x 4 x 2 cluster
c***************************************************************************
        use Global
        implicit none
        logical, parameter :: use_original = .false.
c***************************************************************************
      integer, parameter :: maxnc=300
        integer:: i,i1,i2,ix,iy,iodd,ilow,ind,
     &          ic,ic1,ic2,idim,ik,itest,iex(16),isx(16),isy(16),
     &          isz(16),irwgroup(Nc),e_index,
     &		ig(150),jg(150),kg(150),iv1(3),iv2(3),
     &          j,jmd,jm,j1,j2,jc,jk,jkc,jhigh,
     &          k,n,n1,n2,m,m1,m2,l1,l2,l,ic_max,ic_min, p_r
      real(kind) :: funct(Nc,-nwn:nwn),min_funct(Nc)
        real(kind) :: Kx,Ky,Kz,r1,r2,r3,r4,Rx,Ry,Rz,vx,vy,vz,rv1(3),rv2(3),rv3(3),
     &                Dkt,temp_min,temp_max,ener(Nc),n_fact(Nc)


      interface Inverse
      function Inverse_3x3(a)
c      use global
      integer, parameter :: kind=8
      integer :: a(3,3)
      real(kind) :: Inverse_3x3(3,3)
      end function
      end interface
      
      interface determinant
      function Ideterminant(b)
c      use global
      integer, parameter :: kind=8
      integer :: b(3,3),Ideterminant
      end function
      function Rdeterminant(c)
      use global
      real(kind) :: c(3,3),Rdeterminant
      end function
      end interface
c***************************************************************************
c       Local arrays
c       iex(igroup)      D=3 tetragonal (D4h) point group table(see below)
c       isx(igroup)      D=3 tetragonal (D4h) point group table(see below)
c       isy(igroup)      D=3 tetragonal (D4h) point group table(see below)
c       isz(igroup)      D=3 tetragonal (D4h) point group table(see below)
c       irwgroup(ic)     group operation required to reduce K(ic) to IRW     
c***************************************************************************

        
c       CLUSTER MOMENTUM and DUAL SPACE for the 2D square lattice
 
c                            
c    _     _     _           a2 x a3         
c    a --> g     g1 = 2pi  -----------                
c                          a1*(a2 x a3)       
c                                           
c       Principle translation                      [ square   (-1/2,-1/2)
c       vectors for lattice tiling        Origins  | diamond  ( 0  ,-1/2)
c                                                  [ odd      (-epsilon,-epsilon)
c
c              ^               ^                 ^                  
c Nc=1  a1 = 1 x,       a2 = 1 y,         a3 = 1 z

c
c              ^               ^                 ^          
c Nc=8  a1 = 2 x,       a2 = 2 y,         a3 = 2 z
c
c               ^              ^                 ^  
c Nc=32  a1 = 4 x       a2 = 4 y          a3 = 2 z
c
c               ^     ^          ^     ^          ^
c Nc=36  a1 = 3 x + 3 y, a2 = -3 x + 3 y,  a3 = 2 z
c
c  Group Operations for a tetragonal lattice:
c
c  x --> +/- x
c    \/
c    /\
c  y --> +/- y       
c
c  z --> +/- z
c
c  16 symmetry operations in tetragonal D4h group:
c                            
c  E:                  2C4:            C2:            2C2':
c  x -> x       x ->  y   x -> -y     x ->  x    x ->  x   x -> -x
c  y -> y       y -> -x   y ->  x     y -> -y    y -> -y   y ->  y
c  z -> z       z ->  z   z ->  z     z -> -z    z -> -z   z -> -z
c
c        2C2'':          i:          s_h:             2s_nu:    
c  x ->  y   x -> -y     x -> -x     x ->  x     x ->  x   x -> -x
c  y ->  x   y -> -x     y -> -y     y ->  y     y -> -y   y ->  y
c  z -> -z   z -> -z     z -> -z     z -> -z     z ->  z   z ->  z
c   
c        2S4:                 2s_d:
c  x ->  y   x -> -y     x -> y   x -> -y 
c  y -> -x   y ->  x     y -> x   y -> -x
c  z -> -z   z -> -z     z -> z   z ->  z
c
c
c  pgroup(iop,iex,isx,isy) 
c  iex=0(1) dont (do) exchange x and y
c  isx      sign of x 
c  isy      sign of y
c  isz      sign of z
c	
c  x --> isx*[(1-iex)*x + iex*y]
c  y --> isy*[(1-iex)*y + iex*x]
c  z --> isz*z
c
        iex(1)= 0; isx(1)=  1; isy(1)=  1; isz(1)=  1
        iex(2)= 1; isx(2)=  1; isy(2)= -1; isz(2)=  1
        iex(3)= 1; isx(3)= -1; isy(3)=  1; isz(3)=  1
        iex(4)= 0; isx(4)=  1; isy(4)= -1; isz(4)= -1
        iex(5)= 0; isx(5)=  1; isy(5)= -1; isz(5)= -1
        iex(6)= 0; isx(6)= -1; isy(6)=  1; isz(6)= -1
        iex(7)= 1; isx(7)=  1; isy(7)=  1; isz(7)= -1
        iex(8)= 1; isx(8)= -1; isy(8)= -1; isz(8)= -1
        iex(9)= 0; isx(9)= -1; isy(9)= -1; isz(9)= -1 
        iex(10)=0; isx(10)= 1; isy(10)= 1; isz(10)=-1 
        iex(11)=0; isx(11)= 1; isy(11)= 1; isz(11)= 1 
        iex(12)=0; isx(12)=-1; isy(12)= 1; isz(12)= 1 
        iex(13)=1; isx(13)= 1; isy(13)=-1; isz(13)=-1 
        iex(14)=1; isx(14)=-1; isy(14)= 1; isz(14)=-1
        iex(15)=1; isx(15)= 1; isy(15)= 1; isz(15)= 1
        iex(16)=1; isx(16)=-1; isy(16)=-1; isz(16)= 1


c***********************************************************************        
c       determine the principle translation vectors a1,a2 and a3, and 
c       reciprocal space vectors g1, g2 and g3
c***********************************************************************
        a(3,1)=0; a(3,2)=0
        a(1,3)=0; a(2,3)=0
        If(Nc.eq.1) then                              
          a(1,1)=1
          a(1,2)=0
          a(3,3)=1
          Ncw=1
        else if(Nc.eq.8) then                           
          a(1,1)=2                                 
          a(1,2)=0
	  a(3,3)=2 
          
       else if(Nc.eq.18) then                           
          a(1,1)=3                                 
          a(1,2)=0
	  a(3,3)=2 
       else if(Nc.eq.32) then                          
          a(1,1)=4                                 
          a(1,2)=0
	  a(3,3)=2 
       else if (Nc.eq.36) then
	  a(1,1)=3
	  a(1,2)=3
	  a(3,3)=2 
       else if(Nc.eq.50) then                          
          a(1,1)=5                                 
          a(1,2)=0          
       else if(Nc.eq.72) then                          
          a(1,1)=6                                 
          a(1,2)=0 
	  a(3,3)=2       
      else if(Nc.eq.128) then                          
          a(1,1)=8                                 
          a(1,2)=0 
	  a(3,3)=2 
          
      else if(Nc.eq.98) then                          
          a(1,1)=7                                 
          a(1,2)=0 
	  a(3,3)=2 
      else if(Nc.eq.64) then                          
	  a(1,1)=4
	  a(1,2)=4
	  a(3,3)=2 
       else
          if(myrank.eq.0) write(lud,*) 'Bad Nc, bad, bad'              !0
          stop
        end if    

        a(2,1)=-a(1,2)  ! a2 is perpendicular to a1 and a3
        a(2,2)= a(1,1)

c***********************************************************************        
c       Now set up the K-space tables
c***********************************************************************        

c       Now calculate the principle translations in K-space
c   
c  
c                              _    _
c    _     _    _        2 pi aj x ak          
c    a --> g    g_i =  ------------------      plus cyclic permutations
c                      | ai . (aj x ak)|       
c                                           
c       Principle translation vectors for lattice tiling       
c                                        
c

c       Now calculate the principle translations in K-space
!        gvector(1,1)=twor*pi*a(1,1)/real(a(1,1)*
!     &		a(2,2)-a(2,1)*a(1,2),kind)
!        gvector(1,2)=twor*pi*a(1,2)/real(a(1,1)*
!     &		a(2,2)-a(2,1)*a(1,2),kind)
!        gvector(1,3)=zeror
!        gvector(2,1)=-gvector(1,2)
!        gvector(2,2)= gvector(1,1)
!        gvector(2,3)=zeror
!        gvector(3,1)=zeror
!        gvector(3,2)=zeror
!        gvector(3,3)=twor*pi/real(a(3,3),kind)

      r1=determinant(a)
      
      gvector(1,1)=2*pi*( a(2,2)*a(3,3)-a(3,2)*a(2,3) )/r1  !a2y*a3z-a3y*a2z
      gvector(1,2)=2*pi*( a(2,3)*a(3,1)-a(3,3)*a(2,1) )/r1  !a2z*a3x-a3z*a2x
      gvector(1,3)=2*pi*( a(2,1)*a(3,2)-a(3,1)*a(2,2) )/r1  !a2x*a3y-a3x*a2y
      
      gvector(2,1)=2*pi*( a(3,2)*a(1,3)-a(1,2)*a(3,3) )/r1  !a3y*a1z-a1y*a3z
      gvector(2,2)=2*pi*( a(3,3)*a(1,1)-a(1,3)*a(3,1) )/r1  !a3z*a1x-a1z*a3x
      gvector(2,3)=2*pi*( a(3,1)*a(1,2)-a(1,1)*a(3,2) )/r1  !a3x*a1y-a1x*a3y
        
      gvector(3,1)=2*pi*( a(1,2)*a(2,3)-a(2,2)*a(1,3) )/r1  !a1y*a2z-a2y*a1z
      gvector(3,2)=2*pi*( a(1,3)*a(2,1)-a(2,3)*a(1,1) )/r1  !a1z*a2x-a2z*a1x
      gvector(3,3)=2*pi*( a(1,1)*a(2,2)-a(2,1)*a(1,2) )/r1  !a1x*a2y-a2x*a1y
      
      r1=abs(determinant(gvector))
cc      r1=determinant(gvector)

      if(abs((2*pi)**3/r1-Nc).gt.epsilon) then
        if(myrank.eq.0) write(lud,*) 'This cell volume aint Nc',r1
        stop
      end if
      do i=1,3
      do j=1,3
        r1=sum(a(i,:)*gvector(j,:))/(2*pi)
        if(i.eq.j.and.abs(r1-oner).gt.epsilon) then
          if(myrank.eq.0) write(lud,*) 'bug in g',i,j,r1
          stop
        end if
        if(i.ne.j.and.abs(r1).gt.epsilon) then
          if(myrank.eq.0) write(lud,*) 'bug in g',i,j,r1
          stop
        end if
      end do
      end do
      
      if(iprint.ge.19.and.myrank.eq.0) 
     &  write(lud,"(' gvector=',3(' (',f5.2,1x,f5.2,1x,f5.2,')'))") 
     &  (gvector(l,:),l=1,3)


c*************************************************************************
c       Set the location of the origin in real space, denote the symmetry
c       (high, iodd=0; low, iodd=1) and the number of group operations.
c*************************************************************************
        if(a(1,1).eq.a(1,2)) then  ! The xy part is a diamond (Nc=2,8,18,32,50...)
          write(*,* )'diamond',Nc
          origx= zeror
          origy=-halfr
          iodd=0     ! denote a cluster of high symmetry
          ngroup=16   ! with all ngroup=8 symmetry operations
        else if(a(1,1)*a(1,2).eq.0) then  ! The xy part is a square 
c           write(*,* )'square',Nc
          origx=-halfr
          origy=-halfr
          origz=-halfr
          iodd=0     ! denote a cluster of high symmetry
          ngroup=16  ! with all ngroup=16 symmetry operations
        end if

c****************************************************************************       
c       calculate the table of cluster point locations Rc(i,ic), i=1,x or 2,y
c****************************************************************************       
c 
c  ^ y
c  |
cc-|-c---c-*-c---c---c---c               A point P is within the cluster 
c| | |   |   |   |   |   |               if a vector v from the box origin, 
cc-|-c---c---c---c---c---c               O, to the point satisfies
ca2^ * * * * * * * * |   |               
cc-|-c---c---c---c-*-c---c                    _   _
c| | |   P   |   | * |   |                    v . a1
cc-|-c--/c---c---c-*-c---c               0 < -------   < 1
c| | | /v|   |   | * |   |                   |a1.a1|
cc-|-c/--c---c---c-*-c---c 
c| | /   |   |   | * |   |                      and
cc-|/0---c---c---c-*-c---c                    _   _
c| O--------------->---------> x              v . a2
cc---c---c---c---ca1-c---c               0 < -------   < 1
c|   |   |   |   |   |   |                   |a2.a2|
cc---c---c---c---c---c---c
c                                               and
c                                             _   _
c                                             v . a3
c                                        0 < -------   < 1
c                                            |a3.a3|
c
        if(iprint.ge.19.and.myrank.eq.0) write(lud,*) '  '               !0
        if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) 'ic,Rcx(ic),Rcy(ic),Rcz(ic)'                          !0
        ic=0
        do k=-Nc,Nc
        do j=-Nc,Nc
        do i=-Nc,Nc
          vx=real(i,kind)-origx
          vy=real(j,kind)-origy
          vz=real(k,kind)-origz
          r1=(vx*a(1,1)+vy*a(1,2))/real(a(1,1)**2+a(1,2)**2,kind)  ! a(1,3)=0
          r2=(vx*a(2,1)+vy*a(2,2))/real(a(2,1)**2+a(2,2)**2,kind)  ! a(2,3)=0
          r3= vz*a(3,3)/real(a(3,3)**2,kind)                       ! a(3,1)=a(3,2)=0
          if(r1.lt.oner-epsilon.and.r1.gt.zeror-epsilon.and.
     &      r2.lt.oner-epsilon.and.r2.gt.zeror-epsilon.and.
     &      r3.lt.oner-epsilon.and.r3.gt.zeror-epsilon) then
            ic=ic+1
            Rc(1,ic)=i           ! Rc, like i.e. a(1,1) is integer
            Rc(2,ic)=j
            Rc(3,ic)=k
            if(iprint.ge.19.and.myrank.eq.0) 
     &           write(lud,*) ic,Rc(1,ic),Rc(2,ic),Rc(3,ic)
          end if
        end do
        end do
        end do
c        if(ic.ne.Nc.or.Rc(1,1).ne.0.or.Rc(2,1).ne.0
c     &             .or.Rc(3,1).ne.0) then
c          if(myrank.eq.0) write(lud,*) 'Bug in Rc(ic) loop'             !0
c          stop
c        end if
        icR0=1




c***********************************************************************        
c       create a table mapping points outside the cluster back into it.
c***********************************************************************        

        if(iprint.ge.19.and.myrank.eq.0) write(lud,*) '  '               !0
        if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) 'Equivalency table r'                         !0
        if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) '   ic     igroup   icrequ(ic,ieqiv)'         !0
        do ic=1,Nc
        do igroup=1,ngroup
c         Generate the equivalent points using the ngroup operations
c         x --> isx*[(1-iex)*x + iex*y]
c         y --> isy*[(1-iex)*y + iex*x]
c         z --> isz*z
          n1=isx(igroup)*((1-iex(igroup))*Rc(1,ic)+iex(igroup)*Rc(2,ic))
          m1=isy(igroup)*((1-iex(igroup))*Rc(2,ic)+iex(igroup)*Rc(1,ic))
          l1=isz(igroup)*Rc(3,ic)
          do n=1,2  ! map (n1,m1) back into the cluster.
            vx=real(n1,kind)-origx
            vy=real(m1,kind)-origy
            vz=real(l1,kind)-origz
            r1=(vx*a(1,1)+vy*a(1,2))/                     ! a(1,3)=0
     &         real(a(1,1)*a(1,1)+a(1,2)*a(1,2),kind)
            r2=(vx*a(2,1)+vy*a(2,2))/                     ! a(2,3)=0
     &         real(a(2,1)*a(2,1)+a(2,2)*a(2,2),kind)
            r3=vz*a(3,3)/real(a(3,3)**2,kind)             ! a(1,3)=a(2,3)=0
            if(r1.lt.0.0) then
              n1=n1+a(1,1)
              m1=m1+a(1,2)
            else if(r1.gt.1.0) then
              n1=n1-a(1,1)
              m1=m1-a(1,2)
            end if
            if(r2.lt.0.0) then
              n1=n1+a(2,1)
              m1=m1+a(2,2)
            else if(r2.gt.1.0) then
              n1=n1-a(2,1)
              m1=m1-a(2,2)
            end if
            if(r3.lt.0.0) then
              l1=l1+a(3,3)
            else if(r3.gt.1.0) then
              l1=l1-a(3,3)
            end if
          end do
c         Figure out which point (n1,m1,ll) it is
          itest=0
          do jc=1,Nc
           if(n1.eq.Rc(1,jc).and.m1.eq.Rc(2,jc).and.l1.eq.Rc(3,jc)) then
              icrequ(ic,igroup)=jc
              itest=1
            end if
          end do
          if(itest.ne.1) then
            write(lud,*) 'No equivalent found in icrequ',ic,igroup
            stop
          end if 
          if(iprint.ge.19.and.myrank.eq.0) 
     &        write(lud,*) ic,igroup,icrequ(ic,igroup)                !0
        end do
        end do                               

c***********************************************************************        
c       Now create the tables for the differences of R
c***********************************************************************        
        if(iprint.ge.19.and.myrank.eq.0) write(lud,*) '  '               !0
        if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) 'differences in r'                            !0
        if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) '  ic   jc  icrdiff(ic,jc)'                   !0
        do i=1,Nc
        do j=1,Nc
          Rx=real(Rc(1,i)-Rc(1,j),kind)
          Ry=real(Rc(2,i)-Rc(2,j),kind)
          Rz=real(Rc(3,i)-Rc(3,j),kind)
          itest=0
          do n=-1,1          ! Look in the surrounding 27
          do m=-1,1          ! clusters for Rx and Ry
          do l=-1,1          
            vx=Rx-origx + n*a(1,1) + m*a(2,1) + l*a(3,1)
            vy=Ry-origy + n*a(1,2) + m*a(2,2) + l*a(3,2)
            vz=Rz-origz + n*a(1,3) + m*a(2,3) + l*a(3,3)
            r1=(vx*a(1,1)+vy*a(1,2)+vz*a(1,3))/
     &          real(a(1,1)**2+a(1,2)**2+a(1,3)**2,kind)
            r2=(vx*a(2,1)+vy*a(2,2)+vz*a(2,3))/
     &          real(a(2,1)**2+a(2,2)**2+a(2,3)**2,kind)
            r3=(vx*a(3,1)+vy*a(3,2)+vz*a(3,3))/
     &          real(a(3,1)**2+a(3,2)**2+a(3,3)**2,kind)
            if(r1.lt.oner-epsilon.and.r1.gt.zeror-epsilon.and.
     &        r2.lt.oner-epsilon.and.r2.gt.zeror-epsilon.and.
     &        r3.lt.oner-epsilon.and.r3.gt.zeror-epsilon) then
c             (Rx,Ry) is located within the cluster.  
c             Now we must find which point it is.
              do ic=1,Nc
                r1 = Rx + n*a(1,1) + m*a(2,1) + l*a(3,1)
                r2 = Ry + n*a(1,2) + m*a(2,2) + l*a(3,2)
                r3 = Rz + n*a(1,3) + m*a(2,3) + l*a(3,3)
                if(abs(r1-Rc(1,ic)).lt.epsilon.and.
     &             abs(r2-Rc(2,ic)).lt.epsilon.and.
     &             abs(r3-Rc(3,ic)).lt.epsilon) then
                   icrdiff(i,j)=ic
                   itest=1
                end if
              end do
              if(iprint.ge.19.and.myrank.eq.0) 
     &             write(lud,*) i,j,icrdiff(i,j)
            end if
          end do
          end do
          end do
          if(itest.ne.1) then
            write(lud,*) 'no icrdiff(i,j) found',i,j
            stop
          end if
        end do
        end do

c***********************************************************************        
c       Now set up the K-space tables
c***********************************************************************        
c        
c
c
c
c***********************************************************************        
c       Create the cluster momentum tables Kc(i,ic), i=1,x 2,y 3,z
c***********************************************************************        

c       First find the points in the irreducible wedge.  They will be
c       indexed from 1 to Ncw
        if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) '  ' 
        if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) '  Wedge K points  '                                           !0
        if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) ' ic  i  j  l  Kcx(ic)   Kcy(ic)  Kxz(ic)'     
        ic=0
        Ncw=0
        icK0=1
        icKpi=1

        if(iodd.eq.0) then   ! The cluster is of high symmetry
          do i=0,Nc          ! the wedge is roughly 1/16 of the 1BZ.
          do j=0,i           ! Index j runs only from 0 to i.
          do l=0,Nc
            Kx=i*gvector(1,1) + j*gvector(2,1) + l*gvector(3,1)
            Ky=i*gvector(1,2) + j*gvector(2,2) + l*gvector(3,2)
            Kz=i*gvector(1,3) + j*gvector(2,3) + l*gvector(3,3)
            if(Kx.gt.-pi+epsilon.and.Kx.lt.pi+epsilon.and.
     &        Ky.gt.-pi+epsilon.and.Ky.lt.pi+epsilon.and.
     &        Kz.gt.-pi+epsilon.and.Kz.lt.pi+epsilon) then
c             Kx,Ky and Kz fall within the first Brillouin zone
              ic=ic+1
              Ncw=Ncw+1
              Kc(1,ic)=Kx
              Kc(2,ic)=Ky
              Kc(3,ic)=Kz
              ig(ic)=i
              jg(ic)=j
              kg(ic)=l
              if(abs(Kx).lt.epsilon.and.abs(Ky).lt.epsilon.
     &		and.abs(Kz).lt.epsilon) icK0=ic                             ! (0,0,0) point
              if(abs(Kx-pi).lt.epsilon.and.abs(Ky-pi).lt.epsilon
     &          .and.abs(Kz-pi).lt.epsilon) icKpi=ic     ! (pi,pi,pi) point
              if(iprint.ge.19.and.myrank.eq.0) 
     &            write(lud,5) ic,i,j,l,Kc(1,ic),Kc(2,ic),Kc(3,ic)        
            end if        
          end do
          end do
          end do
        else
          if(myrank.eq.0) write(lud,*) 'woops, iodd=',iodd
          stop
        end if
 5            format(1x,i2,1x,i2,1x,i2,1x,i2,2x,f9.6,2x,f9.6,2x,f9.6)
        
        ntw=nl*Ncw
        if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) 'Nc= ',Nc,' Ncw=',Ncw  

      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) ' Non Wedge K points '    !0
      end if                    
        
c       Now find the remaining points within the 1Bz and outside the 
c       wedge.  These points are indexed from Ncw+1 to Nc       
        do i=-Nc,Nc
        do j=-Nc,Nc
        do l=-Nc,Nc
          Kx=i*gvector(1,1) + j*gvector(2,1) + l*gvector(3,1)
          Ky=i*gvector(1,2) + j*gvector(2,2) + l*gvector(3,2)
          Kz=i*gvector(1,3) + j*gvector(2,3) + l*gvector(3,3)
          itest=0
          do ick=1,Ncw
            if(abs(Kx-Kc(1,ick)).lt.epsilon
     &         .and.
     &         abs(Ky-Kc(2,ick)).lt.epsilon
     &         .and.
     &         abs(Kz-Kc(3,ick)).lt.epsilon)then
c              This K is within the wedge angle
               itest=1
               exit ! exit do loop
            end if
          end do
          if(itest.eq.0.and.Kx.gt.-pi+epsilon.and.Kx.lt.pi+epsilon.and.
     &      Ky.gt.-pi+epsilon.and.Ky.lt.pi+epsilon.and.
     &      Kz.gt.-pi+epsilon.and.Kz.lt.pi+epsilon) then
c           Kx and Ky fall within the first zone and are not IRW
            ic=ic+1
            Kc(1,ic)=Kx
            Kc(2,ic)=Ky
            Kc(3,ic)=Kz
            ig(ic)=i
            jg(ic)=j
            kg(ic)=l
            if(iprint.ge.19.and.myrank.eq.0) 
     &          write(lud,5) ic,i,j,l,Kc(1,ic),Kc(2,ic),Kc(3,ic)
          end if
        end do
        end do
        end do
        
        if(ic.ne.Nc) then
          if(myrank.eq.0) write(lud,*) 'Bug in Kc(ic) loop'              !0
          stop
        end if
        if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) 'icK0=',icK0                                     !0
        if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) 'icKpi=',icKpi                                   !0

        
c***********************************************************************        
c       Create the indirect addressing table ict(i,j) which maps a 
c       point indexed by the principle translations in k-space to ic.
c***********************************************************************        
        do i=-N_sp,N_sp
        do j=-N_sp,N_sp
        do l=-2,2
          Kx=i*gvector(1,1) + j*gvector(2,1) + l*gvector(3,1)
          Ky=i*gvector(1,2) + j*gvector(2,2) + l*gvector(3,2)
          Kz=i*gvector(1,3) + j*gvector(2,3) + l*gvector(3,3)          
          if(Kx.lt.-pi+epsilon) Kx=Kx+2.0*pi
          if(Kx.gt.pi+epsilon) Kx=Kx-2.0*pi
          if(Ky.lt.-pi+epsilon) Ky=Ky+2.0*pi
          if(Ky.gt.pi+epsilon) Ky=Ky-2.0*pi
          if(Kz.lt.-pi+epsilon) Kz=Kz+2.0*pi
          if(Kz.gt.pi+epsilon) Kz=Kz-2.0*pi
          do jc=1,Nc
            if(abs(Kx-Kc(1,jc)).lt.epsilon
     &         .and.
     &         abs(Ky-Kc(2,jc)).lt.epsilon
     &         .and.
     &         abs(Kz-Kc(3,jc)).lt.epsilon) then    
              ict(i,j,l)=jc
            end if
          end do
        end do
        end do
        end do

c***********************************************************************        
c       Create the table ickdiff(ic1,ic2), which takes the difference
c       between two K-points
c***********************************************************************        
        if(iprint.ge.19.and.myrank.eq.0) write(lud,*) '   '               !0
        if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) ' ic1    ic2  ickdiff(ic1,ic2)  '              !0
        do ic1=1,Nc
        do ic2=1,Nc
          Kx=Kc(1,ic1)-Kc(1,ic2)
          Ky=Kc(2,ic1)-Kc(2,ic2)
          Kz=Kc(3,ic1)-Kc(3,ic2)
          if(Kx.lt.-pi+epsilon) Kx=Kx+2.0*pi
          if(Kx.gt.pi+epsilon) Kx=Kx-2.0*pi
          if(Ky.lt.-pi+epsilon) Ky=Ky+2.0*pi
          if(Ky.gt.pi+epsilon) Ky=Ky-2.0*pi
          if(Kz.lt.-pi+epsilon) Kz=Kz+2.0*pi
          if(Kz.gt.pi+epsilon) Kz=Kz-2.0*pi
c         Kdiffx(ic1,ic2)=Kx
c         Kdiffy(ic1,ic2)=Ky
c         Now we need to figure out where this point is!
          do ic=1,nc
            if(abs(Kx-Kc(1,ic)).lt.epsilon
     &         .and.
     &         abs(Ky-Kc(2,ic)).lt.epsilon
     &         .and.
     &         abs(Kz-Kc(3,ic)).lt.epsilon) then    
              ickdiff(ic1,ic2)=ic
            end if
          end do
          if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) ic1,ic2,ickdiff(ic1,ic2)                       !0
        end do
        end do

c***********************************************************************        
c       Create the table ickplus(ic1,ic2), which takes the sum
c       between to K-points
c***********************************************************************        
        if(iprint.ge.19.and.myrank.eq.0) write(lud,*) '   '               !0
        if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) ' ic1    ic2  ickplus(ic1,ic2)  '              !0
        do ic1=1,Nc
        do ic2=1,Nc
          Kx=Kc(1,ic1)+Kc(1,ic2)
          Ky=Kc(2,ic1)+Kc(2,ic2)
          Kz=Kc(3,ic1)+Kc(3,ic2)
          if(Kx.lt.-pi+epsilon) Kx=Kx+2.0*pi
          if(Kx.gt.pi+epsilon) Kx=Kx-2.0*pi
          if(Ky.lt.-pi+epsilon) Ky=Ky+2.0*pi
          if(Ky.gt.pi+epsilon) Ky=Ky-2.0*pi
          if(Kz.lt.-pi+epsilon) Kz=Kz+2.0*pi
          if(Kz.gt.pi+epsilon) Kz=Kz-2.0*pi
c         Kplusx(ic1,ic2)=Kx
c         Kplusy(ic1,ic2)=Ky
c         Now we need to figure out where this point is!
          do ic=1,nc
            if(abs(Kx-Kc(1,ic)).lt.epsilon.
     &		and.abs(Ky-Kc(2,ic)).lt.epsilon
     &          .and.abs(Kz-Kc(3,ic)).lt.epsilon) then
              ickplus(ic1,ic2)=ic
            end if
          end do
          if(iprint.ge.19.and.myrank.eq.0) 
     &       write(lud,*) ic1,ic2,ickplus(ic1,ic2)                       !0
        end do
        end do
        
c***********************************************************************                
c       Now find the points equivalent to ic, ickequ(ic,igroup) where 
c       igroup is the group operation.
c***********************************************************************        
c
c
        if(iprint.ge.19.and.myrank.eq.0) write(lud,*) '   '               !0
        if(iprint.ge.19.and.myrank.eq.0) write(lud,*) 'equivalency in k'  !0
        if(iprint.ge.19.and.myrank.eq.0) 
     &              write(lud,*) '   ic     igroup   ickequ(ic,igroup)'  !0
        do ic=1,Nc
          i=ig(ic)
          j=jg(ic)
          l=kg(ic)
c         Kx and Ky fall within the zone, by construction, now look 
c         for equivalent points
c
c  n --> isx*[(1-iex)*i + iex*j]
c  m --> isy*[(1-iex)*j + iex*i]
c
          do igroup=1,ngroup
            n=isx(igroup)*((1-iex(igroup))*i + iex(igroup)*j)
            m=isy(igroup)*((1-iex(igroup))*j + iex(igroup)*i)
            l=isz(igroup)*l
            Kx=n*gvector(1,1) + m*gvector(2,1) + l*gvector(3,1)
            Ky=n*gvector(1,2) + m*gvector(2,2) + l*gvector(3,2)
            Kz=n*gvector(1,3) + m*gvector(2,3) + l*gvector(3,3)
c           These are the equivalent points.  Now map to 1st BZ
            if(Kx.lt.-pi+epsilon) Kx=Kx+2.0*pi
            if(Kx.gt.pi+epsilon) Kx=Kx-2.0*pi
            if(Ky.lt.-pi+epsilon) Ky=Ky+2.0*pi
            if(Ky.gt.pi+epsilon) Ky=Ky-2.0*pi
            if(Kz.lt.-pi+epsilon) Kz=Kz+2.0*pi
            if(Kz.gt.pi+epsilon) Kz=Kz-2.0*pi
c           Now we need to figure out where this point is!
            itest=0
            do jc=1,Nc
              if(abs(Kx-Kc(1,jc)).lt.epsilon
     &           .and.
     &           abs(Ky-Kc(2,jc)).lt.epsilon
     &           .and.
     &           abs(Kz-Kc(3,jc)).lt.epsilon) then
                ickequ(ic,igroup)=jc
                itest=1
              end if
            end do
            if(iprint.ge.19.and.myrank.eq.0) 
     &           write(lud,*) ic,igroup,ickequ(ic,igroup)                
            if(itest.eq.0) then
!              if(myrank.eq.0) write(lud,*) 'no equivalent K-point found' 
c              if(myrank.eq.0) write(lud,*) 'j1,n2,m2',j1,n2,m2           
!              if(myrank.eq.0) write(lud,*) 'Kx,Ky,Kz',Kx,Ky,Kz              
!              stop
            end if
          end do
        end do

        
c***********************************************************************        
c       Now generate a table which maps any K to an equivalent point
c       in the irreducible wedge (IW), and a table for the degeneracy 
c       of each point in the IW. 
c 
c       The use of the wedge can greatly reduce the storage required.
c***********************************************************************        
c       

        ickdeg=0
c       scan over all of the K points in the 1BZ
        do ic=1,Nc
          i=ig(ic)
          j=jg(ic)
          l=kg(ic)
c         scan over all allowed point group operations, to find the
c         one that maps ic into the IRW
          itest=0
          groupops: do igroup=1,ngroup
            n=isx(igroup)*((1-iex(igroup))*i + iex(igroup)*j)
            m=isy(igroup)*((1-iex(igroup))*j + iex(igroup)*i)
            l=isz(igroup)*l
            Kx=n*gvector(1,1) + m*gvector(2,1) + l*gvector(3,1)
            Ky=n*gvector(1,2) + m*gvector(2,2) + l*gvector(3,2)
            Kz=n*gvector(1,3) + m*gvector(2,3) + l*gvector(3,3)
c           These are the equivalent points.  Now map to 1st BZ
            if(Kx.lt.-pi+epsilon) Kx=Kx+2.0*pi
            if(Kx.gt.pi+epsilon) Kx=Kx-2.0*pi
            if(Ky.lt.-pi+epsilon) Ky=Ky+2.0*pi
            if(Ky.gt.pi+epsilon) Ky=Ky-2.0*pi
            if(Kz.lt.-pi+epsilon) Kz=Kz+2.0*pi
            if(Kz.gt.pi+epsilon) Kz=Kz-2.0*pi
c           See if this point is in the IRW!
            do jc=1,Ncw
              if(abs(Kx-Kc(1,jc)).lt.epsilon
     &           .and.
     &           abs(Ky-Kc(2,jc)).lt.epsilon
     &           .and.
     &           abs(Kz-Kc(3,jc)).lt.epsilon) then
                ickmap(ic)=jc
                ickdeg(jc)=ickdeg(jc)+1
                irwgroup(ic)=igroup
                itest=1
              end if
            end do
            if(itest.eq.1) exit ! at first group op. that takes ic into IRW
          end do groupops
          if(itest.eq.0) then
            if(myrank.eq.0) write(lud,*) 'mapping failed for'
            if(myrank.eq.0) write(lud,*) 'ic',ic    
            if(myrank.eq.0) write(lud,*) 'Kx,Ky,Kz',Kx,Ky,Kz 
            stop
          end if
          
        end do

c       Test the degeneracy table.  Since ickdeg holds the degeneracy
c       of points in the IRW, its sum must equal Nc.
        i=0
        do ic=1,Ncw
          i=i+ickdeg(ic)
        end do
        if(i.ne.Nc) then
          if(myrank.eq.0) write(lud,*) 'bug in ickdeg',i
          stop
        end if
        
c       Test the group table.  irwgroup(ic within wedge)=1, the identity
        do ic=1,Ncw
          if(irwgroup(ic).ne.1) then
            if(myrank.eq.0) write(lud,*) 'irwgroup failed for wedge pt.'
            if(myrank.eq.0) write(lud,*) 'ic',ic
            if(myrank.eq.0) write(lud,*) 'Kx,Ky,Kz',Kx,Ky,Kz
          end if
        end do


!***********************************************************************
!     generate the kt-points in the Brillouin-zone of the
!     super-lattice, i.e. in the Wigner-Seitz cell of the
!     reciprocal space defined by the cluster K-points
!***********************************************************************

	ntot = nover
	Dkt=oner/real(ntot,kind)
	nwsc=0 ! number of k-points in Wigner-Seitz cell
	if (ntot > 0) then
	do i=-ntot+1,ntot
	do j=-ntot+1,ntot
	do k=1,ntot
	    iv1=(/i,j,k/)
	    rv1(:)=Dkt*(iv1(:)-halfr)*pi
	!  Check if this k-point is in Wigner-Seitz cell, i.e. if the closest
	!  K-point is K=0
	    r1=sum(rv1(:)*rv1(:)) ! distance from K=0
	    itest=0
	    do jc=1,Nc
		if (jc == icK0) cycle
		r2=rv1(1)-Kc(1,jc)
		r3=rv1(2)-Kc(2,jc)
		r4=rv1(3)-Kc(3,jc)
		if((r2-pi) > zeror) then; r2=r2-twor*pi; endif
		    if((r2+pi) < zeror) then; r2=r2+twor*pi; endif
			if((r3-pi) > zeror) then; r3=r3-twor*pi; endif
			    if((r3+pi) < zeror) then; r3=r3+twor*pi; endif
				if((r4-pi) > zeror) then; r4=r4-twor*pi; endif
				    if((r4+pi) < zeror) then; r4=r4+twor*pi; endif
			        r2=(r2*r2+r3*r3+r4*r4)
			        if(r2 < r1) then
			            itest=1
			            exit
			        endif
			    
			    enddo ! jc
			    if(itest == 0) then
			        nwsc=nwsc+1
			        kt(1:3,nwsc) = rv1(:)
			        nwsc=nwsc+1 ! -k must also be part of the WS cell
			        kt(1:3,nwsc) = -rv1(:)
			    endif
			enddo ! k
		    enddo ! j
		enddo ! i
	    else ! needed for FSS calculation
		nwsc=1
		kt(:,1)=zeror
	    endif

	    if(myrank == 0 .AND. iprint >= 19) then
		write(lud,*) 'Number of points in WS, nwsc', nwsc
		write(lud,*) 'Nc times nwsc', Nc*nwsc
		write(lud,*) 'Total number of points in BZ', (2*ntot)**3
	    endif



c***********************************************************************
c     Form Epsd(K) for calculating the partial density of states
c***********************************************************************   
      call allocate_PDOS
      estep = (eend-estart)/dble(e_div)
      PDOS(:,:) = 0.d0
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) '    '
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) 
     &  ' ic      Epsbar(ic)        temp_min      temp_max    min_Eps(ic)    max_Eps(ic)'
	Epsbar=zeror
      
        Dkt=0.5d0/real(nover,kind)
        do ic=1,Nc
	temp_max = zeror 
	temp_min = zeror
          do i=-nover+1,nover
          do j=-nover+1,nover
          do l=-nover+1,nover
            kx=Kc(1,ic) + Dkt*((real(i,kind)-0.5)*gvector(1,1) +
     &                       (real(j,kind)-0.5)*gvector(2,1)   +
     &                       (real(l,kind)-0.5)*gvector(3,1))
            ky=Kc(2,ic) + Dkt*((real(i,kind)-0.5)*gvector(1,2) +
     &                       (real(j,kind)-0.5)*gvector(2,2)   +
     &                       (real(l,kind)-0.5)*gvector(3,2))
            kz=Kc(3,ic) + Dkt*((real(i,kind)-0.5)*gvector(1,3) +
     &                       (real(j,kind)-0.5)*gvector(2,3)   +
     &                       (real(l,kind)-0.5)*gvector(3,3))
            Epsbar(ic)=Epsbar(ic) -0.5d0*(cos(kx)+cos(ky)+cos(kz)) -
     &                 tprime*(cos(kx)*cos(ky)+cos(ky)*cos(kz) +
     &			cos(kx)*cos(kz)-1.0) 
           Epsd =-0.5d0*(cos(kx)+cos(ky)+cos(kz)) -
     &                 tprime*(cos(kx)*cos(ky)+cos(ky)*cos(kz) +
     &			cos(kx)*cos(kz)-1.0) 

		e_index = NINT((Epsd-estart)/estep)
		if (e_index >= 0 .and. e_index <= e_div) then
		    PDOS(e_index,ic) = PDOS(e_index,ic) + 1.d0
		else
		    write(*,*) "WARNING: ENERGY EXCEEDS : ", Epsd, estart, eend
		endif
	  if (Epsd.gt.temp_max) temp_max = Epsd
	  if (Epsd.lt.temp_min) temp_min = Epsd   
          end do
          end do
          end do
          Epsbar(ic)=Epsbar(ic)/real((2*nover)**3,kind)
	dEps(ic) = (temp_max - temp_min)/real(nover,kind)
	min_Eps(ic) = temp_min
	max_Eps(ic) = temp_max
       if(myrank.eq.0.and.iprint.ge.19) 
     &	write(lud,"(i6,4x,f12.8,4x,f12.8,4x,f12.8,4x,f12.8,4x,f12.8)") 
     &     ic, Epsbar(ic), temp_min, temp_max,min_Eps(ic), max_Eps(ic)
      	end do


c***********************************************************************
c     Calculate the partial density of states
c
c                     Nc            
c       p_rho(K,w) =  --   SUM_kt delta( eps(K+kt) - w )
c                     N
c                       
c***********************************************************************  
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) '    '
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) 
     &  '  ic         ener(ic)         sum(ener) n_fact' !n_fact is the number of points within each cell ic
	
	n_fact = zeror
	do ic=1,Nc
	do i=0,e_div
	n_fact(ic) = n_fact(ic) + PDOS(i,ic)*estep
	end do
	end do


	do ic=1,Nc
	do i=0,e_div
	PDOS(i,ic) = PDOS(i,ic)/(n_fact(ic)*Nc)
	write(53,*) estart+i*estep, PDOS(i,ic) 
	end do
	write(53,*) '   '
	end do

	do i=0,e_div
        write(53,*) estart+i*estep, sum(PDOS(i,:))
        end do
	close(53)

 	  write(4005,403) 
403      format('#    count         pdos')
	ener = zeror
	do ic =1,Nc
	do i=0,e_div
	write(4005,*) estart+i*estep,PDOS(i,ic)
		      ener(ic) = ener(ic) + PDOS(i,ic)*estep
	end do
c	write(4006,*) ic, ener(ic), sum(ener)
       if(myrank.eq.0.and.iprint.ge.19) 
     &	write(lud,"(i6,4x,f12.8,4x,f12.8,4x,f18.8)") 
     &   ic, ener(ic),  sum(ener), n_fact(ic)
	end do

      p_rho_GcKf = zeror
      do ic=1,Nc
        do n=-nwn+1,nwn-1
	  do i=0,e_div-1
	    if ((wn(n) - (estart+float(i)*estep))*(wn(n) - (estart+float(i-1)*estep)).le.0) then
		p_rho_GcKf(ic,n) = halfr*(PDOS(i,ic) + PDOS(i+1,ic))
	    end if
	  end do
	end do
      end do

c***********************************************************************
c     Form Epsbar_histo(K) for calculating the Histogram
c***********************************************************************     
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) '    '
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) 
     &  ' ic,Epsbar_histo(ic),Epsbar(ic),wn(Epsbar_histo(ic))'
c	Epsbar=zeror
c	do ic=1,Nc
c         do i=1,nwsc
c            rv1(1) = Kc(1,ic)+kt(1,i)
c            rv1(2) = Kc(2,ic)+kt(2,i)
c            rv1(3) = Kc(3,ic)+kt(3,i)
c            Epsbar(ic)= Epsbar(ic)
c     &      -   halfr*(cos(rv1(1))+cos(rv1(2))+cos(rv1(3))) 
c     &            - tprime*(cos(rv1(1))*cos(rv1(2))
c     &                     +cos(rv1(2))*cos(rv1(3))
c     &                     +cos(rv1(1))*cos(rv1(3))-oner)
c                         end do ! i
c                      Epsbar(ic)=Epsbar(ic)/real(nwsc,kind)
c                      do n=-nwn,nwn ! Frequency loop starts
c	           
c	       if (wn(n) - Epsbar(ic).lt. epsilon) Epsbar_histo(ic) = n
c          end do !n
c        if(myrank.eq.0.and.iprint.ge.19) 
c     &	write(lud,"(i6,4x,i8,4x,f14.8,4x,f14.8)") 
c     &     ic,Epsbar_histo(ic),Epsbar(ic),wn(Epsbar_histo(ic))
c	end do !ic

!   finding when wn=eps_bar(ic)
!  1.  define function to find min of:
          funct(:,:)=zeror
      	  do ic=1,Nc
	    do n=-nwn,nwn
              funct(ic,n)=wn(n)-Epsbar(ic)
    	    end do
    	  end do
!  2.  find min of this function:  
          
          do ic=1,Nc   
    	     min_funct(ic)=abs(funct(ic,-nwn))
    	     Epsbar_histo(ic)=-nwn 
    	      do n=-nwn,nwn
    	       if(min_funct(ic) .gt. abs(funct(ic,n))) then 
    	        min_funct(ic)=abs(funct(ic,n))
    	        Epsbar_histo(ic)=n      ! index at which wn=eps_bar
    	       end if	
              end do ! n
        if(myrank.eq.0.and.iprint.ge.19) 
     &	write(lud,"(i6,4x,i8,4x,f14.8,4x,f14.8)") 
     &     ic,Epsbar_histo(ic),Epsbar(ic),wn(Epsbar_histo(ic))
           end do 
  	  
        
c***********************************************************************        
c       Form the Fourier transform tables.  
c       R-->K,  FTCoefs_K_to_R
c       K-->R,  FTCoefs_K_to_R
c
c                    __          the factor of 1/N is required
c                 1  \           so that both G(R)~1/w and G(K)~1/w
c       G(R=0) = --- /  G(K)     for large w
c                 N  --
c                    K
c
c***********************************************************************        
        do ic=1,Nc
        do ik=1,Nc
          r1=0.0
          do idim=1,ndim
            r1=r1+Kc(idim,ik)*Rc(idim,ic)
          end do
          FTCoefs_K_to_R(ik,ic)=exp(-ii*r1)/real(Nc,kind)
          FTCoefs_R_to_K(ik,ic)=exp(+ii*r1)
        end do
        end do

        return
        end subroutine tables_3D

c***************************************************************************
        subroutine tables_3D_cubic                  !3 DIMENSIONAL VERSION
c***************************************************************************
c       This subroutine makes the lookup tables used throughout the
c       3 DIMENSIONAL code. This for a perfect 3D cubic system
c***************************************************************************************
        use Global
        implicit none
c***************************************************************************************
      logical, parameter :: use_original = .false.
      real (kind) :: yy,del,deltat,tc
c***************************************************************************************
      interface transform3D
      function gtransform3D(a,igroup,ixy,ixz,iyz,isx,isy,isz)
c      use global
      integer :: a(3,3),gtransform3D(3,3),igroup,ixy(48),ixz(48),
     &           iyz(48),isx(48),isy(48),isz(48)
      end
      
      function Rtransform3D(Rc1,igroup,ixy,ixz,iyz,isx,isy,isz)
      integer :: Rc1(3),Rtransform3D(3),igroup,ixy(48),ixz(48),iyz(48),
     &     isx(48),isy(48),isz(48)
      end
      
      function Ktransform3D(Kc1,igroup,ixy,ixz,iyz,isx,isy,isz)
c      use global
      integer, parameter :: kind=8
      integer    :: ixy(48),ixz(48),iyz(48),
     &     isx(48),isy(48),isz(48)
      real(kind) :: Kc1(3),Ktransform3D(3)
      end
      end interface 
      
      interface Inverse
      function Inverse_3x3(a)
c      use global
      integer, parameter :: kind=8
      integer :: a(3,3)
      real(kind) :: Inverse_3x3(3,3)
      end function
      end interface
      
      interface determinant
      function Ideterminant(b)
c      use global
      integer, parameter :: kind=8
      integer :: b(3,3),Ideterminant
      end function
      function Rdeterminant(c)
      use global
      real(kind) :: c(3,3),Rdeterminant
      end function
      end interface

c***************************************************************************
c     Locally defined
c***************************************************************************
      logical    :: allowed,denied,denied_lin_alg
      integer, parameter :: maxnc=300
      integer, parameter :: nshellmax=100
      integer :: i,i1,i2,i3,it,ic,ic1,ic2,idim,I_F,Nm,I_FShell,
     &           ik,itest,j,jc,k,kroniker,
     &           l,m,n,p,Lc,iv1(3),iv2(3),
     &           ig(maxnc),jg(maxnc),kg(maxnc),a1(3,3),
     &           ixy(48),ixz(48),iyz(48),isx(48),isy(48),isz(48),
     &           irwgroup(maxnc),group(48),
     &		 p_r,nneigh(nshellmax),nneighl(nshellmax),e_index,
     &     nneighp(nshellmax),ishellmax,nshell,nneigh1(108,6)
      real(kind) :: r1,r2,r3,r4,Dkt,c1(3,3),ener(Nc),
     &              rv1(3),rv2(3),rv3(3),n_fact(Nc),
     &              corner(3,8),Bettsd(4),Bettsf(6),Bettsl(3),
     &              meand,meanf,meanl,cubicity,temp_min,temp_max
c For lin algebra test
      real(kind) :: ainv(3,3),a1ainv(3,3)
      real(kind) :: funct(Nc,-nwn:nwn),min_funct(Nc)
c***************************************************************************************
c Select use of "brute force" vs linear algebra method for 
c checking lattice coincidence. Both methods should agree; lin alg
c method is much faster.
c***************************************************************************************
      logical, parameter :: use_brute_force = .false.
      logical, parameter :: use_lin_alg = .true.

        if(myrank.eq.0.and.iprint.gt.19) then
          write(lud,*) 'start Table'
          flush(lud)
        end if
c***************************************************************************

c  ----------------------------------------------------------------------------
c                       GROUP OPERATIONS
c  ----------------------------------------------------------------------------
c  Only 16A has all 48 operations and has the point group Oh.  The remaining clusters 
c  have lower symmetry than the cubic lattice, with symmetries designated by the
c  Schoenflies symbols (From Ibach and Luth):
c  Cj (j=2,3,4, 6) j-fold rotation axis 
c  Sj j-fold rotation-inversion axis 
c  Dj j 2-fold rotation axes perpendicular to a j-fold principle rotation axis 
c  T  4 three-and 3 two-fold rotation axes, as in a tetrahedron 
c  O  4 three-and 3 four-fold rotation axes, as in a octahedron 
c  Ci a center of inversion 
c  Cs a mirror plane
c  In addition, their are sufixes for mirror planes 
c     h: horizontal=perpendicular to the c  rotation axis, 
c     v: vertical=parallel to the main rotation axis in the plane, 
c     d: diagonal=parallel to the main rotation axis in the plane  
c        bisecting the two-fold rotation axes.
c 
c  The Clusters identified by eetts, have several different symmetries
c
c  C2h   A two fold axis, plus a single mirror plane that is perpendicular to the axis
c  Ci    Contains only the identity and inversion
c  C3i   Triagonal with a 3-fold axis and inversion (GUESS)??
c          * rotations by 2pi/3 about the axix x=y=z
c          * inversion wrt the origin.
c  Oh    The full symmetry of the cubic lattice 
c  D2h   2 2-fold axes perp. to the 2-fold rotation axis, plus a horizonal mirror plane
c  D3h   3 2-fold axes perp. to the 3-fold rotation axis, plus a horizonal mirror plane
c
c  We need to automate the process of finding which of the 48 point group operations are
c  retained by the different clusters.  To do this, we must perform all 48 cubic group
c  operations on the {a} and see if they form an equivalent parallelpiped.  
c
c  Group Operations for a cubic lattice:
c
c  x --> +/- x                  There are 2^D D! operations in D dimensions
c    \/                         or 48 operations in 3D.
c    /\
c  y --> +/- y   
c    \/
c    /\
c  z --> +/- z   
c        .
c  a'=transform3D(a,igroup,ixy,ixz,iyz,isx,isy,isz)
c  ixy=0(1) dont (do) exchange x and y
c  ixz=0(1) dont (do) exchange x and z
c  iyz=0(1) dont (do) exchange y and z
c  isx      sign of x 
c  isy      sign of y
c  isz      sign of z
c
c  x <-- isx*[(1-ixy)*x + ixy*y]  <-- [(1-ixz)*x + ixz*z]<----------------------
c  y <-- isy*[(1-ixy)*y + ixy*x]  <----------------------<-- [(1-iyz)*y + iyz*z]
c  z <-- isz*<------------------- <-- [(1-ixz)*z + ixz*x]<-- [(1-iyz)*z + iyz*y]
c

      ixy( 1)=0;ixz( 1)=0; iyz( 1)=0; isx( 1)=+1; isy( 1)=+1; isz( 1)=+1 !identity  xyz e
      ixy( 2)=1;ixz( 2)=0; iyz( 2)=0; isx( 2)=+1; isy( 2)=+1; isz( 2)=+1 !ref. x=y  yxz o
      ixy( 3)=0;ixz( 3)=1; iyz( 3)=0; isx( 3)=+1; isy( 3)=+1; isz( 3)=+1 !ref. x=z  zyx o
      ixy( 4)=0;ixz( 4)=0; iyz( 4)=1; isx( 4)=+1; isy( 4)=+1; isz( 4)=+1 !ref. y=z  xzy o
      ixy( 5)=1;ixz( 5)=1; iyz( 5)=0; isx( 5)=+1; isy( 5)=+1; isz( 5)=+1 !          yzx e
      ixy( 6)=1;ixz( 6)=0; iyz( 6)=1; isx( 6)=+1; isy( 6)=+1; isz( 6)=+1 !          zxy e
      
      ixy( 7)=0;ixz( 7)=0; iyz( 7)=0; isx( 7)=-1; isy( 7)=+1; isz( 7)=+1 !
      ixy( 8)=1;ixz( 8)=0; iyz( 8)=0; isx( 8)=-1; isy( 8)=+1; isz( 8)=+1 !
      ixy( 9)=0;ixz( 9)=1; iyz( 9)=0; isx( 9)=-1; isy( 9)=+1; isz( 9)=+1 !
      ixy(10)=0;ixz(10)=0; iyz(10)=1; isx(10)=-1; isy(10)=+1; isz(10)=+1 !
      ixy(11)=1;ixz(11)=1; iyz(11)=0; isx(11)=-1; isy(11)=+1; isz(11)=+1 !
      ixy(12)=1;ixz(12)=0; iyz(12)=1; isx(12)=-1; isy(12)=+1; isz(12)=+1 !
     
      ixy(13)=0;ixz(13)=0; iyz(13)=0; isx(13)=+1; isy(13)=-1; isz(13)=+1 !
      ixy(14)=1;ixz(14)=0; iyz(14)=0; isx(14)=+1; isy(14)=-1; isz(14)=+1 !
      ixy(15)=0;ixz(15)=1; iyz(15)=0; isx(15)=+1; isy(15)=-1; isz(15)=+1 !
      ixy(16)=0;ixz(16)=0; iyz(16)=1; isx(16)=+1; isy(16)=-1; isz(16)=+1 !
      ixy(17)=1;ixz(17)=1; iyz(17)=0; isx(17)=+1; isy(17)=-1; isz(17)=+1 !
      ixy(18)=1;ixz(18)=0; iyz(18)=1; isx(18)=+1; isy(18)=-1; isz(18)=+1 !
      
      ixy(19)=0;ixz(19)=0; iyz(19)=0; isx(19)=-1; isy(19)=-1; isz(19)=+1 !
      ixy(20)=1;ixz(20)=0; iyz(20)=0; isx(20)=-1; isy(20)=-1; isz(20)=+1 !
      ixy(21)=0;ixz(21)=1; iyz(21)=0; isx(21)=-1; isy(21)=-1; isz(21)=+1 !
      ixy(22)=0;ixz(22)=0; iyz(22)=1; isx(22)=-1; isy(22)=-1; isz(22)=+1 !
      ixy(23)=1;ixz(23)=1; iyz(23)=0; isx(23)=-1; isy(23)=-1; isz(23)=+1 !
      ixy(24)=1;ixz(24)=0; iyz(24)=1; isx(24)=-1; isy(24)=-1; isz(24)=+1 !
      
      ixy(25)=0;ixz(25)=0; iyz(25)=0; isx(25)=+1; isy(25)=+1; isz(25)=-1 !
      ixy(26)=1;ixz(26)=0; iyz(26)=0; isx(26)=+1; isy(26)=+1; isz(26)=-1 !
      ixy(27)=0;ixz(27)=1; iyz(27)=0; isx(27)=+1; isy(27)=+1; isz(27)=-1 !
      ixy(28)=0;ixz(28)=0; iyz(28)=1; isx(28)=+1; isy(28)=+1; isz(28)=-1 !
      ixy(29)=1;ixz(29)=1; iyz(29)=0; isx(29)=+1; isy(29)=+1; isz(29)=-1 !
      ixy(30)=1;ixz(30)=0; iyz(30)=1; isx(30)=+1; isy(30)=+1; isz(30)=-1 !
      
      ixy(31)=0;ixz(31)=0; iyz(31)=0; isx(31)=-1; isy(31)=+1; isz(31)=-1 !
      ixy(32)=1;ixz(32)=0; iyz(32)=0; isx(32)=-1; isy(32)=+1; isz(32)=-1 !
      ixy(33)=0;ixz(33)=1; iyz(33)=0; isx(33)=-1; isy(33)=+1; isz(33)=-1 !
      ixy(34)=0;ixz(34)=0; iyz(34)=1; isx(34)=-1; isy(34)=+1; isz(34)=-1 !
      ixy(35)=1;ixz(35)=1; iyz(35)=0; isx(35)=-1; isy(35)=+1; isz(35)=-1 !
      ixy(36)=1;ixz(36)=0; iyz(36)=1; isx(36)=-1; isy(36)=+1; isz(36)=-1 !
      
      ixy(37)=0;ixz(37)=0; iyz(37)=0; isx(37)=+1; isy(37)=-1; isz(37)=-1 !
      ixy(38)=1;ixz(38)=0; iyz(38)=0; isx(38)=+1; isy(38)=-1; isz(38)=-1 !
      ixy(39)=0;ixz(39)=1; iyz(39)=0; isx(39)=+1; isy(39)=-1; isz(39)=-1 !
      ixy(40)=0;ixz(40)=0; iyz(40)=1; isx(40)=+1; isy(40)=-1; isz(40)=-1 !
      ixy(41)=1;ixz(41)=1; iyz(41)=0; isx(41)=+1; isy(41)=-1; isz(41)=-1 !
      ixy(42)=1;ixz(42)=0; iyz(42)=1; isx(42)=+1; isy(42)=-1; isz(42)=-1 !
      
      ixy(43)=0;ixz(43)=0; iyz(43)=0; isx(43)=-1; isy(43)=-1; isz(43)=-1 !
      ixy(44)=1;ixz(44)=0; iyz(44)=0; isx(44)=-1; isy(44)=-1; isz(44)=-1 !
      ixy(45)=0;ixz(45)=1; iyz(45)=0; isx(45)=-1; isy(45)=-1; isz(45)=-1 !
      ixy(46)=0;ixz(46)=0; iyz(46)=1; isx(46)=-1; isy(46)=-1; isz(46)=-1 !
      ixy(47)=1;ixz(47)=1; iyz(47)=0; isx(47)=-1; isy(47)=-1; isz(47)=-1 !
      ixy(48)=1;ixz(48)=0; iyz(48)=1; isx(48)=-1; isy(48)=-1; isz(48)=-1 !


c***********************************************************************        
c       determine the principle translation vectors a1,a2 and a3, and 
c       reciprocal space vectors g1, g2 and g3
c***********************************************************************
        a(3,1)=0; a(3,2)=0 
        a(1,3)=0; a(2,3)=0
	a(1,2)=0
        If(Nc.eq.1) then                              
          a(1,1)=1
          a(2,2)=1
          a(3,3)=1
          Ncw=1
        else if(Nc.eq.8) then                           
          a(1,1)=2                                 
          a(2,2)=2
	  a(3,3)=2
        else if(Nc.eq.27) then ! Avoid odd clusters                          
          a(1,1)=3  
	  a(2,2)=3                               
	  a(3,3)=3
        else if(Nc.eq.84) then                           
	  a(1,1)=4
	  a(2,2)=2
	  a(1,2)=3
	  a(3,3)=4 
	  a(3,2)=-2
	  a(2,3)=2 
!Another version of Nc=84
!	  a(1,1)=3
!	  a(2,2)=4
!	  a(1,2)=4
!	  a(3,3)=3 
        else if(Nc.eq.64) then                           
          a(1,1)=4  
	  a(2,2)=4                              
	  a(3,3)=4
        else if(Nc.eq.90) then                           
	  a(1,1)=2
	  a(2,2)=4
	  a(1,2)=4
	  a(3,3)=3 
	  a(2,3)=-3 
	  a(3,2)=3 
        else if(Nc.eq.112) then                           
	  a(1,1)=3
	  a(2,2)=4
	  a(1,2)=4
	  a(3,3)=4 
        else if(Nc.eq.125) then                           
          a(1,1)=5  
	  a(2,2)=5                              
	  a(3,3)=5
        else if(Nc.eq.216) then                           
          a(1,1)=6  
	  a(2,2)=6                              
	  a(3,3)=6
       elseif(abs(determinant(a)).ne.Nc) then
c     Test the a vectors      
        if(myrank.eq.0) write(lud,*) 'This cell volume isnt Nc',
     &                                abs(determinant(a))
        stop
      end if
   
        a(2,1)= -a(1,2)  

c***********************************************************************        
c       Now set up the K-space tables
c***********************************************************************        

c       Now calculate the principle translations in K-space
c   
c  
c                              _    _
c    _     _    _        2 pi aj x ak          
c    a --> g    g_i =  ------------------      plus cyclic permutations
c                      | ai . (aj x ak)|       
c                                           
c       Principle translation vectors for lattice tiling       
c                                        
c
      r1=determinant(a)
      
      gvector(1,1)=2*pi*( a(2,2)*a(3,3)-a(3,2)*a(2,3) )/r1  !a2y*a3z-a3y*a2z
      gvector(1,2)=2*pi*( a(2,3)*a(3,1)-a(3,3)*a(2,1) )/r1  !a2z*a3x-a3z*a2x
      gvector(1,3)=2*pi*( a(2,1)*a(3,2)-a(3,1)*a(2,2) )/r1  !a2x*a3y-a3x*a2y
      
      gvector(2,1)=2*pi*( a(3,2)*a(1,3)-a(1,2)*a(3,3) )/r1  !a3y*a1z-a1y*a3z
      gvector(2,2)=2*pi*( a(3,3)*a(1,1)-a(1,3)*a(3,1) )/r1  !a3z*a1x-a1z*a3x
      gvector(2,3)=2*pi*( a(3,1)*a(1,2)-a(1,1)*a(3,2) )/r1  !a3x*a1y-a1x*a3y
        
      gvector(3,1)=2*pi*( a(1,2)*a(2,3)-a(2,2)*a(1,3) )/r1  !a1y*a2z-a2y*a1z
      gvector(3,2)=2*pi*( a(1,3)*a(2,1)-a(2,3)*a(1,1) )/r1  !a1z*a2x-a2z*a1x
      gvector(3,3)=2*pi*( a(1,1)*a(2,2)-a(2,1)*a(1,2) )/r1  !a1x*a2y-a2x*a1y
      
      r1=abs(determinant(gvector))
cc      r1=determinant(gvector)

      if(abs((2*pi)**3/r1-Nc).gt.epsilon) then
        if(myrank.eq.0) write(lud,*) 'This cell volume is not Nc',r1
        stop
      end if
      do i=1,3
      do j=1,3
        r1=sum(a(i,:)*gvector(j,:))/(2*pi)
        if(i.eq.j.and.abs(r1-oner).gt.epsilon) then
          if(myrank.eq.0) write(lud,*) 'bug in g',i,j,r1
          stop
        end if
        if(i.ne.j.and.abs(r1).gt.epsilon) then
          if(myrank.eq.0) write(lud,*) 'bug in g',i,j,r1
          stop
        end if
      end do
      end do
      
      if(iprint.ge.19.and.myrank.eq.0) 
     &  write(lud,"(' gvector=',3(' (',f5.2,1x,f5.2,1x,f5.2,')'))") 
     &  (gvector(l,:),l=1,3)
!	call initsymm()
!        call findsymmops(a,group,ngroup)

c       Now calculate the principle translations in K-space
!        gvector(1,1)=twor*pi*( a(2,2)*a(3,3)-a(3,2)*a(2,3) )/real(a(1,1)*a(2,2)*a(3,3))
!        gvector(2,2)=twor*pi*( a(3,3)*a(1,1)-a(1,3)*a(3,1) )/real(a(1,1)*a(2,2)*a(3,3))
!        gvector(3,3)=twor*pi*( a(1,1)*a(2,2)-a(2,1)*a(1,2) )/real(a(1,1)*a(2,2)*a(3,3))
!
!	gvector(1,2)=zeror
!        gvector(1,3)=zeror
!        gvector(2,1)=gvector(1,2)
!        gvector(2,3)=zeror
!        gvector(3,1)=zeror
!        gvector(3,2)=zeror

c*************************************************************************
c       Set the location of the origin in real space, denote the symmetry
c       (high, iodd=0; low, iodd=1) and the number of group operations.
c*************************************************************************
        if(a(1,1).eq.a(1,2)) then  ! The xy part is a diamond (Nc=2,8,18,32,50...)
          write(*,* )'diamond',Nc
          origx= zeror
          origy=-halfr
          ngroup=16   ! with all ngroup=8 symmetry operations
        else if(a(1,1)*a(1,2).eq.0) then  ! The xy part is a square 
           write(*,* )'cubic',Nc
          origx=zeror
          origy=zeror
          origz=zeror
          ngroup=48  ! with all ngroup=48 symmetry operations
        end if

c****************************************************************************       
c       calculate the table of cluster point locations Rc(i,ic), i=1,x or 2,y
c****************************************************************************       
c 
c  ^ y
c  |
cc-|-c---c-*-c---c---c---c               A point P is within the cluster 
c| | |   |   |   |   |   |               if a vector v from the box origin, 
cc-|-c---c---c---c---c---c               O, to the point satisfies
ca2^ * * * * * * * * |   |               
cc-|-c---c---c---c-*-c---c                    _   _
c| | |   P   |   | * |   |                    v . a1
cc-|-c--/c---c---c-*-c---c               0 < -------   < 1
c| | | /v|   |   | * |   |                   |a1.a1|
cc-|-c/--c---c---c-*-c---c 
c| | /   |   |   | * |   |                      and
cc-|/0---c---c---c-*-c---c                    _   _
c| O--------------->---------> x              v . a2
cc---c---c---c---ca1-c---c               0 < -------   < 1
c|   |   |   |   |   |   |                   |a2.a2|
cc---c---c---c---c---c---c
c                                               and
c                                             _   _
c                                             v . a3
c                                        0 < -------   < 1
c                                            |a3.a3|
c

c     Check fixed dimensions (PK: Could be allocated)
      if (nc.gt.maxnc) then
        if(myrank.eq.0) write(lud,*) 'Nc>maxnc parameter',nc,maxnc,
     &        'Need to recompile'
        stop
      end if
c     Calculate the cubicity of the clusters      
c
c       c8-------c7    body diagonals   face diagonals
c      /|       /|     --------------   --------------
c     / |      / |       d1 (c1,c7)       f1 (c1,c3)
c    /  |     /  |       d2 (c2,c8)       f2 (c2,c4)
c   c5-------c6  |       d3 (c3,c5)       f3 (c1,c6)
c   |   c4---|---c3      d4 (c4,c6)       f4 (c2,c5)
c a3|  /     |  /                         f5 (c1,c9)
c   | / a2   | /                          f6 (c4,c5)
c   |/       |/
c   c1-------c2
c        a1
c
      corner(:,1) = zeror
      corner(:,2) = a(1,:)
      corner(:,3) = a(1,:)+a(2,:)
      corner(:,4) = a(2,:)
      corner(:,5) = a(3,:)
      corner(:,6) = a(3,:)+a(1,:)
      corner(:,7) = a(3,:)+a(1,:)+a(2,:)
      corner(:,8) = a(3,:)+a(2,:)
      
      Bettsd(1)=sqrt(sum((corner(:,1)-corner(:,7))**2))
      Bettsd(2)=sqrt(sum((corner(:,2)-corner(:,8))**2))
      Bettsd(3)=sqrt(sum((corner(:,3)-corner(:,5))**2))
      Bettsd(4)=sqrt(sum((corner(:,4)-corner(:,6))**2))
      
      Bettsf(1)=sqrt(sum((corner(:,1)-corner(:,3))**2))
      Bettsf(2)=sqrt(sum((corner(:,2)-corner(:,4))**2))
      Bettsf(3)=sqrt(sum((corner(:,1)-corner(:,6))**2))
      Bettsf(4)=sqrt(sum((corner(:,2)-corner(:,5))**2))
      Bettsf(5)=sqrt(sum((corner(:,1)-corner(:,8))**2))
      Bettsf(6)=sqrt(sum((corner(:,4)-corner(:,5))**2))
      
      Bettsl(1)=sqrt(dot_product(corner(:,2),corner(:,2)))
      Bettsl(2)=sqrt(dot_product(corner(:,4),corner(:,4)))
      Bettsl(3)=sqrt(dot_product(corner(:,5),corner(:,5)))
      
      meand=product(Bettsd)**(0.25000000000_kind)
      meanf=product(Bettsf)**(0.16666666667_kind)
      meanl=product(Bettsl)**(0.33333333333_kind)
      r1=1.732050808_kind*meanl/meand
      r2=1.414213562_kind*meanl/meanf
      cubicity=max(r1,oner/r1)*max(r2,oner/r2)
      if(myrank.eq.0) 
     &  write(lud,*) 'cubicity=',cubicity
      

c*********************************************************************************
c     Determine which group operations are allowed    
c*********************************************************************************
      
      ngroup=0; Lc=2*(Nc**.3333+1)
      if(myrank.eq.0.and.iprint.ge.19) then
        write(lud,*) '--------------------------------------------',
     &               '------------------' 
        write(lud,*) '             trans(a)        igroup L  ixy ix',
     &               'z iyz isx isy isz'
        write(lud,*) '--------------------------------------------',
     &               '------------------'
      end if
      do igroup=1,48
         a1=transform3D(a,igroup,ixy,ixz,iyz,isx,isy,isz)

         if (use_brute_force) then

c        look at the lattice points generated by a1 and see if they
c        correspond to those produced by a.
c        First, generate the a1 lattice points
         do i=-Lc,Lc
         do j=-Lc,Lc
         do k=-Lc,Lc
           denied=.true.
           iv1=(/i,j,k/)
           do l=1,3         ! x,ymz components of this lattice point
             rv1(l)=sum(iv1(:)*a1(:,l))
           end do
c          Now generate the a lattice points.  If a and a1 are the 
c          same lattice one of these points MUST correspond to ij
           do n=-4*Lc,4*Lc
           do m=-4*Lc,4*Lc
           do p=-4*Lc,4*Lc
             iv2=(/n,m,p/)
             do l=1,3         ! x,y components of this lattice point
               rv2(l)=sum(iv2(:)*a(:,l))
             end do
             if(sum(abs(rv1-rv2)).lt.epsilon) denied=.false.
           end do
           end do
           end do
           if(denied) exit !no corrs. point was found
         end do
         if(denied) exit
         end do
         if(denied) exit
         end do
         end if
         if (use_lin_alg) then
c PK If lattice points generated by "a1" vectors can be found in 
c lattice points generated by "a" vectors then (a1)x(a^-1) can only 
c contain integers.
            ainv=Inverse_3x3(a)
            a1ainv=matmul(a1,ainv)
            if (all(abs(a1ainv-nint(a1ainv)).le.epsilon)) then 
               denied_lin_alg=.false. ! matrix only contains integers
            else
               denied_lin_alg=.true.
            end if
         end if

c If both algorithms are enabled, check they agree
         if (use_brute_force.and.use_lin_alg) then
            if (denied.neqv.denied_lin_alg) then
               write (*,*) '*** BUG : Denied algorithm disagree'
               write (*,*) 'a11',a1(:,1)
               write (*,*) 'a12',a1(:,2)
               write (*,*) 'a13',a1(:,3)
               
               write (*,*) 'a1',a(:,1)
               write (*,*) 'a2',a(:,2)
               write (*,*) 'a3',a(:,3)
               write (*,*)
               write (*,*) 'ainv',ainv(:,1)
               write (*,*) 'ainv',ainv(:,2)
               write (*,*) 'ainv',ainv(:,3)
               write (*,*)
               
               write (*,*) 'a1ainv',a1ainv(:,1)
               write (*,*) 'a1ainv',a1ainv(:,2)
               write (*,*) 'a1ainv',a1ainv(:,3)
               write (*,*) 'Denied_lin_alg=',denied_lin_alg
               write (*,*) 'Denied_brute_force=',denied
               write (*,*)
               stop
            end if
         end if
         if (use_lin_alg) denied=denied_lin_alg

         if(denied.eqv..false.) then ! a corresponding point was found
            ngroup=ngroup+1
            group(ngroup)=igroup
         end if    

         if(myrank.eq.0.and.iprint.ge.19) then
           allowed=.true.;if (denied) allowed=.false.
           write(lud,"(3('(',i2,1x,i2,1x,i2,')'),2x,i2,2x,l,6(2x,i2))")
     &     a1(1,:),a1(2,:),a1(3,:),igroup,allowed,ixy(igroup),
     &     ixz(igroup),iyz(igroup),isx(igroup),isy(igroup),isz(igroup)
         end if
      end do

     
      if(myrank.eq.0.and.iprint.ge.19) 
     &  write(lud,"('cluster=',a3,' ngroup=',i2,' groups=',48(i2,1x))") 
     &  cluster,ngroup,group(1:ngroup)

c*********************************************************************************
c       calculate the table of cluster point locations Rc(i,ic), i=1,2,3 for x,y,z
c*********************************************************************************       
c 
c  ^ y
c  |
cc-|-c---c-*-c---c---c---c               A point P is within the cluster 
c| | |   |   |   |   |   |               if a vector v from the box origin, 
cc-|-c---c---c---c---c---c               O, to the point satisfies
ca2^ * * * * * * * * |   |               
cc-|-c---c---c---c-*-c---c                        
c| | |   P   |   | * |   |                            _   _
cc-|-c--/c---c---c-*-c---c                   -epsilon  <  v__ai  < 1-epsilon
c| | | /v|   |   | * |   |                        
cc-|-c/--c---c---c-*-c---c 
c| | /   |   |   | * |   |               for all i=1,3 where v__ai is the
cc-|/0---c---c---c-*-c---c               projection of the vector v along
c| O--------------->---------> x         ai
cc---c---c---c---ca1-c---c                              -1
c|   |   |   |   |   |   |               v = A n;  n = A  v ; where A=transpose(a)
cc---c---c---c---c---c---c
c                                        v = i*a1 + j*a2 + k*a3
c
c                                            [ a1 ]
c                                        a = [ a2 ]
c                                            [ a3 ]
c                                        
cc      do i=1,3
cc        a1(:,i)=a(i,:)
cc      end do
	a1=transpose(a)
      c1=Inverse(a1)
      if(iprint.ge.19.and.myrank.eq.0) write(lud,*) '  ' !0
      if(iprint.ge.19.and.myrank.eq.0) 
     &     write(lud,*) 'ic,Rcx(ic),Rcy(ic),Rcz(ic)' !0
      ic=0
      do i=-Nc,Nc
      do j=-Nc,Nc
      do k=-Nc,Nc
         rv1=(/i,j,k/)          ! The coordinates of the point.
         rv2=matmul(c1,rv1)
         if(rv2(1).lt.oner-epsilon.and.rv2(1).gt.zeror-epsilon.and.
     &      rv2(2).lt.oner-epsilon.and.rv2(2).gt.zeror-epsilon.and.
     &      rv2(3).lt.oner-epsilon.and.rv2(3).gt.zeror-epsilon ) then
            ic=ic+1
            Rc(:,ic)=(/i,j,k/) ! Rc, like i.e. a(1,1) is integer
            if(iprint.ge.19.and.myrank.eq.0)                             !0
     &           write(lud,"(1x,i3,3x,3(i2,6x))") ic,Rc(1,ic),Rc(2,ic), !0
     &                                            Rc(3,ic)              !0
            if(sum(abs(Rc(:,ic))).eq.0) icR0=ic
         end if
      end do
      end do
      end do
      if(ic.ne.Nc) then
         if(myrank.eq.0) write(lud,*) 'Bug in Rc(ic) loop' !0
         stop
      end if
      if(iprint.ge.19.and.myrank.eq.0) write(lud,*) 'icR0=',icR0
 

c***********************************************************************        
c     create a table mapping points outside the cluster back into it.
c***********************************************************************        
c     Note we use c1 from above!!        

      if(iprint.ge.19.and.myrank.eq.0) then                  !0
         write(lud,*) '  '                                  !0
         write(lud,*) 'Equivalency table r'                 !0
         write(lud,*) 'ic   igroup   icrequ(ic,ieqiv)'      !0
      end if                                                !0
      do ic=1,Nc
        do igroup=1,ngroup
c       Generate the equivalent points using the ngroup operations
          rv3=transform3D(Rc(:,ic),group(igroup),ixy,ixz,iyz,
     &                    isx,isy,isz)

! Map rv3 back into cluster
          rv2=matmul(c1,rv3)    !project rv3 onto a1,a2,a3
          do while (any(rv2.lt.zeror-epsilon).or.
     &		any(rv2.gt.oner-epsilon))
             do l=1,3
                if(rv2(l).lt.zeror-epsilon) then
                   rv3(:)=rv3(:)+a(l,:)
                else if(rv2(l).gt.oner-epsilon) then
                   rv3(:)=rv3(:)-a(l,:)
                end if
             end do
             rv2=matmul(c1,rv3) !project rv3 onto a1,a2,a3
          end do
c         Figure out which point rv3 is
          itest=0
          do jc=1,Nc
            if(sum(abs(rv3(:)-Rc(:,jc))).lt.epsilon) then ! rv3=Rc(:,ic)
               icrequ(ic,igroup)=jc
                itest=1
               exit
            end if
          end do
          if(itest.eq.0) then
            write(lud,*) 'No equivalent found in icrequ',ic,igroup
            stop
          end if 
          if(iprint.ge.19.and.myrank.eq.0) 
     &      write(lud,"(2x,i3,4x,i3,8x,i4)") ic,igroup,icrequ(ic,igroup) !0
        end do
      end do                               
      

c***********************************************************************        
c     Now create the tables for the differences of R
c***********************************************************************        
c     Note we use c1 from above!!   
      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) '  ' !0
        write(lud,*) 'differences in r' !0
        write(lud,*) '  ic   jc  icrdiff(ic,jc)' !0
      end if
      do i=1,Nc
      do j=1,Nc
         iv1(:)=Rc(:,i)-Rc(:,j)
         itest=0
         do n=-1,1        ! Look in the surrounding 
         do m=-1,1        ! 27 clusters for R by adding and
         do p=-1,1        ! subtracting a1, a2 and a3
            iv2=iv1 + n*a(1,:)+m*a(2,:)+p*a(3,:)
            rv2=matmul(c1,iv2)
            if(rv2(1).lt.oner-epsilon.and.rv2(1).gt.zeror-epsilon.and.
     &         rv2(2).lt.oner-epsilon.and.rv2(2).gt.zeror-epsilon.and.
     &         rv2(3).lt.oner-epsilon.and.rv2(3).gt.zeror-epsilon ) then
c              This point is within the cluster
c              Now we must find which point it is.
               do ic=1,Nc
                  if(sum(abs(iv2(:)-Rc(:,ic))).eq.0) then
                     icrdiff(i,j)=ic
                     itest=itest+1   ! complain if multiple solutions are found
                  end if
               end do
               if(iprint.ge.5.and.myrank.eq.0) 
     &              write(lud,"(2x,i3,4x,i3,6x,i4)") i,j,icrdiff(i,j) !0
            end if
         end do
         end do
         end do
         if(itest.ne.1) then
            write(lud,*) 'no icrdiff(i,j) found',i,j
            stop
         end if
      end do
      end do  
        
c***********************************************************************        
c     Now create the tables for the neighbors of cluster points
c***********************************************************************               
c     Note we use c1 from above, don't overwrite it!!  
c

      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) '  ' !0
        write(lud,*) 'neighbors of ic' !0
        write(lud,*) '  ic   neigh  nneigh(ic,jc)' !0
      end if
      neighbor(1,:)=(/ 1, 0, 0/) ! neighbor in the +x direction
      neighbor(2,:)=(/-1, 0, 0/) ! neighbor in the -x direction
      neighbor(3,:)=(/ 0, 1, 0/) ! neighbor in the +y direction
      neighbor(4,:)=(/ 0,-1, 0/) ! neighbor in the yy direction
      neighbor(5,:)=(/ 0, 0, 1/) ! neighbor in the +z direction
      neighbor(6,:)=(/ 0, 0,-1/) ! neighbor in the -z direction
      do i=1,Nc
      do j=1,2*ndim
         iv1(:)=Rc(:,i)-neighbor(j,:)
         itest=0
         do n=-1,1        ! Look in the surrounding 
         do m=-1,1        ! 27 clusters for R by adding and
         do p=-1,1        ! subtracting a1, a2 and a3
            iv2=iv1 + n*a(1,:)+m*a(2,:)+p*a(3,:)
            rv2=matmul(c1,iv2)
            if(rv2(1).lt.oner-epsilon.and.rv2(1).gt.zeror-epsilon.and.
     &         rv2(2).lt.oner-epsilon.and.rv2(2).gt.zeror-epsilon.and.
     &         rv2(3).lt.oner-epsilon.and.rv2(3).gt.zeror-epsilon ) then
c              This point is within the cluster
c              Now we must find which point it is.
               do ic=1,Nc
                  if(sum(abs(iv2(:)-Rc(:,ic))).eq.0) then
                     nneigh1(i,j)=ic
                     itest=itest+1   ! complain if multiple solutions are found
                  end if
               end do
               if(iprint.ge.19.and.myrank.eq.0) 
     &              write(lud,"(2x,i3,4x,i3,6x,i4)") i,j,nneigh1(i,j) !0
            end if
         end do
         end do
         end do

      end do
      end do  
      
!      if(iprint.ge.19.and.myrank.eq.0) then
!        write(lud,*) 'nearest neighbors of icR0'
!        do j=1,2*ndim
 !         write(lud,"(i2)") nneigh1(icR0,j)
 !       end do
  !    end if

ccccccccccccccccccccccccccccccccccccccccccccccc
c	Test the imperfection

      if(iprint.ge.19.and.myrank.eq.0) then
        write(42,*) '  ' !0
        write(42,*) 'neighbors of ic' !0
        write(42,*) 
     &    'shell  nneigh(shell)  nneighl(shell) nneighp(shell)' !0
      end if

      Nm=abs(determinant(a))

      i1 = 2*(Nm+1)**(1d0/ndim); rv1=0d0; nneigh=0; ishellmax=0
      nneighl(:)=0 
      nneighp(:)=0
      
      do i=-i1,i1
      do j=-i1,i1
      do k=-i1,i1
        iv1=(/ i, j, k/)
        i2=sum(abs(iv1)) ! the lattice neighbor index (according to Betts)
        if(i2.gt.0) then
          if(i2.le.nshellmax) then
            nneighl(i2)=nneighl(i2)+1
          end if
          do n=1,ndim
            rv1(n)=iv1(n) + 0.5d0*sum(a(:,n))
          end do
          rv2=matmul(c1,rv1)
          if(rv2(1).lt.1d0-epsilon.and.rv2(1).gt.0d0-epsilon.and.
     &       rv2(2).lt.1d0-epsilon.and.rv2(2).gt.0d0-epsilon.and.
     &       rv2(3).lt.1d0-epsilon.and.rv2(3).gt.0d0-epsilon ) then
c            This point is within the cluster. Now we need to find
c            the shortest distance between this point and icR0, and
c            the corresponding cluster neighbor index i3.
             i3=i2
             do n=-1,1
             do m=-1,1
             do p=-1,1
               iv2=iv1 + n*a(1,:) + m*a(2,:) + p*a(3,:)
               if(sum(abs(iv2)).lt.i2) i3=sum(abs(iv2))
             end do !p
             end do !m
             end do !n
             nneigh(i3)=nneigh(i3)+1
             if(i3.gt.ishellmax) ishellmax=i3
          end if
        end if
      end do !k
      end do !j
      end do !i
      
      i2=0
      do while (sum(nneighp).lt.Nc-1)
        i2=i2+1
        nneighp(i2)=nneighl(i2)
      nneighp(i2)=Nc-1-sum(nneighp(1:i2-1))
      end do
      
      if(sum(nneigh)+1.ne.Nc.or.sum(nneighp)+1.ne.Nc) then
        if(myrank.eq.0) write(lud,*) 'error in nneigh'
        if(myrank.eq.0) write(lud,*) 'sum(nneigh)=',sum(nneigh)
        if(myrank.eq.0) write(lud,*) 'sum(nneighp)=',sum(nneighp)
        stop
      endif
        
      if(iprint.ge.19.and.myrank.eq.0) then
        do i1=1,ishellmax
          write(lud,"(2x,i3,8x,i3,10x,i3,10x,i3)") 
     &      i1,nneigh(i1),nneighl(i1),nneighp(i1)
        enddo
        write(lud,*) 'ishellmax=',ishellmax
      endif
c     Now determine the perfection.  It seems that Betts defines the
c     ferromagnetic imperfection I_F as the numer of additional spins in
c     the shells beyond j plus the number missing the shells before j,
c     and not counting those in the jth shell.  The integer j seems to
c     be the number which yeilds the smallest I_F

      do j=1,ishellmax
        i3=0                                 ! count of imperfections
        do i1=1,j-1
          i3=i3+abs(nneigh(i1)-nneighl(i1))  ! number in shell i1 for cluster-lattice
        end do
        do i1=j+1,ishellmax                  ! cluster occupancy beyond j is an imperfection
          i3=i3+nneigh(i1)
        end do
        if(j.gt.1) then
          if(i3.lt.I_F) then
             I_F=i3 
             I_FShell=j
           end if
        else  ! j=1
          i_F=i3
          i_FShell=j
        end if
      end do

!	if(iprint.ge.19.and.myrank.eq.0)
 !    &  write(lud,*) 'Ferromagnetic imperfection I_F=',I_F
     
c     Betts proposes a simpler algorithm to calculate the perfection,
c     involving the difference of two quantities.
c
c     F_f = SUM k nneigh(k)
c            k
c
c     F_p = SUM k nneighp(k)
c            k

      i2=0; i3=0
      do i1=1,ishellmax
        i2=i2+i1*nneigh(i1)
        i3=i3+i1*nneighp(i1)
      end do
	if(iprint.ge.19.and.myrank.eq.0)
     &  write(lud,*) 'Betts imperfection I_F=',i2-i3
     

c***********************************************************************        
c     Create the cluster momentum tables Kc(i,ic), i=x,y,z
c***********************************************************************        

c     First find the points in the irreducible wedge.  They will be
c     indexed from 1 to Ncw
      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) '  Wedge K points  '    !0
        write(lud,*) 'ic  i  j  k  Kcx(ic)  Kcy(ic)  Kcz(ic)  ' !0
      end if
      ic=0; Ncw=0; icK0=0; icKpi=0
      do i=-Nc,Nc    !times a1
      do j=-Nc,Nc    !times a2
      do k=-Nc,Nc    !times a3
        iv1=(/i,j,k/)
        do l=1,3         !x,y,z components 
          rv1(l)=sum(iv1(:)*gvector(:,l))
        end do
        if(rv1(1).gt.-pi+epsilon.and.rv1(1).lt.pi+epsilon.and.
     &     rv1(2).gt.-pi+epsilon.and.rv1(2).lt.pi+epsilon.and.
     &     rv1(3).gt.-pi+epsilon.and.rv1(3).lt.pi+epsilon) then
c          The point is in the first SC Brillouin zone.  Now see if
c          the point falls within the IRW.  First see if it is 
c          equivalent to any known point in the IRW (if so it cannot be in IRW)
          Kcops: do igroup=2,ngroup ! exclude the identity
            rv2=transform3D(rv1,group(igroup),ixy,ixz,iyz,isx,isy,isz)
            itest=0
            do jc=1,Ncw
              if(sum(abs(Kc(:,jc)-rv2(:))).lt.epsilon) then ! rv2 = Kc(jc)
                 exit Kcops
              end if
            end do
          end do  Kcops
          if(igroup.gt.ngroup) then ! rv1 is in 1BZ and not equivalent to any point  
             Ncw=Ncw+1        ! in the known IRW.  It must be in the IRW
             ic=ic+1
             Kc(:,ic)=rv1(:)
             ig(ic)=i
             jg(ic)=j
             kg(ic)=k
             if(sum(abs(rv1-pi)).lt.epsilon) icKpi=ic
             if(sum(abs(rv1-zeror)).lt.epsilon) icK0=ic
             if(iprint.ge.19.and.myrank.eq.0) 
     &            write(lud,5) ic,i,j,k,Kc(1,ic),Kc(2,ic),Kc(3,ic) !0
          end if
        end if 
      end do 
      end do 
      end do 

c     Now find the points in the 1BZ but outside the irreducible wedge.  
c     They will be indexed from Ncw+1 to Nc
      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) ' Non Wedge K points '    !0
      end if
      do i=-Nc,Nc              
      do j=-Nc,Nc
      do k=-Nc,Nc
         iv1=(/i,j,k/)
         do l=1,3         !x,y,z components 
            rv1(l)=sum(iv1(:)*gvector(:,l)) ! K point in xyz
         end do
         if(rv1(1).gt.-pi+epsilon.and.rv1(1).lt.pi+epsilon.and.
     &      rv1(2).gt.-pi+epsilon.and.rv1(2).lt.pi+epsilon.and.
     &      rv1(3).gt.-pi+epsilon.and.rv1(3).lt.pi+epsilon) then
c           The point is in the first SC Brillouin zone.  Now see if
c           the point falls within the IRW.  
            do jc=1,Ncw
               if(sqrt(sum((Kc(:,jc)-rv1(:))**2)).lt.epsilon) then
c                 rv1 and Kc(ic) are the same point
                  exit
               end if
            end do
            if(jc.gt.Ncw) then !rv1 is in 1BZ and not in the IRW.
               ic=ic+1
               Kc(:,ic)=rv1(:)
               ig(ic)=i
               jg(ic)=j
               kg(ic)=k
               if(iprint.ge.19.and.myrank.eq.0) 
     &              write(lud,5) ic,i,j,k,Kc(1,ic),Kc(2,ic),Kc(3,ic) !0
            end if
         end if 
      end do 
      end do 
      end do 
        
 5        format(1x,i3,1x,i3,1x,i3,1x,i3,2x,f6.3,3x,f6.3,3x,f6.3)
      
      ntw=nl*Ncw      
      
      if(ic.ne.Nc) then
         if(myrank.eq.0) write(lud,*) 'Bug in Kc(ic) loop' !0
         stop
      end if
      if(iprint.ge.19.and.myrank.eq.0) 
     &  write(lud,"('icK0=',i3,' icKpi=',i3,' Ncw=',i3,' Nc=',i3)") 
     &  icK0,icKpi,Ncw,ic

c***********************************************************************
c     Create the indirect addressing table ict(i,j,k) which maps a 
c     point indexed by the principle translations in k-space to ic.
c***********************************************************************        
      N_sp=sqrt(float(Nc))
      do i=-N_sp,N_sp
      do j=-N_sp,N_sp
      do k=-N_sp,N_sp
         iv1=(/i,j,k/)
         do l=1,ndim        !x,y,z components 
            rv1(l)=sum(iv1(:)*gvector(:,l)) ! K point in xyz
         end do
! Map rv1 into -pi->+pi BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+epsilon) rv1=rv1-twor*pi
         where (rv1.lt.-pi+epsilon) rv1=rv1+twor*pi

         do jc=1,Nc
            if(sum(abs(rv1(:)-Kc(:,jc))).lt.epsilon) then    
               ict(i,j,k)=jc
               exit
            end if
         end do
         if(jc.gt.Nc) then
            write(lud,*) 'ict failed for i,j,k=',i,j,k
            write(lud,*) rv1
            stop
         end if
      end do
      end do
      end do

c***********************************************************************        
c     Create the table ickdiff(ic1,ic2), which takes the difference
c     between K-points
c***********************************************************************        
      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) '   '                           !0
        write(lud,*) 'ic1   ic2  ickdiff(ic1,ic2)  ' !0
      end if
      do ic1=1,Nc
      do ic2=1,Nc
         rv1=Kc(:,ic1)-Kc(:,ic2)
! Map rv1 into -pi->+pi BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+epsilon) rv1=rv1-twor*pi
         where (rv1.lt.-pi+epsilon) rv1=rv1+twor*pi
c        Now we need to figure out where this point is!
         do ic=1,Nc
            if(sum(abs(rv1(:)-Kc(:,ic))).lt.epsilon) then    
               ickdiff(ic1,ic2)=ic
               exit
            end if
         end do
         if(ic.gt.Nc) then
            write(lud,5) 'ickdiff failed for ic1,ic2=',ic1,ic2
            stop
         end if
         if(iprint.ge.19.and.myrank.eq.0) 
     &        write(lud,"(1x,i3,4x,i3,6x,i4)") ic1,ic2,ickdiff(ic1,ic2) !0
      end do
      end do
          
c***********************************************************************        
c     Create the table ickplus(ic1,ic2), which takes the sum
c     between to K-points
c***********************************************************************        
      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) '   '                           !0
        write(lud,*) 'ic1   ic2  ickplus(ic1,ic2)  ' !0
      end if
      do ic1=1,Nc
      do ic2=1,Nc
         rv1=Kc(:,ic1)+Kc(:,ic2)
! Map rv1 into -pi->+pi BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+epsilon) rv1=rv1-twor*pi
         where (rv1.lt.-pi+epsilon) rv1=rv1+twor*pi
c        Now we need to figure out where this point is!
         do ic=1,Nc
            if(sum(abs(rv1(:)-Kc(:,ic))).lt.epsilon) then    
               ickplus(ic1,ic2)=ic
               exit
            end if
         end do
         if(ic.gt.Nc) then
            write(lud,5) 'ickplus failed for ic1,ic2=',ic1,ic2
            stop
         end if
         if(iprint.ge.19.and.myrank.eq.0) 
     &        write(lud,"(1x,i3,4x,i3,6x,i4)") ic1,ic2,ickplus(ic1,ic2) !0
      end do
      end do

c***********************************************************************
c     Now find the points equivalent to ic, ickequ(ic,igroup) where 
c     igroup is the group operation.
c***********************************************************************        
c
c
       if(iprint.ge.19.and.myrank.eq.0) write(lud,*) '   '               !0
        if(iprint.ge.19.and.myrank.eq.0) write(lud,*) 'equivalency in k'  !0
        if(iprint.ge.19.and.myrank.eq.0) 
     &              write(lud,*) '   ic     igroup   ickequ(ic,igroup)'  !0

      do ic=1,Nc
         i=ig(ic)
         j=jg(ic)
         k=kg(ic)
         iv1=(/i,j,k/)
         do l=1,ndim
           rv2(l)=sum(iv1(:)*gvector(:,l))
         end do
c        Kc(ic) falls within the zone, by construction, now look 
c        for equivalent points
c     
         do igroup=1,ngroup
            rv1=transform3D(rv2,group(igroup),ixy,ixz,iyz,isx,isy,isz) 
c           These are the equivalent points.  Now map to 1st BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+epsilon) rv1=rv1-twor*pi
         where (rv1.lt.-pi+epsilon) rv1=rv1+twor*pi
c           Now we need to figure out where this point is.
            do jc=1,Nc
               if(sqrt(sum((rv1(:)-Kc(:,jc))**2)).lt.epsilon) then
                  ickequ(ic,igroup)=jc
                  exit
               end if
            end do
            if(iprint.ge.5.and.myrank.eq.0) 
     &        write(lud,"(2x,i3,4x,i3,6x,i4)") 
     &        ic,igroup,ickequ(ic,igroup)                  !0
            if(jc.gt.Nc) then
               if(myrank.eq.0) then
                 write(lud,*)'no equivalent K-point found' !0
                 write(lud,"('ic,igroup,group(igroup)',3(2x,i2))")
     &             ic,igroup,group(igroup)         !0
                 write(lud,"('K =',3(f6.3,1x))") rv2        !0
                 write(lud,"('K1=',3(f6.3,1x))") rv1       !0
               end if
               stop
            end if
         end do
      end do


c***********************************************************************        
c       Now generate a table which maps any K to an equivalent point
c       in the irreducible wedge (IW), and a table for the degeneracy 
c       of each point in the IW. 
c 
c       The use of the wedge can greatly reduce the storage required.
c       for example for Nc=8, only 4 of the 8 points are in the irreducible 
c       wedge.  Thus, all storage arrays may be reduced in size by a 
c       factor of 2.  As Nc-->oo, the saving factor approaches 48 for
c       clusters of high symmetry.  
c
c***********************************************************************        
c       
      
      ickdeg=0; irwgroup=0
      if(iprint.ge.19.and.myrank.eq.0) then
        write(lud,*) '   '                            !0
        write(lud,*) ' degeneracy in k'              !0
        write(lud,*) ' ic  ickmap(ic)  ickdeg(ic) irwgroup(ic)' !0
      end if
c     scan over all of the K points in the 1BZ
      do ic=1,Nc
         i=ig(ic)
         j=jg(ic)
         k=kg(ic)
         iv1=(/i,j,k/)
         do l=1,ndim
            rv2(l)=sum(iv1(:)*gvector(:,l))
         end do
c        Scan over all allowed point group operations, to find the
c        one that maps ic into the IRW
         groupops: do igroup=1,ngroup
           rv1=transform3D(rv2,group(igroup),ixy,ixz,iyz,isx,isy,isz) 
c          These are the equivalent points.  Now map to 1st BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+epsilon) rv1=rv1-twor*pi
         where (rv1.lt.-pi+epsilon) rv1=rv1+twor*pi
c          See if this point is in the IRW!
           do jc=1,Ncw
              if(sum(abs(rv1(:)-Kc(:,jc))).lt.epsilon) then
                 ickmap(ic)=jc
                 ickdeg(jc)=ickdeg(jc)+1
                 irwgroup(ic)=igroup
                 exit groupops
              end if
           end do
         end do groupops
         if(igroup.gt.ngroup) then
            if(myrank.eq.0) write(lud,*) 'mapping failed for' !0
            if(myrank.eq.0) write(lud,*) 'ic',ic !0
            if(myrank.eq.0) write(lud,"('K=',3(f6.3,1x))") rv2 !0
            stop
         end if
         
      end do !ic

      if(iprint.ge.19.and.myrank.eq.0) then 
        do ic=1,Nc    
          write(lud,"(2x,i3,6x,i3,9x,i4,9x,i4)") 
     &    ic,ickmap(ic),ickdeg(ic),irwgroup(ic)
        end do
      end if

c     Test the degeneracy table.  Since ickdeg holds the degeneracy
c     of points in the IRW, its sum must equal Nc.
      if(sum(ickdeg(1:Ncw)).ne.Nc) then
         if(myrank.eq.0) then
           write(lud,*) 'bug in ickdeg',sum(ickdeg(1:Ncw))
           write(lud,*) 'ickdeg=',ickdeg(1:Ncw)
         end if
         stop
      end if
        
c     Test the group table.  irwgroup(ic within wedge)=1, the identity
      do ic=1,Ncw
         if(irwgroup(ic).ne.1) then
            if(myrank.eq.0) then
              write(lud,*) 'irwgroup failed for wedge pt.'     !0
              write(lud,*) 'ic',ic,'irwgroup(ic)',irwgroup(ic) !0
              write(lud,"('K=(',3(f6.3,1x),')')") Kc(:,ic)      !0
            end if
            stop
         end if
      end do
     
c***********************************************************************
c     generate the kt-points in the Brillouin-zone of the 
c     super-lattice, i.e. in the Wigner-Seitz cell of the  
c     reciprocal space defined by the cluster K-points
c***********************************************************************

      ntot = nover
      Dkt=oner/real(ntot,kind)
      nwsc=0 ! number of k-points in Wigner-Seitz cell 
      if (ntot.gt.0) then
      do i=-ntot+1,ntot
        do j=-ntot+1,ntot
          do k=1,ntot
           iv1=(/i,j,k/)
           rv1(:)=Dkt*(iv1(:)-halfr)*pi
c  Check if this k-point is in Wigner-Seitz cell, i.e. if the closest
c  K-point is K=0 
            r1=sum(rv1(:)*rv1(:)) ! distance from K=0
            itest=0
            do jc=1,Nc
            if (jc.eq.icK0) cycle
             r2=rv1(1)-Kc(1,jc)
             r3=rv1(2)-Kc(2,jc)
             r4=rv1(3)-Kc(3,jc)
             if((r2-pi).gt.zeror) then; r2=r2-twor*pi; endif
             if((r2+pi).lt.zeror) then; r2=r2+twor*pi; endif
             if((r3-pi).gt.zeror) then; r3=r3-twor*pi; endif
             if((r3+pi).lt.zeror) then; r3=r3+twor*pi; endif
             if((r4-pi).gt.zeror) then; r4=r4-twor*pi; endif
             if((r4+pi).lt.zeror) then; r4=r4+twor*pi; endif
             r2=(r2*r2+r3*r3+r4*r4)
             if(r2.lt.r1) then
               itest=1
               exit
             endif
c
            enddo ! jc
            if(itest.eq.0) then 
             nwsc=nwsc+1
             kt(1:3,nwsc) = rv1(:)
             nwsc=nwsc+1 ! -k must also be part of the WS cell
             kt(1:3,nwsc) = -rv1(:)
            endif
          enddo ! k
        enddo ! j
       enddo ! i
       else ! needed for FSS calculation
          nwsc=1
          kt(:,1)=zeror
       endif
c
       if(myrank.eq.0.and.iprint.ge.19) then
          write(lud,*) 'Number of points in WS, nwsc', nwsc
          write(lud,*) 'Nc times nwsc', Nc*nwsc
          write(lud,*) 'Total number of points in BZ', (2*ntot)**3 
       endif

c***********************************************************************
c     Form Epsd(K) for calculating the partial density of states
c***********************************************************************     
      call allocate_PDOS
      estep = (eend-estart)/dble(e_div)
      PDOS(:,:) = 0.d0
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) '    '
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) 
     &  ' ic      Epsbar(ic)        temp_min      temp_max    min_Eps(ic)    max_Eps(ic)'
	Epsbar=zeror
      do ic=1,Nc
	temp_max = - halfr*(cos(Kc(1,ic)+kt(1,1))+
     &			cos(Kc(2,ic)+kt(2,1))+cos(Kc(3,ic)+kt(3,1))) 
     &             - tprime*(cos(Kc(1,ic)+kt(1,1))*cos(Kc(2,ic)+kt(2,1))
     &                      +cos(Kc(2,ic)+kt(2,1))*cos(Kc(3,ic)+kt(3,1))
     &                      +cos(Kc(1,ic)+kt(1,1))*cos(Kc(3,ic)+kt(3,1))-oner)
	temp_min = temp_max
         do i=1,nwsc
            rv1(1) = Kc(1,ic)+kt(1,i)
            rv1(2) = Kc(2,ic)+kt(2,i)
            rv1(3) = Kc(3,ic)+kt(3,i)
            Epsbar(ic)= Epsbar(ic)
c     &            -   halfr*sum(cos(rv1)) 
     &      -   halfr*(cos(rv1(1))+cos(rv1(2))+cos(rv1(3))) 
     &            - tprime*(cos(rv1(1))*cos(rv1(2))
     &                     +cos(rv1(2))*cos(rv1(3))
     &                     +cos(rv1(1))*cos(rv1(3))-oner)
c            Epsd = - halfr*sum(cos(rv1))
             Epsd = - halfr*(cos(rv1(1))+cos(rv1(2))+cos(rv1(3)))
     &            - tprime*(cos(rv1(1))*cos(rv1(2))
     &                     +cos(rv1(2))*cos(rv1(3))
     &                     +cos(rv1(1))*cos(rv1(3))-oner)
		e_index = NINT((Epsd-estart)/estep)
		if (e_index >= 0 .and. e_index <= e_div) then
		    PDOS(e_index,ic) = PDOS(e_index,ic) + 1.d0
		else
		    write(*,*) "WARNING: ENERGY EXCEEDS : ", Epsd, estart, eend
		endif
	  if (Epsd.gt.temp_max) temp_max = Epsd
	  if (Epsd.lt.temp_min) temp_min = Epsd   
         end do ! i
        Epsbar(ic)=Epsbar(ic)/real(nwsc,kind)
	dEps(ic) = (temp_max - temp_min)/real(nover,kind) ! the partial DOS also has nover ticks
	min_Eps(ic) = temp_min
	max_Eps(ic) = temp_max
       if(myrank.eq.0.and.iprint.ge.19) 
     &	write(lud,"(i6,4x,f12.8,4x,f12.8,4x,f12.8,4x,f12.8,4x,f12.8)") 
     &     ic, Epsbar(ic), temp_min, temp_max,min_Eps(ic), max_Eps(ic)
      	end do !ic

c     Check 
c      r1=sum(Epsbar); r2=sum(Epsbar*ickdeg)
c      if(abs(tprime).lt.epsilon.and.abs(r1).gt.epsilon
c     &   .or.
c     &	 abs(tprime).lt.epsilon.and.abs(r2).gt.epsilon) then
c        if(myrank.eq.0.and.iprint.ge.19) then
c          write(lud,*) 'sum Epsbar wrong.',r1,r2
c	end if
c        stop
c      end if
c***********************************************************************
c     Calculate the partial density of states
c
c                     Nc            
c       p_rho(K,w) =  --   SUM_kt delta( eps(K+kt) - w )
c                     N
c                       
c***********************************************************************  

      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) '    '
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) 
     &  '  ic         ener(ic)         sum(ener) n_fact' !n_fact is the number of points within each cell ic
	
	n_fact = zeror
	do ic=1,Nc
	do i=0,e_div
	n_fact(ic) = n_fact(ic) + PDOS(i,ic)*estep
	end do
	end do

	do ic=1,Nc
	do i=0,e_div
	PDOS(i,ic) = PDOS(i,ic)/(n_fact(ic)*Nc)
	write(53,*) estart+i*estep, PDOS(i,ic) 
	end do
	write(53,*) '   '
	end do

	do i=0,e_div
        write(53,*) estart+i*estep, sum(PDOS(i,:))
        end do
	close(53)

 	  write(4005,403) 
403      format('#    count         pdos')
	ener = zeror
	do ic =1,Nc
	do i=0,e_div
	write(4005,*) estart+i*estep,PDOS(i,ic)
		      ener(ic) = ener(ic) + PDOS(i,ic)*estep
	end do
c	write(4006,*) ic, ener(ic), sum(ener)
       if(myrank.eq.0.and.iprint.ge.19) 
     &	write(lud,"(i6,4x,f12.8,4x,f12.8,4x,f18.8)") 
     &   ic, ener(ic),  sum(ener), n_fact(ic)
	end do

      p_rho_GcKf = zeror
      do ic=1,Nc
        do n=-nwn+1,nwn-1
	  do i=0,e_div-1
	    if ((wn(n) - (estart+float(i)*estep))*(wn(n) - (estart+float(i-1)*estep)).le.0) then
		p_rho_GcKf(ic,n) = halfr*(PDOS(i,ic) + PDOS(i+1,ic))
	    end if
	  end do
	end do
      end do


c***********************************************************************
c     Form Epsbar_histo(K) for calculating the Histogram
c***********************************************************************     
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) '    '
      if(myrank.eq.0.and.iprint.ge.19) write(lud,*) 
     &  ' ic,Epsbar_histo(ic),Epsbar(ic),wn(Epsbar_histo(ic))'
c	Epsbar=zeror
c	do ic=1,Nc
c         do i=1,nwsc
c            rv1(1) = Kc(1,ic)+kt(1,i)
c            rv1(2) = Kc(2,ic)+kt(2,i)
c            rv1(3) = Kc(3,ic)+kt(3,i)
c            Epsbar(ic)= Epsbar(ic)
c     &      -   halfr*(cos(rv1(1))+cos(rv1(2))+cos(rv1(3))) 
c     &            - tprime*(cos(rv1(1))*cos(rv1(2))
c     &                     +cos(rv1(2))*cos(rv1(3))
c     &                     +cos(rv1(1))*cos(rv1(3))-oner)
c                         end do ! i
c                      Epsbar(ic)=Epsbar(ic)/real(nwsc,kind)
c                      do n=-nwn,nwn ! Frequency loop starts
c	           
c	       if (wn(n) - Epsbar(ic).lt. epsilon) Epsbar_histo(ic) = n
c          end do !n
c        if(myrank.eq.0.and.iprint.ge.19) 
c     &	write(lud,"(i6,4x,i8,4x,f14.8,4x,f14.8)") 
c     &     ic,Epsbar_histo(ic),Epsbar(ic),wn(Epsbar_histo(ic))
c	end do !ic

!   finding when wn=eps_bar(ic)
!  1.  define function to find min of:
          funct(:,:)=zeror
      	  do ic=1,Nc
	    do n=-nwn,nwn
              funct(ic,n)=wn(n)-Epsbar(ic)
    	    end do
    	  end do
!  2.  find min of this function:  
          
          do ic=1,Nc   
    	     min_funct(ic)=abs(funct(ic,-nwn))
    	     Epsbar_histo(ic)=-nwn 
    	      do n=-nwn,nwn
    	       if(min_funct(ic) .gt. abs(funct(ic,n))) then 
    	        min_funct(ic)=abs(funct(ic,n))
    	        Epsbar_histo(ic)=n      ! index at which wn=eps_bar
    	       end if	
              end do ! n
        if(myrank.eq.0.and.iprint.ge.19) 
     &	write(lud,"(i6,4x,i8,4x,f14.8,4x,f14.8)") 
     &     ic,Epsbar_histo(ic),Epsbar(ic),wn(Epsbar_histo(ic))
           end do 
  	  
    	  



c***********************************************************************        
c       Form the Fourier transform tables.  
c       R-->K,  FTCoefs_K_to_R
c       K-->R,  FTCoefs_K_to_R
c
c                    __          the factor of 1/N is required
c                 1  \           so that both G(R)~1/w and G(K)~1/w
c       G(R=0) = --- /  G(K)     for large w
c                N   --
c                    K
c
c***********************************************************************        
      do ic=1,Nc
         do ik=1,Nc
            r1=zeror
            do idim=1,ndim
               r1=r1+Kc(idim,ik)*Rc(idim,ic)
            end do
            FTCoefs_K_to_R(ik,ic)=exp(-ii*r1)/real(Nc,kind)
            FTCoefs_R_to_K(ik,ic)=exp(+ii*r1)
         end do
      end do

c	do ic = 1,Nc
c	do i = 1,ndim
c	write(235,*) ic,i,Kc(i,ic),Rc(i,ic)
c	end do
c	end do

c     Now test for consistency of Kc and Rc
      do ic=1,Nc
        r1= real(sum(FTCoefs_K_to_R(1:Nc,ic)))-kroniker(ic,icR0)
        r2=aimag(sum(FTCoefs_K_to_R(1:Nc,ic)))
        r3= real(sum(FTCoefs_K_to_R(ic,1:Nc)))-kroniker(ic,icK0)
        r4=aimag(sum(FTCoefs_K_to_R(ic,1:Nc)))
        
        if(abs(r1).gt.epsilon.or.abs(r2).gt.epsilon) then
          write(lud,*) 'Kc and Rc are inconsistent'
          write(lud,*) 'icr,r1,r2=',ic,r1,r2
          stop
        end if
        if(abs(r3).gt.epsilon.or.abs(r4).gt.epsilon) then
          write(lud,*) 'Kc and Rc are inconsistent, K=0'
          write(lud,*) 'ick,r3,r4=',ic,r3,r4
          stop
        end if
      end do
        end subroutine tables_3D_cubic
	!end 
