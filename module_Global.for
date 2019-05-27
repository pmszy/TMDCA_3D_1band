c**************************************************************************** 
c***************************************************************************
c***************************************************************************
c***************************************************************************
c***************************************************************************
c***************************************************************************
c***************************************************************************
	module Global
c****************************************************************************
c       This module should appear in all the main blocks of the code.  
c       It replaces the common block and declares and passes all the 
c       global variables and parameters (those globally accessible).
c       Note that many of the vaiables are allocatable.  They are
c       allocated/deallocated in the calls below for 1 and 2 particle
c       variables.
c****************************************************************************
        use mod_tprof
	implicit none
	save
c**************************************************************************** 
c       set some parameters.
c**************************************************************************** 
c	Character Parameters
c**************************************************************************** 
	character*5 :: cluster
        character(8)  :: date
        character(10) :: time
        character(5)  :: zone
c**************************************************************************** 

c**************************************************************************** 
c       Logical parameters.
c****************************************************************************
	logical, parameter :: iover = .false. !If true, use k_tilde momentum summation; if false, use partial density of states
        logical, parameter :: is_broyden = .false. ! If true use broyden to converge the self-energy
	logical, parameter :: gscript_pole = .false. !If true, use close to critical disorder where Gscript might diverge
        logical, parameter :: ityp = .false. ! if true, TMT for all methods is implemented for error bar calculations.
        logical, parameter :: tmt = .false. ! if true, use the orginal CTMT. Both tmt&tmt_new =off for DCA rho(K)=exp(<Ln(rho(K)>)
        logical, parameter :: tmt_new = .false. ! if true, use one of the new schemes below: see meas.f
        logical, parameter :: tmt1 = .false. ! if true,  new-vers1: rho_typ(K)=exp(<Ln(rho_ii)>)exp(<Ln(rho(K)/rho_ii)>)
        logical, parameter :: tmt2 = .false.      ! if true,  new-vers2: rho_typ(K)=exp(<Ln(rho_ii)>)<(rho(K)/rho_ii)>
        logical, parameter :: tmt3 = .false.    ! if true, rho_typ(K)=(1-a)*exp(<Ln(rho_ii) + a*<rho(K)
        logical, parameter :: tmt4 = .false.    ! if true, rho_typ(K)=exp(<Ln(rho_ii)>)
        logical, parameter :: isw = .false. ! if true, use the inhomogeneous grid else use homogeneous grid
        logical, parameter :: histo_enable = .false. ! Set true if you want to calculate histogram else, always false. It is very expensive
        logical, parameter :: errorbar = .false. ! Set true calculate error bar
!	Interpolation flags. Works only if interpolation is true and case either 1,2, or 3.
        logical, parameter :: RHOTYPICAL = .false. ! Set true calculate error bar
        logical, parameter :: interpolation = .false. !Interpolation works only if case is either 1, 2, or 3. Set below in interpolation_type
        logical, parameter :: spline = .false.      ! if true,  use spline
!        logical, parameter :: spline_interpolation = .false. ! Set true calculate error bar
!        logical, parameter :: star_interpolation = .false. ! Set true calculate error bar
!        logical, parameter :: trilinear_interpolation = .true. ! Set true calculate error bar

c**************************************************************************** 
c       Integer parameters.
c**************************************************************************** 
	integer, parameter :: kind=8
        integer, parameter :: nwnm=400,ntm=100,neta=20,ndim=3,nl=10,e_div=100
        integer, parameter :: ichim=0   ! if 1, measure chim (expensive)
        integer, parameter :: ichit=0   ! if 1, measure chit
        integer, parameter :: ibetts=1  ! If 1, use Betts Lattice
	integer, parameter :: mvlen = 60000 ! Broyden Parameter
	integer, parameter :: mbroylen = 200  ! Broyden Parameter 
	integer, parameter :: lwork=100**2  ! zgetri and zgetrf parameter 
	integer, parameter :: myrank=0,lud=42,iprint=19 
	integer, parameter :: imeas=500000 ! For histogram. Should be larged enough
        integer,dimension(8) :: values
	integer,parameter :: gsflag=0       ! gsflag=0 => P-H symmetry; sflag=1 => no P-H symmetry
        integer, parameter :: tables_type = 1 ! 1=ibetts,2=conventional,3=cubic
        integer, parameter :: interpolation_type = 1 ! 1=trilinear_interpolation,2=star_interpolation,3=spline_interpolation
!#############################################################################
!Remember to set ispline=1 in spectra code while using case 3, else, ispline=0
!#############################################################################
        integer, parameter :: disorder_type = 2! 1=binary,2=box,3=Gaussian(Normal),4=Lorentzian(Cauchy) ,5=Lognormal
c**************************************************************************** 
c       Real parameters.
c**************************************************************************** 
        real(kind), parameter :: tmax=100.0_kind, dpi=dacos(-1.0_kind)
        real(kind), parameter ::  epsilon=0.000001000000000_kind
        real(kind), parameter ::  depsi=0.000001000000000_kind
        real(kind), parameter ::  zeror=0.000000000000000_kind
        real(kind), parameter ::  halfr=0.500000000000000_kind
        real(kind), parameter ::   oner=1.000000000000000_kind
        real(kind), parameter ::   twor=2.000000000000000_kind
        real(kind), parameter :: threer=3.000000000000000_kind
        real(kind), parameter ::  fourr=4.000000000000000_kind
        real(kind), parameter ::     pi=3.141592653589793_kind, hbar=1.0545717260d-34
	real(kind), parameter :: estart=-3.850000000000000_kind, eend=3.850000000000000_kind
	real(kind), parameter :: afact = 0.200000000000000_kind ! Adjust according when using the pole procedure
        real(kind), parameter :: damp=0.005d0 ! damping for the spectra (cannot be larger than -damp)
c**************************************************************************** 
c       Complex parameters.
c****************************************************************************
        complex(kind), parameter :: ii=(0.0_kind,1.0_kind)
        complex(kind), parameter :: eta=(0.0_kind,0.0_kind) ! Broadening parameter
        complex(kind), parameter :: zeroc=(0.0_kind,0.0_kind)
        complex(kind), parameter :: onec=(1.0_kind,0.0_kind)
c**************************************************************************** 
c       Integer Variables Global
c****************************************************************************
        integer :: iseed,iter,niter,nover,nwn,Nc,nt,run,nrun,nf,ned,nwn_c,
     &		meas,nmeas,nslr,isi,iso,isd,isf,isv,isg,nfh,
     &		ncw,a(3,3),ngroup,icR0,icK0,icKpi,ntw,ick,icr,icum,
     &		N_sp,ntot,nwsc,ic00,icpipi,igroup,x_pole,
     &		a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,nproc,threads
c**************************************************************************** 
c       Integer Variables for Broyden
c****************************************************************************
	integer :: vlen,broylen
	integer, allocatable  :: ipiv(:)
c****************************************************************************
c        integer, allocatable  
c****************************************************************************

        integer, allocatable  :: neighbor(:,:),
     &		ickdeg(:),Rc(:,:),ict(:,:,:),
     &		ickmap(:),ick2map(:,:),Epsbar_histo(:),
     &		Rcx(:),Rcy(:),Rcz(:),
     &		icrdiff(:,:),ickdiff(:,:),
     &		ickplus(:,:),ickequ(:,:),icrequ(:,:),ipvt(:)
		   
c**************************************************************************** 
c        Real Variables
c****************************************************************************
        real(kind) :: ed,tprime,tperp,V0,dfac,ifac,test,testl,estep,
     &          delta,wc,gvector(3,3),Epsd,g1x,g1y,g1z,g2x,g2y,g2z,g3x,
     &		g3y,g3z,origx,origy,origz,wtime,ca,w_e,acoef,wmax,dw3,
     &		meas_time,put_time,update_time,geod_time,sumup_time,time1,time2


c**************************************************************************** 
c        Real Variables for spectra
c****************************************************************************
        real(8), allocatable :: g1_sp(:),g2_sp(:),g3_sp(:),Sig_sp_r(:,:,:),
     &                          Sig_sp_i(:,:,:),Derivs_sp_r(:,:,:),
     &                          Derivs_sp_i(:,:,:),Derivs_sp_r3(:,:,:,:),
     &				Derivs_sp_i3(:,:,:,:),Sig_sp_r3(:,:,:,:),Sig_sp_i3(:,:,:,:)

c**************************************************************************** 
c        Real Variables for Broyden
c****************************************************************************
        real(kind)  :: alpha,rms
        real(kind), allocatable  :: ub(:,:),vt_broyden(:,:),fb(:),
     &		df(:),vold(:),a_broyden(:,:),
     &		b_broyden(:,:),d_broyden(:,:),
     &		cm(:),wb(:),vector(:,:)
c**************************************************************************** 
c        Real Allocatables
c***************************************************************************
        real(kind), allocatable  :: V(:),kt(:,:),Kc(:,:),Kcz(:),Epsdbar(:),
     &		Epsbar(:),Kcx(:),Kcy(:),wn(:),dwn(:),weta(:),t(:),
     &		min_Eps(:),max_Eps(:),dEps(:),rhorun(:,:),gammat(:),rho_ij(:,:),rho_ii(:),
     &		chimeas(:),chit(:),Gamma(:,:),gammaf(:,:),gammafs(:,:),d_gamma(:,:),
     &          rho(:,:),rho_logf(:,:),rho_meas(:,:),Gammarun(:,:),
     &          rho_typ(:,:),rhotyprun(:,:),data_rho_local(:,:),logrhorun(:,:),
     &		rho_local(:,:),G_typ_real(:,:),logrho_meas(:,:),Gtyprun_real(:,:),
     &		rho_logfs(:,:),data_rho_typ(:,:),rho_geo(:),rho_geo1(:),
     &		rho_typ_local(:,:),p_rho(:,:),Gamma_data(:,:),
     &		rhof(:,:),rhofs(:,:),drhof(:,:),drho_log(:,:),PDOS(:,:),
     &          logrholoc_meas(:),logrholoc_f(:),logrholoc_run(:),tau(:), taut(:,:),
     &          rholoc_f(:),rholoc_fs(:),rholoc_typt(:),p_rho_GcKf(:,:),
     &		rholoc_meas(:),rholoc_run(:),rholoc_av(:),
     &          lrho_loc(:,:),rhof_rs_loc(:,:),histo_meas(:,:),histo_max(:),histo_min(:),
     &		rho_avrun(:,:),rho_av_f(:,:),rho_av_meas(:,:)

c**************************************************************************** 
c       Complex Variables
c**************************************************************************** 
c       Complex Allocatables
c****************************************************************************
        complex(kind), allocatable :: G0t(:),G0w(:),four(:,:),
     &		gmeas(:,:),GdKK_meas(:,:),V_temp(:,:),Gsc_temp(:,:,:), 
     &		data_G_typ(:,:),G_temp(:,:),Gtyprun(:,:),Gammarun_new(:,:),
     &		G_typf(:,:),Gbarrun(:,:),GcK_data(:,:),Gckf_data(:,:),
     &		gf(:,:),gfs(:,:),gft(:),gfR(:,:),G_typ(:,:),Gammarun_old(:,:),
     &		G(:,:,:),GcKf(:,:),work(:),Gcinv(:,:,:),Ginv(:,:,:),
     &		GcKfsc(:,:),chi(:,:),Gscriptrun(:,:),sigma(:,:),Gcsfsc(:,:),
     &		data(:,:),GdKK(:,:),Gdff(:,:),sigmarun(:,:),sigmaf(:,:),
     &		FTCoefs_K_to_R(:,:),FTCoefs_R_to_K(:,:),Gamma_new(:,:),gammaft(:),
     &		sigma_old(:,:),sigma_SR(:,:),b_data(:),d(:),Gamma_old(:,:),
     &          chicharge_meas(:,:,:),G_av_run(:,:),Gckf_av(:,:),G_av_f(:,:),G_av_meas(:,:),
     &		Gcinv_av(:,:,:),G_chg(:,:,:),Ginv_av(:,:,:),G_av(:,:,:),Gamma_old_av(:,:),
     &		Gamma_new_av(:,:),GcKfsc_av(:,:),Gc_av(:,:),Gdff_av(:,:),gft_av(:),gf_av(:,:),
     &		G_loc1(:,:),sigma_loc(:,:),sigma_temp(:,:),Gloc_temp(:,:),G_loc2(:,:),G_loc(:,:)

c****************************************************************************
c       Define some variables.
c****************************************************************************
c       64 bin precision (kind 8) is explicitly used unless otherwise state
c	throughout this code. To use it on 32 bit system, you have to use kind = 4
c**************************************************************************** 
c       Whenever possible, the following is used to define the
c       nomeclature of subscripts
c**************************************************************************** 
c
c       Gc              cluster greens function
c       K               momentum
c       s               space
c       f               frequency
c       t               euclidean time
c       sc              script (i.e. site excluded)
cc**************************************************************************** 
c       Table of Variables
cc**************************************************************************** 
c       ndim    	Number of spatial dimensions
c       ngroup  	Number of group operations which leave the lattice ivariant
c       nwnm    	Maximum number of frequencies
c       nwn     	number of frequencies used
c       Ncm     	Maximum size of the cluster
c       Nc      	size of the cluster
c       iseed   	seed integer the random number generator
c       iter    	present iteration number
c       niter   	total number of iterations to be done
c       run     	number of runs to be executed per iteration
c       meas    	number of measurements made per run
c       ntr     	number of sweeps before G is recalculated from {s(l)}
c       nslr    	number of sweeps sinse the la	st such recalculation
c       V(p)    	Disorder potential at site p
c       isi     	switch, if isi>=1, increase the # of runs by ifac each iter.
c       isd     	switch, if isd=1 use damping (dfac) in calc. of sigma
c	isf		switch, if isf=1 initalize with an externally provided s(i)
c	isv     	switch, if isv=1 (-1) binary (continuous) disorder 
c       pi      	3.14159
c       ed      	bare orbital energy (measured relative to the chem. pot.)
c       tprime 		next-near-neighbor hopping
c	test		convergence criteria
c       testl   	convergence criteria from last iteration
c       dfac    	damping factor (for sigma, see also isd)
c       ifac    	see isi
c       Kcx     	lookup table for Kcx values
c       Kcy     	lookup table for Kcy values
c       wn      	lookup table of frequencies
c       gf      	all runs accumulator for the single-particle green function
c       gfs     	all runs accumulator for the square single-particle green
c	gft		The damped single-particle green function (only used in put)
c       G       	green function
c       GcKf    	G(k,w) on the cluster
c       GcKfsc  	cluster-excluded G_sc(k,w)
c       Gcsf    	Cluster G(r,w)
c       chi     	temperary array 
c       ii      	sqrt(-1)
c       sigma   	self energy
c       data    	temporary array. All data_ are temporary arrays. E.g., data_rho_local.
c	p_rho		partial density of states
c	G_typ		single run accumulator for the single-particle for the typical green function
c	rho_typ		The K-resolved typical density of states
c	rho_meas	The K-resolved algebraic density of states
c	Gamma		The hybridization rate
c	drhof		Error in ADOS
c	drho_log	Error in TDOS
c	d_gamma		Error in gamma
c	_temp		These are all temporary arrays
c	Gcinv		Inverse of Gscript
c	Ginv		Inverse of G 
c	contains
	end module Global
c************************************************************
        subroutine allocate_arrays
c       allocate some arrays to the desired size.
c************************************************************
        use Global
        implicit none
c        include 'mpif.h'
c************************************************************
        integer info,infot,icount,ifail

        infot=0 
        icount=0
        allocate(ickdeg(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(neighbor(6,3),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(ickmap(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(ick2map(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(icrequ(Nc,e_div),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(icrdiff(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(ickdiff(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(ickplus(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(ickequ(Nc,e_div),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(Rcx(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(Rcy(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(ipvt(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(Rc(ndim,Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(Rcz(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(V(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(Gsc_temp(Nc,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
	allocate (rhof(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
	allocate (rho_ij(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
	allocate (rho_ii(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rhofs(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (drhof(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
	allocate (logrholoc_f(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rholoc_av(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rholoc_fs(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
	allocate (logrholoc_run(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(G_temp(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(V_temp(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(kt(ndim,420**ndim),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(Kc(ndim,Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Kcx(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Kcy(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Kcz(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Epsbar(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Epsbar_histo(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Epsdbar(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(wn(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(dwn(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (weta(neta),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (t(ntm),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (min_Eps(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (max_Eps(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (dEps(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (chimeas(neta),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (chit(ntm),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gamma(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gammarun(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gammarun_old(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gammarun_new(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (gammaf(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (gammafs(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (d_gamma(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gamma_data(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rho(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rho_logf(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rho_meas(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (logrho_meas(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (logrholoc_meas(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rholoc_meas(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rholoc_run(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rholoc_f(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (logrhorun(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rhorun(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rho_local(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rho_logfs(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (data_rho_local(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rho_typ(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (data_rho_typ(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rho_typ_local(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rhotyprun(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (p_rho(e_div,Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (G_typ_real(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (G0t(ntm),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (G0w(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (four(-nwn:nwn,ntm),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (gmeas(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (GdKK_meas(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (data_G_typ(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (G_typ(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (G_typf(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gtyprun(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gtyprun_real(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (GcK_data(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gbarrun(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (gf(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (gfR(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (gfs(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (gft(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (tau(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (taut(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (gammaft(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rholoc_typt(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (G(Nc,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gcinv(Nc,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Ginv(Nc,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (GcKf(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (GcKfsc(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gamma_new(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gamma_old(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gscriptrun(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (chi(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gcsfsc(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (data(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (GdKK(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gdff(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (sigma(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (sigmaf(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (sigmarun(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (sigma_old(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rho_geo(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rho_geo1(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (FTCoefs_K_to_R(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (drho_log(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (FTCoefs_R_to_K(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (ub(mvlen,mbroylen),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", " ub"
        allocate (vt_broyden(mvlen,mbroylen),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at vt_broyden"," info=",info
        allocate (fb(mvlen),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at fb"," info=",info
        allocate (df(mvlen),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at df"," info=",info
        allocate (vold(mvlen),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at vold"," info=",info
        allocate (a_broyden(mbroylen,mbroylen),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (b_broyden(mbroylen,mbroylen),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (d_broyden(mbroylen,mbroylen),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gckf_data(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (cm(mbroylen),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (wb(mbroylen),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (work(lwork),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (ipiv(mbroylen),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (vector(mvlen,2),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (b_data(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (d(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (p_rho_GcKf(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
         allocate (lrho_loc(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rhof_rs_loc(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (histo_meas(Nc,imeas),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (histo_max(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (histo_min(Nc),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount

        allocate (g1_sp(-N_sp:N_sp),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (g2_sp(-N_sp:N_sp),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Sig_sp_r(-N_sp:N_sp,-N_sp:N_sp,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Derivs_sp_r(-N_sp:N_sp,-N_sp:N_sp,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Sig_sp_i(-N_sp:N_sp,-N_sp:N_sp,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Derivs_sp_i(-N_sp:N_sp,-N_sp:N_sp,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (chicharge_meas(Nc,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at chicharge_meas" ," info=",info
        allocate (G_chg(Nc,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Ginv_av(Nc,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gcinv_av(Nc,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (G_av_run(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gckf_av(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (G_av_f(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (G_av_meas(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gamma_new_av(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gamma_old_av(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (G_av(Nc,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (GcKfsc_av(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gc_av(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rho_avrun(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rho_av_f(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (rho_av_meas(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gdff_av(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (gf_av(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (gft_av(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (G_loc1(1:Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (G_loc(1:Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (G_loc2(1:Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (sigma_loc(1:Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (sigma_temp(1:Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate (Gloc_temp(1:Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
	allocate(g3_sp(-N_sp:N_sp),stat=info)
	infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
	allocate(Sig_sp_r3(-N_sp:N_sp,-N_sp:N_sp,-N_sp:N_sp,-nwn:nwn),stat=info)
	infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
	allocate(Sig_sp_i3(-N_sp:N_sp,-N_sp:N_sp,-N_sp:N_sp,-nwn:nwn),stat=info)
	infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
	allocate(Derivs_sp_r3(-N_sp:N_sp,-N_sp:N_sp,-N_sp:N_sp,-nwn:nwn),stat=info)
	infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
	allocate(Derivs_sp_i3(-N_sp:N_sp,-N_sp:N_sp,-N_sp:N_sp,-nwn:nwn),stat=info)
	infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
        allocate(ict(-N_sp:N_sp,-N_sp:N_sp,-N_sp:N_sp),stat=info)
        infot=infot+info
        icount=icount+1; if(info.ge.1) write(42,*) "error at", icount
         if(infot.ne.0.and.myrank.eq.0) then
          write(lud,*) 'allocate_arrays',infot,' failures at step',ifail
          flush(lud)
          stop
        end if

        return
        end subroutine allocate_arrays

c************************************************************
        subroutine deallocate_arrays
c       deallocate the arrays allocated in subroutine allocate_arrays
c************************************************************
        use Global
        implicit none
c        include 'mpif.h'
c************************************************************
        deallocate(ickdeg)
        deallocate(ickmap)
        deallocate(ick2map)
        deallocate(icrequ)
        deallocate(icrdiff)
        deallocate(ickdiff)
        deallocate(ickplus)
        deallocate(ickequ)
        deallocate(Rcx)
        deallocate(Rcy)
        deallocate(Rcz)
        deallocate(V)
        deallocate(V_temp)
        deallocate(Gsc_temp)
        deallocate(G_temp)
        deallocate(kt)
        deallocate(Kc)
        deallocate(Kcx)
        deallocate(Kcy)
        deallocate(Kcz)
        deallocate(Epsbar)
	deallocate(Epsdbar)
        deallocate(Epsbar_histo)
        deallocate(wn)
        deallocate(dwn)
        deallocate(weta)
        deallocate(ipvt)
        deallocate(work)
        deallocate(t)
        deallocate(min_Eps)
        deallocate(max_Eps)
        deallocate(dEps)
        deallocate(chimeas)
        deallocate(chit)
        deallocate(Gamma)
        deallocate(Gammarun)
        deallocate(Gammarun_new)
        deallocate(Gammarun_old)
        deallocate(gammaf)
        deallocate(gammafs)
        deallocate(d_gamma)
        deallocate(Gamma_data)
        deallocate(rho)	
        deallocate(Rc)
        deallocate(rho_logf)
        deallocate(rho_typ)
        deallocate(rhotyprun)
        deallocate(rho_meas)
        deallocate(rho_local)
        deallocate(logrho_meas)
        deallocate(logrholoc_meas)
        deallocate(rholoc_meas)
        deallocate(rholoc_run)
        deallocate(rholoc_av)
        deallocate(logrhorun)
        deallocate(rhorun)
        deallocate(rho_logfs)
        deallocate(data_rho_typ)
	deallocate(rhof)
	deallocate(rhofs)
	deallocate(drhof)
	deallocate(logrholoc_f)
	deallocate(rholoc_f)
	deallocate(rholoc_fs)
	deallocate(logrholoc_run)
        deallocate(rho_typ_local)
        deallocate(p_rho)
        deallocate(G_typ_real)
        deallocate(Gtyprun_real)
        deallocate(Gtyprun)
        deallocate(sigmarun)
        deallocate(sigmaf)
        deallocate(G0t)
        deallocate(G0w)
        deallocate(four)
        deallocate(gmeas)
        deallocate(GdKK_meas)
        deallocate(data_G_typ)
        deallocate(G_typf)
        deallocate(G_typ)
        deallocate(Gckf_data)
        deallocate(data_rho_local)
        deallocate(gf)
        deallocate(gfR)
        deallocate(gfs)
        deallocate(gft)
        deallocate(taut)
        deallocate(tau)
	deallocate(rholoc_typt)
	deallocate(drho_log)
        deallocate(G)
        deallocate(GcKf)
        deallocate(Gbarrun)
        deallocate(GcKfsc)
        deallocate(Gamma_new)
	deallocate(neighbor)
        deallocate(Gamma_old)
        deallocate(Gscriptrun)
        deallocate(chi)
        deallocate(Gcsfsc)
        deallocate(data)
        deallocate(GdKK)
        deallocate(Gdff)
        deallocate(sigma)
        deallocate(sigma_old)
        deallocate(FTCoefs_K_to_R)
        deallocate(FTCoefs_R_to_K)
        deallocate(ub)
        deallocate(vt_broyden)
        deallocate(fb)
        deallocate(df)
        deallocate(vold)
        deallocate(a_broyden)
        deallocate(b_broyden)
        deallocate(d_broyden)
        deallocate(cm)
        deallocate(wb)
        deallocate(ipiv)
        deallocate(vector)
        deallocate(b_data)
        deallocate(d) 
        deallocate(Gcinv)
        deallocate(Ginv)
        deallocate(rho_geo)
        deallocate(rho_geo1)
        deallocate(p_rho_GcKf)
        deallocate(rho_ij)
        deallocate(rho_ii)
        deallocate(rhof_rs_loc)
        deallocate(lrho_loc)
        deallocate(histo_meas)
        deallocate(histo_max)
        deallocate(histo_min)
        deallocate (g1_sp)
        deallocate (g2_sp)
        deallocate (Sig_sp_r)
        deallocate (Derivs_sp_r)
        deallocate (Sig_sp_i)
        deallocate (Derivs_sp_i)
        deallocate (chicharge_meas)
        deallocate (G_chg)
	deallocate(Ginv_av)
	deallocate(Gcinv_av)
        deallocate (G_av_run)
        deallocate (Gckf_av)
        deallocate (G_av_meas)
	deallocate(Gamma_old_av)
	deallocate(Gamma_new_av)
        deallocate (G_av_f)
	deallocate(G_av)
	deallocate(GcKfsc_av)
	deallocate(Gc_av)
	deallocate(rho_avrun)
	deallocate(rho_av_f)
	deallocate(rho_av_meas)
	deallocate(Gdff_av)
	deallocate(gf_av)
	deallocate(gft_av)
	deallocate(G_loc)
	deallocate(G_loc1)
	deallocate(G_loc2)
	deallocate(sigma_loc)
	deallocate(sigma_temp)
	deallocate(Gloc_temp)
        deallocate(g3_sp)
        deallocate(Sig_sp_r3)
        deallocate(Sig_sp_i3)
        deallocate(Derivs_sp_r3)
        deallocate(Derivs_sp_i3)
	deallocate(ict)
        return
        end subroutine deallocate_arrays
c*********************************************************
	subroutine allocate_PDOS

c*********************************************************
	use Global
	implicit none
c*********************************************************
	integer info,infot

	infot=0

	allocate(PDOS(0:e_div,Nc),stat=info)
	infot = infot + info

	if(infot.ne.0.and.myrank.eq.0) then
	     write(lud,*) 'allocate_PDOS',infot,' failures'
	     flush(lud)
 !         stop
        end if
	      
	return
	end subroutine allocate_PDOS
	
c-------------------------------------------------------------------------

	subroutine deallocate_PDOS
c*************************************************************************
	use Global
	implicit none
c*************************************************************************
  
	deallocate(PDOS)
      
	return
	end subroutine deallocate_PDOS
c*****************************************************************************

c*************************************************************************
c	This subroutine is provided for the openmp calls and functions
c*************************************************************************
!  OpenMP runtime library to be used in conjunction with Open64 Compiler Suites.
!
!  Copyright (C) 2003 - 2009 Tsinghua University.
!
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 2.1 of the License, or (at your option) any later version.
!
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
!  
!  Contact information: HPC Institute, Department of Computer Science and Technology,
!  Tsinghua University, Beijing 100084, CHINA, or:
!
!  http://hpc.cs.tsinghua.edu.cn
        module omp_lib_kinds
        integer, parameter :: omp_integer_kind   = 4
        integer, parameter :: omp_logical_kind   = 4
        integer, parameter :: omp_lock_kind      = 8
        integer, parameter :: omp_nest_lock_kind = 8
        end module omp_lib_kinds


        module omp_lib
          use omp_lib_kinds
          integer, parameter :: openmp_version = 199910
              integer (kind=omp_integer_kind) :: ID
              integer (kind=omp_lock_kind) :: lck

          interface
            subroutine omp_destroy_lock (var)
              use omp_lib_kinds
              integer (kind=omp_lock_kind), intent(inout) :: var
            end subroutine omp_destroy_lock
          end interface

          interface
            subroutine omp_destroy_nest_lock (var)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind), intent(inout) :: var
            end subroutine omp_destroy_nest_lock
          end interface

          interface
            function omp_get_dynamic ()
              use omp_lib_kinds
              logical (kind=omp_logical_kind) :: omp_get_dynamic
            end function omp_get_dynamic
          end interface

          interface
            function omp_get_max_threads ()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) :: omp_get_max_threads
            end function omp_get_max_threads
          end interface

          interface
            function omp_get_nested ()
              use omp_lib_kinds
              logical (kind=omp_logical_kind) :: omp_get_nested
            end function omp_get_nested
          end interface

          interface
            function omp_get_num_procs ()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) :: omp_get_num_procs
            end function omp_get_num_procs
          end interface

          interface
            function omp_get_num_threads ()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) :: omp_get_num_threads
            end function omp_get_num_threads
          end interface

          interface
            function omp_get_thread_num ()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) :: omp_get_thread_num
            end function omp_get_thread_num
          end interface

          interface
            function omp_get_wtick ()
              use omp_lib_kinds
              double precision :: omp_get_wtick
            end function omp_get_wtick
          end interface

          interface
            function omp_get_wtime ()
              use omp_lib_kinds
              double precision :: omp_get_wtime
            end function omp_get_wtime
          end interface

          interface
            subroutine omp_init_lock (var)
              use omp_lib_kinds
              integer (kind=omp_lock_kind), intent(out) :: var
            end subroutine omp_init_lock
          end interface

          interface
            subroutine omp_init_nest_lock (var)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind), intent(out) :: var
            end subroutine omp_init_nest_lock
          end interface

          interface
            function omp_in_parallel ()
              use omp_lib_kinds
              logical (kind=omp_logical_kind) :: omp_in_parallel
            end function omp_in_parallel
          end interface

          interface
            subroutine omp_set_dynamic (enable)
              use omp_lib_kinds
              logical (kind=omp_logical_kind), intent(in) :: enable
            end subroutine omp_set_dynamic
          end interface

          interface
            subroutine omp_set_lock (var)
              use omp_lib_kinds
              integer (kind=omp_lock_kind), intent(inout) :: var
            end subroutine omp_set_lock
          end interface

          interface
            subroutine omp_set_nest_lock (var)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind), intent(inout) :: var
            end subroutine omp_set_nest_lock
          end interface

          interface
            subroutine omp_set_nested (enable)
              use omp_lib_kinds
              logical (kind=omp_logical_kind), intent(in) :: enable
            end subroutine omp_set_nested
          end interface

          interface
            subroutine omp_set_num_threads (nthreads)
              use omp_lib_kinds
              integer (kind=omp_integer_kind), intent(in) :: nthreads
            end subroutine omp_set_num_threads
          end interface

          interface
            function omp_test_lock (var)
              use omp_lib_kinds
              logical (kind=omp_logical_kind) :: omp_test_lock
              integer (kind=omp_lock_kind), intent(inout) :: var
            end function omp_test_lock
          end interface

          interface
            function omp_test_nest_lock (var)
              use omp_lib_kinds
              integer (kind=omp_integer_kind) :: omp_test_nest_lock
              integer (kind=omp_nest_lock_kind), intent(inout) :: var
            end function omp_test_nest_lock
          end interface

          interface
            subroutine omp_unset_lock (var)
              use omp_lib_kinds
              integer (kind=omp_lock_kind), intent(inout) :: var
            end subroutine omp_unset_lock
          end interface

          interface
            subroutine omp_unset_nest_lock (var)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind), intent(inout) :: var
            end subroutine omp_unset_nest_lock
          end interface
        end module omp_lib


