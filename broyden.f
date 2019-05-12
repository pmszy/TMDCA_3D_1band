c**************************************************************************** 
      subroutine broyden (vector, vlen, alpha, rms, iter, broylen
     &     ,u,vt,f,df,vold,a,b,d,cm,w,ipiv,mbroylen,mvlen)
! ,ntasks)
c**************************************************************************** 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     History: Original code written by D.D. Johnson (see PRB 38, 12807)
!                  note:  there are a few typos in that paper but 
!                  the code is working!
!              Rewritten by W. A. Shelton for LSMS code 6/21/94
!                  this version is easy to read (no goto!!!! more comments ...)
!                  and is setup for MPP machines (not tested)
!              Rewritten by T. C. Schulthess, ORNL, March 97
!                  this version should work for any code (see comments below)
!
!     Bug fixes:   TCS, 8/5/97 see comments below 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     further comments on how to use this subroutine:
!     (if nothing useful stands here I had no time yet to update these
!     comments, please consult usage in lkkr code version 0.6.3 or later,
!     or call Thomas Schulthess (423) 5768067)
!
!     vector(r,i) -> i=1: old vector (input), scratch (ouput)
!                 -> i=2: new vector (input), mixed vector (output)
!     vlen        -> length of vector
!     alpha       -> linear mixing factor
!     rms         -> RMS difference between old and new vector
!     iter        -> iteration number (if 1, linear mixing, broyden reset)
!     broylen     -> number of iterations that are used for mixing (<=mbroylen)
!     u, vt, f, df, vold, and w are working arrays that need to be saved
!                   between call to this subroutine
!     a, b, d, cm, and ipiv are working arrays that need not be saved
!     mbroylen    -> maximum number of iterations that can be saved
!     mvlen       -> maximum length of vectors
!
!     See declaration for exact dimentsions and types
!
!     There are two options for matrix inversions, a Gaussian
!     elimination routine called invert1 and calls to lapack routines
!     with pivoting (see comments "using invert1" and "using lapack").
!     Obviously only one can be used, comment out the other one.
!
!     When using this subroutine in a parallel code in which only parts
!     of the vectors are known on every node, make sure that the calls
!     to gldsum (global sum) are correct (LKKR and LSMS codes have
!     different calls).
!
!     In a serial code, either comment out the calls to glbsum or
!     provide a dummy subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
!
      integer mbroylen,mvlen ! used for dimensioning
!     
      real*8 vector(mvlen,2)
      integer vlen
      real*8 alpha
      real*8 rms
      integer iter
      integer broylen
      real*8 u(mvlen,mbroylen)
      real*8 vt(mvlen,mbroylen)
      real*8 f(mvlen)
      real*8 df(mvlen)
      real*8 vold(mvlen)
      real*8 a(mbroylen,mbroylen)
      real*8 b(mbroylen,mbroylen)
      real*8 d(mbroylen,mbroylen)
      real*8 cm(mbroylen)
      real*8 w(mbroylen)
      integer ipiv(mbroylen)
!
      integer i,j,k,info,ntasks
      real*8 fac1,fac2,fnorm,dfnorm,w0,work,aij,gmi,cmj,wtmp
!
      integer lastit,lastm1,nn
      real*8 zero,one
      parameter (zero=0.d0,one=1.d0)
      real*8 amix
      save lastit,amix
!
      if (broylen .gt. mbroylen) then
         write(6
     *        ,'('' broyden: broylen='',i5,'' exeeds mbroylen='',i5)'
     *        ) broylen,mbroylen
         stop
      endif
!
      if( iter .eq. 1)then
c      if(.true.)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     first iteration: preform linear mixing, load f and vold, set
!                      different pointers and variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         lastit = 0             ! initialize pointers
         lastm1 = lastit-1
!
         amix = alpha           ! for safety reasons
         do k = 1,vlen
            f(k) = vector(k,2) - vector(k,1)
            vold(k) = vector(k,1)
         enddo
!
         do k = 1,vlen          ! this is the linear mixing
            vector(k,2) = vector(k,1) + amix * f(k)
         enddo
      return
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     iter > 1: this is where the non-linear mixing is done
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         lastit = lastit+1      ! update pointers
         lastm1 = lastit-1
!
         if( iter .gt. broylen) then ! set current lenght of broyden cycle
            nn = broylen
         else
            nn = lastit         !lastm1
         endif
!
         w0=.01d0               ! set weighting factor for the zeroth iteration
!
!---- find: f[i] := vector(2)[i] - vector(1)[i]
!           df   := f[i] - f[i-1]
!
         do k = 1,vlen
            df(k) = vector(k,2) - vector(k,1) - f(k)
         enddo
         do k = 1,vlen
            f(k) = vector(k,2) - vector(k,1)
         enddo
!
!---- find: fnorm  := |f|
!           dfnorm := |df|
!
         dfnorm = zero
         fnorm = zero
         do k = 1,vlen
            dfnorm = dfnorm + df(k)*df(k)
            fnorm  = fnorm  + f(k)*f(k)
         enddo
!	 if(ntasks.gt.1) then
!         call mp_sum_real( dfnorm,1,work)
!         call mp_sum_real( fnorm,1,work)
!
!	 endif
         dfnorm = sqrt( dfnorm )
         fnorm  = sqrt( fnorm )
!
!---- set: vector(2) := alpha*df/|df| + (vector(1) - vold)/|df|
!          vold := vector(1) 
!          vector(1) := df/|df|
!
         fac2 = one/dfnorm
         fac1 = amix*fac2
         do k = 1,vlen
            vector(k,2) = fac1*df(k) + fac2*(vector(k,1) - vold(k))
            vold(k) = vector(k,1)
            vector(k,1) = fac2*df(k)
         enddo
!
!---- store vector(1) and vector(2) in the stacks u and vt restpectively
!          v(i) = v(i+1)
!          u(i) = u(i+1)
!          u(nn) = vector(1)
!          v(nn) = vector(2)
!
         call broy_sav(u,vt,vector,iter-1,broylen,vlen,mvlen)
!
!---- calculate coefficient matrices, a(i,j), and sum cm(i) for corrections:
!         a(i,j<i) = v(j)*v(l)
!         a(i,i) = v(i)*v(i)
!         a(j,i) = a(i,j)
!         cm(i) = sum_(l=1,i) v(l)*f
!
         do j=1,nn - 1          ! off diagonal elements of a(i,j)
           do i = j+1,nn
              aij = zero
              do k = 1,vlen
                 aij = aij + vt(k,j)*vt(k,i)
              enddo
!	      if(ntasks.gt.1) call mp_sum_real(aij,1,work)
              a(i,j) = aij
              a(j,i) = aij
           enddo
        enddo
!
        do i = 1,nn             ! diagonal elements a(i,i) and cm(i)
           aij = zero
           cmj = zero
           do k=1,vlen
              cmj = cmj + vt(k,i)*f(k)
              aij = aij + vt(k,i)*vt(k,i)
           enddo
!	   if(ntasks.gt.1) then
!           call mp_sum_real( aij,1,work)
!           call mp_sum_real( cmj,1,work)
!	   endif
           a(i,i) = aij
           cm(i) = cmj
        enddo
!
!---- shift down weights in stack
!
!     (TCS, bug fixed 8/5/97: replace iter by iter-1 -> see broy_sav)
        if(iter-1 .gt. broylen)then
           do i=1,broylen-1
              w(i)=w(i+1)
           enddo
        endif
        wtmp = zero
        if( rms .gt. 1.0d-09 ) wtmp=2.0*sqrt(0.010d+00/rms)
        if( wtmp .lt. one )    wtmp=1.00d+00
        if(iter .gt. broylen)then
           w(broylen)=wtmp
        else
           w(lastit)=wtmp       !w(lastm1)=wtmp
        endif
!
!---- now calculate the b-matrix:
!        b = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> uses invert1
        do i=1,nn
           do j=1,nn
              d(j,i)= a(j,i)*w(j)*w(i)
              b(j,i)= zero
           enddo
           b(i,i)= one
           d(i,i)= w0**2 + a(i,i)*w(i)*w(i)
        enddo
!
        if(3*mbroylen.gt.mvlen) ! this is very unlikely
     *       stop 'broyden.f: need larger dimension for mvlen'
        call invert1( d, b, nn, 
     *       vector(1,1), vector(mbroylen+1,1), vector(2*mbroylen+1,1)
     *       ,mbroylen)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< uses invert1
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> uses lapack
!        do i=1,nn
!           do j=1,nn
!              b(j,i)= a(j,i)*w(j)*w(i)
!           enddo
!           b(i,i)= w0**2 + a(i,i)*w(i)*w(i)
!        enddo
!        call dgetrf(nn,nn,b,mbroylen,ipiv,info)
!        call dgetri(nn,b,mbroylen, ipiv, d, nn, info )
!        write(6,*) ' optimum lwork', d(1,1)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< uses lapack
!---- mix vectors: 
!        vector(2) = vold + amix*f - sum_i sum_j cm(j)*b(j,i)*w(j)*u(i)*w(i)
!
        do k=1,vlen
           vector(k,2)= vold(k) + amix * f(k)
        enddo
        do i=1,nn
           gmi = zero
           do j=1,nn
              gmi = gmi + cm(j)*b(j,i)*w(j)
           enddo
           do k=1,vlen
              vector(k,2) = vector(k,2) - gmi*u(k,i)*w(i)
           enddo
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      endif
!
      return
      end
!
!******************************************************************************
!
!     ==================================================================
      subroutine broy_sav(fins,fots,vector,itscf,istore,ivsiz,mivsiz)
!     ==================================================================
      integer itscf
      integer ivsiz
      integer mivsiz
      integer istore
!     ==================================================================
      real*8  vector(mivsiz,2)
      real*8  fins(mivsiz,istore)
      real*8  fots(mivsiz,istore)
!      write(6,'('' IN BROY_SAV: istore,itscf '',2i5)') istore,itscf
!     ==================================================================
      if( itscf .le. istore ) then
!     ==================================================================
!     Load the first istore iterations in increasing iteration count
!     ==================================================================
        do i = 1,ivsiz
          fins(i,itscf) = vector(i,2)
        enddo
!     ==================================================================
        do i = 1,ivsiz
          fots(i,itscf) = vector(i,1)
        enddo
!     ==================================================================
      else
!     ==================================================================
!     Re-load so that the ordering is in increasing iteration count
!     ==================================================================
        do j = 1,istore - 1
!          write(6,'('' IN BROY_SAV: j,j+1 '',2i5)') j,j+1
          do i = 1,ivsiz
            fins(i,j) = fins(i,j+1)
          enddo
!     ==================================================================
          do i = 1,ivsiz
            fots(i,j) = fots(i,j+1)
          enddo
!     ==================================================================
        enddo
!     ==================================================================
!     Load current charge densities in the last storage location
!     ==================================================================
        do i = 1,ivsiz
          fins(i,istore) = vector(i,2)
        enddo
!     ==================================================================
        do i = 1,ivsiz
          fots(i,istore) = vector(i,1)
        enddo
!     ==================================================================
      endif
!     ==================================================================
      return
      end
!
!******************************************************************************
!
! temporary for testting, should be lapack routine in future!
!
      subroutine invert1(a,b,m,td,ad,bd,mm)
!     =============================================================
      implicit real*8  (a-h,o-z)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      parameter (mm=5)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      dimension a(mm,mm),b(mm,mm)
      dimension td(mm),ad(mm),bd(mm)
!
! subroutine to preform gaussian elimination
!            no zeros along the diagonal
      save
!
      n=m
      if(n.gt.mm)then
       write(6,'('' invert: matrix a too large'')')
       stop
      endif
!
      do 14 i=1,n
      atmp=a(i,i)
      if(abs(atmp).lt. 1.0d-08)then
        write(6,'('' invert: matrix has zero diagonal'',  
     &            '' element in the '',i4,'' row'')')i
        stop
      endif
  14  continue
!
      if(n.eq.1) go to 605
!
      do 23 i=1,n
!
      do 35 j=1,n
 35      td(j)=a(j,i)/a(i,i)
!
         td(i)=0.0d+00
!
      do 71 k=1,n
         bd(k)=b(i,k)
 71      ad(k)=a(i,k)
!
      do 601 k=1,n
      do 601 j=1,n
         b(j,k)=b(j,k)-(td(j)*bd(k))
 601     a(j,k)=a(j,k)-(td(j)*ad(k))
!
 23   continue
!
      do 603 i=1,n
      do 603 j=1,n
 603     b(j,i)=b(j,i)/a(j,j)
!
      return
!
 605  b(1,1)=1.0d+00/a(1,1)
      return
      end
