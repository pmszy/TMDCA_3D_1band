c**************************************************************************** 
        module mod_tprof
        implicit none

        integer, parameter :: lun6 = 6

        private
        save
c**************************************************************************** 

        logical, public :: timing_on = .FALSE.

        integer, parameter :: maxroutines = 1024
        integer, parameter :: maxlevels = 1024
        integer, parameter :: maxstrlen = 127

        real*8, dimension(maxroutines) :: dictstart, dicttotal
        integer, dimension(maxroutines) :: dictcount
        integer :: nroutine,nlevels

        character(len=maxstrlen), dimension(maxroutines) :: dictname
        character(len=maxstrlen), dimension(maxlevels) :: lastroutine

        public :: profstart, profend, profstat, profinit, dclock

!=======================================================================
        contains
!=======================================================================

        subroutine assert( lcond, msg, ivalue )
        logical, intent(in) :: lcond
        character(len=*), intent(in) :: msg
        integer, intent(in) :: ivalue

        if (.not.lcond) then
          write(lun6,*)  msg, ivalue
          stop 'assertion error '
        endif
        return
        end subroutine assert

!=======================================================================

        real*8 function dclock() 
        integer :: count,count_rate,count_max

        call system_clock(count,count_rate,count_max)
        if (count_rate.ne.0) then
           dclock = real(count)/real(count_rate)
        else
           dclock = 0.0
        endif

        end function dclock

!=======================================================================

        subroutine profinit()

        nroutine = 0

        dictname(:) = ' '
        dictstart(:) = 0.0
        dictcount(:) = 0.0
        dicttotal(:) = 0.0

        nlevels = 0
        lastroutine(:) = ' '

        timing_on = .TRUE.

        end subroutine profinit

!=======================================================================

        subroutine profstart(rname)
        character(len=*),intent(in) :: rname

        character(len=maxstrlen) :: name
        logical :: found,isok
        integer :: i,j,ipos

        name = rname
        nlevels = nlevels + 1

        isok = (1 .le. nlevels).and.(nlevels .le. maxlevels)
        call assert( isok,'** profstart: invalid nlevels ', nlevels )

        lastroutine(nlevels) = name
        found = .false.
        do j=1,nroutine
           i = nroutine - j + 1
           if (dictname(i)(1:1).eq.name(1:1)) then
                found = (dictname(i) .eq. name)
                if (found) then
                   ipos = i
                   exit
                endif
           endif
        enddo

        if (.not.found) then
          nroutine = nroutine + 1
          isok = (nroutine .le. maxroutines)
          call assert(isok,                                                  
     &       '** profstart: nroutine > maxroutines ', nroutine )

          ipos = nroutine
          dictname(ipos) = name
          dictcount(ipos) = 0
          dicttotal(ipos) = 0.0
        endif

        dictstart(ipos) = dclock()
        dictcount(ipos) = dictcount(ipos) + 1

        return
        end subroutine profstart

!=======================================================================

        subroutine profend(rname)
        character(len=*),intent(in) :: rname

        character(len=maxstrlen) :: name
        integer :: i,j,ipos
        logical :: found,isok
        
        real*8 :: tend

        name = rname
        tend = dclock()

        isok = (1.le.nlevels).and.(nlevels.le.maxlevels)
        call assert(isok,                                                    
     &      '** profend: invalid nlevels ', nlevels )

        isok = (name .eq. lastroutine(nlevels))
        if (.not.isok) then
          write(lun6,*) '** profend name != lastroutine(',nlevels,') '
          write(lun6,*) 'name: ', name
          write(lun6,*) 'lastroutine(nlevels): ', lastroutine(nlevels)
          stop '** error ** '
        endif

        found = .false.
        do j=1,nroutine
           i = nroutine - j + 1

           if (dictname(i)(1:1) .eq. name(1:1)) then
                found = (dictname(i) .eq. name)
                if (found) then
                        ipos = i
                        exit
                endif
           endif
        enddo

        if (.not.found) then
                write(lun6,*) '** profend: routine name not found '
                write(lun6,*) 'name: ',name
                stop '** error ** '
        endif

        dicttotal(ipos) = dicttotal(ipos) + (tend - dictstart(ipos));
        nlevels = nlevels - 1;

        return
        end subroutine profend

!=======================================================================

        subroutine profstat(outdev_in)
        implicit none
        integer, optional, intent(in):: outdev_in 
        character(len=maxstrlen) :: fname,fstr
        integer :: i, outdev
        logical, intrinsic :: present

        if (present(outdev_in)) then
           outdev = outdev_in
        else
           outdev = 16
        endif

        fname = 'profstat.dat'
        open(outdev, file=fname, form='formatted',                           
     &       access='sequential',status='unknown')
        rewind(outdev)

        fstr = "(A20,' was called ',i10,' times, total ',f10.2,' secs')"
        do i=1,nroutine
          write(outdev,fstr) dictname(i), dictcount(i), dicttotal(i)
          write(lun6,fstr) dictname(i), dictcount(i), dicttotal(i)
        enddo

        close(outdev)
        return
        end subroutine profstat

        end module mod_tprof
