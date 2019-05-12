c**************************************************************************** 
        subroutine ginit
c**************************************************************************** 
c       This program makes the initial greens functions, ie. those
c       when all s=+1.  It starts with go (the greens function
c	provided by geodgen corresponding to all s=0), and then uses 
c	the greens function updating block upg to form the desired result.
 
 
c************************************************************
	use Global
	implicit none
c************************************************************

        integer i,n,ic


c	Initialize some variables and accumulators
c       initialize the disorder potentials
	if(iter.eq.1.and.isf.eq.0) then
	  if(isv.eq.1) then
	    do i=1,Nc
              V(i)=V0*sign(1.0,ca-ran(iseed))
            end do
	  else
	    do i=1,Nc
	      V(i)=2.0*V0*(ran(iseed)-0.5)
	    end do
	  end if
	else
c	  use the field configuration from the last run or that
c	  was readin from sigma.dat (if isf=1).
	  write(6,*) 'initialized with old potentials'
	end if

c	Initialize some variables and accumulators

        return
        end
