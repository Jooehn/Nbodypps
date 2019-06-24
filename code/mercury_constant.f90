!******************************************************************************
! MODULE: mercury_constant
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Global constant of the Mercury Code, numerical constant that fix
!! some numerical problems.
!
!******************************************************************************

module mercury_constant
  use types_numeriques

implicit none

integer, parameter :: CMAX = 1000 !< CMAX  = maximum number of close-encounter minima monitored simultaneously
integer, parameter :: NMESS = 1000 !< NMESS = maximum number of messages in message.in
integer, parameter :: NFILES = 5000 !< NFILES = maximum number of files that can be open at the same time
real(double_precision), parameter :: HUGE = 9.9d29 !< HUGE  = an implausibly large number
real(double_precision), parameter :: TINY = 4.D-15 !< A small number

!...   convergence criteria for danby
real(double_precision), parameter :: DANBYAC= 1.0d-14
real(double_precision), parameter :: DANBYB = 1.0d-13

!...    loop limits in the Laguerre attempts
integer, parameter :: NLAG1 = 1000
integer, parameter :: NLAG2 = 1000

end module mercury_constant
