!*************************************************************************
!*                                                                       *
!*                                                                       *
!*     BigClust - Stochastic Nonsmooth Optimization based                *
!*     Incremental Clustering Software (version 0.1)                     *
!*                                                                       *
!*     by Napsu Karmitsa 2023 (last modified 3.4.2024).                  *
!*                                                                       *
!*     Precision and parameters for nonsmooth optimization based         *
!*     clustering software.                                              *
!*                                                                       *
!*     The work was financially supported by the Research Council of     *
!*     Finland (Project No. #345804 and #345805).                        *
!*                                                                       *
!*     BigClust software is covered by the MIT license.                  *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     r_precision      ! Precision for reals.
!*     param            ! Parameters
!*     exe_time         ! Execution time.
!*

MODULE r_precision  ! Precision for reals.

  IMPLICIT NONE
  INTEGER, PARAMETER, PUBLIC :: prec = SELECTED_REAL_KIND(2*PRECISION(1.0))  ! old DOUBLE PRECISION
!  INTEGER, PARAMETER, PUBLIC :: prec = SELECTED_REAL_KIND(12) 

END MODULE r_precision


MODULE param  ! Parameters

  USE r_precision, ONLY : prec      ! Precision for reals.
  IMPLICIT NONE

! Intrinsic Functions
  INTRINSIC TINY,HUGE

! Parameters
  INTEGER, PARAMETER, PUBLIC :: maxeps = 20, maxnrs = 2000
  REAL(KIND=prec), PARAMETER, PUBLIC :: &
       zero    = 0.0_prec,    & ! 
       half    = 0.5_prec,    & ! 
       one     = 1.0_prec,    & ! 
       two     = 2.0_prec,    & !
       large   = HUGE(zero),  & !
       small   = TINY(zero)     ! 

END MODULE param


MODULE exe_time  ! Execution time.
  IMPLICIT NONE

  PUBLIC :: getime

CONTAINS
  SUBROUTINE getime(tim)  ! Execution time.
  IMPLICIT NONE
      
! Scalar Arguments
    REAL, INTENT(OUT):: tim  ! Current time, REAL argument.

! Intrinsic Functions
    INTRINSIC CPU_TIME

    CALL CPU_TIME(tim)

  END SUBROUTINE getime

END MODULE exe_time
