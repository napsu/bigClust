!*************************************************************************
!*                                                                       *
!*     BigClust - Stochastic Nonsmooth Optimization based                *
!*     Incremental Clustering Software (version 0.1)                     *
!*                                                                       *
!*     by Napsu Karmitsa 2024 (last modified 3.4.2024).                  *
!*                                                                       *
!*     Computation of the value of the clustering problem and the        *
!*     correbonding subgradients. This file is specially desinged        *
!*     to solve clustering problems with stochastic LMB used within      *
!*     BigClust software.                                                * 
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
!*     obj_fun             ! Computation of clustering and auxiliary 
!*                         ! clustering functions and their subgradients
!*

MODULE obj_fun             ! Computation of the value and the subgradient of the 
                           ! objective function.

    USE r_precision, ONLY : prec  ! Precision for reals.
    IMPLICIT NONE

    PUBLIC :: &
        myf,           &   ! Computation of the value of the objective 
                           ! (building block between SLMB and clustering problem).
        myg,           &   ! Computation of the subgradient of the objective.
                           ! (building block between SLMB and clustering problem).
        auxfunc,       &   ! S Computation of aux clustering problem.
        clusterfunc,   &   ! S Computation of clustering problem.
        fgrad,         &   ! S Computation of the subgradient of the (aux) clustering problem.
        clustfull,     &   ! S Computation of value of clustering function with full batch but partial data
        clustfull2         ! S Computation of value of clustering function with full batch and whole data

CONTAINS
    
    !=============================================
    ! Computation of the value of the objective
    !=============================================
         
    SUBROUTINE myf(n,x,f,iterm)

        USE param, ONLY : large                   ! Parameters.
        USE initclust, ONLY : ns                  ! Switch for auxiliary and real clustering problem.

        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(n), INTENT(IN) :: &
            x  ! Vector of variables.

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(OUT) :: f  ! Value of the function.
        INTEGER, INTENT(IN) :: n           ! Number of variables.
        INTEGER, INTENT(OUT) :: iterm      ! Cause of termination:
                                           !   0  - Everything is ok.
                                           !  -3  - Failure in function calculations

        iterm = 0

        ! Function evaluation
        IF (ns == 1) THEN  ! Auxiliary problem
            CALL auxfunc(x,f)

        ELSE  ! Clustering problem
            CALL clusterfunc(x,f)
        END IF

        ! Error checking.
        IF (f > large) iterm = -3  !
        IF (f < -large) iterm = -3 !

        RETURN
      
    END SUBROUTINE myf

    !=================================================
    ! Computation of the subgradient of the objective
    !=================================================
     
    SUBROUTINE myg(n,x,g,iterm)

        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(n), INTENT(IN) :: x  ! Vector of variables.
        REAL(KIND=prec), DIMENSION(n), INTENT(OUT) :: g ! Subgradient.

        ! Scalar Arguments
        INTEGER, INTENT(IN) :: n                        ! Number of variables.
        INTEGER, INTENT(OUT) :: iterm                   ! Cause of termination:
                                                        !   0  - Everything is ok.
                                                        !  -3  - Failure in subgradient calculations
                                                        !        (assigned by the user).

        iterm = 0

        ! Gradient evaluation.
        CALL fgrad(x,g)

        RETURN

    END SUBROUTINE myg

    !=============================================
    ! Computation of aux clustering problem
    !=============================================

    SUBROUTINE auxfunc(x,fval)

        USE param, ONLY : zero, one, two, large
        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a => a_used, & ! Data matrix
            nft, &         ! Number of features
            nc, &          ! Current number of clusters.
            nrecord => nrec_used, &     ! Number of corts in data.
            dminim, &      ! dminim(nrecord), Distance of a(i) and the nearest centroid.
            xbest, &       ! xbest(nclust*nft), Best solution obtained.
            tnorm, &       ! The total number of norms computed.
            list1          ! list1(i) gives the cluster where point i belongs.
        USE initslmb, ONLY : &     
            nbatch, &      ! Size of the batch
            ib, &          ! The first index of the batch
            batchind       ! Batch indices


        IMPLICIT NONE


        ! Subroutine arguments
        REAL(KIND=prec), DIMENSION(maxdim), INTENT(IN) :: &
            x             ! Vector of variables.
        REAL(KIND=prec), INTENT(OUT) :: &
            fval          ! Value of the auxiliary clustering function.

        ! Local variables
        REAL(KIND=prec) ::  &
            dprev(nc),          &
            f1,f10,f2,f3
        INTEGER :: &
            i,j,k,i1,i2

        INTRINSIC ::  MIN
        
        DO k=1,nc-1
            tnorm=tnorm+one
            dprev(k)=zero
            DO j=1,nft
                dprev(k)=dprev(k)+(x(j)-xbest(j+(k-1)*nft))**2
            END DO
        END DO

        fval=zero
        loop: DO i2=1,nbatch 
            i=batchind(i2+ib) 
            f1=dminim(i)
            f10=4.0_prec*f1
            i1=list1(i)
            IF(dprev(i1) >= f10) THEN
                fval=fval+f1
                CYCLE loop
            END IF
            tnorm=tnorm+one
            f2=zero
            DO j=1,nft
                f2=f2+(a(j,i)-x(j))**2
            END DO
            f3=MIN(f1,f2)
            fval=fval+f3
        END DO loop
        RETURN

    END SUBROUTINE auxfunc

    !=============================================
    ! Computation of clustering problem
    !=============================================

    SUBROUTINE clusterfunc(x,f)

        USE param, ONLY : zero, one, two, large
        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a => a_used, & ! Data matrix.
            nft, &         ! Number of features.
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            nrecord => nrec_used, &     ! Number of recorts in data.
            tnorm          ! The total number of norms computed.
        USE initslmb, ONLY : &
            nbatch, &      ! Size of the batch
            ib, &          ! The first index of the batch
            batchind       ! Batch indices

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec), DIMENSION(maxdim), INTENT(IN) :: &
            x             ! Vector of variables.
        REAL(KIND=prec), INTENT(OUT) :: &
            f             ! Value of the clustering function.
        
        ! Local variables
        REAL(KIND=prec) ::  &
            dcent1(nc,nc), &
            f10,f2,f3

        INTEGER :: &
            i,i2,j,k,kmin

        DO i=1,nc
            dcent1(i,i)=zero
        END DO

        DO i=1,nc
            DO j=i+1,nc
                dcent1(i,j)=zero
                DO k=1,nft
                    dcent1(i,j)=dcent1(i,j)+(x(k+(i-1)*nft)-x(k+(j-1)*nft))**2
                END DO
                dcent1(j,i)=dcent1(i,j)
            END DO
        END DO        

        f=zero
        kmin=1
        DO i2=1,nbatch
            i=batchind(i2+ib) ! take nbatch random indices
    
            f2=large
            f10=large
            loop1: DO k=1,nc
                IF(k >= 2) THEN
                    IF(dcent1(k,kmin) >= f10) THEN
                        CYCLE loop1
                    END IF
                END IF                     
                tnorm=tnorm+one
                f3=zero
                DO j=1,nft
                    f3=f3+(a(j,i)-x(j+(k-1)*nft))**2
                END DO
                IF(f3 < f2) THEN
                    f2=f3
                    f10=4.0_prec*f2
                    kmin=k
                END IF
            END DO loop1
            f=f+f2
        END DO

        RETURN

    END SUBROUTINE clusterfunc

    !==============================================================
    ! Computation of the gradient
    !==============================================================

    SUBROUTINE fgrad(x,grad)

        USE param, ONLY : zero, one, two, large
        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a => a_used, & ! Data matrix.
            nft, &         ! Number of features.
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            nrecord => nrec_used, &     ! Number of recorts in data.
            dminim, &      ! dminim(nrecord), Distance of a(i) and the nearest centroid.
            xbest, &       ! xbest(nc*nft), Best solution obtained.
            list1, &       ! list1(i) gives the cluster where point i belongs.
            tnorm, &       ! The total number of norms computed.
            ns, &          ! Switch for auxiliary and real clustering problem.
            m              ! Number of variables in optimization:
                           !   m = nft    if ns = 1,
                           !   m = nft*nc if ns = 2.
        USE initslmb, ONLY : &
            nbatch, &      ! Size of the batch
            ib, &          ! The first index of the batch
            batchind       ! Batch indices


        IMPLICIT NONE

        REAL(KIND=prec), DIMENSION(maxdim), INTENT(IN) :: &
            x              ! Vector of variables.
        REAL(KIND=prec), DIMENSION(maxdim), INTENT(OUT) :: &
            grad           ! Subgradient.

        ! Array arguments
        REAL(KIND=prec) ::  &
            dcent1(nc,nc), &
            dprev(nc)

        ! Local variables
        REAL(KIND=prec) ::  &
            f1,f10,f12,f3
        INTEGER :: &
            i,j,k,jmin,i1,i2

        DO i=1,m
            grad(i)=zero
        END DO

        IF (ns == 1) THEN
            DO k=1,nc-1
                tnorm=tnorm+one
                dprev(k)=zero
                DO j=1,nft
                    dprev(k)=dprev(k)+(x(j)-xbest(j+(k-1)*nft))**2
                END DO
            END DO
        
            loop1: DO i2=1,nbatch
                i=batchind(i2+ib) ! take nbatch random indices

                f1=dminim(i)
                f10=4.0_prec*f1
                i1=list1(i)
                IF(dprev(i1) >= f10) THEN
                    CYCLE loop1
                END IF

                tnorm=tnorm+one
                f3=zero
                DO j=1,nft
                    f3=f3+(a(j,i)-x(j))**2
                END DO
                IF (f3 < dminim(i)) THEN
                    DO j=1,m
                        grad(j)=grad(j)+two*(x(j)-a(j,i))
                    END DO
                END IF
            END DO loop1

        ELSE
        
            DO i=1,nc
                DO j=i+1,nc
                    dcent1(i,j)=zero
                    DO k=1,nft
                        dcent1(i,j)=dcent1(i,j)+(x(k+(i-1)*nft)-x(k+(j-1)*nft))**2
                    END DO
                    dcent1(j,i)=dcent1(i,j)
                END DO
            END DO        

            jmin=1
            DO i2=1,nbatch
                i=batchind(i2+ib) ! take nbatch random indices

                f10=large
                f12=large
                loop2: DO k=1,nc
                    IF(k >= 2) THEN
                        IF(dcent1(k,jmin) >= f10) THEN
                            CYCLE loop2
                        END IF
                    END IF
                    tnorm=tnorm+one
                    f3=zero
                    DO j=1,nft
                        f3=f3+(a(j,i)-x(j+(k-1)*nft))**2
                    END DO
                    IF (f12 > f3) THEN
                        f12=f3
                        f10=4.0_prec*f12
                        jmin=k
                    END IF
                END DO loop2
                DO j=1,nft
                    grad(j+(jmin-1)*nft)=grad(j+(jmin-1)*nft)+two*(x(j+(jmin-1)*nft)-a(j,i))
                END DO
            END DO
        END IF

    END SUBROUTINE fgrad

    !===================================================
    ! Computation of clustering problem with full batch
    !===================================================

    SUBROUTINE clustfull(x,f)

        USE param, ONLY : zero, one, two, large
        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a => a_used, & ! Data matrix.
            nft, &         ! Number of features.
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            nrecord => nrec_used, &     ! Number of recorts in data.
            tnorm          ! The total number of norms computed.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec), DIMENSION(maxdim), INTENT(IN) :: &
            x             ! Vector of variables.
        REAL(KIND=prec), INTENT(OUT) :: &
            f             ! Value of the clustering function.
        
        ! Local variables 
        REAL(KIND=prec) ::  &
            dcent1(nc,nc), &
            f10,f2,f3

        INTEGER :: &
            i,j,k,kmin

        DO i=1,nc
            dcent1(i,i)=zero
        END DO

        DO i=1,nc
            DO j=i+1,nc
                dcent1(i,j)=zero
                 DO k=1,nft
                    dcent1(i,j)=dcent1(i,j)+(x(k+(i-1)*nft)-x(k+(j-1)*nft))**2
                END DO
                dcent1(j,i)=dcent1(i,j)
            END DO
        END DO        

        f=zero
        kmin=1
        DO i=1,nrecord
    
            f2=large
            f10=large
            loop1: DO k=1,nc
                IF(k >= 2) THEN
                    IF(dcent1(k,kmin) >= f10) THEN
                        CYCLE loop1
                    END IF
                END IF                     
                tnorm=tnorm+one
                f3=zero
                DO j=1,nft
                    f3=f3+(a(j,i)-x(j+(k-1)*nft))**2
                END DO
                IF(f3 < f2) THEN
                    f2=f3
                    f10=4.0_prec*f2
                    kmin=k
                END IF
            END DO loop1
            f=f+f2
        END DO

        RETURN

    END SUBROUTINE clustfull

    !=================================================================
    ! Computation of clustering problem with full batch and full data
    !=================================================================

    SUBROUTINE clustfull2(x,f)

        USE param, ONLY : zero, one, two, large
        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a, &           ! Data matrix.
            nft, &         ! Number of features.
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            nrecord, &     ! Number of recorts in data.
            tnorm          ! The total number of norms computed.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec), DIMENSION(maxdim), INTENT(IN) :: &
            x              ! Vector of variables.
        REAL(KIND=prec), INTENT(OUT) :: &
            f              ! Value of the clustering function.
        
        ! Local variables 
        REAL(KIND=prec) ::  &
            dcent1(nc,nc), &
            f10,f2,f3

        INTEGER :: &
            i,j,k,kmin

        DO i=1,nc
            dcent1(i,i)=zero
        END DO

        DO i=1,nc
            DO j=i+1,nc
                dcent1(i,j)=zero
                 DO k=1,nft
                    dcent1(i,j)=dcent1(i,j)+(x(k+(i-1)*nft)-x(k+(j-1)*nft))**2
                END DO
                dcent1(j,i)=dcent1(i,j)
            END DO
        END DO        

        f=zero
        kmin=1
        DO i=1,nrecord
    
            f2=large
            f10=large
            loop1: DO k=1,nc
                IF(k >= 2) THEN
                    IF(dcent1(k,kmin) >= f10) THEN
                        CYCLE loop1
                    END IF
                END IF                     
                tnorm=tnorm+one
                f3=zero
                DO j=1,nft
                    f3=f3+(a(j,i)-x(j+(k-1)*nft))**2
                END DO
                IF(f3 < f2) THEN
                    f2=f3
                    f10=4.0_prec*f2
                    kmin=k
                END IF
            END DO loop1
            f=f+f2
        END DO

        RETURN

    END SUBROUTINE clustfull2


END MODULE obj_fun
