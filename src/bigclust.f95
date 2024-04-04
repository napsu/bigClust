!*************************************************************************
!*                                                                       *
!*     BigClust - Stochastic Nonsmooth Optimization based                *
!*     Incremental Clustering Software (version 0.1)                     *
!*                                                                       *
!*     by Napsu Karmitsa 2024 (last modified 3.4.2024).                  *
!*                                                                       *
!*     The work was financially supported by the Research Council of     *
!*     Finland (Project No. #345804 and #345805).                        *
!*                                                                       *
!*     BigClust software is covered by the MIT license.                  *
!*                                                                       *
!*************************************************************************
!*
!*
!*     Codes included:
!*
!*     bigclust.f95          - Mainprogram for clustering software
!*                             (this file).
!*     parameters.f95        - Parameters. Inludes modules:
!*                               - r_precision - Precision for reals,
!*                               - param - Parameters,
!*                               - exe_time - Execution time.
!*     initbigclust.f95      - initialization of clustering parameters and SLMB.
!*                             Includes modules:
!*                               - initclust - Initialization of parameters for clustering.
!*                               - initslmb  - Initialization of SLMB.
!*     clusteringmod.f95     - Subroutines for clustering software.
!*     slmb.f95              - SLMB - Stochastic limited memory bundle method.
!*     objfun.f95            - Computation of the cluster function and subgradients values.
!*     subpro.f95            - subprograms for SLMB.
!*
!*     Makefile              - makefile.
!*
!*
!*     To use the software modify initbigclust.f95 as needed.
!*
!*
!*     References:
!*
!*     for BigClust:
!*       N. Karmitsa, V.-P. Eronen, M.M. Mäkelä, A. Airola, T. Pahikkala, "Stochastic limited 
!*       memory bundle method for clustering in big data", 2024
!*
!*     for original LMBM:
!*       N. Haarala, K. Miettinen, M.M. Mäkelä, "Globally Convergent Limited Memory Bundle Method  
!*       for Large-Scale Nonsmooth Optimization", Mathematical Programming, Vol. 109, No. 1,
!*       pp. 181-205, 2007. DOI 10.1007/s10107-006-0728-2.
!*
!*       M. Haarala, K. Miettinen, M.M. Mäkelä, "New Limited Memory Bundle Method for Large-Scale 
!*       Nonsmooth Optimization", Optimization Methods and Software, Vol. 19, No. 6, pp. 673-692, 2004. 
!*       DOI 10.1080/10556780410001689225.
!*
!*     for clustering:
!*       A. Bagirov, N. Karmitsa, S Taheri, "Partitional Clustering via Nonsmooth Optimization",
!*       Springer, 2020.
!*
!*     for NSO:
!*       A. Bagirov, N. Karmitsa, M.M. Mäkelä, "Introduction to nonsmooth optimization: theory, 
!*       practice and software", Springer, 2014.
!*
!*
!******************************************************************************
!*
!*     * PROGRAM bigclust *
!*
!*     Main program for stochastic nonsmooth clustering software with SLMB.
!*
!******************************************************************************

PROGRAM bigclust 

    USE r_precision, ONLY : prec        ! Precision for reals.
    USE param, ONLY : zero, one, large  ! Parameters.

    USE initclust, ONLY : &             ! Initialization of clustering parameters.
        infile, &                       ! Dataset file.
        outfile0, &                     ! Result file with cluster centers.
        outfile1, &                     ! Result file with function values,
                                        !   validity indices, and cpu-time.
        ibig, &                         ! Switch for the size of data, if ibig==1, no computation with whole data are made.
        ipc, &                          ! Printout parameter for cluster centers.
        maxsize, &                      ! Maximum number of candidate points of data set.
        a, &                            ! Data matrix a(nft,nrecord), from input file.
        a_used, &                       ! Partial data matrix a(nft,nrec_used), from input file.
        xbest, &
        dcent, &                        ! Distance (affinity) matrix for cluster centers
        tlimit, &                       ! Time limit, from user.
        tnorm, &                        ! The total number of norms computed.
        nclust, &                       ! Maximum number of clusters, from user.
        nft, &                          ! Number of features in data, from user.
        nrecord, &                      ! Number of instances in data, from user.
        nrec_used, &                    ! Number of data points used at the time.
        nc, &                           ! Current number of clusters, loops from 1 to nclust.
        ns, &                           ! Switch for auxiliary and real clustering problem.
        m, &                            ! Number of variables in optimization:
                                        !   m = nft    if ns = 1,
                                        !   m = nft*nc if ns = 2.
        m1, &                           ! Maximum number of initial solutions.
        ng1,ng2, &                      ! Threshold parameters, ng1,ng2 <= 100
        lcand, &
        init_clustpar, &                ! Furher initialization of parameters.
        def_clustpar                    ! Default values of clustering parameters.
    USE initslmb, ONLY : &              ! For printout purposes:
        nbatch, &                       ! Size of the batch 1 <= nbatch <= nrecords
        maxbatch, &                     ! Maximum number of different batches
        nth, &                          ! Number of iterations after which a new batch is selected
        batchtype, &                    ! Type for selecting the batches:
                                        !     0  - Totally random batches
                                        !     1  - Random batches such that all indices are used at least once      
        batchind
    USE obj_fun, ONLY : clustfull2      ! Computation of full cluster function value
    USE clusteringmod                   ! Subprograms for clustering.
    USE slmb_mod, ONLY : slmb_init      ! SLMB method for optimization.
    USE exe_time, ONLY : getime         ! Execution time.

    IMPLICIT NONE


    REAL(KIND=prec), DIMENSION(:), allocatable :: &
        x, &
        xold, &
        z, &
        x6
    REAL(KIND=prec), DIMENSION(:), allocatable :: &
        x2, &
        x5
    REAL(KIND=prec), DIMENSION(:), allocatable :: &
        amed
    REAL(KIND=prec), DIMENSION(:), allocatable :: &
        fval, &
        fval1, &
        fval2
    REAL(KIND=prec) :: &
        barf, &
        fbarmin1, &
        fbarmin2, &
        db,db2,dn,dn3,sep, & !validity indices
        f, fold, &
        f31, &
        fbarmin, &
        fbest, &
        fclust1, &
        fcurrent, &
        gamma03, &
        gamma04, &
        toler
    REAL :: &
        time1, &
        time3, &
        time4, &
        inittime, inittime2, &
        ns1time, &
        ns2time
    INTEGER :: &
        i,j,k,k1,j1, &
        ibarfmin, & 
        nstart, &
        nstart1, &
        nstart2, &
        counter, &
        nseed
    REAL :: urand    ! Random number
    !character(len=10):: junk  ! Needed e.g. with cord_19_embeddings_2022-06-02.csv
    REAL :: cputotal,usertimetotal

    INTEGER, allocatable :: seed(:)
        
    !==================================================================
    !   Initialization
    !==================================================================
    call random_seed(size = nseed)
    allocate(seed(nseed))
    
    seed = 12345
    call random_seed(put=seed)

    CALL init_clustpar()
    CALL def_clustpar()
    fold = large
    do i=1,nrec_used
        batchind(i)=i
    end do

    allocate(x(nclust*nft),xold(nclust*nft),z(nclust*nft),x6(nclust*nft),amed(nft),fval(m1),fval1(m1),fval2(m1),lcand(nrec_used))
    allocate(x2(maxsize),x5(maxsize))


    OPEN(40,file=outfile0) ! Cluster centers
    OPEN(42,file=outfile1) ! Summary of results
    OPEN(78,file=infile,status='old',form='formatted')


    !==================================================================
    !   Initializing output
    !==================================================================
    WRITE(40, *) 'BigClust: Cluster centers'
    WRITE(40, *)
    WRITE(40, *) 'This file contains the final cluster centers with k =',nclust,' clusters and'
    WRITE(40, *) 'the solutions to all intermediate l-clustering problems with l = 1,...,k-1. '
    if (ipc == 0) then
        WRITE(40, *) 
        WRITE(40, *) 'Writing down the cluster centers is canceled. Select ipc = 1 in initbigclust to enable them.'
    end if
    WRITE(40, *)
    WRITE(42, *) 'Summary of results of BigClust:'
    WRITE(42, *)
    WRITE(42, *) 'Data: ',infile
    WRITE(42, *) 'No. data points:         ',nrecord
    WRITE(42, *) 'No. used data points:    ',nrec_used
    WRITE(42, *) 'No. features:            ',nft
    WRITE(42, *) 'Size of the batch:       ',nbatch
    WRITE(42, *) 'Max no. batches:         ',maxbatch
    WRITE(42, *) 'No. iterations after which a new batch is selected: ', nth
    WRITE(42, *) 'Type of the batch:       ',batchtype
                !     0  - Totally random batches
                !     1  - Random batches such that all indices are used at least once        
    WRITE(42, *) 'Type of data:            ',ibig
                !     0  - Large data
                !     1  - Big data (no computations with full data)
    WRITE(42, *) 'Parameters:              ',m1,ng1,ng2
    WRITE(42, *) 'Seed for random numbers: ',seed(1)
    
    WRITE(42, *) 'Optimization with SLMB.'
    WRITE(42, *)
    WRITE(42, *) 'nclust |', ' f                      |',' D-B       |',' Dunn      |',&
    &' Norms       |', ' CPU-time |', ' Init. time |', ' Valid. time |',&
    &' Aux.Clust. time |', ' Clust. time |'


    !==================================================================
    !   Reading the data
    !==================================================================
    
    tnorm = zero

    !if (ibig == 0) then
        CALL getime(cputotal) 
        !READ(78,*) ! read the first line (needed e.g. with Range-Queries-Aggregates.csv)
        DO i=1,nrecord 
        !    READ(78,'(A)', ADVANCE='NO') junk ! read the first char (needed e.g. with cord_19_embeddings_2022-06-02.csv)
            READ(78,*) (a(k,i),k=1,nft)
        END DO
        CALL getime(usertimetotal) 
        print*,'Data read.',usertimetotal-cputotal
    !else !(ibig == 1)   
        ! We always read the whole data in BigClust version 0.1
    !end if 
    
    fold = large
    CALL getime(time1) 

    ! Read random partial data matrix
    if (nrecord /= nrec_used) then
        do i=1,nrec_used ! 
            call random_number(urand)
            j = i + FLOOR((nrecord+1-i)*urand)  !
            do k=1,nft !
                a_used(k,i)=a(k,j)
            end do
        end do
    else if (ibig == 0) then
        do i=1,nrec_used 
            do k=1,nft
                a_used(k,i) = a(k,i)
            end do
        end do
    else ! ibig == 1
        STOP 'Error: Big data with nrecord == nrec_used.'
    end if    

    
    !================================================================
    !   Computing cluster centers incrementally
    !================================================================
    outerloop: DO nc=1,nclust  ! number of clusters


        PRINT 42,nc
42      FORMAT('Cluster No.:',i10)

        IF(nc > 1) THEN
            if (nrecord*nft > 10000000) then
                m1  = max(m1-50,100)
            else if (nrecord*nft > 1000000) then
                m1  = max(m1-50,200)
            end if
            ng1 = max(ng1-10,10)
            !ng2 = max(ng2-10,20)
            toler=1.0E-04_prec*fclust1/REAL(nc,prec)
            CALL getime(inittime) 
            CALL step2(toler,nstart,x2) !  Step2 computes clusters for each data point and gives inital points
            CALL getime(time4) 
            inittime = time4 - inittime
            fbarmin=large
            ns1time=0.0

            !print*,'Number of aux-problems: ',nstart

            ns=1 ! Auxiliary clustering problem
                
            DO j=1,nstart ! Step 3 in SLMB-Clust algorithm
                DO k=1,nft
                    z(k)=x2(k+(j-1)*nft)
                END DO
                CALL slmb_init(z,x2((j-1)*nft+1:(j-1)*nft+nft),barf,time4) ! Call for SLMB
            
                ns1time=ns1time+time4
                fval(j)=barf
                IF (fbarmin > barf) THEN
                    fbarmin=barf
                    ibarfmin=j
                END IF                
    
            END DO

            nstart2=MAX(1,ng2*nstart/100)
            gamma04=one

            ! Add ifbarmin as center
            nstart1=1
            DO k=1,nft
                x5(k)=x2(k+(ibarfmin-1)*nft)
            END DO
            fval2(1)=fbarmin


            counter=0
            loop871: DO
                
                gamma03=gamma04
                gamma04=gamma04+1.0d-02
                fbarmin1=gamma03*fbarmin
                fbarmin2=gamma04*fbarmin
            
                DO j=1,nstart
                    IF(nstart1 >= nstart2) EXIT loop871 
                    IF (j==ibarfmin) CYCLE 
                    IF ((fval(j) >= fbarmin1).AND.(fval(j) < fbarmin2)) THEN
                        nstart1=nstart1+1
                        DO k=1,nft
                            x5(k+(nstart1-1)*nft)=x2(k+(j-1)*nft)
                        END DO
                        fval2(nstart1)=fval(j)
                    END IF

                END DO
                counter = counter + 1
                if (counter > 10000) then
                    print *,'Warning: loop871 in bigclust failed.',nstart1,nstart2,fbarmin1, fbarmin2
                    EXIT loop871
                end if
            END DO loop871


            nstart=nstart1
            DO i=1,nstart
                fval(i)=fval2(i)
            END DO
            DO i=1,nstart
                DO k=1,nft
                    x2(k+(i-1)*nft)=x5(k+(i-1)*nft)
                END DO
            END DO

            fval1(1)=fval(1)
            nstart2=1
            innerloop: DO j=2,nstart
                DO j1=1,nstart2
                    f31=zero
                    DO k=1,nft
                        f31=f31+(x5(k+(j1-1)*nft)-x2(k+(j-1)*nft))**2
                    END DO
                    IF(f31 <= (1.0E-01_prec*toler)) THEN
                        IF(fval1(j1) >= fval(j)) THEN
                            fval1(j1)=fval(j)
                            DO k=1,nft
                                x5(k+(j1-1)*nft)=x2(k+(j-1)*nft)
                            END DO
                        END IF
                        CYCLE innerloop
                    END IF
                END DO
                nstart2=nstart2+1
                DO k=1,nft
                    x5(k+(nstart2-1)*nft)=x2(k+(j-1)*nft)
                END DO
                fval1(nstart2)=fval(j)
            END DO innerloop
                
            DO i=1,nstart2
                DO k=1,nft
                    x2(k+(i-1)*nft)=x5(k+(i-1)*nft)
                END DO
            END DO 
            nstart=nstart2

            m=nft*nc
            fbest=large
            ns2time=0.0
            !print*,'Number of clustering problems: ',nstart

            DO j=1,nstart
                DO i=1,nft
                    x(i+(nc-1)*nft)=x2(i+(j-1)*nft)
                END DO
                ns=2
                CALL slmb_init(x,x6,fcurrent,time4) ! Call for SLMB
                ! Returns fcurrent in whole data sample not only in batch.

                ns2time=ns2time+time4
                IF (fcurrent < fbest) THEN
                    fbest=fcurrent
                    DO j1=1,m
                        xbest(j1)=x6(j1)
                    END DO
                END IF
            END DO
            f=fbest
            DO j1=1,m
                x(j1)=xbest(j1)
            END DO

            
            !================================================================
            !    Write cluster centers to outfile
            !================================================================
            WRITE(40, *) '____________________________________________________'
            WRITE(40, *)
            WRITE(40, *) '    Number of clusters: ',nc
            WRITE(40, *) '____________________________________________________'
            WRITE(40, *)
            if (ipc == 1) then        
                WRITE(40, *)
    
                DO j=1,nc
                    WRITE(40, *) 'Center of cluster No.',j
                    WRITE(40, *)
                    WRITE(40,49) (x(i+(j-1)*nft),i=1,nft)
                END DO
49              FORMAT(5f16.8)
                WRITE(40, *)
                
            end if

            !================================================================

        ELSE  ! nc=1, Step 0 in SLMB-Clust algorithm
            ns2time=0.0
            CALL step1(f,x) ! computes the centroid and the value of clustering function at the centroid
            fclust1=f 
            

            !================================================================
            !    Write the centroid to outfile
            !================================================================
            WRITE(40, *)
            WRITE(40, *) '____________________________________________________'
            WRITE(40, *)
            WRITE(40, *) '    Number of clusters: ',nc
            WRITE(40, *) '____________________________________________________'
            WRITE(40, *)

            if (ipc == 1) then
                WRITE(40, *) 'Center of cluster No.',nc
                WRITE(40, *)
                WRITE(40,48) (x(i),i=1,nft)
48              FORMAT(5f16.8)
                WRITE(40, *)
            end if
            !================================================================

        END IF

        ! Read random partial data matrix
        if (nrecord /= nrec_used) then
            do i=1,nrec_used ! 
                call random_number(urand)
                j = i + FLOOR((nrecord+1-i)*urand)  !
                do k=1,nft !
                    a_used(k,i)=a(k,j)
                end do
            end do
        end if
    
        DO k=1,nc
            dcent(k,k)=zero
        END DO
        DO k=1,nc
            DO k1=k+1,nc
                dcent(k,k1)=zero
                DO j=1,nft
                    dcent(k,k1)=dcent(k,k1)+(xbest(j+(k-1)*nft)-xbest(j+(k1-1)*nft))**2
                END DO
                dcent(k1,k)=dcent(k,k1)
            END DO
        END DO

        CALL getime(inittime2) 
        CALL check(x,f,db,db2,dn,dn3,sep) ! checks the results and computes validity indices
        CALL getime(time4) 
        
        inittime2 = time4 - inittime2
        if (ibig == 0) then ! We do not compute this with big data
            if (nrecord /= nrec_used) then ! or if entire data set is used
                call clustfull2(x,f)
            end if
            if (fold < f) then
                do i=1,nft
                    xold((nc-1)*nft+i)=a_used(i,1) ! We may have better options, but this is quick.
                end do
                call clustfull2(xold,f)
                fold=f
                x=xold
            else
                fold=f
                xold=x
            end if
        end if


        CALL getime(time3)

        WRITE(40,*)
        WRITE(40,543) f
        WRITE(40,*)
        WRITE(40,142) INT(tnorm)
        time4=time3-time1
        if (time4 > tlimit) then
            print*,' Warning: time limit exceed. Returned final solution has only ',nc,' clusters.'
            EXIT outerloop
        end if
        WRITE(40,*)
        WRITE(40,141) time4
        WRITE(40,*)
        WRITE(40,*)
        WRITE(42,603) nc,f,db,dn,tnorm,time4,inittime,inittime2,ns1time,ns2time
        !print*,'cpu time ',time4
    END DO outerloop
    WRITE(42, *)
    
    WRITE(42, *) 'Data read in ',usertimetotal-cputotal,'seconds.'
    WRITE(42, *)
    WRITE(42, *)

    
    
141 FORMAT('CPU time:                                         ',f12.3)
142 FORMAT('The total number of norms:                  ',i18)
543 FORMAT('The value of the cluster function:',f28.6)
603 FORMAT(i4,f38.4,f12.4,f12.4,f15.0,f10.2,f10.2,f10.2,f10.2,f10.2)
    
    CLOSE(40)
    CLOSE(42)
    CLOSE(78)
        
    deallocate(x,xold,z,x6,amed,fval,fval1,fval2,lcand)
    deallocate(x2,x5)
    deallocate(seed)

    STOP

END PROGRAM bigclust
