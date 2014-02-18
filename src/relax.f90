module relax

use utils, only: get_free_unit
use errors, only: stop_all,warning
use basis_data
use basis, only: FindSqPairInd
use input_data, only: tDyall,tZerothRelaxBoth
use int_management, only: TransposeIntegralArr

implicit none

contains

    !This routine will calculate the orbital relaxation energy, tending to the 
    !Hartree--Fock CBS limit. This calculates the relaxation from the Hartree--Fock
    !determinant only, so is just the single-reference (SR) correction.
    subroutine Calc_SR_CABS_Singles(Energy_CABS_S)
        implicit none
        real(dp), intent(out) :: Energy_CABS_S
        integer :: ierr,i,j,a,lWork,info,OccOrb,OccInd,t
        real(dp) , allocatable :: TempFockMat(:,:),TempFockMat_2(:,:),Work(:),W(:)
        real(dp) , allocatable :: TempMat(:,:),TempMat_2(:,:)
!        real(dp) , allocatable :: TempFockMat_test(:,:)
        integer , allocatable :: OccOrbRotInd(:) 
        logical :: tAisOcc
        character(*), parameter :: t_r='Calc_SR_CABS_Singles'
        logical, parameter :: tDebug=.false.
            
        write(6,*) ""
        write(6,"(A)") "Calculating single-reference orbital relaxation energy..."

        Energy_CABS_S = 0.0_dp

        !We first want to fill a matrix with the fock matrix over the union of the virtual and CABS spaces.
        allocate(TempFockMat(tntmo,tntmo),stat=ierr)
!        allocate(TempFockMat_test(tntmo,tntmo),stat=ierr)
        allocate(TempFockMat_2(tntmo,tntmo),stat=ierr)
        TempFockMat(:,:) = 0.0_dp
        TempFockMat_2(:,:) = 0.0_dp
!        TempFockMat_test(:,:) = 0.0_dp

        do i=1,tntmo
            do j=1,tntmo
                if(min(i,j).gt.tnmo) then
                    TempFockMat(i,j) = FockCABS(i,j)
!                    TempFockMat_test(i,j) = FockCABS(i,j)
                    TempFockMat_2(i,j) = FockCABS(i,j)
                elseif(max(i,j).gt.tnmo) then
                    TempFockMat(i,j) = FockOrbCABS(min(i,j),max(i,j))
!                    TempFockMat_test(i,j) = FockOrbCABS(min(i,j),max(i,j))
                    TempFockMat_2(i,j) = FockOrbCABS(min(i,j),max(i,j))
                else
                    TempFockMat(i,j) = FockOrb(i,j)
!                    TempFockMat_test(i,j) = FockOrb(i,j)
                    TempFockMat_2(i,j) = FockOrb(i,j)
                endif
            enddo
        enddo

        j=0
        do i=1,tnmo
            if(IsOrbOcc(i)) then
                j=j+1
                !Ensure that the TempFockMat matrix is unit for the occupied orbitals
                TempFockMat(i,:) = 0.0_dp
                TempFockMat(:,i) = 0.0_dp
                TempFockMat(i,i) = 5000.0_dp+j !Break the degeneracy of the occupied orbital, so they don't mix.
!                TempFockMat_test(i,:) = 0.0_dp
!                TempFockMat_test(:,i) = 0.0_dp
!                TempFockMat_test(i,i) = 1.0_dp
            endif
        enddo
        
        if(tDebug) then
            write(6,*) "Original Virt:CABS fock Matrix: "
            do i=1,tntmo
                do j=1,tntmo
                    write(6,'(F11.5)',advance='no') TempFockMat(i,j)
                enddo
                write(6,*)
            enddo
        endif

        !OccOrbRotInd gives the index of the unchanged original orbitals in the new basis
        allocate(OccOrbRotInd(tnocc))
        OccOrbRotInd = 0

        !Diagonalise this matrix - ie, rotate the virt + CABS block so that it is diagonal
        allocate(W(tntmo),stat=ierr)
        lWork = -1
        allocate(Work(1))
        call DSYEV('V','U',tntmo,TempFockMat,tntmo,W,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,"Workspace quiery failed")
        lWork = Work(1)
        deallocate(Work)
        allocate(Work(lWork),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Allocation Error")
        call DSYEV('V','U',tntmo,TempFockMat,tntmo,W,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,"Fock matrix diagonalisation failed")
        
        if(tDebug) then
            write(6,*) "U Matrix: "
            do i=1,tntmo
                do j=1,tntmo
                    write(6,'(F10.5)',advance='no') TempFockMat(i,j)
                enddo
                write(6,*)
            enddo
        endif

        do i=1,tnocc
            OccOrb = OccOrbs(i) !Get the index in the original basis
            do j=1,tntmo
                if(abs(TempFockMat(OccOrb,j)).gt.1.D-8) then
                    if(abs(TempFockMat(OccOrb,j)-1.0_dp).lt.1.D-8) then
                        !This is the index in the new basis!
                        OccOrbRotInd(i) = j
                    elseif(abs(TempFockMat(OccOrb,j)+1.0_dp).lt.1.D-8) then
                        !The eigenfunction could also have values = -1, which would also be valid
                        !This is the index in the new basis!
                        OccOrbRotInd(i) = j
                    else
                        write(6,*) "TempFockMat(OccOrb,j): ",TempFockMat(OccOrb,j)
                        write(6,*) "Occupied orbital row (coeffs in new basis): ",i,OccOrb
                        do a=1,tntmo
                            write(6,*) a, TempFockMat(OccOrb,a)
                        enddo
                        call stop_all(t_r,"Occupied orbitals have been rotated from original basis!")
                    endif
                endif
            enddo
        enddo


        allocate(TempMat(tntmo,tntmo),stat=ierr)
        allocate(TempMat_2(tntmo,tntmo),stat=ierr)

!        open(27,file='EigenvectorMat')
!        do i=1,tntmo
!            do j=1,tntmo
!                write(27,"(F14.7)",advance='no') TempFockMat(i,j)
!            enddo
!            write(27,*) ""
!        enddo
!        close(27)
!
!
!        !TEST - check that if the original matrix is rotated, that you get a diagonal matrix...
!        call DGEMM('T','N',tntmo,tntmo,tntmo,1.0_dp,TempFockMat,tntmo,TempFockMat_test,tntmo,0.0_dp,TempMat,tntmo) 
!        call DGEMM('N','N',tntmo,tntmo,tntmo,1.0_dp,TempMat,tntmo,TempFockMat,tntmo,0.0_dp,TempMat_2,tntmo)
!        
!        do i=1,tntmo
!            do j=1,tntmo
!                write(27,"(F14.7)",advance='no') TempMat_2(i,j)
!            enddo
!            write(27,*) ""
!        enddo
!
!        !END TEST


        !Now rotate the full original fock matrix into this new basis
        call DGEMM('T','N',tntmo,tntmo,tntmo,1.0_dp,TempFockMat,tntmo,TempFockMat_2,tntmo,0.0_dp,TempMat,tntmo) 
        call DGEMM('N','N',tntmo,tntmo,tntmo,1.0_dp,TempMat,tntmo,TempFockMat,tntmo,0.0_dp,TempMat_2,tntmo)
        
        !This writes out the rotated basis fock matrix.
        !This should be diagonal in the virtual + CABS blocks
!        do i=1,tntmo
!            do j=1,tntmo
!                write(27,"(F14.7)",advance='no') TempMat_2(i,j)
!            enddo
!            write(27,*) ""
!        enddo


        !TempMat_2 now contains the fock matrix in the new basis
        !Check that the occupied orbitals are still present in their original form (though they may have been moved in position?)
        do i=1,tnocc
            if(abs(TempMat_2(OccOrbRotInd(i),OccOrbRotInd(i))-FockOrb(OccOrbs(i),OccOrbs(i))).gt.1.D-8) then
                write(6,'(A,I5)') 'Original orbital index: ',OccOrbs(i)
                write(6,'(A,I5)') 'Rotated basis orbital index: ',OccOrbRotInd(i)
                write(6,'(A,F17.10)') 'Original fock eigenvalue: ',FockOrb(OccOrbs(i),OccOrbs(i))
                write(6,'(A,F17.10)') 'Rotated fock eigenvalue: ',TempMat_2(OccOrbRotInd(i),OccOrbRotInd(i))
                call stop_all(t_r,"Error in rotation of fock matrix")
            endif
        enddo

        !Now find the energy... 
        do i=1,tnocc
            !i loops over occupied orbitals
            OccInd = OccOrbRotInd(i)

            do a=1,tntmo
                !a loops over rest of space

                !check whether a is an original occupied orbital in the rotated space
                tAisOcc = .false.
                do t=1,tnocc
                    if(OccOrbRotInd(t).eq.a) tAisOcc = .true.
                enddo
                if(tAisOcc) cycle

                Energy_CABS_S = Energy_CABS_S + (TempMat_2(OccInd,a)*TempMat_2(OccInd,a))/(TempMat_2(OccInd,OccInd)-TempMat_2(a,a))

            enddo
        enddo

        Energy_CABS_S = Energy_CABS_S * 2.0_dp  !To account for aa and bb

        write(6,'(A,F17.10)') "Single-reference orbital relaxation energy: ",Energy_CABS_S
        write(6,*) ""

        deallocate(TempFockMat,TempFockMat_2,Work,W,TempMat_2,TempMat,OccOrbRotInd)

    end subroutine Calc_SR_CABS_Singles

    !This routine will provide the [2]_s multireference orbital relaxation correction of Valeev.
    !Part of the X matrix has already been calculated in the reading in of the 2RDM code.
    subroutine Calc_MR_CABS_Singles(Energy_CABS_S)
        implicit none
        real(dp), intent(out) :: Energy_CABS_S
        real(dp), allocatable :: TempArr(:),TempArr_2(:),G_arr_aaaa(:),G_arr_abba(:)
        real(dp), allocatable :: G_arr_abab(:),W(:),Work(:)
        real(dp), allocatable :: CABSGenFockEigens(:,:),RelaxX(:,:)
        real(dp) :: Temp,DDOT,RelaxEnergy
        integer, allocatable :: OccOrbRotInd(:)
        integer :: g_ints_read,y,k,yk_pair,loop,j,jk_pair,x,i,ierr,xk_pair
        integer ::lWork,info
        character(*), parameter :: t_r='Calc_MR_CABS_Singles'
        logical :: tFail
        logical , parameter :: tDebug=.false.

        Energy_CABS_S = 0.0_dp

        write(6,*) ""
        if(tZerothRelaxBoth) then
            write(6,"(A)") "Calculating MR orbital relaxation, using both the Dyall and generalised Fock matrix zeroth-order Hamiltonians"
        elseif(tDyall) then
            write(6,"(A)") "Calculating MR orbital relaxation, using the Dyall zeroth-order Hamiltonian"
        else
            write(6,"(A)") "Calculating MR orbital relaxation, using the generalised Fock matrix zeroth-order Hamiltonian"
        endif
        call flush(6)

        !The first thing to do, is to finish of calculating the X_y^x matrix for the two types of
        !zeroth-order hamiltonian.

        !First, the rest of the cumulant expansion...
        if(tDyall.or.tZerothRelaxBoth) then

            !Need the g integrals to contract the RDMs with.
            g_ints_read = get_free_unit()
            open(g_ints_read,file='G_Dir_ccoo',status='old',form='unformatted',   &
                        access='direct',action='read',recl=reclen_nmo_sq)

            !TempArr_y^j = G_yk^ij 1RDM_i^k
            allocate(TempArr(tnmo_sq),stat=ierr)
            allocate(TempArr_2(tnmo_sq),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc Fail")
            TempArr(:) = 0.0_dp
            TempArr_2(:) = 0.0_dp

            allocate(G_arr_aaaa(tnmo_sq),stat=ierr) !Allocate for all antisymmetric (g_yk)^ij
            allocate(G_arr_abab(tnmo_sq),stat=ierr) !Allocate for all antisymmetric (g_yk)^ij
            allocate(G_arr_abba(tnmo_sq),stat=ierr) !Allocate for all antisymmetric (g_yk)^ij
            if(ierr.ne.0) call stop_all(t_r,"Alloc Fail")

            do y=1,tnmo
                do k=1,tnmo
                    yk_pair = FindSqPairInd(y,k,tntmo)
!                    read(g_ints_read,rec=yk_pair) (G_arr_aaaa(loop),loop=1,tnmo_sq)
                    read(g_ints_read,rec=yk_pair) (G_arr_abab(loop),loop=1,tnmo_sq)

                    !antisym first
                    call TransposeIntegralArr(G_arr_abab,G_arr_abba(:),tnmo)
                    G_arr_abba(:) = -G_arr_abba(:)
                    G_arr_aaaa(:) = G_arr_abab(:) + G_arr_abba(:)


                    !This section deals with the contibution to X_y^x from g_yk^ij 1RDM_i^x 1RDM_j^k
                    !Contract over i to get X_yk^jx (j fast)
                    call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,G_arr_aaaa,tnmo,OneRDM,tnmo,0.0_dp,TempArr_2,tnmo)
                    !Contract over j to get X_yk^xk 
                    call DGEMM('T','N',tnmo,tnmo,tnmo,1.0_dp,TempArr_2,tnmo,OneRDM,tnmo,0.0_dp,TempArr,tnmo)

                    do x=1,tnmo
                        RelaxXDyall(y,x) = RelaxXDyall(y,x) - TempArr(FindSqPairInd(k,x,tnmo))!*2.0_dp 
                    enddo
                    
                    !Contract over i to get X_yk^jx
                    call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,G_arr_abab,tnmo,OneRDM,tnmo,0.0_dp,TempArr_2,tnmo)
                    !Contract over j to get X_yk^xk
                    call DGEMM('T','N',tnmo,tnmo,tnmo,1.0_dp,TempArr_2,tnmo,OneRDM,tnmo,0.0_dp,TempArr,tnmo)

                    do x=1,tnmo
                        RelaxXDyall(y,x) = RelaxXDyall(y,x) - TempArr(FindSqPairInd(k,x,tnmo)) 
                    enddo


                    !This section calculates the contribution to X_y^x from g_yk^ij OneRDM_i^k OneRDM_j^x
                    !contract over i, to give TempArr_2 - jk with j fast.
                    call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,G_arr_aaaa(:),tnmo,OneRDM,tnmo,0.0_dp,TempArr_2,tnmo)
                    G_arr_aaaa(:) = TempArr_2(:)  !This is now (X_yk)^jk
                    call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,G_arr_abba(:),tnmo,OneRDM,tnmo,0.0_dp,TempArr_2,tnmo)
                    G_arr_abba(:) = TempArr_2(:)
                    
                    !Now contract over j, to give TempArr_2 - kx with k fast
                    call DGEMM('T','N',tnmo,tnmo,tnmo,1.0_dp,G_arr_aaaa,tnmo,OneRDM,tnmo,0.0_dp,TempArr_2,tnmo)
                    G_arr_aaaa(:) = TempArr_2(:)
                    call DGEMM('T','N',tnmo,tnmo,tnmo,1.0_dp,G_arr_abba,tnmo,OneRDM,tnmo,0.0_dp,TempArr_2,tnmo)
                    G_arr_abba(:) = TempArr_2(:)

                    !The G_arr_xxxx arrays now contain (X_yk)^kx - we need to sum over all k

                    do x=1,tnmo
                        xk_pair = FindSqPairInd(x,k,tnmo)
                        RelaxXDyall(y,x) = RelaxXDyall(y,x) + G_arr_aaaa(xk_pair) + G_arr_abba(xk_pair)
                    enddo

                enddo
            enddo
            
!            !This section deals with the contibution to X_y^x from g_yk^ij 1RDM_i^x 1RDM_j^k
!            do y=1,tnmo
!                do k=1,tnmo
!                    yk_pair = FindSqPairInd(y,k,tntmo)
!                    read(g_ints_read,rec=yk_pair) (G_arr_abab(loop),loop=1,tnmo_sq)
!
!                    !antisym the aaaa and abab terms.
!                    call TransposeIntegralArr(G_arr_abab(:),G_arr_aaaa(:),tnmo)
!                    G_arr_aaaa(:) = G_arr_abab(:) - G_arr_aaaa(:)
!                    
!
!                    !Contract over i to get X_yk^jx (j fast)
!                    call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,G_arr_aaaa,tnmo,OneRDM,tnmo,0.0_dp,TempArr_2,tnmo)
!                    !Contract over j to get X_yk^xk 
!                    call DGEMM('T','N',tnmo,tnmo,tnmo,1.0_dp,TempArr_2,tnmo,OneRDM,tnmo,0.0_dp,TempArr,tnmo)
!
!                    do x=1,tnmo
!                        RelaxXDyall(y,x) = RelaxXDyall(y,x) - TempArr(FindSqPairInd(k,x,tnmo)) 
!                    enddo
!                    
!                    !Contract over i to get X_yk^jx
!                    call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,G_arr_abab,tnmo,OneRDM,tnmo,0.0_dp,TempArr_2,tnmo)
!                    !Contract over j to get X_yk^xk
!                    call DGEMM('T','N',tnmo,tnmo,tnmo,1.0_dp,TempArr_2,tnmo,OneRDM,tnmo,0.0_dp,TempArr,tnmo)
!
!                    do x=1,tnmo
!                        RelaxXDyall(y,x) = RelaxXDyall(y,x) - TempArr(FindSqPairInd(k,x,tnmo)) 
!                    enddo
!                enddo
!            enddo

            !Now for the first part of X - contract i in f_y^i OneRDM_i^x
            call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,GenFockOrb,tnmo,OneRDM,tnmo,0.0_dp,TempArr,tnmo)

            !Finally, put everything together to get the correct X_y^x for the dyall hamiltonian.
            do y=1,tnmo
                do x=1,tnmo
                    RelaxXDyall(y,x) = (-RelaxXDyall(y,x)/2.0_dp) - TempArr(FindSqPairInd(x,y,tnmo))
                enddo
            enddo

            close(g_ints_read)  !Also need to deallocate stuff
            deallocate(TempArr,TempArr_2,G_arr_aaaa,G_arr_abab,G_arr_abba)
            
!            !TEST - write out X
!            open(27,file='Relax-MRF12')
!            do i=1,tnmo
!                do j=1,tnmo
!                    write(27,*) i,j,RelaxXDyall(i,j)
!                enddo
!            enddo
!            close(27)
            
        endif
        

        if((.not.tDyall).or.tZerothRelaxBoth) then
            !Calculate X-matrix for the non-dyall hamiltonian.
            Temp = DDOT(tnmo_sq,GenFockOrb(1,1),1,OneRDM(1,1),1)
            Temp = Temp * 2.0_dp    !To account for the bb block

            do y=1,tnmo
                do x=1,tnmo
                    RelaxXFock(y,x) = RelaxXFock(y,x) - Temp*OneRDM(y,x)
                enddo
            enddo

        endif
        
        if(tDebug) then
            !TEST - write out X
            open(27,file='Relax-MRF12')
            do i=1,tnmo
                do j=1,tnmo
                    write(27,*) i,j,RelaxXDyall(i,j)
                enddo
            enddo
            close(27)
        endif

        !Now diagonalise the generalised fock matrix in the CABS CABS space.
        !The identity matrix is added for the orbital space.
        !GenFockCABS runs from tnmo+1:tntmo
        allocate(CABSGenFockEigens(tntmo,tntmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Alloc err')
        CABSGenFockEigens(:,:) = 0.0_dp
        CABSGenFockEigens(tnmo+1:tntmo,tnmo+1:tntmo) = GenFockCABS(tnmo+1:tntmo,tnmo+1:tntmo)
        j=0
        do i=1,tnmo
            CABSGenFockEigens(i,i) = 5000.0_dp+j    !This is to break the 'degeneracy' between the occupied orbtials, 
                                                            !to avoid mixing.
        enddo
        
        !OccOrbRotInd gives the index of the unchanged original OBS orbitals in the new rotated basis
        allocate(OccOrbRotInd(tnmo))
        OccOrbRotInd = 0
        
        !Diagonalise this matrix - ie, rotate the CABS block so that it is diagonal
        allocate(W(tntmo),stat=ierr)
        lWork = -1
        allocate(Work(1))
        call DSYEV('V','U',tntmo,CABSGenFockEigens,tntmo,W,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,"Workspace quiery failed")
        lWork = Work(1)
        deallocate(Work)
        allocate(Work(lWork),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Allocation Error")
        call DSYEV('V','U',tntmo,CABSGenFockEigens,tntmo,W,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,"Fock matrix diagonalisation failed")
        deallocate(Work)

        if(tDebug) then
            write(6,*) "U Matrix: "
            do i=1,tntmo
                do j=1,tntmo
                    write(6,'(F10.5)',advance='no') CABSGenFockEigens(i,j)
                enddo
                write(6,*)
            enddo
        endif

        do i=1,tnmo
            do j=1,tntmo
                if(abs(CABSGenFockEigens(i,j)).gt.1.D-8) then
                    if(abs(CABSGenFockEigens(i,j)-1.0_dp).lt.1.D-8) then
                        !This is the index in the new basis!
                        OccOrbRotInd(i) = j
                    elseif(abs(CABSGenFockEigens(i,j)+1.0_dp).lt.1.D-8) then
                        !The eigenfunction could also have values = -1, which would also be valid
                        !This is the index in the new basis!
                        OccOrbRotInd(i) = j
                    else
                        call stop_all(t_r,"OBS orbitals have been rotated from original basis!")
                    endif
                endif
            enddo
        enddo

        if(tDebug) then
            write(6,*) "Original OBS basis      New Basis"
            do i=1,tnmo
                write(6,*) i,OccOrbRotInd(i)
            enddo
        endif

        allocate(RelaxX(tnmo,tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc Err")

        if((.not.tDyall).or.tZerothRelaxBoth) then
            RelaxX(:,:) = RelaxXFock(:,:)
            call FindMRRelaxEnergy(RelaxEnergy,W,RelaxX,CABSGenFockEigens,OccOrbRotInd,tFail)

            if(.not.tFail) then
                write(6,'(A,F17.10)') "Orbital relaxation energy with generalised Fock 0th-order hamiltonian: ",RelaxEnergy
                Energy_CABS_S = RelaxEnergy 
            else
                write(6,'(A)') "Relaxation equations could not be solved for the generalised Fock 0th-order hamiltonian."
                Energy_CABS_S = sqrt(-1.0_dp)
            endif

        endif
        
        !Take the required X matrix, find the coefficients and then the energy.
        if(tDyall.or.tZerothRelaxBoth) then
            RelaxX(:,:) = RelaxXDyall(:,:)
            call FindMRRelaxEnergy(RelaxEnergy,W,RelaxX,CABSGenFockEigens,OccOrbRotInd,tFail)
            
            if(.not.tFail) then
                write(6,'(A,F17.10)') "Orbital relaxation energy with Dyall 0th-order hamiltonian: ",RelaxEnergy
                if(.not.tZerothRelaxBoth) Energy_CABS_S = RelaxEnergy !If tZerothRelaxBoth is on, then it will use the Fock hamiltonian in printed stats. 
            else
                write(6,'(A)') "Relaxation equations could not be solved for the Dyall 0th-order hamiltonian."
                if(.not.tZerothRelaxBoth) Energy_CABS_S = sqrt(-1.0_dp)
            endif

        endif

        write(6,*) ""
        call flush(6)

        deallocate(W,RelaxX,CABSGenFockEigens,OccOrbRotInd)
        if(tDyall.or.tZerothRelaxBoth) deallocate(RelaxXDyall)
        if((.not.tDyall).or.tZerothRelaxBoth) deallocate(RelaxXFock)

    end subroutine Calc_MR_CABS_Singles

    subroutine FindMRRelaxEnergy(RelaxEnergy,W,RelaxX,CABSGenFockEigens,OccOrbRotInd,tDGESVFail)
        use matrixops, only: mat_inv 
        implicit none
        real(dp) , intent(out) :: RelaxEnergy
        real(dp) , intent(in) :: W(tntmo),RelaxX(tnmo,tnmo),CABSGenFockEigens(tntmo,tntmo)
        logical , intent(out) :: tDGESVFail
        integer , intent(in) :: OccOrbRotInd(tnmo)
        logical , allocatable :: IsOrigOrb(:),tRemoveCol(:)
        real(dp) , allocatable :: Mat(:,:),f_tld(:,:),f_tld_2(:,:),f_tld_row(:),Coeffs(:,:),f_tld_row_2(:)
        real(dp) , allocatable :: RotGenFock(:,:),RotOneRDM(:,:),TempMat(:,:),RotX(:,:),Work(:),NonSingCoeffs(:)
        real(dp) , allocatable :: ReturnMat(:),MatInv(:,:),TransformMat(:,:),W_Real(:),W_Imag(:),DiagWork(:)
        real(dp) , allocatable :: NonSingMat(:,:),NonSing_f_tld_row(:),NonSingTransformMat(:,:),NonSingEigenval_Real(:)
        real(dp) , allocatable :: Eigenvec_Left(:,:),Eigenvec_Right(:,:),NonSingTransformMat_Inv(:,:),NonSingEigenval_Imag(:)
        real(sp) , allocatable :: SWork(:)
        real(dp) :: Temp,MaxColEl,MaxRowEl
        integer :: i,j,ierr,info,x,b,betap,y,aprime,a,iter,tnmo_red,ColEl,RowEl,alpha!,k,l
        integer , allocatable :: iPiv(:)
        character(*), parameter :: t_r='FindMRRelaxEnergy'
        logical , parameter :: tDebug=.false.
        logical , parameter :: tIterative=.false.
        logical , parameter :: tFindSpectrum = .false. 
        
        if(tDebug) then
            write(6,*) "X Matrix: "
            do i=1,tnmo
                do j=1,tnmo
                    write(6,'(F10.5)',advance='no') RelaxX(i,j)
                enddo
                write(6,*)
            enddo
        endif

        tDGESVFail = .false.

        !Make a logical array to quickly determine whether an orbital in the new basis is
        !one of the original OBS orbitals, or a rotated CABS orbital
        allocate(IsOrigOrb(tntmo),stat=ierr)
        IsOrigOrb(:)=.false.
        do i=1,tntmo
            do j=1,tnmo
                if(OccOrbRotInd(j).eq.i) then
                    IsOrigOrb(i)=.true.
                    exit
                endif
            enddo
        enddo

        if(tDebug) then
            write(6,*) "Original OBS orbital?"
            do i=1,tntmo
                write(6,*) i,IsOrigOrb(i)
            enddo
        endif

        allocate(f_tld(tnmo,tnxmo),stat=ierr)
        allocate(f_tld_2(tntmo,tntmo),stat=ierr)
        allocate(TempMat(tntmo,tntmo),stat=ierr)
        allocate(Mat(tnmo,tnmo),stat=ierr)
        allocate(MatInv(tnmo,tnmo),stat=ierr)
        allocate(iPiv(tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc Err")

        !Create the f_tld matrix
        call DGEMM('N','N',tnmo,tnxmo,tnmo,-1.0_dp,OneRDM,tnmo,GenFockOrbCABS,tnmo,0.0_dp,f_tld,tnmo)

        !Make this matrix over the whole space to make things easier later.
        f_tld_2(:,:) = 0.0_dp
        do x=1,tnmo
            do b=tnmo+1,tntmo
                f_tld_2(x,b) = f_tld(x,b-tnmo)
                f_tld_2(b,x) = f_tld_2(x,b)
            enddo
        enddo

        deallocate(f_tld)
        allocate(f_tld(tntmo,tntmo))

        if(tDebug) then
            write(6,*) "F_tld - old basis: "
            do i=1,tntmo
                do j=1,tntmo
                    write(6,'(F10.5)',advance='no') f_tld_2(i,j)
                enddo
                write(6,*)
            enddo
        endif

        !We now need to rotate this into the new basis
        call DGEMM('T','N',tntmo,tntmo,tntmo,1.0_dp,CABSGenFockEigens,tntmo,f_tld_2,tntmo,0.0_dp,TempMat,tntmo)
        call DGEMM('N','N',tntmo,tntmo,tntmo,1.0_dp,TempMat,tntmo,CABSGenFockEigens,tntmo,0.0_dp,f_tld,tntmo)
        
        if(tDebug) then
            write(6,*) "F_tld - new basis: "
            do i=1,tntmo
                do j=1,tntmo
                    write(6,'(F10.5)',advance='no') f_tld(i,j)
                enddo
                write(6,*)
            enddo
        endif

        !f_tld now contains the f_tld matrix in the new basis.
        !Also rotate X, fock and OneRDM into the new basis...
        allocate(RotGenFock(tntmo,tntmo))
        allocate(RotOneRDM(tntmo,tntmo))
        allocate(RotX(tntmo,tntmo))

        RotGenFock(:,:) = 0.0_dp
        do i=1,tntmo
            do j=1,tntmo
                if(min(i,j).gt.tnmo) then
                    RotGenFock(i,j) = GenFockCABS(i,j)
                elseif(max(i,j).gt.tnmo) then
                    RotGenFock(i,j) = GenFockOrbCABS(min(i,j),max(i,j))
                elseif(max(i,j).le.tnmo) then
                    RotGenFock(i,j) = GenFockOrb(i,j)
                endif
            enddo
        enddo

        call DGEMM('T','N',tntmo,tntmo,tntmo,1.0_dp,CABSGenFockEigens,tntmo,RotGenFock,tntmo,0.0_dp,TempMat,tntmo)
        call DGEMM('N','N',tntmo,tntmo,tntmo,1.0_dp,TempMat,tntmo,CABSGenFockEigens,tntmo,0.0_dp,RotGenFock,tntmo)

        RotOneRDM(:,:) = 0.0_dp
        RotOneRDM(1:tnmo,1:tnmo) = OneRDM(1:tnmo,1:tnmo)
        call DGEMM('T','N',tntmo,tntmo,tntmo,1.0_dp,CABSGenFockEigens,tntmo,RotOneRDM,tntmo,0.0_dp,TempMat,tntmo)
        call DGEMM('N','N',tntmo,tntmo,tntmo,1.0_dp,TempMat,tntmo,CABSGenFockEigens,tntmo,0.0_dp,RotOneRDM,tntmo)

        RotX(:,:) = 0.0_dp
        RotX(1:tnmo,1:tnmo) = RelaxX(1:tnmo,1:tnmo)
        call DGEMM('T','N',tntmo,tntmo,tntmo,1.0_dp,CABSGenFockEigens,tntmo,RotX,tntmo,0.0_dp,TempMat,tntmo)
        call DGEMM('N','N',tntmo,tntmo,tntmo,1.0_dp,TempMat,tntmo,CABSGenFockEigens,tntmo,0.0_dp,RotX,tntmo)
        
        if(tDebug) then
            write(6,*) "Rotated X Matrix: "
            do i=1,tntmo
                do j=1,tntmo
                    write(6,'(F10.5)',advance='no') RotX(i,j)
                enddo
                write(6,*)
            enddo
        endif

        allocate(f_tld_row(tnmo))
        allocate(f_tld_row_2(tnmo))
        allocate(Coeffs(tntmo,tntmo))
        if(tIterative) then
            allocate(ReturnMat(tnmo),stat=ierr)
            ReturnMat(:) = 0.0_dp
            allocate(Work(tnmo),stat=ierr)
            allocate(SWork(tnmo*(tnmo+1)),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc failed")
        endif
        Coeffs(:,:) = 0.0_dp

        do betap=1,tntmo

            !This only want to loop over cabs orbitals in the rotated basis
            if(IsOrigOrb(betap)) cycle
            
            if(tDebug) then
                write(6,*) "betap = ",betap
            endif

            !Form the relevant matrices
            !Take the transpose, so it can be written as M_xy C_y = f_x
            !i.e. we want x fast
            do y=1,tnmo
                do x=1,tnmo
                    Mat(x,y) = W(betap)*OneRDM(x,y) + RelaxX(y,x)
!                    Mat(x,y) = W(betap)*OneRDM(x,y) + RelaxX(x,y)
                enddo
            enddo
            !Also get the relevant row of f_tld. We need to make sure that we get the correct orbital, since
            !their index will have changed.
            do i=1,tntmo

                if(.not.IsOrigOrb(i)) cycle !We want to extract the OBS orbitals

                do j=1,tnmo
                    if(OccOrbRotInd(j).eq.i) exit
                enddo
                if(j.gt.tnmo) call stop_all(t_r,"Couldn't find orbital")
                !j now indicates the index in the original basis
                if(tDebug) then
                    write(6,'(A,I4,A,I4)') "Filling f_tld with new orbital: ",j, " from rotated orbital: ",i
                endif

                f_tld_row(j) = f_tld(betap,i)
            enddo

!            if(tDebug) then
!                write(6,*) "F_tld vector: "
!                do i=1,tnmo
!                    write(6,*) f_tld_row(i)
!                enddo
!                write(6,*) "Matrix: "
!                do i=1,tnmo
!                    do j=1,tnmo
!                        write(6,'(F10.5)',advance='no') Mat(i,j)
!                    enddo
!                    write(6,*)
!                enddo
!            endif

            if(tFindSpectrum) then


                !Diagonalise the matrix, in order to remove singular components
                !Beware that the matrix is NOT symmetric, and so we need to calculate both left & right eigenvectors.
                allocate(TransformMat(tnmo,tnmo))
                TransformMat(:,:) = Mat(:,:)
                allocate(tRemoveCol(tnmo))
                tRemoveCol(:) = .false.
                allocate(DiagWork(max(1,8*tnmo)))
                allocate(W_Real(tnmo))
                allocate(W_Imag(tnmo))
                allocate(Eigenvec_Left(tnmo,tnmo))
                allocate(Eigenvec_Right(tnmo,tnmo))
                call DGEEV('V','V',tnmo,TransformMat,tnmo,W_Real,W_Imag,Eigenvec_Left,tnmo,Eigenvec_Right,tnmo,DiagWork,max(1,8*tnmo),info)
                if(info.ne.0) call stop_all(t_r,"Diag error")

!                write(6,*) "Eigenvalues:"
                do i=1,tnmo
!                    write(6,*) W_real(i),W_Imag(i)
                    if(abs(W_Imag(i)).gt.1.D-6) then
                        write(6,*) W_real(i),W_Imag(i)
                        write(6,*) "Warning: Complex eigenvalue"
                    endif
                enddo

                !Search for zero eigenvalues
                tnmo_red = tnmo
                tRemoveCol(:) = .false.
                do i=1,tnmo
                    if(SQRT((W_Real(i)**2)+(W_Imag(i)**2)).lt.1.D-8) then
                        !Remove eigen-component
                        write(6,*) "Removing eigenvector: ",W_Real(i),W_Imag(i)
                        tRemoveCol(i) = .true.
                        tnmo_red = tnmo_red - 1
                    endif
                enddo

                !Find new transformation matrix
                allocate(NonSingTransformMat_Inv(tnmo_red,tnmo))
                allocate(NonSingTransformMat(tnmo,tnmo_red))
                allocate(NonSingEigenval_Real(tnmo_red))
                allocate(NonSingEigenval_Imag(tnmo_red))
                j=1
                do i=1,tnmo
                    if(.not.tRemoveCol(i)) then
                        NonSingTransformMat(:,j) = Eigenvec_Right(:,i)
                        NonSingTransformMat_Inv(j,:) = Eigenvec_Left(:,i)
                        NonSingEigenval_Real(j) = W_Real(i)
                        NonSingEigenval_Imag(j) = W_Imag(i)
                        j=j+1
                    endif
                enddo
                if(j.ne.(tnmo_red+1)) call stop_all(t_r,"Error here")

                allocate(NonSing_f_tld_row(tnmo_red))
                NonSing_f_tld_row(:) = 0.0_dp
                !Transform f
                do alpha=1,tnmo_red
                    do i=1,tnmo
                        NonSing_f_tld_row(alpha) = NonSing_f_tld_row(alpha) + f_tld_row(i)*NonSingTransformMat(i,alpha)
                    enddo
                enddo

                !Since our original matrix is now diagonal, we can calculate the coefficients easily.
                allocate(NonSingCoeffs(tnmo_red))
                do alpha=1,tnmo_red
                    NonSingCoeffs(alpha) = NonSing_f_tld_row(alpha)/NonSingEigenval_Real(alpha)
                enddo

                !Now need to transform our NonSingCoeffs back into the original basis. Multiply by NonSingTransformMat^T
                !Put these coefficients back in the C matrix in the new rotated basis
                do i=1,tnmo
                    do alpha=1,tnmo_red
                        Coeffs(betap,OccOrbRotInd(i)) = Coeffs(betap,OccOrbRotInd(i)) + NonSingCoeffs(alpha)*NonSingTransformMat_Inv(alpha,i)
                    enddo
                enddo

                deallocate(tRemoveCol,NonSingCoeffs,NonSingTransformMat,TransformMat,DiagWork,NonSingEigenval_Real)
                deallocate(NonSingTransformMat_Inv,W_Real,W_Imag,Eigenvec_Right,Eigenvec_Left,NonSingEigenval_Imag)
                deallocate(NonSing_f_tld_row)
                
            else
                !Search for zero columns/rows which would make the matrix singular.
                allocate(tRemoveCol(tnmo))
                tRemoveCol(:) = .false.
                tnmo_red = tnmo
                do i=1,tnmo
                    !Run over columns/rows
                    MaxColEl = 0.0_dp
                    MaxRowEl = 0.0_dp
                    do j=1,tnmo
                        !Run over rows - is this column completely zero?
                        if(abs(Mat(j,i)).gt.MaxColEl) MaxColEl = abs(Mat(j,i))
                        if(abs(Mat(i,j)).gt.MaxRowEl) MaxRowEl = abs(Mat(i,j))
                    enddo
                    if((MaxColEl.lt.1.0e-8_dp).or.(MaxRowEl.lt.1.0e-8_dp)) then
                        !Mark this column & row for removal
                        if(tIterative) then
                            call stop_all(t_r,"Iterative inversion not coded up to work with removal of orbitals")
                        endif
                        tRemoveCol(i) = .true.
                        tnmo_red = tnmo_red - 1
                    endif
                enddo

!                write(6,*) "tRemoveCol: "
!                do i=1,tnmo
!                    write(6,*) i,tRemoveCol(i)
!                enddo

                allocate(NonSingMat(tnmo_red,tnmo_red))
                allocate(NonSing_f_tld_row(tnmo_red))
                ColEl = 1
                do i=1,tnmo
                    if(.not.tRemoveCol(i)) then
                        !Copy column (remove the row too)
                        RowEl = 1
                        do j=1,tnmo
                            if(.not.tRemoveCol(j)) then
                                NonSingMat(RowEl,ColEl) = Mat(j,i)
                                RowEl = RowEl + 1
                            endif
                        enddo
                        NonSing_f_tld_row(ColEl) = f_tld_row(i)
                        ColEl = ColEl + 1 
                    else
                        !Removing contribution from this orbital
                        if(abs(f_tld_row(i)).gt.1.D-8) then
                            call warning(t_r,"Apparently significant contribution to relaxation energy from this orbital")
                        endif
                    endif
                enddo

                if(tDebug) then
                    write(6,*) "******************************"
                    write(6,*) "F_tld vector: "
                    do i=1,tnmo_red
                        write(6,*) NonSing_f_tld_row(i)
                    enddo
                    write(6,*) "Matrix: "
                    do i=1,tnmo_red
                        do j=1,tnmo_red
                            write(6,'(F10.5)',advance='no') NonSingMat(i,j)
                        enddo
                        write(6,*)
                    enddo
                    write(6,*) "******************************"
                endif

!                !is A symmetric
!                do k=1,tnmo_red
!                    do l=1,tnmo_red
!                        if(abs(NonSingMat(k,l)-NonSingMat(l,k)).gt.1.D-6) then
!                            write(6,*) "Non symmetric",k,l,abs(NonSingMat(k,l)-NonSingMat(l,k))
!                            write(6,*) NonSingMat(k,l),NonSingMat(l,k)
!                        endif
!                    enddo
!                enddo

                !Now we have to solve the linear equation:
                ! C_y M_yx = f_x   (or M_xy C_y = f_x)
                !Perhaps try DSGESV...?
                if(tIterative) then
                    call DSGESV(tnmo,1,Mat,tnmo,iPiv,f_tld_row,tnmo,ReturnMat,tnmo,Work,SWork,Iter,info)
!                    if(iter.le.0) then
!                        write(6,*) "Iterative solution failed: ",iter
!                    endif
                else
!                    call DGESV(tnmo,1,Mat(1:tnmo,1:tnmo),tnmo,iPiv,f_tld_row(1:tnmo),tnmo,info)
                    call DGESV(tnmo_red,1,NonSingMat(1:tnmo_red,1:tnmo_red),tnmo_red,iPiv,NonSing_f_tld_row(1:tnmo_red),tnmo_red,info)
                endif
                if(info.ne.0) then
                    if(info.gt.0) then
                        write(6,'(A,I5)') "i = ",i
                        write(6,'(A)') "U(i,i) is exactly zero. Factor U is exactly singular, so the solution could not be computed."
                    else
                        write(6,'(A,I5)') "i = ",-i
                        write(6,'(A)') "The i-th argument had an illegal value."
                    endif
                    call warning(t_r,"DGESV error in solving for C")
                    tDGESVFail = .true.
                    exit
                endif
!                if(tDGESVFail) then
!                    write(6,'(A,I4)') "DGESV fail - attempting matrix inversion for betap = ",betap
!                    !Reconstruct the matrices, and attempt to simply invert it.
!                    !Take the transpose, so it can be written as M_xy C_y = f_x
!                    !i.e. we want x fast
!                    do y=1,tnmo
!                        do x=1,tnmo
!                            Mat(x,y) = W(betap)*OneRDM(x,y) + RelaxX(y,x)
!                        enddo
!                    enddo
!                    !Also get the relevant row of f_tld. We need to make sure that we get the correct orbital, since
!                    !their index will have changed.
!                    do i=1,tntmo
!                        if(.not.IsOrigOrb(i)) cycle !We want to extract the OBS orbitals
!                        do j=1,tnmo
!                            if(OccOrbRotInd(j).eq.i) exit
!                        enddo
!                        if(j.gt.tnmo) call stop_all(t_r,"Couldn't find orbital")
!                        !j now indicates the index in the original basis
!                        f_tld_row(j) = f_tld(betap,i)
!                    enddo
!
!                    !Invert Mat...
!                    call mat_inv(Mat,MatInv)
!                    call DGEMM('N','N',tnmo,1,tnmo,1.0_dp,MatInv,tnmo,f_tld_row,tnmo,0.0_dp,f_tld_row_2,tnmo)
!                    f_tld_row(:) = f_tld_row_2(:)
!                    tDGESVFail = .false.
!                endif

                !Put these coefficients back in the C matrix in the new rotated basis
                if(tIterative) then
                    do i=1,tnmo
                        Coeffs(betap,OccOrbRotInd(i)) = ReturnMat(i)
                    enddo
                else
                    ColEl = 1
                    do i=1,tnmo
                        if(tRemoveCol(i)) then
                            Coeffs(betap,OccOrbRotInd(i)) = 0.0_dp
                        else
!                            Coeffs(betap,OccOrbRotInd(i)) = f_tld_row(i)
                            Coeffs(betap,OccOrbRotInd(i)) = NonSing_f_tld_row(ColEl)
                            ColEl = ColEl + 1
                        endif
                    enddo
                endif
            
                deallocate(NonSingMat)
                deallocate(NonSing_f_tld_row)
                deallocate(tRemoveCol)
            
            endif   !Inversion by diagonalisation or not
            
        enddo

        !Sanity check (in the rotated basis)
        if(.not.tDGESVFail) then
            do b=1,tntmo
                if(IsOrigOrb(b)) cycle
                do x=1,tntmo
                    if(.not.IsOrigOrb(x)) cycle
                    Temp = 0.0_dp

                    do y=1,tntmo
                        if(.not.IsOrigOrb(y)) cycle
                        Temp = Temp + Coeffs(b,y)*RotX(y,x)

                        do a=1,tntmo
                            if(IsOrigOrb(a)) cycle

                            Temp = Temp + RotGenFock(b,a)*Coeffs(a,y)*RotOneRDM(y,x) 
                            if((b.ne.a).and.(abs(RotGenFock(b,a)).gt.1.D-8)) then
                                call warning(t_r,"RotGenFock should be diagonal. Errors found in orbital relaxation equations")
                            endif
                        enddo
                    enddo

                    if(abs(Temp-f_tld(b,x)).gt.1.D-8) then
                        write(6,*) Temp,f_tld(b,x),Temp-f_tld(b,x)
                        call warning(t_r,"Errors found in solution to orbital relaxation equations")
                    endif
                enddo
            enddo
        endif


        !Solve for E (in the rotated basis)
        RelaxEnergy = 0.0_dp
        do aprime=1,tntmo
            !Loop over CABS in the new basis
            if(IsOrigOrb(aprime)) cycle 

            do j=1,tntmo
                if(.not.IsOrigOrb(j)) cycle

                do i=1,tntmo
                    if(.not.IsOrigOrb(i)) cycle

                    RelaxEnergy = RelaxEnergy + Coeffs(aprime,j)*RotGenFock(i,aprime)*RotOneRDM(j,i)
                enddo
            enddo
        enddo

        RelaxEnergy = RelaxEnergy * 2.0_dp  !The bb contribution.

        deallocate(iPiv,Mat,f_tld,f_tld_2,f_tld_row,Coeffs,RotGenFock,RotOneRDM,TempMat,RotX,MatInv,f_tld_row_2)
        if(tIterative) then
            deallocate(ReturnMat,Work)
            deallocate(SWork)
        endif

    end subroutine FindMRRelaxEnergy

end module relax
