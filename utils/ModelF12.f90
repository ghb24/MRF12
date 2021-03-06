PROGRAM ReadMO

    IMPLICIT NONE

    logical, parameter :: tCommutatorTrick=.true. 
    Logical, parameter :: tFindTerms=.true.         !Debugging option to write out same spin and opposite spin contributions from ALL terms.
    logical, parameter :: tWriteInts=.true.         !Logging option to write out integrals
    logical, parameter :: tReadRDMs=.true. 
    logical, parameter :: tformattedints=.true.     !Whether to read from formatted, human-readable files, or not.

    INTEGER :: HFOccOrbs    !This is the number of occupied HF orbitals
    INTEGER :: nOrbBasis    !No orbitals in orbital basis
    INTEGER :: nCABSBasis   !No orbitals in CABS
    INTEGER :: nOrbsTot     !Total orbitals in space

    !Units for MO_* 4traf files...
    integer :: unit_fg
    integer :: unit_f12
    integer :: unit_g
    integer :: unit_fock
    integer :: unit_ftf
    integer :: unit_k
    integer :: unit_fpk,unit_ff
    integer :: unit_debug
        
    real(8) , dimension(:,:,:,:) , allocatable :: V,B,Xarr,Phi
    !Full indexed arrays for storing the integrals
    real(8) , dimension(:,:,:,:) , allocatable :: fints,fgints,gints,ffints
    real(8) , dimension(:,:,:,:) , allocatable :: antisym_fg,antisym_g,antisym_f,antisym_ff

    real(8) , dimension(:,:,:,:) , allocatable :: Cumulant(:,:,:,:)

    real(8) , dimension(:,:) , allocatable :: tmat
    real(8) , dimension(:,:) , allocatable :: OneRDM 
    real(8) , dimension(:,:,:,:) , allocatable :: TwoRDM 

    real(8) , dimension(:,:) , allocatable :: fock,mo_k,fpk

    real(8) , dimension(:), allocatable :: OrbEnergies
    real(8) :: NucRep   !Nuclear repulsion energy

    real(8) :: EF12,BF12,XF12,EVSameSpin,EVOppSpin,EBSameSpin,EBOppSpin,EXSameSpin,EXOppSpin
    real(8) :: EBOSSum,EBSSSum,EXOSSum,EXSSSum

    integer :: p,q,x,y,r,s  !*PUBLIC* loop variables - get rid of these when code up properly.

!First want to find the index limits for the different basis sets.
    call find_params()
!    call WriteFCIDUMP()

    call calculateHF()

    call calculateMP2()

    call calculatePermSym()

    write(6,*) ""
    write(6,*) "Finding F12 correction..."
    write(6,*) ""
    write(6,*) ""

    CALL CalculateV()

    !Now, we need to find 2<0 | H | F_pq^xy > (3/8 \delta_px \delta_qy + 1/8 \delta_py \delta_qx)
!    write(6,*) "Applying fixed cusp coalescence conditions to V tensor for first order energy..."
    EF12=0.D0
    do p=1,nOrbBasis*2,1    !this equals x in the first term, and y in the second
        do q=1,nOrbBasis*2,1    !this equals y in the first term, and x in the second
            EF12=EF12+V(p,q,p,q)-V(p,q,q,p)
!            WRITE(6,*) EF12,V(p,q,p,q),V(p,q,q,p) 
        enddo
    enddo
    EF12=EF12/2.D0

    write(6,*) "V F12 correction term gives an energy of: ",EF12
    write(6,*) ""
    call CalcVTermEnergy(V,EVSameSpin,EVOppSpin)
    if(abs((EVOppSpin+EVSameSpin)-EF12).gt.1.D-8) then
        write(6,*) EVOppSpin,EVSameSpin,EF12
        stop 'V correction does not match sum of individual spin components'
    endif
    deallocate(V)

    call CalculateB()
!    call testperm()

!calculate contribution to the energy from B as:
!t_rs^vw B_vw^xy 2RDM_pq^rs t_xy^pq
!= B_rs^xy 2RDM_pq^rs t_xy^pq - B_sr^xy 2RDM_pq^rs t_xy^pq
!= B_rs^xy 2RDM_xy^rs - B_rs^xy 2RDM_yx^rs - B_sr^xy 2RDM_xy^rs + B_sr^xy 2RDM_yx^rs
    BF12=0.D0
    do x=1,nOrbBasis*2
        do y=1,nOrbBasis*2
            do r=1,nOrbBasis*2
                do s=1,nOrbBasis*2
                    BF12=BF12 + B(r,s,x,y)*TwoRDM(x,y,r,s) - &
                                B(r,s,x,y)*TwoRDM(x,y,s,r) - &
                                B(r,s,y,x)*TwoRDM(x,y,r,s) + &
                                B(r,s,y,x)*TwoRDM(x,y,s,r)
                enddo
            enddo
        enddo
    enddo
    BF12=BF12/16.D0

    write(6,*) ""
    if(tFindTerms) then
        write(6,*) "From summing individual energy contibutions in B:"
        write(6,*) "Same spin: ",EBSSSum
        write(6,*) "Opp spin: ",EBOSSum
        write(6,*) "Total: ",EBSSSum+EBOSSum
        write(6,*) ""
        if(abs(BF12-(EBSSSum+EBOSSum)).gt.1.D-8) then
            stop 'B correction does not match sum of individual spin components'
        endif
    endif
    write(6,*) "From full B intermediate: "
    call CalcBTermEnergy(B,EBSameSpin,EBOppSpin)

    write(6,*) "B F12 correction term gives an energy of: ",BF12
    write(6,*) ""
    EF12=EF12+BF12
    deallocate(B)

    call CalculatePhi()
!    unit_debug=get_free_unit()
!    open(unit_debug,file='Phi_Debug',status='unknown')
!    write(6,*) "Non-zero phi-terms are:"
!    do x=1,nOrbBasis*2
!        do y=1,nOrbBasis*2
!            do r=1,nOrbBasis*2
!                do s=1,nOrbBasis*2
!                    if(abs(Phi(x,y,r,s)).gt.1.D-9) then
!                        write(6,*) x,y,r,s,Phi(x,y,r,s)
!                    endif
!                enddo
!            enddo
!        enddo
!    enddo
!    close(unit_debug)

    call CalculateX()

!calculate contribution to the energy from X as:
!t_rs^vw X_vw^xy Phi_pq^rs t_xy^pq
!= X_rs^xy Phi_pq^rs t_xy^pq - X_sr^xy Phi_pq^rs t_xy^pq
!= X_rs^xy Phi_xy^rs - X_rs^xy Phi_yx^rs - X_sr^xy Phi_xy^rs + X_sr^xy Phi_yx^rs
    XF12=0.D0
    do x=1,nOrbBasis*2
        do y=1,nOrbBasis*2
            do r=1,nOrbBasis*2
                do s=1,nOrbBasis*2
                    XF12=XF12 + Xarr(r,s,x,y)*Phi(x,y,r,s) - &
                                Xarr(r,s,x,y)*Phi(y,x,r,s) - &
                                Xarr(s,r,x,y)*Phi(x,y,r,s) + &
                                Xarr(s,r,x,y)*Phi(y,x,r,s)
                enddo
            enddo
        enddo
    enddo
    XF12=XF12/16.D0
    write(6,*) ""
    write(6,*) "X F12 correction term gives an energy of: ",XF12

    write(6,*) ""
    if(tFindTerms) then
        write(6,*) "From summing individual energy contibutions in X:"
        write(6,*) "Same spin: ",EXSSSum
        write(6,*) "Opp spin: ",EXOSSum
        write(6,*) "Total: ",EXSSSum+EXOSSum
        write(6,*) ""
        if(abs(XF12-(EXSSSum+EXOSSum)).gt.1.D-8) then
            stop 'B correction does not match sum of individual spin components'
        endif
    endif
    write(6,*) "From full X intermediate: "
    call CalcXTermEnergy(Xarr,EXSameSpin,EXOppSpin)
    EF12=EF12-XF12

    write(6,"(A)")          "********************************************"
    write(6,"(A)")          "****        Final F12 energy is:        ****"
    write(6,"(A,G25.15,A)") "****    " ,EF12,                "       ****"          
    write(6,"(A)")          "********************************************"

    call CalculateRelaxation()

    contains
    
        
    SUBROUTINE CalculateV()
        integer(2) :: Maxr,Maxs,Maxx,Maxy
        integer :: ierr,Limits(2,4),a_prime,u,t     !Limits(2,4) -> start:end,Index1:Index4
        real(8) , dimension(:,:,:,:) , allocatable :: scratch,scratch2,scratch3
        integer :: gtype,ftype,i,j,k,l
        real(8) :: fintegral,gintegral,EVSS,EVOS

        write(6,*) "Obtaining 2RDM from zeroth order wavefunction..."
        call FindTwoRDM_spin()

        write(6,*) "Calculating V_pq^xy tensor..."

        !Extract 1/2 g_rs^{kappa lambda} r_{kappa lambda}^xy
        !These integrals are the FG integrals from 4traf.
        !Obtain fg_rs^xy
!        write(6,*) "Getting fg_rs^xy tensor from MO_FG file..."
        unit_fg=get_free_unit()
        if(tformattedints) then
            open(unit_fg,file='FGDUMP',status='old',form='formatted')
        else
            open(unit_fg,file='MO_FG',status='old',form='unformatted',access='sequential',action='read')
        endif

        !Fill with spatial orbital integrals 
        allocate(fgints(1:nOrbBasis,1:nOrbBasis,1:nOrbBasis,1:nOrbBasis),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'

        !Fill the fg array (in physical notation). All indices just over orbital basis
        Limits(1,1)=1
        Limits(2,1)=nOrbBasis
        Limits(1,2)=1
        Limits(2,2)=nOrbBasis
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis

        call fill4indarr(unit_fg,fgints,Limits)

!        write(6,'(A)') "Antisymmetrizing fg integrals - converting to spin-orbital representation..."
        allocate(antisym_fg(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        
        !Fill the fg array (in physical notation). All indices just over orbital basis
        Limits(1,1)=1
        Limits(2,1)=nOrbBasis*2
        Limits(1,2)=1
        Limits(2,2)=nOrbBasis*2
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis*2
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis*2

        call antisymcuspints(fgints,antisym_fg,Limits)

!        call CalcVTermRDMEnergy(antisym_fg,EVSS,EVOS)

        deallocate(fgints)

!        call testperm(antisym_fg,Limits)

        !Now we need to calulate the one-body density matrix for the zeroth order wavefunction
        call FindOneRDM_spin()

        !Find the r array, and transform one of the indices with the 1RDM
        !We want r_ua'^xy. The a' space is the CABS only space.
        !These integrals are the F integrals from 4traf.
        !Obtain r_ua'^xy.
        !However, we later want r_tu^xy, so get all indices in that slot
!        write(6,*) "Getting r_ua'^xy tensor from MO_F12 file..."
        unit_f12=get_free_unit()
        if(tformattedints) then
            open(unit_f12,file='F12DUMP',status='old',form='formatted')
        else
            open(unit_f12,file='MO_F12',status='old',form='unformatted',access='sequential',action='read')
        endif

        allocate(fints(1:norbstot,1:norbstot,1:norbbasis,1:norbbasis),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        allocate(antisym_f(1:norbstot*2,1:norbstot*2,1:norbbasis*2,1:norbbasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'

        Limits(1,1)=1
        Limits(2,1)=nOrbsTot
        Limits(1,2)=1
        Limits(2,2)=nOrbsTot
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis

        call fill4indarr(unit_f12,fints,Limits)

!        write(6,*) "R 2 1 1 2: ",fints(2,1,1,2)
!        write(6,*) "R 2 1 2 1: ",fints(2,1,2,1)
        
        Limits(1,1)=1
        Limits(2,1)=nOrbsTot*2
        Limits(1,2)=1
        Limits(2,2)=nOrbsTot*2
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis*2
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis*2
        
!        write(6,*) "antisymmetrizing f integrals"
        call antisymcuspints(fints,antisym_f,Limits)
!        open(24,file="FR-Model",status='unknown')
!        do i=1,nOrbBasis
!            do j=1,nOrbBasis
!                do k=1,nOrbsTot
!                    do l=1,nOrbsTot
!!                        write(24,"(4I4,F15.8)") (i+1)/2,(j+1)/2,(k+1)/2,(l+1)/2,antisym_f(i,j,k,l)
!                        write(24,"(G25.10,2I5)") fints(k,l,i,j),i,j
!                    enddo
!                enddo
!            enddo
!        enddo
!        close(24)
        deallocate(fints)

        !Now contract with 1RDM
!        write(6,*) "Performing contraction of F12 integrals with 1RDM..."
!        write(6,*) "1RDM_t^u r_ua'^xy"

        allocate(scratch(1:nOrbBasis*2,(nOrbBasis*2)+1:(nOrbsTot*2),1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        scratch=0.D0    !This will hold X_ta'^xy

        do t=1,nOrbBasis*2
            do a_prime=(nOrbBasis*2)+1,(nOrbsTot*2)
                do x=1,(nOrbBasis*2)
                    do y=1,(nOrbBasis*2)
                        do u=1,(nOrbBasis*2)
                            scratch(t,a_prime,x,y)=scratch(t,a_prime,x,y)+OneRDM(t,u)*antisym_f(u,a_prime,x,y)
                        enddo
                    enddo
                enddo
            enddo
        enddo

!        write(6,*) "Getting g_rs^ta' tensor from MO_G file..."
        unit_g=get_free_unit()
        if(tformattedints) then
            open(unit_g,file='FCIDUMP',status='old',form='formatted')
        else
            open(unit_g,file='MO_G',status='old',form='unformatted',access='sequential',action='read')
        endif

        allocate(gints(1:nOrbsTot,1:nOrbsTot,1:nOrbsTot,1:nOrbBasis),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        allocate(antisym_g(1:norbsTot*2,1:norbsTot*2,1:norbsTot*2,1:norbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'

        Limits(1,1)=1
        Limits(2,1)=nOrbsTot
        Limits(1,2)=1
        Limits(2,2)=nOrbsTot
        Limits(1,3)=1
        Limits(2,3)=nOrbsTot
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis

        call fill4indarr(unit_g,gints,Limits)

        Limits(1,1)=1
        Limits(2,1)=nOrbsTot*2
        Limits(1,2)=1
        Limits(2,2)=nOrbsTot*2
        Limits(1,3)=1
        Limits(2,3)=nOrbsTot*2
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis*2

!        write(6,*) "antisymmetrizing g integrals..."
        call antisymints(gints,antisym_g,Limits)
        deallocate(gints)

!        write(6,*) "Performing double contraction of g_rs^ta' X_ta'^xy..."

        allocate(scratch2(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        scratch2=0.D0

        do r=1,nOrbBasis*2
            do s=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do t=1,nOrbBasis*2
                            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
!                                scratch2(r,s,x,y)=scratch2(r,s,x,y)+antisym_g(r,s,t,a_prime)*scratch(t,a_prime,x,y)
                                scratch2(r,s,x,y)=scratch2(r,s,x,y)+antisym_g(t,a_prime,r,s)*scratch(t,a_prime,x,y)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        Limits(1,1)=1
        Limits(2,1)=nOrbBasis*2
        Limits(1,2)=1
        Limits(2,2)=nOrbBasis*2
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis*2
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis*2

!        call testperm(scratch2,Limits)

        !scratch2 now is the second tensor in the V_pq^xy expression.
        deallocate(scratch)
!        write(6,*) "calculating 1/2 g_rs^tu r_ru^xy..."
        !Why don't we want to just get this from fg? 
        allocate(scratch(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        scratch=0.D0

        do r=1,nOrbBasis*2
            do s=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do t=1,nOrbBasis*2
                            do u=1,nOrbBasis*2

                                scratch(r,s,x,y)=scratch(r,s,x,y)+(antisym_g(r,s,t,u)*antisym_f(t,u,x,y))/2.D0


!!    !Return 1 if the integral is of spin type <aa|aa>
!!    !Return 2 if the integral is of spin type <ab|ab>
!!    !Return 3 if the integral is of spin type <ab|ba>
!!    !Return 4 if the integral is of spin type <aa|bb> or <aa|ab> or <aa|ba> ...
!!    integer function DetermineSpinType(i,j,k,l)
!!Below are attempts to avoid the need for a spin-orbital representation - needs work.
!                                gtype=DetermineSpinType(r,s,t,u)
!                                if(gtype.eq.1) then
!                                    gintegral=gints(Conv2Spat(r),Conv2Spat(s),Conv2Spat(t),Conv2Spat(u))-gints(Conv2Spat(r),Conv2Spat(s),Conv2Spat(u),Conv2Spat(t))
!                                elseif(gtype.eq.2) then
!                                    gintegral=gints(Conv2Spat(r),Conv2Spat(s),Conv2Spat(t),Conv2Spat(u))
!                                elseif(gtype.eq.3) then
!                                    gintegral=-gints(Conv2Spat(r),Conv2Spat(s),Conv2Spat(u),Conv2Spat(t))
!                                else
!                                    gintegral=0.D0
!                                endif
!                                ftype=DetermineSpinType(t,u,x,y)
!                                if(ftype.eq.1) then
!                                    fintegral=fints(Conv2Spat(t),Conv2Spat(u),Conv2Spat(x),Conv2Spat(y))-gints(Conv2Spat(t),Conv2Spat(u),Conv2Spat(y),Conv2Spat(x))
!                                elseif(ftype.eq.2) then
!                                    fintegral=fints(Conv2Spat(t),Conv2Spat(u),Conv2Spat(x),Conv2Spat(y))
!                                elseif(ftype.eq.3) then
!                                    fintegral=-fints(Conv2Spat(t),Conv2Spat(u),Conv2Spat(y),Conv2Spat(x))
!                                else
!                                    fintegral=0.D0
!                                endif
!                                scratch(r,s,x,y)=scratch(r,s,x,y)+(gintegral*fintegral)/2.D0



                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        allocate(scratch3(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        scratch3=0.D0

        !Now, take scratch and scratch2 away from fgints
        do r=1,nOrbBasis*2
            do s=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        scratch3(r,s,x,y)=antisym_fg(r,s,x,y)-scratch(r,s,x,y)-scratch2(r,s,x,y)
                    enddo
                enddo
            enddo
        enddo

        !Finally, we need to contract scratch3 with the 2RDM, to obtain V_pq^xy

!        write(6,*) "Contracting 2RDM to find V_pq^xy"
        allocate(V(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
        V=0.D0
        if(ierr.ne.0) stop 'error allocating'

        !Calculate individual term energies.
        do p=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do r=1,nOrbBasis*2
                            do s=1,nOrbBasis*2
                                V(p,q,x,y) = V(p,q,x,y) + (antisym_fg(r,s,x,y)*TwoRDM(p,q,r,s))/2.D0
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        call CalcVTermEnergy(V,EVSS,EVOS,1)
        V=0.D0
        do p=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do r=1,nOrbBasis*2
                            do s=1,nOrbBasis*2
                                V(p,q,x,y) = V(p,q,x,y) - (scratch2(r,s,x,y)*TwoRDM(p,q,r,s))/2.D0
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        call CalcVTermEnergy(V,EVSS,EVOS,2)
        V=0.D0
        do p=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do r=1,nOrbBasis*2
                            do s=1,nOrbBasis*2
                                V(p,q,x,y) = V(p,q,x,y) - (scratch(r,s,x,y)*TwoRDM(p,q,r,s))/2.D0
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        call CalcVTermEnergy(V,EVSS,EVOS,3)
        V=0.D0


        do p=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do r=1,nOrbBasis*2      !Eventually, this wants to be over all orbital basis...
                            do s=1,nOrbBasis*2
                                V(p,q,x,y)=V(p,q,x,y)+(scratch3(r,s,x,y)*TwoRDM(p,q,r,s))/2.D0
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        deallocate(scratch3)
        deallocate(scratch,scratch2)

    end subroutine CalculateV

    subroutine CalculatePhi()
        real(8) , allocatable :: Z(:,:,:,:),TwoIndScratch(:,:),TwoIndScratch2(:,:)
        real(8) , allocatable :: FourIndScratch(:,:,:,:),Test(:,:,:,:),Test2(:,:)
        integer :: u,t
        integer :: ierr

!        allocate(Test2(nOrbBasis*2,nOrbBasis*2))
!        Test2 = 0.D0

        !Phi is a little tricky.
        !We have the two commutator permutation operators acting on the electrons in each geminal function.
        !If we calculate the tensor, Z_pq^rs, then Phi_pq^rs = Z_pq^rs - Z_qp^rs - Z_pq^sr + Z_qp^sr
        !NOTE: Does this not just equal 2[ Z_pq^rs - Z_pq^sr ] : perhaps not - test this.

        allocate(Z(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        allocate(Phi(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        Z=0.D0
        Phi=0.D0

        allocate(TwoIndScratch(nOrbBasis*2,nOrbBasis*2),stat=ierr)
        allocate(TwoIndScratch2(nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        TwoIndScratch=0.D0
        TwoIndScratch2=0.D0

        do t=1,nOrbBasis*2
            do s=1,nOrbBasis*2
                do u=1,nOrbBasis*2
                    TwoIndScratch(t,s)=TwoIndScratch(t,s)+fock(t,u)*OneRDM(u,s)  !f_t^u 1RDM_u^s
                enddo
            enddo
        enddo

        do q=1,nOrbBasis*2
            do s=1,nOrbBasis*2
                do t=1,nOrbBasis*2
                    TwoIndScratch2(q,s)=TwoIndScratch2(q,s)+OneRDM(q,t)*TwoIndScratch(t,s)  !1RDM_q^t f_t^u 1RDM_u^s
                enddo
            enddo
        enddo

        do p=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
                        Z(p,q,r,s)=Z(p,q,r,s)+OneRDM(p,r)*TwoIndScratch2(q,s)
                    enddo
                enddo
            enddo
        enddo
        deallocate(TwoIndScratch2)
        !Calculated Term 1

        allocate(FourIndScratch(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        FourIndScratch=0.D0

        !Now for the cumulants...
        allocate(Cumulant(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        Cumulant=0.D0

        write(6,*) "Calculating 2-cumulant..."
        do u=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
                        Cumulant(u,q,r,s)=TwoRDM(u,q,r,s) + OneRDM(u,s)*OneRDM(q,r) - OneRDM(u,r)*OneRDM(q,s)
                    enddo
                enddo
            enddo
        enddo


        
        do p=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do r=1,nOrbBasis*2
                    do t=1,nOrbBasis*2
                        do u=1,nOrbBasis*2
                            FourIndScratch(p,q,r,t)=FourIndScratch(p,q,r,t)+fock(u,t)*Cumulant(p,q,r,u)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        do p=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
                        do t=1,nOrbBasis*2
                            Z(p,q,r,s)=Z(p,q,r,s)+(OneRDM(t,s)*FourIndScratch(p,q,r,t))/2.D0
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !Included Term 2

        FourIndScratch=0.D0

        do t=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
                        do u=1,nOrbBasis*2
                            FourIndScratch(t,q,r,s)=FourIndScratch(t,q,r,s)+fock(t,u)*Cumulant(u,q,r,s)
                        enddo
                    enddo
                enddo
            enddo
        enddo

!        !BooBug
!        allocate(Test(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2))
!        Test = 0.D0

        do p=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
                        do t=1,nOrbBasis*2
                            Z(p,q,r,s)=Z(p,q,r,s)+(OneRDM(p,t)*FourIndScratch(t,q,r,s))/2.D0
!                            test(p,q,r,s) = test(p,q,r,s) +(OneRDM(p,t)*FourIndScratch(t,q,r,s))/2.D0
                        enddo
                    enddo
                enddo
            enddo
        enddo

!        deallocate(Test)

        !Final term!
        TwoIndScratch=0.D0
        do s=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do t=1,nOrbBasis*2
                    do u=1,nOrbBasis*2
                        TwoIndScratch(q,s)=TwoIndScratch(q,s)+fock(t,u)*Cumulant(u,q,t,s)
                    enddo
                enddo
            enddo
        enddo

        do p=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
                        Z(p,q,r,s)=Z(p,q,r,s) - (OneRDM(p,r)*TwoIndScratch(q,s))
                    enddo
                enddo
            enddo
        enddo

        write(6,*) "Z calculated. Applying permutation operators to obtain Phi"



        do p=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
                        Phi(p,q,r,s)=Z(p,q,r,s)-Z(q,p,r,s)-Z(p,q,s,r)+Z(q,p,s,r)
                    enddo
                enddo
            enddo
        enddo

        write(6,*) "Phi tensor obtained."

!        do p=1,nOrbBasis*2
!            do q=1,nOrbBasis*2
!                do r=1,nOrbBasis*2
!                    do s=1,nOrbBasis*2
!                        if(abs(Phi(p,q,r,s)-Phi(r,s,p,q)).gt.1.D-8) then
!                            write(6,*) "Phi NOT symmetric",Phi(p,q,r,s),Phi(r,s,p,q)
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo
!        
!        open(24,file="FR-Model_aaaa",status='unknown')
!        open(25,file="FR-Model_abab",status='unknown')
!        open(26,file="FR-Model_abba",status='unknown')
!        do p=1,nOrbBasis*2,2
!            do q=1,nOrbBasis*2,2
!                do r=1,nOrbBasis*2,2
!                    do s=1,nOrbbasis*2,2
!                        write(24,"(G25.10,4I5)") Phi(p,q,r,s),(p+1)/2,(q+1)/2,(r+1)/2,(s+1)/2
!                    enddo
!                enddo
!            enddo
!        enddo
!        do p=1,nOrbBasis*2,2
!            do q=2,nOrbBasis*2,2
!                do r=1,nOrbBasis*2,2
!                    do s=2,nOrbbasis*2,2
!                        write(25,"(G25.10,4I5)") Phi(p,q,r,s),(p+1)/2,(q+1)/2,(r+1)/2,(s+1)/2
!                    enddo
!                enddo
!            enddo
!        enddo
!        do p=1,nOrbBasis*2,2
!            do q=2,nOrbBasis*2,2
!                do r=2,nOrbBasis*2,2
!                    do s=1,nOrbbasis*2,2
!                        write(26,"(G25.10,4I5)") Phi(p,q,r,s),(p+1)/2,(q+1)/2,(r+1)/2,(s+1)/2
!                    enddo
!                enddo
!            enddo
!        enddo
!        do p=1,nOrbBasis*2,2
!            do q=2,nOrbBasis*2,2
!                do r=2,nOrbBasis*2,2
!                    do s=2,nOrbBasis*2,2
!                        if(abs(Phi(p,q,r,s)).gt.1.D-8) stop 'Spin Error'
!                    enddo
!                enddo
!            enddo
!        enddo


    end subroutine CalculatePhi

    subroutine CalculateX()
        integer :: ierr,uu,t,w,u,a_prime
        real(8) , allocatable :: scratch(:,:,:,:)
        real(8) , allocatable :: TermWise(:,:,:,:)
        real(8) :: EXOS,EXSS
!        open(24,file="FR-Model",status='unknown')

        write(6,*) "Calculating X_uw^xy tensor..."
        allocate(Xarr(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        Xarr=0.D0

        if(tFindTerms) then
            EXOSSum=0.D0
            EXSSSum=0.D0
            allocate(TermWise(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
            if(ierr.ne.0) stop 'error allocating'
            TermWise=0.D0
        endif


        !Since the first term is simply the antisymmetrised ff integrals, we can copy them across
        Xarr(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2)=  &
            antisym_ff(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2)

!        write(6,*) "1 1 1 1: ", Xarr(1,2,1,2)
!        write(6,*) "1 2 1 2: ", Xarr(1,4,1,4)

        deallocate(antisym_ff)  !we no longer need these integrals.

        if(tFindTerms) then
            TermWise=Xarr
            call CalcXTermEnergy(TermWise,EXSS,EXOS,1)
            EXSSSum=EXSSSum+EXSS
            EXOSSum=EXOSSum+EXOS
            TermWise=0.D0
        endif

        allocate(scratch(1:nOrbBasis*2,(nOrbBasis*2)+1:nOrbsTot*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        scratch=0.D0

        do t=1,nOrbBasis*2
            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do u=1,nOrbBasis*2
                            scratch(t,a_prime,x,y)=scratch(t,a_prime,x,y)+OneRDM(t,u)*antisym_f(u,a_prime,x,y)
                        enddo
                    enddo
                enddo
            enddo
        enddo

!        write(6,*) "Double Contraction 1"

        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do t=1,nOrbBasis*2
                            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
                                Xarr(u,w,x,y)=Xarr(u,w,x,y)-antisym_f(t,a_prime,u,w)*scratch(t,a_prime,x,y)
                                if(tFindTerms) then
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)-antisym_f(t,a_prime,u,w)*scratch(t,a_prime,x,y)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        deallocate(scratch)
        if(tFindTerms) then
            call CalcXTermEnergy(TermWise,EXSS,EXOS,2)
            EXSSSum=EXSSSum+EXSS
            EXOSSum=EXOSSum+EXOS
            TermWise=0.D0
        endif
        
!        write(6,*) "Double Contraction 2"

        !Now for the final term - double contraction over orbital space.
        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do t=1,nOrbBasis*2
                            do uu=1,nOrbBasis*2
                                Xarr(u,w,x,y)=Xarr(u,w,x,y)-(antisym_f(u,w,t,uu)*antisym_f(t,uu,x,y))/2.D0
                                if(tFindTerms) then
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)-(antisym_f(u,w,t,uu)*antisym_f(t,uu,x,y))/2.D0
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if(tFindTerms) then
            call CalcXTermEnergy(TermWise,EXSS,EXOS,3)
            EXSSSum=EXSSSum+EXSS
            EXOSSum=EXOSSum+EXOS
            deallocate(TermWise)
        endif

        write(6,*) "Calculated all terms in X"

    end subroutine CalculateX


    subroutine CalculateB()
    integer :: ierr,TwoIndLims(2),Limits(2,4),i,j,k,l
    real(8) , allocatable :: OneIndOrbSumFockR(:,:,:,:),ftfints(:,:,:,:),scratch(:,:,:,:),scratch2(:,:,:,:)
    real(8) , allocatable :: TermWise(:,:,:,:),TermWise2(:,:,:,:)
    real(8) , allocatable :: TermWise3(:,:,:,:),TermWise4(:,:,:,:)
    real(8) , allocatable :: Test(:,:,:,:),Test2(:,:,:,:),Test3(:,:,:,:)
    integer :: p_prime,u,w,SpinType,a_prime,q_prime,b_prime,r_prime,c_prime,count
    real(8) :: EBSS,EBOS,ETMP,Val

        write(6,*) "Calculating B_uw^xy tensor..."
        allocate(B(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        B=0.D0
        if(tFindTerms) then
            allocate(TermWise(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
            if(ierr.ne.0) stop 'error allocating'
            allocate(TermWise2(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
            if(ierr.ne.0) stop 'error allocating'
            TermWise=0.D0
            TermWise2=0.D0
            EBSSSum=0.D0
            EBOSSum=0.D0
        endif

        !The antisymmetrised f_PQ^xy integrals are already stored in antisym_f
        !For three of the terms, the *FIRST* index wants to be transformed into a 
        !fock basis. Read in the fock basis.
        !Fill with spin orbital integrals 
        allocate(fock(1:nOrbsTot*2,1:nOrbsTot*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        TwoIndLims(1)=1
        TwoIndLims(2)=nOrbsTot*2

        if(.not.tReadRDMs) then
            unit_fock=get_free_unit()
            if(tformattedints) then
                open(unit_fock,file='MO_F',status='old',form='formatted')
            else
                open(unit_fock,file='MO_F',status='old',form='unformatted',access='sequential',action='read')
            endif

            !Fill the fg array (in physical notation). All indices just over orbital basis
            write(6,*) "Reading fock matrix from disk."
            call fill2indarr_spin(unit_fock,fock,TwoIndLims)
            !Now we need the exchange matrix only
            !Now get the k matrix
            unit_k=get_free_unit()
            if(tformattedints) then
                open(unit_k,file='MO_K',status='old',form='formatted')
            else
                open(unit_k,file='MO_K',status='old',form='unformatted',access='sequential',action='read')
            endif

            !Fill with spin orbital integrals 
            allocate(mo_k(1:nOrbsTot*2,1:nOrbsTot*2),stat=ierr)
            if(ierr.ne.0) stop 'error allocating'

            !Fill the fg array (in physical notation). All indices just over orbital basis
            TwoIndLims(1)=1
            TwoIndLims(2)=nOrbsTot*2
            write(6,*) "Reading exchange matrix from disk."
            call fill2indarr_spin(unit_k,mo_k,TwoIndLims)
!            write(6,*) "K matrix read in"
        else
            write(6,*) "Calculating *generalised* fock matrices..."
            call CalcTMat(fock,TwoIndLims)
            call CalcGenFockMats(fock,TwoIndLims)
        endif

!        write(6,*) "Fock matrix read in - is this correct in the spin-orbital formulation?"

        write(19,*) "Fock matrix: "
        do p=1,nOrbsTot*2
            do q=1,nOrbsTot*2
                write(19,"(G15.8)",advance='no')    fock(p,q)
            enddo
            write(19,"(A)") " " 
        enddo
!
!        write(20,*) "ONERDM matrix: "
!        do p=1,nOrbBasis*2
!            do q=1,nOrbBasis*2
!                write(20,"(G15.8)",advance='no')    OneRDM(p,q)
!            enddo
!            write(20,"(A)") " " 
!        enddo

        !For the first two terms, transform the first index according to the fock matrix
        !summing over the ORBITAL space
!        write(6,*) "Calculating first-index fock-transformed F12 integrals."
        allocate(OneIndOrbSumFockR(nOrbsTot*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        OneIndOrbSumFockR=0.D0
        do p_prime=1,nOrbsTot*2     !Allow the transformed index to run over whole of space.
            do q=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do p=1,nOrbBasis*2  !Only contracted over the orbital space.
                            OneIndOrbSumFockR(p_prime,q,x,y)=OneIndOrbSumFockR(p_prime,q,x,y)+fock(p_prime,p)*antisym_f(p,q,x,y)
                        enddo
                    enddo
                enddo
            enddo
        enddo

!        allocate(Test(nOrbsTot*2,nOrbBasis*2),stat=ierr)
!        if(ierr.ne.0) stop 'error allocating'
!        Test=0.D0
!        allocate(Test2(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,(nOrbBasis*2+1):nOrbsTot*2),stat=ierr)
!        if(ierr.ne.0) stop 'error allocating'
!        Test2=0.D0
!!        allocate(Test3(nOrbBasis*2,(nOrbBasis*2+1):nOrbsTot*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
!!        if(ierr.ne.0) stop 'error allocating'
!!        Test3=0.D0
!
!        do p_prime=1,nOrbsTot*2
!            do q=1,nOrbBasis*2
!                do p=1,nOrbBasis*2
!                    Test(p_prime,q) = Test(p_prime,q) + fock(p_prime,p)*OneRDM(p,q)
!                enddo
!            enddo
!        enddo
!
!        do u=1,nOrbBasis*2
!            do w=1,nOrbBasis*2
!                do q=1,nOrbBasis*2
!                    do a_prime=(nOrbBasis*2+1),nOrbsTot*2
!                        do p_prime=1,nOrbsTot*2
!                            Test2(u,w,q,a_prime) = Test2(u,w,q,a_prime) + antisym_f(p_prime,a_prime,u,w)*Test(p_prime,q)
!                        enddo
!                    enddo
!                enddo
!            enddo
!        enddo

!
!        do q=1,nOrbBasis*2
!            do b_prime=nOrbBasis*2+1,nOrbsTot*2
!                do u=1,nOrbBasis*2
!                    do w=1,nOrbBasis*2
!                        do x=1,nOrbBasis*2
!                            do y=1,nOrbBasis*2
!                                Test3(q,b_prime,u,w) = Test3(q,b_prime,u,w) + antisym_f(q,b_prime,x,y) * TwoRDM(x,y,u,w)
!                            enddo
!                        enddo
!                    enddo
!                enddo
!            enddo
!        enddo
!
!        do u=1,nOrbBasis*2
!            do w=1,nOrbBasis*2
!                do q=1,nOrbBasis*2
!                    do a_prime=nOrbBasis*2+1,nOrbsTot*2
!                        do p=1,nOrbBasis*2
!                            Test(u,w,q,a_prime) = Test(u,w,q,a_prime) + antisym_f(p,a_prime,u,w) * OneRDM(p,q)
!!                            if(u.eq.1.and.w.eq.3.and.q.eq.1.and.a_prime.eq.5) then
!!                                write(6,*) Test(u,w,q,a_prime),antisym_f(p,a_prime,u,w),OneRDM(p,q)
!!                            endif
!                        enddo
!                    enddo
!                enddo
!            enddo
!        enddo
!
!        do u=1,nOrbBasis*2
!            do w=1,nOrbBasis*2
!                do q=1,nOrbBasis*2
!                    do b_prime=(nOrbBasis*2+1),nOrbsTot*2
!                        do a_prime=(nOrbBasis*2+1),nOrbsTot*2
!                            Test2(u,w,q,b_prime) = Test2(u,w,q,b_prime) + Test(u,w,q,a_prime)*fock(a_prime,b_prime)
!!                            if(u.eq.1.and.w.eq.3.and.q.eq.1.and.b_prime.eq.5) then
!!                                write(6,*) Test2(u,w,q,b_prime),Test(u,w,q,a_prime),fock(a_prime,b_prime)
!!                            endif
!                        enddo
!                    enddo
!                enddo
!            enddo
!        enddo
!
!
!        open(24,file="FR-Model",status='unknown')
!
!        do u=1,nOrbsTot*2,2
!            do w=2,nOrbsTot*2,2
!                do q=1,nOrbBasis*2,2
!                    do a_prime=(nOrbBasis*2+1),nOrbsTot*2,2
!                do p_prime=2,nOrbsTot*2,2
!                        write(24,"(G25.10,2I5)") Test2(u,w,q,a_prime),(u+1)/2,(w+1)/2
!                        write(24,"(G25.10,2I5)") antisym_f(p_prime,a_prime,u,w),(u+1)/2,(w+1)/2
!                 write(24,"(G25.10,2I5)") mo_k(u,w),(u+1)/2,(w+1)/2
!                    enddo
!                enddo
!            enddo
!        enddo

!        do q=1,nOrbBasis*2,2
!            do p_prime=1,nOrbsTot*2,2
!                write(24,"(G25.10,2I5)") Test(p_prime,q),p_prime,q
!            enddo
!        enddo
!
!        ETMP = 0.D0
!        do u=1,nOrbBasis*2,2
!            do w=1,nOrbBasis*2,2
!                do q=1,nOrbBasis*2,2
!                    do b_prime=(nOrbBasis*2+1),nOrbsTot*2,2
!                        ETMP = ETMP - (Test3(q,b_prime,u,w)*Test2(u,w,q,b_prime))/2.D0
!!                        write(24,"(G25.10,2I5)") Antisym_f(p,a_prime,u,w),(u+1)/2,(w+1)/2
!                        write(24,"(G25.10,2I5)") Test2(u,w,q,b_prime),(u+1)/2,(w+1)/2
!                    enddo
!                enddo
!                do q=1,nOrbBasis*2,2
!                    do b_prime=(nOrbBasis*2+1),nOrbsTot*2,2
!!                        ETMP = ETMP - (Test3(q+1,b_prime-1,u,w)*Test2(u,w,q,b_prime))/2.D0
!!                        write(24,"(G25.10,2I5)") Antisym_f(p,a_prime,u,w),(u+1)/2,(w+1)/2
!                        write(24,"(G25.10,2I5,A)") Test3(q,b_prime,u,w),(u+1)/2,(w+1)/2," RTLD"
!                    enddo
!                enddo
!
!            enddo
!        enddo
!        write(6,*) "ETMP: ",ETMP

!        do u=1,nOrbBasis*2,2
!            do w=1,nOrbBasis*2,2
!                do q=1,nOrbBasis*2,2
!                    do a_prime=nOrbBasis*2+1,nOrbsTot*2,2
!                        write(24,"(G25.10,2I5)") Test(a_prime,q,u,w),(u+1)/2,(w+1)/2
!                    enddo
!                enddo
!            enddo
!        enddo
!        do i=1,nOrbBasis*2,2
!            do j=1,nOrbBasis*2,2
!        count=0
!        do q=1,nOrbBasis*2,2
!            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2,2
!                count=count+1
!                write(24,"(I5,G25.10,4I5)") count,antisym_f(a_prime,q,i,j),(i+1)/2,(j+1)/2,(a_prime+1)/2,(q+1)/2
!            enddo
!        enddo
!        enddo
!    enddo
        
!    do i=1,nOrbBasis*2,2
!        do j=2,nOrbBasis*2,2
!        count=0
!        do q=2,nOrbBasis*2,2
!            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2,2
!                count=count+1
!                write(24,"(I5,G25.10,4I5)") count,antisym_f(a_prime,q,i,j),(i+1)/2,j/2,(a_prime+1)/2,q/2
!            enddo
!        enddo
!        enddo
!    enddo
!    do i=2,nOrbBasis*2,2
!        do j=1,nOrbBasis*2,2
!        count=0
!        do q=1,nOrbBasis*2,2
!            do a_prime=(nOrbBasis*2)+2,nOrbsTot*2,2
!                count=count+1
!                write(24,"(I5,G25.10,4I5)") count,antisym_f(a_prime,q,i,j),i/2,(j+1)/2,a_prime/2,(q+1)/2
!            enddo
!        enddo
!        enddo
!    enddo
!    do i=1,nOrbBasis*2,2
!        do j=2,nOrbBasis*2,2
!        count=0
!        do q=2,nOrbBasis*2,2
!            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2,2
!                count=count+1
!                write(24,"(I5,G25.10,4I5)") count,antisym_f(a_prime,q,i,j),i,j,a_prime,q
!            enddo
!        enddo
!        enddo
!    enddo
!        write(24,*) "a a a a"
!       do i=1,nOrbBasis*2,2
!           do j=1,nOrbBasis*2,2
!               do k=1,nOrbBasis*2,2
!                    do l=1,nOrbBasis*2,2
!                        if(abs(Test(i,j,k,l)).gt.1.D-9) then
!                            write(24,"(4I4,F15.8)") (i+1)/2,(j+1)/2,(k+1)/2,(l+1)/2,Test(i,j,k,l)
!                        endif
!!                        write(24,"(4I4,F15.8)") (i+1)/2,(j+1)/2,(k+1)/2,(l+1)/2,antisym_f(i,j,k,l)
!                        write(24,"(4I4,F15.8)") i,j,k,l,antisym_f(i,j,k,l)
!                    enddo
!                enddo
!            enddo
!        enddo
!        write(24,*) "************ a a a b"
!        do i=1,nOrbBasis*2,2
!            do j=1,nOrbBasis*2,2
!                do k=1,nOrbBasis*2,2
!                    do l=2,nOrbBasis*2,2
!                        if(abs(Test(i,j,k,l)).gt.1.D-9) then
!                            write(24,"(4I4,F15.8)") (i+1)/2,(j+1)/2,(k+1)/2,(l)/2,Test(i,j,k,l)
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo
!        write(24,*) "************ a a b a"
!        do i=1,nOrbBasis*2,2
!            do j=1,nOrbBasis*2,2
!                do k=2,nOrbBasis*2,2
!                    do l=1,nOrbBasis*2,2
!                        if(abs(Test(i,j,k,l)).gt.1.D-9) then
!                            write(24,"(4I4,F15.8)") (i+1)/2,(j+1)/2,(k)/2,(l+1)/2,Test(i,j,k,l)
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo
!        write(24,*) "************ a b a a"
!        do i=1,nOrbBasis*2,2
!            do j=2,nOrbBasis*2,2
!                do k=1,nOrbBasis*2,2
!                    do l=1,nOrbBasis*2,2
!                        if(abs(Test(i,j,k,l)).gt.1.D-9) then
!                            write(24,"(4I4,F15.8)") (i+1)/2,(j)/2,(k+1)/2,(l+1)/2,Test(i,j,k,l)
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo
!        write(24,*) "************ b a a a"
!        do i=2,nOrbBasis*2,2
!            do j=1,nOrbBasis*2,2
!                do k=1,nOrbBasis*2,2
!                    do l=1,nOrbBasis*2,2
!                        if(abs(Test(i,j,k,l)).gt.1.D-9) then
!                            write(24,"(4I4,F15.8)") (i)/2,(j+1)/2,(k+1)/2,(l+1)/2,Test(i,j,k,l)
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo
!        write(24,*) "************ b b a a"
!        do i=2,nOrbBasis*2,2
!            do j=2,nOrbBasis*2,2
!                do k=1,nOrbBasis*2,2
!                    do l=1,nOrbBasis*2,2
!                        if(abs(Test(i,j,k,l)).gt.1.D-9) then
!                            write(24,"(4I4,F15.8)") (i)/2,(j)/2,(k+1)/2,(l+1)/2,Test(i,j,k,l)
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo
!        write(24,*) "************ b a b a"
!        do i=2,nOrbBasis*2,2
!            do j=1,nOrbBasis*2,2
!                do k=2,nOrbBasis*2,2
!                    do l=1,nOrbBasis*2,2
!                        if(abs(Test(i,j,k,l)).gt.1.D-9) then
!                            write(24,"(4I4,F15.8)") (i)/2,(j+1)/2,(k)/2,(l+1)/2,Test(i,j,k,l)
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo
!        write(24,*) "************ b a a b"
!        do i=2,nOrbBasis*2,2
!            do j=1,nOrbBasis*2,2
!                do k=1,nOrbBasis*2,2
!                    do l=2,nOrbBasis*2,2
!                        if(abs(Test(i,j,k,l)).gt.1.D-9) then
!                            write(24,"(4I4,F15.8)") (i)/2,(j+1)/2,(k+1)/2,(l)/2,Test(i,j,k,l)
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo
!        write(24,*) "************ a b b a"
!        do i=1,nOrbBasis*2,2
!            do j=2,nOrbBasis*2,2
!                do k=2,nOrbBasis*2,2
!                    do l=1,nOrbBasis*2,2
!                        if(abs(Test(i,j,k,l)).gt.1.D-9) then
!                            write(24,"(4I4,F15.8)") (i+1)/2,(j)/2,(k)/2,(l+1)/2,Test(i,j,k,l)
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo
!        write(24,*) "************ a b a b"
!        do i=1,nOrbBasis*2,2
!            do j=2,nOrbBasis*2,2
!                do k=1,nOrbBasis*2,2
!                    do l=2,nOrbBasis*2,2
!                        if(abs(Test(i,j,k,l)).gt.1.D-9) then
!                            write(24,"(4I4,F15.8)") (i+1)/2,(j)/2,(k+1)/2,(l)/2,Test(i,j,k,l)
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo
!
!        close(24)

        !Now calculate the term -r_uw^rq f_r^p r_pq^xy - put it straight into the B array
!        write(6,*) "Performing double-contraction to obtain B-term: -r_uw^rq f_r^p r_pq^xy..."
        !TERM 1
        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        !Sum over r & q
                        do r=1,nOrbBasis*2
                            do q=1,nOrbBasis*2
                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(u,w,r,q)*OneIndOrbSumFockR(r,q,x,y)
                                if(tFindTerms) then
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)-antisym_f(u,w,r,q)*OneIndOrbSumFockR(r,q,x,y)
                                endif
                            enddo
                        enddo
                        !Seperately, calculate TERM 2!
                        !These can easily be combined later into one term.
                        do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
                            do q=1,nOrbBasis*2
!                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(u,w,a_prime,q)*OneIndOrbSumFockR(a_prime,q,x,y)
!Use perm sym:
                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(a_prime,q,u,w)*OneIndOrbSumFockR(a_prime,q,x,y)
                                if(tFindTerms) then
                                    TermWise2(u,w,x,y)=TermWise2(u,w,x,y)-antisym_f(a_prime,q,u,w)*OneIndOrbSumFockR(a_prime,q,x,y)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        deallocate(OneIndOrbSumFockR)
        if(tFindTerms) then
            call CalcBTermEnergy(TermWise,EBSS,EBOS,1)
            EBSSSum=EBSSSum+EBSS
            EBOSSum=EBOSSum+EBOS
            call CalcBTermEnergy(TermWise2,EBSS,EBOS,2)
            EBSSSum=EBSSSum+EBSS
            EBOSSum=EBOSSum+EBOS
            TermWise=0.D0
            TermWise2=0.D0
        endif
        write(6,*) "Summed in TERMS 1 and 2"



        !First find f_p^a' r_a'q^xy
!        write(6,*) "Calculating first-index fock-transformed F12 integrals, transformed over the CABS space."
        allocate(OneIndOrbSumFockR(nOrbBasis*2,nOrbsTot*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        OneIndOrbSumFockR=0.D0
        do p=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do a_prime=(nOrbBasis*2)+1,nOrbsTot*2  !Only contracted over the CABS space.
                            OneIndOrbSumFockR(p,q,x,y)=OneIndOrbSumFockR(p,q,x,y)+fock(p,a_prime)*antisym_f(a_prime,q,x,y)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        !Now can use this to find TERM 3: -r_uw^pq f_p^a' r_a'q^xy
        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do p=1,nOrbBasis*2
                            do q=1,nOrbBasis*2
                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(u,w,p,q)*OneIndOrbSumFockR(p,q,x,y)
                                if(tFindTerms) then
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)-antisym_f(u,w,p,q)*OneIndOrbSumFockR(p,q,x,y)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if(tFindTerms) then
            call CalcBTermEnergy(TermWise,EBSS,EBOS,3)
            EBSSSum=EBSSSum+EBSS
            EBOSSum=EBOSSum+EBOS
            TermWise=0.D0
        endif
        write(6,*) "Summed in TERM 3"

        !-r_uw^pa' 1RDM_p^q f_q^p' r_p'a'^xy
        !Unfortunately, we have to recalculate OneIndOrbSumFockR, since it is now contracted over the
        !formally complete space - not the CABS space.
        OneIndOrbSumFockR=0.D0
        do p=1,nOrbBasis*2
            do q_prime=(nOrbBasis*2)+1,nOrbsTot*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do p_prime=1,nOrbsTot*2  !Now contract over the formally complete space.
                            OneIndOrbSumFockR(p,q_prime,x,y)=OneIndOrbSumFockR(p,q_prime,x,y)+  &
                                fock(p,p_prime)*antisym_f(p_prime,q_prime,x,y)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        !Allocate scratch space for 1RDM contracted bit
        allocate(scratch(1:nOrbBasis*2,(nOrbBasis*2)+1:nOrbsTot*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        scratch=0.D0
        do p=1,nOrbBasis*2
            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do q=1,nOrbBasis*2
                            scratch(p,a_prime,x,y)=scratch(p,a_prime,x,y)+OneRDM(p,q)*OneIndOrbSumFockR(q,a_prime,x,y)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        !Find TERM7:
        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do p=1,nOrbBasis*2
                            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
!                                B(u,w,x,y)=B(u,w,x,y)-(antisym_f(u,w,p,a_prime)*scratch(p,a_prime,x,y))
!Use perm sym:
                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(p,a_prime,u,w)*scratch(p,a_prime,x,y)
                                if(tFindTerms) then
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)-antisym_f(p,a_prime,u,w)*scratch(p,a_prime,x,y)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        deallocate(scratch)
        deallocate(OneIndOrbSumFockR)
        if(tFindTerms) then
            call CalcBTermEnergy(TermWise,EBSS,EBOS,7)
            EBSSSum=EBSSSum+EBSS
            EBOSSum=EBOSSum+EBOS
            TermWise=0.D0
        endif
        write(6,*) "Summed in TERM 7"
        
        !Now calculate TERM 5
        allocate(OneIndOrbSumFockR((nOrbBasis*2)+1:nOrbsTot*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        OneIndOrbSumFockR=0.D0
        do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
            do q=1,nOrbBasis*2
                do y=1,nOrbBasis*2
                    do x=1,nOrbBasis*2
                        do b_prime=(nOrbBasis*2)+1,nOrbsTot*2  !Only contracted over the CABS space.
                            OneIndOrbSumFockR(a_prime,q,y,x)=OneIndOrbSumFockR(a_prime,q,y,x)+  &
                                fock(a_prime,b_prime)*antisym_f(b_prime,q,y,x)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        !Now contract with 1RDM
        allocate(scratch(1:nOrbBasis*2,(nOrbBasis*2)+1:nOrbsTot*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        scratch=0.D0
        do p=1,nOrbBasis*2
            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do q=1,nOrbBasis*2
                            scratch(p,a_prime,x,y)=scratch(p,a_prime,x,y)+OneRDM(p,q)*OneIndOrbSumFockR(a_prime,q,y,x)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        !Now add the TERM 5 contribution to B
        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do p=1,nOrbBasis*2
                            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
!                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(u,w,p,a_prime)*scratch(p,a_prime,x,y)
!Use perm sym:
                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(p,a_prime,u,w)*scratch(p,a_prime,x,y)
                                if(tFindTerms) then
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)-antisym_f(p,a_prime,u,w)*scratch(p,a_prime,x,y)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        deallocate(scratch,OneIndOrbSumFockR)
        if(tFindTerms) then
            call CalcBTermEnergy(TermWise,EBSS,EBOS,5)
            EBSSSum=EBSSSum+EBSS
            EBOSSum=EBOSSum+EBOS
            TermWise=0.D0
        endif
        write(6,*) "Summed in TERM 5"
        
!!I believe that term 5 is incorrect - recalculate it in a different way...
!        allocate(OneIndOrbSumFockR(nOrbBasis*2,(nOrbBasis*2+1):nOrbsTot*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
!        if(ierr.ne.0) stop 'error allocating'
!        OneIndOrbSumFockR=0.D0
!        do q=1,nOrbBasis*2
!            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
!                do x=1,nOrbBasis*2
!                    do y=1,nOrbBasis*2
!                        do b_prime=(nOrbBasis*2)+1,nOrbsTot*2  !Only contracted over the CABS space.
!                            OneIndOrbSumFockR(q,a_prime,x,y)=OneIndOrbSumFockR(q,a_prime,x,y)+  &
!                                fock(a_prime,b_prime)*antisym_f(q,b_prime,x,y)
!                        enddo
!                    enddo
!                enddo
!            enddo
!        enddo
!
!        !Now contract with 1RDM
!        allocate(scratch(1:nOrbBasis*2,(nOrbBasis*2)+1:nOrbsTot*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
!        if(ierr.ne.0) stop 'error allocating'
!        scratch=0.D0
!        do p=1,nOrbBasis*2
!            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
!                do x=1,nOrbBasis*2
!                    do y=1,nOrbBasis*2
!                        do q=1,nOrbBasis*2
!                            scratch(p,a_prime,x,y)=scratch(p,a_prime,x,y)+OneRDM(p,q)*OneIndOrbSumFockR(q,a_prime,x,y)
!                        enddo
!                    enddo
!                enddo
!            enddo
!        enddo
!
!        !Now add the TERM 5 contribution to B
!        do u=1,nOrbBasis*2
!            do w=1,nOrbBasis*2
!                do x=1,nOrbBasis*2
!                    do y=1,nOrbBasis*2
!                        do p=1,nOrbBasis*2
!                            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
!!                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(u,w,p,a_prime)*scratch(p,a_prime,x,y)
!!Use perm sym:
!!                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(p,a_prime,u,w)*scratch(p,a_prime,x,y)
!                                if(tFindTerms) then
!                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)-antisym_f(p,a_prime,u,w)*scratch(p,a_prime,x,y)
!                                endif
!                            enddo
!                        enddo
!                    enddo
!                enddo
!            enddo
!        enddo
!        deallocate(scratch,OneIndOrbSumFockR)
!        if(tFindTerms) then
!            call CalcBTermEnergy(TermWise,EBSS,EBOS,99)
!!            EBSSSum=EBSSSum+EBSS
!!            EBOSSum=EBOSSum+EBOS
!            TermWise=0.D0
!        endif
!        write(6,*) "Summed in TERM 5*"

        !find: f_p'^p 1RDM_p^q r_qa'^xy for terms 4 & 6
        allocate(scratch(1:nOrbsTot*2,(nOrbBasis*2)+1:nOrbsTot*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        allocate(scratch2(nOrbBasis*2,(nOrbBasis*2)+1:nOrbsTot*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        scratch=0.D0
        scratch2=0.D0
        !Split calculation of scratch up into two contractions...

        do p=1,nOrbBasis*2
            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do q=1,nOrbBasis*2
                            scratch2(p,a_prime,x,y)=scratch2(p,a_prime,x,y)+(OneRDM(p,q)*antisym_f(q,a_prime,x,y))
                        enddo
                    enddo
                enddo
            enddo
        enddo

        do p_prime=1,nOrbsTot*2
            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do p=1,nOrbBasis*2
                            scratch(p_prime,a_prime,x,y)=scratch(p_prime,a_prime,x,y)+(fock(p_prime,p)*scratch2(p,a_prime,x,y))
                        enddo
                    enddo
                enddo
            enddo
        enddo

        deallocate(scratch2)

!        do p_prime=1,nOrbsTot*2
!            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
!                do x=1,nOrbBasis*2
!                    do y=1,nOrbBasis*2
!                        do q=1,nOrbBasis*2
!                            do p=1,nOrbBasis*2
!                                scratch(p_prime,a_prime,x,y)=scratch(p_prime,a_prime,x,y)+(fock(p_prime,p)*OneRDM(p,q)*antisym_f(q,a_prime,x,y))
!                            enddo
!                        enddo
!                    enddo
!                enddo
!            enddo
!        enddo

        !Now just contract with integrals for TERM 6:
        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do p_prime=1,nOrbsTot*2
                            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
!                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(u,w,p_prime,a_prime)*scratch(p_prime,a_prime,x,y)
!Use perm sym:
                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(p_prime,a_prime,u,w)*scratch(p_prime,a_prime,x,y)
                                if(tFindTerms) then
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)-antisym_f(p_prime,a_prime,u,w)*scratch(p_prime,a_prime,x,y)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if(tFindTerms) then
            call CalcBTermEnergy(TermWise,EBSS,EBOS,6)
            EBSSSum=EBSSSum+EBSS
            EBOSSum=EBOSSum+EBOS
            TermWise=0.D0
        endif
        write(6,*) "Summed in TERM 6"

        !Contract with 1RDM and integrals for TERM 4
        allocate(scratch2(nOrbBasis*2,(nOrbBasis*2)+1:nOrbsTot*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        scratch2=0.D0
        do p=1,nOrbBasis*2
            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do q=1,nOrbBasis*2
                            scratch2(p,a_prime,x,y)=scratch2(p,a_prime,x,y)+OneRDM(p,q)*scratch(q,a_prime,x,y)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        deallocate(scratch)

        !Final contraction for TERM 4
        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do p=1,nOrbBasis*2
                            do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
!                                B(u,w,x,y)=B(u,w,x,y)+antisym_f(u,w,p,a_prime)*scratch2(p,a_prime,x,y)
!Use perm sym:
                                B(u,w,x,y)=B(u,w,x,y)+antisym_f(p,a_prime,u,w)*scratch2(p,a_prime,x,y)
                                if(tFindTerms) then
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)+antisym_f(p,a_prime,u,w)*scratch2(p,a_prime,x,y)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        deallocate(scratch2)
        if(tFindTerms) then
            call CalcBTermEnergy(TermWise,EBSS,EBOS,4)
            EBSSSum=EBSSSum+EBSS
            EBOSSum=EBOSSum+EBOS
            TermWise=0.D0
        endif
        write(6,*) "Summed in TERM 4"


        !*********************************************************************************************
        if(tCommutatorTrick) then

            write(6,*) "Using commutator trick..."

            !calculate tau term - 1/2 FTF_vw^xy . This has geminal indices in both the bra and the ket.
            !Therefore, we need to operate on the integrals by (3/8 + 1/8P_xy)(3/8 + 1/8P_vw.
            !Since there is only one term, do this on-the-fly with the spatial orbital integrals, rather
            !than explicitly transforming it into the spin-orbital representation first.
            !First, read these integrals in. They are only needed once, so deallocate after.
!            write(6,*) "Getting FTF_vw^xy tensor from MO_FTF file..."
            unit_ftf=get_free_unit()
            if(tformattedints) then
                open(unit_ftf,file='FTFDUMP',status='old',form='formatted')
            else
                open(unit_ftf,file='MO_FTF',status='old',form='unformatted',access='sequential',action='read')
            endif

            !Fill with spatial orbital integrals 
            allocate(ftfints(1:nOrbBasis,1:nOrbBasis,1:nOrbBasis,1:nOrbBasis),stat=ierr)
            if(ierr.ne.0) stop 'error allocating'

            !Fill the ftf array (in physical notation). All indices just over orbital basis
            Limits(1,1)=1
            Limits(2,1)=nOrbBasis
            Limits(1,2)=1
            Limits(2,2)=nOrbBasis
            Limits(1,3)=1
            Limits(2,3)=nOrbBasis
            Limits(1,4)=1
            Limits(2,4)=nOrbBasis

            call fill4indarr(unit_ftf,ftfints,Limits)

            !TERM 8 (Multiply contribution by 4, since it seems as though the DALTON integrals need to be.
            do u=1,nOrbBasis*2
                do w=1,nOrbBasis*2
                    do x=1,nOrbBasis*2
                        do y=1,nOrbBasis*2
                            SpinType=DetermineSpinType(u,w,x,y)
                            if(SpinType.eq.1) then
                                !<aa|aa>
!                                if(abs(ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(x),Conv2Spat(y))).gt.1.D-8) then
!                                write(6,*) Conv2Spat(u),Conv2Spat(w),Conv2Spat(x),Conv2Spat(y),ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(x),Conv2Spat(y))
!                                endif
                                B(u,w,x,y)=B(u,w,x,y)+(ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(x),Conv2Spat(y))- &
                                                       ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(y),Conv2Spat(x)))/16.D0
                            elseif(SpinType.eq.2) then
                                !<ab|ab>
                                B(u,w,x,y)=B(u,w,x,y)+(ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(x),Conv2Spat(y))) &
                                                            *(5.D0/32.D0) + &
                                                       (ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(y),Conv2Spat(x))) &
                                                            *(3.D0/32.D0)
                            elseif(SpinType.eq.3) then
                                !<ab|ba>
                                B(u,w,x,y)=B(u,w,x,y)-(ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(y),Conv2Spat(x))) &
                                                            *(5.D0/32.D0) - &
                                                       (ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(x),Conv2Spat(y))) &
                                                            *(3.D0/32.D0)
                            endif
                            if(tFindTerms) then
                                if(SpinType.eq.1) then
                                    !<aa|aa>
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)+    &
                                    (ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(x),Conv2Spat(y))- &
                                       ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(y),Conv2Spat(x)))/16.D0
                                elseif(SpinType.eq.2) then
                                    !<ab|ab>
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)+    &
                                    (ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(x),Conv2Spat(y))) &
                                                     *(5.D0/32.D0) + &
                                           (ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(y),Conv2Spat(x))) &
                                                              *(3.D0/32.D0)
                                elseif(SpinType.eq.3) then
                                    !<ab|ba>
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)-    &
                                        (ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(y),Conv2Spat(x))) &
                                                                *(5.D0/32.D0) - &
                                                           (ftfints(Conv2Spat(u),Conv2Spat(w),Conv2Spat(x),Conv2Spat(y))) &
                                                                *(3.D0/32.D0)
                                endif
!                                if(abs(TermWise(u,w,x,y)).gt.1.D-7) then
!                                    write(6,*) u,w,x,y,TermWise(u,w,x,y)
!                                endif
                            endif

                        enddo
                    enddo
                enddo
            enddo
            deallocate(ftfints)
            if(tFindTerms) then
                call CalcBTermEnergy(TermWise,EBSS,EBOS,8)
                EBSSSum=EBSSSum+EBSS
                EBOSSum=EBOSSum+EBOS
                TermWise=0.D0
            endif
            write(6,*) "Summed in TERM 8"

            if((.not.tReadRDMs).and.(.not.tformattedints)) then

                !Fill with spin orbital integrals 
                allocate(fpk(1:nOrbsTot*2,1:nOrbsTot*2),stat=ierr)
                if(ierr.ne.0) stop 'error allocating'

                !Now get the f+k matrix
                unit_fpk=get_free_unit()
                open(unit_fpk,file='MO_FpK',status='old',form='unformatted',access='sequential',action='read')

                !Fill the fg array (in physical notation). All indices just over orbital basis
                TwoIndLims(1)=1
                TwoIndLims(2)=nOrbsTot*2
                write(6,*) "Reading fock+exchange matrix from disk."
                call fill2indarr_spin(unit_fpk,fpk,TwoIndLims)
!                write(6,*) "Fock+K matrix read in"
            else
                !Need to calculate *generalised* f+k matrix
                !This should already have been generated earlier.
!                call GetfpkMat()
            endif

            !Now we also need to read in the r^2 matrix from MO_FF. Since both set of indices here are over
            !the geminal functions, we need to put in the cusp conditions to both of them.
!            write(6,*) "Getting rr_vw^xy tensor from MO_FG file (both geminal functions)..."
            unit_ff=get_free_unit()
            if(tformattedints) then
                open(unit_ff,file='FFDUMP',status='old',form='formatted')
            else
                open(unit_ff,file='MO_FF',status='old',form='unformatted',access='sequential',action='read')
            endif

            !Fill with spatial orbital integrals 
            allocate(ffints(1:nOrbsTot,1:nOrbBasis,1:nOrbBasis,1:nOrbBasis),stat=ierr)
            if(ierr.ne.0) stop 'error allocating'

            !Fill the fg array (in physical notation). All indices just over orbital basis
            Limits(1,1)=1
            Limits(2,1)=nOrbsTot
            Limits(1,2)=1
            Limits(2,2)=nOrbBasis
            Limits(1,3)=1
            Limits(2,3)=nOrbBasis
            Limits(1,4)=1
            Limits(2,4)=nOrbBasis

            call fill4indarr(unit_ff,ffints,Limits)

!            write(6,'(A)') "Antisymmetrizing ff integrals - converting to spin-orbital representation, "
!            write(6,'(A)') "now with both indices geminal..."
            allocate(antisym_ff(1:nOrbsTot*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
            if(ierr.ne.0) stop 'error allocating'
            
            !Fill the ff array (in physical notation). All indices just over *geminal* basis
            Limits(1,1)=1
            Limits(2,1)=nOrbsTot*2
            Limits(1,2)=1
            Limits(2,2)=nOrbBasis*2
            Limits(1,3)=1
            Limits(2,3)=nOrbBasis*2
            Limits(1,4)=1
            Limits(2,4)=nOrbBasis*2

    !        call antisymints(fgints,antisym_fg,Limits)
            call antisym2cuspints(ffints,antisym_ff,Limits)
            deallocate(ffints)

!        allocate(Test(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
!        Test=0.D0
!        do u=1,nOrbBasis*2
!            do w=1,nOrbBasis*2
!                do x=1,nOrbBasis*2
!                    do y=1,nOrbBasis*2
!                        do p_prime=1,nOrbsTot*2
!                            Test(u,w,x,y) = Test(u,w,x,y) + fpk(u,p_prime)*antisym_ff(p_prime,w,x,y)
!                        enddo
!                    enddo
!                enddo
!            enddo
!        enddo
!
!            !BooBug
!        open(24,file="FR-Model",status='unknown')
!
!        do x=1,nOrbBasis*2,2
!            do y=2,nOrbBasis*2,2
!                Val=0.D0
!                    do w=2,nOrbBasis*2,2
!                do u=1,nOrbBasis*2,2
!                        Val = Val - Test(u,w,x,y)*TwoRDM(x,y,w,u)/2.D0
!                        write(24,"(2G20.10,4I5)") Test(u,w,x,y),TwoRDM(x,y,w,u),(x+1)/2,y/2,w/2,(u+1)/2
!                    enddo
!                enddo
!!                write(24,"(G20.10,2I5)") Val,(x+1)/2,(y)/2
!            enddo
!        enddo



!                do w=2,nOrbBasis*2,2
!                    do u=1,nOrbBasis*2,2
!!                        write(24,"(G25.10,2I5)") Test2(u,w,q,a_prime),(u+1)/2,(w+1)/2
!!                        write(24,"(G25.10,2I5)") antisym_f(p_prime,a_prime,u,w),(u+1)/2,(w+1)/2
!!                        write(24,"(G25.10,2I5)") mo_k(u,w),(u+1)/2,(w+1)/2
!                        write(24,"(g25.10,4I5)") Test(u,w,x,y),(u+1)/2,(w+1)/2,(x+1)/2,(y+1)/2
!                    enddo
!                enddo
!            enddo
!        enddo



            if(tFindTerms) then
                allocate(TermWise3(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
                if(ierr.ne.0) stop 'error allocating'
                allocate(TermWise4(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
                if(ierr.ne.0) stop 'error allocating'
                TermWise3=0.D0
                TermWise4=0.D0
            endif

            do u=1,nOrbBasis*2
                do w=1,nOrbBasis*2
                    do x=1,nOrbBasis*2
                        do y=1,nOrbBasis*2
                            do p_prime=1,nOrbsTot*2
    !                            B(u,w,x,y)=B(u,w,x,y)+((fpk(u,p_prime)*antisym_ff(p_prime,w,x,y))/2.D0)+((antisym_ff(u,w,p_prime,y)*fpk(p_prime,x))/2.D0)
    !Perm sym: Term 9 + Term 11 + Term 10 + Term 12
                                B(u,w,x,y)=B(u,w,x,y)+((fpk(u,p_prime)*antisym_ff(p_prime,w,x,y))/2.D0)+    &
                                ((antisym_ff(p_prime,y,u,w)*fpk(p_prime,x))/2.D0)+  &
                                ((antisym_ff(p_prime,u,y,x)*fpk(w,p_prime))/2.D0)+  &
                                ((antisym_ff(p_prime,x,w,u)*fpk(p_prime,y))/2.D0)
                                if(tFindTerms) then
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)+((fpk(u,p_prime)*antisym_ff(p_prime,w,x,y))/2.D0)    !Term 9
                                    TermWise2(u,w,x,y)=TermWise2(u,w,x,y)+((antisym_ff(p_prime,y,u,w)*fpk(p_prime,x))/2.D0)  !Term 11
                                    TermWise3(u,w,x,y)=TermWise3(u,w,x,y)+((antisym_ff(p_prime,u,y,x)*fpk(w,p_prime))/2.D0)  !Term 10
                                    TermWise4(u,w,x,y)=TermWise4(u,w,x,y)+((antisym_ff(p_prime,x,w,u)*fpk(p_prime,y))/2.D0)  !Term 12
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            if(tFindTerms) then
                call CalcBTermEnergy(TermWise,EBSS,EBOS,9)
                EBSSSum=EBSSSum+EBSS
                EBOSSum=EBOSSum+EBOS
                call CalcBTermEnergy(TermWise3,EBSS,EBOS,10)
                EBSSSum=EBSSSum+EBSS
                EBOSSum=EBOSSum+EBOS
                call CalcBTermEnergy(TermWise2,EBSS,EBOS,11)
                EBSSSum=EBSSSum+EBSS
                EBOSSum=EBOSSum+EBOS
                call CalcBTermEnergy(TermWise4,EBSS,EBOS,12)
                EBSSSum=EBSSSum+EBSS
                EBOSSum=EBOSSum+EBOS
                TermWise=0.D0
                TermWise2=0.D0
                TermWise3=0.D0
                TermWise4=0.D0
            endif
            write(6,*) "SUMMED IN TERMS 9, 10, 11 and 12"

!New contractions to calculate terms 1-3 + 13 without double contraction over the whole space
        allocate(OneIndOrbSumFockR(nOrbsTot*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        OneIndOrbSumFockR=0.D0
        do p_prime=1,nOrbsTot*2     !Allow the transformed index to run over whole of space.
            do q=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do p=1,nOrbBasis*2  !Only contracted over the orbital space.
                            OneIndOrbSumFockR(p_prime,q,x,y)=OneIndOrbSumFockR(p_prime,q,x,y)+fpk(p_prime,p)*antisym_f(p,q,x,y)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        !Now calculate the term -r_uw^rq f_r^p r_pq^xy - put it straight into the B array
!        write(6,*) "Performing double-contraction to obtain B-term: -r_uw^rq f_r^p r_pq^xy..."
        !TERM 1
        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        !Sum over r & q
                        do r=1,nOrbBasis*2
                            do q=1,nOrbBasis*2
!                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(u,w,r,q)*OneIndOrbSumFockR(r,q,x,y)
                                if(tFindTerms) then
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)-antisym_f(u,w,r,q)*OneIndOrbSumFockR(r,q,x,y)
                                endif
                            enddo
                        enddo
                        !Seperately, calculate TERM 2!
                        !These can easily be combined later into one term.
                        do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
                            do q=1,nOrbBasis*2
!                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(u,w,a_prime,q)*OneIndOrbSumFockR(a_prime,q,x,y)
!Use perm sym:
!                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(a_prime,q,u,w)*OneIndOrbSumFockR(a_prime,q,x,y)
                                if(tFindTerms) then
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)-antisym_f(a_prime,q,u,w)*OneIndOrbSumFockR(a_prime,q,x,y)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        deallocate(OneIndOrbSumFockR)
        write(6,*) "Summed in TERMS 1* and 2*"



        !First find f_p^a' r_a'q^xy
!        write(6,*) "Calculating first-index fock-transformed F12 integrals, transformed over the CABS space."
        allocate(OneIndOrbSumFockR(nOrbBasis*2,nOrbsTot*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        OneIndOrbSumFockR=0.D0
        do p=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do a_prime=(nOrbBasis*2)+1,nOrbsTot*2  !Only contracted over the CABS space.
                            OneIndOrbSumFockR(p,q,x,y)=OneIndOrbSumFockR(p,q,x,y)+fpk(p,a_prime)*antisym_f(a_prime,q,x,y)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        !Now can use this to find TERM 3: -r_uw^pq f_p^a' r_a'q^xy
        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do p=1,nOrbBasis*2
                            do q=1,nOrbBasis*2
!                                B(u,w,x,y)=B(u,w,x,y)-antisym_f(u,w,p,q)*OneIndOrbSumFockR(p,q,x,y)
                                if(tFindTerms) then
                                    TermWise(u,w,x,y)=TermWise(u,w,x,y)-antisym_f(u,w,p,q)*OneIndOrbSumFockR(p,q,x,y)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        write(6,*) "Summed in TERM 3*"
            deallocate(fpk)


        allocate(scratch(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,(nOrbBasis*2)+1:nOrbsTot*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        scratch=0.D0

        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do r=1,nOrbBasis*2
                    do b_prime=(nOrbBasis*2)+1,nOrbsTot*2
                        do p_prime=1,nOrbsTot*2
                            scratch(u,w,r,b_prime)=scratch(u,w,r,b_prime)+mo_k(p_prime,r)*antisym_f(p_prime,b_prime,u,w)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do r=1,nOrbBasis*2
                            do b_prime=(nOrbBasis*2)+1,nOrbsTot*2
                                TermWise(u,w,x,y)=TermWise(u,w,x,y)-scratch(u,w,r,b_prime)*antisym_f(r,b_prime,x,y)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        deallocate(scratch)

        allocate(scratch(nOrbBasis*2,nOrbBasis*2,(nOrbBasis*2)+1:nOrbsTot*2,(nOrbBasis*2)+1:nOrbsTot*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        scratch=0.D0

        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do c_prime=(nOrbBasis*2)+1,nOrbsTot*2
                    do b_prime=(nOrbBasis*2)+1,nOrbsTot*2
                        do p_prime=1,nOrbsTot*2
                            scratch(u,w,c_prime,b_prime)=scratch(u,w,c_prime,b_prime)+mo_k(p_prime,c_prime)*antisym_f(p_prime,b_prime,u,w)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do c_prime=(nOrbBasis*2)+1,nOrbsTot*2
                            do b_prime=(nOrbBasis*2)+1,nOrbsTot*2
                                TermWise(u,w,x,y)=TermWise(u,w,x,y)-scratch(u,w,c_prime,b_prime)*antisym_f(c_prime,b_prime,x,y)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        deallocate(scratch)

        allocate(scratch(nOrbBasis*2,nOrbBasis*2,(nOrbBasis*2)+1:nOrbsTot*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        scratch=0.D0

        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do c_prime=(nOrbBasis*2)+1,nOrbsTot*2
                    do q=1,nOrbBasis*2
                        do a_prime=(nOrbBasis*2)+1,nOrbsTot*2
                            scratch(u,w,c_prime,q)=scratch(u,w,c_prime,q)+mo_k(a_prime,c_prime)*antisym_f(a_prime,q,u,w)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        do u=1,nOrbBasis*2
            do w=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do c_prime=(nOrbBasis*2)+1,nOrbsTot*2
                            do q=1,nOrbBasis*2
                                TermWise(u,w,x,y)=TermWise(u,w,x,y)-scratch(u,w,c_prime,q)*antisym_f(c_prime,q,x,y)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        deallocate(scratch)

        if(tFindTerms) then
            call CalcBTermEnergy(TermWise,EBSS,EBOS,14)
!            EBSSSum=EBSSSum+EBSS
!            EBOSSum=EBOSSum+EBOS
            TermWise=0.D0
        endif
        write(6,*) "Summed in TERMS 13*"












            allocate(scratch(1:nOrbsTot*2,1:nOrbsTot*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
            if(ierr.ne.0) stop 'error allocating'
            scratch=0.D0

            !Find k_p'^r' r_r'q'^xy = X_p'q'^xy
            do p_prime=1,nOrbsTot*2
                do q_prime=1,nOrbsTot*2
                    do x=1,nOrbBasis*2
                        do y=1,nOrbBasis*2
                            do r_prime=1,nOrbsTot*2
                                scratch(p_prime,q_prime,x,y)=scratch(p_prime,q_prime,x,y)+mo_k(p_prime,r_prime)*antisym_f(r_prime,q_prime,x,y)
                            enddo
                        enddo
                    enddo
                enddo
            enddo

!            write(6,*) "Double contraction over whole space coming up - this will be slow...!"

            !Now double contraction over whole space of scratch with r_uw^p'q'  THIS WILL BE REALLY REALLY SLOW!!!
            do u=1,nOrbBasis*2
                do w=1,nOrbBasis*2
                    do x=1,nOrbBasis*2
                        do y=1,nOrbBasis*2
                            do p_prime=1,nOrbsTot*2
                                do q_prime=1,nOrbsTot*2
                                    B(u,w,x,y)=B(u,w,x,y)-antisym_f(p_prime,q_prime,u,w)*scratch(p_prime,q_prime,x,y)
                                    if(tFindTerms) then
                                        TermWise(u,w,x,y)=TermWise(u,w,x,y)-antisym_f(p_prime,q_prime,u,w)*scratch(p_prime,q_prime,x,y)
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            deallocate(scratch,mo_k)
            if(tFindTerms) then
                call CalcBTermEnergy(TermWise,EBSS,EBOS,13)
                EBSSSum=EBSSSum+EBSS
                EBOSSum=EBOSSum+EBOS
                TermWise=0.D0
                deallocate(TermWise,TermWise2,TermWise3,TermWise4)
            endif
            write(6,*) "SUMMED IN TERM 13"

        ELSE    !Don't use commutator trick - only one term, rather than the 4 terms above
                !However, this will have very slow CABS RI convergence.

            write(6,*) "NOT using commutator trick - this will have slow CABS RI convergence."

            !Still want to read in the r^2 matrix from MO_FF for later though. 
            !Since both set of indices here are over
            !the geminal functions, we need to put in the cusp conditions to both of them.
            write(6,*) "Getting rr_vw^xy tensor from MO_FG file (both geminal functions)..."
            unit_ff=get_free_unit()
            if(tformattedints) then
                open(unit_ff,file='FFDUMP',status='old',form='formatted')
            else
                open(unit_ff,file='MO_FF',status='old',form='unformatted',access='sequential',action='read')
            endif

            !Fill with spatial orbital integrals 
            allocate(ffints(1:nOrbsTot,1:nOrbBasis,1:nOrbBasis,1:nOrbBasis),stat=ierr)
            if(ierr.ne.0) stop 'error allocating'

            !Fill the fg array (in physical notation). All indices just over orbital basis
            Limits(1,1)=1
            Limits(2,1)=nOrbsTot
            Limits(1,2)=1
            Limits(2,2)=nOrbBasis
            Limits(1,3)=1
            Limits(2,3)=nOrbBasis
            Limits(1,4)=1
            Limits(2,4)=nOrbBasis

            call fill4indarr(unit_ff,ffints,Limits)

            write(6,'(A)') "Antisymmetrizing ff integrals - converting to spin-orbital representation, "
            write(6,'(A)') "now with both indices geminal..."
            allocate(antisym_ff(1:nOrbsTot*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
            if(ierr.ne.0) stop 'error allocating'
            
            !Fill the ff array (in physical notation). All indices just over *geminal* basis
            Limits(1,1)=1
            Limits(2,1)=nOrbsTot*2
            Limits(1,2)=1
            Limits(2,2)=nOrbBasis*2
            Limits(1,3)=1
            Limits(2,3)=nOrbBasis*2
            Limits(1,4)=1
            Limits(2,4)=nOrbBasis*2

    !        call antisymints(fgints,antisym_fg,Limits)
            call antisym2cuspints(ffints,antisym_ff,Limits)
            deallocate(ffints)

            !Now to obtain: r_uw^p'q' f_p'^r' r_r'q'^xy
            allocate(scratch(1:nOrbstot*2,1:nOrbsTot*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
            if(ierr.ne.0) stop 'error allocating'
            scratch=0.D0

            do p_prime=1,nOrbsTot*2
                do q_prime=1,nOrbsTot*2
                    do x=1,nOrbBasis*2
                        do y=1,nOrbBasis*2
                            do r_prime=1,nOrbsTot*2
                                scratch(p_prime,q_prime,x,y)=scratch(p_prime,q_prime,x,y)+fock(p_prime,r_prime)*antisym_f(r_prime,q_prime,x,y)
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            !Now contract with r

            do u=1,nOrbBasis*2
                do w=1,nOrbBasis*2
                    do x=1,nOrbBasis*2
                        do y=1,nOrbBasis*2
                            do p_prime=1,nOrbsTot*2
                                do q_prime=1,nOrbsTot*2
                                    B(u,w,x,y)=B(u,w,x,y)+antisym_f(p_prime,q_prime,u,w)*scratch(p_prime,q_prime,x,y)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            deallocate(scratch)

        ENDIF

    end subroutine CalculateB
    
    
    subroutine antisym2cuspints(ints,antisym_ints,Lims)
        !** The ints are in spatial coordinates and antisym_ints returned in spin coords. **
        ! This routine assumes that BOTH sets of indices are geminal indices.
        !Here, assume that <vw|f|xy> is is a spinless orbital integral (i.e. spin integrated out)
        !Therefore <vw||xy> for the different spins is:
        !<aa||aa> = 1/16[ <vw|xy> - <vw|yw> ]
        !<ab||ab> = 5/32 <vw|xy> + 3/32 <vw|yx>
        !<ab||ba> = -5/32 <vw|yx> - 3/32 <vw|xy>
        !I believe that if the geminal indices are in the bra, then the same will hold.
        integer, intent(in) :: Lims(2,4)
        real(8), intent(in) :: ints(Lims(1,1):Lims(2,1)/2,Lims(1,2):Lims(2,2)/2,Lims(1,3):Lims(2,3)/2,Lims(1,4):Lims(2,4)/2)
        real(8), intent(out) :: antisym_ints(Lims(1,1):Lims(2,1),Lims(1,2):Lims(2,2),Lims(1,3):Lims(2,3),Lims(1,4):Lims(2,4))
        integer :: i,j,k,l,SpinType

        antisym_ints=0.D0

!        write(6,*) "Antisymmetrizing integrals into spin-orbital notation: "
!        write(6,*) "i: ",Lims(1,1)," -> ",Lims(2,1)
!        write(6,*) "j: ",Lims(1,2)," -> ",Lims(2,2)
!        write(6,*) "k: ",Lims(1,3)," -> ",Lims(2,3)
!        write(6,*) "l: ",Lims(1,4)," -> ",Lims(2,4)

        if(Lims(2,3).ne.Lims(2,4)) stop 'Limits wrong for antisymmetrization'

        !do <ij||kl> = <ij|kl> - <ij|lk>
        do i=Lims(1,1),Lims(2,1)
            do j=Lims(1,2),Lims(2,2)
                do k=Lims(1,3),Lims(2,3)
                    do l=Lims(1,4),Lims(2,4)
                        SpinType=DetermineSpinType(i,j,k,l)
!                        write(6,*) i,j,k,l,SpinType

                        if(spinType.eq.1) then
                            !<aa||aa> = 1/16[ <vw|xy> - <vw|yw> ]
                            antisym_ints(i,j,k,l) = (ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(k),Conv2Spat(l))- &
                                                     ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(l),Conv2Spat(k)))/16.D0
                        elseif(spinType.eq.2) then
                            !<ab||ab> = 5/32 <vw|xy> + 3/32 <vw|yx>
                           antisym_ints(i,j,k,l)=(5.D0/32.D0)*ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(k),Conv2Spat(l))+ &
                                              (3.D0/32.D0)*ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(l),Conv2Spat(k))
                        elseif(spinType.eq.3) then
                            !<ab||ba> = -5/32 <vw|yx> - 3/32 <vw|xy>
                          antisym_ints(i,j,k,l)=-(3.D0/32.D0)*ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(k),Conv2Spat(l))- &
                                              (5.D0/32.D0)*ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(l),Conv2Spat(k))
                        else
                            antisym_ints(i,j,k,l)=0.D0
                        endif
!                        antisym_ints(i,j,k,l)=ints(i,j,k,l)-ints(i,j,l,k)
                    enddo
                enddo
            enddo
        enddo

    end subroutine antisym2cuspints
    
    subroutine antisymcuspints(ints,antisym_ints,Lims)
        !** The ints are in spatial coordinates and antisym_ints returned in spin coords. **
        !Here, assume that <ij|f|xy> is is a spinless orbital integral (i.e. spin integrated out)
        !Therefore <ij||xy> for the different spins is:
        !<aa||aa> = 1/4[ <ij|xy> - <ij|yx> ]
        !<ab||ab> = 3/8 <ij|xy> + 1/8 <ij|yx>
        !<ab||ba> = -3/8 <ij|yx> - 1/8 <ij|xy>
        !I believe that if the geminal indices are in the bra, then the same will hold.
        integer, intent(in) :: Lims(2,4)
        real(8), intent(in) :: ints((Lims(1,1)+1)/2:Lims(2,1)/2,(Lims(1,2)+1)/2:Lims(2,2)/2,(Lims(1,3)+1)/2:Lims(2,3)/2,(Lims(1,4)+1)/2:Lims(2,4)/2)
        real(8), intent(out) :: antisym_ints(Lims(1,1):Lims(2,1),Lims(1,2):Lims(2,2),Lims(1,3):Lims(2,3),Lims(1,4):Lims(2,4))
        integer :: i,j,k,l,SpinType

        antisym_ints=0.D0

!        write(6,*) "Antisymmetrizing integrals into spin-orbital notation: "
!        write(6,*) "i: ",Lims(1,1)," -> ",Lims(2,1)
!        write(6,*) "j: ",Lims(1,2)," -> ",Lims(2,2)
!        write(6,*) "k: ",Lims(1,3)," -> ",Lims(2,3)
!        write(6,*) "l: ",Lims(1,4)," -> ",Lims(2,4)

        if(Lims(2,3).ne.Lims(2,4)) stop 'Limits wrong for antisymmetrization'

        do i=Lims(1,1),Lims(2,1)
            do j=Lims(1,2),Lims(2,2)
                do k=Lims(1,3),Lims(2,3)
                    do l=Lims(1,4),Lims(2,4)
                        SpinType=DetermineSpinType(i,j,k,l)
!                        write(6,*) i,j,k,l,SpinType

                        if(spinType.eq.1) then
                            !<aa||aa> = 1/4 [ <ij|kl> - <ij|lk>
                            antisym_ints(i,j,k,l) = (ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(k),Conv2Spat(l))- &
                                                     ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(l),Conv2Spat(k)))/4.D0
                        elseif(spinType.eq.2) then
                            !<ab||ab> = 3/8 <ij|kl> + 1/8 <ij|lk>
                           antisym_ints(i,j,k,l)=(3.D0/8.D0)*ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(k),Conv2Spat(l))+ &
                                              (1.D0/8.D0)*ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(l),Conv2Spat(k))
                        elseif(spinType.eq.3) then
                            !<ab||ba> = -3/8 <ij|lk> - 1/8 <ij|kl>
                          antisym_ints(i,j,k,l)=-(1.D0/8.D0)*ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(k),Conv2Spat(l))- &
                                              (3.D0/8.D0)*ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(l),Conv2Spat(k))
                        else
                            antisym_ints(i,j,k,l)=0.D0
                        endif
!                        antisym_ints(i,j,k,l)=ints(i,j,k,l)-ints(i,j,l,k)
                    enddo
                enddo
            enddo
        enddo

    end subroutine antisymcuspints

    integer function Conv2Spat(i)
        integer , intent(in) :: i

        if(mod(i,2).eq.1) then
            !alpha orbital
            Conv2Spat=(i+1)/2
        else
            Conv2Spat=i/2
        endif
    end function Conv2Spat

    !Return 1 if the integral is of spin type <aa|aa>
    !Return 2 if the integral is of spin type <ab|ab>
    !Return 3 if the integral is of spin type <ab|ba>
    !Return 4 if the integral is of spin type <aa|bb> or <aa|ab> or <aa|ba> ...
    integer function DetermineSpinType(i,j,k,l)
        integer , intent(in) :: i,j,k,l
        logical :: ispin,jspin,kspin,lspin

        if(mod(i,2).eq.1) then
            Ispin=.true.        !Odd indices -> alpha orbital -> Logical TRUE
        else
            Ispin=.false.       !Even indices -> beta orbital -> Logical FALSE
        endif

        if(mod(j,2).eq.1) then
            Jspin=.true.
        else
            Jspin=.false.
        endif

        if(mod(k,2).eq.1) then
            Kspin=.true.
        else
            Kspin=.false.
        endif

        if(mod(l,2).eq.1) then
            Lspin=.true.
        else
            Lspin=.false.
        endif

        if((Ispin.eqv.Jspin).and.(Ispin.eqv.Kspin).and.(Ispin.eqv.Lspin)) then
            DetermineSpinType=1
        elseif((Ispin.eqv.Kspin).and.(Jspin.eqv.Lspin).and.(Ispin.neqv.Jspin)) then
            DetermineSpinType=2
        elseif((Ispin.eqv.Lspin).and.(Jspin.eqv.Kspin).and.(Ispin.neqv.Jspin)) then
            DetermineSpinType=3
        else
            DetermineSpinType=4
        endif
    end function DetermineSpinType

    subroutine antisymints(ints,antisym_ints,Lims)
        integer, intent(in) :: Lims(2,4)
        real(8), intent(in) :: ints(Lims(1,1):Lims(2,1)/2,Lims(1,2):Lims(2,2)/2,Lims(1,3):Lims(2,3)/2,Lims(1,4):Lims(2,4)/2)
        real(8), intent(out) :: antisym_ints(Lims(1,1):Lims(2,1),Lims(1,2):Lims(2,2),Lims(1,3):Lims(2,3),Lims(1,4):Lims(2,4))
        integer :: i,j,k,l,SpinType

        antisym_ints=0.D0

        write(6,*) "Antisymmetrizing integrals in spin-orbital notation: "
!        write(6,*) "i: ",Lims(1,1)," -> ",Lims(2,1)
!        write(6,*) "j: ",Lims(1,2)," -> ",Lims(2,2)
!        write(6,*) "k: ",Lims(1,3)," -> ",Lims(2,3)
!        write(6,*) "l: ",Lims(1,4)," -> ",Lims(2,4)

        if(Lims(2,3).ne.Lims(2,4)) then
            !since k and l have different limits, actually do: <ij||ab> = <ij|ab> - <ji|ab>

            do i=Lims(1,1),Lims(2,1)
                do j=Lims(1,2),Lims(2,2)
                    do k=Lims(1,3),Lims(2,3)
                        do l=Lims(1,4),Lims(2,4)
                            SpinType=DetermineSpinType(i,j,k,l)
                            if(spinType.eq.1) then
                                !<aa|aa>
                                antisym_ints(i,j,k,l)=ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(k),Conv2Spat(l)) &
                                    -ints(Conv2Spat(j),Conv2Spat(i),Conv2Spat(k),Conv2Spat(l))
                            elseif(SpinType.eq.2) then
                                !<ab||ab>
                                antisym_ints(i,j,k,l)=ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(k),Conv2Spat(l))
                            elseif(SpinType.eq.3) then
                                !<ab||ba>
                                antisym_ints(i,j,k,l)=-ints(Conv2Spat(j),Conv2Spat(i),Conv2Spat(k),Conv2Spat(l))
                            else
                                antisym_ints(i,j,k,l)=0.D0
                            endif
                        enddo
                    enddo
                enddo
            enddo

        else
            !do <ij||kl> = <ij|kl> - <ij|lk>

            do i=Lims(1,1),Lims(2,1)
                do j=Lims(1,2),Lims(2,2)
                    do k=Lims(1,3),Lims(2,3)
                        do l=Lims(1,4),Lims(2,4)
                            SpinType=DetermineSpinType(i,j,k,l)
                            if(spinType.eq.1) then
                                !<aa|aa>
                                antisym_ints(i,j,k,l)=ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(k),Conv2Spat(l)) &
                                    -ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(l),Conv2Spat(k))
                            elseif(SpinType.eq.2) then
                                !<ab||ab>
                                antisym_ints(i,j,k,l)=ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(k),Conv2Spat(l))
                            elseif(SpinType.eq.3) then
                                !<ab||ba>
                                antisym_ints(i,j,k,l)=-ints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(l),Conv2Spat(k))
                            else
                                antisym_ints(i,j,k,l)=0.D0
                            endif
!                            antisym_ints(i,j,k,l)=ints(i,j,k,l)-ints(i,j,l,k)
                        enddo
                    enddo
                enddo
            enddo
        endif

    end subroutine antisymints

    subroutine FindTwoRDM_spin()
        integer :: ierr,i,j,k,l,TwoRDM_read,ios
        real*8 :: TwoRDMVal
        logical :: exists

        allocate(TwoRDM(nOrbBasis*2,nOrbBasis*2,nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        TwoRDM(:,:,:,:) = 0.D0

        if(tReadRDMs) then

            write(6,*) "Reading in 2RDM from disk..."
            inquire(file='TwoRDM',exist=exists)
            if(.not.exists) then
                stop "Cannot find file 'TwoRDM' to read in One-body density matrix"
            endif

            TwoRDM_read=get_free_unit()
            open(TwoRDM_read,file='TwoRDM',status='old',action='read')
            rewind(TwoRDM_read)

            do while(.true.)
                read(TwoRDM_read,"(4I6,G25.17)",iostat=ios) i,j,k,l,TwoRDMVal
                if(ios.gt.0) stop "Error reading TwoRDM"
                if(ios.lt.0) exit
                !Read in in spin-orbital notation.

                if((i.gt.nOrbBasis*2).or.(j.gt.nOrbBasis*2).or.(k.gt.nOrbBasis*2).or.(l.gt.nOrbBasis*2)) then
                    stop 'Bounds incorrect for 2RDM'
                endif
    
                !Symmetric permutations
                if((abs(TwoRDM(i,j,k,l)-TwoRDMVal).gt.1.D-7).and.(abs(TwoRDM(i,j,k,l)).gt.1.D-8)) then
                    write(6,*) i,j,k,l,TwoRDM(i,j,k,l),TwoRDMVal
                    stop 'Error in filling 2RDM 1'
                else
!                    write(6,*) TwoRDM(i,j,k,l),TwoRDMVal,i,j,k,l
                    TwoRDM(i,j,k,l)=TwoRDMVal
                endif
                if((abs(TwoRDM(j,i,l,k)-TwoRDMVal).gt.1.D-7).and.(abs(TwoRDM(j,i,l,k)).gt.1.D-8)) then
!                    write(6,*) TwoRDM(j,i,l,k),TwoRDMVal,j,i,l,k
                    stop 'Error in filling 2RDM 2'
                else
!                    write(6,*) TwoRDM(j,i,l,k),TwoRDMVal,j,i,l,k
                    TwoRDM(j,i,l,k)=TwoRDMVal
                endif
                if((abs(TwoRDM(k,l,i,j)-TwoRDMVal).gt.1.D-7).and.(abs(TwoRDM(k,l,i,j)).gt.1.D-8)) then
                    stop 'Error in filling 2RDM 3'
                else
!                    write(6,*) TwoRDM(k,l,i,j),TwoRDMVal,k,l,i,j
                    TwoRDM(k,l,i,j)=TwoRDMVal
                endif
                if((abs(TwoRDM(l,k,j,i)-TwoRDMVal).gt.1.D-7).and.(abs(TwoRDM(l,k,j,i)).gt.1.D-8)) then
                    stop 'Error in filling 2RDM 4'
                else
!                    write(6,*) TwoRDM(l,k,j,i),TwoRDMVal,l,k,j,i
                    TwoRDM(l,k,j,i)=TwoRDMVal
                endif

                !Antisymmetric permutations
                if((abs(TwoRDM(j,i,k,l)+TwoRDMVal).gt.1.D-7).and.(abs(TwoRDM(j,i,k,l)).gt.1.D-8)) then
                    stop 'Error in filling 2RDM 5'
                else
!                    write(6,*) TwoRDM(j,i,k,l),TwoRDMVal,j,i,k,l
                    TwoRDM(j,i,k,l)=-1.D0 * TwoRDMVal
                endif
                if((abs(TwoRDM(i,j,l,k)+TwoRDMVal).gt.1.D-7).and.(abs(TwoRDM(i,j,l,k)).gt.1.D-8)) then
                    stop 'Error in filling 2RDM 6'
                else
!                    write(6,*) TwoRDM(i,j,l,k),TwoRDMVal,i,j,l,k
                    TwoRDM(i,j,l,k)=-1.D0 * TwoRDMVal
                endif
                if((abs(TwoRDM(k,l,j,i)+TwoRDMVal).gt.1.D-7).and.(abs(TwoRDM(k,l,j,i)).gt.1.D-8)) then
                    stop 'Error in filling 2RDM 7'
                else
!                    write(6,*) TwoRDM(k,l,j,i),TwoRDMVal,k,l,j,i
                    TwoRDM(k,l,j,i)=-1.D0 * TwoRDMVal
                endif
                if((abs(TwoRDM(l,k,i,j)+TwoRDMVal).gt.1.D-7).and.(abs(TwoRDM(l,k,i,j)).gt.1.D-8)) then
                    stop 'Error in filling 2RDM 8'
                else
!                    write(6,*) TwoRDM(l,k,i,j),TwoRDMVal,l,k,i,j
                    TwoRDM(l,k,i,j)=-1.D0 * TwoRDMVal
                endif
            enddo
            close(TwoRDM_read)

        else
            write(6,*) "Choosing spin 2RDM from the HF determinant..."
            write(6,*) "Therefore, 2RDM_pq^rs = 1/2\delta_pr \delta_qs - 1/2\delta_ps \delta_qr..."

            do p=1,HFOccOrbs*2
                do q=1,HFOccOrbs*2
                    do r=1,HFOccOrbs*2
                        do s=1,HFOccOrbs*2
                            if((p.eq.r).and.(q.eq.s)) then
                                TwoRDM(p,q,r,s)=TwoRDM(p,q,r,s)+1.0
                            endif 
                            if((r.eq.q).and.(p.eq.s)) then
    !                            if(mod(p,2).eq.mod(q,2)) then
                                TwoRDM(p,q,r,s)=TwoRDM(p,q,r,s)-1.0
    !                            endif
                            endif
                        enddo
                    enddo
                enddo
            enddo
        endif
        
    end subroutine FindTwoRDM_spin
        
    !Routine to find 1RDM
    subroutine findOneRDM_spin()
        integer :: ierr,i,j,OneRDM_read,ios
        real*8 :: OneRDMVal
        logical :: exists
            
        allocate(OneRDM(nOrbBasis*2,nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        OneRDM=0.D0

        if(tReadRDMs) then

            inquire(file='OneRDM',exist=exists)
            if(.not.exists) then
                stop "Cannot find file 'OneRDM' to read in One-body density matrix"
            endif
            OneRDM_read=get_free_unit()
            open(OneRDM_read,file='OneRDM',status='old',action='read')
            rewind(OneRDM_read)

            do while(.true.)
                read(OneRDM_read,"(2I6,G25.17)",iostat=ios) i,j,OneRDMVal
                if(ios.gt.0) stop "Error reading OneRDM"
                if(ios.lt.0) exit
                !Read in in spin-orbital notation.
                OneRDM(i,j)=OneRDMVal
                OneRDM(j,i)=OneRDMVal
            enddo
            close(OneRDM_read)

        else

            write(6,*) "Constructing One-body density matrix..."

            write(6,*) "Choosing 1RDM to be single determinant: 1RDM_p^q = \delta_pq"
    !        write(6,*) "Therefore, 1RDM_p^q = \delta_pq..."

            do i=1,HFOccOrbs*2
                OneRDM(i,i)=1.D0    !Since we are working in spatial orbitals, shouldn't this be two?
            enddo

        endif

    end subroutine findOneRDM_spin
    
    !Fill a four-index array with integrals from a file.
    !The basis indicators, indicate the basis which the array wants to run over.
    subroutine fill4indarrspin(fu,arr,Lim)
        
        integer, intent(in) :: Lim(2,4)  !Limits(1,:) = minimum indices of array dimensions ; (2,:) = max.
        real(8), intent(out) :: arr(Lim(1,1):Lim(2,1),Lim(1,2):Lim(2,2),Lim(1,3):Lim(2,3),Lim(1,4):Lim(2,4))
        integer , intent(in) :: fu
        integer(8) :: MaxLength,length
        integer(2) , allocatable :: Indices(:,:)
        real(8) , allocatable :: Buf(:)
        integer :: iBuffer,i,j,k,l,ii,jj,kk,ll,Mini,Minj,Mink,Minl,Maxi,Maxj,Maxk,Maxl
        integer(2) :: MaxFilei,MaxFilej,MaxFilek,MaxFilel
        integer :: iSpin,jSpin,kSpin,lSpin,MaxiSpin,MaxjSpin,MaxkSpin,MaxlSpin,MiniSpin,MinjSpin,MinkSpin,MinlSpin
        
        CALL FindMaxIndices(fu,MaxFilei,MaxFilej,MaxFilek,MaxFilel)
!        write(6,"(A,4I5)") "Index upper limits for the file are (physical notation): ",MaxFilei,MaxFilej,MaxFilek,MaxFilel

        MaxFilei=MaxFilei*2 !Now storing spin orbitals but written on disk as spatial.
        MaxFilej=MaxFilej*2 !Now storing spin orbitals but written on disk as spatial.
        MaxFilek=MaxFilek*2 !Now storing spin orbitals but written on disk as spatial.
        MaxFilel=MaxFilel*2 !Now storing spin orbitals but written on disk as spatial.
        
        call check_limits(MaxFilei,MaxFilej,MaxFilek,MaxFilel,Lim)

        Mini=(Lim(1,1)+1)/2
        Maxi=Lim(2,1)/2
        Minj=(Lim(1,2)+1)/2
        Maxj=Lim(2,2)/2
        Mink=(Lim(1,3)+1)/2
        Maxk=Lim(2,3)/2
        Minl=(Lim(1,4)+1)/2
        Maxl=Lim(2,4)/2

        write(6,*) "In spatial orbital notation, limits of these integrals which will be stored are: "
        write(6,*) "i: ",Mini," -> ",Maxi
        write(6,*) "j: ",Minj," -> ",Maxj
        write(6,*) "k: ",Mink," -> ",Maxk
        write(6,*) "l: ",Minl," -> ",Maxl

        arr=0.D0

        rewind(fu)
        read(fu) maxlength
        allocate(Indices(4,MaxLength),Buf(MaxLength))

        iBuffer=0
        do while(.true.)

            iBuffer=iBuffer+1
!            WRITE(6,*) "Reading buffer: ",iBuffer

            read(fu,err=12) length,Indices(1:4,1:length),Buf(1:length)
            IF(length.le.0) goto 13
!            WRITE(6,"(I12,A)") length," integrals read. Writing Integrals..."

            do i=1,length
                !Its read in in *chemical notation*...immediately switch to physical
                ii=Indices(1,i)
                kk=Indices(2,i)
                jj=Indices(3,i)
                ll=Indices(4,i)
!                WRITE(luWrite,'(1X,G20.12,4I3)') Buf(i),Indices(1,i),Indices(2,i),Indices(3,i),Indices(4,i)

                IF((ii.le.Maxi).and.(ii.ge.Mini).and.(jj.le.Maxj).and.(jj.ge.Minj).and.(kk.le.Maxk) &
                    .and.(kk.ge.Minj).and.(ll.le.Maxl).and.(ll.ge.Minl)) THEN

                    do iSpin=ii*2-1,ii*2
                        do jSpin=jj*2-1,jj*2
                            do kSpin=kk*2-1,kk*2
                                do lSpin=ll*2-1,ll*2
                                    if(AllowedSpin(iSpin,jSpin,kSpin,lSpin)) then
                                        arr(iSpin,jSpin,kSpin,lSpin)=Buf(i)
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo

                ENDIF
                !Now for rest of permutations...
  
                !<ji|lk>
                IF((jj.le.Maxi).and.(jj.ge.Mini).and.(ii.le.Maxj).and.(ii.ge.Minj).and.(ll.le.Maxk) &
                    .and.(ll.ge.Minj).and.(kk.le.Maxl).and.(kk.ge.Minl)) THEN

                    do iSpin=jj*2-1,jj*2
                        do jSpin=ii*2-1,ii*2
                            do kSpin=ll*2-1,ll*2
                                do lSpin=kk*2-1,kk*2
                                    if(AllowedSpin(iSpin,jSpin,kSpin,lSpin)) then
                                        arr(iSpin,jSpin,kSpin,lSpin)=Buf(i)
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo

!                    arr(jj,ii,ll,kk)=Buf(i)
                endif

                !<kl|ij>
                IF((kk.le.Maxi).and.(kk.ge.Mini).and.(ll.le.Maxj).and.(ll.ge.Minj).and.(ii.le.Maxk) &
                    .and.(ii.ge.Minj).and.(jj.le.Maxl).and.(jj.ge.Minl)) THEN

                    do iSpin=kk*2-1,kk*2
                        do jSpin=ll*2-1,ll*2
                            do kSpin=ii*2-1,ii*2
                                do lSpin=jj*2-1,jj*2
                                    if(AllowedSpin(iSpin,jSpin,kSpin,lSpin)) then
                                        arr(iSpin,jSpin,kSpin,lSpin)=Buf(i)
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo

!                    arr(kk,ll,ii,jj)=Buf(i)
                endif

                !<lk|ji>
                IF((ll.le.Maxi).and.(ll.ge.Mini).and.(kk.le.Maxj).and.(kk.ge.Minj).and.(jj.le.Maxk) &
                    .and.(jj.ge.Minj).and.(ii.le.Maxl).and.(ii.ge.Minl)) THEN

                    do iSpin=ll*2-1,ll*2
                        do jSpin=kk*2-1,kk*2
                            do kSpin=jj*2-1,jj*2
                                do lSpin=ii*2-1,ii*2
                                    if(AllowedSpin(iSpin,jSpin,kSpin,lSpin)) then
                                        arr(iSpin,jSpin,kSpin,lSpin)=Buf(i)
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo

!                    arr(ll,kk,jj,ii)=Buf(i)
                endif

                !<kj|il>
                IF((kk.le.Maxi).and.(kk.ge.Mini).and.(jj.le.Maxj).and.(jj.ge.Minj).and.(ii.le.Maxk) &
                    .and.(ii.ge.Minj).and.(ll.le.Maxl).and.(ll.ge.Minl)) THEN

                    do iSpin=kk*2-1,kk*2
                        do jSpin=jj*2-1,jj*2
                            do kSpin=ii*2-1,ii*2
                                do lSpin=ll*2-1,ll*2
                                    if(AllowedSpin(iSpin,jSpin,kSpin,lSpin)) then
                                        arr(iSpin,jSpin,kSpin,lSpin)=Buf(i)
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo

!                    arr(kk,jj,ii,ll)=Buf(i)
                endif

                !<li|jk>
                IF((ll.le.Maxi).and.(ll.ge.Mini).and.(ii.le.Maxj).and.(ii.ge.Minj).and.(jj.le.Maxk) &
                    .and.(jj.ge.Minj).and.(kk.le.Maxl).and.(kk.ge.Minl)) THEN

                    do iSpin=ll*2-1,ll*2
                        do jSpin=ii*2-1,ii*2
                            do kSpin=jj*2-1,jj*2
                                do lSpin=kk*2-1,kk*2
                                    if(AllowedSpin(iSpin,jSpin,kSpin,lSpin)) then
                                        arr(iSpin,jSpin,kSpin,lSpin)=Buf(i)
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo

!                    arr(ll,ii,jj,kk)=Buf(i)
                endif

                !<il|kj>
                IF((ii.le.Maxi).and.(ii.ge.Mini).and.(ll.le.Maxj).and.(ll.ge.Minj).and.(kk.le.Maxk) &
                    .and.(kk.ge.Minj).and.(jj.le.Maxl).and.(jj.ge.Minl)) THEN

                    do iSpin=ii*2-1,ii*2
                        do jSpin=ll*2-1,ll*2
                            do kSpin=kk*2-1,kk*2
                                do lSpin=jj*2-1,jj*2
                                    if(AllowedSpin(iSpin,jSpin,kSpin,lSpin)) then
                                        arr(iSpin,jSpin,kSpin,lSpin)=Buf(i)
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo

!                    arr(ii,ll,kk,jj)=Buf(i)
                endif

                !<jk|li>
                IF((jj.le.Maxi).and.(jj.ge.Mini).and.(kk.le.Maxj).and.(kk.ge.Minj).and.(ll.le.Maxk) &
                    .and.(ll.ge.Minj).and.(ii.le.Maxl).and.(ii.ge.Minl)) THEN

                    do iSpin=jj*2-1,jj*2
                        do jSpin=kk*2-1,kk*2
                            do kSpin=ll*2-1,ll*2
                                do lSpin=ii*2-1,ii*2
                                    if(AllowedSpin(iSpin,jSpin,kSpin,lSpin)) then
                                        arr(iSpin,jSpin,kSpin,lSpin)=Buf(i)
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo

!                    arr(jj,kk,ll,ii)=Buf(i)
                endif
            enddo

!            WRITE(6,*) "Finished integral buffer."

        enddo

12      STOP 'Error reading block'

13      CONTINUE

!        WRITE(6,*) "All Integrals read"

        DEALLOCATE(Indices)
        DEALLOCATE(Buf)

        MiniSpin=Mini*2-1
        MinjSpin=Minj*2-1
        MinkSpin=Mink*2-1
        MinlSpin=Minl*2-1
        MaxiSpin=Maxi*2
        MaxjSpin=Maxj*2
        MaxkSpin=Maxk*2
        MaxlSpin=Maxl*2

        write(6,*) "Checking all permutations stored. Spin-orbital limits are: "
        write(6,*) "i: ",MiniSpin," -> ",MaxiSpin
        write(6,*) "j: ",MinjSpin," -> ",MaxjSpin
        write(6,*) "k: ",MinkSpin," -> ",MaxkSpin
        write(6,*) "l: ",MinlSpin," -> ",MaxlSpin

        !check
        do i=MiniSpin,MaxiSpin
            do j=MinjSpin,MaxjSpin
                do k=MinkSpin,MaxkSpin
                    do l=MinlSpin,MaxlSpin
                        !Fill the array with all permutations of the integrals
                        if(arr(i,j,k,l).gt.1.D-12) then

                            !<ji|lk>
                            if((j.le.MaxiSpin).and.(j.ge.MiniSpin).and.(i.le.MaxjSpin).and.(i.ge.MinjSpin).and.(l.le.MaxkSpin) &
                                .and.(l.ge.MinkSpin).and.(k.le.MaxlSpin).and.(k.ge.MinlSpin)) then
                                if((abs(arr(j,i,l,k)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif

                            !<kl|ij>
                            if((k.le.MaxiSpin).and.(k.ge.MiniSpin).and.(l.le.MaxjSpin).and.(l.ge.MinjSpin).and.(i.le.MaxkSpin) &
                                .and.(i.ge.MinkSpin).and.(j.le.MaxlSpin).and.(j.ge.MinlSpin)) then
                                if((abs(arr(k,l,i,j)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif

                            !<lk|ji>
                            if((l.le.MaxiSpin).and.(l.ge.MiniSpin).and.(k.le.MaxjSpin).and.(k.ge.MinjSpin).and.(j.le.MaxkSpin) &
                                .and.(j.ge.MinkSpin).and.(i.le.MaxlSpin).and.(i.ge.MinlSpin)) then
                                if((abs(arr(l,k,j,i)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif

                            !<kj|il>
                            if((k.le.MaxiSpin).and.(k.ge.MiniSpin).and.(j.le.MaxjSpin).and.(j.ge.MinjSpin).and.(i.le.MaxkSpin) &
                                .and.(i.ge.MinkSpin).and.(l.le.MaxlSpin).and.(l.ge.MinlSpin)) then
                                if((abs(arr(k,j,i,l)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif
                            
                            !<li|jk>
                            if((l.le.MaxiSpin).and.(l.ge.MiniSpin).and.(i.le.MaxjSpin).and.(i.ge.MinjSpin).and.(j.le.MaxkSpin) &
                                .and.(j.ge.MinkSpin).and.(k.le.MaxlSpin).and.(k.ge.MinlSpin)) then
                                if((abs(arr(l,i,j,k)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif

                            !<il|kj>
                            if((i.le.MaxiSpin).and.(i.ge.MiniSpin).and.(l.le.MaxjSpin).and.(l.ge.MinjSpin).and.(k.le.MaxkSpin) &
                                .and.(k.ge.MinkSpin).and.(j.le.MaxlSpin).and.(j.ge.MinlSpin)) then
                                if((abs(arr(i,l,k,j)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif

                            !<jk|li>
                            if((j.le.MaxiSpin).and.(j.ge.MiniSpin).and.(k.le.MaxjSpin).and.(k.ge.MinjSpin).and.(l.le.MaxkSpin) &
                                .and.(l.ge.MinkSpin).and.(i.le.MaxlSpin).and.(i.ge.MinlSpin)) then
                                if((abs(arr(j,k,l,i)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif

                        endif
                    enddo
                enddo
            enddo
        enddo

    end subroutine fill4indarrspin

    logical function AllowedSpin(i,j,k,l)
        integer , intent(in) :: i,j,k,l !spin orbital indices of integral

        if(((mod(i,2).eq.0).and.(mod(k,2).eq.0)).or.((mod(i,2).eq.1).and.(mod(k,2).eq.1))) then
            AllowedSpin=.true.
        else
            AllowedSpin=.false.
            return
        endif

        if(((mod(j,2).eq.0).and.(mod(l,2).eq.0)).or.((mod(j,2).eq.1).and.(mod(l,2).eq.1))) then
            return
        else
            AllowedSpin=.false.
            return
        endif

    end function AllowedSpin
    
    !Fill a two-index, one-electron array with integrals from a file.
    !The basis indicators, indicate the basis which the array wants to run over.
    !We do not assume any permutational symmetry at all.
    !Upper and lower limits must be the same, and the matrix square.
    !Reading in as SPIN orbitals.
    subroutine fill2indarr_spin(fu,arr,TwoIndLims)
        integer , intent(in) :: fu,TwoIndLims(2) !Just minimum and maximum index of a square array.
        real(8) , intent(out) :: arr(TwoIndLims(1):TwoIndLims(2),TwoIndLims(1):TwoIndLims(2))
        integer(2) :: MaxFilei,MaxFilej
        integer(8) :: MaxLength,Length
        integer(2) , allocatable :: Indices(:,:)
        real(8) , allocatable :: Buf(:)
        real(8) :: xout(2)
        integer :: iBuffer,ii,jj,i,ierr

        if(TwoIndLims(2).lt.TwoIndLims(1)) stop 'Error with 1 electron integral limits'
        call FindMax2Indices(fu,MaxFilei,MaxFilej)
!        write(6,"(A,2I5)") "Index upper limits for the file are: ",MaxFilei,MaxFilej

        if((TwoIndLims(2).gt.MaxFilei*2).or.(TwoIndLims(2).gt.MaxFilej*2)) then
            stop 'Indices on disk not sufficient to read in 1e- matrix.'
        endif

        arr=0.D0
        rewind(fu)
        if(.not.tformattedints) then
            read(fu) maxlength,xout(2)
        else
            maxlength = 1
            length = 1
        endif
        allocate(Indices(2,MaxLength),Buf(MaxLength))
        iBuffer=0
        do while(.true.)

            iBuffer=iBuffer+1
!            WRITE(6,*) "Reading buffer: ",iBuffer

            if(tformattedints) then
                read(fu,*,iostat=ierr) Buf(1),Indices(1,1),Indices(2,1)
                if(ierr.lt.0) then
                    goto 13
                elseif(ierr.gt.0) then
                    goto 12
                endif
            else
                read(fu,err=12) length,Indices(1:2,1:length),Buf(1:length)
                IF(length.le.0) goto 13
            endif
!            WRITE(6,"(I12,A)") length," integrals read. Writing Integrals..."

            do i=1,length
                !Its read in in *chemical notation*...immediately switch to physical
                ii=Indices(1,i)
                jj=Indices(2,i)
!                WRITE(luWrite,'(1X,G20.12,4I3)') Buf(i),Indices(1,i),Indices(2,i),Indices(3,i),Indices(4,i)

                IF((ii*2.le.TwoIndLims(2)).and.(ii*2.ge.TwoIndLims(1)).and.(jj*2.le.TwoIndLims(2)).and.(jj*2.ge.TwoIndLims(1))) then
                    !beta orbitals
                    arr((ii*2)-1,(jj*2)-1)=Buf(i)
                    arr((jj*2)-1,(ii*2)-1)=Buf(i)
                    !alpha
                    arr((ii*2),(jj*2))=Buf(i)
                    arr((jj*2),(ii*2))=Buf(i)

                    !test
!                    arr((ii*2),(jj*2)-1)=Buf(i)
!                    arr((jj*2),(ii*2)-1)=Buf(i)
!                    arr((jj*2)-1,(ii*2))=Buf(i)
!                    arr((ii*2)-1,(jj*2))=Buf(i)
                    
                ENDIF
            enddo
!            WRITE(6,*) "Finished integral buffer."
        enddo

12      STOP 'Error reading block'

13      CONTINUE

!        do i=1,TwoIndLims(2)
!            do j=1,TwoIndLims(2)
!                write(6,*) i,j,arr(i,j)
!            enddo
!        enddo

!        WRITE(6,*) "All Integrals read"

        DEALLOCATE(Indices)
        DEALLOCATE(Buf)

    end subroutine fill2indarr_spin

    !Fill a two-index, one-electron array with integrals from a file.
    !The basis indicators, indicate the basis which the array wants to run over.
    !We do not assume any permutational symmetry at all.
    !Upper and lower limits must be the same, and the matrix square.
    !Reading in as SPATIAL orbitals.
    subroutine fill2indarr(fu,arr,TwoIndLims)
        integer , intent(in) :: fu,TwoIndLims(2) !Just minimum and maximum index of a square array.
        real(8) , intent(out) :: arr(TwoIndLims(1):TwoIndLims(2),TwoIndLims(1):TwoIndLims(2))
        integer(2) :: MaxFilei,MaxFilej
        integer(8) :: MaxLength,Length
        integer(2) , allocatable :: Indices(:,:)
        real(8) , allocatable :: Buf(:)
        real(8) :: xout(2)
        integer :: iBuffer,ii,jj,i

        if(TwoIndLims(2).lt.TwoIndLims(1)) stop 'Error with 1 electron integral limits'
        call FindMax2Indices(fu,MaxFilei,MaxFilej)
!        write(6,"(A,2I5)") "Index upper limits for the file are: ",MaxFilei,MaxFilej

        if((TwoIndLims(2).gt.MaxFilei).or.(TwoIndLims(2).gt.MaxFilej)) then
            stop 'Indices on disk not sufficient to read in 1e- matrix.'
        endif

        arr=0.D0
        rewind(fu)
        read(fu) maxlength,xout(2)
!        write(6,*) "Maxlength = ",maxlength
        allocate(Indices(2,MaxLength),Buf(MaxLength))
        iBuffer=0
        do while(.true.)

            iBuffer=iBuffer+1
!            WRITE(6,*) "Reading buffer: ",iBuffer

            read(fu,err=12) length,Indices(1:2,1:length),Buf(1:length)
            IF(length.le.0) goto 13
!            WRITE(6,"(I12,A)") length," integrals read. Writing Integrals..."

            do i=1,length
                !Its read in in *chemical notation*...immediately switch to physical
                ii=Indices(1,i)
                jj=Indices(2,i)
!                WRITE(luWrite,'(1X,G20.12,4I3)') Buf(i),Indices(1,i),Indices(2,i),Indices(3,i),Indices(4,i)

                IF((ii.le.TwoIndLims(2)).and.(ii.ge.TwoIndLims(1)).and.(jj.le.TwoIndLims(2)).and.(jj.ge.TwoIndLims(1))) then
                    arr(ii,jj)=Buf(i)
                    arr(jj,ii)=Buf(i)
                ENDIF
            enddo
!            WRITE(6,*) "Finished integral buffer."
        enddo

12      STOP 'Error reading block'

13      CONTINUE

!        WRITE(6,*) "All Integrals read"

        DEALLOCATE(Indices)
        DEALLOCATE(Buf)

    end subroutine fill2indarr

    !Fill a four-index array with integrals from a file.
    !The basis indicators, indicate the basis which the array wants to run over.
    !We assume an 8-fold permutational symmetry with these integrals.
    subroutine fill4indarr(fu,arr,Lim)
        
        integer, intent(in) :: Lim(2,4)  !Limits(1,:) = minimum indices of array dimensions ; (2,:) = max.
        real(8), intent(out) :: arr(Lim(1,1):Lim(2,1),Lim(1,2):Lim(2,2),Lim(1,3):Lim(2,3),Lim(1,4):Lim(2,4))
        integer , intent(in) :: fu
        integer(8) :: MaxLength,length
        integer(2) , allocatable :: Indices(:,:)
        real(8) , allocatable :: Buf(:)
        integer :: iBuffer,i,j,k,l,ii,jj,kk,ll,Mini,Minj,Mink,Minl,Maxi,Maxj,Maxk,Maxl,ierr
        integer(2) :: MaxFilei,MaxFilej,MaxFilek,MaxFilel
        
        CALL FindMaxIndices(fu,MaxFilei,MaxFilej,MaxFilek,MaxFilel)
!        write(6,"(A,4I5)") "Index upper limits for the file are (physical notation): ",MaxFilei,MaxFilej,MaxFilek,MaxFilel
        
        call check_limits(MaxFilei,MaxFilej,MaxFilek,MaxFilel,Lim)

        Mini=Lim(1,1)
        Maxi=Lim(2,1)
        Minj=Lim(1,2)
        Maxj=Lim(2,2)
        Mink=Lim(1,3)
        Maxk=Lim(2,3)
        Minl=Lim(1,4)
        Maxl=Lim(2,4)

        arr=0.D0

        rewind(fu)

        if(.not.tformattedints) then
            read(fu) maxlength
        else
            maxlength = 1
            length = maxlength
        endif
!        write(6,*) "maxlength = ",maxlength
        allocate(Indices(4,MaxLength),Buf(MaxLength))

        iBuffer=0
        do while(.true.)

            iBuffer=iBuffer+1
!            WRITE(6,*) "Reading buffer: ",iBuffer

            if(tformattedints) then
                read(fu,*,iostat=ierr) Buf(1),Indices(1:4,1)
                if(ierr.lt.0) then
                    goto 13 !eof
                elseif(ierr.gt.0) then
                    goto 12
                endif
                if((Indices(1,1).eq.0).or.(Indices(2,1).eq.0).or.(Indices(3,1).eq.0).or.(Indices(4,1).eq.0)) then
                    cycle   !We only want to pick up 4 index integrals
                endif
            else
                read(fu,err=12) length,Indices(1:4,1:length),Buf(1:length)
                IF(length.le.0) goto 13
            endif
!            WRITE(6,"(I12,A)") length," integrals read. Writing Integrals..."

            do i=1,length
                !Its read in in *chemical notation*...immediately switch to physical
                ii=Indices(1,i)
                kk=Indices(2,i)
                jj=Indices(3,i)
                ll=Indices(4,i)
!                WRITE(luWrite,'(1X,G20.12,4I3)') Buf(i),Indices(1,i),Indices(2,i),Indices(3,i),Indices(4,i)

                IF((ii.le.Maxi).and.(ii.ge.Mini).and.(jj.le.Maxj).and.(jj.ge.Minj).and.(kk.le.Maxk) &
                    .and.(kk.ge.Minj).and.(ll.le.Maxl).and.(ll.ge.Minl)) THEN
                    arr(ii,jj,kk,ll)=Buf(i)
                ENDIF
                !Now for rest of permutations...
  
                !<ji|lk>
                IF((jj.le.Maxi).and.(jj.ge.Mini).and.(ii.le.Maxj).and.(ii.ge.Minj).and.(ll.le.Maxk) &
                    .and.(ll.ge.Minj).and.(kk.le.Maxl).and.(kk.ge.Minl)) THEN
                    arr(jj,ii,ll,kk)=Buf(i)
                endif

                !<kl|ij>
                IF((kk.le.Maxi).and.(kk.ge.Mini).and.(ll.le.Maxj).and.(ll.ge.Minj).and.(ii.le.Maxk) &
                    .and.(ii.ge.Minj).and.(jj.le.Maxl).and.(jj.ge.Minl)) THEN
                    arr(kk,ll,ii,jj)=Buf(i)
                endif

                !<lk|ji>
                IF((ll.le.Maxi).and.(ll.ge.Mini).and.(kk.le.Maxj).and.(kk.ge.Minj).and.(jj.le.Maxk) &
                    .and.(jj.ge.Minj).and.(ii.le.Maxl).and.(ii.ge.Minl)) THEN
                    arr(ll,kk,jj,ii)=Buf(i)
                endif

                !<kj|il>
                IF((kk.le.Maxi).and.(kk.ge.Mini).and.(jj.le.Maxj).and.(jj.ge.Minj).and.(ii.le.Maxk) &
                    .and.(ii.ge.Minj).and.(ll.le.Maxl).and.(ll.ge.Minl)) THEN
                    arr(kk,jj,ii,ll)=Buf(i)
                endif

                !<li|jk>
                IF((ll.le.Maxi).and.(ll.ge.Mini).and.(ii.le.Maxj).and.(ii.ge.Minj).and.(jj.le.Maxk) &
                    .and.(jj.ge.Minj).and.(kk.le.Maxl).and.(kk.ge.Minl)) THEN
                    arr(ll,ii,jj,kk)=Buf(i)
                endif

                !<il|kj>
                IF((ii.le.Maxi).and.(ii.ge.Mini).and.(ll.le.Maxj).and.(ll.ge.Minj).and.(kk.le.Maxk) &
                    .and.(kk.ge.Minj).and.(jj.le.Maxl).and.(jj.ge.Minl)) THEN
                    arr(ii,ll,kk,jj)=Buf(i)
                endif

                !<jk|li>
                IF((jj.le.Maxi).and.(jj.ge.Mini).and.(kk.le.Maxj).and.(kk.ge.Minj).and.(ll.le.Maxk) &
                    .and.(ll.ge.Minj).and.(ii.le.Maxl).and.(ii.ge.Minl)) THEN
                    arr(jj,kk,ll,ii)=Buf(i)
                endif
            enddo

!            WRITE(6,*) "Finished integral buffer."

        enddo

12      STOP 'Error reading block'

13      CONTINUE

!        WRITE(6,*) "All Integrals read"

        DEALLOCATE(Indices)
        DEALLOCATE(Buf)

        !check
        do i=Mini,Maxi
            do j=Minj,Maxj
                do k=Mink,Maxk
                    do l=Minl,Maxl
                        !Fill the array with all permutations of the integrals
                        if(arr(i,j,k,l).gt.1.D-12) then

                            !<ji|lk>
                            if((j.le.Maxi).and.(j.ge.Mini).and.(i.le.Maxj).and.(i.ge.Minj).and.(l.le.Maxk) &
                                .and.(l.ge.Mink).and.(k.le.Maxl).and.(k.ge.Minl)) then
                                if((abs(arr(j,i,l,k)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif

                            !<kl|ij>
                            if((k.le.Maxi).and.(k.ge.Mini).and.(l.le.Maxj).and.(l.ge.Minj).and.(i.le.Maxk) &
                                .and.(i.ge.Mink).and.(j.le.Maxl).and.(j.ge.Minl)) then
                                if((abs(arr(k,l,i,j)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif

                            !<lk|ji>
                            if((l.le.Maxi).and.(l.ge.Mini).and.(k.le.Maxj).and.(k.ge.Minj).and.(j.le.Maxk) &
                                .and.(j.ge.Mink).and.(i.le.Maxl).and.(i.ge.Minl)) then
                                if((abs(arr(l,k,j,i)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif

                            !<kj|il>
                            if((k.le.Maxi).and.(k.ge.Mini).and.(j.le.Maxj).and.(j.ge.Minj).and.(i.le.Maxk) &
                                .and.(i.ge.Mink).and.(l.le.Maxl).and.(l.ge.Minl)) then
                                if((abs(arr(k,j,i,l)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif
                            
                            !<li|jk>
                            if((l.le.Maxi).and.(l.ge.Mini).and.(i.le.Maxj).and.(i.ge.Minj).and.(j.le.Maxk) &
                                .and.(j.ge.Mink).and.(k.le.Maxl).and.(k.ge.Minl)) then
                                if((abs(arr(l,i,j,k)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif

                            !<il|kj>
                            if((i.le.Maxi).and.(i.ge.Mini).and.(l.le.Maxj).and.(l.ge.Minj).and.(k.le.Maxk) &
                                .and.(k.ge.Mink).and.(j.le.Maxl).and.(j.ge.Minl)) then
                                if((abs(arr(i,l,k,j)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif

                            !<jk|li>
                            if((j.le.Maxi).and.(j.ge.Mini).and.(k.le.Maxj).and.(k.ge.Minj).and.(l.le.Maxk) &
                                .and.(l.ge.Mink).and.(i.le.Maxl).and.(i.ge.Minl)) then
                                if((abs(arr(j,k,l,i)-arr(i,j,k,l)).gt.1.D-8)) then
                                    stop 'error in reading in integrals'
                                endif
                            endif

                        endif
                    enddo
                enddo
            enddo
        enddo

    end subroutine fill4indarr
    
    subroutine check_limits(Maxi,Maxj,Maxk,Maxl,Lims)
        integer(2) , intent(in) :: Maxi,Maxj,Maxk,Maxl
        integer , intent(in) :: Lims(2,4)
        logical :: LimitsOK

        LimitsOK=.false.

        !<ij|kl>
        if((Lims(2,1).gt.Maxi).or.(Lims(2,2).gt.Maxj).or.(Lims(2,3).gt.Maxk).or.(Lims(2,4).gt.Maxl)) then
            LimitsOK=.false.
        else
            return
        endif
        !<ji|lk>
        if((Lims(2,2).gt.Maxi).or.(Lims(2,1).gt.Maxj).or.(Lims(2,4).gt.Maxk).or.(Lims(2,3).gt.Maxl)) then
            LimitsOK=.false.
        else
            return
        endif
        !<kl|ij>
        if((Lims(2,3).gt.Maxi).or.(Lims(2,4).gt.Maxj).or.(Lims(2,1).gt.Maxk).or.(Lims(2,2).gt.Maxl)) then
            LimitsOK=.false.
        else
            return
        endif
        !<lk|ji>
        if((Lims(2,4).gt.Maxi).or.(Lims(2,3).gt.Maxj).or.(Lims(2,2).gt.Maxk).or.(Lims(2,1).gt.Maxl)) then
            LimitsOK=.false.
        else
            return
        endif
        !<kj|il>
        if((Lims(2,3).gt.Maxi).or.(Lims(2,2).gt.Maxj).or.(Lims(2,1).gt.Maxk).or.(Lims(2,4).gt.Maxl)) then
            LimitsOK=.false.
        else
            return
        endif
        !<li|jk>
        if((Lims(2,4).gt.Maxi).or.(Lims(2,1).gt.Maxj).or.(Lims(2,2).gt.Maxk).or.(Lims(2,3).gt.Maxl)) then
            LimitsOK=.false.
        else
            return
        endif
        !<il|kj>
        if((Lims(2,1).gt.Maxi).or.(Lims(2,4).gt.Maxj).or.(Lims(2,3).gt.Maxk).or.(Lims(2,2).gt.Maxl)) then
            LimitsOK=.false.
        else
            return
        endif
        !<jk|li>
        if((Lims(2,2).gt.Maxi).or.(Lims(2,3).gt.Maxj).or.(Lims(2,4).gt.Maxk).or.(Lims(2,1).gt.Maxl)) then
            LimitsOK=.false.
        else
            return
        endif

        if(.not.LimitsOK) then
            stop 'Index limits in file not sufficient'
        endif

    end subroutine check_limits
    
    subroutine CalculateHF()
        integer :: Limits(2,4),ierr,i,j
        real(8) :: HFEnergy

!        write(6,*) "Getting g_rs^xy tensor from MO_G file..."
        unit_g=get_free_unit()
        if(tformattedints) then
            open(unit_g,file='FCIDUMP',status='old',form='formatted')
        else
            open(unit_g,file='MO_G',status='old',form='unformatted',access='sequential',action='read')
        endif
        
        allocate(gints(1:nOrbBasis,1:nOrbBasis,1:nOrbBasis,1:nOrbBasis),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'

        !Fill the fg array (in physical notation). All indices just over orbital basis
        Limits(1,1)=1
        Limits(2,1)=nOrbBasis
        Limits(1,2)=1
        Limits(2,2)=nOrbBasis
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis

        call fill4indarr(unit_g,gints,Limits)
!        call fill4indarrspin(unit_g,gints,Limits)

!        write(6,*) "Antisymmetrizing g integrals"
        allocate(antisym_g(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        
        Limits(1,1)=1
        Limits(2,1)=nOrbBasis*2
        Limits(1,2)=1
        Limits(2,2)=nOrbBasis*2
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis*2
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis*2

        call antisymints(gints,antisym_g,Limits)

        deallocate(gints)

        !Now calculate HF

        HFEnergy=0.D0
        do i=1,HFOccOrbs*2
            HFEnergy=HFEnergy+OrbEnergies(i)
        enddo

        do i=1,HFOccOrbs*2
            do j=1,HFOccOrbs*2
                HFEnergy=HFEnergy-antisym_g(i,j,i,j)/2.D0
            enddo
        enddo

        HFEnergy=HFEnergy+NucRep

        write(6,*) "*** HF Energy is: ",HFEnergy
        write(6,*) ""
        close(unit_g)
        deallocate(antisym_g)

    end subroutine CalculateHF 

    subroutine CalculateRelaxation()
        integer :: Limits(2,4),ierr,i,j,k,l,SpinType
        real(8) , allocatable :: RelaxX(:,:)
        real(8) :: MP2Energy

!        write(6,*) "Getting g_rs^xy tensor from MO_G file..."
        unit_g=get_free_unit()
        if(tformattedints) then
            open(unit_g,file='FCIDUMP',status='old',form='formatted')
        else
            open(unit_g,file='MO_G',status='old',form='unformatted',access='sequential',action='read')
        endif
        
        if(allocated(gints)) deallocate(gints)
        allocate(gints(1:nOrbBasis,1:nOrbBasis,1:nOrbBasis,1:nOrbBasis),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'

        !Fill the fg array (in physical notation). All indices just over orbital basis
        Limits(1,1)=1
        Limits(2,1)=nOrbBasis
        Limits(1,2)=1
        Limits(2,2)=nOrbBasis
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis

        call fill4indarr(unit_g,gints,Limits)
!        call fill4indarrspin(unit_g,gints,Limits)

!        write(6,*) "Antisymmetrizing g integrals"
        if(allocated(antisym_g)) deallocate(antisym_g)
        allocate(antisym_g(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'

        do i=1,nOrbBasis*2
            do j=1,nOrbBasis*2
                do k=1,nOrbBasis*2
                    do l=1,nOrbBasis*2
                        SpinType=DetermineSpinType(i,j,k,l)
                        if(spinType.eq.1) then
                            !<aa|aa>
                            antisym_g(i,j,k,l)=gints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(k),Conv2Spat(l))! &
!                                -ints(Conv2Spat(j),Conv2Spat(i),Conv2Spat(k),Conv2Spat(l))
                        elseif(SpinType.eq.2) then
                            !<ab||ab>
                            antisym_g(i,j,k,l)=gints(Conv2Spat(i),Conv2Spat(j),Conv2Spat(k),Conv2Spat(l))
!                        elseif(SpinType.eq.3) then
!                            !<ab||ba>
!                            antisym_ints(i,j,k,l)=-ints(Conv2Spat(j),Conv2Spat(i),Conv2Spat(k),Conv2Spat(l))
                        else
                            antisym_g(i,j,k,l)=0.D0
                        endif
                    enddo
                enddo
            enddo
        enddo

!        Limits(1,1)=1
!        Limits(2,1)=nOrbBasis*2
!        Limits(1,2)=1
!        Limits(2,2)=nOrbBasis*2
!        Limits(1,3)=1
!        Limits(2,3)=nOrbBasis*2
!        Limits(1,4)=1
!        Limits(2,4)=nOrbBasis*2
!        call antisymints(gints,antisym_g,Limits)


        deallocate(gints)

        allocate(RelaxX(nOrbBasis*2,nOrbBasis*2))
        RelaxX(:,:) = 0.D0

        !This is for Dyall hamiltonian
        do y=1,nOrbBasis*2
            do x=1,nOrbBasis*2
                do k=1,nOrbBasis*2
                    do i=1,nOrbBasis*2
                        do j=1,nOrbBasis*2
                            RelaxX(y,x) = RelaxX(y,x) - antisym_g(y,k,i,j)*Cumulant(i,j,x,k)/2.D0
!                            RelaxX(y,x) = RelaxX(y,x) + antisym_g(y,k,i,j)*TwoRDM(i,j,x,k)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        do y=1,nOrbBasis*2
            do x=1,nOrbBasis*2
                do i=1,nOrbBasis*2
                    RelaxX(y,x) = RelaxX(y,x) - fock(y,i)*OneRDM(i,x)
                enddo
            enddo
        enddo

!        !This is for Fock hamiltonian
!        do y=1,nOrbBasis*2
!            do x=1,nOrbBasis*2
!                do i=1,nOrbBasis*2
!                    do j=1,nOrbBasis*2
!                        RelaxX(y,x) = RelaxX(y,x) + fock(j,i)*Cumulant(y,i,x,j)
!!                        RelaxX(y,x) = RelaxX(y,x) + fock(j,i)*TwoRDM(y,i,x,j)
!                    enddo
!                enddo
!            enddo
!        enddo
!
!        do y=1,nOrbBasis*2
!            do x=1,nOrbBasis*2
!                do i=1,nOrbBasis*2
!                    do j=1,nOrbBasis*2
!                        RelaxX(y,x) = RelaxX(y,x) - fock(j,i)*OneRDM(i,x)*OneRDM(y,j)
!!                        RelaxX(y,x) = RelaxX(y,x) - fock(j,i)*OneRDM(y,x)*OneRDM(i,j)
!                    enddo
!                enddo
!            enddo
!        enddo

        open(85,file='Relax-Model')
        do y=1,nOrbBasis*2,2
            do x=1,nOrbBasis*2,2
                write(85,*) (y+1)/2,(x+1)/2,RelaxX(y,x)
            enddo
        enddo
        close(85)

    end subroutine CalculateRelaxation

    subroutine CalculateMP2()
        integer :: Limits(2,4),ierr,i,j,a,b
        real(8) :: MP2Energy

!        write(6,*) "Getting g_rs^xy tensor from MO_G file..."
        unit_g=get_free_unit()
        if(tformattedints) then
            open(unit_g,file='FCIDUMP',status='old',form='formatted')
        else
            open(unit_g,file='MO_G',status='old',form='unformatted',access='sequential',action='read')
        endif
        
        allocate(gints(1:nOrbBasis,1:nOrbBasis,1:nOrbBasis,1:nOrbBasis),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'

        !Fill the fg array (in physical notation). All indices just over orbital basis
        Limits(1,1)=1
        Limits(2,1)=nOrbBasis
        Limits(1,2)=1
        Limits(2,2)=nOrbBasis
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis

        call fill4indarr(unit_g,gints,Limits)
!        call fill4indarrspin(unit_g,gints,Limits)

!        write(6,*) "Antisymmetrizing g integrals"
        allocate(antisym_g(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        Limits(1,1)=1
        Limits(2,1)=nOrbBasis*2
        Limits(1,2)=1
        Limits(2,2)=nOrbBasis*2
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis*2
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis*2

        call antisymints(gints,antisym_g,Limits)

        deallocate(gints)

        !Now calculate MP2

        MP2Energy=0.D0
        do i=1,HFOccOrbs*2
            do j=1,HFOccOrbs*2
                do a=HFOccOrbs*2+1,nOrbBasis*2
                    do b=HFOccOrbs*2+1,nOrbBasis*2
                        MP2Energy=MP2Energy+(antisym_g(i,j,a,b)**2)/(OrbEnergies(i)+OrbEnergies(j)-OrbEnergies(a)-OrbEnergies(b))
!                        if((i.eq.j).or.(a.eq.b)) then
!                            write(6,*) "contrib: ",i,j,a,b,(antisym_g(i,j,a,b)**2)/(OrbEnergies(i)+OrbEnergies(j)-OrbEnergies(a)-OrbEnergies(b))
!                        endif
                    enddo
                enddo
            enddo
        enddo

        write(6,*) "*** MP2 Energy is: ",MP2Energy/4.D0
        write(6,*) ""
        close(unit_g)
        deallocate(antisym_g)

    end subroutine CalculateMP2


    SUBROUTINE WriteFCIDUMP()
        IMPLICIT NONE
        INTEGER :: i,iBuffer,nBuffers
        INTEGER(8) :: Length,MaxLength
        INTEGER(2) :: Maxi,Maxj,Maxk,Maxl

        return  !Ignore this subroutine

        OPEN(8,FILE='MO_G',STATUS='old',FORM='UNFORMATTED',access='sequential')
        rewind(8)
        OPEN(10,FILE='FCIDUMP',STATUS='unknown',form='formatted')
        read(8) maxlength
        WRITE(6,*) "Max Buffer Length: ",maxlength

        CALL FindNoBuffers(8,nBuffers)

        WRITE(6,*) "Number of integral buffers in the file: ",nBuffers

        CALL FindMaxIndices(8,Maxi,Maxj,Maxk,Maxl)

        WRITE(6,*) "Max indices (chemical notation) are: ",Maxi,Maxj,Maxk,Maxl

        Maxi=nOrbBasis
        Maxj=nOrbBasis
        Maxk=nOrbBasis
        Maxl=nOrbBasis

        CALL Write4IndIntegrals(8,10,Maxi,Maxj,Maxk,Maxl)

        CLOSE(8)
        CLOSE(10)

    END SUBROUTINE WriteFCIDUMP
    
    SUBROUTINE Write4IndIntegrals(luRead,luWrite,Maxi,Maxj,Maxk,Maxl)
        IMPLICIT NONE
        INTEGER , INTENT(IN) :: luRead,luWrite
        INTEGER(2) , INTENT(IN) :: Maxi,Maxj,Maxk,Maxl
        INTEGER(2) , ALLOCATABLE :: Indices(:,:)
        REAL*8 , ALLOCATABLE :: Buf(:)
        INTEGER(8) :: MaxLength,Length
        INTEGER :: iBuffer,i

        rewind(luRead)
        read(luRead) MaxLength

        ALLOCATE(Indices(4,MaxLength),Buf(MaxLength))

        iBuffer=0
        do while(.true.)

            iBuffer=iBuffer+1
            WRITE(6,*) "Reading buffer: ",iBuffer

            read(luRead,err=12) length,Indices(1:4,1:length),Buf(1:length)
            IF(length.le.0) goto 13
            WRITE(6,"(I12,A)") length," integrals read. Writing Integrals..."

            do i=1,length
                !Chemical notation...
                IF((Indices(1,i).le.Maxi).and.(Indices(2,i).le.Maxj).and.(Indices(3,i).le.Maxk).and.(Indices(4,i).le.Maxl)) THEN
                    WRITE(luWrite,'(1X,G20.12,4I3)') Buf(i),Indices(1,i),Indices(2,i),Indices(3,i),Indices(4,i)
                ENDIF
            enddo

            WRITE(6,*) "Finished integral buffer."

        enddo

12      STOP 'Error reading block'

13      CONTINUE

        WRITE(6,*) "All Integrals read"

        DEALLOCATE(Indices)
        DEALLOCATE(Buf)

    END SUBROUTINE Write4IndIntegrals

    SUBROUTINE FindNoBuffers(iUnit,nBuffers)
        IMPLICIT NONE
        INTEGER , INTENT(IN) :: iUnit
        INTEGER , INTENT(OUT) :: nBuffers
        INTEGER(8) :: MaxLength, length
        REWIND(iUnit)
        READ(iUnit) MaxLength

        nBuffers=0
        do while(.true.)
            read(8,err=10) length
            IF(length.gt.0) THEN
                nBuffers=nBuffers+1
            ELSE
                EXIT
            ENDIF
        enddo
        REWIND(8)
        RETURN
10      STOP 'Error in finding number of buffers'
    END SUBROUTINE FindNoBuffers

!The maximum indices in a particular file, where the indices correspond to a *1 electron integral*!
    SUBROUTINE FindMax2Indices(iUnit,Maxi,Maxj)
        IMPLICIT NONE
        INTEGER(2) , INTENT(OUT) :: Maxi,Maxj
        INTEGER , INTENT(IN) :: iUnit
        INTEGER(8) :: MaxLength,Length
        INTEGER(2) , ALLOCATABLE :: Indices(:,:)
        INTEGER :: i,iBuf,j,ierr
        REAL(8) :: xout(2),integral

        rewind(iUnit)
        Maxi=0
        Maxj=0
        if(tformattedints) then
            do while(.true.)
                read(iunit,*,iostat=ierr) integral,i,j
                if(ierr.lt.0) then
                    exit
                elseif(ierr.gt.0) then
                    goto 11
                endif
                if(i.gt.Maxi) Maxi = i
                if(j.gt.Maxj) Maxj = j
            enddo
        else
            read(iUnit) MaxLength,xout(2)
            ALLOCATE(Indices(2,MaxLength))
            iBuf=0
            do while(.true.)
                iBuf=iBuf+1
    !            WRITE(6,*) "Attempting to read buffer: ",iBuf
                read(iUnit,err=11) length,Indices(1:2,1:length)
                IF(length.gt.0) THEN
                    do i=1,length
                        IF(Indices(1,i).gt.Maxi) Maxi=Indices(1,i)
                        IF(Indices(2,i).gt.Maxj) Maxj=Indices(2,i)
                    enddo
                ELSE
                    EXIT
                ENDIF
            enddo
            DEALLOCATE(Indices)
        endif
        REWIND(iUnit)
        RETURN
11      STOP 'Error when finding largest indices'
    END SUBROUTINE FindMax2Indices

!The maximum indices in a particular file, where the indices are returned in PHYSICAL NOTATION!
    SUBROUTINE FindMaxIndices(iUnit,Maxi,Maxj,Maxk,Maxl)
        IMPLICIT NONE
        INTEGER(2) , INTENT(OUT) :: Maxi,Maxj,Maxk,Maxl
        INTEGER , INTENT(IN) :: iUnit
        INTEGER(8) :: MaxLength,Length
        INTEGER(2) , ALLOCATABLE :: Indices(:,:)
        INTEGER :: i,iBuf,j,k,l,ierr
        real(8) :: Integral

        rewind(iUnit)
        Maxi=0
        Maxj=0
        Maxk=0
        Maxl=0
        if(tformattedints) then
            do while(.true.)
                read(iUnit,*,iostat=ierr) integral,i,k,j,l
                if(ierr.lt.0) then
                    exit    !eof
                elseif(ierr.gt.0) then
                    goto 11
                endif
                if(i.gt.Maxi) Maxi = i
                if(j.gt.Maxj) Maxj = j
                if(k.gt.Maxk) Maxk = k
                if(l.gt.Maxl) Maxl = l
            enddo
        else
            read(iUnit) MaxLength
            ALLOCATE(Indices(4,MaxLength))
            iBuf=0
            do while(.true.)
                iBuf=iBuf+1
    !            WRITE(6,*) "Attempting to read buffer: ",iBuf
                read(iUnit,err=11) length,Indices(1:4,1:length)
                IF(length.gt.0) THEN
                    do i=1,length
                        IF(Indices(1,i).gt.Maxi) Maxi=Indices(1,i)
                        IF(Indices(2,i).gt.Maxk) Maxk=Indices(2,i)
                        IF(Indices(3,i).gt.Maxj) Maxj=Indices(3,i)
                        IF(Indices(4,i).gt.Maxl) Maxl=Indices(4,i)
                    enddo
                ELSE
                    EXIT
                ENDIF
            enddo
            DEALLOCATE(Indices)
        endif
        REWIND(iUnit)
        RETURN
11      STOP 'Error when finding largest indices'
    END SUBROUTINE FindMaxIndices
    
!The minimum indices in a particular file, where the indices are returned in PHYSICAL NOTATION!
    SUBROUTINE FindMinIndices(iUnit,Mini,Minj,Mink,Minl)
        IMPLICIT NONE
        INTEGER(2) , INTENT(OUT) :: Mini,Minj,Mink,Minl
        INTEGER , INTENT(IN) :: iUnit
        INTEGER(8) :: MaxLength,Length
        INTEGER(2) , ALLOCATABLE :: Indices(:,:)
        INTEGER :: i,iBuf

        rewind(iUnit)
        read(iUnit) MaxLength
        ALLOCATE(Indices(4,MaxLength))
        Mini=30000
        Minj=30000
        Mink=30000
        Minl=30000
        iBuf=30000
        do while(.true.)
            iBuf=iBuf+1
            WRITE(6,*) "Attempting to read buffer: ",iBuf
            read(iUnit,err=11) length,Indices(1:4,1:length)
            IF(length.gt.0) THEN
                do i=1,length
                    IF(Indices(1,i).lt.Mini) Mini=Indices(1,i)
                    IF(Indices(2,i).lt.Mink) Mink=Indices(2,i)
                    IF(Indices(3,i).lt.Minj) Minj=Indices(3,i)
                    IF(Indices(4,i).lt.Minl) Minl=Indices(4,i)
                enddo
            ELSE
                EXIT
            ENDIF
        enddo
        DEALLOCATE(Indices)
        REWIND(iUnit)
        RETURN
11      STOP 'Error when finding largest indices'
    END SUBROUTINE FindMinIndices
    
    subroutine find_params()
        INTEGER :: i,j,symorbs(8),total,unit_orb_info,ierr
        character(100) :: junk,junk2
        
        unit_orb_info=get_free_unit()
        open(unit_orb_info,file='OrbitalInfo',status='old',form='formatted',action='read')
        do i=1,5
            !there should be 5 lines in this file, all with the same format
            symorbs=0
            read(unit_orb_info,*) junk,symorbs(:),junk2,total
            
            do j=2,8
                if(symorbs(j).ne.0) stop 'cannot cope with symmetry - turn off symmetry before proceeding'
            enddo
            if(symorbs(1).ne.total) stop 'Not all orbitals accounted for...'

            if(i.eq.1) then
                HFOccOrbs=total
                write(6,*) total," occupied HF orbitals detected."
            elseif(i.eq.2) then
                nOrbBasis=total
                write(6,*) total," orbitals in orbital basis."
            elseif(i.eq.3) then
                nCABSBasis=total
                write(6,*) total, " orbitals in CABS basis."
            elseif(i.eq.5) then
                nOrbsTot=total
                write(6,*) total, " total orbitals in basis."
                if(nOrbsTot.ne.(nOrbBasis+nCABSBasis)) stop 'total orbitals not equal to sum of orbital and CABS basis'
            endif
        enddo

        allocate(OrbEnergies(nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        OrbEnergies=0.D0
        do i=1,nOrbBasis
            read(unit_orb_info,*) junk,OrbEnergies(i*2-1)
            OrbEnergies(i*2)=OrbEnergies(i*2-1)
        enddo

        write(6,*) "Spin orbital energies are: "
        do i=1,nOrbBasis*2
            write(6,*) i,OrbEnergies(i)
        enddo

        read(unit_orb_info,*) junk,junk2,NucRep
        write(6,*) "Nuclear repulsion energy: ",NucRep

        close(unit_orb_info)

    end subroutine find_params
    
    subroutine CalculatePermSym()
        integer :: ierr,Limits(2,4),PermSym,i,j,k,l
        logical :: tBraKetSym,tElecPerm,tPermElecOrbs1,tPermElecOrbs2

        
        write(6,*) "Testing 1-geminal pair spin-orb integral permutational symmetry..."
        unit_fg=get_free_unit()
        if(tformattedints) then
            open(unit_fg,file='FGDUMP',status='old',form='formatted')
        else
            open(unit_fg,file='MO_FG',status='old',form='unformatted',access='sequential',action='read')
        endif
        !Fill with spatial orbital integrals 
        allocate(fgints(1:nOrbBasis,1:nOrbBasis,1:nOrbBasis,1:nOrbBasis),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        !Fill the fg array (in physical notation). All indices just over orbital basis
        Limits(1,1)=1
        Limits(2,1)=nOrbBasis
        Limits(1,2)=1
        Limits(2,2)=nOrbBasis
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis

        call fill4indarr(unit_fg,fgints,Limits)
        close(unit_fg)

!        write(6,'(A)') "Antisymmetrizing fg integrals - converting to spin-orbital representation..."
        allocate(antisym_fg(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        
        !Fill the fg array (in physical notation). All indices just over orbital basis
        Limits(1,1)=1
        Limits(2,1)=nOrbBasis*2
        Limits(1,2)=1
        Limits(2,2)=nOrbBasis*2
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis*2
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis*2

        call antisymcuspints(fgints,antisym_fg,Limits)
        deallocate(fgints)

        tBraKetSym=.true.
        tElecPerm=.true.
        tPermElecOrbs1=.true.
        tPermElecOrbs2=.true.
        do i=1,nOrbBasis*2
            do j=1,nOrbBasis*2
                do k=1,nOrbBasis*2
                    do l=1,nOrbBasis*2
                        !<ij|kl> = <kl|ij> ?
                        if(abs(antisym_fg(i,j,k,l)-antisym_fg(k,l,i,j)).gt.1.D-8) then
                            tBraKetSym=.false.
                        endif
                        if(abs(antisym_fg(i,j,k,l)-antisym_fg(j,i,l,k)).gt.1.D-8) then
                            tElecPerm=.false.
                        endif
                        if(abs(antisym_fg(i,j,k,l)-antisym_fg(k,j,i,l)).gt.1.D-8) then
                            tPermElecOrbs1=.false.
                        endif
                        if(abs(antisym_fg(i,j,k,l)-antisym_fg(i,l,k,j)).gt.1.D-8) then
                            tPermElecOrbs2=.false.
                        endif
                    enddo
                enddo
            enddo
        enddo
        deallocate(antisym_fg)
        if(tPermElecOrbs1.neqv.tPermElecOrbs2) then
            write(6,*) "Permutation of orbitals of a each electron does not have same symmetry! Electrons not equivalent!"
            stop 'Permutation of orbitals of a each electron does not have same symmetry! Electrons not equivalent!'
        endif
        write(6,*) "For 1 geminal pair 4-index spin integrals, permutational symmetry is: "
        PermSym=0
        if(tBraKetSym) then
            write(6,"(A)") "Bra-Ket symmetry -> <ij|kl> = <kl|ij>"
            PermSym=PermSym+1
        else
            write(6,*) "No bra-ket permutational symmetry"
        endif
        if(tElecPerm) then
            write(6,"(A)") "Electron exchange symmetry -> <ij|kl> = <ji|lk>"
            PermSym=PermSym+1
        else
            write(6,*) "No electron exchange symmetry"
        endif
        if(tPermElecOrbs1) then
            write(6,"(A)") "Permutation of orbitals within either electron allowed -> <ij|kl> = <kj|il> = <il|kj>"
            PermSym=PermSym+2
        else
            write(6,*) "No permutation of orbitals within either electron allowed"
        endif
        if(PermSym.eq.0) then
            write(6,"(A)") "No permutational symmetry of the four-index, one-geminal pair antisymetrized spin integrals allowed"
        else
            write(6,"(A,I4)") "Permutational symmetry of the four-index, one-geminal pair antisymmetrized spin integrals is order: ",PermSym*2
        endif
        write(6,*) ""

        write(6,*) "Now testing 2-geminal geminal permutational symmetry..."
        unit_ff=get_free_unit()
        if(tformattedints) then
            open(unit_ff,file='FFDUMP',status='old',form='formatted')
        else
            open(unit_ff,file='MO_FF',status='old',form='unformatted',access='sequential',action='read')
        endif

        !Fill with spatial orbital integrals 
        allocate(ffints(1:nOrbsTot,1:nOrbBasis,1:nOrbBasis,1:nOrbBasis),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'

        !Fill the fg array (in physical notation). All indices just over orbital basis
        Limits(1,1)=1
        Limits(2,1)=nOrbsTot
        Limits(1,2)=1
        Limits(2,2)=nOrbBasis
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis

        call fill4indarr(unit_ff,ffints,Limits)
        close(unit_ff)

        allocate(antisym_ff(1:nOrbsTot*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        
        !Fill the ff array (in physical notation). All indices just over *geminal* basis
        Limits(1,1)=1
        Limits(2,1)=nOrbsTot*2
        Limits(1,2)=1
        Limits(2,2)=nOrbBasis*2
        Limits(1,3)=1
        Limits(2,3)=nOrbBasis*2
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis*2

        call antisym2cuspints(ffints,antisym_ff,Limits)
        deallocate(ffints)
        
        tBraKetSym=.true.
        tElecPerm=.true.
        tPermElecOrbs1=.true.
        tPermElecOrbs2=.true.
        do i=1,nOrbBasis*2
            do j=1,nOrbBasis*2
                do k=1,nOrbBasis*2
                    do l=1,nOrbBasis*2
                        !<ij|kl> = <kl|ij> ?
                        if(abs(antisym_ff(i,j,k,l)-antisym_ff(k,l,i,j)).gt.1.D-8) then
                            tBraKetSym=.false.
                        endif
                        if(abs(antisym_ff(i,j,k,l)-antisym_ff(j,i,l,k)).gt.1.D-8) then
                            tElecPerm=.false.
                        endif
                        if(abs(antisym_ff(i,j,k,l)-antisym_ff(k,j,i,l)).gt.1.D-8) then
                            tPermElecOrbs1=.false.
                        endif
                        if(abs(antisym_ff(i,j,k,l)-antisym_ff(i,l,k,j)).gt.1.D-8) then
                            tPermElecOrbs2=.false.
                        endif
                    enddo
                enddo
            enddo
        enddo
        deallocate(antisym_ff)
        if(tPermElecOrbs1.neqv.tPermElecOrbs2) then
            write(6,*) "Permutation of orbitals of a each electron does not have same symmetry! Electrons not equivalent!"
            stop 'Permutation of orbitals of a each electron does not have same symmetry! Electrons not equivalent!'
        endif
        write(6,*) "For 2 geminal pair 4-index spin antisymmetrized integrals, permutational symmetry is: "
        PermSym=0
        if(tBraKetSym) then
            write(6,"(A)") "Bra-Ket symmetry -> <ij|kl> = <kl|ij>"
            PermSym=PermSym+1
        else
            write(6,*) "No bra-ket permutational symmetry"
        endif
        if(tElecPerm) then
            write(6,"(A)") "Electron exchange symmetry -> <ij|kl> = <ji|lk>"
            PermSym=PermSym+1
        else
            write(6,*) "No electron exchange symmetry"
        endif
        if(tPermElecOrbs1) then
            write(6,"(A)") "Permutation of orbitals within either electron allowed -> <ij|kl> = <kj|il> = <il|kj>"
            PermSym=PermSym+2
        else
            write(6,*) "No permutation of orbitals within either electron allowed"
        endif
        if(PermSym.eq.0) then
            write(6,"(A)") "No permutational symmetry of the four-index, two-geminal pair antisymetrized spin integrals allowed"
        else
            write(6,"(A,I4)") "Permutational symmetry of the four-index, two-geminal pair antisymmetrized spin integrals is order: ",PermSym*2
        endif
        write(6,*) ""

    end subroutine CalculatePermSym


    subroutine testperm(TestTensor,Limits,tFail)
        integer :: i,j,k,l,Limits(2,4)
        real(8) :: TestTensor(Limits(1,1):Limits(2,1),Limits(1,2):Limits(2,2),Limits(1,3):Limits(2,3),Limits(1,4):Limits(2,4))
        logical, intent(in) :: tFail
        logical :: tPerm1,tPerm2,tPerm3

        tPerm1=.true.
        tPerm2=.true.
        tPerm3=.true.
            
        do i=1,nOrbsTot*2
            do j=1,nOrbsTot*2
                do k=1,nOrbsTot*2
                    do l=1,nOrbsTot*2
                        if(i.lt.Limits(1,1)) cycle
                        if(j.lt.Limits(1,2)) cycle
                        if(k.lt.Limits(1,3)) cycle
                        if(l.lt.Limits(1,4)) cycle
                        if(i.gt.Limits(2,1)) cycle
                        if(j.gt.Limits(2,2)) cycle
                        if(k.gt.Limits(2,3)) cycle
                        if(l.gt.Limits(2,4)) cycle
                        !<ij|kl> = <lk|ji>
                        if(abs(TestTensor(i,j,k,l)-TestTensor(l,k,j,i)).gt.1.D-8) then
                            if(tFail) then
                                write(6,*) i,j,k,l
                                write(6,*) TestTensor(i,j,k,l),TestTensor(l,k,j,i)
                                STOP 'Perm sym does not apply -> ijkl - lkji'
                            else
                                if(tPerm1) then
                                    write(6,*) i,j,k,l
                                    write(6,*) TestTensor(i,j,k,l),TestTensor(l,k,j,i)
                                    write(6,*) 'Perm sym does not apply -> ijkl - lkji'
                                endif
                                tPerm1=.false.
                            endif
                        endif
                        !<ij|kl> = <ji|lk>
                        if(abs(TestTensor(i,j,k,l)-TestTensor(j,i,l,k)).gt.1.D-8) then
                            if(tFail) then
                                write(6,*) i,j,k,l
                                write(6,*) TestTensor(i,j,k,l),TestTensor(j,i,l,k)
                                STOP 'Perm sym does not apply -> ijkl - jilk'
                            else
                                if(tPerm2) then
                                    write(6,*) i,j,k,l
                                    write(6,*) TestTensor(i,j,k,l),TestTensor(j,i,l,k)
                                    write(6,*) 'Perm sym does not apply -> ijkl - jilk'
                                endif
                                tPerm2=.false.
                            endif
                        endif
                        !<ij|kl> = <kl|ij>
                        if(abs(TestTensor(i,j,k,l)-TestTensor(k,l,i,j)).gt.1.D-8) then
                            if(tFail) then
                                write(6,*) i,j,k,l
                                write(6,*) TestTensor(i,j,k,l),TestTensor(k,l,i,j)
                                STOP 'Perm sym does not apply -> ijkl - klij'
                            else
                                if(tPerm3) then
                                    write(6,*) i,j,k,l
                                    write(6,*) TestTensor(i,j,k,l),TestTensor(k,l,i,j)
                                    write(6,*) 'Perm sym does not apply -> ijkl - klij'
                                endif
                                tPerm3=.false.
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo

        write(6,*) "PERMUTATIONAL SYMMETRY APPLIES!!!"
        call flush(6)

    end subroutine testperm


    function get_free_unit() result(free_unit)

        ! Returns:
        !    The first free file unit above 10 and less than or equal to
        !    the paramater max_unit (currently set to 200).

        integer, parameter :: max_unit = 100
        integer :: free_unit
        integer :: i
        logical :: t_open, t_exist

        do i = 10, max_unit
            inquire(unit=i, opened=t_open, exist=t_exist)
            if (.not.t_open .and. t_exist) then
                free_unit = i
                exit
            end if
        end do
        if (i == max_unit+1) stop 'Cannot find a free unit below max_unit.'

    end function get_free_unit
    
    subroutine CalcXTermEnergy(XTerm,EXSameSpinLoc,EXOppSpinLoc,TermNo)
        implicit none
        real(8), intent(in) :: XTerm(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2)
        real(8), intent(out) :: EXSameSpinLoc,EXOppSpinLoc
        integer, optional, intent(in) :: TermNo
        real(8) :: temp, temp2

        if(present(TermNo)) then
            write(6,"(A,I4)") "Calculating individual contribution to X intermediate energy from term: ",TermNo
        endif

        !All Xterms the same.
        EXSameSpinLoc=0.D0
        EXOppSpinLoc=0.D0
        do x=1,nOrbBasis*2,2
            do y=1,nOrbBasis*2,2
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
                        EXSameSpinLoc=EXSameSpinLoc + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
                                                      XTerm(r,s,x,y)*Phi(y,x,r,s) - &
                                                      XTerm(r,s,x,y)*Phi(x,y,s,r) + &
                                                      XTerm(r,s,x,y)*Phi(y,x,s,r)
                    enddo
                enddo
            enddo
        enddo
        do x=2,nOrbBasis*2,2
            do y=2,nOrbBasis*2,2
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
                        EXSameSpinLoc=EXSameSpinLoc + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
                                                      XTerm(r,s,x,y)*Phi(y,x,r,s) - &
                                                      XTerm(r,s,x,y)*Phi(x,y,s,r) + &
                                                      XTerm(r,s,x,y)*Phi(y,x,s,r)
                    enddo
                enddo
            enddo
        enddo
        EXSameSpinLoc=EXSameSpinLoc/16.D0  !We divide by 16 here!
        write(6,*) "Same spin: ",EXSameSpinLoc

        do x=1,nOrbBasis*2,2
            do y=2,nOrbBasis*2,2
!                temp2=0.D0
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
!                        temp=0.D0
                        EXOppSpinLoc=EXOppSpinLoc + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
                                                    XTerm(r,s,x,y)*Phi(y,x,r,s) - &
                                                    XTerm(r,s,x,y)*Phi(x,y,s,r) + &
                                                    XTerm(r,s,x,y)*Phi(y,x,s,r)
!                        temp=temp + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
!                                                    XTerm(r,s,x,y)*Phi(y,x,r,s) - &
!                                                    XTerm(r,s,x,y)*Phi(x,y,s,r) + &
!                                                    XTerm(r,s,x,y)*Phi(y,x,s,r)
!                        temp2=temp2 + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
!                                                    XTerm(r,s,x,y)*Phi(y,x,r,s) - &
!                                                    XTerm(r,s,x,y)*Phi(x,y,s,r) + &
!                                                    XTerm(r,s,x,y)*Phi(y,x,s,r)
!                        if(abs(temp).gt.1.D-8) then
!                        write(6,*) x,y,temp/8.D0, XTerm(r,s,x,y),Phi(x,y,r,s),-XTerm(r,s,x,y),Phi(y,x,r,s),-XTerm(r,s,x,y),Phi(x,y,s,r),XTerm(r,s,x,y),Phi(y,x,s,r)
!                        endif
                    enddo
                enddo
!                write(6,*) x,y,temp2/8.D0
            enddo
        enddo
!        write(6,*) "**********************************************************"

        do x=2,nOrbBasis*2,2
            do y=1,nOrbBasis*2,2
!                temp2=0.D0
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
!                        temp=0.D0
                        EXOppSpinLoc=EXOppSpinLoc + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
                                                    XTerm(r,s,x,y)*Phi(y,x,r,s) - &
                                                    XTerm(r,s,x,y)*Phi(x,y,s,r) + &
                                                    XTerm(r,s,x,y)*Phi(y,x,s,r)
!                        temp=temp + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
!                                                    XTerm(r,s,x,y)*Phi(y,x,r,s) - &
!                                                    XTerm(r,s,x,y)*Phi(x,y,s,r) + &
!                                                    XTerm(r,s,x,y)*Phi(y,x,s,r)
!                        temp2=temp2 + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
!                                                    XTerm(r,s,x,y)*Phi(y,x,r,s) - &
!                                                    XTerm(r,s,x,y)*Phi(x,y,s,r) + &
!                                                    XTerm(r,s,x,y)*Phi(y,x,s,r)
!                        if(abs(temp).gt.1.D-8) then
!                        write(6,*) x,y,temp/8.D0, XTerm(r,s,x,y),Phi(x,y,r,s),-XTerm(r,s,x,y),Phi(y,x,r,s),-XTerm(r,s,x,y),Phi(x,y,s,r),XTerm(r,s,x,y),Phi(y,x,s,r)
!                        endif
                    enddo
                enddo
!                write(6,*) x,y,temp2/8.D0
            enddo
        enddo
        EXOppSpinLoc=EXOppSpinLoc/16.D0  !We divide by 16 here!
        write(6,*) "Opp spin: ",EXOppSpinLoc
        write(6,*) "Total term energy: ",EXSameSpinLoc+EXOppSpinLoc


!*****************************************************************************************************
!The way its done eventually

!        EXSameSpinLoc=0.D0
!        EXOppSpinLoc=0.D0
!        do x=1,nOrbBasis*2,2
!            do y=1,nOrbBasis*2,2
!                do r=1,nOrbBasis*2
!                    do s=1,nOrbBasis*2
!                        EXSameSpinLoc=EXSameSpinLoc + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
!                                                      XTerm(r,s,x,y)*Phi(y,x,r,s) - &
!                                                      XTerm(s,r,x,y)*Phi(x,y,r,s) + &
!                                                      XTerm(s,r,x,y)*Phi(y,x,r,s)
!                    enddo
!                enddo
!            enddo
!        enddo
!        do x=2,nOrbBasis*2,2
!            do y=2,nOrbBasis*2,2
!                do r=1,nOrbBasis*2
!                    do s=1,nOrbBasis*2
!                        EXSameSpinLoc=EXSameSpinLoc + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
!                                                      XTerm(r,s,x,y)*Phi(y,x,r,s) - &
!                                                      XTerm(s,r,x,y)*Phi(x,y,r,s) + &
!                                                      XTerm(s,r,x,y)*Phi(y,x,r,s)
!                    enddo
!                enddo
!            enddo
!        enddo
!        EXSameSpinLoc=EXSameSpinLoc/16.D0  !We divide by 16 here!
!        write(6,*) "Same spin: ",EXSameSpinLoc
!
!        do x=1,nOrbBasis*2,2
!            do y=2,nOrbBasis*2,2
!                temp2=0.D0
!                do r=1,nOrbBasis*2
!                    do s=1,nOrbBasis*2
!                        temp=0.D0
!                        EXOppSpinLoc=EXOppSpinLoc + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
!                                                    XTerm(r,s,x,y)*Phi(y,x,r,s) - &
!                                                    XTerm(s,r,x,y)*Phi(x,y,r,s) + &
!                                                    XTerm(s,r,x,y)*Phi(y,x,r,s)
!                        temp=temp + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
!                                                    XTerm(r,s,x,y)*Phi(y,x,r,s) - &
!                                                    XTerm(s,r,x,y)*Phi(x,y,r,s) + &
!                                                    XTerm(s,r,x,y)*Phi(y,x,r,s)
!                        temp2=temp2 + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
!                                                    XTerm(r,s,x,y)*Phi(y,x,r,s) - &
!                                                    XTerm(s,r,x,y)*Phi(x,y,r,s) + &
!                                                    XTerm(s,r,x,y)*Phi(y,x,r,s)
!                        write(6,*) x,y,temp/8.D0, XTerm(r,s,x,y),Phi(x,y,r,s),-XTerm(r,s,x,y),Phi(y,x,r,s),-XTerm(s,r,x,y),Phi(x,y,r,s),XTerm(s,r,x,y),Phi(y,x,r,s)
!                    enddo
!                enddo
!                write(6,*) x,y,temp2/8.D0
!            enddo
!        enddo
!        do x=2,nOrbBasis*2,2
!            do y=1,nOrbBasis*2,2
!                do r=1,nOrbBasis*2
!                    do s=1,nOrbBasis*2
!                        EXOppSpinLoc=EXOppSpinLoc + XTerm(r,s,x,y)*Phi(x,y,r,s) - &
!                                                    XTerm(r,s,x,y)*Phi(y,x,r,s) - &
!                                                    XTerm(s,r,x,y)*Phi(x,y,r,s) + &
!                                                    XTerm(s,r,x,y)*Phi(y,x,r,s)
!                    enddo
!                enddo
!            enddo
!        enddo
!        EXOppSpinLoc=EXOppSpinLoc/16.D0  !We divide by 16 here!
!        write(6,*) "Opp spin: ",EXOppSpinLoc
!        write(6,*) "Total term energy: ",EXSameSpinLoc+EXOppSpinLoc

    end subroutine CalcXTermEnergy
    
    subroutine CalcBTermEnergy(BTerm,EBSameSpinLoc,EBOppSpinLoc,TermNo)
        implicit none
        real(8), intent(in) :: BTerm(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2)
        real(8), intent(out) :: EBSameSpinLoc,EBOppSpinLoc
        integer, optional, intent(in) :: TermNo

        if(present(TermNo)) then
            write(6,"(A,I4)") "Calculating individual contribution to B intermediate energy from term: ",TermNo
        endif

        EBSameSpinLoc=0.D0
        EBOppSpinLoc=0.D0

        do x=1,nOrbBasis*2,2
            do y=1,nOrbBasis*2,2
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
                        EBSameSpinLoc=EBSameSpinLoc + BTerm(r,s,x,y)*TwoRDM(x,y,r,s) - &
                                                      BTerm(r,s,x,y)*TwoRDM(x,y,s,r) - &
                                                      BTerm(r,s,y,x)*TwoRDM(x,y,r,s) + &
                                                      BTerm(r,s,y,x)*TwoRDM(x,y,s,r)
                    enddo
                enddo
            enddo
        enddo
        do x=2,nOrbBasis*2,2
            do y=2,nOrbBasis*2,2
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
                        EBSameSpinLoc=EBSameSpinLoc + BTerm(r,s,x,y)*TwoRDM(x,y,r,s) - &
                                                      BTerm(r,s,x,y)*TwoRDM(x,y,s,r) - &
                                                      BTerm(r,s,y,x)*TwoRDM(x,y,r,s) + &
                                                      BTerm(r,s,y,x)*TwoRDM(x,y,s,r)
                    enddo
                enddo
            enddo
        enddo
        EBSameSpinLoc=EBSameSpinLoc/16.D0  !We divide by 16 here!
        write(6,*) "Same spin: ",EBSameSpinLoc

        do x=1,nOrbBasis*2,2
            do y=2,nOrbBasis*2,2
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
                        EBOppSpinLoc=EBOppSpinLoc + BTerm(r,s,x,y)*TwoRDM(x,y,r,s) - &
                                                    BTerm(r,s,x,y)*TwoRDM(x,y,s,r) - &
                                                    BTerm(r,s,y,x)*TwoRDM(x,y,r,s) + &
                                                    BTerm(r,s,y,x)*TwoRDM(x,y,s,r)
                    enddo
                enddo
            enddo
        enddo
        do x=2,nOrbBasis*2,2
            do y=1,nOrbBasis*2,2
                do r=1,nOrbBasis*2
                    do s=1,nOrbBasis*2
                        EBOppSpinLoc=EBOppSpinLoc + BTerm(r,s,x,y)*TwoRDM(x,y,r,s) - &
                                                    BTerm(r,s,x,y)*TwoRDM(x,y,s,r) - &
                                                    BTerm(r,s,y,x)*TwoRDM(x,y,r,s) + &
                                                    BTerm(r,s,y,x)*TwoRDM(x,y,s,r)
                    enddo
                enddo
            enddo
        enddo
        EBOppSpinLoc=EBOppSpinLoc/16.D0  !We divide by 16 here!
        write(6,*) "Opp spin: ",EBOppSpinLoc
        write(6,*) "Total term energy: ",EBSameSpinLoc+EBOppSpinLoc

    end subroutine CalcBTermEnergy
    
    subroutine CalcVTermRDMEnergy(VTerm,EVSameSpinLoc,EVOppSpinLoc)
        implicit none
        real(8), intent(in) :: VTerm(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2)
        real(8), intent(out) :: EVSameSpinLoc,EVOppSpinLoc
        real(8), allocatable :: Temp(:,:,:,:)

        allocate(Temp(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2))
        Temp=0.D0
        
        do p=1,nOrbBasis*2
            do q=1,nOrbBasis*2
                do x=1,nOrbBasis*2
                    do y=1,nOrbBasis*2
                        do r=1,nOrbBasis*2      !Eventually, this wants to be over all orbital basis...
                            do s=1,nOrbBasis*2
                                Temp(p,q,x,y)=Temp(p,q,x,y)+(VTerm(r,s,x,y)*TwoRDM(p,q,r,s))/2.D0
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

!    write(6,*) "Decomposing into spin contributions gives: "
        EVSameSpinLoc=0.D0
        EVOppSpinLoc=0.D0

        do p=1,nOrbBasis*2,2
            do q=1,nOrbBasis*2,2
                EVSameSpinLoc=EVSameSpinLoc+Temp(p,q,p,q)-Temp(p,q,q,p)
!                write(6,*) p,q,Temp(p,q,p,q)-Temp(p,q,q,p)
            enddo
        enddo
        do p=2,nOrbBasis*2,2
            do q=2,nOrbBasis*2,2
                EVSameSpinLoc=EVSameSpinLoc+Temp(p,q,p,q)-Temp(p,q,q,p)
!                write(6,*) p,q,Temp(p,q,p,q)-Temp(p,q,q,p)
            enddo
        enddo
        EVSameSpinLoc=EVSameSpinLoc/2.D0  !We divide by two here!

        write(6,*) "Same spin: ",EVSameSpinLoc
        do p=1,nOrbBasis*2,2
            do q=2,nOrbBasis*2,2
                EVOppSpinLoc=EVOppSpinLoc+Temp(p,q,p,q)-Temp(p,q,q,p)
!                write(6,*) p,q,Temp(p,q,p,q)-Temp(p,q,q,p)
            enddo
        enddo
        do p=2,nOrbBasis*2,2
            do q=1,nOrbBasis*2,2
                EVOppSpinLoc=EVOppSpinLoc+Temp(p,q,p,q)-Temp(p,q,q,p)
!                write(6,*) p,q,Temp(p,q,p,q)-Temp(p,q,q,p)
            enddo
        enddo
        EVOppSpinLoc=EVOppSpinLoc/2.D0

        write(6,*) "Opp spin: ",EVOppSpinLoc
        write(6,*) "Total term energy: ",EVSameSpinLoc+EVOppSpinLoc
        deallocate(temp)
    end subroutine CalcVTermRDMEnergy
    
    subroutine CalcVTermEnergy(VTerm,EVSameSpinLoc,EVOppSpinLoc,TermNo)
        implicit none
        real(8), intent(in) :: VTerm(1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2,1:nOrbBasis*2)
        real(8), intent(out) :: EVSameSpinLoc,EVOppSpinLoc
        integer, intent(in), optional :: TermNo

        if(present(TermNo)) then
            write(6,*) "Calculating contribution to V from term: ",TermNo
        endif

!    write(6,*) "Decomposing into spin contributions gives: "
        EVSameSpinLoc=0.D0
        EVOppSpinLoc=0.D0

        do p=1,nOrbBasis*2,2
            do q=1,nOrbBasis*2,2
                EVSameSpinLoc=EVSameSpinLoc+VTerm(p,q,p,q)-VTerm(p,q,q,p)
            enddo
        enddo
        do p=2,nOrbBasis*2,2
            do q=2,nOrbBasis*2,2
                EVSameSpinLoc=EVSameSpinLoc+VTerm(p,q,p,q)-VTerm(p,q,q,p)
            enddo
        enddo
        EVSameSpinLoc=EVSameSpinLoc/2.D0  !We divide by two here!

        write(6,*) "Same spin: ",EVSameSpinLoc
        do p=1,nOrbBasis*2,2
            do q=2,nOrbBasis*2,2
                EVOppSpinLoc=EVOppSpinLoc+VTerm(p,q,p,q)-VTerm(p,q,q,p)
            enddo
        enddo
        do p=2,nOrbBasis*2,2
            do q=1,nOrbBasis*2,2
                EVOppSpinLoc=EVOppSpinLoc+VTerm(p,q,p,q)-VTerm(p,q,q,p)
            enddo
        enddo
        EVOppSpinLoc=EVOppSpinLoc/2.D0

        write(6,*) "Opp spin: ",EVOppSpinLoc
        write(6,*) "Total term energy: ",EVSameSpinLoc+EVOppSpinLoc
    end subroutine CalcVTermEnergy
        
    !If tReadRDMs is true, then we have already read in the RDMs, and now need to calculate
    !generalised fock matrices, f+k and k
    subroutine CalcGenFockMats(fock,TwoIndLims)
        implicit none
        integer, intent(in) :: TwoIndLims(2)
        real*8, intent(out) :: fock(TwoIndLims(1):TwoIndLims(2),TwoIndLims(1):TwoIndLims(2))
        integer :: k,l,p,q,Limits(2,4),ierr,unit_g
        real*8 , allocatable :: gints(:,:,:,:)

        write(6,*) "Calculating generalised fock, fock+exchange and exchange matrices..."

        fock(:,:)=0.D0
        do k=1,norbstot*2
            do l=1,norbstot*2
                fock(k,l)=tmat(k,l)
                do p=1,norbbasis*2
                    do q=1,norbbasis*2
                        fock(k,l)=fock(k,l)+antisym_g(k,p,l,q)*OneRDM(q,p)
                    enddo
                enddo
            enddo
        enddo

        allocate(fpk(1:nOrbsTot*2,1:nOrbsTot*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        allocate(mo_k(1:nOrbsTot*2,1:nOrbsTot*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        fpk(:,:)=0.D0
        mo_k(:,:)=0.D0

        unit_g=get_free_unit()
        if(tformattedints) then
            open(unit_g,file='FCIDUMP',status='old',form='formatted')
        else
            open(unit_g,file='MO_G',status='old',form='unformatted',access='sequential',action='read')
        endif

        allocate(gints(1:nOrbsTot,1:nOrbsTot,1:nOrbsTot,1:nOrbBasis),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'

        Limits(1,1)=1
        Limits(2,1)=nOrbsTot
        Limits(1,2)=1
        Limits(2,2)=nOrbsTot
        Limits(1,3)=1
        Limits(2,3)=nOrbsTot
        Limits(1,4)=1
        Limits(2,4)=nOrbBasis

        call fill4indarr(unit_g,gints,Limits)

        do k=1,norbstot*2
            do l=1,norbstot*2
                fpk(k,l)=tmat(k,l)  !This will be zero if spins aren't the same
                do p=1,norbbasis*2
                    do q=1,norbbasis*2
                        if((mod(k,2).eq.mod(p,2)).and.(mod(l,2).eq.mod(q,2)).and.(mod(k,2).eq.mod(l,2))) then
                            ! aaaa or bbbb - all spins the same.
                            if(mod(k,2).eq.0) then
                                fpk(k,l) = fpk(k,l) + (gints(k/2,p/2,l/2,q/2) * OneRDM(q,p))
!                                mo_k(k,l) = mo_k(k,l) + (gints(k/2,p/2,q/2,l/2) * OneRDM(q,p))
                                mo_k(k,l) = mo_k(k,l) + (gints(k/2,l/2,q/2,p/2) * OneRDM(q,p))
                            else
                                fpk(k,l) = fpk(k,l) + (gints((k+1)/2,(p+1)/2,(l+1)/2,(q+1)/2) * OneRDM(q,p))
!                                mo_k(k,l) = mo_k(k,l) + (gints((k+1)/2,(p+1)/2,(q+1)/2,(l+1)/2) * OneRDM(q,p))
                                mo_k(k,l) = mo_k(k,l) + (gints((k+1)/2,(l+1)/2,(q+1)/2,(p+1)/2) * OneRDM(q,p))
                            endif
                        elseif((mod(k,2).eq.mod(l,2)).and.(mod(p,2).eq.mod(q,2)).and.(mod(k,2).ne.mod(p,2))) then
                            ! abab or baba
                            if(mod(k,2).eq.0) then
                                fpk(k,l) = fpk(k,l) + (gints(k/2,(p+1)/2,l/2,(q+1)/2) * OneRDM(q,p))
                            else
                                fpk(k,l) = fpk(k,l) + (gints((k+1)/2,p/2,(l+1)/2,q/2) * OneRDM(q,p))
                            endif
                        elseif((mod(k,2).eq.mod(q,2)).and.(mod(p,2).eq.mod(l,2)).and.(mod(k,2).ne.mod(p,2))) then
                            !abba or baab
                            if(mod(k,2).eq.0) then
!                                mo_k(k,l) = mo_k(k,l) + (gints(k/2,(p+1)/2,q/2,(l+1)/2) * OneRDM(q,p))
                                mo_k(k,l) = mo_k(k,l) + (gints(k/2,(l+1)/2,q/2,(p+1)/2) * OneRDM(q,p))
                            else
!                                mo_k(k,l) = mo_k(k,l) + (gints((k+1)/2,p/2,(q+1)/2,l/2) * OneRDM(q,p))
                                mo_k(k,l) = mo_k(k,l) + (gints((k+1)/2,l/2,(q+1)/2,p/2) * OneRDM(q,p))
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo

        deallocate(gints)
        close(unit_g)

    end subroutine CalcGenFockMats


    subroutine CalcTMat(fock,TwoIndLims)
        implicit none
        integer, intent(in) :: TwoIndLims(2)
        real*8, intent(inout) :: fock(TwoIndLims(1):TwoIndLims(2),TwoIndLims(1):TwoIndLims(2))
        integer :: i,n,m,unit_fock,ierr

        !Temporarily store the real fock matrix 
        unit_fock=get_free_unit()
        if(tformattedints) then
            open(unit_fock,file='MO_F',status='old',form='formatted')
        else
            open(unit_fock,file='MO_F',status='old',form='unformatted',access='sequential',action='read')
        endif

        write(6,"(A)") "Reading fock matrix from disk to calculate one-electron hamiltonian matrix elements..."
        call fill2indarr_spin(unit_fock,fock,TwoIndLims)
        close(unit_fock)

        allocate(TMAT(1:nOrbsTot*2,1:nOrbsTot*2),stat=ierr)
        if(ierr.ne.0) stop 'error allocating'
        tmat(:,:)=0.D0

        do n=1,nOrbsTot*2
            do m=1,nOrbsTot*2
                tmat(n,m) = fock(n,m)
                do i=1,HFOccOrbs*2
                    tmat(n,m) = tmat(n,m) - antisym_g(n,i,m,i)
                enddo
            enddo
        enddo
    end subroutine CalcTMat

END PROGRAM ReadMO




