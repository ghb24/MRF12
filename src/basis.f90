module basis

use basis_data
use input_data, only: tDaltonFormat,tUEG,tComplexOrbs,gamma_length,tReadBin
use input_data, only: tExplicitCommTerm 
use const

implicit none

contains

    subroutine init_basis()

        call check_format()

        if(tDaltonFormat) then
            !read basis for dalton/4traf created input
            call read_dalt_basis()
        else
            !Read basis for slater-type geminal basis from 
            !FCIDUMP formatted files
            call read_slat_basis()
        endif

    end subroutine init_basis

    !Routine to tell us whether to flip orbital contributions in the RDMs or not.
    subroutine flip_orbitals()
        use utils, only: get_free_unit
        use errors, only: stop_all
        implicit none
        logical :: exists
        integer :: haaunit
        real(dp), allocatable :: haa(:,:)
        logical, allocatable :: tCheckedOrb(:)
        integer :: i,j,ios
        real(dp) :: HaaVal
        character(len=*), parameter :: t_r="flip_orbitals"

        allocate(OrbPhaseFlip(tnmo))
        OrbPhaseFlip(:) = 1

        inquire(file='HAA',exist=exists)
        if(exists) then
            !This is a NOCI calculation. Determine orbitals which need to be flipped.

            haaunit = get_free_unit()
            open(unit=haaunit,file='HAA',status='old',form='formatted',action='read')
            rewind(haaunit)

            allocate(haa(tnmo,tnmo))    !Store one-electron integrals from qchem
            haa(:,:) = 0.0_dp
            allocate(tCheckedOrb(tnmo)) !Whether the phase of this orbital has been fixed or not.
            tCheckedOrb(:) = .false.

            do while(.true.)
                read(haaunit,*,iostat=ios) i,j,HaaVal
                if(ios.gt.0) call stop_all(t_r,"Error reading HAA")
                if(ios.lt.0) exit   !EOF
                haa(i,j) = HaaVal
                haa(j,i) = HaaVal
            enddo
            close(haaunit)

            do i=1,tnmo
                do j=1,tnmo
                    if(i.eq.j) cycle
                    if(abs(abs(tmat(i,j))-abs(haa(i,j))).gt.1.D-7) then
                        write(6,*) i,j,tmat(i,j),haa(i,j)
                        call stop_all(t_r,"Magnitude of one-electron integrals is different between qchem and dalton. Check reference")
                    endif
                enddo
            enddo

!            tCheckedOrb(1) = .true.
!            do i=1,tnmo
!                if(abs(tmat(i,1)).gt.1.D-7) then
!                    !i is coupled to 1
!
!                    tCheckedOrb(i)=.true.
!
!                    if(abs(tmat(i,1)-haa(i,1)).gt.1.D-6) then
!                        !Need to flip phase of i
!                        OrbPhaseFlip(i) = -1
!                    endif
!                endif
!            enddo
            !If tCheckedOrb(x) is true, then it is a consistent phase w.r.t. 1


            do i=1,tnmo
                do j=1,tnmo

                    if(abs(tmat(i,j)).gt.1.D-7) then
                        !Coupled
                    
                        if(abs(tmat(i,j)-(haa(i,j)*OrbPhaseFlip(i)*OrbPhaseFlip(j))).gt.1.D-7) then
                            !Need to flip an orbital sign
                            if(tCheckedOrb(i).and.tCheckedOrb(j)) then
                                !Signs wrong, but have already fixed signs of both orbitals
                                call stop_all(t_r,"Frustrated system")
                            endif

                            if(tCheckedOrb(i).eqv.tCheckedOrb(j)) then
                                !Neither orbital has been checked. We can flip either to use as reference
                                OrbPhaseFlip(i) = -1
                            else
                                if(tCheckedOrb(i)) then
                                    !We have checked i previously, therefore flip j
                                    OrbPhaseFlip(j) = -1
                                else
                                    OrbPhaseFlip(i) = -1
                                endif
                            endif
                        endif

                        tCheckedOrb(i) = .true.
                        tCheckedOrb(j) = .true.

                    endif
                enddo
            enddo

            !Sanity check at the end
            do i=1,tnmo
                do j=1,tnmo
                    if((abs(tmat(i,j)).gt.1.D-7).and.abs(tmat(i,j)-(haa(i,j)*OrbPhaseFlip(i)*OrbPhaseFlip(j))).gt.5.D-7) then
                        call stop_all(t_r,"Inconsistent orbital signs")
                    endif
                    if(.not.tCheckedOrb(i)) then
                        write(6,*) "Orbital: ",i
                        call stop_all(t_r,"Have not checked an orbital for sign consistency")
                    endif
                    if(.not.tCheckedOrb(j)) then
                        write(6,*) "Orbital: ",j
                        call stop_all(t_r,"Have not checked an orbital for sign consistency")
                    endif
                enddo
            enddo

            write(6,"(A)") "Reference orbitals correct, and rotated consistently"

            deallocate(haa,tCheckedOrb)

        endif

    end subroutine flip_orbitals


    !Routine to work out wheter to expect the 4traf-type integral files
    !or not.
    subroutine check_format()
        use errors, only: stop_all
        implicit none
        logical :: exists,exists2,exists3
        character(*), parameter :: t_r="check_format"

        !Search for the file FGDUMP to indicate whether we are 
        !reading in formatted STG functions or not.
        inquire(file='MO_FG',exist=exists)
        if(exists) then
            write(6,"(A)") "Calculating in GTG-basis from DALTON/4traf integral files..."
            tComplexOrbs = .false.
            tDaltonFormat = .true.
            write(6,"(A)") "Assuming that all orbitals and integrals are strictly real..."
        else
            write(6,"(A)") "Calculating in STG-basis from FCIDUMP-type formatted integral files..."
            tDaltonFormat = .false.
            inquire(file='FGDUMPBIN',exist=exists2)
            if(exists2) then
                write(6,"(A)") "Reading in integrals from binary files..."
                tReadBin=.true.
            else
                tReadBin=.false.
                inquire(file='FGDUMP',exist=exists3)
                if(.not.exists3) then
                    call stop_all(t_r,"Cannot find integral files")
                endif
                write(6,"(A)") "Reading in integrals from formatted files..."
            endif
        endif
    end subroutine check_format

!read basis information from the integral files from XXXDUMP files.
    subroutine read_slat_basis()
        use utils, only: get_free_unit
        use errors, only: stop_all, warning
        implicit none
        integer :: nOrb,nAux,nElec,Ms2,iSym,PropBitLen,nProp(3)
        integer :: gdumpunit,i,j,k,l,OrbNum,Pos,ios,ierr,KinMatUnit
        real(dp) :: Z,GAM,KinEneElem
        logical :: UHF,exists,UEG,tRemoveSym,KinExist
        integer(i8) :: OrbSym(1000)
        character(*), parameter :: t_r='read_slat_basis'
        namelist /FCI/ nOrb,nAux,nElec,Ms2,OrbSym,iSym,UHF,UEG,GAM,PropBitLen,nProp

        write(6,'(A)') "Reading in basis data from GDUMP header"
        write(6,'(A)') "No symmetry information will be used in this calculation"

        !defaults, in case they are not in the FCIDUMP header
        UHF = .false.
        iSym = 1
        nOrb = 0
        nAux = 0
        nElec = 0
        Ms2 = 0
        OrbSym(:) = 0
        PropBitLen = 0
        gam = 0.0_dp
        UEG = .false.

        if(tReadBin) then
            !Binary integral files. The header is written to a seperate formatted file GDUMPHEAD
            inquire(file='GDUMPHEAD',exist=exists)
            if(.not.exists) call stop_all(t_r,"No GDUMPHEAD file found")
        else
            inquire(file='GDUMP',exist=exists)
            if(.not.exists) call stop_all(t_r,"No GDUMP file found")
        endif

        gdumpunit = get_free_unit()
        if(tReadBin) then
            open(unit=gdumpunit,file='GDUMPHEAD',status='old',form='formatted',action='read')
        else
            open(unit=gdumpunit,file='GDUMP',status='old',form='formatted',action='read')
        endif
        rewind(gdumpunit)

        !Read FCIDUMP header
        read(gdumpunit,fci)

        if(UHF.or.(Ms2.ne.0).or.(iSym.ne.1)) call stop_all(t_r,"Cannot deal with open-shell systems currently")
        if(nElec.le.0) call stop_all(t_r,"No electron number read in")
        if(nOrb.le.0) call stop_all(t_r,"No basis number read in")
        if(gam.eq.0.0_dp) call stop_all(t_r,"No gamma lengthscale specified")

        !Store global basis variables (see basis_data for definitions)
        NEl = nElec
        gamma_length = gam
        tnocc = NEl/2
        tnmo = nOrb
        tnbf = nOrb !(I assume all linear dependancy issues are non-existant/sorted in integral program)
        tnxmo = nAux
        tnxbf = nAux
        tntbf = nOrb + nAux
        tntmo = nOrb + nAux
        if(UEG) then
            !UEG system
            tUEG = .true.
            tComplexOrbs = .true.
            write(6,"(A)") "UEG system detected."
            write(6,"(A)") "Complex orbitals, but real integrals assumed..."
        else
            !Solid-state system
            write(6,"(A)") "Real orbitals (and integrals) assumed. If this is not true, then tComplexOrbs should be set in the code"
            tUEG = .false.
        endif

        allocate(OrbSyms(tntmo))
        tRemoveSym=.false.
        do i=1,tntmo
            if(OrbSym(i).eq.0) then
                write(6,"(A,I6)") "Orbital: ",i
                call warning(t_r,"Undefined symmetry label for orbital - Turning off symmetry")
                tRemoveSym=.true.
                exit
            endif
            OrbSyms(i) = OrbSym(i)  !Transfer it to the global array
        enddo
        if(tRemoveSym) OrbSyms(:) = 1

        !Turn off symetry to start - we need to think of a way to put in
        !symmetry nicely (probably just copy neci, but will need to change
        !the way that daltons symmetry is specified.
        nsym = 1
        allocate(nfrz(1),nocc(1),nmo(1),nbf(1),nxmo(1),nxbf(1),ntbf(1),ntmo(1))

        !Initially, do not contract over frozen orbitals
        !This worries me slightly, since it means that now my CABS space is
        !not complete....
        !This may be ok, since the plane waves from the PAW core electrons
        !go over whole space, and we are just missing bessel functions in core.
        !Are we attempting to freeze orbitals
        inquire(file='FROZENORBS',exist=exists)
        if(exists) then
            write(6,"(A)") "FROZENORBS file detected..."
            write(6,"(A)") "Freezing orbitals will not be allowed in this calculation, and all electrons will be correlated"
        endif
        tnfrz = 0
        nfrz(1) = tnfrz 
        nocc(1) = tnocc
        nmo(1) = tnmo
        nbf(1) = tnbf
        nxmo(1) = tnxmo
        nxbf(1) = tnxbf
        ntbf(1) = tntbf
        ntmo(1) = tntmo

        if(tReadBin) then
            !If reading in binary integrals, close this file (which just contains the header), and reopen the full integral file.
            close(gdumpunit)
            inquire(file='GDUMPBIN',exist=exists)
            if(.not.exists) call stop_all(t_r,"Cannot find GDUMPBIN file")
            open(unit=gdumpunit,file='GDUMPBIN',status='old',form='unformatted',action='read')
        endif
        
        allocate(OrbEnergies(tntmo))
        OrbEnergies(:) = 0.0_dp
        !If reading in FCIDUMP-type integrals, the KE integrals are written out explicitly.
        !There is no need to necessarily calculate them later, apart from self-consistency.
        !We will read them in here.
        allocate(tmat(tntmo,tntmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc error")
        tmat(:,:) = 0.0_dp  !This should be diagonal for the UEG

        inquire(file='KinEne',exist=Kinexist)
        if(Kinexist) then
            !Load KE eigenvalues (hack for the UEG)
            KinMatUnit = get_free_unit() 
            open(unit=KinMatUnit,file='KinEne',status='old')
            rewind(KinMatUnit)
            allocate(Kinmat(tntmo,tntmo),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc error")
            Kinmat(:,:) = 0.0_dp  !This should be diagonal for the UEG

            do while(.true.)
                read(KinMatUnit,*,iostat=ios) KinEneElem,i,j
                Kinmat(i,j) = KinEneElem
                if(ios.gt.0) call stop_all(t_r,"Error read KinEne")
                if(ios.lt.0) exit
            enddo
        else
            if(tExplicitCommTerm) call stop_all(t_r,"Cannot calculate RI'd FTF term with KE matrix")
        endif

        !Pair spaces:
        tnmo_sq = tnmo*tnmo
        taux_nmo = tnmo*tnxmo
        tntmo_sq = tntmo*tntmo
        tntmo_nmo = tntmo*tnmo
        tntmo_nxmo = tntmo*tnxmo

        !Determine record lengths, assuming the record length is in chunks of 4 bytes
        reclen_nmo_sq = tnmo_sq*2
        reclen_ntmo_sq = tntmo_sq*2

        ref_ene = 0.0_dp    !Set reference energy to zero. This is not read in, and will have to be calculated later.
        pot_nuc = 0.0_dp

        !However, the pot_nuc energy can be read in.
        do while(.true.)
            if(tReadBin) then
                read(gdumpunit,iostat=ios) Z,i,k,j,l    !Chemical notation. Spatial orbitals.
            else
                read(gdumpunit,*,iostat=ios) Z,i,k,j,l    !Chemical notation. Spatial orbitals.
            endif
            if(ios.gt.0) call stop_all(t_r,"Error read GDUMP")
            if(ios.lt.0) exit

            if(i.eq.0) then
                !Core energy
                pot_nuc = Z
            elseif(k.eq.0) then
                !HF eigenvalues
                OrbEnergies(i) = Z
            elseif(j.eq.0) then
                !1e KE integrals 
                !For the UEG, these should be diagonal.
                tmat(i,k) = Z
                tmat(k,i) = Z
                if((i.ne.k).and.(abs(Z).gt.1.D-8).and.(tUEG)) then
                    call stop_all(t_r,"UEG basis is not a KE eigenbasis")
                endif
            endif
        enddo
        close(gdumpunit)
        
        if(pot_nuc.eq.0.0_dp) then
            write(6,"(A)") "No nuclear potential energy read in."
        else
            write(6,"(A,F20.10)") "Nuclear potential energy read in as: ",pot_nuc
        endif

!Check that the orbitals are ordered in increasing energy
        do i=2,tntmo
            if(OrbEnergies(i).lt.OrbEnergies(i-1)) call stop_all(t_r,"Assumed orbital ordering by energy is not correct")
        enddo

        allocate(OccOrbs(tnocc)) !This stores the indices of the occupied orbitals.
        allocate(IsOrbOcc(tnmo))    !This is a logical array indicating whether an orbital is occupied or not
        IsOrbOcc(:) = .false.

        !If we assume simple energy ordering of the orbitals, this is simple to calculate.
        do i=1,tnocc
            OccOrbs(i) = i
            IsOrbOcc(i) = .true.
        enddo

        !OccOrbNum will take the index of an occupied orbital, and tell you how many orbitals are occupied before it.
        !OccOrbs should be the inverse of this.
        !Basically, OccOrbNum is a mapping from the occupied orbitals, to a number 1 -> tnocc
        allocate(OccOrbNum(tnmo))
        OccOrbNum = 0
        OrbNum=1
        do i=1,tnmo
            if(IsOrbOcc(i)) then
                OccOrbNum(i)=OrbNum
                OrbNum=OrbNum+1
            endif
        enddo
        if(OrbNum.ne.tnocc+1) then
            call stop_all(t_r,"Error in setting up OccOrbNum")
        endif
        do i=1,tnmo
            if(IsOrbOcc(i)) then
                if(OccOrbs(OccOrbNum(i)).ne.i) then
                    call stop_all(t_r,"Error in setting up OccOrbNum")
                endif
            else
                if(OccOrbNum(i).ne.0) then
                    call stop_all(t_r,"Error in setting up OccOrbNum")
                endif
            endif
        enddo

        !Allocate arrays for frozen orbital data
        allocate(IsOrbFrz(tnmo))   !Logical to tell if frozen
        IsOrbFrz(:) = .false.
        allocate(FrzOrbs(tnfrz))    !Index of frozen orbitals
        FrzOrbs(:) = 0
        nfrz(:) = 0

        !Create orbital mapping for when reading in frozen density matrices
        !This means that orbital i from NECI corresponds to orbital OrbMapping(i) in this code.
        allocate(OrbMapping(tnmo))
        OrbMapping(:) = 0
        Pos = 0
        do i=1,tnmo
            if(.not.IsOrbFrz(i)) then
                Pos = Pos + 1
                OrbMapping(Pos) = i
            endif
        enddo

        !write out details of the orbital space
        write(6,'("HF occupied Orbs: ",8i4,"  total: ",i4)') nocc(1:nsym),(0, i=nsym+1,nsym),tnocc
        if(tnfrz.gt.0) then
            write(6,'("Frozen Orbs: ",8i4,"  total: ",i4)') nFrz(1:nsym),(0, i=nsym+1,nsym),tnfrz
        endif
        write(6,'("Orbital basis: ",8i4,"  total: ",i4)') nmo(1:nsym),(0, i=nsym+1,nsym),tnmo

            
        if(tnxmo.eq.0) then
            call warning(t_r,"Proceeding without auxiliary &
     &basis. This is highly ill-advised and will result in *very* slow orbital convergence") 
            call stop_all(t_r,"Cannot currently cope without an auxiliary basis set.")
        endif
        write(6,'("CABS space: ",8i4,"  total: ",i4)') nxmo(1:nsym),(0, i=nsym+1,nsym),tnxmo
        write(6,'("Total basis: ",8i4,"  total: ",i4)') ntmo(1:nsym),(0, i=nsym+1,nsym),tntmo
        if(tntmo.ne.(tnmo+tnxmo)) call stop_all(t_r,"Error in counting orbitals")
        write(6,*) ""
        write(6,"(A)") "Orb  Sym   Energy "
        do i=1,tnmo
            if(IsOrbFrz(i)) then
                write(6,"(I10,I5,F20.10,A)") i,OrbSyms(i),OrbEnergies(i),"  FROZEN"
            else
                write(6,"(I10,I5,F20.10)") i,OrbSyms(i),OrbEnergies(i)
            endif
        enddo
        write(6,*) ""

        !Create arrays indicating starting orbital of each symmetry
        allocate(StartSym_nmo(nsym))
        allocate(EndSym_nmo(nsym))
        StartSym_nmo(:) = 0
        EndSym_nmo(:) = 0
        StartSym_nmo(nsym) = tnmo - nmo(nsym) + 1
        do i=nsym-1,1,-1
            StartSym_nmo(i) = StartSym_nmo(i+1) - nmo(i)
        enddo
        if(StartSym_nmo(1).ne.1) then
            write(6,'("Starting indices: ",8I5)') StartSym_nmo(1:nsym)
            call stop_all(t_r,"Error setting up orb sym start arrays")
        endif
        do i=1,nsym
            EndSym_nmo(i) = StartSym_nmo(i) + nmo(i) - 1
        enddo
        if(EndSym_nmo(nsym).ne.tnmo) then
            write(6,'("Ending indices: ",8I5)') EndSym_nmo(1:nsym)
            call stop_all(t_r,"Error setting up orb sym end arrays")
        endif
        
        allocate(StartSym_nxmo(nsym))
        allocate(EndSym_nxmo(nsym))
        StartSym_nxmo(:) = 0
        EndSym_nxmo(:) = 0
        StartSym_nxmo(nsym) = tntmo - nxmo(nsym) + 1
        do i=nsym-1,1,-1
            StartSym_nxmo(i) = StartSym_nxmo(i+1) - nxmo(i)
        enddo
        if(StartSym_nxmo(1).ne.tnmo+1) then
            write(6,'("Starting indices CABS: ",8I5)') StartSym_nxmo(1:nsym)
            call stop_all(t_r,"Error setting up aux sym start arrays")
        endif
        do i=1,nsym
            EndSym_nxmo(i) = StartSym_nxmo(i) + nxmo(i) - 1
        enddo
        if(EndSym_nxmo(nsym).ne.tntmo) then
            write(6,'("Ending indices CABS:   ",8I5)') EndSym_nxmo(1:nsym)
            call stop_all(t_r,"Error setting up aux sym end arrays")
        endif

        allocate(AuxOrbSyms(tnmo+1:tntmo))
        AuxOrbSyms(:) = 0
        do j=1,nsym
            do i=StartSym_nxmo(j),EndSym_nxmo(j)
                AuxOrbSyms(i) = j
            enddo
        enddo

    end subroutine read_slat_basis

!read daltons basis information from SIRIFC
!note: for reasons I don't fully understand, read in integers from
!dalton must be declared as integer(i4)'s. These integers will have a _dalt
!appended to the variable name.
    subroutine read_dalt_basis()
        use dalt_funcs, only: find_label
        use utils, only: get_free_unit
        use errors, only: stop_all, warning
        implicit none
        integer :: sirifcunit, unit_err, i, idum, nmo_i, loop, j, nsym_temp, auxunit, Occ, OrbitalNum
        integer :: k, MinEIndex, unit_frz, OrbNum, j_OccOrb, MappedOrb, Pos, MO_unit
        logical :: sirexist,tmoeqao,auxexist,exists
        integer(i4) :: nsym_dalt,norbt_dalt,ncmot_dalt,nocc_dalt(8),nlambda_dalt(8),norb_dalt(8)
        integer(i4) :: nbast_dalt 
        real(dp), allocatable :: mocoeffs(:)
        character(*), parameter :: t_r='read_dalt_basis'

        !We want to look in the file SIRIFC
        sirifcunit = get_free_unit() 
        inquire(file='SIRIFC',exist=sirexist)
        if(.not.sirexist) then
            call stop_all(t_r,"Cannot find SIRIFC file")
        endif

        open(unit=sirifcunit,file='SIRIFC',status='old',form='unformatted',access='sequential')
        rewind(sirifcunit)

        !find relevant part in file.
        unit_err = 6  !write errors to stdout
        call find_label('TRCCINT ',sirifcunit,unit_err)

        !Get data on orbitals, as well as nuclear potential and reference energy (HF energy?)
        read(sirifcunit) nsym_dalt,norbt_dalt,nbast_dalt,ncmot_dalt,nocc_dalt(1:nsym_dalt),nlambda_dalt(1:nsym_dalt),  &
     &          norb_dalt(1:nsym_dalt),pot_nuc,ref_ene

!        write(6,*) nsym_dalt,norbt_dalt,ncmot_dalt,nocc_dalt(:),nlambda_dalt(:)

        !Transfer data into non-dalton format, allocate and put into the global variables
        nsym=nsym_dalt
        if(nsym.gt.8) call stop_all(t_r,"Error in calculating nsym")
        write(6,"(A,I4)") "Number of symmetry elements detected in basis: ",nsym
        write(6,"(A,F20.10)") "Nuclear potential energy read in as: ",pot_nuc
        write(6,"(A,F20.10)") "Reference energy (HF energy) read in as: ",ref_ene
        write(6,*) ""
        write(6,"(A)") "Reading in orbital basis set from SIRIFC, assuming it is canonical..."

        !Allocate nsym dependant arrays.
        allocate(nocc(nsym),nmo(nsym),nbf(nsym),nxmo(nsym),nxbf(nsym),ntbf(nsym),nfrz(nsym),ntmo(nsym))

        tmoeqao=.true.  !Logical to tell us if no. ao's = no. mo's
        do i=1,nsym
            nocc(i) = nocc_dalt(i)
            nmo(i) = nlambda_dalt(i)
            nbf(i) = norb_dalt(i)
            if(nmo(i).ne.nbf(i)) tmoeqao=.false.
        enddo

        !calculate the sum of orbitals
        tnocc = sum(nocc(1:nsym))
        tnmo = sum(nmo(1:nsym))
        tnbf = sum(nbf(1:nsym))

        !Assume that all occupied orbitals doubly occupied - this will need to change for open-shell
        NEl = tnocc*2

        tnmo_sq = tnmo*tnmo
        
        !Do we have frozen orbitals
        inquire(file='FROZENORBS',exist=exists)
        if(exists) then
            write(6,"(A)") "FROZENORBS file detected"
            unit_frz=get_free_unit()
            open(unit_frz,file='FROZENORBS',status='old',form='formatted',action='read')
            read(unit_frz,*) tnfrz
            write(6,"(A,I4,A)") "Freezing ",tnfrz," core orbitals..."
            close(unit_frz)
        else
            tnfrz = 0
        endif

        allocate(OccOrbs(tnocc)) !This stores the indices of the occupied orbitals.
        allocate(IsOrbOcc(tnmo))    !This is a logical array indicating whether an orbital is occupied or not
        IsOrbOcc(:)=.false.
        !This will have to be updated once we attempt to take into account open-shell systems.
        OrbitalNum = 0
        Occ = 1
        do i=1,nsym
            !Orbitals are organised in symmetry, then energy (from DALTON anyway...same as MOLPRO) 
            !Assume that lowest energy orbitals are occupied ones.
            if(i.ne.1) then
                !Add in orbital indices from previous symmetry block
                OrbitalNum = OrbitalNum + nmo(i-1)
            endif
            do j=1,nocc(i)
                OccOrbs(Occ)=OrbitalNum+j
                IsOrbOcc(OrbitalNum+j)=.true.
                Occ = Occ + 1
            enddo
        enddo

        !OccOrbNum will take the index of an occupied orbital, and tell you how many orbitals are occupied before it.
        !OccOrbs should be the inverse of this.
        !Basically, OccOrbNum is a mapping from the occupied orbitals, to a number 1 -> tnocc
        allocate(OccOrbNum(tnmo))
        OccOrbNum = 0
        OrbNum=1
        do i=1,tnmo
            if(IsOrbOcc(i)) then
                OccOrbNum(i)=OrbNum
                OrbNum=OrbNum+1
            endif
        enddo
        if(OrbNum.ne.tnocc+1) then
            call stop_all(t_r,"Error in setting up OccOrbNum")
        endif
        do i=1,tnmo
            if(IsOrbOcc(i)) then
                if(OccOrbs(OccOrbNum(i)).ne.i) then
                    call stop_all(t_r,"Error in setting up OccOrbNum")
                endif
            else
                if(OccOrbNum(i).ne.0) then
                    call stop_all(t_r,"Error in setting up OccOrbNum")
                endif
            endif
        enddo

!        do i=1,tnmo
!            write(6,*) OccOrbNum(i)
!        enddo
        
        !allocate memory for reading in of individual orbital energies and symmetries.
        if(norbt_dalt.ne.tnmo) then
            call stop_all(t_r,"Mismatch between number of orbitals in orbital basis")
        endif
        allocate(OrbEnergies(norbt_dalt))
        allocate(OrbSyms(norbt_dalt))

        read(sirifcunit) OrbEnergies(1:norbt_dalt),OrbSyms(1:norbt_dalt)

        !Allocate arrays for frozen orbital data
        allocate(IsOrbFrz(tnmo))   !Logical to tell if frozen
        IsOrbFrz(:) = .false.
        allocate(FrzOrbs(tnfrz))    !Index of frozen orbitals
        FrzOrbs(:) = 0
        nfrz(:) = 0

        do i=1,tnfrz
            !Find the index of each of the frozen orbs
            !Find minimum energy orbitals

            !Find initial MinEIndex
            !We want to find any as yet unfrozen orbital
            do j=1,tnocc
                j_OccOrb = OccOrbs(j)
                if(.not.IsOrbFrz(j_OccOrb)) then
                    MinEIndex=j_OccOrb
                    exit
                endif
            enddo
            if(j.eq.tnocc+1) call stop_all(t_r,"All orbitals already frozen!")

!            write(6,*) "Initial choice of MinEIndex: ",MinEIndex

            do j=1,tnocc
                j_OccOrb = OccOrbs(j)
                if(OrbEnergies(j_OccOrb).lt.OrbEnergies(MinEIndex)) then
!                    write(6,*) j_OccOrb, " lower energy than ",MinEIndex
                    !Lower energy, but have we already considered it?
                    do k=1,i-1
                        if(j_OccOrb.eq.FrzOrbs(k)) exit
                    enddo
                    if(k.eq.i) then
                        !We completed the loop without finding that the orbital had already been chosen
                        MinEIndex=j_OccOrb
                    endif
                endif
            enddo

            FrzOrbs(i) = MinEIndex
            IsOrbFrz(MinEIndex) = .true.
            nFrz(OrbSyms(MinEIndex)) = nFrz(OrbSyms(MinEIndex)) + 1
!            write(6,"(A,I6)") "Freezing orbital: ",MinEIndex
        enddo

        !Create orbital mapping for when reading in frozen density matrices
        !This means that orbital i from NECI corresponds to orbital OrbMapping(i) in this code.
        allocate(OrbMapping(tnmo))
        OrbMapping(:) = 0

        Pos = 0
        do i=1,tnmo

            if(.not.IsOrbFrz(i)) then
                Pos = Pos + 1
                OrbMapping(Pos) = i
            endif
        enddo
!        write(6,*) "Mapping : "
!        do i=1,tnmo
!            write(6,*) i,OrbMapping(i)
!        enddo

        !write out details of the orbital space
        write(6,'("HF occupied Orbs: ",8i4,"  total: ",i4)') nocc(1:nsym),(0, i=nsym+1,nsym),tnocc
        if(tnfrz.gt.0) then
            write(6,'("Frozen Orbs: ",8i4,"  total: ",i4)') nFrz(1:nsym),(0, i=nsym+1,nsym),tnfrz
        endif
        if(tmoeqao) then
            !We only consider the nmo basis, since we assume we are only interested in the MO orbitals 
            !(which should have already been transformed.)
            write(6,'("Orbital basis: ",8i4,"  total: ",i4)') nmo(1:nsym),(0, i=nsym+1,nsym),tnmo
        else
            write(6,'("AO Orbital basis: ",8i4,"  total: ",i4)') nbf(1:nsym),(0, i=nsym+1,nsym),tnbf
            write(6,'("MO Orbital basis: ",8i4,"  total: ",i4)') nmo(1:nsym),(0, i=nsym+1,nsym),tnmo
            write(6,"(A)") "WARNING: AO basis and MO basis differ in orbital symmetries."
            write(6,"(A)") "This must be due to linear dependencies in the AO basis."
            write(6,"(A)") "The AO basis will be ignored, but beware of bugs in code."
        endif

        write(6,*) ""
        write(6,"(A)") "Orb  Sym   Energy "
        do i=1,norbt_dalt
            if(IsOrbFrz(i)) then
                write(6,"(I10,I5,F20.10,A)") i,OrbSyms(i),OrbEnergies(i),"  FROZEN"
            else
                write(6,"(I10,I5,F20.10)") i,OrbSyms(i),OrbEnergies(i)
            endif
        enddo
        write(6,*) ""


!        read(sirifcunit)    !overread this record
!        write(6,*) "ncmot: ",ncmot_dalt
        allocate(mocoeffs(ncmot_dalt))
        mocoeffs(:) = 0.0_dp
        read(sirifcunit) mocoeffs(1:ncmot_dalt)      !This *would* read in the AO/MO transformation matrix, 
                                                !but I don't think we have any need for it. 
                                                !It can be read in in the future if needs be.

        MO_unit=get_free_unit()
        open(unit=MO_unit,file='MO_Coeffs')
        write(MO_unit,*) mocoeffs(1:ncmot_dalt)
        close(MO_unit)
        deallocate(mocoeffs)

        close(sirifcunit)

        !Now read in the auxiliary basis
        write(6,"(A)") "Reading in auxiliary CABS basis set from AUXBAS..."
        !We want to look in the file AUXBAS
        auxunit = get_free_unit() 
        inquire(file='AUXBAS',exist=auxexist)
        if(.not.auxexist) then
            call warning(t_r,"Cannot find AUXBAS file. Proceeding without auxiliary    &
     &          basis. This is highly ill-advised and will result in *very* slow orbital convergence") 

            !TODO: Set all auxiliary basis variables to zero
            !zero variables
            tnxbf = 0
            tnxmo = 0
            nxmo(:) = 0
            nxbf(:) = 0

            ntbf(1:nsym) = nbf(1:nsym)  !The total number of orbitals 
            tntbf = tnbf 

        else
            open(unit=auxunit,file='AUXBAS',status='old',form='formatted',access='sequential')
            rewind(auxunit)

            read(auxunit,"(I5)") nsym_temp  !Read in the number of symmetries again
            if(nsym_temp.ne.nsym) then
                call stop_all(t_r,"nsym inconsistent")
            endif

            !zero variables
            tnxbf = 0
            tntbf = 0
            tnxmo = 0
            tntmo = 0

            do i=1,nsym
                ! read in: symmetry number, number of orbitals in orbital basis in this sym, 
                ! no. orbs in CABS in this sym, no. orbs in full space in this sym
                read(auxunit,*) idum,nmo_i,nxmo(i),ntbf(i)

                !checks...
                if(idum.ne.i) call stop_all(t_r,"symmetry loop not consistent")
                if (nmo_i.ne.nmo(i)) then
                    if(nmo_i.eq.nbf(i)) then
                        write(6,*) "Number of orbitals matches number of AOs, but not MOs - is this ok??"
                    endif
                    call stop_all(t_r,"Orbitals removed from mo basis")
                endif

                !Calculate number of AO functions in CABS as difference between total number and number in mo basis.
                nxbf(i) = ntbf(i) - nbf(i)
                ntmo(i) = nmo(i) + nxmo(i)

                tnxbf = tnxbf + nxbf(i)     !Sum the total number of AO functions in CABS
                tnxmo = tnxmo + nxmo(i)     !Find the total number of MO functions in CABS
                tntbf = tntbf + ntbf(i)     !Find the total number of functions which make up the space.
                tntmo = tntmo + ntmo(i)     !Find the total number of MO functions which make up the space.

                !Loop over the CMO coefficients, since we are not reading these in now.
                loop=nxmo(i)*(2+ntbf(i)/4)
                if(mod(ntbf(i),4).ne.0) loop=loop+nxmo(i)
                do j=1,loop
                    read(auxunit,*)
                enddo

            enddo
            close(auxunit,status='keep')

        endif

        if(tnxmo.eq.tnxbf) then
            !We only consider the nmo basis, since we assume we are only interested in the MO orbitals 
            !(which should have already been transformed.)
            write(6,'("CABS space: ",8i4,"  total: ",i4)') nxmo(1:nsym),(0, i=nsym+1,nsym),tnxmo
            write(6,'("Total basis: ",8i4,"  total: ",i4)') ntmo(1:nsym),(0, i=nsym+1,nsym),tntmo
        else
            write(6,'("CABS AO space: ",8i4,"  total: ",i4)') nxbf(1:nsym),(0, i=nsym+1,nsym),tnxbf
            write(6,'("CABS MO space: ",8i4,"  total: ",i4)') nxmo(1:nsym),(0, i=nsym+1,nsym),tnxmo
            write(6,'("Total basis: ",8i4,"  total: ",i4)') ntmo(1:nsym),(0, i=nsym+1,nsym),tntmo
            write(6,"(A)") "WARNING: AO basis and MO CABS basis differ in number."
            write(6,"(A)") "This must be due to linear dependencies in the AO basis."
            write(6,"(A)") "The AO basis will be ignored, but beware of bugs in code."
        endif

        !TODO: WE ALSO WANT TO FIND OUT THE SYMMETRIES OF THE CABS ORBITALS!
        !write out auxiliary CABS basis...

        if(tntmo.ne.(tnmo+tnxmo)) call stop_all(t_r,"Error in counting orbitals")

        taux_nmo = tnmo*tnxmo
        tntmo_sq = tntmo*tntmo
        tntmo_nmo = tntmo*tnmo
        tntmo_nxmo = tntmo*tnxmo

        !Determine record lengths, assuming the record length is in chunks of 4 bytes
        reclen_nmo_sq = tnmo_sq*2
        reclen_ntmo_sq = tntmo_sq*2


        !Create arrays indicating starting orbital of each symmetry
        allocate(StartSym_nmo(nsym))
        allocate(EndSym_nmo(nsym))
        StartSym_nmo(:) = 0
        EndSym_nmo(:) = 0
        StartSym_nmo(nsym) = tnmo - nmo(nsym) + 1
        do i=nsym-1,1,-1
            StartSym_nmo(i) = StartSym_nmo(i+1) - nmo(i)
        enddo
        if(StartSym_nmo(1).ne.1) then
            write(6,'("Starting indices: ",8I5)') StartSym_nmo(1:nsym)
            call stop_all(t_r,"Error setting up orb sym start arrays")
        endif
        do i=1,nsym
            EndSym_nmo(i) = StartSym_nmo(i) + nmo(i) - 1
        enddo
        if(EndSym_nmo(nsym).ne.tnmo) then
            write(6,'("Ending indices: ",8I5)') EndSym_nmo(1:nsym)
            call stop_all(t_r,"Error setting up orb sym end arrays")
        endif
        
        allocate(StartSym_nxmo(nsym))
        allocate(EndSym_nxmo(nsym))
        StartSym_nxmo(:) = 0
        EndSym_nxmo(:) = 0
        StartSym_nxmo(nsym) = tntmo - nxmo(nsym) + 1
        do i=nsym-1,1,-1
            StartSym_nxmo(i) = StartSym_nxmo(i+1) - nxmo(i)
        enddo
        if(StartSym_nxmo(1).ne.tnmo+1) then
            write(6,'("Starting indices CABS: ",8I5)') StartSym_nxmo(1:nsym)
            call stop_all(t_r,"Error setting up aux sym start arrays")
        endif
        do i=1,nsym
            EndSym_nxmo(i) = StartSym_nxmo(i) + nxmo(i) - 1
        enddo
        if(EndSym_nxmo(nsym).ne.tntmo) then
            write(6,'("Ending indices CABS:   ",8I5)') EndSym_nxmo(1:nsym)
            call stop_all(t_r,"Error setting up aux sym end arrays")
        endif

        allocate(AuxOrbSyms(tnmo+1:tntmo))
        AuxOrbSyms(:) = 0
        do j=1,nsym
            do i=StartSym_nxmo(j),EndSym_nxmo(j)
                AuxOrbSyms(i) = j
            enddo
        enddo

!        write(6,'("Starting indices OBS : ",8I5)') StartSym_nmo(1:nsym)
!        write(6,'("Ending indices OBS   : ",8I5)') EndSym_nmo(1:nsym)
!        write(6,'("Starting indices CABS: ",8I5)') StartSym_nxmo(1:nsym)
!        write(6,'("Ending indices CABS  : ",8I5)') EndSym_nxmo(1:nsym)


    end subroutine read_dalt_basis
        
    pure integer function FourIndSym(i,j,k,l)
        implicit none
        integer, intent(in) :: i,j,k,l

        FourIndSym=ieor(ieor(ieor(OrbSyms(i)-1,OrbSyms(j)-1),OrbSyms(k)-1),OrbSyms(l)-1)

    end function FourIndSym

    pure integer function FindTriPairInd(i,j)
        !Finds unique pair index for two electrons
        implicit none
        integer, intent(in) :: i,j

        if(i.gt.j) then
            FindTriPairInd=(i*(i-1))/2+j
        else
            FindTriPairInd=(j*(j-1))/2+i
        endif
    end function FindTriPairInd
    
    !Index the pair (a_prime,p) with a_prime fast!
    !a_prime is over CABS space, and p over orbital
    pure integer function FindCabsOrbPairInd(a_prime,p)
        implicit none
        integer, intent(in) :: a_prime,p

        FindCabsOrbPairInd = tnxmo*(p-1) + (a_prime-tnmo)

    end function FindCabsOrbPairInd
                        

    !Always right hand index is 'fast' index
    pure integer function FindSqPairInd(i,j,nBas)
        !Finds unique square pair index for two orbitals
        implicit none
        integer, intent(in) :: i,j,nBas

        FindSqPairInd=((i-1)*nBas)+j

    end function FindSqPairInd
    
    !Gets an unique index for an integral <i,j|k,l>
    pure integer function UMatInd(i,j,k,l)
        integer, intent(in) :: i,j,k,l
        integer :: a,b

        !Combine indices I and K, ensuring I>K
        if(i.gt.k) then
            a=(i*(i-1))/2+k
        else
            a=(k*(k-1))/2+i
        endif

        !combine indices j and l, ensuring J>K
        if(j.gt.l) then
            b=(j*(j-1))/2+l
        else
            b=(l*(l-1))/2+j
        endif

        !combine (ik) and (jl) in a unique way  (k > l or if k = l then i > j)
        if(a.gt.b) then
            UmatInd=(a*(a-1))/2+b
        else
            UmatInd=(b*(b-1))/2+a
        endif

    end function UMatInd


end module basis
