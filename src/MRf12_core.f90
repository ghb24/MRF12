program MRF12

use input_data
use const
use timing
use errors, only: stop_all
implicit none
real :: start_time
type(timer), save :: init_basis_timer
type(timer), save :: LoadFockOrb_timer
type(timer), save :: StoreGInts_timer
type(timer), save :: ReadOneRDM_timer
type(timer), save :: CalcOneHamInts_timer
type(timer), save :: CalcGenFock_timer
type(timer), save :: ReadTwoRDM_timer
type(timer), save :: CalcGenKandFpK_timer
type(timer), save :: CreateFCIDUMP_timer
type(timer), save :: CalcHFEnergy_timer 
type(timer), save :: CheckRDMProperties_timer 
type(timer), save :: CalcMP2Energy_timer 
type(timer), save :: StoreRInts_timer 
type(timer), save :: CreateR_TLD_timer 
type(timer), save :: CalcPhi_timer 
type(timer), save :: CreateR_TLD_Phi_timer 
type(timer), save :: CreateR_TLD_Both_timer
type(timer), save :: CalcF12Correction_timer 
type(timer), save :: Calc_SR_CABS_Singles_timer
type(timer), save :: Calc_MR_CABS_Singles_timer

call init_calc()

call run_calc()

call end_calc()

contains

    subroutine init_calc()
        use int_management
        use basis , only: init_basis,flip_orbitals
        use calc , only: CalcHFEnergy, CalcMP2Energy,CheckRDMProperties
        ! Initialise calculation.
        use report, only: environment_report
        logical :: exists,exists_2
        character(*), parameter :: t_r='init_calc'
        real(dp) :: DiskSpace

        write(6,"(A)") "***  Starting F12 calculation...  ***"

        call cpu_time(start_time) !record start time
        call environment_report() !write details of code

        call name_timers()

        !If want to read input in the future, do so here.
      
        call set_timer(init_basis_timer)
        call init_basis()
        call halt_timer(init_basis_timer)
        
        call CalcDiskReq(DiskSpace)
        write(6,'(F10.3,A)') DiskSpace," Gb total disk space approximately required to complete calculation."
        write(6,*) ""
        call flush(6)

        !read in fock matrix into seperate blocks as spatial orbitals.
        call set_timer(LoadFockOrb_timer)
        call LoadFockOrb()
        call halt_timer(LoadFockOrb_timer)

        call set_timer(StoreGInts_timer)
        call StoreGInts()
        call halt_timer(StoreGInts_timer)
        
        !calulate the one electron integrals
        call set_timer(CalcOneHamInts_timer)
        call CalcOneHamInts()
        call halt_timer(CalcOneHamInts_timer)
        
        !Routine to tell us whether to flip orbital contributions in the RDMs or not.
        call flip_orbitals()
        
        !Determine whether to use HF density matrix or not.
        !Also, determine whether to integrate 2RDM to get 1RDM,
        !and whether to read in the RDM in spin orbital notation,
        !or whether to read in the RDMs in spatial orbitals from
        !separate spin-type files.
        call DetermineRDMType()
        
        !If tHFRDMs is true, then it will just store the HF OneRDM
        call set_timer(ReadOneRDM_timer)
        call ReadOneRDM()
        call halt_timer(ReadOneRDM_timer)
        
        !Calculate the *generalised* fock matrices in seperate blocks
        call set_timer(CalcGenFock_timer)
        call CalcGenFock()
        call halt_timer(CalcGenFock_timer)
        
        call set_timer(ReadTwoRDM_timer)
        call ReadTwoRDM()
        call halt_timer(ReadTwoRDM_timer)

        !Calculate the *generalised* exchange and f+k matrices
        call set_timer(CalcGenKandFpK_timer)
        call CalcGenKandFpK()
        call halt_timer(CalcGenKandFpK_timer)

        if(tFCIDUMP) then
            write(6,"(A)",advance='no') "Dumping FCIDUMP G integrals... "
            call set_timer(CreateFCIDUMP_timer)
            call CreateFCIDUMP()
            call halt_timer(CreateFCIDUMP_timer)
            write(6,"(A)") "FCI Dumped"
        endif

        call set_timer(CalcHFEnergy_timer)
        call CalcHFEnergy()
        call halt_timer(CalcHFEnergy_timer)
        call set_timer(CheckRDMProperties_timer)
        call CheckRDMProperties()
        call halt_timer(CheckRDMProperties_timer)
        call set_timer(CalcMP2Energy_timer)
        call CalcMP2Energy()
        call halt_timer(CalcMP2Energy_timer)
        
    end subroutine init_calc

    subroutine run_calc()
        use int_management, only: StoreRInts,CreateR_TLD,CalcPhi
        use int_management, only: CreateR_TLD_Phi,CreateR_TLD_Both
        use calc, only : CalcF12Correction
        use calc_combo, only : CalcF12Correction_Combo
        use relax, only : Calc_SR_CABS_Singles,Calc_MR_CABS_Singles
        use basis_data, only: FockOrb,VarE,ref_ene,EMP2
        real(dp) :: Energy_F12,EF12,TotMP2,E_CABS_S
        character(*), parameter :: t_r='run_calc'

        call set_timer(StoreRInts_timer)
        call StoreRInts()
        call halt_timer(StoreRInts_timer)
        call set_timer(CalcPhi_timer)
        call CalcPhi() 
        call halt_timer(CalcPhi_timer)

        write(6,"(A)") "Creating R_TLD files..."
        call flush(6)
        call set_timer(CreateR_TLD_Both_timer)
        call CreateR_TLD_Both
        call halt_timer(CreateR_TLD_Both_timer)
        write(6,"(A)") "done."
        call flush(6)

!        call set_timer(CreateR_TLD_timer)
!        call CreateR_TLD() 
!        call halt_timer(CreateR_TLD_timer)
!        call set_timer(CreateR_TLD_Phi_timer)
!        call CreateR_TLD_Phi() 
!        call halt_timer(CreateR_TLD_Phi_timer)
        
        ! Run calculation according to input options.
        
        if(tHFRDMs) then
            call set_timer(Calc_SR_CABS_Singles_timer)
            call Calc_SR_CABS_Singles(E_CABS_S)
            call halt_timer(Calc_SR_CABS_Singles_timer)
        else
            call set_timer(Calc_SR_CABS_Singles_timer)
            call Calc_MR_CABS_Singles(E_CABS_S)
            call halt_timer(Calc_SR_CABS_Singles_timer)
        endif

        Energy_F12=0.D0
        call set_timer(CalcF12Correction_timer)
        if(tCalcCombo) then
            call CalcF12Correction_Combo(Energy_F12)
        else
            call CalcF12Correction(Energy_F12)
        endif
        call halt_timer(CalcF12Correction_timer)

        EF12 = Energy_F12
        TotMP2 = EMP2 + ref_ene

        write(6,"(A)") ""
        write(6,"(A)")          "**************************************************************"
        write(6,"(A)")          "****                                                      ****"
        write(6,"(A,F20.12,A)") "****   Hartree--Fock energy:       ",ref_ene,          "   ****"
        write(6,"(A)")          "****                                                      ****"
        write(6,"(A,F20.12,A)") "****   MP2 energy:                 ",EMP2   ,          "   ****"
        write(6,"(A)")          "****                                                      ****"
        write(6,"(A,F20.12,A)") "****   F12 correction:             ",EF12,             "   ****"
        write(6,"(A)")          "****                                                      ****"
        if(tHFRDMs) then
        write(6,"(A,F20.12,A)") "****   Orbital relaxation (CABS_S):",E_CABS_S,         "   ****"
        else
        if((.not.tDyall).or.tZerothRelaxBoth) then
        write(6,"(A,F20.12,A)") "****   Orbital relaxation (Fock):  ",E_CABS_S,         "   ****"
        else
        write(6,"(A,F20.12,A)") "****   Orbital relaxation (Dyall): ",E_CABS_S,         "   ****"
        endif
        endif
        write(6,"(A)")          "****                                                      ****"
        write(6,"(A,F20.12,A)") "****   Total CBS correction:       ",EF12+E_CABS_S,    "   ****"
        write(6,"(A)")          "****                                                      ****"
        if(tHFRDMs) then
        write(6,"(A,F20.12,A)") "****   MP2-F12 total energy:       ",TotMP2+EF12+E_CABS_S,"   ****"
        else
        write(6,"(A,F20.12,A)") "****   MR-F12 total energy:        ",VarE+EF12+E_CABS_S,"   ****"
        endif
        write(6,"(A)")          "****                                                      ****"
        write(6,"(A)")          "**************************************************************"

    end subroutine run_calc

    subroutine end_calc()

        ! Finalise calculation and clean up

        use report, only : end_report
        use utils, only: get_free_unit
        use basis_data, only : reclen_nmo_sq
        real :: end_time
        integer :: unit_2RDM_aaaa,unit_2RDM_abab,unit_2RDM_abba

        !Deallocate allocated memory (fock etc.)
!        call DeallocMem()

        if(.not.tHFRDMs) then
            !delete the new RDM files
            unit_2RDM_aaaa=get_free_unit()
            open(unit_2RDM_aaaa,file='TwoRDM_Dir_aaaa',status='old',form='unformatted',    &
                access='direct',recl=reclen_nmo_sq)
            close(unit_2RDM_aaaa,status='delete')
            unit_2RDM_abab=get_free_unit()
            open(unit_2RDM_abab,file='TwoRDM_Dir_abab',status='old',form='unformatted',    &
                access='direct',recl=reclen_nmo_sq)
            close(unit_2RDM_abab,status='delete')
            unit_2RDM_abba=get_free_unit()
            open(unit_2RDM_abba,file='TwoRDM_Dir_abba',status='old',form='unformatted',    &
                access='direct',recl=reclen_nmo_sq)
            close(unit_2RDM_abba,status='delete')
        endif

        call end_timing()
        call print_timing_report()
        call cpu_time(end_time)
        call end_report(end_time-start_time)

    end subroutine end_calc
        
!Determine whether to use HF density matrix or not.
!Also, determine whether to integrate 2RDM to get 1RDM,
!and whether to read in the RDM in spin orbital notation,
!or whether to read in the RDMs in spatial orbitals from
!separate spin-type files.
    subroutine DetermineRDMType()
        implicit none
        logical :: exists,exists_2,exists_3
        character(*), parameter :: t_r="DetermineRDMType"

        write(6,"(A)") ""

        inquire(file='TwoRDM_aaaa',exist=exists)
        if(exists) then
            !We are reading in spatial orbitals from separate spin files
            tSpatRDMs = .true.
            tHFRDMs = .false.
            inquire(file='TwoRDM',exist=exists_2)
            if(exists_2) then
                call stop_all(t_r,"Both TwoRDM_aaaa and TwoRDM files found. Unknown spin/spatial orbital format")
            endif
            write(6,"(A)") "RDM files detected in spatial orbital format..."

            inquire(file='OneRDM',exist=exists_3)
            if(.not.exists_3) then
                write(6,"(A)") "'TwoRDM' file detected, but not 'OneRDM'"
                write(6,"(A)") "OneRDM values will be calculated from integration of TwoRDM."
                tIntegrateTwoRDM=.true.
            else
                tIntegrateTwoRDM=.false.
                write(6,"(A)") "OneRDM file detected..."
            endif
        else
            inquire(file='TwoRDM',exist=exists_2)
            if(exists_2) then
                !We are reading in one file in spin-orbital notation
                write(6,"(A)") "RDM files detected in spin orbital format..."
                tSpatRDMs = .false.
                tHFRDMs = .false.

                inquire(file='OneRDM',exist=exists_3)
                if(.not.exists_3) then
                    write(6,"(A)") "'TwoRDM' file detected, but not 'OneRDM'"
                    write(6,"(A)") "OneRDM values will be calculated from integration of TwoRDM."
                    tIntegrateTwoRDM=.true.
                else
                    tIntegrateTwoRDM=.false.
                    write(6,"(A)") "OneRDM file detected..."
                endif
            else
                !There are no RDM files to read in
                write(6,"(A)") "No RDM files detected..."
                write(6,"(A)") "A Hartree--Fock density matrix will be assumed"
                tHFRDMs = .true.
            endif
        endif

    end subroutine DetermineRDMType

    !This routine calculates the disk space required to complete the calculation (in Gb)
    subroutine CalcDiskReq(Memory)
        use basis_data, only: tnmo_sq,tntmo_sq,tnmo,tnfrz,tntmo
        implicit none
        real(dp) , intent(out) :: Memory
        logical :: exists
        real(dp) :: tnmo_sq_loc,tntmo_sq_loc,tnmo_loc,tnfrz_loc,tntmo_loc

        tnmo_sq_loc = real(tnmo_sq,dp)
        tntmo_sq_loc = real(tntmo_sq,dp)
        tnmo_loc = real(tnmo,dp)
        tntmo_loc = real(tntmo,dp)
        tnfrz_loc = real(tnfrz,dp)

        write(6,'(A,F10.3)') 'G and R = ',tnmo_sq_loc*tntmo_sq_loc*2.0_dp*2.0_dp*8.0_dp*9.31322575D-10
        Memory = tnmo_sq_loc*tntmo_sq_loc*2.0_dp*2.0_dp*8.0_dp*9.31322575D-10   !This is the G and R direct files (two each)

        write(6,'(A,F10.3)') 'Phi files = ',tnmo_sq_loc*tnmo_sq_loc*3.0_dp*8.0_dp*9.31322575D-10
        Memory = Memory + tnmo_sq_loc*tnmo_sq_loc*3.0_dp*8.0_dp*9.31322575D-10  !This is for the three Phi files

        if(tCalcCombo) then
            write(6,'(A,F10.3)') 'rtld files = ',tnmo_sq_loc*tntmo_sq_loc*2.0_dp*8.0_dp*9.31322575D-10
            Memory = Memory + tnmo_sq_loc*tntmo_sq_loc*2.0_dp*8.0_dp*9.31322575D-10 !This is for the three rtld files
        else
            write(6,'(A,F10.3)') 'rtld files = ',tnmo_sq_loc*tntmo_sq_loc*3.0_dp*8.0_dp*9.31322575D-10
            Memory = Memory + tnmo_sq_loc*tntmo_sq_loc*3.0_dp*8.0_dp*9.31322575D-10 !This is for the three rtld files
        endif

        if(tCalcCombo) then
            write(6,'(A,F10.3)') 'phi_tld files = ',(tnmo_loc-tnfrz_loc)*(tnmo_loc-tnfrz_loc)*tntmo_loc*tntmo_loc*2.0_dp*8.0_dp*9.31322575D-10
            Memory = Memory + (tnmo_loc-tnfrz_loc)*(tnmo_loc-tnfrz_loc)*tntmo_loc*tntmo_loc*2.0_dp*8.0_dp*9.31322575D-10 !This is for the three phi_tld files
        else
            write(6,'(A,F10.3)') 'phi_tld files = ',(tnmo_loc-tnfrz_loc)*(tnmo_loc-tnfrz_loc)*tntmo_loc*tntmo_loc*3.0_dp*8.0_dp*9.31322575D-10
            Memory = Memory + (tnmo_loc-tnfrz_loc)*(tnmo_loc-tnfrz_loc)*tntmo_loc*tntmo_loc*3.0_dp*8.0_dp*9.31322575D-10 !This is for the three phi_tld files
        endif


        write(6,'(A,F10.3)') 'If 2RDM = ',tnmo_sq_loc*tnmo_sq_loc*3.0_dp*8.0_dp*9.31322575D-10
        inquire(file='TwoRDM',exist=exists)
        if(exists) then
            Memory = Memory + tnmo_sq_loc*tnmo_sq_loc*3.0_dp*8.0_dp*9.31322575D-10  !This is for the 2RDM direct access files.
        endif

    end subroutine CalcDiskReq
    
    subroutine name_timers()

        init_basis_timer%timer_name='init_basis'
        LoadFockOrb_timer%timer_name='LoadFockOrb'
        StoreGInts_timer%timer_name='StoreGInts'
        ReadOneRDM_timer%timer_name='ReadOneRDM'
        CalcOneHamInts_timer%timer_name='CalcOneHamInts'
        CalcGenFock_timer%timer_name='CalcGenFock'
        ReadTwoRDM_timer%timer_name='ReadTwoRDM'
        CalcGenKandFpK_timer%timer_name='CalcGenKandFpK'
        CreateFCIDUMP_timer%timer_name='CreateFCIDUMP'
        CalcHFEnergy_timer%timer_name='CalcHFEnergy'
        CheckRDMProperties_timer%timer_name='CheckRDMProperties'
        CalcMP2Energy_timer%timer_name='CalcMP2Energy'
        StoreRInts_timer%timer_name='StoreRInts'
        CreateR_TLD_timer%timer_name='CreateR_TLD'
        CalcPhi_timer%timer_name='CalcPhi'
        CreateR_TLD_Phi_timer%timer_name='CreateR_TLD_Phi'
        CalcF12Correction_timer%timer_name='CalcF12Correction'
        CreateR_TLD_Both_timer%timer_name='CreateR_TLD_Both'
        Calc_SR_CABS_Singles_timer%timer_name='Calc_SR_CABS_Singles'
        Calc_MR_CABS_Singles_timer%timer_name='Calc_SR_CABS_Singles'

    end subroutine name_timers

end program MRF12
