module input_data

use const

implicit none

save

logical , parameter :: tFCIDUMP = .true.
logical , parameter :: tRDUMP = .false.
logical :: tHFRDMs    !false if reading in RDMs 

!Checks permutational symmetry of 2RDM.
!2RDM should be antisymmetric wrt to swapping of {ij} or {kl}
logical , parameter :: tCheckRDMPermSym = .true.
logical :: tIntegrateTwoRDM
logical , parameter :: tWritePairEnergies = .false.
!If this is true, calculate the orbital relaxation energies for both zeroth order hamiltonians
logical , parameter :: tZerothRelaxBoth = .true.
!Which form of the zeroth order hamiltonian to use for the MR orbital relaxation 
logical , parameter :: tDyall = .true.  

!This logical turns on optimisations in the code
!the F12 correction goes through different code, and only one permutation of each term is calculated.
!In addition, the abba energy is assumed to be the same as the abab energy, and not seperately calculated.
logical , parameter :: tCalcCombo = .true. 

!Mainly debugging option. Will calculate the troublesome non-truncating B term (rfr), but not include it
!in the energy - just write it out for comparison
!Must use HF density matrices (MP2F12) and tCalcCombo to calculate.
logical , parameter :: tExplicitCommTerm = .false.

!This logical indicates whether we are reading in from DALTON/4traf created files or not.
logical :: tDaltonFormat 
!If we are not reading in the dalton integrals, this flag will tell us whether the integrals are in binary format or not.
logical :: tReadBin

!This defines the format of the RDM files which are being read in, and whether they represent spatial or spin orbital notation
logical :: tSpatRDMs

!This indicates whether the calculation is on the uniform electron gas
logical :: tUEG

!This indicates whether the orbitals are complex or not.
!If they are complex, then <ij|kl> /= <il|kj>, and will have to be considered separately
logical :: tComplexOrbs
real(dp) :: gamma_length    !This is the gamma parameter in inverse bohr for the F12.
!This is only needed for STG, since there, the FTF term is a simple change from the FF term.

end module input_data
