!     
! File:   global_variables.F90
! Author: fellipe
!
! Created on May 20, 2015, 9:15 AM
!

MODULE global_variables
IMPLICIT NONE

CHARACTER(LEN=100) :: input_file, stress_file 
INTEGER :: num_poin_appro_rdf = 1, n_interval = 64
INTEGER :: natoms , nframes , number_of_bins_gr , number_of_bins_sq , status_open , correl_leng , numb_perio
INTEGER :: radial_distribution_flag , static_factor_flag , static_factor_option , Fs_flag , Fs_Tgap, Pc_flag, cluster_flag
INTEGER :: Fc_flag , Fc_Tgap , MSD_flag , MSD_Tgap , contacs_number_flag = 0 , create_contact_group_flag = 0
INTEGER :: MSD_by_contact_number_flag , Nc_by_contact_number_flag, Sq_Tgap , ave_number_Sq , gr_Tgap 
INTEGER :: nframes_stress , tp , flag_rheology , stress_acf_flag , stress_acf_Tgap, stress_acf_option
INTEGER :: G_primes_flag,  G_primes_option
INTEGER :: n_types
INTEGER :: n_molecules, n_parti_per_molecule, max_contacts
INTEGER, DIMENSION(:), ALLOCATABLE :: index_part, index
INTEGER, DIMENSION(:,:), ALLOCATABLE  :: Nc, particle_type
REAL(8) :: vetor_length , radius , cutoff , rho , sigma_LJ = 3.405 , delta_t , dump_interval , fase , stress_ampli
REAL(8) :: box_lower_limit , box_upper_limit , L , Fs_q , Fc_q  , Kb , Temp
REAL(8) :: delta_t_stress , dump_interval_stress , S , strain_ampli
REAL(8) :: dist_cluster
REAL(8), DIMENSION(:), ALLOCATABLE  ::  q , distance_rdf , S_qx_ave , S_qy_ave , S_qz_ave , S_q_ave , Ls 
REAL(8), DIMENSION(:), ALLOCATABLE  :: stress_fit , strain , time_stress,  correlation_stress , distance_rdf_fitted
REAL(8), DIMENSION(:), ALLOCATABLE  :: Fsx , Fsy , Fsz , MSD , Dt , S_qx , S_qy , S_qz , S_q , bond_length 
REAL(8), DIMENSION(:,:), ALLOCATABLE  :: rdf , x , rdf_smooth , rdf_fitted , stress , stress_acf_ave   
REAL(8), DIMENSION(:,:,:), ALLOCATABLE  :: vxyz_frames,rxyz, rx, ry, rz, rx2, ry2, rz2, rdf_types
REAL(8) , DIMENSION(:,:) , ALLOCATABLE :: Pc, group_size
REAL(8), DIMENSION(13) :: count_part
REAL(8), DIMENSION(:), ALLOCATABLE :: rdf_ave
REAL(8), DIMENSION(:), ALLOCATABLE :: mass
REAL(8), DIMENSION(4) :: mm  !max 4 types of particles 
REAL(8) :: total_mass 
INTEGER :: n_bins_angle
INTEGER, DIMENSION(:), ALLOCATABLE :: n_atoms_rdf
REAL(8), DIMENSION(:), ALLOCATABLE :: rho_rdf

!Density profile variables
INTEGER :: density_profile_flag , number_of_bins_dens_profile , dens_pro_Tgap , axis_dens

!Orientation profile variables
INTEGER :: orien_profile_flag , number_of_bins_orien_profile , orien_pro_Tgap , axis_orien, n_molecules_orien, &
           n_parti_per_molecule_orien, n_components

REAL(8) :: Lx, Ly, Lz

END MODULE global_variables
