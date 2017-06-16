/*******************************************************
*      Copyright (C) 1995, 1998 Greg Landrum
*
*  This file is part of yaehmop.
*
*   This is free software.
* 
*  Permission is granted to modify, or otherwise fold, spindle, and mutilate this
*    code provided all copyright notices are left intact.
*
*  This code may be distributed to your heart's content, in whatever form,
*    provided no fee is charged for the distribution, all copyright notices are
*    left intact, and the source is distributed (without fee) along with any
*    binaries to anyone who requests it.
*
*  There are, of course, no warranties at all on this program.
*
********************************************************************/

/********

  the prototypes for all functions

  Created by greg Landrum October 1995

********/

/***
  Edit History:
  
  11.04.98 gL:
   prototypes for fatal_bug and nonfatal_bug modified
  
  August '98: WG
    - intercell_COOP_check added

***/


/* this is stolen from steve gifford */
#ifndef PROTO
# if defined(_NO_PROTO) || defined(_alpha) || defined(MIPSEL)
#  define PROTO(x) ()
# else /* defined(_NO_PROTO) || defined(__alpha) || defined(MIPSEL) */
#  define PROTO(x) x
# endif /* defined(_NO_PROTO) || defined(__alpha) || defined(MIPSEL) */
#endif /* PROTO */

extern void main PROTO((int, char **));

extern real eval_COOP PROTO((COOP_type *,detail_type *,cell_type *,long,
			      avg_prop_info_type *,hermetian_matrix_type,\
			      K_orb_ptr_type *,long *));
extern void gen_COOP PROTO((detail_type *,cell_type *,long,
			      avg_prop_info_type *,hermetian_matrix_type, K_orb_ptr_type *,
			      long *));
extern void gen_avg_COOPs PROTO((detail_type *,cell_type *,long,
			      avg_prop_info_type *,hermetian_matrix_type, K_orb_ptr_type *,
			      long *));
extern void intercell_COOP_check PROTO((COOP_type *));	
extern void gen_total_DOS PROTO((detail_type *,cell_type *,long,avg_prop_info_type *,
			   K_orb_ptr_type *));
extern real orb_contribution PROTO((long,long,long,avg_prop_info_type *,long));
extern real FMO_contribution PROTO((long,long,long,avg_prop_info_type *,long));
extern real atom_contribution PROTO((long,long,long,long,avg_prop_info_type *,long,long *));
extern void gen_projected_DOS PROTO((detail_type *,cell_type *,long,avg_prop_info_type *,
				 K_orb_ptr_type *,long *));
extern void build_k_hamil_FAT PROTO((cell_type *,hermetian_matrix_type,hermetian_matrix_type,
				     hermetian_matrix_type,long));
extern void build_k_hamil_THIN PROTO((cell_type *,hermetian_matrix_type,hermetian_matrix_type,
				     hermetian_matrix_type,long));
extern void build_k_overlap_FAT PROTO((cell_type *,k_point_type *,hermetian_matrix_type,
				     hermetian_matrix_type,long));
extern void build_k_overlap_THIN PROTO((cell_type *,detail_type *,
					k_point_type *,hermetian_matrix_type,
					hermetian_matrix_type,long));
extern void build_all_K_overlaps PROTO((cell_type *,detail_type *,
					hermetian_matrix_type,
				        hermetian_matrix_type,
					long,long,long *));

extern void R_space_Hamiltonian PROTO((cell_type *,detail_type *,hermetian_matrix_type,
				       hermetian_matrix_type,long,long *));
extern void full_R_space_Hamiltonian PROTO((cell_type *,detail_type *,hermetian_matrix_type,
				       hermetian_matrix_type,long,long *,char));
extern void zero_overlaps PROTO((real *,detail_type *,long,long,char,long *));
extern void calc_R_overlap PROTO((real *,cell_type *,detail_type *,long,point_type,char,long *));
extern void R_space_overlap_matrix PROTO((cell_type *,detail_type *,hermetian_matrix_type,
					  long,long,long *,long));
extern long find_atom PROTO((atom_type *,long,long));
extern void eval_Zmat_locs PROTO((atom_type *,long,long,char));
extern void calc_avg_occups PROTO((detail_type *,cell_type *,long,K_orb_ptr_type *,
				   avg_prop_info_type *,prop_type *,real *));
extern void print_avg_occups PROTO((detail_type *,cell_type *,long,K_orb_ptr_type *,
				   avg_prop_info_type *,prop_type,real *));
extern void calc_avg_FMO_occups PROTO((detail_type *,long,K_orb_ptr_type *,
				       avg_prop_info_type *,real *));
extern void calc_avg_OP PROTO((detail_type *,cell_type *,long,K_orb_ptr_type *,
			       avg_prop_info_type *,hermetian_matrix_type,
			       prop_type));
extern void find_crystal_occupations PROTO((detail_type *,real,long,K_orb_ptr_type *,real *));
extern void store_avg_prop_info PROTO((detail_type *,long,eigenset_type,hermetian_matrix_type,
				       long,real *,avg_prop_info_type *));
extern int sort_energies_helper PROTO((const void *,const void *));
extern void sort_avg_prop_info PROTO((detail_type *,long,avg_prop_info_type *,
				      K_orb_ptr_type *));
extern void gen_symm_lines PROTO((band_info_type *));
extern void construct_band_structure PROTO((cell_type *,detail_type *,
					    hermetian_matrix_type,hermetian_matrix_type,
					    hermetian_matrix_type,hermetian_matrix_type,
					    complex *,complex *,
					    eigenset_type,real *,real *,real *,complex *,
					    long, long *));
extern void eval_charge_matrix PROTO((cell_type *,eigenset_type,hermetian_matrix_type,
				      long,long *,real *,real *));
extern void check_a_cell PROTO((atom_type *,point_type,long,real,char *));
extern void check_nn_contacts PROTO((cell_type *,detail_type *details));
extern void build_distance_matrix PROTO((cell_type *,detail_type *details));
extern void dump_distance_mats PROTO((cell_type *,detail_type *details));
extern void fill_distance_matrix PROTO((cell_type *,long,float *,point_type *));
extern void display_lattice_parms PROTO((cell_type *));
extern long factorial PROTO((long));
extern void free_atomic_energy PROTO((cell_type *,real *,real *));
extern void AO_occupations PROTO((cell_type *,long,real *,long *,real *));
extern void eval_electrostatics PROTO((cell_type *, long,eigenset_type,real *,
				       real *,long *,real *,real *,real *,
				       real *,real *));
extern void read_geom_frag PROTO((FILE *,geom_frag_type *));
extern void write_atom_parms PROTO((detail_type *,atom_type *,long,char));
extern void write_atom_coords PROTO((atom_type *,long,char,char));
extern void fill_atomic_parms PROTO((atom_type *,long,FILE *));
extern void parse_printing_options PROTO((FILE *,detail_type *,cell_type *));
extern void read_inputfile PROTO((cell_type *,detail_type *,char *,long *,long **,FILE *));
extern void fatal PROTO((char *));
extern void fatal_bug PROTO((char *,char *,long));
extern void nonfatal_bug PROTO((char *,char *,long));
extern void error PROTO((char *));
extern void upcase PROTO((char *));
extern char *safe_strcpy PROTO((char *,char *));
extern void center_text_string PROTO((char *src, char *dest, long dest_len));
extern void left_just_text_string PROTO((char *src, char *dest, long dest_len));
extern void map_orb_num_to_name PROTO((char *,long,long *,long,atom_type *,long));
extern void debugmat PROTO((real *,long,long,real));
extern void printmat PROTO((real *,long,long,FILE *,real,char,long));
extern void print_labelled_mat PROTO(( real *, long, long, FILE *, real,
			atom_type *,long, long *,long,char,char,long));
extern void print_sym_mat PROTO((real *,long,long,FILE *,char *, char *,long));
extern long skipcomments PROTO((FILE *,char *,char));
extern void find_atoms_orbs PROTO((long,long,long,long *,long *,long *));
extern long overlap_tab_from_vect PROTO((point_type *,cell_type *));
extern void check_for_errors PROTO((cell_type *, detail_type *,long));
extern void build_orbital_lookup_table PROTO((cell_type *,long *,long **));
extern void print_MOs PROTO((detail_type *,long,eigenset_type,long,atom_type *,
			     long,long,long *));
extern double d_sign PROTO((double, double));
extern void parse_integer_string PROTO((char *,long **, long *));
extern void dump_hermetian_mat PROTO((long, real *, long));
extern void dump_sparse_mat PROTO((FILE *, real *, long, real));
extern void loop_over_k_points PROTO((cell_type *,detail_type *,
				      hermetian_matrix_type,hermetian_matrix_type,
				      hermetian_matrix_type,hermetian_matrix_type,
				      complex *, complex *,
				      eigenset_type,real *,real *,real *,complex *,
				      prop_type *,avg_prop_info_type *,long,long *));

extern void sparsify_hermetian_matrix PROTO((real, hermetian_matrix_type, long));
extern void sparsify_matrix PROTO((real, real *, real *, long));
extern void allocate_matrices PROTO((cell_type *,detail_type *,
				     hermetian_matrix_type *,hermetian_matrix_type *,
				     hermetian_matrix_type *,hermetian_matrix_type *,
				     complex **, complex **,
				     eigenset_type *,real **,real **,real **,
				     complex **,
				     prop_type *,avg_prop_info_type **,long long,
				     long *, long *, K_orb_ptr_type **));
extern void mov PROTO((real *,real *,real *,real *,
		       long,long,real,long,long,
		       long,long,atom_type *));
extern void calc_occupations PROTO((detail_type *,real,long,real *,
				    eigenset_type));
extern void reduced_mulliken PROTO((long,long,long *,real *,real *));
extern void FMO_reduced_mulliken PROTO((detail_type *,long,long,real *,real *));
extern void eval_mulliken PROTO((cell_type *,eigenset_type,hermetian_matrix_type,
				 long,real *,long *,real *,real *,real *));
extern void modified_mulliken PROTO((cell_type *,eigenset_type,hermetian_matrix_type,
				 long,real *,long *,real *,real *,real *,real *));
extern void read_NEW3file PROTO((cell_type *,detail_type *,FILE *));
extern void find_princ_axes PROTO((atom_type *,point_type *,real[3][3],real [3],
				   long));

extern void vector_diff PROTO((point_type *,point_type *,point_type *));
extern void normalize_vector PROTO((point_type *,point_type *));
extern real dot_prod PROTO((point_type *,point_type *));  
extern void scale_vector PROTO((point_type *,real));
extern void cross_prod PROTO((point_type *,point_type *,point_type *));
extern void mult_matrices PROTO((real *,real *,real *,long));
extern void auto_walsh PROTO((real *,long,real,real));
extern void walsh_update PROTO((cell_type *,detail_type *,long,char));
extern void walsh_output PROTO((detail_type *,cell_type *,long,eigenset_type,
				hermetian_matrix_type,hermetian_matrix_type,
				prop_type,long *,long));
extern void eval_xtal_coord_locs PROTO((cell_type *,char));
extern void update_zetas PROTO((cell_type *,real *,real,long *,char));
extern void init_FMO_file PROTO((detail_type *,long,real));
extern void build_FMO_overlap PROTO((detail_type *,long,long,hermetian_matrix_type,long *));
extern void build_FMO_hamil PROTO((detail_type *,long,long,hermetian_matrix_type,long *));
extern void diagonalize_FMO PROTO((detail_type *,real *,real *,real *,complex *,
				   complex *,complex *));
extern void gen_FMO_tform_matrices PROTO((detail_type *));
extern void tform_wavefuncs_to_FMO_basis PROTO((detail_type *,long,long,eigenset_type,long *));
extern void tform_matrix_to_FMO_basis PROTO((detail_type *,long,long,real *,real *,
					     real *, real *,complex_matrix_type,long *));
extern void tform_hermetian_matrix_to_FMO_basis PROTO((detail_type *,long,long,hermetian_matrix_type,
						       real *, real *,hermetian_matrix_type,long *));


extern void charge_to_num_electrons PROTO((cell_type *));
extern void update_chg_it_parms PROTO((detail_type *,cell_type *,real *,long *,long,long *));
extern void fill_chg_it_parms PROTO((atom_type *,long,long,FILE *));
extern void parse_charge_iteration PROTO((FILE *,detail_type *,cell_type *));
extern void update_muller_it_parms PROTO((detail_type *,cell_type *,
					  real *,long *,long,long *));
extern void calc_muller_init_parms PROTO((atom_type *));
extern void calc_muller_parms PROTO((atom_type *atom, real s_occup,real p_occup,real d_occup,
		       real *s_E,real *s_Z,real *p_E,real *p_Z,real *d_E,
		       real *d_Z));

extern void parse_muller_parms PROTO((FILE *,detail_type *,cell_type *));
extern void parse_equiv_atoms PROTO((FILE *,detail_type *,cell_type *));

extern void handle_sigint PROTO(());

/*********
  the fortran stuff
**********/
extern long abfns PROTO((real *,real *,real *,real *,real *,long *,long *,
			 long *,long *, long *, long *));
extern void lovlap PROTO((real *,real *,real *,real *,real *,real *,
			 long *,long *,long *,long *,long *,long *));			  
extern void cboris PROTO((long *,long *,real *,real *,real *,real *,
                          real *,real *,long *));
  
#ifndef SYM_OPS_DEFINED
#include "symmetry.h"
#endif
extern void name_sym_element PROTO((sym_op_type *,FILE *,long, long));
extern sym_op_type *make_new_sym_op PROTO((void));
extern void gen_sym_ops PROTO((sym_op_type **,long *));
extern void find_off_axis_sym_ops
  PROTO((detail_type *,cell_type *,point_type *,point_type *,
	 point_type [3],sym_op_type *,long *));
extern void find_sym_ops PROTO((detail_type *,cell_type *));
extern void find_walsh_sym_ops PROTO((cell_type *,detail_type *));
extern void find_MO_symmetries PROTO((long,detail_type *,cell_type *,
				      eigenset_type, hermetian_matrix_type,
				      long *));
extern void compare_molecules PROTO((atom_type *,point_type *,point_type *,
				     long,long *,char *,real));
extern void transform_p_orbs PROTO((real *,real [T_MAT_DIM][T_MAT_DIM]));
extern void transform_d_orbs PROTO((real *,real [D_T_MAT_DIM][D_T_MAT_DIM]));
extern void transform_orbitals PROTO((atom_type *,real *,sym_op_type *));
extern void full_transform PROTO((atom_type *,point_type,real [3][3],long));
extern void transform_atoms PROTO((atom_type *,real [T_MAT_DIM][T_MAT_DIM],long));
extern void transform_one_point
  PROTO((point_type *,real [T_MAT_DIM][T_MAT_DIM]));
extern void translate_atoms PROTO((point_type *, point_type, long));
extern void transform_atomic_locs PROTO((point_type *,real [T_MAT_DIM][T_MAT_DIM],long));
extern void transform_3x3_transpose PROTO((point_type *,real [3][3],long));
				

extern void postprocess_FMO PROTO((cell_type *,detail_type *,
				      hermetian_matrix_type,hermetian_matrix_type,
				      hermetian_matrix_type,hermetian_matrix_type,
				      complex *, complex *,
				      eigenset_type,real *,real *,real *,complex *,
				      prop_type *,avg_prop_info_type *,long,long *));

extern void postprocess_FCO PROTO((cell_type *,detail_type *,
				      hermetian_matrix_type,hermetian_matrix_type,
				      hermetian_matrix_type,hermetian_matrix_type,
				      complex *, complex *,
				      eigenset_type,real *,real *,real *,complex *,
				      prop_type *,avg_prop_info_type *,long,long *));

extern void postprocess_results PROTO((cell_type *,detail_type *,
				      hermetian_matrix_type,hermetian_matrix_type,
				      hermetian_matrix_type,hermetian_matrix_type,
				      complex *, complex *,
				      eigenset_type,real *,real *,real *,complex *,
				      prop_type *,avg_prop_info_type *,long,long *));
extern void process_geom_frags PROTO((cell_type *));
extern void insert_geom_frags PROTO((cell_type *));


extern void compare_crystal_basis PROTO((cell_type *,point_type *,point_type *,
					 point_type *,long,long *,char *,real));
extern void compare_crystal_lattice PROTO((cell_type *,point_type *,point_type *,
					   char *, real, real *));
extern void reduce_kpoints
  PROTO((detail_type *details,cell_type *cell,
	 point_type *raw_points,long num_raw_points,
	 real *multiplicity,long *num_points,long orthogonal_axes));

extern long check_for_orthogonal_basis PROTO((point_type vects[3],long dim,real tol));
extern long atoms_are_equiv  PROTO((cell_type *cell,point_type *loc1,
				   point_type *loc2,long *which_cell,
				   real symm_tol,point_type *));

extern void calc_reciprocal_lattice PROTO((cell_type *cell));
extern void gen_k_point_mesh PROTO((point_type **points,long num_per_vect[3],
				    long *num_generated, char include_high_symm_p,
				    long dim,real offset));
extern void automagic_k_points PROTO((detail_type *details,cell_type *cell));

extern long *my_malloc PROTO((long long));
extern long *my_calloc PROTO((long long, long));
extern long *my_realloc PROTO((long *, long));

#ifdef INCLUDE_NETCDF_SUPPORT
extern void netCDF_handle_error PROTO((long));
extern void netCDF_init_file PROTO((detail_type *details,cell_type *cell,long num_orbs,long));
extern void netCDF_close_file PROTO((detail_type *details));
extern void netCDF_write_Es PROTO((detail_type *details,long num_orbs,
				   avg_prop_info_type *prop_info));

#endif





#ifdef USE_LAPACK
#ifndef F2C_INCLUDE
typedef long int integer;
typedef double doublereal;
#endif
extern long zhegv_ PROTO((integer *itype, char *jobz, char *uplo, integer *n,
		 doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
		 doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork,
		 integer *info));
extern long zheev_ PROTO((char *jobz, char *uplo, integer *n, doublecomplex 
	*a, integer *lda, doublereal *w, doublecomplex *work, integer *lwork, 
	doublereal *rwork, integer *info));
#endif		 

#ifdef USING_THE_MAC
extern FILE *choose_mac_file(char *file_name, char);
#endif

