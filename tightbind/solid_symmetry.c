/*******************************************************
*      Copyright (C) 1997 Greg Landrum
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



/****************************************************************************
*
*     this file contains everything needed to deal with finding symmetry
*      and applying symmetry operations to crystals
*
*  created:  greg landrum  November 1997
*
*****************************************************************************/
#include "bind.h"
#include "symmetry.h"

/****************************************************************************
*
*                   Procedure   reduce_kpoints
*
* Arguments:  details: pointer to detail_type
*                cell: pointer to cell_type
*          raw_points: pointer to point_type
*      num_raw_points: int
*        multiplicity: pointer to real
*          num_points: pointer to int
*       nonortho_axes: int
*
* Returns: none
*
* Action: uses the symmetry operations of 'cell 
*   to reduce the number 
*   of k points in 'raw_points.  the multiplicity of each polong is returned
*   in 'multiplicity.  If a polong has an initial multiplicity of -1, it
*   will be ignored here.  
*  The final number of points is returned in 'num_points.
*
*  If 'nonortho_axes is non-zero then the points are transformed
*   into the cartesian basis before the comparisons are done.
*
*  NOTE:  This code is actually not correct in the general case where the
*   three lattice vectors do not polong along the cartesian axes.  Luckily, 
*   since the code which finds the symmetry elements is too stupid to find
*   any elements that don't lie along axes, this will work.
*   ah... the joys of stupid code.  :-)
*
****************************************************************************/
void reduce_kpoints(detail_type *details,cell_type *cell,
		    point_type *raw_points,long num_raw_points,
		    real *multiplicity,long *num_points,
		    long nonortho_axes)
{
  long p1,p2;
  long i,j;
  point_type *new_points,*orig_points;
  sym_op_type *this_op;

  new_points = (point_type *)my_calloc(num_raw_points,sizeof(point_type));
  if(!new_points) fatal("reduce_kpoints can't allocate new_points");

  /* initialize the multiplicities of nonredundant points to 1 */
  for(p1=0;p1<num_raw_points;p1++){
    if(multiplicity[p1] != -1) multiplicity[p1] = 1;
  }

  /* transform the points to cartesians if necessary */
  if( nonortho_axes ){
    orig_points = (point_type *)my_calloc(num_raw_points,sizeof(point_type));
    if(!new_points) fatal("reduce_orthogonal_kpoints can't allocate orig_points");
    bcopy((char *)raw_points,(char *)orig_points,
	  num_raw_points*sizeof(point_type));	  
    for(p1=0;p1<num_raw_points;p1++){
      if(multiplicity[p1] != -1){
	raw_points[p1].x = orig_points[p1].x * cell->recip_vects[0].x +
	   orig_points[p1].y * cell->recip_vects[1].x +
	     orig_points[p1].z * cell->recip_vects[2].x; 
	raw_points[p1].y = orig_points[p1].x * cell->recip_vects[0].y +
	   orig_points[p1].y * cell->recip_vects[1].y +
	     orig_points[p1].z * cell->recip_vects[2].y; 
	raw_points[p1].z = orig_points[p1].x * cell->recip_vects[0].z +
	   orig_points[p1].y * cell->recip_vects[1].z +
	     orig_points[p1].z * cell->recip_vects[2].z;
      }
    }
  }

  /* until I get off my butt and fix it, sym_ops_present is a global... ew */
  this_op = sym_ops_present;
  while(this_op){
    if( this_op->type != Identity && !this_op->redundant ){
      bcopy((char *)raw_points,(char *)new_points,
	    num_raw_points*sizeof(point_type));
      for(p1 = 0;p1<num_raw_points;p1++){
	if( multiplicity[p1] > 0 ){
	  if( this_op->type == Rotation || this_op->type == Improper_Rotation ){
	    /* do all the rotations */
	    for(i=1;i<this_op->order;i++){
	      transform_one_point(&(new_points[p1]),this_op->t_mat);
	      /* loop over the other points and see if we can find a match */
	      for(p2 = p1+1;p2<num_raw_points;p2++){
		if( multiplicity[p2] > 0 &&
		   POINTS_ARE_THE_SAME(&(new_points[p1]),&(raw_points[p2]),
				       details->symm_tol)){
		  multiplicity[p1] += multiplicity[p2];
		  multiplicity[p2] = 0;
		}
	      }
	    }
	  } else {
	    /* it's not a rotation, so we just need to do one transformation */
	    transform_one_point(&(new_points[p1]),this_op->t_mat);
	    /* loop over the other points and see if we can find a match */
	    for(p2 = p1+1;p2<num_raw_points;p2++){
	      if( multiplicity[p2] > 0 &&
		 POINTS_ARE_THE_SAME(&(new_points[p1]),&(raw_points[p2]),
				     details->symm_tol)){
		multiplicity[p1] += multiplicity[p2];
		multiplicity[p2] = 0;
	      }
	    }
	  }
	}
      }
    }
    this_op = this_op->next;
  }
  /* copy back in the original coordinates if needed */
  if( nonortho_axes ){
    bcopy((char *)orig_points,(char *)raw_points,
	  num_raw_points*sizeof(point_type));
  }
  /* run through, eliminate -k if k is present, and count the points */
  *num_points = 0;
  for(p1 = 0; p1<num_raw_points; p1++){
    if( multiplicity[p1] > 0 ){
      for(p2 = p1+1; p2<num_raw_points; p2++){
	if( multiplicity[p2] > 0 ){
	  new_points[0].x = -raw_points[p2].x;
	  new_points[0].y = -raw_points[p2].y;
	  new_points[0].z = -raw_points[p2].z;
	  if( POINTS_ARE_THE_SAME(&(raw_points[p1]),&(new_points[0]),
				  details->symm_tol)){
	    multiplicity[p1] += multiplicity[p2];
	    multiplicity[p2] = 0;
	  }
	}
      }
      (*num_points)++;
    }
  }
  free(new_points);
  if(nonortho_axes) free(orig_points);

  /* the last bit is to reduce the weights of points at zone edges/corners */
  if(details->use_high_symm_p){
    for(p1=0;p1<num_raw_points;p1++){
      if( multiplicity[p1] > 0 ){
	if( fabs(raw_points[p1].x) == 0.5 ) multiplicity[p1] /= 2.0;
	if( fabs(raw_points[p1].y) == 0.5 ) multiplicity[p1] /= 2.0;
	if( fabs(raw_points[p1].z) == 0.5 ) multiplicity[p1] /= 2.0;
      }
    }
  }
}

/****************************************************************************
*
*                   Function check_for_orthogonal_basis
*
* Arguments:  vects: point_type[3]
*               dim:  int
*               tol: real
*
* Returns: int
*
* Action: checks for an orthogonal basis.
*   the return value is determined by setting bits in the return
*   value for each set of non-orthogonal vectors as follows:
*    bit 0: 0 and 1
*    bit 1: 0 and 2
*    bit 2: 1 and 2
*
****************************************************************************/
long check_for_orthogonal_basis(point_type vects[3],long dim,real tol)
{
  long result=0;

  if( dim == 0 ){
    NONFATAL_BUG("check_for_orthogonal_basis called for a molecule.");
    return(0);
  }
  if( dim >= 2){
    if(fabs(dot_prod(&(vects[0]),&(vects[1]))) > tol){
      result = result | 1;
    }
    if( dim == 3 ){
      if(fabs(dot_prod(&(vects[0]),&(vects[2]))) > tol){
	result = result | 2;
      }
      if(fabs(dot_prod(&(vects[1]),&(vects[2]))) > tol){
	result = result | 4;
      }
    }
  }
  return(result);
}

/****************************************************************************
*
*                   Function make_atoms_equiv
*
* Arguments:  cell: pointer to cell_type
*         loc1,loc2: pointers to point_type
*        which_cell: pointer to int
*          symm_tol: real
*        cell_dim: pointer to point_type
*
* Returns: int
*
* Action: determines if the locations pointed to by 'loc1 and 'loc2
*    can be made equivalent by translation vectors  of 'cell.
*
****************************************************************************/


/****************************************************************************
*
*                   Function atoms_are_equiv
*
* Arguments:  cell: pointer to cell_type
*         loc1,loc2: pointers to point_type
*        which_cell: pointer to int
*          symm_tol: real
*        cell_dim: pointer to point_type
*
* Returns: int
*
* Action: determines if the locations pointed to by 'loc1 and 'loc2
*    can be made equivalent by translation vectors  of 'cell.
*
****************************************************************************/
long atoms_are_equiv(cell_type *cell,point_type *loc1,point_type *loc2,
		    long *which_cell,real symm_tol,point_type *cell_dim)
{
  long a_loc,b_loc,c_loc;
  long i,j,k,itab,jtab;
  point_type temploc;
  long result;
  

  /* start simply, just compare the positions */
  if( POINTS_ARE_THE_SAME(loc1,loc2,symm_tol) ){
    which_cell[0] = 0;which_cell[1] = 0;which_cell[2] = 0;
    return(1);
  }

  /*******

    loop over translation vectors... this could be more efficient,
      but we are not going to be spending much time here anyway.

  *******/
  if( cell->dim >= 1 ){
    for(i=-1;i<2;i++){
      temploc.x = loc2->x + i*cell_dim[0].x;
      temploc.y = loc2->y + i*cell_dim[0].y;
      temploc.z = loc2->z + i*cell_dim[0].z;
      if( POINTS_ARE_THE_SAME(loc1,&temploc,symm_tol) ){
	which_cell[0] = i;which_cell[1] = 0;which_cell[2] = 0;
	return(1);
      }
      if( cell->dim > 1 ){
	for(j=-1;j<2;j++){
	  temploc.x = loc2->x + i*cell_dim[0].x+ j*cell_dim[1].x;
	  temploc.y = loc2->y + i*cell_dim[0].y+ j*cell_dim[1].y;
	  temploc.z = loc2->z + i*cell_dim[0].z+ j*cell_dim[1].z;
	  if( POINTS_ARE_THE_SAME(loc1,&temploc,symm_tol) ){
	    which_cell[0] = i;which_cell[1] = j;which_cell[2] = 0;
	    return(1);
	  }
	  if( cell->dim == 3 ){
	    for(k=-1;k<2;k++){
	      temploc.x = loc2->x + i*cell_dim[0].x +
		j*cell_dim[1].x+ k*cell_dim[2].x; 
	      temploc.y = loc2->y + i*cell_dim[0].y +
		j*cell_dim[1].y+ k*cell_dim[2].y;
	      temploc.z = loc2->z + i*cell_dim[0].z +
		j*cell_dim[1].z+ k*cell_dim[2].z;
	      if( POINTS_ARE_THE_SAME(loc1,&temploc,symm_tol) ){
		which_cell[0] = i;which_cell[1] = j;which_cell[2] = k;
		return(1);
	      }
	    }
	  }
	}
      }
    }
  }

  return(0);
}


/****************************************************************************
*
*                   Procedure compare_crystal_lattice
*
* Arguments:  cell: pointer to cell_type
*       vects1,vects2: pointers to point_type
*           present: pointer to char
*          symm_tol: real
*         mapping: pointer to real
*
* Returns: none
*
* Action:  checks to see if the two sets of lattice vectors
*    are the same
*  
*   if they are then 'present is  set to be nonzero
*
*   the mapping between vectors is passed back in 'mapping, which should
*   be at least nine elements long.
*
*****************************************************************************/
void compare_crystal_lattice(cell_type *cell,point_type *vects1,point_type *vects2,
			   char *present,real symm_tol,real *mapping)
{
  long i,j,k,l,m;
  char found;
  long which_cell[3];
  point_type temp_vect;
  real dota,dotb,dotc;

  /* initialize present to 0 so that we can bomb out at any time */
  *present = 0;

  memset(mapping,0,9*sizeof(real));

  /*********

    loop over each of the translation vectors

  *********/
  for(i=0;i<cell->dim;i++){
    found = 0;

    /********

      start off the easy way, by just checking against
      plus and minus all the other lattice vectors....
      if we are lucky (i.e. the lattice vectors are orthogonal to
      each other), then this is all we need to do.  If not, we have
      to do a little bit more work...

      once again, this could be more efficient, but this is not
      going to be where the program spends its time, so it's not
      worth making the code overly complicated... in this case
      the brute force approach is (I think) easier to understand.

    *********/
    for(j=0;j<cell->dim && !found;j++){
      temp_vect.x = vects1[j].x;
      temp_vect.y = vects1[j].y;
      temp_vect.z = vects1[j].z;
      if(POINTS_ARE_THE_SAME(&(temp_vect),&(vects2[i]),symm_tol)){
	found = 1;
	mapping[i*3+j] = 1;
      }
      if( !found ){
	temp_vect.x = -vects1[j].x;
	temp_vect.y = -vects1[j].y;
	temp_vect.z = -vects1[j].z;
	if(POINTS_ARE_THE_SAME(&(temp_vect),&(vects2[i]),symm_tol)){
	  found = 1;
	  mapping[i*3+j] = -1;
	}
      }
    }
    if( !found ){
      /********

	okay, we did not find it... sooooo that means we have
	to explore sums and differences of lattice vectors... ick.

	We'll use dot products here to eliminate vectors from
	consideration... a zero dot product indicates no possibility
	of contribution.

      ********/
      dota = dot_prod(&vects1[0],&vects2[i]);
      if( cell->dim > 1 ){
	dotb = dot_prod(&vects1[1],&vects2[i]);
	if(cell->dim > 2 ){
	  dotc = dot_prod(&vects1[2],&vects2[i]);
	}else dotc = 0;
      }else dotb = 0;

      if( fabs(dota) > symm_tol){
	if( fabs(dotb) > symm_tol ){
	  /* this does all possible combinations of +/-a +/- b */
	  for(k=-1;k<2 && !found ;k++){
	    for(l=-1;l<2 && !found ;l++){
	      temp_vect.x = k*vects1[0].x + l*vects1[1].x;
	      temp_vect.y = k*vects1[0].y + l*vects1[1].y;
	      temp_vect.z = k*vects1[0].z + l*vects1[1].z;
	      if(POINTS_ARE_THE_SAME(&(vects2[i]),&(temp_vect),symm_tol)){
		found = 1;
		mapping[i*3] = k;
		mapping[i*3+1] = l;
	      }
	    }
	  }
	}
	if( !found && fabs(dotc) > symm_tol ){
	  /* this does all possible combinations of +/-a +/- c */
	  for(k=-1;k<2 && !found ;k++){
	    for(l=-1;l<2 && !found ;l++){
	      temp_vect.x = k*vects1[0].x + l*vects1[2].x;
	      temp_vect.y = k*vects1[0].y + l*vects1[2].y;
	      temp_vect.z = k*vects1[0].z + l*vects1[2].z;
	      if(POINTS_ARE_THE_SAME(&(vects2[i]),&(temp_vect),symm_tol)){
		found = 1;
		mapping[i*3] = k;
		mapping[i*3+2] = l;
	      }
	    }
	  }
	}
      }
      if( !found && fabs(dotb) > symm_tol && fabs(dotc) > symm_tol ){
	/* this does all possible combinations of +/-b +/- c */
	for(k=-1;k<2 && !found ;k++){
	  for(l=-1;l<2 && !found ;l++){
	    temp_vect.x = k*vects1[1].x + l*vects1[2].x;
	    temp_vect.y = k*vects1[1].y + l*vects1[2].y;
	    temp_vect.z = k*vects1[1].z + l*vects1[2].z;
	    if(POINTS_ARE_THE_SAME(&(vects2[i]),&(temp_vect),symm_tol)){
	      found = 1;
	      mapping[i*3+1] = k;
	      mapping[i*3+2] = l;
	    }
	  }
	}
      }
    }
    /* if we didn't find the vector, return failure now */
    if( !found ) return;
  }
  *present = 1;
}


/****************************************************************************
*
*                   Procedure compare_crystal_basis
*
* Arguments:  cell: pointer to cell_type
*       locs1,locs2: pointers to point_type
*         num_atoms: int
*       equiv_atoms: pointer to int
*           present: pointer to char
*          symm_tol: real
*
* Returns: none
*
* Action:  checks to see if the two sets of atomic positions
*   in 'locs1 and 'locs2 are the same using the translation vectors
*   defined in 'cell.
*  
*   As the name implies, this just checks the basis, it is also
*    important to make sure the lattice is invariant, you have
*    to use compare_crystal_lattice to do this.
*
*   if they are then 'present is  set to be nonzero
*
*  This also constructs the list of equivalent atoms: 'equiv_atoms
*
*****************************************************************************/
void compare_crystal_basis(cell_type *cell,point_type *locs1,point_type *locs2,
			   point_type *cell_dim,
			   long num_atoms,long *equiv_atoms,char *present,real symm_tol)
{
  long i,j;
  char found;
  long which_cell[3];

  /* initialize present to 0 so that we can bomb out at any time */
  *present = 0;


  /*********

    loop over each of the atoms in locs1, checking to see if it is present in locs2

  *********/
  for(i=0;i<num_atoms;i++){
    found = 0;
    /* don't do dummy atoms */
    if( cell->atoms[i].at_number >= 0 ){
      for(j=0;j<num_atoms && !found;j++){
	if( cell->atoms[j].at_number >= 0 ){
	  /* check to see if they are the same atom type */
	  if( cell->atoms[i].at_number == cell->atoms[j].at_number ){
	    /* now check to see if they are in the same location */
	    if( atoms_are_equiv(cell,&locs1[i],&locs2[j],which_cell,symm_tol,
				cell_dim)){
	      /* they are */
	      found = 1;
	      equiv_atoms[j] = i;
	    }
	  }
	}
      }
    } else{
      /* it's a dummy atom, pretend we found it and put a -1 in the equiv_atoms array */
      found = 1;
      equiv_atoms[i] = -1;
    }
    
    /* if we didn't find this atom, we might as well go ahead and return */
    if( !found ) return;
  }
  /* we found every atom, so set present to 1 and return */
  *present = 1;
}

