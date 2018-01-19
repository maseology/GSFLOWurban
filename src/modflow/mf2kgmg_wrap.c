

#include "mf2kgmg.h"

/*
 * Allocates GMG data
 */
void mf2kgmg_allocate_ (int* GMGID, int* NCOL, int* NROW, int* NLAY,
                      int* IPREC, int* ISM, int* ISC,
					  double* RELAX, int* ISIZ, int* IERR) {
   MF2KGMG_ALLOCATE(GMGID, NCOL, NROW, NLAY,
                      IPREC, ISM, ISC,
                      RELAX, ISIZ, IERR);
}

void mf2kgmg_free_ (int* gmgid) {
	MF2KGMG_FREE(gmgid);
}

///*
// *  MF2KGMG_FREE: Deallocates  vectors and operators for cell-centered
// *  finite-diference problem.
// */
//void MF2KGMG_FREE(int* GMGID) 
//{
//  MF2KGMG_operator* GMG_ptr=(MF2KGMG_operator*)*GMGID;
//
//  GEN_free(&GMG_ptr->CCFDMG);
//  GEN_free(&GMG_ptr->PCG);
//  r_free(&GMG_ptr->r);
//  r_free(&GMG_ptr->z);
//
//  /* Dealocate CCFD outside normal contex */
//  free(GMG_ptr->CCFD_ptr->DD);
//  free(GMG_ptr->CCFD_ptr);
//
//  return;
//}

/*
 * Local functions for assembling CCFD matrix.
 */
/* Single precision CCFD assembly method. */
void MF2KGMG_SCCFD_ASSEMBLE(MF2KGMG_operator* GMG_ptr, 
                            double* BIGR0, float* RHS, float* HCOF,
                            float* HNOFLO, double* HNEW);

/* Double precision CCFD assembly method. */
void MF2KGMG_DCCFD_ASSEMBLE(MF2KGMG_operator* GMG_ptr, 
                            double* BIGR0, double* RHS, double* HCOF,
                            double* HNOFLO, double* HNEW);

/* Assemble GMG data:
 * Arguments that are pointer to void are either singel precision
 * or double precision and are resolved at run time.
 */

void mf2kgmg_assemble_(int* GMGID, double* BIGR0, void* CR, void* CC, void* CV,
                      void* HCOF, double* HNEW, void* RHS,
                      void* HNOFLO, int* IBOUND, int* IERR)
{
   MF2KGMG_ASSEMBLE(GMGID, BIGR0, CR, CC, CV,
                      HCOF, HNEW, RHS,
                      HNOFLO, IBOUND, IERR);
}

/*
 *  MF2KGMG_EVAL: Computes head change returning l2-norm of residual (BIGR)
 *  and number of iterations (ITER).
 */
void mf2kgmg_eval_(int* GMGID, int* ITER, double* BIGR, double* DRCLOSE,
                  int* IITER, int* IOUTGMG, int* IOUT)

{
   MF2KGMG_EVAL(GMGID, ITER, BIGR, DRCLOSE,
                  IITER, IOUTGMG, IOUT);
}

/*
 *  MF2KGMG_UPDATE: Adds damped head change to current approximation.
 */
void mf2kgmg_update_(int* GMGID, double* HNEW, double* DDAMP)
{
   MF2KGMG_UPDATE(GMGID, HNEW, DDAMP);
}


/* Calculate l2-norm of residual and return location
 * of max residual.
 */
void mf2kgmg_bigr_(int* GMGID, double* BIGR, int* IBIGR, int* JBIGR, int* KBIGR)
{
   MF2KGMG_BIGR(GMGID, BIGR, IBIGR, JBIGR, KBIGR);
}

/*
 * Calculate max of head change and return location of max head change.
 * Absalute value of BIGH is max-norm of head change.
 */
void mf2kgmg_bigh_(int* GMGID, double* BIGH, int* IBIGH, int* JBIGH, int* KBIGH)
{
 MF2KGMG_BIGH(GMGID, BIGH, IBIGH, JBIGH, KBIGH);
}

void mf2kgmg_sccfd_assemble_(MF2KGMG_operator* GMG_ptr,
                            double* BIGR0, float* RHS, float* HCOF,
                            float* HNOFLO, double* HNEW)
{
   MF2KGMG_SCCFD_ASSEMBLE(GMG_ptr,
                            BIGR0, RHS, HCOF,
                            HNOFLO, HNEW);
}

void mfwkgmg_dccfd_assemble_(MF2KGMG_operator* GMG_ptr,
                            double* BIGR0, double* RHS, double* HCOF,
                            double* HNOFLO, double* HNEW)
{
   MF2KGMG_DCCFD_ASSEMBLE(GMG_ptr,
                            BIGR0, RHS, HCOF,
                            HNOFLO, HNEW);
}

void resprint(int* IOUT, int* I, double* RES, double* CFAC) {
	resprint_(IOUT, I, RES, CFAC);
}
