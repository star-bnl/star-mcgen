/*:>--------------------------------------------------------------------
**: FILE:       astore_xdf.c.template
**: HISTORY:
**:             00jan93-v000a-cet- Created by stic Version
**:  Id: idl.y,v 1.8 1996/10/15 18:33:35 ward Exp  
**:<------------------------------------------------------------------*/
#include "astore_xdf.h"

long type_of_call astore_xdf_(
  TABLE_HEAD_ST         *gtable_h,       EG_GENER_ST           *gtable ,
  TABLE_HEAD_ST         *etable_h,       EG_EVENT_ST           *etable ,
  TABLE_HEAD_ST         *ttable_h,       EG_TRACK_ST           *ttable ,
  TABLE_HEAD_ST         *vtable_h,      EG_VERTEX_ST           *vtable )
{
/*:>--------------------------------------------------------------------
**: ROUTINE:    astore_xdf_
**: DESCRIPTION: Physics Analysis Module ANSI C template.
**:             This is an ANSI C Physics Analysis Module template
**:             automatically generated by stic from astore_xdf.idl.
**:             Please edit comments and code.
**: AUTHOR:     cet - C.E.Tull, cetull@lbl.gov
**: ARGUMENTS:
**:       IN:
**:             gtable    - Generator table
**:             etable    - Event table
**:             ttable    - Track table
**:             vtable    - Vertex table
**:    INOUT:
**:      OUT:
**: RETURNS:    STAF Condition Value
**:		** KILLS PROGRAM **
**:>------------------------------------------------------------------*/

   static FILE *f=NULL;
   static XDR x;
   static DS_DATASET_T *d=NULL;		/* "root" dataset */
   static DS_DATASET_T *g=NULL;		/* generator table */
   static DS_DATASET_T *e=NULL;		/* event table */
   static DS_DATASET_T *t=NULL;		/* track table */
   static DS_DATASET_T *v=NULL;		/* vertex table */

   if( NULL == f ){
      if( (NULL == (f = fopen("evgen.xdf", "wb")))
      ||  !dsNewDataset(&d,"evgen")
      ||  !dsAddTable(d,"gener", EG_GENER_SPEC, gtable_h->maxlen
		, &gtable)
      ||  !dsFindEntry(&g,d,"gener")
      ||  !dsAddTable(d,"event", EG_EVENT_SPEC, etable_h->maxlen
		, &etable)
      ||  !dsFindEntry(&e,d,"event")
      ||  !dsAddTable(d,"track", EG_TRACK_SPEC, ttable_h->maxlen
		, &ttable)
      ||  !dsFindEntry(&t,d,"track")
      ||  !dsAddTable(d,"vertex", EG_VERTEX_SPEC,vtable_h->maxlen
		, &vtable)
      ||  !dsFindEntry(&v,d,"vertex")
      ){
	 exit(0);			/* KILL PROGRAM */
	 return STAFCV_BAD;
      }
   }

   xdrstdio_create(&x,f,XDR_ENCODE);

   if( !dsSetTableRowCount(g,gtable_h->nok)
   ||  !dsSetTableRowCount(e,etable_h->nok)
   ||  !dsSetTableRowCount(t,ttable_h->nok)
   ||  !dsSetTableRowCount(v,vtable_h->nok)
   ){
      exit(0);			/* KILL PROGRAM */
      return STAFCV_BAD;
   }

   xdr_dataset(&x,&d);

/* Successful completion of analysis module... */
   return STAFCV_OK;
}
