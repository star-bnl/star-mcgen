/* .h */
/* This file was made by the idl compiler "stic". Do not edit.
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:

        description: crs unpacking program

 */
/* cr_star.h */
#ifndef CR_STAR_H
#define CR_STAR_H
/*----------------------------------------------- INCLUDES   --*/
#include "PAM.h"
/*----------------------------------------------- MACROS     --*/
#define CR_STAR_RANK 0
/*----------------------------------------------- TYPEDEFS   --*/
typedef STAFCV_T (*CR_STAR_FT)
(
);
/*----------------------------------------------- PROTOTYPES --*/
extern CC_P STAFCV_T cr_star_ (
);
#ifdef __cplusplus
extern CC_P STAFCV_T cr_star_load_ami(amiBroker *broker);
#endif /* __cplusplus */
#endif /* CR_STAR_H */
