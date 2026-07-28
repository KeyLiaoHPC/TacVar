/* Minimal stub when libtirpc-devel is unavailable. Enough for non-RPC lmbench bins. */
#ifndef TACVAR_STUB_RPC_TYPES_H
#define TACVAR_STUB_RPC_TYPES_H
#include <sys/types.h>
typedef int bool_t;
typedef unsigned long u_long;
typedef unsigned int u_int;
typedef char *caddr_t;
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif
#endif
