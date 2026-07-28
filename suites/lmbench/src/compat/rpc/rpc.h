/* Minimal stub when libtirpc-devel is unavailable. Do not build lat_rpc against this. */
#ifndef TACVAR_STUB_RPC_H
#define TACVAR_STUB_RPC_H
#include <rpc/types.h>
typedef void *CLIENT;
typedef void *SVCXPRT;
typedef bool_t (*xdrproc_t)();
typedef struct { long tv_sec; long tv_usec; } rpc_timeout_stub;
#define TIMEOUT ((rpc_timeout_stub){25,0})
#endif
