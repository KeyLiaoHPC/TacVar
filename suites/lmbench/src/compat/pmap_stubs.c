/* Weak stubs for portmap symbols when linking without libtirpc-devel. */
#include <sys/types.h>
#include <rpc/types.h>

int pmap_unset(unsigned long prog, unsigned long vers)
{
    (void)prog; (void)vers; return 0;
}
int pmap_set(unsigned long prog, unsigned long vers, unsigned long prot, unsigned short port)
{
    (void)prog; (void)vers; (void)prot; (void)port; return 0;
}
unsigned short pmap_getport(void *addr, unsigned long prog, unsigned long vers, unsigned int prot)
{
    (void)addr; (void)prog; (void)vers; (void)prot; return 0;
}
void *get_myaddress(void *addr) { (void)addr; return 0; }
