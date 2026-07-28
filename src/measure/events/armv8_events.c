/**
 * @file armv8_events.c
 * @brief ARMv8 architectural PMU event aliases.
 */
#include "events.h"
#include <stdlib.h>
#include <string.h>
#include <strings.h>

typedef struct {
    const char *name;
    uint32_t code;
} tacvar_arm_event_t;

static const tacvar_arm_event_t g_arm_events[] = {
    { "SW_INCR",            0x00 },
    { "L1I_CACHE_REFILL",   0x01 },
    { "L1I_TLB_REFILL",     0x02 },
    { "L1D_CACHE_REFILL",   0x03 },
    { "L1D_CACHE",          0x04 },
    { "L1D_TLB_REFILL",     0x05 },
    { "LD_RETIRED",         0x06 },
    { "ST_RETIRED",         0x07 },
    { "INST_RETIRED",       0x08 },
    { "INSTRUCTIONS",       0x08 },
    { "EXC_TAKEN",          0x09 },
    { "EXC_RETURN",         0x0A },
    { "CID_WRITE_RETIRED",  0x0B },
    { "PC_WRITE_RETIRED",   0x0C },
    { "BR_IMMED_RETIRED",   0x0D },
    { "BR_RETURN_RETIRED",  0x0E },
    { "UNALIGNED_LDST_RETIRED", 0x0F },
    { "BR_MIS_PRED",        0x10 },
    { "CPU_CYCLES",         0x11 },
    { "CYCLES",             0x11 },
    { "BR_PRED",            0x12 },
    { "MEM_ACCESS",         0x13 },
    { "L1I_CACHE",          0x14 },
    { "L1D_CACHE_WB",       0x15 },
    { "L2D_CACHE",          0x16 },
    { "L2D_CACHE_REFILL",   0x17 },
    { "L2D_CACHE_WB",       0x18 },
    { "BUS_ACCESS",         0x19 },
    { "MEMORY_ERROR",       0x1A },
    { "INST_SPEC",          0x1B },
    { "TTBR_WRITE_RETIRED", 0x1C },
    { "BUS_CYCLES",         0x1D },
    { "CHAIN",              0x1E },
    { "BR_RETIRED",         0x21 },
    { "BR_MIS_PRED_RETIRED", 0x22 },
    { NULL, 0 }
};

uint32_t tacvar_armv8_parse_event(const char *name)
{
    int i;
    char *end = NULL;
    unsigned long v;

    if (!name || !*name)
        return 0xFFFFFFFFu;
    if (name[0] == '0' && (name[1] == 'x' || name[1] == 'X')) {
        v = strtoul(name, &end, 16);
        if (end && *end == '\0')
            return (uint32_t)v;
        return 0xFFFFFFFFu;
    }
    for (i = 0; g_arm_events[i].name; i++) {
        if (strcasecmp(name, g_arm_events[i].name) == 0)
            return g_arm_events[i].code;
    }
    return 0xFFFFFFFFu;
}
