/**
 * @file x86_events.c
 * @brief Portable x86 PERFEVTSEL aliases for kmod/rdpmc backend.
 *
 * Codes include USR|OS|EN bits (0x430000) where applicable.
 * Raw hex strings (0x...) are accepted as full PERFEVTSEL values.
 */
#include "events.h"
#include <stdlib.h>
#include <string.h>
#include <strings.h>

typedef struct {
    const char *name;
    uint64_t code;
} tacvar_x86_event_t;

/* USR=1 OS=1 EN=1 → bits 16,17,22 = 0x430000 */
#define X86_EV(event, umask) (0x430000ULL | ((uint64_t)(umask) << 8) | (uint64_t)(event))

static const tacvar_x86_event_t g_x86_events[] = {
    { "CPU_CLK_UNHALTED.THREAD",     X86_EV(0x3C, 0x00) },
    { "CPU_CLK_THREAD_UNHALTED:THREAD_P", X86_EV(0x3C, 0x00) },
    { "INST_RETIRED.ANY_P",          X86_EV(0xC0, 0x00) },
    { "INSTRUCTIONS",                X86_EV(0xC0, 0x00) },
    { "BR_INST_RETIRED.ALL_BRANCHES", X86_EV(0xC4, 0x00) },
    { "BR_MISP_RETIRED.ALL_BRANCHES", X86_EV(0xC5, 0x00) },
    { "LONGEST_LAT_CACHE.REFERENCE", X86_EV(0x2E, 0x4F) },
    { "LONGEST_LAT_CACHE.MISS",      X86_EV(0x2E, 0x41) },
    { "MEM_INST_RETIRED.ALL_LOADS",  X86_EV(0xD0, 0x81) },
    { "MEM_INST_RETIRED.ALL_STORES", X86_EV(0xD0, 0x82) },
    { NULL, 0 }
};

uint64_t tacvar_x86_parse_event(const char *name)
{
    int i;
    char *end = NULL;
    unsigned long long v;

    if (!name || !*name)
        return 0;
    if (name[0] == '0' && (name[1] == 'x' || name[1] == 'X')) {
        v = strtoull(name, &end, 16);
        if (end && *end == '\0' && v != 0)
            return (uint64_t)v;
        return 0;
    }
    for (i = 0; g_x86_events[i].name; i++) {
        if (strcasecmp(name, g_x86_events[i].name) == 0)
            return g_x86_events[i].code;
    }
    return 0;
}
