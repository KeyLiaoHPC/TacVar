/**
 * @file events.h
 * @brief Event name → code helpers for asm counter backends.
 */
#ifndef TACVAR_EVENTS_H
#define TACVAR_EVENTS_H

#include <stdint.h>

/** Parse x86 event alias or hex; returns 0 on failure. */
uint64_t tacvar_x86_parse_event(const char *name);

/** Parse ARMv8 PMU event alias or hex; returns 0xFFFFFFFF on failure. */
uint32_t tacvar_armv8_parse_event(const char *name);

#endif
