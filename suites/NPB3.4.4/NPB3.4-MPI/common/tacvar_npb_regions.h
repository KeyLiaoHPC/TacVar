/**
 * @file tacvar_npb_regions.h
 * @brief Register logical NPB-MPI timer regions into region_info.csv.
 */
#ifndef TACVAR_NPB_REGIONS_H
#define TACVAR_NPB_REGIONS_H

#ifdef __cplusplus
extern "C" {
#endif

/** Register all logical regions for the named benchmark (e.g. "mg"). */
void tacvar_npb_register_regions(const char *benchmark);

#ifdef __cplusplus
}
#endif
#endif
