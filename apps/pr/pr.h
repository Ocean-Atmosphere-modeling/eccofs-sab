/* 
** svn $Id 
******************************************************************************* 
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               ** 
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    ** 
******************************************************************************* 
** 
** Experimental System for Predicting Shelf and Slope Optics (ONR MAB MURI)  */

#define NLM_DRIVER

/* Basic physics options */
#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV                       /* momentum mixing on s-surfaces */
#define TS_DIF2
#define MIX_GEO_TS              /* tracer mixing on constant z surfaces */
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CRAIG_BANNER
#define CHARNOK
#define WIND_MINUS_CURRENT

#define ATM_PRESS
#define PRESS_COMPENSATE

#define  WTYPE_GRID
#define  SINGLE_PRECISION
#define OMEGA_IMPLICIT
#define COOL_SKIN
#define  LIMIT_STFLX_COOLING


/* Basic numerics options */
#define UV_U3HADVECTION
#define UV_C4VADVECTION
#define DJ_GRADPS
#define CURVGRID
#define MASKING

/* Outputs */
#undef OUT_DOUBLE
#define AVERAGES
#define AVERAGES_FLUXES
#define AVERAGES_AKV           
#define AVERAGES_AKT           
#define AVERAGES_AKS 

/* Surface and bottom boundary conditions */
#define UV_QDRAG
#define BULK_FLUXES    
#define GLS_MIXING
#define EMINUSP     /* evap from latent heat and combine with NCEP rain */
#define LONGWAVE_OUT /* define to read downward longwave, compute outgoing */
#define KANTHA_CLAYSON
#define N2S2_HORAVG
#undef SSH_TIDES       /* Activated tides for initial, simple case */
#undef ADD_FSOBC         /* Tide data is added to OBC from HYCOM */
#undef UV_TIDES        /* Reactivated when tidal data acquired     */
#undef ADD_M2OBC         /* Tide data is added to OBC from HYCOM */
#undef TIDE_GENERATING_FORCES
#define SOLAR_SOURCE   /* solar shortwave distributed over water column */
#define ANA_BSFLUX
#define ANA_BTFLUX

#define NO_LBC_ATT
