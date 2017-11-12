#ifndef _ZERO_FLIP_H
#define _ZERO_FLIP_H
//
//
//
//
//
//
//
//--------------------------DEBUG MACRO---------------------
// #ifndef _ZERO_NDEBUG_FLIP
// #define _ZERO_NDEBUG_FLIP
// #endif

// #ifndef _CARTESIAN_LS_DEBUG_FLIP
// #define _CARTESIAN_LS_DEBUG_FLIP
// #endif

/*#ifndef _SINGLE_STEP_DEBUG_FLIP
#define _SINGLE_STEP_DEBUG_FLIP
#endif*/
//----------------------------------------------------------
//----------------------------------------------------------

#ifndef _MESHTYPE_ARK
#define _MESHTYPE_ARK "Quad"
#endif

#ifndef _FLUX_SCHEME_ARK
#define _FLUX_SCHEME_ARK "UW"							//UW = up wind,CD = central difference;
#endif

#ifndef _BC_ARK
#define _BC_ARK "_4_5" //BB = bounce back,NEE = non-equilibrium extrapolation,DS = diffusive scattering
#endif

#ifndef _QMODEL_ARK
#define _QMODEL_ARK "D2V16"
#endif

#ifndef _MESHFILE_NAME_ARK
#define _MESHFILE_NAME_ARK "_Car_Blend_ShockTube"
#endif
//----------------Boundary Condition Macro------------------
/*#ifndef _CARTESIAN_MESH_FLIP
#define _CARTESIAN_MESH_FLIP
#endif*/

/*#ifndef _BB_BOUNDARY_SCHEME_FLIP
#define _BB_BOUNDARY_SCHEME_FLIP
#endif*/

/*#ifndef _NEE_BOUNDARY_SCHEME_FLIP
#define _NEE_BOUNDARY_SCHEME_FLIP
#endif*/

#ifndef _P_INLET_4_BCS_FLIP
#define _P_INLET_4_BCS_FLIP
#endif

#ifndef _P_OUTLET_5_BCS_FLIP
#define _P_OUTLET_5_BCS_FLIP
#endif

#ifndef _PERIODIC_12_8_BCs_FLIP
#define _PERIODIC_12_8_BCs_FLIP
#endif

/*#ifndef _Wall_3_BCs_FLIP
#define _Wall_3_BCs_FLIP
#endif*/

//-------------------------------Isothermal-------------------------------
#ifndef _ARK_ISOTHERMAL_FLIP
#define _ARK_ISOTHERMAL_FLIP
#endif
//
#ifndef _ARK_LIMITER_FLIP
#define _ARK_LIMITER_FLIP
#endif
//----------------------------------------------------------------------------

/*#ifndef _OUTPUT_L2NORM_ERROR_FLIP
#define _OUTPUT_L2NORM_ERROR_FLIP
#endif*/

#ifndef _PRINT_ERROR_MSG_FLIP
#define _PRINT_ERROR_MSG_FLIP  cout<<"File : "<<__FILE__<<"  Line : "\
<<__LINE__<<"  fun : "<<__func__<<'\n';
#endif
//
//
//
//
//
//
//
//
//
#endif