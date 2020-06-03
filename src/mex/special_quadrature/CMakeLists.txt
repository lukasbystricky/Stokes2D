include_directories(
  ${CMAKE_SOURCE_DIR}/mex/common
)

## MEX functions
matlab_add_mex(
	NAME mex_SQ_slp
	SRC ${CMAKE_SOURCE_DIR}/mex/common/SQ_common.cpp mex_SQ_slp.cpp
)

matlab_add_mex(
	NAME mex_SQ_dlp
	SRC ${CMAKE_SOURCE_DIR}/mex/common/SQ_common.cpp mex_SQ_dlp.cpp
)
