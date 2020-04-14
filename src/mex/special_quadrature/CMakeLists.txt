include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)


## MEX functions
matlab_add_mex(
	NAME mex_SQ_slp
	SRC mex_SQ_slp.cpp
)

matlab_add_mex(
	NAME mex_SQ_dlp
	SRC mex_SQ_dlp.cpp
)
