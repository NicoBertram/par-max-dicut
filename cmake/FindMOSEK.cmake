find_path(MOSEK_INCLUDE_DIRS
    NAMES mosek.h fusion.h
    HINTS ${MOSEK_DIR} $ENV{MOSEK_HOME}
    PATH_SUFFIXES h)

find_library(MOSEK_LIBRARY
    NAMES libmosek64.so
    HINTS ${MOSEK_DIR} $ENV{MOSEK_HOME}
    PATH_SUFFIXES bin)

find_library(FUSION_LIBRARY
    NAMES libfusion64.so
    HINTS ${MOSEK_DIR} $ENV{MOSEK_HOME}
    PATH_SUFFIXES bin)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MOSEK DEFAULT_MSG MOSEK_LIBRARY FUSION_LIBRARY)
