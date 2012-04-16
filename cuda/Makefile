################################################################################
#
# Build script for project
#
################################################################################

verbose = 1
#maxregisters = 20
SMVERSION_template = 1
#usesm20 = 1

# Add source files here
EXECUTABLE	   := entbody 
CUFILES_sm_13  := main.cu
CCFILES		   := 
NVCCFLAGS += -g -G
CFFLAGS += -g -G
################################################################################
# Rules and targets
ROOTDIR := /a/CUDA_SDK_31/C/common
BINDIR := ./build
include common.mk

