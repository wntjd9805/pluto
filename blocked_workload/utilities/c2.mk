UTIL_MK_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
UTIL_DIR ?= $(UTIL_MK_DIR)
-include $(UTIL_DIR)/common.mk
