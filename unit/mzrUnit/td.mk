# The name of this module.
MODULE_NAME := mzr

# Other modules (really libs made by this build) that this
# module requires for linking.
REQUIRED_MODULES := domUtils

# Libraries that are required for linking this module,
# but are not made by this build.
EXTRA_LIBS :=

DOT := $(DOT)/$(MODULE_NAME)Unit

include $(DOT)/cc/td.mk
include $(DOT)/dbg-o/td.mk
include $(DOT)/opt-o/td.mk
include $(DOT)/prf-o/td.mk

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~

DOT := $(call dotdot,$(DOT))
