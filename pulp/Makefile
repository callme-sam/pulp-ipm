SRC_DIR = ../src
INC_DIR = ../inc

APP = pulp-ipm
APP_SRCS += pulp-main.c $(wildcard $(SRC_DIR)/*.c)

ifdef USE_CLUSTER
APP_CFLAGS += -DCLUSTER -DNUM_CLUSTER=$(USE_CLUSTER)
ifdef NUM_CORES
APP_CFLAGS += -DNUM_CORES=$(NUM_CORES)
else
APP_CFLAGS += -DNUM_CORES=1
endif
endif

APP_CFLAGS += -Os -g -I$(INC_DIR) -lm -MMD -MP -DTARGET_EXECUTION
APP_LDFLAGS += -Os -g -lm

include $(RULES_DIR)/pmsis_rules.mk