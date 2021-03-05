DIRS := $(wildcard */.)

all: $(DIRS)
clean: $(DIRS)
semiclean: $(DIRS)
atomicBarriers: $(DIRS)
bindCores: $(DIRS)
bindThreads: $(DIRS)
atomicBarriersBindCores: $(DIRS)
atomicBarriersBindThreads: $(DIRS)
$(DIRS):
	+$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: all clean semiclean atomicBarriers bindCores bindThreads atomicBarriersBindCores atomicBarriersBindThreads $(DIRS)
