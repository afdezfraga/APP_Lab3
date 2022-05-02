TARGETS = stencil

EXTRAITEMS=heat.svg

LDFLAGS=-lm

MPICC=mpicc

all: $(TARGETS)

$(TARGETS): % : %.c
	$(MPICC) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGETS) *.o

cleanall: clean
	rm -f $(EXTRAITEMS)