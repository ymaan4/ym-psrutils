#all: extract_spulse fetch_chunk normByRMS
all: normByRMS apply_rfifindMask
IDIR=include/
VPATH=src/:src/ext/
ODIR=src/obj/
LDIR=lib
BINDIR=bin/
MYBIN=/home/maan/pulsar_softwares/bin/

CC=gcc
CFLAGS=-I$(IDIR) -Wno-unused-result -O3 -march=native
LIBS=-lm 

_DEPS = header.h apply_rfifindMask.h utils.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = pack_unpack.o  scaledata.o  strings_equal.o  read_block.o  send_stuff.o swap_bytes.o  nsamples.o  read_header.o sizeof_file.o median.o dmedian.o read_block_orig.o med_rms.o utils.o bcast_header.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	@mkdir -p $(BINDIR)
	$(CC) -c -o $@ $< $(CFLAGS)

apply_rfifindMask: $(OBJ) $(ODIR)/apply_rfifindMask.o
	$(CC) -o $(BINDIR)/apply_rfifindMask $^ $(CFLAGS) $(LIBS)

normByRMS: $(OBJ) $(ODIR)/normByRMS.o
	$(CC) -o $(BINDIR)/normByRMS $^ $(CFLAGS) $(LIBS)

.PHONY: install
install:$(BINDIR)/apply_rfifindMask
	install $(BINDIR)/apply_rfifindMask $(MYBIN)/

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o *~ core

