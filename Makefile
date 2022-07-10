PROJNAME      = egm96
SRCS          = src
PLIBNAME      = lib$(PROJNAME).a
PROGNAME1     = $(PROJNAME)
PROGNAME2     = $(PROJNAME)_gen
PROGS         = $(PROGNAME1) $(PROGNAME2)
CC            = gcc
CFLAGS        = -Wall
LDFLAGS       = -lm -s
AR            = ar
PREFIX        = /usr/local
INSTALL       = install
LN            = ln -fs
GROFF2PDF     = groff -m man -T pdf
RM            = rm -f

.PHONY: all clean install

all: $(PROGS)

clean:
	$(RM) $(SRCS)/*.o $(PLIBNAME) $(PROGS) man_*.pdf

$(SRCS)/$(PROJNAME).o: $(SRCS)/$(PROJNAME).c
	$(CC) $(CFLAGS) -c $^ -o $@

$(PLIBNAME): $(SRCS)/$(PROJNAME).o
	$(AR) rcs $@ $^

$(PROGNAME1): $(SRCS)/$(PROGNAME1)_cli.c $(PLIBNAME)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(PROGNAME2): $(SRCS)/$(PROGNAME2).c $(PLIBNAME)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

manual: man_$(PROGNAME1).pdf man_$(PROGNAME2).pdf

man_$(PROGNAME1).pdf: man/man1/$(PROGNAME1).1
	$(GROFF2PDF) $^ > $@

man_$(PROGNAME2).pdf: man/man1/$(PROGNAME2).1
	$(GROFF2PDF) $^ > $@

install: $(PROGS)
	$(INSTALL) -d $(PREFIX)/bin
	$(INSTALL) -m 0755 $(PROGS) $(PREFIX)/bin/
	$(INSTALL) -d $(PREFIX)/share/man/man1
	$(INSTALL) -m 0644 man/man1/*.1 $(PREFIX)/share/man/man1
	$(INSTALL) -d $(PREFIX)/share/doc/$(PROJNAME)
	$(INSTALL) -m 0644 LICENSE README.md doc/* $(PREFIX)/share/doc/$(PROJNAME)
