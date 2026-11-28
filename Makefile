# liblzr (static + shared library) and the lzr command-line tool.
#
#   make                  build/c/liblzr.a, build/c/liblzr.so, build/c/lzr
#   make check            C API tests + CLI / tar round trips (tests/cli/run.sh)
#   make install          PREFIX=/usr/local, DESTDIR supported
#   make clean
#
# The Python extension (setup.py) compiles the same core sources.

PREFIX  ?= /usr/local
BINDIR  ?= $(PREFIX)/bin
LIBDIR  ?= $(PREFIX)/lib
INCDIR  ?= $(PREFIX)/include
MANDIR  ?= $(PREFIX)/share/man
PKGDIR  ?= $(LIBDIR)/pkgconfig

CC      ?= cc
CFLAGS  ?= -O3
WARN    = -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers
ALL_CFLAGS = -std=gnu11 $(WARN) $(CFLAGS) -Iinclude -pthread
LDLIBS  = -lm -pthread

VERSION := $(shell sed -n 's/^\#define LZR_LIB_VERSION_\(MAJOR\|MINOR\|PATCH\) *\([0-9]*\).*/\2/p' include/lzr.h | paste -sd.)
MAJOR   := $(firstword $(subst ., ,$(VERSION)))

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SHLIB      = liblzr.$(VERSION).dylib
  SHLIB_LINK = liblzr.dylib
  SHLIB_SO   = liblzr.$(MAJOR).dylib
  SHFLAGS    = -dynamiclib -install_name $(LIBDIR)/$(SHLIB_SO) -compatibility_version $(MAJOR) -current_version $(VERSION)
else
  SHLIB      = liblzr.so.$(VERSION)
  SHLIB_LINK = liblzr.so
  SHLIB_SO   = liblzr.so.$(MAJOR)
  SHFLAGS    = -shared -Wl,-soname,$(SHLIB_SO)
endif

B       = build/c
CORE    = api enc dec opt crc32 huf fse
HDRS    = $(wildcard src/lzrans/*.h) include/lzr.h
STATIC_OBJS = $(CORE:%=$(B)/static/%.o)
SHARED_OBJS = $(CORE:%=$(B)/shared/%.o)

all: lib $(B)/lzr
lib: $(B)/liblzr.a $(B)/$(SHLIB)

$(B)/static/%.o: src/lzrans/%.c $(HDRS)
	@mkdir -p $(@D)
	$(CC) $(ALL_CFLAGS) -c -o $@ $<

$(B)/shared/%.o: src/lzrans/%.c $(HDRS)
	@mkdir -p $(@D)
	$(CC) $(ALL_CFLAGS) -fPIC -fvisibility=hidden -c -o $@ $<

$(B)/liblzr.a: $(STATIC_OBJS)
	rm -f $@
	$(AR) rcs $@ $^

$(B)/$(SHLIB): $(SHARED_OBJS)
	$(CC) $(ALL_CFLAGS) $(LDFLAGS) $(SHFLAGS) -o $@ $^ $(LDLIBS)
	ln -sf $(SHLIB) $(B)/$(SHLIB_SO)
	ln -sf $(SHLIB) $(B)/$(SHLIB_LINK)

# The CLI links the static library: one self-contained binary.
$(B)/lzr: cli/lzr.c include/lzr.h $(B)/liblzr.a
	$(CC) $(ALL_CFLAGS) $(LDFLAGS) -o $@ cli/lzr.c $(B)/liblzr.a $(LDLIBS)
	ln -sf lzr $(B)/unlzr
	ln -sf lzr $(B)/lzrcat

$(B)/test_api: tests/c/test_api.c include/lzr.h $(B)/$(SHLIB)
	$(CC) $(ALL_CFLAGS) $(LDFLAGS) -o $@ tests/c/test_api.c -L$(B) -llzr -Wl,-rpath,$(abspath $(B)) $(LDLIBS)

$(B)/lzr.pc: Makefile include/lzr.h
	@mkdir -p $(@D)
	printf 'prefix=%s\nlibdir=%s\nincludedir=%s\n\nName: liblzr\nDescription: lzrans (LZR 1.0) compression library\nVersion: %s\nLibs: -L$${libdir} -llzr\nLibs.private: -lm -pthread\nCflags: -I$${includedir}\n' \
	  '$(PREFIX)' '$(LIBDIR)' '$(INCDIR)' '$(VERSION)' > $@

check: all $(B)/test_api
	$(B)/test_api
	sh tests/cli/run.sh $(B)/lzr

install: all $(B)/lzr.pc
	install -d $(DESTDIR)$(BINDIR) $(DESTDIR)$(LIBDIR) $(DESTDIR)$(INCDIR) $(DESTDIR)$(PKGDIR) $(DESTDIR)$(MANDIR)/man1
	install -m 755 $(B)/lzr $(DESTDIR)$(BINDIR)/lzr
	ln -sf lzr $(DESTDIR)$(BINDIR)/unlzr
	ln -sf lzr $(DESTDIR)$(BINDIR)/lzrcat
	install -m 644 $(B)/liblzr.a $(DESTDIR)$(LIBDIR)/liblzr.a
	install -m 755 $(B)/$(SHLIB) $(DESTDIR)$(LIBDIR)/$(SHLIB)
	ln -sf $(SHLIB) $(DESTDIR)$(LIBDIR)/$(SHLIB_SO)
	ln -sf $(SHLIB) $(DESTDIR)$(LIBDIR)/$(SHLIB_LINK)
	install -m 644 include/lzr.h $(DESTDIR)$(INCDIR)/lzr.h
	install -m 644 $(B)/lzr.pc $(DESTDIR)$(PKGDIR)/lzr.pc
	install -m 644 docs/lzr.1 $(DESTDIR)$(MANDIR)/man1/lzr.1

uninstall:
	rm -f $(DESTDIR)$(BINDIR)/lzr $(DESTDIR)$(BINDIR)/unlzr $(DESTDIR)$(BINDIR)/lzrcat \
	  $(DESTDIR)$(LIBDIR)/liblzr.a $(DESTDIR)$(LIBDIR)/$(SHLIB) $(DESTDIR)$(LIBDIR)/$(SHLIB_SO) \
	  $(DESTDIR)$(LIBDIR)/$(SHLIB_LINK) $(DESTDIR)$(INCDIR)/lzr.h $(DESTDIR)$(PKGDIR)/lzr.pc \
	  $(DESTDIR)$(MANDIR)/man1/lzr.1

clean:
	rm -rf $(B)

.PHONY: all lib check install uninstall clean
