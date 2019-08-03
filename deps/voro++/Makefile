# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Tell make that these are phony targets
.PHONY: all help clean install uninstall

include config.mk

# Build all of the executable files
all:
	$(MAKE) -C src
	$(MAKE) -C examples

# Build the help files (with Doxygen)
help:
	$(MAKE) -C src help

# Clean up the executable files
clean:
	$(MAKE) -C src clean
	$(MAKE) -C examples clean

# Install the executable, man page, and shared library
install:
	$(MAKE) -C src
	$(INSTALL) -d $(IFLAGS_EXEC) $(INSTALLDIR)/bin
	$(INSTALL) -d $(IFLAGS_EXEC) $(INSTALLDIR)/lib
	$(INSTALL) -d $(IFLAGS_EXEC) $(INSTALLDIR)/man
	$(INSTALL) -d $(IFLAGS_EXEC) $(INSTALLDIR)/man/man1
	$(INSTALL) -d $(IFLAGS_EXEC) $(INSTALLDIR)/include
	$(INSTALL) -d $(IFLAGS_EXEC) $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS_EXEC) src/voro++ $(INSTALLDIR)/bin
	$(INSTALL) $(IFLAGS) man/voro++.1 $(INSTALLDIR)/man/man1
	$(INSTALL) $(IFLAGS) src/libvoro++.a $(INSTALLDIR)/lib
	$(INSTALL) $(IFLAGS) src/voro++.hh $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS) src/c_loops.hh $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS) src/cell.hh $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS) src/common.hh $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS) src/config.hh $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS) src/container.hh $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS) src/container_prd.hh $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS) src/rad_option.hh $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS) src/pre_container.hh $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS) src/unitcell.hh $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS) src/v_base.hh $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS) src/v_compute.hh $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS) src/wall.hh $(INSTALLDIR)/include/voro++
	$(INSTALL) $(IFLAGS) src/worklist.hh $(INSTALLDIR)/include/voro++

# Uninstall the executable, man page, and shared library
uninstall:
	rm -f $(INSTALLDIR)/bin/voro++
	rm -f $(INSTALLDIR)/man/man1/voro++.1
	rm -f $(INSTALLDIR)/lib/libvoro++.a
	rm -f $(INSTALLDIR)/include/voro++/voro++.hh
	rm -f $(INSTALLDIR)/include/voro++/c_loops.hh
	rm -f $(INSTALLDIR)/include/voro++/cell.hh
	rm -f $(INSTALLDIR)/include/voro++/common.hh
	rm -f $(INSTALLDIR)/include/voro++/config.hh
	rm -f $(INSTALLDIR)/include/voro++/container.hh
	rm -f $(INSTALLDIR)/include/voro++/container_prd.hh
	rm -f $(INSTALLDIR)/include/voro++/pre_container.hh
	rm -f $(INSTALLDIR)/include/voro++/rad_option.hh
	rm -f $(INSTALLDIR)/include/voro++/unitcell.hh
	rm -f $(INSTALLDIR)/include/voro++/v_base.hh
	rm -f $(INSTALLDIR)/include/voro++/v_compute.hh
	rm -f $(INSTALLDIR)/include/voro++/wall.hh
	rm -f $(INSTALLDIR)/include/voro++/worklist.hh
	rmdir $(INSTALLDIR)/include/voro++
