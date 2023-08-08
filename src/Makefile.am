ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -Wall -I$(srcdir) -I$(srcdir)/include $(BOOST_CPPFLAGS)
AM_CXXFLAGS = $(ALL_CXXFLAGS) -std=c++0x $(JELLYFISH2_0_CFLAGS)
AM_LDFLAGS = $(JELLYFISH2_0_LIBS)

noinst_HEADERS = $(YAGGO_SOURCES)
bin_PROGRAMS =
dist_bin_SCRIPTS = find_path.py create_exon_fasta.py add_read_counts.py generate_gtf.py majority_vote.py nifflr.sh quantification.py fastqToFasta.pl
data_DATA =
BUILT_SOURCES = $(YAGGO_SOURCES)
CLEANFILES =
DISTCLEANFILES = $(BUILT_SOURCES)

# Yaggo automatic rules with silencing
V_YAGGO = $(V_YAGGO_$(V))
V_YAGGO_ = $(V_YAGGO_$(AM_DEFAULT_VERBOSITY))
V_YAGGO_0 = @echo "  YAGGO   " $@;
.yaggo.hpp:
#	$(V_YAGGO)$(YAGGO) --license $(srcdir)/header-license -o $@ $<
	$(V_YAGGO)$(YAGGO) -o $@ $<

YAGGO_SOURCES = # Append all file to be built by yaggo


# # Print the value of a variable
# print-%:
# 	@echo -n $($*)

