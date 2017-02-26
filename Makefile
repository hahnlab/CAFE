#
# Compiler flags
#
CC     = gcc
VPATH = cafe:libcommon:libtree:tests
INCLUDE:=-I cafe -I libtree -I libcommon
CFLAGS = -Wall -std=c11 $(INCLUDE)
CXXFLAGS = -Wall $(INCLUDE)
LINKFLAGS = -lpthread -lm
TESTLIBS=-lCppUTest -lCppUTestExt 

#
# Project files
#
CSRCS=cafe_family.c cafe_main.c cafe_report.c cafe_tree.c cafe_shell.c birthdeath.c chooseln_cache.c phylogeny.c tree.c fminsearch.c grpcmp.c histogram.c  matrix_exponential.c regexpress.c utils_string.c gmatrix.c hashtable.c mathfunc.c memalloc.c utils.c
CXXSRCS=branch_cutting.cpp cafe_commands.cpp conditional_distribution.cpp \
        error_model.cpp Globals.cpp lambda.cpp log_buffer.cpp reports.cpp likelihood_ratio.cpp pvalue.cpp simerror.cpp
TESTSRCS=command_tests.cpp error_model_tests.cpp family_tests.cpp lambda_tests.cpp test.cpp
COBJS=$(CSRCS:.c=.o)
CXXOBJS=$(CXXSRCS:.cpp=.o)
SRCS=$(CSRCS) $(CXXSRCS)
OBJS=$(COBJS) $(CXXOBJS)
EXE  = cafe

#
# Debug build settings
#
DBGDIR = debug
DBGEXE = $(DBGDIR)/$(EXE)
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))
DBGCFLAGS = -g -O0 -DDEBUG -DVERBOSE

#
# Release build settings
#
RELDIR = release
RELEXE = $(RELDIR)/$(EXE)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))
RELCFLAGS = -O3 -DNDEBUG

#
# Test build settings
#
TESTDIR = test
TESTEXE = $(TESTDIR)/runtests
TESTOBJS=$(addprefix $(TESTDIR)/, $(TESTSRCS:.cpp=.o)) $(addprefix $(TESTDIR)/, $(OBJS))
TESTCFLAGS = -g -rdynamic -O0 -DDEBUG
TESTCPPFLAGS = -g -rdynamic -O0 -DDEBUG

.PHONY: all clean debug prep release remake test

# Default build
all: prep release

#
# Debug rules
#
debug: $(DBGEXE)

$(DBGEXE): $(DBGOBJS) main.o
	$(CXX) $(CXXFLAGS) $(DBGCFLAGS) -o $(DBGEXE) $^ $(LINKFLAGS) 

$(DBGDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(DBGCFLAGS) -o $@ $<

$(DBGDIR)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(DBGCFLAGS) -o $@ $<

$(DBGDIR)/main.o: main.cpp
	$(CXX) -c $(CXXFLAGS) $(DBGCFLAGS) -o $@ $<

#
# Test rules
#
test: $(TESTEXE)

$(TESTEXE) : $(TESTOBJS)
	$(CXX) $(CXXFLAGS) $(TESTCFLAGS) -o $(TESTEXE) $^ $(LINKFLAGS) $(TESTLIBS)

$(TESTDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(TESTCFLAGS) -o $@ $<

$(TESTDIR)/%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $(TESTCPPFLAGS) -o $@ $<

#
# Release rules
#
release: $(RELEXE)

$(RELEXE): $(RELOBJS) main.o
	$(CXX) $(CXXFLAGS) $(RELCFLAGS) -o $(RELEXE) $^ $(LINKFLAGS) 

$(RELDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(RELCFLAGS) -o $@ $<

$(RELDIR)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(RELCFLAGS) -o $@ $<

$(RELDIR)/main.o: main.cpp
	$(CXX) -c $(CXXFLAGS) $(DBGCFLAGS) -o $@ $<

#
# Other rules
#
prep:
	@mkdir -p $(DBGDIR) $(RELDIR) $(TESTDIR)

remake: clean all

clean:
	rm -f $(RELEXE) $(RELOBJS) $(DBGEXE) $(DBGOBJS) $(TESTEXE) $(TESTOBJS)
