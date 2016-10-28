CC=gcc
OPTS= -O2 -ffast-math -funroll-loops -Wall -g -std=c11 # -pthread
CPP_OPTS= -O2 -ffast-math -funroll-loops -Wall -g 
INCLUDE:=-I . -I ../libtree -I ../libcommon
LINKS:=-lpthread -lm
LIBS:=../lib/libtree.a ../lib/libcommon.a
SRCS=cafe_family.c cafe_main.c cafe_report.c cafe_tree.c cafe_shell.c cafe_commands.cpp lambda.cpp reports.cpp
OBJS=cafe_family.o cafe_main.o cafe_report.o cafe_tree.o cafe_shell.o cafe_commands.o lambda.o reports.o
EXEC=cafe
COMMON=../libcommon/
TREE=../libtree/
ARFILE=../lib/libcafe.a

all: $(EXEC)

$(EXEC): $(LIBS) $(OBJS)
	ar -r $(ARFILE) $(OBJS)
	g++ main.cpp $(INCLUDE) $(CPP_OPTS) $(ARFILE) $(LIBS) $(LINKS) -o $@

$(LIBS):
	make -C $(COMMON)
	make -C $(TREE)

$(OBJS):
.c.o: $(LIBS)
	$(CC) -c $(INCLUDE) $(OPTS) -D__DEBUG__  $< -o $@
.cpp.o: $(LIBS)
	$(CC) -c $(INCLUDE) $(CPP_OPTS) -D__DEBUG__  $< -o $@

# what have I done wrong?
lint: $(SRCS)
	lint $(SRCS) 

# clean out the dross
clean:
	rm -rf cafe ../lib/* *~ *.o *.bak core cafe.dSYM

