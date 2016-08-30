CC=gcc
OPTS= -O2 -ffast-math -funroll-loops -Wall -g -std=c99 # -pthread
INCLUDE:=-I . -I ../libtree -I ../libcommon
LINKS:=-lpthread -lm
LIBS:=../lib/libtree.a ../lib/libcommon.a
SRCS=cafe_cgi.c cafe_family.c cafe_main.c cafe_report.c cafe_tree.c cafe_shell.c
OBJS=$(SRCS:.c=.o)
EXEC=cafe
COMMON=../libcommon/
TREE=../libtree/
ARFILE=../lib/libcafe.a

all: $(EXEC)

$(EXEC): $(LIBS) $(OBJS)
	ar -r $(ARFILE) $(OBJS)
	$(CC) main.c $(INCLUDE) $(OPTS) $(ARFILE) $(LIBS) $(LINKS) -o $@

$(LIBS):
	make -C $(COMMON)
	make -C $(TREE)

$(OBJS):
.c.o: $(LIBS)
	$(CC) -c $(INCLUDE) $(OPTS) -D__DEBUG__  $< -o $@

# what have I done wrong?
lint: $(SRCS)
	lint $(SRCS) 

# clean out the dross
clean:
	rm -rf cafe ../lib/* *~ *.o *.bak core cafe.dSYM

