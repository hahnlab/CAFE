OBJS:=$(patsubst %.c,../lib/%.o,$(wildcard *.c))
CC=gcc
OPTS= -O2 -I . -I ../libcommon -ffast-math -funroll-loops  # -pthread
ARFILE=../lib/libtree.a

all: $(OBJS)
	ar -r $(ARFILE) $(OBJS) 

$(OBJS):../lib/%.o:%.c
	$(CC) -c -Wall $(OPTS) -D__DEBUG__ -o $@ $^

clean:
	@rm -rf $(OBJS) $(ARFILE) 
