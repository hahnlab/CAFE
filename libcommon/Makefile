OBJS:=$(patsubst %.c,../lib/%.o,$(wildcard *.c))
CC=gcc
OPTS=-g -O2 -I. -ffast-math -funroll-loops  # -pthread
LINKS=-lpthread -lm 
ARFILE=../lib/libcommon.a

all: $(OBJS)
	ar -r $(ARFILE) $(OBJS)

$(OBJS):../lib/%.o:%.c
	$(CC) -c -Wall $(OPTS) -D__DEBUG__ -o $@ $^

clean:
	@rm -rf $(OBJS) $(ARFILE)
