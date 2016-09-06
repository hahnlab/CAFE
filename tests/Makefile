test:
	g++ -g -I../libcommon -I../cafe -I../libtree -c test.c
	g++ -g -L../lib -o runtests test.o -lcafe -ltree -lcommon -lCppUTest -lCppUTestExt -pthread
