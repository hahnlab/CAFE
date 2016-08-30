test:
	g++ -I${GTEST_DIR}/include -I../libcommon -I../cafe -I../libtree -c test.c
	g++ -L../lib -o runtests test.o -lcafe -ltree -lcommon ${GTEST_DIR}/libgtest.a -pthread
