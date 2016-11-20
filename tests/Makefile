test: ../cafe/cafe
	g++ -g -I../libcommon -I../cafe -I../libtree -c test.cpp lambda_tests.cpp family_tests.cpp command_tests.cpp
	g++ -g -L../lib -o runtests test.o lambda_tests.o family_tests.o command_tests.o -lcafe -ltree -lcommon -lCppUTest -lCppUTestExt -pthread
