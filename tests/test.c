#include "gtest/gtest.h"

extern "C"
{
  #include <utils_string.h>
  #include <cafe_shell.h>
}

class SplitTest : public ::testing::Test {
	char c[10];
protected:
	virtual void SetUp() {
		strcpy(c, "a b");
		pArray = string_pchar_space_split(c);
	}

	pArrayList pArray;
};

char c[100];

// Tests factorial of 0.
TEST_F(SplitTest, SplitsOnSpaces) {
  EXPECT_EQ(2, pArray->size);
  EXPECT_STREQ("a", (char *)(pArray->array[0]));
  EXPECT_STREQ("b", (char *)(pArray->array[1]));
}

TEST(CommandEchoTest, Echoes) {
	EXPECT_EQ(0, cafe_shell_dispatch_command(c));
}

int main(int argc, char **argv) {
	cafe_shell_init();
	strcpy(c, "echo fish");
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
