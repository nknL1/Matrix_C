all: s21_matrix.a test gcov_report
CC = -Wall -Werror -Wextra -std=c11
DEF = -o testing -lcheck -lm
MATR = s21_matrix.c
TEST = s21_test_matrix.c
clang_formats:
	clang-format -i *.c
	clang-format -i *.h

s21_matrix.a:
	gcc $(CC) -c $(MATR) -lm
	ar rcs s21_matrix.a *.o

leaks:
	gcc $(MATR) $(TEST) $(DEF)
	valgrind --tool=memcheck --leak-check=yes ./testing
	rm -rf testing

test:
	gcc $(MATR) $(TEST) $(DEF)
	./testing
	rm -rf testing

clean:
	rm -rf *.o *.a *.gch *.gcno *.gcda *.info report testing

gcov_report:
	gcc --coverage $(MATR) $(TEST) $(DEF)
	./testing
	lcov -t "test" -o test_report.info -c -d . 
	genhtml -o report test_report.info
