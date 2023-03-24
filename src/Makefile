#CFLAGS=-Wall -Werror -Wextra
CFLAGS = -Wall -Werror -Wextra -std=c11

all: clean s21_math.a test gcov_report

ubuntu: clean s21_math.a test_ubuntu gcov_report

test: main.c s21_math.c
	gcc -o test $(CFLAGS) --coverage main.c s21_math.c -lcheck -lpthread -lm  -I /usr/local/include -L/usr/local/lib
	
test_ubuntu: main.c s21_math.c
	gcc -o test $(CFLAGS) --coverage main.c s21_math.c -lcheck -lm -lpthread -lsubunit
	
s21_math.a: s21_math.o
	ar rc s21_math.a s21_math.o
	ranlib s21_math.a
	

s21_math.o: s21_math.c
	gcc $(CFLAGS) -c s21_math.c
	
gcov_report: test
	./test
	gcovr --html-details rep.html
	

check:
	cp ../materials/linters/.clang-format ../src/.clang-format
	clang-format -i *.c
	clang-format -i *.h
	

clean:
	rm -rf *.o *.gcda *.a *.html *.gcno *.css
	rm -f test
	rm -f test_ubuntu
	
git:
	git add *.c
	git add *.h
	git add Makefile
	#git commit -m "add from make"
	#git push
