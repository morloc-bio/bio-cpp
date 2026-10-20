.PHONY: all test clean

all: test

test:
	morloc make -o test test.loc
	./test test | tee test.out
	tail -n 1 test.out | grep -qx true
	./test testReference | tee test-ref.out
	tail -n 1 test-ref.out | grep -qx true

clean:
	rm -rf test test-build test.out test-ref.out
