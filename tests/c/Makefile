# Sanitizer stress test for the C core:  make -C tests/c run
SRC = ../../src/lzrans
CFLAGS = -O1 -g -std=c11 -Wall -Wextra -fsanitize=address,undefined -fno-sanitize-recover=all
fuzz_roundtrip: fuzz_roundtrip.c $(SRC)/frame.c $(SRC)/lz77.c $(SRC)/rans.c $(SRC)/crc32.c
	$(CC) $(CFLAGS) -o $@ $^
run: fuzz_roundtrip
	./fuzz_roundtrip 300
clean:
	rm -f fuzz_roundtrip
.PHONY: run clean
