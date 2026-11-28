# Sanitizer stress test for the C core:  make -C tests/c run
SRC = ../../src/lzrans
CFLAGS = -O1 -g -std=gnu11 -Wall -Wextra -pthread -fsanitize=address,undefined -fno-sanitize-recover=all -DLZR_PARANOID
fuzz_roundtrip: $(wildcard $(SRC)/*.h) fuzz_roundtrip.c $(SRC)/enc.c $(SRC)/dec.c $(SRC)/opt.c $(SRC)/crc32.c $(SRC)/huf.c $(SRC)/fse.c
	$(CC) $(CFLAGS) -o $@ $(filter %.c,$^) -lm
run: fuzz_roundtrip
	./fuzz_roundtrip 300
# Data-race check of the threaded match finder
tsan: fuzz_roundtrip.c $(SRC)/enc.c $(SRC)/dec.c $(SRC)/opt.c $(SRC)/crc32.c $(SRC)/huf.c $(SRC)/fse.c
	$(CC) -O1 -g -std=gnu11 -pthread -fsanitize=thread -o fuzz_tsan $^ -lm && ./fuzz_tsan 40
# Public API tests (test_api.c) under ASan/UBSan
test_api: $(wildcard $(SRC)/*.h) ../../include/lzr.h test_api.c $(SRC)/api.c $(SRC)/enc.c $(SRC)/dec.c $(SRC)/opt.c $(SRC)/crc32.c $(SRC)/huf.c $(SRC)/fse.c
	$(CC) $(CFLAGS) -o $@ $(filter %.c,$^) -lm
api: test_api
	./test_api
clean:
	rm -f fuzz_roundtrip fuzz_tsan test_api
.PHONY: run tsan api clean
