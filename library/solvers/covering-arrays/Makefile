CC ?= gcc
CFLAGS ?= -O2 -Wall -Wextra -std=c11 -march=native
ca_sa: ca_sa.c
	$(CC) $(CFLAGS) -o ca_sa ca_sa.c -lm
clean:
	rm -f ca_sa
.PHONY: clean
