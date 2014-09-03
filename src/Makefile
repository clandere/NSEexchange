# Makefile for CES library

CC = gcc
LIBS = -L/usr/lib/ -lm -lgsl -L/usr/lib/R/lib -lR
INCS = -I/usr/share/R/include -I /data/apps/R/3.0.1/lib64/R/include
#Flags for optimization
CFLAGS = -c -mtune=generic -O3 -Wall -fPIC $(INCS)
#Flags for debugging
#CFLAGS = -c -g -mtune=generic -Wall -fPIC $(INCS) 

OBJS = shared/iofuns_CES.o shared/seqfuns_CES.o shared/synfuns_CES.o

all: calc_eta_NSE_R_wrapper.o simulate_sequence.o calc_denom_exact.o calc_eta_moments.o $(OBJS)
	gcc -std=gnu99 -shared -o exchange_NSE.so calc_denom_exact.o calc_eta_moments.o simulate_sequence.o calc_eta_NSE_R_wrapper.o $(OBJS) $(LIBS)
	
calc_denom_exact.o: calc_denom_exact.c
	gcc -std=gnu99 $(CFLAGS) -pipe -c calc_denom_exact.c -o calc_denom_exact.o
	
calc_eta_moments.o: calc_eta_moments.c
	gcc -std=gnu99 $(CFLAGS) -pipe -c calc_eta_moments.c -o calc_eta_moments.o

simulate_sequence.o: simulate_sequence.c
	gcc -std=gnu99 $(CFLAGS) -pipe -c simulate_sequence.c -o simulate_sequence.o

calc_eta_NSE_R_wrapper.o: calc_eta_NSE_R_wrapper.c
	gcc -std=gnu99 $(CFLAGS) -pipe -c calc_eta_NSE_R_wrapper.c -o calc_eta_NSE_R_wrapper.o

shared/iofuns.o: shared/iofuns_CES.c shared/iofuns_CES.h shared/CES_struct.h
	$(CC) $(CFLAGS) -o shared/iofuns_CES.o shared/iofuns_CES.c $(LIBS)

shared/seqfuns_CES.o: shared/seqfuns_CES.c shared/seqfuns_CES.h shared/CES_struct.h
	$(CC) $(CFLAGS) -o shared/seqfuns_CES.o shared/seqfuns_CES.c $(LIBS)

shared/synfuns_CES.o: shared/synfuns_CES.c shared/synfuns_CES.h shared/CES_struct.h
	$(CC) $(CFLAGS) -o shared/synfuns_CES.o shared/synfuns_CES.c $(LIBS)
