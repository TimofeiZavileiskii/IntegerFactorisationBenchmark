CC = g++
CFLAGS = -Wall -g0 -Os -pg

all: factorise_integers

factorise_integers: IntegerFactorisationBenchmark.o PollardsRho.o TrialDivision.o TrialDivision.o PollardsP1.o ECM.o
	$(CC) $(CFLAGS) -o factorise_integers IntegerFactorisationBenchmark.o PollardsRho.o TrialDivision.o PollardsP1.o ECM.o Utils.o -lgmp


IntegerFactorisationBenchmark.o: IntegerFactorisationBenchmark.cpp PollardsRho.h TrialDivision.h
	$(CC) $(CFLAGS) -c IntegerFactorisationBenchmark.cpp -lgmp


PollardsRho.o: PollardsRho.h PollardsRho.cpp
	$(CC) $(CFLAGS) -c PollardsRho.cpp -lgmp


PollardsP1.o: PollardsP1.cpp PollardsP1.h Utils.o Utils.h
	$(CC) $(CFLAGS) -c PollardsP1.cpp -lgmp


ECM.o: ECM.cpp ECM.h Utils.h Utils.o
	$(CC) $(CFLAGS) -c ECM.cpp -lgmp


TrialDivision.o: TrialDivision.h TrialDivision.cpp
	$(CC) $(CFLAGS) -c TrialDivision.cpp -lgmp


Utils.o: Utils.h Utils.cpp
	$(CC) $(CFLAGS) -c Utils.cpp -lgmp