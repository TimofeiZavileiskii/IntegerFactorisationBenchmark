CC = g++
CFLAGS = -Wall -g0 -Os
INCLUDES = -lgmp 
GPU_FLAG = -arch=sm_75

all: generate_benchmark factorise_integers test


generate_benchmark: GenerateBenchmark.o GeneratePrimes.o
	$(CC) $(CFLAGS) -o generate_benchmark GenerateBenchmark.o GeneratePrimes.o $(INCLUDES)


test: Test.o GeneratePrimes.o
	$(CC) $(CFLAGS) -o test Test.o GeneratePrimes.o $(INCLUDES)


Test.o: Test.cpp EllipticCurves.h
	$(CC) $(CFLAGS) -c Test.cpp $(INCLUDES)


GenerateBenchmark.o: GeneratePrimes.o
	$(CC) $(CFLAGS) -c GenerateBenchmark.cpp $(INCLUDES)


GeneratePrimes.o: GeneratePrimes.cpp GeneratePrimes.h
	$(CC) $(CFLAGS) -c GeneratePrimes.cpp $(INCLUDES)


factorise_integers: IntegerFactorisationBenchmark.o PollardsRho.o TrialDivision.o TrialDivision.o PollardsP1.o ECM.o TrialDivisionCuda.o Utils.o
	$(CC) $(CFLAGS) -o factorise_integers IntegerFactorisationBenchmark.o PollardsRho.o TrialDivision.o PollardsP1.o ECM.o Utils.o $(INCLUDES)


IntegerFactorisationBenchmark.o: IntegerFactorisationBenchmark.cpp PollardsRho.h TrialDivision.h
	$(CC) $(CFLAGS) -c IntegerFactorisationBenchmark.cpp $(INCLUDES)


PollardsRho.o: PollardsRho.h PollardsRho.cpp Utils.o Utils.h
	$(CC) $(CFLAGS) -c PollardsRho.cpp $(INCLUDES)


PollardsP1.o: PollardsP1.cpp PollardsP1.h Utils.o Utils.h
	$(CC) $(CFLAGS) -c PollardsP1.cpp $(INCLUDES)


ECM.o: ECM.cpp ECM.h Utils.h Utils.o
	$(CC) $(CFLAGS) -c ECM.cpp $(INCLUDES)


TrialDivision.o: TrialDivision.h TrialDivision.cpp
	$(CC) $(CFLAGS) -c TrialDivision.cpp $(INCLUDES)


TrialDivisionCuda.o: TrialDivisionCuda.h TrialDivisionCuda.cu
	nvcc $(GPU_FLAG) -c TrialDivisionCuda.cu


Utils.o: Utils.h Utils.cpp
	$(CC) $(CFLAGS) -c Utils.cpp $(INCLUDES)
