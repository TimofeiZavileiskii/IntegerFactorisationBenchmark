CC = g++
CFLAGS = -Wall -Ofast -g -I ./extern/cxxopts/include -L/usr/local/cuda/lib64 -flto
INCLUDES = -lgmp -lcudart -lpari
GPU_FLAG = -arch=sm_75 -lgmp -O3


all: generate_benchmark factorise_integers test


generate_benchmark: GenerateBenchmark.o GeneratePrimes.o GeneratePrimes.cpp GeneratePrimes.h GenerateBenchmark.cpp
	$(CC) $(CFLAGS) -o generate_benchmark GenerateBenchmark.o GeneratePrimes.o $(INCLUDES)


test: Test.o GeneratePrimes.o Test.cpp EllipticCurves.h Utils.h
	$(CC) $(CFLAGS) -o test Test.o GeneratePrimes.o Utils.o $(INCLUDES)


Test.o: Test.cpp EllipticCurves.h Utils.h
	$(CC) $(CFLAGS) -c Test.cpp $(INCLUDES)


GenerateBenchmark.o: GeneratePrimes.o
	$(CC) $(CFLAGS) -c GenerateBenchmark.cpp $(INCLUDES)


GeneratePrimes.o: GeneratePrimes.cpp GeneratePrimes.h
	$(CC) $(CFLAGS) -c GeneratePrimes.cpp $(INCLUDES)


factorise_integers: IntegerFactorisationBenchmark.o PollardsRho.o TrialDivision.o TrialDivision.o PollardsP1.o ECM.o TrialDivisionCuda.o Utils.o PollardsRhoCuda.o FactorisationStats.o ECMCuda.o
	$(CC) $(CFLAGS) -o factorise_integers IntegerFactorisationBenchmark.o PollardsRho.o TrialDivision.o PollardsP1.o ECM.o Utils.o PollardsRhoCuda.o TrialDivisionCuda.o FactorisationStats.o ECMCuda.o $(INCLUDES)


IntegerFactorisationBenchmark.o: IntegerFactorisationBenchmark.cpp PollardsRho.h TrialDivision.h PollardsRhoCuda.h
	$(CC) $(CFLAGS) -c IntegerFactorisationBenchmark.cpp $(INCLUDES)


PollardsRho.o: PollardsRho.h PollardsRho.cpp Utils.o Utils.h
	$(CC) $(CFLAGS) -c PollardsRho.cpp $(INCLUDES)


PollardsRhoCuda.o: PollardsRhoCuda.cu
	nvcc $(GPU_FLAG) -c PollardsRhoCuda.cu


PollardsP1.o: PollardsP1.cpp PollardsP1.h Utils.o Utils.h
	$(CC) $(CFLAGS) -c PollardsP1.cpp $(INCLUDES)


ECM.o: ECM.cpp ECM.h Utils.h Utils.o
	$(CC) $(CFLAGS) -c ECM.cpp $(INCLUDES)


ECMCuda.o: ECMCuda.cu ECMCuda.h
	nvcc $(GPU_FLAG) -c ECMCuda.cu


TrialDivision.o: TrialDivision.h TrialDivision.cpp
	$(CC) $(CFLAGS) -c TrialDivision.cpp $(INCLUDES)


TrialDivisionCuda.o: TrialDivisionCuda.h TrialDivisionCuda.cu
	nvcc $(GPU_FLAG) -c TrialDivisionCuda.cu


FactorisationStats.o: FactorisationStats.cpp FactorisationStats.h
	$(CC) $(CFLAGS) -c FactorisationStats.cpp $(INCLUDES)


Utils.o: Utils.h Utils.cpp
	$(CC) $(CFLAGS) -c Utils.cpp $(INCLUDES)


clean:
	rm -f ./*.o factorise_integers test generate_benchmark
