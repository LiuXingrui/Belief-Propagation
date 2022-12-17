WARN =-Wall -Wsign-compare -Wextra
CPP=g++ -O -g ${WARN} 

hpcc:
	ssh -X xliu261@cluster.hpcc.ucr.edu


MBP:BP.h BP.cpp modified_BP.h modified_BP.cpp modified_BP_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o MBP BP.h BP.cpp modified_BP.h modified_BP.cpp modified_BP_main.cpp `pkg-config --libs itpp`


MLSR:BP.h BP.cpp modified_BP.h modified_BP.cpp ML_suc_rate_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o MLSR BP.h BP.cpp modified_BP.h modified_BP.cpp ML_suc_rate_main.cpp `pkg-config --libs itpp`

check: checkerboard.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
	${CPP}  `pkg-config --cflags itpp` -o checkerboard checkerboard.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

test: test.cpp  BP.h BP.cpp
	${CPP}  `pkg-config --cflags itpp` -o test test.cpp BP.h BP.cpp  `pkg-config --libs itpp`

test2: test.cpp 
	${CPP}  `pkg-config --cflags itpp` -o test test.cpp   `pkg-config --libs itpp`

qtest1: quantum_test1.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
	${CPP}  `pkg-config --cflags itpp` -o qtest1 quantum_test1.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

cla_BP:BP.h BP.cpp classical_main.cpp
	${CPP}  `pkg-config --cflags itpp` -o classical_BP BP.h BP.cpp classical_main.cpp  `pkg-config --libs itpp`


quan_BP:BP.h BP.cpp quantum_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o quantum_BP2 BP.h BP.cpp quantum_main.cpp  `pkg-config --libs itpp`



