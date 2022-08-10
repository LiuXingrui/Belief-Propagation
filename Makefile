WARN =-Wall -Wsign-compare -Wextra
CPP=g++ -O -g ${WARN} 
 

quan_BP:BP.h BP.cpp quantum_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o quantum_BP BP.h BP.cpp quantum_main.cpp  `pkg-config --libs itpp`


cla_BP::BP.h BP.cpp classical_main.cpp
	${CPP}  `pkg-config --cflags itpp` -o classical_BP BP.h BP.cpp classical_main.cpp  `pkg-config --libs itpp`

