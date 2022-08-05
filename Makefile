WARN =-Wall -Wsign-compare -Wextra
CPP=g++ -O -g ${WARN} 

hpcc:
	ssh -X xliu261@cluster.hpcc.ucr.edu



check: checkerboard.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
	${CPP}  `pkg-config --cflags itpp` -o checkerboard checkerboard.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

test: test.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
	${CPP}  `pkg-config --cflags itpp` -o test test.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

test2: test.cpp 
	${CPP}  `pkg-config --cflags itpp` -o test test.cpp   `pkg-config --libs itpp`
 

gen_Hx:gen_surface.cpp
	${CPP}  -o gen_surface_star gen_surface.cpp

hypergraph: gen_HPC.cpp BP_basic_func1.cpp
	${CPP}  `pkg-config --cflags itpp` -o HPC gen_HPC.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`
 
qtest1: quantum_test1.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
	${CPP}  `pkg-config --cflags itpp` -o qtest1 quantum_test1.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

cla_BP:BP.h BP.cpp classical_main.cpp
	${CPP}  `pkg-config --cflags itpp` -o classical_BP BP.h BP.cpp classical_main.cpp  `pkg-config --libs itpp`

quan_BP:quantum_BP.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
		${CPP}  `pkg-config --cflags itpp` -o qBP quantum_BP.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

ana_quan_BP:trapping_sets_ana.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
		${CPP}  `pkg-config --cflags itpp` -o aqBP trapping_sets_ana.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

ana_quan_BP_try_again:trapping_sets_ana.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
		${CPP}  `pkg-config --cflags itpp` -o aqBP_try_again trapping_sets_ana.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

ana_quan_BP_print_mes:trapping_sets_ana.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
		${CPP}  `pkg-config --cflags itpp` -o aqBP_print_mes trapping_sets_ana.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

xquan_BP:quantum_BP_x_only.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
		${CPP}  `pkg-config --cflags itpp` -o xqBP quantum_BP_x_only.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

hquan_BP:quantum_BP_x_only_for_hpcc.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
		${CPP} -std=c++11 -fopenmp  `pkg-config --cflags itpp` -o hqBP quantum_BP_x_only.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

xquan_BP2:quantum_BP_x_only.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
		${CPP}  `pkg-config --cflags itpp` -o xqBP2 quantum_BP_x_only.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

gen_cyc: gen_cyclic.cpp  BP_basic_func1.cpp
	${CPP}  `pkg-config --cflags itpp` -o gen_cyc gen_cyclic.cpp  BP_basic_func1.cpp  `pkg-config --libs itpp`

tran: HT.cpp  BP_basic_func1.cpp
	${CPP}  `pkg-config --cflags itpp` -o tran HT.cpp  BP_basic_func1.cpp  `pkg-config --libs itpp`

display: BP_quantum_func1.cpp display_mat.cpp BP_basic_func1.cpp
	${CPP}  `pkg-config --cflags itpp` -o display  BP_quantum_func1.cpp display_mat.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`
