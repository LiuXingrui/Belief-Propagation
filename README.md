# Belief Propagation

BP.cpp, BP.h contains the functions used in the programs. classical_main.cpp is for classical BP, quantum_main.cpp is for quantum BP. 
To compile, type
```bash
make quan_BP
```  
and 
```bash
make cla_BP
```  

## Classicl BP
To run the program:  
```bash  
./classical_BP  <H_file> <data_file> <p>   <number of codewords>  <max  iterations> <debug>
```
###  Paremeters
H_file: the file stored the parity checks  
data_file: the file to store the output data  
p: bit error rate  
debug:  debug=0 for parallel schedule, debug=1 for serial  schedule  

### Data-file structure
  The data stored in data_file are  
  column: parameter  
  0:n   
  1:logical-error rate  
  2:p  
  3:average iterations for successful decoding
  4:number_of_suc_decoding   
  5:bit_error_rate_after_decoding  
  
 
 ## Quantum BP
 To run the program:  
```bash  
./quantum_BP <Hx_file> <Hz_file>  <pavg/wt> <range> <max number of failed decoding/number of cws to decode><lmax> <data_file><debug><channel><alpha><decode_pavg><decode_prange><lambda>
```

  quantum_BP will decode CSS code with error probability distributed between pmin=pavg*(1-range) and pmax=pavg*(1-range), the decoder will use  dec_pmin=dec_pavg*(1-dec_range) and dec_pmax=dec_pavg*(1-dec_range), so we can use different error model for the error channel and decoder.  If range>1, pmin will be set to 0.
  
  If <max number of failed decoding/number of cws to decode> is less than 0, then that is number of codewords to decode. If it is greater than 0, that is this max number of failed decoding, when the program produces such number of failed decodings, it stopped.  For example, -100 means the program will decode 100 codewords. 100 means the program will decode until decoding failed for 100 codewords.  
  
  if <pavg/wt> is greater or equal than 1, it is explained as the weight of the errors and pavg/dec_pavg will be setted to be wt/n, for example, if  <pavg/wt> =2, the program will only generate wt=2 errors.   
  
   channel is for choosing which error channel to use.
  
   channel==0: bsc, only decode z-errors(use Hx)  
   channel==1: depolarizing, pavg is the depolarizing rate  (have not added this function yet).
  

debug is for bitwise option.   
  debug%2==1 : if reach maximum iteratios, try another random probability distribution with the same pmin and pmax, decode it again.  Do this for 5 times.
  (debug/2)%2==1 : parallel  schedule  
   (debug/4)%2==1: print the real_e, output_e and residual e when reach max iterations and failed (if use OSD, printing will be after the OSD), only works for the checkerboard codes.  
   (debug/8)%2==1: use OSD after BP fails.  
   (debug/16)%2==1: use LLR_avg.   (don't use this now, it won't work.)     
   Tips:
   0001=1: try again  
   0010=2: para  
   0011=3: try again and para  
   0100=4: print mes  
   0101=5: print and try again  
   0110=6: print and  para
   0111=7: print and para and try again  
   1001=9: OSD and try again  
   1011=11: para and OSD and try again  

   alpha is for testing larger step size.  
   lambda won't work now, just set it=1 is ok.  
    
   
   ### Data-file structure
   The data stored in data_file are  n d fail_rate pavg/wt range avg_iter_for_suc num_of_suc_dec num_of_cws syn_fail max_fail syn_fail_rate max_fail_rate:  
   column:parameter
   0:n  
   1:d  (have not calculate it in this prog yet, so outputs for d are setted to -1)  
   2: logical-error rate  
   3:pavg/wt    
   4: r 
   5:average iterations for successful decoding
   6:number of successful decodings
   7: number of codewords
   8: number of syn_fail  (the failure such that output_e gives the right syndrome but failed)
   9:number of max_fail  (failure that reached maximum iterations)
   10:syn_fail_rate  (fail rate for such 2 failures)
   11: max_fail_rate  
   12: decode_p  
   13:decode_prange    
   14:alpha  
   15: number of OSD_suc  
   16:lambda  
   
   ## Parity check format

An example for the parity check matrix file:
  
3 2  
1 2   
2 3  
Which is the matrix:  
  1 1 0  
  0 1 1  
  
  
  another example:  
  5 3  
  1 2 3  
  2 3 4  
  3 4 5  
 is:  
 1 1 1 0 0  
 0 1 1 1 0  
 0 0 1 1 1  
 
 The file stored parity check matrix use a similar format as MayKay's alist format, the only difference is there is no transposed matrix stored here. The first line is n n-k, then shows the position of 1 in corresponding rows.(The first element is 1 rather than 0.)
 
M1, M2, M3... are some matrices from MayKay's website.  
Toric code example: toric_Hx_n10_cyc11 is Hx for d=10 toric code.  

## Modified BP
```bash
make MBP
./MBP  Hx_file Hz_file p  max_failed_cws/max_num_of_cws lmax data_file options channel num_of_row_reduction debug
```

## Estimate Maximum-likelyhood(ML) decoding success rate
```bash
make MLSR
./MLSR  Hx_file Hz_file p  num_large_wt_errors debug data_file
./MLSR ../S_BP/checkerboard_Hx_d3  ../S_BP/checkerboard_Hz_d3 0.8 300 0 data/ML_data.data
```
