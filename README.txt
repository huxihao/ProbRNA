 
     T h e   R e a d M e   F i l e   F o r   P r o b R N A   P r o j e c t

                                 Created by Xihao Hu 
		                                  on Dec 13, 2013

 The source code is supporting the paper:
 > Hu et. al. "Computational identification of protein binding sites on
               RNA using high-throughput RNA structure-probing data"

 Working Environment:
    * CPU ~8*2.40GHz, RAM > 4G (suggested to run 8 threads to speed up)
    * Python 2.7.2 
    * R version 2.10.1
         - randomForest
         - ROCR
    * A work path with 10G free spaces, here we use ./work as an example
    
-> Please download the full PARS data set from "PARS_DATA_SET.rar (14.2M)" at

    http://yiplab.cse.cuhk.edu.hk/probrna/

 and then uncompress it to the work folder.

-> Test the functions of basic modules

	python src/fit_pars.py WorkPath=work Action=test
	python src/learn_domain.py WorkPath=work Action=test
	python src/validate_kim.py WorkPath=work Action=test

-> If running the above commands is successful from your shell, you can use the
 "Produce.sh" script to get some results with a small range of parameters:

	sh Produce.sh all &> log.txt &

-> Wait ~1 day to finish with a single processor

 The file "PredAll_MixPoiLin_2_PARS_V1_S1.csv" in the work folder saves the 
 fitting outcomes, where the 8 columns are: 

 <index>:	unique index for each RNA sequence
 <tag>:		bases with negative values were excluded from statistical analysis
 <seq>:		each bases on the RNA sequence
 <v1>:		the V1 read count from the PARS data set
 <s1>:		the S1 read count from the PARS data set
 <gene>:	gene name
 <pbv>:		probability of V1 read count belong to the cluster with high values 
 <pbs>:		probability of S1 read count belong to the cluster with high values 

 The <v1> and <s1> columns are the input data, and <pbv> and <pbs> are the 
 output data of the Mixture of Poisson Linear Model.

 The table format was designed to be suitable for R to read. Interested users
 can refer to "fit_pars_all.R" and "fit_pars_3m.R" for the optimization process.

-> You may need further change the parameters defined in the head of Produce.sh 
 to suggested values to reproduce *all* the results in mentioned paper:
 
#!/bin/sh
TH=8     # number of threads (suggest 8)
MK=20    # maximum of K, fitting window size (suggest 20)
MW=200   # maximum of w, scanning window size (suggest 200)

 Then, run again by: sh Produce.sh all &> log.txt &
 
-> It will take a few days using 8 processors. All results will be saved to the
 current folder by names of case*.txt. They are corresponding to the paper by:
 
 case1  -> Fig 1 and Fig S2
 case2  -> Fig 1 and Fig S2
 case3a -> Fig S3
 case3b -> Fig 3
 case4  -> Fig 4
 case5  -> Figure of AUCs from combined features
 case6a -> Tab S3
 case6b -> Fig 5
 case6b -> Table of nucleotide composition in Puf3p
 case7b -> Fig S6

 You can find the expected result files from the ./result folder.

 Any question can be sent to "Xihao Hu" <huxihao@gmail.com> 

//End of read me
