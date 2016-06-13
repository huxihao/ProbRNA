#!/bin/sh
TH=1     # number of threads (suggest 8)
MK=10     # maximum of K, fitting window size (suggest 20)
MW=200    # maximum of w, scanning window size (suggest 200)

if [ $# -eq 0 ] 
then
echo "Please input the case that you want to reproduce. eg: Produce.sh case1"
echo "Here we use "$TH" threads for speeding up. You may change it for your equipment."
exit 1
fi

# Compare PL, MP, MPL
if [ $1 = "case1" ] # case 1 is for synthetic data only and so computed seperately
then
python src/fit_pars.py WorkPath=work Action=compare_3m Paras="$MK"_1000
mv work/log.txt case1.txt
fi

# Fit all PARS data with MP, MPL, MPLC-same, MPLC-oppo
if [ $1 = "all" -o $1 = "case2" ]
then
python src/fit_pars.py WorkPath=work Action=combine_fit Paras="$TH"_AllModel_"$MK"
mv work/log.txt case2.txt
fi

# Evaluate fitted PARS data on zipcode regions
if [ $1 = "all" -o $1 = "case3" ]
then
python src/learn_domain.py WorkPath=work Action=figure3 Paras="$TH"_40_"$MK"
mv work/log.txt case3a.txt
python src/learn_domain.py WorkPath=work Action=figure3 Paras="$TH"_100_"$MK"
mv work/log.txt case3b.txt
fi

# Evaluate individual features on zipcode regions
if [ $1 = "all" -o $1 = "case4" ]
then
python src/learn_domain.py WorkPath=work Action=figure4 Paras="$TH"_2_"$MW"
mv work/log.txt case4.txt
fi

# Evaluate combined features on zipcode regions
if [ $1 = "all" -o $1 = "case5" ]
then
python src/learn_domain.py WorkPath=work Action=figure6 Paras="$TH"_2_"$MW"
mv work/log.txt case5.txt
fi

# Kim data set one
if [ $1 = "all" -o $1 = "case6" ]
then
python src/validate_kim.py WorkPath=work Action=data_sizes Paras=gPAR_0_All
mv work/log.txt case6a.txt
python src/validate_kim.py WorkPath=work Action=choose_list Paras=gPAR_1000_All
python src/learn_domain.py WorkPath=work Action=figure5_permute_kim Paras="$TH"_2_"$MW"
mv work/log.txt case6b.txt
fi

# Kim data set two
if [ $1 = "all" -o $1 = "case7" ]
then
python src/validate_kim.py WorkPath=work Action=data_sizes Paras=Puf3p_0_All
mv work/log.txt case7a.txt
python src/validate_kim.py WorkPath=work Action=choose_list Paras=Puf3p_0_All
python src/learn_domain.py WorkPath=work Action=figure5_permute_kim Paras="$TH"_2_"$MW"
mv work/log.txt case7b.txt
fi

echo "End successfully!"

