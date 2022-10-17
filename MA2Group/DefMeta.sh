#!/bin/sh    
#DefMeta.sh    
conda activate myenv_x86
#merge SeqDatabase_2, MAData_1, MAData_3, MAData_4 
cat MA2_data.1.fasta SeqDatabase_2.fasta >SeqDatabase_6.1.fasta
cat MA2_data.3.fasta SeqDatabase_6.1.fasta > SeqDatabase_6.2.fasta
cat MA2_data.13.fasta SeqDatabase_6.2.fasta > SeqDatabase_6.0.fasta


      
#step 1 - Remove duplicates from the SeqDatabase_6.0.fasta  
awk '{gsub(/\r/,"")} />/&&list~substr($1,2){flag=0} />/&&list!~substr($1,2){n++;list=list", "substr($1,2);flag=1}flag==1{print}'<SeqDatabase_6.0.fasta >SeqDatabase_6.1.fasta
     
    
#step2 - generate metadata for the SeqDatabase_6.1 fasta file 
awk '/>/{print substr($0,2)}'<SeqDatabase_6.1.fasta| awk 'BEGIN{print "ACC,GenSp,SubSp,Strain1,Strain2,Gene,comments"}$4=="16S"{$4=",,"$4}$4=="strain"{$4=","$4}{$2=","$2; $3=substr($3,2)","}/strain/&&!/subsp/{$4=","$4;gsub(/ strain/,",strain")}!/strain/&&!/subsp/{$4=","$4;gsub(/ 16S ribosomal/,", 16S ribosomal")}/strain/&&/subsp/{$5=$5","} !/strain/&&/subsp/{$6=","$6; gsub(/ 16S ribosomal/,", 16S ribosomal")}{gsub(/ 16S ribosomal/,",16S ribosomal");gsub(/ strain/,",strain");gsub(/ partial/,"partial");gsub(/ complete/,"complete")}{print}' >SeqData6_meta.csv
                                                      
#step 3 - remove the excess information from the header of the SeqDatabase file
awk '/>/{print $1} !/>/&&$1!=""{print}'<SeqDatabase_6.1.fasta >SeqDatabase_6.fasta
                                                                                      
#step 4 - do a muscle alignment
muscle -super5 SeqDatabase_6.fasta -output SeqDatabase_6.aln.fasta

