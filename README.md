# CancerImage
Idea of use deeplearing Algorithm for detecting ealy cancer,use cfDNA methylation 。 
1.first use fastp merge R1 and R2，  
2、merged reads  mapping to the reference genome   
3、recognize the single reads methylation pattern，every 20 or more cpgs site per windows,this is the parameters of this programe.  
4、each Windows will have 100+ reads depth ,and select a median depth for the whole genome depth as the jeight of the matrix and the width is the number of cpg sites and channels number is the selected final DMR number   
5、use CNN algorithm to train and test the for cancer early detection   
6、evalute the performance of this algorithm.
