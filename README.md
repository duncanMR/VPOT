# VPOT
VPOT - Variant Prioritisation Ordering Tool. (using Python3)

VPOT is a Python tool written to allow prioritisation of variants in ANNOVAR annotated VCF files. VPOT provides three functions for the purpose of speeding up variant discovery.
* 1 - priority tool
* 2 - gene filter
* 3 - samples filtering
                                                                                                                                
 Entry :                                                                                                                               
 
       python3 VPOT.py  -  will return a help screen                                                                                 
       
       python3 VPOT.py priority <location for output file+prefix> <file of input VCF files> - will create a default Prioritisation Parameter File, PPF
       
       python3 VPOT.py prioirty <location for output file+prefix> <file of input VCF files> <PPF>
       
       python3 VPOT.py genef <location for output file+prefix> <VPOT prioritiy output> <gene list>
       
       python3 VPOT.py samplef <location for output file+prefix> <VPOT prioritiy output> <sample selection file>
                                                                                                                                                                                         
       <file of input VCF files> format (one sample per line):                                                                                                                    
          location of VCF file<tab>sample id                                                                                                                                      
          eg: test_input/data/sample1.vcf	S01                                                                                                                                    
                                                                                                                                                                                  
       <gene list> format  (one gene per line):                                                                                                                                   
          eg: ACTC1                                                                                                                                                               
                                                                                                                                                                                  
       <sample selection file > format (based on pedigree ped file format):                                                                                                       
          eg  FAM_1	PATIENT1  PATIENT2  PATIENT3	1	2
              FAM_1	PATIENT2  ND1	      ND2     	1	1
                                                                                                                           
             file of samples_IDs with corresponding include/exclude setting based on the affected column in the ped file                                                          
                     -  PATIENT1  2  = affected status (variant IS in sample PATIENT)                                                                                                                  
                     -  PATIENT2  1  = unaffected status (variant IS NOT in sample PATIENT)                                                                                                              
                     - combination of these values will determine if a variant is maintained or not                                                                               
                     - for above case, a variant is maintain if it is found in PATIENT1 and not in PATIENT2.                                                                             
                     - Note: if there are more samples than the ones stated, then they do not influence the variant selection.                                                
                                                                                                                                                                                  
## SETTING UP THE VPOT PARAMETER FILE FOR OPTION PRIORITY                                                                                                                                  

### 1. Setting up PF population filter in parameter file

 To provide population frequency threshold for selection of variants.

 FORMAT :    PF	<population frequency annotation name in VCF>	Value

 Example :   PF	ExAC_ALL	0.01

 This will tell VPOT to use the ExAC_ALL annotation values as a variant filtering criteria. VPOT will return variants that are <= to the value given, in this case 0.01.
 Multiple PF lines can be provide if you want to filter based on a combination of population frequency datasets. Note it is a AND logical approach, so the return variant would have
 met all the PF criteria. 

### 2. Setting up PD predictors in parameter file 

 To provide a point value to the various categories a predictor might return for a variant VPOT allow the user to determine the point assigned to a specific score or prediction.
 VPOT allows for as many breakdown/levels of differentiation as the user wants or the predictor needs.

 FORMAT :   
 
 PD	<predictor annotation name in VCF>	<(A)lpha/(N)umeric prediction type>	[prediction value/score]	[VPOT value for prediction value/score] .....[repeat as many values as you need] 
 see example below or in parameter files in the test data folder. 

 Multiple PD lines can be provided if you want to use multiple predictors 

 ALPHA PREDICTORS 

 For predictors that return alphanumeric prediction categories the PD line will list each category and its assigned VPOT value.
 Example :
 MutationTaster predicts an alteration as one of four possible types:

 disease causing - i.e. probably deleterious - D 
 disease causing automatic - i.e. known to be deleterious, see section dbSNP / TGP / ClinVar / HGMD for details- A
 polymorphism - i.e. probably harmless - P 
 polymorphism automatic - i.e. known to be harmless, see section dbSNP / TGP / ClinVar / HGMD for details N

 PD	MutationTaster_pred	A	N	0	P	0	D	1	A	2	
 
 where N = 0
       P = 0
       D = 1
       A = 2

 NUMERIC PREDICTORS 

 For predictors that return a numeric score the PD line will list scoring threshold and its assigned VPOT value.
 Example :
 CADD predicts an alteration with a score range from 0-30+ :

 CADD > 20 has been recommend by some papers as a good pathogenicity threshold
 CADD between 10 and 20 - median value for all possible canonical splice site changes and non-synonymous variants, as stated by CADD 
 CADD < 10 - probably harmless 

 Example :   PD	CADD_phred	N	10	0	20	1	20	2
 
 where x < 10 = 0
 
       10 < x < 20 = 1
       
       x > 20 = 2

 so the parameter line will contain
 S1|V1|S2|V2|S3|V3|.......|Sn|Vn|Sn+1|Vn+1
 
 where S = the predictor score
 
       V = VPOT value
       
      Sn = Sn+1
      
 and the method is :
 
        x < S1 then variant is assigned V1 
        
   S1 < x < S2 then variant is assigned V2 
   
   S2 < x < S3 then variant is assigned V3 
   
  ....... 
  
 Sn-1 < x < Sn then variant is assigned Vn 
 
        x > Sn+1 then variant is assigned Vn+1 

### 3. Setting up VT annotation in parameter file 

 To provide a point value to certain variant types that might not be well covered by predictors, eg STOPGAIN/SPLICING. 
 This allow highlighting of vertain variant types.

 FORMAT :    VT	<annotation field name in VCF that holds the >	<the variant type>	Value

 Example :   VT	VARIANT_TYPE	exonic_stopgain_	50

 For variant type "exonic_stopgain" an extra 50 will be added to the variant's cumulative value.
 Multiple VT lines can be provided to provide different stratification of variants 


### 4. Setting up GN gene symbol in parameter file 


 To provide the annotation field that contains the gene the variant is located in. 

 FORMAT :    GN	<annotation field name in VCF that holds the gene name> 

 Example :   GN	Gene Symbol

 ONLY one field should be provided for this parameter option.

## SETTING UP THE GENE SELECTION FILE FOR OPTION GENEF                                                                                                                                
                                            
 The gene selection is based on a text file with a single gene name each line.
                                                                         
## SETTING UP THE SAMPLE SELECTION FILE FOR OPTION SAMPLEF                                                                                                                                
                                            
 The sample selection is based on the pedigree ped file format, where the affected column is used to determine the selection of a variant.                                        
                                            
 If you are selecting variants based on one sample, eg all variants for a sample from a multi-sample VCF, then you can just specify that single sample in the 
 sample slection file, eg only want all sample PATIENT1 variants :
 
             FAM_1	PATIENT1	PATIENT2	PATIENT3	1	2       
                                            
 If you want variants that appear in PATIENT1 but PATIENT2, a segregation request then you will need a line for each sample 
 
             FAM_1	PATIENT1	PATIENT2	PATIENT3	1	2
             
             FAM_1	PATIENT2	ND1	ND2	1	1
             
                     -  PATIENT1  2  = affected (variant IS in sample PATIENT1)
                     
                     -  PATIENT2  1  = unaffected (variant IS NOT in sample PATIENT2)
                     
 so a combination of these values will determine if a variant is maintained or not                                                                              
 for the above case, a variant is maintain if it is found in PATIENT1 and not in PATIENT2.                                                                             
 Note: if there are more samples than the ones stated, then they do not influence the variant selection.                                                                          

## see TUTORIAL in the test_data directory for test examples of each function and also test files.
