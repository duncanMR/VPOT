###########################################################################################################
# #
###########################################################################################################
#
import sys, re, glob, os, subprocess, time #
import numpy as np #
from shutil import copyfile #
#
tab='\t' # 
nl='\n' #
#
###########################################################################################################
# Define global variables
##########################################################################################################
#
def init(): ##
    ######################################
    #additions for gene panel function:
	global cancer_type #
	cancer_type="" #
	global final_output_file_list #
	final_output_file_list=[] #
	global panel_gene_list #
	panel_gene_list=[] #
	global panel_name_list #
	panel_name_list=[] #
	global n_panels #
	n_panels=0 #
    ######################################
	global VPOT_option #
	VPOT_option=0 #
	global input_file #
	input_file="" #
	global output_dir #
	output_dir="" #
	global final_output_file #
	final_output_file="" #
	global ufunction #
	ufunction="" #
	global parameter_file #
	parameter_file="" #
	global gene_list #
	gene_list="" #
	global Population #
	Population = ["ExAC","1000g"] #
	global Exonic #
	Exonic = ["ExonicFunc","VARIANT_TYPE"] #
	global Inheritance_model #
	Inheritance_model = ["DN","AR","AD","CH"] #
	global chromosome_list #
	chromosome_list = ["0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","M","MT"] #
	global MAFval #
	MAFval="0.01" #
	global PF #
	PF="PF" #
	global PD #
	PD="PD" #
	global VT #
	VT="VT" #
	global GN #
	GN="GN" #
	global VS #
	VS="VS" #
	global QC #
	QC="QC" #
	global QC_PASS #
	QC_PASS=False #
	global startlen #
	startlen=0 #
	global PD_type #
	PD_type=2 #
	global PD_low #
	PD_low=3 #
	global PD_mid #
	PD_mid=5 #
	global PD_high #
	PD_high=7 #
	global Maxval # allow for max of 28 characters values for an alpha predictor.
	Maxval=28 # allow for max of 28 characters values for an alpha predictor.
	global CV #
	CV="CV" #
	global Non_alt_GT_types #
	Non_alt_GT_types = ["0","."] #
	#
	global Maxcoverage # read coverage of sample must be equal or greater then this to included 
	Maxcoverage=0 # read coverage of sample must be equal or greater then this to included 
	#Maxcoverage=10 # read coverage of sample must be equal or greater then this to included 
	global Hete_Balance # Allele balance for heteozygous call, must be equal or greater then this to accept as hete, else is homozygous ref 
	Hete_Balance=0 # Allele balance for heteozygous call, must be equal or greater then this to accept as hete, else is homozygous ref 
	global Genotype_Quality # Genotype quality score must be equal or greater then this to include 
	Genotype_Quality=0 #  Genotype quality score must be equal or greater then this to include 
	global VariantScoreThreshold # Variant score threshold - do not include variant in VPOL less than or equal to this value
	VariantScoreThreshold=0 # Variant score threshold - do not include variant in VPOL less than or equal to this value 
	global VariantPercentageThreshold # Variant percentage threshold - do not include variant in VPOL below or equal to this percentage value
	VariantPercentageThreshold=0 # Variant percentage threshold - do not include variant in VPOL below or equal to this percentage value 
	global sample_coverage # VCF format codes to look out for 
	global sample_coverage_loc # location of VCF format codes for sample 
	sample_coverage = ["GT","NR","NV","DP","AD","GQ"] # VCF format codes to look out for 
	sample_coverage_loc = [-1,-1,-1,-1,-1,-1] # location of VCF format codes for sample 
#	sample_coverage = ["GT","NR","NV","DP","AD"] # VCF format codes to look out for 
#	sample_coverage_loc = [-1,-1,-1,-1,-1] # location of VCF format codes for sample 
	#
	global GT_val #
	global NR_val #
	global NV_val #
	global DP_val #
	global AD_val #
	global GQ_val #
	GT_val=0 #
	NR_val=1 #
	NV_val=2 #
	DP_val=3 #
	AD_val=4 #
	GQ_val=5 #
#
	global pop_array #
	global pred_array #
	global PF_array #
	global PD_array #
	global VT_array #
	global GN_value #
	global GENE_NAME #
	global Gene_ref #
	global Sample_ids #
	pop_array=[] #
	pred_array=[] #
	PF_array=[] #
	PD_array=[] #
	VT_array=[] #
	GN_value="" #
	GENE_NAME="GENE_NAME" #
	Gene_ref="NONE" #
	Sample_ids=[] #
# 
	global input_type_gVCF #
	input_type_gVCF=False #
	global txt_start #
	global txt_end# 
	txt_start=-1 #
	txt_end=-1 # 
	global sample_loc #
	global sample_ID #
	global INFO_loc #
	global FORMAT_loc #
	sample_loc=-1 #
	sample_ID="" #
	INFO_loc=-1 #
	FORMAT_loc=-1 #
#GT_loc=-1 #
	global header_ln #
	global blank_variant_ln #
	global inh_model #
	inh_model="NONE" #
	header_ln="" #
	blank_variant_ln="" #
#
	global temp_output_file #
	global working_file1 #
	global working_file2 #
	global working_file3 #
	global working_file4 #
	global full_file1 #
	global full_file2 #
	global sort_file1 #
	global sort_file2 #
	global sort_file3 #
	temp_output_file="" #
	working_file1="" #
	working_file2="" #
	working_file3="" #
	working_file4="" #
	full_file1="" #
	full_file2="" #
	sort_file1="" #
	sort_file2="" #
	sort_file3="" #
# output files
###########################################################################################################
#
###########################################################################################################
# determine if input field is a numeric field                                                             #
###########################################################################################################
def is_number(s):
	try:
		float(s)
#		print "float :",s #
		return True 
	except ValueError :
#		print "not float :",s #
		return False
#
###########################################################################################################
#
###########################################################################################################
def update_existing_variants(work_file,value): #
#
#	print "update_existing_variants(): " #
##
	with open(work_file,'r',encoding="utf-8") as samples_file, open(full_file2,'w',encoding="utf-8") as final_file : # 
#ed		outline1 = header_ln+nl #
#		print "header line )",outline1 #
#		print "test blank new variant line )",outline2 #
#		print "blank variant line splited )",variant_line #
#ed		final_file.write(outline1) # write new header line with new sample ID to the new combined samples output file
		#
		for line2 in samples_file: # work each line of current samples variant file 
#			print "Line in combined samples file: ",line2 #
			line2_array=re.split('\t|\n|\r',line2) # split into file location and sample id
			if (line2_array[0] != "#CHROM"): # not header line
				if (value == "1"): # 1st file then no need for any addition
					outline3='\t'.join(line2_array[:-1])+nl #
				else : # not first - then just do the zero
					outline3='\t'.join(line2_array[:-1])+tab+value+nl #
#				print "new variant : ", outline3 #
				final_file.write(outline3) #
			#
#				
	
#
###########################################################################################################
#
###########################################################################################################
def incorporate_this_src_into_full_file(): #
#
#	print "incorporate_this_src_into_full_file():" #
# working_file1 = new sample's variants
# full_file1 = current combined samples' variants
# full_file2 = the new combined samples' variants file - will contain new sample's variants 
#
	new_sample_not_finished=True #
	current_sample_not_finished=True #
#
#	COMMAND="sort -V -k1,5 "+working_file1+" | grep -v \"CHROM\" > "+sort_file1 #  
#	COMMAND="sort -V -k1,5 "+working_file1+" > "+sort_file1 #  
#	COMMAND="sort -V -u -k1,5 "+working_file1+" > "+sort_file1 #  
#	COMMAND="sort -u -k1,2V -k3,5 "+working_file1+" > "+sort_file1 #  
	COMMAND="sort -u -k1,1V -k2,2g -k4,5 "+working_file1+" > "+sort_file1 #  
	subprocess.call(COMMAND, shell=True) #
	copyfile(sort_file1,working_file1) # copy back
#	sample_output_file=output_dir+"sample_"+sample_ID+"_"+suffix+".txt" #
#	copyfile(sort_file1,sample_output_file) # copy back
#
#	COMMAND="sort -V -k1,5 "+full_file1+" | grep -v \"CHROM\" > "+sort_file2 #  
#	COMMAND="sort -V -k1,5 "+full_file1+" > "+sort_file2 #  
#	COMMAND="sort -k1,2V -k3,5 "+full_file1+" > "+sort_file2 #  
	COMMAND="sort -k1,1V -k2,2g -k4,5 "+full_file1+" > "+sort_file2 #  
	subprocess.call(COMMAND, shell=True) #
	copyfile(sort_file2,full_file1) # copy back
#
	#	print working_file1 #
	with open(working_file1,'r',encoding="utf-8") as new_sample, open(full_file1,'r',encoding="utf-8") as current_sample, open(full_file2,'w',encoding="utf-8") as final_file : # 
#ed		outline1 = header_ln+nl #
		variant_line=re.split('\t|\n',blank_variant_ln) # split up the blank variant line 
#ed		final_file.write(outline1) # write new header line with new sample ID to the new combined samples output file
		#
		S1=new_sample.readline() # 
		C1=current_sample.readline() #
#		print "S1:", S1 #
#		print "C1:", C1 # 
		while (new_sample_not_finished and current_sample_not_finished): # while we are working the two files
			line1_array=re.split('\t|\n|\r',S1) # split the variant line
			line2_array=re.split('\t|\n|\r',C1) # split the variant line
			new_sample_len=len(line1_array) #
			current_sample_len=len(line2_array) #
#
#			print ("new variant - S1 : ",line1_array[:5]) #
#			print ("variant already in file C1 : ",line2_array[:5]) # same variant
			if (len(line1_array[0]) > 3): # have chr prefix
				wrk_chr=line1_array[0]
				new_chr=wrk_chr[3:] #
				wrk_chr=line2_array[0] #
				existing_chr=wrk_chr[3:] #
			else:
				new_chr=line1_array[0] #
				existing_chr=line2_array[0] #
				
			if (new_chr == "X"):
				new_chr=23
			elif (new_chr == "Y"):
				new_chr=24
			elif ((new_chr == "M") or (new_chr == "MT")):
				new_chr=25
			else:
				new_chr=int(new_chr)
#
			if (existing_chr == "X"):
				existing_chr=23
			elif (existing_chr == "Y"):
				existing_chr=24
			elif ((existing_chr == "M") or (existing_chr == "MT")):
				existing_chr=25
			else:
				existing_chr=int(existing_chr)
#
#			print ("new chr - ",str(new_chr))
#			print ("existing chr - ",str(existing_chr))
#
#			if (line1_array[:5] == line2_array[:5]): # same variant
#			if (line1_array[0] == line2_array[0]): # same chr
			if (new_chr == existing_chr): # same chr
				if (int(line1_array[1]) == int(line2_array[1])) :  # same POS (left flank)
					if (line1_array[3] == line2_array[3]): # same REF
						if (line1_array[4] == line2_array[4]): # same ALT - so the same variants
#							outline3='\t'.join(line2_array[:-1])+tab+"1"+nl # 
#							print ("this sample genotype : ",line1_array[new_sample_len-2])
							outline3='\t'.join(line2_array[:-1])+tab+line1_array[new_sample_len-2]+nl # add in the new sample to the existing variant detail
							final_file.write(outline3) # 
#							print ("adding new sample to current variant 1 -",outline3) #
							S1=new_sample.readline() #
							C1=current_sample.readline() #
						elif (line1_array[4] < line2_array[4]): # not same variant and new variant is before current variant
#							new_sample_len=len(line1_array) #
							variant_line[:new_sample_len-2]=line1_array[:new_sample_len-2] # one less to exclude the genotype and \n 
#							variant_line[:new_sample_len-1]=line1_array[:new_sample_len-1] # one less to exclude the genotype and \n 
							outline3='\t'.join(variant_line[:-1])+tab+line1_array[new_sample_len-2]+nl # replace the last 0 with a 1
#							outline3='\t'.join(variant_line[:-1])+tab+"1"+nl # replace the last 0 with a 1
							final_file.write(outline3) #
#							print ("saving new variant with new sample 1 -",outline3) #
							S1=new_sample.readline() #
						else: # (line1_array[4] > line2_array[4]): # this variant is after current
							outline3='\t'.join(line2_array[:-1])+tab+"0"+nl # add a 0 for current variant for sample
							final_file.write(outline3) # 
#							print ("saving current variant with a 0 sample 1 -",outline3) #
							C1=current_sample.readline() #
					elif (line1_array[3] < line2_array[3]): # different ref - earlier
#						new_sample_len=len(line1_array) #
						variant_line[:new_sample_len-2]=line1_array[:new_sample_len-2] # one less to exclude the genotype and \n 
#						variant_line[:new_sample_len-1]=line1_array[:new_sample_len-1] # one less to exclude the \n 
#						outline3='\t'.join(variant_line[:-1])+tab+"1"+nl # replace the last 0 with a 1
						outline3='\t'.join(variant_line[:-1])+tab+line1_array[new_sample_len-2]+nl # add in the new sample 
						final_file.write(outline3) #
#						print ("saving new variant and new sample 2 -",outline3) #
						S1=new_sample.readline() #
					else: # (line1_array[3] > line2_array[3]) #  saved variant ref > 
						outline3='\t'.join(line2_array[:-1])+tab+"0"+nl # 
						final_file.write(outline3) # 
#						print ("saving current variant already in file with 0 for new sample 2 -",outline3) #
						C1=current_sample.readline() #
				elif (int(line1_array[1]) < int(line2_array[1])) :  # earlier position 
#					new_sample_len=len(line1_array) #
					variant_line[:new_sample_len-2]=line1_array[:new_sample_len-2] # one less to exclude the genotype and \n 
#					variant_line[:new_sample_len-1]=line1_array[:new_sample_len-1] # one less to exclude the \n 
					outline3='\t'.join(variant_line[:-1])+tab+line1_array[new_sample_len-2]+nl # add the sample genotype
#					outline3='\t'.join(variant_line[:-1])+tab+"1"+nl # replace the last 0 with a 1
					final_file.write(outline3) #
#					print ("saving new variant and new sample 3-",outline3) #
					S1=new_sample.readline() #
				else: # (int(line1_array[1]) > int(line2_array[1])) :  # 
					outline3='\t'.join(line2_array[:-1])+tab+"0"+nl # 
					final_file.write(outline3) # 
#					print ("saving current variant already in file with 0 for new sample 3 -",outline3) #
					C1=current_sample.readline() #
#			elif (line1_array[0] < line2_array[0]): # different chromosome
			elif (new_chr < existing_chr): # different chromosome
#				print ("new < exist")
#				new_sample_len=len(line1_array) #
				variant_line[:new_sample_len-2]=line1_array[:new_sample_len-2] # one less to exclude the genotype and \n 
#				variant_line[:new_sample_len-1]=line1_array[:new_sample_len-1] # one less to exclude the \n 
				outline3='\t'.join(variant_line[:-1])+tab+line1_array[new_sample_len-2]+nl # add the sample genotype
#				outline3='\t'.join(variant_line[:-1])+tab+"1"+nl # replace the last 0 with a 1
				final_file.write(outline3) #
#				print ("saving new variant and new sample 4 -",outline3) #
				S1=new_sample.readline() #
			else: # (line1_array[0] > line2_array[0]): #
#				print ("new > exist")
				outline3='\t'.join(line2_array[:-1])+tab+"0"+nl # 
				final_file.write(outline3) # 
#				print ("saving current variant already in file with 0 for new sample 5 -",outline3) #
				C1=current_sample.readline() #
#
			if not S1 : # no more lines in new sample
				new_sample_not_finished=False #
#				print ("new sample file finished") #
#
			if not C1 : # no more lines in old combined
				current_sample_not_finished=False #
#				print ("current sample file finished") #
		#
		while (new_sample_not_finished): # there are still line to work
			line1_array=re.split('\t|\n|\r',S1) # split the variant line
#			print "sample file left", line1_array #
			new_sample_len=len(line1_array) #
			variant_line[:new_sample_len-2]=line1_array[:new_sample_len-2] # one less to exclude the genotype and \n 
#			variant_line[:new_sample_len-1]=line1_array[:new_sample_len-1] # one less to exclude the \n 
			outline3='\t'.join(variant_line[:-1])+tab+line1_array[new_sample_len-2]+nl # add the sample genotype
#			outline3='\t'.join(variant_line[:-1])+tab+"1"+nl # replace the last 0 with a 1
			final_file.write(outline3) #
#			print ("saving new variant and new sample 6 -",outline3) #
			S1=new_sample.readline() #
			if not S1 : # no more lines in new sample
				new_sample_not_finished=False #

#
		while (current_sample_not_finished): # there are still line to work
			line2_array=re.split('\t|\n|\r',C1) # split the variant line
#			print "curr file left", line2_array #
			outline3='\t'.join(line2_array[:-1])+tab+"0"+nl # 
			final_file.write(outline3) # 
#			print ("saving current variant already in file with 0 for new sample 6 -",outline3) #
			C1=current_sample.readline() #
			if not C1 : # no more lines in new sample
				current_sample_not_finished=False #
				
# all done - now sort the file 
#
#	COMMAND="grep \"CHROM\" "+full_file2+" > "+sort_file1 # save the header line  
#	subprocess.call(COMMAND, shell=True) #
#
#	COMMAND="sort -V -k1,5 "+full_file2+" > "+sort_file1 #  
#	COMMAND="sort -k1,2V -k3,5 "+full_file2+" > "+sort_file1 #  
	COMMAND="sort -k1,1V -k2,2g -k4,5 "+full_file2+" > "+sort_file1 #  
	subprocess.call(COMMAND, shell=True) #
#
	copyfile(sort_file1,full_file2) # copy back
#
##################################################################
