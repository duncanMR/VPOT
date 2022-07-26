###########################################################################################################
# #
###########################################################################################################
#./VPOT.py genepanelf 
import sys, re, glob, os, subprocess, time #
import numpy as np #
import VPOT_conf #
import pandas as pd
from shutil import copyfile #
#
#

###########################################################################################################
# Define global variables
##########################################################################################################
suffix="" #no suffix needed
#
supplied_args=0 #
#
tab='\t' # 
nl='\n' #
#
GENE_loc=-1 #
#
info_msg1_1="VPOT genef : Invalid number of inputs, must have four :" 
info_msg1_2="VPOT genef : 1) output destination directory + prefix" #
info_msg1_3="VPOT genef : 2) Input file - output from VPOT prioritisation process" 
info_msg1_4="VPOT genef : 3) Gene list in csv format" #
info_msg1_5="VPOT genef : 4) Cancer type (must be column name in genelist)" #
#
############################################################################################################
#
###########################################################################################################
#
###########################################################################################################

def initial_setup():
#
	global supplied_args 
	#
	print ("initial_setup():") #
	print (suffix) #
#	print sys.argv #
	supplied_args=len(sys.argv) #
	print (supplied_args) #
#
	if (supplied_args != 6 ):  # arg [0] is the python program
		print (info_msg1_1+nl+info_msg1_2+nl+info_msg1_3+nl+info_msg1_4+nl+info_msg1_5) #
		return 1 #
	else :
		VPOT_conf.output_dir=sys.argv[2] #
		VPOT_conf.input_file=sys.argv[3] #
		VPOT_conf.gene_list=sys.argv[4] #
		VPOT_conf.cancer_type=sys.argv[5] #

		#extract list of genes by panel from input csv
		df = pd.read_csv(VPOT_conf.gene_list,sep=",")
		df.fillna("MiscCPG", inplace=True) #If panel is unspecified, lump into misc category
		grp_obj = df.groupby(VPOT_conf.cancer_type) #returns groupby object
		VPOT_conf.panel_name_list = grp_obj.groups.keys() #e.g. BreastCancerHighRisk
		#make list of output files, one for each panel
		VPOT_conf.final_output_file_list = [VPOT_conf.output_dir + (panel) + "_output.txt" for panel in VPOT_conf.panel_name_list]
		#make one list of genes for each panel
		VPOT_conf.panel_gene_list=[list(grp_obj.get_group(panel)["Gene"]) for panel in VPOT_conf.panel_name_list]
		print ("output : ",VPOT_conf.final_output_file_list) #
	
	return 0 #
#
#
###########################################################################################################
#
###########################################################################################################
def filter_the_variants(): #
#
# input file for filtering is the output from the VPOT process. This means the variant gene name is in a column named GENE_NAME.
	global GENE_loc #
#
#	print ("filter_the_variants(): ") #
	#
#	print (VPOT_conf.input_file,VPOT_conf.final_output_file) # 
	for n in range(len(VPOT_conf.panel_name_list)):
		print(n)
		with open(VPOT_conf.input_file,'r',encoding="utf-8") as variants_file, open(VPOT_conf.final_output_file_list[n],'w',encoding="utf-8") as filtered_file : # 
			for line1 in variants_file: # work each line of new sample vcf file 
				write_it=False # initialise score 
				line_parts=re.split('\t|\n|\r',line1) # split the variant up
				#			print "line part 0 : ",line_parts[0] #
				if ("#CHROM" != line_parts[2]): #
					#				print src_line1 #
					#				
					write_it=filter_variants_by_GN(line_parts, n) # check get priority score
					#
				else : # save the header line	
					write_it=True # initialise score 
					for i, content in enumerate(line_parts): # return the value and index number of the GENE_NAME item in the line array 
						#					print "content-",content,"/",i				#
						if (content == VPOT_conf.GENE_NAME) : # 
							GENE_loc=i #save sample location
							#						print "INFO_loc: ",INFO_loc #
							#	
				if (write_it): #
					filtered_file.write(line1) # write the line to final output file 
###########################################################################################################
#
###########################################################################################################
def filter_variants_by_GN(INFO_details, n): #
#
#	print "filter_variants_by_GN(INFO_details): " #
#
	val=False #
	for gene_id in VPOT_conf.panel_gene_list[n]: # work each line of new sample vcf file 
		gene_id1=gene_id.rstrip() #
#		print ("Gene ID",gene_id1,"in",INFO1[i+1]  # move to pred_array slot
		if ( gene_id1 in INFO_details[GENE_loc] ) :		# is this the gene we are looking for   
			gene_detail=re.split(',',INFO_details[GENE_loc]) # split into annotations
			#print INFO_details[GENE_loc],"/",gene_detail,"/",gene_id1 #
			for k in range(len(gene_detail)): #
#				print gene_detail[k],"/",gene_detail,"/",gene_id1 #
				if ( gene_id1 == gene_detail[k] ) :		# is this the gene we are looking for   
#					print "match: ",gene_detail[k],"/",gene_detail,"/",gene_id1 #
					val=True #
					break #
#		print "genes go around-",gene_id1 #
		if ( val ) :		# we have found the gene   
			break # then get out and move to next variant
						
#
	return val #
#	
##
###########################################################################################################
#
###########################################################################################################
def main(): #
##
	VPOT_conf.init() #
	#
	print ("Gene Filter - Main") #
	#
	if (initial_setup() != 0): #
#		print "no good" #
		return #
	#
# Now filter the input file by gene list 
	filter_the_variants() #
	#
