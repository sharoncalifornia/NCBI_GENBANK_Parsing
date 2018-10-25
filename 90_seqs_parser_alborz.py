from Bio import SeqIO
from Bio import Entrez
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
import os
import re

TRUE = 1
FALSE = 0
begin_count =FALSE
hot_spring_flag=FALSE
new_record_flag=FALSE

Entrez.email="sharonzh10@yahoo.com"
temp_gbk_filename = "90_seqs_temp_gb.gbk"
temp_gbk_file_handle = open(temp_gbk_filename, 'w')

geo_gbk_filename="90_seqs_geo.gbk"
if os.path.isfile(geo_gbk_filename):
        os.remove(geo_gbk_filename)


debug_filename="90_geo_seqs_ID.txt"
if os.path.isfile(debug_filename):
        os.remove(debug_filename)


file_in_fmt7 = "output_90_seqs_7.out"

#file_write_out="geo_03_05.txt"
#if os.path.isfile(file_write_out):
#   os.remove(file_write_out)

#fout=open(file_write_out,"a")
#fout.write("organism\t"+"source\t"+"query_id\t"+"denovo_ID\t"+"Hit_ID\t"+"Hit_Def\t"+"similarity\t"+"sequence length\t"+"mismatch\t"+"gap open\t"+"e-value\t"+"bit score\t"+"sequence_annotation\n")
#fout.close()
geo_thermal_temp_line=""
temp_line=""

def FetchGenBankRecord(Pubmed_ID):
    if os.path.isfile(temp_gbk_filename):
        os.remove(temp_gbk_filename)
        
        net_handle = Entrez.efetch(db="nucleotide", id=Pubmed_ID, rettype="gb", retmode="text")
        temp_gbk_file_handle=open(temp_gbk_filename, "a")
        temp_gbk_file_handle.write(net_handle.read())
        temp_gbk_file_handle.close()
        net_handle.close()

    return 1

def SearchHotSpring(hot_spring_flag):
    temp_gbk_file_handle = open(temp_gbk_filename, 'r')
    hot_spring_flag =0
 
    for line in temp_gbk_file_handle:
        if line.find("hot spring")>1 or line.find("warm spring")>1 or line.find("geothermal spring")>1 or line.find("sulfidic spring")>1 or line.find("sulfurous spring")>1 or line.find("sulfide spring")>1 or line.find("sulphur spring")>1 or line.find("sulfur spring")>1 or line.find("saline spring")>1:
            hot_spring_flag = 1
            geo_thermal_temp_line=line
            genbank_file_handle=open(geo_gbk_filename,"a")
            genbank_file_handle.write("string matched:\t"+geo_thermal_temp_line+'\n')
            genbank_file_handle.close()    
    return hot_spring_flag

def inputGenbankInfo():
    
    temp_gbk_file_handle = open(temp_gbk_filename, 'r')
    
    for seq_record in SeqIO.parse(temp_gbk_file_handle, "gb"):
        
        debug_file_handle.write(seq_record.annotations["organism"]+'\t')
        debug_file_handle.write(seq_record.annotations["source"]+'\n')

        
    temp_gbk_file_handle.close()
        
    return

def WriteFileHotSpring():
    fout=open(file_write_out,"a")
    fout.write(query_ID+'\t')
    fout.write(denovo_ID_long+'\t')
    fout.write(hit_ID+'\t')
    fout.write(hit_def+'\t')
    fout.close()
    return


file_fmt7_handle = open(file_in_fmt7,'r')

for line_fmt7 in file_fmt7_handle:
        if line_fmt7.find("2 hits found")>0:
                new_record_flag=TRUE
        if line_fmt7.find("gi")>0 :
                temp_line = line_fmt7
                x=re.split("\t",temp_line)
                if float(x[2])>90.0:
                        y=re.split('[|]+',x[1])
                        FetchGenBankRecord(y[1])
                        hot_spring_flag=SearchHotSpring(hot_spring_flag)
                        if (hot_spring_flag ==1 and  new_record_flag == TRUE):
                                debug_file_handle = open(debug_filename, 'a')
                                inputGenbankInfo()

                                debug_file_handle.write(line_fmt7)
                                debug_file_handle.close()
                                net_handle = Entrez.efetch(db="nucleotide", id=y[1], rettype="gb", retmode="text")
                                genbank_file_handle=open(geo_gbk_filename,"a")
                                genbank_file_handle.write("ID"+x[0]+"\t"+"gi_string"+"\t"+y[1]+'\n')
                                genbank_file_handle.write("string matched:\t"+geo_thermal_temp_line+'\n')
                                genbank_file_handle.write(temp_line)
                                genbank_file_handle.write(net_handle.read())
                                genbank_file_handle.close()
                                net_handle.close()
                                hot_spring_flag=FALSE
                                new_record_flag=FALSE
                                geo_thermal_temp_line=""
                temp_line=""

file_fmt7_handle.close()





