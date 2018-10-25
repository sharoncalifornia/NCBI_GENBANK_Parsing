from Bio import SeqIO
from Bio import Entrez
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
import os
import re

TRUE = 1
FALSE = 0
begin_count =FALSE
hot_spring_hit=0

Entrez.email="sharonzh10@yahoo.com"
file_in_fmt6 = "output_43_seqs_fmt6.out"
read_geo_filename = "43_geo_seqs_ID.txt"
write_filename = "43_seqs_geo_3hits.txt"
write_file_handle = open(write_filename, 'w')

gbk_filename="43_seqs_geo_3hits.gbk"
if os.path.isfile(gbk_filename):
        os.remove(gbk_filename)

temp_gbk_filename="temp.gbk"

def FetchGenBankRecord(Pubmed_ID):
    flag = 0
    
    if os.path.isfile(temp_gbk_filename):
        os.remove(temp_gbk_filename)
        
    net_handle = Entrez.efetch(db="nucleotide", id=Pubmed_ID, rettype="gb", retmode="text")
    
    temp_gbk_file_handle = open(temp_gbk_filename, 'w')
    temp_gbk_file_handle.write(net_handle.read())
    temp_gbk_file_handle.close()
    flag =SearchHotSpring(flag)
    
    net_handle = Entrez.efetch(db="nucleotide", id=Pubmed_ID, rettype="gb", retmode="text")
    gbk_file_handle=open(gbk_filename, "a")
    gbk_file_handle.write(net_handle.read())
    gbk_file_handle.close()
    net_handle.close()
    print ("inside FetchGenBankRecord")

    return flag


def SearchHotSpring(hot_spring_flag):
    
        temp_gbk_file_handle = open(temp_gbk_filename, 'r')
 
        for line in temp_gbk_file_handle:
                if line.find("hot spring")>1 or line.find("warm spring")>1 or line.find("geothermal spring")>1 or line.find("sulfidic spring")>1 or line.find("sulfurous spring")>1 or line.find("sulfide spring")>1 or line.find("sulphur spring")>1 or line.find("sulfur spring")>1 or line.find("saline spring")>1:
                    hot_spring_flag = 1
                    gbk_file_handle=open(gbk_filename,"a")
                    gbk_file_handle.write("string matched:\t"+line+'\n')
                    gbk_file_handle.close()

        temp_gbk_file_handle.close()

        return hot_spring_flag



def inputGenbankInfo():
    
    temp_gbk_file_handle = open(temp_gbk_filename, 'r')
    
    for seq_record in SeqIO.parse(temp_gbk_file_handle, "gb"):
        
        debug_file_handle.write(seq_record.annotations["organism"]+'\t')
        debug_file_handle.write(seq_record.annotations["source"]+'\n')

        
    temp_gbk_file_handle.close()
        
    return

file_geo_handle = open(read_geo_filename, 'r')

for line_geo in file_geo_handle:
    
    if line_geo.find("gi|")>0:
        x= re.split("\t",line_geo)
        print "inside geo_file           "+ x[0]
        gbk_file_handle=open(gbk_filename, "a")
        gbk_file_handle.write("sequence_id----------------->"+x[0]+'\n')
        gbk_file_handle.close()
        file_fmt6_handle = open(file_in_fmt6,'r')
        for line_fmt6 in file_fmt6_handle:
           
                y = re.split("\t",line_fmt6)

                if (y[0]==(x[0])) & (len(y[0])==len(x[0])):
    
                    print (y[1])
                    genbank_id = re.split('[|]+',y[1])
                    print (genbank_id[1])
                    gbk_file_handle=open(gbk_filename, "a")
                    gbk_file_handle.write(line_fmt6)
                    gbk_file_handle.close()
                    hot_spring_hit = FetchGenBankRecord(genbank_id[1])
                    
        hot_spring_hit=0
        file_fmt6_handle.close()
file_geo_handle.close()





