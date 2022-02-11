import os
import re

#Demultiplexes using a barcode file and stores new files in pradout directory
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool1_3.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq1/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool2_3.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq2/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool3_3.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq3/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool4_3.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq4/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool5_3.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq5/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool6_3.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq6/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool7_3.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq7/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool8_3.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq8/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool1_4.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq9/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool2_4.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq10/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool3_4.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq11/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool4_4.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq12/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool5_4.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq13/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool6_4.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq14/");
os.system("process_radtags -i fastq -o stacks/ -b ./sample_tables/pool7_4.csv -e ecoRI --renz_1 mspI -t 145 -r -c -q -p ./fastq/fastq15/");
#declare String#
sub = ""
#declare array "list"#
arr = []
i = 0
#os.getcwd()+#

for file in os.listdir('./stacks'):
    if file.endswith('.fq'):
        os.system("ustacks -t fastq -i %s -m 3 -M 4 -p 20 -f ./stacks/%s -o stacks" %(i,file) )
        #Remove the last 3 characters of the  file name and store them in a temp variable file2#
        file2 = file[0:-3]
        #concatenate #
        sub =  " " +  sub  + " -s ./stacks/"+ file2
        arr.append(file2)
	print('File name: %s'%file2)
       	i = i + 1

print("cstacks -P ./stacks -M popmap_2019_no_GB_HP -n 4 -p 20")
os.system("cstacks -P ./stacks -M popmap_2019_no_GB_HP -n 4 -p 20")

print("sstacks -P ./stacks -M popmap_2019_no_GB_HP -p 20")
os.system("sstacks -P ./stacks -M popmap_2019_no_GB_HP -p 20")

os.system("tsv2bam -P ./stacks -M popmap_2019_no_GB_HP -t 10")

os.system("gstacks -P ./stacks -M popmap_2019_no_GB_HP -t 10")


#Executes populations program exporting statistics in different formats and move them to output directory#
#-r 0.80 for parameter optimization based on Paris et al., 2017
os.system("populations -P ./stacks --popmap popmap_2019_no_GB_HP -t 10 -p 17 -r 0.80") 
