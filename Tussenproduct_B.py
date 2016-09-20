"""
Author: Jasper van Dalum
Version: 1.1
"""

import os
import sys

def menu():
    print('==========OPTIONS==========')
    print('1. Download Files')
    print('2. BLAST')
    print('3. E-value sort')
    print('4. Get Proteïns')
    print('5. Alles')
    print('6. Exit')
    m_input = input('Input: ')
    return str(m_input)
    
def download():
    cmd1 = 'wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000189315.1_Devil_ref_v7.0/GCF_000189315.1_Devil_ref_v7.0_rna.fna.gz'
    cmd2 = 'gunzip GCF_000189315.1_Devil_ref_v7.0_rna.fna'
    cmd3 = 'mv GCF_000189315.1_Devil_ref_v7.0_rna.fna TasDev.fa'
    cmd = [cmd1, cmd2, cmd3]
    for i in cmd:
        os.system(i)
    print('Files zijn gedownload.')

def BLAST():
    cmd1 = 'formatdb -i TasDev.fa -p F -o'
    cmd2 = 'blastall -p blastn -m 8 -d TasDev.fa -i sequenties_groep_04.fa -o out_blast.txt'
    cmd = [cmd1, cmd2]
    for i in cmd:
        os.system(i)
    print('Done. De output staat in out_blast.txt.')

def eValue():
    os.system("""awk '{if ($11=="0.0") print $2}' out_blast.txt >> E0.txt""")
    print('De files met een lage E-Value staan in E0.txt')

def getProtein():
    sys.os("""cat TasDev.fa | egrep -f E0.txt | sed 's/ /\n/g' | egrep 'gene:' | sed 's/gene://g' > protein.txt""")
    print('De namen van de eiwitten met een lage E-Value staan in protein.txt')

def main():
    while 1 == 1:
        m_input = menu()
        if m_input == '1':
            download()
        if m_input == '2':
            BLAST()
        if m_input == '3':
            eValue()
        if m_input == '4':
            getProtein()
        if m_input == '5':
            download()
            BLAST()
            eValue()
            getProtein()
        if m_input == '6':
            sys.exit()
            
main()

