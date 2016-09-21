"""
Author: Jasper van Dalum
Version: 1.1
"""

import os
import sys

def menu():
    print('\n\n\n==========OPTIONS==========')
    print('1. Download Files')
    print('2. BLAST')
    print('3. E-value sort')
    print('4. Get Proteïns')
    print('5. Print Results')
    print('6. Alles')
    print('7. Exit')
    m_input = input('Input: ')
    return str(m_input)
    
def download():
    cmd1 = 'wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000189315.1_Devil_ref_v7.0/GCF_000189315.1_Devil_ref_v7.0_protein.faa.gz'
    cmd2 = 'gunzip GCF_000189315.1_Devil_ref_v7.0_protein.faa'
    cmd3 = 'mv GCF_000189315.1_Devil_ref_v7.0_protein.faa TasDev.fa'
    #Niet mogelijk om van deze site te downloaden, needs fix.
    #cmd4 = 'wget http://www.ncbi.nlm.nih.gov/genomes/Genome2BE/genome2srv.cgi?action=GetFeatures4Grid&amp;download=on&amp;type=Proteins&amp;genome_id=3066&genome_assembly_id=34269&mode=2&is_locus=1&is_locus_tag=0&optcols=1,0,,0,0&amp;filterText=replicon.gi%20=%20-1'
    cmd = [cmd1, cmd2, cmd3]
    for i in cmd:
        os.system(i)
    print('Files zijn gedownload.')

def BLAST():
    cmd1 = 'formatdb -i TasDev.fa -p T -o'
    cmd2 = 'blastall -p blastx -m 8 -d TasDev.fa -i sequenties_groep_04.fa -o out_blast.txt'
    cmd = [cmd1, cmd2]
    for i in cmd:
        os.system(i)
    print('Done. De output staat in out_blast.txt.')

def TopHits():
    os.system('sort -k1,1 -k12,12nr -k11,11n  out_blast.txt | sort -u -k1,1 --merge > E0.txt')
    os.system("cat E0.txt | tr '\t' '!'| awk -F '!' '{print $2}' > tophits.txt")
    print('De files met een lage E-Value staan in E0.txt')

def getProtein():
    #moet greppen protein code uit protein tabel (ProteinTable.txt)
    os.system('cat ProteinTable.txt | egrep -f tophits.txt | awk "{print $1}" > protein.txt')
    os.system("""cat protein.txt | tr '\t' '!' | awk -F '!' '{print $2 "\t" $6 "\t" $7 "\t" $8 "\t" $10}' > results.txt""")
    print('De namen van de eiwitten met een lage E-Value staan in protein.txt')

def printResults():
    inFile = open('results.txt')
    inhoud = inFile.read()
    print('Replicon Acc \t GeneID \t Locus \t proteïn product \t Protein name ')
    print(inhoud)

def main():
    while 1 == 1:
        m_input = menu()
        if m_input == '1':
            download()
        if m_input == '2':
            BLAST()
        if m_input == '3':
            TopHits()
        if m_input == '4':
            getProtein()
        if m_input == '5':
            printResults()
        if m_input == '6':
            download()
            BLAST()
            TopHits()
            getProtein()
        if m_input == '7':
            sys.exit()
            
main()

