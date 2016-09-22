"""
Author: Jasper van Dalum, Paul de Raadt, Brenda, Duncan, Joery 
Version: 2.0
"""

import os
import sys
from subprocess import call

def menu():
    print('\n\n\n==========OPTIONS==========')
    print('1. Download Files')
    print('2. BLAST')
    print('3. Gather results')
    print('4. Write results to file')
    print('5. All')
    print('6. Exit')
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
    print('Download: 1/1')

def BLAST():
    cmd1 = 'formatdb -i TasDev.fa -p T -o'
    cmd2 = 'blastall -p blastx -m 8 -d TasDev.fa -i sequenties_groep_04.fa -o out_blast.txt'
    cmd = [cmd1, cmd2]
    for i in cmd:
        os.system(i)
    print('BLAST: 1/1')

def getResults():
    os.system('sort -k1,1 -k12,12nr -k11,11n  out_blast.txt | sort -u -k1,1 --merge > E0.txt')
    os.system("cat E0.txt | tr '\t' '!'| awk -F '!' '{print $2}' > tophits.txt")
    print('Gather results: 1/3')
    getProtein()
    print('Gather results: 2/3')
    getPathways()
    print('Gather results: 3/3')

def getProtein():
    os.system('cat ProteinTable.txt | egrep -f tophits.txt | awk "{print $1}" > protein.txt')
    os.system("""cat protein.txt | tr '\t' '!' | awk -F '!' '{print $2 "\t" $6 "\t" $7 "\t" $8 "\t" $10}' > ProteinResults.txt""")

def getPathways():
    #os.system('sudo apt install lynx')
    os.system("cat protein.txt | awk '{print $6}' > linklist.txt")
    inFile = open('linklist.txt', 'r')
    inhoud = inFile.readlines()
    new = ''
    for regel in inhoud:
        newregel = ('http://www.kegg.jp/dbget-bin/www_bget?shr:' + regel)
        new += newregel
    inFile.close()
    inFile = open('linklist.txt', 'w')
    inFile.write(new)
    inFile.close()
    inFile = open('linklist.txt', 'r')
    inhoud = inFile.readlines()
    for i in inhoud:
        cmd = i
        cmd = cmd.replace('\n', '')
        os.system("lynx %s -dump -nolist | tr -d '\n' | tr ')' '\n' | sed 's/Module/\\n/g' | sed 's/Brite/\\n/g' |egrep Pathway >> PathwayResults.txt"%cmd)

def writeResults():
    inFile = open('ProteinResults.txt', 'r')
    ProtRes = inFile.readlines()
    inFile.close()
    inFile = open('PathwayResults.txt', 'r')
    PathRes = inFile.readlines()
    inFile.close()
    inFile = open('Results.txt', 'w')
    for i in range(len(ProtRes)-1):
        inFile.write(ProtRes[i] + PathRes[i])
    inFile.close()
    print('Writing results: 1/1')

def main():
    while 1 == 1:
        m_input = menu()
        if m_input == '1':
            download()
        if m_input == '2':
            BLAST()
        if m_input == '3':
            getResults()
        if m_input == '4':
            writeResults()
        if m_input == '5':
            download()
            BLAST()
            getResults()
            writeResults()
        if m_input == '6':
            sys.exit()
            
main()

