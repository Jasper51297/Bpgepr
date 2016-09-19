"""
Author: Jasper van Dalum
Version: 1.1
"""

import subprocess
import time
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
    cmd1 = 'wget ftp://ftp.ensembl.org/pub/release-85/fasta/sarcophilus_harrisii/cds/Sarcophilus_harrisii.DEVIL7.0.cds.all.fa.gz'
    cmd2 = 'gunzip Sarcophilus_harrisii.DEVIL7.0.cds.all.fa'
    cmd3 = 'mv Sarcophilus_harrisii.DEVIL7.0.cds.all.fa TasDev.fa'
    cmd = [cmd1, cmd2, cmd3]
    for i in cmd:
        time.sleep(2)
        subprocess.Popen(i.split(' '), stdout=subprocess.PIPE )
    print('Files zijn gedownload.')

def BLAST():
    cmd1 = 'formatdb -i TasDev.fa -p F -o'
    cmd2 = 'blastall -p blastn -m 8 -d TasDev.fa -i sequenties_groep_04.fa -o out_blast.txt'
    cmd = [cmd1, cmd2]
    for i in cmd:
        subprocess.Popen(i.split(' '), stdout=subprocess.PIPE )
        time.sleep(4)
    print('Done. De output staat in out_blast.txt.')

def eValue():
    subprocess.Popen("""awk '{if ($11=="0.0") print $2}' out_blast.txt >> E0.txt""", shell = True)
    print('De files met een lage E-Value staan in E0.txt')
    time.sleep(1)

def getProtein():
    subprocess.Popen("""cat TasDev.fa | egrep -f E0.txt | sed 's/ /\n/g' | egrep 'gene:' | sed 's/gene://g' > protein.txt""", shell = True)
    print('De namen van de eiwitten met een lage E-Value staan in protein.txt')
    time.sleep(1)

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

