#!/usr/bin/python
"""
Author: Jasper van Dalum, Paul de Raadt, Brenda van den Berg, Duncan, Joery 
Version: 2.0
"""
#GeneID (KEY) | Gennaam | Genstart | Genstop | AccessionCode | EiwitNaam | PathwayName | PathwayDesc | EC | geneSeq | protSeq | exonStart | exonStop

import os
import psycopg2
import sys
from subprocess import call

def menu():
    print '\n\n\n==========OPTIONS=========='
    print '1. Download Files'
    print '2. BLAST'
    print '3. Gather results'
    print '4. Create Database'
    print '5. Exit'
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
    print 'Download: 1/1'

def BLAST():
    cmd1 = 'formatdb -i TasDev.fa -p T -o'
    cmd2 = 'blastall -p blastx -m 8 -d TasDev.fa -i sequenties_groep_04.fa -o out_blast.txt'
    cmd = [cmd1, cmd2]
    for i in cmd:
        os.system(i)
    print 'BLAST: 1/1'

def getResults():
    os.system("sort -k1,1 -k12,12nr -k11,11n  out_blast.txt | sort -u -k1,1 --merge | tr '\t' '!'| awk -F '!' '{print $2}' > tophits.txt")
    results = getProtein()
    for geneid in results.keys():
        lst, link = results[geneid], ('http://www.kegg.jp/dbget-bin/www_bget?shr:' + str(geneid))
        os.system("lynx %s -dump -nolist > kegg"%link)
        pathwayname, pathwaydesc = getPathways()
        EC = getEC()
        geneSeq, protSeq = getSeq(lst[3])
        exonStart, exonStop = getExon(geneid)
        varlist = [pathwayname, pathwaydesc, EC, geneSeq, protSeq, exonStart, exonStop]
        for var in varlist:
            lst.append(var)
        results[geneid] = lst
    return results
        
def getProtein():
    results = {}
    os.system('cat ProteinTable.txt | egrep -f tophits.txt | awk "{print $1}" > protein.txt')
    proteins = os.popen("""cat protein.txt | tr '\t' '!' | awk -F '!' '{print $6 "\t" $7 "\t" $3 "\t" $4 "\t" $8 "\t" $10}'""").readlines()
    for line in proteins:
        line = line.replace('\r', '')
        line = line.replace('\n', '')
        line = line.split('\t')
        key = line[0]
        del line[0]
        results[key] = line
    return results

def getPathways():
    pathway = os.popen("cat kegg | tr -d '\n' | tr ')' '\n' | sed 's/Module/\\n/g' | sed 's/Brite/\\n/g' |egrep Pathway | sed 's/Pathway//g' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'").read()
    pathways = pathway.split('   ')
    a, pathwayname, pathwaydesc = 0, '', ''
    for pathway in pathways:
        if a%2 == 0:
            pathwayname += pathways[a]
            pathwayname += '\t'
        else:
            pathwaydesc += pathways[a]
            pathwaydesc += '\t'
        a += 1
    pathwaydesc, pathwayname = pathwaydesc.strip(), pathwayname.strip()
    return pathwayname, pathwaydesc

def getEC():
    EC = os.popen("cat kegg|egrep EC: | awk -F '[\[\]]' '{print $2}'").read()
    EC = EC.replace('EC:', '')
    EC = EC.replace('\n', '')
    EC = EC.replace(' ', '\t')
    return EC

def getSeq(AccCode):
    geneSeq = os.popen("cat kegg | tr -d '\n' | sed 's/DB search/\\n!/g' | sed 's/NT seq/\\n/g' | egrep ! | tr -d ! | tr -d ' ' | tr -d '\n'").read()
    protSeq = os.popen("cat TasDev.fa | tr -d '\n' | tr '>' '\n' | egrep %s | sed 's/]/\\n!/g' | egrep ! | tr -d ! | tr -d '\n'"%AccCode).read()
    return geneSeq, protSeq

def getExon(geneid):
    exons =  os.popen("cat ref_Devil_ref_v7.0_scaffolds.gff3 | egrep %s | egrep exon | awk '{print $4, $5}'"%geneid).readlines()
    exonStart, exonStop = '', ''
    for pos in exons:
        pos = pos.replace('\n', '\t')
        poslist = pos.split(' ')
        exonStart += poslist[0]
        exonStart += '\t'
        exonStop += poslist[1]
    exonStop = exonStop.strip()
    exonStart = exonStart.strip()
    return exonStart, exonStop

def login():
    conn_string = "host='127.0.0.1' dbname='bpgepr' user='user' password='password'"
    print "Connecting to database\n	->%s" % (conn_string)
    conn = psycopg2.connect(conn_string)
    print "Connected!\n"
    conn.autocommit = True
    cursor = conn.cursor()
    return cursor, conn

def deleteTables(cursor):
    tablelist = ['EIWIT_PATHWAY_04', 'PATHWAY_04', 'EXON_04', 'EIWIT_04', 'GEN_04']
    for table in tablelist:
        cursor.execute("DROP TABLE IF EXISTS %s"%table)

def createTables(cursor):
    GEN_04 = "CREATE TABLE GEN_04 (GenID VARCHAR(25) NOT NULL UNIQUE PRIMARY KEY, Naam VARCHAR(25), Chromosoom SMALLINT, Start_gen INT, Stop_gen INT, Sequentie TEXT);"
    EXON_04 = "CREATE TABLE EXON_04 (GenID VARCHAR(25) NOT NULL REFERENCES GEN_04(GenID), Exon_start INT NOT NULL, Exon_stop INT NOT NULL, PRIMARY KEY(GenID, Exon_start));"
    EIWIT_04 = "CREATE TABLE EIWIT_04 (Accession_code VARCHAR(25) UNIQUE PRIMARY KEY, GenID VARCHAR(25) REFERENCES GEN_04(GenID), Naam TEXT, EC_nummer VARCHAR(25), Sequentie TEXT);"
    PATHWAY_04 = "CREATE TABLE PATHWAY_04 (PathwayID VARCHAR(25) PRIMARY KEY, Omschrijving TEXT);"
    EIWIT_PATHWAY_04 ="CREATE TABLE EIWIT_PATHWAY_04 (Accession_code VARCHAR(25) REFERENCES EIWIT_04(Accession_code), PathwayID VARCHAR(25) REFERENCES PATHWAY_04(PathwayID), PRIMARY KEY(Accession_code, PathwayID));"
    commands = [GEN_04, EXON_04, EIWIT_04, PATHWAY_04, EIWIT_PATHWAY_04]
    for sql in commands:
        cursor.execute(sql)
        
def main():
    os.system("""command -v lynx >/dev/null 2>&1 || { echo >&2 "I require Lynx but it's not installed.  ."; sudo apt install lynx; }""")
    while 1 == 1:
        m_input = menu()
        if m_input == '1':
            download()
        if m_input == '2':
            BLAST()
        if m_input == '3':
            results = getResults()
            print results
            outFile = open('Results.txt', 'w')
            results = str(results)
            outFile.write(results)
            outFile.close()
        if m_input == '4':
            cursor, conn = login()
            deleteTables(cursor)
            createTables(cursor)
            
            conn.close()
        if m_input == '5':
            sys.exit()
            
main()
