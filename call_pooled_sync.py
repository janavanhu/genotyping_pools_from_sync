import numpy
import numpy as np
import random
from scipy import stats
import sys
import os


def read_syncAll(name,NAME,minT,maxTall,maxTpop,nPops,minPops,minAF,minCount):
    ws = open(name, "r")
    WS = open(NAME, "w")
    genome={}
    count =0 #count variable positions
    count6 =0 #count all positions fulfilling coverage requirements
    c=0
    for line in ws:
        c+=1
        l =line.split('\t')
        pos = str(l[0])+ '__'+ str(l[1])
        a=0
        b1=0
        A=[]
        B=[]
        C=[]
        D=[]
        E=[]
        G_All=0
        for p in range(3,nPops+3):                       #get reads for each pop
                    G=list(map(int,l[p].split(':')))    #convert into list
                    G_All+= sum(G)                      #count number of reads across pops
                    a+=1
                    if sum(G) >=minT and sum(G) <= maxTpop: #          #only count population reads within coverage thresholds
                            A.append(int(G[0]))             #read numbers for each nucleotide
                            B.append(int(G[1]))
                            C.append(int(G[2]))
                            D.append(int(G[3]))
                            E.append(sum(G))                #number of reads per pop                    
                            b1+=1
                    else:
                            A.append(int(0))             #read numbers for each nucleotide
                            B.append(int(0))
                            C.append(int(0))
                            D.append(int(0))
                            E.append(int(0))                #number of reads per pop
        L = [A,B,C,D]
        b=0
        X=0
        Y=0
        if G_All>=minT and G_All<= maxTall and b1>=minPops:#
            WS.write(line) 
            count6+=1
    print('number of sites in input file: ',c)
    print('number of sites fulfilling coverae requirements (variable and non-variable): ',count6)
    WS.close()
    ws.close()
    return()




def read_syncVariant(name,NAME,minT,maxTall,maxTpop,nPops,minPops,minAF,minCount):
    ws = open(name, "r")
    WS = open(NAME, "w")
    genome={}
    count =0 #count variable positions
    for line in ws:
        l =line.split('\t')
        pos = str(l[0])+ '__'+ str(l[1])
        a=0
        b1=0
        A=[]
        B=[]
        C=[]
        D=[]
        E=[]
        G_All=0
        for p in range(3,nPops+3):                       #get reads for each pop
                    G=list(map(int,l[p].split(':')))    #convert into list
                    G_All+= sum(G)                      #count number of reads across pops
                    a+=1
                    if sum(G) >=minT and sum(G) <= maxTpop: #          #only count population reads within coverage thresholds
                            A.append(int(G[0]))             #read numbers for each nucleotide
                            B.append(int(G[1]))
                            C.append(int(G[2]))
                            D.append(int(G[3]))
                            E.append(sum(G))                #number of reads per pop                    
                            b1+=1
                    else:
                            A.append(int(0))             #read numbers for each nucleotide
                            B.append(int(0))
                            C.append(int(0))
                            D.append(int(0))
                            E.append(int(0))                #number of reads per pop
        L = [A,B,C,D]
        b=0
        X=0
        Y=0
        if G_All>=minT and G_All<= maxTall and b1>=minPops:# 
            count2=0
            nucs =['A','T','C','G']
            minF=[]
            minC= []            
            for p in L:
                count3=0                #track population for read counts
                if sum(p)>=minCount:           #sum of all reads for nucleotide across pops
                    X+=1
                    minC.append(nucs[count2])
                    for x in p:
                            if E[count3]>0:
                                    if x/E[count3]>=minAF:       #Y counts if nucleotide has freq >= minAF in at least 1 pop
                                        Y+=1
                                        minF.append(nucs[count2])               #keep track of nucleotides with minimum freq at pos
                                        break
                            count3+=1
                count2+=1
            count2=0
        if X>=2 and G_All >=minT and Y>=2 and G_All<= maxTall :       #X>=2 to have at least two different nucleotides with minimum of minCount reads >SNP, G_All read coverage across pops, Y indicates number of nucleotides observed with freq >=minAF in at least 1 pop
                    WS.write(line)      #
                    count+=1
                    if pos not in genome.keys():
                        genome[pos] = {}
                        for nuc in nucs:
                                if nuc in minF and nuc in minC:
                                        genome[pos][nuc]=[]             #only have SNPs with minimum coverage/count in gemome dic
                    genome[pos]['ref']=str(l[2])
                    b =3
                    d=0  #identifier of pops
                    for x in Single_pops:
                                geno = l[b].split(':')
                                C=0
                                c=0
                                for nuc in nucs:
                                        if nuc in genome[pos].keys() and E[d]>0:
                                                if int(geno[c])>=1 and float(geno[c])/E[d]>=minAF: # and sum(list(map(int,geno)))>=50 and sum(list(map(int,geno)))<=500
                                                        genome[pos][nuc].append(Single_pops_names[d]) # float(geno[c] for AFcount             #only append counts for SNPs with minimum coverage in gemome dic
                                                        C+=1
                                        c+=1
                                
                                if sum(list(map(int,geno)))>=minT and C>=1 and sum(list(map(int,geno)))<=maxTpop:#
                                        x[pos]=[]                                               #use single pop at frequencey if coverage >=minT and count for a nuc considered in genome_dic                                
                                d+=1
                                b+=1
                        

    print('number of sites fulfilling coverage requirements (variable only): ',count)
    WS.close()
    ws.close()
    return(genome)



def get_AF_genepop(genome,minAF):
    pseudo_indivs=1/minAF
    print(pseudo_indivs)
    ws = open(NAME, "r")
    C=0
    g_list=[]
    p_count=0
    for line in ws:
                    p_count+=1
                    l =line.split('\t')
                    if len(l)>1:
                            pos = str(l[0])+ '__'+ str(l[1])
                            b =3
                            pos_N =['A','T','C','G']
                            if pos in genome.keys() :
                                g_list.append(pos)
                                count2 =0
                                for x in Single_pops:
                                    if pos in x.keys():
                                        genos = ['01','02','03','04']
                                        total = 0
                                        a = 0
                                        geno=l[b].split(':')
                                        for t in geno[0:4]:
                                            if pos_N[a] in genome[pos].keys() and int(t)>=1:        #only take read numbers of SNPs into account that fulfilled above criteria and are in genome dic
                                                if Single_pops_names[count2] in genome[pos][pos_N[a]]:
                                                        total+=int(t)
                                            a+=1
                                        for count in range(0,4):
                                            if pos_N[count] in genome[pos].keys() and int(geno[count])>= 1:     #only use reads of nucleotides considered (genome dic) in AF calculation
                                                if Single_pops_names[count2] in genome[pos][pos_N[count]]:
                                                        AF=int(geno[count])/float(total)                    
                                                        AF2 = round(AF,3)
                                                        g = [genos[count]] * int(round(AF2*pseudo_indivs,0))                 #append nuc code as many times as AF*pseudo_indivs (total number of indivs in genepop file)
                                                        random.shuffle(g)
                                                        for element in g:
                                                                x[pos].append(element)
                                        if len(x[pos])==0:
                                                del x[pos]
                                                C+=1
                                        elif len(x[pos])!=pseudo_indivs and len(x[pos])!=0:                               #if number uneven remove or add most common
                                               COUNT= x[pos].count(x[pos][0])
                                               IT = x[pos][0]
                                               for item in set(x[pos]):
                                                   if x[pos].count(item)>COUNT:
                                                       COUNT = x[pos].count(item)
                                                       IT =item
                                               if len(x[pos])>pseudo_indivs:
                                                    dif = int(len(x[pos])-pseudo_indivs)
                                                    for more in range(0,dif):
                                                        x[pos].remove(IT)
                                               else:
                                                    dif = int(pseudo_indivs - len(x[pos]))
                                                    for more in range(0,dif):
                                                        x[pos].append(IT)
                                               
                                    count2+=1
                                    b+=1
    print(len(g_list))
    return(g_list)





def write_genepop(g_list,nPops,name,minAF):
    pseudo_indivs=1/minAF
    ws = open(name, "w")
    A=str(sys.argv[:])+'\n'
    for i in g_list:
             A= A+str(i)+ ',\t'#
    ws.write(A+'\n')#
    a=0
    b=0
    for i in Single_pops:#
                ws.write('Pop\n')
                for count in range(0,int(pseudo_indivs)):
                    b+=1
                    geno = str(Single_pops_names[a])+',\t'#
                    for e in g_list:
                        if e in i.keys():
                            geno = geno + i[e][count]+'\t'
                        else:
                            geno = geno + 'NA\t'
                    ws.write(geno+'\n')
                a+=1

    





##########################
name = str(sys.argv[1])
minT= int(sys.argv[2])
maxTall= int(sys.argv[3])
maxTpop= int(sys.argv[4])
nPops= int(sys.argv[5])
minPops= int(sys.argv[6])
minAF= float(sys.argv[7])
minCount= int(sys.argv[8])

##
Single_pops = [{} for _ in range(nPops)]#
Single_pops_names = []
for i in range(nPops):
    Single_pops_names.append(str(i))


##new syncs to write
NAME= str(name[:-5])+'_All_sites.sync'
NAME2= str(name[:-5])+'_All_variable_sites.sync'
NAME3= str(name[:-5])+'.genepop'

read_syncAll(name,NAME,minT,maxTall,maxTpop,nPops,minPops,minAF,minCount)
genome = read_syncVariant(NAME,NAME2,minT,maxTall,maxTpop,nPops,minPops,minAF,minCount)
g_list=get_AF_genepop(genome,minAF)
write_genepop(g_list,nPops,NAME3,minAF)

#

