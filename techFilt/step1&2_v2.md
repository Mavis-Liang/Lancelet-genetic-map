```python
import shutil,sys,os,time
import numpy as np
```

### 1. 导入文件，并形成四个list


```python
with open("Extracted.vcf","r") as h:  #h = gzip.GzipFile("52_bbe_copy.vcf.gz") 
    i=0
    line_list=[]
    tmp_list=[]
    chrom=[]
    pos=[]
    ref=[]
    alt=[]
    line_list=h.readlines()[1:]
    for line in line_list:
        tmp_list=line.split()
        chrom.append(tmp_list[0])
        pos.append(tmp_list[1])
        ref.append((tmp_list[2]).split(","))
        alt.append((tmp_list[3]).split(","))
        i+=1
pos=[int(x) for x in pos]
```

### 2.把所有变异分为五类：普通SNP，普通INDEL，复合，长，带星号的，写成新的一列“CnL”。


```python
i=0
CnL=[]  # "0.1" for normal SNP and "0.2" for normal INDEL, "1.1" for complex, "1.2" for long, "1.3" for *. (既是complex又是long的，则优先判断为complex)
star=0
norm_snp=0
norm_indel=0
long=0
compx=0
while i<len(chrom):
    if "*" in alt[i]:
        CnL.append(1.3)
        star+=1
        
    #These are INDELs.
    elif len(max(ref[i]+alt[i], key=len, default='1'))>1:
        ##complex
        if "*" in alt[i+1] or (len(ref[i][0])<len(max(alt[i], key=len)) and len(min(alt[i],key=len))==1):
            CnL.append(1.1)
            compx+=1
        ##long
        elif len(max(ref[i]+alt[i], key=len, default='1'))>4:
            CnL.append(1.2)
            long+=1
        ##normal
        else:
            CnL.append(0.2)
            norm_indel+=1
        
#These are normal SNPs.
    else: 
        CnL.append(0.1)
        norm_snp+=1
    i+=1
    
print("normal SNP: ",norm_snp,", normal INDEL: ",norm_indel,", Complex Indel: ",compx,", long INDEL: ",long,", star: ",star)
```

    normal SNP:  26746475 , normal INDEL:  3661658 , Complex Indel:  2516666 , long INDEL:  2249629 , star:  4381171


### 3.往下删除10bp。第一层循环：根据当前位点是否为complex/long；第二层循环：根据其下游的位点是否在10bp内。 


```python
i=0
j=0
maxPos=0
newPos=0
while i<len(CnL):
    if CnL[i]==0.1 or CnL[i]==0.2 or CnL[i]==1.3:
        pass
    
    elif CnL[i]==1.1 or CnL[i]==1.2:    #Complex/Long INDEL deletion.
        maxPos=pos[i]+len(ref[i][0])-1
        j=i+1
        while CnL[j]==1.3:
            newPos=pos[j]+len(ref[j][0])-1
            if newPos>maxPos:
                maxPos=newPos    #用最远的那个
            
            j+=1
            
            if j>len(CnL)-1:
                break
        
        while chrom[i]==chrom[j] and (CnL[j]==0.1 or CnL[j]==0.2) and pos[j]<=maxPos+10:
            CnL[j]=1.4   # For the normal SNP or INDEL which is within 10 bp downstream of complex/long, we assign the number "1.4".
            j+=1
            if j>len(CnL)-1:
                break
    i+=1            
```

### 4.相同方式，倒着数10bp


```python
i=len(CnL)-1
j=0
while i>0:
    if CnL[i]==0.1 or CnL[i]==0.2 or CnL[i]==1.3 or CnL[i]==1.4:
        pass
    elif CnL[i]==1.1 or CnL[i]==1.2:
        j=i-1
        while chrom[i]==chrom[j] and (CnL[j]==0.1 or CnL[j]==0.2 or CnL[j]==1.4) and pos[j]>=pos[i]-10:
            CnL[j]=1.5   # For the normal SNP or INDEL which is within 10 bp upstream of complex/long, we assign the number "1.5".
            j-=1
            if j<0:
                break    #防止编号j溢出label
    i-=1         
```

### 5. 第一步的统计


```python
kept=0
down10=0
up10=0
for snp in CnL:
    if snp==0.1 or snp==0.2:
        kept+=1
    elif snp==1.4:
        down10+=1
    elif snp==1.5:
        up10+=1
        
print("In step 1, we remove ",len(CnL)-kept," variants (including",down10+up10,"are broadened to 10bp around), ",kept," kept.")
```

    In step 1, we remove  17524292  variants (including 8376826 are broadened to 10bp around),  22031307  kept.


## 二、 用第一步剩下的位点，删除高密度。


```python
#之前是相当于每10bp取一个SNP，现在是相当于把高密度（>10bp）的区域都删了（留下第一个）。
#另起一列标签“HighDens”，高密度区域写成“1”，其它则为“0”。

HighDens=[]
for i in range(0,len(pos)):
    HighDens.append(0)

i=0
j=i+1
start = time.time()
while i<len(pos)-2:
    if CnL[i]==0.1 or CnL[i]==0.2 or HighDens[i]==1:
        while CnL[j]==1.1 or CnL[j]==1.2 or CnL[j]==1.3 or CnL[j]==1.4 or CnL[j]==1.5:    ##Skip which are already marked as removed.
            j+=1
            if j>=len(pos)-2:
                break
        if chrom[j]==chrom[i] and pos[j]<pos[i]+len(ref[i][0])+9:
            HighDens[j]=1
    else:
        pass
    i=j
    j=i+1
    if i>=len(pos)-2:
        break

end = time.time()

t=end-start

print("Runtime is ：",t) 
```

    Runtime is ： 31.42001247406006


### 第二步的统计


```python
dense=0
for i in range(0,len(pos)):
    if HighDens[i]==1:
        dense+=1
print("In step 2, we remove ",dense)
```

    In step 2, we remove  12823311


### 部分结果展示


```python
#查看 step 2 部分筛选结果
for i in range(0,500):
     print(chrom[i],pos[i],CnL[i],HighDens[i],ref[i],alt[i])
```

    scaffold1 488 0.1 0 ['T'] ['C']
    scaffold1 506 0.1 0 ['C'] ['T']
    scaffold1 753 0.1 0 ['T'] ['G']
    scaffold1 900 0.1 0 ['A'] ['C']
    scaffold1 1029 1.2 0 ['A'] ['AGGGACCT']
    scaffold1 1032 1.4 0 ['G'] ['GT']
    scaffold1 1085 1.1 0 ['GGGTAGGTAGGTAGGTAGTTCGGTAGGTAAGTA'] ['G']
    scaffold1 1097 1.3 0 ['A'] ['C', '*']
    scaffold1 1106 1.3 0 ['GGTAGGTAA'] ['*', 'G']
    scaffold1 1114 1.3 0 ['A'] ['*', 'G']
    scaffold1 1155 0.1 0 ['G'] ['T']
    scaffold1 1157 0.1 1 ['A'] ['C']
    scaffold1 1174 0.1 0 ['G'] ['A']
    scaffold1 1197 0.1 0 ['G'] ['C']
    scaffold1 1213 0.1 0 ['A'] ['T']
    scaffold1 1217 0.1 1 ['G'] ['T']
    scaffold1 1222 0.1 1 ['A'] ['T']
    scaffold1 1248 0.1 0 ['C'] ['T']
    scaffold1 1252 0.1 1 ['G'] ['A']
    scaffold1 1255 0.1 1 ['A'] ['C']
    scaffold1 1256 0.1 1 ['C'] ['G']
    scaffold1 1257 0.1 1 ['T'] ['C']
    scaffold1 1260 0.1 1 ['C'] ['G']
    scaffold1 1266 0.1 1 ['C'] ['A']
    scaffold1 1268 0.1 1 ['C'] ['T']
    scaffold1 1269 0.1 1 ['A'] ['T']
    scaffold1 1271 0.1 1 ['T'] ['G']
    scaffold1 1282 0.1 0 ['G'] ['C']
    scaffold1 1295 0.1 0 ['G'] ['T']
    scaffold1 1307 0.1 0 ['T'] ['A']
    scaffold1 1312 0.1 1 ['A'] ['G']
    scaffold1 1319 0.1 1 ['T'] ['C']
    scaffold1 1325 0.1 1 ['G'] ['T']
    scaffold1 1330 0.1 1 ['G'] ['A']
    scaffold1 1343 0.1 0 ['A'] ['C']
    scaffold1 1349 0.1 1 ['T'] ['C']
    scaffold1 1355 0.1 1 ['A'] ['T']
    scaffold1 1361 0.1 1 ['T'] ['A']
    scaffold1 1363 0.1 1 ['G'] ['A']
    scaffold1 1381 0.1 0 ['T'] ['C']
    scaffold1 1387 0.1 1 ['A'] ['G']
    scaffold1 1388 0.1 1 ['T'] ['C']
    scaffold1 1436 0.1 0 ['G'] ['C']
    scaffold1 1471 0.1 0 ['A'] ['G']
    scaffold1 1478 0.1 1 ['G'] ['T']
    scaffold1 1485 0.1 1 ['G'] ['A']
    scaffold1 1495 0.2 0 ['CAT'] ['C']
    scaffold1 1500 0.1 1 ['T'] ['C', 'G']
    scaffold1 1501 0.1 1 ['T'] ['C']
    scaffold1 1505 0.1 1 ['C'] ['A']
    scaffold1 1508 0.1 1 ['T'] ['G', 'C']
    scaffold1 1511 0.1 1 ['C'] ['T']
    scaffold1 1514 0.1 1 ['T'] ['C']
    scaffold1 1527 0.1 0 ['G'] ['A']
    scaffold1 1529 0.1 1 ['A'] ['T']
    scaffold1 1533 0.1 1 ['C'] ['G']
    scaffold1 1538 0.1 1 ['A'] ['G']
    scaffold1 1541 0.1 1 ['T'] ['C']
    scaffold1 1542 0.1 1 ['T'] ['G']
    scaffold1 1545 0.1 1 ['C'] ['G']
    scaffold1 1548 0.1 1 ['C'] ['A']
    scaffold1 1550 0.1 1 ['G'] ['C']
    scaffold1 1551 0.1 1 ['G'] ['A']
    scaffold1 1556 0.1 1 ['A'] ['G']
    scaffold1 1564 0.1 1 ['G'] ['A']
    scaffold1 1565 0.1 1 ['C'] ['T']
    scaffold1 1571 0.1 1 ['A'] ['T']
    scaffold1 1574 0.1 1 ['A'] ['G']
    scaffold1 1579 0.1 1 ['A'] ['G']
    scaffold1 1582 0.1 1 ['T'] ['G']
    scaffold1 1588 0.1 1 ['C'] ['T']
    scaffold1 1594 0.1 1 ['T'] ['C']
    scaffold1 1600 0.2 1 ['ACAT'] ['A']
    scaffold1 1609 0.1 1 ['C'] ['T']
    scaffold1 1621 0.1 0 ['G'] ['A']
    scaffold1 1625 0.1 1 ['T'] ['C']
    scaffold1 1630 1.5 0 ['C'] ['A']
    scaffold1 1636 1.1 0 ['TTGACTACGGGGTTCGCCCA'] ['T']
    scaffold1 1638 1.3 0 ['G'] ['*', 'GAAAGGA']
    scaffold1 1640 1.3 0 ['C'] ['*', 'CCATG']
    scaffold1 1641 1.3 0 ['TACGGGG'] ['*', 'T']
    scaffold1 1642 1.3 0 ['A'] ['*', 'AGTG']
    scaffold1 1649 1.3 0 ['TCGCCCATG'] ['*', 'T']
    scaffold1 1657 1.3 0 ['G'] ['*', 'A']
    scaffold1 1677 1.2 0 ['T'] ['TTCAACCCTGCA']
    scaffold1 1679 1.5 0 ['G'] ['C', 'A']
    scaffold1 1680 1.5 0 ['A'] ['AAT', 'AAC']
    scaffold1 1681 1.1 0 ['C'] ['A', 'CCCTCAAAAGGTTAAAG']
    scaffold1 1682 1.5 0 ['G'] ['C', 'A']
    scaffold1 1683 1.1 0 ['C'] ['CCTGCATGACCCTCAAAAGGTTAAAGGT', 'CCTGCATGACCCTCAAAAGGTTAACGGT', 'A']
    scaffold1 1684 1.4 0 ['C'] ['A', 'G']
    scaffold1 1685 1.4 0 ['A'] ['G']
    scaffold1 1686 1.4 0 ['T'] ['C']
    scaffold1 1712 0.1 0 ['T'] ['G']
    scaffold1 1763 0.1 0 ['C'] ['G']
    scaffold1 1777 0.1 0 ['A'] ['G']
    scaffold1 1847 0.1 0 ['G'] ['C']
    scaffold1 1967 0.1 0 ['G'] ['T']
    scaffold1 2028 0.1 0 ['T'] ['A']
    scaffold1 2073 1.2 0 ['A'] ['AGGTCTCGGACCGGAGTTTTTTTCACGATTTTGGATCCGGCGGACACCCT']
    scaffold1 2102 0.1 0 ['C'] ['A']
    scaffold1 2124 0.1 0 ['A'] ['T']
    scaffold1 2205 1.2 0 ['GTACATT'] ['G']
    scaffold1 2452 0.1 0 ['C'] ['T']
    scaffold1 2460 0.1 1 ['G'] ['A']
    scaffold1 2473 0.1 0 ['G'] ['A']
    scaffold1 2480 0.2 1 ['TCTC'] ['T']
    scaffold1 2499 0.2 0 ['AG'] ['A']
    scaffold1 2504 0.1 1 ['T'] ['A']
    scaffold1 2512 1.5 0 ['CT'] ['C']
    scaffold1 2516 1.5 0 ['T'] ['A']
    scaffold1 2522 1.2 0 ['T'] ['TGAACCCGC']
    scaffold1 2550 1.2 0 ['C'] ['CTTACCTCT']
    scaffold1 2585 0.1 0 ['T'] ['C']
    scaffold1 2589 0.1 1 ['T'] ['C']
    scaffold1 2591 0.2 1 ['GGAA'] ['G']
    scaffold1 2610 0.1 0 ['A'] ['C']
    scaffold1 2617 0.1 1 ['T'] ['G']
    scaffold1 2628 0.1 0 ['C'] ['T']
    scaffold1 2632 0.1 1 ['C'] ['T']
    scaffold1 2633 0.1 1 ['A'] ['G']
    scaffold1 2797 0.1 0 ['G'] ['C']
    scaffold1 2798 0.1 1 ['T'] ['A']
    scaffold1 2803 0.1 1 ['A'] ['T']
    scaffold1 2809 1.5 0 ['A'] ['G']
    scaffold1 2810 1.5 0 ['A'] ['T']
    scaffold1 2815 1.5 0 ['C'] ['G']
    scaffold1 2816 1.2 0 ['AACTT'] ['A']
    scaffold1 2824 1.4 0 ['T'] ['TC']
    scaffold1 2826 1.4 0 ['GC'] ['G']
    scaffold1 2828 1.4 0 ['C'] ['A']
    scaffold1 2831 0.2 0 ['GA'] ['G']
    scaffold1 2838 0.1 1 ['C'] ['T']
    scaffold1 2842 0.1 1 ['A'] ['G']
    scaffold1 2843 0.1 1 ['C'] ['A']
    scaffold1 2846 0.2 1 ['T'] ['TTA']
    scaffold1 2853 0.2 1 ['T'] ['TCA']
    scaffold1 2854 0.2 1 ['TG'] ['T']
    scaffold1 2859 0.1 1 ['G'] ['A']
    scaffold1 2866 0.1 1 ['G'] ['A']
    scaffold1 2871 0.1 1 ['C'] ['A']
    scaffold1 2876 0.2 1 ['AC'] ['A']
    scaffold1 2878 0.2 1 ['ACT'] ['A']
    scaffold1 2882 0.2 1 ['CT'] ['C']
    scaffold1 2913 0.1 0 ['C'] ['T']
    scaffold1 2988 0.1 0 ['T'] ['C']
    scaffold1 2994 0.1 1 ['C'] ['T']
    scaffold1 3003 0.1 1 ['T'] ['G']
    scaffold1 3233 1.2 0 ['C'] ['CTACTTAG']
    scaffold1 3343 0.1 0 ['T'] ['A']
    scaffold1 3350 0.1 1 ['G'] ['A']
    scaffold1 3359 0.1 1 ['C'] ['T']
    scaffold1 3410 0.1 0 ['A'] ['G']
    scaffold1 3413 0.1 1 ['C'] ['T']
    scaffold1 3422 0.1 1 ['C'] ['T']
    scaffold1 3452 1.2 0 ['TGTACG'] ['T']
    scaffold1 3461 1.4 0 ['T'] ['A']
    scaffold1 3464 1.5 0 ['T'] ['G']
    scaffold1 3473 1.2 0 ['CAATCGG'] ['C']
    scaffold1 3481 1.2 0 ['TAGACAAGAAG'] ['T']
    scaffold1 3492 1.2 0 ['TAAAAAGCAGCGG'] ['T']
    scaffold1 3505 1.2 0 ['TCCCAAC'] ['T']
    scaffold1 3512 1.2 0 ['TGGAAAGA'] ['T']
    scaffold1 3522 1.4 0 ['G'] ['GA']
    scaffold1 3597 1.2 0 ['CACGTG'] ['C']
    scaffold1 3644 0.1 0 ['C'] ['T']
    scaffold1 3656 1.2 0 ['C'] ['CCTTAGGCCATGTTGATTTGATTATATGGATGACATCCGCGCGCGCATGAA', 'CCTTAGGCCATGTTGATTTGATTATATGGATGACATCCGCGCGCGCATGAATTTTCGGCCG']
    scaffold1 3657 1.5 0 ['T'] ['C']
    scaffold1 3660 1.1 0 ['A'] ['AGGCCATGTTGATTTGATTATATGGATGACAT', 'AGG', 'T', 'AGGCCATGTTGATTTGATTATATGGATGACATCCGCGCGCGCATGAATTTTCGG', 'AG']
    scaffold1 3661 1.5 0 ['C'] ['G']
    scaffold1 3662 1.1 0 ['C'] ['CATGTTGATTTGATTATATGGATGA', 'G']
    scaffold1 3665 1.2 0 ['T'] ['TCCGCGCGCGCA']
    scaffold1 3667 1.2 0 ['G'] ['GAATTTTCGGCCGTTTC']
    scaffold1 3721 0.1 0 ['C'] ['T']
    scaffold1 3742 1.5 0 ['T'] ['C']
    scaffold1 3745 1.2 0 ['C'] ['CATATGAATATATGGCATTAACAATGATGTTTATTTTTATGCTCTTTTATACATGTATTCACTGTTTTTCTCTTTG']
    scaffold1 3758 0.1 0 ['C'] ['T']
    scaffold1 3806 0.1 0 ['C'] ['T']
    scaffold1 3815 0.1 1 ['G'] ['C']
    scaffold1 3824 1.5 0 ['C'] ['T']
    scaffold1 3833 1.5 0 ['G'] ['GT']
    scaffold1 3834 1.1 0 ['GCC'] ['G']
    scaffold1 3835 1.3 0 ['C'] ['*', 'CAAAT']
    scaffold1 3838 1.5 0 ['G'] ['C']
    scaffold1 3839 1.5 0 ['T'] ['TGGC']
    scaffold1 3840 1.5 0 ['C'] ['G']
    scaffold1 3841 1.1 0 ['TCA'] ['T']
    scaffold1 3843 1.3 0 ['AT'] ['*', 'A']
    scaffold1 3845 1.2 0 ['ATGAC'] ['A']
    scaffold1 3861 0.2 0 ['AT'] ['A']
    scaffold1 3863 1.5 0 ['G'] ['C']
    scaffold1 3869 1.5 0 ['A'] ['ACC']
    scaffold1 3873 1.1 0 ['AAT'] ['A', 'ATCAT']
    scaffold1 3877 1.5 0 ['T'] ['A']
    scaffold1 3878 1.1 0 ['ATTGT'] ['A']
    scaffold1 3881 1.3 0 ['G'] ['A', '*']
    scaffold1 3882 1.3 0 ['T'] ['A', '*']
    scaffold1 3894 0.1 0 ['T'] ['A']
    scaffold1 3895 0.1 1 ['G'] ['A']
    scaffold1 3910 1.5 0 ['C'] ['A']
    scaffold1 3912 1.5 0 ['G'] ['GGAT']
    scaffold1 3913 1.2 0 ['AGTGT'] ['A']
    scaffold1 3920 1.5 0 ['T'] ['G']
    scaffold1 3926 1.5 0 ['T'] ['G']
    scaffold1 3928 1.2 0 ['GTTGT'] ['G']
    scaffold1 3948 0.1 0 ['G'] ['C']
    scaffold1 3960 0.1 0 ['A'] ['G']
    scaffold1 3973 0.1 0 ['G'] ['A']
    scaffold1 3989 0.1 0 ['C'] ['A']
    scaffold1 3993 0.1 1 ['C'] ['A']
    scaffold1 4011 1.5 0 ['A'] ['G']
    scaffold1 4015 1.2 0 ['C'] ['CTGCCCCATTAGCCTGGGTACCATCCGGATTTTCGCTCCGCAACGCTCGAAAAAATCACTTTCTAGCATTGCGGAGTGAAAATCCGAATGGTACCCAGGCTA', 'CTGCCCCATTAGCCGGGGTACCATCCGGATTTTTGCTCCGCAACGCTCGAAAAAATCACTTTCTAGCATTGCGGAGTGAAAATCCGAATGGTACCCAGGCTA', 'CTGCCCCATTAGCCGGGGTACCATCCGGATTTTTGCTCCGCAACGCTCGAAAAAATCACTTTCTAGCATTGCGGAGTGAAAATCCGAATGGTTCGCAGGCTA', 'CTGCCCCATTAGCCTGGGTACCATCCGGATTTTCGCTCCGCAACGCTCGAAAAAATCACTTTCTAGCATTGCGGAGTGAAAATCCGAATGGTTCGCAGGCTA', 'CTGCCCCATTAGCCGGGGTACCATCCGGATTTTTGCTCCGCAACGCTCGAAAAAATCACTTTCTAGCATTGCGGAGAGAAAATCCGAATGGTACCCAGGCTA']
    scaffold1 4046 1.5 0 ['C'] ['T']
    scaffold1 4047 1.2 0 ['ATACAAATG'] ['A']
    scaffold1 4071 0.2 0 ['CCA'] ['C']
    scaffold1 4094 0.1 0 ['T'] ['G']
    scaffold1 4095 0.1 1 ['A'] ['G']
    scaffold1 4105 0.1 0 ['T'] ['C']
    scaffold1 4110 0.1 1 ['A'] ['G']
    scaffold1 4210 0.1 0 ['C'] ['T']
    scaffold1 4216 0.1 1 ['T'] ['C']
    scaffold1 4260 0.2 0 ['C'] ['CTA']
    scaffold1 4278 1.5 0 ['A'] ['G']
    scaffold1 4280 1.5 0 ['ACT'] ['A']
    scaffold1 4283 1.2 0 ['AACAG'] ['A']
    scaffold1 4294 1.4 0 ['G'] ['C']
    scaffold1 4308 0.1 0 ['C'] ['T']
    scaffold1 4319 1.2 0 ['A'] ['ACTCACACTTACTATCCTGGTTTATCCAGTCT']
    scaffold1 4320 1.4 0 ['T'] ['C']
    scaffold1 4325 1.4 0 ['C'] ['T']
    scaffold1 4328 1.4 0 ['T'] ['A']
    scaffold1 4337 0.1 0 ['T'] ['C']
    scaffold1 4338 0.1 1 ['T'] ['C']
    scaffold1 4343 0.1 1 ['C'] ['T']
    scaffold1 4354 0.1 0 ['T'] ['C']
    scaffold1 4403 1.1 0 ['TTAA'] ['T']
    scaffold1 4405 1.3 0 ['A'] ['*', 'T']
    scaffold1 4411 1.4 0 ['G'] ['T']
    scaffold1 4417 0.1 0 ['G'] ['T']
    scaffold1 4420 0.2 1 ['CT'] ['C']
    scaffold1 4428 0.1 1 ['A'] ['G']
    scaffold1 4436 0.1 1 ['T'] ['G']
    scaffold1 4441 0.2 1 ['CTA'] ['C']
    scaffold1 4457 0.1 0 ['G'] ['A']
    scaffold1 4476 0.1 0 ['A'] ['G']
    scaffold1 4490 0.1 0 ['T'] ['C']
    scaffold1 4497 1.5 0 ['A'] ['G']
    scaffold1 4498 1.5 0 ['T'] ['C']
    scaffold1 4500 1.5 0 ['G'] ['T']
    scaffold1 4503 1.1 0 ['CA'] ['C', 'CATTTTTCTA']
    scaffold1 4506 1.4 0 ['TGGC'] ['T']
    scaffold1 4511 1.4 0 ['A'] ['C']
    scaffold1 4512 1.4 0 ['A'] ['C']
    scaffold1 4521 0.1 0 ['A'] ['T']
    scaffold1 4526 0.1 1 ['T'] ['G']
    scaffold1 4531 0.1 1 ['T'] ['A']
    scaffold1 4532 0.1 1 ['G'] ['A']
    scaffold1 4559 0.1 0 ['C'] ['T']
    scaffold1 4599 0.1 0 ['T'] ['C']
    scaffold1 4604 0.1 1 ['A'] ['C']
    scaffold1 4607 0.1 1 ['C'] ['T']
    scaffold1 4608 0.1 1 ['A'] ['G']
    scaffold1 4613 0.1 1 ['T'] ['G']
    scaffold1 4616 0.1 1 ['C'] ['T']
    scaffold1 4618 0.1 1 ['G'] ['C']
    scaffold1 4630 0.1 0 ['A'] ['G']
    scaffold1 4634 0.1 1 ['T'] ['C']
    scaffold1 4636 0.1 1 ['A'] ['G']
    scaffold1 4639 0.1 1 ['C'] ['A']
    scaffold1 4645 1.5 0 ['A'] ['G']
    scaffold1 4647 1.5 0 ['A'] ['G']
    scaffold1 4650 1.5 0 ['A'] ['C']
    scaffold1 4652 1.2 0 ['ATTAAGCTG'] ['A']
    scaffold1 4665 1.4 0 ['C'] ['G']
    scaffold1 4667 1.4 0 ['A'] ['T']
    scaffold1 4670 1.4 0 ['G'] ['T']
    scaffold1 4680 0.1 0 ['G'] ['C']
    scaffold1 4691 1.1 0 ['GACTT'] ['G']
    scaffold1 4693 1.3 0 ['C'] ['*', 'G']
    scaffold1 4694 1.3 0 ['T'] ['*', 'C']
    scaffold1 4697 1.1 0 ['GGTGAA'] ['G']
    scaffold1 4699 1.3 0 ['T'] ['*', 'A']
    scaffold1 4700 1.3 0 ['G'] ['*', 'A']
    scaffold1 4706 1.4 0 ['T'] ['C']
    scaffold1 4724 1.2 0 ['GTAAA'] ['G']
    scaffold1 4731 1.5 0 ['C'] ['CT']
    scaffold1 4732 1.5 0 ['A'] ['G']
    scaffold1 4733 1.5 0 ['T'] ['A']
    scaffold1 4735 1.5 0 ['G'] ['GCA']
    scaffold1 4737 1.5 0 ['G'] ['A', 'C']
    scaffold1 4739 1.5 0 ['C'] ['CTT']
    scaffold1 4740 1.2 0 ['T'] ['TGTTA']
    scaffold1 4741 1.1 0 ['CAACCAG'] ['C']
    scaffold1 4744 1.3 0 ['C'] ['*', 'T']
    scaffold1 4748 1.1 0 ['TGGTAGA'] ['T']
    scaffold1 4751 1.3 0 ['T'] ['*', 'G']
    scaffold1 4756 1.5 0 ['TAG'] ['T']
    scaffold1 4762 1.5 0 ['T'] ['TG']
    scaffold1 4765 1.2 0 ['A'] ['ATGACCG']
    scaffold1 4785 0.1 0 ['C'] ['A']
    scaffold1 4789 1.5 0 ['C'] ['A']
    scaffold1 4795 1.5 0 ['C'] ['A']
    scaffold1 4796 1.1 0 ['CTTTTG'] ['C', 'ATTTTG']
    scaffold1 4799 1.3 0 ['T'] ['*', 'TTC']
    scaffold1 4801 1.3 0 ['G'] ['*', 'T']
    scaffold1 4807 1.4 0 ['T'] ['TTA']
    scaffold1 4815 0.1 0 ['T'] ['C']
    scaffold1 4818 0.1 1 ['T'] ['C']
    scaffold1 4833 0.1 0 ['T'] ['C']
    scaffold1 4905 0.1 0 ['T'] ['C']
    scaffold1 4908 0.1 1 ['C'] ['A']
    scaffold1 4921 0.1 0 ['A'] ['G']
    scaffold1 4969 0.1 0 ['C'] ['T']
    scaffold1 4980 0.1 0 ['A'] ['T']
    scaffold1 5002 1.5 0 ['GT'] ['G']
    scaffold1 5004 1.1 0 ['TTC'] ['T']
    scaffold1 5006 1.3 0 ['CCCCCTTGCAGAATGTTTCAAGCATGTTATTGAATGTTAGTGTGTAAGAAGTTTTG'] ['C', '*']
    scaffold1 5008 1.3 0 ['C'] ['*', 'CAT']
    scaffold1 5011 1.3 0 ['T'] ['*', 'A']
    scaffold1 5016 1.3 0 ['G'] ['*', 'C']
    scaffold1 5034 1.3 0 ['ATT'] ['*', 'A']
    scaffold1 5037 1.3 0 ['G'] ['*', 'GCA']
    scaffold1 5038 1.3 0 ['A'] ['*', 'C']
    scaffold1 5065 1.2 0 ['CAGTTAAGTTGTAAAACAGTACTGGAGTAGAAAGGTGAGATATCCT'] ['C']
    scaffold1 5114 1.5 0 ['GT'] ['G']
    scaffold1 5117 1.2 0 ['ATTGGGGTGAC'] ['A']
    scaffold1 5130 1.2 0 ['ACTGTTTTACCAAGTTCC'] ['A']
    scaffold1 5154 1.2 0 ['ACCCATGCAGAATGTTTTC'] ['A']
    scaffold1 5182 1.4 0 ['A'] ['G']
    scaffold1 5193 0.1 0 ['A'] ['T']
    scaffold1 5308 0.1 0 ['C'] ['T']
    scaffold1 5330 0.1 0 ['T'] ['A']
    scaffold1 5339 0.1 1 ['G'] ['A']
    scaffold1 5342 0.1 1 ['C'] ['T']
    scaffold1 5352 0.2 0 ['CA'] ['C']
    scaffold1 5365 0.1 0 ['T'] ['C']
    scaffold1 5401 0.1 0 ['C'] ['T']
    scaffold1 5404 0.1 1 ['G'] ['T']
    scaffold1 5428 0.1 0 ['T'] ['C']
    scaffold1 5437 0.1 1 ['A'] ['G']
    scaffold1 5455 0.1 0 ['G'] ['A']
    scaffold1 5476 1.5 0 ['T'] ['A']
    scaffold1 5480 1.2 0 ['T'] ['TTTAC']
    scaffold1 5485 1.5 0 ['T'] ['TC']
    scaffold1 5488 1.2 0 ['TCATA'] ['T']
    scaffold1 5507 0.1 0 ['T'] ['G']
    scaffold1 5508 0.1 1 ['G'] ['T']
    scaffold1 5514 0.1 1 ['T'] ['C']
    scaffold1 5525 0.2 0 ['T'] ['TAC']
    scaffold1 5530 1.5 0 ['G'] ['GTA']
    scaffold1 5531 1.5 0 ['C'] ['CAT']
    scaffold1 5532 1.5 0 ['G'] ['GCTT']
    scaffold1 5536 1.5 0 ['A'] ['G']
    scaffold1 5539 1.5 0 ['T'] ['TTGG']
    scaffold1 5540 1.2 0 ['CAGAG'] ['C']
    scaffold1 5546 1.2 0 ['A'] ['AATTCT']
    scaffold1 5547 1.2 0 ['TCAGC'] ['T']
    scaffold1 5554 1.2 0 ['CAAGTAAA'] ['C']
    scaffold1 5563 1.4 0 ['G'] ['C']
    scaffold1 5568 1.4 0 ['CAG'] ['C']
    scaffold1 5571 1.4 0 ['C'] ['CTGG']
    scaffold1 5573 0.1 0 ['T'] ['C']
    scaffold1 5576 0.1 1 ['T'] ['A']
    scaffold1 5579 0.1 1 ['A'] ['G']
    scaffold1 5594 0.1 0 ['C'] ['T']
    scaffold1 5637 0.1 0 ['G'] ['A']
    scaffold1 5701 0.1 0 ['G'] ['A']
    scaffold1 5718 0.1 0 ['T'] ['C']
    scaffold1 5755 0.2 0 ['C'] ['CTG']
    scaffold1 5839 0.2 0 ['ATC'] ['A']
    scaffold1 5858 0.1 0 ['A'] ['T']
    scaffold1 5917 1.1 0 ['TCTC'] ['T', 'CCTC']
    scaffold1 5920 1.3 0 ['C'] ['*', 'T']
    scaffold1 5936 0.1 0 ['T'] ['C']
    scaffold1 5950 0.1 0 ['T'] ['G']
    scaffold1 5967 0.1 0 ['T'] ['C']
    scaffold1 5974 0.1 1 ['C'] ['T']
    scaffold1 5982 0.1 1 ['T'] ['C']
    scaffold1 6001 0.1 0 ['G'] ['A']
    scaffold1 6010 1.5 0 ['C'] ['T']
    scaffold1 6012 1.5 0 ['C'] ['A']
    scaffold1 6014 1.5 0 ['G'] ['T']
    scaffold1 6017 1.1 0 ['G'] ['A', 'GTT']
    scaffold1 6025 1.4 0 ['A'] ['C']
    scaffold1 6048 0.1 0 ['G'] ['A']
    scaffold1 6078 0.1 0 ['T'] ['G']
    scaffold1 6102 0.1 0 ['A'] ['T']
    scaffold1 6119 0.1 0 ['C'] ['A']
    scaffold1 6173 0.1 0 ['A'] ['G']
    scaffold1 6191 0.1 0 ['T'] ['A']
    scaffold1 6209 0.1 0 ['A'] ['G']
    scaffold1 6228 0.1 0 ['A'] ['G']
    scaffold1 6232 0.1 1 ['A'] ['T', 'G']
    scaffold1 6266 0.1 0 ['G'] ['T']
    scaffold1 6267 0.1 1 ['G'] ['A']
    scaffold1 6273 0.1 1 ['G'] ['A']
    scaffold1 6297 0.1 0 ['T'] ['C']
    scaffold1 6299 0.1 1 ['A'] ['T']
    scaffold1 6309 0.1 0 ['C'] ['G']
    scaffold1 6337 0.1 0 ['G'] ['A']
    scaffold1 6341 0.1 1 ['C'] ['T']
    scaffold1 6358 1.1 0 ['CAAG'] ['C', 'AAAG']
    scaffold1 6361 1.3 0 ['G'] ['C', '*']
    scaffold1 6403 0.1 0 ['G'] ['T']
    scaffold1 6415 0.2 0 ['CA'] ['C']
    scaffold1 6428 0.1 0 ['G'] ['C']
    scaffold1 6523 0.1 0 ['G'] ['C']
    scaffold1 6534 0.1 0 ['T'] ['C']
    scaffold1 6558 0.1 0 ['G'] ['A', 'T']
    scaffold1 6560 0.1 1 ['T'] ['C']
    scaffold1 6571 1.5 0 ['GT'] ['G']
    scaffold1 6574 1.1 0 ['G'] ['A', 'GTAA']
    scaffold1 6575 1.5 0 ['T'] ['A']
    scaffold1 6576 1.5 0 ['C'] ['A', 'T']
    scaffold1 6578 1.5 0 ['G'] ['C']
    scaffold1 6579 1.1 0 ['T'] ['A', 'TG']
    scaffold1 6581 1.4 0 ['G'] ['C']
    scaffold1 6584 1.4 0 ['C'] ['A']
    scaffold1 6585 1.4 0 ['C'] ['T']
    scaffold1 6587 1.4 0 ['C'] ['A']
    scaffold1 6981 1.2 0 ['C'] ['CAGTGATTG']
    scaffold1 7051 0.1 0 ['T'] ['G']
    scaffold1 7155 0.1 0 ['G'] ['T']
    scaffold1 7225 0.1 0 ['A'] ['C']
    scaffold1 7233 1.5 0 ['C'] ['CTT']
    scaffold1 7238 1.1 0 ['TG'] ['T']
    scaffold1 7239 1.3 0 ['G'] ['*', 'T']
    scaffold1 7244 1.4 0 ['G'] ['T']
    scaffold1 7259 1.5 0 ['G'] ['A']
    scaffold1 7261 1.5 0 ['A'] ['G']
    scaffold1 7264 1.5 0 ['A'] ['ATG']
    scaffold1 7265 1.2 0 ['C'] ['CAGTTT']
    scaffold1 7267 1.5 0 ['A'] ['AG']
    scaffold1 7270 1.2 0 ['T'] ['TGATA']
    scaffold1 7273 1.5 0 ['T'] ['TC']
    scaffold1 7275 1.2 0 ['T'] ['TTAAGGG']
    scaffold1 7277 1.4 0 ['TC'] ['T']
    scaffold1 7282 1.4 0 ['C'] ['A']
    scaffold1 7284 1.4 0 ['C'] ['T']
    scaffold1 7286 0.1 0 ['G'] ['C']
    scaffold1 7288 0.1 1 ['G'] ['A']
    scaffold1 7306 0.1 0 ['G'] ['A']
    scaffold1 7307 0.1 1 ['C'] ['A']
    scaffold1 7325 0.1 0 ['C'] ['T']
    scaffold1 7336 1.1 0 ['AAG'] ['A']
    scaffold1 7337 1.3 0 ['A'] ['*', 'G']
    scaffold1 7344 1.4 0 ['G'] ['A']
    scaffold1 7345 1.4 0 ['T'] ['A']
    scaffold1 7368 0.1 0 ['A'] ['T']
    scaffold1 7373 0.1 1 ['A'] ['G']
    scaffold1 7378 1.5 0 ['T'] ['C']
    scaffold1 7384 1.1 0 ['CGTTTG'] ['CCAAAGTTTG', 'C']
    scaffold1 7389 1.3 0 ['G'] ['A', '*']
    scaffold1 7390 1.4 0 ['C'] ['G']
    scaffold1 7407 0.2 0 ['AG'] ['A']
    scaffold1 7409 0.2 1 ['AG'] ['A']
    scaffold1 7420 0.1 0 ['A'] ['C']
    scaffold1 7435 1.2 0 ['A'] ['ATTAATC']
    scaffold1 7439 1.4 0 ['AC'] ['A']
    scaffold1 7442 1.4 0 ['C'] ['G']
    scaffold1 7475 0.1 0 ['G'] ['A']
    scaffold1 7480 0.1 1 ['A'] ['G']
    scaffold1 7500 1.2 0 ['TCTTTACATATC'] ['T']
    scaffold1 7520 1.4 0 ['T'] ['A']
    scaffold1 7566 0.1 0 ['G'] ['T']
    scaffold1 7621 0.1 0 ['G'] ['A']
    scaffold1 7635 0.1 0 ['T'] ['C']
    scaffold1 7642 0.1 1 ['A'] ['G']
    scaffold1 7644 0.1 1 ['T'] ['A']
    scaffold1 7645 0.1 1 ['G'] ['A']
    scaffold1 7653 0.1 1 ['C'] ['T']
    scaffold1 7660 0.1 1 ['T'] ['G']
    scaffold1 7664 0.1 1 ['T'] ['A', 'C']
    scaffold1 7666 0.1 1 ['T'] ['C']
    scaffold1 7680 0.1 0 ['G'] ['A']
    scaffold1 7691 0.1 0 ['T'] ['G']
    scaffold1 7702 0.1 0 ['T'] ['C']
    scaffold1 7715 0.1 0 ['C'] ['T']
    scaffold1 7751 0.1 0 ['T'] ['C']
    scaffold1 7758 0.1 1 ['G'] ['A']
    scaffold1 7762 0.1 1 ['T'] ['G']
    scaffold1 7767 0.1 1 ['G'] ['C']
    scaffold1 7772 1.5 0 ['TTG'] ['T']
    scaffold1 7777 1.5 0 ['C'] ['G']
    scaffold1 7778 1.1 0 ['G'] ['C', 'GGC']
    scaffold1 7836 0.1 0 ['T'] ['C']
    scaffold1 7846 0.1 0 ['T'] ['C']
    scaffold1 7849 0.2 1 ['C'] ['CTGA']
    scaffold1 7853 0.1 1 ['C'] ['T']
    scaffold1 7883 0.1 0 ['T'] ['C']
    scaffold1 7887 0.1 1 ['T'] ['C']
    scaffold1 7919 0.1 0 ['C'] ['A']
    scaffold1 7937 0.1 0 ['T'] ['C']
    scaffold1 7939 0.1 1 ['T'] ['C']
    scaffold1 7962 0.1 0 ['A'] ['C']
    scaffold1 7969 0.1 1 ['G'] ['A']
    scaffold1 7993 0.1 0 ['G'] ['T']
    scaffold1 7996 0.1 1 ['A'] ['G']
    scaffold1 7999 0.2 1 ['G'] ['GA']
    scaffold1 8012 1.2 0 ['C'] ['CAAAAAAAAAA', 'CAAAAAAAAAAAA', 'CAAAAAAAA', 'CAAAAAAAAA']

