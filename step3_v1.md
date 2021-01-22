```python
import shutil,sys,os,time
```

## 一、读取52_bbe_copy.vcf、INDEL_filter.vcf、SNP_filter.vcf、NO_QD_INDEL.vcf、NO_QD_SNP.vcf


```python
with open("Extracted.vcf","r") as h:  #with open("52_bbe_copy.vcf" ,"r") as h:
    count=0
    line_list=[]
    tmp_list=[]
    total_chrom=[]
    total_pos=[]
    line_list=h.readlines()[1:]
    for line in line_list:
            tmp_list=line.split()
            total_chrom.append(tmp_list[0])
            total_pos.append(tmp_list[1])
            count+=1
total_pos=[int(x) for x in total_pos]
print("Total variants: ",count)
```

    Total variants:  39555599



```python
with open("../VariantFiltration/52_bbe_copy_INDEL_filter.vcf","r") as h:
    count=0
    line_list=[]
    tmp_list=[]
    indel_chrom=[]
    indel_pos=[]
    indelQD=[]
    line_list=h.readlines()
    for line in line_list:
        if not line.startswith('#'):
            tmp_list=line.split()
            indel_chrom.append(tmp_list[0])
            indel_pos.append(tmp_list[1])
            indelQD.append(tmp_list[6])
            count+=1
indel_pos=[int(x) for x in indel_pos]
print("Counts of indel defined by GATK: ",count)
```

    Counts of indel defined by GATK:  8411651



```python
with open("../VariantFiltration/52_bbe_copy_SNP_filter.vcf","r") as h:
    count=0
    line_list=[]
    tmp_list=[]
    snp_chrom=[]
    snp_pos=[]
    snpQD=[]
    line_list=h.readlines()
    for line in line_list:
        if not line.startswith('#'):
            tmp_list=line.split()
            snp_chrom.append(tmp_list[0])
            snp_pos.append(tmp_list[1])
            snpQD.append(tmp_list[6])
            count+=1
snp_pos=[int(x) for x in snp_pos]
print("Counts of SNP defined by GATK: ",count)
```

    Counts of SNP defined by GATK:  29444171



```python
with open("../VariantFiltration/NO_QD_INDEL.vcf","r") as h:
    count=0
    line_list=[]
    tmp_list=[]
    noQD_indel_chrom=[]
    noQD_indel_pos=[]
    line_list=h.readlines()
    for line in line_list:
        if not line.startswith('#'):
            tmp_list=line.split()
            noQD_indel_chrom.append(tmp_list[0])
            noQD_indel_pos.append(tmp_list[1])
            count+=1
noQD_indel_pos=[int(x) for x in noQD_indel_pos]
print("Counts of INDEL without QD-value: ",count)

with open("../VariantFiltration/NO_QD_SNP.vcf","r") as h:
    count=0
    line_list=[]
    tmp_list=[]
    noQD_snp_chrom=[]
    noQD_snp_pos=[]
    line_list=h.readlines()
    for line in line_list:
        if not line.startswith('#'):
            tmp_list=line.split()
            noQD_snp_chrom.append(tmp_list[0])
            noQD_snp_pos.append(tmp_list[1])
            count+=1
noQD_snp_pos=[int(x) for x in noQD_snp_pos]
print("Counts of SNP without QD-value: ",count)
```

    Counts of INDEL without QD-value:  279
    Counts of SNP without QD-value:  321


## 二、一些统计


```python
indel_filter=0
for QD in indelQD:
    if QD=="Filter":
        indel_filter+=1
print("Of all the indels, %s are labeled as 'PASS', while %s are labeled as 'Filter'."%(len(indelQD)-indel_filter,indel_filter))
snp_filter=0
for QD in snpQD:
    if QD=="Filter":
        snp_filter+=1
print("Of all the SNPs, %s are labeled as 'PASS', while %s are labeled as 'Filter'."%(len(snpQD)-snp_filter,snp_filter))
```

    Of all the indels, 8088933 are labeled 'PASS', while 322718 are labeled 'Filter'.
    Of all the indels, 23716734 are labeled 'PASS', while 5727437 are labeled 'Filter'.


## 三、转化成一列标记


```python
#"1" for "Filter" SNP; "2" for "Filter" INDEL;
#"0.1" for 'PASS' and with-QD-value SNP; "0.2" for 'PASS' and with-OD-value INDEL;
#"3" for "PASS" but without-QD-value SNP; "4" for "PASS" but without-QD-value INDEL;
#"99" for variants undetected by GATK.
```


```python
QD=[]
for p in range(len(total_pos)):
    QD.append(99)

i=0    
for j in range(len(snp_pos)):
    while not (total_chrom[i]==snp_chrom[j] and total_pos[i]==snp_pos[j]): 
        i+=1
    if snpQD[j]=="PASS":
        QD[i]=0.1
    elif snpQD[j]=="Filter":
        QD[i]=1
    i+=1

i=0
for j in range(len(indel_pos)):
    while not (total_chrom[i]==indel_chrom[j] and total_pos[i]==indel_pos[j]): 
        i+=1
    if indelQD[j]=="PASS":
        QD[i]=0.2
    elif indelQD[j]=="Filter":
        QD[i]=2
    i+=1

i=0    
for j in range(len(noQD_indel_pos)):
    while not (total_chrom[i]==noQD_indel_chrom[j] and total_pos[i]==noQD_indel_pos[j]): 
        i+=1
    QD[i]=4
    i+=1
    
i=0    
for j in range(len(noQD_snp_pos)):
    while not (total_chrom[i]==noQD_snp_chrom[j] and total_pos[i]==noQD_snp_pos[j]): 
        i+=1
    QD[i]=3
    i+=1 
```

## 四、转换结果


```python
passSNP=passINDEL=yi=er=san=si=na=0
for qd in QD:
    if qd==0.1:
        passSNP+=1
    elif qd==0.2:
        passINDEL+=1
    elif qd==1:
        yi+=1
    elif qd==2:
        er+=1
    elif qd==3:
        san+=1
    elif qd==4:
        si+=1
    else:
        na+=1
print("%s SNP left, %s indel left; %s SNP labeled 'filter', %s indel labeled 'filter'; %s SNP pass \n but no QD, %s indel pass but no QD; %s variants cannot be detected by GATK."%(passSNP,passINDEL,yi,er,san,si,na))
```

    23716413 SNP left, 8088654 indel left; 5727437 SNP labeled 'filter', 322718 indel labeled 'filter'; 321 SNP pass 
     but no QD, 279 indel pass but no QD; 1699777 variants cannot be detected by GATK.



```python
QD[:20]
```




    [0.1,
     0.1,
     1,
     0.1,
     2,
     2,
     0.2,
     0.1,
     0.2,
     0.1,
     0.1,
     0.1,
     0.1,
     0.1,
     1,
     0.1,
     1,
     0.1,
     0.1,
     0.1]


