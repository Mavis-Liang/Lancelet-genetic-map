```python
import shutil,sys,os,gzip
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```

### 1. 导入文件，并形成四个list


```python
with open("Extractedaa","r") as h:  #h = gzip.GzipFile("52_bbe_copy.vcf.gz") 
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

### 2.把所有变异分为四类：普通SNP/INDEL，复合/长插入，复合/长缺失，带星号的，写成新的一列“label”。


```python
i=0
label=[]  # "0" for normal SNP and INDEL, "1" for complex/long insertions, "2" for complex/long deletions, "3" for *.

while i<len(chrom):
    if "*" in alt[i]:
        label.append(3)
        
    #These are INDELs.
    elif len(max(ref[i]+alt[i], key=len, default='1'))>1:
        
        #Insertion
        if len(ref[i][0])<len(max(alt[i], key=len)):
            #SNP overlapping INDEL.
            if len(min(alt[i],key=len))==1:
                label.append(1)
            #Long Insertion
            elif len(max(alt[i],key=len))>4:
                label.append(1)
            else:
                label.append(0)
        
        #Deletion
        else:
            if "*" in alt[i+1]:
                label.append(2)
            elif len(ref[i][0])>4:
                label.append(2)
            else:
                label.append(0)
        
#These are normal SNPs.
    else: 
        label.append(0)
    i+=1
```

### 3. 第一层循环：根据当前位点是否为complex/long；第二层循环：根据其下游的位点是否在10bp内。


```python
i=0
j=0
while i<len(label):
    if label[i]==0 or label[i]==3:
        pass
    elif label[i]==1:    # Complex/Long INDEL insertion.
        j=i+1
        while chrom[i]==chrom[j] and label[j]==0 and pos[j]<=pos[i]+10:
            label[j]=4   # For the normal SNP or INDEL which is within 10 bp downstream of complex/long, we assign the number "4".
            j+=1
            if j>len(label)-1:
                break    #防止编号j溢出label
    elif label[i]==2:    #Complex/Long INDEL deletion.
        j=i+1
        while label[j]==3:
            j+=1
            if j>len(label)-1:
                break
        while chrom[i]==chrom[j] and label[j]==0 and pos[j]<=pos[i]+len(ref[i][0])+10:
            label[j]=4
            j+=1
            if j>len(label)-1:
                break
    i+=1            
```

### 4.相同方式，倒着数10bp


```python
i=len(label)-1
j=0
while i>0:
    if label[i]==0 or label[i]==3 or label[i]==4:
        pass
    elif label[i]==1 or label[i]==2:    # This time, we can combine the insertion and deletion, since their starting positions are both the same as their origninal positions.
        j=i-1
        while chrom[i]==chrom[j] and (label[j]==0 or label[j]==4) and pos[j]>=pos[i]-10:
            label[j]=5   # For the normal SNP or INDEL which is within 10 bp upstream of complex/long, we assign the number "5".
            j-=1
            if j<0:
                break    #防止编号j溢出label
    i-=1         
```

### 5.统计


```python
i=0
for snp in label:
    if snp==0:
        i+=1
print("In step 1, we remove ",len(label)-i," variants, ",i," kept.")
```

    In step 1, we remove  919375  variants,  3036184  kept.

