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

### 2.把所有变异分为五类：普通SNP，普通INDEL，复合，长，带星号的，写成新的一列“label”。


```python
i=0
label=[]  # "0.1" for normal SNP and "0.2" for normal INDEL, "1.1" for complex, "1.2" for long, "1.3" for *. (既是complex又是long的，则优先判断为complex)
star=0
norm_snp=0
norm_indel=0
long=0
compx=0
while i<len(chrom):
    if "*" in alt[i]:
        label.append(1.3)
        star+=1
        
    #These are INDELs.
    elif len(max(ref[i]+alt[i], key=len, default='1'))>1:
        ##complex
        if "*" in alt[i+1] or (len(ref[i][0])<len(max(alt[i], key=len)) and len(min(alt[i],key=len))==1):
            label.append(1.1)
            compx+=1
        ##long
        elif len(max(ref[i]+alt[i], key=len, default='1'))>4:
            label.append(1.2)
            long+=1
        ##normal
        else:
            label.append(0.2)
            norm_indel+=1
        
#These are normal SNPs.
    else: 
        label.append(0.1)
        norm_snp+=1
    i+=1
    
print("normal SNP: ",norm_snp,", normal INDEL: ",norm_indel,", Complex Indel: ",compx,", long INDEL: ",long,", star: ",star)
```

    normal SNP:  26746475 , normal INDEL:  3661658 , Complex Indel:  2516666 , long INDEL:  2249629 , star:  4381171


### 3. 第一层循环：根据当前位点是否为complex/long；第二层循环：根据其下游的位点是否在10bp内。


```python
i=0
j=0
maxPos=0
newPos=0
while i<len(label):
    if label[i]==0.1 or label[i]==0.2 or label[i]==1.3:
        pass
    
    elif label[i]==1.1 or label[i]==1.2:    #Complex/Long INDEL deletion.
        maxPos=pos[i]+len(ref[i][0])-1
        j=i+1
        while label[j]==1.3:
            newPos=pos[j]+len(ref[j][0])-1
            if newPos>maxPos:
                maxPos=newPos    #用最远的那个
            
            j+=1
            
            if j>len(label)-1:
                break
        
        while chrom[i]==chrom[j] and (label[j]==0.1 or label[j]==0.2) and pos[j]<=maxPos+10:
            label[j]=1.4   # For the normal SNP or INDEL which is within 10 bp downstream of complex/long, we assign the number "1.4".
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
    if label[i]==0.1 or label[i]==0.2 or label[i]==1.3 or label[i]==1.4:
        pass
    elif label[i]==1.1 or label[i]==1.2:
        j=i-1
        while chrom[i]==chrom[j] and (label[j]==0.1 or label[j]==0.2 or label[j]==1.4) and pos[j]>=pos[i]-10:
            label[j]=1.5   # For the normal SNP or INDEL which is within 10 bp upstream of complex/long, we assign the number "1.5".
            j-=1
            if j<0:
                break    #防止编号j溢出label
    i-=1         
```

### 5. 统计


```python
kept=0
down10=0
up10=0
for snp in label:
    if snp==0.1 or snp==0.2:
        kept+=1
    elif snp==1.4:
        down10+=1
    elif snp==1.5:
        up10+=1
        
print("In step 1, we remove ",len(label)-kept," variants (including",down10+up10,"are broadened to 10bp around), ",kept," kept.")
```

    In step 1, we remove  17524292  variants (including 8376826 are broadened to 10bp around),  22031307  kept.



```python
#查看 step 1 部分筛选结果
for i in range(0,500):
    print(chrom[i],pos[i],label[i],ref[i],alt[i])
```

## 二、 删除高密度


```python
#之前是相当于每10bp取一个SNP，现在是相当于把高密度（>10bp）的区域都删了（留下第一个）。
i=0
j=i+1
start = time.time()
while i<len(pos)-2:
    if label[i]==0.1 or label[i]==0.2 or label[i]==2:
        while label[j]==1.1 or label[j]==1.2 or label[j]==1.3 or label[j]==1.4 or label[j]==1.5:    ##Skip which are already marked as removed.
            j+=1
            if j>=len(pos)-2:
                break
        if chrom[j]==chrom[i] and pos[j]<pos[i]+len(ref[i][0])+9:
            label[j]=2
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

    Runtime is ： 34.0363347530365


### Statistics for Step 2


```python
dense=0
kept2=0
pos_kept=[]
chrom_kept=[]
for i in range(0,len(pos)):
    if label[i]==0.1 or label[i]==0.2:
        kept2+=1
        pos_kept.append(pos[i])   #顺便存入文件，用R包画密度图。
        chrom_kept.append(chrom[i])
    elif label[i]==2:
        dense+=1
print("In step 2, we remove ",dense," variants, ",kept2," kept.")


np_chrom_kept=np.array(chrom_kept)
np_pos_kept=np.array(pos_kept)
index=np.arange(1,len(pos_kept)+1)
l=np.vstack((index,np_chrom_kept,np_pos_kept)).T
np.savetxt("step1&2kept_addIndex.txt",l,fmt='%s',delimiter='\t')
```

    In step 2, we remove  12823311  variants,  9207996  kept.



```python
#查看 step 2 部分筛选结果
for i in range(0,500):
     print(chrom[i],pos[i],label[i],ref[i],alt[i])
```
