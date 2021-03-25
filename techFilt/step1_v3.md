```python
import shutil,sys,os,gzip
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```

### 1. 导入文件，并形成四个list


```python
with open("52_bbe_copy.vcf","r") as h:  #h = gzip.GzipFile("52_bbe_copy.vcf.gz") 
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
```

### 2.把complex和long拎出来


```python
i=0
var=0   #Counting all the variants
indel=0
ins=0  #Counting insertions.
delete=0   #Counting deletions.
missing=0   # counting *
comp=0  #Counting complex.
long=0
indel_len=[]  # The lengths of the indels.
chrom_long_comp=[]
begin_long_comp=[]  #the beginning positions of long or complex.
end_long_comp=[]   #the ending positions of long or complex.
chrom_missing=[]  #the chrom of the *. they are within the deletions, so we just need to rule out the positions of themselves. 
pos_missing=[]  


while i<len(chrom):
    var+=1
    if "*" in alt[i]:
        chrom_missing.append(chrom[i])
        pos_missing.append(pos[i])
        missing+=1
        
    #These are indels.
    elif len(max(ref[i]+alt[i], key=len, default='1'))>1:
        indel+=1
        
        #Insertion
        if len(ref[i][0])<len(max(alt[i], key=len)):
            ins+=1
            indel_len.append(len(max(alt[i], key=len)))
            #SNPs overlapping Indels.
            if len(min(alt[i],key=len))==1:
                chrom_long_comp.append(chrom[i])
                begin_long_comp.append(pos[i])
                end_long_comp.append(pos[i])
                comp+=1
            #Long insertions
            if len(max(alt[i],key=len))>50:
                chrom_long_comp.append(chrom[i])
                begin_long_comp.append(pos[i])
                end_long_comp.append(pos[i])
                long+=1
        
        #Deletion
        else:
            indel_len.append(len(ref[i][0]))
            delete+=1
            if "*" in alt[i+1]:
                chrom_long_comp.append(chrom[i])
                begin_long_comp.append(pos[i])
                end_long_comp.append(int(pos[i])+len(ref[i][0]))   ###For the deletions, the ending positions should begin after adding lengths to the original positions.
                comp+=1
            if len(ref[0])>50:
                chrom_long_comp.append(chrom[i])
                begin_long_comp.append(pos[i])
                end_long_comp.append(int(pos[i])+len(ref[i][0]))
                long+=1
        
#These are snps.
    else: 
        pass
    i+=1
```

### 3. 一些简单的统计


```python
#Statistics
print("complex indel: ",comp)
print("long indel: ",long)
print("insertion: ",ins)
print("delete: ",delete)
print("total indel: ",indel) #这里的indel是不算带*的indel的。（带*的indel外面已经有一层更大的indel算在里面了）
print("total variants: ",var)
print("max indel: ",max(indel_len))
print("ratio of complex+long and variants:",(comp+long)/var)
print("*: ",missing)
```

    complex indel:  2506099
    long indel:  116599
    insertion:  4474261
    delete:  3953692
    total indel:  8427953
    total variants:  39555599
    max indel:  741
    ratio of complex+long and variants: 0.06630409009859767
    *:  4381171


### 4. 画图


```python
n, bins, patches = plt.hist(x=indel_len,bins=[0.0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100], color='#0504aa',alpha=0.7, rwidth=0.85)
 
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Length of Indels')
plt.ylabel('Counts')
plt.title("Distribution of the Length of Indels")
plt.text(30, 1700000, r'max=741bp, 116,599>50bp')
maxfreq = n.max()
# 设置y轴的上限
plt.ylim(ymax=2700000)
#fig=plt.gcf()
#fig.savefig('./666.jpg')
```




    (0.0, 2700000.0)




![png](output_8_1.png)


### 5. 将起点向上数10bp，将终点向下数10bp；先分别在矩阵里做，然后分别展平成一维array，然后连接起来；scaffold名称也同步进行这些操作；带*的位点则不用这样（到时直接将原位点去除）。
##### PS. 如果不做上下10bp的话，将删除的位点有6,887,270个。


```python
#Transform the lists into numpy arrays.
np_begin_long_comp=np.array(list(map(int,begin_long_comp))).reshape((len(begin_long_comp),1))
np_end_long_comp=np.array(list(map(int,end_long_comp))).reshape((len(end_long_comp),1))
np_chrom_long_comp=np.array(chrom_long_comp).reshape((len(chrom_long_comp),1))

np_pos_missing=np.array(list(map(int,pos_missing)))
np_chrom_missing=np.array(chrom_missing)
```


```python
#Expand the location to +-10bp around.
bp_up_matrix=np.zeros((11,len(np_begin_long_comp)))
for i in range(0,10):
    bp_up_matrix[i] = np.array([i-10]*len(np_begin_long_comp))
bp_up_matrix=bp_up_matrix.astype(int)

bp_down_matrix=np.zeros((11,len(np_end_long_comp)))
for i in range(0,10):
    bp_down_matrix[i] = np.array([i+10]*len(np_begin_long_comp))
bp_down_matrix=bp_down_matrix.astype(int)
```


```python
#Broadcasting
begin_up_matrix=bp_up_matrix.T+np_begin_long_comp
end_down_matrix=bp_down_matrix.T+np_end_long_comp
chrom_matrix=np.tile(chrom_long_comp,(1,11))
```


```python
#flattening
begin_up_array=begin_up_matrix.flatten()
end_down_array=end_down_matrix.flatten()
pos_all=np.append(begin_up_array,np.append(end_down_array,np_pos_missing))
chrom_array=chrom_matrix.flatten()
chrom_all=np.append(chrom_array,np.append(chrom_array,np_chrom_missing))
```


```python
pos_all.shape
```


    (62080527,)


```python
#Write to txt file.
chrom_pos=np.vstack((chrom_all,pos_all)).T
np.savetxt("chrom_pos_Broad.txt",chrom_pos,fmt='%s',delimiter='\t')
```

### 6.用vcftools求chrom_pos_Broad.txt与原文件的交集，即需要删除的位点。第一步删去了8,655,746个(21.9%)。（后面再做成软删除）
