```python
import shutil,sys,os
```

## 1.把fasta中的小写位点写下来


```python
#把fasta中的小写位点写成两个list
with open("Branchiostoma.belcheri_v18h27.r3_ref_genome.softmasked.fa","r") as f:
    chrom=[]
    pos=[]
    linelist=f.read().splitlines()
    ####用i数scaffold，j数行（每换一个scaffold归零一次），z数位点（每换一行归零一次）。
    i=0
    while i<len(linelist):
        j=1
        while not linelist[i+j].startswith('>'):
            for z in range(len(linelist[i+j])):
                if linelist[i+j][z].isupper():
                    pass
                else:
                    chrom.append(linelist[i][1:])
                    pos.append(z+1+(j-1)*50)
            j+=1
            if i+j>=len(linelist):
                break
        i+=j
z=list(zip(chrom,pos))
```

## 2.与vcf中的位点配对筛选


```python
with open("Extracted.vcf","r") as h:  #with open("52_bbe_copy.vcf" ,"r") as h:
    line_list=[]
    tmp_list=[]
    vcf_chrom=[]
    vcf_pos=[]
    line_list=h.readlines()[1:]
    for line in line_list:
            tmp_list=line.split()
            vcf_chrom.append(tmp_list[0])
            vcf_pos.append(tmp_list[1])
vcf_pos=[int(x) for x in vcf_pos]
```


```python
#求交集
z_vcf=list(zip(vcf_chrom,vcf_pos))
intersect=list(set(z_vcf)&set(z))
sorted_intersect=sorted(intersect,key=lambda x:(x[0],x[1]))
```


```python
#写入lcr的标记
lcr=[]
for i in range(len(vcf_pos)):
    lcr.append(0)
i=0    
for j in range(len(sorted_intersect)):
    while not z_vcf[i]==sorted_intersect[j]: 
        i+=1
    lcr[i]=1
    i+=1
```

## 3.统计


```python
print("lcr in fasta: ",len(pos))
count=0
for i in lcr:
    if i==1:
        count+=1
print("lcr in vcf: ",count)
```

    lcr in fasta:  114764956
    lcr in vcf:  11442773



```python
#查看前20个
for i in range(20):
    print(sorted_intersect[i],z_vcf[i],lcr[i])
```

    ('scaffold1', 488) ('scaffold1', 488) 1
    ('scaffold1', 506) ('scaffold1', 506) 1
    ('scaffold1', 753) ('scaffold1', 753) 1
    ('scaffold1', 900) ('scaffold1', 900) 1
    ('scaffold1', 1029) ('scaffold1', 1029) 1
    ('scaffold1', 1032) ('scaffold1', 1032) 1
    ('scaffold1', 1097) ('scaffold1', 1085) 0
    ('scaffold1', 1106) ('scaffold1', 1097) 1
    ('scaffold1', 1114) ('scaffold1', 1106) 1
    ('scaffold1', 1155) ('scaffold1', 1114) 1
    ('scaffold1', 1157) ('scaffold1', 1155) 1
    ('scaffold1', 1174) ('scaffold1', 1157) 1
    ('scaffold1', 1197) ('scaffold1', 1174) 1
    ('scaffold1', 1213) ('scaffold1', 1197) 1
    ('scaffold1', 1217) ('scaffold1', 1213) 1
    ('scaffold1', 1222) ('scaffold1', 1217) 1
    ('scaffold1', 1248) ('scaffold1', 1222) 1
    ('scaffold1', 1252) ('scaffold1', 1248) 1
    ('scaffold1', 1255) ('scaffold1', 1252) 1
    ('scaffold1', 1256) ('scaffold1', 1255) 1

