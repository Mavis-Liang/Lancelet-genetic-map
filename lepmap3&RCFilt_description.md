

[toc]

### 0. 预备知识
①两点测交（Two-point Analysis）：在减数分裂中期，同源染色体联会，会发生交换recombination（如下图c）。距离越远的两个位点，重组的频率越高。因此，如果能计算出交换的次数，则能对应到一个相对的遗传距离。（当交换偶数次时，也不能表现出重组型）
<img src="https://img-blog.csdnimg.cn/20210312172328865.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80Nzc3MjEyMw==,size_16,color_FFFFFF,t_70" alt="在这里插入图片描述" style="zoom:33%;" />

②重组率r（recombination fraction）：每两个marker间都有一个r，由发生重组的个体比上所有个体得到。r=0，没有发生重组，两个marker完全连锁；r>0.5，相当于两个marker自由组合，即两个marker完全独立，完全不连锁（在不同的染色体上）；r>0且r<0.5时，两个marker不完全连锁，这个时候，r越大，它们之间的距离越远，可以根据几种函数关系估计距离。r在雌性和雄性中可能不同，在软件中可以分别设置参数。特点：样本越大、r越小，估计越准确。

③LOD score：r的检验统计量LR服从卡方分布。而LOD是LR的变型。r越大，LOD越小，连锁关系越小。一般规定LOD大于3时，两个marker间存在连锁关系（即在同一Linkage Group中，后续步骤就是在LG内的排序）。实际情况下，也会把临界LOD设为4、5、6等小于10的值。

计算公式：
$$
L(r|n)=\frac{n!}{n_{NR}!n_{R}}(1-r)^{NR}r^R
$$

这是一个二项式，$n_{NR}$代表没有发生重组的个体，$n_R$代表发生了重组的个体。为了检验$H_0$: r=0.5，$H_1$: r<0.5，设卡方统计量LR（log-likelihood ratio)：

$$
LR=-2log[\frac{L(r=0.5|n)}{L(\hat{r}|n)}]
$$

这里的$\hat{r}$指①式的最大似然估计值（MLE）

$$
LOD=0.217LR
$$

LOD和r的关系：
<img src="https://img-blog.csdnimg.cn/20210314223557691.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80Nzc3MjEyMw==,size_16,color_FFFFFF,t_70" alt="在这里插入图片描述" style="zoom:75%;" />

LOD和$n_{R}$的关系：
<img src="https://img-blog.csdnimg.cn/20210314223557691.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80Nzc3MjEyMw==,size_16,color_FFFFFF,t_70" alt="在这里插入图片描述" style="zoom:75%;" />

### 1.Lep-map3下载及安装
用`mkdir`新建一个文件夹，然后到官网下载lep-map3到这个目录中：[lep-map3下载地址](https://sourceforge.net/projects/lep-map3/files/)。用`unzip`解压，得到四个文件，`README`，`COPYING`，`src`，和`bin`。

直接在该目录下输入`java -cp bin/ ParentCall2`可以调用相应模块（ParentCall2），则说明安装成功。

### 2.上软件前的准备
制作一个pedigree.txt的新文件。第一行是自己定的家族名字，第二行是每个个体的名字，第三、四行是第二行个体对应的母亲和父亲的名字，第5行是各个个体的性别（1是雄性，2是雌性，0为未知），第六行全写为0即可。第一第二列是“CHR”和“POS”。各列之间用tab隔开。（PS。我用的家系中，父本缺失，这种情况下，应该在pedigree中补全家族信息：在第二行加个“male”，在第三行F1对应写上“male”。原本的vcf不变。）
![在这里插入图片描述](https://img-blog.csdnimg.cn/20210311134618273.PNG?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80Nzc3MjEyMw==,size_16,color_FFFFFF,t_70#pic_center =500x)
可以通过`java -cp bin/ ParentCall2 data=../pedigree.txt`检验这六行有没有问题（譬如，error504可能是某行列数不正确、error503可能是没用tab空开等等）



### 3.ParentCall2
这一步是根据父本母本的基因型，汇总所有个体的likelihood信息，以便于后续LOD值的计算。
`zcat ../filename.vcf.gz | java -cp bin/ ParentCall2 data=../pedigree.txt vcfFile=- |gzip >data.call.gz`

该命令中，先解压vcf，使其通过标准流流入模块中，data参数是之前的pedigree.txt文件，removeNonInformative=1指去掉vcf中的注释行（不算#开头那一行个体名称），最后把输出文件压缩。

执行后输出下列结果。
<img src="https://img-blog.csdnimg.cn/20210311180143403.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80Nzc3MjEyMw==,size_16,color_FFFFFF,t_70" alt="在这里插入图片描述" style="zoom:33%;" />



得到的结果文件：
![在这里插入图片描述](https://img-blog.csdnimg.cn/20210311174759529.PNG?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80Nzc3MjEyMw==,size_16,color_FFFFFF,t_70#pic_center =500x)


### 4.SeparateChromosomes2
这一步通过计算两两snp之间的LOD值，然后与我们自定义的LOD值比对，来判断这两个位点是否在同一个连锁群当中。若大于设定的LOD，则在同一连锁群。


若是single family，则不需要使用filtering模块。直接调用SeparateChromosome2。值得注意的是，lep-map3文章作者提供了另一套流程（SeparateIdenticals,  JoinSingles2Identicals,  SeparateChromosomes2），其中还涉及到数据模拟。但新版Lep-map3强调可以不用这么麻烦。

``java -cp bin/ SeparateChromosomes2 data=data.call lodLimit=5 >map5.nofilt.txt
``

结果如下：
![<img src="https://img-blog.csdnimg.cn/20210311175757728.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80Nzc3MjEyMw==,size_16,color_FFFFFF,t_70" alt="(https://img-blog.csdnimg.cn/20210311175644438.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80Nzc3MjEyMw==,size_16,color_FFFFFF,t_70)" style="zoom:30%;" />

得到一个txt文件，里面有一列数字，每一行对应到原来的snp，数字代表分到的LG，0代表离群值。
<img src="https://img-blog.csdnimg.cn/20210311192933996.png" alt="在这里插入图片描述" style="zoom:45%;" />

### 5.分window重复上述步骤

①用脚本fileSplit.ipynb，按scaffold不同将文件分割成多个1000个SNP，不足1000SNP的scaffold也为一个文件。得到2043个文件。

②逐个通过软件：
``
#!/bin/bash
for i in `ls ../../newLanc/bioFilt_split`;do
java -cp bin/ ParentCall2 data=../../newLanc/pedigree.txt vcfFile=../../newLanc/bioFilt_split/$i >./data12.call/data12_$i.call;
java -cp bin/ SeparateChromosomes2 data=./data12.call/data12_$i.call lodLimit=12 >./result12/$i.txt
done
``

``
for i in `ls ../../newLanc/bioFilt_split`;do
java -cp bin/ ParentCall2 data=../../newLanc/pedigree.txt vcfFile=../../newLanc/bioFilt_split/$i >./data13.call/data13_$i.call;
java -cp bin/ SeparateChromosomes2 data=./data13.call/data13_$i.call lodLimit=13 >./result13/$i.txt
done
``

``
for i in `ls ../../newLanc/bioFilt_split`;do
java -cp bin/ ParentCall2 data=../../newLanc/pedigree.txt vcfFile=../../newLanc/bioFilt_split/$i >./data14.call/data14_$i.call;
java -cp bin/ SeparateChromosomes2 data=./data14.call/data14_$i.call lodLimit=14 >./result14/$i.txt
done
echo "The job is done!"
``

③用recomFilt_v1.ipynb进行统计和软删除。



|        | >80%线 | >85%线 |
| :----: | :----: | :----: |
| LOD=12 | 56.9%  | 53.5%  |
| LOD=13 | 51.4%  | 46.9%  |
| LOD=14 | 41.8%  | 33.6%  |

*表格内百分比为将要保留的整段1000snp占所有的1000snp的比例
