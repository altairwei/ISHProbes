找找评估 RNA 二级结构，dimer 的库，写一个脚本自动搜索最好的探针。
用 Biopython 读取 fasta 文件，评估与排序所有 k-mer 的二级结构评分。
如果有可能的话，生成两到三个探针，并评估相互之间的 dimer 结构。
有没有公司能直接合成 DIG-labeled 200bp 探针？不行就合成带有 T7 启动子和终止子的 DNA （单链或双链？）


## 思路1：OligoMiner 设计小麦全转录组探针库

将小麦全转录组序列构建成库（基因序列），然后 OligoMiner 设计探针库。

问题是如何处理可变剪切？~~利用 GFF 文件，从 FASTA 中将序列一个个提取出来。~~ 

- 不能用 gene ，因为不包含 UTR 区域，还包含 intron。
- 不能用 mRNA，因为重复的区域太多了，除非只保留一个。
- 不能用 exon 因为它太短了。
- CDS + three_prime_UTR 也不行，CDS 也是重复的。
- 只有整个基因组了。

先把所有的 mRNA 提取出来，选择长度最长的。

使用人员在探针库中找到合适的探针，人工检查它们是否位于想要的区域内。

安装：beliveau-lab/OligoMiner 与 altairwei/gfftools

指定基因组区域，提取存在的探针，以 FASTA 文件输出。

在选定 probes 后，可以用 ViennaRNA Package 模拟一下探针的热动力学特性。


