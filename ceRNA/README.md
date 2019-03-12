# ceRNA-analysis
## miRNA 靶mRNA lncRNA预测

## ceRNA 网络构建
1.差异mRNA lncRNA与数据库比对，得到mRNA-miRNA pairs lncRNA-miRNA pairs
<br>2.具有共同MRE（miRNA response element）的mRNA lnRNA为mRNA-lncRNA pairs
<br>3.筛选具有超过3个共同MRE的mRNA-lncRNA pairs
<br>4.筛选mRNA lncRNA 相关性cor>0.9 padj<0.05的mRNA-lncRNA pairs
<br>5.超几何检验筛选 FDR<0.01的 mRNA-lncRNA pairs
<br>6.最终获得mRNA-lncRNA pairs

## 下游分析
1. 共表达网络，寻找module，进行功能富集分析
2. Topological 筛选 degree Betweenness closeness
    * Betweenness for each node is defined as the number of the shortest paths that pass through the node
    * Closeness of a node is a measure of centrality in a network, calculated as the sum of the length of the shortest
paths between the node and all other nodes in the network graph.         
3. dysregulation ceRNA

## 详细技术
### 超几何检验
使用公式(R代码)

    1-phyper(x-1,K,N-K,M,lower.tail = TRUE, log.p = FALSE)
    x: lncRNA-mRNA pairs common miRNA number
    K,M: miRNA number of interact whith lncRNA/mRNA in lncRNA-mRNA pairs（分别对每一对基因进行计算）
    N: total miRNA number in this analysis(not K+M)
示例

## ceRNA筛选条件
lncRNA mRNA正相关；lncRNA mRNA共同结合miRNA超过3个； cor>0.9 padj<0.05; 超几何检验FDR<0.01; miRNA lncRNA负相关； miRNA mRNA负相关
## 功能注释
ceRNA网络内所有mRNA GO KEGG 富集分析
<br>hubgene分析:筛选degree前10%为hubgene

## 问题
预测结果如何进行超几何检验
