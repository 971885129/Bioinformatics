# microArray-data-analysis

## affy illumina agilent 三种芯片间的区别
* 芯片的单通道 多通道

* 芯片取log2的原因
raw gene expression levels measured from microarray fluorescence intensity typically have a skewed log-normal distribution resulting from a multiplicative error during the amplification process
The log transformation allows to normalize the data distribution and use classical parametric statistics such as the t-test for analysis.

## affy芯片处理

## limma处理芯片
* 标准化

      neqc算法
      COMbat去除批次效应


## oligo处理芯片

