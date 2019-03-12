# GTEx-data-download
## 1.进入GTEx Portal
* 国内网站进入后无内容，用谷歌浏览器可打开
* Database/Data Download 下“Annotations”下载注释；“RNA-Seq Data”下载表达量文件
## 2.下载RNAseq counts数据
* GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz
* 该文件每一行为基因名，列为样本名（所有组织的样本全在一个表格中）
## 3.下载注释文件
* GTEx_v7_Annotations_SampleAttributesDS.txt
* 由于需要找出某一组织的样本，故需要对表达量文件的列名进行注释，提取目标组织的样本表达量
* 注：注释文件内不仅包含RNAseq的注释，还包括WES等的注释信息
