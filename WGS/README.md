# WGS-data-analysis-script-journal
全基因组测序流程搭建日志
## 所使用软件
FastaQC
BWA
Samtools

## 分析所需文件准备
#### 注意：建议参考基因组使用hg38，GRCh38的到的结果中染色体编号没有chr，与其他库的文件无法匹配
* 比对参考基因组

        /media/sdb/user/wxm/database/Ensembl/release-94/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa
        
* samtools为参考基因组建立索引（GATK需要）

        nohup time samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa > faidx.o 2> faidx.e
* .dict文件

        gatk CreateSequenceDictionary -R Homo_sapiens.GRCh38.dna.primary_assembly.fa
* BWA 索引

        nohup time bwa index /media/sdb/user/wxm/database/Ensembl/release-94/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa >log.o 2>log.e &
* GATK要求参考基因组有.dict文件

*
*
*


## 测试数据下载

        新发现fastq-dump --split-files SRR1770413 这样使用可直接下载对应文件，并拆分为fastq文件

    #递归下载目标路径下的文件（韩国人基因组数据）
    nohup wget -r -nd ftp://ftp.kobic.re.kr/pub/KPGP/2015_release_candidate/WGS/KPGP-00001/ >log.o 2>log.e &
    #千人基因组原始测序数据下载地址
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/
    #千人基因组按人种划分数据（无原始数据）
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/

## 测序数据质控

        fastqc -o L1 -t 1 /media/sdb/user/wxm/test/WGS/data/KPGP-00001_L1_R1.fq.gz /media/sdb/user/wxm/test/WGS/data/KPGP-                      00001_L1_R2.fq.gz
        
## 数据预处理
![image](https://github.com/971885129/WGS-data-analysis-script-journal/blob/master/image/Data%20pre-processing%20for%20variant%20discovery.png)
* bwa的使用

        语法
        bwa mem -t 50 -M <ref> <reads.fq1> <reads.fq2>
        mem 适用于reads长度70~1M bp长度的序列
        -t 调用线程数
        -M Mark shorter split hits as secondary (for Picard compatibility)
        -R read group header line such as '@RG\tID:foo\tSM:bar'

bwa mem 的-R设置 read group信息，非常重要,后面BQSR步骤会用到。它是read数据的组别标识，其中的ID，PL和SM信息在正式的项目中是不能缺少的(如果样本包含多个测序文库的话，LB信息也不要省略)，另外由于考虑到与GATK的兼容关系，PL（测序平台）信息不能随意指定，必须是：ILLUMINA，SLX，SOLEXA，SOLID，454，LS454，COMPLETE，PACBIO，IONTORRENT，CAPILLARY，HELICOS或UNKNOWN这12个中的一个。


* samtools的使用

        samtools view
        -b: out bam(不设置该参数输出sam)
        -h: print header for the SAM output
        -F: filtering flag(数字4代表没有比对到参考序列)
        -q: minimum mapping quality
        -S: input is SAM(默认下是输入bam)
 
        samtools rmdup
        -s: single-end reads
        -S: Paired-end reads作为single-end reads处理
        
* GATK安装

        wget https://github.com/broadinstitute/gatk/releases/download/4.0.11.0/gatk-4.0.11.0.zip

* GATK bundle下载

        wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz
        wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
        wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/af-only-gnomad.hg38.vcf.gz
        wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/af-only-gnomad.hg38.vcf.gz.tbi
        wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/GetPileupSummaries/small_exac_common_3.hg38.vcf.gz
        wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/GetPileupSummaries/small_exac_common_3.hg38.vcf.gz.tbi
* 建立索引

        gatk IndexFeatureFile --feature-file dbsnp_146.hg38.vcf
        gatk IndexFeatureFile --feature-file Mills_and_1000G_gold_standard.indels.hg38.vcf
        gatk IndexFeatureFile --feature-file 1000G_phase1.snps.high_confidence.hg38.vcf
* GATK 要求两个文件来访问参考文件和对参考文件进行安全检查

        .fai文件
        samtools faidx genome.fa
        .dict文件
        gatk CreateSequenceDictionary -R genome.fa 

* GATK BaseRecalibrator工具

        time gatk BaseRecalibrator \
                -R /media/sda/database/UCSC/hg38/Sequence/BWAIndex/genome.fa \
                -I KPGP-00001_L1.sorted.rmdup.bam \
                --known-sites /media/sdb/user/wxm/database/GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf \
                --known-sites /media/sdb/user/wxm/database/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
                --known-sites /media/sdb/user/wxm/database/GATK_bundle/dbsnp_146.hg38.vcf \
                -O KPGP-00001_L1.sorted.rmdup.base_recal.table && echo "** KPGP-00001_L1.sorted.markdup.recal_data.table done **"        
* GATK ApplyBQSR工具

        time gatk ApplyBQSR \
                --bqsr-recal-file KPGP-00001_L1.sorted.rmdup.base_recal.table \
                 -R /media/sda/database/UCSC/hg38/Sequence/BWAIndex/genome.fa \
                -I KPGP-00001_L1.sorted.rmdup.bam \
                -O KPGP-00001_L1.sorted.rmdup.BQSR.bam && echo "** ApplyBQSR done **"
* 比对
        
        bwa mem -t 50 -M /media/sda/database/UCSC/hg38/Sequence/BWAIndex/genome.fa -R '@RG\tID:L1\tPL:illumina\tLB:library\tSM:L1' /media/sdb/user/wxm/test/WGS/data/KPGP-00001_L1_R1.fq.gz /media/sdb/user/wxm/test/WGS/data/KPGP-00001_L1_R2.fq.gz > KPGP-00001_L1.sam

* 比对后对数据进行合并（同一个样本不同lane的bam合并）
* 排序 & 去除multiple mapping、PCR duplication及低质量比对结果

       #排序
       samtools sort -@ 15 -O BAM -o KPGP-00001_L1.sorted.bam  KPGP-00001_L1.sam &&
       #融合merge（需增加-h参数）
       
       #转换为sam格式并过滤没有匹配上和低质量reads
       samtools view -h -F 4 -q 5 KPGP-00001_L1.sorted.bam |
       #转换为bam
       samtools view -bS >KPGP-00001_L1.tmp.bam
       #去除RNA duplication reads
       samtools rmdup KPGP-00001_L1.tmp.bam KPGP-00001_L1.sorted.rmdup.bam &&
       #建立索引
       samtools index KPGP-00001_L1.sorted.rmdup.bam
* 插入突变位点重新比对(新版GATK已经删除该流程)

* 碱基质量校正（base quality score recalibrate, BQSR）

        time gatk BaseRecalibrator \
                -R /media/sda/database/UCSC/hg38/Sequence/BWAIndex/genome.fa \
                -I KPGP-00001_L1.sorted.rmdup.bam \
                --known-sites /media/sdb/user/wxm/database/GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf \
                --known-sites /media/sdb/user/wxm/database/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
                --known-sites /media/sdb/user/wxm/database/GATK_bundle/dbsnp_146.hg38.vcf \
                -O KPGP-00001_L1.sorted.rmdup.base_recal.table && echo "** KPGP-00001_L1.sorted.markdup.recal_data.table done **"   
                
        time gatk ApplyBQSR \
                --bqsr-recal-file KPGP-00001_L1.sorted.rmdup.base_recal.table \
                 -R /media/sda/database/UCSC/hg38/Sequence/BWAIndex/genome.fa \
                -I KPGP-00001_L1.sorted.rmdup.bam \
                -O KPGP-00001_L1.sorted.rmdup.BQSR.bam && echo "** ApplyBQSR done **"
        #建立索引
        time samtools index KPGP-00001_L1.sorted.rmdup.BQSR.bam && echo "** ${sample}.sorted.markdup.BQSR.bam index done **"
* 变异检测前确定bam文件是否符合GATK要求，运行

        gatk ValidateSamFile -I KPGP-00001_L1.sorted.rmdup.BQSR.bam，如果显示 no error，则可以用HaplotypeCaller call SNP/Indel

## call SNP/InDel

### Germline Variation
使用GATK 的 HaplotypeCaller
![image](https://github.com/971885129/WGS-data-analysis-script-journal/blob/master/image/Germline%20Variation.png)
### Somatic Variation(肿瘤)
* 一个问题，肿瘤的variation中是否去除了Germline Variation

                已去除
* 测试数据下载地址（全外显子数据）

        ftp://gsapubftp-anonymous@ftp.broadinstitute.org/tutorials/datasets/tutorial_11136.tar.gz
![image](https://github.com/971885129/WGS-data-analysis-script-journal/blob/master/image/Somatic%20Variation.png)

#### 1.根据normal 样本得到 panel of normal
* 首先对每个normal 样本，运行Mutect2

        gatk Mutect2 \   
            -R /media/sda/database/UCSC/hg38/Sequence/BWAIndex/genome.fa \   
            -I normal1.bam \   
            -tumor normal1_sample_name \   
            --germline-resource /media/sdb/user/wxm/database/GATK_bundle/af-only-gnomad.hg38.vcf.gz \   
            -O normal1_for_pon.vcf.gz
* 然后使用CreateSomaticPanelOfNormals命令创建panel of normal

        gatk CreateSomaticPanelOfNormals \
            -vcfs normal1_for_pon_vcf.gz \   
            -vcfs normal2_for_pon_vcf.gz \   
            -vcfs normal3_for_pon_vcf.gz \   
            -O pon.vcf.gz

 This effectively omits common germline variant alleles from the PoN
<br>If you don't have a PoN, you may get more false positives that are sequencing artifacts. If possible, you can try to find publicly available data to use in your PoN. Have a look at this article for more information. Also note, having a matched normal and germline resource will help a lot with filtering out potential germline variants.

#### 2. normal_vs_turmor 得到体细胞突变

            gatk Mutect2 \
                -R /media/sda/database/UCSC/hg38/Sequence/BWAIndex/genome.fa \
                -I tumor.bam \
                -tumor tumor_sample_name \
                -I normal.bam \
                -normal normal_sample_name \
                --germline-resource /media/sdb/user/wxm/database/GATK_bundle/af-only-gnomad.hg38.vcf.gz \
                --af-of-alleles-not-in-resource 0.00003125 \
                --panel-of-normals pon.vcf.gz \
                -O somatic.vcf.gz

af-of-alleles-not-in-resource指定germline-resource 变异位点的频率，低于该频率的位点认为是一个不可靠的生殖细胞突变位点。
#### 3.过滤VCF文件
* 第一步，运行GetPileupSummaries

        gatk GetPileupSummaries \
            -I tumor.bam \
            -V /media/sdb/user/wxm/database/GATK_bundle/small_exac_common_3.hg38.vcf.gz \
            -L /media/sdb/user/wxm/database/GATK_bundle/small_exac_common_3.hg38.vcf.gz \
            -O pileups.table
* 第二步，运行CalculateContamination

        gatk  CalculateContamination \
            -I pileups.table \
            -O contamination.table
* 第三步，运行FilterMutectCalls

        gatk FilterMutectCalls \
            -V somatic.vcf.gz \
            -contamination-table contamination.table \
            -O filtere
## CNV
* Control-Freec 软件
* 该软件为预编译版本，由于编译所用glibc和gcc版本较高，而更新glibc会导致系统不稳定。故下载11.3版本，仅需更新gcc即可
* gcc更新

        #打印 libstdc++.so.6的输出信息中限定GLIBC库的信息
        strings /usr/lib64/libstdc++.so.6 | grep GLIBC
        #查看目前的libstdc++.so.6 是连接到哪一个？
        ll  /usr/lib64/libstdc++.so.6
        lrwxrwxrwx. 1 root root 19 4月  19 2017 /usr/lib64/libstdc++.so.6 -> libstdc++.so.6.0.19
        #尝试拷贝现有的
        cp /media/sdb/user/wxm/software/anaconda/install/lib/libstdc++.so.6.0.24 /usr/lib64/
        /bin/rm /usr/lib64/libstdc++.so.6
        ln -s /usr/lib64/libstdc++.so.6.0.24 /usr/lib64/libstdc++.so.6
        strings /usr/lib64/libstdc++.so.6 | grep GLIBC
* Control-Freec下载        
        
        wget https://github.com/BoevaLab/FREEC/archive/v11.3.tar.gz
        tar -xzf v11.3.tar.gz
        cd FREEC-11.5
        
* call CNV

        #配置导入文件
        #执行脚本
        freec -config config.txt


* 配置文件

## SV
### BreakDancer
* 1.软件下载（cpp版本）

        http://breakdancer.sourceforge.net/morecpp.html
* 2.安装

        cd breakdancer
        mkdir build
        cd build
        cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/biosoftware/breakdancer/build
        make
        make install
* 3.安装bam2cfg.pl所调用perl module

        安装module包括：Statistics::Descriptive（及其依赖module: List::MoreUtils、Exporter::Tiny、List::MoreUtils::XS）
        问题：安装module位置不是当前perl调用module路径，安装到了root下的路径
        解决： 1.更改安装路径（未实现）
              2.将安装的module拷贝到调用路径：/usr/lib64/perl5/

### Pindel

### Meerkat

### SvABA（Broad institute 2018 用于中、小SV的鉴定）

## LOH

## 融合基因

##
