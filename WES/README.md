# WES-data-analysis

#### CNV
* 测试用外显子捕获区文件下载

      参考博客：http://www.bioinfo-scrounger.com/archives/642
      exon.list（GATK-style）：http://www.bioinfo-scrounger.com/data/exon.list
* 转换为control-FREEC可识别格式

      编辑exon.list文件，改为可导入Control-FRECC的格式（替换“：” “-”为“\t”）
      python脚本对其去重、排序（位置E:\2019.1.18_WES_CNV）
* 构建config文件

      [general]
      chrLenFile = /media/sdb/user/wxm/database/UCSC/hg38/WholeGenomeFasta/genome.fa.fai
      window = 0
      #step = 10000
      ploidy = 2
      chrFiles = /media/sdb/user/wxm/database/UCSC/hg38/Chromosomes
      minMappabilityPerWindow = 0.7
      sex = XY
      breakPointType = 4
      maxThreads = 30
      sambamba = /media/sdb/user/wxm/software/sambamba/sambamba
      SambambaThreads = 30 

      [sample]
      mateFile = /media/sdb/user/wxm/test/WES/1.Alignment/TCRBOA6-T-WEX.sorted.rmdup.BQSR.bam
      inputFormat = BAM
      mateOrientation = FR

      [control]
      mateFile = /media/sdb/user/wxm/test/WES/1.Alignment/TCRBOA6-N-WEX.sorted.rmdup.BQSR.bam
      inputFormat = BAM
      mateOrientation = FR

      [target]
      captureRegions = /media/sdb/user/wxm/database/CCDS/exon_nodup_sort

* 运行Control-FRECC

      /biosoftware/FREEC-11.3/src/freec -conf /media/sdb/user/wxm/test/WES/3.CNV/config_CNV.txt

