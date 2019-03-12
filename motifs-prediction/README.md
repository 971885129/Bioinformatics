# motifs-prediction software

## Clover
* https://zlab.bu.edu/clover/
* 下载地址https://zlab.bu.edu/~mfrith/downloads/clover-2015-08-28.tar.gz
* 解压 & 安装

      gunzip clover-2015-08-28.tar.gz
      tar -xvf clover-2015-08-28.tar
      cd clover-2015-08-28/
      make
      ln -s /opt/Clover/clover-2015-08-28/clover /bin/clover
完成

## rsat系列软件安装
* 下载地址http://download.rsat.eu/
* 解压& 安装

      cd ../rsat/
      perl perl-scripts/configure_rsat.pl
      source RSAT_config.bashrc
      make -f makefiles/init_rsat.mk init
      make -f makefiles/install_rsat.mk perl_modules_list
      make -f makefiles/install_rsat.mk perl_modules_check
      make -f makefiles/install_rsat.mk perl_modules_install


## PoSSuM
* 下载地址https://bibiserv.cebitec.uni-bielefeld.de/resources/download/possumsearch/PoSSuM-2_0-linux-gnu.x86_64.tar.gz
* 解压 & 安装

      gunzip PoSSuM-2_0-linux-gnu.x86_64.tar.gz
      tar -xf PoSSuM-2_0-linux-gnu.x86_64.tar
      cd PoSSuM-2_0/
      将bin目录加入环境变量
完成
