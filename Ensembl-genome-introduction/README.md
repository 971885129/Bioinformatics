# Ensembl-genome-introduction
### Ensemble基因组数据的形式包含以下2种：
* masked/unmasked
<br>dna_sm- Repeats soft-masked (converts repeat nucleotidesto lowercase)
<br>dna_rm- Repeats masked (converts repeats to to N's)
<br>dna- No masking
* toplevel / primary assembly
<br>toplevel- Includes haplotype information (notsure how aligners deal with this)
<br>primary_assembly– contains all toplevel sequenceregions excluding haplotypes and patches. This is best used for performingsequence similarity searches where patch and haplotype sequences would confuseanalysis.
* 根据README中的介绍，primary_assembly 和 toplevel相比不包含haplotype, 更适合用于比对，对于mask/un mask 通常选择softmask或者unmasked, 一般不用rm的。
