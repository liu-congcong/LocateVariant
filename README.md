# LocateVariant

Locate variants in genes:

* Variant: variant id
* Chromosome: chromosome id
* Position: position (1-based)
* Transcript: transcript id | - (intergenic region)
* Type: exon | intron | promoter | cds | 5'utr | 3'utr
* Start: start position of the type (1-based)
* End: end position of the type (1-based)
* **Offset: offset will be added soon**

|Variant|Chromosome|Position|Transcript|Type|Start|End|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|tmp_1_3592_G_C|1|3592|transcript:AT1G01010.1|promoter|1631|3630|
|tmp_1_3651_C_G|1|3651|transcript:AT1G01010.1|5'utr|3631|3759|
|tmp_1_3651_C_G|1|3651|transcript:AT1G01010.1|exon|3631|3913|
|ENSVATH10438726|1|3755|transcript:AT1G01010.1|5'utr|3631|3759|
|ENSVATH04500126|1|3767|transcript:AT1G01010.1|cds|3760|3913|
|tmp_1_3779_T_A|1|3779|transcript:AT1G01010.1|exon|3631|3913|
|tmp_1_3779_T_A|1|3779|transcript:AT1G01010.1|cds|3760|3913|
|ENSVATH04500127|1|3944|transcript:AT1G01010.1|intron|3914|3995|
|tmp_1_3965_T_C|1|3965|transcript:AT1G01010.1|intron|3914|3995|
|ENSVATH04500212|1|11760|transcript:AT1G01030.2|3'utr|11649|11863|
|ENSVATH04500212|1|11760|transcript:AT1G01030.1|3'utr|11649|11863|

## Installation

```shell
git clone https://github.com/liu-congcong/LocateVariant.git
cd LocateVariant
gcc *.c -o LocateVariant
```

## Usage

```shell
LocateVariant -vcf VCF -gff GFF > OUTPUT
```

## Change logs

* 2021/06/19: Add 5'utr and 3'utr support.
