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

**Locate variant:**

|Variant|Chromosome|Position|Transcript|Type|Start|End|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|tmp_1_3592_G_C|1|3592|transcript:AT1G01010.1|promoter|1631|3630|
|tmp_1_3651_C_G|1|3651|transcript:AT1G01010.1|5'utr|3631|3759|
|tmp_1_3651_C_G|1|3651|transcript:AT1G01010.1|exon|3631|3913|
|...|...|...|...|...|...|...|
|tmp_1_3779_T_A|1|3779|transcript:AT1G01010.1|cds|3760|3913|
|ENSVATH04500127|1|3944|transcript:AT1G01010.1|intron|3914|3995|
|...|...|...|...|...|...|...|
|ENSVATH04500212|1|11760|transcript:AT1G01030.2|3'utr|11649|11863|
|ENSVATH04500212|1|11760|transcript:AT1G01030.1|3'utr|11649|11863|

**Print gene structure:**

|Transcript|Strand|Type|Chromosome|Start|End|
|:-:|:-:|:-:|:-:|:-:|:-:|
|transcript:AT1G01010.1|+|promoter|1|1631|3630|
|transcript:AT1G01010.1|+|5'utr|1|3631|3759|
|transcript:AT1G01010.1|+|exon|1|3631|3913|
|transcript:AT1G01010.1|+|cds|1|3760|3913|
|transcript:AT1G01010.1|+|intron|1|3914|3995|
|transcript:AT1G01010.1|+|cds|1|3996|4276|
|transcript:AT1G01010.1|+|exon|1|3996|4276|
|transcript:AT1G01010.1|+|intron|1|4277|4485|
|transcript:AT1G01010.1|+|exon|1|4486|4605|
|transcript:AT1G01010.1|+|cds|1|4486|4605|
|transcript:AT1G01010.1|+|intron|1|4606|4705|
|transcript:AT1G01010.1|+|exon|1|4706|5095|
|transcript:AT1G01010.1|+|cds|1|4706|5095|
|transcript:AT1G01010.1|+|intron|1|5096|5173|
|transcript:AT1G01010.1|+|exon|1|5174|5326|
|transcript:AT1G01010.1|+|cds|1|5174|5326|
|transcript:AT1G01010.1|+|intron|1|5327|5438|
|transcript:AT1G01010.1|+|cds|1|5439|5630|
|transcript:AT1G01010.1|+|exon|1|5439|5899|
|transcript:AT1G01010.1|+|3'utr|1|5631|5899|

## Installation

```shell
git clone https://github.com/liu-congcong/LocateVariant.git
cd LocateVariant
gcc *.c -o LocateVariant
```

## Usage

```shell
# Locate variant. #
LocateVariant -vcf VCF -gff GFF > OUTPUT

# Print gene structure. #
LocateVariant -gff GFF -print_gene_structure > OUTPUT
```

## Change logs

* 2021/06/19: Add 5'utr and 3'utr support.
* 2021/06/25: Add a function to print gene structure.
