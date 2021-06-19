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

<img width="686" alt="image" src="https://user-images.githubusercontent.com/34596618/122646816-bb4fec00-d153-11eb-91a5-da4b333c4942.png">

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
