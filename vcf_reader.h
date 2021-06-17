#ifndef __VCF_READER_H__

#define __VCF_READER_H__

#define LINE 10240

#define FILE_NAME 10240

typedef struct Variant
{
    char *variant;
    unsigned long position;
    struct Variant *next;
} Variant;

typedef struct ChromosomeVariant
{
    char *chromosome;
    unsigned long variant_number;
    Variant *variant;
    Variant *temp_variant;
    struct ChromosomeVariant *next;
} ChromosomeVariant;

int read_vcf_file(char *, ChromosomeVariant ***, size_t);

#endif