#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include "hash.h"
#include "gff_reader.h"
#include "vcf_reader.h"

static int print_help(char *program)
{
    printf("Usage:\n");
    printf("%s -gff GFF -vcf VCF.\n", program);
    printf("%s -gff GFF -print_gene_structure.\n", program);
    exit(EXIT_SUCCESS);
}

static int locate_variant(char *gff, char *vcf)
{
    unsigned long hash_size = 1572869ul;

    ChromosomeTranscript **chromosome_transcript = NULL;
    ChromosomeVariant **chromosome_variant = NULL;

    read_gff_file(gff, &chromosome_transcript, hash_size);
    read_vcf_file(vcf, &chromosome_variant, hash_size);

    char type[9];
    printf("Variant\tChromosome\tPosition\tTranscript\tType\tStart\tEnd\n");
    for (unsigned long hash_index = 0; hash_index < hash_size; hash_index++)
    {
        for (ChromosomeVariant *chromosome_variant_node = chromosome_variant[hash_index]; chromosome_variant_node; chromosome_variant_node = chromosome_variant_node->next)
        {
            Variant *variant_node = chromosome_variant_node->variant; // Head node of variants on a chromosome.
            char *chromosome = chromosome_variant_node->chromosome;
            ChromosomeTranscript *chromosome_transcript_node = find_chromosome_transcript(chromosome_transcript, hash_size, chromosome);
            if (chromosome_transcript_node)
            {
                Transcript *Transcript_Node = chromosome_transcript_node->transcript; // Head node of transcripts on a chromosome.
                unsigned long transcript_number = chromosome_transcript_node->transcript_number;

                /* For each variant. */
                for (unsigned long variant_index = 0; variant_index < chromosome_variant_node->variant_number; variant_index++)
                {
                    bool flag = true; // Intergeneic region.
                    unsigned long transcript_index = 0;
                    while (transcript_index < transcript_number)
                    {
                        Transcript *transcript_node = Transcript_Node + transcript_index;
                        if (variant_node->position < transcript_node->start)
                        {
                            if (flag)
                            {
                                printf("%s\t%s\t%lu\t-\t-\t-\t-\n", variant_node->variant, chromosome, variant_node->position);
                            }
                            break;
                        }
                        else if (variant_node->position <= transcript_node->end) // Overlap.
                        {
                            Element *Element_Node = transcript_node->element;
                            for (unsigned long element_index = 0; element_index < transcript_node->element_number; element_index++)
                            {
                                Element *element_node = Element_Node + element_index;
                                if (variant_node->position >= element_node->positions[0] && variant_node->position <= element_node->positions[1])
                                {
                                    switch (element_node->type)
                                    {
                                    case 'e':
                                        strcpy(type, "exon");
                                        break;
                                    case 'i':
                                        strcpy(type, "intron");
                                        break;
                                    case 'c':
                                        strcpy(type, "cds");
                                        break;
                                    case 'p':
                                        strcpy(type, "promoter");
                                        break;
                                    case '5':
                                        strcpy(type, "5'utr");
                                        break;
                                    case '3':
                                        strcpy(type, "3'utr");
                                        break;
                                    default:
                                        break;
                                    }
                                    printf("%s\t%s\t%lu\t%s\t%s\t%lu\t%lu\n", variant_node->variant, chromosome, variant_node->position, transcript_node->transcript, type, element_node->positions[0], element_node->positions[1]);
                                }
                            }
                            transcript_index++;
                            flag = false;
                        }
                        else if (transcript_index)
                        {
                            transcript_index++;
                        }
                        else
                        {
                            Transcript_Node++;
                            transcript_number--;
                        }
                    }
                    variant_node++;
                }
            }
            else
            {
                for (unsigned long variant_index = 0; variant_index < chromosome_variant_node->variant_number; variant_index++)
                {
                    printf("%s\t%s\t%lu\t-\t-\t-\t-\n", variant_node->variant, chromosome, variant_node->position);
                    variant_node++;
                }
            }
        }
    }
    free_chromosome_transcript_hash(chromosome_transcript, hash_size);
    return 0;
}

static int print_gene_structure(char *gff)
{
    unsigned long hash_size = 1572869ul;

    ChromosomeTranscript **chromosome_transcript = NULL;

    read_gff_file(gff, &chromosome_transcript, hash_size);

    char type[9];

    printf("Transcript\tStrand\tType\tChromosome\tStart\tEnd\n");
    for (unsigned long hash_index = 0; hash_index < hash_size; hash_index++)
    {
        for (ChromosomeTranscript *chromosome_transcript_node = chromosome_transcript[hash_index]; chromosome_transcript_node; chromosome_transcript_node = chromosome_transcript_node->next)
        {
            // chromosome: chromosome_transcript_node->chromosome
            for (unsigned long transcript_index = 0; transcript_index < chromosome_transcript_node->transcript_number; transcript_index++)
            {
                Transcript *transcript_node = chromosome_transcript_node->transcript + transcript_index;
                for (unsigned long element_index = 0; element_index < transcript_node->element_number; element_index++)
                {
                    Element *element = transcript_node->element + element_index;
                    switch (element->type)
                    {
                    case 'p':
                        strcpy(type, "promoter");
                        break;
                    case '5':
                        strcpy(type, "5'utr");
                        break;
                    case '3':
                        strcpy(type, "3'utr");
                        break;
                    case 'e':
                        strcpy(type, "exon");
                        break;
                    case 'c':
                        strcpy(type, "cds");
                        break;
                    case 'i':
                        strcpy(type, "intron");
                        break;
                    default:
                        break;
                    }
                    printf("%s\t%c\t%s\t%s\t%lu\t%lu\n", transcript_node->transcript, transcript_node->strand, type, chromosome_transcript_node->chromosome, element->positions[0], element->positions[1]);
                }
            }
        }
    }
    free_chromosome_transcript_hash(chromosome_transcript, hash_size);
    return 0;
}

int main(int argc, char *argv[])
{
    char gff[FILE_NAME];
    char vcf[FILE_NAME];

    bool gff_flag = false;
    bool vcf_flag = false;
    bool gene_structure_flag = false;

    for (int arg_index = 1; arg_index < argc; arg_index++)
    {
        if (!strcmp(argv[arg_index], "-gff"))
        {
            assert(arg_index + 1 < argc);
            strncpy(gff, argv[arg_index + 1], FILE_NAME);
            gff[FILE_NAME - 1] = 0;
            gff_flag = true;
        }
        else if (!strcmp(argv[arg_index], "-vcf"))
        {
            assert(arg_index + 1 < argc);
            strncpy(vcf, argv[arg_index + 1], FILE_NAME);
            vcf[FILE_NAME - 1] = 0;
            vcf_flag = true;
        }
        else if (!strcmp(argv[arg_index], "-print_gene_structure"))
        {
            gene_structure_flag = true;
        }
    }

    if (gff_flag && vcf_flag)
    {
        locate_variant(gff, vcf);
    }
    else if (gff_flag && gene_structure_flag)
    {
        print_gene_structure(gff);
    }
    else
    {
        print_help(argv[0]);
    }
    return 0;
}