#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "hash.h"
#include "gff_reader.h"
#include "vcf_reader.h"

static int print_help()
{
    printf("Usage:\nLocateVariant -gff GFF -vcf VCF.\n");
    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
    char gff_file[FILE_NAME];
    char vcf_file[FILE_NAME];
    int nesscessary_parameters = 0;
    for (int arg_index = 1; arg_index < argc - 1; arg_index += 2)
    {
        if (!strcmp(argv[arg_index], "-gff"))
        {
            strncpy(gff_file, argv[arg_index + 1], FILE_NAME - 1);
            gff_file[FILE_NAME - 1] = 0;
            nesscessary_parameters++;
        }
        else if (!strcmp(argv[arg_index], "-vcf"))
        {
            strncpy(vcf_file, argv[arg_index + 1], FILE_NAME - 1);
            vcf_file[FILE_NAME - 1] = 0;
            nesscessary_parameters++;
        };
    };
    if (nesscessary_parameters != 2)
    {
        print_help();
    };
    unsigned long hash_size = 1572869;

    ChromosomeTranscript **chromosome_transcript = NULL;
    ChromosomeVariant **chromosome_variant = NULL;

    read_gff_file(gff_file, &chromosome_transcript, hash_size);
    read_vcf_file(vcf_file, &chromosome_variant, hash_size);

    char type[9];

    printf("Variant\tChromosome\tPosition\tTranscript\tType\tStart\tEnd\n");
    for (unsigned long hash_index = 0; hash_index < hash_size; hash_index++)
    {
        for (ChromosomeVariant *chromosome_variant_node = chromosome_variant[hash_index]; chromosome_variant_node; chromosome_variant_node = chromosome_variant_node->next)
        {
            Variant *Variant_Node = chromosome_variant_node->variant;
            char *chromosome = chromosome_variant_node->chromosome;
            ChromosomeTranscript *chromosome_transcript_node = find_chromosome_transcript(chromosome_transcript, hash_size, chromosome);
            if (chromosome_transcript_node)
            {
                Transcript *Transcript_Node = chromosome_transcript_node->transcript;
                unsigned long transcript_offset = 0;
                for (unsigned long variant_index = 0; variant_index < chromosome_variant_node->variant_number; variant_index++)
                {
                    bool flag = true; // Intergeneic region.
                    Variant *variant_node = Variant_Node + variant_index;
                    for (unsigned long transcript_index = 0; transcript_index < chromosome_transcript_node->transcript_number; transcript_index++)
                    {
                        Transcript *transcript_node = Transcript_Node + transcript_offset + transcript_index;
                        if (variant_node->position < transcript_node->start)
                        {
                            if (flag)
                            {
                                printf("%s\t%s\t%lu\t-\t-\t-\t-\n", variant_node->variant, chromosome, variant_node->position);
                            };
                            break;
                        }
                        else if (variant_node->position <= transcript_node->end) // Overlap.
                        {
                            // printf("%s\t%s\t%lu\t%s\t%lu-%lu\n", variant_node->variant, chromosome, variant_node->position, transcript_node->transcript, transcript_node->start, transcript_node->end);
                            Element *Element_Node = transcript_node->element;
                            for (unsigned long element_index = 0; element_index < transcript_node->exon_number * 2 + transcript_node->cds_number - 1; element_index++)
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
                                    default:
                                        break;
                                    };
                                    printf("%s\t%s\t%lu\t%s\t%s\t%lu\t%lu\n", variant_node->variant, chromosome, variant_node->position, transcript_node->transcript, type, element_node->positions[0], element_node->positions[1]);
                                };
                            };
                            flag = false;
                        }
                        else if (!transcript_index)
                        {
                            transcript_offset++;
                        };
                    };
                };
            }
            else
            {
                for (unsigned long variant_index = 0; variant_index < chromosome_variant_node->variant_number; variant_index++)
                {
                    printf("%s\t%s\t%lu\t-\t-\t-\t-\n", (Variant_Node + variant_index)->variant, chromosome, (Variant_Node + variant_index)->position);
                };
            };
        };
    };
    free_chromosome_transcript_hash(chromosome_transcript, hash_size);
    //free_chromosome_transcript_hash//
    return 0;
}