#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "hash.h"
#include "vcf_reader.h"

static int create_chromosome_variant_hash(ChromosomeVariant ***chromosome_variant_hash, unsigned long hash_size)
{
    *chromosome_variant_hash = malloc(sizeof(ChromosomeVariant *) * hash_size);
    memset(*chromosome_variant_hash, 0, sizeof(ChromosomeVariant *) * hash_size);
    return 0;
}

static int add2chromosome_variant_hash(ChromosomeVariant **hash, unsigned long hash_size, char *chromosome, unsigned long position, char *variant)
{
    bool flag = false;
    unsigned long hash_value = ElfHash(chromosome) % hash_size;
    ChromosomeVariant *chromosome_variant_node = hash[hash_value];
    while (chromosome_variant_node)
    {
        if (!strcmp(chromosome_variant_node->chromosome, chromosome))
        {
            flag = true;
            break;
        }
        else if (chromosome_variant_node->next)
        {
            chromosome_variant_node = chromosome_variant_node->next;
        }
        else
        {
            break;
        };
    };
    if (!flag)
    {
        ChromosomeVariant *new_chromosome_variant_node = malloc(sizeof(ChromosomeVariant));
        new_chromosome_variant_node->chromosome = malloc(sizeof(char) * (strlen(chromosome) + 1));
        strcpy(new_chromosome_variant_node->chromosome, chromosome);
        new_chromosome_variant_node->variant_number = 0;
        new_chromosome_variant_node->variant = NULL;
        new_chromosome_variant_node->temp_variant = NULL;
        new_chromosome_variant_node->next = NULL;
        if (chromosome_variant_node)
        {
            chromosome_variant_node->next = new_chromosome_variant_node;
        }
        else
        {
            hash[hash_value] = new_chromosome_variant_node;
        };
        chromosome_variant_node = new_chromosome_variant_node;
    };
    chromosome_variant_node->variant_number++;

    Variant *new_variant_node = malloc(sizeof(Variant));
    new_variant_node->variant = malloc(sizeof(char) * (strlen(variant) + 1));
    strcpy(new_variant_node->variant, variant);
    new_variant_node->position = position;
    new_variant_node->next = NULL;
    if (chromosome_variant_node->temp_variant)
    {
        chromosome_variant_node->temp_variant->next = new_variant_node;
    }
    else
    {
        chromosome_variant_node->variant = new_variant_node;
    };
    chromosome_variant_node->temp_variant = new_variant_node;
    return 0;
}

static int compare(const void *x1, const void *x2)
{
    return ((Variant *)x1)->position < ((Variant *)x2)->position ? -1 : 1;
}

static int sort_chromosome_variant_hash(ChromosomeVariant **hash, unsigned long hash_size)
{
    for (unsigned long hash_index = 0; hash_index < hash_size; hash_index++)
    {
        for (ChromosomeVariant *chromosome_variant_node = hash[hash_index]; chromosome_variant_node; chromosome_variant_node = chromosome_variant_node->next)
        {
            Variant *new_variant_node = malloc(sizeof(Variant) * chromosome_variant_node->variant_number);
            unsigned long variant_index = 0;
            while (chromosome_variant_node->variant)
            {
                Variant *variant_node = chromosome_variant_node->variant->next;
                memcpy(new_variant_node + variant_index, chromosome_variant_node->variant, sizeof(Variant));
                (new_variant_node + variant_index++)->next = NULL;
                free(chromosome_variant_node->variant);
                chromosome_variant_node->variant = variant_node;
            };
            qsort(new_variant_node, chromosome_variant_node->variant_number, sizeof(Variant), compare);
            chromosome_variant_node->variant = new_variant_node;
        };
    };
    return 0;
}

int read_vcf_file(char *file, ChromosomeVariant ***chromosome_variant_hash, unsigned long hash_size)
{
    char buffer[LINE];
    char chromosome[LINE]; // 0
    unsigned long position; // 1
    char variant[LINE]; // 2

    create_chromosome_variant_hash(chromosome_variant_hash, hash_size);

    FILE *open_file = fopen(file, "r");
    while (fgets(buffer, LINE, open_file))
    {
        bool new_line = true;
        if ((buffer[0] != '#') && new_line)
        {
            const char *format = "%s\t%lu\t%[^\t]";
            sscanf(buffer, format, chromosome, &position, variant);
            add2chromosome_variant_hash(*chromosome_variant_hash, hash_size, chromosome, position, variant);
        };
        new_line = buffer[strlen(buffer) - 1] == '\n' ? true : false;
    };
    fclose(open_file);
    sort_chromosome_variant_hash(*chromosome_variant_hash, hash_size);
    return 0;
}