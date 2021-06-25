unsigned long ElfHash(char *string)
{
    /*
    53ul, 97ul, 193ul, 389ul, 769ul, 1543ul, 3079ul, 6151ul, 12289ul
    24593ul, 49157ul, 98317ul, 196613ul, 393241ul, 786433ul, 1572869ul,
    3145739ul, 6291469ul, 12582917ul, 25165843ul, 50331653ul, 100663319ul,
    201326611ul, 402653189ul, 805306457ul, 1610612741ul, 3221225473ul, 4294967291ul
    */
    unsigned long hash = 0;
    unsigned long high;
    while (*string)
    {
        hash = (hash << 4) + *string++; // 将hash值左移4位后加上字符ascii
        high = hash & 0xF0000000;
        if (high)
        {
            hash ^= high >> 24;
        }
        hash &= ~high;
    }
    return hash;
}