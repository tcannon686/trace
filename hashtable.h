
typedef struct hashtable_entry hashtable_entry_t;

typedef union
{
    void *pointer;
    int integer;
    vecc_t real;
    vector_t vector;
} hashtable_value_t;

typedef struct hashtable_entry {
    char key[32];
    hashtable_value_t value;
    hashtable_entry_t *next_ptr;
} hashtable_entry_t;

typedef struct {
    hashtable_entry_t *entries;
    int capacity;
    int count;
    double max_load_factor;
} hashtable_t;

int HashTableHash(char *key);
hashtable_t *HashTableNew(int capacity, double max_load_factor);
hashtable_t *HashTableNewDefault();
hashtable_entry_t *HashTableGet(hashtable_t *table, char *key);
hashtable_entry_t *HashTableGetOrInsert(hashtable_t *table, char *key);
void HashTableSet(hashtable_entry_t *entry, hashtable_value_t value);
hashtable_t *HashTableCopy(hashtable_t *table, int new_capacity);
void HashTableFree(hashtable_t *table);
void HashTableRehash(hashtable_t *table, int new_capacity);
