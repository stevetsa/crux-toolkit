/*************************************************************************//**
 * \file hash.c
 * AUTHOR: David Crawshaw, Chris Park
 * $Revision: 1.15 $
 * \brief: Object for hashing.
 ****************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hash.h"
//#include "carp.h"
#include "parse_arguments.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"

// TODO why does hash include parameter and not the other way around?

// Table is sized by primes to minimise clustering.
static const unsigned int sizes[] = {
    53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
    196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
    50331653, 100663319, 201326611, 402653189, 805306457, 1610612741
};

static const unsigned int sizes_count = sizeof(sizes) / sizeof(sizes[0]);
static const FLOAT_T load_factor = 0.65;

/**
 * \struct record
 * \brief record for each value/key pair
 */
struct record {
  unsigned int hash; ///< hash algorithm code
  char* key; ///< the key for the record
  void* value; ///< value the record contains
  int count; ///< count for duplicate adds(if adding two of the same key count 1 is added)
};

/**
 * \struct hash
 * \brief hash table, contains the records
 */
struct hash {
  RECORD_T* records; ///< record holds key & values
  unsigned int records_count; ///< number of records
  unsigned int size_index; ///< index into the size array, thus can get the size of the hash table
};

/**
 *\struct hash_iterator
 *\brief An object that iterates over the keys in a hash 
 */
struct hash_iterator {
  HASH_T* hash; ///< the hash to iterate over
  int hash_idx;   ///< current hash key to return
  int hash_total; ///< total number of hash keys
};

// Function definition, description found below 
BOOLEAN_T add_hash_when_grow(
  HASH_T* h, ///< Hash object to add -in/out
  char *key, ///< key of the record to add -in
  void *value, ///< value to add to be hashed if needed -in
  int count  ///< the count of the record to grow -in
  );

/**
 * Increase the size of hash table
 *\returns the TRUE if successfully inceased hash table size, else FALSE
 */
static int hash_grow(
  HASH_T* h ///< Hash object to grow -out
  )
{
  unsigned int i;
  RECORD_T* old_recs;
  unsigned int old_recs_length;
  
  old_recs_length = sizes[h->size_index];
  old_recs = h->records;
  
  if (h->size_index == sizes_count - 1){
    return FALSE;
  }
  // increase larger hash table
  if ((h->records = mycalloc(sizes[++h->size_index],
                             sizeof(RECORD_T))) == NULL) {
    h->records = old_recs;
    return FALSE;
  }
  
  h->records_count = 0;
  
  // rehash table
  for (i=0; i < old_recs_length; i++){
    if (old_recs[i].hash && old_recs[i].key){
      add_hash_when_grow(h, old_recs[i].key, old_recs[i].value, old_recs[i].count);
    }
  }
  
  free(old_recs);
  
  return TRUE;
}

/**
 * hash algorithm 
 * \returns the array slot
 */
static unsigned int strhash(
  const char *str ///< string into the hash function -in
  )
{
  int c;
  int hash = 5381;
  while ((c = *str++)){
    hash = hash * 33 + c;
  }
  return hash == 0 ? 1 : hash;
}

/**
 * Create new hashtable with capacity.
 *\returns the hash table
 */
HASH_T* new_hash(
  unsigned int capacity ///< The capacity of the new hash table
  ) 
{
  HASH_T* h;
  unsigned int i, sind = 0;
  
  capacity /= load_factor;
  
  for (i=0; i < sizes_count; i++) 
    if (sizes[i] > capacity) { sind = i; break; }
  
  if ((h = malloc(sizeof(HASH_T))) == NULL) return NULL;
  if ((h->records = mycalloc(sizes[sind], sizeof(RECORD_T))) == NULL) {
    free(h);
    return NULL;
  }
  
  h->records_count = 0;
  h->size_index = sind;
  
  return h;
}

/**
 * free hashtable
 */
void free_hash(
  HASH_T* h ///< Hash object to free -in
  )
{
  unsigned int idx = 0;
  unsigned int size = sizes[h->size_index];
  // free up all key & values strings
  for(; idx < size; ++idx){
    if(h->records[idx].key != NULL){
      free(h->records[idx].key);
    }
    if(h->records[idx].value != NULL){
      free(h->records[idx].value);
    }    
  }
  
  free(h->records);
  free(h);
}

/**
 * add key and value to hash table.
 * If key exists, free current value and allocate and set new one
 * If key not found, allocate key, and allocate and set value
 * Does not copy value (for use with void pointers).
 *\returns TRUE if successfully adds to new record, else FALSE
 */
BOOLEAN_T add_or_update_hash(
  HASH_T* h, ///< Hash object to add to -in/out
  const char *key, ///< key of the record to add or update -in
  void *value ///< value to associate with the key -in
  )
{
  RECORD_T* recs;  
  int rc;
  unsigned int off, ind, size, code;

  if (key == NULL || *key == '\0') return FALSE;
  
  code = strhash(key);
  recs = h->records;
  size = sizes[h->size_index];
  
  ind = code % size;
  off = 0;
  
  // probe down until reaching open slot
  // Quadratic probing used
  while (recs[ind].key){     
    // if find duplicate key, thus identical item
    if ((code == recs[ind].hash) && recs[ind].key &&
        strcmp(key, recs[ind].key) == 0){
      // free existing value
       free(recs[ind].value); 
      // set new value
      recs[ind].value = value;
      return TRUE;
    }
    else{
      // continue to search
      ind = (code + (int)pow(++off,2)) % size;
    }
  }

  // key not found, add it
  // first check size
  if (h->records_count > sizes[h->size_index] * load_factor) {
    rc = hash_grow(h);
    if (rc) return FALSE;
  }


  recs[ind].hash = code;
  recs[ind].key = my_copy_string(key);
  recs[ind].value = value;
  recs[ind].count = 1;
  
  h->records_count++;
  
  return TRUE;
}



/**
 * add key and value to hash table.
 * If key exists, free current value and allocate and set new one
 * If key not found, allocate key, and allocate and set value
 * Copies the value
 *\returns TRUE if successfully adds to new record, else FALSE
 */
BOOLEAN_T add_or_update_hash_copy(
  HASH_T* h, ///< Hash object to add to -in/out
  const char *key, ///< key of the record to add or update -in
  void *value ///< value to associate with the key -in
  )
{
  char* new_value = my_copy_string(value);
  return add_or_update_hash(h, key, new_value);
}


/**
 * add key and value to hash table.
 * Must add a heap allocated key, value may be NULL
 * If finds duplicate key, just increase count by 1
 *\returns TRUE if successfully adds to new record, else FALSE
 */
BOOLEAN_T add_hash(
  HASH_T* h, ///< Hash object to add -in/out
  char *key, ///< key of the record to add -in
  void *value ///< value to add to be hashed if needed -in
  )
{
    RECORD_T* recs;
    int rc;
    unsigned int off, ind, size, code;

    if (key == NULL || *key == '\0') return FALSE;
    if (h->records_count > sizes[h->size_index] * load_factor) {
        rc = hash_grow(h);
        if (rc) return FALSE;
    }

    code = strhash(key);
    recs = h->records;
    size = sizes[h->size_index];

    ind = code % size;
    off = 0;

    // probe down until reach open slot
    // Quadratic probing used
    while (recs[ind].key){     
      // if find duplicate key, thus identical item
      if ((code == recs[ind].hash) && recs[ind].key &&
          strcmp(key, recs[ind].key) == 0){
        // increment count
        ++recs[ind].count;        
        free(key);
        free(value);
        return TRUE;
      }
      else{
        // continue to search
        ind = (code + (int)pow(++off,2)) % size;
      }
    }
    
    recs[ind].hash = code;
    recs[ind].key = key;
    recs[ind].value = value;
    recs[ind].count = 1;

    h->records_count++;
    
    return TRUE;
}

/**
 * add key and value to hash table.
 *\returns TRUE if successfully adds to new record, else FALSE
 */
BOOLEAN_T add_hash_when_grow(
  HASH_T* h, ///< Hash object to add -in/out
  char *key, ///< key of the record to add -in
  void *value, ///< value to add to be hashed if needed -in
  int count  ///< the count of the record to grow -in
  )
{
    RECORD_T* recs;
    int rc;
    unsigned int off, ind, size, code;

    if (key == NULL || *key == '\0') return FALSE;
    if (h->records_count > sizes[h->size_index] * load_factor) {
        rc = hash_grow(h);
        if (rc) return FALSE;
    }

    code = strhash(key);
    recs = h->records;
    size = sizes[h->size_index];

    ind = code % size;
    off = 0;

    // probe down until reach open slot
    while (recs[ind].key){
      ind = (code + (int)pow(++off,2)) % size;
    }
    
    recs[ind].hash = code;
    recs[ind].key = key;
    recs[ind].value = value;
    recs[ind].count = count;
    
    h->records_count++;
    
    return TRUE;
}

/**
 * Updates the value for the key
 * Must already have a existing value for the key
 * Copies the value, thus no need to pass in a heap allocated value
 *\returns TRUE if successfully updates hash value, else FALSE
 */
BOOLEAN_T update_hash_value(
  HASH_T* h, ///< Hash object to add -in/out
  const char *key, ///< key of the record to update -in
  void *value ///< value to add to be hash -in
  )
{
  RECORD_T* recs;  
  unsigned int off, ind, size, code;
  
  if (key == NULL || *key == '\0') return FALSE;
  
  code = strhash(key);
  recs = h->records;
  size = sizes[h->size_index];
  
  ind = code % size;
  off = 0;
  
  // probe down until reach open slot
  // Quadratic probing used
  while (recs[ind].key){     
    // if find duplicate key, thus identical item
    if ((code == recs[ind].hash) && recs[ind].key &&
        strcmp(key, recs[ind].key) == 0){
      // free existing value
      free(recs[ind].value); 
      // set new value
      recs[ind].value = my_copy_string(value); 
      return TRUE;
    }
    else{
      // continue to search
      ind = (code + (int)pow(++off,2)) % size;
    }
  }
  //carp(CARP_ERROR, "Failed to find key %s in hash table", key);
  return FALSE;
}

/**
 * Get the value of the record for the hash key
 *\return the value for the hash record, returns NULL if can't find key
 */
void* get_hash_value(
  HASH_T* h, ///< working hash object -in
  const char *key  ///< the key of the record to retrieve -in
  )
{
  RECORD_T* recs;
  unsigned int off, ind, size;
  unsigned int code = strhash(key);
  
  recs = h->records;
  size = sizes[h->size_index];
  ind = code % size;
  off = 0;
  
  // search on hash which remains even if a record has been removed,
  // so remove_hash() does not need to move any collision records
  while (recs[ind].hash) {
    if ((code == recs[ind].hash) && recs[ind].key &&
        strcmp(key, recs[ind].key) == 0)
      return recs[ind].value;
    ind = (code + (int)pow(++off,2)) % size;
  }
  
  return NULL;
}

/**
 * Get a pointer to the variable pointing to the value of the record 
 * for the hash key
 * BAD!! This is not the right thing to do.  It is so parse_arguments can
 * directly insert values into the hash without using the key.  Soon I should
 * change parse_arguments so that it takes a reference to the hash and uses
 * the option name to insert the value into the hash
 *\return a pointer to variable pointing to the value for the hash record, returns NULL if can't find key
 */
void** get_hash_value_ref(
  HASH_T* h, ///< working hash object -in
  char *key  ///< the key of the record to retrieve -in
  )
{
  RECORD_T* recs;
  unsigned int off, ind, size;
  unsigned int code = strhash(key);
  
  recs = h->records;
  size = sizes[h->size_index];
  ind = code % size;
  off = 0;
  
  // search on hash which remains even if a record has been removed,
  // so remove_hash() does not need to move any collision records
  while (recs[ind].hash) {
    if ((code == recs[ind].hash) && recs[ind].key &&
        strcmp(key, recs[ind].key) == 0)
      return &recs[ind].value;
    ind = (code + (int)pow(++off,2)) % size;
  }
  
  return NULL;
}

/**
 * Get the count of the record for the hash key
 *\return the count for the hash record, returns -1 if can't find key
 */
int get_hash_count(
  HASH_T* h, ///< working hash object -in
  char *key  ///< the key of the record to retrieve -in
  )
{
  RECORD_T* recs;
  unsigned int off, ind, size;
  unsigned int code = strhash(key);
  
  recs = h->records;
  size = sizes[h->size_index];
  ind = code % size;
  off = 0;
  
  // search on hash which remains even if a record has been removed,
  // so remove_hash() does not need to move any collision records
  while (recs[ind].hash) {
    if ((code == recs[ind].hash) && recs[ind].key &&
        strcmp(key, recs[ind].key) == 0)
      return recs[ind].count;
    ind = (code + (int)pow(++off,2)) % size;
  }
  
  return -1;
}

/**
 * Remove key from table, returning value.
 *\returns the value of removing record, returns NULL if can't find key
 */
void* remove_hash(
  HASH_T* h, ///< working hash object -in
  char *key  ///< the key of the record to remove -in
  )
{
  unsigned int code = strhash(key);
  RECORD_T* recs;
  void * value;
  unsigned int off, ind, size;
  
  recs = h->records;
  size = sizes[h->size_index];
  ind = code % size;
  off = 0;
  
  while (recs[ind].hash) {
    if ((code == recs[ind].hash) && recs[ind].key &&
        strcmp(key, recs[ind].key) == 0) {
      // do not erase hash, so probes for collisions succeed
      value = recs[ind].value;
      // free key
      free(recs[ind].key);
      recs[ind].key = 0;
      recs[ind].value = 0;
      h->records_count--;
      return value;
    }
    ind = (code + (int)pow(++off, 2)) % size;
  }
  
  return NULL;
}

/**
 *\returns total number of keys in the hashtable.
 */  
unsigned int hash_size(
  HASH_T* h ///< working hash object -in
  )
{
  return h->records_count;
}


/**
 * hash_iterator routines!
 */

/**
 *\returns a new memory allocated hash iterator
 */
HASH_ITERATOR_T* new_hash_iterator(
  HASH_T* hash ///< the hash collection to iterate -out
  ){
  if (hash == NULL){
    carp(CARP_FATAL, "Null hash collection passed to hash iterator");
  }
  
  // allocate a new hash iterator
  HASH_ITERATOR_T* hash_iterator = 
    (HASH_ITERATOR_T*)mycalloc(1, sizeof(HASH_ITERATOR_T));
  
  // set items
  hash_iterator->hash = hash;
  hash_iterator->hash_idx = 0;
  hash_iterator->hash_total = sizes[hash->size_index];

  return hash_iterator;
}

/**
 * Does the hash_iterator have another hash object to return?
 * \returns TRUE, if hash iterator has a next hash, else FALSE
 */
BOOLEAN_T hash_iterator_has_next(
  HASH_ITERATOR_T* hash_iterator ///< the working  hash iterator -in
  )
{
  HASH_T* hash = hash_iterator->hash;
  while (hash_iterator->hash_idx < hash_iterator->hash_total && 
         hash->records[hash_iterator->hash_idx].key == NULL){
    hash_iterator->hash_idx++;
  }
  return (hash_iterator->hash_idx < hash_iterator->hash_total);
}

/**
 * \returns the next the hash struct
 */
char* hash_iterator_next(
  HASH_ITERATOR_T* hash_iterator ///< the working hash iterator -in
  )
{
  HASH_T* hash = hash_iterator->hash;
  return hash->records[hash_iterator->hash_idx++].key;
}

/**
 * free the memory allocated iterator
 */
void free_hash_iterator(
  HASH_ITERATOR_T* hash_iterator ///< the hash iterator to free
  )
{
  if (hash_iterator != NULL){
    free(hash_iterator);
  }
}



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */