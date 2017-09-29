#include "complex_scalar.h"

#ifndef VECTOR_H
#define VECTOR_H

#define VECTOR_INIT(vec, capacity) vector_init(&vec, capacity)
#define VECTOR_RESIZE(vec, capacity) vector_resize(&vec, capacity)
#define VECTOR_ADD(vec, item) vector_add(&vec, (void *) item)
#define VECTOR_ADD_DUPLICATE(vec, item) vector_add_duplicate(&vec, (void *) item)
#define VECTOR_SET(vec, id, item) vector_set(&vec, id, (void *) item)
#define VECTOR_GET(vec, type, id) (type) vector_get(&vec, id)
#define VECTOR_DELETE(vec, id) vector_delete(&vec, id)
#define VECTOR_CAPACITY(vec) vector_length(&vec)
#define VECTOR_TOTAL(vec) vector_total(&vec)
#define VECTOR_ITEMS(vec) vector_items(&vec)
#define VECTOR_FREE(vec) vector_free(&vec)
#define VECTOR_SHOW(vec) vector_show(&vec)
#define VECTOR_COPY(vec1, vec2) vector_copy(&vec1, &vec2)
#define VECTOR_DUPLICATE(vec1, vec2) vector_duplicate(&vec1, &vec2)

typedef struct vector {
	void ** items;
	int capacity;
	int total;
} vector;

void vector_init(vector * v, int capacity);
int vector_capacity(vector * v);
int  vector_total(vector * v);
void ** vector_items(vector * v);
static void vector_resize(vector * v, int capacity);
void vector_add(vector * v, void * item);
void vector_add_duplicate(vector * v, void * item);
void vector_set(vector * v, int index, void * item);
void * vector_get(vector * v, int index);
void vector_delete(vector * v, int index);
void vector_free(vector * v);
void vector_show(vector * v);
void vector_copy(vector * v, vector * newV);
void vector_duplicate(vector * v, vector * newV);

#endif
