#include "vector.h"

int vector_capacity(vector * v){
	return v->capacity;
}

int vector_total(vector * v){
	return v->total;
}

void ** vector_items(vector * v){
	return v->items;
}

void vector_init(vector * v, int capacity){
	v->capacity = capacity;
	v->total = 0;
	v->items = malloc(sizeof(void *) * v->capacity);
}

static void vector_resize(vector * v, int capacity){
	#ifdef DEBUG_ON
	printf("vector_resize: %d to %d\n", v->capacity, capacity);
	#endif

	v->items = realloc(v->items, sizeof(void *) * capacity);
	if(v->items){
		v->capacity = capacity;
	}
}

void vector_add(vector * v, void * item){
	if(v->capacity == v->total){
		vector_resize(v, v->capacity + 1);
	}
	v->items[v->total++] = item;
}

void vector_add_duplicate(vector * v, void * item){
	complex * item_dup = (complex *)malloc(sizeof(complex));
	if(v->capacity == v->total){
		vector_resize(v, v->capacity + 1);
	}
	complex_copy(item_dup, (complex *)item);
	v->items[v->total++] = (void *)item_dup;
}

void vector_set(vector * v, int index, void * item){
	if(index >= 0 && index < v->total){
		complex_copy((complex *)v->items[index], (complex *)item);
	}
}

void * vector_get(vector * v, int index){
	if(index >= 0 && index < v->total){
		return v->items[index];
	}
	return NULL;
}

void vector_delete(vector * v, int index){
	if(index >= 0 && index < v->total){
		for(int i = index; i < v->total - 1; i++){
			v->items[i] = v->items[i + 1];
		}
		v->total--;
		if(v->total > 0){
			vector_resize(v, v->capacity - 1);
		}else if(v->total == 0){
			vector_free(v);
		}
	}else{
		return;
	}
}

void vector_free(vector * v){
	v->capacity = 0;
	v->total = 0;
	free(v->items);
	v->items = NULL;
}

void vector_show(vector * v){
	for(int i=0; i<vector_total(v); i++){
		complex_show((complex *)vector_get(v, i));
	}
	printf("\n");
}

void vector_copy(vector * v, vector * newV){
	vector_init(newV, vector_capacity(v));

	for(int i=0; i<v->total; i++){
		vector_add(newV, vector_get(v, i));
	}
}

void vector_duplicate(vector * v, vector * newV){
	vector_init(newV, vector_capacity(v));

	for(int i=0; i<v->total; i++){
		vector_add_duplicate(newV, vector_get(v, i));
	}
}