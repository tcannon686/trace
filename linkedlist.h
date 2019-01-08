

struct list_element;
typedef struct list_element list_element_t;

typedef union
{
	void *data_ptr;
	long data_long;
	int data_int;
	short data_short;
	char data_char;
	double data_double;
	float data_float;
} list_data_t;

typedef struct list_element
{
	list_data_t data;
	list_element_t *next_ptr;
	list_element_t *last_ptr;
} list_element_t;

typedef struct
{
	list_element_t *base_ptr;
	list_element_t *top_ptr;
	int count;
	void (*free)(void *);
	int (*cmp)(void *, void *);
} list_t;

/*
 * Creates a new list. free is a function pointer that frees element data, or
 * NULL if the data should not be freed. cmp is a function for comparing data,
 * or NULL if not needed.
 */
list_t *ListNew(void (*free)(void *), int (*cmp)(void *, void *));

/*
 * Iterates through the list and frees each element and its data.
 */
void ListFree(list_t *list_ptr);

void ListAppend(list_t *list_ptr, list_data_t data);

/* Helper macro for appending to a list. */
#define ListAppendType(type, list_ptr, data) { list_data_t _d; _d.type = data; ListAppend(list_ptr, _d); }
#define ListAppendPointer(list_ptr, data) ListAppendType(data_ptr, list_ptr, data)
#define ListAppendLong(list_ptr, data) ListAppendType(data_long, list_ptr, data)
#define ListAppendInt(list_ptr, data) ListAppendType(data_int, list_ptr, data)
#define ListAppendShort(list_ptr, data) ListAppendType(data_short, list_ptr, data)
#define ListAppendChar(list_ptr, data) ListAppendType(data_char, list_ptr, data)
#define ListAppendDouble(list_ptr, data) ListAppendType(data_double, list_ptr, data)
#define ListAppendFloat(list_ptr, data) ListAppendType(data_float, list_ptr, data)
list_data_t ListGet(list_t *list_ptr, int index);
list_data_t ListGetFirst(list_t *list_ptr);
list_data_t ListGetLast(list_t *list_ptr);
/*
 * Removes an element, and frees its data if should_free is true.
 */
void ListRemove(list_t *list_ptr, int index, int should_free);
void ListRemoveFirst(list_t *list_ptr, int should_free);
void ListRemoveLast(list_t *list_ptr, int should_free);
int ListIndexOf(list_t *list_ptr, list_data_t data);

#define ListSize(list_ptr)				list_ptr->count

/*
 * Replaces a for loop to iterate through list_ptr, storing the data of each
 * element in data_ptr.
 */
#define ListIterate(list_ptr, data_ptr)			for(\
	list_element_t *_next_ptr, *_current_ptr = list_ptr->base_ptr;\
	_current_ptr != NULL ? (_next_ptr = _current_ptr->next_ptr,\
		*((list_data_t *)(data_ptr)) = _current_ptr->data, 1) : 0;\
	_current_ptr = _next_ptr)


