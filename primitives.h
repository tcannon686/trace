

typedef struct sphere
{
    vector_t origin;
    float radius;
} sphere_t;

void InitPrimitives(render_settings_t *rs);
void CommandsSetPrimitives(hashtable_t *table);
