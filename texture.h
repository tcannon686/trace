

vector_t Texture2dSampleNearestNeighbor(texture2d_t *tex_ptr, vector_t co);
vector_t Texture2dSampleBilinear(texture2d_t *tex_ptr, vector_t co);
vector_t Texture2dSample(texture2d_t *tex_ptr, vector_t co);
void CommandsSetTexture(hashtable_t *table_cmds);
void Texture2dFree(void *ptr);
void ImageFree(void *ptr);

