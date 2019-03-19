

vector_t PhongShader(
	void *shader_data,
    hit_t *hit_ptr,
    render_params_t *rp_ptr,
    int samples, int max_samples,
    int iteration, int max_iterations);

vector_t NormalShader(
    void *shader_data,
    hit_t *hit_ptr,
    render_params_t *rp_ptr,
    int samples, int max_samples,
    int iteration, int max_iterations);

vector_t DepthShader(
    void *shader_data,
    hit_t *hit_ptr,
    render_params_t *rp_ptr,
    int samples, int max_samples,
    int iteration, int max_iterations);

vector_t PositionShader(
    void *shader_data,
    hit_t *hit_ptr,
    render_params_t *rp_ptr,
    int samples, int max_samples,
    int iteration, int max_iterations);
