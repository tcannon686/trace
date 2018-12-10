

#ifdef INCLUDE_GUI
int CmdRenderWindow(render_settings_t *rs);
#endif

int CmdQuit(render_settings_t *rs);
int CmdVertex(render_settings_t *rs);
int CmdMakeFace(render_settings_t *rs);
int CmdNormal(render_settings_t *rs);

int CmdMatIndex(render_settings_t *rs);
int CmdMatSetNumber(render_settings_t *rs);
int CmdMatSetInteger(render_settings_t *rs);
int CmdMatSetVector(render_settings_t *rs);
int CmdMatShader(render_settings_t *rs);

int CmdLightPosition(render_settings_t *rs);
int CmdLightSetVector(render_settings_t *rs);
int CmdLightSetNumber(render_settings_t *rs);
int CmdLightSetInteger(render_settings_t *rs);

int CmdMakeMaterial(render_settings_t *rs);
int CmdMakeLight(render_settings_t *rs);

int CmdInFile(render_settings_t *rs);
int CmdOutFile(render_settings_t *rs);
int CmdInStdIn(render_settings_t *rs);

int CmdSky(render_settings_t *rs);

int CmdCamFov(render_settings_t *rs);

int CmdRenderSamples(render_settings_t *rs);
int CmdRenderShadowSamples(render_settings_t *rs);
int CmdRenderSectionSize(render_settings_t *rs);
int CmdRenderThreads(render_settings_t *rs);
int CmdRenderIterations(render_settings_t *rs);
int CmdRender(render_settings_t *rs);

int CmdPrintTriangles(render_settings_t *rs);

int CmdTransformPop(render_settings_t *rs);
int CmdTransformPush(render_settings_t *rs);
int CmdTransformTranslate(render_settings_t *rs);
int CmdTransformScale(render_settings_t *rs);
int CmdTransformRotate(render_settings_t *rs);


void CommandsSetStandard(hashtable_t *table_cmds);


