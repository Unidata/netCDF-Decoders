#ifdef __cplusplus
extern "C" void user_makeparamtable(char *userfile);
extern "C" void user_printparamtable(void);
extern "C" char *user_pname(int param);
extern "C" char *user_plongname(int param);
extern "C" int  user_gribpcode(char *pname);
extern "C" char *user_gribunits(int param);
#elif defined(__STDC__)
extern void user_makeparamtable(char *userfile);
extern void user_printparamtable(void);
extern char *user_pname(int param);
extern char *user_plongname(int param);
extern int  user_gribpcode(char *pname);
extern char *user_gribunits(int param);
#else
extern void user_makeparamtable(/* char *userfile */);
extern void user_printparamtable(/* void */);
extern char *user_pname(/* int param */);
extern char *user_plongname(/* int param */);
extern int  user_gribpcode(/* char *pname */);
extern char *user_gribunits(/* int param */);
#endif
