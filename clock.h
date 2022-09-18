#ifdef _WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif

void start_chrono();

double stop_chrono();