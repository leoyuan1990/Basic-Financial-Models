#ifndef _MERSENNETWISTER_H_
#define _MERSENNETWISTER_H_

#include <cstdlib>
#include <cmath>

void			init_genrand(unsigned long s);
void			init_by_array(unsigned long *init_key, unsigned long key_length);
unsigned long	genrand_int32(void);
long			genrand_int31(void);
double			genrand_real1(void);
double			genrand_real2(void);
double			genrand_real3(void);
double			genrand_res53(void);
double			gasdev();

#endif