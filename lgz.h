#ifndef __LGZ_H__
#define __LGZ_H__
#include<stdint.h>

typedef struct lgz
{
	uint8_t *num;
	int length;
}lgz;

typedef struct lgz_div
{
	lgz *quot;
	lgz *rem;
}lgz_div_t;

lgz *new_lgz();
lgz_div_t *new_lgz_div_t();
lgz *lgz_make_copy(lgz *num);
void lgz_set(lgz *b, lgz *a);
void lgz_resize(lgz *num, int length);
void lgz_free(lgz *num);
void lgz_div_t_free(lgz_div_t *num);
void lgz_print(lgz *num);
void lgz_memset(lgz *num, uint8_t bit);
void lgz_enter(lgz *num, uint64_t input);
int64_t lgz_int64(lgz *num);
void lgz_enter_str(lgz *num, char *str, int length);
void lgz_enter_num(lgz *num, char *str);
int lgz_eqz(lgz *a);
int lgz_eqo(lgz *a);

uint8_t lgz_sign(lgz *num);
int lgz_mag(lgz *num);
int lgz_leq(lgz *a, lgz *b);
lgz *lgz_min(lgz *a, lgz *b);
lgz *lgz_max(lgz *a, lgz *b);
void lgz_lshift(lgz *a, int num);
void lgz_rshift(lgz *a, int num);
void lgz_add(lgz *c, lgz *a, lgz *b);
void lgz_sub(lgz *c, lgz *a, lgz *b);
void lgz_mul(lgz *c, lgz *a, lgz *b);
void lgz_div(lgz_div_t *c, lgz *a, lgz *b);
void lgz_add32(lgz *c, lgz *a, int32_t num);
void lgz_mul32(lgz *c, lgz *a, int32_t num);
void lgz_div32(lgz_div_t *c, lgz *a, int32_t num);
void lgz_mod32(lgz *c, lgz *a, int32_t num);
void lgz_mod(lgz *c, lgz *a, lgz *b);
void lgz_powm(lgz *c, lgz *a, uint64_t p, lgz *n);
void lgz_sqrt(lgz *b, lgz *a);
void lgz_gcd(lgz *c, lgz *a, lgz *b);

#endif //__LGZ_H__
