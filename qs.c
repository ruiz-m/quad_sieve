#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<stdint.h>
#include<string.h>

#include"lgz.h"
#include"heap.h"

#define MEM 4

typedef struct quad
{
	uint32_t pr;
	uint32_t sol[2];
}quad;

quad *quadros;
lgz **mem;
lgz_div_t *c;



void add(int8_t *r1, int8_t *r2, int32_t size)
{
	for(int32_t i=0; i<size; ++i)
	{
		r1[i] ^= r2[i];
	}
}

void qs_rref(int8_t *b, int8_t **A, int32_t row, int32_t col)
{
	if(col < row)
	{
		return;
	}

	//get ref	
	int32_t i = 0;
	int32_t j = 0;
	while(i < row && j < col)
	{
		while(j < col && A[i][j] == 0)
		{
			for(int32_t k=i+1; k<row; ++k)
			{
				if(A[k][j] == 1)
				{
					int8_t *tmp = A[i];
					A[i] = A[k];
					A[k] = tmp;
					goto yabo;
				}
			}
			++j;
		}

		yabo:
		if(j < col)
		{
			for(int32_t k=i+1; k<row; ++k)
			{
				if(A[k][j] == 1)
				{
					add(A[k], A[i], col);
				}
			}
		}
		++i;
		++j;
	}

	//get rref
	i=row-1;
	for(; i>=0; --i)
	{
		j=0;
		for(; j<col; ++j)
		{
			if(A[i][j] == 1)
			{
				for(int32_t k=0; k<i; ++k)
				{
					if(A[k][j] == 1)
					{
						add(A[k], A[i], col);
					}
				}
				break;
			}
		}
	}
	
	int32_t prev_j = 0;
	for(i=0; i<row; ++i)
	{
		for(j=0; j<col; ++j)
		{
			if(A[i][j] == 1)
			{
				for(int32_t k=prev_j; k<j; ++k)
				{
					b[k] = 1;
				}
				prev_j = j+1;
				for(int32_t k=j+1; k<col; ++k)
				{
					b[j] ^= A[i][k];
				}
				break;
			}
		}
	}
	for(int32_t k=prev_j; k<col; ++k)
	{
		b[k] = 1;
	}
}

uint32_t qs_powm(int64_t a, uint32_t p, uint32_t n)
{
	uint32_t res = 1;
	a = a%n;

	if(a == 0)
	{
		return 0;
	}
	if(p == 0)
	{
		return 1;
	}

	while(p > 0)
	{
		if(p&1)
		{
			res = (res*a)%n;
		}

		p >>= 1;
		a = (a*a)%n;
	}
	return res;
}

int32_t qs_legendre(lgz *a, uint32_t p)
{
	if(p == 2)
	{
		return a->num[0]&1;
	}
	
	lgz_memset(mem[0], 0);
	lgz_mod32(mem[0], a, p);
	int64_t amod = lgz_int64(mem[0]);
	//amod can fit in uint32_t type

	return qs_powm(amod, (p-1)/2, p); 
}

uint32_t find_nquad(uint32_t p)
{
	uint32_t yabo = (p-1)/2;
	for(uint32_t i=2; i<p; ++i)
	{
		if(qs_powm(i, yabo, p) == p-1)
		{
			return i;
		}
	}
}

//only works if |a| = 2^n in F_n 
uint32_t find_order2(uint64_t a, uint32_t n)
{
	uint32_t pow = 0;
	while(a != 1)
	{
		a = (a*a)%n;
		++pow;
	}
	
	return pow;
}

uint32_t qs_tonelli_shanks(uint32_t a, uint32_t p)
{
	uint32_t q = p-1;
	uint32_t s = 0;

	while((q&1) == 0)
	{
		q >>= 1;
		++s;
	}

	uint32_t z = find_nquad(p);
	uint32_t nres = qs_powm(z, q, p);
	uint64_t t = qs_powm(a, q, p);
	uint64_t r = qs_powm(a, (q+1)/2, p);
	uint32_t prev_order = s;

	while(t != 1)
	{
		uint32_t curr_order = find_order2(t, p);
		uint32_t b = qs_powm(nres, (1<<(prev_order-curr_order-1)), p);
		nres = ((uint64_t)b*b)%p;
		r *= b;
		t *= nres;
		r %= p;
		t %= p;

		prev_order = curr_order;
	}

	return r;
}

void qs_get_quad(uint64_t B, lgz *n)
{
	quadros = (quad*)malloc(B*sizeof(quad));
	uint32_t count = 0;
	uint32_t sub_count = 0x10000;
	
	uint32_t *pr = (uint32_t*)malloc(0x10000*sizeof(uint32_t));
	FILE *fp = fopen("./prime/prime32.txt", "r");
	
	while(count < B)
	{
		if(feof(fp))
		{
			break;
		}
		if(sub_count == 0x10000)
		{
			sub_count = 0;
			fread(pr, sizeof(uint32_t), 0x10000, fp);
		}
		
		if(qs_legendre(n, pr[sub_count]) == 1)
		{
			quadros[count].pr = pr[sub_count];
			++count;	
		}
		++sub_count;
	}
	
	free(pr);
	fclose(fp);
}

void qs_polyQ(lgz *out, int64_t x, lgz *sqrtn, lgz *n)
{
	lgz_enter(mem[2], x);
	lgz_add(mem[2], mem[2], sqrtn);
	lgz_mul(out, mem[2], mem[2]);
	lgz_sub(out, out, n);
	if(lgz_sign(out) == 255)
	{
		lgz_mul32(out, out, -1);
	}
}

int32_t I64log2(uint64_t x)
{
	int32_t res = -1;
	while(x != 0)
	{
		x >>= 1;
		++res;
	}
	return res;
}

float qs_ln(lgz *n)
{
	int mag = lgz_mag(n);
	uint8_t last = n->num[mag-1];

	int log2 = I64log2(last);
	log2 += 8*(mag-1);
	return (float)(log2*0.69315);
}

void qs_factor(lgz *n)
{
	//find B 
	float ln_n = qs_ln(n);
	float ex = exp(sqrt(ln_n*log(ln_n)));
	uint32_t B = ceil(pow(ex, sqrt(2)/4.0));
	uint64_t count = 0;

	lgz *a = new_lgz();
	lgz *b = new_lgz();
	lgz *sqrtn = new_lgz();
	printf("compute square root\n");
	lgz_sqrt(sqrtn, n);
	printf("sqrtn=");lgz_print(sqrtn);

	uint32_t min = lgz_mag(sqrtn)+1;
	uint32_t col = B+1;

	int8_t *scr = (int8_t*)calloc(col, sizeof(int8_t));
	uint32_t *queue = (uint32_t*)malloc(B*sizeof(uint32_t));
	int64_t *smooth_x = (int64_t*)malloc(col*sizeof(int64_t));
	int8_t **A = (int8_t**)malloc(B*sizeof(int8_t*));
	for(int64_t i=0; i<B; ++i)
	{
		A[i] = (int8_t*)calloc(col,sizeof(int8_t));
	}
	heap *h1 = new_heap(2*B);
	heap *h2 = new_heap(2*B);
	
	printf("B=%llu\n", B);
	printf("col=%u\n", col);
	printf("min=%u\n", min);

	qs_get_quad(B, n);
	
	getchar();
	for(int64_t i=0; i<B; ++i)
	{
		lgz_memset(mem[0], 0);
		lgz_mod32(mem[0], n, quadros[i].pr);
		int64_t n_modp = lgz_int64(mem[0]);
		
		lgz_memset(mem[0], 0);
		lgz_mod32(mem[0], sqrtn, quadros[i].pr);
		int64_t sqrtn_modp = lgz_int64(mem[0]);
		
		int32_t sol = qs_tonelli_shanks(n_modp, quadros[i].pr);
		int32_t sol1 = sol-sqrtn_modp;
		int32_t sol2 = quadros[i].pr-sol-sqrtn_modp;

		while(sol1 < 0)
		{
			sol1 += quadros[i].pr;
		}
		while(sol2 < 0)
		{
			sol2 += quadros[i].pr;
		}
		quadros[i].sol[0] = sol1;
		quadros[i].sol[1] = sol2;
	
		if(quadros[i].pr == 2)
		{
			hInsert(h1, sol1, i);
			hInsert(h2, sol1-2, i);
		}
		else
		{
			hInsert(h1, sol1, i);
			hInsert(h1, sol2, i);
			hInsert(h2, sol1-(int64_t)quadros[i].pr, i);
			hInsert(h2, sol2-(int64_t)quadros[i].pr, i);
		}
		//printf("%u, ", quadros[i].pr);
	}
	//printf("\n");
	/*hPrint(h1);
	hPrint(h2);
	getchar();*/

	while(count < col)
	{
		uint32_t fc = 0;
		int64_t x = peekX(h1);

		//hPrint(h1);
		while(peekX(h1) == x)
		{
			uint32_t qInd = peekQind(h1);
			queue[fc] = qInd;
			++fc;

			hDelete(h1);
			hInsert(h1, x+quadros[qInd].pr, qInd);
		}
		//printf("x=%lld\nfc=%u\n", x, fc);
		//getchar();
		if(fc >= min)
		{
			//printf("x=%lld\nfc=%u\n", x, fc);
			//getchar();
			memset(scr, 0, col);
			lgz_memset(mem[0], 0);
			qs_polyQ(mem[0], x, sqrtn, n);
			
			for(uint32_t j=0; j<fc; ++j)
			{
				uint32_t qInd = queue[j];
				uint32_t p = quadros[queue[j]].pr;

				while(!lgz_eqo(mem[0]))
				{
					lgz_memset(c->quot, 0);
					lgz_div32(c, mem[0], p);
					if(!lgz_eqz(c->rem))
					{
						break;
					}
					lgz_set(mem[0], c->quot);
					++scr[qInd];
					scr[qInd] &= 1;
				}
				if(lgz_eqo(mem[0]))
				{
					for(int32_t k=0; k<B; ++k)
					{
						A[k][count] = scr[k];
					}
					smooth_x[count] = x;
					++count;
					if(count%0x10 == 0)
					{
						printf("count=%u\n", count);
					}
					if(count == col)
					{
						goto end;
					}
				}
			}
		}

		fc = 0;
		x = peekX(h2);
		//hPrint(h2);
		while(peekX(h2) == x)
		{
			uint32_t qInd = peekQind(h2);
			queue[fc] = qInd;
			++fc;

			hDelete(h2);
			hInsert(h2, x-quadros[qInd].pr, qInd);
		}
		//printf("x=%lld\nfc=%u\n", x, fc);
		//getchar();
		if(fc >= min)
		{
			//printf("x=%lld\nfc=%u\n", x, fc);
			//getchar();
			memset(scr, 0, col);
			lgz_memset(mem[0], 0);
			qs_polyQ(mem[0], x, sqrtn, n);

			for(uint32_t j=0; j<fc; ++j)
			{
				uint32_t qInd = queue[j];
				uint32_t p = quadros[queue[j]].pr;

				while(!lgz_eqo(mem[0]))
				{
					lgz_memset(c->quot, 0);
					lgz_div32(c, mem[0], p);
					if(!lgz_eqz(c->rem))
					{
						break;
					}
					lgz_set(mem[0], c->quot);
					++scr[qInd];
					scr[qInd] &= 1;
				}
				if(lgz_eqo(mem[0]))
				{
					for(int32_t k=0; k<B; ++k)
					{
						A[k][count] = scr[k];
					}
					smooth_x[count] = x;
					++count;
					if(count%0x10 == 0)
					{
						printf("count=%u\n", count);
					}
				}
			}
		}
	}
	end:
	
	/*printf("smooth x\n");
	for(int i=0; i<col; ++i)
	{
		printf("%lld, ", smooth_x[i]);
	}
	printf("\n");*/

	/*for(int i=0; i<B; ++i)
	{
		for(int j=0; j<col; ++j)
		{
			printf("%d ", A[i][j]);
		}
		printf("\n");
	}
	printf("\n");*/
	
	memset(scr, 0, B);
	qs_rref(scr, A, B, col);

	/*for(int i=0; i<B; ++i)
	{
		for(int j=0; j<col; ++j)
		{
			printf("%d ", A[i][j]);
		}
		printf("\n");
	}
	printf("\n");*/
	/*for(int i=0; i<col; ++i)
	{
		printf("%d ", scr[i]);
	}
	printf("\n");*/
	
	lgz_enter(mem[3], 1);
	for(int i=0; i<col; ++i)
	{
		if(scr[i])
		{
			lgz_memset(mem[1], 0);
			lgz_memset(mem[0], 0);
			qs_polyQ(mem[1], smooth_x[i], sqrtn, n);
			lgz_mul(mem[0], mem[3], mem[1]);
			lgz_set(mem[3], mem[0]);
		}
	}
	lgz_sqrt(a, mem[0]);
	printf("\na=");
	lgz_print(a);

	lgz_enter(mem[3], 1);	
	for(int i=0; i<col; ++i)
	{
		if(scr[i])
		{
			lgz_memset(mem[0], 0);
			lgz_memset(mem[1], 0);
			lgz_add32(mem[0], sqrtn, smooth_x[i]);
			if(lgz_sign(mem[0]) == 255)
			{
				lgz_mul32(mem[0], mem[0], -1);
			}
			lgz_mul(mem[1], mem[3], mem[0]);
			lgz_set(mem[3], mem[1]);
		}
	}
	lgz_set(b, mem[3]);
	printf("b="); lgz_print(b);

	lgz_memset(mem[0], 0);
	lgz_add(mem[0], a, b);
	lgz_gcd(mem[1], mem[0], n);
	printf("gcd="); lgz_print(mem[1]);

	lgz_free(sqrtn);
	lgz_free(a);
	lgz_free(b);
	free(scr);
	hFree(h1);
	hFree(h2);
	free(smooth_x);
	for(int64_t i=0; i<B; ++i)
	{
		free(A[i]);
	}
	free(A);
}

void mem_alloc()
{
	mem = (lgz**)malloc(MEM*sizeof(lgz*));
	for(int32_t i=0; i<MEM; ++i)
	{
		mem[i] = new_lgz();
	}
	

	c = new_lgz_div_t();
}

void mem_free()
{
	for(int32_t i=0; i<MEM; ++i)
	{
		lgz_free(mem[i]);
	}
	lgz_div_t_free(c);

	free(mem);
	
	free(quadros);
}

int main(int argc, char **argv)
{
	if(argc != 2)
	{
		return 1;
	}

	lgz *num = new_lgz();
	int32_t u = 0;
	while(argv[1][u] != 0)
	{
		int8_t c = argv[1][u] - '0';
		lgz_mul32(num, num, 10);
		lgz_add32(num, num, c);
		++u;
	}
	lgz_print(num);
	
	mem_alloc();

	qs_factor(num);
	
	lgz_free(num);
	mem_free();

	return 0;
}
