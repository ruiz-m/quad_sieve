#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include"lgz.h"

#define DEFAULT 8

lgz *new_lgz()
{
	lgz *num = (lgz*)malloc(sizeof(lgz));
	num->num = (uint8_t*)calloc(DEFAULT, 1);
	num->length = DEFAULT;
	return num;
}

lgz_div_t *new_lgz_div_t()
{
	lgz_div_t *tup = (lgz_div_t*)malloc(sizeof(lgz_div_t));
	tup->quot = new_lgz();
	tup->rem = new_lgz();
	return tup;
}

//not in lgz.h
lgz *new_lgz_size(int size)
{
	lgz *num = (lgz*)malloc(sizeof(lgz));
	
	int actual_size = size/4;
	if(size % 4 != 0)
	{
		++actual_size;
	}
	if(actual_size == 1)
	{
		actual_size = 2;
	}

	num->num = (uint8_t*)calloc(actual_size, 4);
	num->length = 4*actual_size;
	return num;
}

//not in lgz.h
void lgz_size_clear(lgz *num, int size)
{
	int actual_size = size/4;
	if(size % 4 != 0)
	{
		++actual_size;
	}
	if(actual_size == 1)
	{
		actual_size = 2;
	}

	uint8_t *cnum = (uint8_t*)calloc(actual_size,4);
	uint8_t *temp = num->num;
	num->num = cnum;
	num->length = 4*actual_size;

	free(temp);
}

lgz *lgz_make_copy(lgz *num)
{
	lgz *copy = (lgz*)malloc(sizeof(lgz));
	copy->num = (uint8_t*)malloc(num->length);
	copy->length = num->length;
	
	for(int i=0;i<num->length;++i)
	{
		copy->num[i] = num->num[i];
	}
	return copy;
}

void lgz_set(lgz *b, lgz *a)
{
	if(b->length < a->length)
	{
		uint8_t *cnum = (uint8_t*)malloc(a->length);
		
		uint8_t *temp = b->num;
		b->num = cnum;
		b->length = a->length;

		free(temp);
	}
	
	lgz_memset(b,lgz_sign(a));

	for(int i=0;i<a->length;++i)
	{
		b->num[i] = a->num[i];
	}
}

void lgz_resize(lgz *num, int size)
{
	if(size <= num->length)
	{
		return;
	}

	int orig_length = num->length;
	uint8_t last = num->num[num->length-1];
	
	int actual_size = size/4;
	if(size % 4 != 0)
	{
		++actual_size;
	}
	if(actual_size == 1)
	{
		actual_size = 2;
	}

	num->length = 4*actual_size;

	uint8_t *cnum = (uint8_t*)malloc(num->length);
	uint8_t *temp = num->num;
	num->num = cnum;
	
	if(last >= 128)
	{
		lgz_memset(num, 255);
	}
	else
	{
		lgz_memset(num, 0);
	}
	
	for(int i=0;i<orig_length;++i)
	{
		num->num[i] = temp[i];
	}
	
	free(temp);
}

void lgz_free(lgz *num)
{
	if(num != NULL)
	{
		free(num->num);
		free(num);
	}
}

void lgz_div_t_free(lgz_div_t *num)
{
	if(num != NULL)
	{
		lgz_free(num->quot);
		lgz_free(num->rem);
		free(num);
	}
}

void lgz_print(lgz *num)
{	
	int length = num->length;
	for(int i=0;i<length;++i)
	{
		if(num->num[length-i-1] < 16)
		{
			printf("0");
		}
		printf("%x", num->num[length-i-1]);
	}
	printf("\n");
}

void lgz_memset(lgz *num, uint8_t bit)
{
	memset(num->num, bit, num->length);
}

int lgz_enter(lgz *num, int64_t input)
{
	if(input < 0)
	{
		lgz_memset(num, 255);
	}
	else
	{
		lgz_memset(num, 0);
	}

	int i = 0;
	while(input != 0)
	{
		num->num[i] = input;
		input = input >> 8;
		++i;
	}
	return 0;
}

int64_t lgz_long(lgz *num)
{
	int64_t res = 0;
	for(int i=num->length-1;i>=0;--i)
	{
		res = (res << 8)+num->num[i];
	}
	return res;
}

int lgz_enter_str(lgz *num, char *str, int length)
{
	if(length > num->length)
	{
		lgz_size_clear(num, length+2);
	}
	else
	{
		lgz_memset(num, 0);
	}

	for(int i=0;i<length;++i)
	{
		num->num[i] = (uint8_t)str[length-i-1];
	}
}

int lgz_enter_num(lgz *num, char *str)
{
	int length = strlen(str);
	if(length > num->length)
	{
		lgz_size_clear(num, ((length*9)/20));
	}
	else
	{
		lgz_memset(num, 0);
	}
	
	for(int i=0;i<length-1;++i)
	{
		lgz_add32(num, num, (str[i]-'0'));
		lgz_mul32(num, num, 10);
	}
	lgz_add32(num, num, (str[length-1]-'0'));
}

int lgz_eqz(lgz *a)
{
	for(int i=0;i<a->length;++i)
	{
		if(a->num[i] != 0)
		{
			return 0;
		}
	}
	return 1;
}

uint8_t lgz_sign(lgz *num)
{
	if(num->num[num->length-1] < 128)
		return 0;
	if(num->num[num->length-1] >= 128)
		return 255;

	return -1;
}

int lgz_mag(lgz *num)
{
	uint8_t sign = lgz_sign(num);
	int i = num->length-1;
	while(i >= 0 && num->num[i] == sign)
	{
		--i;
	}
	return ++i; 
}

int lgz_leq(lgz *a, lgz *b)
{
	int signa = lgz_sign(a);
	int signb = lgz_sign(b);
	int lengtha = a->length;
	int lengthb = b->length;
	
	//printf("a="); lgz_print(a); printf("length=%d\n", a->length);
	//printf("b="); lgz_print(b); printf("length=%d\n", b->length);

	if(signa < signb)
		return 0;

	if(signb < signa)
		return 1;

	int a_mag = lgz_mag(a);
	int b_mag = lgz_mag(b);
	if(a_mag < b_mag)
	{
		if(signa == 255)
		{
			return 0;	
		}
		else
		{
			return 1;
		}
	}
	else if(a_mag > b_mag)
	{
		if(signa == 255)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	
	//a_mag == b_mag
	for(int i=0;i<a_mag;++i)
	{
		if(a->num[a_mag-i-1] < b->num[b_mag-i-1])
		{
			return 1;
		}
		if(b->num[b_mag-i-1] < a->num[a_mag-i-1])
		{
			return 0;
		}
	}
	return 2;
}

lgz *lgz_min(lgz *a, lgz *b)
{
	int res = lgz_leq(a, b);
	if(res == 1)
		return a;
	if(res == 0)
		return b;

	return a;
}

lgz *lgz_max(lgz *a, lgz *b)
{
	int res = lgz_leq(b, a);
	if(res == 1)
		return a;
	if(res == 0)
		return b;

	return a;
}

void lgz_lshift(lgz *a, int num)
{
	int j = a->length-1;
	int n = j-num;
	for(int i=n; i>=0; --i)
	{
		a->num[j] = a->num[i];
		--j;
	}
	for(; j>=0; --j)
	{
		a->num[j] = 0;
	}
}

void lgz_rshift(lgz *a, int num)
{
	uint8_t sign = lgz_sign(a);
	int length = a->length;
	int j = 0;
	for(int i=num; i<length; ++i)
	{
		a->num[j] = a->num[i];
		++j;	
	}
	for(;j<length;++j)
	{
		a->num[j] = sign;
	}
}

void lgz_add(lgz *c, lgz *a, lgz *b)
{
	lgz *a_oper = a;
	lgz *b_oper = b;
	
	int sum = 0;
	int carry = 0;
	uint8_t rem;
	int signa;
	int signb;
	
	if(b->length < a->length)
	{
		a_oper = b;
		b_oper = a;
	}
	
	signa = lgz_sign(a_oper);
	signb = lgz_sign(b_oper);
	
	//this is a problem if c = a or c = b	
	if(c->length < b_oper->length)
	{
		lgz_resize(c, b_oper->length);
	}
	
	for(int i=0;i<a_oper->length;++i)
	{
		sum = a_oper->num[i] + b_oper->num[i] + carry;
		c->num[i] = sum;
		carry = sum >> 8;
	}
	for(int i=a_oper->length;i<b_oper->length;++i)
	{
		sum = b_oper->num[i] + signa + carry;
		c->num[i] = sum;
		carry = sum >> 8;
	}
	rem = c->num[b_oper->length-1];

	if(c->length == b_oper->length)
	{
		if(rem < 128 && signa == 255 && signb == 255)
		{
			c->num[b_oper->length-1] = 255;
			lgz_resize(c, (c->length << 1));
			c->num[b_oper->length-1] = rem;
		}
		else if(rem >= 128 && signa == 0 && signb == 0)
		{
			c->num[b_oper->length-1] = 0;
			lgz_resize(c, (c->length << 1));
			c->num[b_oper->length-1] = rem;
		}
	}
}

void lgz_sub(lgz *c, lgz *a, lgz *b)
{	
	int sum = 0;
	int carry = 0;
	int signb = 0;
	int nsignb = 0;
	
	if(a->length < b->length)
	{
		lgz_resize(a, b->length);
	}
	
	signb = lgz_sign(b);
	nsignb = ~signb+1;
	
	if(c->length != a->length)
	{
		lgz_size_clear(c, a->length);
	}

	for(int i=0; i<b->length; ++i)
	{
		sum = a->num[i] + (~b->num[i]+1) + carry;
		c->num[i] = sum;
		carry = sum >> 8;
	}
	for(int i=b->length; i<a->length; ++i)
	{
		sum = a->num[i] + nsignb + carry;
		c->num[i] = sum;
		carry = sum >> 8;
	}
}

void lgz_mul(lgz *c, lgz *a, lgz *b)
{
	lgz *a_oper = a;
	lgz *b_oper = b;
	
	int lengthc = 0;
	int prod = 0;
	int carry = 0;
	int tmp = 0;

	if(b->length < a->length)
	{
		a_oper = b;
		b_oper = a;
	}
	
	int mix = lgz_mag(a_oper)+lgz_mag(b_oper);
	if(c->length < mix)
	{
		lgz_size_clear(c, mix);
	}

	lengthc = c->length;
	
	for(int i=0;i<a_oper->length;++i)
	{
		tmp = i;
		for(int j=0;j<b_oper->length && tmp<lengthc;++j)
		{
			prod = c->num[tmp] + (b_oper->num[j]*a_oper->num[i]) + carry;
			c->num[tmp] = (uint8_t)prod;
			carry = prod >> 8;
			++tmp;
		}
		while(carry != 0 && tmp < lengthc)
		{
			c->num[tmp] = (uint8_t)carry;
			carry >>= 8;
			++tmp;
		}
		carry = 0;
	}
}

void lgz_add32(lgz *c, lgz *a, int num)
{
	int64_t sum = 0;
	int64_t carry = 0;
	int signa = lgz_sign(a);
	if(signa == 255)
	{
		signa = ~(0);
	}

	if(c->length < a->length)
	{
		lgz_resize(c, a->length+4);
	}

	uint32_t *a_32 = (uint32_t*)a->num;
	uint32_t *c_32 = (uint32_t*)c->num;
	
	sum = (int64_t)a_32[0] + num + carry;
	c_32[0] = sum;
	carry = sum >> 32;

	for(int i=1;i<(a->length/4);++i)
	{
		sum = (int64_t)a_32[i] + signa + carry;
		c_32[i] = sum;
		carry = sum >> 32;
	}
}

void lgz_mul32(lgz *c, lgz *a, int num)
{
	int64_t prod = 0;
	uint32_t carry = 0;		
	
	int mix = lgz_mag(a)+4;
	if(c->length < mix)
	{
		lgz_resize(c,mix);
	}
	
	uint32_t *a_32 = (uint32_t*)a->num;
	uint32_t *c_32 = (uint32_t*)c->num;
	
	for(int i=0; i<(a->length/4); ++i)
	{
		prod = (int64_t)num*a_32[i];
		c_32[i] = prod + carry;
		carry = prod >> 32;
	}
}

void lgz_div(lgz_div_t *c, lgz *a, lgz *b)
{
	lgz_set(c->rem, a);
	
	lgz *quot = c->quot;
	lgz *rem = c->rem;

	int rem_mag = lgz_mag(rem);
	int b_mag = lgz_mag(b);
	int sig = b->num[b_mag-1];
	int yo = rem->num[rem_mag-1];
	int signb = lgz_sign(b);
	int nsignb = ~signb+1;
	int track = rem_mag-b_mag;
	int q = 0;
	int sum = 0;
	int carry = 0;
	
	//printf("b_mag=%d\n", b_mag);
	//printf("sig=%d\n", sig);
	
	if(quot->length <= track)
	{
		lgz_resize(quot, track+2);
	}


	//getchar();
	while(track >= 0)
	{
		q = yo/sig;
		//printf("track=%d\nyo=%x\nq=%x\n",track,yo, q);
		if(q == 0)
		{
			yo = (yo << 8) + rem->num[track+b_mag-2];
			--track;
			//getchar();
			continue;
		}
		
		int k = track;
		carry = 0;
		for(int j=0;j<b->length && k < rem->length;++j)
		{
			sum = rem->num[k] + q*(~b->num[j]+1) + carry;
			rem->num[k] = (uint8_t)sum;
			carry = sum >> 8;
			++k;
		}
		while(k < rem->length)
		{
			sum = rem->num[k] + q*nsignb + carry;
			rem->num[k] = (uint8_t)sum;
			carry = sum >> 8;
			++k;
		}
		//lgz_print(rem);
		
		//getchar();

		while(lgz_sign(rem) == 255)
		{
			k = track;
			carry = 0;
			for(int j=0;j<b->length && k < rem->length;++j)
			{
				sum = rem->num[k] + b->num[j]+carry;
				rem->num[k] = (uint8_t)sum;
				carry = sum >> 8;
				++k;
			}
			while(k < rem->length)
			{
				sum = rem->num[k] + signb + carry;
				rem->num[k] = (uint8_t)sum;
				carry = sum >> 8;
				++k;
			}
			--q;
			//lgz_print(rem);
			//getchar();
		}
		quot->num[track] = q;
		//printf("quot=");
		//lgz_print(quot);
		//getchar();
		yo = (rem->num[track+b_mag-1] << 8) + rem->num[track+b_mag-2];
		--track;
	}
}


void lgz_mod(lgz *c, lgz *a, lgz *b)
{
	lgz_set(c, a);
	lgz *rem = c;

	int rem_mag = lgz_mag(rem);
	int b_mag = lgz_mag(b);
	int sig = b->num[b_mag-1];
	int yo = rem->num[rem_mag-1];
	int signb = lgz_sign(b);
	int nsignb = ~signb+1;
	int track = rem_mag-b_mag;
	int q = 0;
	int sum = 0;
	int carry = 0;
	
	while(track >= 0)
	{
		q = yo/sig;
		if(q == 0)
		{
			yo = (yo << 8) + rem->num[track+b_mag-2];
			--track;
			continue;
		}
		
		int k = track;
		carry = 0;
		for(int j=0;j<b->length && k < rem->length;++j)
		{
			sum = rem->num[k] + q*(~b->num[j]+1) + carry;
			rem->num[k] = (uint8_t)sum;
			carry = sum >> 8;
			++k;
		}
		while(k < rem->length)
		{
			sum = rem->num[k] + q*nsignb + carry;
			rem->num[k] = (uint8_t)sum;
			carry = sum >> 8;
			++k;
		}
		while(lgz_sign(rem) == 255)
		{
			k = track;
			carry = 0;
			for(int j=0;j<b->length && k < rem->length;++j)
			{
				sum = rem->num[k] + b->num[j]+carry;
				rem->num[k] = (uint8_t)sum;
				carry = sum >> 8;
				++k;
			}
			while(k < rem->length)
			{
				sum = rem->num[k] + signb + carry;
				rem->num[k] = (uint8_t)sum;
				carry = sum >> 8;
				++k;
			}
			--q;
		}
		yo = (rem->num[track+b_mag-1] << 8) + rem->num[track+b_mag-2];
		--track;
	}
}

uint64_t zsqrt(uint64_t num)
{
	uint64_t res = 0;
	uint64_t bit = ((uint64_t)1) << 62;

	while(bit > num)
	{
		bit >>= 2;
	}
	
	while(bit != 0)
	{
		if(num >= res + bit)
		{
			num -= res + bit;
			res = (res >> 1) + bit;
		}
		else
		{
			res >>= 1;
		}
		bit >>= 2;
	}
	return res;
}

void lgz_powm(lgz *c, lgz *a, uint64_t p, lgz *n)
{
	if(p == 1)
	{
		lgz_mod(c, a, n);
		return;
	}
	if(p == 2)
	{
		lgz_memset(c, 0);
		lgz_mul(c, a, a);
		return;
	}
	uint64_t sqrtp = zsqrt(p);
	uint64_t rem = p - (sqrtp*sqrtp);
	lgz *amod = new_lgz();
	lgz *d = new_lgz();
	
	lgz_mod(amod, a, n);	

	if(sqrtp == 2)
	{
		lgz_powm(d, amod, sqrtp, n);
		lgz_memset(c, 0);
		lgz_mod(c, d, n);
		lgz_powm(d, c, sqrtp, n);
		lgz_memset(c, 0);
		lgz_mod(c, d, n);
	}
	else
	{
		lgz_memset(c, 0);
		lgz_powm(d, amod, sqrtp, n);
		lgz_powm(c, d, sqrtp, n);
	}

	for(uint64_t i=0; i<rem; ++i)
	{
		lgz_memset(d, 0);
		lgz_mul(d, amod, c);
		
		lgz_memset(c, 0);
		lgz_mod(c, d, n);
	}

	lgz_free(d);
	lgz_free(amod);
}

void lgz_sqrt(lgz *r, lgz *a)
{
	int n = lgz_mag(a);
	int m = ((n+1) >> 1)-1;	
	
	if(r->length < m)
	{
		lgz_resize(r, m);
	}

	lgz *rsq = new_lgz_size(n);
	lgz *esq = new_lgz_size(n);
	lgz *tmp1 = new_lgz_size(n);
	lgz *tmp2 = new_lgz_size(n);
	lgz *tmp3 = new_lgz_size(n);
	lgz *prev_tmp1 = new_lgz_size(n);

	esq->num[2*m] = 1;
	int c = 0;
	//printf("a="); lgz_print(a);	
	while(m >= 0)
	{
		c = 1;
		lgz_add(tmp1,rsq,esq); //tmp1 = r^2 + e^2
		//printf("tmp1="); lgz_print(tmp1);
		lgz_set(tmp2,r); //tmp2 = 2*e*r
		lgz_lshift(tmp2,m);
		lgz_mul32(tmp2,tmp2,2);
		//printf("tmp2="); lgz_print(tmp1);
		lgz_add(tmp1,tmp1,tmp2);//tmp1 = tmp1 + tmp2
		//printf("tmp1="); lgz_print(tmp1);
		if(lgz_leq(tmp1, a) != 0)
		{
			do
			{
				++c;
				//printf("C=%d\n",c);
				//(r+c*e)*(r+c*e) = r^2 + 2*c*e*r + c^2*e^2
				lgz_set(prev_tmp1,tmp1);
				lgz_mul32(tmp1,esq,c*c); //tmp1 = c^2*e^2
				//printf("tmp1="); lgz_print(tmp1);
				lgz_mul32(tmp3,tmp2,c); //tmp3 = tmp2*c = c*2*e*r
				//printf("tmp3="); lgz_print(tmp3);
				lgz_add(tmp1,tmp1,rsq);//tmp1 = tmp1+r^2 = c^2*e^2 + r^2
				//printf("tmp1="); lgz_print(tmp1);
				lgz_add(tmp1,tmp1,tmp3); //tmp1 = tmp1 + tmp3 = (r+c*e)^2
				//printf("tmp1="); lgz_print(tmp1);
				//getchar();
			}while(lgz_leq(tmp1, a));
			--c;
			lgz_set(rsq,prev_tmp1);
			r->num[m] = c;
		}
		lgz_rshift(esq,2);
		--m;
	}

	lgz_free(rsq);
	lgz_free(esq);
	lgz_free(tmp1);
	lgz_free(tmp2);
	lgz_free(tmp3);
	lgz_free(prev_tmp1);
}

void lgz_gcd(lgz *c, lgz *a, lgz *b)
{
	lgz *d = new_lgz();
	lgz_div_t *qr = new_lgz_div_t();
	
	lgz_set(c, a);
	lgz_set(qr->rem, b);
	
	while(!lgz_eqz(qr->rem))
	{
		lgz_set(d, c);
		lgz_set(c, qr->rem);
		lgz_div(qr, d, c);
	}

	lgz_free(d);
	lgz_div_t_free(qr);
}
