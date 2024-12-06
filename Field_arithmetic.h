#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

uint64_t PRM[] = {535425013,174332635,444665496,192778653,388389189,518147849,304619691,363717891,15281728,0};

uint64_t* Subtraction(uint64_t *a, uint64_t *b, int size) 
{
    uint64_t* c = (uint64_t*)calloc(size, sizeof(uint64_t));
    if (c == NULL) 
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    int64_t borrow = 0;
    for (int i = 0; i < size; i++) 
    {
        int64_t temp = (int64_t)a[i] - (int64_t)b[i] - borrow;
        
        if (temp < 0) 
        {
            temp += 0x20000000; 
            borrow = 1;          
        } else 
        {
            borrow = 0;          
        }
        
        c[i] = temp & 0x1FFFFFFF; 
    }
    return c;
}

uint64_t* Multiplication(uint64_t* a, uint64_t* b, int deg) 
{
    int s = 2*deg;
    uint64_t* c = (uint64_t *)calloc(s+1,sizeof(uint64_t));
    if(c==NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    uint64_t* d = (uint64_t *)calloc(s+2,sizeof(uint64_t));
    if(d==NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    for (int i = 0; i <= deg; i++) {
        for (int j = 0; j <= deg; j++) 
        {
            c[i + j] += a[i] * b[j];
        }
    }

    for (int i = 0; i < s; i++) 
    {
        d[i] = c[i] & 0x1FFFFFFF;
        c[i + 1] += (c[i] >> 29);
    }
    d[s] = c[s] & 0x1FFFFFFF;
    d[s+1] = (c[s] >> 29);

    free(c);
    return d;
}

int rem_greater_prm(uint64_t* rem, uint64_t* prm)
{
    for(int i = 9; i>=0; i--)
    {
        if(rem[i]>prm[i])
        {
            return 1;
        }
        else if(rem[i]<prm[i])
        {
            return 0;
        }
    }
    return 0;
}

uint64_t* Barrett_reduction(uint64_t x[])
{
    uint64_t T[] = {450887704,490307913,387807083,403879883,291135210,307268612,110539282,24605042,70628772,35};
    uint64_t stp1[10],stp3[10];  
    for(int i = 0; i<10; i++)
    {
        stp1[i] = x[i+8];
    }
    uint64_t* stp2 = Multiplication(stp1, T,9);
    for(int i = 0; i<10; i++)
    {
        stp3[i] = stp2[i+10];
    }
    uint64_t* temp = Multiplication(stp3, PRM, 9);
    uint64_t* rem = Subtraction(x, temp, 10);
    while(rem_greater_prm(rem, PRM)==1)
    {
        rem = Subtraction(rem,PRM,10);
    }
    return rem;
}

int sum_grt_prm(uint64_t* sum, uint64_t* prm)
{
    for(int i = 8; i>=0; i--)
    {
        if(sum[i]>prm[i])
        {
            return 1;
        }
        else if(sum[i]<prm[i])
        {
            return 0;
        }
    }
    return 0;
}
 
 //Addition of two 256 bit numbers in Zp*
uint64_t * Addition(uint64_t* a, uint64_t* b)
{
    uint64_t* s = (uint64_t *)calloc(9,sizeof(uint64_t));
    if(s == NULL)
    {
        fprintf(stderr,"Memory allocation field.\n");
        exit(1);
    }

    for(int i = 0; i < 9; i++)
    {
        s[i] = a[i] + b[i];
    }
    for(int i = 0; i < 8; i++)
    {
        s[i+1] = s[i+1] + (s[i] >> 29);
        s[i] = s[i] & 0x1FFFFFFF;
    }
    if(sum_grt_prm(s,PRM) == 1)
    {
        return Subtraction(s,PRM,9);
    }
    else{ return s;}
}