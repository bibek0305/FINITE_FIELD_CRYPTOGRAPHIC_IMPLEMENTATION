# include"Field_arithmetic.h"
typedef struct 
{
    uint64_t* x;
    uint64_t* y;
}Point;


uint64_t* Exponent_LTR(uint64_t *gen, uint64_t *exp )
{
    uint64_t* h = (uint64_t *)calloc(9, sizeof(uint64_t));
    h[0] = 1;
    for(int i = 8; i >= 0; i--)
    {
        unsigned char nbits[29];
        int l=0;
        for(int j = 28; j >= 0; j--)
        {
            nbits[l++] = (exp[i] >> j) & 1;
        }
        for(int k = 0; k < 29; k++)
        {
            h = Barrett_reduction(Multiplication(h, h, 8));
            if(nbits[k] == 1)
            {
                h = Barrett_reduction(Multiplication(gen, h, 8));
            }
        }
    }
    return h;
}

uint64_t *Exponent_RTL(uint64_t *gen, uint64_t*exp)
{
    uint64_t *h = (uint64_t *) calloc(9,sizeof(uint64_t));
    h[0]=1;
    for(int i=0; i<=8; i++)
    {
        unsigned char nbits[29];
        int l=0;
        for(int j =0; j<=28; j++)
        {
            nbits[l++] = (exp[i]>>j) & 1;
        }
        for(int m = 0; m<=28; m++)
        {
            if(nbits[m] == 1)
            {
                h = Barrett_reduction(Multiplication(h, gen, 8));
            }
            gen = Barrett_reduction(Multiplication(gen,gen,8));
        }
    }
    return h;
}

uint64_t* Exponent_Montgomery(uint64_t *gen, uint64_t *exp )
{
    uint64_t* S = (uint64_t *)calloc(9, sizeof(uint64_t));
    S[0] = 1;
    uint64_t *R;
    R = gen;
    for(int i = 8; i >= 0; i--)
    {
        unsigned char nbits[29];
        int l=0;
        for(int j = 28; j >= 0; j--)
        {
            nbits[l++] = (exp[i] >> j) & 1;
        }
        for(int k = 0; k < 29; k++)
        {
            if(nbits[k] == 0)
            {
                R = Barrett_reduction(Multiplication(S, R, 8));
                S = Barrett_reduction(Multiplication(S, S, 8));
            }
            else
            {
                S = Barrett_reduction(Multiplication(S, R, 8));
                R = Barrett_reduction(Multiplication(R, R, 8));
            }
        }
    }
    return S;
}

// Inverse of an element in Zp*
uint64_t *Inverse(uint64_t * z)
{
    uint64_t b[] = {2, 0, 0, 0, 0, 0, 0, 0, 0};
    return Exponent_RTL(z,Subtraction(PRM,b,9));
}

uint64_t* Mult_zp(uint64_t* a, uint64_t* b)
{
    return Barrett_reduction(Multiplication(a,b,8));
}

uint64_t* Neg_val(uint64_t* a)
{
    return Subtraction(PRM, a, 9);
}

int Is_not_equal(uint64_t* a, uint64_t* b)
{
    int flag = 0;
    for(int i = 0; i < 9; i++)
    {
        if(a[i] != b[i])
        {
            flag = 1;
            break;
        }
    }
    return flag;
} 

int Is_not_zero(uint64_t* a)
{
    int flag = 0;
    for(int i = 0; i < 9; i++)
    {
        if(a[i] != 0)
        {
            flag = 1;
            break;
        }
    }
    return flag;
}

Point Doubling(Point a)
{
    uint64_t temp[] = {8,0,0,0,0,0,0,0,0};
    uint64_t* A = Neg_val(temp);
    Point Q;
    if(Is_not_zero(a.y) == 0)
    {
        uint64_t k[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        Q.x = k;
        Q.y = k;
        return Q;
    }
    else
    {
        uint64_t r[9] = {3, 0, 0, 0, 0, 0, 0, 0, 0};
        uint64_t s[9] = {2, 0, 0, 0, 0, 0, 0, 0, 0};
        uint64_t *temp1, *temp2, *temp3, *temp4, *temp5;
        temp1 = Mult_zp( Mult_zp( a.x, a.x ), r);//3x1^2
        temp1 = Addition(temp1, A);//(3x1^2+A)
        temp2 = Mult_zp(s, a.y);//2y1
        temp3 = Mult_zp(temp1, Inverse(temp2));//(3x1^2+A)/2y1
        temp4 = Mult_zp(temp3, temp3);//((3x1^2+A)/2y1)^2
        Q.x = Addition(temp4, Neg_val(Mult_zp(s, a.x)));//((3x1^2+A)/2y1)^2-2x1
        temp5 = Addition(a.x, Neg_val(Q.x));//(x1-x3)
        Q.y = Addition( Mult_zp(temp3, temp5), Neg_val(a.y));//((3x1^2+A)/2y1)(x1-x3)-y1
        return Q;
    }
}

Point Elliptic_add(Point a, Point b)//a(x1,y1) & b(x2,y2)
{
    Point P;
    if((Is_not_zero(a.x)==0 && Is_not_zero(a.y)==0) || (Is_not_zero(b.x)==0 && Is_not_zero(b.y)==0))
    //Either of two points is identiy element (0,0)
    {
        if (Is_not_zero(a.x)==0 && Is_not_zero(a.y)==0)
        {
            P.x = b.x;
            P.y = b.y;
        }
        else
        {
            P.x = a.x;
            P.y = a.y;
        }
        return P;
    }
    else
    {
        if(Is_not_equal(a.x , b.x) == 1)
        {
            uint64_t *temp1, *temp2, *temp3, *temp4, *temp5;
            temp1 = Addition(b.y, Neg_val(a.y));//(y2-y1)
            temp2 = Addition(b.x, Neg_val(a.x));//(x2-x1)
            temp3 = Mult_zp(temp1, Inverse(temp2));//(y2-y1)/(x2-x1)
            temp4 = Mult_zp(temp3, temp3);//((y2-y1)/(x2-x1))^2
            P.x = Addition(Addition(temp4, Neg_val(a.x)), Neg_val(b.x));//((y2-y1)/(x2-x1))^2-x1-x2
            temp5 = Addition(a.x, Neg_val(P.x));//(x1-x3)
            P.y = Addition(Mult_zp(temp3, temp5), Neg_val(a.y));//((y2-y1)/(x2-x1))(x1-x3)-y1
            return P;
        }
        else if((Is_not_equal(a.x, b.x) == 0) && (Is_not_equal(a.y, b.y) == 1))
        {
            uint64_t z[9] = {0};
            P.x = z;
            P.y = z;
            return P;
        }
        else
        {
            Doubling(a);
        }
    }
}

Point Scalar_multiplication_RTL(Point gen, uint64_t* exp)
{
    Point H;
    uint64_t array[9] = {0};
    H.x = array;
    H.y = array;
    for(int i = 0; i< 9; i++)
    {
        unsigned char nbits[29];
        int l=0;
        for(int j =0; j< 29; j++)
        {
            nbits[l++] = (exp[i]>>j) & 1;
        }
        for(int m = 0; m<=28; m++)
        {
            if(nbits[m] == 1)
            {
                H = Elliptic_add(H, gen);
            }
            gen = Doubling(gen);
        }
    }
    return H;
}

Point Scalar_multiplication_LTR(Point gen, uint64_t *exp )
{
    Point H;
    uint64_t array[9] = {0};
    H.x = array;
    H.y = array;
    for(int i = 8; i >= 0; i--)
    {
        unsigned char nbits[29];
        int l=0;
        for(int j = 28; j >= 0; j--)
        {
            nbits[l++] = (exp[i] >> j) & 1;
        }
        for(int k = 0; k < 29; k++)
        {
            H = Doubling(H);
            if(nbits[k] == 1)
            {
                H = Elliptic_add(H, gen);
            }
        }
    }
    return H;
}

// int main() {
//     // Define example points on the curve
//     uint64_t x1[9] = {67761129,
//  194562371,
//  350961085,
//  359779665,
//  514199853,
//  89720943,
//  299790580,
//  440825741,
//  13008848};
//     uint64_t y1[9] = {132682569,
//  200877299,
//  22644068,
//  402358525,
//  452581711,
//  277320122,
//  122613,
//  404562878,
//  8634712
// };

//     uint64_t x2[9] = {189922301,
//  234269113,
//  407949731,
//  400047267,
//  183356059,
//  369030301,
//  114949782,
//  488582683,
//  5929631};
//      uint64_t y2[9] = {320179514,
//  202838025,
//  491063075,
//  532933308,
//  260827449,
//  398024043,
//  66459634,
//  201223031,
//  3902662};

//     Point P1 = {x1, y1};
//     Point P2 = {x2, y2};

//     Point result_add = Elliptic_add(P1, P2);
//     printf("Elliptic_add Result:\n");
//     printf("x: ");
//     for (int i = 0; i < 9; i++) printf("%lu ", result_add.x[i]);
//     printf("\ny: ");
//     for (int i = 0; i < 9; i++) printf("%lu ", result_add.y[i]);
//     printf("\n");

//     return 0;  
// }

int main() {
   // Initialize generator point
    uint64_t gen_x[9] = {493191603,
 215291332,
 272113980,
 485734153,
 462269611,
 14322031,
 128322386,
 220900844,
 14192347};
    uint64_t gen_y[9] = {228135575,
 106391292,
 518512409,
 154429979,
 213301539,
 406308432,
 323086405,
 48352955,
 14936292};
    Point generator;
    generator.x = gen_x;
    generator.y = gen_y;

    // Example exponent
    uint64_t exp[9] = {100, 0, 0, 0, 0, 0, 0, 0, 0};

    // Perform scalar multiplication
    Point result = Scalar_multiplication_RTL(generator, exp);

    // Output result
    printf("Resulting Point:\n");
    printf("X: ");
    for (int i = 0; i < 9; i++) {
        printf("%lu ", result.x[i]);
    }
    printf("\nY: ");
    for (int i = 0; i < 9; i++) {
        printf("%lu ", result.y[i]);
    }
    printf("\n");

    return 0;
}
