#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <openssl/bn.h>

void status(int, int, void *) {}

void print_prime(char *prime)
{
        int i;
        for(i = 0; i < strlen(prime) && prime[i] == '0'; i++);
        for(;i < strlen(prime);i++)
                printf("%c", prime[i]);
        printf("\n");
}

int main(int argc, char *argv[])
{
        char *prime;
        BIGNUM *num_tmp;
        long int num_bits = 0;
        
        if(argc >= 2 && argv[1]) 
                num_bits = atol(argv[1]);
        else
                num_bits = 1024;
        printf("Prime generator by (c) 2004 Paolo Ardoino < paolo.ardoino@gmail.com >\n usage: ./genprimes [num_bits]\nGenerating %ld bits primes.\nWait…\n",num_bits);
        num_tmp = BN_new();
        for (int i = 0;i < 10;++ i) {
                BN_generate_prime(num_tmp,num_bits,1,NULL,NULL,status,NULL);
                prime = (char *)malloc(BN_num_bytes(num_tmp));
                prime = BN_bn2dec(num_tmp);
                print_prime(prime);
                free(prime);
		printf("--------------------------------\n");
        }
        BN_free(num_tmp);
}
