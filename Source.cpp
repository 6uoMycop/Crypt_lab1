#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdint.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define KEY_LEN  4
#define KEY_SEED 2,10,7,13

int readFile(uint8_t** pBuf, const char* pFileName)
{
    int iLen = 0;

    FILE* F = fopen(pFileName, "rb");
    if (F == NULL)
    {
        return -1;
    }

    fseek(F, 0L, SEEK_END);
    iLen = ftell(F);
    fseek(F, 0L, SEEK_SET);

    *pBuf = (unsigned char*)calloc(iLen, 1);
    fread(*pBuf, 1, iLen, F);
    
    return (iLen);
}

int comp(const uint8_t* a, const uint8_t* b)
{
    return *a - *b;
}

int getAlphabet(uint8_t** pAlphabet, float** pFreq, int** pQuantLetters, int iLen, uint8_t* pPlainText)
{

    *pQuantLetters = (int*)malloc(iLen * sizeof(int));

    *pAlphabet = (uint8_t*)malloc(iLen * sizeof(uint8_t));
    memcpy(*pAlphabet, pPlainText, iLen);
    qsort(*pAlphabet, iLen, sizeof(uint8_t), (int(*) (const void*, const void*)) comp);

    int ind = 0;
    int tmp = 1;
    for (int i = 1; i < iLen + 1; i++)
    {
        if ((*pAlphabet)[ind] == (*pAlphabet)[i])
        {
            (*pAlphabet)[i] = (uint8_t)0;
            tmp++;
        }
        else
        {
            (*pQuantLetters)[ind] = tmp;
            tmp = 1;

            ind++;
            (*pAlphabet)[ind] = (*pAlphabet)[i];
        }
    }

    int iNumAlphLetters = ind;
    *pAlphabet = (uint8_t*)realloc(*pAlphabet, iNumAlphLetters * sizeof(uint8_t));

    tmp = 0;
    for (int i = 0; i < iNumAlphLetters; i++)
    {
        tmp += (*pQuantLetters)[i];
    }

    *pFreq = (float*)malloc(iNumAlphLetters * sizeof(float));
    for (int i = 0; i < iNumAlphLetters; i++)
    {
        (*pFreq)[i] = (float)(*pQuantLetters)[i] / (float)tmp;
    }
    
    return (iNumAlphLetters);
}

void writeRez(uint8_t* pAlphabet, int* pQuantLetters, int iLen)
{
    FILE* A = NULL;
    FILE* F = NULL;

    A = fopen("alphabet.txt",    "wb");
    if (A == NULL)
    {
        return;
    }
    F = fopen("frequencies.txt", "w");
    if (F == NULL)
    {
        fclose(A);
        return;
    }

    fwrite(pAlphabet, 1, iLen, A);

    for (int i = 0; i < iLen; i++)
    {
        fprintf(F, "%i\n", pQuantLetters[i]);
    }

    fclose(A);
    fclose(F);
}

void vigenereEncrypt(uint8_t** pCypherText, uint8_t* pPlainText, int iLenPlain, uint8_t* pAlphabet, int iNumLetters)
{
    int iKey[KEY_LEN] = { KEY_SEED };
    
    *pCypherText = (uint8_t*)malloc(iLenPlain * sizeof(uint8_t));

    int indRez;
    int indAlph;
    for (int i = 0; i < iLenPlain; i++)
    {
        for (int j = 0; j < iNumLetters; j++)
        {
            if (pPlainText[i] == pAlphabet[j])
            {
                indAlph = j;
                break;
            }
        }
        indRez = (iKey[i % KEY_LEN] + indAlph) % iNumLetters;
        (*pCypherText)[i] = pAlphabet[indRez];
    }

    FILE* V = fopen("vigenere.txt", "wb");
    fwrite(*pCypherText, 1, iLenPlain, V);
    fclose(V);
}

int check(uint8_t a, uint8_t m)
{
    for (uint8_t c; m; )
    {
        c = a % m;
        a = m;
        m = c;
    }

    if (a == 1 || a == -1)
    {
        return 1;
    }
    return 0;
}

void initAffine(uint8_t* a, uint8_t* b, int m)
{
    uint8_t* group = (uint8_t*)malloc(m * sizeof(uint8_t));
    int groupSize = 0;

    for (int i = 1; i < m; i++)
    {
        if(check(i, m))
        {
            group[groupSize] = i;
            groupSize++;
        }
    }

    srand(time(NULL));
    *a = group[rand() % groupSize];
    *b = rand() % m;

    printf("a = %u; b = %u\n", *a, *b);

    free(group);
}

int findIndexOfSymbol(uint8_t symbol, uint8_t* pAlphabet, int iNumLetters)
{
    for (int i = 0; i < iNumLetters; i++)
    {
        if (pAlphabet[i] == symbol)
        {
            return i;
        }
    }
}
void affineEncrypt(uint8_t** pCypherText, uint8_t* pPlainText, int iLenPlain, uint8_t* pAlphabet, uint8_t a, uint8_t b, int m)
{
    *pCypherText = (uint8_t*)malloc(iLenPlain * sizeof(uint8_t));

    for (int i = 0; i < iLenPlain; i++)
    {
        (*pCypherText)[i] = pAlphabet[((a * findIndexOfSymbol(pPlainText[i], pAlphabet, m) + b) % m)];
    }

    FILE* A = fopen("affine.txt", "wb");
    fwrite(*pCypherText, 1, iLenPlain, A);
    fclose(A);
}

int main(int argc, char* argv[])
{
    uint8_t* pPlainText     = NULL;
    uint8_t* pAlphabet      = NULL;
    uint8_t* pCypherTextV   = NULL;
    uint8_t* pCypherTextA   = NULL;
    float*   pFreq          = NULL;
    int*     pQuantLetters  = NULL;
    uint8_t  a;
    uint8_t  b;

    int iLenPlain           = 0;
    int iNumAlphabetLetters = 0;

    
    if (argc != 2)
    {
        printf("Wrong args. Terminate.\n");
        return (-1);
    }

    iLenPlain = readFile(&pPlainText, argv[1]);

    if (pPlainText == NULL)
    {
        printf("Wrong file. Terminate.\n");
        return (-1);
    }

    iNumAlphabetLetters = getAlphabet(&pAlphabet, &pFreq, &pQuantLetters, iLenPlain, pPlainText);

    writeRez(pAlphabet, pQuantLetters, iNumAlphabetLetters);

    vigenereEncrypt(&pCypherTextV, pPlainText, iLenPlain, pAlphabet, iNumAlphabetLetters);

    initAffine(&a, &b, iNumAlphabetLetters);

    affineEncrypt(&pCypherTextA, pPlainText, iLenPlain, pAlphabet, a, b, iNumAlphabetLetters);

    free(pPlainText);
    free(pAlphabet);
    free(pCypherTextV);
    free(pCypherTextA);
    free(pFreq);
    free(pQuantLetters);

    system("pause");
    return (0);
}