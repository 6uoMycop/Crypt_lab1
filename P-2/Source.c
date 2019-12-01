#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdint.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h> 

#define PRINT_CONST 30
#define MAX_KEY_LEN 10

int comp(const uint8_t* a, const uint8_t* b)
{
    return *a - *b;
}

int compInt(const int* a, const int* b)
{
    return *b - *a;
}

int compFloat(const void* a, const void* b)
{
    float fa = *(const float*)a;
    float fb = *(const float*)b;
    return (fa > fb) - (fa < fb);
}

int compLetters(const void* a, const void* b)
{
    struct letter
    {
        uint8_t symbol;
        int quant;
    };

    return ((struct letter*)b)->quant - ((struct letter*)a)->quant;
}

struct nGram
{
    uint8_t* pNGram;     // n-gram
    int      iN;         // n
    int      iNum;       // number of instances
    int*     pPositions; // positions of n-grams
};

int compNGrams(const void* a, const void* b)
{
    return ((struct nGram*)b)->iN - ((struct nGram*)a)->iN;
}

int readText(uint8_t** pBuf, const char* pFileName)
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

void sortAlphabet(uint8_t** pAlphabet, int iNumLetters, int** pQuantities)
{
    struct letter
    {
        uint8_t symbol;
        int quant;
    }*tmpAlphabet = NULL;

    tmpAlphabet = (struct letter*)malloc(iNumLetters * sizeof(struct letter));

    for (int i = 0; i < iNumLetters; i++)
    {
        tmpAlphabet[i].symbol = (*pAlphabet)[i];
        tmpAlphabet[i].quant = (*pQuantities)[i];
    }

    qsort(tmpAlphabet, iNumLetters, sizeof(struct letter), compLetters);

    for (int i = 0; i < iNumLetters; i++)
    {
        (*pAlphabet)[i] = tmpAlphabet[i].symbol;
        (*pQuantities)[i] = tmpAlphabet[i].quant;
    }

    free(tmpAlphabet);
}

int readAlphabet(uint8_t** pAlphabet, int** pQuant)
{
    FILE* A = fopen("alphabet.txt", "rb");
    FILE* F = fopen("frequencies.txt", "r");

    char* pBuf = NULL;
    int iNum = 0;
    //uint8_t c;

    fseek(F, 0L, SEEK_END);
    iNum = ftell(F);
    fseek(F, 0L, SEEK_SET);

    pBuf = (char*)calloc(iNum + 1, sizeof(char));
    fread(pBuf, 1, iNum, F);
    pBuf[iNum] = '\0';


    *pQuant = (int*)calloc(iNum, sizeof(int));

    char* ptr = strtok(pBuf, "\n");
    for (int i = 0; ptr != NULL; i++)
    {
        (*pQuant)[i] = atoi(ptr);
        ptr = strtok(NULL, "\n");
    }

    fseek(A, 0L, SEEK_END);
    iNum = ftell(A);
    fseek(A, 0L, SEEK_SET);

    *pAlphabet = (uint8_t*)calloc(iNum, 1);
    fread(*pAlphabet, 1, iNum, A);

    sortAlphabet(pAlphabet, iNum, pQuant);

    fclose(A);
    fclose(F);
    free(pBuf);

    return (iNum);
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

void calculateAlphabet(uint8_t* pBlock, int iBlockLen, uint8_t** pAlphabet, int* iNumLetters, int** pQuantities)
{
    *pAlphabet =   (uint8_t*)malloc(iBlockLen * sizeof(uint8_t));
    *pQuantities = (int*)malloc(iBlockLen * sizeof(int));

    memcpy(*pAlphabet, pBlock, iBlockLen);
    qsort(*pAlphabet, iBlockLen, sizeof(uint8_t), (int(*) (const void*, const void*)) comp);

    int ind = 0;
    int tmp = 1;
    for (int i = 1; i < iBlockLen + 1; i++)
    {
        if ((*pAlphabet)[ind] == (*pAlphabet)[i])
        {
            (*pAlphabet)[i] = (uint8_t)0;
            tmp++;
        }
        else
        {
            (*pQuantities)[ind] = tmp;
            tmp = 1;

            ind++;
            (*pAlphabet)[ind] = (*pAlphabet)[i];
        }
    }
    sortAlphabet(pAlphabet, ind, pQuantities);

    *iNumLetters = ind;
}

float findDelta(int iNumLetters, float* frequencies)
{
    float sum = 0;
    for (int i = 0; i < iNumLetters; i++)
    {
        sum += frequencies[i] * frequencies[i];
    }
    return sum;
}

float findIc(uint8_t* pCypherText, int iNumLetters, int ind, int iKeyLen)
{
    uint8_t* pBlockAlphabet = NULL;
    uint8_t* pBlock = NULL;
    int* pQuantBlockLetters = NULL;
    int numLetters = 0;
    int n = 0;
    int sum = 0;

    pBlock = (uint8_t*)calloc(iNumLetters, 1);
    for (int i = 0; i < iNumLetters; i++)
    {
        if (ind + i * iKeyLen < iNumLetters)
        {
            pBlock[n] = pCypherText[ind + i * iKeyLen];
            n++;
        }
    }

    if (n == 0)
    {
        free(pBlock);
        return -0.0;
    }

    calculateAlphabet(pBlock, n, &pBlockAlphabet, &numLetters, &pQuantBlockLetters);
    
    for (int i = numLetters - 1; i > 0; i--)
    {
        if (pQuantBlockLetters == 0)
        {
            iNumLetters--;
        }
    }

    for (int i = 0; i < numLetters; i++)
    {
        sum += pQuantBlockLetters[i] * (pQuantBlockLetters[i] - 1);
    }

    free(pBlockAlphabet);
    free(pQuantBlockLetters);


    //int sum = 0;
    //
    //for (int i = 0; i < iNumLetters; i++)
    //{
    //    sum += pQuantities[i] * (pQuantities[i] - 1);
    //}
    //
    //return (float)sum / (float)(iLen * (iLen - 1));

    free(pBlock);
    return ((float)sum / (float)(n * (n - 1)));
}

void friedman1(uint8_t* pCypherText, int iLen, uint8_t* pAlphabet, int iNumLetters, int* pQuantities)
{
    int numOfParts;
    float delta;
    float sum;
    float* frequencies = NULL;
    float* ICs =         NULL;
    float* subs =        NULL;

    frequencies = (float*)malloc(iNumLetters * sizeof(float));
    ICs =         (float*)malloc(iLen        * sizeof(float));
    subs =        (float*)calloc(iLen,         sizeof(float));

    for (int i = 0; i < iNumLetters; i++)
    {
        frequencies[i] = (float)pQuantities[i] / (float)iLen;
    }

    delta = findDelta(iNumLetters, frequencies);

    printf("Delta: %f\n\n", delta);
    for(int d = 2; d < ((iLen < MAX_KEY_LEN) ? iLen : MAX_KEY_LEN); d++)
    {
        printf("Supposed key length: %8i\t\t", d);

        sum = 0;
        for (int j = 0; j < d; j++)
        {
            ICs[j] = findIc(pCypherText, iNumLetters, j, d);
            printf("%.6f ", ICs[j]);
            sum += ICs[j];
        }

        for (int j = 0; j < d; j++)
        {
            float tmpSub = delta - ICs[j];
            if (tmpSub < 0.0)
            {
                tmpSub *= -1.0;
            }
            if (subs[d] == 0.0 || subs[d] > tmpSub)
            {
                subs[d] = tmpSub;
            }
        }
        printf("\n");
    }

    float minSub = delta - subs[2];
    if (minSub < 0.0)
    {
        minSub *= -1.0;
    }
    int minInd = 2;
    for (int d = 3; d < ((iLen < MAX_KEY_LEN) ? iLen : MAX_KEY_LEN); d++)
    {
        if (subs[d] < 0.0)
        {
            subs[d] *= -1.0;
        }

        if (subs[d] < minSub)
        {
            minSub = subs[d];
            minInd = d;
        }
    }

    printf("Appropriate key length: %i\nDelta = %.6f\nDifference = %.6f\n", minInd, delta, minSub);

    // 
    // cleanup:
    // 
    free(frequencies);
    free(ICs);
    free(subs);
}

void swap(int* a, int i, int j)
{
    int s = a[i];
    a[i] = a[j];
    a[j] = s;
}

void NextSet(int* a, int n)
{
    int j = n - 2;
    while (j != -1 && a[j] >= a[j + 1])
    {
        j--;
    }
    if (j == -1)
    {
        for (int i = 0; i < n; i++)
        {
            a[i] = i;
        }
        return; // больше перестановок нет. по кругу
    }
    int k = n - 1;
    while (a[j] >= a[k])
    {
        k--;
    }
    swap(a, j, k);
    int l = j + 1, r = n - 1; // сортируем оставшуюся часть последовательности
    while (l < r)
    {
        swap(a, l++, r--);
    }
}

void freqAnalysis(uint8_t* pCypherText, int iLen, uint8_t* pAlphabet, int iNumLetters, int* pQuantities)
{
    uint8_t* pCypherAlphabet = NULL;
    int* pQuantCypherLetters = NULL;
    int numCypherLetters = 0;
    uint8_t* pPlainText = NULL;
    uint8_t* pTmpChangeLetters = NULL;

    struct freqSet
    {
        int quant;                 // frequency of letter
                                   // 
        int size;                  // number of values
        uint8_t* set;              // values
                                   // 
        int* permutation;          // 
    };
    struct freqSet* sets = NULL;


    pPlainText = (uint8_t*)malloc(iLen);
    pTmpChangeLetters = (uint8_t*)malloc(iNumLetters);

    calculateAlphabet(pCypherText, iLen, &pCypherAlphabet, &numCypherLetters, &pQuantCypherLetters);

    // analyse frequencies
    int ind = 0;
    int setSize = 0;
    int allSize = 0;
    uint8_t* tmpPtr = NULL;
    for (int i = 1; i < iNumLetters + 1; i++)
    {
        setSize++;
        tmpPtr = (uint8_t*)realloc(tmpPtr, setSize);
        tmpPtr[setSize - 1] = pAlphabet[i - 1];
        if (pQuantities[ind] != pQuantities[i])
        //if (pQuantities[ind] == pQuantities[i])
        //{
        //    setSize++;
        //    tmpPtr = (uint8_t*)realloc(tmpPtr, setSize);
        //    tmpPtr[setSize - 1] = pAlphabet[ind];
        //}
        //else
        {
            allSize++;
            sets = (struct freqSet*)realloc(sets, allSize * sizeof(struct freqSet));
            sets[allSize - 1].quant = pQuantities[ind];
            sets[allSize - 1].set = tmpPtr;
            sets[allSize - 1].size = setSize;
            
            sets[allSize - 1].permutation = (int*)calloc(setSize, sizeof(int));
            for (int j = 0; j < setSize; j++)
            {
                sets[allSize - 1].permutation[j] = j;
            }
            
            tmpPtr = NULL;
            setSize = 0;
            ind = i;
        }
    }



    // match frequencies
    char userInput[32] = { '\0' };
    //int* changeIndexes = NULL;
    int tmpIndex = 0;

    


    //changeIndexes = (int*)calloc(allSize, sizeof(int));
    int changeIndex = -1;

    while (1)
    {
        memset(pTmpChangeLetters, 0, iNumLetters);

        for (int i = 0; i < allSize; i++)
        {
            printf((sets[i].size == 1) ? ("") : ((changeIndex == i) ? "[%3i] ": "%5i "), i);
        }
        printf("\n");
        for (int i = 0; i < allSize; i++)
        {
            printf((sets[i].size == 1) ? ("") : ("%5i "), sets[i].size);
        }
        printf("\n");


        ind = 0;
        for (int i = 0; i < allSize; i++)
        {
            if (changeIndex == i)
            {
                NextSet(sets[i].permutation, sets[i].size);
                for (int k = 0; k < sets[i].size; k++)
                    printf("%i ", sets[i].permutation[k]);
                printf("\n");
            }
            for (int j = 0; j < sets[i].size; j++)
            {
                pTmpChangeLetters[ind] = sets[i].set[sets[i].permutation[j]];
                ind++;
            }
        }


        for (int i = 0; i < iLen && i < 100; i++)
        {
            pPlainText[i] = pTmpChangeLetters[findIndexOfSymbol(pCypherText[i], pCypherAlphabet, iNumLetters)];
            printf("%c", pPlainText[i]);
        }

        // change letters
        //for (int i = 0; i < iLen; i++)
        //{
        //    pPlainText[i] = pAlphabet[findIndexOfSymbol(pCypherText[i], pTmpChangeLetters, iNumLetters)];
        //}

        
        printf("\n . . .\n\nContinue? (Esc - no / Enter - print decrypted text to file / ={Number from 0 to %i} - rotate position / Else - continue rotating the same position)\n> ", allSize);
        userInput[0] = _getch();
        if (userInput[0] == 27)
        {
            break;
        }
        else if (userInput[0] == 13)
        {
            FILE* outF = fopen("freqAnalysisRes.txt", "wb");
            fwrite(pPlainText, 1, iLen, outF);
            fclose(outF);
            changeIndex = -1;
        }
        else if (userInput[0] == 61)
        {
            scanf("%s", &userInput[1]);
            tmpIndex = atoi(&userInput[1]);
            if (tmpIndex > -1 && tmpIndex < allSize)
            {
                if (userInput[0] == '=')
                {
                    changeIndex = tmpIndex;
                }
            }
        }
    }

    // equation solve
    for (int i = 0; i < iNumLetters; i++)
    {
        for (int j = 0; j < iNumLetters; j++)
        {

        }
    }



    cleanup:
    free(pCypherAlphabet);
    free(pQuantCypherLetters);
    free(pPlainText);
    free(tmpPtr);
    free(pTmpChangeLetters);
    for (int i = 0; i < allSize; i++)
    {
        free(sets[i].permutation);
    }
    free(sets);
}

float findMIc(uint8_t* pPart1, int iPartNum1, uint8_t* pPart2, int iPartNum2, uint8_t* pAlphabet, int iNumLetters, int* pQuantities)
{
    uint8_t* pBlockAlphabet1 = NULL;
    uint8_t* pBlockAlphabet2 = NULL;
    int* pQuantBlockLetters1 = NULL;
    int* pQuantBlockLetters2 = NULL;
    int* pQuantBlock1 = NULL;
    int* pQuantBlock2 = NULL;
    int numLetters1 = 0;
    int numLetters2 = 0;
    int sum = 0;

    pQuantBlock1 = (int*)calloc(iNumLetters, sizeof(int));
    pQuantBlock2 = (int*)calloc(iNumLetters, sizeof(int));

    calculateAlphabet(pPart1, iPartNum1, &pBlockAlphabet1, &numLetters1, &pQuantBlockLetters1);
    calculateAlphabet(pPart2, iPartNum2, &pBlockAlphabet2, &numLetters2, &pQuantBlockLetters2);

    for (int i = 0; i < numLetters1; i++)
    {
        for (int j = 0; j < iNumLetters; j++)
        {
            if (pBlockAlphabet1[i] == pAlphabet[j])
            {
                pQuantBlock1[j] = pQuantBlockLetters1[i];
                break;
            }
        }
    }
    for (int i = 0; i < numLetters2; i++)
    {
        for (int j = 0; j < iNumLetters; j++)
        {
            if (pBlockAlphabet2[i] == pAlphabet[j])
            {
                pQuantBlock2[j] = pQuantBlockLetters2[i];
                break;
            }
        }
    }
    
    for (int i = 0; i < iNumLetters; i++)
    {
        sum += pQuantBlock1[i] * pQuantBlock2[i];
    }


    free(pBlockAlphabet1);
    free(pBlockAlphabet2);
    free(pQuantBlockLetters1);
    free(pQuantBlockLetters2);
    free(pQuantBlock1);
    free(pQuantBlock2);
    return (float)sum / (float)(iPartNum1 * iPartNum2);
}

int findK0(int** k0, uint8_t* pPart, int iPartNum, uint8_t* pAlphabet, int iNumLetters, int* pQuantities)
{
    uint8_t* pBlockAlphabet = NULL;
    int* pQuantBlockLetters = NULL;
    //int* pQuantBlock = NULL;
    int numLetters = 0;

    //pQuantBlock = (int*)calloc(iNumLetters, sizeof(int));

    calculateAlphabet(pPart, iPartNum, &pBlockAlphabet, &numLetters, &pQuantBlockLetters);

    int index = 0;
    for (int i = 0; i < iPartNum; i++)
    {
        if (pQuantBlockLetters[0] == pQuantBlockLetters[i])
        {
            *k0 = (int*)realloc(*k0, (index + 1) * sizeof(int));
            (*k0)[index] = findIndexOfSymbol(pBlockAlphabet[0], pAlphabet, iNumLetters) - findIndexOfSymbol(pAlphabet[i], pAlphabet, iNumLetters);
            while ((*k0)[index] < 0)
            {
                (*k0)[index] += iNumLetters;
            }
            (*k0)[index] %= iNumLetters;
            index++;
        }
    }

    free(pBlockAlphabet);
    free(pQuantBlockLetters);
    //free(pQuantBlock);
    return index;
}

void friedman2(uint8_t* pCypherText, int iLen, uint8_t* pAlphabet, int iNumLetters, int* pQuantities)
{
    float     delta;
    float     sum;
    int       iMaxNumOfSymbols;
    float*    frequencies =        NULL;
    float*    MICs =               NULL;
    uint8_t** cypheredParts =      NULL; // blocks of bytes which have been encrypted on the same key's byte
    int*      numOfSymbolsInPart = NULL; // number of bytes in a part
    int*      k_0 =                NULL; // possible k_0's
    int*      offsets =            NULL;

    numOfSymbolsInPart = (int*)calloc(iLen, sizeof(int));
    offsets = (int*)calloc(iLen, sizeof(int));
    MICs = (float*)malloc(iLen * iNumLetters * sizeof(float));

    cypheredParts = (uint8_t**)malloc(iLen * sizeof(uint8_t*));
    for (int i = 0; i < iLen; i++)
    {
        cypheredParts[i] = (uint8_t*)calloc(iLen, 1);
    }

    frequencies = (float*)malloc(iNumLetters * sizeof(float));
    for (int i = 0; i < iNumLetters; i++)
    {
        frequencies[i] = (float)pQuantities[i] / (float)iLen;
    }

    delta = findDelta(iNumLetters, frequencies);
    printf("Delta: %f\n\n", delta);

    uint8_t* tmpYj = NULL;
    tmpYj = (uint8_t*)malloc(iLen);
    for (int d = 2; d < ((iLen < MAX_KEY_LEN) ? iLen : MAX_KEY_LEN); d++)
    {
        int* tmpKey = (int*)calloc(d, sizeof(int));
        printf("Supposed key length: %8i\n", d);

        // divide cyphertext to blocks of bytes which have been encrypted on the same key's byte
        for (int i = 0; i < d; i++)
        {
            memset(cypheredParts[i], 0, iLen);
        }
        memset(numOfSymbolsInPart, 0, iLen);
        for (int i = 0; i < iLen; i++)
        {
            cypheredParts[i % d][numOfSymbolsInPart[i % d]] = pCypherText[i];
            numOfSymbolsInPart[i % d]++;
        }

        for (int i = 0; i < d; i++)         // for all pairs Yi, Yj
        {                                   //
            for (int j = i + 1; j < d; j++) // 
            {
                int tmpMIind = 0;
                memset(MICs, 0, iLen * iNumLetters * sizeof(float));

                for (int r = 1; r < iNumLetters - 1; r++) // for all possible rotations (r = 1..|A|-1)
                {
                    for (int p = 0; p < numOfSymbolsInPart[j]; p++) // encrypt Yj by ROTr()
                    {
                        tmpYj[p] = pAlphabet[(findIndexOfSymbol(cypheredParts[j][p], pAlphabet, iNumLetters) + r) % iNumLetters];
                    }
                    // calculate MIc(Yi,ROTr(Yj))
                    MICs[tmpMIind] = findMIc(cypheredParts[i], numOfSymbolsInPart[i], tmpYj, numOfSymbolsInPart[j], pAlphabet, iNumLetters, pQuantities);
                    tmpMIind++;
                }

                // find minimal difference
                float diff = delta - MICs[0];
                float tmpDiff;
                int ind = 0;
                if (diff < 0)
                {
                    diff *= -1.0;
                }
                //printf("\n%.6f ", MICs[0]);
                for (int p = 1; p < tmpMIind; p++)
                {
                    //printf("%.6f ", MICs[p]);
                    if (MICs[p] == 0)
                    {
                        continue;
                    }
                    tmpDiff = delta - MICs[p];
                    if (tmpDiff < 0)
                    {
                        tmpDiff *= -1.0;
                    }
                    if(tmpDiff < diff)
                    {
                        diff = tmpDiff;
                        ind = p;
                    }
                }
                printf("\n");
                offsets[j] = ind;
                printf("Difference delta-MIc=%.6f; k_%i-k_%i= %i\n", diff, i, j, ind);
            }
            break; // check only differences k_0 - k_j, j=1..d-1
        }

        // find k_0 using frequency analysis
        free(k_0);
        k_0 = NULL;
        int sizeK0 = findK0(&k_0, cypheredParts[0], numOfSymbolsInPart[0], pAlphabet, iNumLetters, pQuantities);
        if (k_0 == NULL)
        {
            continue;
        }

        for (int i = 0; i < sizeK0; i++)
        {
            tmpKey[0] = k_0[i];
            printf("Supposed key: %i ", tmpKey[0]);
            for (int j = 1; j < d; j++) // generate key
            {
                tmpKey[j] = tmpKey[0] + offsets[j];
                printf("%i ", tmpKey[j]);
            }
            printf("\n");
        }
        free(tmpKey);
        printf("\n");
    }

    //printf("Appropriate key length: %i\nDelta = %.6f\nDifference = %.6f\n", minInd, delta, minSub);

    // 
    // cleanup:
    // 
    free(frequencies);
    free(MICs);
    free(tmpYj);

    for (int i = 0; i < iLen; i++)
    {
        free(cypheredParts[i]);
    }
    free(cypheredParts);
    free(numOfSymbolsInPart);
    free(offsets);
    free(k_0);
}

void insertNGram(struct nGram** pNGrams, int *iSizeNGrams, uint8_t* pTmpNGram, int iN, int iPos0, int iPos1)
{
    if (*iSizeNGrams != 0)
    {
        for (int i = 0; i < *iSizeNGrams; i++)
        {
            // if found:
            if ((*pNGrams)[i].iN == iN && memcmp((*pNGrams)[i].pNGram, pTmpNGram, iN) == 0)
            {
                (*pNGrams)[i].pPositions = (int*)realloc((*pNGrams)[i].pPositions, ((*pNGrams)[i].iNum + 1) * sizeof(int));
                (*pNGrams)[i].pPositions[(*pNGrams)[i].iNum] = iPos1;
                (*pNGrams)[i].iNum++;
                return;
            }
        }
    }
    // if not found:
    *pNGrams = (struct nGram*)realloc(*pNGrams, ((*iSizeNGrams) + 1) * sizeof(struct nGram));
    (*pNGrams)[(*iSizeNGrams)].pPositions = (int*)malloc(2 * sizeof(int));
    (*pNGrams)[(*iSizeNGrams)].pPositions[0] = iPos0;
    (*pNGrams)[(*iSizeNGrams)].pPositions[1] = iPos1;
    (*pNGrams)[(*iSizeNGrams)].pNGram = (uint8_t*)malloc(iN);
    memcpy((*pNGrams)[(*iSizeNGrams)].pNGram, pTmpNGram, iN);
    (*pNGrams)[(*iSizeNGrams)].iN = iN;
    (*pNGrams)[(*iSizeNGrams)].iNum = 2;
    (*iSizeNGrams)++;

    //qsort(*pNGrams, *iSizeNGrams, sizeof(struct nGram), compNGrams);
}

int GCD(int a, int b)
{
    while (a > 0 && b > 0)
    {
        if (a > b)
        {
            a %= b;
        }
        else
        {
            b %= a;
        }
    }
    return a + b;
}

void kasiski(uint8_t* pCypherText, int iLen)
{
    struct nGram* pNGrams =    NULL;
    int*          pGCDs = NULL;
    int           sizeNGrams = 0;

    uint8_t*      tmpNGram =   NULL;
    int           tmpN =       0;

    tmpNGram = (uint8_t*)malloc(iLen);

    for (int i = 0; i < iLen; i++)
    {
        for (int j = i + 2; j < iLen; j++)
        {

            if (pCypherText[i] == pCypherText[j]) // found match
            {
                memset(tmpNGram, 0, iLen);
                tmpN = 0;

                // check if it is n-gram, n > 1
                for (int p = i, q = j; pCypherText[p] == pCypherText[q] && p < j && q < iLen; p++, q++)
                {
                    tmpNGram[tmpN] = pCypherText[p];
                    tmpN++;
                    if (tmpN > 1) // save n-gram
                    {
                        insertNGram(&pNGrams, &sizeNGrams, tmpNGram, tmpN, i, j);
                    }
                }
            }
        }
    }
    
    // sort positions and delete repetitions
    for (int p = 0; p < sizeNGrams; p++)
    {
        qsort(pNGrams[p].pPositions, pNGrams[p].iNum, sizeof(int), compInt);
        int ind = 0;
        for (int i = 1; i < pNGrams[p].iNum + 1; i++)
        {
            if (pNGrams[p].pPositions[ind] == pNGrams[p].pPositions[i])
            {
                pNGrams[p].pPositions[i] = 0;
            }
            else
            {
                ind++;
                pNGrams[p].pPositions[ind] = pNGrams[p].pPositions[i];
            }
        }
        pNGrams[p].iNum = ind;
    }

    // calculate GCDs
    pGCDs = (int*)malloc((sizeNGrams + 1) * sizeof(int));
    for (int n = 0; n < sizeNGrams; n++)
    {
        pGCDs[n] = GCD(pNGrams[n].pPositions[0], pNGrams[n].pPositions[1]);
        for (int i = 2; i < pNGrams[n].iNum; i++)
        {
            pGCDs[n] = GCD(pGCDs[n], pNGrams[n].pPositions[i]);
        }
    }

    printf("n-grams:\n\n");
    for (int i = 0; i < sizeNGrams; i++)
    {
        printf("n-gram: \"");
        for (int j = 0; j < pNGrams[i].iN; j++)
        {
            printf("%c", pNGrams[i].pNGram[j]);
        }
        printf("\"\nn=%i\nGCD:%i\nNumber of instances: %i\nPositions: ", pNGrams[i].iN, pGCDs[i], pNGrams[i].iNum);
        for (int j = 0; j < pNGrams[i].iNum; j++)
        {
            printf("%6i ", pNGrams[i].pPositions[j]);
        }
        printf("\n\n");
    }

    printf("\n\n");

    qsort(pGCDs, sizeNGrams, sizeof(int), compInt);

    int tmpNum = 0;
    for (int i = 0; i < sizeNGrams; i++)
    {
        tmpNum++;
        if (pGCDs[i] != pGCDs[i + 1])
        {
            printf((pGCDs[i] == 1) ? "" : "GCD (supposed key length) %4i - %4i times\n", pGCDs[i], tmpNum);
            tmpNum = 0;
        }
    }
    printf("\n\n");

    // cleanup
    for (int i = 0; i < sizeNGrams; i++)
    {
        free(pNGrams[i].pPositions);
        free(pNGrams[i].pNGram);
    }
    free(pNGrams);
    free(tmpNGram);
    free(pGCDs);
}

int main()
{
    uint8_t* pAlphabet      = NULL;
    uint8_t* pCypherTextV   = NULL;
    uint8_t* pCypherTextA   = NULL;
    int*     pQuantLetters  = NULL;

    int iLenTextV           = 0;
    int iLenTextA           = 0;
    int iNumAlphabetLetters = 0;
    
    iNumAlphabetLetters = readAlphabet(&pAlphabet, &pQuantLetters);

    iLenTextV = readText(&pCypherTextV, "vigenere.txt");
    iLenTextA = readText(&pCypherTextA, "affine.txt");
    
    //freqAnalysis(pCypherTextA, iLenTextA, pAlphabet, iNumAlphabetLetters, pQuantLetters);
    //friedman1(pCypherTextV, iLenTextV, pAlphabet, iNumAlphabetLetters, pQuantLetters);
    //friedman2(pCypherTextV, iLenTextV, pAlphabet, iNumAlphabetLetters, pQuantLetters);
    kasiski(pCypherTextV, iLenTextV);


    free(pAlphabet);
    free(pCypherTextV);
    free(pCypherTextA);
    free(pQuantLetters);

    system("pause");
    return 0;
}
