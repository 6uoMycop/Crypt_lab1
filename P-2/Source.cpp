#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdint.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h> 

#define KEY_LEN 4
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

/*
void freqAnalysis(uint8_t* pCypherText, int iLen, uint8_t* pAlphabet, int iNumLetters, int* pQuantities)
{
    int       iMaxNumOfSymbols =            (iLen % KEY_LEN == 0) ? (iLen / KEY_LEN) : (iLen / KEY_LEN + 1); // max number of symbols in part
    uint8_t** cypheredParts =               NULL; // blocks of bytes which have been encrypted on the same key's byte
    int*      numOfSymbolsInPart =          NULL; // number of bytes in a part
    int*      tmpFreqOfBlockSymbols =       NULL; // quantity of certain symbols in a block
    int*      tmpFreqOfBlockSymbolsSorted = NULL; // sorted quantity of certain symbols in a block
    uint8_t** rezSymbols =                  NULL; // array of symbols with maximal frequency
    int*      numOfRezSymbols =             NULL; // quantity of candidats for certain block


    cypheredParts =               (uint8_t**)malloc(KEY_LEN * sizeof(uint8_t*));
    rezSymbols =                  (uint8_t**)malloc(KEY_LEN * sizeof(uint8_t*));
    numOfSymbolsInPart =          (int*)calloc(KEY_LEN, sizeof(int));
    tmpFreqOfBlockSymbols =       (int*)malloc(iNumLetters * sizeof(int));
    tmpFreqOfBlockSymbolsSorted = (int*)malloc(iNumLetters * sizeof(int));
    numOfRezSymbols =             (int*)calloc(KEY_LEN, sizeof(int));

    for (int i = 0; i < KEY_LEN; i++)
    {
        cypheredParts[i] = (uint8_t*)calloc(iMaxNumOfSymbols, 1);
        rezSymbols[i]    = NULL;
    }

    // divide cyphertext to blocks of bytes which have been encrypted on the same key's byte
    for (int j = 0; j < iLen; j += KEY_LEN, numOfSymbolsInPart[i]++)
    {
        cypheredParts[i][numOfSymbolsInPart[i]] = pCypherText[j];
    }
    
    // get quantities of symbols
    uint8_t* tmpAlph = NULL;
    int* tmpQuant = NULL;
    
    calculateAlphabet(cypheredParts[]);
    
    free(tmpAlph);
    free(tmpQuant);
    
    
    ////
        //// find the most frequent symbol(s) in the block
        ////
        //
        //memset(tmpFreqOfBlockSymbols, 0, iNumLetters * sizeof(int));
        //for (int j = 0; j < numOfSymbolsInPart[i]; j++)
        //{
        //    tmpFreqOfBlockSymbols[findIndexOfSymbol(cypheredParts[i][j], pAlphabet, iNumLetters)]++;
        //}
        //
        //// find the maximal quantity
        //int maxQuant = 1;
        //memcpy(tmpFreqOfBlockSymbolsSorted, tmpFreqOfBlockSymbols, iNumLetters * sizeof(int));
        //qsort(tmpFreqOfBlockSymbolsSorted, iNumLetters, sizeof(int), (int(*) (const void*, const void*)) compInt);
        //while (tmpFreqOfBlockSymbolsSorted[maxQuant - 1] == tmpFreqOfBlockSymbolsSorted[maxQuant] && maxQuant < iNumLetters)
        //{
        //    maxQuant++;
        //}
        //
        //// find all the symbols with maximal quantity
        //for (int j = 0; j < iNumLetters; j++)
        //{
        //    int tmpCntr = -1;
        //    int tmpInd = 0;
        //    for (int k = 0; k < iNumLetters; k++)
        //    {
        //        if (tmpFreqOfBlockSymbols[k] == 0)
        //        {
        //            continue;
        //        }
        //        if (tmpFreqOfBlockSymbols[j] == tmpFreqOfBlockSymbols[k])
        //        {
        //            tmpCntr++;
        //            tmpInd = j;
        //        }
        //    }
        //    if (tmpCntr == maxQuant - 1)
        //    {
        //        if (rezSymbols[i] == NULL)
        //        {
        //            rezSymbols[i] = (uint8_t*)malloc(1);
        //        }
        //        else
        //        {
        //            rezSymbols[i] = (uint8_t*)realloc(rezSymbols[i], numOfRezSymbols[i] + 1);
        //        }
        //        rezSymbols[i][numOfRezSymbols[i]] = pAlphabet[tmpInd];
        //        numOfRezSymbols[i]++;
        //    }
        //}
    




    
    //
    // find the most frequent symbol(s) in alphabet
    //
    uint8_t* mostFrequentSymbols = NULL;
    int numOfMostFrequent = 0;

    // find the maximal quantity
    int maxQuantity = 1;
    memcpy(tmpFreqOfBlockSymbolsSorted, pQuantities, iNumLetters * sizeof(int));
    qsort(tmpFreqOfBlockSymbolsSorted, iNumLetters, sizeof(int), (int(*) (const void*, const void*)) compInt);
    maxQuantity = tmpFreqOfBlockSymbolsSorted[0];

    // find all the symbols with maximal quantity
    for (int i = 0, tmpQuant = 0; i < iNumLetters; i++)
    {
        if (pQuantities[i] == maxQuantity)
        {
            if (mostFrequentSymbols == NULL)
            {
                mostFrequentSymbols = (uint8_t*)malloc(1);
            }
            else
            {
                mostFrequentSymbols = (uint8_t*)realloc(mostFrequentSymbols, numOfMostFrequent + 1);
            }
            mostFrequentSymbols[numOfMostFrequent] = pAlphabet[i];
            numOfMostFrequent++;
        }
    }

    // 
    // matching of symbols and find key
    // 
    //char userInput[32];
    int indexInMostFrequentSymbols = 0; // indicates current symbol to choose in alphabet
    int indexesInRezSymbols[4] = { 0 }; // indicate current symbol in every block of resulting blocks

    

    for (int p = 0; p < numOfRezSymbols[0]; p++)
    {
        for (int q = 0; q < numOfRezSymbols[1]; q++)
        {
            for (int r = 0; r < numOfRezSymbols[2]; r++)
            {
                for (int s = 0; s < numOfRezSymbols[3]; s++)
                {
                    for (int t = 0; t < numOfMostFrequent; t++)
                    {
                        int indArr[KEY_LEN] = { p, q, r, s };
                        uint8_t tmpKey[KEY_LEN] = { 0 };

                        printf("- - - Finding key - - -\n\n");

                        for (int i = 0; i < KEY_LEN; i++)
                        {
                            int tmpSymbX = findIndexOfSymbol(rezSymbols[i][indArr[i]], pAlphabet, iNumLetters);
                            int tmpSymbY = findIndexOfSymbol(mostFrequentSymbols[t], pAlphabet, iNumLetters);
                            printf("\nSearching key byte N%i:\n\n", i);
                            tmpKey[i] = (tmpSymbX - tmpSymbY + iNumLetters) % iNumLetters;

                            if (numOfRezSymbols[i] == 1)
                            {
                                if (numOfMostFrequent == 1) // case 0 - can match automaticaly
                                {
                                    printf("Auto match: %c (%u) -> %u => K[%i] = %u\n",
                                        rezSymbols[i][0],
                                        rezSymbols[i][0],
                                        mostFrequentSymbols[0],
                                        i,
                                        tmpKey[i]);
                                }
                                else // case 1.0 - one multi choice
                                {
                                    printf("There are %i most frequent symbols in language:\n", numOfMostFrequent);
                                    for (int j = 0; j < numOfMostFrequent; j++)
                                    {
                                        printf("%c (%u)\t", mostFrequentSymbols[j], mostFrequentSymbols[j]);
                                    }
                                    printf("\nVariative match: %c (%u) -> %u => K[%i] = %u\n",
                                        rezSymbols[i][0],
                                        rezSymbols[i][0],
                                        mostFrequentSymbols[indexInMostFrequentSymbols],
                                        i,
                                        tmpKey[i]);
                                }
                            }
                            else
                            {
                                if (numOfMostFrequent == 1) // case 1.1 - one multi choice
                                {
                                    printf("There are %i most frequent symbols in resulting set:\n", numOfRezSymbols[i]);
                                    for (int j = 0; j < numOfRezSymbols[i]; j++)
                                    {
                                        printf("%c (%u)\t", rezSymbols[i][j], rezSymbols[i][j]);
                                    }
                                    printf("\nVariative match: %c (%u) -> %u => K[%i] = %u\n",
                                        rezSymbols[i][indexesInRezSymbols[i]],
                                        rezSymbols[i][indexesInRezSymbols[i]],
                                        mostFrequentSymbols[0],
                                        i,
                                        tmpKey[i]);
                                }
                                else // case 2 - two multi choices
                                {
                                    printf("There are %i most frequent symbols in language:\n", numOfMostFrequent);
                                    for (int j = 0; j < numOfMostFrequent; j++)
                                    {
                                        printf("%c (%u)\t", mostFrequentSymbols[j], mostFrequentSymbols[j]);
                                    }
                                    printf("and there are %i most frequent symbols in resulting set:\n", numOfRezSymbols[i]);
                                    for (int j = 0; j < numOfRezSymbols[i]; j++)
                                    {
                                        printf("%c (%u)\t", rezSymbols[i][j], rezSymbols[i][j]);
                                    }
                                    printf("\nVariative match: %c (%u) -> %u => K[%i] = %u\n",
                                        rezSymbols[i][indexesInRezSymbols[i]],
                                        rezSymbols[i][indexesInRezSymbols[i]],
                                        mostFrequentSymbols[indexInMostFrequentSymbols],
                                        i,
                                        tmpKey[i]);
                                }
                            }
                        }
                        printf("\n- - - - - - -\n\nKey: %u %u %u %u\n\nDecrypted text:\n", tmpKey[0], tmpKey[1], tmpKey[2], tmpKey[3]);
                        int indAlph = 0;
                         for (int i = 0; i < PRINT_CONST && i < iLen; i++)
                        {
                            for (int j = 0; j < iNumLetters; j++)
                            {
                                if (pCypherText[i] == pAlphabet[j])
                                {
                                    indAlph = j;
                                    break;
                                }
                            }
                            int ind = (indAlph - (int)tmpKey[(i % KEY_LEN)]);
                            while (ind < 0)
                            {
                                ind += iNumLetters;
                            }
                            printf("%c", pAlphabet[(ind % iNumLetters)]);
                        }
                        printf("\n. . . . . . . . . \n\nContinue? (No - Esc / Print decrypted text to file - Enter / Yes - anything)\n> ");
                        int x = _getch();
                        if (x == 27)
                        {
                            goto cleanup;
                        }
                        else if (x == 13)
                        {
                            FILE* outF = fopen("freqAnalysisRes.txt", "wb");
                            indAlph = 0;
                            for (int i = 0; i < iLen; i++)
                            {
                                for (int j = 0; j < iNumLetters; j++)
                                {
                                    if (pCypherText[i] == pAlphabet[j])
                                    {
                                        indAlph = j;
                                        break;
                                    }
                                }
                                int ind = (indAlph - (int)tmpKey[(i % KEY_LEN)]);
                                while (ind < 0)
                                {
                                    ind += iNumLetters;
                                }
                                fwrite(&pAlphabet[(ind % iNumLetters)], 1, 1, outF);
                            }
                        }

                    }
                }
            }
        }
    }

    

    printf("\n- - -Nothing to be done - - -\n\n");

    //
    // cleanup
    //
    cleanup:
    for (int i = 0; i < KEY_LEN; i++)
    {
        free(cypheredParts[i]);
        free(rezSymbols[i]);
    }
    free(cypheredParts);
    free(numOfSymbolsInPart);
    free(tmpFreqOfBlockSymbols);
    free(tmpFreqOfBlockSymbolsSorted);
    free(rezSymbols);
    free(numOfRezSymbols);
    //free(mostFrequentSymbols);
}
*/

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
        int* permutation =  NULL;  // 
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

float findMIc(uint8_t* pPart1, int iPartNum1, uint8_t* pPart2, int iPartNum2)
{

}

void friedman2(uint8_t* pCypherText, int iLen, uint8_t* pAlphabet, int iNumLetters, int* pQuantities)
{
    int numOfParts;
    float delta;
    float sum;
    float* frequencies = NULL;
    float* MICs = NULL;
    float* subs = NULL;

    frequencies = (float*)malloc(iNumLetters * sizeof(float));
    MICs = (float*)malloc(iLen * sizeof(float));
    subs = (float*)calloc(iLen, sizeof(float));

    for (int i = 0; i < iNumLetters; i++)
    {
        frequencies[i] = (float)pQuantities[i] / (float)iLen;
    }

    delta = findDelta(iNumLetters, frequencies);

    printf("Delta: %f\n\n", delta);
    for (int d = 2; d < ((iLen < MAX_KEY_LEN) ? iLen : MAX_KEY_LEN); d++)
    {
        printf("Supposed key length: %8i\t\t", d);

        for (int j = 0; j < d; j++)
        {
            //MICs[j] = findMIc(pCypherText, iNumLetters, j, d);
        }

        printf("\n");
    }

    //float minSub = delta - subs[2];
    //if (minSub < 0.0)
    //{
    //    minSub *= -1.0;
    //}
    //int minInd = 2;
    //for (int d = 3; d < ((iLen < MAX_KEY_LEN) ? iLen : MAX_KEY_LEN); d++)
    //{
    //    if (subs[d] < 0.0)
    //    {
    //        subs[d] *= -1.0;
    //    }
    //
    //    if (subs[d] < minSub)
    //    {
    //        minSub = subs[d];
    //        minInd = d;
    //    }
    //}

    //printf("Appropriate key length: %i\nDelta = %.6f\nDifference = %.6f\n", minInd, delta, minSub);

    // 
    // cleanup:
    // 
    free(frequencies);
    free(MICs);
    free(subs);
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
    friedman1(pCypherTextV, iLenTextV, pAlphabet, iNumAlphabetLetters, pQuantLetters);




    free(pAlphabet);
    free(pCypherTextV);
    free(pCypherTextA);
    free(pQuantLetters);

    system("pause");
    return 0;
}

