#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdint.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h> 

#define KEY_LEN 4
#define PRINT_CONST 30

int compInt(const int* a, const int* b)
{
    return *b - *a;
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

int readAlphabet(uint8_t** pAlphabet, int** pQuant)
{
    FILE* A = fopen("alphabet.txt", "rb");
    FILE* F = fopen("frequencies.txt", "r");

    char* pBuf = NULL;
    int iNum = 0;
    uint8_t c;

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

    for (int i = 0; i < KEY_LEN; i++)
    {
        // divide cyphertext to blocks of bytes which have been encrypted on the same key's byte
        for (int j = i; j < iLen; j += KEY_LEN, numOfSymbolsInPart[i]++)
        {
            cypheredParts[i][numOfSymbolsInPart[i]] = pCypherText[j];
        }

        //
        // find the most frequent symbol(s) in the block
        //
        
        memset(tmpFreqOfBlockSymbols, 0, iNumLetters * sizeof(int));
        for (int j = 0; j < numOfSymbolsInPart[i]; j++)
        {
            tmpFreqOfBlockSymbols[findIndexOfSymbol(cypheredParts[i][j], pAlphabet, iNumLetters)]++;
        }

        // find the maximal quantity
        int maxQuant = 1;
        memcpy(tmpFreqOfBlockSymbolsSorted, tmpFreqOfBlockSymbols, iNumLetters * sizeof(int));
        qsort(tmpFreqOfBlockSymbolsSorted, iNumLetters, sizeof(int), (int(*) (const void*, const void*)) compInt);
        while (tmpFreqOfBlockSymbolsSorted[maxQuant - 1] == tmpFreqOfBlockSymbolsSorted[maxQuant] && maxQuant < iNumLetters)
        {
            maxQuant++;
        }

        // find all the symbols with maximal quantity
        for (int j = 0; j < iNumLetters; j++)
        {
            int tmpCntr = -1;
            int tmpInd = 0;
            for (int k = 0; k < iNumLetters; k++)
            {
                if (tmpFreqOfBlockSymbols[k] == 0)
                {
                    continue;
                }
                if (tmpFreqOfBlockSymbols[j] == tmpFreqOfBlockSymbols[k])
                {
                    tmpCntr++;
                    tmpInd = j;
                }
            }
            if (tmpCntr == maxQuant - 1)
            {
                if (rezSymbols[i] == NULL)
                {
                    rezSymbols[i] = (uint8_t*)malloc(1);
                }
                else
                {
                    rezSymbols[i] = (uint8_t*)realloc(rezSymbols[i], numOfRezSymbols[i] + 1);
                }
                rezSymbols[i][numOfRezSymbols[i]] = pAlphabet[tmpInd];
                numOfRezSymbols[i]++;
            }
        }
    }

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
    char userInput[32];
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
    /*
    for (int n = 0; ; n++)
    {
        uint8_t tmpKey[KEY_LEN] = { 0 };
        for (int i = 0; i < KEY_LEN; i++)
        {
            fixed[i] = false;
        }

        fixed[(n % 3)] = true;

        printf("- - - Finding key - - -\n\n");

        for (int i = 0; i < KEY_LEN; i++)
        {
            int tmpSymbX = findIndexOfSymbol(rezSymbols[i][indexesInRezSymbols[0]],           pAlphabet, iNumLetters);
            int tmpSymbY = findIndexOfSymbol(mostFrequentSymbols[indexInMostFrequentSymbols], pAlphabet, iNumLetters);
            printf("\nSearching key byte N%i:\n\n", i);
            if (numOfRezSymbols[i] == 1)
            {
                if (numOfMostFrequent == 1) // case 0 - can match automaticaly
                {
                    tmpKey[i] = (tmpSymbX - tmpSymbY + iNumLetters) % iNumLetters;
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
                    tmpKey[i] = (tmpSymbX - tmpSymbY + iNumLetters) % iNumLetters;
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
                    tmpKey[i] = (tmpSymbX - tmpSymbY + iNumLetters) % iNumLetters;
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
                    tmpKey[i] = (tmpSymbX - tmpSymbY + iNumLetters) % iNumLetters;
                    printf("\nVariative match: %c (%u) -> %u => K[%i] = %u\n",
                        rezSymbols[i][indexesInRezSymbols[i]],
                        rezSymbols[i][indexesInRezSymbols[i]],
                        mostFrequentSymbols[indexInMostFrequentSymbols],
                        i,
                        tmpKey[i]);
                }
            }
        }
        printf("Key: %u %u %u %u\nDecrypted text:\n", tmpKey[0], tmpKey[1], tmpKey[2], tmpKey[3]);
        int indAlph = 0;
        for (int i = 0; i < iLen; i++)
        {
            for (int j = 0; j < iNumLetters; j++)
            {
                if (pAlphabet[i] == pCypherText[j])
                {
                    indAlph = j;
                    break;
                }
            }
            printf("%c", pAlphabet[(tmpKey[i % KEY_LEN] - indAlph + iNumLetters) % iNumLetters]);
        }
        printf("\nContinue? (y/n)\n> ");
        scanf("%s", &userInput);
        if (userInput[0] != 'y')
        {
            goto cleanup;
        }

        // prepare indexes for next iteration
        indexInMostFrequentSymbols = (indexInMostFrequentSymbols + 1) % numOfMostFrequent;
        for (int i = 0; i < KEY_LEN; i++)
        {
            if (fixed[i])
            {
                indexesInRezSymbols[i] = (indexesInRezSymbols[i] + 1) % numOfRezSymbols[i];
            }
        }
    }*/
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
    free(mostFrequentSymbols);
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

    freqAnalysis(pCypherTextV, iLenTextV, pAlphabet, iNumAlphabetLetters, pQuantLetters);





    free(pAlphabet);
    free(pCypherTextV);
    free(pCypherTextA);
    free(pQuantLetters);

    system("pause");
    return 0;
}

