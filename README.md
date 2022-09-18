# Star-Alignment
Multiple sequence alignment using star-alignment algorithm with improving the score of alignment.


### Some Examples :

Input Format :
``` python
>Number of sequences
>sequence1
>sequence2
>    .
>    .
>    .
>sequenceN

```

Output Format :
``` python
Alignment Score : int
Alignment :
    sequence1
    sequence2
        .
        .
        .
    sequenceN
```

---
Input 1 :
```python
4 
TYIMREAQYESAQ
TCIVMREAYE
YIMQEVQQER
WRYIAMREQYES
```
Output 1 :
```python
51
-TYI-MREAQYESAQ
-TCIVMREA-YE---
--YI-MQEVQQER--
WRYIAMRE-QYES--
```
---

---
Input 2 :
``` python
5
TAGCTACCAGGA
CAGCTACCAGG
TAGCTACCAGT
CAGCTATCGCGGC
CAGCTACCAGGA
```
Output 2 :
``` python
240
TAGCTA-C-CAGGA
CAGCTA-C-CAGG-
TAGCTA-C-CA-GT
CAGCTATCGC-GGC
CAGCTA-C-CAGGA
```
---
