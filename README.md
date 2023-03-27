# NMF
# simple NMF "hello world" in C
Non-negative matrix factorization (NMF or NNMF)

credits https://github.com/0xfffffff7/NMF

compile/execute with <code>$ gcc nnmf.c -Wall -lm -o nnmf && ./nnmf</code>

output example

```
V =
        50.0         50.0         50.0         50.0 
       100.0         20.0         30.0         90.0 
       100.0         50.0         40.0        100.0 
        90.0         10.0         10.0         80.0 
        10.0        100.0         90.0         20.0 

W x H =
        50.2         51.8         48.2         49.7 
        98.5         22.9         26.2         92.0 
       102.3         45.2         45.8         97.0 
        88.1          9.1         13.5         81.5 
        15.1        100.5         88.4         20.8 

W =
       864.1        560.0 
       320.0       1235.1 
       706.1       1246.1 
        86.7       1125.1 
      1744.8         20.8 

H * 1000 =
         7.7         57.6         50.5         11.1 
        77.7          3.6          8.1         71.6 
```
