# BinaryGoppaCode

## GF2 operation
+ Addition
+ Multiplication
+ Square
+ Division
+ GCD
----

----
## 다항식 곱셈 검증을 위한 방법

SageMath 프로그램을 이용하여 검증한다.

설치하지 않아도 온라인으로 이용할 수 있다.

https://sagecell.sagemath.org/


    bin.<z> = GF(2)[]
    poly = z^13+z^4+z^3+z+1
    S = bin.quotient(poly, 'X')

    #example
    A = z^5+z^4+z^3+z^2+z^0
    B = z^9+z^8+z^7+z^6+z^3+z^0
    C = z^12+z^10+z^9+z^8+z^5+z^1+z^0
    S(A*B) == S(C)


----
### 
나중에 openssl bn 사용 예정