<br><br><br><br><br><br><br>
<center>
<font face="黑体" size=10>
计算机科学和机器学习中的
<br>
代数学、拓扑学、微积分以及
</br>
最优化理论
</font>
</center>

<br><br><br><br><br>

<center>
<big><big><big><br>作者：<u>Zhenyu Wu</u></big></big></big>
<big><big><big><br>版本号：<u>V0.0.1</u></big></big></big>
<big><big><big><br>时间：<u>2020年02月22日</u></big></big></big>

<big><big><big><br>e-mail：<u>807698462@qq.com</u></big></big></big>
<big><big><big><br>GitHub：<u>https://github.com/wzy6642/Math-in-Machine-Learning-and-Deep-Learning</u></big></big></big>
</center>
<br><br><br><br><br><br><br>

# 第一章
<big><big><big>
群、环、域
</big></big></big>


&emsp;&emsp;在接下来的四章中，我们将介绍最基本的四种代数结构（群、环、域以及向量空间，这些都是近世代数的内容），其中着重介绍的是向量空间。同时呢，也会对线性代数的基本概念进行介绍，其中包括向量空间、子空间、线性组合、线性无关、基、商空间、线性映射、矩阵、基变换、直和、线性形式、对偶空间、超平面、线性变换等。

&emsp;&emsp;我们首先引入笛卡尔积的概念：假设我们有集合 $A=\{a_1,a_2,a_3\}, B=\{b_1,b_2,b_3\}$，那么 $A$ 与 $B$ 的笛卡尔积记作 $A \circ B$，计算结果为有序数对，即 $A \circ B=\{(a,b)|a \in A,b \in B\}$ 。不难发现当 $A,B$ 的元素均为实数时，笛卡尔积表示平面直角坐标系。当然我们也可以定义某一种具体的运算方式，例如 $A+B$ （除此之外，我们还可以定义 "$-$" 以及 "$\times$" 运算），此时 $A+B=\{a+b|a \in A,b \in B\}$。当然，集合的元素不仅仅局限在实数范围内，更一般的，我们可以用字符串作为其元素，例如：$A=\{东，西\}，B=\{南，北\}$ 那么 $A \circ B = \{(东，南)，(东，北)，(西，南)，(西，北)\}$。那么经过笛卡尔积的运结果包含多少个元素呢？不难发现，当被作用的集合均为有限集时，最终计算结果的元素个数就是各个集合元素个数的乘积。

## 1.1 群、子群、陪集
&emsp;&emsp;实数组成的集合$R$有两种运算操作：加法 $+$：$R + R \rightarrow R$ ，以及乘法 $\times$：$R \times R \rightarrow R$，实数集之间的加法和乘法运算本质上是一个阿贝尔群。接下来我们回忆一下群的定义。

**定义1.1：** 两个集合 $G$ 通过二元运算符操作 $\cdot$ 便可得到一个群<font color=#ea4335>（注：这里的 $\cdot$ 可取 $\times$ 或 $+$，下同）</font>，例如：$G \times G \rightarrow G$ <font color=#ea4335>（当且仅当<u>三个 $G$ 都相同</u>的情况下才叫二元运算，且本章仅对二元运算进行讨论）</font>，该运算操作可以将集合中的每一组元素 $a,b\in G$ 进行有效的结合，从而得到 $a\cdot b\in G$。我们先对一些计算的结果是否为二元运算进行讨论：
> 1、我们有 $Z \cdot Z \rightarrow Z$，其中 $Z$ 表示整数集
>> (1)、当 $\cdot$ 表示 $+$ 时，是二元运算；<br>
>> (2)、当 $\cdot$ 表示 $-$ 时，是二元运算；<br>
>> (3)、当 $\cdot$ 表示 $\times$ 时，是二元运算；<br>
>> (4)、当 $\cdot$ 表示 $\div$ 时，不是二元运算，因为分母不能为0，且两个整数相除的计算结果不一定为整数。

> 2、我们有 $Z \cdot Z \rightarrow Z$，其中 $Z=Z-\{0\}$ 表示非零整数集
>> (1)、当 $\cdot$ 表示 $+$ 时，不是二元运算，因为相反数的和为0；<br>
>> (2)、当 $\cdot$ 表示 $-$ 时，不是二元运算，因为相同数的差为0；<br>
>> (3)、当 $\cdot$ 表示 $\times$ 时，是二元运算；<br>
>> (4)、当 $\cdot$ 表示 $\div$ 时，不是二元运算，因为商不一定为整数。

> 3、我们有 $R^+ \cdot R^+ \rightarrow R^+$，其中 $R^+$ 表示正实数集
>> (1)、当 $\cdot$ 表示 $+$ 时，是二元运算；<br>
>> (2)、当 $\cdot$ 表示 $-$ 时，不是二元运算，因为计算结果可能为0或者负数；<br>
>> (3)、当 $\cdot$ 表示 $\times$ 时，是二元运算；<br>
>> (4)、当 $\cdot$ 表示 $\div$ 时，是二元运算。

> 4、我们有 $Q \cdot Q \rightarrow Q$，其中 $Q=Q-\{0\}$ 表示非零有理数集
>> (1)、当 $\cdot$ 表示 $+$ 时，不是二元运算，因为相反数的和为0；<br>
>> (2)、当 $\cdot$ 表示 $-$ 时，不是二元运算，因为相同数的差为0；<br>
>> (3)、当 $\cdot$ 表示 $\times$ 时，是二元运算；<br>
>> (4)、当 $\cdot$ 表示 $\div$ 时，是二元运算。

> 5、我们有 $V \cdot V \rightarrow V$，其中 $V$ 表示元素为实数的向量，即 $V=\{a_i|a_i \in R\}$，由于向量加法的运算方式为：$\{a_1,a_2,a_3\}+\{b_1,b_2,b_3\}=\{a_1+b_1,a_2+b_2,a_3+b_3\}$，即对应位置元素相加，不难发现其计算结果的元素仍然为实数，所以该运算过程为二元运算。

> 6、我们有 $M_n \cdot M_n \rightarrow M_n$，其中 $M_n$ 表示实数组成的 $n$ 阶方阵，不难发现 $n$ 阶实数方阵相加的结果仍为 $n$ 阶实数方阵， $n$ 阶实数方阵相减的结果仍为 $n$ 阶实数方阵， $n$ 阶实数方阵相乘的结果仍为 $n$ 阶实数方阵，所以上述运算均为二元运算。

&emsp;&emsp;我们通常将二元运算称为 $\cdot$ 运算在 $A$ 上封闭。那么如果 $A$ 是有限集如何判断封闭性呢？其实只需要判断是否满足下图即可：

<div align=center><img src="figure/有限封闭.jpg" width="40%"/></div>
<center>
图1.1 判断有限集是否封闭
</center>

左上角表示二元运算符，第一行表示集合 $A$，第一列也表示集合 $A$，剩余部分表示二元运算结果，可见运算结果和集合包含的元素是一致的。所以，判断是否封闭，我们只需要判断各部分间元素是否均相同即可。

&emsp;&emsp;群具有四个性质：封闭性（$\forall a,b \in G \Rightarrow a \cdot b \in G$）、结合律、单位元素、逆元素。我们定义 $a,b,c$ 为集合 $G$ 中的元素，即 $a,b,c\in G$，$e$ 为集合 $G$ 的单位元素，且 $e\in G$ ，集合 $G$ 中的每一个元素均可逆，那么有下述等式恒成立（必须满足二元运算）：

（等式1：结合律）$a \cdot(b \cdot c)=(a \cdot b) \cdot c$

（等式2：单位元）$a \cdot e=e \cdot a=a$

（等式3：逆元素）对于集合$G$中的每一个元素$a$，$a \in G$，都有$a^{-1} \in G$，满足$a \cdot a^{-1}=a^{-1} \cdot a=e$

&emsp;&emsp;交换律：对于群 $G$ 中任意的两个元素 $a,b$ ，即 $a,b \in G$，若 $a \cdot b=b \cdot a$，那么我们称群 $G$ 是可交换的，也称之为<font color=#ea4335>阿贝尔群</font>。若图1.1中蓝色部分组成的方阵中的元素关于主对角线对称，即 $a_{ij}=a_{ji}$，此时就满足交换律。

&emsp;&emsp;两个实数集 $M$ 使用二元运算符 $\cdot$ 进行运算：$M \times M \rightarrow M$，若元素仅仅满足结合律（<u>当元素个数较少时，我们可以遍历每一种情况判断是否相等</u>）和单位元，那么我们将其称作<font color=#ea4335>独异点</font>。例如，两个由自然数组成的集合 $N=\{0,1,\cdots,n,\cdots \}$ 进行加法运算，构成可交换独异点，但由于其不满足等式3，所以不能被称为群。我们接下来给出几个群的例子供大家学习。

<div align=center><img src="figure/群的判断.jpg" width="100%"/></div>
<center>
图1.2 群的判断
</center>

**示例1.1：**

&emsp;&emsp;1. 将两个集合 $Z=\{\cdots,-n,\cdots,-1,0,1,\cdots,n,\cdots\}$ 进行相加便构成一个单位元为0的阿贝尔群。但是两个 $Z^*$ 进行相乘并不能得到群，其中 $Z^*=Z-\{0\}$ 。分析过程如下：

（1） 相加

&emsp;&emsp;a、封闭性：$\forall a,b \in Z$，我们都有 $a+b \in Z$，所以满足封闭性。

&emsp;&emsp;b、结合律：$\forall a,b,c \in Z$，我们都有 $(a+b)+c=a+(b+c)$，所以满足结合律。

&emsp;&emsp;c、单位元：$\forall a \in Z$，我们都有 $a+0=0+a=a$，所以满足单位元。

&emsp;&emsp;d、逆元素：$\forall a \in Z$，我们都有 $a+b=0$，且 $b \in Z$，所以满足逆元素。

&emsp;&emsp;e、交换律：$\forall a,b \in Z$，我们都有 $a+b=b+a$，所以满足交换律。

&emsp;&emsp;综上，$Z$ 是单位元为0的阿贝尔群。

（2） 相乘

&emsp;&emsp;a、封闭性：$\forall a,b \in Z^*$，我们都有 $a \times b \in Z^*$，所以满足封闭性。

&emsp;&emsp;b、结合律：$\forall a,b,c \in Z^*$，我们都有 $(a \times b) \times c=a \times (b \times c)$，所以满足结合律。

&emsp;&emsp;c、单位元：$\forall a \in Z^*$，我们都有 $a \times 1=1 \times a=a$，所以满足单位元。

&emsp;&emsp;d、逆元素：$\forall a \in Z^*$，我们不一定有 $b \times a=1,b \in Z$，所以不满足逆元素。

&emsp;&emsp;综上，$Z$ 不是群，仅为独异点。需要注意的是 $Z,Z^*$ 的元素均为无限的。

&emsp;&emsp;2. 通过将两个有理数组成的集合 $Q$（集合中的元素均可写为 $p/q$ 的形式，其中 $p,q\in Z$ 且 $q\neq0$）进行相加，可以得到一个单位元为0的阿贝尔群。将两个集合 $Q^*$ 进行相乘也可得到一个单位元为1的阿贝尔群，其中 $Q^*=Q-\{0\}$。分析过程如下：

（1） 相加

&emsp;&emsp;a、封闭性：$\forall a,b \in Q$，我们都有 $a+b \in Q$，所以满足封闭性。

&emsp;&emsp;b、结合律：$\forall a,b,c \in Q$，我们都有 $(a+b)+c=a+(b+c)$，所以满足结合律。

&emsp;&emsp;c、单位元：$\forall a \in Q$，我们都有 $a+0=0+a=a$，所以满足单位元。

&emsp;&emsp;d、逆元素：$\forall a \in Q$，我们都有 $a+b=0,b \in Q$，所以满足逆元素。

&emsp;&emsp;e、交换律：$a,b \in Q$，我们都有 $a+b=b+a$，所以满足交换律。

&emsp;&emsp;综上，$Q$ 是单位元为0的阿贝尔群。

（2） 相乘

&emsp;&emsp;a、封闭性：$\forall a,b \in Q^*$，我们都有 $a \times b \in Q^*$，所以满足封闭性。

&emsp;&emsp;b、结合律：$\forall a,b,c \in Q^*$，我们都有 $(a \times b) \times c=a \times (b \times c)$，所以满足结合律。

&emsp;&emsp;c、单位元：$\forall a \in Q^*$，我们都有 $a \times 1=1 \times a=a$，所以满足单位元。

&emsp;&emsp;d、逆元素：$\forall a \in Q^*$，我们都有 $a \times b=1,b \in Q^*$，所以满足逆元素。

&emsp;&emsp;e、交换律：$a,b \in Q^*$，我们都有 $a \times b=b \times a$，所以满足交换律。

&emsp;&emsp;综上，$Q^*$ 是单位元为1的阿贝尔群。需要注意的是 $Q,Q^*$ 的元素均为无限的。

&emsp;&emsp;我们发现上述的群都满足交换律，那么有没有不满足交换律的群呢？我们在这里给出一个例子加以说明。 $M^+$ 为由无穷个 $n \times n$ 的可逆方阵（行列式不为0，且元素均为实数）构成的集合，那么 $M^+ \times M^+$ 是不是群呢?分析过程如下：

&emsp;&emsp;a、封闭性：我们取 $\forall A,B \in M^+$，都有 $A \times B \in M^+$，所以满足封闭性。

&emsp;&emsp;b、结合律：$\forall A,B,C \in M^+$，都有 $(A \times B) \times C=A \times (B \times C)$，所以满足结合律（矩阵乘法满足结合律）。

&emsp;&emsp;c、单位元：$\forall A \in M^+$，都有 $A \times I_n = I_n \times A=A$，所以满足单位元，其中 $I_n$ 为 $n$ 阶单位阵。

&emsp;&emsp;d、逆元素：$\forall A \in M^+$，都有 $A \times B=I_n$，其中 $B=A^{-1}$，所以满足逆元素。

&emsp;&emsp;e、交换律：$\forall A,B \in M^+$，一般而言 $A \times B \neq B \times A$，所以不满足交换律。

&emsp;&emsp;综上，$M^+$ 是一个单位元为 $I_n$ 的群。

&emsp;&emsp;3. 给定一个非空集合 $S$，若有作用方式 $f$ 可以使两个相同集合 $S$ 之间满足双射关系（也可称为 $S$ 的排列），即 $f:S\rightarrow S$，此时，通过函数与函数之间的运算便可构成一个群（例如，将函数 $f$ 和函数 $g$ 通过复合运算得到计算结果 $f \circ g$，其中 $f$ 和 $g$ 均可使集合 $S$ 到其自身之间满足一一映射），当集合 $S$ 中的元素个数超过两个时，所构成的群并不是一个阿贝尔群。集合 $S=\{1,\cdots,n\}$ 所构成的置换群通常被记作 $S_n$，也被称为 $n$ 个元素构成的对称群。举例如下：

设 $S=\{1, 2, 3\}$，我们定义 $f$ 的映射关系为：
$$
f:S \rightarrow S,f=\{(1,2),(2,3),(3,1)\} \\
f^{-1}:S \rightarrow S,f^{-1}=\{(1,3),(2,1),(3,2)\}
$$
同时定义 $g$ 的映射关系为：
$$
g:S\rightarrow S,g=\{(1,3),(2,1),(3,2)\} \\
g^{-1}:S\rightarrow S,g^{-1}=\{(1,2),(2,3),(3,1)\}
$$
$h$ 的映射关系为：
$$
h:S\rightarrow S,h=\{(1,2),(2,1),(3,3)\} \\
h^{-1}:S\rightarrow S,h^{-1}=\{(1,2),(2,1),(3,3)\}
$$

则，我们有 $f \circ g$ 的计算结果：
$$
f \circ g:S\rightarrow S,f \circ g=\{(1,1),(2,2),(3,3)\}
$$
可知，$f \circ g$ 的复合结果依然为满射，且其逆过程为：
$$
(f \circ g)^{-1}:S\rightarrow S,(f \circ g)^{-1}=\{(1,1),(2,2),(3,3)\}
$$
接下来我们利用群的相关定义进行验证，其中我们将 $f \circ g \circ h$ 记作 $\varphi$，那么其逆过程记为 $\varphi^{-1}$：

&emsp;&emsp;封闭性：$\forall \ a \in S$，都有 $\varphi(a) \in S$，所以满足封闭性；

&emsp;&emsp;结合律：$\forall \ a \in S$，都有 $(f \circ g) \circ h=f \circ (g \circ h)$，所以满足结合律；

&emsp;&emsp;单位元：不难发现 $f \circ g$ 构成了恒同映射，即 $\forall a \in S,f \circ g(a)=a$，所以 $\forall a \in S,f \circ g \circ h(a)=h \circ f \circ g(a)=h(a)$；

&emsp;&emsp;逆映射：对于每一个$a \in S$，均有 $f \circ g \circ h \circ (f \circ g \circ h(a))^{-1}=(f \circ g \circ h)^{-1} \circ f \circ g \circ h(a)=f \circ g(a)$

&emsp;&emsp;故而 $f \circ g \circ h(a)$ 的复合结果为群，且为对称群。这部分知识来自于[知乎专栏](https://zhuanlan.zhihu.com/p/85203425)。

&emsp;&emsp;4. 对于任意的正整数 $p \in N$，定义在 $Z$ 上的同余关系记作 $m \equiv n \ (mod\ p)$，具体定义如下：

$$
m \equiv n \ (mod\ p) \ \Leftrightarrow \ m-n=kp \ ,存在k \in Z
$$
其中，$\equiv$ 表示同余符号，即 $m \ mod \ p \equiv n \ mod \ p$ ，读者很容易证明这是一个恒等关系，此外，将同余号两边同时进行相加或相乘，相等关系不变，即若 $m_1 \equiv n_1 \ (mod\ p)$ 且 $m_2 \equiv n_2 \ (mod\ p)$ ，则 $m_1+m_2 \equiv n_1+n_2 \ (mod\ p)$ ， $m_1m_2 \equiv n_1n_2 \ (mod\ p)$ 。我们在这里给出一个算例：
$$
5 \equiv 3 \ (mod \ 2)\\
11 \equiv 7 \ (mod \ 2)\\ \Rightarrow
加运算：16 \equiv 10 \ (mod \ 2)\\\Rightarrow
乘运算：55 \equiv 21 \ (mod \ 2)
$$
我们将一组等价类对 $p$ 取余的相加和相乘操作用如下记号进行表示：
$$
[m]+[n]=[m+n]\\
[m] \cdot [n]=[mn]
$$
读者很容易证明将一组对 $p$ 取余的同余类进行相加可以得到单位元为 $[0]$ 的阿贝尔群，我们将这个群记作 $Z/pZ$ 。分析过程如下：

根据上面的定义，我们可以用 $[n]$ 表示对 $p$ 取余为 $n$ 的同余类，那么该集合可以表示为 $[n]=\{a\ | \ a=kp+n,k \in Z\}$。

&emsp;&emsp;封闭性：$\forall \ [a],[b] \in Z$，都有 $[a+b] \in Z$，所以满足封闭性；

&emsp;&emsp;结合律：$\forall \ [a],[b],[c] \in Z$，都有 $[a+b]+[c]=[a]+[b+c]$，所以满足结合律；

&emsp;&emsp;单位元：$\forall \ [a] \in Z$，都有 $[a]+[0]=[0]+[a]=[a]$，所以满足单位元；

&emsp;&emsp;逆元素：$\forall \ [a] \in Z$，都有 $[a]+[-a]=[-a]+[a]=[0]$，且 $[-a] \in Z$，所以满足逆元素；

&emsp;&emsp;交换律：$\forall \ [a],[b] \in Z$，都有 $[a]+[b]=[b]+[a]$，所以满足交换律。

&emsp;&emsp;综上，$Z/pZ=\{[0],[1],[2],\cdots,[p-1]\}$ 是一个单位元为 $[0]$ 的阿贝尔群。

&emsp;&emsp;5. 将一组系数为实数或复数的 $n \times n$ 的可逆矩阵进行相乘可以得到一个单位元为单位矩阵 $I_n$ 的群，这个群被称为一般线性群，对于系数为实数的记作 $GL(n,R)$ 对于系数为复数的记作 $GL(n,C)$ 。假设我们有可逆矩阵 $A,B$ 以及实数 $\lambda,\beta$，则有如下计算过程（和2相结合分析）：
$$
结合律：(\lambda A \times \beta )B = \lambda (A \times \beta B)\\
单位元：\lambda A \times \beta B \times I = I \times \lambda A \times \beta B = \lambda A \times \beta B\\
逆元素：(\lambda A \times \beta B)^{-1} \times (\lambda A \times \beta B) = I
$$
注：$\because$ 矩阵可逆的充要条件之一是它的行列式不等于0，$\therefore$ 两个可逆矩阵相乘得到矩阵仍然是可逆矩阵。

&emsp;&emsp;6. 将一组系数为实数或复数的 $n \times n$ 的可逆矩阵 $A$ 进行相乘，其中矩阵的行列式为1，即 $det(A)=1$ ，可以得到一个单位元为单位矩阵 $I_n$ 的群，这个群被称为特殊线性群，对于系数为实数的记作 $SL(n,R)$ 对于系数为复数的记作 $SL(n,C)$ 。证明方式如5。

&emsp;&emsp;7. 将一组系数为实数的 $n \times n$ 的矩阵 $Q$ 进行相乘，可以得到一个单位元为单位矩阵 $I_n$ 的群，其中矩阵 $Q$ 满足 $QQ^T=Q^TQ=I_n$ 。我们有 $Q^{-1}=Q^T$ ，这个群被称为正交群，记作 $O(n)$。

&emsp;&emsp;8. 将一组系数为实数的 $n \times n$ 的矩阵 $Q$ 进行相乘，可以得到一个单位元为单位矩阵 $I_n$ 的群，其中矩阵 $Q$ 满足 $QQ^T=Q^TQ=I_n \ 且 \ det(Q)=1$ 。就像示例7中一样，我们有 $Q^{-1}=Q^T$  ，这个群被称为特殊正交群或旋转群，记作 $SO(n)$。

&emsp;&emsp;在示例5~8中，除了 $SO(2)$ 为阿贝尔群以外，当 $n \geq 2$ 时均为非阿贝尔群。我们通常将集合相加后得到的阿贝尔群用 $G$ 进行表示，此时元素 $a \in G$ 的逆元 $a^{-1}$ 可以表示为 $-a$ 。群的单位元（幺元）是独一无二的，我们可以得到一些更一般的结论。

**命题1.1：** 若有二元运算符 $\cdot$ : 使得 $M \times M \rightarrow M$ 的计算结果是一个群，且 $e' \in M$ 是左单位元，$e'' \in M$ 是右单位元，也即：
$$
G2l:\ \ 对于任意的\  a \in M \ \ 都有 \ \ e' \cdot a = a \\
G2r:\ \ 对于任意的\  a \in M \ \ 都有 \ \ a \cdot e'' = a
$$
那么我们有 $e'=e''$ 。

证明过程如下：若我们令等式 $G2l$ 中 $a=e''$，那么我们有：
$$
e' \cdot e'' = e''
$$
若我们令等式 $G2r$ 中 $a=e'$，我们有：
$$
e' \cdot e'' = e'
$$
那么，我们得到如下等式：
$$
e' = e' \cdot e'' = e''
$$
综上，我们便得到了 $e'=e''$ （存在且唯一）。

&emsp;&emsp;**命题1.1**说明了独异点的幺元是唯一的，而所有的群都是独异点，所以群的幺元都是唯一的。此外，群中的每一个元素都有其对应的逆元，接下来我们给出一个命题：

**命题1.2：** 在独异点 $M$ 中有幺元 $e$ ，若某元素 $a \in M$ 有左逆元 $a' \in M$ 和右逆元 $a'' \in M$，也即：
$$
G3l:\ \ a' \cdot a = e\\
G3r:\ \ a \cdot a'' = e
$$
则有 $a'=a''$。

证明过程如下：结合公式 $G3l$ 以及 $e$ 为群的幺元，我们可以得到
$$
(a' \cdot a) \cdot a'' = e \cdot a'' = a''
$$
同样的，结合公式 $G3r$ 以及 $e$ 为群的幺元，我们可以得到
$$
a' \cdot (a \cdot a'') = a' \cdot e = a'
$$
由于 $M$ 是独异点，所以二元运算符 $\cdot$ 符合结合律，故而有
$$
a' = a' \cdot (a \cdot a'') = (a' \cdot a) \cdot a'' = a''
$$
得证 $a'=a''$ （存在且唯一）。

**注意：** 群的单位元（等式2）以及群的逆元素（等式3）的证明可以被弱化为仅要求 $G2r$ （右单位元存在）和 $G3r$ （对于群中每一个元素均存在右逆元）存在（或者是 $G2l$ 和 $G3l$ 存在）。通过证明 $G2l$ 和 $G3l$ 成立来证明等式2（公理2）以及等式3（公理3）成立是一个行之有效的方法。

**定义1.2：** 若群 $G$ 由有限的 $n$ 个元素组成，我们称群 $G$ 为 $n$ 阶群。若群 $G$ 的元素个数是无穷的，我们称群 $G$ 为无穷阶群。若群为有限群，那么其阶数我们使用符号 $|G|$ 进行表示。除此之外，我们还可对群中某个元素的阶进行分析。在一个群 $G$ 中，使得满足 $a \cdot a \cdot \ \cdots\  \cdot a=e$ （共计 $n$ 个元素 $a$ 做二元运算）的最小正整数 $n$ 叫做 $a$ 的阶，若这样的 $n$ 不存在，称 $a$ 是无穷阶，或者叫 $a$ 的阶是无穷。

&emsp;&emsp;例如 $G=Z$，我们有二元运算 $+$，使得 $G + G \rightarrow G$，不难发现 $G$ 是幺元为0的阿贝尔群。那么，当 $a=e$ 时满足 $e=e$，所以 $e$ 的阶为1，对于其它 $a \neq e$ 的元素呢，显然不论多少个 $a$ 进行相加都得不到 $e$，所以我们说这些元素阶都是无穷的。

&emsp;&emsp;例如 $G=Q-\{0\}$，我们有二元运算 $\times$，使得 $G \times G \rightarrow G$，不难发现 $G$ 是幺元为1的阿贝尔群。那么，当 $a=e$ 时满足 $e=e$，所以 $e$ 的阶为1，对于其它 $a \neq e$ 的元素呢，显然-1的偶次幂为1，所以其阶为2，除了1、-1之外，不论多少个 $a$ 进行相乘都得不到 $e$，所以我们说这些元素阶都是无穷的。

&emsp;&emsp;我们目前接触到的群都为无穷阶群，那么有没有有限群呢？我们在这里给出一个例子：集合 $G=\{1,\frac{-1+ \sqrt{3} \imath}{2},\frac{-1- \sqrt{3} \imath}{2}\}$，我们对其做二元运算 $G \times G \rightarrow G$，可以得到单位元为1的阿贝尔群，验证过程如下：

&emsp;&emsp;a、封闭性：1与任何数相乘都为其自身，所以我们只需要验证两个不为1的元素的乘积结果是否是 $G$ 的元素，$\frac{-1+ \sqrt{3} \imath}{2} \times \frac{-1- \sqrt{3} \imath}{2}=1 \in G$，$\frac{-1+ \sqrt{3} \imath}{2} \times \frac{-1+ \sqrt{3} \imath}{2}=\frac{-1- \sqrt{3} \imath}{2} \in G$，$\frac{-1- \sqrt{3} \imath}{2} \times \frac{-1- \sqrt{3} \imath}{2}=\frac{-1+ \sqrt{3} \imath}{2} \in G$，所以满足封闭性；

&emsp;&emsp;b、结合律：由于此处的乘法为普通代数乘法，所以满足结合律；

&emsp;&emsp;c、单位元：$\forall a,b \in G,都有a \times 1=1 \times a=a$，所以满足单位元；

&emsp;&emsp;d、逆元素：$\frac{-1+ \sqrt{3} \imath}{2} \times  \frac{-1- \sqrt{3} \imath}{2}=1$，所以其互逆，而1的逆元就是1，所以满足逆元素；

&emsp;&emsp;e、交换律：普通代数乘法显然满足交换律。

&emsp;&emsp;综上所述，$G$ 为单位元为1的阿贝尔群，且群的阶为3，即 $|G|=3$。我们再来分析一下各元素的阶，不难发现1的阶是1，$\frac{-1+ \sqrt{3} \imath}{2}$ 的阶是3，$\frac{-1- \sqrt{3} \imath}{2}$ 的阶也为3。我们在这里给出群 $G$ 的乘法表：

|            $\times$             |                1                | $\frac{-1+ \sqrt{3} \imath}{2}$ | $\frac{-1- \sqrt{3} \imath}{2}$ |
|:-------------------------------:|:-------------------------------:|:-------------------------------:|:-------------------------------:|
|                1                |                1                | $\frac{-1+ \sqrt{3} \imath}{2}$ | $\frac{-1- \sqrt{3} \imath}{2}$ |
| $\frac{-1+ \sqrt{3} \imath}{2}$ | $\frac{-1+ \sqrt{3} \imath}{2}$ | $\frac{-1- \sqrt{3} \imath}{2}$ |                1                |
| $\frac{-1- \sqrt{3} \imath}{2}$ | $\frac{-1- \sqrt{3} \imath}{2}$ |                1                | $\frac{-1+ \sqrt{3} \imath}{2}$ |
不难发现其乘积结果部分组成矩阵的元素关于主对角线对称，矩阵所对应行列式的值为0，所以不可逆。我们接下来给出几个有限群的例子加以学习：这部分知识来自于[百度文库](https://wenku.baidu.com/view/56ef1e1d4b35eefdc8d33373.html)。

<div align=center><img src="figure/有限群.jpg" width="100%"/></div>

对于给定的群 $G$ ，对于任意的两个子集 $R,S \subseteq G$，我们令
$$
RS=\{r \cdot s\ |\ r \in R,\ s \in S\}
$$
特殊的，对于任意的 $g \in G$，如果 $R=\{g\}$，我们记作
$$
gS=\{g \cdot s\ |\ s \in S\}
$$
同样的，如果 $S=\{g\}$，我们记作
$$
Rg=\{r \cdot g\ |\ r \in R\}
$$
&emsp;&emsp;从现在开始，我们将乘法运算符进行省略，将 $g_1 \cdot g_2$ 写作 $g_1g_2$。

**定义1.3：** $G$ 为一个群，对于任意的 $g \in G$，我们令 $L_g$ 表示用 $g$ 左平移，具体计算方式为对于任意的 $a \in G$ 有 $L_g(a)=ga$。同样的，我们令 $R_g$ 表示用 $g$ 右平移，计算方式为对于任意的 $a \in G$ 有 $R_g(a)=ag$。

&emsp;&emsp;我们会经常用到下面这些简单的结论。

**命题1.3：** 给定群 $G$，其左平移 $L_g$ 和右平移 $R_g$ 得到的结果均满足双射。我们在这里仅给出左平移 $L_g$ 的证明过程，右平移 $R_g$ 的证明方式类似，证明方式如下：

&emsp;&emsp;若 $L_g(a)=L_g(b)$，那么有 $ga=gb$，我们在等式两边同乘 $g^{-1}$，便可得到 $a=b$，所以 $L_g$ 满足单射。对于任意的 $b \in G$，我们有 $L_g(g^{-1}b)=gg^{-1}b=b$，所以 $L_g$ 满足满射。因此，$L_g$ 满足双射。

&emsp;&emsp;我们在这里给出单射、满射以及双射的图解表示：
<div align=center><img src="figure/单射、满射、双射.jpg" width="80%"/></div>
<center>
图1.3 单射、满射、双射的示意图
</center>

&emsp;&emsp;特殊的，若映射过程的定义域和值域一样，我们将这个过程称之为变换，即有作用方式 $\varphi$，使得 $A \rightarrow A$。为了方便起见，我们在这里将映射后的结果记为 $\bar{A}$。

&emsp;&emsp;若我们定义 $A$ 上的二元运算为 $\cdot$，$\bar{A}$ 上的二元运算为 $\bar{\cdot}$ （此处的 $A$ 和 $\bar{A} 不一定相同$），对于 $A$ 中的元素 $a$ 经过作用方法 $\varphi$ 后得到 $\bar{a}$。那么，若我们定义 $A$ 上有两个元素进行二元运算 $a \cdot b$，经过函数映射后其计算结果为 $\varphi(a \cdot b)= \overline{a \cdot b}$，而对于每一个元素又都有 $\varphi(a)=\bar{a} \in \bar{A}$，$\varphi(b)=\bar{b} \in \bar{A}$，如果满足 $a \cdot b \rightarrow \bar{a} \bar{\cdot} \bar{b}$ （或者写成 $\overline{a \cdot b}=\bar{a} \bar{\cdot} \bar{b}$），称 $\varphi$ 是 $A$ 到 $\bar{A}$ 的同态映射。我们在这里给出一些例子方便大家理解：

<div align=center><img src="figure/同态.jpeg" width="60%"/></div>
<center>
图1.4 同态映射的判断方式
</center>
其中“象”指函数对变量的映射结果。

&emsp;&emsp;1、我们令 $A=R$，$\bar{A}=R^+$，其中 $A$ 上的二元运算为 $+$，$\bar{A}$ 上的二元运算为 $\times$，作用方式 $\varphi$ 为 $x \rightarrow e^x$，即求某一个变量的指数函数值。那么，$\forall x,y \in A,\overline{x+y}=e^{x+y}=e^{x}e^{y}=\bar{x}\bar{y}=\bar{x} \times \bar{y}$，故而为同态映射。这里的作用方式为双射。

&emsp;&emsp;2、我们有 $A=Z$，$A$ 上的二元运算为 $+$，$\bar{A}=\{1,-1\}$，$\bar{A}$ 上的二元运算为 $\times$，作用方式 $\varphi_1: \forall a \in A, a \rightarrow 1$。不难发现，$\forall a,b \in A$ 有 $a \cdot b = a+b \in A$，所以 $\overline{a+b}=\varphi_1(a+b)=1$。而对于 $\bar{a} \bar{\cdot} \bar{b}=\varphi_1(a) \times \varphi_1(b)=1 \times 1=1$，综上，$\overline{a+b}=\bar{a} \bar{\cdot} \bar{b}$，所以为同态映射。这里的作用方式为非满射，非单射。

&emsp;&emsp;3、在2中若把映射方式改为 $\varphi_2: \forall a \in A, a \rightarrow -1$，此时不难发现 $\overline{a+b}=\varphi_2(a+b)=-1$，$\bar{a} \bar{\cdot} \bar{b}=\varphi_2(a) \times \varphi_2(b)=(-1) \times (-1)=1$，故而不为同态映射。

&emsp;&emsp;4、在2中若把映射方式改为 $\varphi_3: \forall a为奇数 \in A, a \rightarrow -1，\forall a为偶数 \in A, a \rightarrow 1$。此时为同态映射。这里的作用方式为满射。

&emsp;&emsp;需要注意的是 $A$ 与 $\bar{A}$ 有同态映射且作用方式为满射我们才称 $A$ 与 $\bar{A}$ <font color=#ea4335>同态</font>，例如例4。当映射方式为双射时，此时我们称之为<font color=#ea4335>同构</font>，例如例1，我们将其记为 $A \cong \bar{A}$。我们先给出同态的一些性质：

&emsp;&emsp;假定对于代数运算 $\cdot$ 与 $\bar{\cdot}$ 来说，$A$ 与 $\bar{A}$ 同态，那么（1）若 $\cdot$ 满足结合律，则 $\bar{\cdot}$ 也满足结合律；（2）若 $\cdot$ 满足交换律，则 $\bar{\cdot}$ 也满足交换律。

&emsp;&emsp;我们再对同构进行分析，我们先给出一个例子：有 $A=\{1,2,3\}$，$\bar{A}=\{4,5,6\}$，$A$ 上的二元运算方式 $\cdot$ 如表1，$\bar{A}$ 上的二元运算方式 $\bar{\cdot}$ 如表2，$\varphi$ 是 $A \rightarrow \bar{A}$ 的映射，映射关系如表3：

<div align=center><img src="figure/同构.jpg" width="90%" /></div>

不难发现 $\varphi$ 是一个映射，且 $A$ 中的元素均有象且唯一。由于 $\varphi(a \cdot b)=\varphi(3)=6$，且 $\varphi(a) \bar{\cdot} \varphi(b)=6$，故而满足同态映射。$\varphi$ 满足双射。综合分析可知，$A \cong \bar{A}$，即满足同构。

&emsp;&emsp;那么，如果我们将 $\bar{\cdot}$ 的计算方式修改为 $\forall x,y \in \bar{A},x \bar{\cdot} y=5$，此时 $A$ 与 $\bar{A}$ 还满足同构吗？答案是肯定的，虽然此时的 $\varphi$ 不能满足同态映射，但是我们能找到 $\varphi'$ 使其满足同态映射，比如 $1 \rightarrow 4, 2 \rightarrow 6, 3 \rightarrow 5$，那么 $\varphi'$ 就可使 $A \cong \bar{A}$，所以还是满足同构的。综上，我们只要能找到一个（**存在**）映射方式（<u>这个映射方式不一定是给定的</u>），使得 $A \cong \bar{A}$，我们便说它们之间是同构的。这个存在的映射方式（同构映射）需要满足以下4点：

（1） $\varphi'$ 是 $A$ 到 $\bar{A}$ 的映射；

（2） $\varphi'$ 是同态映射；

（3） $\varphi'$ 是满射；

（4） $\varphi'$ 是单射。

&emsp;&emsp;特殊的，若 $\bar{A}=A$ 且 $A \cong \bar{A}$，那么我们将映射 $\varphi$ 称之为 $A$ 的自同构。举例如下：

<div align=center><img src="figure/自同构.jpg" width="90%" /></div>

**定义1.4：** 给定一个群 $G$ ，$G$ 的子集 $H$ 是其子群的充要条件是：

（1） $G$ 的幺元 $e$ 也是 $H$ 的元素（$e \in H$）;

（2） 对于所有的 $h_1,h_2 \in H$，都有 $h_1h_2 \in H$;

（3） 对于所有的 $h \in H$，都有 $h^{-1} \in H$。

&emsp;&emsp;定义的证明：我们先证明充分性 $\Leftarrow$，封闭性——由于 $\forall \  h_1,h_2 \in H,h_1h_2 \in H$，所以满足封闭性。结合律——乘法满足结合律。单位元——根据 $\forall \ h_1 \in H,h_1^{-1} \in H$，且 $\forall \ h_1,h_2 \in H,h_1h_2 \in H$，所以 $h_1h_1^{-1} \in H$，我们将 $h_1h_1^{-1}$ 记作 $e$，那么 $e \in H$。逆元素——由于 $\forall \ h,h^{-1} \in H,h^{-1}h=e$，所以逆元存在；我们再证明必要性 $\Rightarrow$，若 $H$ 为 $G$ 的子群，那么 $H$ 必须满足封闭性，所以能推出（2）。我们将群 $H$ 中的单位元记为 $e'$ （因为此处我们还不能确定子群的单位元是否等于群的单位元），$\forall \ a \in H,e'a=a$，同时由于 $H$ 是 $G$ 的子集，所以 $e',a \in G$，我们又知在 $G$ 中 $ea=a$，所以 $e'=e$，<u>即 $G$ 中的单位元就是 $H$ 中的单位元</u>。由于 $H$ 是一个群，那么 $\forall \ h \in H,hh^{-1}=e\ 且 h^{-1} \in H$，由于 $H$ 是 $G$ 的子集，所以 $h,h^{-1} \in G$，而在 $G$ 中显然有 $hh'=e$，所以 $h^{-1}=h'$，<u>即任意元素在群中的逆元素和子群中的逆元素相等</u>。

&emsp;&emsp;命题1.4的证明过程我们留作练习。

**命题1.4：** 给定一个群 $G$，其子集 $H \subseteq G$ 是群 $G$ 的子群 $\Leftrightarrow$ $H$ 非空且对于任意的 $h_1,h_2 \in H$，都有 $h_1h_2^{-1} \in H$。

&emsp;&emsp;命题的证明：我们先证明必要性 $\Rightarrow$，$\forall \ h_1,h_2 \in H$，由定义1.4的（3）知 $h_2^{-1} \in H$，再结合定义1.4的（2）知 $h_1h_2^{-1} \in H$；我们再证明充分性 $\Leftarrow$，$\forall \ h_1 \in H$，都有 $h_1h_1^{-1}=e \in H$，那么此时 $h_1h_1^{-1},h_1 \in H$，所以 $h_1h_1^{-1}h_1^{-1}=eh_1^{-1}=h_1^{-1} \in H$。由上一步知 $\forall \ h_1,h_2 \in H$，有 $h_1,h_2^{-1} \in H$，所以有 $h_1(h_2^{-1})^{-1}=h_1h_2 \in H$。

&emsp;&emsp;若群 $G$ 是有限群，那么可以使用下述判断方法。

**命题1.5：** 给定有限群 $G$，其子集 $H \subseteq G$ 是群 $G$ 的子群 $\Leftrightarrow$ （1） $e \in H$; （2） 两个 $H$ 做乘积运算后得到的结果是封闭的。

证明：我们先证明必要性 $\Rightarrow$，由于 $H$ 是 $G$ 的子群，所以其肯定满足封闭性且具有单位元 $e$；我们再证明充分性 $\Leftarrow$，封闭性——命题1.5的（2）就是封闭性的定义。结合律——乘积运算满足结合律。单位元——命题1.5的（1）就是单位元。逆元素—— $\forall \ a \in H$，根据封闭性可知 $a,a^1,a^2,\cdots,a^n \in H$，因为 $H$ 的元素个数有限，所以 $\exists \ i,j \in N \ (j>i) \ ,使得 \ a^i=a^j$，这里是为什么呢？若某个元素的 $n$ 次幂之间都不相等，那必然造成集合元素个数是无穷的，这与元素个数有限相悖，所以要使集合的元素个数有限，那元素的 $n$ 次幂的计算结果之间一定存在"周期性"。不难发现 $a^j=a^ia^{j-i}$，由于 $a^i=a^j$，那么 $a^{j-i}=1$，也即 $aa^{j-i-1}=1$。若 $j-i>1$，那么 $a^{-1}=a^{j-i-1} \in H$。若 $j-i=1$，那么 $a=1$，故而 $a^{-1}=a \in H$，因此 $\forall \ a \in H,都有  \ a^{-1} \in H,且 \ aa^{-1}=e \in H$。得证逆元素存在。综上所述，命题1.5成立。这部分知识来自于[百度文库](https://wenku.baidu.com/view/3e3b660a79563c1ec5da7184.html?sxts=1583120406111)。

**示例1.2：**

&emsp;&emsp;1. 对于任意的整数 $n \in Z$，集合 $nZ=\{nk\ |\ k \in Z\}$ 是群 $Z$ 的子群。

&emsp;&emsp;2. 对于 $n \times n$ 的可逆矩阵而言，若其满足 $GL^{+}(n,R)=\{A \in GL(n,R)\ |\ det(A)>0\}$，此时 $GL^{+}(n,R)$ 是群 $GL(n,R)$ 的子群。

&emsp;&emsp;3. 群 $SL(n,R)$ 是群 $GL(n,R)$ 的子群。

&emsp;&emsp;4. 群 $O(n)$ 是群 $GL(n,R)$ 的子群。

&emsp;&emsp;5. 群 $SO(n)$ 是群 $O(n)$ 的子群，同时也是群 $SL(n,R)$ 的子群。

&emsp;&emsp;6. 不难发现，每一个 $2 \times 2$ 的旋转矩阵 $R \in SO(2)$ 都可以被写作
$$
R=
\left(
 \begin{matrix}
   \cos \theta & -\sin \theta\\
   \sin \theta & \cos \theta
 \end{matrix}
\right),\ \ \ 其中\ 0 \leq \theta < 2\pi
$$
&emsp;&emsp;  我们在这里给出二阶矩阵的逆矩阵的计算方法
$$
\left(
 \begin{matrix}
   a & b\\
   c & d
 \end{matrix}
\right)^{-1}=\frac{1}{ad-bc}
\left(
 \begin{matrix}
   d & -b\\
   -c & a
 \end{matrix}
\right)
$$

&emsp;&emsp;  在下例中，$SO(2)$ 可以被看作是 $SO(3)$ 的子群
$$
R=
\left(
 \begin{matrix}
   \cos \theta & -\sin \theta\\
   \sin \theta & \cos \theta
 \end{matrix}
\right),\ \
Q=
\left(
 \begin{matrix}
   \cos \theta & -\sin \theta & 0\\
   \sin \theta & \cos \theta & 0\\
   0 & 0 & 1
 \end{matrix}
\right)
$$
&emsp;&emsp; 我们在这里给出旋转矩阵的定义，首先观察下图：
<div align=center><img src="figure/旋转矩阵.jpg" width = 400 height = 400 /></div>
<center>
图1.4 绕原点二维旋转
</center>

&emsp;&emsp;首先我们需要明确的是，二维旋转是围绕坐标原点旋转，如图1.4所示。图中点 $V$ 绕原点逆时针转过 $\theta$ 角到达 $V'$ 点处。假设点 $V$ 的坐标为 $(x,y)$，那么 $V'$ 点的坐标为 $(x',y')$，其中点 $V$ 到原点的距离为 $r$，且射线 $OV$ 与 $x$ 轴的夹角为 $\phi$。那么我们能得到如下等式：
$$
x = r \cos \phi, \ \ \ \ y = r \sin \phi\\
x' = r \cos (\theta + \phi), \ \ \ \ y' = r \sin (\theta + \phi)
$$
我们对 $x'$ 和 $y'$ 的等式进行展开可得：
$$
x' = r \cos \theta \cos \phi - r \sin \theta \sin \phi\\
y' = r \sin \theta \cos \phi + r \cos \theta \sin \phi
$$
我们将 $x$ 和 $y$ 的表达式代入上式可得：
$$
x' = x \cos \theta - y \sin \theta, \ \ \ \ y' = x \sin \theta + y \cos \theta
$$
我们将上述结果用矩阵进行表示可得：
$$
\left[
 \begin{matrix}
   x'\\
   y'
 \end{matrix}
\right]=
\left[
 \begin{matrix}
   \cos \theta & -\sin \theta\\
   \sin \theta & \cos \theta
 \end{matrix}
\right]
\left[
 \begin{matrix}
   x\\
   y
 \end{matrix}
\right]
$$
不难发现此处的系数矩阵便是我们的二阶旋转矩阵。那么我们尝试考虑一下三维空间内绕 $x$ 轴旋转的情况是什么样子的呢，我们先给出一个直观的感受。
<div align=center><img src="figure/三维旋转矩阵.jpg" width = 400 height = 400 /></div>
<center>
图1.5 绕x轴三维旋转
</center>

如图1.5所示，我们对 $OY$ 和 $OZ$ 绕 $x$ 轴旋转 $\theta$ 角度，分别到达 $OY'$ 以及 $OZ'$ 的位置，那么 $Y'$ 的坐标为 $(0,r \cos \theta,r \sin \theta)$，$Z'$ 的坐标为 $(0,-r \sin \theta,r \cos \theta)$，为了简单起见，我们取 $r=1$，那么 $Y$ 的坐标为 $(0,1,0)$，$Z$ 的坐标为 $(0,0,1)$，而 $X$ 的坐标始终不变。所以我们有如下等式：
$$
x' = x\\
y' = y \cos \theta - z \sin \theta \\
z' = y \sin \theta + z \cos \theta
$$
其中，$y',z'$ 的计算方式参见二维情况，我们将上述等式组利用矩阵进行表示如下：
$$
\left[
 \begin{matrix}
   x'\\
   y'\\
   z'
 \end{matrix}
\right]=
\left[
 \begin{matrix}
   1 & 0 & 0\\
   0 & \cos \theta & -\sin \theta\\
   0 & \sin \theta & \cos \theta
 \end{matrix}
\right]
\left[
 \begin{matrix}
   x\\
   y\\
   z
 \end{matrix}
\right]
$$
不难发现此形式和三阶旋转矩阵相同，只是示例6中的三阶旋转矩阵是绕 $z$ 轴旋转得到的。那么，旋转矩阵有什么性质呢？我们以二阶为例进行分析。我们将二阶旋转矩阵记为 $R$ ，不难发现其行列式计算结果为1，且矩阵的转置等于矩阵的逆，也即 $R$ 为正交阵。

&emsp;&emsp;7. 形如下式这种 $2 \times 2$ 的上三角矩阵
$$
\left(
 \begin{matrix}
   a & b\\
   0 & c
 \end{matrix}
\right)\ \ a,b,c \in R\ \ \ \ a,c \neq 0
$$
&emsp;&emsp; 是群 $GL(2,R)$ 的子群。

&emsp;&emsp;8. 集合 $V$ 由4个矩阵组成，这些矩阵的具体形式如下
$$
\left(
 \begin{matrix}
   \pm 1 & 0\\
   0 & \pm 1
 \end{matrix}
\right)
$$
集合 $V$ 是群 $GL(2,R)$ 的子群，被称为克莱因四元群。

**定义1.5：** 若 $H$ 为 $G$ 的一个子群，并且对于任意的 $g \in G$，形如 $gH$ 的计算方式称为 $H$ 在 $G$ 中的左陪集，形如 $Hg$ 的计算方式称为 $H$ 在 $G$ 中的右陪集。$H$ 的左陪集（右陪集亦同）中包含一种等价关系 $\sim$（$a \sim b \Leftrightarrow ab^{-1} \in H$，$a \sim b \Leftrightarrow b^{-1}a \in H$） ，定义如下：对于所有的 $g_1,g_2 \in G$ 有
$$
g_1 \sim g_2 \Leftrightarrow g_1H=g_2H, \ \ 同理\\
g_1 \sim g_2 \Leftrightarrow Hg_1=Hg_2
$$

证明：右陪集——我们将 $H$ 表示为 $H=\{h_1,h_2,h_3,\cdots\}$，那么 $Hg_2=\{h_1g_2,h_2g_2,h_3g_2,\cdots\}$。若 $g_1 \sim g_2$，则 $g_1g_2^{-1} \in H$，即 $\exists \ h \in H,g_1g_2^{-1}=h \Rightarrow g_1=hg_2 \in Hg_2$。$\forall \ hg_2 \in Hg_2,h \in H$，$hg_2g_2^{-1}=he=h \in H$，即 $hg_2 \sim g_2$。左陪集—— $H=\{h_1,h_2,h_3,\cdots\}$，那么 $g_2H=\{g_2h_1,g_2h_2,g_2h_3,\cdots\}$。若 $g_1 \sim g_2$，则 $g_2^{-1}g_1 \in H$，即 $\exists \ h \in H,g_2^{-1}g_1=h \Rightarrow g_1=g_2h \in g_2H$。$\forall \ g_2h \in g_2H, h \in H$，$g_2^{-1}g_2h=eh=h \in H$，即 $g_2h \sim g_2$。

显然，$\sim$ 是一种等价关系。我们这里先说明<font color=#ea4335>关系</font>的定义：对于 $R：A \times A \rightarrow D$，其中 $D=\{对，错\}$，$R$ 为Relation的首字母，若 $R(a,b)=对$，那么我们说 $(a,b)$ 满足关系 $R$，记为 $aRb$。例如：$A=\{1,2\}$，$R$ 表示 $>$，那么 $>(1,2)=\{错\}，>(2,1)=\{对\}$，所以 $(2,1)$ 满足关系 $>$，记为 $2>1$。那么什么是<font color=#ea4335>等价关系</font>呢？等价关系首先应该满足关系，除此之外还要满足：（1）反身性：$\forall a \in A,a \sim a$;（2）对称性：$\forall a,b \in A$，若 $a \sim b$，则 $b \sim a$;（3）传递性：$\forall a,b,c \in A$，若 $a \sim b,b \sim c$，则 $a \sim c$。例如：相等、三角形相似、三角形全等都是等价关系。我们给出一个例子加以说明：

&emsp;&emsp;例：我们有 $A=Z$，关系定义如下：当 $a=b \ (mod\ n)$ 时，$R(a,b) \rightarrow 对$，否则，$R(a,b) \rightarrow 错$，其中 $n$ 为正整数。我们来验证 $R$ 满足等价关系：

（1）：$R(a,a) \Rightarrow a\  mod \ n=a\  mod \ n$ 恒成立，所以满足反身性；

（2）：$aRb \Rightarrow a \ mod \ n=b \ mod \ n \Rightarrow b \ mod \ n=a \ mod \ n$
$\Rightarrow bRa$，所以满足对称性；

（3）：$aRb,bRc \Rightarrow a\ mod \ n=b \ mod \ n,\ b \ mod \ n=$
$c \ mod \ n \Rightarrow a \ mod \ n=c\ mod \ n$，所以满足传递性。

我们在这里将满足同余关系的所有元素可以归为一类，将其称为剩余类，例如余数为0的记作 $[0]$，余数为1的记作 $[1]$，依此类推，那么我们可以将整个集合 $A$ 划分为 $n$ 个互不相交的类，我们将该过程称为集合的分类。

&emsp;&emsp;有了这些知识以后呢，我们便可对陪集中的等价关系进行证明。自反性—— $H$ 是 $G$ 的子群，所以 $G$ 的单位元等于 $H$ 的单位元，元素在 $G$ 中的逆元等于其在 $H$ 中的逆元。故而有 $\forall \ a \in H,aa^{-1} \in H$；对称性——$ab^{-1} \in H \rightarrow (ab^{-1})^{-1} \in H \rightarrow ba^{-1}=(ab^{-1})^{-1} \in H$，这里可以根据 $ab^{-1}ba^{-1}=e$ 得出；传递性——$ab^{-1} \in H,bc^{-1} \in H \rightarrow ab^{-1}bc^{-1}=ac^{-1} \in H$。

&emsp;&emsp;现在，我们引入如下结论：

**命题1.6：** 给定一个群 $G$ 以及 $G$ 的任意一个子群 $H$，我们都有 $\forall \ g_1,g_2 \in G,\ g_1H=g_2H \Leftrightarrow g_2^{-1}g_1H=H \Leftrightarrow g_2^{-1}g_1 \in H$，证明过程如下：

&emsp;&emsp;若我们利用双射 $L_{g_2^{-1}}$ 同时对 $g_1H$ 和 $g_2H$ 进行作用，我们能分别得到 $L_{g_2^{-1}}(g_1H)=g_2^{-1}g_1H$，$L_{g_2^{-1}}(g_2H)=g_2^{-1}g_2H=H$，所以 $g_1H=g_2H \Leftrightarrow g_2^{-1}g_1H=H$，若 $g_2^{-1}g_1H=H$，由于 $1 \in H$，所以我们根据封闭性能得到 $g_2^{-1}g_1 \in H$。相反的，如果 $g_2^{-1}g_1 \in H$，由于 $H$ 是一个群，所以其左平移 $L_{g_2^{-1}g_1}$ 是 $H$ 的双射，所以 $g_2^{-1}g_1H=H$。综上，$g_2^{-1}g_1H=H \Leftrightarrow g_2^{-1}g_1 \in H$。

&emsp;&emsp;因此元素 $g \in G$ 的等价类是陪集 $gH$（或陪集 $Hg$）。子群不同的左（右）陪集就是将整个群划分为若干互不相交的集合。由于 $L_g$ 是 $H$ 到 $gH$ 的双射，所以 $H$ 的所有陪集 $gH$ 具有相同的基数（元素的个数）且与 $H$ 的元素个数相等。映射 $L_{g^{-1}} \circ R_g$ 是左陪集 $gH$ 与右陪集 $Hg$ 之间的映射，由于 $\forall \ h \in H,L_{g^{-1}}(R_g(gh))=L_{g^{-1}}(ghg)=g^{-1}ghg=hg$，所以映射为双射，故而它们有相同的基数。综上，左（右）陪集元素个数与 $H$ 的元素个数相等。我们从特定的群 $G$ 中取特定的陪集 $gH$，可以得到如下结论：

**命题1.7：** （拉格朗日定理）对于任意的有限群 $G$ 以及其任意子群 $H$，$H$ 的阶数可以除尽 $G$ 的阶数。

**定义1.6：** 给定一个有限群 $G$ 以及 $G$ 的一个子群 $H$，如果 $n=|G|,\ h=|H|$，那么比值 $\frac{n}{h}$ 可以记为 $(G:H)$，并将其称为 $H$ 在 $G$ 中的指数。

&emsp;&emsp;指数 $(G:H)$ 是 $G$ 的子群 $H$ 中左陪集（或右陪集）的个数，命题1.7可以被表述为 $|G|=(G:H)|H|$。$H$ 中所有左陪集组成的集合（一般来说，这个集合不是一个群）可以记作 $G/H$，一个陪集便是 $G/H$ 中的一个元素。我们这里给出一个有限群的例子加以说明：

&emsp;&emsp;我们令 $G$ 表示6阶群，$H$ 为 $G$ 的子群且 $H$ 的阶数为2
<div align=center><img src="figure/陪集.jpg" width="70%" /></div>

&emsp;&emsp;我们先分析 $H$ 的左陪集，$\forall \ g \in G$，有
<div align=center><img src="figure/左陪集.jpg" width="100%" /></div>

不难发现 $e \sim a$，$b \sim f$，$c \sim d$。所以 $H$ 的左陪集有3个，分别为 $(e,a)$，$(b,f)$，$(c,d)$，每个左陪集均和 $H$ 的元素个数相同，且左陪集个数等于群 $G$ 的阶除以 $H$ 的阶。$G/H=((e,a),(b,f),(c,d))$。

&emsp;&emsp;我们再分析 $H$ 的右陪集，$\forall \ g \in G$，有
<div align=center><img src="figure/右陪集.jpg" width="100%" /></div>

不难发现 $e \sim a$，$b \sim d$，$c \sim f$。所以 $H$ 的右陪集有3个，分别为 $(e,a)$，$(b,d)$，$(c,f)$，每个右陪集均和 $H$ 的元素个数相同，且右陪集个数等于群 $G$ 的阶除以 $H$ 的阶。$G/H=((e,a),(b,d),(c,f))$。

&emsp;&emsp;需要注意的是虽然左陪集的元素个数等于右陪集的元素个数，但是左陪集和右陪集不一定相同哦。

**示例1.3：**

&emsp;&emsp;1. $n$ 是任意的一个正整数，$Z$ （此处的二元运算为 $+$，即 $Z+Z \rightarrow Z$）的子群 $nZ$。0的陪集是集合 $\{0\}$，并且任何非零整数 $m \in Z$ 的陪集是：
$$
m+nZ=\{m+nk|k \in Z\}
$$
&emsp;&emsp;通过用 $n$ 除 $m$，对于 $0 \leq r \leq n-1$ 存在 $r$ 满足 $m=nq+r$。同时我们会发现 $r$ 是陪集 $m+nZ$ 的最小正元素。这意味着 $Z$ 的子群 $nZ$ 的陪集与以 $n$ 为模的余数组成的集合 $\{0,1,\cdots,n-1\}$ 之间满足双射关系，或者说 $Z$ 的子群 $nZ$ 的陪集与 $Z/nZ$ 之间满足双射关系。

&emsp;&emsp;2. $SL(n,R)$ 是 $GL(n,R)$ 的子群，$SL(n,R)$ 的陪集是矩阵组成的集合：
$$
A \ SL(n,R)=\{AB \ | \ B \in SL(n,R)\}, \ A \in GL(n,R)
$$
&emsp;&emsp;由于 $A$ 可逆，所以其行列式不为零 $det(A) \neq 0$，并且当 $det(A)>0$ 时我们可以将矩阵 $A$ 写作 $A=(det(A))^{1/n}((det(A))^{-1/n}A)$，当 $det(A)<0$ 时我们可以将矩阵 $A$ 写作 $A=(-det(A))^{1/n}((-det(A))^{-1/n}A)$。但是，当 $det(A)>0$ 时 $(det(A))^{-1/n}A \in SL(n,R)$；当 $det(A)<0$ 时，$-(-det(A))^{-1/n}A \in SL(n,R)$，所以陪集 $A \ SL(n,R)$ 包含矩阵：
$$
若  \ det(A)>0，(det(A))^{1/n}I_n\\
若  \ det(A)<0，-(-det(A))^{1/n}I_n
$$
且 $SL(n,R)$ 的陪集与 $R$ 之间具有双射关系。

&emsp;&emsp;3. $SO(n)$ 是 $GL^{+}(n,R)$ 的子群，$SO(n)$ 的陪集是矩阵组成的集合：
$$
A \ SO(n)=\{AQ \ | \ Q \in SO(n)\}，A \in GL^{+}(n,R)
$$
我们若将矩阵用极坐标进行表示，可以发现 $SO(n)$ 的陪集和 $n$ 阶正的对称的正定矩阵之间满足双射关系，且这些对称矩阵的特征值均为正。

&emsp;&emsp;4. 群 $SO(2)$ 是群 $SO(3)$ 的子群，$SO(2)$ 的陪集是由矩阵组成的集合：
$$
Q \ SO(2)=\{QR \ | \ R \in SO(2)\}，Q \in SO(3)
$$
群 $SO(3)$ 的作用我们可以表述为将球体 $R^3$ 表面 $S^2$ 上的一点 $x \in S^2$ 进行旋转，旋转过程满足 $\forall  \ Q \in SO(3), \ x \rightarrow Qx$，其中 $S^2=\{(x,y,z) \in R^3 \ | \ x^2+y^2+z^2=1\}$。我们用 $N=(0,0,1)$ 表示球体 $S^2$ 的北极点，不难发现，当 $N$ 固定时，$SO(2)$ 恰好是 $SO(3)$ 的子群。这导致陪集 $Q \ SO(2)$ 的所有的旋转 $Q \ R$ 都将点 $N$ 映射到点 $Q \ N \in S^2$ ，并且可以发现 $SO(2)$ 的陪集与 $S^2$ 上的点满足双射关系。此映射关系满足满射与 $SO(3)$ 对 $S^2$ 的作用满足传递性有关。$\forall \ x \in S^2$，均有旋转矩阵 $Q \in SO(3)$，使得 $QN=x$。

&emsp;&emsp;通过下式我们可以定义左陪集（或者右陪集）之间的乘积运算：
$$
(g_1H)(g_2H)=(g_1g_2)H
$$
但是这个运算性质并不是普遍存在的，除非子群 $H$ 拥有特殊的性质才成立。在示例1.3中，可以对例1进行陪集间的乘积运算，但是对例2，3中的陪集并不能做乘积运算。那么哪些陪集可以做乘积运算呢？若子群 $H$ 的核满足同态，那么其左陪集可以进行乘积运算，我们在这里给出同态的定义。

**定义1.7：** 给定两个任意的群 $G$ 和 $G'$，我们定义映射关系 $\varphi$，使得 $G \rightarrow G'$ 是同态 $\Leftrightarrow$ $\forall \ g_1,g_2 \in G, \ \varphi(g_1g_2)=\varphi(g_1)\varphi(g_2)$。当 $g_1=g_2=e \in G$ 时，我们可以发现 $\varphi(e)=e'$，当 $g_1=g,g_2=g^{-1}$ 时，我们可以发现 $\varphi(g^{-1})=(\varphi(g))^{-1}$。

**示例1.4：**

&emsp;&emsp;1. 映射 $\varphi:Z \rightarrow Z/nZ$，其中 $\forall \ m \in Z, \ \varphi(m)=m \ mod \ n$ 是一个同态。

&emsp;&emsp;2. 映射 $det:GL(n,R) \rightarrow R$ 是一个同态，因为 $\forall \ A,B \in GL(n,R),det(AB)=det(A)det(B)$。同样的映射 $det:O(n) \rightarrow R$ 也是一个同态。

&emsp;&emsp;如果 $\varphi:G \rightarrow G'$ 和 $\psi:G' \rightarrow G''$ 是群同态，那么 $\psi \circ \varphi:G \rightarrow G''$ 也是一个同态。如果 $\varphi:G \rightarrow G'$ 是群的同态，并且 $H \subseteq G,H' \subseteq G'$ 是两个子群，那么有 $Im \ H=\varphi(H)=\{\varphi(g) \ | \ g \in H\}$ 是 $G'$ 的子群，并且 $\varphi^{-1}(H')=\{g \in G \ | \ \varphi(g) \in H'\}$ 是 $G$ 的子群。特殊的，当 $H'=\{e'\}$ 时，计算结果为 $\varphi$ 的核，记作 $Ker \ \varphi$。

**定义1.8：** 若映射 $\varphi:G \rightarrow G'$ 是群的同态，且 $H \subseteq G$ 是 $G$ 的一个子群，那么 $G'$ 的子群可由下式计算 $Im \ H=\varphi(H)=\{\varphi(g) \ | \  g \in H\}$，其中 $Im \ H$ 称作 $H$ 在 $\varphi$ 之下的象，同时也是 $G'$ 的子群。$Ker \ \varphi=\{g \in G \ | \ \varphi(g)=e'\}$ 称作 $\varphi$ 的核。

**示例1.5：**

&emsp;&emsp;1. 同态 $\varphi:Z \rightarrow Z/nZ$ 的核是 $nZ$。

&emsp;&emsp;2. 同态 $det:GL(n,R) \rightarrow R$ 的核是 $SL(n,R)$，同样的，同态 $det:O(n) \rightarrow R$ 的核是 $SO(n)$。

&emsp;&emsp;我们将经常对满足单射的群的同态进行分析。

**命题1.8：** 如果映射 $\varphi:G \rightarrow G'$ 是群的同态，那么 $\varphi:G \rightarrow G'$ 是一个单射  $\Leftrightarrow$ $Ker \ \varphi=\{e\}$ （我们也可以写作 $Ker \ \varphi=(0)$）。

&emsp;&emsp;证明：假设 $\varphi$ 满足单射。若 $\varphi(e)=e'$， $\varphi(g)=e'$，则 $\varphi(g)=\varphi(e)$。因为 $\varphi$ 满足单射，所以 $g=e$，所以 $Ker \ \varphi=\{e\}$。相反的，假定 $Ker \ \varphi=\{e\}$，若 $\varphi(g_1)=\varphi(g_2)$，那么等式两端同乘 $(\varphi(g_1))^{-1}$，我们有 $e'=(\varphi(g_1))^{-1}\varphi(g_1)=(\varphi(g_1))^{-1}\varphi(g_2)$，由于 $\varphi$ 是一个同态，所以 $(\varphi(g_1))^{-1}=\varphi(g_1^{-1})$，所以
$$
e'=(\varphi(g_1))^{-1}\varphi(g_2)=\varphi(g_1^{-1})\varphi(g_2)=\varphi(g_1^{-1}g_2)
$$
上式表明 $g_1^{-1}g_2 \in Ker \ \varphi$，但由于 $Ker \ \varphi=\{e\}$ 所以我们有 $g_1^{-1}g_2=e$，故而 $g_2=g_1$，验证了 $\varphi$ 满足单射。

**定义1.9：** 若存在同态 $\psi:G' \rightarrow G$，我们说群的同态 $\varphi:G \rightarrow G'$ 是一个同构。也就是说
$$
\psi \circ \varphi=id_G \ \ 并且  \ \ \varphi \circ \psi=id_{G'}\ \ \ \ \ \ \ \ \  (\sharp)
$$
如果 $\varphi$ 是同构我们说群 $G$ 与 $G'$ 是同构的。当 $G'=G$ ，我们便将其称之为自同构。

&emsp;&emsp;命题1.2的证明过程显示了如果一个群同态 $\varphi:G \rightarrow G'$ 是一个同构，那么存在唯一的同态 $\psi:G' \rightarrow G$ 满足条件 $\sharp$。这个同态被记为 $\varphi^{-1}$。

&emsp;&emsp;左平移 $L_g$ 和右平移 $R_g$ 是 $G$ 的自同构。

&emsp;&emsp;假设 $\varphi:G \rightarrow G'$ 是一个双射同态，并且 $\varphi^{-1}$ 是 $\varphi$ 的逆映射，那我们对于 $\forall \ a,b \in G$，都有
$$
\varphi(\varphi^{-1}(a)\varphi^{-1}(b))=\varphi(\varphi^{-1}(a))\varphi(\varphi^{-1}(b))=ab\\
\varphi^{-1}(ab)=\varphi^{-1}(a)\varphi^{-1}(b)
$$
这证明了 $\varphi^{-1}$ 是一个同态，因此我们给出下述结论。

**命题1.9：** 双射群的同态 $\varphi:G \rightarrow G'$ 是一个同构。

&emsp;&emsp;我们先对性质 $(\ast)\ \ \ \forall \ g \in G,\ gH=Hg$ 进行分析，在等式两端同时乘 $g^{-1}$，便有 $\forall \ g \in G,\ gHg^{-1}=H$，并且呢 $(\ast\ast)\ \ \ \forall \ g \in G,\ gHg^{-1} \subseteq G$。这是因为 $\forall \ g \in G,\ gHg^{-1} \subseteq H$ 意味着 $H \subseteq g^{-1}Hg$。

**命题1.10：** 令 $\varphi:G \rightarrow G'$ 是一个群同态，那么 $H=Ker \ \varphi$ 满足性质($\ast$)和($\ast\ast$)。

&emsp;&emsp;证明：我们有
$$
\varphi(ghg^{-1})=\varphi(g)\varphi(h)\varphi(g^{-1})=\varphi(g)e'\varphi(g)^{-1}=\varphi(g)\varphi(g)^{-1}=e'
$$
$\forall \ h \in H=Ker \ \varphi$ 以及 $\forall \ g \in G$。所以，通过定义 $H=Ker \ \varphi$，我们有 $gHg^{-1} \subseteq H$。
