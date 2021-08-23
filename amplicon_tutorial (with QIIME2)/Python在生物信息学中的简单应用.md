# Python在生物信息学中的简单应用

Yufei Zeng, zengyf93@qq.com  
THU， 10/05/2020


-----
学习前的准备
1. 选择合适的工具  
IDE——Integrated development environment
* Pycharm
* Atom + wsl
* Jupyter notebook
2. 良好的编程习惯
* 注释清晰明白
* 学会版本管理 (git的使用)
* 多记笔记
* 尽可能优化流程
3. 学会看bug
* bug的第一行通常是问题所在
* 注意编码、缩进、数据类型等基础bug
4. 如何面对未知的新任务？  
(任务导向是学习新知识的最好帮手)
* Google (stackoverflow)
* Bing
* 百度

### 1. Python的数据类型

###### 1.1 数值（init）字符串（string）
 数值与字符串是生信操作中最常面对数据类型
1.1.1 变量赋值(var)  

数值型(init)  
```
a = 1
b = '30'
type(a)
type(b)
```
字符串型(string)
```
a = '1'
a = 'AATCGG'
```
1.1.2 运算符
数值型运算
```
3+2
4/2
5%2

5/2
5*1.0/2
10*1.0/3

3**2
```

字符型运算
```
# join the string
  ‘ATGC’+'NNNT'

# multipy the string
  'ATGC'*3

# Now we operate the string var
  a, b = 'ATCTCC','TCTAAA'
  c = a + b

# slice
  a[1],a[1:3],a[-1]
  c[1:]
  c[1::3]

# reverse
  c[::-1]

# index
  a.index('C')

# replace (or delete)
  a.replace('ATC','NNN')
  a.replace('ATC','')

# string length
  len(a)
  len(c)

# boolean
 a in c  
 'AT' in c  
 ```



###### 1.2 列表（list）与字典(dict)
 列表与字典都是python中常用的数据存储结构（快速小型的数据集）
1.2.1 列表
```
# make a empty list
 a = []

# set indivdual element of list  
 a[0] = 2
 a[2] = 2


 a = [1,2,3]
 b = ['4']

# they are different
 a+b
 a.append(b) #  列表嵌套

#  slice
 c = a+b
 c[1]
 c[1:3]
 c[-1]
 c[1:]

# index
 c.index('2')
 c.Index(2)

# Delete
del c[3]

# list length
 len(c)

# boolean
 a in c
 '2' in c
```
 for循环
for i in range:
    do something
for i in range：
    for n in range:
        do something



1.2.2 字典
字典
```
# make a dict
  a = {}
  a['key'] = 'value'

# make a dict of standard gentic
 a['GCT']='Ala'
 a['TGT']='Cys'
 a['Glu']='Cys'
 print(a)

# get information from dict  
 a.items()
 a.keys()
 a.values()
```

### 2. Python的常见操作
###### 2.1 如何打开，读取及写入文件
2.1.1 基础命令  
```
# open file
  f = open('path','r','buffer')
```
||mode的选项|
|-|-|
|w|写入|
|a|追加|
|r|读取|
* 如果文件不存在，open()函数就会抛出一个IOError的错误，并且给出错误码和详细的信息告诉你文件不存在
```
# 一次性读取
f.read()

# 一次性读取（分行）
f.readlines()

# 逐行读取
f.readline()
```
文件使用完毕后必须关闭，因为文件对象会占用操作系统的资源
```
f.close()
```
但是代码调试时往往运行过程中就会报错，这时候with语句可以保证文件的关闭
```
with open('/path', 'r') as f:
    print(f.read())
```
如果需要写入内容到新文件：
```
f = open('try.txt', 'w')
f.write('Hello, world!')
f.close()
```
* 如果文件已存在，会直接覆盖（相当于删掉后新写入一个文件）。如果希望追加到文件末尾,可以mode选择'a'（append）以追加

###### 2.2 自动化操作的基础：循环与判断语句、异常处理
2.2.1 条件判断
```
# 增加elif或者else
Cys = ['TGT','TGC']
seq = 'TAA'
if seq in Cys:
    print('The seq included cysteine')
elif seq == 'TAA' or 'TAG' or 'TGA':
    print('''It's a stop codon.''')
else:
    print('No cycteine found')

```
* 注意if判断是从上往下进行的，设置多个条件的时候互相之间不要嵌套  

2.2.2 for循环  
* for循环的对象: string, list, enumerator，**file**.
先介绍range
```
range(0,5)
range(1,10,3)

# python3无法直接表达range序列，需要list一下
list()
```
操作一下
```
seq = 'TTTAAAGCG'
for i in range(0,len(seq),3):
    print(seq[i:i+3])
for i in seq:
    print(i)
```

2.2.3 while循环
```
sum = 0
n = 99
while n > 0:
    sum = sum + n
    n = n - 2
print(sum)
```
2.2.4 其他控制语句
用break终止循环,当n = 11时，条件满足，执行break语句
```
n = 1
while n <= 100:
    if n > 10:
        break
    print(n)
    n = n + 1
```
next可以在循环中多执行一次循环；
用contiune跳过某些循环。
下面这个例子结合二者实现某些特定的功能
```
f = open('path','r')
for line in f:
    if line[0] != '>':
        continue
    print(line)
    print(next(f))

f.close()
```
* 尽可能地简化你的代码，不要添加过多的控制语句。  

2.2.5 异常处理
```
try:
    fh = open("testfile", "w")
    fh.write("这是一个测试文件，用于测试异常!!")
except IOError:
    print "Error: 没有找到文件或读取文件失败"
else:
    print "内容写入文件成功"
    fh.close()

# 还可以使用try-finally结构，即使报错也继续往下执行
```


###### 2.3 简单脚本演示：
2.3.1 样品名中出现无法识别的特殊字符，如何批量重命名样品名？  
2.3.2 如何快速生成数据的mapping文件？  
2.3.3 如何根据样品名将总样本按照样本名分开成单独样本？  
2.3.4 如何得到一段序列的碱基互补序列？  

### 3. 课堂作业
如何将fastq文件转化成为fasta文件？  
如何将环状基因组拆分成1K长度为单位的序列集，并统计每段长度上的GC含量？
