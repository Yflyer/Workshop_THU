# Linux, & QIIME2
Yufei Zeng
yfzeng0827@hotmail.com
github.com/Yflyer/
THU， 04/15/2022

-----


## 1. Install R and Rstudio
###### 1.1 Information for local use
所有的分析都将在Rstudio内进行，请大家现在自己PC上安装对应的R （R4版本） 和 Rstudio 
* R的安装与使用的地址都最好不要含有中文
### *About Rstudio*
RStudio是一款R语言具有调试、可视化等功能的IDE，解决了R自带的环境操作不便这个问题

###### 1.2 

Please,  
(1) connect to THU vpn in advance  
(2) install a remote terminal (e.g., Xshell)

```
IP：166.111.42.42
User: test0
password: ********
```
*Why we use Linux?*
* High stability
* Open source
* Perfect For Programing



## 2. Install Anaconda or Miniconda  

\# optional: show current path
<br>`pwd`<br>

\# download conda<br>
`wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh`
<br> or <br>
`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`  
conda有python2和python3的版本，推荐python3版本的  
*wget命令添加 --O 选项，可以更改文件存放路径及重命名，如wget --O /mnt/e/Anaconda3。sh

\# install conda  
`bash Anaconda3-2019.10-Linux-x86_64.sh`  
or  
`bash Miniconda3-latest-Linux-x86_64.sh`  
* 如果conda.sh的文件是直接拷贝的，上述命令文件名加上路径信息即可，如/mnt/e/Anaconda3-2019.10-Linux-x86_64.sh  
* wsl用户可能需要添加sudo前置，sudo即为管理员权限

\# optional: chown the conda file  
wsl用户sudo后可能需要conda文件夹权限，根据debug提示输入:  
`sudo chown $UID:$GID ~/.conda`  

\# restart the terminal to use conda

\# common conda mangement command
```
conda info #显示conda版本，env path有关的信息
conda update conda #conda 升级  
rm -rf anaconda  # conda卸载（在安装目录下）  
conda env list # 列出所有conda环境  
```  