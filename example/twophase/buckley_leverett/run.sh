 ###
 # @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 # @Date: 2024-05-15 11:19:06
 # @LastEditTime: 2024-08-20 14:27:56
 # @FilePath: /cpgrid/example/twophase/buckley_leverett/run.sh
 # @Description: 
 # 
### 
mpirun -np 1 ./twophase -options_file optionsfile

# run: sh run.sh > log.txt
#      sh run.sh > .log 
#-options_file optionsfile