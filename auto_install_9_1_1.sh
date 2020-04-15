#!/usr/bin/expect -f

set timeout -1

spawn ./install_9_1_1.lisp --build-dir=/home/zeli/local/src/ --install-dir=/home/zeli/local/opt/ --module-dir=/home/zeli/local/modules/modulefiles/ --debug -j 2

expect " - Abort Installation.\r"

send -- "2\r"

expect "*(yes or no)*"

send -- "no\r"

expect eof 
