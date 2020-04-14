#!/usr/bin/expect -f

set timeout -1

spawn ./install_9_1_1.lisp

expect "q - Abort Installation.\r"

send -- "9\r"

expect "export MODULEPATH=*\r"

send -- "no\r"

expect eof 
