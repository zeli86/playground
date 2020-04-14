#!/usr/bin/expect -f

set timeout -1

spawn ./install_9_1_1.lisp

expect " - Abort Installation.\r"

send -- "a\r"

expect "export MODULEPATH=*\r"

send -- "yes\r"

expect eof 
