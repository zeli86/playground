#!/usr/bin/expect -f

set timeout -1

spawn ./install_9.1.1.lisp

expect "q - Abort Installation.\r"

send -- "a\r"

expect eof 
