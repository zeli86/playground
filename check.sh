#!/bin/bash
cppcheck -Iinclude/ -Isource/oct/ -Isource/realtime/ -Isource/groundstates/ -Isource/stationary/ -Isource/lib/ -j 6 . 2> err.txt
