#!/usr/bin/env python3
'''
SHAKE128 executable
'''

import sys
from keccak import shake128_hex


def main():
    '''
    Main module.
    '''
    res_hex = shake128_hex(MSG, BYTES_NUM)
    print(res_hex)


if __name__ == "__main__":
    try:
        BYTES_NUM = int(sys.argv[1])
    except:
        BYTES_NUM = 32
    try:
        MSG = sys.stdin.buffer.read()
    except:
        MSG = b''
    main()
