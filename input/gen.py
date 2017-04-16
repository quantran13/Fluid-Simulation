#!/usr/bin/python3

import sys

def printUsage(progName):
    print("Usage: {} <n> <output file>\n".format(progName))
    exit(0)

def check(x, n):
    return x >= n/4 and x <= 3*n/4

def main():
    if len(sys.argv) != 3:
        printUsage(sys.argv[0])

    n = int(sys.argv[1])
    fileName = sys.argv[2]
    fout = open(fileName, "w")

    fout.write(str(n)+"\n")
    fout.write(str(1)+"\n")
    fout.write(str(1)+"\n")

    for i in range(n):
        for j in range(n):
            for k in range(n):
                if check(i, n) and check(j, n) and check(k, n):
                    print(1)
                else:
                    print(0)

    fout.close()

main()

