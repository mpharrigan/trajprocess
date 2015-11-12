#!/bin/bash


if [ ! -e mock3-reference.tar.bz2 ]
then
    wget http://web.stanford.edu/~harrigan/mock3-reference.tar.bz2
    tar -xf mock3-reference.tar.bz2
fi

if [ ! -e mock4-reference.tar.bz2 ]
then
    wget http://web.stanford.edu/~harrigan/mock4-reference.tar.bz2
    tar -xf mock4-reference.tar.bz2
fi
