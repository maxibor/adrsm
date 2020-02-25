#!/bin/bash

mkdir -p $PREFIX/bin
cp adrsm $PREFIX/bin/
mkdir -p $PREFIX/bin/lib
cp -r lib/* $PREFIX/bin/lib/
mkdir -p $PREFIX/data/quality
cp data/quality/*.p $PREFIX/data/quality/


