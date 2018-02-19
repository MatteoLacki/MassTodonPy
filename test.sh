#!/bin/bash
for i in `seq 1 100`;
do
    echo $i
    masstodon_example_call;
done
