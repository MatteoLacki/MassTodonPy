#!/bin/bash
for i in `seq 1 10`;
do
    echo $i
    masstodon_example_call;
done
